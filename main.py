"""逐网架主线：1 初始化 → 2 首轮射线 → 3 最大违反量切割 → 4 边界补充 → 5 全局查漏。"""
import argparse
from pathlib import Path
from time import perf_counter

import numpy as np
from gurobipy import GRB
from threadpoolctl import threadpool_limits

from Network.four_bus_five_corridor import FourBus
from Network.case33bw import Case33
from model import GridPhysics, MasterProblem, SubProblem, RemainingRegionModel, PLANNING_TOL
from monitor import RunMonitor, RegionTimeout, SynchronizedReplay, _union
from region import RegionState, GEOMETRY_TOL, contains, polytope_volume, initial_polytope
from plot import union_volume
from vertify import validate_socp_region

# 运行设置：DIMENSION 同时控制计算、扫描及原生回放的维数。
NETWORK = Case33                  # FourBus / Case33
DIMENSION = 2                     # 2 / 3；Case33: (18,25) / (18,25,30)
BUDGET = 20000.                   # FourBus 建设预算；Case33 用其 switch_budget
CASE_TIME_LIMIT = 50.             # 构域总时限，含初始化和记录，不含事后扫描
SOLVER_THREADS = 20                # 构域及单进程扫描中，每个求解器的线程数
REGION_TAU = .005                 # 全局覆盖的径向精度；不是可行性容差
CUT_THRESHOLD = .02               # 连续小割的面积/体积比例
CUT_PATIENCE = 3
RAY_THRESHOLD = 1e-2              # 射线收益上界不足本网架认证面积/体积的 1% 时跳过；首轮不筛选
POINT_TOL = 1e-2                  # 同一几何点的最大坐标差，kW
DIVISIONS = 80                    # FourBus 每轴扫描格数
SCAN_DIVISIONS = {2: 100, 3: 60}   # Case33 每轴扫描格数，按维数选择
SCAN_WORKERS = 20                 # 扫描进程数；并行时每个求解器用 1 个线程
FORCE_RESCAN = False              # True 强制重扫；False 在已有功率范围内复用
SHOW_UI = True
ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT/'results'/'mainline'
SCAN_OUTPUT = ROOT/'results'/'scans'


def recording_path(case, output=OUTPUT, load_nodes=None, dimension=DIMENSION):
    nodes = ({'fourbus': (1, 2, 3), 'case33': (18, 25, 30)}[case][:dimension]
             if load_nodes is None else load_nodes)
    return Path(output)/f"{case}_{'_'.join(map(str, nodes))}.json.gz"


def union_measure(polytopes, d):
    """二维面积 / 三维体积；重叠网架只计一次。"""
    return _union(polytopes).area if d == 2 else union_volume(polytopes)


def register_power(power, powers, point_tol):
    """以 kW 最大坐标差匹配固定代表；不舍入坐标，也不沿近点链移动代表。"""
    for index, representative in enumerate(powers):
        if np.max(np.abs(power-representative)) <= point_tol:
            return index
    powers.append(np.asarray(power, dtype=float).copy())
    return len(powers)-1


def stage_candidates(region, x):
    """只检查当前 N_x 的顶点；其他网架的认证边界不参与局部选点。"""
    outer = region.records[tuple(x)]['outer']
    return outer[~contains(outer, region.inner_equations(x))]


def ray_gain(inner, point):
    """射线可行端点在线段上，新增测度不超过把目标点直接收入凸包的增量。"""
    return polytope_volume(np.vstack([inner, point]))-polytope_volume(inner)


def coverage_halfspaces(region):
    """从大到小保留认证域；被单个已保留凸包包含的排除项无需重复建模。"""
    keys = [key for key, row in region.records.items() if len(row['inner'])]
    keys.sort(key=lambda key: polytope_volume(region.records[key]['inner']), reverse=True)
    kept, equations = [], []
    for key in keys:
        if any(contains(region.records[key]['inner'], eq).all() for eq in equations):
            continue
        kept.append(key)
        equations.append(region.inner_equations(key))
    return kept, equations


def ray_support(equations, budget, x, anchor, power, *, threads, time_limit, numeric_focus=0):
    """固定 x，从同网架认证锚点朝不可行顶点求最远可行点。"""
    deadline = perf_counter()+time_limit
    problem = MasterProblem(equations, budget=budget,
        fixed_plan=equations.network.decode_plan(x), threads=threads)
    with problem.model as model:
        ray_fraction = model.addVar(ub=1., name='ray_fraction')
        for i, a, q in zip(equations.network.load_nodes, anchor, power):
            model.addConstr(problem.loads[i]/equations.network.base == float(a/equations.network.base)
                +ray_fraction*float((q-a)/equations.network.base), name=f'ray_power[{i}]')
        model.setObjective(ray_fraction, GRB.MAXIMIZE)
        model.Params.Aggregate = 0
        model.Params.ScaleFlag = 0
        model.Params.NumericFocus = numeric_focus
        model.Params.BarQCPConvTol = 1e-9
        model.Params.TimeLimit = max(0., deadline-perf_counter())
        model.optimize()
        if model.Status == GRB.TIME_LIMIT:
            raise TimeoutError('射线支撑达到总时限')
        if model.Status != GRB.OPTIMAL:
            raise RuntimeError(f'Ray support: status={model.Status}')
        if model.MaxVio > PLANNING_TOL:
            raise RuntimeError(f'Ray support: MaxVio={model.MaxVio:g}')
        return dict(x=x.copy(), p=problem.power.X, state=problem.state.X,
                    ray_fraction=float(ray_fraction.X), feasible=True)


def build_sequential_region(network, *, budget, monitor, seconds=CASE_TIME_LIMIT, threads=SOLVER_THREADS,
                            threshold=CUT_THRESHOLD, patience=CUT_PATIENCE, tau=REGION_TAU,
                            numeric_focus=0, point_tol=POINT_TOL,
                            ray_threshold=RAY_THRESHOLD):
    """新网架先遍历 G 做射线；每阶段小割后补边界，随即由 physical 搜索查漏。"""
    d = len(network.load_nodes)
    bounds = np.full(d, network.power_limit)
    region = RegionState(bounds, network.power_limit, tau)
    monitor.begin(network, 'socp', budget, region, seconds)
    equations = GridPhysics(network, 'socp')
    oracle = SubProblem(equations, threads=threads, numeric_focus=numeric_focus)
    cache = {}
    powers, applied = [], set()
    initialized = set()
    counts = dict(initial=0, sp=0, cuts=0, ray=0, global_search=0)
    status, certified, coverage = 'time_limit', False, None
    x = None
    stage = 0
    monitor._emit('settings', phase='初始化', threshold=threshold, patience=patience,
                  tau=tau, numeric_focus=numeric_focus, point_tol=point_tol, ray_threshold=ray_threshold,
                  validation_note='独立 SOCP 扫描结果在此显示')

    def check(point):
        # 每个点都读取最新凸包：本批前面的端点认证后，中间点立即免求解。
        if contains([point], region.inner_equations(x))[0]:
            return None
        index = register_power(point*bounds, powers, point_tol)
        power = powers[index]
        key = tuple(x), index
        if key in applied or contains([power/bounds], region.inner_equations(x))[0]:
            return None
        if key not in cache:
            monitor.sp_start(x, power, oracle.calls+1)
            cache[key] = oracle.solve(x, power, time_limit=monitor.remaining(seconds), score_only=True)
            monitor.sp_end(cache[key])
        answer = cache[key]
        if answer['feasible']:
            if region.add_point(x, power/bounds):
                monitor.updated(region, 'feasible', x=x, power=power, checked=answer)
            return None
        return index

    def add_ray(power, phase, initial_sweep=False):
        inner = region.records[tuple(x)]['inner']
        gain = ray_gain(inner, power/bounds)
        gain_threshold = ray_threshold*polytope_volume(inner)
        # 低维认证集先补足维度；满维凸包才使用面积/体积收益筛选。
        if not initial_sweep and gain_threshold > 0. and gain < gain_threshold:
            monitor._emit('ray_skip', phase='射线筛选', ray=None, ray_target=power,
                          ray_gain_bound=gain*np.prod(bounds), ray_gain_threshold=gain_threshold*np.prod(bounds))
            return
        anchor = inner.mean(axis=0)*bounds
        monitor._emit('ray_start', phase=phase,
            ray_gain_bound=gain*np.prod(bounds), ray_gain_threshold=gain_threshold*np.prod(bounds),
            ray=dict(scheme=monitor._scheme(x), anchor=anchor, target=power, p=None))
        answer = ray_support(equations, budget, x, anchor, power,
                             threads=threads, time_limit=monitor.remaining(seconds), numeric_focus=numeric_focus)
        counts['ray'] += 1
        region.add_point(x, answer['p']/bounds)
        monitor._emit('ray_end', phase=phase,
            ray=dict(scheme=monitor._scheme(x), anchor=anchor, target=power, p=answer['p']),
            ray_fraction=answer['ray_fraction'],
            **monitor._geometry(region))

    try:
        # 1. d 个轴向 MP2 和一个总负荷 MP2；后者给出初始网架 A。
        for index, direction in enumerate([*np.eye(d), np.ones(d)]):
            monitor.initializing(index)
            problem = MasterProblem(equations, budget=budget, direction=direction, threads=threads)
            with problem.model:
                counts['initial'] += 1
                answer = problem.solve(time_limit=monitor.remaining(seconds))
            if answer is None:
                region.tighten_bounds(np.zeros(d), 0.)
                status, certified = 'empty', True
                break
            axis_bounds = region.axis_bounds.copy()
            if index < d:
                axis_bounds[index] = min(axis_bounds[index], answer['bound'])
            region.tighten_bounds(axis_bounds, answer['bound'] if index == d else region.total_bound)
            monitor._emit('initial_bounds', phase='初始化', **monitor._geometry(region))
            if index == d:
                x = answer['x']
                region.add_scheme(x, network.decode_plan(x), network.cost_offset+network.cost@x)
                region.add_point(x, answer['p']/bounds)
                monitor.seed(answer, region)

        while x is not None:
            # 2. 新网架先遍历 G 顶点做射线认证；每个网架仅执行一次。
            monitor.remaining(seconds)
            stage += 1
            small_cuts = 0
            monitor._emit('scheme_start', phase='单网架构域', stage=stage,
                active_scheme=monitor._scheme(x), small_cuts=0, area_ratio=None,
                ray=None, sp_point=None, global_point=None, **monitor._geometry(region))
            if tuple(x) not in initialized:
                monitor._emit('initial_sweep_start', phase='首轮射线认证')
                for point in initial_polytope(bounds, region.total_bound, region.axis_bounds):
                    if not contains([point], region.inner_equations(x))[0]:
                        add_ray(point*bounds, '首轮射线认证', initial_sweep=True)
                initialized.add(tuple(x))
                monitor._emit('initial_sweep_end', phase='首轮射线认证', ray=None)
            # 3. 只评分当前 N_x 的顶点；收完可行点后，仅为最大 eta 点生成共享割。
            reason = 'covered'
            while True:
                monitor.remaining(seconds)
                points = stage_candidates(region, x)
                if not len(points):
                    break
                monitor.selecting(len(points), False)
                pending = set()
                for point in points:
                    index = check(point)
                    if index is not None:
                        pending.add(index)
                if not len(stage_candidates(region, x)):
                    break
                pending = [index for index in pending
                           if not contains([powers[index]/bounds], region.inner_equations(x))[0]]
                if not pending:
                    reason = 'point_resolution'
                    break
                index = max(pending, key=lambda i: (cache[(tuple(x), i)]['eta'],
                                                    float(powers[i].sum()), tuple(powers[i])))
                power, checked = powers[index], cache[(tuple(x), index)]
                monitor._emit('cut_start', phase='获选点生成割',
                              sp_point=dict(scheme=monitor._scheme(x), p=power), eta=checked['eta'])
                checked['cut'] = oracle.generate_cut(x, power, checked['cone_normals'],
                                                     time_limit=monitor.remaining(seconds))

                # 冻结 G'；μ(N\G')=μ(N∪G')−μ(G')，二维用面积、三维用体积。
                inner = [row['inner'] for row in region.records.values()]
                inner_area = union_measure(inner, d)
                before = region.records[tuple(x)]['outer']
                region.apply_cut(checked['cut'])
                applied.add((tuple(x), index))
                after = region.records[tuple(x)]['outer']
                denominator = polytope_volume(before) if stage == 1 else inner_area
                removed_area = (polytope_volume(before)-polytope_volume(after) if stage == 1 else
                                union_measure([before, *inner], d)-union_measure([after, *inner], d))
                area_ratio = removed_area/denominator if denominator > 0. else None
                small_cuts = small_cuts+1 if area_ratio is not None and area_ratio < threshold else 0
                counts['cuts'] += 1
                monitor.updated(region, 'cut', x=x, power=power, checked=checked)
                monitor._emit('cut_measure', phase='单网架构域', area_ratio=area_ratio,
                    small_cuts=small_cuts)

                # 加割后沿当前方向补可行点；连续小割只结束局部阶段。
                add_ray(power, '射线补点')
                if small_cuts >= patience:
                    reason = 'area_stagnation'
                    break

            # 4. 补一轮最新 N_x 顶点；实际可行射线点加入 N'_x，不能跨网架取凸包。
            monitor._emit('boundary_start', phase='边界补充', ray=None, sp_point=None)
            targets = {}
            for point in stage_candidates(region, x):
                index = register_power(point*bounds, powers, point_tol)
                targets.setdefault(index, point)
            while targets:
                inner = region.records[tuple(x)]['inner']
                index = max(targets, key=lambda i: ray_gain(inner, targets[i]))
                point = targets.pop(index)
                if not contains([point], region.inner_equations(x))[0]:
                    add_ray(point*bounds, '边界补充')
            monitor._emit('boundary_end', phase='边界补充', ray=None)
            monitor._emit('scheme_end', phase='网架阶段结束', stage_reason=reason, sp_point=None)

            # 5. 在 G 去掉容许扩边的 G' 后查漏；x、p、y 自由，见证带完整物理可行解。
            monitor.global_start(0)
            _, inner_halfspaces = coverage_halfspaces(region)
            problem = RemainingRegionModel(equations, budget, bounds, region.total_bound,
                region.cuts, inner_halfspaces, tau, axis_bounds=region.axis_bounds,
                mode='physical', threads=threads)
            with problem.model:
                counts['global_search'] += 1
                answer = problem.solve(GEOMETRY_TOL, time_limit=monitor.remaining(seconds))
            coverage = answer['bound']
            if not answer['complete']:
                x = answer['x']
                region.add_scheme(x, network.decode_plan(x), network.cost_offset+network.cost@x)
                region.add_point(x, answer['p']/bounds)
            monitor.global_end(answer, region)
            if answer['complete']:
                status, certified = 'certified', True
                break
    except (RegionTimeout, TimeoutError):
        status = 'time_limit'

    counts['sp'] = oracle.calls
    result = dict(status=status, certified=certified, axis_bounds=region.axis_bounds,
        coverage_bound=coverage, counts=counts, timing=monitor.timing(), **region.finish(certified))
    monitor.finish(result, region)
    return result


def scan_path(network, budget, output=SCAN_OUTPUT):
    """同算例、同有序节点、同预算共用一份扫描。"""
    return Path(output)/network.name/'_'.join(map(str, network.load_nodes))/f'budget_{budget:g}.npz'


def save_scan(path, reference):
    """只存功率上限与格点状态；误差指标每次重新计算。"""
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, bounds=reference['bounds'],
                        states=np.asarray(reference['states'], dtype=np.int8))


def scan_reference(network, budget, divisions, bounds, path, *, threads, workers, progress,
                   force_rescan=False):
    """1 读取已有范围；2 覆盖则复用原网格；3 扩大范围或强制要求时重扫。"""
    bounds = np.asarray(bounds)
    if path.exists():
        with np.load(path) as data:
            reference = dict(bounds=data['bounds'], states=data['states'])
        if not force_rescan and np.all(bounds <= reference['bounds']):
            print(f"复用扫描: {path}，上限 {reference['bounds']} kW，网格 {reference['states'].shape}", flush=True)
            return reference
        bounds = np.maximum(bounds, reference['bounds'])
    reference = validate_socp_region(network, budget, divisions, np.ceil(bounds),
                                     threads=threads, workers=workers, progress=progress)
    save_scan(path, reference)
    print(f"保存扫描: {path}，上限 {reference['bounds']} kW，网格 {reference['states'].shape}", flush=True)
    return reference


def run(network, *, budget=BUDGET, divisions=DIVISIONS, show_ui=SHOW_UI, output=None,
        tau=REGION_TAU, time_limit=CASE_TIME_LIMIT, threads=SOLVER_THREADS,
        scan=True, reference=None, scan_workers=SCAN_WORKERS,
        force_rescan=FORCE_RESCAN, scan_output=SCAN_OUTPUT):
    """构域与独立扫描分开计时；只保存一份原生增量回放。"""
    monitor = RunMonitor(output=output, algorithm='逐网架主线')

    def calculate():
        with threadpool_limits(limits=1):
            result = build_sequential_region(network, budget=budget, monitor=monitor,
                seconds=time_limit, threads=threads, tau=tau,
                numeric_focus=3 if isinstance(network, Case33) else 0)
        if scan or reference is not None:
            path = scan_path(network, budget, scan_output)
            if reference is not None and not force_rescan:
                recorded = RunMonitor()
                recorded.load_recording(reference)
                assert recorded.state['network'] == network.name and recorded.state['budget'] == budget
                assert tuple(recorded.state['load_nodes']) == network.load_nodes
                save_scan(path, recorded.state['validation'])
            with threadpool_limits(limits=1):
                reference_grid = scan_reference(network, budget, divisions, result['axis_bounds'], path,
                    threads=threads, workers=scan_workers, progress=monitor.scanning, force_rescan=force_rescan)
                monitor.validation(reference_grid, result)
        return result

    result = monitor.execute(calculate, show_ui=show_ui)
    print(f"{network.name} {len(network.load_nodes)}D: {result['status']}, "
          f"{result['timing']['total_seconds']:.3f}s, {result['counts']}", flush=True)
    if 'validation' in monitor.state:
        for key, metrics in monitor.state['validation']['metrics'].items():
            print(f"{key}: 遗漏 {metrics['mr_percent']}%, 多余 {metrics['fr_percent']}%", flush=True)
    return result


def main(case=None, load_nodes=None, divisions=None, *, dimension=None, seconds=CASE_TIME_LIMIT,
         show=SHOW_UI, scan=True, reference=None, output=OUTPUT):
    case = {FourBus: 'fourbus', Case33: 'case33'}[NETWORK] if case is None else case
    dimension = DIMENSION if dimension is None else dimension
    load_nodes = ({'fourbus': (1, 2, 3), 'case33': (18, 25, 30)}[case][:dimension]
                  if load_nodes is None else tuple(load_nodes))
    network_type, budget = {'fourbus': (FourBus, BUDGET), 'case33': (Case33, Case33.switch_budget)}[case]
    divisions = (SCAN_DIVISIONS[len(load_nodes)] if case == 'case33' else DIVISIONS) if divisions is None else divisions
    return run(network_type(load_nodes=load_nodes), budget=budget, divisions=divisions,
        show_ui=show, output=recording_path(case, output, load_nodes), tau=REGION_TAU,
        time_limit=seconds, threads=SOLVER_THREADS, scan=scan, reference=reference, scan_workers=SCAN_WORKERS,
        force_rescan=FORCE_RESCAN, scan_output=SCAN_OUTPUT)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='二维/三维逐网架规划域与原生过程回放')
    parser.add_argument('--case', choices=('fourbus', 'case33', 'both'),
                        default={FourBus: 'fourbus', Case33: 'case33'}[NETWORK])
    parser.add_argument('--dimension', type=int, choices=(2, 3), default=DIMENSION)
    parser.add_argument('--load-nodes', type=lambda value: tuple(map(int, value.split(','))))
    parser.add_argument('--seconds', type=float, default=CASE_TIME_LIMIT)
    parser.add_argument('--divisions', type=int)
    parser.add_argument('--no-ui', action='store_true', help='无窗口计算，仍保存回放')
    parser.add_argument('--no-scan', action='store_true', help='只构域，不启动独立扫描')
    parser.add_argument('--reference', type=Path, help='复用同案例、同负荷节点及预算的 SOCP 扫描')
    parser.add_argument('--replay', action='store_true', help='只回放，不调用求解器')
    args = parser.parse_args()
    cases = ('fourbus', 'case33') if args.case == 'both' else (args.case,)
    if args.replay:
        SynchronizedReplay([recording_path(case, load_nodes=args.load_nodes, dimension=args.dimension)
                            for case in cases]).root.mainloop()
    else:
        for case in cases:
            main(case, args.load_nodes, args.divisions, dimension=args.dimension, seconds=args.seconds,
                 show=SHOW_UI and not args.no_ui, scan=not args.no_scan, reference=args.reference)
