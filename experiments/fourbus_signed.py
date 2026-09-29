"""固定功率因数的正负接入实验；默认 FourBus，供 Case33 复用。

python experiments/fourbus_signed.py
python experiments/fourbus_signed.py --dimension 2 --seconds 20
python experiments/fourbus_signed.py --replay results/fourbus_signed/mode_1/1_2/recording.json.gz
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from itertools import product
import gzip
import json
from pathlib import Path
import sys
from time import perf_counter
import webbrowser

import gurobipy as gp
import numpy as np
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from Network.four_bus_five_corridor import FourBus
from model import GridPhysics, MasterProblem, RemainingRegionModel, PLANNING_TOL, new_model
from experiments.port_subproblem import PortSubProblem, solve_conic
from main import (stage_candidates, register_power, ray_gain,
                  coverage_halfspaces, union_measure)
from monitor import RunMonitor, RegionTimeout, _plain, _merge
from region import (RegionState, initial_polytope, contains, halfspaces,
                    polytope_volume, GEOMETRY_TOL)

# 1. 实验参数：构域与扫描分别计时；各符号分区共用构域总时限。
MODE = 1
TIME_LIMITS = {2: 20., 3: 60.}
BUDGET = 20000.
SOLVER_THREADS = 4
LOAD_PF = .95
PV_PF = 1.
PV_Q_SIGN = 1.                     # +1 注入无功；-1 吸收无功；PF=1 时均为零
REGION_TAU = .005
CUT_THRESHOLD = .02
CUT_PATIENCE = 3
POINT_TOL = 1e-2
RAY_THRESHOLD = 1e-4               # 射线最大几何收益不足本网架认证测度的 0.01% 时跳过
RAY_CONE_MARGIN = 1e-6             # 射线锥约束的保守内缩量，标幺；不改变全局物理域
SCAN_DIVISIONS = {2: 160, 3: 80}
SCAN_WORKERS = 20
FORCE_RESCAN = False
OUTPUT = ROOT/'results'/'fourbus_signed'


def voltage_flow_bounds(network):
    """由压降等式和电流锥推导有效界；两端电压均有限，反送也有界。"""
    voltage = np.r_[network.vmax, 1.]
    sending = voltage[network.senders[network.type_corridor]]
    receiving = voltage[network.receivers[network.type_corridor]]
    ellmax = (np.sqrt(sending)+np.sqrt(receiving))**2/(network.r**2+network.reactance**2)
    return ellmax, np.sqrt(np.maximum(sending, receiving)*ellmax)


def port_bounds(network):
    """|节点净功率| 不超过相邻走廊各最大两端有功界之和。"""
    _, voltage_bound = voltage_flow_bounds(network)
    capacity = np.minimum(network.capacity, voltage_bound)
    return network.base*np.array([sum(capacity[block].max()
        for c, block in zip(network.corridors, network.type_slices) if node in c.endpoints)
        for node in network.load_nodes])


class PortPhysics(GridPhysics):
    """固定 sign 后用非负幅值装配原方程；真实接入功率 p=sign*u。"""

    def __init__(self, network, sign):
        self.sign = np.asarray(sign)
        network.q_ratio = np.where(self.sign > 0, np.tan(np.arccos(LOAD_PF)),
                                   PV_Q_SIGN*np.tan(np.arccos(PV_PF)))
        network.power_limit = float(port_bounds(network).sum())
        super().__init__(network, 'socp')

    def _build_variable_bounds(self):
        super()._build_variable_bounds()
        net = self.network
        lower_q = np.minimum(net.q_ratio*self.sign*port_bounds(net)/net.base, 0.).sum()
        pmax = net.capacity.copy()
        qmax = np.full(net.n_types, min(net.source_qmax, net.source_smax)-lower_q)
        ellmax = 2*pmax/net.r
        if not np.all(np.isfinite(net.capacity)):
            voltage_ell, _ = voltage_flow_bounds(net)
            ellmax = np.minimum(voltage_ell, (qmax-net.fixed_q.sum()/net.base)/net.reactance)
            voltage = np.r_[net.vmax, 1.]
            end_voltage = np.maximum(voltage[net.senders[net.type_corridor]],
                                     voltage[net.receivers[net.type_corridor]])
            pmax = np.minimum(pmax, np.sqrt(end_voltage*ellmax))
            qmax = np.minimum(qmax, pmax)
        self.pmax, self.qmax, self.ellmax = (dict(zip(self.keys, a)) for a in (pmax, qmax, ellmax))
        self.pmin, self.qmin = dict(zip(self.keys, -pmax)), dict(zip(self.keys, -qmax))
        self.y_lb_global = np.r_[-pmax, -qmax, np.zeros(net.n_types+net.n+2*net.n_corridors)]
        drop = list(self.drop_max.values())
        self.y_ub_global = np.r_[pmax, qmax, ellmax, net.vmax, drop, drop]

    def add_operation(self, model, x, p, eta=0.):
        # 全部实验求解首次建模即采用同一数值设置，不作失败后的参数切换。
        model.Params.BarHomogeneous, model.Params.Aggregate = 1, 0
        model.Params.ScaleFlag = 1 if self.network.name == 'case33bw' else 0
        if isinstance(eta, gp.Var):
            net = self.network
            eta.UB = float(max(np.max(np.abs(net.fixed_p)), np.max(np.abs(net.fixed_q)),
                np.max(np.maximum(1., np.abs(net.q_ratio))*port_bounds(net)))/net.base)
        signed = {i: int(s)*p[i] for i, s in zip(self.network.load_nodes, self.sign)}
        operation = super().add_operation(model, x, signed, eta)
        for e in self.outgoing[self.network.root]:
            for k in self.types[e]:
                model.addConstr(-operation.P[e, k]+self.r[e, k]*operation.ell[e, k]
                                <= self.pmax[e, k]*x[e, k])
                model.addConstr(-operation.Q[e, k]+self.reactance[e, k]*operation.ell[e, k]
                                <= self.qmax[e, k]*x[e, k])
        return operation


def port_ray(equations, x, anchor, power, time_limit):
    """固定网架的近边界射线点；原约束残差达标后才收入认证域。"""
    deadline = perf_counter()+time_limit
    with new_model('port_ray', SOLVER_THREADS) as model:
        choice = dict(zip(equations.keys, x))
        ray_fraction = model.addVar(ub=1., name='ray_fraction')
        p = {i: float(a)+ray_fraction*float(q-a)
             for i, a, q in zip(equations.network.load_nodes, anchor, power)}
        operation = equations.add_operation(model, choice, p)
        model.setObjective(ray_fraction, gp.GRB.MAXIMIZE)
        solution, violation, _ = solve_conic(model, equations, operation, x, {},
            SOLVER_THREADS, deadline, 1e-9, RAY_CONE_MARGIN)
        if violation > PLANNING_TOL:
            raise RuntimeError(f'Ray: MaxVio={violation:g}, anchor={list(anchor)}, p={list(power)}')
        return dict(x=x.copy(), p=anchor+solution[ray_fraction.index]*(power-anchor),
                    state=solution[[v.index for v in operation.state.tolist()]],
                    ray_fraction=float(solution[ray_fraction.index]), feasible=True)


def build_partition(network, sign, mode, seconds, *, budget=BUDGET):
    """2. 一个符号区内逐网架生长；缓存、凸包、割均不跨符号区。"""
    d = len(sign)
    bounds = port_bounds(network) if mode else np.full(d, network.power_limit)
    region = RegionState(bounds, float(bounds.sum()), REGION_TAU)
    monitor = RunMonitor(algorithm='固定功率因数分区')
    monitor.begin(network, 'socp', budget, region, seconds)
    equations = PortPhysics(network, sign) if mode else GridPhysics(network, 'socp')
    oracle = PortSubProblem(equations, threads=SOLVER_THREADS)
    cache, powers, applied, initialized = {}, [], set(), set()
    counts = dict(initial=0, sp=0, cuts=0, ray=0, global_search=0)
    certified, coverage, stage, x = False, None, 0, None

    def register(answer):
        selection = answer['x']
        added = region.add_scheme(selection, network.decode_plan(selection), network.cost_offset+network.cost@selection)
        region.add_point(selection, answer['p']/bounds)
        if added:
            if np.any(network.fixed_p) or np.any(network.fixed_q):
                # Case33 有固定背景负荷：原点也必须实际认证。
                power = np.zeros(d)
                monitor._emit('origin_start', phase='零接入认证', **monitor._geometry(region))
                monitor.sp_start(selection, power, oracle.calls+1)
                checked = oracle.solve(selection, power, time_limit=monitor.remaining(seconds), score_only=True)
                monitor.sp_end(checked)
                if checked['feasible']:
                    region.add_point(selection, power)
            else:
                # FourBus 无固定负荷，零潮流和单位电压解析可行。
                region.add_point(selection, np.zeros(d))
        return selection

    def add_ray(power, phase, initial=False):
        inner = region.records[tuple(x)]['inner']
        gain, threshold = ray_gain(inner, power/bounds), RAY_THRESHOLD*polytope_volume(inner)
        if not initial and threshold > 0. and gain < threshold:
            monitor._emit('ray_skip', phase='射线筛选', ray=None, ray_target=power)
            return
        anchor = inner.mean(axis=0)*bounds
        ray = dict(scheme=monitor._scheme(x), anchor=anchor, target=power, p=None)
        monitor._emit('ray_start', phase=phase, ray=ray)
        answer = port_ray(equations, x, anchor, power, monitor.remaining(seconds))
        counts['ray'] += 1
        region.add_point(x, answer['p']/bounds)
        monitor._emit('ray_end', phase=phase, ray={**ray, 'p': answer['p']},
                      **monitor._geometry(region))

    try:
        # 2.1 轴向及总幅值支撑给初始外包络，最后一个方向给初始网架。
        for j, direction in enumerate([*np.eye(d), np.ones(d)]):
            monitor.initializing(j)
            problem = MasterProblem(equations, budget=budget, direction=direction, threads=SOLVER_THREADS)
            with problem.model:
                counts['initial'] += 1
                answer = problem.solve(time_limit=monitor.remaining(seconds))
            axis_bounds = region.axis_bounds.copy()
            if j < d:
                axis_bounds[j] = min(axis_bounds[j], answer['bound'])
            region.tighten_bounds(axis_bounds, answer['bound'] if j == d else region.total_bound)
            x = register(answer)
            monitor.seed(answer, region)

        while True:
            monitor.remaining(seconds)
            stage += 1
            small_cuts = 0
            monitor._emit('scheme_start', phase='单网架构域', stage=stage,
                active_scheme=monitor._scheme(x), ray=None, sp_point=None, global_point=None,
                **monitor._geometry(region))
            # 2.2 新网架先射线认证本符号区 G 的全部顶点，不按收益筛选。
            if tuple(x) not in initialized:
                for point in initial_polytope(bounds, region.total_bound, region.axis_bounds):
                    if not contains([point], region.inner_equations(x))[0]:
                        add_ray(point*bounds, '首轮射线认证', initial=True)
                initialized.add(tuple(x))
            # 2.3 只评分 N_x 顶点；只为最大 eta 点解一次取割 LP。
            while True:
                monitor.remaining(seconds)
                pending = set()
                for point in stage_candidates(region, x):
                    if contains([point], region.inner_equations(x))[0]:
                        continue
                    # 在评价前统一朝认证凸包内移最多 POINT_TOL/4 kW；只认证实际评价点。
                    power = point*bounds
                    anchor = region.records[tuple(x)]['inner'].mean(axis=0)*bounds
                    power += min(.5, POINT_TOL/(4*np.max(np.abs(anchor-power))))*(anchor-power)
                    index = register_power(power, powers, POINT_TOL)
                    power, key = powers[index], (tuple(x), index)
                    if key in applied or contains([power/bounds], region.inner_equations(x))[0]:
                        continue
                    if key not in cache:
                        monitor.sp_start(x, power, oracle.calls+1)
                        cache[key] = oracle.solve(x, power, time_limit=monitor.remaining(seconds), score_only=True)
                        monitor.sp_end(cache[key])
                    if cache[key]['feasible']:
                        region.add_point(x, power/bounds)
                        monitor.updated(region, 'feasible', x=x, power=power, checked=cache[key])
                    else:
                        pending.add(index)
                pending = [i for i in pending if not contains([powers[i]/bounds], region.inner_equations(x))[0]]
                if not pending:
                    break
                index = max(pending, key=lambda i: cache[(tuple(x), i)]['eta'])
                power, checked = powers[index], cache[(tuple(x), index)]
                checked['cut'] = oracle.generate_cut(x, power, checked['cone_normals'],
                                                     time_limit=monitor.remaining(seconds))
                inner = [r['inner'] for r in region.records.values()]
                before = region.records[tuple(x)]['outer']
                region.apply_cut(checked['cut'])
                after = region.records[tuple(x)]['outer']
                denominator = polytope_volume(before) if stage == 1 else union_measure(inner, d)
                removed = (polytope_volume(before)-polytope_volume(after) if stage == 1 else
                           union_measure([before, *inner], d)-union_measure([after, *inner], d))
                ratio = removed/denominator
                small_cuts = small_cuts+1 if ratio < CUT_THRESHOLD else 0
                applied.add((tuple(x), index))
                counts['cuts'] += 1
                monitor.updated(region, 'cut', x=x, power=power, checked=checked)
                monitor._emit('cut_measure', phase='单网架构域', area_ratio=ratio, small_cuts=small_cuts)
                add_ray(power, '射线补点')
                if small_cuts >= CUT_PATIENCE:
                    break
            # 2.4 最新外顶点补充射线，已进入凸包的点不再求解。
            targets = list(stage_candidates(region, x))
            while targets:
                inner = region.records[tuple(x)]['inner']
                j = max(range(len(targets)), key=lambda k: ray_gain(inner, targets[k]))
                point = targets.pop(j)
                if not contains([point], region.inner_equations(x))[0]:
                    add_ray(point*bounds, '边界补充')
            # 2.5 完整物理查漏：全部合法网架开放，覆盖上界才是停止证书。
            monitor.global_start(0)
            _, faces = coverage_halfspaces(region)
            problem = RemainingRegionModel(equations, budget, bounds, region.total_bound,
                region.cuts, faces, REGION_TAU, axis_bounds=region.axis_bounds,
                mode='physical', threads=SOLVER_THREADS)
            with problem.model:
                counts['global_search'] += 1
                answer = problem.solve(GEOMETRY_TOL, time_limit=monitor.remaining(seconds))
            coverage = answer['bound']
            certified = answer['complete']
            if not certified:
                x = register(answer)
            monitor.global_end(answer, region)
            if certified:
                break
    except (TimeoutError, RegionTimeout):
        pass  # 总时限是用户指定的正常停止条件；仅保存已经取得的证据。
    counts['sp'] = oracle.calls
    result = dict(status='certified' if certified else 'time_limit', certified=certified,
        coverage_bound=coverage, counts=counts, timing=monitor.timing(), **region.finish(certified))
    monitor.finish(result, region)
    return monitor, result


def signed_values(value, sign, prefix, key=''):
    """3. 写文件时一次性把幅值、割系数、网架标签转成带符号坐标。"""
    if value is None:
        return None
    if key in ('scheme', 'active_scheme'):
        return prefix+value
    if key in ('p', 'anchor', 'target', 'vertices', 'ray_target'):
        return (np.asarray(value)*sign).tolist()
    if key == 'cut':
        cut = np.asarray(value).copy()
        cut[1:1+len(sign)] *= sign
        return cut.tolist()
    if isinstance(value, dict):
        return {k: signed_values(v, sign, prefix, k) for k, v in value.items()}
    if isinstance(value, list):
        return [signed_values(v, sign, prefix) for v in value]
    return value


def combine_history(partitions, bounds, signs):
    history, envelopes = [], {}
    d = len(bounds)
    for sign in signs:
        label = ''.join('+' if s > 0 else '-' for s in sign)
        envelopes[label] = [np.asarray(list(product((0., 1.), repeat=d)))*bounds*sign]
    history.append(dict(elapsed=0., patch=dict(event='start', phase='初始化',
        schemes={}, cut_history={}, active_scheme=None, ray=None, sp_point=None, global_point=None,
        global_outer=[p.tolist() for polys in envelopes.values() for p in polys])))
    for sign, offset, monitor, result in partitions:
        label = ''.join('+' if s > 0 else '-' for s in sign)
        prefix = label+':'
        local = {}
        for item in monitor.history:
            _merge(local, item['patch'])
            patch = signed_values(item['patch'], sign, prefix)
            patch['sign'] = list(sign)
            patch['partition'] = label
            if 'schemes' in item['patch']:
                patch['schemes'] = {prefix+k: {**v,
                    'inner': (np.asarray(v['inner']).reshape(-1, d)*sign).tolist(),
                    'outer': (np.asarray(v['outer']).reshape(-1, d)*sign).tolist()}
                    for k, v in item['patch']['schemes'].items()}
            if 'cut_history' in patch:
                patch['cut_history'] = {prefix+k: v for k, v in patch['cut_history'].items()}
            if item['patch']['event'] == 'region_end':
                envelopes[label] = [np.asarray(row['vertices'])*sign for row in result['outer']]
            else:
                envelopes[label] = [initial_polytope(np.asarray(local['bounds']), local['total_bound'],
                                     np.asarray(local['axis_bounds']))*local['bounds']*sign]
            patch['global_outer'] = [p.tolist() for polys in envelopes.values() for p in polys]
            history.append(dict(elapsed=offset+item['elapsed'], patch=patch))
    return history


def scan_line(args, *, network_type=FourBus, budget=BUDGET):
    """4. 独立扫描；同符号、同网架的两个可行端点认证中间格点。"""
    nodes, mode, lower, upper, divisions, index = args
    d = len(nodes)
    states = np.empty(divisions, dtype=np.int8)
    coordinates = lower+(np.arange(divisions)[:, None]*np.eye(d)[-1]
                         +np.r_[index, 0.]+.5)*(upper-lower)/divisions
    with threadpool_limits(limits=1):
        for last_sign in ((1, -1) if mode else (1,)):
            selection = (coordinates[:, -1] >= 0) if last_sign == 1 else (coordinates[:, -1] < 0)
            indices = np.flatnonzero(selection)
            if not len(indices):
                continue
            sign = np.where(coordinates[indices[0]] >= 0., 1, -1)
            network = network_type(load_nodes=nodes)
            equations = PortPhysics(network, sign) if mode else GridPhysics(network, 'socp')
            problem = MasterProblem(equations, power=np.zeros(d), budget=budget, threads=1)
            problem.model.setObjective(0.)  # 参考扫描只判断存在可行网架，不求最小投资。
            with problem.model as model:
                model.ModelName = f'scan_{index}_{last_sign}'
                model.Params.BarQCPConvTol = 1e-10 if network_type.__name__ == 'Case33' else 1e-8
                model.Params.NumericFocus = 0
                model.Params.BarHomogeneous = 1
                model.Params.Aggregate = 1 if network_type.__name__ == 'Case33' else 0
                model.Params.ScaleFlag = 1 if network_type.__name__ == 'Case33' else 0
                model.update()
                fixed = [c for c in model.getConstrs() if c.ConstrName.startswith('fixed_power[')]
                checked, intervals = {}, [(0, len(indices)-1)]
                while intervals:
                    first, last = intervals.pop()
                    for k in (first, last):
                        j = indices[k]
                        if j not in checked:
                            model.ModelName = f'scan_{index}_{j}'
                            model.setAttr('RHS', fixed, coordinates[j]*sign)
                            if network_type is FourBus:
                                model.reset()
                            answer = problem.solve()
                            checked[j] = None if answer is None else answer['x']
                            states[j] = -1 if answer is None else 1
                    a, b = checked[indices[first]], checked[indices[last]]
                    if a is not None and b is not None and np.array_equal(a, b):
                        states[indices[first:last+1]] = 1
                    elif last-first > 1:
                        middle = (first+last)//2
                        intervals.extend(((first, middle), (middle, last)))
    return index, states


def scan_reference(nodes, mode, lower, upper, divisions, force_rescan, *, network_type=FourBus, budget=BUDGET):
    case = network_type.__name__.lower()
    path = ROOT/'results'/'scans'/f'{case}_signed'/f'mode_{mode}'/'_'.join(map(str, nodes))/(
        f'pf_{LOAD_PF:g}_{PV_PF:g}_{PV_Q_SIGN:+g}_budget_{budget:g}.npz')
    lower, upper = np.floor(lower), np.ceil(upper)
    if path.exists():
        with np.load(path) as saved:
            reference = dict(axis_lower=saved['axis_lower'], bounds=saved['bounds'], states=saved['states'])
        if not force_rescan and np.all(lower >= reference['axis_lower']) and np.all(upper <= reference['bounds']):
            print(f'复用扫描：{path}', flush=True)
            return reference
        lower, upper = np.minimum(lower, reference['axis_lower']), np.maximum(upper, reference['bounds'])
    started = perf_counter()
    shape = (divisions,)*len(nodes)
    states = np.empty(shape, dtype=np.int8)
    jobs = ((nodes, mode, lower, upper, divisions, index) for index in np.ndindex(shape[:-1]))
    last_print = started
    with ProcessPoolExecutor(max_workers=SCAN_WORKERS) as pool:
        for count, (index, line) in enumerate(pool.map(partial(scan_line, network_type=network_type, budget=budget), jobs), 1):
            states[index] = line
            if perf_counter()-last_print > 15.:
                print(f'扫描 {count*divisions}/{divisions**len(nodes)}，{perf_counter()-started:.1f}s', flush=True)
                last_print = perf_counter()
    reference = dict(axis_lower=lower, bounds=upper, states=states)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **reference)
    print(f'扫描完成：{perf_counter()-started:.2f}s，{path}', flush=True)
    return reference


def metrics(reference, result):
    states = np.asarray(reference['states'])
    lower, upper = np.asarray(reference['axis_lower']), np.asarray(reference['bounds'])
    indices = np.indices(states.shape).reshape(states.ndim, -1).T
    points = lower+(indices+.5)*(upper-lower)/np.asarray(states.shape)
    truth = states.ravel() == 1
    answer = {}
    for key in ('inner', 'outer'):
        inside = np.zeros(len(points), dtype=bool)
        for row in result[key]:
            inside |= contains(points, halfspaces(row['vertices']))
        missed, extra = np.count_nonzero(truth & ~inside), np.count_nonzero(inside & ~truth)
        answer[key] = dict(mr_percent=100.*missed/truth.sum(),
            fr_percent=100.*extra/inside.sum() if inside.any() else None,
            missed_cells=int(missed), extra_cells=int(extra),
            reference_cells=int(truth.sum()), computed_cells=int(inside.sum()))
    return answer


def run(dimension, *, seconds, mode=MODE, divisions=None, force_rescan=FORCE_RESCAN,
        network_type=FourBus, load_nodes=(1, 2, 3), budget=BUDGET, output=OUTPUT):
    """5. 顺序运行符号区，保存一份完整轨迹和无需求解器的网页回放。"""
    from experiments.signed_replay import export_replay
    nodes = tuple(load_nodes[:dimension])
    signs = [np.asarray(s) for s in product((1, -1), repeat=dimension)] if mode else [np.ones(dimension, int)]
    network = network_type(load_nodes=nodes)
    bounds = port_bounds(network) if mode else np.full(dimension, network.power_limit)
    partitions = []
    started = perf_counter()
    deadline = started+seconds
    with threadpool_limits(limits=1):
        for j, sign in enumerate(signs):
            offset = perf_counter()-started
            allowance = (deadline-perf_counter())/(len(signs)-j)
            monitor, result = build_partition(network_type(load_nodes=nodes), sign, mode, allowance, budget=budget)
            partitions.append((sign, offset, monitor, result))
            print(f'{dimension}D {sign.tolist()}: {result["status"]}, '
                  f'{result["timing"]["total_seconds"]:.2f}s, {result["counts"]}', flush=True)
    elapsed = perf_counter()-started
    result = dict(status='certified' if all(r['certified'] for _, _, _, r in partitions) else 'time_limit',
                  certified=all(r['certified'] for _, _, _, r in partitions),
                  timing=dict(total_seconds=elapsed), inner=[], outer=[],
                  counts={key: sum(r['counts'][key] for _, _, _, r in partitions)
                          for key in ('initial', 'sp', 'cuts', 'ray', 'global_search')})
    for sign, _, _, row in partitions:
        for key in ('inner', 'outer'):
            result[key].extend([{**p, 'sign': sign, 'vertices': np.asarray(p['vertices'])*sign} for p in row[key]])
    history = combine_history(partitions, bounds, signs)
    vertices = np.vstack([p['vertices'] for p in result['outer']])
    lower, upper = vertices.min(axis=0), vertices.max(axis=0)
    output = Path(output)/f'mode_{mode}'/'_'.join(map(str, nodes))
    output.mkdir(parents=True, exist_ok=True)
    data = _plain(dict(version=1, settings=dict(mode=mode, load_nodes=nodes, budget=budget,
        load_pf=LOAD_PF, pv_pf=PV_PF, pv_q_sign=PV_Q_SIGN, time_limit=seconds, threads=SOLVER_THREADS,
        tau=REGION_TAU, ray_threshold=RAY_THRESHOLD, ray_cone_margin=RAY_CONE_MARGIN,
        point_tol=POINT_TOL, axis_lower=lower, bounds=upper),
        history=history, result=result, reference=None))
    recording = output/'recording.json.gz'
    with gzip.open(recording, 'wt', encoding='utf-8') as stream:
        json.dump(data, stream, ensure_ascii=False, allow_nan=False, separators=(',', ':'))
    export_replay(data, output/'replay.html')
    reference = scan_reference(nodes, mode, lower, upper,
        SCAN_DIVISIONS[dimension] if divisions is None else divisions, force_rescan,
        network_type=network_type, budget=budget)
    data['reference'] = _plain(reference)
    data['result']['metrics'] = _plain(metrics(reference, result))
    with gzip.open(recording, 'wt', encoding='utf-8') as stream:
        json.dump(data, stream, ensure_ascii=False, allow_nan=False, separators=(',', ':'))
    export_replay(data, output/'replay.html')
    print(json.dumps(data['result']['metrics'], ensure_ascii=False), flush=True)
    return data


def main(*, network_type=FourBus, load_nodes=(1, 2, 3), budget=BUDGET, time_limits=TIME_LIMITS, output=OUTPUT):
    parser = argparse.ArgumentParser(description=f'{network_type.__name__} 固定功率因数正负接入与回放')
    parser.add_argument('--dimension', type=int, choices=(2, 3))
    parser.add_argument('--seconds', type=float)
    parser.add_argument('--mode', type=int, choices=(0, 1), default=MODE)
    parser.add_argument('--divisions', type=int)
    parser.add_argument('--force-rescan', action='store_true', default=FORCE_RESCAN)
    parser.add_argument('--replay', type=Path)
    args = parser.parse_args()
    if args.replay:
        from experiments.signed_replay import export_replay
        with gzip.open(args.replay, 'rt', encoding='utf-8') as stream:
            data = json.load(stream)
        path = args.replay.parent/'replay.html'
        export_replay(data, path)
        webbrowser.open(path.resolve().as_uri())
        return
    for dimension in ((args.dimension,) if args.dimension else (2, 3)):
        run(dimension, seconds=time_limits[dimension] if args.seconds is None else args.seconds,
            mode=args.mode, divisions=args.divisions, force_rescan=args.force_rescan,
            network_type=network_type, load_nodes=load_nodes, budget=budget, output=output)


if __name__ == '__main__':
    main()
