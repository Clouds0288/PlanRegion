"""配电网可规划域：算法主流程、MP1/MP2/SP 协作、计时与结果存取。

python main.py --recompute --ui
导入本模块只加载定义，不求解、不启动网页。界面与终端共用只读进度事件。
"""
from dataclasses import dataclass
from contextlib import contextmanager
from io import BytesIO
import json
import warnings

from argparse import ArgumentParser, BooleanOptionalAction
from hashlib import sha256
from pathlib import Path
from time import perf_counter

import numpy as np
from threadpoolctl import threadpool_limits

from Network.case33bw import Case33
from Network.four_bus_five_corridor import FourBus
from model import PlanningEquations, PlanningModel, PlanningSP, RemainingRegionModel
from region import GEOMETRY_TOL, RegionState, sample_region, contains
from vertify import METHODS, METHOD_NAMES, validate_ac_region, comparison_labels, validation_summary
from plot import RunMonitor, json_value, print_summary

# 常用设置：直接修改这里，或用命令行参数覆盖。
NETWORK = 'case33'                   # 'case33' 或 'fourbus'
CANDIDATE_COUNT = 16                 # case33 支持 4 / 8 / 16
BUDGETS = None                       # None 使用网架默认预算；也可填列表
DIVISIONS = 8
REGION_TAU = .002                   # SOCP 连续内外域的相对径向精度，与 DIVISIONS 无关。
RECOMPUTE = True                     # False 只读取已有 result.npz
SHOW_UI = True                      # True 自动打开本地实时监视页面
OUTPUT = None                       # None 按网架自动选择结果目录
ROOT = Path(__file__).resolve().parent



def emit(observer, event, **data):
    """不向观察者暴露求解器；数组由监视器在需要时复制。"""
    if observer is not None:
        observer(event, **data)


@contextmanager
def numerical_threads(observer=None):
    """线程限制是计时设置；Windows 的 DLL 探测失败不应中断模型求解。"""
    limiter, error = None, None
    for attempt in range(1, 4):
        try:
            limiter = threadpool_limits(limits=1)
            break
        except OSError as exc:
            if str(exc) not in ('GetModuleFileNameEx failed', 'EnumProcessModulesEx failed'):
                raise
            error = str(exc)
    status = dict(requested_threads=1, applied=limiter is not None, attempts=attempt)
    if limiter is None:
        status['error'] = error
        message = (f'数值库线程检测重试 3 次仍失败（{error}）；使用库当前线程数继续。'
                   '模型与精度不变，本次耗时不能按单线程基准比较。')
        warnings.warn(message, RuntimeWarning, stacklevel=2)
        emit(observer, 'runtime_warning', message=message, thread_control=status)
    try:
        yield status
    finally:
        if limiter is not None:
            limiter.restore_original_limits()


def joint_benders(equations, *, power=None, budget=np.inf, direction=None,
                  cuts=(), start=None, radial_gap_kw=1e-3, observer=None,
                  min_total=None, fixed_x=None, incumbent=None):
    """MP1/MP2 与 SP 的联合割迭代；observer 只接收阶段信息。"""
    emit(observer, 'query_model', message='建立联合割主问题与连续子问题', reused_cuts=len(cuts))
    problem = PlanningModel(equations, power=power, budget=budget,
                            direction=direction, cuts_only=True, min_total=min_total, fixed_x=fixed_x)
    oracle, generated = PlanningSP(equations), []
    minimizing = power is not None or min_total is not None
    mode = 'MP1' if minimizing else 'MP2'
    bound = -np.inf if minimizing else np.inf
    with problem.model:
        for cut in cuts:
            problem.add_cut(cut)
        if start is not None:
            problem.x.Start = start
        if incumbent is not None:
            incumbent_cost = problem.use_incumbent(incumbent, budget=budget, power=power,
                                                    min_total=min_total)
        iteration = 0
        while True:
            iteration += 1
            emit(observer, 'mp_start', iteration=iteration, mode=mode, message=f'求解主问题 {mode}')
            answer = problem.solve()
            if answer is None:
                emit(observer, 'query_end', status='infeasible', message='主问题已证不可行')
                return None, generated
            bound = max(bound, answer['bound']) if minimizing else min(bound, answer['bound'])
            answer['bound'] = bound
            if incumbent is not None and bound >= incumbent_cost-1e-7:
                answer = dict(incumbent, objective=incumbent_cost, bound=bound, status='optimal')
                emit(observer, 'query_end', status='optimal', objective=incumbent_cost, bound=bound,
                     message='MP1 全局费用下界达到已有可行方案费用，最优性已认证')
                return answer, generated
            if answer['x'] is None:
                emit(observer, 'query_end', status='unknown', bound=bound,
                     message=f"无整数候选，保持未确定：{answer['termination']}")
                return answer, generated
            point = answer['p']
            if not minimizing and point.sum() > 0.:
                point = point * max(0., 1. - radial_gap_kw / point.sum())
            emit(observer, 'sp_start', point=point, choice=equations.choice(answer['x']),
                 bound=bound, objective=answer['objective'], message='验证当前选型的运行约束 SP')
            checked = oracle.solve(answer['x'], point)
            emit(observer, 'sp_end', feasible=checked['feasible'], eta=checked['eta'],
                 message=f"SP {'通过' if checked['feasible'] else '未获证'}：{checked['termination']}")
            if checked['feasible']:
                value = answer['objective'] if minimizing else point.sum()
                gap = value-bound if minimizing else bound-value
                tolerance = 1e-7 if minimizing else radial_gap_kw+1e-5
                answer.update(p=point, state=checked['state'], feasible=True, objective=value,
                              status='optimal' if gap <= tolerance else 'feasible')
                emit(observer, 'query_end', status=answer['status'], objective=value, bound=bound,
                     message='获得运行可行证书')
                return answer, generated
            if checked['cut'] is None:
                answer.update(status='unknown', termination=checked['termination'],
                              feasible=False, objective=None)
                emit(observer, 'query_end', status='unknown', message='没有可靠新割，保持未确定')
                return answer, generated
            generated.append(checked['cut'])
            problem.add_cut(checked['cut'])
            emit(observer, 'cut', cut=checked['cut'], selection=answer['x'], point=point,
                 pool_size=len(cuts)+len(generated), message='加入全局有效联合割，继续求解 MP')


class ContinuousRegion:
    def __init__(self, network, method, budget, bounds, query, *, tau=.002,
                 observer=None, cuts=(), start=None, reuse=None):
        self.started = perf_counter()
        self.network, self.method, self.budget = network, method, float(budget)
        self.bounds, self.query = np.asarray(bounds, dtype=float), query
        self.tau = 0. if method == 'linear' else float(tau)
        if not 0 <= self.tau < 1:
            raise ValueError('tau must be in [0, 1)')
        self.observer, self.start = observer, start
        self.equations = PlanningEquations(network, method)
        fingerprint = sha256(method.encode())
        fingerprint.update(str(network.power_limit).encode())
        for array in (self.bounds, self.equations.c, self.equations.F, self.equations.G,
                      self.equations.upper, self.equations.box_constant, self.equations.box_selection,
                      self.equations.cost, self.equations.senders, self.equations.receivers,
                      *self.equations.groups):
            array = np.asarray(array)
            fingerprint.update(str((array.shape, array.dtype)).encode())
            fingerprint.update(array.tobytes())
        self.model_key = fingerprint.hexdigest()
        self.reuse = reuse
        if reuse is not None:
            previous_budget = np.inf if reuse['budget'] is None else reuse['budget']
            if reuse.get('model_key') != self.model_key or previous_budget > self.budget:
                raise ValueError('只能继承同一物理模型、评价箱和不更高预算的证书')
            cuts = [*cuts, *reuse['cuts']]
        unique_cuts = {np.asarray(c, dtype=float).tobytes(): np.asarray(c, dtype=float) for c in cuts}
        self.region = RegionState(self.bounds, network.power_limit, self.tau, unique_cuts.values())
        self.oracle = PlanningSP(self.equations)
        self.counts = dict(mp1=0, mp2=0, residual=0, sp=0, cuts=0, vertices=0,
                           sp_skipped=0, residual_points=0, gap_support_sp=0, covered_support_sp=0,
                           reused_schemes=0, reused_vertices=0,
                           reused_cuts=len(self.region.cuts))
        self.times = dict(mp_seconds=0., sp_seconds=0., geometry_seconds=0.)
        self.maximum = None
        self.clock = None
        self.skipped = set()

    def emit(self, event, **data):
        if self.observer is not None:
            self.observer(event, **data)

    def observe(self, event, **data):
        now = perf_counter()
        if self.clock is not None and event in ('sp_start', 'sp_end', 'query_end'):
            kind, previous = self.clock
            self.times[kind+'_seconds'] += now-previous
            self.clock = None
        if event == 'mp_start':
            self.counts[data.get('mode', 'MP2').lower()] += 1
            self.clock = ('mp', now)
        elif event == 'sp_start':
            self.counts['sp'] += 1
            self.clock = ('sp', now)
        elif event == 'cut':
            self.counts['cuts'] += 1
        self.emit(event, **data)

    
    def geometry(self, **extra):
        before = perf_counter()
        self.emit('geometry', geometry=self.region.geometry(), tau=self.tau,
                  counts_algorithm=self.counts.copy(), pool_size=len(self.region.cuts), **extra)
        self.times['geometry_seconds'] += perf_counter()-before

    
    def add_cut(self, cut, x, point):
        self.counts['cuts'] += 1
        self.emit('cut', cut=cut, selection=x, point=point, pool_size=len(self.region.cuts)+1,
                  message='SP 联合割收紧所有方案的连续外域')
        before = perf_counter()
        self.region.apply_cut(cut)
        self.times['geometry_seconds'] += perf_counter()-before

    def check(self, x, z, *, origin='vertex'):
        point = z*self.bounds
        self.counts['sp'] += 1
        self.counts['vertices' if origin == 'vertex' else 'residual_points'] += 1
        if origin == 'gap-support':
            self.counts['gap_support_sp'] += 1
        self.emit('point', point=point, query=self.counts['sp'], mode='SP-'+origin,
                  message='补齐未覆盖片区所需的同方案支撑证书' if origin == 'gap-support'
                  else '检查全局搜索发现的未覆盖负荷点' if origin == 'residual'
                  else '检查并集尚未覆盖的连续外域顶点')
        self.emit('sp_start', point=point, choice=self.equations.choice(x),
                  message='固定方案，SP 认证连续域候选')
        started = perf_counter()
        checked = self.oracle.solve(x, point)
        self.times['sp_seconds'] += perf_counter()-started
        self.emit('sp_end', feasible=checked['feasible'], eta=checked['eta'],
                  message=f"SP {'通过' if checked['feasible'] else '未通过'}：{checked['termination']}")
        return checked

    
    def skip_covered(self, x, point, owner):
        key = (tuple(x), np.asarray(point).tobytes())
        if key in self.skipped:
            return
        self.skipped.add(key)
        row = self.region.records[owner]
        self.counts['sp_skipped'] += 1
        self.emit('sp_skip', point=point*self.bounds, choice=self.equations.choice(x),
                  covered_by=row['choice'], covered_by_x=row['x'], covered_cost=row['cost'],
                  mode='covered', counts_algorithm=self.counts.copy(),
                  message=f"已由方案 {row['choice']} 的认证内域覆盖，跳过 SP；不归入当前方案内域")

    def refine(self, x, seed=None, *, witness=None, max_checks=96):
        """每批只探索并集外的候选；局部顶点耗尽后仍须全局检查内部空隙。"""
        x = np.asarray(x, dtype=int)
        if self.region.add_scheme(x, self.equations.choice(x).tolist(), self.equations.cost@x):
            self.geometry(message='发现建设方案，建立连续候选多面体')
        row = self.region.records[tuple(x)]
        if seed is not None:
            self.region.add_point(x, np.asarray(seed)/self.bounds)
            self.geometry(message='登记 MP 已有可行证书，无需再次调用 SP')
        pending = None if witness is None else (1-self.tau)*np.asarray(witness)/self.bounds
        for _ in range(max_checks):
            if not len(row['outer']):
                break
            before = perf_counter()
            origin = 'vertex'
            point = None
            if pending is not None:
                owner = self.region.covering_schemes([pending], preferred=x)[0]
                if owner is not None:
                    pending = None
                else:
                    support = self.region.witness_support(x, pending)
                    for candidate in support:
                        if len(row['inner']) and contains([candidate], self.region.inner_equations(x))[0]:
                            continue
                        point, origin = candidate, 'gap-support'
                        support_owner = self.region.covering_schemes([candidate], preferred=x)[0]
                        if support_owner is not None:
                            self.counts['covered_support_sp'] += 1
                        self.emit('gap_support', point=candidate*self.bounds,
                                  uncovered_witness=pending*self.bounds,
                                  support_points=np.asarray(support)*self.bounds,
                                  covered_by=None if support_owner is None else self.region.records[support_owner]['choice'],
                                  choice=self.equations.choice(x),
                                  message='全局搜索已发现内部未覆盖片区，补齐其同方案单纯形证书')
                        break
                    if point is None:
                        point, origin = pending, 'residual'
                        pending = None
            if point is None:
                origin = 'vertex'
                point = self.region.next_point(x, skipped=lambda p, owner: self.skip_covered(x, p, owner))
            self.times['geometry_seconds'] += perf_counter()-before
            if point is None:
                self.geometry(message='该方案候选顶点已被全局内域覆盖；转向全局搜索，继续检查内部空隙')
                return True
            checked = self.check(x, point, origin=origin)
            if checked['feasible']:
                before = perf_counter()
                self.region.add_point(x, point)
                self.times['geometry_seconds'] += perf_counter()-before
                self.geometry(message='同方案可行证书扩展连续内域')
            elif checked['cut'] is not None:
                old = row['outer'].copy()
                self.add_cut(checked['cut'], x, point*self.bounds)
                if pending is not None:
                    cut = checked['cut']
                    if cut[0]+cut[1:4]@(pending*self.bounds/(1-self.tau))+cut[4:]@x < -1e-12:
                        pending = None  # 新割已排除原见证，不再认证其过期支撑点。
                self.geometry(message='有效割完成连续多面体裁剪')
                if np.array_equal(old, row['outer']):
                    self.emit('unknown', message='割未产生可靠几何进展，保留未确定区域')
                    return False
            else:
                self.emit('unknown', message='顶点既无可行证书也无有效割，保留未确定区域')
                return False
        self.geometry(message='本批未覆盖候选处理完成，转向全局剩余域搜索' if len(row['outer'])
                      else '有效割已排除该方案的整个候选域')
        return True

    
    def residual(self):
        self.counts['residual'] += 1
        self.emit('mp_start', mode='residual', iteration=self.counts['residual'],
                  message='全局搜索尚未覆盖的连续区域（包含未发现方案）')
        state = self.region
        problem = RemainingRegionModel(self.equations, self.budget, self.bounds, state.total_bound,
                                       state.cuts, state.inner_halfspaces(), self.tau)
        with problem.model:
            answer = problem.solve(GEOMETRY_TOL)
        self.times['mp_seconds'] += answer.pop('solve_seconds')
        return answer

    
    def finish(self, status, coverage, maximum):
        before = perf_counter()
        domain = self.region.finish(status == 'certified')
        self.times['geometry_seconds'] += perf_counter()-before
        elapsed = perf_counter()-self.started
        self.times['total_seconds'] = elapsed
        self.times['other_seconds'] = max(0., elapsed-sum(v for k, v in self.times.items() if k != 'total_seconds'))
        result = dict(method=self.method, budget=None if np.isinf(self.budget) else self.budget,
                      model_key=self.model_key, strategy='union',
                      reuse_from_budget=None if self.reuse is None else self.reuse['budget'],
                      timing_mode='independent' if self.reuse is None else 'incremental',
                      status=status, tau=self.tau, geometry_tolerance=GEOMETRY_TOL,
                      coverage_bound=coverage, max_total=maximum, max_total_bound=self.region.total_bound,
                      max_point=None if self.maximum is None else self.maximum['p'].tolist(),
                      max_choice=None if self.maximum is None else self.equations.choice(self.maximum['x']).tolist(),
                      max_cost=None if self.maximum is None else float(self.equations.cost@self.maximum['x']),
                      counts=self.counts.copy(), timing=self.times.copy(), **domain)
        self.emit('region_end', region=result, geometry=result['geometry'], global_outer=result['outer'],
                  coverage_complete=status == 'certified', counts_algorithm=self.counts.copy(),
                  message='连续内外域已获全局覆盖证书' if status == 'certified' else '连续域存在未确定部分')
        return result

    def solve(self):
        self.geometry(message='初始化连续外域，开始 MP2 → MP1 → SP 协作')
        answer, new = self.query(self.equations, budget=self.budget, cuts=self.region.cuts,
                                 start=self.start, observer=self.observe)
        self.region.cuts.extend(new)
        if answer is None:
            self.region.total_bound = 0.
            return self.finish('certified', None, 0.)
        if not answer['feasible']:
            return self.finish('unknown', None, None)
        self.region.total_bound = min(self.region.total_bound, answer['bound'])
        self.maximum = answer
        maximum = float(answer['p'].sum())
        self.emit('maximum', max_total=maximum, max_total_bound=self.region.total_bound,
                  point=answer['p'], choice=self.equations.choice(answer['x']),
                  message='MP2 获得最大总负荷的可行下界与全局上界')
        # MP1 固定总量而非节点分配；避免把最大总量点当成整个三维域。
        cheapest, new = self.query(self.equations, min_total=max(0., maximum-1e-6),
                                   budget=self.budget, cuts=self.region.cuts, start=answer['x'],
                                   incumbent=answer, observer=self.observe)
        self.region.cuts.extend(new)
        if self.reuse is not None:
            before = perf_counter()
            for certificate in self.reuse['certificates']:
                x = np.asarray(certificate['x'], dtype=int)
                cost = float(self.equations.cost@x)
                points = np.asarray(certificate['inner']).reshape(-1, 3)
                if cost > self.budget+1e-8 or not len(points):
                    continue
                self.region.add_scheme(x, self.equations.choice(x).tolist(), cost)
                self.region.add_point(x, points)
                self.counts['reused_schemes'] += 1
                self.counts['reused_vertices'] += len(points)
            self.times['geometry_seconds'] += perf_counter()-before
            self.geometry(reuse_from_budget=self.reuse['budget'],
                          message='继承更低预算下同一模型的认证内域；只探索当前预算新增部分')
        seeds = [answer]
        if cheapest is not None and cheapest['feasible']:
            seeds.insert(0, cheapest)
        for seed in seeds:
            if not self.refine(seed['x'], seed['p']):
                return self.finish('unknown', None, maximum)
        while True:
            witness = self.residual()
            self.emit('coverage', coverage_bound=witness['bound'], coverage_complete=witness['complete'],
                      message='已证明全部整数方案的连续外域被覆盖' if witness['complete'] else '仍存在未覆盖的连续候选区域')
            if witness['complete']:
                return self.finish('certified', witness['bound'], maximum)
            if witness['x'] is None:
                return self.finish('unknown', witness['bound'], maximum)
            self.emit('candidate', point=witness['p'], choice=self.equations.choice(witness['x']),
                      message='剩余区域搜索发现候选建设方案')
            before = self.region.progress
            if not self.refine(witness['x'], witness=witness['p']):
                return self.finish('unknown', witness['bound'], maximum)
            after = self.region.progress
            if before == after:
                self.emit('unknown', message='剩余区域搜索与几何容差不一致，停止并保留未确定')
                return self.finish('unknown', witness['bound'], maximum)


def build_continuous_region(network, method, budget, bounds, *, tau=REGION_TAU,
                            observer=None, budget_index=0, budgets=None, reuse=None):
    """单预算连续构域；混合方法只传递有效割，SOCP 内域重新认证。"""
    if method not in ('linear', 'socp', 'hybrid'):
        raise ValueError('continuous method must be linear, socp or hybrid')
    started = perf_counter()
    phases = ('linear', 'socp') if method == 'hybrid' else (method,)
    cuts, phase_results = [], []
    for phase_number, phase in enumerate(phases, 1):
        emit(observer, 'phase_start', method=method, phase=phase, phase_number=phase_number,
             phase_count=len(phases), budget_index=budget_index, budget=budget,
             budgets=[budget] if budgets is None else budgets, bounds=bounds, load_nodes=network.load_nodes,
             representation='continuous', message=f'{method} · 预算 {budget:g} · {phase} 连续构域')
        solver = ContinuousRegion(network, phase, budget, bounds, joint_benders,
                                  tau=tau, observer=observer, cuts=cuts,
                                  reuse=None if reuse is None else reuse.get(phase))
        region = solver.solve()
        if reuse is not None:
            reuse[phase] = region
        phase_results.append(region)
        cuts = [np.asarray(c) for c in region['cuts']]
        emit(observer, 'phase_end', message=f'{phase} 连续构域：{region["status"]}')
    region = dict(region, method=method, budget_index=budget_index)
    if method == 'hybrid':
        region['phase_timing'] = [r['timing'] for r in phase_results]
        region['timing'] = {k: sum(r['timing'][k] for r in phase_results) for k in region['timing']}
        region['counts'] = {k: sum(r['counts'][k] for r in phase_results) for k in region['counts']}
        region['linear_status'] = phase_results[0]['status']
    region['timing']['total_seconds'] = perf_counter()-started
    return region


@dataclass  # 仅保存三态网格与复现配置，不维护方案表或逐查询日志。
class BenchmarkResult:  # 指标与颜色由唯一基础数据推导。
    states: np.ndarray  # (方法,预算,N,N,N)，-1 域外、0 未确定、1 域内。
    metadata: dict  # 网架、预算、评价箱、源码指纹和各方法总耗时。

    @property  # 将 JSON 中的无限预算恢复为数值。
    def budgets(self):  # 返回计算使用的预算序列。
        return tuple(np.inf if b is None else b for b in self.metadata['budgets'])  # null 表示无限预算。

    @property  # 坐标上界只在配置中保存一次。
    def bounds(self):  # 返回公共评价箱，单位 kW。
        return np.asarray(self.metadata['bounds'])  # 与三维坐标顺序一致。

    @property  # 负荷节点直接读取配置。
    def load_nodes(self):  # 返回三条坐标轴对应的实际节点。
        return tuple(self.metadata['load_nodes'])  # 不另维护节点副本。

    @property  # 步长由评价箱和分辨率推导。
    def spacing(self):  # 返回三个方向的单元宽度。
        return self.bounds/self.metadata['divisions']  # 等体积网格的长度尺度。

    
    @property
    def labels(self):
        return comparison_labels(self.states)

    
    @property
    def summary(self):
        return validation_summary(self.states, self.metadata)

    def save(self, folder):  # 一次实验只写一个原始结果文件。
        folder = Path(folder)  # 统一路径表示。
        folder.mkdir(parents=True,exist_ok=True)  # 创建本次实验目录。
        temporary = folder/'result.npz.tmp'  # 完整写入后再替换正式文件。
        with temporary.open('wb') as stream:  # 中断写入不破坏上一次完整结果。
            np.savez_compressed(stream,states=self.states,metadata=json.dumps(self.metadata,ensure_ascii=False))  # 仅存三态标签与唯一配置。
        temporary.replace(folder/'result.npz')  # 原子发布本次结果。

    @classmethod  # 无需建立空实例即可读取结果。
    def load(cls, folder):  # 只读取当前简化格式，不保留旧方案表兼容分支。
        with np.load(BytesIO((Path(folder)/'result.npz').read_bytes()),allow_pickle=False) as data:  # 释放文件句柄后解码。
            return cls(data['states'],json.loads(str(data['metadata'])))  # 指标和图形均重新从基础数据生成。


def evaluation_bounds(network, observer=None):
    """用三个无预算 LP 轴向全局上界确定公共评价箱。"""
    equations, bounds = PlanningEquations(network, 'linear'), []
    for axis, direction in enumerate(np.eye(3), 1):
        emit(observer, 'preparation', message=f'公共评价箱：求解轴向上界 {axis}/3')
        problem = PlanningModel(equations, direction=direction)
        with problem.model:
            answer = problem.solve()
        if answer is None or answer['bound'] is None or not np.isfinite(answer['bound']):
            raise RuntimeError('No finite planning bound for the common evaluation box')
        bounds.append(min(network.power_limit, answer['bound']))
    return np.ceil(np.asarray(bounds)/10.)*10.


def run(network=None, *, budgets=None, divisions=DIVISIONS, recompute=RECOMPUTE,
        show_ui=SHOW_UI, output=None, plots=True, open_browser=True, keep_ui=False,
        tau=REGION_TAU, reuse_budgets=False):
    """先构造连续域，再独立采样校核；完整算法事件用于同一 HTML 的回放。"""

    network = Case33(candidate_count=CANDIDATE_COUNT) if network is None else network
    if not np.isfinite(tau) or not 0 <= tau < 1:
        raise ValueError('tau 必须是 [0, 1) 内的有限数')
    budgets = np.asarray(network.budgets if budgets is None else budgets, dtype=float)
    if (not isinstance(divisions, (int, np.integer)) or divisions < 1
            or budgets.ndim != 1 or not len(budgets) or np.any(np.isnan(budgets))
            or np.any(budgets < 0) or np.any(budgets[1:] <= budgets[:-1])):
        raise ValueError('divisions 必须是正整数；budgets 必须是严格递增的非负预算列表')
    output = Path(output) if output is not None else ROOT/'results'/network.name/f'planning_{len(network.line_options)}'
    output = output.resolve()
    with RunMonitor(show_ui=show_ui, output=output, open_browser=open_browser, record=recompute) as monitor:
        monitor('start', network=network.name, cost_unit=network.cost_unit, message='准备计算' if recompute else '读取已有结果')
        if recompute:
            with numerical_threads(monitor) as thread_control:
                prepared = perf_counter()
                bounds = evaluation_bounds(network, monitor)
                metadata = dict(network=network.name, planning=True, cost_unit=network.cost_unit,
                                load_nodes=list(network.load_nodes), candidate_count=len(network.line_options),
                                budgets=[None if np.isinf(b) else float(b) for b in budgets],
                                divisions=int(divisions), bounds=bounds.tolist(),
                                preparation_seconds=perf_counter()-prepared, seconds={}, show_ui=show_ui,
                                mode='continuous', tau=tau, continuous=[], reuse_budgets=bool(reuse_budgets),
                                thread_control=thread_control,
                                strategy='union')
                metadata['hashes'] = {name: sha256((ROOT/name).read_bytes()).hexdigest()
                                      for name in ('main.py', 'model.py', 'region.py', 'vertify.py',
                                                   'plot.py', 'live_view.html', *network.sources)}
                result = BenchmarkResult(np.zeros((len(METHODS), len(budgets))+(divisions,)*3, dtype=np.int8), metadata)
                points = (np.indices((divisions,)*3).reshape(3, -1).T+.5)*bounds/divisions
                for number, method in enumerate(('linear', 'socp', 'hybrid'), 1):
                    monitor('method_start', method=method, method_number=number, method_count=4,
                            message=f'开始方法 {number}/4：{dict(zip(METHODS, METHOD_NAMES))[method]}')
                    method_seconds = 0.
                    reuse = {} if reuse_budgets else None
                    for j, budget in enumerate(budgets):
                        region = build_continuous_region(network, method, budget, bounds, tau=tau,
                                                         observer=monitor, budget_index=j, budgets=budgets, reuse=reuse)
                        method_seconds += region['timing']['total_seconds']
                        result.states[METHODS.index(method), j] = sample_region(region, points, bounds).reshape((divisions,)*3)
                        metadata['continuous'].append(json_value(region))
                        monitor('budget_end', results=metadata['continuous'], message=f'预算 {budget:g} 完成，保存连续域与实测时间')
                    metadata['seconds'][method] = method_seconds
                    monitor('method_end', seconds=metadata['seconds'][method], message=f'{method} 方法完成')
                monitor('method_start', method='ac', method_number=4, method_count=4,
                        representation='validation', message='独立 AC 校核（采样不参与连续构域）')
                started = perf_counter()
                result.states[METHODS.index('ac')] = validate_ac_region(network, budgets, divisions, bounds, observer=monitor)
                metadata['seconds']['ac'] = perf_counter()-started
                monitor('method_end', seconds=metadata['seconds']['ac'], message='独立 AC 校核完成')
                metrics = {(row['method'], row['budget']): row for row in result.summary}
                for region in metadata['continuous']:
                    region['validation'] = metrics[(region['method'], region['budget'])]
                metadata['wall_seconds_before_export'] = perf_counter()-monitor.started
            monitor('saving', message='保存 result.npz')
            result.save(output)
            monitor('completed', results=metadata['continuous'], metadata=metadata,
                    message='连续域计算与独立校核完成，可回放全部算法事件')
        else:
            result = BenchmarkResult.load(output)
            if (output/'replay.json').exists():
                monitor.load_recording(output/'replay.json')
            else:
                monitor('completed', network=result.metadata['network'], results=result.metadata.get('continuous', []),
                        bounds=result.bounds, budgets=result.budgets, load_nodes=result.load_nodes,
                        message='已载入结果；该历史文件没有算法回放记录')
        print_summary(result)
        if plots or show_ui:
            exporting = perf_counter()
            monitor.save_snapshot()
            print(f'回放导出耗时：{perf_counter()-exporting:.3f}s', flush=True)
        if show_ui:
            if keep_ui:
                try:
                    input('计算已完成，网页可继续查看；按 Enter 或 Ctrl+C 关闭监视服务。\n')
                except (EOFError, KeyboardInterrupt):
                    pass
    return result


def main(argv=None):
    import sys
    # Windows 重定向到管道时可能默认 cp1252，统一中文日志及 --help 的编码。
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, 'reconfigure'):
            stream.reconfigure(encoding='utf-8')
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--network', choices=('case33', 'fourbus'), default=NETWORK)
    parser.add_argument('--candidates', type=int, choices=(4, 8, 16), default=CANDIDATE_COUNT)
    parser.add_argument('--divisions', type=int, default=DIVISIONS)
    parser.add_argument('--tau', type=float, default=REGION_TAU, help='SOCP 连续域径向精度，默认 0.002')
    parser.add_argument('--reuse-budgets', action=BooleanOptionalAction, default=False,
                        help='递增预算继承同模型证书；--no-reuse-budgets 用于独立预算计时')
    parser.add_argument('--budgets', type=float, nargs='+', default=BUDGETS, help='递增预算，inf 表示无限')
    parser.add_argument('--ui', action=BooleanOptionalAction, default=SHOW_UI, help='显示实时本地网页')
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--recompute', dest='recompute', action='store_true')
    mode.add_argument('--load', dest='recompute', action='store_false')
    parser.set_defaults(recompute=RECOMPUTE)
    parser.add_argument('--output', type=Path, default=OUTPUT)
    parser.add_argument('--no-plots', action='store_true')
    parser.add_argument('--no-open-browser', action='store_true')
    parser.add_argument('--no-hold', action='store_true', help='计算完成后直接退出网页服务；仍保存静态快照')
    args = parser.parse_args(argv)
    network = FourBus() if args.network == 'fourbus' else Case33(candidate_count=args.candidates)
    try:
        run(network, budgets=args.budgets, divisions=args.divisions, recompute=args.recompute,
            show_ui=args.ui, output=args.output, plots=not args.no_plots,
            open_browser=not args.no_open_browser, keep_ui=not args.no_hold, tau=args.tau,
            reuse_budgets=args.reuse_budgets)
    except KeyboardInterrupt:
        return 130
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
