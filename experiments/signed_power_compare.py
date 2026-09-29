"""二维正负净功率：独立比较主线调度与 KKT 认证并集排除，不改生产入口。

python experiments/signed_power_compare.py --case both --seconds 1000 --divisions 80
python experiments/signed_power_compare.py --replay results/signed_power/fourbus/comparison.json.gz
"""
import argparse
import gzip
from itertools import product
import json
from pathlib import Path
import sys
from time import perf_counter

import gurobipy as gp
from gurobipy import GRB
import numpy as np
from scipy.spatial import ConvexHull
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from Network.four_bus_five_corridor import FourBus
from Network.case33bw import Case33
from model import GridPhysics, MasterProblem, SubProblem, PLANNING_TOL
from region import RegionState, GEOMETRY_TOL, clip_polytope, contains, halfspaces
from experiments.fourbus_outer_kkt import KktViolation, UncoveredKktViolation, certify_load
from experiments.fourbus_outer_volume import budget_schemes

DIRECTIONS = np.array([[1., 0.], [-1., 0.], [0., 1.], [0., -1.],
                       [1., 1.], [-1., -1.], [1., -1.], [-1., 1.]])
CORNERS = np.array(list(product((0., 1.), repeat=2)))


class SignedPhysics(GridPhysics):
    """只替换实验中的有效变量界；共用全部原始 DistFlow/SOCP 方程和 SP。"""

    def __init__(self, network, power_lower, power_upper):
        self.power_lower = np.asarray(power_lower, dtype=float)
        self.power_upper = np.asarray(power_upper, dtype=float)
        super().__init__(network, 'socp')

    def _build_variable_bounds(self):
        super()._build_variable_bounds()
        net = self.network
        self.source_pmax = min(net.source_pmax, net.source_smax)
        self.source_qmax = min(net.source_qmax, net.source_smax)
        lower_p, lower_q = net.loads(self.power_lower)
        available_p = self.source_pmax-np.minimum(lower_p, 0.).sum()
        available_q = self.source_qmax-np.minimum(lower_q, 0.).sum()
        loss_max = self.source_pmax-lower_p.sum()
        pmax = np.minimum(net.capacity, available_p)
        qmax = np.full(net.n_types, available_q)
        ellmax = np.minimum(loss_max, 2*pmax)/net.r
        self.pmax, self.qmax, self.ellmax = (dict(zip(self.keys, a)) for a in (pmax, qmax, ellmax))
        self.pmin, self.qmin = dict(zip(self.keys, -pmax)), dict(zip(self.keys, -qmax))
        drop = list(self.drop_max.values())
        self.y_lb_global = np.r_[-pmax, -qmax, np.zeros(net.n_types+net.n+2*net.n_corridors)]
        self.y_ub_global = np.r_[pmax, qmax, ellmax, net.vmax, drop, drop]

    def add_operation(self, model, x, p, eta=0.):
        operation = super().add_operation(model, x, p, eta)
        net = self.network
        roots = [(e, k) for e in self.outgoing[net.root] for k in self.types[e]]
        model.addConstr(gp.quicksum(operation.P[key] for key in roots) >= -self.source_pmax)
        model.addConstr(gp.quicksum(operation.Q[key] for key in roots) >= -self.source_qmax)
        for e, k in roots:
            model.addConstr(-operation.P[e, k]+self.r[e, k]*operation.ell[e, k]
                            <= self.pmax[e, k]*x[e, k])
            model.addConstr(-operation.Q[e, k]+self.reactance[e, k]*operation.ell[e, k]
                            <= self.qmax[e, k]*x[e, k])
        return operation


class SignedMaster(MasterProblem):
    def __init__(self, equations, **kwargs):
        super().__init__(equations, **kwargs)
        self.power.LB, self.power.UB = equations.power_lower, equations.power_upper


def expanded_faces(points, padding):
    shifts = (2*CORNERS-1)*padding
    return ConvexHull((points[:, None, :]+shifts[None, :, :]).reshape(-1, 2)).equations


class SignedRegion(RegionState):
    """xi=(p-axis_lower)/bounds；原始 kW 割仅在几何入口转换一次。"""

    def __init__(self, axis_lower, axis_bounds, cuts):
        self.axis_lower = np.asarray(axis_lower)
        super().__init__(np.asarray(axis_bounds)-axis_lower, np.inf, 0., cuts)
        self.axis_bounds = np.asarray(axis_bounds)

    def normalize(self, power):
        return (np.asarray(power)-self.axis_lower)/self.bounds

    def physical(self, point):
        return self.axis_lower+np.asarray(point)*self.bounds

    def inset(self, x, point, padding, witness=None):
        """朝同网架认证凸包移至多 padding/4；搜索/停止仍检查未移动的域。"""
        inner = self.records[tuple(x)]['inner']
        anchor = inner.mean(axis=0) if len(inner) else self.normalize(witness['p'])
        distance = np.max(np.abs(anchor-point))
        fraction = min(.5, padding/(4*distance))
        return point+fraction*(anchor-point)

    def outer_for(self, x):
        outer = CORNERS.copy()
        for cut in self.cuts:
            outer = clip_polytope(outer, cut[0]+cut[1:3]@self.axis_lower+cut[3:]@x,
                                  cut[1:3]*self.bounds)
        return outer

    def add_scheme(self, x, choice, cost):
        key = tuple(x)
        if key in self.records:
            return False
        self.records[key] = dict(x=np.asarray(x), choice=choice, cost=float(cost),
            inner=np.empty((0, 2)), outer=self.outer_for(x), inner_equations=None, inner_box=None)
        return True

    def apply_cut(self, cut):
        self.cuts.append(np.asarray(cut))
        for row in self.records.values():
            row['outer'] = clip_polytope(row['outer'],
                cut[0]+cut[1:3]@self.axis_lower+cut[3:]@row['x'], cut[1:3]*self.bounds)

    def covered(self, points, padding, own=None):
        points = np.asarray(points).reshape(-1, 2)
        covered = np.zeros(len(points), dtype=bool)
        rows = self.records.values() if own is None else [self.records[tuple(own)]]
        for row in rows:
            if len(row['inner']):
                covered |= contains(points, expanded_faces(row['inner'], padding))
        return covered

    def snapshot(self):
        return [dict(x=row['x'].copy(), inner=self.physical(row['inner']),
                     outer=self.physical(row['outer'])) for row in self.records.values()]


class SignedKkt(KktViolation):
    """原 KKT 的负荷下界与排除面作仿射迁移；不枚举网架或负荷点。"""

    def __init__(self, equations, budget, region, *, padding=1e-5, threads=4):
        self.axis_lower = region.axis_lower
        super().__init__(equations, budget, region.axis_bounds, equations.network.power_limit, threads=threads)
        self.problem.power.LB = self.axis_lower
        net = equations.network
        loads = net.loads(region.physical(CORNERS))
        self.violation.UB = float(max(np.max(np.abs(row)) for row in loads))
        for cut in region.cuts:
            self.problem.add_cut(cut)
            self._add_outer_row(-cut[1:3], cut[0], cut[3:])
        for k, row in enumerate(region.records.values()):
            if not len(row['inner']):
                continue
            faces = expanded_faces(row['inner'], padding)
            select = self.model.addVars(len(faces), vtype=GRB.BINARY, name=f'outside[{k}]')
            self.model.addConstr(select.sum() == 1.)
            for f, face in enumerate(faces):
                big_m = max(0., -float(face[2]+np.minimum(face[:2], 0.).sum()))
                outer_p = -face[:2]/region.bounds
                outer_rhs = float(face[2]+outer_p@region.axis_lower)+big_m
                self.model.addConstr(gp.LinExpr(outer_p, self.problem.power.tolist())
                                    <= outer_rhs-big_m*select[f])
                self._add_outer_row(outer_p, outer_rhs, np.zeros(net.n_types),
                                    selector=select[f], selector_rhs=-big_m)
        self._bound_outer_duals()
        relaxed = self.template.sp_eta < 0.
        self.padding_bound = float(padding*np.max(np.abs(self.template.sp_p[relaxed])@region.bounds
                                  /(-self.template.sp_eta[relaxed]))+PLANNING_TOL)
        self.model.update()

    def _add_outer_row(self, outer_p, outer_rhs, outer_x, **kwargs):
        index = len(self.outer_dual)
        if index < 4 and index % 2 == 0:
            outer_rhs = -float(self.axis_lower[index//2])
        super()._add_outer_row(outer_p, outer_rhs, outer_x, **kwargs)

    solve = UncoveredKktViolation.solve


class Trace:
    """实验的单份最简轨迹；只保存阶段、候选、几何变化，扫描不生成回放帧。"""

    def __init__(self, name, seconds):
        self.name, self.seconds = name, seconds
        self.started = perf_counter()
        self.events = []
        self.counts = dict(initialization=0, support=0, sp=0, cuts=0, global_search=0, certification=0)
        self.last_print = 0.

    def remaining(self, limit=1000.):
        remaining = self.seconds-(perf_counter()-self.started)
        if remaining <= 0.:
            raise TimeoutError('Experiment time limit')
        return min(limit, remaining)

    def event(self, phase, region=None, *, x=None, power=None, cut=None, eta=None, bound=None):
        elapsed = perf_counter()-self.started
        self.events.append(dict(phase=phase, seconds=elapsed, x=x, p=power, cut=cut, eta=eta,
                                bound=bound, records=None if region is None else region.snapshot()))
        if phase in ('initialization', 'global_search', 'finish') or elapsed-self.last_print > 10.:
            print(f'{self.name}: {phase}, t={elapsed:.2f}s, SP={self.counts["sp"]}, '
                  f'cuts={self.counts["cuts"]}, p={power}, bound={bound}', flush=True)
            self.last_print = elapsed


def initialize(equations, budget, threads, trace):
    answers, cuts = [], []
    for direction in DIRECTIONS:
        problem = SignedMaster(equations, budget=budget, direction=direction, threads=threads)
        with problem.model:
            answer = problem.solve(time_limit=trace.remaining())
        answers.append(answer)
        cuts.append(np.r_[answer['bound'], -direction, np.zeros(equations.network.n_types)])
        trace.counts['initialization'] += 1
        trace.event('initialization', x=answer['x'], power=answer['p'])
    lower = -np.array([answers[1]['bound'], answers[3]['bound']])
    upper = np.array([answers[0]['bound'], answers[2]['bound']])
    region = SignedRegion(lower, upper, cuts)
    for answer in answers:
        register(equations, region, answer['x'])
        region.add_point(answer['x'], region.normalize(answer['p']))
    trace.event('initialized', region)
    return region


def register(equations, region, x):
    net = equations.network
    return region.add_scheme(x, net.decode_plan(x), net.cost_offset+net.cost@x)


def physical_search(equations, budget, region, padding, threads, trace):
    problem = SignedMaster(equations, budget=budget, cuts=region.cuts, threads=threads)
    model = problem.model
    problem.power.LB, problem.power.UB = region.axis_lower, region.axis_bounds
    gamma = model.addVar(lb=-4., ub=4., name='uncovered_distance')
    for k, row in enumerate(region.records.values()):
        if not len(row['inner']):
            continue
        faces = expanded_faces(row['inner'], padding)
        select = model.addVars(len(faces), vtype=GRB.BINARY)
        model.addConstr(select.sum() == 1.)
        for f, face in enumerate(faces):
            big_m = 4.-face[2]-np.minimum(face[:2], 0.).sum()
            coefficients = face[:2]/region.bounds
            expression = gp.LinExpr(coefficients, problem.power.tolist())
            expression += float(face[2]-coefficients@region.axis_lower)
            model.addConstr(gamma <= expression+big_m*(1-select[f]))
    model.setObjective(gamma, GRB.MAXIMIZE)
    model.Params.BestObjStop = 1e-6
    with model:
        model.Params.TimeLimit = trace.remaining()
        model.optimize()
        if model.Status == GRB.INFEASIBLE:
            return dict(complete=True, bound=0., x=None, p=None)
        if model.Status == GRB.TIME_LIMIT:
            raise TimeoutError('Physical uncovered search')
        assert model.Status in (GRB.OPTIMAL, GRB.USER_OBJ_LIMIT) and model.MaxVio <= PLANNING_TOL
        bound = float(model.ObjBound)
        return dict(complete=bound <= GEOMETRY_TOL, bound=bound,
                    x=np.rint(problem.x.X).astype(int), p=problem.power.X)


def run_mainline(equations, budget, *, seconds, threads, padding=.005):
    """1 八向初始化；2 顶点评分；3 每 32 次 SP 的物理全局搜索；4 更新；5 结束。"""
    trace = Trace('mainline', seconds)
    region = initialize(equations, budget, threads, trace)
    sp = SubProblem(equations, threads=threads)
    cache, witness = {}, None
    sp_since_global, bound, certified = 0, None, False
    try:
        while True:
            trace.remaining()
            candidates = None
            # 2. 先完成活动见证的同网架支撑，保持主线的局部扩域顺序。
            if witness is not None:
                x, point = witness['x'], region.normalize(witness['p'])
                if region.covered([point], padding)[0]:
                    region.add_point(x, point)
                    witness = None
                else:
                    support = region.witness_support(x, point)
                    support = [q for q in support if not region.covered([q], padding, own=x)[0]]
                    point = max(support, key=lambda q: (float(region.physical(q).sum()), tuple(q)))
                    candidates = [(x, region.inset(x, point, padding, witness))]
            supporting = candidates is not None
            if not supporting:
                candidates = [(row['x'], region.inset(row['x'], point, padding)) for row in region.records.values()
                              for point, covered in zip(row['outer'], region.covered(row['outer'], padding))
                              if not covered]
                # 3. 顶点覆盖不代表并集无孔洞；完整物理 G 自由选择所有 x,p,y。
                if not candidates or sp_since_global >= 32:
                    witness = physical_search(equations, budget, region, padding, threads, trace)
                    trace.counts['global_search'] += 1
                    bound, certified = witness['bound'], witness['complete']
                    if not certified:
                        register(equations, region, witness['x'])
                    trace.event('global_search', region, x=witness['x'], power=witness['p'], bound=bound)
                    if certified:
                        break
                    sp_since_global = 0
                    continue
            scored = []
            # 4. 固定 x,p 求 SP；先加入本批所有证书，再选择最大 eta 的一条割。
            for x, point in candidates:
                power = region.physical(point)
                key = tuple(x), tuple(power)
                if key not in cache:
                    cache[key] = sp.solve(x, power, time_limit=trace.remaining())
                    sp_since_global += 1
                    trace.counts['sp'] = sp.calls
                    trace.event('sp', region, x=x, power=power, eta=cache[key]['eta'])
                checked = cache[key]
                scored.append((x, point, checked))
                if checked['feasible']:
                    region.add_point(x, point)
                    trace.event('certified', region, x=x, power=power)
            pending = [row for row in scored if not row[2]['feasible'] and
                       (supporting or not region.covered([row[1]], padding)[0])]
            if pending:
                x, point, checked = max(pending, key=lambda row: (row[2]['eta'], float(row[1]@region.bounds)))
                region.apply_cut(checked['cut'])
                trace.counts['cuts'] += 1
                trace.event('cut', region, x=x, power=region.physical(point), cut=checked['cut'], eta=checked['eta'])
    except TimeoutError:
        pass  # 时间预算是实验的显式停止条件；不将未完成写成认证完成。
    return finish(region, trace, certified, bound, padding=padding)


def run_kkt(equations, budget, *, seconds, threads, epsilon=.01):
    """1 初始化；2 G 排除认证并集；3 SP 加割；4 固定 p 换网架认证；5 扩域。"""
    trace = Trace('kkt', seconds)
    region = initialize(equations, budget, threads, trace)
    sp = SubProblem(equations, threads=threads)
    supported = set()

    def support(x):
        for direction in DIRECTIONS:
            problem = SignedMaster(equations, budget=budget, direction=direction,
                                   fixed_plan=equations.network.decode_plan(x), threads=threads)
            with problem.model:
                answer = problem.solve(time_limit=trace.remaining())
            region.add_point(x, region.normalize(answer['p']))
            trace.counts['support'] += 1
            trace.event('support', region, x=x, power=answer['p'])
        supported.add(tuple(x))

    oracle = SignedMaster(equations, budget=budget, cuts_only=True, threads=threads)
    eta = oracle.model.addVar(name='violation')
    equations.add_operation(oracle.model, oracle.choices, oracle.loads, eta)
    oracle.model.setObjective(eta, GRB.MINIMIZE)
    certified, bound, padding_bound = False, None, None
    search = None
    try:
        for x in list(region.records):
            support(np.asarray(x))
        search = SignedKkt(equations, budget, region, threads=threads)
        padding_bound = search.padding_bound
        bound = float(search.violation.UB)
        while True:
            # 2. 界足够小才能停；候选的保守下界等于零不跳过其 SP。
            answer = search.solve(trace.remaining(30.), epsilon)
            trace.counts['global_search'] += 1
            bound = max(answer['bound'], padding_bound)
            certified = bound <= epsilon
            trace.event('global_search', region, x=answer['x'], power=answer['p'], bound=bound)
            if certified:
                break
            trace.remaining()
            if answer['x'] is None:
                continue
            x, power = answer['x'], answer['p']
            register(equations, region, x)
            # 3. 原 SP 和原对偶取割，不改变 eta 容差。
            checked = sp.solve(x, power, time_limit=trace.remaining())
            trace.counts['sp'] = sp.calls
            trace.event('sp', region, x=x, power=power, eta=checked['eta'])
            if checked['feasible']:
                feasible, owner = True, x
            else:
                region.apply_cut(checked['cut'])
                search.add_cut(checked['cut'])
                trace.counts['cuts'] += 1
                trace.event('cut', region, x=x, power=power, cut=checked['cut'], eta=checked['eta'])
                # 4. 一个方案失败不删除 p；全部合法网架均参与固定 p 认证。
                found = certify_load(oracle, eta, power, trace.remaining())
                trace.counts['certification'] += 1
                if found['status'] == 'time_limit':
                    raise TimeoutError('Load certification')
                feasible, owner = found['feasible'], found['x']
                trace.event('certification', region, x=owner, power=power, eta=found['eta'])
            if feasible:
                register(equations, region, owner)
                region.add_point(owner, region.normalize(power))
                trace.event('certified', region, x=owner, power=power)
                if tuple(owner) not in supported:
                    support(owner)
                search.model.dispose()
                search = SignedKkt(equations, budget, region, threads=threads)
                search.violation.UB = bound
    except TimeoutError:
        pass
    finally:
        oracle.model.dispose()
        if search is not None:
            search.model.dispose()
    return finish(region, trace, certified, bound, epsilon=epsilon, padding_bound=padding_bound)


def finish(region, trace, certified, bound, **settings):
    trace.event('finish', region, bound=bound)
    result = dict(status='certified' if certified else 'time_limit', certified=certified,
                  seconds=perf_counter()-trace.started, counts=trace.counts, bound=bound,
                  axis_lower=region.axis_lower, axis_bounds=region.axis_bounds, cuts=region.cuts,
                  inner=[dict(vertices=region.physical(row['inner'])) for row in region.records.values()
                         if len(row['inner'])], schemes=len(region.records), events=trace.events, **settings)
    return result, region


def postprocess(equations, budget, result, region, threads):
    """事后投影：枚举仅用于画外域，绝不反馈选点。主线加上已证明的覆盖包络。"""
    from shapely.geometry import Polygon
    from shapely.ops import unary_union
    polytopes = [region.physical(region.outer_for(x)) for x in budget_schemes(equations, budget, threads)]
    outer = unary_union([Polygon(poly[ConvexHull(poly).vertices]) for poly in polytopes
                         if len(poly) >= 3 and np.linalg.matrix_rank(poly-poly[0]) == 2])
    if result['certified'] and 'padding' in result:
        envelopes = []
        for row in region.records.values():
            if len(row['inner']):
                shifts = (2*CORNERS-1)*(result['padding']+GEOMETRY_TOL)
                points = (row['inner'][:, None, :]+shifts[None, :, :]).reshape(-1, 2)
                envelopes.append(Polygon(region.physical(points[ConvexHull(points).vertices])))
        outer = outer.intersection(unary_union(envelopes))
    result['outer'] = polygon_records(outer)


def polygon_records(geometry):
    if geometry.is_empty:
        return []
    if geometry.geom_type == 'Polygon':
        return [dict(vertices=np.asarray(geometry.exterior.coords),
                     holes=[np.asarray(ring.coords) for ring in geometry.interiors])]
    return [row for item in geometry.geoms for row in polygon_records(item)]


def scan(equations, budget, lower, upper, divisions, threads):
    """6. 两方法结束后逐点完整 SOCP；没有正象限单调性或算法区域筛选。"""
    started = perf_counter()
    states = np.empty((divisions, divisions), dtype=np.int8)
    problem = SignedMaster(equations, budget=budget, threads=threads)
    with problem.model:
        problem.model.setObjective(0.)
        problem.model.Params.TimeLimit = 120.
        for index in np.ndindex(states.shape):
            power = lower+(np.asarray(index)+.5)*(upper-lower)/divisions
            problem.power.LB = problem.power.UB = power
            problem.model.optimize()
            assert problem.model.Status in (GRB.OPTIMAL, GRB.INFEASIBLE), (index, problem.model.Status)
            if problem.model.Status == GRB.OPTIMAL:
                assert problem.model.MaxVio <= PLANNING_TOL, (index, problem.model.MaxVio)
            states[index] = 1 if problem.model.Status == GRB.OPTIMAL else -1
            if index[1] == divisions-1 and (index[0]+1) % 5 == 0:
                print(f'SOCP scan {index[0]+1}/{divisions}, {perf_counter()-started:.1f}s', flush=True)
    return dict(axis_lower=lower, axis_bounds=upper, states=states, seconds=perf_counter()-started)


def metrics(result, reference):
    from shapely import intersects_xy
    from shapely.geometry import Polygon
    from shapely.ops import unary_union
    states = np.asarray(reference['states'])
    lower, upper = np.asarray(reference['axis_lower']), np.asarray(reference['axis_bounds'])
    grid = np.indices(states.shape).reshape(2, -1).T
    points = lower+(grid+.5)*(upper-lower)/len(states)
    truth = states.ravel() == 1
    result['metrics'] = {}
    for key in ('inner', 'outer'):
        geometry = unary_union([Polygon(np.asarray(row['vertices'])[ConvexHull(row['vertices']).vertices])
            if key == 'inner' else Polygon(row['vertices'], row.get('holes', []))
            for row in result[key] if len(row['vertices']) >= 3])
        predicted = intersects_xy(geometry, points[:, 0], points[:, 1])
        missed, excess = np.count_nonzero(truth & ~predicted), np.count_nonzero(~truth & predicted)
        result['metrics'][key] = dict(missed=int(missed), excess=int(excess),
            mr_percent=100.*missed/np.count_nonzero(truth),
            fr_percent=100.*excess/max(1, np.count_nonzero(predicted)), area=float(geometry.area))


def save(data, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, 'wt', encoding='utf-8') as stream:
        json.dump(data, stream, ensure_ascii=False, separators=(',', ':'),
                  default=lambda value: value.tolist() if isinstance(value, np.ndarray) else value.item())


def draw_panel(ax, result, reference, load_nodes):
    from matplotlib.colors import ListedColormap
    from shapely.geometry import Polygon
    from shapely.ops import unary_union
    lower, upper = np.asarray(reference['axis_lower']), np.asarray(reference['axis_bounds'])
    states = np.asarray(reference['states'])
    ax.pcolormesh(np.linspace(lower[0], upper[0], len(states)+1),
                 np.linspace(lower[1], upper[1], len(states)+1), (states.T == 1).astype(int),
                 cmap=ListedColormap(['white', '#d6e3ed']), vmin=0, vmax=1, rasterized=True)
    for key, color, linestyle in [('outer', '#b86640', '--'), ('inner', '#247b72', '-')]:
        polygons = []
        for row in result[key]:
            poly = np.asarray(row['vertices'])
            if len(poly) < 3:
                ax.plot(poly[:, 0], poly[:, 1], color=color, marker='.', lw=.7)
                continue
            if 'holes' not in row:
                poly = poly[ConvexHull(poly).vertices]
            polygons.append(Polygon(poly, row.get('holes', [])))
        for row in polygon_records(unary_union(polygons)):
            poly = np.asarray(row['vertices'])
            poly = np.vstack([poly, poly[:1]])
            ax.plot(poly[:, 0], poly[:, 1], color=color, ls=linestyle, lw=1.25)
            for hole in row.get('holes', []):
                hole = np.asarray(hole)
                ax.plot(hole[:, 0], hole[:, 1], color=color, ls=linestyle, lw=1.25)
    ax.axhline(0., color='.65', lw=.5, zorder=1.5)
    ax.axvline(0., color='.65', lw=.5, zorder=1.5)
    margin = .025*(upper-lower)
    ax.set(xlim=(lower[0]-margin[0], upper[0]+margin[0]),
           ylim=(lower[1]-margin[1], upper[1]+margin[1]),
           xlabel=rf'$p_{{{load_nodes[0]}}}$ (kW)', ylabel=rf'$p_{{{load_nodes[1]}}}$ (kW)')


def plot_comparison(data, output):
    """定量网格：每个算法一面板，同扫描、同坐标；PNG 预览及可编辑 SVG/PDF。"""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
                         'font.size': 9, 'svg.fonttype': 'none',
                         'pdf.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False})
    rows = 2 if data['case'] == 'case33' else 1
    fig, axes = plt.subplots(rows, 2, figsize=(7.2, 3.35 if rows == 1 else 6.1), squeeze=False)
    for index, name in enumerate(('mainline', 'kkt')):
        ax = axes[0, index]
        draw_panel(ax, data['results'][name], data['reference'], data['load_nodes'])
        ax.text(-.13, 1.04, chr(97+index), transform=ax.transAxes, weight='bold', fontsize=11)
        ax.set_title('Mainline (signed)' if name == 'mainline' else 'KKT (signed)', fontsize=10)
        if rows == 2:
            ax = axes[1, index]
            draw_panel(ax, data['results'][name], data['reference'], data['load_nodes'])
            upper = np.asarray(data['reference']['axis_bounds'])
            ax.set(xlim=(-2*upper[0], 1.05*upper[0]), ylim=(-.25*upper[1], 1.05*upper[1]))
            ax.text(-.13, 1.04, chr(99+index), transform=ax.transAxes, weight='bold', fontsize=11)
    handles = [Line2D([], [], color='#bacfdf', lw=5, label='SOCP grid reference'),
               Line2D([], [], color='#247b72', label='Certified inner union'),
               Line2D([], [], color='#b86640', ls='--', label='Outer union')]
    fig.legend(handles=handles, loc='lower center', ncol=3, frameon=False, fontsize=8)
    fig.subplots_adjust(left=.12, right=.98, bottom=.23 if rows == 1 else .12,
                        top=.90 if rows == 1 else .94, wspace=.36, hspace=.38)
    fig.savefig(output/'comparison.png', dpi=600)
    fig.savefig(output/'comparison.svg')
    fig.savefig(output/'comparison.pdf')
    plt.close(fig)


def replay(path):
    """两个原生面板同步回放，橙点 SP、红点 G；参考扫描不占回放帧。"""
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Slider
    with gzip.open(path, 'rt', encoding='utf-8') as stream:
        data = json.load(stream)
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    fig.canvas.manager.set_window_title(f'PlanRegion | {data["case"]} | signed power')
    fig.subplots_adjust(bottom=.22, wspace=.28)
    names = ('mainline', 'kkt')
    count = max(len(data['results'][name]['events']) for name in names)
    slider = Slider(fig.add_axes([.15, .07, .70, .04]), 'Frame', 0, count-1, valinit=0, valstep=1)

    def update(value):
        for ax, name in zip(axes, names):
            result = data['results'][name]
            events = result['events']
            frame = min(int(value), len(events)-1)
            event = events[frame]
            records = next((row['records'] for row in reversed(events[:frame+1]) if row['records'] is not None), [])
            current = {key: [dict(vertices=row[key]) for row in records if len(row[key])]
                       for key in ('inner', 'outer')}
            if event['phase'] == 'finish':
                current = {key: result[key] for key in ('inner', 'outer')}
            ax.clear()
            draw_panel(ax, current, data['reference'], data['load_nodes'])
            if event['p'] is not None:
                power = event['p']
                ax.scatter(*power, color='#c53b31' if event['phase'] == 'global_search' else '#ec982e', s=35, zorder=5)
            if event['cut'] is not None:
                cut, x = np.asarray(event['cut']), np.asarray(event['x'])
                lower, upper = np.asarray(data['reference']['axis_lower']), np.asarray(data['reference']['axis_bounds'])
                p1 = np.linspace(lower[0], upper[0], 200)
                if abs(cut[2]) > 1e-14:
                    ax.plot(p1, -(cut[0]+cut[3:]@x+cut[1]*p1)/cut[2], color='#a74726', lw=1.)
                else:
                    ax.axvline(-(cut[0]+cut[3:]@x)/cut[1], color='#a74726', lw=1.)
            ax.set_title(f'{name}: {frame+1}/{len(events)} | {event["phase"]} | {event["seconds"]:.2f}s')
        fig.canvas.draw_idle()
    slider.on_changed(update)
    update(0)
    plt.show()


def run_case(case, *, seconds=1000., epsilon=.01, divisions=80, threads=4, output=ROOT/'results/signed_power'):
    network = FourBus(load_nodes=(1, 2)) if case == 'fourbus' else Case33(load_nodes=(18, 25))
    budget = 20000. if case == 'fourbus' else network.switch_budget
    study_limit = min(network.source_pmax, network.source_smax)*network.base
    equations = SignedPhysics(network, np.full(2, -study_limit), np.full(2, study_limit))
    output = Path(output)/case
    data = dict(case=case, load_nodes=network.load_nodes, budget=budget, power_lower=equations.power_lower,
                power_upper=equations.power_upper, q_ratio=network.q_ratio, voltage_min=np.sqrt(network.vmin),
                voltage_max=np.sqrt(network.vmax), seconds_limit=seconds, threads=threads, results={})
    regions = {}
    with threadpool_limits(limits=1):
        for name, algorithm in [('mainline', run_mainline), ('kkt', run_kkt)]:
            options = dict(epsilon=epsilon) if name == 'kkt' else {}
            result, region = algorithm(equations, budget, seconds=seconds, threads=threads, **options)
            data['results'][name], regions[name] = result, region
            postprocess(equations, budget, result, region, threads)
            save(data, output/'comparison.json.gz')
            print(case, name, result['status'], result['seconds'], result['counts'], flush=True)
        lower = np.min([r.axis_lower for r in regions.values()], axis=0)
        upper = np.max([r.axis_bounds for r in regions.values()], axis=0)
        data['reference'] = scan(equations, budget, lower, upper, divisions, threads)
        for name, result in data['results'].items():
            metrics(result, data['reference'])
            print(case, name, result['metrics'], flush=True)
        save(data, output/'comparison.json.gz')
        plot_comparison(data, output)
    return data


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=('fourbus', 'case33', 'both'), default='both')
    parser.add_argument('--seconds', type=float, default=1000.)
    parser.add_argument('--epsilon', type=float, default=.01)
    parser.add_argument('--divisions', type=int, default=80)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--output', type=Path, default=ROOT/'results/signed_power')
    parser.add_argument('--replay', type=Path)
    args = vars(parser.parse_args())
    replay_path, case = args.pop('replay'), args.pop('case')
    if replay_path:
        replay(replay_path)
    else:
        for name in ('fourbus', 'case33') if case == 'both' else (case,):
            run_case(name, **args)
