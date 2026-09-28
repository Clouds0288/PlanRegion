"""二维 FourBus：原双线性 G 与负荷 LP 互补重写的独立对照。

python experiments/fourbus_outer_kkt.py
python experiments/fourbus_outer_kkt.py --plot-only
python experiments/fourbus_outer_kkt.py --coverage --epsilon .01 --seconds 100 --repeats 1

默认两方法保持整数网架、原始 SOCP SP 和 R_k<=epsilon 的停止证书。
--coverage 对比原 KKT 与认证并集排除版本，后者使用负荷并集的违反量上界。
主线、物理约束与取割实现均不改动。
正式对照采用 4 线程；互补乘子使用解析有效界，数值终止状态不接收。
"""
import argparse
from fractions import Fraction
import gzip
from itertools import combinations
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
from experiments.fourbus_outer import GlobalViolation
from experiments.fourbus_outer_scan import point_query
from experiments.fourbus_outer_volume import budget_schemes
from Network.four_bus_five_corridor import FourBus
from model import GridPhysics, MasterProblem, PLANNING_TOL, SubProblem
from region import clip_polytope, contains, halfspaces, initial_polytope, polytope_vertices

OUTPUT = ROOT/'results/fourbus_outer_kkt'
COVERAGE_PAD = 1e-5


class KktViolation(GlobalViolation):
    """固定 x/原 SP 对偶后的负荷 LP，以完整 KKT 替换 p*乘子。"""

    def __init__(self, equations, budget, axis_bounds, total_bound, *, threads=4):
        super().__init__(equations, budget, axis_bounds, total_bound, threads=threads)
        model, template = self.model, self.template
        model.update()
        model.remove([row for row in model.getQConstrs() if row.QCName == 'dual_objective'])
        model.ModelName = 'G_load_LP_KKT'
        model.Params.NonConvex = 0
        dimension = len(axis_bounds)
        size = equations.network.n_types

        # 1. 负荷 LP：outer_p @ p <= outer_rhs + outer_x @ x。
        # 原 MasterProblem 已含原始可行性；这里新增对偶与互补条件。
        self.outer_dual, self.outer_active = [], []
        self.outer_p, self.outer_x = [], []
        self.selector_terms = []
        theta = template.sp_p.T @ self.row_dual
        self.load_stationarity = [model.addConstr(-theta[j].item() == 0.,
                                   name=f'load_stationarity[{j}]') for j in range(dimension)]
        self.outer_choice_dual = model.addMVar(size, lb=-GRB.INFINITY, name='outer_choice_dual')
        self.choice_balance = [model.addConstr(self.outer_choice_dual[j].item() == 0.,
                               name=f'outer_choice_balance[{j}]') for j in range(size)]

        # 2. x*(outer_x.T@lambda) 用二元指示约束表达；不设乘子大 M。
        self.outer_choice_term = model.addMVar(size, lb=-GRB.INFINITY, name='outer_choice_term')
        for j, choice in enumerate(self.problem.choices.values()):
            term = self.outer_choice_term[j].item()
            model.addGenConstrIndicator(choice, 0, term == 0.)
            model.addGenConstrIndicator(choice, 1, term == self.outer_choice_dual[j].item())
        dual_value = (self.choice_term.sum()+self.outer_choice_term.sum()
                      -template.sp_rhs@self.row_dual-template.cone_constant@self.cone_dual).item()
        self.dual_objective = model.addConstr(self.violation <= dual_value, name='kkt_objective')

        # 3. 完整登记负荷上下界、原网络总量界及 MP2 总量界。
        for j, row in enumerate(np.eye(dimension)):
            self._add_outer_row(-row, 0., np.zeros(size))
            self._add_outer_row(row, float(axis_bounds[j]), np.zeros(size))
        self._add_outer_row(np.ones(dimension)/equations.network.base,
                            equations.network.power_limit/equations.network.base, np.zeros(size))
        self._add_outer_row(np.ones(dimension), total_bound, np.zeros(size))
        self._bound_outer_duals()
        model.update()

    def _add_outer_row(self, outer_p, outer_rhs, outer_x, *, selector=None, selector_rhs=0.):
        """新增一行的 lambda/互补选择，同步驻点与替换目标。"""
        model = self.model
        index = len(self.outer_dual)
        dual = model.addVar(name=f'outer_dual[{index}]')
        active = model.addVar(vtype=GRB.BINARY, name=f'outer_active[{index}]')
        self.outer_dual.append(dual)
        self.outer_active.append(active)
        self.outer_p.append(np.asarray(outer_p))
        self.outer_x.append(np.asarray(outer_x))
        slack = float(outer_rhs)+gp.LinExpr(outer_x, self.problem.x.tolist())
        slack -= gp.LinExpr(outer_p, self.problem.power.tolist())
        if selector is not None:
            slack += selector_rhs*selector
            term = model.addVar(name=f'selector_term[{index}]')
            model.addGenConstrIndicator(selector, 0, term == 0.)
            model.addGenConstrIndicator(selector, 1, term == dual)
            model.chgCoeff(self.dual_objective, term, -float(selector_rhs))
            self.selector_terms.append((index, term))
        model.addGenConstrIndicator(active, 0, dual == 0.)
        model.addGenConstrIndicator(active, 1, slack == 0.)
        for row, coefficient in zip(self.load_stationarity, outer_p):
            model.chgCoeff(row, dual, float(coefficient))
        for row, coefficient in zip(self.choice_balance, outer_x):
            model.chgCoeff(row, dual, -float(coefficient))
        model.chgCoeff(self.dual_objective, dual, -float(outer_rhs))

    def _bound_outer_duals(self):
        """二维 LP 有最优基本对偶解；遍历二行基给 lambda 推导有限上界。"""
        template = self.template
        relaxed = template.sp_eta < 0.
        assert template.sp_p.shape[1] == 2 and not np.any(template.sp_p[~relaxed])
        theta = sorted({tuple(Fraction(float(value))/Fraction(float(-eta)) for value in row)
                        for row, eta in zip(template.sp_p[relaxed], template.sp_eta[relaxed])})
        rows = [tuple(Fraction(float(value)) for value in row) for row in self.outer_p]
        upper = [Fraction(0) for _ in rows]
        for i, j in combinations(range(len(rows)), 2):
            a, b = rows[i]
            c, d = rows[j]
            determinant = a*d-b*c
            if determinant == 0:
                continue
            for first, second in theta:
                upper[i] = max(upper[i], (d*first-c*second)/determinant)
                upper[j] = max(upper[j], (a*second-b*first)/determinant)
        limits = np.nextafter(np.array([float(value) for value in upper]), np.inf)
        for dual, value in zip(self.outer_dual, limits):
            dual.UB = value
        for index, term in self.selector_terms:
            term.UB = limits[index]
        coefficients = np.asarray(self.outer_x).T
        lower = [sum(min(Fraction(float(c)), 0)*u for c, u in zip(row, upper)) for row in coefficients]
        higher = [sum(max(Fraction(float(c)), 0)*u for c, u in zip(row, upper)) for row in coefficients]
        lower = np.nextafter(np.array([float(value) for value in lower]), -np.inf)
        higher = np.nextafter(np.array([float(value) for value in higher]), np.inf)
        self.outer_choice_dual.LB, self.outer_choice_dual.UB = lower, higher
        self.outer_choice_term.LB, self.outer_choice_term.UB = np.minimum(0., lower), np.maximum(0., higher)

    def add_cut(self, cut):
        """alpha+beta@p+delta@x>=0 => -beta@p<=alpha+delta@x。"""
        self.problem.add_cut(cut)
        dimension = len(self.problem.loads)
        self._add_outer_row(-cut[1:1+dimension], cut[0], cut[1+dimension:])
        self._bound_outer_duals()


def exclusion_halfspaces(points):
    """认证凸包与 δ 方盒的 Minkowski 和；不把扩边点登记为可行。"""
    shifts = COVERAGE_PAD*np.array([[-1., -1.], [-1., 1.], [1., -1.], [1., 1.]])
    expanded = (np.asarray(points)[:, None, :]+shifts[None, :, :]).reshape(-1, 2)
    return ConvexHull(expanded).equations


class UncoveredKktViolation(KktViolation):
    """在认证并集的补集上求 G；固定面选择后，对完整负荷 LP 写 KKT。"""

    def __init__(self, equations, budget, axis_bounds, total_bound, inner, cuts=(), *, threads=4):
        super().__init__(equations, budget, axis_bounds, total_bound, threads=threads)
        model = self.model
        model.ModelName = 'G_uncovered_load_LP_KKT'
        size = equations.network.n_types
        for cut in cuts:
            cut = np.asarray(cut)
            self.problem.add_cut(cut)
            self._add_outer_row(-cut[1:3], cut[0], cut[3:])
        for k, points in enumerate(inner):
            faces = exclusion_halfspaces(points)
            select = model.addVars(len(faces), vtype=GRB.BINARY, name=f'outside[{k}]')
            model.addConstr(select.sum() == 1.)
            for f, face in enumerate(faces):
                # F@xi+g >= 0 是外侧；M 由 xi∈[0,1]^2 的精确下界推出。
                big_m = max(0., -float(face[2]+np.minimum(face[:2], 0.).sum()))
                outer_p = -face[:2]/np.asarray(axis_bounds)
                outer_rhs = float(face[2])+big_m
                model.addConstr(gp.LinExpr(outer_p, self.problem.power.tolist())
                                <= outer_rhs-big_m*select[f])
                self._add_outer_row(outer_p, outer_rhs, np.zeros(size),
                                    selector=select[f], selector_rhs=-big_m)
        self._bound_outer_duals()
        model.update()

    def solve(self, time_limit, epsilon):
        """空剩余域给 0 界；超时只有有效全局上界能认证完成。"""
        model = self.model
        model.Params.TimeLimit = time_limit
        model.Params.BestBdStop = epsilon
        model.optimize()
        if model.Status == GRB.INFEASIBLE:
            return dict(objective=None, bound=0., status=int(model.Status), x=None, p=None)
        if model.Status not in (GRB.OPTIMAL, GRB.TIME_LIMIT, GRB.USER_OBJ_LIMIT):
            raise RuntimeError(f'Uncovered G status={model.Status}')
        if model.SolCount and model.MaxVio > PLANNING_TOL:
            raise RuntimeError(f'Uncovered G MaxVio={model.MaxVio:g}')
        bound = min(float(model.ObjBound), self.violation.UB)
        x = np.rint(self.problem.x.X).astype(int) if model.SolCount else None
        power = self.problem.power.X if model.SolCount else None
        self.violation.UB = bound
        return dict(objective=self.candidate_lower(x, power) if model.SolCount else None,
                    bound=bound, status=int(model.Status), x=x, p=power)


def certify_load(problem, eta, power, time_limit):
    """V(p)=min_(x,y) eta；全部合法网架参与，原物理容差接受可行点。"""
    model = problem.model
    problem.power.LB = problem.power.UB = power
    model.Params.TimeLimit = time_limit
    model.optimize()
    if model.Status == GRB.TIME_LIMIT:
        return dict(status='time_limit', feasible=False, eta=None, bound=float(model.ObjBound), x=None)
    if model.Status != GRB.OPTIMAL or model.MaxVio > PLANNING_TOL:
        raise RuntimeError(f'Load certification status={model.Status}')
    value = float(eta.X)
    feasible = max(0., value)+model.MaxVio <= PLANNING_TOL
    if not feasible and model.ObjBound <= PLANNING_TOL:
        raise RuntimeError(f'Load certification unresolved: eta={value}, bound={model.ObjBound}')
    return dict(status='optimal', feasible=feasible, eta=value, bound=float(model.ObjBound),
                x=np.rint(problem.x.X).astype(int))


def run_coverage(initialization, *, epsilon=.01, seconds=100., time_limit=5., threads=4):
    """1 初始化认证；2 剩余 G；3 SP 加割；4 固定 p 换网架认证；5 更新并集。"""
    started = perf_counter()
    deadline = started+seconds
    equations = GridPhysics(FourBus(load_nodes=(1, 2)), 'socp')
    bounds = np.asarray(initialization['axis_bounds'])
    budget = initialization['budget']
    result = dict(initialization, method='coverage', epsilon=epsilon, cuts=[], trace=[],
                  status='time_limit', certified=False, bound_scope='uncovered_pairs',
                  counts=dict(global_search=0, sp=0, certification=0, support=0))
    records = {}

    def add_point(x, power):
        key = tuple(x)
        fresh = key not in records
        if fresh:
            # 本算例无固定负荷：零潮流、v=1 使原点在每个合法网架下可行。
            assert not np.any(equations.network.fixed_p) and not np.any(equations.network.fixed_q)
            records[key] = np.zeros((1, 2))
        records[key] = polytope_vertices(np.vstack([records[key], np.asarray(power)/bounds]))
        return fresh

    def support(x):
        # 每个新网架只做一次三方向支撑，避免仅有线段导致逐点平移搜索。
        for direction in [*np.eye(2), np.ones(2)]:
            problem = MasterProblem(equations, budget=budget, direction=direction,
                                    fixed_plan=equations.network.decode_plan(x), threads=threads)
            with problem.model:
                answer = problem.solve(time_limit=max(0., deadline-perf_counter()))
            result['counts']['support'] += 1
            add_point(x, answer['p'])

    # 1. 复用 MP2 的可行点，各网架分别补三方向，不跨网架取凸包。
    for seed in initialization['initial']:
        add_point(seed['x'], seed['p'])
    for key in list(records):
        support(np.asarray(key))
    oracle = MasterProblem(equations, budget=budget, cuts_only=True, threads=threads)
    eta = oracle.model.addVar(name='violation')
    equations.add_operation(oracle.model, oracle.choices, oracle.loads, eta)
    oracle.model.setObjective(eta, GRB.MINIMIZE)
    sp = SubProblem(equations, threads=threads)
    search = UncoveredKktViolation(equations, budget, bounds, initialization['total_bound'],
                                  list(records.values()), threads=threads)
    template = search.template
    relaxed = template.sp_eta < 0.
    assert not np.any(template.sp_p[~relaxed])
    padding_bound = float(COVERAGE_PAD*np.max(np.abs(template.sp_p[relaxed])@bounds
                                             /(-template.sp_eta[relaxed]))+PLANNING_TOL)
    result.update(bound=float(search.violation.UB), padding_bound=padding_bound,
                  union_bound=float(search.violation.UB))
    with oracle.model:
        while perf_counter() < deadline:
            # 2. G 排除所有已认证网架的并集，而不只是候选 x 自身的内域。
            answer = search.solve(min(time_limit, deadline-perf_counter()), epsilon)
            result['counts']['global_search'] += 1
            result.update(bound=answer['bound'], union_bound=max(answer['bound'], padding_bound))
            record = dict(step=result['counts']['global_search'], objective=answer['objective'],
                          bound=answer['bound'], status=answer['status'], schemes=len(records))
            result['trace'].append(record)
            print(f"coverage G{record['step']}: lower={answer['objective']} upper={answer['bound']:.8g}, "
                  f"schemes={len(records)}", flush=True)
            if result['union_bound'] <= epsilon:
                result.update(status='union_residual_certified', certified=True)
                break
            if answer['objective'] is None or answer['objective'] <= PLANNING_TOL:
                result['status'] = 'global_search_unresolved'
                break
            if perf_counter() >= deadline:
                break
            x, power = answer['x'], answer['p']
            assert not any(contains(power/bounds, halfspaces(points))[0] for points in records.values())
            record.update(x=x.tolist(), p=power.tolist())
            # 3. 固定 G 的搭配检查原 SP；联合割照常适用于所有合法 x。
            checked = sp.solve(x, power, time_limit=min(30., deadline-perf_counter()))
            result['counts']['sp'] += 1
            assert not checked['feasible']
            assert answer['objective'] <= checked['eta']+2e-7 <= answer['bound']+4e-7
            cut = checked['cut']
            result['cuts'].append(cut.tolist())
            record.update(eta=checked['eta'])
            search.add_cut(cut)
            if perf_counter() >= deadline:
                break
            # 4. 固定同一负荷，x/y 全部放开；不能以刚才的 x 不可行删除 p。
            certified = certify_load(oracle, eta, power, min(30., deadline-perf_counter()))
            result['counts']['certification'] += 1
            record.update(certification_eta=certified['eta'], certification_bound=certified['bound'])
            if certified['status'] == 'time_limit':
                result['status'] = 'certification_time_limit'
                break
            # 5. 成功认证扩大所属网架凸包；失败仅保留已加的有效割。
            if certified['feasible']:
                record['certified_x'] = certified['x'].tolist()
                fresh = add_point(certified['x'], power)
                if fresh:
                    support(certified['x'])
                search.model.dispose()
                search = UncoveredKktViolation(equations, budget, bounds, initialization['total_bound'],
                                              list(records.values()), result['cuts'], threads=threads)
                search.violation.UB = result['bound']
    search.model.dispose()
    result.update(inner=[dict(x=np.asarray(key, dtype=int).tolist(), vertices=(points*bounds).tolist())
                         for key, points in records.items()],
                  seconds=perf_counter()-started,
                  total_seconds=perf_counter()-started+initialization['initial_seconds'])
    print(f"coverage: {result['status']}, cuts={len(result['cuts'])}, "
          f"total={result['total_seconds']:.4f}s", flush=True)
    return result


def initialize(budget, threads):
    """共用两个轴方向及总负荷方向的完整 SOCP MP2。"""
    started = perf_counter()
    equations = GridPhysics(FourBus(load_nodes=(1, 2)), 'socp')
    initial = []
    for direction in [*np.eye(2), np.ones(2)]:
        problem = MasterProblem(equations, budget=budget, direction=direction, threads=threads)
        with problem.model:
            answer = problem.solve(time_limit=60.)
        initial.append(dict(direction=direction.tolist(), x=answer['x'].tolist(),
                            p=answer['p'].tolist(), bound=float(answer['bound'])))
    return dict(budget=budget, initial=initial, initial_seconds=perf_counter()-started,
                axis_bounds=[item['bound'] for item in initial[:-1]], total_bound=initial[-1]['bound'])


def run_trial(method, initialization, *, epsilon=.1, seconds=300., time_limit=5., threads=4):
    """1 初始化；2 全局选点；3 判上界；4 原 SP；5 加割再搜索。"""
    started = perf_counter()
    deadline = started+seconds
    equations = GridPhysics(FourBus(load_nodes=(1, 2)), 'socp')
    result = dict(initialization, method=method, epsilon=epsilon, cuts=[], trace=[],
                  status='time_limit', certified=False, counts=dict(global_search=0, sp=0))
    # 1. 同一初始外域、零条割、相同完整 SP；只改变 G 的表达方式。
    search_type = {'bilinear': GlobalViolation, 'kkt': KktViolation}[method]
    search = search_type(equations, initialization['budget'], initialization['axis_bounds'],
                         initialization['total_bound'], threads=threads)
    sp = SubProblem(equations, threads=threads)
    add_cut = search.add_cut if method == 'kkt' else search.problem.add_cut
    search.model.update()
    result['bound'] = float(search.violation.UB)
    with search.model:
        while perf_counter() < deadline:
            # 2. G 搜索所有预算内 x,p，不使用已知方案列表或扫描结果。
            answer = search.solve(min(time_limit, deadline-perf_counter()), epsilon)
            result['counts']['global_search'] += 1
            result['bound'] = answer['bound']
            record = dict(step=result['counts']['global_search'], objective=answer['objective'],
                          bound=answer['bound'], status=answer['status'], global_seconds=answer['seconds'])
            result['trace'].append(record)
            print(f"{method} G{record['step']}: lower={answer['objective']} upper={answer['bound']:.8g}",
                  flush=True)
            # 3. 只使用全局上界终止；无有效反例时如实保留未认证状态。
            if answer['bound'] <= epsilon:
                result.update(status='residual_certified', certified=True)
                break
            if answer['objective'] is None or answer['objective'] <= PLANNING_TOL:
                result['status'] = 'global_search_unresolved'
                break
            if perf_counter() >= deadline:
                break
            x, power = answer['x'], answer['p']
            record.update(x=x.tolist(), p=power.tolist())
            # 4. 固定相同 x,p，原 SP 仅调 y/eta；其结果不来自 KKT 模型。
            checked = sp.solve(x, power, time_limit=min(30., deadline-perf_counter()))
            result['counts']['sp'] += 1
            assert not checked['feasible']
            assert answer['objective'] <= checked['eta']+2e-7 <= answer['bound']+4e-7
            cut = checked['cut']
            record.update(eta=checked['eta'],
                          cut_at_candidate=float(cut[0]+cut[1:3]@power+cut[3:]@x))
            # 5. 两方法加入同一种有效联合割；KKT 同步更新外域的最优性条件。
            add_cut(cut)
            result['cuts'].append(cut.tolist())
    result.update(seconds=perf_counter()-started,
                  total_seconds=perf_counter()-started+initialization['initial_seconds'])
    print(f"{method}: {result['status']}, cuts={len(result['cuts'])}, "
          f"total={result['total_seconds']:.4f}s", flush=True)
    return result


def reference_scan(initialization, divisions, threads):
    """独立逐点完整 SOCP：固定 p，自由选择全部合法 x/y，不含算法割。"""
    started = perf_counter()
    equations = GridPhysics(FourBus(load_nodes=(1, 2)), 'socp')
    bounds = np.asarray(initialization['axis_bounds'])
    axes = [(np.arange(divisions)+.5)*bound/divisions for bound in bounds]
    power = np.stack(np.meshgrid(*axes, indexing='ij'), axis=-1).reshape(-1, 2)
    states = np.empty(len(power), dtype=np.int8)
    problem = MasterProblem(equations, budget=initialization['budget'], threads=threads)
    problem.model.setObjective(0.)
    with problem.model:
        for i, p in enumerate(power):
            states[i] = 1 if point_query(problem, p) else -1
            if (i+1) % (divisions*10) == 0 or i+1 == len(power):
                print(f'SOCP reference {i+1}/{len(power)}', flush=True)
    return dict(bounds=bounds.tolist(), divisions=divisions,
                states=states.reshape(divisions, divisions).tolist(), seconds=perf_counter()-started)


def compare_reference(result, reference, schemes):
    """外域为全部预算内网架的并集；此处枚举仅用于事后测量和绘图。"""
    from shapely.geometry import MultiPoint
    from shapely.ops import unary_union
    bounds = np.asarray(result['axis_bounds'])
    divisions = reference['divisions']
    axes = [(np.arange(divisions)+.5)*bound/divisions for bound in bounds]
    power = np.stack(np.meshgrid(*axes, indexing='ij'), axis=-1).reshape(-1, 2)
    feasible = np.asarray(reference['states']).ravel() == 1
    initial = (power.sum(axis=1) <= result['total_bound']+PLANNING_TOL)
    retained = np.zeros(len(power), dtype=bool)
    cuts = np.asarray(result['cuts']).reshape(-1, 1+2+schemes.shape[1])
    base = initial_polytope(bounds, result['total_bound'], bounds)
    polytopes = []
    for x in schemes:
        constants = cuts[:, 0]+cuts[:, 3:]@x
        retained |= initial & np.all(constants[:, None]+cuts[:, 1:3]@power.T >= -PLANNING_TOL, axis=0)
        vertices = base.copy()
        for cut in cuts:
            vertices = clip_polytope(vertices, cut[0]+cut[3:]@x, cut[1:3]*bounds)
            if not len(vertices):
                break
        if len(vertices):
            polytopes.append(MultiPoint(vertices*bounds).convex_hull)
    outer = unary_union(polytopes)
    result.update(fr_percent=float(100*np.count_nonzero(retained & ~feasible)/np.count_nonzero(retained)),
                  mr_percent=float(100*np.count_nonzero(~retained & feasible)/np.count_nonzero(feasible)),
                  outer_area=float(outer.area), outer_wkt=outer.wkt)
    assert not np.any(~retained & feasible), 'An outer region excluded an independent SOCP feasible point'


def plot_comparison(data, output):
    """两张同尺度区域图：离散 SOCP 参考、最终外域及扫描多余点。"""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from shapely import contains_xy, from_wkt
    plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
                         'font.size': 8, 'axes.labelsize': 8, 'legend.fontsize': 7,
                         'svg.fonttype': 'none', 'pdf.fonttype': 42,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.linewidth': .7})
    reference = data['reference']
    bounds, n = np.asarray(reference['bounds']), reference['divisions']
    edges = [np.linspace(0., bound, n+1) for bound in bounds]
    centers = np.stack(np.meshgrid(*[(edge[:-1]+edge[1:])/2 for edge in edges], indexing='ij'), axis=-1)
    feasible = np.asarray(reference['states']) == 1
    for letter, method in zip('abc', data['runs']):
        runs = data['runs'][method]
        result = sorted(runs, key=lambda row: row['total_seconds'])[len(runs)//2]
        outer = from_wkt(result['outer_wkt'])
        retained = contains_xy(outer, centers[..., 0], centers[..., 1])
        colors = np.zeros(feasible.shape, dtype=int)
        colors[feasible] = 1
        colors[retained & ~feasible] = 2
        colors[~retained & feasible] = 3
        fig, ax = plt.subplots(figsize=(3.5039370079, 3.1496062992))  # 89 x 80 mm
        ax.pcolormesh(*edges, colors.T, cmap=ListedColormap(['white', '#BCCDD5', '#E7C397', '#B87483']),
                      vmin=0, vmax=3, shading='flat', rasterized=True)
        polygons = [outer] if outer.geom_type == 'Polygon' else list(outer.geoms)
        for polygon in polygons:
            if polygon.geom_type == 'Polygon':
                for ring in [polygon.exterior, *polygon.interiors]:
                    coordinates = np.asarray(ring.coords)
                    ax.plot(coordinates[:, 0], coordinates[:, 1], color='#A35D20', lw=1.15)
        ax.set(xlim=(0, bounds[0]*1.01), ylim=(0, bounds[1]*1.03),
               xlabel=r'Load $p_1$ (kW)', ylabel=r'Load $p_2$ (kW)')
        ax.text(-.15, 1.04, letter, transform=ax.transAxes, fontweight='bold', fontsize=10)
        ax.text(0., 1.04, {'bilinear': 'Bilinear', 'kkt': 'KKT reformulation',
                          'coverage': 'KKT with certified-union exclusion'}[method],
                transform=ax.transAxes, fontsize=8)
        fig.legend(handles=[Patch(facecolor='#BCCDD5', label='SOCP scan'),
                            Patch(facecolor='#E7C397', label='Excess'),
                            Line2D([], [], color='#A35D20', lw=1.15, label='Outer boundary')],
                   loc='lower center', bbox_to_anchor=(.53, .015), ncol=3, frameon=False,
                   handlelength=1.3, columnspacing=.8, handletextpad=.5)
        fig.subplots_adjust(left=.17, right=.98, bottom=.24, top=.91)
        fig.savefig(output/f'{method}_vs_reference.png', dpi=600)
        fig.savefig(output/f'{method}_vs_reference.svg')
        fig.savefig(output/f'{method}_vs_reference.pdf')
        plt.close(fig)


def run_comparison(output=OUTPUT, *, repeats=3, divisions=80, threads=4,
                   epsilon=.1, seconds=300., time_limit=5., coverage=False):
    """同条件交替运行；一份原始结果、共同扫描和两张矢量图。"""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    with threadpool_limits(limits=1):
        initialization = initialize(20000., threads)
        methods = ('kkt', 'coverage') if coverage else ('bilinear', 'kkt')
        data = dict(initialization=initialization,
                    settings=dict(repeats=repeats, divisions=divisions, threads=threads, kkt_dual_bounds=True,
                                  epsilon=epsilon, seconds=seconds, time_limit=time_limit,
                                  coverage=coverage, coverage_pad=COVERAGE_PAD if coverage else None),
                    runs={method: [] for method in methods})
        for repeat in range(repeats):
            order = methods if repeat % 2 == 0 else methods[::-1]
            for method in order:
                if method == 'coverage':
                    result = run_coverage(initialization, epsilon=epsilon, seconds=seconds,
                                          time_limit=time_limit, threads=threads)
                else:
                    result = run_trial(method, initialization, epsilon=epsilon, seconds=seconds,
                                       time_limit=time_limit, threads=threads)
                data['runs'][method].append(result)
                with gzip.open(output/'comparison.json.gz', 'wt', encoding='utf-8') as stream:
                    json.dump(data, stream, ensure_ascii=False)
        data['reference'] = reference_scan(initialization, divisions, threads)
        equations = GridPhysics(FourBus(load_nodes=(1, 2)), 'socp')
        schemes = budget_schemes(equations, initialization['budget'], threads)
        for method, runs in data['runs'].items():
            for result in runs:
                compare_reference(result, data['reference'], schemes)
            print(json.dumps(dict(method=method, times=[row['total_seconds'] for row in runs],
                                  cuts=[len(row['cuts']) for row in runs],
                                  fr_percent=[row['fr_percent'] for row in runs],
                                  mr_percent=[row['mr_percent'] for row in runs])), flush=True)
        with gzip.open(output/'comparison.json.gz', 'wt', encoding='utf-8') as stream:
            json.dump(data, stream, ensure_ascii=False)
    plot_comparison(data, output)
    return data


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUTPUT)
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--divisions', type=int, default=80)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--epsilon', type=float, default=.1)
    parser.add_argument('--seconds', type=float, default=300.)
    parser.add_argument('--time-limit', type=float, default=5.)
    parser.add_argument('--plot-only', action='store_true')
    parser.add_argument('--coverage', action='store_true')
    arguments = vars(parser.parse_args())
    if arguments.pop('plot_only'):
        with gzip.open(arguments['output']/'comparison.json.gz', 'rt', encoding='utf-8') as stream:
            plot_comparison(json.load(stream), arguments['output'])
    else:
        run_comparison(**arguments)
