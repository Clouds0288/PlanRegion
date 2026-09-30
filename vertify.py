"""独立 AC 潮流、建设方案认证与网格校核；求解失败直接报错。"""
from time import perf_counter
from pathlib import Path
from functools import partial
from concurrent.futures import ProcessPoolExecutor
from contextlib import nullcontext

import gurobipy as gp
from gurobipy import GRB
import numpy as np
from threadpoolctl import threadpool_limits

from model import DEFAULT_SOLVER_THREADS, GridPhysics, MasterProblem, PortPhysics, LOAD_PF, PV_PF, PV_Q_SIGN
from Network.four_bus_five_corridor import FourBus

AC_TOL = 1e-9
FIXED_POINT_TOL = 1e-12
GLOBAL_AC_TOL = 1e-7
AC_TIME_LIMIT = 10.
AC_ITERATIONS = 160


class ACPowerFlow:
    """独立完整 AC；仅读取网架数据，不使用 LP/SOCP 的方程矩阵或乘子。

    适用于非负 P/Q 负荷、正阻抗、根电压固定为 1 p.u. 的径向网络。
    不读取 GridPhysics 的矩阵；以支路递推独立实现 AC 电流等式。
    """

    def __init__(self, network, *, threads=DEFAULT_SOLVER_THREADS):
        self.threads = threads
        self.network = network
        self.model = None

    def state(self, power, ell):
        p, q = self.network.loads(power)
        return self._state(p, q, ell)

    def _state(self, p, q, ell):
        c = self.network
        ell = np.asarray(ell).reshape(-1, c.n)
        P, Q = p+ell*c.r, q+ell*c.reactance
        for i in reversed(c.order):
            if c.parent[i] >= 0:
                P[:, c.parent[i]] += P[:, i]
                Q[:, c.parent[i]] += Q[:, i]
        v, u = np.ones_like(P), np.ones_like(P)
        for i in c.order:
            if c.parent[i] >= 0:
                u[:, i] = v[:, c.parent[i]]
            # 完整支路压降：v=u-2(rP+χQ)+(r²+χ²)*ell；χ 对应 reactance。
            v[:, i] = (u[:, i]-2*(c.r[i]*P[:, i]+c.reactance[i]*Q[:, i])
                       +(c.r[i]**2+c.reactance[i]**2)*ell[:, i])
        return P, Q, v, u

    def violation(self, P, Q, v):
        c = self.network
        ps, qs = P[:, c.roots].sum(axis=1), Q[:, c.roots].sum(axis=1)
        # 所有项目均写为“违反量”：≤0 才满足全部限值，未启用上界产生 -inf。
        return np.maximum.reduce([
            np.max(c.vmin-v, axis=1), np.max(v-c.vmax, axis=1),
            np.max(P-c.capacity, axis=1), ps-c.source_pmax, qs-c.source_qmax,
            np.hypot(ps, qs)-c.source_smax])

    def classify(self, power, return_currents=False):
        """1 可行，-1 已证不可行；迭代未收敛直接报错。"""
        c = self.network
        power = np.asarray(power).reshape(-1, len(c.load_nodes))
        p, q = c.loads(power)
        ell = np.zeros((len(power), c.n))
        status = np.zeros(len(power), dtype=np.int8)
        active = np.arange(len(power))
        for _ in range(AC_ITERATIONS):
            if not len(active):
                break
            P, Q, v, u = self._state(p[active], q[active], ell[active])
            # ell 从零单调递增：功率是下界、电压是上界；仅这些越限能提前拒绝。
            ps, qs = P[:, c.roots].sum(axis=1), Q[:, c.roots].sum(axis=1)
            bad = (np.any(v < c.vmin-AC_TOL, axis=1)
                   | np.any(P > c.capacity+AC_TOL, axis=1)
                   | (ps > c.source_pmax+AC_TOL) | (qs > c.source_qmax+AC_TOL)
                   | (np.hypot(ps, qs) > c.source_smax+AC_TOL) | np.any(u <= 0, axis=1))
            residual = np.max(np.abs(P*P+Q*Q-u*ell[active]), axis=1)
            good = ~bad & (residual <= FIXED_POINT_TOL) & np.all(v <= c.vmax+AC_TOL, axis=1)
            status[active[bad]], status[active[good]] = -1, 1
            keep = ~(bad | good)
            ell[active[keep]] = (P[keep]**2+Q[keep]**2)/u[keep]
            active = active[keep]
        if len(active):
            raise RuntimeError(f'AC power flow did not converge in {AC_ITERATIONS} iterations: p={power[active].tolist()}')
        return (status, ell) if return_currents else status


    def _build_global(self, environment):
        """显式非凸 AC 等式模型，仅供主动调用的独立交叉核验。"""
        c = self.network
        m = gp.Model('independent_AC_reference', env=environment)
        m.Params.OutputFlag = 0
        m.Params.Threads = self.threads
        m.Params.NonConvex = 2
        m.Params.FeasibilityTol = m.Params.OptimalityTol = AC_TOL
        m.Params.DualReductions = 0
        P = m.addVars(c.n, lb=0, ub=c.capacity.tolist())
        Q = m.addVars(c.n, lb=0, ub=min(c.source_qmax, c.source_smax))
        v = m.addVars(c.n, lb=c.vmin.tolist(), ub=c.vmax.tolist())
        bound = min(c.source_smax**2, c.source_pmax**2+c.source_qmax**2)/c.vmin.min()
        ell = m.addVars(c.n, lb=0, ub=bound)
        bp, bq = [], []
        for i in range(c.n):
            up = 1. if c.parent[i] < 0 else v[int(c.parent[i])]
            # P_i-ΣP_child-r_i*ell_i=dP_i；Q_i-ΣQ_child-χ_i*ell_i=dQ_i（标幺）。
            bp.append(m.addConstr(P[i]-gp.quicksum(P[int(j)] for j in c.children[i])-c.r[i]*ell[i] == 0))
            bq.append(m.addConstr(Q[i]-gp.quicksum(Q[int(j)] for j in c.children[i])-c.reactance[i]*ell[i] == 0))
            m.addConstr(v[i] == up-2*(c.r[i]*P[i]+c.reactance[i]*Q[i])
                        +(c.r[i]**2+c.reactance[i]**2)*ell[i])
            m.addQConstr(P[i]*P[i]+Q[i]*Q[i] == up*ell[i])
        ps, qs = (gp.quicksum(values[int(i)] for i in c.roots) for values in (P, Q))
        if np.isfinite(c.source_pmax):
            m.addConstr(ps <= c.source_pmax)
        if np.isfinite(c.source_qmax):
            m.addConstr(qs <= c.source_qmax)
        if np.isfinite(c.source_smax):
            m.addQConstr(ps*ps+qs*qs <= c.source_smax**2)
        m.setObjective(0.)
        m.update()
        self.model = m, ell, bp, bq

    def global_status(self, power, environment, time_limit=AC_TIME_LIMIT):
        if self.model is None:
            self._build_global(environment)
        m, ell, bp, bq = self.model
        m.Params.TimeLimit = max(0.,time_limit)
        p, q = self.network.loads(power)
        m.setAttr('RHS', bp, p[0])
        m.setAttr('RHS', bq, q[0])
        m.optimize()
        if m.Status == GRB.INFEASIBLE:
            return -1
        if m.Status != GRB.OPTIMAL:
            raise RuntimeError(f'Global AC: status={m.Status}, p={list(power)}')
        current = np.array([ell[i].X for i in ell])
        P, Q, v, u = self.state(power, current)
        if (np.max(np.abs(P*P+Q*Q-u*current)) > GLOBAL_AC_TOL
                or self.violation(P, Q, v).max() > GLOBAL_AC_TOL):
            raise RuntimeError(f'Global AC certificate exceeds {GLOBAL_AC_TOL:g}: p={list(power)}')
        return 1

    def close(self):
        if self.model is not None:
            self.model[0].dispose()

def ac_planning_query(equations, power, *, threads=DEFAULT_SOLVER_THREADS):
    """SOCP 搜索建设方案；完整 AC 等式负责认证。"""
    problem = MasterProblem(equations, power=power, threads=threads)
    bound = -np.inf
    with problem.model:
        while True:
            answer = problem.solve()
            if answer is None:
                return None
            answer.pop('state')  # SOCP 状态不能作为独立 AC 证书返回。
            bound = max(bound, answer['bound'])
            answer['bound'] = bound
            oracle = ACPowerFlow(equations.network.tree(answer['x']), threads=threads)
            status = int(oracle.classify(power)[0])
            if status == 1:
                answer.update(feasible=True, status='optimal' if answer['objective']-bound <= 1e-7 else 'feasible')
                return answer
            problem.exclude(answer['x'])


def classify_orthant(states, index, status):
    """调用者须验证非负负荷、正阻抗和 vmax≥根电压；不跨方案构造凸包。"""
    if status == 1:
        states[tuple(slice(0,int(i)+1) for i in index)] = 1
    elif status == -1:
        states[tuple(slice(int(i),None) for i in index)] = -1


def split_grid_box(lower, upper):
    axis = int(np.argmax(upper-lower))
    middle = (lower[axis]+upper[axis])//2
    left, right = upper.copy(), lower.copy()
    left[axis], right[axis] = middle, middle+1
    return ((lower,left),(right,upper))


def validate_ac_region(network, budgets, divisions, bounds, *, threads=DEFAULT_SOLVER_THREADS):
    """连续构域完成后的独立 AC 网格校核；网格不参与规划域构造。"""
    assert (np.all(network.r > 0.) and np.all(network.reactance >= 0.)
            and np.all(network.original_p >= 0.) and np.all(network.original_q >= 0.)
            and np.all(network.q_ratio >= 0.) and np.all(network.vmax >= 1.))
    bounds, budgets = np.asarray(bounds), np.asarray(budgets)
    states = np.zeros((len(budgets),)+(divisions,)*3, dtype=np.int8)
    equations = GridPhysics(network, 'socp')
    visited = set()
    pending = [(j, np.zeros(3, dtype=int), np.full(3, divisions-1, dtype=int))
               for j in reversed(range(len(budgets)))]
    while pending:
        j, lower, upper = pending.pop()
        block = (j,)+tuple(slice(int(a), int(b)+1) for a, b in zip(lower, upper))
        if np.all(states[block] != 0):
            continue
        for index in (upper, lower):
            key = tuple(map(int, index))
            if states[(j,)+key] != 0 or key in visited:
                continue
            power = (index+.5)*bounds/divisions
            answer = ac_planning_query(equations, power, threads=threads)
            visited.add(key)
            labels = np.zeros(len(budgets), dtype=np.int8)
            if answer is None:
                labels[:] = -1
            else:
                labels[budgets < answer['bound']-1e-7] = -1
                labels[budgets >= answer['objective']] = 1
            for k, status in enumerate(labels):
                classify_orthant(states[k], index, status)
            if np.all(states[block] != 0):
                break
        if np.any(states[block] == 0) and np.any(upper > lower):
            pending.extend((j, a, b) for a, b in split_grid_box(lower, upper))
    if np.any(states == 0):
        raise RuntimeError(f'AC grid has {np.count_nonzero(states == 0)} unclassified cells')
    return states


def _scan_line(args):
    """一个扫描行复用完整模型；只改变固定负荷等式的右端。"""
    network, budget, divisions, bounds, index, threads = args
    d = len(network.load_nodes)
    with threadpool_limits(limits=1):
        problem = MasterProblem(GridPhysics(network, 'socp'), power=np.zeros(d), budget=budget, threads=threads)
        states = np.empty(divisions, dtype=np.int8)
        with problem.model as model:
            model.Params.NumericFocus = 3
            model.Params.Presolve = 0 if network.name == 'four_bus_five_corridor' else -1
            model.Params.BarHomogeneous = 1
            model.Params.Aggregate = 0
            model.Params.ScaleFlag = 0
            model.update()
            fixed = [row for row in model.getConstrs() if row.ConstrName.startswith('fixed_power[')]
            for j in range(divisions):
                power = (np.asarray((*index, j))+.5)*bounds/divisions
                model.setAttr('RHS', fixed, power)
                states[j] = 1 if problem.solve() is not None else -1
    return index, states


def validate_socp_region(network, budget, divisions, bounds, *, threads=DEFAULT_SOLVER_THREADS,
                         workers=1, progress=lambda completed, total: None):
    """独立 d 维扫描：逐点固定 p，完整 SOCP 的全部 x、y 自由。"""
    started = perf_counter()
    bounds = np.asarray(bounds, dtype=float)
    d = len(network.load_nodes)
    states = np.empty((divisions,)*d, dtype=np.int8)
    total = divisions**d
    progress(0, total)
    jobs = ((network, budget, divisions, bounds, index, 1 if workers > 1 else threads)
            for index in np.ndindex((divisions,)*(d-1)))
    with ProcessPoolExecutor(max_workers=workers) if workers > 1 else nullcontext() as pool:
        lines = pool.map(_scan_line, jobs) if workers > 1 else map(_scan_line, jobs)
        for completed, (index, line) in enumerate(lines, 1):
            states[index] = line
            progress(completed*divisions, total)
    return dict(bounds=bounds, states=states, scan_seconds=perf_counter()-started)


def scan_line(args, *, network_type=FourBus, budget=20000.):
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


def scan_path(network, budget, output=Path(__file__).resolve().parent/'results'/'scans', mode=1):
    """按算例、模式、节点、功率因数和预算隔离扫描。"""
    case = type(network).__name__.lower()
    return Path(output)/f'{case}_signed'/f'mode_{mode}'/'_'.join(map(str, network.load_nodes))/(
        f'pf_{LOAD_PF:g}_{PV_PF:g}_{PV_Q_SIGN:+g}_budget_{budget:g}.npz')


def save_scan(path, reference):
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, axis_lower=reference['axis_lower'], bounds=reference['bounds'],
                        states=np.asarray(reference['states'], dtype=np.int8))


def scan_reference(network, budget, divisions, bounds, path, *, axis_lower, mode,
                   threads, workers, progress, force_rescan=False):
    """复用包含请求范围的原网格；扩界或强制要求时独立重扫。"""
    lower, upper = np.floor(axis_lower), np.ceil(bounds)
    if path.exists():
        with np.load(path) as saved:
            reference = dict(axis_lower=saved['axis_lower'], bounds=saved['bounds'], states=saved['states'])
        if not force_rescan and np.all(lower >= reference['axis_lower']) and np.all(upper <= reference['bounds']):
            print(f'Reuse scan: {path}', flush=True)
            return reference
        lower, upper = np.minimum(lower, reference['axis_lower']), np.maximum(upper, reference['bounds'])
    started = perf_counter()
    shape = (divisions,)*len(network.load_nodes)
    states = np.empty(shape, dtype=np.int8)
    jobs = ((network.load_nodes, mode, lower, upper, divisions, index) for index in np.ndindex(shape[:-1]))
    last_print = started
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for count, (index, line) in enumerate(pool.map(partial(scan_line, network_type=type(network), budget=budget), jobs), 1):
            states[index] = line
            progress(count*divisions, states.size)
            if perf_counter()-last_print > 15.:
                print(f'Scan {count*divisions}/{states.size}, {perf_counter()-started:.1f}s', flush=True)
                last_print = perf_counter()
    reference = dict(axis_lower=lower, bounds=upper, states=states)
    save_scan(path, reference)
    print(f'Scan completed: {perf_counter()-started:.2f}s, {path}', flush=True)
    return reference
