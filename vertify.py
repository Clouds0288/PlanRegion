"""独立 AC 校核：全部合法拓扑、固定格架缓存、增量补扫与结果对比。"""
import argparse
from contextlib import contextmanager, nullcontext
from concurrent.futures import ProcessPoolExecutor, as_completed
from copy import deepcopy
from dataclasses import asdict
from functools import partial
from itertools import product
import hashlib
import json
import os
from pathlib import Path
import shutil
from time import perf_counter
from uuid import uuid4

import gurobipy as gp
from gurobipy import GRB
import numpy as np
from threadpoolctl import threadpool_limits

from model import DEFAULT_SOLVER_THREADS, GridPhysics, MasterProblem, PortPhysics
from model import PLANNING_TOL, LOAD_PF, PV_PF, PV_Q_SIGN

AC_TOL = 1e-9
FIXED_POINT_TOL = 1e-12
AC_ITERATIONS = 160


class ACPowerFlow:
    """支路递推的完整 AC 潮流状态；仅读取运行树数据，不使用 LP/SOCP 的方程矩阵或乘子。
    state(power, ell) 由节点负荷与支路电流平方递推 (P,Q,v,u)，供扫描复核 AC 见证与运行界。"""

    def __init__(self, network):
        self.network = network

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


def budget_schemes(equations, budget, threads=4):
    """只检查建设预算、锁定线路、道路与根向树；不以零负荷筛除合法拓扑。"""
    net = equations.network
    options = []
    for corridor, allowed in zip(net.corridors, net.road_allowed):
        if not corridor.switchable:
            selected = net.initial_plan[corridor.id]
            if not allowed and selected is not None:
                return np.empty((0, net.n_types), dtype=int)
            options.append((selected,))
        else:
            options.append((None, *(kind.id for kind in corridor.types)) if allowed else (None,))
    schemes = []
    for selected in product(*options):
        x = net.encode_plan(dict(zip((c.id for c in net.corridors), selected)))
        if net.cost_offset+net.cost@x > budget:
            continue
        try:
            net.tree(x)
        except ValueError:
            continue  # 该组合不是覆盖必接节点的根向树。
        schemes.append(x)
    return np.asarray(schemes, dtype=int).reshape(-1, net.n_types)


def signed_ac_witness(network, x, power, *, mode=1):
    """独立树递推只提供可行证书；不使用正负荷单调性拒绝反送功率点。"""
    tree = network.tree(x)
    oracle = ACPowerFlow(tree)
    power = np.asarray(power, dtype=float).reshape(-1, len(network.load_nodes))
    ratio = np.where(power >= 0., np.tan(np.arccos(LOAD_PF)),
                     PV_Q_SIGN*np.tan(np.arccos(PV_PF))) if mode else network.q_ratio
    p = (tree.fixed_p+power@tree.E.T)/tree.base
    q = (tree.fixed_q+(power*ratio)@tree.E.T)/tree.base
    ell = np.zeros((len(power), tree.n))
    feasible = np.zeros(len(power), dtype=bool)
    residual = np.full(len(power), np.inf)
    active = np.arange(len(power))
    for _ in range(AC_ITERATIONS):
        if not len(active):
            break
        P, Q, v, u = oracle._state(p[active], q[active], ell[active])
        residual[active] = np.max(np.abs(P*P+Q*Q-u*ell[active]), axis=1)
        converged = residual[active] <= FIXED_POINT_TOL
        ps, qs = P[:, tree.roots].sum(axis=1), Q[:, tree.roots].sum(axis=1)
        valid = (np.all(v >= tree.vmin-AC_TOL, axis=1)
                 & np.all(v <= tree.vmax+AC_TOL, axis=1)
                 & np.all(ell[active] <= tree.ell_limit+AC_TOL, axis=1)
                 & np.all(P <= tree.capacity+AC_TOL, axis=1)
                 & np.all(-P+tree.r*ell[active] <= tree.capacity+AC_TOL, axis=1)
                 & (ps <= tree.source_pmax+AC_TOL) & (qs <= tree.source_qmax+AC_TOL)
                 & (np.hypot(ps, qs) <= tree.source_smax+AC_TOL))
        feasible[active[converged & valid]] = True
        # 无正电压或迭代越界仅放弃该见证；后续完整等式模型负责未决点。
        keep = ~converged & np.all(u > 0., axis=1) & np.all(v > 0., axis=1)
        ell[active[keep]] = (P[keep]**2+Q[keep]**2)/u[keep]
        active = active[keep]
    return dict(feasible=feasible, ell=ell, residual=residual)


def ac_interval_possible(network, x, power, *, mode=1):
    """外扩的 AC 必要区间；False 为本树不可行，True 不构成可行证书。"""
    tree = network.tree(x)
    oracle = ACPowerFlow(tree)
    power = np.asarray(power, dtype=float).reshape(-1, len(network.load_nodes))
    ratio = np.where(power >= 0., np.tan(np.arccos(LOAD_PF)),
                     PV_Q_SIGN*np.tan(np.arccos(PV_PF))) if mode else network.q_ratio
    p = (tree.fixed_p+power@tree.E.T)/tree.base
    q = (tree.fixed_q+(power*ratio)@tree.E.T)/tree.base
    lower = np.zeros((len(power), tree.n))
    # Kirchhoff 电流定律和三角不等式，使用实际节点电压下限。
    upper = (np.hypot(p, q)/np.sqrt(tree.vmin-AC_TOL)@tree.D.T)**2+AC_TOL
    upper = np.minimum(upper, tree.ell_limit+AC_TOL)
    possible = np.ones(len(power), dtype=bool)
    active = np.arange(len(power))
    for _ in range(32):
        if not len(active):
            break
        Plo, Qlo, vhi, uhi = oracle._state(p[active], q[active], lower[active])
        Phi, Qhi, vlo, ulo = oracle._state(p[active], q[active], upper[active])
        keep = (np.all(vhi >= tree.vmin-AC_TOL, axis=1)
                & np.all(vlo <= tree.vmax+AC_TOL, axis=1)
                & np.all(lower[active] <= upper[active]+AC_TOL, axis=1)
                & (Plo[:, tree.roots].sum(axis=1) <= tree.source_pmax+AC_TOL)
                & (Qlo[:, tree.roots].sum(axis=1) <= tree.source_qmax+AC_TOL))
        possible[active[~keep]] = False
        active = active[keep]
        if not len(active):
            break
        Plo, Qlo, Phi, Qhi = (a[keep] for a in (Plo, Qlo, Phi, Qhi))
        ulo = np.maximum(ulo[keep], np.r_[tree.vmin, 1.][tree.parent]-AC_TOL)
        uhi = np.minimum(uhi[keep], np.r_[tree.vmax, 1.][tree.parent]+AC_TOL)
        pmin = np.where((Plo <= 0.) & (Phi >= 0.), 0., np.minimum(Plo**2, Phi**2))
        qmin = np.where((Qlo <= 0.) & (Qhi >= 0.), 0., np.minimum(Qlo**2, Qhi**2))
        lower[active] = np.maximum(lower[active], (pmin+qmin)/uhi-AC_TOL)
        upper[active] = np.minimum(upper[active],
            (np.maximum(Plo**2, Phi**2)+np.maximum(Qlo**2, Qhi**2))/ulo+AC_TOL)
    return possible


AC_CACHE_METHOD = 'ac_socp_grid_v4'
SCAN_OUTPUT = Path(__file__).resolve().parent/'results'/'scan'
SCAN_FIELDS = ('states', 'witness_x', 'residual', 'socp_states', 'socp_witness_x', 'socp_residual')
SAVE_SECONDS = 30.   # 扫描中两次落盘的最短间隔：避免紧接着替换刚写入的缓存（Windows 上该文件可能仍被扫描占用）
SCAN_TOL = 1e-6      # 参考扫描接受可行解的质量门槛（Gurobi 默认可行性容差）；不可行判定仍用 PLANNING_TOL。
                     # 只决定原先给不出标签的点，已有标签不变，故不进入缓存身份


def ac_network(network):
    """冻结构域前的网络数据，AC 与 SOCP 保留相同逐线路限额。"""
    return deepcopy(network)


def reference_box(network, budget=None, *, mode=1, output=SCAN_OUTPUT):
    """同限流 SOCP 的坐标全局上界包含 AC 域，供两套扫描共用。"""
    network = ac_network(network)
    budget = network.budgets[0] if budget is None else budget
    path = scan_path(network, budget, output, mode)/'bounds.npz'
    identity = ac_identity(network, budget, mode)
    with _scan_lock(path):
        if path.exists():
            with np.load(path) as saved:
                if (str(saved['identity']) == identity and str(saved['method']) == 'socp_global_bound_v1'
                        and str(saved.get('power_unit', '')) == 'kW'):
                    return saved['axis_lower'], saved['bounds']
            raise ValueError(f'AC bound identity differs: {path}')
        d = len(network.load_nodes)
        lower, upper, supports = np.zeros(d), np.zeros(d), []
        with threadpool_limits(limits=1):
            for sign in product((-1, 1), repeat=d) if mode else [(1,)*d]:
                for axis in range(d):
                    problem = scan_problem(network, budget, sign, mode=mode, direction=np.eye(d)[axis])
                    with problem.model as model:
                        model.Params.TimeLimit = 60.
                        model.Params.Presolve = 2
                        model.optimize()
                        if model.Status == GRB.INFEASIBLE:
                            continue
                        if model.Status not in (GRB.OPTIMAL, GRB.NODE_LIMIT, GRB.TIME_LIMIT) or not np.isfinite(model.ObjBound):
                            raise RuntimeError(f'No finite SOCP coordinate bound: sign={sign}, axis={axis}, status={model.Status}')
                        bound = float(model.ObjBound)*problem.objective_scale
                        supports.append(dict(sign=sign, axis=axis, bound_kw=bound, status=int(model.Status)))
                        if sign[axis] > 0:
                            upper[axis] = max(upper[axis], bound)
                        else:
                            lower[axis] = min(lower[axis], -bound)
        lower, upper = np.floor(lower-1e-4) if mode else lower, np.ceil(upper+1e-4)
        if np.any(upper <= lower):
            raise ValueError('Empty AC reference bounds')
        temporary = path.with_suffix('.tmp')
        with temporary.open('wb') as stream:
            np.savez_compressed(stream, method='socp_global_bound_v1', identity=identity, power_unit='kW',
                axis_lower=lower, bounds=upper, supports=json.dumps(supports))
        temporary.replace(path)
        return lower, upper


def ac_identity(network, budget, mode=1):
    network = ac_network(network)
    payload = asdict(network)
    # 正负模式逐点重建功率因数，构域过程中临时写入的这两个值不是 AC 条件。
    if mode:
        payload.pop('q_ratio')
        payload.pop('power_limit')
    payload.update(method=AC_CACHE_METHOD, budget=float(budget), mode=int(mode),
                   power_factors=[LOAD_PF, PV_PF, PV_Q_SIGN],
                   tolerances=[AC_TOL, FIXED_POINT_TOL, PLANNING_TOL],
                   operating_arrays={key: getattr(network, key) for key in
                                     ('r', 'reactance', 'capacity', 'ell_limit', 'cost')})
    encoded = json.dumps(payload, sort_keys=True, default=lambda a: a.tolist(), separators=(',', ':'))
    return hashlib.sha256(encoded.encode()).hexdigest()


def scan_path(network, budget, output=SCAN_OUTPUT, mode=1):
    return Path(output)/network.name/'_'.join(map(str, network.load_nodes))/ac_identity(network, budget, mode)


def region_path(path, reference):
    grid = np.concatenate([np.asarray(reference[key], float) for key in ('origin', 'step', 'start')]
                          + [np.array(reference['states'].shape, float)])
    return Path(path)/('region_'+hashlib.sha256(grid.tobytes()).hexdigest()[:16]+'.npz')


@contextmanager
def _scan_lock(path):
    """OS 文件锁随进程退出释放，禁止两个扫描覆盖同一份缓存。"""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.with_suffix('.lock').open('a+b') as stream:
        stream.seek(0)
        if not stream.read(1):
            stream.write(b'0')
            stream.flush()
        stream.seek(0)
        try:
            if os.name == 'nt':
                import msvcrt
                msvcrt.locking(stream.fileno(), msvcrt.LK_NBLCK, 1)
            else:
                import fcntl
                fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as exc:
            raise RuntimeError(f'AC cache is in use: {path}') from exc
        try:
            yield
        finally:
            stream.seek(0)
            if os.name == 'nt':
                msvcrt.locking(stream.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                fcntl.flock(stream.fileno(), fcntl.LOCK_UN)


def save_scan(path, reference):
    """一个功率区域同时保存 AC、SOCP 标签与证书。"""
    if reference.get('method') != AC_CACHE_METHOD:
        raise ValueError('Only identified AC references can enter the AC cache')
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name+'.'+uuid4().hex+'.tmp')
    fields = ('axis_lower', 'bounds', 'origin', 'step', 'start')+SCAN_FIELDS
    try:
        with temporary.open('wb') as stream:
            np.savez_compressed(stream, **{key: reference[key] for key in fields},
                method=AC_CACHE_METHOD, metadata=json.dumps(reference['metadata'], ensure_ascii=False))
            stream.flush()
            os.fsync(stream.fileno())
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


class ACReferenceMismatch(ValueError):
    """参考来源或物理配置不同；区别于同一配置的数据损坏或证书冲突。"""


def load_scan(path, identity=None):
    with np.load(path, allow_pickle=False) as saved:
        if str(saved.get('method', '')) != AC_CACHE_METHOD or 'metadata' not in saved:
            raise ACReferenceMismatch(f'Not an identified AC cache: {path}')
        answer = {key: saved[key].copy() for key in
                  ('axis_lower', 'bounds', 'origin', 'step', 'start')+SCAN_FIELDS}
        answer.update(method=AC_CACHE_METHOD, metadata=json.loads(str(saved['metadata'])))
    if identity is not None and answer['metadata']['identity'] != identity:
        raise ACReferenceMismatch(f'AC cache physical identity differs: {path}')
    shape = np.array(answer['states'].shape)
    d = len(shape)
    if (d not in (2, 3) or np.any(shape <= 0) or not np.isin(answer['states'], [-1, 0, 1]).all()
            or any(answer[key].shape != (d,) for key in ('axis_lower', 'bounds', 'origin', 'step', 'start'))
            or not np.isfinite(answer['step']).all() or np.any(answer['step'] <= 0)
            or not np.array_equal(answer['start'], np.rint(answer['start']))
            or answer['residual'].shape != tuple(shape)
            or answer['witness_x'].shape != (*shape, answer['metadata']['n_types'])):
        raise ValueError(f'Malformed AC grid: {path}')
    if (not np.allclose(answer['axis_lower'], answer['origin']+answer['start']*answer['step'], atol=1e-8, rtol=0)
            or not np.allclose(answer['bounds'], answer['axis_lower']+shape*answer['step'], atol=1e-8, rtol=0)):
        raise ValueError(f'AC grid coordinates disagree: {path}')
    return answer


def _empty_ac_grid(network, budget, mode, origin, step, start, shape):
    lower = origin+start*step
    return dict(method=AC_CACHE_METHOD, axis_lower=lower, bounds=lower+np.array(shape)*step,
        origin=origin, step=step, start=start, states=np.zeros(shape, dtype=np.int8),
        witness_x=np.zeros((*shape, network.n_types), dtype=np.int8), residual=np.full(shape, np.nan),
        socp_states=np.zeros(shape, dtype=np.int8),
        socp_witness_x=np.zeros((*shape, network.n_types), dtype=np.int8), socp_residual=np.full(shape, np.nan),
        metadata=dict(identity=ac_identity(network, budget, mode), network=network.name,
            network_fingerprint=network.fingerprint, load_nodes=list(network.load_nodes),
            budget=float(budget), mode=int(mode), n_types=network.n_types,
            power_factors=[LOAD_PF, PV_PF, PV_Q_SIGN],
            current_limit_a=(np.sqrt(network.ell_limit)*network.base/(np.sqrt(3)*network.voltage_kv)).tolist(),
            imported_sources=[], errors=[]))


def _reuse_ac_points(source, target):
    """复制同配置、同坐标的两套标签；两种扫描分别复用已完成点。"""
    if source['metadata']['identity'] != target['metadata']['identity']:
        return 0
    indices = np.argwhere((source['states'] != 0) | (source['socp_states'] != 0))
    power = source['axis_lower']+(indices+.5)*source['step']
    mapped = (power-target['axis_lower'])/target['step']-.5
    rounded = np.rint(mapped)
    good = (np.all(np.abs(mapped-rounded) <= 1e-8, axis=1)
            & np.all(rounded >= 0, axis=1) & np.all(rounded < target['states'].shape, axis=1))
    old, new = tuple(indices[good].T), tuple(rounded[good].astype(int).T)
    count = 0
    for prefix in ('', 'socp_'):
        values, previous = source[prefix+'states'][old], target[prefix+'states'][new]
        if np.any((values != 0) & (previous != 0) & (previous != values)):
            raise ValueError('Conflicting scan certificates at identical power coordinates')
        transfer = (values != 0) & ((previous == 0)
            | (~np.isfinite(target[prefix+'residual'][new]) & np.isfinite(source[prefix+'residual'][old])))
        count += int(np.count_nonzero((values != 0) & (previous == 0)))
        selected_old, selected_new = tuple(a[transfer] for a in old), tuple(a[transfer] for a in new)
        for key in ('states', 'witness_x', 'residual'):
            target[prefix+key][selected_new] = source[prefix+key][selected_old]
    target['metadata']['imported_sources'] = sorted(set(target['metadata']['imported_sources'])
        | set(source['metadata']['imported_sources']))
    return count


def import_ac_reference(network, budget, source, path=None, *, mode=1):
    """导入显式指定的同物理配置区域；AC 与 SOCP 数据始终成对保存。"""
    path = scan_path(network, budget, mode=mode) if path is None else Path(path)
    incoming = load_scan(source, ac_identity(network, budget, mode))
    incoming['metadata']['imported_sources'].append(str(Path(source).resolve()))
    destination = region_path(path, incoming)
    with _scan_lock(path/'scan'):
        count = sum(np.count_nonzero(incoming[key]) for key in ('states', 'socp_states'))
        if destination.exists():
            current = load_scan(destination)
            count = _reuse_ac_points(incoming, current)
            incoming = current
        save_scan(destination, incoming)
    return int(count)


def scan_problem(network, budget, sign, *, mode=1, ac=False, power=None, direction=None):
    """共用原始约束装配；AC 只改变电流关系，SOCP 保留锥不等式。"""
    equations = PortPhysics(deepcopy(network), sign) if mode else GridPhysics(deepcopy(network), 'socp')
    problem = MasterProblem(equations, power=power, budget=budget, direction=direction, threads=1)
    model = problem.model
    model.Params.NumericFocus = 3 if ac and network.name == 'case33bw' else 1
    model.Params.Presolve = (0 if ac else 2) if network.name == 'case33bw' else -1
    if power is not None:
        model.setObjective(0.)
        if not ac:
            eta = model.addVar(name='eta')
            for constraint in model.getConstrs():
                if constraint.ConstrName.startswith(('active_balance[', 'reactive_balance[', 'voltage_drop[')):
                    expression = model.getRow(constraint)-constraint.RHS
                    model.addConstr(expression <= eta, name=constraint.ConstrName+'_upper')
                    model.addConstr(expression >= -eta, name=constraint.ConstrName+'_lower')
                    model.remove(constraint)
            model.setObjective(eta)
    if ac:
        model.Params.NonConvex = 2
        model.update()
        for constraint in model.getQConstrs():
            if constraint.QCName.startswith('current_cone['):
                constraint.QCSense = '='
    model.update()
    return problem


def scan_line(args, *, network, budget, mode=1, ac=False):
    """独立逐点全拓扑求解；没有构域割、拓扑筛选或数值失败重试。求解失败、超时或可行解质量超过 SCAN_TOL 的点
    保持未决（0），不当作不可行。"""
    indices, coordinates = args
    states = np.zeros(len(indices), dtype=np.int8)
    witness_x = np.zeros((len(indices), network.n_types), dtype=np.int8)
    residual = np.full(len(indices), np.nan)
    signs = np.where(coordinates >= 0, 1, -1)
    with threadpool_limits(limits=1):
        for sign in np.unique(signs, axis=0):
            selected = np.flatnonzero(np.all(signs == sign, axis=1))
            problem = scan_problem(network, budget, sign, mode=mode, ac=ac, power=np.zeros(len(sign)))
            with problem.model as model:
                fixed = [c for c in model.getConstrs() if c.ConstrName.startswith('fixed_power[')]
                for j in selected:
                    model.ModelName = f'{"ac" if ac else "socp"}_point_{coordinates[j].tolist()}'
                    model.setAttr('RHS', fixed, abs(coordinates[j]) if mode else coordinates[j])
                    # 1. 求解；数值失败或超时的点保持未决
                    try:
                        answer = problem.solve(time_limit=60., tolerance=SCAN_TOL)
                    except (RuntimeError, TimeoutError):
                        continue
                    # 2. 不可行：已证不可行，或 SOCP 的 η 下界超过 PLANNING_TOL
                    if answer is None or (not ac and model.ObjBound > PLANNING_TOL):
                        states[j] = -1
                        continue
                    # 3. 可行见证的质量：SOCP 为 MaxVio+η，AC 再用潮流方程复核残差与运行界
                    residual[j], violation = model.MaxVio, 0.
                    if not ac:
                        residual[j] += model.ObjVal
                    else:
                        tree = problem.equations.network.tree(answer['x'])
                        ell = answer['state'][problem.equations.ell_slice][tree.type_indices]
                        P, Q, v, u = ACPowerFlow(tree).state(coordinates[j], ell)
                        residual[j] = float(np.max(np.abs(P*P+Q*Q-u*ell)))
                        violation = float(max(np.max(tree.vmin-v), np.max(v-tree.vmax), np.max(ell-tree.ell_limit),
                            np.max(P-tree.capacity), np.max(-P+tree.r*ell-tree.capacity),
                            P[:, tree.roots].sum()-tree.source_pmax, Q[:, tree.roots].sum()-tree.source_qmax,
                            np.hypot(P[:, tree.roots].sum(), Q[:, tree.roots].sum())-tree.source_smax))
                    if max(residual[j], violation) > SCAN_TOL:
                        continue
                    states[j], witness_x[j] = 1, answer['x']
    return dict(indices=indices, states=states, witness_x=witness_x, residual=residual,
                global_calls=len(indices), errors=[])


def ac_scan_line(args, *, network, budget, schemes, mode=1):
    """独立树递推提供 AC 见证；必要区间排除后未决点求完整 AC 等式。"""
    indices, coordinates = args
    states = np.zeros(len(indices), dtype=np.int8)
    witness_x = np.zeros((len(indices), network.n_types), dtype=np.int8)
    residual = np.full(len(indices), np.nan)
    possible = np.ones(len(indices), bool) if schemes is None else np.zeros(len(indices), bool)
    with threadpool_limits(limits=1):
        for x in (() if schemes is None else schemes):
            remaining = np.flatnonzero(states == 0)
            if not len(remaining):
                break
            answer = signed_ac_witness(network, x, coordinates[remaining], mode=mode)
            good = remaining[answer['feasible']]
            states[good], witness_x[good], residual[good] = 1, x, answer['residual'][answer['feasible']]
            remaining = remaining[~answer['feasible']]
            if len(remaining):
                possible[remaining] |= ac_interval_possible(network, x, coordinates[remaining], mode=mode)
    states[(states == 0) & ~possible] = -1
    remaining = np.flatnonzero(states == 0)
    if len(remaining):
        answer = scan_line((indices[remaining], coordinates[remaining]), network=network, budget=budget, mode=mode, ac=True)
        states[remaining], witness_x[remaining], residual[remaining] = answer['states'], answer['witness_x'], answer['residual']
    return dict(indices=indices, states=states, witness_x=witness_x, residual=residual,
                global_calls=len(remaining), errors=[])


def scan_ac_reference(network, budget, reference, path=None, *, workers=20, mode=1,
                      progress=lambda completed, total: None, force_rescan=False):
    """固定格架的区域缓存，同时补齐独立 AC 与 SOCP；每 SAVE_SECONDS 及结束或中断时保存。"""
    network = ac_network(network)
    path = scan_path(network, budget, mode=mode) if path is None else Path(path)
    lower, upper = np.asarray(reference['axis_lower'], float), np.asarray(reference['bounds'], float)
    shape, refinement = tuple(reference['shape']), reference.get('refinement', 1)
    started = perf_counter()
    identity = ac_identity(network, budget, mode)
    with _scan_lock(path/'scan'):
        sources = sorted(path.glob('region_*.npz'), key=lambda p: p.stat().st_mtime_ns)
        previous = load_scan(sources[-1], identity) if sources else None
        origin, step, start = lower.copy(), (upper-lower)/np.array(shape), np.zeros(len(shape), int)
        if previous is not None:
            lower, upper = np.minimum(lower, previous['axis_lower']), np.maximum(upper, previous['bounds'])
            origin, step = previous['origin'].copy(), previous['step'].copy()
            if refinement > 1:
                origin += .5*step*(1-1/refinement)
                step /= refinement
            start = np.floor((lower-origin)/step+1e-8).astype(int)
            end = np.ceil((upper-origin)/step-1e-8).astype(int)
            if not mode:
                start = np.maximum(start, np.ceil(-origin/step-.5-1e-8).astype(int))
            shape = tuple(end-start)
        answer = _empty_ac_grid(network, budget, mode, origin, step, start, shape)
        states = answer['states']
        destination = region_path(path, answer)
        if force_rescan:
            if destination.exists():
                shutil.copy2(destination, path/('before_'+uuid4().hex+'.npz'))
        else:
            for source in sources:
                _reuse_ac_points(load_scan(source, identity), answer)
        ac_pending, socp_pending = (np.flatnonzero(answer[key].ravel() == 0) for key in ('states', 'socp_states'))
        pending = np.union1d(ac_pending, socp_pending)
        reused = states.size-len(pending)
        progress(2*states.size-len(ac_pending)-len(socp_pending), 2*states.size)
        if len(pending):
            save_scan(destination, answer)
            schemes = budget_schemes(GridPhysics(deepcopy(network), 'socp'), budget, threads=1) if len(ac_pending) else None
            functions = {'': partial(ac_scan_line, network=network, budget=budget, schemes=schemes, mode=mode),
                         'socp_': partial(scan_line, network=network, budget=budget, mode=mode)}
            jobs = []
            for prefix, missing in (('', ac_pending), ('socp_', socp_pending)):
                for offset in range(0, len(missing), 128):
                    indices = missing[offset:offset+128]
                    cells = np.column_stack(np.unravel_index(indices, shape))
                    jobs.append((prefix, (indices, answer['axis_lower']+(cells+.5)*step)))
            completed = 2*answer['states'].size-len(ac_pending)-len(socp_pending)
            last_print = last_save = perf_counter()
            try:
                with ProcessPoolExecutor(max_workers=workers) if workers > 1 else nullcontext() as pool:
                    if pool:
                        futures = {pool.submit(functions[prefix], job): prefix for prefix, job in jobs}
                        lines = ((futures[future], future.result()) for future in as_completed(futures))
                    else:
                        lines = ((prefix, functions[prefix](job)) for prefix, job in jobs)
                    for prefix, line in lines:
                        indices = line['indices']
                        for key in ('states', 'residual'):
                            answer[prefix+key].ravel()[indices] = line[key]
                        answer[prefix+'witness_x'].reshape(-1, network.n_types)[indices] = line['witness_x']
                        completed += len(indices)
                        if perf_counter()-last_save > SAVE_SECONDS:
                            save_scan(destination, answer)
                            last_save = perf_counter()
                        progress(completed, 2*answer['states'].size)
                        if perf_counter()-last_print > 15:
                            print(f'AC/SOCP {completed}/{2*answer["states"].size} labels, reused paired points={reused}', flush=True)
                            last_print = perf_counter()
            finally:
                save_scan(destination, answer)
        else:
            save_scan(destination, answer)
        answer.update(scan_seconds=perf_counter()-started, reused_points=reused, computed_points=len(pending),
            ac_computed_points=len(ac_pending), socp_computed_points=len(socp_pending), cache_path=str(destination))
        print(f'AC/SOCP completed: {answer["states"].size} points, reused={reused}, '
              f'AC computed={len(ac_pending)}, SOCP computed={len(socp_pending)}', flush=True)
        return answer




def export_comparison(network, budget, result, ac_reference, output):
    """导出同坐标的实验、AC、SOCP 标签及三组遗漏率/多余率。"""
    from monitor import RunMonitor
    from region import contains, halfspaces
    states = np.asarray(ac_reference['states'])
    socp_states = np.asarray(ac_reference['socp_states'])
    if ac_reference.get('method') != AC_CACHE_METHOD:
        raise ValueError('Comparison requires an identified AC reference')
    mode = ac_reference['metadata'].get('mode', 1)
    if ac_reference['metadata']['identity'] != ac_identity(network, budget, mode):
        raise ACReferenceMismatch('Comparison AC physical identity differs')
    monitor = RunMonitor()
    monitor.validation(ac_reference, result)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    lower, upper = np.asarray(ac_reference['axis_lower']), np.asarray(ac_reference['bounds'])
    indices = np.indices(states.shape).reshape(states.ndim, -1).T
    power = lower+(indices+.5)*(upper-lower)/np.array(states.shape)
    labels = {}
    with threadpool_limits(limits=1):
        for key in ('inner', 'outer'):
            labels[key] = np.zeros(len(power), dtype=bool)
            for row in result[key]:
                labels[key] |= contains(power, halfspaces(row['vertices']))
    np.savez_compressed(output/'comparison.npz', comparison_method='paired_scan_v2', power=power,
        ac_states=states.ravel(), socp_states=socp_states.ravel(), **labels, load_nodes=network.load_nodes, budget=budget,
        current_limit_a=np.sqrt(network.ell_limit)*network.base/(np.sqrt(3)*network.voltage_kv),
        axis_lower=lower, bounds=upper, shape=states.shape, ac_identity=ac_reference['metadata']['identity'])
    with (output/'comparison.csv').open('w', encoding='utf-8', newline='') as stream:
        np.savetxt(stream, np.column_stack([power, states.ravel(), socp_states.ravel(), labels['inner'], labels['outer']]),
            delimiter=',', comments='', fmt=['%.9f']*states.ndim+['%d']*4,
            header=','.join([*(f'p_{node}_kw' for node in network.load_nodes), 'ac_state', 'socp_state', 'inner', 'outer']))
    summary = dict(comparison_method='paired_scan_v2', network=network.name,
        load_nodes=list(network.load_nodes), budget=budget, shape=list(states.shape), cells=int(states.size),
        axis_lower_kw=lower.tolist(), bounds_kw=upper.tolist(), power_unit='kW',
        ac_identity=ac_reference['metadata']['identity'], ac_cache=ac_reference.get('cache_path'),
        current_limit_a=ac_reference['metadata'].get('current_limit_a'),
        socp_feasible=int(np.count_nonzero(socp_states == 1)), ac_feasible=int(np.count_nonzero(states == 1)),
        socp_undecided=int(np.count_nonzero(socp_states == 0)), ac_undecided=int(np.count_nonzero(states == 0)),
        result_status=result['status'],
        result_certified=result['certified'], reused_points=ac_reference.get('reused_points', 0),
        computed_points=ac_reference.get('computed_points', 0),
        metrics=monitor.validation_state['validation']['metrics'],
        comparisons=monitor.validation_state['validation']['comparisons'])
    (output/'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
    return summary


def main():
    from Network.case33bw import Case33, CURRENT_LIMIT
    from Network.four_bus_five_corridor import FourBus
    from monitor import RunMonitor
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('recording', type=Path, help='为已有构域记录补算同限流 AC / SOCP 校验')
    parser.add_argument('--workers', type=int, default=20)
    parser.add_argument('--divisions', type=int)
    parser.add_argument('--refine', type=int, default=1, help='按整数倍加密已有格点，旧坐标仍复用')
    parser.add_argument('--force-rescan', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    monitor = RunMonitor(output=args.recording)
    monitor.load_recording(args.recording)
    state, result = monitor.state, monitor.state['result']
    network = (Case33(load_nodes=tuple(state['load_nodes']), current_limit=CURRENT_LIMIT)
               if state['network'] == 'case33bw' else FourBus(load_nodes=tuple(state['load_nodes'])))
    budget, mode = state['budget'], state['mode']
    path = scan_path(network, budget, mode=mode)
    lower, upper = reference_box(network, budget, mode=mode)
    lower, upper = np.minimum(lower, result['axis_lower']), np.maximum(upper, result['axis_bounds'])
    divisions = args.divisions or (160 if len(network.load_nodes) == 2 else 80)
    reference = scan_ac_reference(network, budget,
        dict(axis_lower=lower, bounds=upper, shape=(divisions,)*len(lower), refinement=args.refine),
        path, workers=args.workers, mode=mode, progress=monitor.scanning, force_rescan=args.force_rescan)
    monitor.validation_state.pop('error', None)
    monitor.validation(reference, result)
    monitor.save()
    output = args.output or args.recording.parent/(args.recording.name.removesuffix('.json.gz')+'_comparison')
    print(json.dumps(export_comparison(network, budget, result, reference, output), ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
