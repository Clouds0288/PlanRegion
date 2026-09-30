"""Case33 带符号 AC 事后校验；直接读取主线记录，保留独立 SOCP 参考。"""
import argparse
import gzip
import json
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path
from time import perf_counter

import numpy as np
from gurobipy import GRB
from threadpoolctl import threadpool_limits

from Network.case33bw import Case33
from experiments.fourbus_outer_volume import budget_schemes
from model import GridPhysics, MasterProblem, PortPhysics, PLANNING_TOL, LOAD_PF, PV_PF, PV_Q_SIGN
from monitor import RunMonitor
from region import contains, halfspaces
from vertify import (ACPowerFlow, AC_TOL, AC_ITERATIONS, FIXED_POINT_TOL,
                     scan_path, scan_reference)


def signed_ac_witness(network, x, power):
    """独立树递推只提供可行证书；不使用正负荷单调性拒绝反送功率点。"""
    tree = network.tree(x)
    oracle = ACPowerFlow(tree, threads=1)
    power = np.asarray(power, dtype=float).reshape(-1, len(network.load_nodes))
    ratio = np.where(power >= 0., np.tan(np.arccos(LOAD_PF)),
                     PV_Q_SIGN*np.tan(np.arccos(PV_PF)))
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


def ac_interval_possible(network, x, power):
    """外扩的 AC 必要区间；False 为本树不可行，True 不构成可行证书。"""
    tree = network.tree(x)
    oracle = ACPowerFlow(tree, threads=1)
    power = np.asarray(power, dtype=float).reshape(-1, len(network.load_nodes))
    ratio = np.where(power >= 0., np.tan(np.arccos(LOAD_PF)),
                     PV_Q_SIGN*np.tan(np.arccos(PV_PF)))
    p = (tree.fixed_p+power@tree.E.T)/tree.base
    q = (tree.fixed_q+(power*ratio)@tree.E.T)/tree.base
    lower = np.zeros((len(power), tree.n))
    # Kirchhoff 电流定律和三角不等式，使用实际节点电压下限。
    upper = (np.hypot(p, q)/np.sqrt(tree.vmin-AC_TOL)@tree.D.T)**2+AC_TOL
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


def ac_scan_line(args, *, network, budget, schemes):
    """逐点 AC 证书；SOCP 不可行是有效排除，其余未决点解完整非凸等式。"""
    lower, upper, divisions, index, socp_states = args
    coordinates = lower+(np.arange(divisions)[:, None]*np.eye(len(lower))[-1]
                         +np.r_[index, 0.]+.5)*(upper-lower)/divisions
    states = np.where(socp_states == -1, -1, 0).astype(np.int8)
    with threadpool_limits(limits=1):
        possible = np.zeros(divisions, dtype=bool) if len(schemes) else np.ones(divisions, dtype=bool)
        for x in schemes:
            remaining = np.flatnonzero(states == 0)
            if not len(remaining):
                break
            answer = signed_ac_witness(network, x, coordinates[remaining])
            states[remaining[answer['feasible']]] = 1
            remaining = remaining[~answer['feasible']]
            if len(remaining):
                possible[remaining] |= ac_interval_possible(network, x, coordinates[remaining])
        states[(states == 0) & ~possible] = -1
        global_calls = 0
        for last_sign in (1, -1):
            selected = np.flatnonzero((states == 0) & ((coordinates[:, -1] >= 0) == (last_sign == 1)))
            if not len(selected):
                continue
            sign = np.where(coordinates[selected[0]] >= 0., 1, -1)
            equations = PortPhysics(Case33(load_nodes=network.load_nodes), sign)
            problem = MasterProblem(equations, power=np.zeros(len(lower)), budget=budget, threads=1)
            with problem.model as model:
                model.setObjective(0.)
                model.Params.NonConvex = 2
                model.Params.NumericFocus = 1
                model.Params.Aggregate = 1
                model.update()
                for constraint in model.getQConstrs():
                    if constraint.QCName.startswith('current_cone['):
                        constraint.QCSense = '='
                fixed = [row for row in model.getConstrs() if row.ConstrName.startswith('fixed_power[')]
                for j in selected:
                    model.ModelName = f'ac_scan_{index}_{j}'
                    model.setAttr('RHS', fixed, coordinates[j]*sign)
                    answer = problem.solve(time_limit=60.)
                    global_calls += 1
                    if answer is None:
                        states[j] = -1
                        continue
                    tree = equations.network.tree(answer['x'])
                    ell = answer['state'][equations.ell_slice][tree.type_indices]
                    P, Q, v, u = ACPowerFlow(tree).state(coordinates[j], ell)
                    residual = float(np.max(np.abs(P*P+Q*Q-u*ell)))
                    violation = float(max(np.max(tree.vmin-v), np.max(v-tree.vmax),
                        np.max(P-tree.capacity), np.max(-P+tree.r*ell-tree.capacity),
                        P[:, tree.roots].sum()-tree.source_pmax,
                        Q[:, tree.roots].sum()-tree.source_qmax,
                        np.hypot(P[:, tree.roots].sum(), Q[:, tree.roots].sum())-tree.source_smax))
                    if residual > PLANNING_TOL or violation > PLANNING_TOL:
                        raise RuntimeError(f'{model.ModelName}: independent AC residual={residual}, violation={violation}')
                    states[j] = 1
    return index, states, global_calls


def _save_ac_scan(path, network, budget, reference, states):
    """连同模型身份和对应 SOCP 网格保存，避免跨节点/预算或参考复用。"""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix('.tmp')
    with temporary.open('wb') as stream:
        np.savez_compressed(stream, axis_lower=reference['axis_lower'], bounds=reference['bounds'],
            states=states, socp_states=reference['states'], method='ac_equality_v1',
            network=network.name, load_nodes=network.load_nodes, budget=budget,
            power_factors=[LOAD_PF, PV_PF, PV_Q_SIGN])
    temporary.replace(path)


def scan_ac_reference(network, budget, reference, path, *, workers=20):
    """在完整 SOCP 同一坐标网格上生成 AC 电流等式参考，可从已完成行续算。"""
    lower, upper = np.asarray(reference['axis_lower']), np.asarray(reference['bounds'])
    socp_states = np.asarray(reference['states'], dtype=np.int8)
    if not np.isin(socp_states, [-1, 1]).all():
        raise ValueError('AC comparison requires a complete SOCP reference')
    shape = socp_states.shape
    states = np.zeros(shape, dtype=np.int8)
    checkpoint = path.with_suffix('.partial.npz')
    for saved_path in (path, checkpoint):
        if saved_path.exists():
            with np.load(saved_path) as saved:
                if (set(('method', 'network', 'load_nodes', 'budget', 'power_factors', 'socp_states')) <= set(saved.files)
                        and str(saved['method']) == 'ac_equality_v1' and str(saved['network']) == network.name
                        and np.array_equal(saved['load_nodes'], network.load_nodes) and float(saved['budget']) == budget
                        and np.array_equal(saved['power_factors'], [LOAD_PF, PV_PF, PV_Q_SIGN])
                        and np.array_equal(saved['socp_states'], socp_states)
                        and saved['states'].shape == shape and np.isin(saved['states'], [-1, 0, 1]).all()
                        and np.array_equal(saved['axis_lower'], lower) and np.array_equal(saved['bounds'], upper)):
                    states[:] = saved['states']
                    break
    if np.all(states != 0):
        print(f'Reuse AC scan: {path}', flush=True)
        return dict(axis_lower=lower, bounds=upper, states=states, method='ac_equality')
    schemes = budget_schemes(GridPhysics(network, 'socp'), budget, threads=1)
    jobs = ((lower, upper, shape[0], index, socp_states[index])
            for index in np.ndindex(shape[:-1]) if np.any(states[index] == 0))
    started = last_print = perf_counter()
    completed, global_calls = int(np.count_nonzero(states)), 0
    try:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            function = partial(ac_scan_line, network=network, budget=budget, schemes=schemes)
            for index, line, calls in pool.map(function, jobs):
                states[index] = line
                completed += shape[-1]
                global_calls += calls
                if perf_counter()-last_print > 15.:
                    _save_ac_scan(checkpoint, network, budget, reference, states)
                    print(f'AC {completed}/{states.size}, global calls {global_calls}, {perf_counter()-started:.1f}s', flush=True)
                    last_print = perf_counter()
    except BaseException:
        _save_ac_scan(checkpoint, network, budget, reference, states)
        raise
    answer = dict(axis_lower=lower, bounds=upper, states=states, method='ac_equality')
    _save_ac_scan(path, network, budget, reference, states)
    checkpoint.unlink(missing_ok=True)
    print(f'AC completed: {states.size} cells, feasible={np.count_nonzero(states == 1)}, {perf_counter()-started:.1f}s', flush=True)
    return answer


def export_comparison(network, budget, result, reference, ac_reference, output):
    """统一真实 kW 网格：AC、SOCP、内域、外域逐行对比。"""
    output.mkdir(parents=True, exist_ok=True)
    states = np.asarray(reference['states'])
    if (np.shape(ac_reference['states']) != states.shape
            or not np.array_equal(ac_reference['axis_lower'], reference['axis_lower'])
            or not np.array_equal(ac_reference['bounds'], reference['bounds'])
            or not np.isin(states, [-1, 1]).all()
            or not np.isin(ac_reference['states'], [-1, 1]).all()):
        raise ValueError('Comparison requires complete references on exactly the same grid')
    indices = np.indices(states.shape).reshape(states.ndim, -1).T
    lower, upper = np.asarray(reference['axis_lower']), np.asarray(reference['bounds'])
    power = lower+(indices+.5)*(upper-lower)/np.asarray(states.shape)
    labels = {}
    with threadpool_limits(limits=1):
        for key in ('inner', 'outer'):
            labels[key] = np.zeros(len(power), dtype=bool)
            for row in result[key]:
                labels[key] |= contains(power, halfspaces(row['vertices']))
    ac_states = np.asarray(ac_reference['states']).ravel()
    np.savez_compressed(output/'comparison.npz', power=power, socp_states=states.ravel(),
        ac_states=ac_states, **labels, load_nodes=network.load_nodes, budget=budget,
        axis_lower=lower, bounds=upper, shape=states.shape)
    with gzip.open(output/'comparison.csv.gz', 'wt', encoding='utf-8', newline='') as stream:
        header = ','.join([*(f'p_{node}_kw' for node in network.load_nodes), 'socp_state', 'ac_state', 'inner', 'outer'])
        np.savetxt(stream, np.column_stack([power, states.ravel(), ac_states, labels['inner'], labels['outer']]),
                   delimiter=',', header=header, comments='', fmt=['%.9f']*states.ndim+['%d']*4)
    metrics = {}
    for name, grid in (('socp', reference), ('ac', ac_reference)):
        monitor = RunMonitor()
        monitor.validation(grid, result)
        metrics[name] = monitor.validation_state['validation']['metrics']
    summary = dict(network=network.name, mode=1, load_nodes=list(network.load_nodes), budget=budget,
        shape=list(states.shape), cells=int(states.size), axis_lower=lower.tolist(), bounds=upper.tolist(),
        power_unit='kW', result_status=result['status'], result_certified=result['certified'],
        socp_feasible=int(np.count_nonzero(states == 1)), ac_feasible=int(np.count_nonzero(ac_states == 1)),
        socp_only=int(np.count_nonzero((states.ravel() == 1) & (ac_states == -1))), metrics=metrics)
    (output/'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('recording', type=Path)
    parser.add_argument('--workers', type=int, default=20)
    parser.add_argument('--divisions', type=int, default=80)
    parser.add_argument('--output', type=Path, default=Path('results/validation_case33_3d'))
    args = parser.parse_args()
    monitor = RunMonitor(output=args.recording)
    monitor.load_recording(args.recording)
    state, result = monitor.state, monitor.state['result']
    if state['network'] != 'case33bw' or state['mode'] != 1:
        raise ValueError('This validator currently supports signed Case33 recordings')
    network, budget = Case33(load_nodes=tuple(state['load_nodes'])), state['budget']
    reference = scan_reference(network, budget, args.divisions, result['axis_bounds'], scan_path(network, budget),
        axis_lower=result['axis_lower'], mode=1, threads=1, workers=args.workers, progress=monitor.scanning)
    monitor.validation_state.pop('error', None)
    monitor.validation(reference, result)
    monitor.save()
    ac_reference = scan_ac_reference(network, budget, reference, args.output/'ac_reference.npz', workers=args.workers)
    summary = export_comparison(network, budget, result, reference, ac_reference, args.output)
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)


if __name__ == '__main__':
    main()
