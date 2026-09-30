"""Adaptive support-face experiment; no enumeration of construction schemes.

Run from any directory::

    python experiments/test_support_face_certification_fourbus_2d.py --compare-global
    python experiments/test_support_face_certification_fourbus_2d.py --record
    python monitor.py results/support_face_test/fourbus_2d_physical/monitor.json.gz

Only FourBus(load_nodes=(1, 2)), nonnegative loads and the original complete SOCP
are used. All construction variables remain binary. A complete physical MISOCP
searches for an uncovered physical witness before a new scheme is constructed.
Its global upper bound is required for global certification; there is no candidate SP.
``--no-outer-cuts`` disables conditional support halfspaces for an ablation.
Scanning is independent and never used to select points or stop. The optional
physical global query is only a post-run comparison. Production code is unchanged.
"""
from __future__ import annotations

import argparse
import csv
import json
from itertools import product
from pathlib import Path
import sys
from time import perf_counter

import gurobipy as gp
from gurobipy import GRB
import numpy as np
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Network.four_bus_five_corridor import FourBus
from main import REGION_TAU
from model import GridPhysics, MasterProblem, PLANNING_TOL, RemainingRegionModel
from monitor import RunMonitor
from region import (GEOMETRY_TOL, clip_polytope, contains, halfspaces,
                    initial_polytope, polytope_vertices, polytope_volume, RegionState)
from vertify import validate_socp_region

BOUND_PAD = 1e-10
DEFAULT_OUTPUT = ROOT / 'results' / 'support_face_test' / 'fourbus_2d_physical'
CERTIFIED = {'GEOMETRY_CERTIFIED', 'SUPPORT_CERTIFIED'}


class SupportReplay:
    """Write actual execution events to the existing native monitor, without feedback."""

    def __init__(self, network, budget, tau, epsilon_geom, output, *, bounds=None, monitor=None):
        self.network, self.epsilon_geom = network, epsilon_geom
        bounds = np.full(len(network.load_nodes), network.power_limit) if bounds is None else np.asarray(bounds)
        self.region = RegionState(bounds, network.power_limit, tau)
        self.monitor = (RunMonitor(output=Path(output).resolve(), algorithm='支持面认证 · 按需发现')
                        if monitor is None else monitor)
        self.monitor.begin(network, 'socp', budget, self.region, np.inf)
        self.support_calls = 0
        self.monitor._emit('settings', phase='初始化', tau=tau,
            validation_note='构域后显示独立 SOCP 扫描；全局见证由完整 physical 模型求得并复核')

    def frame(self, x, inner, outer, event, phase, **values):
        self.region.add_scheme(x, self.network.decode_plan(x),
                               self.network.cost_offset+self.network.cost @ x)
        row = self.region.records[tuple(x)]
        row.update(inner=np.asarray(inner).copy(), outer=np.asarray(outer).copy(), inner_equations=None)
        label = self.monitor._scheme(x)
        self.monitor._emit(event, phase=phase, active_scheme=label,
            stage=len(self.region.records), **self.monitor._geometry(self.region), **values)

    def support_start(self, x, normal):
        label = self.monitor._scheme(x)
        self.monitor._emit('support_start', phase='网架支撑', support_scheme=label,
            support_direction=np.asarray(normal), active_scheme=label,
            sp_point=None, global_point=None, seed_point=None, eta=None, feasible=None)

    def support_end(self, x, inner, outer, normal, answer, status, hit=False):
        self.support_calls += int(not hit)
        label = self.monitor._scheme(x)
        point = None if answer['point'] is None else dict(
            scheme=label, p=answer['point']*self.region.bounds)
        title = {'INITIAL': '初始支持点', 'VIOLATED': '发现边界补点',
                 'SUPPORT_CERTIFIED': '支持面通过', 'UNRESOLVED': '支持面未决'}[status]
        phase = f'{title} · LB {answer["lb"]:.6g} / UB {answer["ub"]:.6g}'
        self.frame(x, inner, outer, 'support_end', phase, sp_point=point,
            support_calls=self.support_calls, support_normal=normal,
            support_lb=answer['lb'], support_ub=answer['ub'], support_status=status,
            support_cache_hit=hit)

    def support_cut(self, x, inner, outer, normal, ub):
        cut = conditional_support_cut(x, normal, ub, self.region.bounds, self.region.axis_bounds)
        self.region.cuts.append(cut)
        label = self.monitor._scheme(x)
        self.frame(x, inner, outer, 'cut', '支持上界裁剪',
            cut_history={str(len(self.region.cuts)): dict(cut=cut, scheme=label, kind='support')})

    def global_start(self):
        self.monitor._emit('residual_start', phase='完整物理查漏（MISOCP）', active_scheme=None,
            global_search=self.monitor.state['global_search']+1, global_point=None,
            sp_point=None, seed_point=None, eta=None, feasible=None)

    def global_end(self, answer):
        label = None if answer['x'] is None else self.monitor.scheme_ids.get(tuple(answer['x']), '候选')
        point = None if label is None else dict(scheme=label, x=answer['x'], p=answer['p'])
        phase = '完整物理域已覆盖' if answer['complete'] else '完整物理见证：原约束复核通过'
        if answer['bound'] is not None:
            phase += f' · 上界 {answer["bound"]:.6g}'
        self.monitor._emit('residual_end', phase=phase, global_point=point,
            coverage_bound=answer['bound'], coverage_complete=answer['complete'])

    def finish(self, summary):
        certified = summary['totals']['global_certified']
        inner = [dict(choice=row['choice'], cost=row['cost'], vertices=row['inner'])
                 for row in summary['schemes']]
        outer = ([dict(vertices=expanded_polygon(row['inner'], self.region.bounds,
                    self.region.tau, self.epsilon_geom)) for row in summary['schemes']]
                 if certified else [dict(vertices=initial_polytope(self.region.bounds,
                    self.region.total_bound, self.region.axis_bounds)*self.region.bounds)])
        result = dict(status=summary['totals']['status'], certified=certified,
            coverage_bound=summary['coverage']['bound'], inner=inner, outer=outer,
            axis_bounds=self.region.axis_bounds, counts=summary['totals'])
        self.monitor.finish(result, self.region)
        self.monitor.save()


def geometry_margins(inner_faces, outer_vertices, tau):
    """All faces at once, in xi=p/bounds; no solver or area-based decision."""
    A, b = inner_faces[:, :-1], inner_faces[:, -1]
    if not len(outer_vertices):
        # An unexpectedly empty floating-point polygon is not an infeasibility proof.
        return np.full(len(b), np.inf), np.full(len(b), np.inf)
    scores = outer_vertices @ A.T
    face_upper = scores.max(axis=0)
    return face_upper, (1-tau)*face_upper+b


def normal_key(normal):
    normal = np.asarray(normal, dtype=float)
    length = np.linalg.norm(normal)
    if not np.isfinite(length) or length == 0:
        raise ValueError('Support requires a nonzero finite normal')
    # Never flip a sign: h(a) and h(-a) are different support problems.
    return tuple(np.round(normal/length, 10))


def cached_support(entry, normal, box_upper):
    """Reuse nearby unit directions conservatively on 0<=xi<=box_upper."""
    length = np.linalg.norm(normal)
    direction = np.asarray(normal)/length
    difference = direction-np.asarray(entry['normal'])
    ub = length*(entry['ub']+np.maximum(difference, 0.) @ box_upper)
    point = entry['point']
    lb = -np.inf if point is None else float(np.asarray(normal) @ point)
    return dict(entry, lb=lb, ub=float(ub), normal=np.asarray(normal))


def classify_support(lb, ub, point, normal, offset, tau, epsilon_geom):
    """Only an upper bound certifies safety; only an audited point proves a miss."""
    rho_lb, rho_ub = (1-tau)*lb+offset, (1-tau)*ub+offset
    violated = point is not None and (1-tau)*float(normal @ point)+offset > epsilon_geom
    if violated and rho_ub <= epsilon_geom:
        raise RuntimeError('Inconsistent support upper bound and certified incumbent')
    if rho_ub <= epsilon_geom:
        return 'SUPPORT_CERTIFIED', rho_lb, rho_ub
    if violated:
        return 'VIOLATED', rho_lb, rho_ub
    return 'UNRESOLVED', rho_lb, rho_ub


def audit_incumbent(problem, x):
    """Substitute the actual solution in every original constraint; no extra solve.

    Includes topology/budget, bounds, integrality, quadratic rows, and the original
    Lorentz norms. The latter is an additional conservative physical check.
    """
    model = problem.model
    variables = model.getVars()
    values = np.asarray(model.getAttr('X', variables))
    if not np.all(np.isfinite(values)):
        return np.inf
    lower = np.asarray(model.getAttr('LB', variables))
    upper = np.asarray(model.getAttr('UB', variables))
    integer = np.asarray(model.getAttr('VType', variables)) != GRB.CONTINUOUS
    residual = [0., float(model.MaxVio), float(np.max(lower-values)),
                float(np.max(values-upper)), float(np.max(np.abs(problem.x.X-x)))]
    if np.any(integer):
        residual.append(float(np.max(np.abs(values[integer]-np.rint(values[integer])))))
    rows = model.getConstrs()
    errors = model.getA() @ values-np.asarray(model.getAttr('RHS', rows))
    senses = np.asarray(model.getAttr('Sense', rows))
    residual.append(float(np.max(np.where(senses == '=', np.abs(errors),
                                         np.where(senses == '<', errors, -errors)))))
    for row in model.getQConstrs():
        error = model.getQCRow(row).getValue()-row.QCRHS
        residual.append(abs(error) if row.QCSense == '=' else
                        error if row.QCSense == '<' else -error)
    if problem.operation is not None:
        for head, tail in problem.operation.cones:
            residual.append(float(np.linalg.norm([item.getValue() for item in tail])
                                  -head.getValue()))
    return float(max(residual))


class SupportOracle:
    """One complete, fixed-binary-scheme model with normalized support objectives."""

    def __init__(self, equations, x, bounds, budget, threads, time_limit):
        self.x, self.bounds = np.asarray(x, int), np.asarray(bounds, float)
        self.scheme_key = tuple(self.x)
        self.time_limit = time_limit
        self.problem = MasterProblem(equations, budget=budget,
            fixed_plan=equations.network.decode_plan(x), threads=threads)
        model = self.problem.model
        # Fixed settings from the first solve; never change precision after a failure.
        model.Params.NumericFocus, model.Params.BarHomogeneous = 3, 1
        model.Params.Aggregate, model.Params.ScaleFlag = 0, 0
        self.support_cache = {}
        self.calls = self.cache_hits = self.cache_misses = 0
        self.solve_seconds = 0.
        self.history = []

    def close(self):
        self.problem.model.dispose()

    def solve(self, normal):
        """Set max sum(a_i*p_i/bounds_i) directly: ObjVal/ObjBound need no /base."""
        problem, model = self.problem, self.problem.model
        normal = np.asarray(normal, float)/np.linalg.norm(normal)
        objective = gp.quicksum(float(a/b)*p for a, b, p in
                               zip(normal, self.bounds, problem.loads.values()))
        model.setObjective(objective, GRB.MAXIMIZE)
        model.Params.TimeLimit = self.time_limit
        self.calls += 1
        started = perf_counter()
        model.optimize()
        seconds = perf_counter()-started
        self.solve_seconds += seconds
        status = int(model.Status)
        reliable_status = status in (GRB.OPTIMAL, GRB.TIME_LIMIT, GRB.NODE_LIMIT,
                                    GRB.ITERATION_LIMIT, GRB.WORK_LIMIT, GRB.USER_OBJ_LIMIT)
        objective_value, bound = -np.inf, np.inf
        if model.SolCount:
            objective_value = float(model.ObjVal)
        if reliable_status:
            try:
                bound = float(model.ObjBound)
            except (gp.GurobiError, AttributeError):
                pass
        if not np.isfinite(bound) or abs(bound) >= GRB.INFINITY:
            bound = np.inf
        ub = float(np.nextafter(bound+BOUND_PAD*(1+abs(bound)), np.inf))
        point = state = None
        residual = np.inf
        if reliable_status and model.SolCount:
            residual = audit_incumbent(problem, self.x)
            if residual <= PLANNING_TOL:
                point, state = problem.power.X/self.bounds, problem.state.X.copy()
        lb = -np.inf if point is None else float(normal @ point)
        if lb > ub:
            raise RuntimeError('Audited support incumbent exceeds the solver upper bound')
        entry = dict(lb=lb, ub=ub, point=point, state=state, status=status,
                     solve_seconds=seconds, normal=normal, objective=objective_value,
                     bound=bound, max_violation=residual,
                     incumbent_certified=point is not None)
        self.history.append(entry)
        return entry

    def get(self, normal):
        key = (self.scheme_key, normal_key(normal))
        hit = key in self.support_cache
        if hit:
            self.cache_hits += 1
        else:
            self.cache_misses += 1
            self.support_cache[key] = self.solve(normal)
        # The common box is xi in [0,1]^d, independent of later outer clipping.
        return cached_support(self.support_cache[key], normal, np.ones(len(self.bounds))), hit


def zero_certificate(equations):
    """FourBus only: zero P/Q/ell/slacks, unit voltages, zero fixed loads."""
    net = equations.network
    if (type(net) is not FourBus or net.load_nodes not in ((1, 2), (1, 2, 3)) or
            np.any(net.fixed_p) or np.any(net.fixed_q) or
            np.any(net.vmin > 1.) or np.any(net.vmax < 1.) or
            min(net.source_pmax, net.source_qmax, net.source_smax) < 0.):
        raise ValueError('The analytic zero-load certificate applies only to this FourBus case')
    state = np.zeros(equations.slack_slice.stop)
    state[equations.v_slice] = 1.
    return dict(p=np.zeros(len(net.load_nodes)), state=state, max_violation=0.)


def initial_bounds(equations, budget, threads, time_limit, *, replay=None, bounds=None, deadline=np.inf):
    """d+1 global MP2 bounds in kW, plus the total-power solution as first seed."""
    d = len(equations.network.load_nodes)
    bounds = np.full(d, equations.network.power_limit) if bounds is None else np.asarray(bounds)
    answers = []
    for direction in np.vstack((np.eye(d), np.ones(d))):
        if replay is not None:
            replay.monitor.initializing(len(answers))
        problem = MasterProblem(equations, budget=budget, direction=direction, threads=threads)
        with problem.model:
            answer = problem.solve(time_limit=max(0., min(time_limit, deadline-perf_counter())))
        if answer is None:
            raise RuntimeError('The complete FourBus model unexpectedly has no feasible scheme')
        answers.append(answer)
        if replay is not None:
            axis = replay.region.axis_bounds.copy()
            total = replay.region.total_bound
            if len(answers) <= d:
                axis[len(answers)-1] = min(bounds[len(answers)-1], answer['bound']+BOUND_PAD*bounds.max())
            else:
                total = min(total, answer['bound']+BOUND_PAD*bounds[0])
            replay.region.tighten_bounds(axis, total)
            replay.monitor._emit('initial_bounds', phase='初始化',
                **replay.monitor._geometry(replay.region))
    axis_bounds = np.minimum(bounds, [item['bound']+BOUND_PAD*bounds.max() for item in answers[:d]])
    total_bound = min(equations.network.power_limit, answers[d]['bound']+BOUND_PAD*bounds.max())
    return bounds, axis_bounds, total_bound, answers[d]


def certify_scheme(equations, x, scheme, bounds, axis_bounds, total_bound, *, budget,
                   tau, epsilon_geom, threads, support_seconds, max_iterations,
                   outer_cuts=True, verbose=True, replay=None, deadline=np.inf):
    """Geometry first; every new feasible violation immediately invalidates old faces."""
    started = perf_counter()
    net = equations.network
    cost = float(net.cost_offset+net.cost @ x)
    origin = zero_certificate(equations)
    d = len(bounds)
    inner = np.zeros((1, d))
    outer = initial_polytope(bounds, total_bound, axis_bounds)
    oracle = SupportOracle(equations, x, bounds, budget, threads, support_seconds)
    face_history, iterations, violating_points = [], [], []
    certificates = [origin]
    unique_faces = set()
    certified, reason = False, 'ITERATION_LIMIT'
    initial_seconds = 0.
    if replay is not None:
        replay.frame(x, inner, outer, 'scheme_start', '开始当前方案支持面认证',
                     sp_point=None, global_point=None, seed_point=None)
    try:
        for normal in np.vstack((np.eye(d), np.ones(d))):
            oracle.time_limit = max(0., min(support_seconds, deadline-perf_counter()))
            if replay is not None:
                replay.support_start(x, normal)
            answer, _ = oracle.get(normal)
            if answer['point'] is not None:
                inner = polytope_vertices(np.vstack((inner, answer['point'])))
                certificates.append(dict(p=answer['point']*bounds, state=answer['state'],
                                         max_violation=answer['max_violation']))
            if replay is not None:
                replay.support_end(x, inner, outer, normal, answer, 'INITIAL')
            if outer_cuts and np.isfinite(answer['ub']):
                outer = clip_polytope(outer, answer['ub'], -normal)
                if replay is not None:
                    replay.support_cut(x, inner, outer, normal, answer['ub'])
        initial_inner, initial_outer = inner.copy(), outer.copy()
        initial_seconds = perf_counter()-started
        if verbose:
            print(f'\n{"="*56}\nSCHEME {scheme:02d} | cost {cost:.0f}\n{"="*56}', flush=True)
            print(f'initial inner vertices: {len(inner)} | initial outer vertices: {len(outer)}', flush=True)
        for iteration in range(1, max_iterations+1):
            if perf_counter() >= deadline:
                reason = 'TIME_LIMIT'
                break
            faces = halfspaces(inner)
            face_upper, margins = geometry_margins(faces, outer, tau)
            geometry_ok = margins <= epsilon_geom
            if replay is not None:
                replay.frame(x, inner, outer, 'geometry_check',
                    f'几何认证：{geometry_ok.sum()}/{len(faces)} 面通过',
                    geometry_certified=int(geometry_ok.sum()), faces=len(faces),
                    geometry_margins=margins, sp_point=None)
            rows = [dict(scheme=scheme, iteration=iteration, face=f,
                         **{f'normal_{j+1}': face[j] for j in range(d)}, offset=face[-1],
                         geometry_upper=face_upper[f], geometry_margin=margins[f],
                         remaining_margin=margins[f], lb=None, ub=None, rho_lb=None,
                         rho_ub=None, cache_hit=None, solve_seconds=0.,
                         status='GEOMETRY_CERTIFIED' if geometry_ok[f] else 'DEFERRED')
                    for f, face in enumerate(faces)]
            unique_faces.update(tuple(np.round(face, 10)) for face in faces)
            calls_before, hits_before = oracle.calls, oracle.cache_hits
            rebuilt, outer_updated = False, False
            pending = np.flatnonzero(~geometry_ok)
            pending = pending[np.argsort(-margins[pending], kind='stable')]
            if verbose:
                print(f'iteration {iteration:02d} | faces {len(faces)} | geometry '
                      f'{geometry_ok.sum()}/{len(faces)} ({geometry_ok.mean():.2%}) | '
                      f'inner/outer vertices {len(inner)}/{len(outer)}', flush=True)
            for f in pending:
                normal, offset = faces[f, :-1], faces[f, -1]
                oracle.time_limit = max(0., min(support_seconds, deadline-perf_counter()))
                if replay is not None:
                    replay.support_start(x, normal)
                answer, hit = oracle.get(normal)
                status, rho_lb, rho_ub = classify_support(
                    answer['lb'], answer['ub'], answer['point'], normal, offset, tau, epsilon_geom)
                if replay is not None:
                    replay.support_end(x, inner, outer, normal, answer, status, hit)
                rows[f].update(status=status, lb=answer['lb'], ub=answer['ub'],
                    rho_lb=rho_lb, rho_ub=rho_ub, cache_hit=hit,
                    solve_seconds=0. if hit else answer['solve_seconds'],
                    remaining_margin=min(margins[f], rho_ub), solver_status=answer['status'],
                    max_violation=answer['max_violation'])
                if verbose:
                    print(f'  normal={normal} b={offset:+.8g} geom={margins[f]:+.6g}\n'
                          f'  LB={answer["lb"]:.10g} UB={answer["ub"]:.10g} '
                          f'allowed={(epsilon_geom-offset)/(1-tau):.10g} '
                          f'{status}{" (cache)" if hit else ""}', flush=True)
                if outer_cuts and np.isfinite(answer['ub']):
                    updated = clip_polytope(outer, answer['ub'], -normal)
                    outer_updated = not np.array_equal(updated, outer)
                    outer = updated
                    if replay is not None and not hit:
                        replay.support_cut(x, inner, outer, normal, answer['ub'])
                if status == 'VIOLATED':
                    updated = polytope_vertices(np.vstack((inner, answer['point'])))
                    if np.array_equal(updated, inner):
                        rows[f]['status'] = 'UNRESOLVED'
                        rows[f]['note'] = 'Violated point lost at geometry resolution'
                    else:
                        inner, rebuilt = updated, True
                        violating_points.append(answer['point']*bounds)
                        certificates.append(dict(p=answer['point']*bounds, state=answer['state'],
                                                 max_violation=answer['max_violation']))
                        if replay is not None:
                            replay.frame(x, inner, outer, 'feasible', '加入认证边界点，重建内域')
                        if verbose:
                            print(f'  new point {answer["point"]*bounds} kW; rebuild inner hull', flush=True)
                if rebuilt or outer_updated:
                    # Recompute ALL faces/margins, including after a useful safe outer cut.
                    for row in rows:
                        if row['status'] == 'DEFERRED':
                            row['status'] = 'STALE_AFTER_REBUILD' if rebuilt else 'DEFERRED_OUTER_UPDATE'
                    break
            face_history.extend(rows)
            statuses = [row['status'] for row in rows]
            iterations.append(dict(iteration=iteration, number_of_inner_faces=len(faces),
                number_geometry_certified=int(geometry_ok.sum()),
                number_support_requested=oracle.calls-calls_before+oracle.cache_hits-hits_before,
                number_support_certified=statuses.count('SUPPORT_CERTIFIED'),
                number_violated=statuses.count('VIOLATED'), number_unresolved=statuses.count('UNRESOLVED'),
                geometry_certification_ratio=float(geometry_ok.mean()),
                support_calls=oracle.calls, support_cache_hits=oracle.cache_hits,
                support_cache_misses=oracle.cache_misses,
                inner_area=polytope_volume(inner)*np.prod(bounds),
                outer_area=polytope_volume(outer)*np.prod(bounds),
                max_geometry_margin=float(margins.max()),
                max_remaining_margin=max(row['remaining_margin'] for row in rows),
                support_solve_seconds=oracle.solve_seconds, total_seconds=perf_counter()-started))
            if rebuilt or outer_updated:
                continue
            certified = all(status in CERTIFIED for status in statuses)
            reason = 'CERTIFIED' if certified else 'UNRESOLVED'
            break
        # A read-only final assessment of the CURRENT hull; stale face lists never certify it.
        checked = perf_counter()
        faces = halfspaces(inner)
        upper, margins = geometry_margins(faces, outer, tau)
        remaining = margins.copy()
        final_statuses = []
        for f, face in enumerate(faces):
            status = 'GEOMETRY_CERTIFIED' if margins[f] <= epsilon_geom else 'UNRESOLVED'
            key = (tuple(x), normal_key(face[:-1]))
            if status == 'UNRESOLVED' and key in oracle.support_cache:
                entry = cached_support(oracle.support_cache[key], face[:-1], np.ones(d))
                status, _, rho_ub = classify_support(entry['lb'], entry['ub'], entry['point'],
                    face[:-1], face[-1], tau, epsilon_geom)
                remaining[f] = min(remaining[f], rho_ub)
            final_statuses.append(status)
        terminal_seconds = perf_counter()-checked
        certified = certified and all(status in CERTIFIED for status in final_statuses)
        if replay is not None:
            replay.frame(x, inner, outer, 'scheme_end',
                '当前方案全部面通过，继续全局查漏' if certified else '当前方案尚未完成认证',
                scheme_certified=certified, sp_point=None)
        runtime = perf_counter()-started
        return dict(scheme=scheme, x=np.asarray(x), choice=net.decode_plan(x), cost=cost,
            iterations=len(iterations), faces_seen=len(face_history), unique_faces=len(unique_faces),
            geometry_certified=sum(row['status'] == 'GEOMETRY_CERTIFIED' for row in face_history),
            geometry_certification_ratio=sum(row['status'] == 'GEOMETRY_CERTIFIED'
                                             for row in face_history)/max(1, len(face_history)),
            final_geometry_certification_ratio=float(np.mean(margins <= epsilon_geom)),
            support_calls=oracle.calls, initial_support_calls=d+1, refinement_support_calls=oracle.calls-d-1,
            support_cache_hits=oracle.cache_hits, support_cache_misses=oracle.cache_misses,
            support_certified=sum(row['status'] == 'SUPPORT_CERTIFIED' for row in face_history),
            violating_points_added=len(violating_points),
            unresolved_faces=final_statuses.count('UNRESOLVED'),
            final_inner_area=polytope_volume(inner)*np.prod(bounds),
            final_outer_area=polytope_volume(outer)*np.prod(bounds),
            max_geometry_margin=float(margins.max()), max_remaining_margin=float(remaining.max()),
            certified=certified, status=reason, runtime_seconds=runtime,
            initial_seconds=initial_seconds, refinement_seconds=runtime-initial_seconds,
            terminal_recheck_seconds=terminal_seconds, support_solve_seconds=oracle.solve_seconds,
            initial_inner_vertices=initial_inner*bounds, initial_outer_vertices=initial_outer*bounds,
            inner=inner*bounds, outer=outer*bounds, violating_points=np.asarray(violating_points).reshape(-1, d),
            face_history=face_history, iteration_history=iterations, certificates=certificates,
            support_history=oracle.history, final_faces=faces, final_face_statuses=final_statuses)
    finally:
        oracle.close()


def conditional_support_cut(x, normal, ub, bounds, axis_bounds):
    """Lift h_x(normal)<=ub to a valid joint cut without excluding any scheme.

    Hamming distance is zero only at the supported binary x. For every different
    binary scheme, M makes this inequality redundant on the common power box.
    The returned layout is the original alpha + beta@p + delta@x >= 0 (p in kW).
    """
    x, normal = np.asarray(x, int), np.asarray(normal, float)
    if not np.isfinite(ub):
        raise ValueError('A conditional support cut requires a finite upper bound')
    box_upper = np.maximum(normal, 0.) @ (axis_bounds/bounds)
    big_m = max(0., float(box_upper-ub))
    return np.r_[ub+big_m*x.sum(), -normal/bounds, big_m*(1-2*x)]


def discover_schemes(equations, seed, bounds, axis_bounds, total_bound, *, budget,
                     tau, epsilon_geom, threads, support_seconds, global_seconds,
                     max_iterations, max_discovery_iterations, outer_cuts, verbose, replay=None, deadline=np.inf):
    """Construct schemes only when an uncovered, physically feasible point needs one.

    All legal binary schemes and original physical constraints remain in every
    query. There are no candidate SPs, scheme enumeration or no-good cuts.
    Local certificates alone never terminate the global algorithm.
    """
    started = perf_counter()
    rows, cuts, history = [], [], []
    x = seed['x']
    complete, bound, error = False, None, None
    status = 'DISCOVERY_LIMIT'
    conditional_support_cuts = 0
    global_solve_seconds = 0.
    try:
        for iteration in range(max_discovery_iterations+1):
            if perf_counter() >= deadline:
                raise TimeoutError('Partition time limit')
            if x is not None:
                row = certify_scheme(equations, x, len(rows)+1, bounds, axis_bounds,
                    total_bound, budget=budget, tau=tau, epsilon_geom=epsilon_geom,
                    threads=threads, support_seconds=support_seconds,
                    max_iterations=max_iterations, outer_cuts=outer_cuts, verbose=verbose, replay=replay,
                    deadline=deadline)
                rows.append(row)
                if not row['certified']:
                    status = row['status']
                    break
                if outer_cuts:
                    for entry in row['support_history']:
                        if np.isfinite(entry['ub']):
                            cuts.append(conditional_support_cut(x, entry['normal'],
                                entry['ub'], bounds, axis_bounds))
                            conditional_support_cuts += 1
                x = None
            if iteration == max_discovery_iterations:
                break
            step = dict(iteration=iteration+1, schemes=len(rows), cuts=len(cuts),
                        complete=False, status='SEARCHING')
            history.append(step)
            query_started = perf_counter()
            if replay is not None:
                replay.global_start()
            problem = RemainingRegionModel(equations, budget, bounds, total_bound,
                cuts, [row['final_faces'] for row in rows], tau,
                axis_bounds=axis_bounds, threads=threads)
            with problem.model:
                try:
                    answer = problem.solve(epsilon_geom, time_limit=max(0., min(global_seconds, deadline-perf_counter())))
                    step.update(bound=answer['bound'], complete=answer['complete'])
                    if not answer['complete']:
                        if not answer['feasible'] or problem.problem.operation is None:
                            raise RuntimeError('Physical search returned no physical witness')
                        residual = audit_incumbent(problem.problem, answer['x'])
                        step.update(max_violation=residual, x=answer['x'], p=answer['p'])
                        if residual > PLANNING_TOL:
                            raise RuntimeError(f'Physical witness audit: {residual:g} > {PLANNING_TOL:g}')
                        step.update(candidate_feasible=True, state=problem.problem.state.X.copy())
                finally:
                    step['solver_status'] = int(problem.model.Status)
                    seconds = perf_counter()-query_started
                    global_solve_seconds += seconds
                    step['global_seconds'] = seconds
            if replay is not None:
                replay.global_end(answer)
            bound = answer['bound']
            if answer['complete']:
                complete, status = True, 'CERTIFIED'
                step['status'] = status
                break
            if any(np.array_equal(answer['x'], row['x']) for row in rows):
                # A real miss in a supposedly finished scheme must not be discarded.
                status = 'UNRESOLVED'
                error = 'Uncovered feasible witness conflicts with a local scheme certificate'
                step['status'] = status
                break
            else:
                x = answer['x']
                step['status'] = 'NEW_SCHEME'
            if verbose:
                print(f'discovery {iteration+1:03d} | {step["status"]} | '
                      f'constructed schemes {len(rows)} | bound {bound:+.6g}', flush=True)
    except (TimeoutError, RuntimeError) as exc:
        status = 'TIME_LIMIT' if isinstance(exc, TimeoutError) else 'UNRESOLVED'
        error = str(exc)
        if history and not history[-1]['complete']:
            history[-1].update(status=status, error=error)
    return dict(rows=rows, cuts=cuts, history=history, complete=complete,
                bound=bound, status=status, error=error,
                global_milp_calls=0, global_physical_calls=len(history), discovery_sp_calls=0,
                discovery_cut_lp_calls=0,
                conditional_support_cuts=conditional_support_cuts,
                global_solve_seconds=global_solve_seconds, discovery_sp_seconds=0.,
                total_seconds=perf_counter()-started)


def expanded_polygon(inner, bounds, tau, epsilon_geom):
    """All tolerance-expanded faces intersected with the nonnegative common box."""
    polygon = np.asarray(list(product((0., 1.), repeat=len(bounds))))
    for face in halfspaces(np.asarray(inner)/bounds):
        polygon = clip_polytope(polygon, epsilon_geom-face[-1], -(1-tau)*face[:-1])
    return polygon*bounds


def compare_global(equations, rows, budget, bounds, axis_bounds, total_bound,
                   tau, epsilon_geom, threads, time_limit):
    """A terminal query on exactly the same final union, not an end-to-end speedup."""
    all_faces = [halfspaces(row['inner']/bounds) for row in rows]
    # Keep a representative of each nested exclusion polygon, as in mainline.
    order = sorted(range(len(rows)), key=lambda i: -rows[i]['final_inner_area'])
    kept = []
    for i in order:
        if not any(contains(rows[i]['inner']/bounds, all_faces[j]).all() for j in kept):
            kept.append(i)
    started = perf_counter()
    problem = RemainingRegionModel(equations, budget, bounds, total_bound, [],
        [all_faces[i] for i in kept], tau, axis_bounds=axis_bounds, threads=threads)
    answer = dict(complete=False, bound=None)
    error = None
    with problem.model:
        try:
            answer = problem.solve(epsilon_geom, time_limit=time_limit)
        except (TimeoutError, RuntimeError) as exc:
            error = str(exc)
            try:
                answer['bound'] = float(problem.model.ObjBound)/problem.distance_scale
            except (gp.GurobiError, AttributeError):
                pass
        status = int(problem.model.Status)
    return dict(seconds=perf_counter()-started, global_mip_calls=1,
                status=status, error=error, input_schemes=len(rows), exclusion_schemes=len(kept),
                scope='terminal query on the same final inner union; not full region construction',
                **answer)


def scan_comparison(network, budget, rows, bounds, axis_bounds, divisions,
                    tau, epsilon_geom, threads, workers):
    last = [-1]
    def progress(completed, total):
        percent = int(100*completed/total)
        if percent//10 != last[0]:
            print(f'independent SOCP scan: {completed}/{total} ({percent}%)', flush=True)
            last[0] = percent//10
    reference = validate_socp_region(network, budget, divisions, axis_bounds,
                                    threads=threads, workers=workers, progress=progress)
    index = np.indices((divisions, divisions)).reshape(2, -1).T
    power = (index+.5)*axis_bounds/divisions
    inner, envelope = np.zeros(len(power), bool), np.zeros(len(power), bool)
    for row in rows:
        faces = halfspaces(row['inner']/bounds)
        inner |= contains(power/bounds, faces, tolerance=epsilon_geom)
        envelope |= contains((1-tau)*power/bounds, faces, tolerance=epsilon_geom)
    feasible = reference['states'].reshape(-1) == 1
    missing, false_inner = feasible & ~envelope, ~feasible & inner
    return dict(divisions=divisions, total_points=len(power), feasible_points=int(feasible.sum()),
                inner_points=int(inner.sum()), missing_points=int(missing.sum()),
                false_inner_points=int(false_inner.sum()),
                missing_rate=float(missing.sum()/feasible.sum()) if feasible.any() else None,
                false_inner_rate=float(false_inner.sum()/inner.sum()) if inner.any() else None,
                scan_seconds=reference['scan_seconds'], bounds=axis_bounds,
                states=reference['states'], missing_power=power[missing], false_inner_power=power[false_inner])


def save_plots(rows, bounds, axis_bounds, tau, epsilon_geom, output, reference):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon as Patch
    from matplotlib.lines import Line2D
    from scipy.spatial import ConvexHull
    from shapely.geometry import MultiPoint
    from shapely.ops import unary_union

    def polygon(ax, vertices, **style):
        if len(vertices) >= 3 and polytope_volume(vertices) > 0:
            ax.add_patch(Patch(vertices[ConvexHull(vertices).vertices], **style))
        elif len(vertices):
            ax.plot(vertices[:, 0], vertices[:, 1], color=style.get('edgecolor', '#0072B2'))

    def finish(fig, ax, name):
        ax.set(xlabel='Node 1 load (kW)', ylabel='Node 2 load (kW)',
               xlim=(-.8, axis_bounds[0]*1.05), ylim=(-.8, axis_bounds[1]*1.05))
        ax.set_aspect('equal', adjustable='box')
        ax.spines[['top', 'right']].set_visible(False)
        ax.legend(loc='upper right', frameon=False, fontsize=8)
        fig.tight_layout()
        for suffix in ('png', 'pdf'):
            fig.savefig(output/f'{name}.{suffix}', dpi=200, bbox_inches='tight')
        plt.close(fig)

    with plt.rc_context({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'pdf.fonttype': 42, 'axes.linewidth': .8}):
        for row in rows:
            fig, ax = plt.subplots(figsize=(6.4, 5.1))
            polygon(ax, row['outer'], facecolor='#eeeeee', edgecolor='#888888', linewidth=1.3, label='Safe outer')
            polygon(ax, expanded_polygon(row['inner'], bounds, tau, epsilon_geom),
                    facecolor='none', edgecolor='#D55E00', linewidth=1.4, linestyle='--', label='Allowed envelope')
            polygon(ax, row['inner'], facecolor='#56B4E966', edgecolor='#0072B2', linewidth=1.3, label='Certified inner')
            points = row['violating_points']
            if len(points):
                ax.scatter(*points.T, s=24, color='#CC79A7', marker='x', label='Added violations', zorder=5)
            checked = [f for f in row['face_history'] if f['ub'] is not None and np.isfinite(f['ub'])]
            if checked:
                face = checked[-1]
                normal = np.array([face['normal_1'], face['normal_2']])/bounds
                power = np.linspace(0., axis_bounds[0], 200)
                if abs(normal[1]) > 1e-12:
                    other = (face['ub']-normal[0]*power)/normal[1]
                    ax.plot(power, other, ':', color='#009E73', lw=1.2, label='Last checked support UB')
                else:
                    ax.axvline(face['ub']/normal[0], ls=':', color='#009E73', lw=1.2, label='Last checked support UB')
            ax.set_title(f'S{row["scheme"]:02d} | cost {row["cost"]:,.0f} | {row["status"]}\n'
                         f'Geometry visits {row["geometry_certification_ratio"]:.1%}; '
                         f'SOCP {row["support_calls"]} (initial 3)')
            finish(fig, ax, f'scheme_{row["scheme"]:02d}')
        fig, ax = plt.subplots(figsize=(7.0, 5.4))
        colors = plt.get_cmap('tab20')
        for i, row in enumerate(rows):
            polygon(ax, row['inner'], facecolor=(*colors(i % 20)[:3], .13),
                    edgecolor=colors(i % 20), linewidth=.8)
        envelope = unary_union([MultiPoint(expanded_polygon(row['inner'], bounds, tau, epsilon_geom)).convex_hull
                                for row in rows])
        for part in getattr(envelope, 'geoms', [envelope]):
            if hasattr(part, 'exterior'):
                ax.plot(*part.exterior.xy, color='#D55E00', lw=1.8, ls='--')
                for ring in part.interiors:
                    ax.plot(*ring.xy, color='#D55E00', lw=1.8, ls='--')
        legend = [Line2D([], [], color='#0072B2', label='Individual certified inner regions'),
                  Line2D([], [], color='#D55E00', ls='--', label='Union of allowed envelopes')]
        if reference is not None:
            n = reference['divisions']
            u, v = ((np.arange(n)+.5)*axis_bounds[j]/n for j in range(2))
            ax.contour(u, v, (reference['states'] == 1).T, levels=[.5], colors=['#222222'], linewidths=[1.0])
            legend.append(Line2D([], [], color='#222222', label=f'Independent SOCP grid ({n} x {n})'))
        ax.set_title(f'FourBus (1, 2) | {sum(r["certified"] for r in rows)}/{len(rows)} discovered schemes certified\n'
                     f'Radial tolerance {tau:g}; no convex hull across schemes')
        ax.legend(handles=legend, loc='upper right', frameon=False, fontsize=8)
        # Explicit handles survive finish's standard legend call.
        for handle in legend:
            ax.add_line(handle)
        finish(fig, ax, 'union')


def json_value(value):
    if isinstance(value, np.ndarray):
        return json_value(value.tolist())
    if isinstance(value, dict):
        return {str(k): json_value(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_value(v) for v in value]
    if isinstance(value, np.generic):
        return json_value(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def save_reports(summary, output):
    output.mkdir(parents=True, exist_ok=True)
    (output/'summary.json').write_text(json.dumps(json_value(summary), indent=2,
                                                ensure_ascii=False, allow_nan=False), encoding='utf-8')
    rows = summary['schemes']
    columns = ('scheme', 'cost', 'iterations', 'faces_seen', 'unique_faces', 'geometry_certified',
        'geometry_certification_ratio', 'final_geometry_certification_ratio', 'support_calls',
        'initial_support_calls', 'refinement_support_calls', 'support_cache_hits', 'support_cache_misses',
        'violating_points_added', 'unresolved_faces', 'final_inner_area', 'final_outer_area',
        'max_geometry_margin', 'max_remaining_margin', 'certified', 'status',
        'support_solve_seconds', 'runtime_seconds')
    with (output/'scheme_summary.csv').open('w', encoding='utf-8-sig', newline='') as stream:
        writer = csv.DictWriter(stream, columns, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(json_value(rows))
    history = [face for row in rows for face in row['face_history']]
    columns = list(dict.fromkeys(key for row in history for key in row))
    with (output/'face_history.csv').open('w', encoding='utf-8-sig', newline='') as stream:
        writer = csv.DictWriter(stream, columns)
        writer.writeheader()
        writer.writerows(json_value(history))


def run_experiment(*, budget=20000., tau=REGION_TAU, epsilon_geom=GEOMETRY_TOL,
                   threads=4, support_seconds=30., max_iterations=100, divisions=80,
                   workers=4, outer_cuts=True, comparison=False, global_seconds=30.,
                   scan=True, plots=True, output=DEFAULT_OUTPUT, verbose=True,
                   max_discovery_iterations=200, record=False):
    if (not np.isfinite(budget) or budget < 0 or not 0 <= tau < 1 or
            not np.isfinite(epsilon_geom) or epsilon_geom < 0 or
            min(threads, workers, max_iterations, max_discovery_iterations) < 1 or divisions < 2 or
            min(support_seconds, global_seconds) <= 0):
        raise ValueError('Invalid budget, tolerance, time limit, or positive count')
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    wall = perf_counter()
    print('Adaptive support certification: no scheme enumeration; '
          'global coverage requires the complete physical MISOCP certificate.', flush=True)
    network = FourBus(load_nodes=(1, 2))
    equations = GridPhysics(network, 'socp')
    replay = SupportReplay(network, budget, tau, epsilon_geom, output/'monitor.json.gz') if record else None
    with threadpool_limits(limits=1):
        started = perf_counter()
        bounds, axis_bounds, total_bound, seed = initial_bounds(equations, budget, threads, support_seconds,
                                                               replay=replay)
        initialization_seconds = perf_counter()-started
        discovery = discover_schemes(equations, seed, bounds, axis_bounds, total_bound,
            budget=budget, tau=tau, epsilon_geom=epsilon_geom, threads=threads,
            support_seconds=support_seconds, global_seconds=global_seconds,
            max_iterations=max_iterations, max_discovery_iterations=max_discovery_iterations,
            outer_cuts=outer_cuts, verbose=verbose, replay=replay)
        rows = discovery['rows']
        support_seconds_total = sum(row['runtime_seconds'] for row in rows)
        totals = dict(schemes=len(rows), certified_schemes=sum(row['certified'] for row in rows),
            faces_seen=sum(row['faces_seen'] for row in rows),
            unique_scheme_faces=sum(row['unique_faces'] for row in rows),
            geometry_certified=sum(row['geometry_certified'] for row in rows),
            support_calls=sum(row['support_calls'] for row in rows),
            initial_support_calls=sum(row['initial_support_calls'] for row in rows),
            refinement_support_calls=sum(row['refinement_support_calls'] for row in rows),
            support_cache_hits=sum(row['support_cache_hits'] for row in rows),
            violating_points_added=sum(row['violating_points_added'] for row in rows),
            global_certified=discovery['complete'] and all(row['certified'] for row in rows),
            status=discovery['status'], enumeration_mip_calls=0, initial_global_mp2_calls=3,
            enumeration_seconds=0., global_initialization_seconds=initialization_seconds,
            global_milp_calls=discovery['global_milp_calls'],
            global_physical_calls=discovery['global_physical_calls'],
            discovery_sp_calls=discovery['discovery_sp_calls'],
            discovery_cut_lp_calls=discovery['discovery_cut_lp_calls'],
            conditional_support_cuts=discovery['conditional_support_cuts'],
            global_solve_seconds=discovery['global_solve_seconds'],
            discovery_sp_seconds=discovery['discovery_sp_seconds'],
            support_certification_seconds=support_seconds_total,
            refinement_seconds=sum(row['refinement_seconds'] for row in rows),
            terminal_recheck_seconds=sum(row['terminal_recheck_seconds'] for row in rows),
            support_solve_seconds=sum(row['support_solve_seconds'] for row in rows),
            total_algorithm_seconds=initialization_seconds+discovery['total_seconds'])
        totals['geometry_certification_ratio'] = (totals['geometry_certified']/totals['faces_seen']
                                                  if totals['faces_seen'] else None)
        totals['support_calls_per_face_visit'] = (totals['support_calls']/totals['faces_seen']
                                                  if totals['faces_seen'] else None)
        summary = dict(schema='support-face-fourbus-physical-v3', settings=dict(
            strategy='adaptive_physical_discovery',
            budget=budget, load_nodes=(1, 2), tau=tau, epsilon_geom=epsilon_geom,
            planning_tol=PLANNING_TOL, bound_pad=BOUND_PAD, bounds=bounds,
            axis_bounds=axis_bounds, total_bound=total_bound, threads=threads,
            support_seconds=support_seconds, max_iterations=max_iterations, outer_cuts=outer_cuts,
            global_seconds=global_seconds, max_discovery_iterations=max_discovery_iterations,
            divisions=divisions, scan_workers=workers, gurobi_version=gp.gurobi.version()),
            totals=totals, schemes=rows, cuts=discovery['cuts'],
            discovery_history=discovery['history'], coverage=dict(
                complete=discovery['complete'], bound=discovery['bound'],
                status=discovery['status'], error=discovery['error'], mode='physical'),
            comparison=None, reference=None)
        if replay is not None:
            replay.finish(summary)
            summary['recording'] = str(replay.monitor.output)
        # Persist certificates before independent optional work, including failed scans.
        save_reports(summary, output)
        if comparison:
            print('Comparing existing physical global search on the final inner union...', flush=True)
            summary['comparison'] = compare_global(equations, rows, budget, bounds, axis_bounds,
                total_bound, tau, epsilon_geom, threads, global_seconds)
            save_reports(summary, output)
        if scan:
            summary['reference'] = scan_comparison(network, budget, rows, bounds, axis_bounds,
                divisions, tau, epsilon_geom, threads, workers)
            save_reports(summary, output)
            if replay is not None:
                replay.monitor.validation(summary['reference'], replay.monitor.state['result'])
                replay.monitor.save()
        if plots:
            save_plots(rows, bounds, axis_bounds, tau, epsilon_geom, output, summary['reference'])
        totals['total_wall_seconds'] = perf_counter()-wall
        save_reports(summary, output)
    print('\n================ FINAL SUMMARY ================')
    for key, value in totals.items():
        print(f'{key:32s}: {value}')
    if totals['global_certified']:
        print('The complete physical region over ALL legal x is covered by the discovered envelopes.')
        print('Therefore union_(all legal x) D_x is contained in union_(discovered x) E_tau(I_x).')
    elif summary['coverage']['error']:
        print('coverage unresolved             :', summary['coverage']['error'])
    if summary['reference'] is not None:
        print('scan missing rate               :', summary['reference']['missing_rate'])
        print('scan false-inner rate           :', summary['reference']['false_inner_rate'])
    if summary['comparison'] is not None:
        print('existing global-search seconds  :', summary['comparison']['seconds'])
        print('existing global complete        :', summary['comparison']['complete'])
    print('================================================', flush=True)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument('--budget', type=float, default=20000.)
    parser.add_argument('--tau', type=float, default=REGION_TAU)
    parser.add_argument('--epsilon-geom', type=float, default=GEOMETRY_TOL)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--support-seconds', type=float, default=30.)
    parser.add_argument('--max-iterations', type=int, default=100)
    parser.add_argument('--max-discovery-iterations', type=int, default=200)
    parser.add_argument('--divisions', type=int, default=80)
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--no-outer-cuts', action='store_true')
    parser.add_argument('--compare-global', action='store_true')
    parser.add_argument('--global-seconds', type=float, default=30.)
    parser.add_argument('--no-scan', action='store_true')
    parser.add_argument('--no-plots', action='store_true')
    parser.add_argument('--quiet', action='store_true')
    parser.add_argument('--record', action='store_true', help='Save native monitor.py replay frames')
    args = parser.parse_args()
    run_experiment(output=args.output, budget=args.budget, tau=args.tau, epsilon_geom=args.epsilon_geom,
        threads=args.threads, support_seconds=args.support_seconds, max_iterations=args.max_iterations,
        divisions=args.divisions, workers=args.workers, outer_cuts=not args.no_outer_cuts,
        comparison=args.compare_global, global_seconds=args.global_seconds,
        max_discovery_iterations=args.max_discovery_iterations,
        record=args.record,
        scan=not args.no_scan, plots=not args.no_plots, verbose=not args.quiet)


if __name__ == '__main__':
    main()
