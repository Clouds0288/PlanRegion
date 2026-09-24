"""Case33 load-space boundary experiment. No production algorithms are edited.

python -m tests.benchmark_boundary_search --stage validate --output results/...
python -m tests.benchmark_boundary_search --stage run --variant C --budget 2 ...
All power interfaces are kW; geometry and theta are dimensionless.
"""
import argparse
import hashlib
import json
import platform
import sys
from itertools import product
from pathlib import Path
from time import perf_counter
from unittest.mock import patch

import gurobipy as gp
import numpy as np
from gurobipy import GRB
from threadpoolctl import threadpool_limits

import main as application
import model
from Network.case33bw import Case33
from region import GEOMETRY_TOL, contains, halfspaces
from tests.benchmark_math_acceleration import SOURCE_FILES, fingerprints, audit
from tests.planning_checks import margin
from vertify import ACPowerFlow

ROOT = Path(__file__).resolve().parents[1]
TAU = .002
BOUNDS = np.array([1520., 6130., 1490.])
DIRECTIONS = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1],
                       [1, 1, 1], [.5, 1, .5], [1, .5, .5]])


def json_value(value):
    if isinstance(value, np.ndarray):
        return json_value(value.tolist())
    if isinstance(value, np.generic):
        return json_value(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    if isinstance(value, dict):
        return {str(k): json_value(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_value(v) for v in value]
    return value


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(json_value(value), ensure_ascii=False, indent=2,
                               allow_nan=False), encoding='utf-8')


def maximal(points):
    """Exact coordinate dominance; never discard an outer box using tolerance."""
    kept = []
    for point in np.asarray(points).reshape(-1, 3):
        if any(np.all(point <= old) for old in kept):
            continue
        kept = [old for old in kept if not np.all(old <= point)]
        kept.append(point.copy())
    return np.asarray(kept).reshape(-1, 3)


def inner_scales(outer_vertices, inner_points):
    """Largest inner-box ray scale, capped at one, for every outer vertex."""
    inner_scale = np.zeros(len(outer_vertices))
    for j, point in enumerate(outer_vertices):
        active = point > 0.
        if not np.any(active):
            inner_scale[j] = float(len(inner_points) > 0)
        elif len(inner_points):
            inner_scale[j] = min(1., np.max(np.min(inner_points[:, active]/point[active], axis=1)))
    return inner_scale


def coverage(outer_vertices, inner_points):
    coverage_gap = 1.-inner_scales(outer_vertices, inner_points)
    residual = np.full(len(outer_vertices), np.inf)
    for j, point in enumerate(outer_vertices):
        if len(inner_points):
            residual[j] = np.min(np.max((1-TAU)*point-inner_points, axis=1))
    return coverage_gap, residual


def clip_boxes(outer_vertices, direction, theta_upper, bounds=BOUNDS):
    """Global disjunction OR_{h:direction_h>0} xi_h <= theta_upper*direction_h/b_h."""
    active = np.flatnonzero(direction > 0.)
    cutoff = theta_upper*direction/bounds+GEOMETRY_TOL
    updated, children = [], []
    for point in outer_vertices:
        if np.any(point[active] <= cutoff[active]):
            updated.append(point)
        else:
            for h in active:
                clipped = point.copy()
                clipped[h] = cutoff[h]
                children.append(clipped)
    # Existing vertices retain their age; all newly generated vertices go last.
    return maximal(updated+children)


class RayOracle:
    """Shared global cuts; every ray owns its own MP and original eta SP."""

    def __init__(self, network, budget, threads=20):
        self.network, self.budget, self.threads = network, budget, threads
        self.linear = model.PlanningEquations(network, 'linear')
        self.equations = model.PlanningEquations(network, 'socp')
        self.sp = model.PlanningSP(self.equations, threads=threads)
        self.cuts, self.certificates, self.calls = [], [], []
        self.timing = dict(build=0., mp=0., sp=0.)

    def query(self, direction, deadline, *, certificates=(), early=False):
        started = perf_counter()
        direction = np.asarray(direction)
        active = direction > 0.
        theta_lower = None
        for row in certificates:
            value = float(np.min(row['p'][active]/direction[active]))
            theta_lower = value if theta_lower is None else max(theta_lower, value)
        theta_upper = min(self.network.power_limit/direction.sum(), float(np.min(BOUNDS[active]/direction[active])))
        # In boundary search this ray starts at a maximal vertex of a certified
        # outer polyblock. No feasible point can strictly dominate that vertex.
        if early:
            theta_upper = min(theta_upper, 1.)
        tick = perf_counter()
        problem = model.PlanningModel(self.linear, budget=self.budget, cuts=self.cuts, threads=self.threads)
        m = problem.model
        theta = m.addVar(ub=theta_upper, name='theta')
        problem.power.UB = BOUNDS
        m.addConstrs((problem.loads[i]/self.network.base == theta*float(value/self.network.base)
                     for i, value in zip(self.network.load_nodes, direction)), name='ray')
        m.setObjective(theta, GRB.MAXIMIZE)
        self.timing['build'] += perf_counter()-tick
        count = len(self.calls)
        status = 'time_limit'
        with m:
            while perf_counter() < deadline:
                m.Params.TimeLimit = max(0., min(20., deadline-perf_counter()))
                tick = perf_counter()
                m.optimize()
                self.timing['mp'] += perf_counter()-tick
                entry = dict(direction=direction.copy(), status=m.Status,
                             seconds=perf_counter()-tick, bound=float(m.ObjBound),
                             candidate=None, feasible=False, cut=False)
                self.calls.append(entry)
                if m.Status == GRB.INFEASIBLE:
                    # Origin was established before normal boundary queries.
                    raise AssertionError('Ray MP infeasible despite a certified origin')
                if np.isfinite(m.ObjBound):
                    theta_upper = min(theta_upper, float(m.ObjBound)+GEOMETRY_TOL)
                if theta_lower is not None and theta_upper+GEOMETRY_TOL < theta_lower:
                    raise AssertionError('Global ray bound cuts a known feasible point')
                if not m.SolCount or perf_counter() >= deadline:
                    break
                entry['candidate'] = float(theta.X)
                x = np.rint(problem.x.X).astype(int)
                # Use the exact ray point, not the MP's rounded constraint residual.
                power = float(theta.X)*direction
                entry.update(x=x.copy(), p=power.copy())
                tick = perf_counter()
                checked = self.sp.solve(x, power, time_limit=max(0., min(10., deadline-perf_counter())))
                entry['sp_seconds'] = perf_counter()-tick
                self.timing['sp'] += entry['sp_seconds']
                entry.update(feasible=checked['feasible'], cut=checked['cut'] is not None)
                if checked['feasible']:
                    row = dict(x=x, p=power, state=checked['state'], theta=float(theta.X))
                    self.certificates.append(row)
                    theta_lower = row['theta'] if theta_lower is None else max(theta_lower, row['theta'])
                    if theta_upper-theta_lower <= TAU*max(theta_upper, 1e-12) or early and theta_lower >= 1-TAU:
                        status = 'bounded'
                        break
                    # Continue MIP bound improvement without imposing a cost cap.
                elif checked['cut'] is not None:
                    self.cuts.append(checked['cut'])
                    problem.add_cut(checked['cut'])
                else:
                    status = 'unknown'
                    break
                if early and theta_upper < 1-TAU:
                    status = 'clipped'
                    break
        return dict(direction=direction, theta_lower=theta_lower, theta_upper=theta_upper,
                    status=status, seconds=perf_counter()-started, mp=len(self.calls)-count)


def direct_ray(network, budget, direction, seconds, threads):
    started = perf_counter()
    equations = model.PlanningEquations(network, 'socp')
    problem = model.PlanningModel(equations, budget=budget, threads=threads)
    m = problem.model
    theta = m.addVar(name='theta')
    problem.power.UB = BOUNDS
    m.addConstrs((problem.loads[i]/network.base == theta*float(value/network.base)
                 for i, value in zip(network.load_nodes, direction)), name='ray')
    m.setObjective(theta, GRB.MAXIMIZE)
    with m:
        m.Params.TimeLimit = max(0., seconds-(perf_counter()-started))
        m.optimize()
        feasible = m.SolCount > 0 and m.MaxVio <= model.PLANNING_TOL
        row = None
        if feasible:
            row = dict(x=np.rint(problem.x.X).astype(int), p=problem.power.X,
                       state=problem.state.X, theta=float(theta.X))
        return dict(theta_lower=None if row is None else row['theta'],
                    theta_upper=float(m.ObjBound), solver_status=m.Status,
                    seconds=perf_counter()-started, certificate=row)


def origin(oracle, deadline):
    problem = model.PlanningModel(oracle.equations, power=np.zeros(3), budget=oracle.budget,
                                  threads=oracle.threads)
    with problem.model:
        answer = problem.solve(time_limit=max(0., min(20., deadline-perf_counter())))
    if answer is None or not answer['feasible']:
        return None
    return dict(x=answer['x'], p=answer['p'], state=answer['state'])


def boundary_search(network, budget, variant, seconds, threads, events):
    started = perf_counter()
    deadline = started+seconds
    oracle = RayOracle(network, budget, threads)
    tick = perf_counter()
    first = origin(oracle, deadline)
    origin_seconds = perf_counter()-tick
    certificates = [] if first is None else [first]
    outer_vertices = np.ones((1, 3))
    inner_points = np.array([r['p']/BOUNDS for r in certificates]).reshape(-1, 3)
    logical_cuts, queries, history = [], [], []
    geometry_seconds = selection_seconds = 0.
    peak = 1
    pending = set()
    status = 'unknown'

    def update(direction, early):
        nonlocal outer_vertices, inner_points, geometry_seconds, peak
        query = oracle.query(direction, min(deadline, perf_counter()+20.), certificates=certificates, early=early)
        queries.append(query)
        tick = perf_counter()
        certificates.extend(oracle.certificates)
        oracle.certificates.clear()
        inner_points = maximal([r['p']/BOUNDS for r in certificates])
        outer_vertices = clip_boxes(outer_vertices, direction, query['theta_upper'])
        logical_cuts.append(dict(direction=direction.copy(), theta_upper=query['theta_upper']))
        peak = max(peak, len(outer_vertices))
        geometry_seconds += perf_counter()-tick
        snapshot()

    def snapshot():
        tick = perf_counter()
        coverage_gap, residual = coverage(outer_vertices, inner_points)
        record = dict(elapsed=perf_counter()-started, max_coverage_gap=float(coverage_gap.max()),
                      max_covering_residual=float(residual.max()), outer_count=len(outer_vertices),
                      inner_count=len(inner_points), cuts=len(oracle.cuts), sp=oracle.sp.calls,
                      queries=len(queries), pending=len(pending), inner_points=inner_points.copy())
        history.append(record)
        events.write(json.dumps(json_value({k: v for k, v in record.items() if k != 'inner_points'}))+'\n')
        events.flush()
        return perf_counter()-tick

    snapshot()
    if first is not None:
        for direction in np.eye(3)*BOUNDS:
            if perf_counter() >= deadline:
                break
            update(direction, False)
        while perf_counter() < deadline:
            tick = perf_counter()
            coverage_gap, residual = coverage(outer_vertices, inner_points)
            available = [j for j in range(len(outer_vertices))
                         if residual[j] > GEOMETRY_TOL and tuple(outer_vertices[j]) not in pending]
            if np.all(residual <= GEOMETRY_TOL):
                status = 'certified'
                break
            if not available:
                pending.clear()
                available = np.flatnonzero(residual > GEOMETRY_TOL).tolist()
            index = available[0] if variant == 'B' else max(available, key=lambda j: coverage_gap[j])
            point = outer_vertices[index].copy()
            selection_seconds += perf_counter()-tick
            before = (len(oracle.cuts), len(certificates), outer_vertices.copy())
            update(point*BOUNDS, True)
            if np.array_equal(before[2], outer_vertices) and before[1] == len(certificates):
                pending.add(tuple(point))
            else:
                pending.intersection_update(map(tuple, outer_vertices))
    elapsed = perf_counter()-started
    coverage_gap, residual = coverage(outer_vertices, inner_points)
    if first is not None and np.all(residual <= GEOMETRY_TOL):
        status = 'certified'
    elif elapsed >= seconds:
        status = 'time_limit'
    snapshot()
    return dict(status=status, seconds=elapsed, inner_points=inner_points, outer_vertices=outer_vertices,
                certificates=certificates, joint_cuts=oracle.cuts, logical_cuts=logical_cuts,
                queries=queries, calls=oracle.calls, history=history,
                max_coverage_gap=float(coverage_gap.max()), max_covering_residual=float(residual.max()),
                counts=dict(sp=oracle.sp.calls, mp=1+len(oracle.calls), cuts=len(oracle.cuts),
                            queries=len(queries), schemes=len({tuple(r['x']) for r in certificates}),
                            peak_outer_vertices=peak, pending=len(pending)),
                timing=dict(**oracle.timing, origin=origin_seconds, geometry=geometry_seconds, selection=selection_seconds))


def check_certificates(network, budget, certificates, cuts=(), logical_cuts=()):
    tick = perf_counter()
    equations = model.PlanningEquations(network, 'socp')
    minimum_margin, minimum_cut = None, None
    failures, unknown, checked = [], 0, 0
    rng = np.random.default_rng(20260923)
    for j, row in enumerate(certificates):
        x, power = np.asarray(row['x']), np.asarray(row['p'])
        value = margin(equations, x, power, np.asarray(row['state']))
        minimum_margin = value if minimum_margin is None else min(value, minimum_margin)
        tree = network.tree(x)
        if value < -model.PLANNING_TOL or tree.cost > budget+1e-8:
            failures.append(dict(certificate=j, kind='raw_or_budget', margin=value, cost=tree.cost))
        points = np.vstack([power, power*rng.random((3, 3))])
        labels = ACPowerFlow(tree, threads=1).classify(points)
        checked += len(labels)
        unknown += int(np.sum(labels == 0))
        if np.any(labels == -1):
            failures.append(dict(certificate=j, kind='AC', labels=labels))
        if len(cuts):
            value = float(np.min(np.asarray(cuts)@np.r_[1., power, x]))
            minimum_cut = value if minimum_cut is None else min(minimum_cut, value)
            if value < -model.PLANNING_TOL:
                failures.append(dict(certificate=j, kind='joint_cut', value=value))
        for cut in logical_cuts:
            direction = np.asarray(cut['direction'])
            active = direction > 0.
            if np.all(power[active]/BOUNDS[active] >
                      cut['theta_upper']*direction[active]/BOUNDS[active]+GEOMETRY_TOL):
                failures.append(dict(certificate=j, kind='logical_cut'))
    return dict(failures=failures, ac_unknown=unknown, ac_points_checked=checked,
                certificates=len(certificates), minimum_original_margin=minimum_margin,
                minimum_joint_cut=minimum_cut, seconds=perf_counter()-tick)


def baseline(network, budget, variant, seconds, threads, events):
    certificates, cuts, history = [], [], []
    OriginalModel, OriginalSP, OriginalRemaining = model.PlanningModel, model.PlanningSP, model.RemainingRegionModel
    counts = dict(mp=0, sp=0, cuts=0)
    timing = dict(mp=0., sp=0., remaining_build=0., remaining=0.)

    class MeasuredModel(OriginalModel):
        def solve(self, *args, **kwargs):
            tick = perf_counter()
            result = super().solve(*args, **kwargs)
            timing['mp'] += perf_counter()-tick
            counts['mp'] += 1
            if result is not None and result['feasible']:
                certificates.append(dict(x=result['x'], p=result['p'], state=result['state']))
            return result

    class MeasuredSP(OriginalSP):
        def solve(self, x, power, *args, **kwargs):
            tick = perf_counter()
            result = super().solve(x, power, *args, **kwargs)
            timing['sp'] += perf_counter()-tick
            counts['sp'] += 1
            if result['feasible']:
                certificates.append(dict(x=x.copy(), p=power.copy(), state=result['state']))
            return result

    class MeasuredRemaining(OriginalRemaining):
        def __init__(self, equations, *args, **kwargs):
            tick = perf_counter()
            self.linear_outer = variant == 'A-L' and kwargs.get('mode') == 'light'
            if self.linear_outer:
                equations = model.PlanningEquations(network, 'linear')
                kwargs['mode'] = 'physical'
            super().__init__(equations, *args, **kwargs)
            timing['remaining_build'] += perf_counter()-tick

        def solve(self, *args, **kwargs):
            tick = perf_counter()
            result = super().solve(*args, **kwargs)
            timing['remaining'] += perf_counter()-tick
            if self.linear_outer:
                result['feasible'] = False
            elif result['feasible']:
                certificates.append(dict(x=result['x'], p=result['p'], state=self.problem.state.X))
            return result

    started = perf_counter()
    def progress(event, **data):
        if event == 'cut':
            cuts.append(data['cut'].copy())
        if event in ('feasible', 'region_end', 'phase_start'):
            history.append(dict(elapsed=perf_counter()-started,
                                inner=[r['inner'].copy() for r in data['records'] if len(r['inner'])]))
        if event in ('region_end', 'query_end', 'residual_start'):
            events.write(json.dumps(json_value(dict(event=event, elapsed=perf_counter()-started,
                counts=data['counts_algorithm'], schemes=len(data['records']))))+'\n')
            events.flush()

    with patch.object(application, 'PlanningModel', MeasuredModel), patch.object(application, 'PlanningSP', MeasuredSP), patch.object(application, 'RemainingRegionModel', MeasuredRemaining):
        result = application.build_continuous_region(network, 'socp', budget, BOUNDS, tau=TAU,
            time_limit=seconds, threads=threads, progress=progress)
    elapsed = perf_counter()-started
    counts.update(result['counts'], schemes=len(result['inner']))
    return dict(status=result['status'], seconds=elapsed, result=result, certificates=certificates,
                joint_cuts=cuts, counts=counts, timing=timing, history=history)


def grid_progress(data, variant, seconds):
    grid = (np.array(list(product(range(32), repeat=3)))+.5)/32
    records = []
    for at in sorted(set([v for v in (30., 60., 120., 300.) if v <= seconds]+[seconds])):
        snapshots = [h for h in data['history'] if h['elapsed'] <= at]
        mask = np.zeros(len(grid), dtype=bool)
        if snapshots:
            snapshot = snapshots[-1]
            if variant in ('B', 'C'):
                for point in np.asarray(snapshot['inner_points']).reshape(-1, 3):
                    mask |= np.all(grid <= point, axis=1)
            else:
                for poly in snapshot['inner']:
                    mask |= contains(grid, halfspaces(np.asarray(poly)))
        records.append(dict(seconds=at, inner_grid_points=int(mask.sum()),
                            grid_size=len(grid), sample_semantics='certified_inner_membership'))
    return records


def run(args):
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    network = Case33(upgrade_count=8)
    prepared = json.loads((ROOT/'results/case33_math_acceleration/20260923/bounds_n8.json').read_text(encoding='utf-8'))
    before = fingerprints()
    assert before == prepared['source_hashes'], 'Common bounds require identical source hashes'
    metadata = dict(schema_version='case33_boundary_v1', variant=args.variant, budget=args.budget,
                    seed=args.seed, time_limit=args.seconds, threads=args.threads, bounds=BOUNDS,
                    tau=TAU, source_hashes=before, network_fingerprint=network.fingerprint,
                    python=platform.python_version(), gurobi=gp.gurobi.version(), machine=platform.platform(),
                    experiment_hash=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    write_json(output/'manifest.json', metadata)
    (output/'TEST_PLAN.md').write_bytes((ROOT/'docs/case33_boundary_test_plan.md').read_bytes())
    original_new_model = model.new_model
    def seeded(name, threads):
        problem = original_new_model(name, threads)
        problem.Params.Seed = args.seed
        return problem
    with patch.object(model, 'new_model', seeded), threadpool_limits(limits=1), (output/'events.jsonl').open('w', encoding='utf-8') as events:
        if args.stage == 'validate':
            data = []
            for index, direction in enumerate(DIRECTIONS*BOUNDS):
                oracle = RayOracle(network, 2., args.threads)
                first = origin(oracle, perf_counter()+20.)
                assert first is not None, 'Background origin was not certified'
                result = oracle.query(direction, perf_counter()+args.seconds, certificates=[first])
                direct = direct_ray(network, 2., direction, args.seconds, args.threads)
                lower = [v for v in (result['theta_lower'], direct['theta_lower']) if v is not None]
                consistent = not lower or max(lower) <= min(result['theta_upper'], direct['theta_upper'])+GEOMETRY_TOL
                audit_result = check_certificates(network, 2., [first]+oracle.certificates+
                    ([] if direct['certificate'] is None else [direct['certificate']]), oracle.cuts)
                row = dict(index=index, direction=direction, ray=result, direct=direct, consistent=consistent,
                           audit=audit_result, calls=oracle.calls, joint_cuts=oracle.cuts,
                           certificates=[first]+oracle.certificates)
                data.append(row)
                write_json(output/'validation.json', data)
                print(json.dumps(json_value(dict(direction=index, ray=result, direct={k:v for k,v in direct.items() if k!='certificate'},
                                                  consistent=consistent, audit=audit_result))), flush=True)
                assert consistent and not audit_result['failures'], 'Validation failed; performance stages must stop'
            return
        function = baseline if args.variant in ('A', 'A-L') else boundary_search
        data = function(network, args.budget, args.variant, args.seconds, args.threads, events)
        data.update(metadata)
        write_json(output/'result.json', data)
        data['audit'] = check_certificates(network, args.budget, data['certificates'], data['joint_cuts'], data.get('logical_cuts', []))
        if args.variant in ('A', 'A-L'):
            data['inner_vertex_audit'] = audit(network, data['result'],
                [('socp', r['x'], r['p'], r['state']) for r in data['certificates']], data['joint_cuts'])
        data['grid_progress'] = grid_progress(data, args.variant, args.seconds)
        data['source_unchanged'] = before == fingerprints()
        write_json(output/'result.json', data)
        print(json.dumps(json_value({k:data[k] for k in ('variant', 'budget', 'seed', 'status', 'seconds', 'counts', 'grid_progress', 'audit')})), flush=True)
        assert not data['audit']['failures'], 'Audit failed; performance comparison must stop'
        if 'inner_vertex_audit' in data:
            assert not data['inner_vertex_audit']['ac_failures'], 'Inner vertex AC audit failed'


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--stage', choices=('run', 'validate'), default='run')
    parser.add_argument('--variant', choices=('A', 'A-L', 'B', 'C'), default='C')
    parser.add_argument('--budget', type=float, default=2.)
    parser.add_argument('--seconds', type=float, default=120.)
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--threads', type=int, default=20)
    parser.add_argument('--output', required=True)
    run(parser.parse_args())
