"""Independent 2-D planning experiment; all model powers remain 3-D and in kW.

python -m tests.case33_two_stage --mode probe --output results/case33_two_stage/probe
"""
import argparse
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from time import perf_counter

import gurobipy as gp
from gurobipy import GRB
import numpy as np
from shapely.geometry import Polygon, box, mapping
from shapely.ops import unary_union
from scipy.spatial import ConvexHull
from threadpoolctl import threadpool_limits

from Network.case33bw import Case33
from model import PlanningEquations, PlanningModel, PLANNING_TOL, evaluation_bounds
from tests.planning_checks import margin
from vertify import ACPowerFlow

ROOT = Path(__file__).resolve().parents[1]
SOURCE_FILES = ('main.py', 'model.py', 'region.py', 'vertify.py',
                'Network/__init__.py', 'Network/case33bw.py', 'Network/data/case33bw.m')
BOUND_PAD_KW = 1e-4


def serial(value):
    if isinstance(value, np.ndarray):
        return serial(value.tolist())
    if isinstance(value, np.generic):
        return serial(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    if isinstance(value, dict):
        return {str(k): serial(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [serial(v) for v in value]
    return value


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(serial(value), ensure_ascii=False, indent=2,
                               allow_nan=False), encoding='utf-8')


def fingerprints():
    return {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
            for name in SOURCE_FILES}


def add_tree_relaxation(problem):
    """Extended rooted-tree relaxation; preserves every feasible integer tree."""
    equations, model = problem.equations, problem.model
    net = equations.network
    assert np.all(net.required) and np.all(net.fixed_p >= 0.) and np.all(net.fixed_q >= 0.)
    arcs = [arc for ends in equations.ends.values() for arc in (ends, ends[::-1])]
    incoming = {i: [arc for arc in arcs if arc[1] == i] for i in (net.root, *net.nodes)}
    outgoing = {i: [arc for arc in arcs if arc[0] == i] for i in (net.root, *net.nodes)}
    parent_arc = model.addVars(arcs, ub=1., name='parent_arc')
    for arc in incoming[net.root]:
        parent_arc[arc].UB = 0.
    for e, (i, j) in equations.ends.items():
        model.addConstr(parent_arc[i, j]+parent_arc[j, i] ==
                        gp.quicksum(problem.choices[e, k] for k in equations.types[e]))
    for i in net.nodes:
        model.addConstr(gp.quicksum(parent_arc[arc] for arc in incoming[i]) == 1.)
        problem.operation.v[i].UB = min(1., equations.vmax[i])
    commodity_flow = model.addVars([(i, j, k) for i, j in arcs for k in net.nodes],
                                   ub=1., name='commodity_flow')
    for k in net.nodes:
        for i, j in arcs:
            model.addConstr(commodity_flow[i, j, k] <= parent_arc[i, j])
        for i in net.nodes:
            model.addConstr(gp.quicksum(commodity_flow[a, b, k] for a, b in incoming[i])-
                            gp.quicksum(commodity_flow[a, b, k] for a, b in outgoing[i]) == float(i == k))
    for e, (i, j) in equations.ends.items():
        for k in equations.types[e]:
            for flow, loss, upper in ((problem.operation.P[e, k], equations.r[e, k], equations.pmax[e, k]),
                     (problem.operation.Q[e, k], equations.reactance[e, k], equations.qmax[e, k])):
                model.addConstr(flow <= upper*parent_arc[i, j])
                model.addConstr(flow >= loss*problem.operation.ell[e, k]-upper*parent_arc[j, i])


class SliceOracle:
    """Full integer SOCP oracle, using original model constraints unchanged."""

    def __init__(self, budget=2., axes=(0, 1), threads=20, seed=0, strengthen=True):
        self.network = Case33(upgrade_count=8)
        self.equations = PlanningEquations(self.network, 'socp')
        self.problem = PlanningModel(self.equations, budget=budget, threads=threads)
        if strengthen:
            add_tree_relaxation(self.problem)
        self.axes = list(axes)
        self.fixed_axis = next(i for i in range(3) if i not in axes)
        self.budget = budget
        self.square_kw = self.network.power_limit
        self.bounds = evaluation_bounds(self.network, threads=threads)
        self.problem.model.Params.Seed = seed
        self.problem.model.Params.NumericFocus = 2
        self.certificates = []
        self.events = []
        self.trees = {}

    def solve(self, weights=(0., 1.), threshold=0., power=None, time_limit=120.,
              gap_kw=1., axis_only=False, label='query', upper=None, fixed_first=False):
        start = perf_counter()
        problem, net = self.problem, self.network
        model = problem.model
        weights = np.asarray(weights, dtype=float)
        problem.power.LB = np.zeros(3)
        problem.power.UB = self.bounds
        problem.power[self.fixed_axis].UB = 0.
        problem.power[self.axes[0]].LB = threshold
        if fixed_first:
            problem.power[self.axes[0]].UB = threshold
        if upper is not None:
            problem.power[self.axes[1]].UB = min(self.bounds[self.axes[1]], upper)
        if axis_only:
            for axis, weight in zip(self.axes, weights):
                if weight == 0.:
                    problem.power[axis].UB = 0.
        if power is not None:
            problem.power.LB = problem.power.UB = np.asarray(power)
        model.setObjective(gp.quicksum(float(w)*problem.power[a].item()/net.base
                                      for a, w in zip(self.axes, weights)), GRB.MAXIMIZE)
        model.Params.TimeLimit = time_limit
        model.Params.MIPGapAbs = gap_kw / net.base
        model.Params.SolutionLimit = 1 if power is not None else GRB.MAXINT
        warm_calls = 0
        if self.certificates:
            # Complete a few nearby, already known plans as continuous SOCPs.
            # Their solutions are starts only; no construction is fixed or
            # excluded in the global planning model.
            first_target = threshold if power is None else power[self.axes[0]]
            order = sorted(self.certificates, key=lambda c:
                abs(c['p'][self.axes[0]]-first_target) if threshold > 0. or power is not None
                else -float(weights @ c['p'][self.axes]))
            plans = []
            for candidate in order:
                key = tuple(candidate['x'])
                if key not in plans:
                    plans.append(key)
                if len(plans) == 4:
                    break
            model.update()
            x_indices = [v.index for v in problem.x.tolist()]
            p_indices = [v.index for v in problem.power.tolist()]
            state_indices = [v.index for v in problem.state.tolist()]
            for key in plans:
                fixed = model.copy()
                variables = fixed.getVars()
                for index, value in zip(x_indices, key):
                    variables[index].LB = variables[index].UB = value
                fixed.Params.Threads = 1
                fixed.Params.TimeLimit = 3.
                fixed.Params.SolutionLimit = GRB.MAXINT
                fixed.optimize()
                warm_calls += 1
                if fixed.SolCount and fixed.MaxVio <= PLANNING_TOL:
                    values = np.asarray(fixed.getAttr('X', variables))
                    x = np.asarray(key)
                    p, state = values[p_indices], values[state_indices]
                    raw_margin = margin(self.equations, x, p, state)
                    if key not in self.trees:
                        self.trees[key] = ACPowerFlow(net.tree(x))
                    if raw_margin >= -PLANNING_TOL and self.trees[key].classify(p)[0] == 1:
                        self.certificates.append(dict(x=x, p=p, state=state,
                            max_vio=float(fixed.MaxVio), raw_margin=raw_margin, origin='fixed_plan_start'))
                fixed.dispose()
        if self.certificates:
            eligible = [c for c in self.certificates
                        if np.all(c['p'] >= problem.power.LB-1e-7)
                        and np.all(c['p'] <= problem.power.UB+1e-7)]
            if axis_only:
                eligible = [c for c in eligible if all(c['p'][a] <= 1e-7
                    for a, w in zip(self.axes, weights) if w == 0.)]
            if power is not None:
                eligible = [c for c in eligible if np.max(np.abs(c['p']-power)) <= 1e-7]
            if eligible:
                chosen = max(eligible, key=lambda c: weights @ c['p'][self.axes])
                problem.x.Start = chosen['x']
                problem.power.Start = chosen['p']
                problem.state.Start = chosen['state']
        model.Params.TimeLimit = max(.01, time_limit-(perf_counter()-start))
        model.optimize()
        answer = dict(label=label, weights=weights, threshold=threshold,
                      supplied_upper=upper, fixed_first=fixed_first,
                      status=int(model.Status), bound=None, objective=None,
                      feasible=False, x=None, p=None, state=None, max_vio=None,
                      ac_status=None, raw_margin=None, nodes=float(model.NodeCount),
                      warm_calls=warm_calls,
                      seconds=perf_counter()-start)
        if model.Status == GRB.INFEASIBLE:
            answer['infeasible'] = True
        else:
            answer['infeasible'] = False
            if np.isfinite(model.ObjBound):
                answer['bound'] = float(model.ObjBound * net.base + BOUND_PAD_KW)
            if model.SolCount:
                x = np.rint(problem.x.X).astype(int)
                p, state = problem.power.X.copy(), problem.state.X.copy()
                answer.update(x=x, p=p, state=state, max_vio=float(model.MaxVio),
                              objective=float(weights @ p[self.axes]),
                              raw_margin=margin(self.equations, x, p, state))
                key = tuple(x)
                if key not in self.trees:
                    self.trees[key] = ACPowerFlow(net.tree(x))
                answer['ac_status'] = int(self.trees[key].classify(p)[0])
                answer['feasible'] = (model.MaxVio <= PLANNING_TOL
                                      and answer['raw_margin'] >= -PLANNING_TOL
                                      and answer['ac_status'] == 1
                                      and net.cost @ x <= self.budget+1e-8)
                if answer['feasible']:
                    self.certificates.append(dict(x=x, p=p, state=state,
                        max_vio=answer['max_vio'], raw_margin=answer['raw_margin']))
        self.events.append(answer)
        print(json.dumps(serial({k: answer[k] for k in ('label', 'status', 'bound',
              'objective', 'feasible', 'max_vio', 'ac_status', 'seconds')})), flush=True)
        return answer

    def close(self):
        self.problem.model.dispose()


@dataclass
class Interval:
    left: float
    right: float
    upper: float
    depth: int = 0


def clip_polygon(vertices, weights, bound):
    """Sutherland-Hodgman clipping; outer inequalities use their padded bound."""
    if len(vertices) == 0:
        return np.empty((0, 2))
    result = []
    for first, second in zip(vertices, np.roll(vertices, -1, axis=0)):
        one, two = np.dot(weights, first)-bound, np.dot(weights, second)-bound
        if one <= 0.:
            result.append(first)
        if (one <= 0.) != (two <= 0.):
            result.append(first+(second-first)*one/(one-two))
    return np.array(result).reshape(-1, 2)


def maximal(points):
    points = np.asarray(points).reshape(-1, 2)
    keep = []
    for point in points:
        if any(np.all(point <= old) for old in keep):
            continue
        keep = [old for old in keep if not np.all(old <= point)]
        keep.append(point)
    return np.array(keep).reshape(-1, 2)


def scheme_facets(certificates, axes=(0, 1)):
    """Positive-normal facets of downward convex inner regions, grouped by x."""
    groups = {}
    for certificate in certificates:
        groups.setdefault(tuple(certificate['x']), []).append(np.asarray(certificate['p'])[list(axes)])
    inner_facets = []
    for points in groups.values():
        points = np.asarray(points)
        cloud = np.unique(np.vstack((points, points*[1., 0.], points*[0., 1.], [0., 0.])), axis=0)
        if np.min(np.max(cloud, axis=0)) <= 0.:
            inner_facets.append(np.array([[1., 0., -cloud[:, 0].max()], [0., 1., -cloud[:, 1].max()]]))
        else:
            rows = ConvexHull(cloud).equations
            rows = rows[np.all(rows[:, :2] >= -1e-12, axis=1)]
            rows[:, :2] = np.maximum(rows[:, :2], 0.)
            rows[:, 2] = -np.max(cloud @ rows[:, :2].T, axis=0)
            inner_facets.append(rows)
    return inner_facets


def convex_distance(points, facets):
    """Infinity distance from points to one downward convex polygon."""
    distances = []
    for point in np.asarray(points):
        requirements = []
        smallest = int(np.argmin(point))
        largest = 1-smallest
        corner = np.maximum(point-point[smallest], 0.)
        for row in facets:
            weights, bound = row[:2], -row[2]
            if weights @ point <= bound:
                requirements.append(0.)
            elif weights @ corner <= bound:
                requirements.append((weights @ point-bound)/weights.sum())
            else:
                requirements.append(point[largest]-bound/weights[largest])
        distances.append(max(requirements, default=0.))
    return np.asarray(distances)


def interval_gap(node, supports, inner_points, inner_facets=None, outer_vertices=None):
    """Conservative directed L-infinity distance in kW for a complete strip."""
    upper = node.upper
    for cut in supports:
        wx, wy = cut['weights']
        if wy > 0.:
            upper = min(upper, (cut['bound']-wx*node.left)/wy)
    if upper < 0.:
        return 0.
    corner = np.array([node.right, upper])
    gap = float(np.min(np.max(np.maximum(corner-inner_points, 0.), axis=1)))
    if inner_facets is not None:
        vertices = clip_polygon(outer_vertices, [-1., 0.], -node.left)
        vertices = clip_polygon(vertices, [1., 0.], node.right)
        vertices = clip_polygon(vertices, [0., 1.], node.upper)
        if not len(vertices):
            return 0.
        # One complete node must fit in the expansion of ONE convex inner set.
        # Taking min over schemes separately for each vertex is not valid.
        convex_gap = min(float(convex_distance(vertices, facets).max()) for facets in inner_facets)
        gap = min(gap, convex_gap+1e-5)
    return gap


def region_from_intervals(stage1_outer, intervals):
    pieces = [stage1_outer.intersection(box(n.left, 0., n.right, max(0., n.upper)))
              for n in intervals if n.upper >= 0. and n.right > n.left]
    return unary_union(pieces)


def hull_direction(vertices, inner_points):
    """Most violated positive-normal facet of the convexified inner evidence."""
    cloud = np.vstack((inner_points, inner_points*[1., 0.], inner_points*[0., 1.], [0., 0.]))
    hull = ConvexHull(cloud)
    best = (-np.inf, None)
    for row in hull.equations:
        if np.min(row[:2]) < -1e-9 or row[:2].sum() <= 1e-12:
            continue
        weights = np.maximum(row[:2], 0.)
        weights /= weights.sum()
        gap = float(np.max(vertices @ weights)-np.max(inner_points @ weights))
        if gap > best[0]:
            best = gap, weights
    return best


def two_stage(output, epsilon_kw=20., time_limit=120., max_seconds=7200., resume=False):
    start = perf_counter()
    previous_result = json.loads((output/'result.json').read_text(encoding='utf-8')) if resume else None
    if not resume and (output/'result.json').exists():
        raise FileExistsError('Use a new output directory or --resume to preserve existing results.')
    oracle = SliceOracle()
    square_kw, axes = oracle.square_kw, oracle.axes
    vertices = np.array([[0., 0.], [square_kw, 0.], [square_kw, square_kw], [0., square_kw]])
    snapshots = [dict(stage=0, step=0, vertices=vertices.copy(), seconds=0.)]
    supports, logical_cuts, intervals, history = [], [], [], []
    metadata = dict(schema_version=1, network='case33bw', load_nodes=oracle.network.load_nodes,
        axes=axes, fixed_axis=oracle.fixed_axis, fixed_power_kw=0., budget=oracle.budget,
        upgrade_count=8, square_kw=square_kw, bounds=oracle.bounds,
        epsilon_kw=epsilon_kw, source_hashes=fingerprints(), algorithm='objective_space_branch_and_bound',
        threads=20, seed=0, physical_tolerance=PLANNING_TOL)
    for weights, bound, label in (([1., 1.], square_kw, 'source_total'),
             ([1., 0.], oracle.bounds[axes[0]], 'linear_axis_1'),
             ([0., 1.], oracle.bounds[axes[1]], 'linear_axis_2')):
        supports.append(dict(weights=weights, bound=float(bound), label=label, objective=None))
        vertices = clip_polygon(vertices, weights, bound)
        snapshots.append(dict(stage=1, step=len(supports), vertices=vertices.copy(),
                              cut=supports[-1], seconds=perf_counter()-start))

    def checkpoint(status, gap=None):
        inner_points = maximal([c['p'][axes] for c in oracle.certificates])
        result = dict(metadata=metadata, status=status, supports=supports,
            logical_cuts=logical_cuts, snapshots=snapshots, history=history,
            stage1_outer=vertices, inner_points=inner_points,
            intervals=[vars(n) for n in intervals], max_gap_kw=gap,
            certificates=oracle.certificates, events=oracle.events,
            seconds=perf_counter()-start)
        if intervals:
            result['stage2_outer'] = mapping(region_from_intervals(Polygon(vertices), intervals))
        save(output/'result.json', result)
        return result

    if previous_result is None:
        for weights in ((0., 1.), (1., 0.)):
            answer = oracle.solve(weights, time_limit=time_limit, gap_kw=1., axis_only=True,
                                  label='support_axis')
            if answer['bound'] is not None:
                supports.append({k: answer[k] for k in ('weights', 'bound', 'objective', 'label')})
                vertices = clip_polygon(vertices, weights, answer['bound'])
                snapshots.append(dict(stage=1, step=len(supports), vertices=vertices.copy(),
                                      cut=supports[-1], seconds=perf_counter()-start))
            checkpoint('stage1')
        if len(oracle.certificates) < 2:
            oracle.close()
            return checkpoint('insufficient_inner_certificates')
        previous = None
        for _ in range(12):
            inner_points = maximal([c['p'][axes] for c in oracle.certificates])
            support_gap, weights = hull_direction(vertices, inner_points)
            if support_gap <= epsilon_kw/2:
                break
            repeated = previous is not None and np.linalg.norm(weights-previous) < 1e-5
            answer = oracle.solve(weights, time_limit=time_limit*(2 if repeated else 1),
                                  gap_kw=1., label='support_adaptive')
            previous = weights.copy()
            if answer['bound'] is not None:
                supports.append({k: answer[k] for k in ('weights', 'bound', 'objective', 'label')})
                vertices = clip_polygon(vertices, weights, answer['bound'])
                snapshots.append(dict(stage=1, step=len(supports), vertices=vertices.copy(),
                                      cut=supports[-1], seconds=perf_counter()-start))
            checkpoint('stage1')
            if perf_counter()-start >= max_seconds/2:
                break
        stage1_seconds = perf_counter()-start
        metadata['stage1_seconds'] = stage1_seconds
        metadata['stage1_support_gap_kw'] = hull_direction(vertices,
            maximal([c['p'][axes] for c in oracle.certificates]))[0]
    else:
        assert previous_result['metadata']['source_hashes'] == fingerprints()
        metadata = previous_result['metadata']
        assert epsilon_kw == metadata['epsilon_kw']
        metadata['resumed_with_valid_interval_bounds'] = True
        interrupted = metadata.pop('interrupted_query_seconds', 0.)
        start -= previous_result['seconds']+interrupted
        metadata['interrupted_query_total_seconds'] = metadata.get('interrupted_query_total_seconds', 0.)+interrupted
        supports, logical_cuts = previous_result['supports'], previous_result['logical_cuts']
        snapshots, history = previous_result['snapshots'], previous_result['history']
        vertices = np.asarray(previous_result['stage1_outer'])
        intervals = [Interval(**n) for n in previous_result['intervals']]
        oracle.certificates = [{**c, **{k: np.asarray(c[k]) for k in ('x', 'p', 'state')}}
                               for c in previous_result['certificates']]
        oracle.events = previous_result['events']
        stage1_seconds = metadata['stage1_seconds']
    save(output/'metadata.json', metadata)
    metadata['tree_relaxation'] = True
    metadata['inner_certification'] = 'same_plan_convex_hulls'
    # Valid supports tighten all subsequent planning queries, without restricting x.
    for cut in supports:
        oracle.problem.model.addConstr(gp.quicksum(float(w)*oracle.problem.power[a].item()
            for a, w in zip(axes, cut['weights'])) <= float(cut['bound']))
    if not intervals:
        intervals = [Interval(0., float(vertices[:, 0].max()), float(vertices[:, 1].max()))]
    step = len(history)
    while perf_counter()-start < max_seconds:
        inner_points = maximal([c['p'][axes] for c in oracle.certificates])
        inner_facets = scheme_facets(oracle.certificates, axes)
        gaps = [interval_gap(n, supports, inner_points, inner_facets, vertices) for n in intervals]
        if not gaps or max(gaps) <= epsilon_kw:
            result = checkpoint('certified', max(gaps, default=0.))
            break
        selected = int(np.argmax(gaps))
        node = intervals[selected]
        threshold = (node.left+node.right)/2.
        # For a narrow unresolved strip improve its left-end bound instead of
        # infinitely bisecting a vertical solver gap.
        refine_bound = node.right-node.left <= epsilon_kw/2.
        if refine_bound:
            threshold = node.left
        answer = oracle.solve(threshold=threshold, time_limit=time_limit*(2 if refine_bound else 1),
                              gap_kw=epsilon_kw/8, label=f'branch_{step}',
                              upper=node.upper, fixed_first=True)
        upper = -1. if answer['infeasible'] else answer['bound']
        if upper is None:
            checkpoint('unresolved_bound', max(gaps))
            continue
        logical_cuts.append(dict(threshold=threshold, bound=upper,
                                 infeasible=answer['infeasible'], step=step))
        if refine_bound:
            intervals[selected].upper = min(node.upper, upper)
        else:
            intervals[selected:selected+1] = [Interval(node.left, threshold, node.upper, node.depth+1),
                Interval(threshold, node.right, min(node.upper, upper), node.depth+1)]
        for other in intervals:
            if other.left >= threshold:
                other.upper = min(other.upper, upper)
        intervals = [n for n in intervals if n.upper >= 0.]
        new_inner = maximal([c['p'][axes] for c in oracle.certificates])
        new_facets = scheme_facets(oracle.certificates, axes)
        new_gap = max((interval_gap(n, supports, new_inner, new_facets, vertices) for n in intervals), default=0.)
        history.append(dict(step=step, selected=vars(node), selected_gap_kw=gaps[selected],
            max_gap_kw=new_gap, threshold=threshold, bound=upper, query_seconds=answer['seconds'],
            active_intervals=len(intervals), seconds=perf_counter()-start))
        snapshots.append(dict(stage=2, step=step, intervals=[vars(n).copy() for n in intervals],
                              cut=logical_cuts[-1], max_gap_kw=new_gap, seconds=perf_counter()-start))
        print(f'BB {step}: max certified gap {new_gap:.3f} kW; {len(intervals)} strips', flush=True)
        checkpoint('stage2', new_gap)
        step += 1
    else:
        result = checkpoint('time_limit', max(gaps, default=None))
    metadata['stage2_seconds'] = perf_counter()-start-stage1_seconds
    result = checkpoint(result['status'], result['max_gap_kw'])
    oracle.close()
    return result


def reference_scan(output, spacing_kw=5., time_limit=120., threads=20, resume=False,
                   start_index=None, stop_index=None):
    """Independent regular-grid reference; no stage-1/2 cuts or plans imported.

    Column global bounds certify infeasible grid points. Every feasible grid
    point is checked by independent AC using a witness construction. This is
    a pointwise grid classification with monotonic acceleration, not N^2 MIPs.
    """
    start = perf_counter()
    previous = json.loads((output/'scan.json').read_text(encoding='utf-8')) if resume else None
    if not resume and (output/'scan.json').exists():
        raise FileExistsError('Use a new reference output or --resume.')
    oracle = SliceOracle(seed=17, threads=threads)
    net, axes = oracle.network, oracle.axes
    coordinates = np.arange(0., oracle.square_kw+1e-8, spacing_kw)
    if coordinates[-1] < oracle.square_kw:
        coordinates = np.r_[coordinates, oracle.square_kw]
    states = np.zeros((len(coordinates), len(coordinates)), dtype=np.int8)
    # These bounds are recomputed in this independent run.
    states[coordinates > oracle.bounds[axes[0]], :] = -1
    states[:, coordinates > oracle.bounds[axes[1]]] = -1
    columns = []
    previous_upper = float(oracle.bounds[axes[1]])
    first_index = 0
    if previous is not None:
        assert previous['source_hashes'] == fingerprints() and previous['spacing_kw'] == spacing_kw
        grid = np.load(output/'scan_grid.npz')
        coordinates, states = grid['coordinates'], grid['states']
        columns = previous['columns']
        first_index = columns[-1]['index']+1
        previous_upper = columns[-1]['upper']
        start -= previous['seconds']+previous.get('interrupted_query_seconds', 0.)
        oracle.events = previous['events']
        oracle.certificates = [{**c, **{k: np.asarray(c[k]) for k in ('x', 'p', 'state')}}
                               for c in previous['certificates']]
        for c in oracle.certificates:
            key = tuple(c['x'])
            if key not in oracle.trees:
                oracle.trees[key] = ACPowerFlow(net.tree(c['x']))
    if start_index is not None:
        first_index = max(first_index, start_index)
    for i, threshold in enumerate(coordinates):
        if i < first_index:
            continue
        if stop_index is not None and i >= stop_index:
            break
        if threshold > oracle.bounds[axes[0]]:
            break
        top = coordinates[(coordinates <= previous_upper) & (states[i] != -1)][-1]
        trial = np.zeros(3)
        trial[axes] = threshold, top
        top_witness = next((key for key, tree in oracle.trees.items()
                            if tree.classify(trial)[0] == 1), None)
        # Thresholds increase in this independent scan. A preceding global
        # vertical bound remains valid for every subsequent column.
        if top_witness is not None:
            answer = dict(infeasible=False, bound=previous_upper)
        else:
            oracle.problem.model.addConstr(oracle.problem.power[axes[1]] <= previous_upper)
            answer = oracle.solve(threshold=threshold, time_limit=time_limit, gap_kw=spacing_kw*.95,
                                  label=f'scan_column_{i}')
        if answer['infeasible']:
            states[i:, :] = -1
            columns.append(dict(index=i, threshold=threshold, lower=None, upper=-1.))
            break
        upper = min(previous_upper, answer['bound'] if answer['bound'] is not None else previous_upper)
        previous_upper = upper
        states[i:, coordinates > upper] = -1
        # Independent AC evaluates every point below a candidate boundary in
        # batches. Failed checks under one plan never certify global infeasibility.
        eligible = [c for c in oracle.certificates if c['p'][axes[0]] >= threshold-1e-8]
        if top_witness is not None:
            eligible.append(dict(x=np.asarray(top_witness), p=trial))
        lower = None
        witness_x = None
        for certificate in sorted(eligible, key=lambda c: c['p'][axes[1]], reverse=True):
            candidates = np.where((states[i] == 0) & (coordinates <= certificate['p'][axes[1]]+1e-8))[0]
            if not len(candidates):
                continue
            power = np.zeros((len(candidates), 3))
            power[:, axes[0]], power[:, axes[1]] = threshold, coordinates[candidates]
            ac = oracle.trees[tuple(certificate['x'])].classify(power)
            states[i, candidates[ac == 1]] = 1
            if np.any(states[i] == 1):
                candidate_lower = float(coordinates[np.where(states[i] == 1)[0][-1]])
                if lower is None or candidate_lower > lower:
                    lower, witness_x = candidate_lower, certificate['x']
        # Resolve the usually single grid point between the incumbent and bound
        # by a fixed-load full planning query, with no imported stage cuts.
        unresolved = np.where(states[i] == 0)[0]
        for j in unresolved[::-1]:
            power = np.zeros(3)
            power[axes] = threshold, coordinates[j]
            point = oracle.solve(weights=(0., 0.), power=power,
                                 time_limit=time_limit, label=f'scan_point_{i}_{j}')
            if point['infeasible']:
                states[i:, j:] = -1
            elif point['feasible']:
                candidates = np.where((states[i] == 0) & (coordinates <= coordinates[j]))[0]
                batch = np.zeros((len(candidates), 3))
                batch[:, axes[0]], batch[:, axes[1]] = threshold, coordinates[candidates]
                ac = oracle.trees[tuple(point['x'])].classify(batch)
                states[i, candidates[ac == 1]] = 1
                lower = float(coordinates[np.where(states[i] == 1)[0][-1]])
                witness_x = point['x']
                break
        columns.append(dict(index=i, threshold=threshold, lower=lower, upper=upper, witness_x=witness_x,
                            cached_AC_witness=top_witness is not None))
        np.savez_compressed(output/'scan_grid.npz', coordinates=coordinates, states=states)
        save(output/'scan.json', dict(schema_version=1, spacing_kw=spacing_kw, columns=columns,
            seconds=perf_counter()-start, events=oracle.events, certificates=oracle.certificates,
            counts={str(k): int(np.count_nonzero(states == k)) for k in (-1, 0, 1)},
            source_hashes=fingerprints(), independent=True, seed=17, threads=threads, tree_relaxation=True))
        print(f'SCAN {i}: p1={threshold:.1f}; feasible <= {lower}; global upper {upper:.3f}', flush=True)
    np.savez_compressed(output/'scan_grid.npz', coordinates=coordinates, states=states)
    result = dict(schema_version=1, spacing_kw=spacing_kw, columns=columns,
        seconds=perf_counter()-start, events=oracle.events, certificates=oracle.certificates,
        counts={str(k): int(np.count_nonzero(states == k)) for k in (-1, 0, 1)},
        source_hashes=fingerprints(), independent=True, seed=17, threads=threads, tree_relaxation=True,
        status='classified_grid' if np.all(states != 0) else 'unknown_grid_points',
        start_index=first_index, stop_index=stop_index)
    save(output/'scan.json', result)
    oracle.close()
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=('probe', 'tree_probe', 'run', 'scan'), default='probe')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--time-limit', type=float, default=60.)
    parser.add_argument('--epsilon-kw', type=float, default=20.)
    parser.add_argument('--spacing-kw', type=float, default=5.)
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--threads', type=int, default=20)
    parser.add_argument('--start-index', type=int)
    parser.add_argument('--stop-index', type=int)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    original = fingerprints()
    with threadpool_limits(limits=1):
        if args.mode == 'tree_probe':
            oracle = SliceOracle(strengthen=True)
            previous = json.loads((args.output.parent/'run/result.json').read_text(encoding='utf-8'))
            oracle.certificates = [{**c, **{k: np.asarray(c[k]) for k in ('x', 'p', 'state')}}
                                   for c in previous['certificates']]
            oracle.solve(threshold=848.4519119781603, upper=1376.7965564353835, fixed_first=True,
                         time_limit=args.time_limit, gap_kw=2.5, label='tree_relaxation_probe')
            save(args.output/'probe.json', oracle.events)
            oracle.close()
        elif args.mode == 'probe':
            oracle = SliceOracle()
            for weights in ((0., 1.), (1., 0.), (.75, .25), (.5, .5)):
                oracle.solve(weights, time_limit=args.time_limit, gap_kw=2.,
                             axis_only=sum(w > 0 for w in weights) == 1,
                             label='probe_'+str(weights))
                save(args.output/'probe.json', oracle.events)
            oracle.close()
        elif args.mode == 'run':
            two_stage(args.output, epsilon_kw=args.epsilon_kw, time_limit=args.time_limit, resume=args.resume)
        else:
            reference_scan(args.output, spacing_kw=args.spacing_kw, time_limit=args.time_limit,
                           threads=args.threads, resume=args.resume,
                           start_index=args.start_index, stop_index=args.stop_index)
    assert original == fingerprints(), 'Production source changed during the experiment.'


if __name__ == '__main__':
    main()
