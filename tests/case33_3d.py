"""Standalone three-dimensional, two-stage planning experiment (kW coordinates)."""
import argparse
from dataclasses import dataclass
from itertools import product
import json
from pathlib import Path
from time import perf_counter

import gurobipy as gp
from gurobipy import GRB
import numpy as np
from scipy.spatial import ConvexHull
from threadpoolctl import threadpool_limits

from Network.case33bw import Case33
from model import PlanningEquations, PlanningModel, PLANNING_TOL, evaluation_bounds
from region import clip_polytope, polytope_vertices, polytope_volume
from vertify import ACPowerFlow
from tests.planning_checks import margin
from tests.case33_two_stage import add_tree_relaxation, save, serial, fingerprints

PAD_KW = 1e-4
MASKS = np.asarray(list(product((0., 1.), repeat=3)))


class PlanningOracle3D:
    """Original physical constraints; support or minimum load-distance objective."""

    def __init__(self, budget, threads=8, seed=0, strengthen=True):
        self.network = Case33(upgrade_count=8)
        self.equations = PlanningEquations(self.network, 'socp')
        self.problem = PlanningModel(self.equations, budget=budget, threads=threads)
        if strengthen:
            add_tree_relaxation(self.problem)
        self.budget = budget
        self.bounds = evaluation_bounds(self.network, threads=threads)
        self.square_kw = self.network.power_limit
        model = self.problem.model
        self.distance_kw = model.addVar(ub=0., name='load_distance_kw')
        self.projected_deficit_kw = model.addVar(ub=0., name='projected_deficit_kw')
        self.projected_rows = []
        self.distance_rows = [model.addConstr(
            (self.problem.power[i].item()+self.distance_kw)/self.network.base >= 0.) for i in range(3)]
        model.Params.NumericFocus = 2
        model.Params.Seed = seed
        model.Params.MIPGap = 0.
        model.update()
        self.x_indices = [v.index for v in self.problem.x.tolist()]
        self.p_indices = [v.index for v in self.problem.power.tolist()]
        self.state_indices = [v.index for v in self.problem.state.tolist()]
        self.certificates, self.events, self.trees = [], [], {}
        self.certificate_keys = set()

    def accept(self, values, max_vio, origin):
        x = np.rint(values[self.x_indices]).astype(int)
        p, state = values[self.p_indices], values[self.state_indices]
        key = tuple(x)
        if key not in self.trees:
            self.trees[key] = ACPowerFlow(self.network.tree(x))
        raw_margin = margin(self.equations, x, p, state)
        ac_status = int(self.trees[key].classify(p)[0])
        feasible = (max_vio <= PLANNING_TOL and raw_margin >= -PLANNING_TOL
                    and ac_status == 1 and self.network.cost @ x <= self.budget+1e-8)
        certificate = dict(x=x, p=p, state=state, max_vio=max_vio,
                           raw_margin=raw_margin, ac_status=ac_status, origin=origin)
        stamp = key+tuple(np.round(p, 6))
        if feasible and stamp not in self.certificate_keys:
            self.certificates.append(certificate)
            self.certificate_keys.add(stamp)
        return feasible, certificate

    def solve(self, weights=None, target=None, fixed=None, time_limit=60.,
              gap_kw=.5, label='query', axis_only=False, warm_count=4, stop_distance_kw=None,
              cut_normals=None):
        start = perf_counter()
        problem, net = self.problem, self.network
        model = problem.model
        mode = 'projected_distance' if cut_normals is not None else ('distance' if target is not None else 'support')
        weights = np.zeros(3) if weights is None else np.asarray(weights, dtype=float)
        target = None if target is None else np.asarray(target, dtype=float)
        problem.power.LB, problem.power.UB = np.zeros(3), self.bounds
        self.distance_kw.UB = self.square_kw if mode == 'distance' else 0.
        self.projected_deficit_kw.UB = self.square_kw if mode == 'projected_distance' else 0.
        model.remove(self.projected_rows)
        self.projected_rows = []
        for i, row in enumerate(self.distance_rows):
            row.RHS = target[i]/net.base if mode == 'distance' else 0.
        if cut_normals is not None:
            cut_normals = np.asarray(cut_normals)
            for normal in cut_normals:
                self.projected_rows.append(model.addConstr((gp.quicksum(float(normal[i])*problem.power[i].item()
                    for i in range(3))+self.projected_deficit_kw)/net.base >= float(normal @ target)/net.base))
        if axis_only:
            problem.power.UB = np.where(weights > 0., self.bounds, 0.)
        for i, value in (fixed or {}).items():
            problem.power[i].LB = problem.power[i].UB = value
        deficit = self.projected_deficit_kw if cut_normals is not None else self.distance_kw
        objective = (deficit/net.base if target is not None else
                     gp.quicksum(float(weights[i])*problem.power[i].item()/net.base for i in range(3)))
        model.setObjective(objective, GRB.MINIMIZE if target is not None else GRB.MAXIMIZE)
        model.Params.MIPGapAbs = gap_kw/net.base
        model.Params.SolutionLimit = GRB.MAXINT
        model.update()
        score = (lambda c: np.max(np.maximum(target-c['p'], 0.))) if target is not None else (
            lambda c: -float(weights @ c['p']))
        if fixed:
            previous_score = score
            score = lambda c: (sum(abs(c['p'][i]-v) for i, v in fixed.items()), previous_score(c))
        plans = []
        for certificate in sorted(self.certificates, key=score):
            key = tuple(certificate['x'])
            if key not in plans:
                plans.append(key)
            if len(plans) >= warm_count:
                break
        best_values, best_objective = None, np.inf if target is not None else -np.inf
        best_certificate = None
        for key in plans:
            warm = model.copy()
            variables = warm.getVars()
            for index, value in zip(self.x_indices, key):
                variables[index].LB = variables[index].UB = value
            warm.Params.Threads, warm.Params.TimeLimit = 1, 2.
            warm.optimize()
            if warm.SolCount:
                values = np.asarray(warm.getAttr('X', variables))
                accepted, warm_certificate = self.accept(values, float(warm.MaxVio), 'fixed_plan_start')
                obj = float(warm.ObjVal*net.base)
                better = obj < best_objective if target is not None else obj > best_objective
                if accepted and better:
                    best_values, best_objective = values, obj
                    best_certificate = warm_certificate
            warm.dispose()
            if stop_distance_kw is not None and best_values is not None and np.max(target-best_certificate['p']) <= stop_distance_kw:
                break
        if stop_distance_kw is not None and best_values is not None and np.max(target-best_certificate['p']) <= stop_distance_kw:
            # The accepted warm point is a feasible certificate only. No global
            # lower bound is claimed, and this branch cannot generate a cut.
            answer = dict(mode=mode, label=label, weights=weights, target=target, fixed=fixed, cut_normals=cut_normals,
                          status=None, global_solve=False, bound=None, objective=best_objective,
                          feasible=True, infeasible=False, seconds=perf_counter()-start, nodes=0.,
                          **best_certificate)
            self.events.append({k: v for k, v in answer.items() if k != 'state'})
            print(json.dumps(serial({k: answer[k] for k in ('label','global_solve','objective','seconds')})), flush=True)
            return answer
        if best_values is not None:
            model.setAttr('Start', model.getVars(), best_values)
        model.Params.TimeLimit = max(.1, time_limit-(perf_counter()-start))
        model.optimize()
        answer = dict(mode=mode, label=label, weights=weights, target=target, fixed=fixed, cut_normals=cut_normals,
                      status=int(model.Status), global_solve=True, bound=None, objective=None, feasible=False,
                      infeasible=model.Status == GRB.INFEASIBLE, x=None, p=None, state=None,
                      seconds=perf_counter()-start, nodes=float(model.NodeCount))
        if not answer['infeasible']:
            if np.isfinite(model.ObjBound):
                answer['bound'] = float(model.ObjBound*net.base +
                                        (-PAD_KW if target is not None else PAD_KW))
            if model.SolCount:
                values = np.asarray(model.getAttr('X', model.getVars()))
                feasible, certificate = self.accept(values, float(model.MaxVio), label)
                answer.update(certificate, feasible=feasible, objective=float(model.ObjVal*net.base))
        self.events.append({k: v for k, v in answer.items() if k != 'state'})
        print(json.dumps(serial({k: answer[k] for k in ('label', 'status', 'bound', 'objective',
                                                        'feasible', 'seconds')})), flush=True)
        return answer

    def add_support(self, weights, bound):
        self.problem.model.addConstr(gp.quicksum(float(w)*self.problem.power[i].item()
                                                for i, w in enumerate(weights)) <= bound)

    def close(self):
        self.problem.model.dispose()


def downward_facets(points):
    """Conservative nonnegative-normal facets of one downward convex hull."""
    points = np.maximum(np.asarray(points).reshape(-1, 3)-1e-5, 0.)
    cloud = np.unique((points[:, None, :]*MASKS).reshape(-1, 3), axis=0)
    active = np.max(cloud, axis=0) > 1e-8
    rows = []
    if active.sum() >= 2:
        equations = ConvexHull(cloud[:, active]).equations
        for equation in equations:
            if np.min(equation[:-1]) < -1e-9 or equation[:-1].sum() <= 1e-10:
                continue
            weights = np.zeros(3)
            weights[active] = np.maximum(equation[:-1], 0.)
            weights /= weights.sum()
            rows.append(np.r_[weights, -max(0., np.max(cloud @ weights)-1e-6)])
    for axis in np.where(~active)[0] if active.sum() >= 2 else range(3):
        weights = np.eye(3)[axis]
        rows.append(np.r_[weights, -cloud[:, axis].max()])
    rows = np.asarray(rows).reshape(-1, 4)
    _, indices = np.unique(np.round(rows, 9), axis=0, return_index=True)
    return rows[np.sort(indices)]


def scheme_hulls(certificates):
    groups = {}
    for certificate in certificates:
        groups.setdefault(tuple(certificate['x']), []).append(certificate['p'])
    return [(key, downward_facets(points)) for key, points in groups.items()]


def hull_distance(points, facets):
    """Exact L-infinity distance to a downward convex set, including axis faces."""
    points = np.asarray(points).reshape(-1, 3)
    weights, bound = facets[:, :3], -facets[:, 3]
    remaining = np.broadcast_to(weights, (len(points), *weights.shape)).copy()
    distance = np.zeros(remaining.shape[:2])
    for _ in range(3):
        total = remaining.sum(axis=2)
        numerator = np.einsum('pfi,pi->pf', remaining, points)-bound
        distance = np.maximum(0., np.divide(numerator, total, out=distance.copy(), where=total > 0.))
        remaining *= points[:, None, :] > distance[:, :, None]
    return distance.max(axis=1)


@dataclass
class Cell:
    vertices: np.ndarray
    depth: int = 0
    gap: float = np.inf
    scheme: tuple = ()


def cell_gap(vertices, hulls):
    """min_scheme max_vertex distance bounds every point of a convex cell."""
    distances = np.array([hull_distance(vertices, facets) for _, facets in hulls])
    worst = distances.max(axis=1)
    best = int(worst.argmin())
    point_distances = distances.min(axis=0)
    return float(worst[best]+1e-5), hulls[best][0], point_distances


def select_cell(cells, hulls):
    """Lazy exact maximum: old bounds stay valid while the inner union grows."""
    while True:
        chosen = int(np.argmax([cell.gap for cell in cells]))
        cell = cells[chosen]
        cell.gap, cell.scheme, point_distances = cell_gap(cell.vertices, hulls)
        if cell.gap >= max(other.gap for other in cells)-1e-9:
            return chosen, point_distances


def subtract_orthant(vertices, threshold):
    """Disjoint-interior decomposition of P intersect (OR_i p_i <= threshold_i)."""
    return subtract_disjunction(vertices, threshold, np.eye(3))


def subtract_disjunction(vertices, threshold, cut_normals):
    """Retain OR_k normal[k] @ p <= threshold[k], with disjoint interiors."""
    remaining, pieces = np.asarray(vertices), []
    for normal, bound in zip(cut_normals, threshold):
        retained = clip_polytope(remaining, bound, -normal)
        if len(retained):
            pieces.append(retained)
        remaining = clip_polytope(remaining, -bound, normal)
        if not len(remaining):
            break
    return pieces


def apply_orthant(cells, threshold):
    """A globally proved orthant exclusion tightens every intersected node."""
    return apply_disjunction(cells, threshold, np.eye(3))


def apply_disjunction(cells, threshold, cut_normals):
    result = []
    for cell in cells:
        if np.any((cell.vertices @ cut_normals.T).max(axis=0) <= threshold+1e-9):
            result.append(cell)
        else:
            result.extend(Cell(v, cell.depth+1, cell.gap, cell.scheme)
                          for v in subtract_disjunction(cell.vertices, threshold, cut_normals))
    return result


def disjunction_directions(target, hulls):
    candidates = []
    for _, rows in hulls:
        residuals = rows[:, :3] @ target+rows[:, 3]
        face = int(residuals.argmax())
        candidates.append((residuals[face], rows[face, :3]))
    directions = []
    for _, normal in sorted(candidates, key=lambda pair: pair[0]):
        if all(np.linalg.norm(normal-old) > .04 for old in directions):
            directions.append(normal)
        if len(directions) == 3:
            break
    return np.asarray(directions)


def cover_partition(vertices, facets, epsilon_kw):
    """Partition along a convex inner set's distance neighborhood, losing no volume."""
    normals = (facets[:,None,:3]*MASKS[None,1:,:]).reshape(-1,3)
    scale = normals.sum(axis=1)
    bounds = np.repeat(-facets[:,3],len(MASKS)-1)+epsilon_kw*scale
    active = scale > 1e-12
    normals,bounds = normals[active]/scale[active,None],bounds[active]/scale[active]
    residuals = vertices @ normals.T-bounds
    if np.any(residuals.min(axis=0) > 1e-9):
        return None
    indices = np.where(residuals.max(axis=0) > 1e-9)[0]
    indices = indices[np.argsort(-residuals.max(axis=0)[indices])]
    remaining,pieces = vertices,[]
    for index in indices:
        normal,bound = normals[index],bounds[index]
        if np.max(remaining @ normal-bound) <= 1e-9:
            continue
        outside = clip_polytope(remaining,-bound,normal)
        remaining = clip_polytope(remaining,bound,-normal)
        if len(outside):
            pieces.append(outside)
        if not len(remaining):
            return None
    # A tangent face gives the unchanged outside node again; it is not progress.
    if np.linalg.matrix_rank(remaining-remaining[0], tol=1e-7) < np.linalg.matrix_rank(vertices-vertices[0], tol=1e-7):
        return None
    return pieces+[remaining] if pieces else None


def support_direction(vertices, certificates, excluded=()):
    facets = downward_facets([c['p'] for c in certificates])
    if len(excluded):
        facets = facets[[all(np.max(np.abs(row[:3]-old)) > .02 for old in excluded) for row in facets]]
    if not len(facets):
        return 0., None
    gaps = vertices @ facets[:, :3].T+facets[:, 3]
    index = int(np.max(gaps, axis=0).argmax())
    return float(gaps[:, index].max()), facets[index, :3]


def construct(output, budget, epsilon_kw=20., threads=8, time_limit=60.,
              max_seconds=21600., resume=False, oblique=False):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    previous = json.loads((output/'result.json').read_text(encoding='utf-8')) if resume else None
    if not resume and (output/'result.json').exists():
        raise FileExistsError(output/'result.json')
    start = perf_counter()
    oracle = PlanningOracle3D(budget, threads=threads)
    supports, branch_cuts, history, cells = [], [], [], []
    stage1 = MASKS*oracle.square_kw
    stage1_seconds, stage1_count = 0., 0
    if previous:
        assert previous['source_hashes'] == fingerprints()
        assert previous['budget'] == budget and previous['epsilon_kw'] == epsilon_kw
        start -= previous['seconds']
        supports, branch_cuts, history = previous['supports'], previous['branch_cuts'], previous['history']
        oracle.events = previous['events']
        oracle.certificates = [{**c, **{k: np.asarray(c[k]) for k in ('x', 'p', 'state')}}
                               for c in previous['certificates']]
        stage1 = np.asarray(previous['stage1_vertices'])
        stage1_seconds, stage1_count = previous['stage1_seconds'], previous['stage1_count']
        cells = [Cell(np.asarray(c['vertices']), c['depth'], c['gap'], tuple(c['scheme'])) for c in previous['cells']]
        for cut in supports:
            oracle.add_support(cut['weights'], cut['bound'])

    def checkpoint(status, gap):
        value = dict(schema_version=3, algorithm_revision=7, budget=budget, epsilon_kw=epsilon_kw,
            oblique=oblique,
            square_kw=oracle.square_kw, bounds=oracle.bounds, load_nodes=oracle.network.load_nodes,
            status=status, max_gap_kw=gap, seconds=perf_counter()-start,
            stage1_seconds=stage1_seconds, stage1_count=stage1_count,
            supports=supports, branch_cuts=branch_cuts, stage1_vertices=stage1,
            cells=[dict(vertices=c.vertices, depth=c.depth, gap=c.gap, scheme=c.scheme) for c in cells],
            history=history, certificates=oracle.certificates, events=oracle.events,
            source_hashes=fingerprints(), threads=threads,
            stage1_volume=polytope_volume(stage1),
            stage2_volume=sum(polytope_volume(c.vertices) for c in cells))
        save(output/'result.json', value)
        return value

    def insert_support(weights, bound, label, event=None):
        nonlocal stage1
        supports.append(dict(weights=np.asarray(weights), bound=bound, label=label, event=event))
        stage1 = clip_polytope(stage1, bound, -np.asarray(weights))
        oracle.add_support(weights, bound)

    if not previous:
        insert_support(np.ones(3), oracle.square_kw, 'total_load')
        for i in range(3):
            insert_support(np.eye(3)[i], oracle.bounds[i], f'necessary_axis_{i}')
        for i in range(3):
            answer = oracle.solve(np.eye(3)[i], axis_only=True, time_limit=time_limit,
                                  label=f'axis_{i}')
            if answer['bound'] is not None:
                insert_support(np.eye(3)[i], answer['bound'], f'axis_{i}', len(oracle.events)-1)
    if not previous or previous['status'] == 'stage1':
        assert oracle.certificates, 'No feasible origin or load certificate.'
        excluded = [np.asarray(event['weights']) for event in oracle.events
                    if event['mode'] == 'support' and event['bound'] is not None
                    and event['objective'] is not None and event['bound']-event['objective'] > epsilon_kw/2]
        for iteration in range(max(0, stage1_count-3), 40):
            gap, weights = support_direction(stage1, oracle.certificates, excluded)
            if gap <= epsilon_kw/2:
                break
            answer = oracle.solve(weights, time_limit=time_limit, label=f'support_{iteration}')
            if answer['bound'] is not None:
                insert_support(weights, answer['bound'], f'support_{iteration}', len(oracle.events)-1)
            if answer['objective'] is None or answer['bound'] is None or answer['bound']-answer['objective'] > epsilon_kw/2:
                excluded.append(weights)
            stage1_seconds = perf_counter()-start
            stage1_count = len(oracle.events)
            checkpoint('stage1', gap)
        stage1_seconds, stage1_count = perf_counter()-start, len(oracle.events)
        cells = [Cell(stage1)]
    if not cells:
        cells = [Cell(stage1)]
    if previous:
        for cut in branch_cuts:
            cells = apply_disjunction(cells, np.asarray(cut['threshold']), np.asarray(cut.get('cut_normals', np.eye(3))))
    hulls = scheme_hulls(oracle.certificates)
    revision = -1
    iteration = len(history)
    while perf_counter()-start < max_seconds:
        if revision != len(oracle.certificates):
            hulls = scheme_hulls(oracle.certificates)
            revision = len(oracle.certificates)
        chosen, point_distances = select_cell(cells, hulls)
        cell = cells[chosen]
        if cell.gap <= epsilon_kw:
            result = checkpoint('certified', cell.gap)
            oracle.close()
            return result
        record = dict(iteration=iteration, selected_gap_kw=cell.gap,
                      cell_count=len(cells), certificate_count=len(oracle.certificates))
        if point_distances.max() > epsilon_kw:
            target = cell.vertices[int(point_distances.argmax())]
            cut_normals = disjunction_directions(target, hulls) if oblique else None
            answer = oracle.solve(target=target, time_limit=time_limit,
                                  label=f'branch_{len(branch_cuts)}', warm_count=8,
                                  gap_kw=epsilon_kw/5, stop_distance_kw=epsilon_kw-.001,
                                  cut_normals=cut_normals)
            if cut_normals is not None and answer['global_solve'] and (answer['bound'] is None or answer['bound'] <= .01):
                cut_normals = None
                answer = oracle.solve(target=target, time_limit=time_limit, label=f'branch_{len(branch_cuts)}_orthant',
                    warm_count=8, gap_kw=epsilon_kw/5, stop_distance_kw=epsilon_kw-.001)
            if answer['bound'] is not None and answer['bound'] > .01:
                normals = np.eye(3) if cut_normals is None else cut_normals
                threshold = normals @ target-answer['bound']
                branch_cuts.append(dict(threshold=threshold, target=target,
                                        bound=answer['bound'], event=len(oracle.events)-1, cut_normals=normals))
                cells = apply_disjunction(cells, threshold, normals)
                record.update(action='orthant_cut' if cut_normals is None else 'oblique_cut', cut=len(branch_cuts)-1, target=target)
            else:
                record.update(action='inner_refinement', target=target)
                if not answer['feasible']:
                    history.append(record)
                    result = checkpoint('unresolved_query', cell.gap)
                    oracle.close()
                    return result
        else:
            choices = sorted(hulls,key=lambda pair: -np.sum(hull_distance(cell.vertices,pair[1]) <= epsilon_kw-2e-5))
            children = None
            for key,facets in choices[:3]:
                children = cover_partition(cell.vertices,facets,epsilon_kw-2e-5)
                if children is not None:
                    break
            if children is not None:
                following = [Cell(v,cell.depth+1,cell.gap,cell.scheme) for v in children]
                gap = float(hull_distance(children[-1],facets).max()+1e-5)
                if gap < following[-1].gap:
                    following[-1].gap,following[-1].scheme = gap,key
                record.update(action='coverage_split',piece_count=len(children))
            else:
                axis = int(np.ptp(cell.vertices, axis=0).argmax())
                middle = (cell.vertices[:, axis].min()+cell.vertices[:, axis].max())/2
                children = [clip_polytope(cell.vertices, middle, -np.eye(3)[axis]),
                            clip_polytope(cell.vertices, -middle, np.eye(3)[axis])]
                following = [Cell(v,cell.depth+1,cell.gap,cell.scheme) for v in children if len(v)]
                record.update(action='geometric_split',axis=axis,threshold=middle)
            cells[chosen:chosen+1] = following
        history.append(record)
        iteration += 1
        if iteration % 10 == 0 or record['action'] not in ('geometric_split','coverage_split'):
            checkpoint('stage2', max(c.gap for c in cells))
            print(f'GEOMETRY budget={budget:g} step={iteration} cells={len(cells)} gap={cell.gap:.3f}', flush=True)
    result = checkpoint('time_limit', max(c.gap for c in cells))
    oracle.close()
    return result


def scan_slice(output, budget, fixed_value, spacing_kw=10., threads=4,
               time_limit=90., resume=False, start_index=0, stop_index=None):
    """Independent pointwise AC grid reference, with global column bounds."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    previous = json.loads((output/'scan.json').read_text(encoding='utf-8')) if resume else None
    if not resume and (output/'scan.json').exists():
        raise FileExistsError(output/'scan.json')
    start = perf_counter()
    oracle = PlanningOracle3D(budget, threads=threads, seed=71)
    coordinates = np.r_[np.arange(0., oracle.square_kw, spacing_kw), oracle.square_kw]
    states = np.zeros((len(coordinates), len(coordinates)), dtype=np.int8)
    states[coordinates > oracle.bounds[0], :] = -1
    states[:, coordinates > oracle.bounds[1]] = -1
    states[coordinates[:, None]+coordinates[None, :]+fixed_value > oracle.square_kw] = -1
    columns, upper = [], float(oracle.bounds[1])
    if previous:
        assert previous['source_hashes'] == fingerprints()
        assert previous['fixed_value'] == fixed_value and previous['budget'] == budget
        states = np.load(output/'scan_grid.npz')['states']
        columns, upper = previous['columns'], previous['columns'][-1]['upper']
        start_index = max(start_index, columns[-1]['index']+1)
        oracle.certificates = [{**c, **{k: np.asarray(c[k]) for k in ('x', 'p', 'state')}}
                               for c in previous['certificates']]
        for c in oracle.certificates:
            oracle.trees[tuple(c['x'])] = ACPowerFlow(oracle.network.tree(c['x']))
        oracle.events = previous['events']
        start -= previous['seconds']

    def checkpoint():
        np.savez_compressed(output/'scan_grid.npz', coordinates=coordinates, states=states)
        result = dict(schema_version=3, budget=budget, fixed_value=fixed_value, spacing_kw=spacing_kw,
                      independent=True, source_hashes=fingerprints(), columns=columns,
                      events=oracle.events, certificates=oracle.certificates,
                      seconds=perf_counter()-start, threads=threads,
                      counts={str(k): int(np.sum(states == k)) for k in (-1, 0, 1)},
                      status='classified_grid' if np.all(states != 0) else 'unknown_grid_points')
        save(output/'scan.json', result)
        return result

    for i, threshold in enumerate(coordinates):
        if i < start_index:
            continue
        if stop_index is not None and i >= stop_index:
            break
        if threshold > oracle.bounds[0]:
            break
        states[i:, coordinates > upper] = -1
        candidates = np.where(states[i] == 0)[0]
        if not len(candidates):
            break
        # Earlier column proofs remain valid as p18 increases at fixed p33.
        upper = min(upper, float(coordinates[candidates[-1]]+spacing_kw))
        oracle.add_support((0.,1.,0.), upper)
        trial = np.array([threshold, coordinates[candidates[-1]], fixed_value])
        witness = next((key for key, tree in oracle.trees.items() if tree.classify(trial)[0] == 1), None)
        if witness is None:
            answer = oracle.solve(weights=(0, 1, 0), fixed={0: threshold, 2: fixed_value},
                                  time_limit=time_limit, gap_kw=spacing_kw*.7, label=f'column_{i}')
            if answer['infeasible']:
                states[i:, :] = -1
                columns.append(dict(index=i, upper=-1., lower=None, witness_x=None))
                checkpoint()
                break
            if answer['bound'] is not None:
                upper = min(upper, answer['bound'])
            states[i:, coordinates > upper] = -1
            eligible = sorted(oracle.certificates,
                              key=lambda c: -c['p'][1] if c['p'][0] >= threshold-1e-7 and c['p'][2] >= fixed_value-1e-7 else np.inf)
            witness = next((tuple(c['x']) for c in eligible if c['p'][0] >= threshold-1e-7 and
                            c['p'][2] >= fixed_value-1e-7 and oracle.trees[tuple(c['x'])].classify(
                                [threshold, min(c['p'][1], upper), fixed_value])[0] == 1), None)
        lower = None
        # Scan-owned plans only: AC-test the column to obtain the strongest
        # feasible bracket, even when its previous certificate was at smaller p18.
        candidates = np.where(states[i] == 0)[0]
        batch = np.column_stack((np.full(len(candidates), threshold), coordinates[candidates],
                                 np.full(len(candidates), fixed_value)))
        best_grid_index = -1
        for key, tree in oracle.trees.items():
            feasible_indices = np.where(tree.classify(batch) == 1)[0]
            if len(feasible_indices) and feasible_indices[-1] > best_grid_index:
                best_grid_index = int(feasible_indices[-1])
                witness = key
        if witness is not None:
            feasible = oracle.trees[witness].classify(batch) == 1
            states[i, candidates[feasible]] = 1
        attempted = set()
        while True:
            pending = [j for j in np.where(states[i] == 0)[0] if j not in attempted]
            if not pending:
                break
            j = pending[len(pending)//2]
            attempted.add(j)
            answer = oracle.solve(fixed={0: threshold, 1: coordinates[j], 2: fixed_value},
                                  time_limit=time_limit, label=f'point_{i}_{j}')
            if answer['infeasible']:
                states[i:, j:] = -1
            elif answer['feasible']:
                witness = tuple(answer['x'])
                candidates = np.where(states[i] == 0)[0]
                candidates = candidates[coordinates[candidates] <= coordinates[j]]
                batch = np.column_stack((np.full(len(candidates), threshold), coordinates[candidates],
                                         np.full(len(candidates), fixed_value)))
                ac = oracle.trees[witness].classify(batch)
                states[i, candidates[ac == 1]] = 1
        if np.any(states[i] == 1):
            lower = float(coordinates[np.where(states[i] == 1)[0][-1]])
        columns.append(dict(index=i, threshold=threshold, lower=lower, upper=upper, witness_x=witness))
        checkpoint()
        print(f'SLICE budget={budget:g} p33={fixed_value:g} p18={threshold:g} lower={lower} upper={upper:.3f}', flush=True)
    result = checkpoint()
    oracle.close()
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=('probe', 'run', 'scan'), required=True)
    parser.add_argument('--budget', type=float, default=2.)
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--time-limit', type=float, default=60.)
    parser.add_argument('--max-seconds', type=float, default=21600.)
    parser.add_argument('--epsilon-kw', type=float, default=20.)
    parser.add_argument('--fixed-value', type=float, default=0.)
    parser.add_argument('--spacing-kw', type=float, default=10.)
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--oblique', action='store_true')
    args = parser.parse_args()
    with threadpool_limits(limits=1):
        if args.mode == 'run':
            construct(args.output, args.budget, args.epsilon_kw, args.threads, args.time_limit,
                      args.max_seconds, args.resume, args.oblique)
        elif args.mode == 'scan':
            scan_slice(args.output, args.budget, args.fixed_value, args.spacing_kw, args.threads,
                       args.time_limit, args.resume)
        else:
            oracle = PlanningOracle3D(args.budget, threads=args.threads)
            for axis in range(3):
                oracle.solve(np.eye(3)[axis], axis_only=True, time_limit=args.time_limit,
                             label=f'axis_{axis}')
            oracle.solve(target=np.array([900., 3000., 900.]), time_limit=args.time_limit, label='distance_probe')
            save(args.output/'probe.json', dict(budget=args.budget, events=oracle.events,
                 certificates=oracle.certificates, source_hashes=fingerprints()))
            oracle.close()


if __name__ == '__main__':
    main()
