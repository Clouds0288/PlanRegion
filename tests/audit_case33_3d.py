"""Offline physical, whole-volume and independent-slice certificate audit."""
import argparse
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.spatial import ConvexHull
from shapely.geometry import Polygon, MultiPoint, box
from shapely.ops import unary_union
from shapely import covers,points as shapely_points
from threadpoolctl import threadpool_limits

from Network.case33bw import Case33
from model import PlanningEquations, PLANNING_TOL
from region import clip_polytope, polytope_volume
from vertify import ACPowerFlow
from tests.planning_checks import margin
from tests.case33_two_stage import save, fingerprints
from tests.case33_3d import MASKS, scheme_hulls, downward_facets, cell_gap, hull_distance, subtract_disjunction
from tests.audit_case33_two_stage import directed_distance


def final_cells(result):
    """Apply every globally valid stage-2 exclusion also to already covered cells."""
    cells = [np.asarray(c['vertices']) for c in result['cells']]
    for cut in result['branch_cuts']:
        threshold = np.asarray(cut['threshold'])
        normals = np.asarray(cut.get('cut_normals', np.eye(3)))
        refined = []
        for vertices in cells:
            if np.any((vertices @ normals.T).max(axis=0) <= threshold+1e-8):
                refined.append(vertices)
            else:
                refined.extend(subtract_disjunction(vertices, threshold, normals))
        cells = refined
    return cells


def section(vertices, fixed_value):
    vertices = clip_polytope(vertices, fixed_value, [0, 0, -1])
    vertices = clip_polytope(vertices, -fixed_value, [0, 0, 1])
    return MultiPoint(vertices[:, :2]).convex_hull if len(vertices) else Polygon()


def actual_region(grid):
    coordinates, states = grid['coordinates'], grid['states']
    return unary_union([box(0, 0, coordinates[i], coordinates[np.where(row == 1)[0][-1]])
                        for i, row in enumerate(states) if np.any(row == 1) and coordinates[i] > 0])


def infeasibility_mask(scan, coordinates, bounds, total_bound):
    """Reconstruct every rejected grid point from independent global evidence."""
    first, second = coordinates[:,None], coordinates[None,:]
    proved = (first > bounds[0]) | (second > bounds[1]) | (first+second+scan['fixed_value'] > total_bound)
    for event in scan['events']:
        if not event.get('global_solve',True):
            continue
        if scan.get('reused'):
            if event['label'].startswith('scan_point_'):
                i,j = map(int,event['label'].split('_')[-2:])
                p = np.array([i,j])*scan['provenance']['source_spacing_kw']
                if event['infeasible']:
                    proved |= (first >= p[0]) & (second >= p[1])
                continue
            threshold = event['threshold']
        else:
            fixed = event['fixed']
            assert abs(fixed['2']-scan['fixed_value']) < 1e-8
            threshold = fixed['0']
            if '1' in fixed:
                if event['infeasible']:
                    proved |= (first >= threshold) & (second >= fixed['1'])
                continue
        if event['infeasible']:
            proved |= first >= threshold
        elif event['bound'] is not None:
            proved |= (first >= threshold) & (second > event['bound'])
    for check in scan.get('original_model_checks',[]):
        if check['classification'] == -1 and check['answer']['infeasible']:
            assert check['source_hashes'] == scan['source_hashes']
            p = check['p']
            proved |= (first >= p[0]) & (second >= p[1])
    return proved


def audit(directory, require_slices=True):
    directory = Path(directory)
    result = json.loads((directory/'result.json').read_text(encoding='utf-8'))
    start = perf_counter()
    net = Case33(upgrade_count=8)
    equations = PlanningEquations(net, 'socp')
    trees = {}
    min_margin, min_support_slack, min_cut_slack = np.inf, np.inf, np.inf
    ac_pass = 0
    for c in result['certificates']:
        x, p, state = (np.asarray(c[k]) for k in ('x', 'p', 'state'))
        assert net.cost @ x <= result['budget']+1e-8
        assert np.all((x == 0) | (x == 1))
        tree = net.tree(x)
        assert tree.n == net.n
        key = tuple(x)
        if key not in trees:
            trees[key] = ACPowerFlow(tree)
        ac_pass += int(trees[key].classify(p)[0] == 1)
        min_margin = min(min_margin, margin(equations, x, p, state))
        for cut in result['supports']:
            min_support_slack = min(min_support_slack, cut['bound']-np.asarray(cut['weights']) @ p)
        for cut in result['branch_cuts']:
            normals = np.asarray(cut.get('cut_normals', np.eye(3)))
            min_cut_slack = min(min_cut_slack, np.max(np.asarray(cut['threshold'])-normals @ p))
    hulls = scheme_hulls(result['certificates'])
    max_gap, recorded_gap_violation = 0., 0.
    for cell in result['cells']:
        gap, _, _ = cell_gap(np.asarray(cell['vertices']), hulls)
        max_gap = max(max_gap, gap)
        recorded_gap_violation = max(recorded_gap_violation, gap-cell['gap'])
    for cut in result['branch_cuts']:
        event = result['events'][cut['event']]
        assert event.get('global_solve', True) and event['mode'] in ('distance', 'projected_distance')
        normals = np.asarray(cut.get('cut_normals', np.eye(3)))
        if event['mode'] == 'projected_distance':
            np.testing.assert_allclose(normals, event['cut_normals'], atol=1e-12)
        else:
            np.testing.assert_allclose(normals, np.eye(3), atol=1e-12)
        np.testing.assert_allclose(cut['threshold'], normals @ np.asarray(event['target'])-event['bound'], atol=1e-8)
    for cut in result['supports']:
        if cut['event'] is not None:
            event = result['events'][cut['event']]
            assert event.get('global_solve', True) and event['mode'] == 'support'
            assert abs(cut['bound']-event['bound']) < 1e-8
    cells = final_cells(result)
    upper_volume = sum(polytope_volume(v) for v in cells)
    # One common translation preserves disjoint interiors. Every translated
    # nonnegative point is in the downward inner union by the distance proof.
    inner_volume_lower = 0.
    for vertices in cells:
        translated = vertices-max_gap
        for axis in range(3):
            translated = clip_polytope(translated, 0., np.eye(3)[axis])
        inner_volume_lower += polytope_volume(translated)
    convex_facets = downward_facets([c['p'] for c in result['certificates']])
    convex_inner = MASKS*result['square_kw']
    for row in convex_facets:
        convex_inner = clip_polytope(convex_inner, -row[3], -row[:3])
    retained_convex_volume = 0.
    for vertices in cells:
        crossing = convex_facets[np.max(vertices @ convex_facets[:, :3].T+convex_facets[:,3], axis=0) > 1e-9]
        clipped = vertices
        for row in crossing:
            clipped = clip_polytope(clipped, -row[3], -row[:3])
            if not len(clipped):
                break
        retained_convex_volume += polytope_volume(clipped)
    nonconvex_volume_lower = max(0., polytope_volume(convex_inner)-retained_convex_volume)
    convex_support_gap = float(hull_distance(result['stage1_vertices'], convex_facets).max())
    slices, slice_pass = [], True
    for fixed_value in (0., 500., 1000.):
        folder = directory/f'slice_{fixed_value:g}'
        if not (folder/'scan.json').exists():
            if require_slices:
                slice_pass = False
            continue
        scan = json.loads((folder/'scan.json').read_text(encoding='utf-8'))
        grid = np.load(folder/'scan_grid.npz')
        coordinates, states = grid['coordinates'], grid['states']
        reference = actual_region(grid)
        stage1 = section(np.asarray(result['stage1_vertices']), fixed_value)
        final = unary_union([section(v, fixed_value) for v in cells])
        outside = reference.difference(final.buffer(1e-5)).area
        reference_points = coordinates[np.argwhere(states == 1)]
        reference_points_outside = int(np.sum(~covers(final.buffer(1e-5),shapely_points(reference_points))))
        ac_count, ac_fail = 0, 0
        for column in scan['columns']:
            i = column['index']
            indices = np.where(states[i] == 1)[0]
            if not len(indices):
                continue
            key = tuple(column['witness_x'])
            assert net.cost @ np.asarray(key) <= result['budget']+1e-8
            if key not in trees:
                trees[key] = ACPowerFlow(net.tree(key))
            batch = np.column_stack((np.full(len(indices), coordinates[i]), coordinates[indices],
                                     np.full(len(indices), fixed_value)))
            ac_count += len(indices)
            ac_fail += int(np.sum(trees[key].classify(batch) != 1))
        counts = {str(k): int(np.sum(states == k)) for k in (-1, 0, 1)}
        proved = infeasibility_mask(scan, coordinates, np.asarray(result['bounds']),result['square_kw'])
        infeasible_unproved = int(np.sum((states == -1) & ~proved))
        proof_conflicts = int(np.sum((states == 1) & proved))
        grid_tops = [(coordinates[i], coordinates[np.where(row == 1)[0][-1]])
                     for i, row in enumerate(states) if np.any(row == 1)]
        slice_gap = (directed_distance(final, grid_tops) if grid_tops else
                     (0. if final.is_empty else None))
        passed = (counts['0'] == 0 and ac_count == counts['1'] and ac_fail == 0 and outside < 1e-3
                  and infeasible_unproved == 0 and proof_conflicts == 0
                  and reference_points_outside == 0
                  and scan['independent'] and scan['source_hashes'] == result['source_hashes'])
        slice_pass &= passed
        slices.append(dict(fixed_value=fixed_value, grid_counts=counts, actual_area=reference.area,
                           stage1_area=stage1.area, stage2_area=final.area, actual_outside=outside,
                           outer_to_reference_gap_kw=slice_gap,
                           infeasible_unproved=infeasible_unproved,proof_conflicts=proof_conflicts,
                           reference_points_outside=reference_points_outside,
                           ac_replayed=ac_count, ac_fail=ac_fail, passed=passed, seconds=scan['seconds']))
    checks = dict(production_unchanged=result['source_hashes'] == fingerprints(),
                  certified_status=result['status'] == 'certified',
                  physical_margin=min_margin >= -PLANNING_TOL,
                  all_endpoints_ac=ac_pass == len(result['certificates']),
                  support_cut_validity=min_support_slack >= -1e-3,
                  orthant_cut_validity=min_cut_slack >= -1e-3,
                  whole_volume_gap=max_gap <= result['epsilon_kw'],
                  largest_gap_nonincreasing=bool(np.all(np.diff([h['selected_gap_kw'] for h in result['history']]) <= .001)),
                  reference_slices=slice_pass)
    data = dict(budget=result['budget'], checks=checks, passed=all(checks.values()),
                audit_scope='full' if require_slices else 'construction_with_available_slices',
                source_hashes=result['source_hashes'], min_margin=min_margin,
                min_support_slack_kw=min_support_slack, min_cut_slack_kw=min_cut_slack,
                certificate_count=len(result['certificates']), ac_pass=ac_pass,
                scheme_count=len(hulls), max_gap_kw=max_gap,
                recorded_gap_violation_kw=recorded_gap_violation,
                stage1_volume=result['stage1_volume'], stage2_raw_volume=result['stage2_volume'],
                stage2_volume=upper_volume, inner_volume_lower=inner_volume_lower,
                nonconvex_volume_lower=nonconvex_volume_lower,
                convex_support_gap_upper_kw=convex_support_gap,
                removed_stage2_volume=result['stage1_volume']-upper_volume,
                final_cell_count=len(cells), construction_seconds=result['seconds'],
                stage1_seconds=result['stage1_seconds'], slices=slices,
                audit_seconds=perf_counter()-start)
    save(directory/'audit.json', data)
    save(directory/'final_geometry.json', dict(cells=cells, max_gap_kw=max_gap,
                                              inner_volume_lower=inner_volume_lower,
                                              outer_volume=upper_volume))
    print(json.dumps(data, default=lambda x: x.item() if isinstance(x, np.generic) else x), flush=True)
    return data


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', type=Path)
    parser.add_argument('--without-slices', action='store_true')
    args = parser.parse_args()
    with threadpool_limits(limits=1):
        answer = audit(args.directory, not args.without_slices)
    raise SystemExit(0 if answer['passed'] else 1)
