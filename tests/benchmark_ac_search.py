"""Independent AC audit of the saved physical-search experiment, without re-solving it.

Finite budgets enumerate every affordable design. Infinite-budget positive labels
need an actual AC witness; negative labels need a global all-design search. The
all-upgraded design is only a fast source of witnesses, never an exclusion rule.
"""
from argparse import ArgumentParser
from hashlib import sha256
from pathlib import Path
from time import perf_counter
import json

import numpy as np

from Network.case33bw import Case33
from main import numerical_threads
from model import PlanningEquations
from plot import RunMonitor, save_benchmark_report
from io import StringIO
from region import sample_region, classify_orthant, split_grid_box
from vertify import (ACPowerFlow, AC_TOL, FIXED_POINT_TOL, ac_planning_query,
                     disagreement_interval, validate_power_flow)
from tests.benchmark_physical_search import (ROOT, VARIANTS, LABELS, fingerprints,
                                             write_json, case_name)

DEFAULT = ROOT/'results/case33bw/latest'


def affordable_designs(network, budget):
    """独立参考只枚举预算内方案；按费用剪枝，32 候选预算 2 仅有 498 种。"""
    if not np.isfinite(budget):
        raise ValueError('Unlimited budgets require global search, not enumeration')
    result = []
    def visit(index, cost, chosen):
        if index == len(network.line_options):
            choice = tuple(chosen)
            result.append((cost, choice, ACPowerFlow(network.design(choice))))
            return
        for option, price in enumerate(network.line_options[index].cost):
            if price < 0.:
                raise ValueError('Budget pruning requires nonnegative investment costs')
            if cost+price <= budget:
                visit(index+1, cost+price, chosen+[option])
    visit(0, 0., [])
    return sorted(result, key=lambda row: (row[0], row[1]))


def union_labels(designs, points):
    points = np.atleast_2d(points)
    labels = np.full(len(points), -1, dtype=np.int8)
    calls = 0
    for _, _, oracle in designs:
        pending = np.flatnonzero(labels != 1)
        if not len(pending):
            break
        status = oracle.classify(points[pending])
        calls += len(pending)
        labels[pending[status == 1]] = 1
        labels[pending[status == 0]] = 0
    return labels, calls


def fixed_rays(oracle, rays, steps=24):
    lower, upper = np.zeros(len(rays)), np.ones(len(rays))
    if oracle.classify(np.zeros((1, 3)))[0] != 1:
        raise ValueError('Ray audit requires a feasible zero independent-load point')
    top = oracle.classify(rays)
    lower[top == 1] = 1.
    unresolved = np.zeros(len(rays), dtype=bool)
    for _ in range(steps):
        active = np.flatnonzero((lower < upper) & ~unresolved)
        midpoint = (lower[active]+upper[active])/2
        status = oracle.classify(rays[active]*midpoint[:, None])
        lower[active[status == 1]] = midpoint[status == 1]
        upper[active[status == -1]] = midpoint[status == -1]
        unresolved[active[status == 0]] = True
    return lower, upper, int(unresolved.sum())


def certify_infinite_grid(network, points, divisions, folder):
    started = perf_counter()
    seed = ACPowerFlow(network.design([1]*len(network.projects))).classify(points)
    # Only positive AC certificates are reusable across the union.
    states = np.where(seed == 1, 1, 0).astype(np.int8).reshape((divisions,)*3)
    equation = PlanningEquations(network, 'socp')
    visited, queries, exclusions = set(), [], 0
    pending = [(np.zeros(3, dtype=int), np.full(3, divisions-1, dtype=int))]

    def observer(event, **data):
        nonlocal exclusions
        if event == 'ac_exclude':
            exclusions += 1

    while pending:
        lower, upper = pending.pop()
        block = tuple(slice(int(a), int(b)+1) for a, b in zip(lower, upper))
        if np.all(states[block] != 0):
            continue
        for index in (upper, lower):
            key = tuple(map(int, index))
            if states[key] or key in visited:
                continue
            visited.add(key)
            flat_index = np.ravel_multi_index(key, states.shape)
            point = points[flat_index]
            tick = perf_counter()
            answer = ac_planning_query(equation, point, observer=observer)
            status = -1 if answer is None else (1 if answer['feasible'] else 0)
            classify_orthant(states, index, status)
            queries.append(dict(index=key, point=point, status=status,
                                seconds=perf_counter()-tick,
                                bound=None if answer is None else answer['bound'],
                                choice=None if answer is None or answer['x'] is None
                                else equation.choice(answer['x'])))
            if len(queries) % 25 == 0:
                write_json(folder/'progress.json', dict(stage='infinite_grid',
                    queries=len(queries), unresolved=int((states == 0).sum()),
                    seconds=perf_counter()-started))
                print('INFINITE GRID', network.name, len(network.projects), len(queries),
                      'remaining', int((states == 0).sum()), flush=True)
            if np.all(states[block] != 0):
                break
        if np.any(states[block] == 0) and np.any(upper > lower):
            pending.extend(split_grid_box(lower, upper))
    return states.ravel(), dict(seconds=perf_counter()-started, queries=queries,
        all_upgraded_ac_witnesses=int((seed == 1).sum()), ac_design_exclusions=exclusions,
        unknown=int((states == 0).sum()),
        rule='Seed positive AC witnesses only; all negative labels require all-design global search.')


def build_ac_reference(network, records, divisions, folder):
    started = perf_counter()
    bounds = np.asarray(records[0]['bounds'])
    for record in records:
        np.testing.assert_array_equal(record['bounds'], bounds)
    points = (np.indices((divisions,)*3).reshape(3, -1).T+.5)*bounds/divisions
    designs = affordable_designs(network, 2.)
    grid = {b: np.full(len(points), -1, dtype=np.int8) for b in (0., 1., 2.)}
    for index, (cost, choice, oracle) in enumerate(designs):
        eligible = [b for b in grid if cost <= b]
        active = np.flatnonzero(np.any(np.vstack([grid[b] != 1 for b in eligible]), axis=0))
        labels = oracle.classify(points[active])
        for b in eligible:
            known = grid[b][active] == 1
            grid[b][active[(labels == 0) & ~known]] = 0
            grid[b][active[labels == 1]] = 1
        if (index+1) % 30 == 0:
            print('AC DESIGNS', len(network.projects), index+1, '/', len(designs), flush=True)
            write_json(folder/'progress.json', dict(stage='finite_designs',
                candidate_count=len(network.projects), completed=index+1, total=len(designs),
                seconds=perf_counter()-started))
    finite_seconds = perf_counter()-started
    grid[None], infinite = certify_infinite_grid(network, points, divisions, folder)
    reference = dict(candidate_count=len(network.projects), bounds=bounds, divisions=divisions,
        points=len(points), finite_designs=len(designs), finite_seconds=finite_seconds,
        infinite_grid=infinite, seconds=perf_counter()-started, source_hashes=fingerprints(),
        budgets=[dict(budget=b, designs=sum(c <= b for c, _, _ in designs) if b is not None else None,
            feasible_cells=int((grid[b] == 1).sum()), infeasible_cells=int((grid[b] == -1).sum()),
            unknown_cells=int((grid[b] == 0).sum())) for b in (0., 1., 2., None)])
    np.savez_compressed(folder/'reference_grid.npz', points=points,
        **{f'b{"inf" if b is None else int(b)}': grid[b] for b in grid})
    write_json(folder/'reference.json', reference)
    return reference, grid, designs, points


def converged_diagnostics(oracle, points):
    points = np.atleast_2d(points)
    ell = np.zeros((len(points), oracle.network.n))
    converged = False
    for _ in range(500):
        P, Q, v, u = oracle.state(points, ell)
        residual = np.max(np.abs(P*P+Q*Q-u*ell), axis=1)
        if np.max(residual) <= FIXED_POINT_TOL:
            converged = True
            break
        if np.any(u <= 0):
            break
        ell = (P*P+Q*Q)/u
    c = oracle.network
    return [dict(converged=converged, equation_residual=float(residual[i]),
                 maximum_violation_pu=float(oracle.violation(P, Q, v)[i]),
                 undervoltage_pu=float(max(0., np.max(np.sqrt(c.vmin)-np.sqrt(np.maximum(v[i], 0.))))),
                 minimum_voltage_pu=float(np.sqrt(np.maximum(v[i], 0.)).min()))
            for i in range(len(points))]


def audit_ac_case(network, record, reference, ac_labels, designs, points):
    started = perf_counter()
    budget, result = record['budget'], record['region']
    bounds = np.asarray(record['bounds'])
    boundary, failures, total, accepted, unknown = [], [], 0, 0, 0
    for row in result['certificates']:
        powers = np.asarray(row['inner']).reshape(-1, 3)*bounds
        if not len(powers):
            continue
        oracle = ACPowerFlow(network.design(row['choice']))
        labels = oracle.classify(powers)
        total += len(powers)
        accepted += int((labels == 1).sum())
        unknown += int((labels == 0).sum())
        for i in np.flatnonzero(labels != 1):
            point = powers[i]
            lo, hi, ray_unknown = fixed_rays(oracle, point[None, :])
            detail = dict(point=point, choice=row['choice'], fixed_design_status=int(labels[i]),
                diagnostic=converged_diagnostics(oracle, point)[0],
                inward_retreat_kw=float((1-lo[0])*point.sum()),
                retreat_bracket_kw=float((hi[0]-lo[0])*point.sum()), retreat_unknown=ray_unknown)
            failures.append(detail)
        boundary.append(dict(choice=row['choice'], points=len(powers),
                             feasible=int((labels == 1).sum()), unknown=int((labels == 0).sum())))
    if failures:
        fail_points = np.array([f['point'] for f in failures])
        if budget is not None:
            union, _ = union_labels([d for d in designs if d[0] <= budget], fail_points)
        else:
            union = ACPowerFlow(network.design([1]*len(network.projects))).classify(fail_points)
            equation = PlanningEquations(network, 'socp')
            for i in np.flatnonzero(union != 1):
                answer = ac_planning_query(equation, fail_points[i])
                union[i] = -1 if answer is None else (1 if answer['feasible'] else 0)
        for failure, status in zip(failures, union):
            failure['union_status'] = int(status)
    maximum = None
    if result['max_point'] is not None:
        oracle = ACPowerFlow(network.design(result['max_choice']))
        maximum = dict(point=result['max_point'], choice=result['max_choice'],
            total_kw=result['max_total'], ac_status=int(oracle.classify(result['max_point'])[0]),
            diagnostic=converged_diagnostics(oracle, result['max_point'])[0])
    boundary_seconds = perf_counter()-started
    sampled = sample_region(result, points, bounds)
    uncertainty = disagreement_interval(sampled, ac_labels)
    inside, ac_inside = sampled == 1, ac_labels == 1
    extra = int(np.count_nonzero(inside & (ac_labels == -1)))
    missed = int(np.count_nonzero((sampled == -1) & ac_inside))
    retained_missed = int(np.count_nonzero(~inside & ac_inside))
    volume_error = disagreement_interval(np.where(inside, 1, -1), ac_labels)
    return dict(candidate_count=record['candidate_count'], budget=budget, variant=record['variant'],
        construction_status=record['status'], construction_seconds=record['solve_seconds'],
        strict_vertices=dict(total=total, feasible=accepted, infeasible=total-accepted-unknown,
            unknown=unknown, details=boundary, failures=failures,
            union_infeasible=sum(f['union_status'] == -1 for f in failures),
            union_unknown=sum(f['union_status'] == 0 for f in failures),
            maximum_inward_retreat_kw=max((f['inward_retreat_kw'] for f in failures), default=0.)),
        maximum=maximum, grid=dict(**uncertainty, points=len(points),
            ac_inside=int(ac_inside.sum()), inner_inside=int(inside.sum()),
            definite_extra=extra, definite_missed=missed, ac_not_retained=retained_missed,
            ac_retained_percent=100*int(np.count_nonzero(inside & ac_inside))/int(ac_inside.sum()) if ac_inside.any() else None,
            retained_inner_fr_percent=100*extra/int(inside.sum()) if inside.any() else None,
            ac_in_unknown=int(np.count_nonzero((sampled == 0) & ac_inside)),
            method_unknown=int((sampled == 0).sum()), ac_unknown=int((ac_labels == 0).sum())),
        volume_error=volume_error, boundary_ac_seconds=boundary_seconds, comparison_seconds=perf_counter()-started), sampled


def render_ac_report(output, summary):
    originals = [json.loads(p.read_text(encoding='utf-8')) for p in sorted(output.glob('n*_*/result.json'))]
    originals.sort(key=lambda r: (r['candidate_count'], np.inf if r['budget'] is None else r['budget'],
                                 VARIANTS.index(r['variant'])))
    write_json(output/'summary.json', dict(records=originals,
        source_unchanged=summary['source_unchanged'], total=len(originals),
        certified=sum(r['status'] == 'certified' for r in originals),
        ac_validation=dict(divisions=summary['protocol']['divisions'],
                           summary_file='ac_validation/summary.json', findings=summary['findings'])))
    save_benchmark_report(output/'ac_validation', originals, LABELS, ac=summary, case_root=output)
    save_benchmark_report(output, originals, LABELS, ac=summary)


def attach_ac_to_replay(folder, original, check, divisions, ac_labels):
    """补充最终校验指标；原构域事件、时间和历史帧原样保留。"""
    original['region']['validation'] = dict(check['volume_error'], divisions=divisions,
        strict_vertices=check['strict_vertices'], maximum=check['maximum'], metric_scope='retained_inner_union')
    original['ac_audit'] = check
    write_json(folder/'result.json', original)
    with RunMonitor(record=False, output=folder, stream=StringIO()) as monitor:
        monitor.load_recording(folder/'replay.json')
        monitor.state['results'] = [original['region']]
        monitor.state['states'] = [ac_labels.tolist()]
        monitor.state['divisions'] = divisions
        monitor.state['ac_validation'] = dict(divisions=divisions, **check['volume_error'])
        monitor.save_snapshot()


def run_ac_benchmark(args):
    destination = args.output/'ac_validation'
    destination.mkdir(exist_ok=True)
    if args.report_only:
        summary = json.loads((destination/'summary.json').read_text(encoding='utf-8'))
        render_ac_report(args.output, summary)
        return
    before = fingerprints()
    original_protocol = json.loads((args.output/'protocol.json').read_text(encoding='utf-8'))
    if before != original_protocol['source_hashes']:
        raise ValueError('Production sources differ from saved benchmark; refuse mixed-version comparison')
    protocol = dict(divisions=args.divisions, ac_tolerance=AC_TOL,
        fixed_point_tolerance=FIXED_POINT_TOL, finite_budgets=[0, 1, 2],
        source_hashes=before,
        purpose='AC checks of saved domains only; sampling does not replace continuous-domain certification')
    write_json(destination/'protocol.json', protocol)
    records, references, bounds = [], [], {}
    with numerical_threads() as threads:
        if not threads['applied']:
            raise RuntimeError('Single-thread timing required')
        for count in args.counts:
            folder = destination/f'n{count}'
            folder.mkdir(exist_ok=True)
            originals = [json.loads(p.read_text(encoding='utf-8')) for p in sorted(args.output.glob(f'n{count}_*/result.json'))]
            if not originals:
                raise ValueError(f'No construction records for {count} candidates')
            for original in originals:
                assert original['source_hashes'] == before
            network = Case33(candidate_count=count)
            bounds[str(count)] = originals[0]['bounds']
            if (folder/'reference.json').exists():
                reference = json.loads((folder/'reference.json').read_text(encoding='utf-8'))
                assert reference['divisions'] == args.divisions and reference['source_hashes'] == before
                data = np.load(folder/'reference_grid.npz')
                points = data['points']
                grid = {b: data[f'b{"inf" if b is None else int(b)}'] for b in (0., 1., 2., None)}
                designs = affordable_designs(network, 2.)
            else:
                reference, grid, designs, points = build_ac_reference(network, originals, args.divisions, folder)
            references.append(reference)
            comparison = {f'ac_b{"inf" if b is None else int(b)}': labels for b, labels in grid.items()}
            for original in originals:
                record, labels = audit_ac_case(network, original, reference, grid[original['budget']], designs, points)
                records.append(record)
                name = case_name(count, np.inf if original['budget'] is None else original['budget'], original['variant'])
                write_json(folder/f'{name}.json', record)
                attach_ac_to_replay(args.output/name, original, record, args.divisions, grid[original['budget']])
                comparison[f'{original["variant"]}_b{"inf" if original["budget"] is None else int(original["budget"])}'] = labels
                print('AUDITED', name, 'vertices', record['strict_vertices']['feasible'], '/', record['strict_vertices']['total'],
                      'AC coverage', None if record['grid']['ac_retained_percent'] is None else
                      round(record['grid']['ac_retained_percent'], 5), flush=True)
            np.savez_compressed(folder/'comparison_grid.npz', **comparison)
            # Check a second independent AC implementation on base and upgraded networks.
            cross = []
            for choice in ([0]*count, [1]*count):
                cross.append(dict(choice=choice, **validate_power_flow(network.design(choice))))
            write_json(folder/'nodal_cross_checks.json', cross)
    records.sort(key=lambda r: (r['candidate_count'], np.inf if r['budget'] is None else r['budget'], VARIANTS.index(r['variant'])))
    certified = [r for r in records if r['construction_status'] == 'certified']
    vertex_total = sum(r['strict_vertices']['total'] for r in records)
    vertex_accepted = sum(r['strict_vertices']['feasible'] for r in records)
    definite_extra = sum(r['grid']['definite_extra'] for r in records)
    retained_missed = sum(r['grid']['ac_not_retained'] for r in certified)
    findings = [f'共检查 {len(records)} 组结果、{vertex_total} 个保留内域顶点；严格 AC 直接通过 {vertex_accepted} 个。'
                f'全部组网格中合计多算 {definite_extra} 点；已完成构域组的保留内域合计漏算 {retained_missed} 点（按案例累计）。',
                '比较结论须同时考虑构域是否完成、AC 误差和构域耗时；未完成案例的采样覆盖率只代表当前进度。']
    summary = dict(protocol=protocol, records=records, references=references, bounds=bounds,
        findings=findings, source_unchanged=before == fingerprints(),
        input_hashes={str(p.relative_to(args.output)): sha256(p.read_bytes()).hexdigest()
                      for p in args.output.glob('n*_*/result.json')})
    write_json(destination/'summary.json', summary)
    render_ac_report(args.output, summary)
    print('COMPLETE', destination/'report.html', flush=True)


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=DEFAULT)
    parser.add_argument('--counts', type=int, nargs='+', default=[8, 16])
    parser.add_argument('--divisions', type=int, default=32)
    parser.add_argument('--report-only', action='store_true')
    run_ac_benchmark(parser.parse_args())
