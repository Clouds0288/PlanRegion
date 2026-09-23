"""Production mainline benchmark with a frozen joint-cut baseline.

Run: python -m tests.benchmark_physical_search --counts 8 16
Each variant runs sequentially in a fresh process, with the same tolerances,
geometry and per-solve limits. Audits are outside solver timing.
"""
from argparse import ArgumentParser
from hashlib import sha256
from pathlib import Path
from time import perf_counter, sleep
import json
import subprocess
import sys
import traceback

import numpy as np
from gurobipy import GRB

import main as production
from threadpoolctl import threadpool_limits
from model import evaluation_bounds
from region import ContinuousRegion
from model import PlanningEquations, PlanningModel, PlanningSP, PLANNING_TOL
from plot import json_value, save_benchmark_report, region_view
from region import GEOMETRY_TOL
from plot import sample_region
from tests.audit_continuous import independent_feasible
from tests.reference import _dispatch_equations, dispatch_state

ROOT = Path(__file__).resolve().parents[1]
SOURCE_FILES = ('main.py', 'model.py', 'region.py', 'vertify.py', 'plot.py',
                'live_view.html', 'Network/__init__.py', 'Network/case33bw.py',
                'Network/data/case33bw.m')
VARIANTS = ('baseline', 'direct_mp', 'direct_all', 'auto')
LABELS = dict(baseline='旧联合割基准', direct_mp='完整 MP1/MP2＋轻量剩余域',
              direct_all='完整 MP1/MP2＋物理剩余域', auto='主线 auto')


def fingerprints():
    return {name: sha256((ROOT/name).read_bytes()).hexdigest() for name in SOURCE_FILES}


def write_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(json_value(value), ensure_ascii=False, indent=2,
                                    allow_nan=False), encoding='utf-8')
    # Windows readers can briefly deny replacement of a progress snapshot.
    for attempt in range(10):
        try:
            temporary.replace(path)
            break
        except PermissionError:
            if attempt == 9:
                raise
            sleep(.025)


def joint_benders(equations, *, power=None, budget=np.inf,
                  cuts=(), start=None, radial_gap_kw=1e-3,
                  min_total=None, fixed_plan=None, incumbent=None, deadline=np.inf, threads=1, oracle=None):
    """供独立对照使用的 MP1/MP2 与 SP 联合割迭代。"""
    problem = PlanningModel(equations, power=power, budget=budget,
                            cuts_only=True, min_total=min_total, fixed_plan=fixed_plan, threads=threads)
    oracle, generated = oracle or PlanningSP(equations, threads=threads), []
    minimizing = power is not None or min_total is not None
    bound = -np.inf if minimizing else np.inf
    with problem.model:
        for cut in cuts:
            problem.add_cut(cut)
        if start is not None:
            problem.x_vector.Start = start
        if incumbent is not None:
            incumbent_cost = problem.use_incumbent(incumbent)
        while True:
            if perf_counter() >= deadline:
                return dict(feasible=False, x=None, p=None, state=None, objective=None,
                            bound=bound, status='unknown', termination='case_time_limit'), generated
            answer = problem.solve(time_limit=min(20., max(0., deadline-perf_counter())))
            if answer is None:
                return None, generated
            bound = max(bound, answer['bound']) if minimizing else min(bound, answer['bound'])
            answer['bound'] = bound
            if incumbent is not None and bound >= incumbent_cost-1e-7:
                answer = dict(incumbent, objective=incumbent_cost, bound=bound, status='optimal')
                return answer, generated
            if answer['x'] is None:
                return answer, generated
            point = answer['p']
            if not minimizing and point.sum() > 0.:
                point = point * max(0., 1. - radial_gap_kw / point.sum())
            checked = oracle.solve(answer['x'], point, time_limit=min(
                5. if equations.method == 'linear' else 10., max(0., deadline-perf_counter())))
            if checked['feasible']:
                value = answer['objective'] if minimizing else point.sum()
                gap = value-bound if minimizing else bound-value
                tolerance = 1e-7 if minimizing else radial_gap_kw+1e-5
                answer.update(p=point, state=checked['state'], feasible=True, objective=value,
                              status='optimal' if gap <= tolerance else 'feasible')
                return answer, generated
            if checked['cut'] is None:
                answer.update(status='unknown', feasible=False, objective=None)
                return answer, generated
            generated.append(checked['cut'])
            problem.add_cut(checked['cut'])



def independent_boundary_certificate(network, equations, point):
    """Alternate independent primal witness, evaluated at the unchanged point.

    Near an axis, phase I and minimum-current objectives may stall. A radial
    objective can find a better state; its objective alone NEVER certifies the
    tested point. Recheck all original linear/conic residuals at that point.
    """
    if np.sum(point) <= 0:
        return dict(passed=False, point=np.asarray(point).tolist())
    try:
        current = dispatch_state(network, point)
    except RuntimeError as error:
        return dict(passed=False, point=np.asarray(point).tolist(), error=str(error))
    residual = equations.c+equations.F@point+equations.G@current
    margins = [float(residual[:equations.linear_count].min())]
    offset = equations.linear_count
    for size in equations.sizes:
        margins.append(float(residual[offset]-np.linalg.norm(residual[offset+1:offset+size])))
        offset += size
    return dict(passed=bool(np.isfinite(residual).all() and min(margins) >= -PLANNING_TOL),
                point=np.asarray(point).tolist(), min_margin=min(margins),
                tolerance=PLANNING_TOL, current=current.tolist(), point_changed=False)


def audit_region(network, budget, bounds, result):
    """Independent fixed-scheme equations for every retained inner vertex."""
    started = perf_counter()
    checked, failures, numerical_rechecks = 0, [], []
    for certificate in region_view(result, bounds)['geometry']:
        points = np.asarray(certificate['inner']['vertices']).reshape(-1, 3)
        if not len(points):
            continue
        design = network.design(certificate['choice'])
        equations = _dispatch_equations(design)
        valid = independent_feasible(equations, 'socp', points)
        for index in np.flatnonzero(~valid):
            recheck = independent_boundary_certificate(design, equations, points[index])
            numerical_rechecks.append(dict(choice=design.x.tolist(), **recheck))
            valid[index] = recheck['passed']
        checked += len(points)
        if design.cost > budget+1e-8 or not valid.all():
            failures.append(dict(choice=certificate['choice'], cost=design.cost,
                                 failed_points=points[~valid].tolist()))
    equation = PlanningEquations(network, 'socp')
    probe = PlanningModel(equation, budget=budget, threads=1)
    with probe.model:
        reference = probe.solve(time_limit=20.)
    maximum_error = None if reference is None or result['max_total'] is None else abs(reference['objective']-result['max_total'])
    maximum_passed = bool(reference and reference['feasible'] and
                          (maximum_error is None or maximum_error <= .003))
    # A small deterministic set of global support queries tests the returned
    # envelopes; it is not used as a substitute for the global coverage bound.
    boundary_checks, boundary_failures, support_unknown = 0, [], 0
    for normal in np.vstack([np.eye(3), np.ones(3), [[1, 6, 1], [3, 1, 4], [1, 1, 5]]]):
        query = PlanningModel(equation, budget=budget, threads=1)
        with query.model:
            query.model.setObjective(normal@query.power/network.base, GRB.MAXIMIZE)
            answer = query.solve(time_limit=20.)
        if answer is None or not answer['feasible']:
            support_unknown += 1
            continue
        label = int(sample_region(result, [answer['p']], bounds)[0])
        boundary_checks += 1
        if label == -1:
            boundary_failures.append(answer['p'].tolist())
    coverage_passed = result['status'] != 'certified' or (
        result['coverage_bound'] is None or result['coverage_bound'] <= GEOMETRY_TOL)
    return dict(passed=not failures and not boundary_failures and maximum_passed and coverage_passed,
                inner_vertices_checked=checked, inner_failures=failures,
                numerical_rechecks=numerical_rechecks,
                reference_maximum=reference, maximum_error_kw=maximum_error,
                boundary_checks=boundary_checks, boundary_failures=boundary_failures,
                support_queries_unknown=support_unknown, coverage_passed=coverage_passed,
                seconds=perf_counter()-started,
                scope='Independent SOCP inner-vertex checks and seven global support points; AC not tested.')


def case_name(count, budget, variant):
    return f'n{count}_b{"inf" if np.isinf(budget) else f"{budget:g}"}_{variant}'


def baseline_query(equations, *, cuts, **kwargs):
    """旧联合割接口仅在基准测试中适配。"""
    answer, generated = joint_benders(equations, cuts=cuts, **kwargs)
    cuts.extend(generated)
    return answer


def worker(args):
    count, budget, variant = args.counts[0], args.budgets[0], args.variants[0]
    output = args.output/case_name(count, budget, variant)
    output.mkdir(parents=True, exist_ok=True)
    before_hashes = fingerprints()
    with threadpool_limits(limits=1):
        network = production.Case33(candidate_count=count)
        prepare = perf_counter()
        bounds = evaluation_bounds(network, threads=1)
        preparation_seconds = perf_counter()-prepare
        mode = dict(baseline='light', direct_mp='light', direct_all='physical', auto='auto')[variant]
        solver = ContinuousRegion(network, 'socp', budget, bounds,
            query=baseline_query if variant == 'baseline' else None,
            residual_mode=mode, time_limit=args.limit, threads=1)
        region = production.solve_region(solver)
        solve_seconds = region['timing']['total_seconds']
        region['budget_index'] = 0
        record = dict(candidate_count=count, budget=budget, variant=variant,
                      status=region['status'], solve_seconds=solve_seconds,
                      termination='case_time_limit' if region['status']=='time_limit' else 'returned',
                      preparation_seconds=preparation_seconds, case_time_limit=args.limit,
                      threads=1, bounds=bounds, region=region,
                      source_hashes=before_hashes, source_unchanged=before_hashes == fingerprints(),
                      test_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
        write_json(output/'result.json', record)
        record['audit'] = audit_region(network, budget, bounds, region)
        record['source_unchanged'] &= before_hashes == fingerprints()
        write_json(output/'result.json', record)
    print(json.dumps(json_value(dict(case=output.name, status=record['status'],
              seconds=solve_seconds, sp=region['counts']['sp'],
              max_total=region['max_total'], audit=record['audit']['passed'])), ensure_ascii=False), flush=True)


def report(output):
    records = [json.loads(p.read_text(encoding='utf-8')) for p in sorted(output.glob('n*_b*/result.json'))]
    records.sort(key=lambda r: (r['candidate_count'], np.inf if r['budget'] is None else r['budget'],
                               VARIANTS.index(r['variant'])))
    summary = dict(records=records, source_unchanged=all(r['source_unchanged'] for r in records),
                   certified=sum(r['status'] == 'certified' for r in records), total=len(records))
    write_json(output/'summary.json', summary)
    save_benchmark_report(output, records, LABELS)
    return summary


def supervisor(args):
    args.output.mkdir(parents=True, exist_ok=True)
    for saved in args.output.glob('n*_b*/result.json'):
        record = json.loads(saved.read_text(encoding='utf-8'))
        if record['source_hashes'] != fingerprints():
            raise ValueError('Saved results use different production sources; choose a new --output directory')
    write_json(args.output/'protocol.json', dict(counts=args.counts, budgets=args.budgets,
        variants=args.variants, time_limit=args.limit, source_hashes=fingerprints(),
        tau=production.REGION_TAU, physical_tolerance=PLANNING_TOL, geometry_tolerance=GEOMETRY_TOL,
        repeats=1, order='Variant order rotates by case; cases run sequentially.'))
    index = 0
    for count in args.counts:
        for budget in args.budgets:
            order = args.variants[index % len(args.variants):]+args.variants[:index % len(args.variants)]
            index += 1
            for variant in order:
                folder = args.output/case_name(count, budget, variant)
                if (folder/'result.json').exists():
                    existing = json.loads((folder/'result.json').read_text(encoding='utf-8'))
                    if 'audit' in existing:
                        continue
                folder.mkdir(parents=True, exist_ok=True)
                command = [sys.executable, '-u', '-X', 'utf8', '-m', 'tests.benchmark_physical_search', '--worker',
                           '--counts', str(count), '--budgets', str(budget), '--variants', variant,
                           '--limit', str(args.limit), '--output', str(args.output)]
                print('START', folder.name, flush=True)
                with (folder/'worker.log').open('w', encoding='utf-8') as log:
                    process = subprocess.Popen(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
                    try:
                        code = process.wait(timeout=args.limit+240.)
                    except subprocess.TimeoutExpired:
                        process.kill()
                        code = process.wait()
                        write_json(folder/'error.json', dict(error='worker watchdog', exit_code=code))
                if code:
                    print('FAILED', folder.name, 'exit', code, flush=True)
                elif (folder/'result.json').exists():
                    row = json.loads((folder/'result.json').read_text(encoding='utf-8'))
                    print('DONE', folder.name, row['status'], round(row['solve_seconds'], 3),
                          'SP', row['region']['counts']['sp'], 'audit', row.get('audit', {}).get('passed'), flush=True)
                report(args.output)
    report(args.output)


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--counts', nargs='+', type=int, default=[8, 16])
    parser.add_argument('--budgets', nargs='+', type=float, default=[0., 1., 2., np.inf])
    parser.add_argument('--variants', nargs='+', choices=VARIANTS, default=['auto'])
    parser.add_argument('--limit', type=float, default=300.)
    parser.add_argument('--output', type=Path, default=ROOT/'results/case33bw/latest')
    parser.add_argument('--worker', action='store_true')
    parser.add_argument('--report-only', action='store_true')
    arguments = parser.parse_args()
    arguments.output = arguments.output.resolve()
    if arguments.report_only:
        report(arguments.output)
    elif arguments.worker:
        try:
            worker(arguments)
        except Exception:
            traceback.print_exc()
            raise
    else:
        supervisor(arguments)
