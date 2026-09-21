"""Isolated MP1/MP2/residual ablation; never edits production modules.

Run: python -m tests.benchmark_physical_search --counts 8 16
Each variant runs sequentially in a fresh process, with the same tolerances,
geometry, event recording and per-solve limits. Audits are outside solver timing.
"""
from argparse import ArgumentParser
from contextlib import nullcontext
from hashlib import sha256
from io import StringIO
from pathlib import Path
from time import perf_counter
from unittest.mock import patch
import html
import json
import subprocess
import sys
import traceback

import numpy as np

import main as production
import model as models
from model import PlanningModel, PlanningSP, RemainingRegionModel, PLANNING_TOL
from plot import RunMonitor, json_value
from region import GEOMETRY_TOL, sample_region
from tests.audit_continuous import independent_feasible
from tests.reference import _dispatch_equations

ROOT = Path(__file__).resolve().parents[1]
SOURCE_FILES = ('main.py', 'model.py', 'region.py', 'vertify.py', 'plot.py',
                'live_view.html', 'Network/__init__.py', 'Network/case33bw.py',
                'Network/data/case33bw.m')
VARIANTS = ('baseline', 'direct_mp', 'direct_all')
LABELS = dict(baseline='原流程', direct_mp='完整 MP1/MP2', direct_all='完整 MP1/MP2 + 物理剩余域')


def fingerprints():
    return {name: sha256((ROOT/name).read_bytes()).hexdigest() for name in SOURCE_FILES}


def write_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(json_value(value), ensure_ascii=False, indent=2,
                                    allow_nan=False), encoding='utf-8')
    temporary.replace(path)


class BenchmarkTimeout(RuntimeError):
    pass


class ExperimentRegion(production.ContinuousRegion):
    """Only test-side query/residual adapters; inherit production geometry unchanged."""

    def __init__(self, network, budget, bounds, variant, limit, monitor):
        self.variant, self.limit, self.monitor = variant, limit, monitor
        self.deadline, self.expired = np.inf, False
        self.phase_seconds = dict(MP2=0., MP1=0., residual=0.)
        self.query_results = []
        self.reused_certificates = 0
        self.physical_repairs = 0
        self.last_residual = None
        self.cached_witness = None
        self.last_checkpoint = 0.
        self.last_stage = 'preparing'
        super().__init__(network, 'socp', budget, bounds, self.query_adapter,
                         tau=production.REGION_TAU, observer=self.record)

    def record(self, event, **data):
        self.last_stage = event
        self.monitor(event, **data)
        if perf_counter()-self.last_checkpoint >= 2. or event in ('maximum', 'query_end', 'coverage'):
            write_json(self.monitor.output/'progress.json', dict(
                variant=self.variant, elapsed=perf_counter()-self.started, event=event,
                counts=self.counts, phases=self.phase_seconds, queries=self.query_results,
                max_total=None if self.maximum is None else float(sum(self.maximum['p'])),
                residual=self.last_residual))
            self.last_checkpoint = perf_counter()
        if not self.expired and perf_counter() >= self.deadline:
            raise BenchmarkTimeout('case wall-time budget reached')

    def query_adapter(self, equations, **kwargs):
        phase = 'MP1' if kwargs.get('power') is not None or kwargs.get('min_total') is not None else 'MP2'
        before = perf_counter()
        answer = None
        try:
            if self.variant == 'baseline':
                answer, cuts = production.joint_benders(equations, **kwargs)
            else:
                answer, cuts = self.direct_query(equations, **kwargs)
            return answer, cuts
        finally:
            seconds = perf_counter()-before
            self.phase_seconds[phase] += seconds
            self.query_results.append(dict(
                phase=phase, seconds=seconds,
                status='unfinished' if answer is None else answer['status'],
                objective=None if answer is None else answer['objective'],
                bound=None if answer is None else answer['bound'],
                feasible=False if answer is None else answer['feasible']))

    def direct_query(self, equations, *, power=None, budget=np.inf, direction=None,
                     cuts=(), start=None, radial_gap_kw=1e-3, observer=None,
                     min_total=None, fixed_x=None, incumbent=None):
        minimizing = power is not None or min_total is not None
        mode = 'MP1' if minimizing else 'MP2'
        production.emit(observer, 'query_model', message='测试：完整物理约束查询', reused_cuts=len(cuts))
        production.emit(observer, 'mp_start', mode=mode, iteration=1,
                        message=f'测试：完整物理模型 {mode}')
        problem = PlanningModel(equations, power=power, budget=budget, direction=direction,
                                min_total=min_total, fixed_x=fixed_x)
        with problem.model:
            for cut in cuts:
                problem.add_cut(cut)
            if start is not None:
                problem.x.Start = start
            if incumbent is not None:
                problem.use_incumbent(incumbent, budget=budget, power=power, min_total=min_total)
                problem.x.Start = incumbent['x']
                problem.state.Start = incumbent['state']
            answer = problem.solve(time_limit=min(20., max(0., self.deadline-perf_counter())))
        if answer is None:
            production.emit(observer, 'query_end', status='infeasible', message='完整物理模型已证不可行')
            return None, []
        if answer['x'] is not None and not minimizing and answer['p'].sum() > 0.:
            # Match the baseline's 0.001 kW retreat before MP1; certify the new
            # state against original equations rather than assume monotonicity.
            point = answer['p']*max(0., 1.-radial_gap_kw/answer['p'].sum())
            state = equations.restore(answer['x'], point, answer['state'])
            if equations.margin(answer['x'], point, state) >= -PLANNING_TOL:
                answer.update(p=point, state=state, feasible=True, objective=float(point.sum()))
        if answer['x'] is not None and not answer['feasible']:
            self.physical_repairs += 1
            production.emit(observer, 'sp_start', point=answer['p'],
                            choice=equations.choice(answer['x']), message='数值证书未通过，保留必要 SP 复核')
            checked = PlanningSP(equations).solve(answer['x'], answer['p'])
            production.emit(observer, 'sp_end', feasible=checked['feasible'], eta=checked['eta'])
            if checked['feasible']:
                answer.update(state=checked['state'], feasible=True)
        if answer['feasible']:
            value = answer['objective']
            gap = value-answer['bound'] if minimizing else answer['bound']-value
            answer['status'] = 'optimal' if gap <= (1e-7 if minimizing else radial_gap_kw+1e-5) else 'feasible'
            self.reused_certificates += 1
        else:
            answer['status'] = 'unknown'
        production.emit(observer, 'query_end', status=answer['status'], objective=answer['objective'],
                        bound=answer['bound'], message='测试：复用完整物理模型的原始证书')
        return answer, []

    def residual(self):
        self.counts['residual'] += 1
        self.emit('mp_start', mode='residual', iteration=self.counts['residual'],
                  message='测试：全局剩余域搜索')
        started = perf_counter()
        state = self.region

        def full_model(*args, **kwargs):
            kwargs['cuts_only'] = False
            return PlanningModel(*args, **kwargs)

        # A constructor-only substitution, confined to this worker process.
        # All union exclusion constraints remain those of the original class.
        context = patch.object(models, 'PlanningModel', full_model) if self.variant == 'direct_all' else nullcontext()
        with context:
            problem = RemainingRegionModel(self.equations, self.budget, self.bounds, state.total_bound,
                                           state.cuts, state.inner_halfspaces(), self.tau)
        with problem.model:
            answer = problem.solve(GEOMETRY_TOL, time_limit=min(120., max(0., self.deadline-perf_counter())))
            if self.variant == 'direct_all' and answer['x'] is not None:
                certificate = self.equations.restore(answer['x'], answer['p'], problem.problem.state.X)
                margin = self.equations.margin(answer['x'], answer['p'], certificate)
                if margin >= -PLANNING_TOL:
                    self.cached_witness = (answer['x'].copy(), answer['p'].copy())
                    self.reused_certificates += 1
                else:
                    # The incumbent is still a search candidate, but is not
                    # inserted into the inner union without a certificate.
                    self.physical_repairs += 1
                answer['physical_margin'] = margin
            answer['solver_status'] = int(problem.model.Status)
            answer['binary_variables'] = problem.model.NumBinVars
            answer['linear_constraints'] = problem.model.NumConstrs
            answer['quadratic_constraints'] = problem.model.NumQConstrs
        self.times['mp_seconds'] += answer.pop('solve_seconds')
        self.phase_seconds['residual'] += perf_counter()-started
        self.last_residual = dict(answer)
        return answer

    def refine(self, x, seed=None, *, witness=None, max_checks=96):
        if witness is not None and self.cached_witness is not None:
            known_x, known_p = self.cached_witness
            if np.array_equal(known_x, x) and np.array_equal(known_p, witness):
                seed = known_p
                self.cached_witness = None
        return super().refine(x, seed=seed, witness=witness, max_checks=max_checks)


def audit_region(network, budget, bounds, result):
    """Independent fixed-scheme equations for every retained inner vertex."""
    started = perf_counter()
    checked, failures = 0, []
    for certificate in result['certificates']:
        points = np.asarray(certificate['inner']).reshape(-1, 3)*bounds
        if not len(points):
            continue
        design = network.design(certificate['choice'])
        valid = independent_feasible(_dispatch_equations(design), 'socp', points)
        checked += len(points)
        if design.cost > budget+1e-8 or not valid.all():
            failures.append(dict(choice=certificate['choice'], cost=design.cost,
                                 failed_points=points[~valid].tolist()))
    equation = models.PlanningEquations(network, 'socp')
    probe = PlanningModel(equation, budget=budget)
    with probe.model:
        reference = probe.solve(time_limit=20.)
    maximum_error = None if reference is None or result['max_total'] is None else abs(reference['objective']-result['max_total'])
    maximum_passed = bool(reference and reference['feasible'] and
                          (maximum_error is None or maximum_error <= .003))
    # A small deterministic set of global radial queries tests the returned
    # envelopes; it is not used as a substitute for the global coverage bound.
    boundary_checks, boundary_failures, radial_unknown = 0, [], 0
    for direction in np.vstack([np.eye(3), np.ones(3), [[1, 6, 1], [3, 1, 4], [1, 1, 5]]]):
        query = PlanningModel(equation, budget=budget, direction=direction)
        with query.model:
            answer = query.solve(time_limit=20.)
        if answer is None or not answer['feasible']:
            radial_unknown += 1
            continue
        label = int(sample_region(result, [answer['p']], bounds)[0])
        boundary_checks += 1
        if label == -1:
            boundary_failures.append(answer['p'].tolist())
    coverage_passed = result['status'] != 'certified' or (
        result['coverage_bound'] is None or result['coverage_bound'] <= GEOMETRY_TOL)
    return dict(passed=not failures and not boundary_failures and maximum_passed and coverage_passed,
                inner_vertices_checked=checked, inner_failures=failures,
                reference_maximum=reference, maximum_error_kw=maximum_error,
                boundary_checks=boundary_checks, boundary_failures=boundary_failures,
                radial_queries_unknown=radial_unknown, coverage_passed=coverage_passed,
                seconds=perf_counter()-started,
                scope='Independent SOCP inner-vertex checks and seven global rays; AC not tested.')


def case_name(count, budget, variant):
    return f'n{count}_b{"inf" if np.isinf(budget) else f"{budget:g}"}_{variant}'


def worker(args):
    count, budget, variant = args.counts[0], args.budgets[0], args.variants[0]
    output = args.output/case_name(count, budget, variant)
    output.mkdir(parents=True, exist_ok=True)
    before_hashes = fingerprints()
    with production.numerical_threads() as threads:
        if not threads['applied']:
            raise RuntimeError('Single-thread timing unavailable; benchmark refused')
        network = production.Case33(candidate_count=count)
        prepare = perf_counter()
        bounds = production.evaluation_bounds(network)
        preparation_seconds = perf_counter()-prepare
        with RunMonitor(record=True, output=output, open_browser=False, stream=StringIO()) as monitor:
            monitor('start', network=network.name, cost_unit=network.cost_unit,
                    candidate_count=count, message=LABELS[variant])
            monitor('method_start', method='socp', method_number=1, method_count=1)
            monitor('phase_start', method='socp', phase='socp', phase_number=1, phase_count=1,
                    budget_index=0, budget=budget, budgets=[budget], bounds=bounds,
                    load_nodes=network.load_nodes, representation='continuous', message=LABELS[variant])
            solver = ExperimentRegion(network, budget, bounds, variant, args.limit, monitor)
            solver.deadline = solver.started+args.limit
            termination = 'returned'
            try:
                region = solver.solve()
            except BenchmarkTimeout:
                solver.expired = True
                termination = 'case_time_limit'
                maximum = None if solver.maximum is None else float(sum(solver.maximum['p']))
                bound = None if solver.last_residual is None else solver.last_residual['bound']
                region = solver.finish('timeout', bound, maximum)
            solve_seconds = perf_counter()-solver.started
            solver.expired = True
            monitor('method_end', seconds=solve_seconds)
            monitor('completed', message=f'{LABELS[variant]}：{region["status"]}', region=region)
            monitor.save_snapshot()
            record = dict(candidate_count=count, budget=budget, variant=variant,
                          status=region['status'], termination=termination,
                          solve_seconds=solve_seconds, preparation_seconds=preparation_seconds,
                          case_time_limit=args.limit, thread_control=threads, bounds=bounds,
                          phases=solver.phase_seconds, queries=solver.query_results,
                          certificates_reused=solver.reused_certificates,
                          numerical_repairs=solver.physical_repairs,
                          last_residual=solver.last_residual, region=region,
                          source_hashes=before_hashes, source_unchanged=before_hashes == fingerprints(),
                          test_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
            write_json(output/'result.json', record)
        record['audit'] = audit_region(network, budget, bounds, region)
        record['source_unchanged'] &= before_hashes == fingerprints()
        write_json(output/'result.json', record)
    print(json.dumps(json_value(dict(case=output.name, status=record['status'],
              seconds=solve_seconds, phases=record['phases'], sp=region['counts']['sp'],
              max_total=region['max_total'], audit=record['audit']['passed'])), ensure_ascii=False), flush=True)


def report(output):
    records = [json.loads(p.read_text(encoding='utf-8')) for p in sorted(output.glob('n*_b*/result.json'))]
    groups, comparisons = {}, []
    for row in records:
        groups.setdefault((row['candidate_count'], row['budget']), {})[row['variant']] = row
    for key, group in groups.items():
        old = group.get('baseline')
        for name in ('direct_mp', 'direct_all'):
            new = group.get(name)
            if not old or not new:
                continue
            a, b = old['region'], new['region']
            conflict_count = None
            if a['inner'] and b['inner']:
                bounds = np.asarray(old['bounds'])
                points = (np.indices((20,)*3).reshape(3, -1).T+.5)*bounds/20
                conflict_count = int(np.count_nonzero(sample_region(a, points, bounds)*sample_region(b, points, bounds) == -1))
            comparisons.append(dict(candidate_count=key[0], budget=key[1], variant=name,
                both_certified=old['status'] == new['status'] == 'certified',
                speed_ratio=old['solve_seconds']/new['solve_seconds'],
                ratio_is_lower_bound=old['status'] == 'timeout' and new['status'] == 'certified',
                definite_grid_conflicts=conflict_count,
                maximum_difference_kw=None if a['max_total'] is None or b['max_total'] is None else abs(a['max_total']-b['max_total']),
                volume_intervals_overlap=max(a['inner_volume'], b['inner_volume']) <= min(a['outer_volume'], b['outer_volume'])*(1+1e-8)))
    compact = [dict(candidate_count=r['candidate_count'], budget=r['budget'], variant=r['variant'],
                    status=r['status'], seconds=r['solve_seconds'], phases=r['phases'],
                    counts=r['region']['counts'], maximum=r['region']['max_total'],
                    maximum_bound=r['region']['max_total_bound'], coverage_bound=r['region']['coverage_bound'],
                    volume_gap=r['region']['volume_gap'], audit=r.get('audit'),
                    source_unchanged=r['source_unchanged']) for r in records]
    summary = dict(protocol='Fresh sequential processes; SOCP; tau=.002; original tolerance; full replay; audits excluded',
                   records=compact, comparisons=comparisons,
                   source_unchanged=all(r['source_unchanged'] for r in records))
    write_json(output/'summary.json', summary)
    rows = []
    for row in records:
        region, phases = row['region'], row['phases']
        name = case_name(row['candidate_count'], np.inf if row['budget'] is None else row['budget'], row['variant'])
        audit = row.get('audit', {})
        maximum_text = '—' if region['max_total'] is None else format(region['max_total'], '.3f')
        rows.append(f'<tr><td>{row["candidate_count"]}</td><td>{"无限" if row["budget"] is None else row["budget"]}</td>'
            f'<td>{html.escape(LABELS[row["variant"]])}</td><td>{row["status"]}</td>'
            f'<td>{row["solve_seconds"]:.3f}</td><td>{phases["MP2"]:.3f}</td><td>{phases["MP1"]:.3f}</td>'
            f'<td>{phases["residual"]:.3f}</td><td>{region["counts"]["sp"]}</td>'
            f'<td>{maximum_text}</td>'
            f'<td>{"通过" if audit.get("passed") else "待查"}</td><td><a href="{name}/live_view.html">回放</a></td></tr>')
    document = '''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>8 / 16 候选线路 · 独立提速测试</title><style>
body{font:15px/1.6 system-ui,sans-serif;margin:32px;color:#202630;background:#fafbfc}h1{font-size:24px}
.table{overflow:auto}table{border-collapse:collapse;background:#fff;width:100%}th,td{padding:9px 12px;text-align:left;border-bottom:1px solid #dde2e8;white-space:nowrap}th{background:#edf1f6}a{color:#145ca3}
</style><h1>8 / 16 候选线路 · SOCP 独立提速测试</h1>
<p>正式源代码保持不变。三组使用相同物理方程、τ = 0.002、原始可行性容差、单线程和完整回放。总耗时包含构域和记录，不包含公共评价箱准备及独立审核。</p>
<p>certified：获得全局覆盖证书；unknown：求解器或几何过程未完成；timeout：达到整组时间预算。超时不能当作完整域结果。各组单次 MP 上限 20 秒，剩余搜索上限 120 秒。</p>
<div class="table"><table><thead><tr><th>候选线路</th><th>预算</th><th>方案</th><th>完成状态</th><th>总秒数</th><th>MP2 秒</th><th>MP1 秒</th><th>剩余搜索秒</th><th>SP 次数</th><th>最大负荷 kW</th><th>审核</th><th>过程</th></tr></thead><tbody>'''
    document += ''.join(rows)+'</tbody></table></div><p>审核包含独立消元方程对所有保留内域顶点的复核，以及七条全局射线边界检查。全域完整性仍由全局覆盖界判定。本轮不测试 AC。</p><p><a href="summary.json">原始汇总数据</a></p></html>'
    (output/'report.html').write_text(document, encoding='utf-8')
    return summary


def supervisor(args):
    args.output.mkdir(parents=True, exist_ok=True)
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
    parser.add_argument('--variants', nargs='+', choices=VARIANTS, default=list(VARIANTS))
    parser.add_argument('--limit', type=float, default=300.)
    parser.add_argument('--output', type=Path, default=ROOT/'results/case33bw/physical_search_benchmark')
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
