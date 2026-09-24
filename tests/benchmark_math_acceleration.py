"""Isolated Case33 acceleration experiments; production modules are never edited.

Examples (run from the repository root):
  python tests/benchmark_math_acceleration.py --prepare --count 8
  python tests/benchmark_math_acceleration.py --variant baseline --budget 2
  python tests/benchmark_math_acceleration.py --variant combined --budget 2
"""
import argparse
import hashlib
import inspect
import json
import platform
import sys
import textwrap
from itertools import product
from pathlib import Path
from time import perf_counter

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import gurobipy as gp
import numpy as np
from gurobipy import GRB
from threadpoolctl import threadpool_limits

import main as application
import model
from Network.case33bw import Case33
from region import GEOMETRY_TOL, RegionState, contains, halfspaces, initial_polytope
from tests.planning_checks import margin
from vertify import ACPowerFlow

OUTPUT = ROOT/'results'/'case33_math_acceleration'/'20260923'
SOURCE_FILES = ('main.py', 'model.py', 'region.py', 'vertify.py',
                'Network/__init__.py', 'Network/case33bw.py',
                'Network/data/case33bw.m', 'tests/planning_checks.py')
BOX_CORNERS = np.array(list(product((0., 1.), repeat=3)))


def fingerprints():
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            for name in SOURCE_FILES}


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False,
        default=lambda v: v.tolist() if isinstance(v, np.ndarray) else v.item()), encoding='utf-8')


def tight_remaining(self, equations, budget, bounds, total_bound, cuts,
                    inner_halfspaces, tau, *, mode='light', threads=20):
    """The original disjunction with an analytic upper bound and face-specific M."""
    self.mode = mode
    self.problem = problem = model.PlanningModel(equations, budget=budget,
        cuts_only=mode == 'light', threads=threads)
    m = self.model = problem.model
    problem.power.UB = bounds
    m.addConstr(problem.power.sum() <= total_bound)
    self.distance_scale = scale = 1000.
    for cut in cuts:
        cut = np.asarray(cut)
        magnitude = max(np.max(np.abs(np.r_[cut[0], cut[1:4]*bounds, cut[4:]])), 1e-20)
        problem.add_cut(cut*scale/magnitude)
    box = initial_polytope(bounds, total_bound)
    extrema = [box@((1-tau)*eq[:, :3]).T+eq[:, 3] for eq in inner_halfspaces]
    # Outward rounding is more conservative than the exact extrema.
    padding = 1e-10
    upper = min((float(values.max())+padding for values in extrema), default=4.)
    delta = m.addVar(lb=-4.*scale, ub=upper*scale, name='uncovered_distance')
    for k, (eq, values) in enumerate(zip(inner_halfspaces, extrema)):
        select = m.addVars(len(eq), vtype=GRB.BINARY, name=f'outside_{k}')
        m.addConstr(select.sum() == 1)
        lower = values.min(axis=0)-padding
        for f, face in enumerate(eq):
            a = face[:3]*(1-tau)
            big_m = max(0., upper-lower[f])
            expression = gp.quicksum(float(a[j]/bounds[j])*problem.power[j].item() for j in range(3))
            m.addConstr(delta <= scale*(expression+float(face[3])+float(big_m)*(1-select[f])))
    m.setObjective(delta, GRB.MAXIMIZE)
    self.tight_upper = upper


class DownwardRegion(RegionState):
    def add_point(self, x, point):
        points = np.asarray(point).reshape(-1, 3)
        corners = (points[:, None, :]*BOX_CORNERS).reshape(-1, 3)
        return super().add_point(x, corners)


def single_linear_solve():
    """Copy only this experimental method; keep the original SOCP route verbatim."""
    source = textwrap.dedent(inspect.getsource(model.PlanningSP.solve))
    old = '        model.optimize()\n        if model.Status != GRB.OPTIMAL or model.ObjVal <= 0.:'
    new = ('        if equations.method != \'linear\' or model.Status != GRB.OPTIMAL:\n'
           '            model.optimize()\n'
           '        if model.Status != GRB.OPTIMAL or model.ObjVal <= 0.:')
    assert source.count(old) == 1
    namespace = dict(vars(model))
    exec(compile(source.replace(old, new), '<experimental-linear-SP>', 'exec'), namespace)
    return namespace['solve']


def audit(network, result, certificates, cuts):
    """Outside measured time: original equations, all inner vertices, independent AC."""
    started = perf_counter()
    original_margins = [margin(model.PlanningEquations(network, method), x, power, state)
                        for method, x, power, state in certificates]
    checked, ac_bad, ac_unknown, budget_bad = 0, [], [], []
    valid_samples = []
    for row in result['inner']:
        x = network.encode_plan(row['choice'])
        tree = network.tree(x)
        if tree.cost > result['budget']+1e-8:
            budget_bad.append(tree.cost)
        points = np.asarray(row['vertices']).reshape(-1, 3)
        reference = ACPowerFlow(tree, threads=1)
        labels, currents = reference.classify(points, return_currents=True)
        checked += len(points)
        for index in np.flatnonzero(labels != 1):
            P, Q, v, u = reference.state(points[index], currents[index])
            detail = dict(point=points[index], choice=row['choice'],
                          limit_violation=float(reference.violation(P, Q, v).max()),
                          equality_residual=float(np.max(np.abs(P*P+Q*Q-u*currents[index]))))
            (ac_bad if labels[index] == -1 else ac_unknown).append(detail)
        valid_samples.extend((x, point) for point in points[labels == 1])
    minimum_cut = None
    if valid_samples and cuts:
        values = np.array([np.r_[1., point, x] for x, point in valid_samples])
        minimum_cut = float(np.min(values@np.asarray(cuts).T))
    return dict(original_state_count=len(certificates),
                minimum_original_margin=min(original_margins, default=None),
                original_state_failures=sum(v < -model.PLANNING_TOL for v in original_margins),
                inner_vertices_checked=checked, ac_failures=ac_bad, ac_unknown=ac_unknown,
                budget_failures=budget_bad, minimum_cut_on_ac_feasible_points=minimum_cut,
                coverage_bound_accepted=(result['status'] != 'certified' or
                    result['coverage_bound'] is None or result['coverage_bound'] <= GEOMETRY_TOL),
                seconds=perf_counter()-started)


def prepare(args):
    output = Path(args.output)
    network = Case33(upgrade_count=args.count)
    before = fingerprints()
    started = perf_counter()
    bounds = model.evaluation_bounds(network, threads=args.threads)
    data = dict(candidate_count=args.count, bounds=bounds, threads=args.threads,
                seconds=perf_counter()-started, source_hashes=before,
                network_fingerprint=network.fingerprint,
                gurobi_version=gp.gurobi.version(), python=platform.python_version(),
                machine=platform.platform(), source_unchanged=before == fingerprints())
    write_json(output/f'bounds_n{args.count}.json', data)
    for name in SOURCE_FILES:
        path = output/'source'/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes((ROOT/name).read_bytes())
    print(json.dumps(data, default=lambda value: value.tolist()), flush=True)


def run(args):
    output = Path(args.output)
    before = fingerprints()
    prepared = json.loads((output/f'bounds_n{args.count}.json').read_text(encoding='utf-8'))
    assert before == prepared['source_hashes'], 'Production sources changed after bounds preparation'
    network = Case33(upgrade_count=args.count)
    assert network.fingerprint == prepared['network_fingerprint']
    assert (np.all(network.r > 0.) and np.all(network.reactance >= 0.)
            and np.all(network.original_p >= 0.) and np.all(network.original_q >= 0.)
            and np.all(network.q_ratio >= 0.) and np.all(network.vmax >= 1.))
    bounds = np.asarray(prepared['bounds'])
    label = f'n{args.count}_b{args.budget:g}_{args.method}_{args.variant}_t{args.threads}_r{args.repeat}'
    directory = output/label
    directory.mkdir(parents=True, exist_ok=True)
    events, mp_calls, sp_calls, remaining_calls, cuts, certificates = [], [], [], [], [], []
    started = perf_counter()
    last_checkpoint = started
    sp_method = single_linear_solve() if args.variant in ('linear_single', 'combined') else model.PlanningSP.solve

    class TimedModel(model.PlanningModel):
        def solve(self, *values, **options):
            tick = perf_counter()
            answer = super().solve(*values, **options)
            mp_calls.append(dict(seconds=perf_counter()-tick, runtime=self.model.Runtime,
                status=self.model.Status, nodes=self.model.NodeCount,
                answer_status=None if answer is None else answer['status'],
                objective=None if answer is None else answer['objective'],
                bound=None if answer is None else answer['bound']))
            if answer is not None and answer['feasible']:
                certificates.append((self.equations.method, answer['x'], answer['p'], answer['state']))
            return answer

    class TimedSP(model.PlanningSP):
        def solve(self, x, power, *values, **options):
            tick = perf_counter()
            answer = sp_method(self, x, power, *values, **options)
            sp_calls.append(dict(seconds=perf_counter()-tick, feasible=answer['feasible'], cut=answer['cut'] is not None))
            if answer['feasible']:
                certificates.append((self.equations.method, x.copy(), np.asarray(power).copy(), answer['state']))
            return answer

    class TimedRemaining(model.RemainingRegionModel):
        def __init__(self, equations, budget, bounds, total_bound, cuts, inner_halfspaces, tau, **options):
            tick = perf_counter()
            self.linear_outer = (args.variant.startswith('linear_outer') and
                equations.method == 'socp' and options.get('mode') == 'light')
            if self.linear_outer:
                equations = model.PlanningEquations(equations.network, 'linear')
                options = dict(options, mode='physical')
            constructor = tight_remaining if args.variant in ('tight_m', 'combined', 'linear_outer_combined') else model.RemainingRegionModel.__init__
            constructor(self, equations, budget, bounds, total_bound, cuts, inner_halfspaces, tau, **options)
            self.build_seconds = perf_counter()-tick
            self.input = dict(method=equations.method, budget=budget, bounds=bounds,
                total_bound=total_bound, cuts=cuts, inner_halfspaces=inner_halfspaces, tau=tau, **options)

        def solve(self, *values, **options):
            tick = perf_counter()
            answer = super().solve(*values, **options)
            if self.linear_outer:
                answer['feasible'] = False  # An LP witness is not an SOCP certificate.
            remaining_calls.append(dict(build_seconds=self.build_seconds,
                seconds=perf_counter()-tick, runtime=self.model.Runtime,
                status=self.model.Status, nodes=self.model.NodeCount,
                binary_variables=self.model.NumBinVars, linear_rows=self.model.NumConstrs,
                cones=self.model.NumQConstrs, upper=getattr(self, 'tight_upper', 4.),
                complete=answer['complete'], bound=answer['bound']))
            write_json(directory/'last_residual_input.json', self.input)
            return answer

    application.PlanningModel = TimedModel
    application.PlanningSP = TimedSP
    application.RemainingRegionModel = TimedRemaining
    if args.variant in ('downward', 'combined', 'linear_outer_combined'):
        application.RegionState = DownwardRegion

    def progress(event, **data):
        nonlocal last_checkpoint
        now = perf_counter()
        if event == 'cut':
            cuts.append(data['cut'].copy())
        if event in ('phase_start', 'mp_start', 'query_end', 'residual_start', 'region_end') or now-last_checkpoint > 10.:
            record = dict(event=event, elapsed=now-started, message=data.get('message'),
                          counts=data['counts_algorithm'], schemes=len(data['records']))
            events.append(record)
            write_json(directory/'progress.json', dict(case=label, latest=record, events=events))
            last_checkpoint = now
            if event in ('phase_start', 'query_end', 'region_end'):
                print(json.dumps(dict(case=label, **record), ensure_ascii=False), flush=True)

    print(json.dumps(dict(case=label, event='start', time_limit=args.time_limit)), flush=True)
    result = application.build_continuous_region(network, args.method, args.budget, bounds,
        threads=args.threads, tau=.002, time_limit=args.time_limit,
        residual_mode=args.residual_mode, progress=progress)
    seconds = perf_counter()-started
    data = dict(case=label, variant=args.variant, candidate_count=args.count, budget=args.budget,
        method=args.method, threads=args.threads, time_limit=args.time_limit, bounds=bounds,
        source_hashes=before, source_unchanged=before == fingerprints(),
        measured_seconds=seconds, mp_calls=mp_calls, sp_calls=sp_calls,
        remaining_calls=remaining_calls, cuts=cuts, result=result)
    write_json(directory/'result.json', data)
    if args.method != 'linear':
        data['audit'] = audit(network, result, certificates, cuts)
    write_json(directory/'result.json', data)
    print(json.dumps(dict(case=label, event='done', status=result['status'], seconds=seconds,
        counts=result['counts'], maximum=result['max_total'], upper=result['max_total_bound'],
        inner_schemes=len(result['inner']), residual_calls=len(remaining_calls),
        audit=None if 'audit' not in data else {k: v for k, v in data['audit'].items()
            if k not in ('ac_failures', 'ac_unknown')})), flush=True)


def replay(args):
    """Paired repetitions on identical SP inputs and saved residual models."""
    from unittest.mock import patch
    output = Path(args.output)
    network = Case33(upgrade_count=args.count)
    equations = model.PlanningEquations(network, 'linear')
    x = network.encode_plan(network.initial_plan)
    powers = np.random.default_rng(20260923).random((64, 3))*[1500., 6000., 1500.]
    implementations = {'baseline': model.PlanningSP.solve, 'linear_single': single_linear_solve()}
    implementations['baseline'](model.PlanningSP(equations, threads=args.threads), x, powers[0])
    runs, references, mismatches, maximum_cut_difference = [], {}, [], 0.
    original_optimize = gp.Model.optimize
    for repetition in range(5):
        order = ('baseline', 'linear_single') if repetition % 2 == 0 else ('linear_single', 'baseline')
        for variant in order:
            optimize_calls = 0

            def optimize(problem, *values, **options):
                nonlocal optimize_calls
                optimize_calls += 1
                return original_optimize(problem, *values, **options)

            oracle = model.PlanningSP(equations, threads=args.threads)
            started = perf_counter()
            with patch.object(gp.Model, 'optimize', optimize):
                answers = [implementations[variant](oracle, x, p) for p in powers]
            runs.append(dict(repetition=repetition+1, variant=variant,
                seconds=perf_counter()-started, optimize_calls=optimize_calls,
                feasible=sum(a['feasible'] for a in answers),
                cuts=sum(a['cut'] is not None for a in answers)))
            for index, answer in enumerate(answers):
                if index not in references:
                    references[index] = answer
                other = references[index]
                if (answer['feasible'] != other['feasible'] or
                    (answer['cut'] is None) != (other['cut'] is None)):
                    mismatches.append([repetition, variant, index])
                if answer['cut'] is not None and other['cut'] is not None:
                    maximum_cut_difference = max(maximum_cut_difference,
                        float(np.max(np.abs(answer['cut']-other['cut']))))
    linear = dict(samples=powers, repetitions=5, runs=runs, mismatches=mismatches,
        maximum_cut_difference=maximum_cut_difference,
        medians={v: float(np.median([r['seconds'] for r in runs if r['variant'] == v])) for v in implementations})
    write_json(output/'linear_sp_replay.json', linear)
    print(json.dumps(dict(experiment='linear_sp_replay', **{k: v for k, v in linear.items() if k not in ('samples', 'runs')})), flush=True)

    residual = []
    for budget in ('0', '2', 'inf'):
        path = output/f'n{args.count}_b{budget}_socp_baseline_t{args.threads}_r1'/'last_residual_input.json'
        if not path.exists():
            continue
        data = json.loads(path.read_text(encoding='utf-8'))
        eq = model.PlanningEquations(network, data['method'])
        bounds = np.asarray(data['bounds'])
        inners = [np.asarray(a) for a in data['inner_halfspaces']]
        box = initial_polytope(bounds, data['total_bound'])
        extrema = [box@((1-data['tau'])*p[:, :3]).T+p[:, 3] for p in inners]
        upper = min((float(a.max())+1e-10 for a in extrema), default=4.)
        minimum_inactive_slack = min((float(np.min(a.min(axis=0)+
            np.maximum(0., upper-a.min(axis=0)+1e-10)-upper)) for a in extrema), default=0.)

        class TightRemaining(model.RemainingRegionModel):
            __init__ = tight_remaining

        for repetition in range(3):
            order = ('baseline', 'tight_m') if repetition % 2 == 0 else ('tight_m', 'baseline')
            for variant in order:
                started = perf_counter()
                cls = TightRemaining if variant == 'tight_m' else model.RemainingRegionModel
                problem = cls(eq, data['budget'], bounds, data['total_bound'], data['cuts'],
                    inners, data['tau'], mode=data['mode'], threads=args.threads)
                with problem.model:
                    answer = problem.solve(GEOMETRY_TOL, time_limit=10.)
                    row = dict(budget=budget, repetition=repetition+1, variant=variant,
                        seconds=perf_counter()-started, runtime=problem.model.Runtime,
                        nodes=problem.model.NodeCount, status=problem.model.Status,
                        complete=answer['complete'], bound=answer['bound'], witness=answer['x'] is not None,
                        analytic_gamma_upper=upper, minimum_inactive_slack=minimum_inactive_slack)
                residual.append(row)
        print(json.dumps(dict(experiment='residual_replay', budget=budget,
            medians={v:float(np.median([r['seconds'] for r in residual if r['budget']==budget and r['variant']==v]))
                     for v in ('baseline', 'tight_m')})), flush=True)
    write_json(output/'residual_replay.json', residual)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--prepare', action='store_true')
    parser.add_argument('--replay', action='store_true')
    parser.add_argument('--count', type=int, default=8)
    parser.add_argument('--budget', type=float, default=2.)
    parser.add_argument('--variant', choices=('baseline', 'tight_m', 'downward', 'linear_single',
        'combined', 'linear_outer', 'linear_outer_combined'), default='baseline')
    parser.add_argument('--method', choices=('socp', 'linear', 'hybrid'), default='socp')
    parser.add_argument('--threads', type=int, default=20)
    parser.add_argument('--time-limit', type=float, default=120.)
    parser.add_argument('--residual-mode', choices=('auto', 'light', 'physical'), default='auto')
    parser.add_argument('--repeat', type=int, default=1)
    parser.add_argument('--output', default=str(OUTPUT))
    arguments = parser.parse_args()
    with threadpool_limits(limits=1):
        (prepare if arguments.prepare else replay if arguments.replay else run)(arguments)
