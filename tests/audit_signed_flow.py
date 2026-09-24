"""Independent physics, before/after results and paired timings for signed flows."""
from argparse import ArgumentParser
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, product
import json
from pathlib import Path
import sys
from time import perf_counter

import numpy as np
from threadpoolctl import threadpool_limits

from Network.case33bw import Case33
from Network.four_bus_five_corridor import FourBus
import model
from tests.reference import dispatch_support, fixed_topology


def query_options(name, network):
    powers = ([10., 10., 10.], [25., 15., 20.], [30., 30., 30.]) if name == 'fourbus' else (
        [100., 800., 150.], [250., 1800., 350.], [300., 3000., 500.])
    return [dict(budget=b) for b in network.budgets]+[dict(power=p) for p in powers]


def direct_checks(before, after):
    assert len(before['queries']) == len(after['queries']) == 42
    differences = []
    for old, new in zip(before['queries'], after['queries']):
        assert (old['case'], old['method'], old['query']) == (new['case'], new['method'], new['query'])
        assert old['status'] == new['status'] == 'optimal'
        assert old['feasible'] and new['feasible']
        assert new['margin'] >= -model.PLANNING_TOL
        delta = abs(old['objective']-new['objective'])
        assert delta < (1e-7 if new['query'].startswith('power') else .002)
        if new['method'] == 'socp':
            assert new['ac_status'] == 1
        differences.append(delta)
    return dict(queries=42, socp_ac_passed=21, maximum_objective_difference=max(differences))


def all_fourbus_trees():
    network = FourBus()
    trees, checks, maximum = 0, 0, 0.
    equations = {method: model.PlanningEquations(network, method) for method in ('linear', 'socp')}
    for corridors in combinations(network.corridors, network.n):
        reached = {network.root}
        for _ in range(network.n):
            for corridor in corridors:
                if reached.intersection(corridor.endpoints):
                    reached.update(corridor.endpoints)
        if len(reached) != network.n+1:
            continue
        trees += 1
        for types in product(*(c.types for c in corridors)):
            plan = dict.fromkeys(c.id for c in network.corridors)
            plan.update((c.id, t.id) for c, t in zip(corridors, types))
            for method, e in equations.items():
                problem = model.PlanningModel(e, fixed_plan=plan, threads=1)
                with problem.model:
                    answer = problem.solve(time_limit=60.)
                reference = dispatch_support(network.tree(network.encode_plan(plan)), method, np.ones(3))
                delta = abs(answer['objective']-reference['value'])
                assert answer['status'] == 'optimal' and answer['feasible']
                assert delta < .002, (method, plan, delta)
                maximum = max(maximum, delta)
                checks += 1
    assert trees == 8 and checks == 432
    return dict(trees=trees, designs=checks//2, lp_socp_checks=checks,
                maximum_boundary_difference_kw=maximum)


def region_checks(before, after, networks):
    assert len(before['regions']) == len(after['regions']) == 9
    rows, queries = [], 0
    for name, network in networks:
        old_rows = [row for row in before['regions'] if row['case'] == name]
        new_rows = [row for row in after['regions'] if row['case'] == name]
        references = {}
        for old, new in zip(old_rows, new_rows):
            assert old['method'] == new['method']
            for row in (old, new):
                method = 'socp' if row['method'] == 'hybrid' else row['method']
                key = method, tuple(row['bounds']), row['budget']
                if key not in references:
                    e = model.PlanningEquations(network, method)
                    points = (np.indices((3,)*3).reshape(3, -1).T+.5)*row['bounds']/3
                    labels = []
                    for point in points:
                        problem = model.PlanningModel(e, power=point, budget=row['budget'], threads=1)
                        with problem.model:
                            answer = problem.solve(time_limit=60.)
                        assert answer is None or answer['status'] == 'optimal'
                        labels.append(-1 if answer is None else 1)
                        queries += 1
                    references[key] = np.array(labels)
                expected = references[key]
                actual = np.array(row['grid_labels'])
                known = actual != 0
                np.testing.assert_array_equal(actual[known], expected[known])
            rows.append(dict(case=name, method=new['method'], before_status=old['status'],
                             after_status=new['status'], before_seconds=old['seconds'],
                             after_seconds=new['seconds'],
                             before_bounds=old['bounds'], after_bounds=new['bounds'],
                             classified=int(np.count_nonzero(new['grid_labels'])), samples=len(expected)))
    return dict(independent_point_queries=queries, regions=rows)


def paired_timings(old_model, old_networks, networks, repeats):
    rows = []
    for name, network in networks:
        for method in ('linear', 'socp'):
            modules = dict(before=old_model, after=model)
            equations = {key: module.PlanningEquations(old_networks[name] if key == 'before' else network, method)
                         for key, module in modules.items()}
            samples = {key: [] for key in modules}
            sizes = {}
            for iteration in range(repeats+1):
                order = ('before', 'after') if iteration % 2 == 0 else ('after', 'before')
                for key in order:
                    module, e = modules[key], equations[key]
                    build, solve, nodes = 0., 0., 0.
                    for options in query_options(name, network):
                        started = perf_counter()
                        problem = module.PlanningModel(e, threads=1, **options)
                        with problem.model:
                            problem.model.update()
                            build += perf_counter()-started
                            sizes[key] = dict(binary=problem.model.NumBinVars, variables=problem.model.NumVars,
                                              constraints=problem.model.NumConstrs, cones=problem.model.NumQConstrs,
                                              physics_variables=problem.state.shape[0],
                                              connectivity=sum(v.VarName.startswith(('connectivity', 'f['))
                                                               for v in problem.model.getVars()))
                            started = perf_counter()
                            answer = problem.solve(time_limit=60.)
                            solve += perf_counter()-started
                            nodes += problem.model.NodeCount
                            assert answer['status'] == 'optimal' and answer['feasible']
                    if iteration:
                        samples[key].append(dict(build_seconds=build, solve_seconds=solve,
                                                 total_seconds=build+solve, nodes=nodes))
            median = {key: float(np.median([row['total_seconds'] for row in values]))
                      for key, values in samples.items()}
            rows.append(dict(case=name, method=method, queries_per_repeat=7, repeats=repeats,
                             samples=samples, median_seconds=median,
                             after_over_before=median['after']/median['before'], sizes=sizes))
    return rows


def run(directory, repeats):
    before = json.loads((directory/'before.json').read_text(encoding='utf-8'))
    after = json.loads((directory/'after.json').read_text(encoding='utf-8'))
    if after.get('case33_topology_scope') != 'initial_tree':
        raise ValueError('Paired comparisons require initial-tree Case33 results')
    source = directory/'source_before'
    for name, path in (('network_before', source/'Network/__init__.py'),
                       ('network_before.case33bw', source/'Network/case33bw.py'),
                       ('network_before.four_bus_five_corridor', source/'Network/four_bus_five_corridor.py')):
        specification = spec_from_file_location(name, path)
        module = module_from_spec(specification)
        sys.modules[name] = module
        specification.loader.exec_module(module)
    old_case = sys.modules['network_before.case33bw'].Case33
    old_networks = {'fourbus': sys.modules['network_before.four_bus_five_corridor'].FourBus(),
                    'case33_8': old_case(candidate_count=8), 'case33_16': old_case(candidate_count=16)}
    spec = spec_from_file_location('model_before_signed_flow', directory/'source_before/model.py')
    old_model = module_from_spec(spec)
    spec.loader.exec_module(old_model)
    networks = [('fourbus', FourBus()), ('case33_8', fixed_topology(Case33(upgrade_count=8))),
                ('case33_16', fixed_topology(Case33(upgrade_count=16)))]
    report = dict(threads=1, direct=direct_checks(before, after))
    with threadpool_limits(limits=1):
        report['case33_topology_scope'] = 'initial_tree'
        report['fourbus_physics'] = all_fourbus_trees()
        print(json.dumps(report), flush=True)
        report['region_audit'] = region_checks(before, after, networks)
        print(json.dumps(report['region_audit']), flush=True)
        report['paired_timings'] = paired_timings(old_model, old_networks, networks, repeats)
    (directory/'audit.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report['paired_timings']), flush=True)


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, default=Path('results/signed_flow'))
    parser.add_argument('--repeats', type=int, default=3)
    args = parser.parse_args()
    run(args.directory, args.repeats)
