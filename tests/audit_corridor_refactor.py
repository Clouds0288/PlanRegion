"""Compare saved runs on identical topology scopes and verify independent physics."""
from argparse import ArgumentParser
import json
from pathlib import Path

import numpy as np
from threadpoolctl import threadpool_limits

from Network.case33bw import Case33
from Network.four_bus_five_corridor import FourBus
from model import PlanningEquations, PlanningModel
from tests.reference import dispatch_support, fixed_topology


def audit(folder):
    before, after = [json.loads((folder/name).read_text(encoding='utf-8'))
                     for name in ('before.json', 'after.json')]
    if after.get('case33_topology_scope') != 'initial_tree':
        raise ValueError('Before/after comparison requires the same initial-tree scope for Case33')
    networks = {'fourbus': FourBus(), 'case33_8': fixed_topology(Case33(upgrade_count=8)),
                'case33_16': fixed_topology(Case33(upgrade_count=16))}
    report = dict(passed=True, topology_scope='initial_tree', physics_checks=[], queries=[], regions=[])
    for name, network in networks.items():
        for method in ('linear', 'socp'):
            x = network.encode_plan(network.initial_plan)
            reference = dispatch_support(network.tree(x), method, np.ones(3))
            problem = PlanningModel(PlanningEquations(network, method), fixed_plan=network.initial_plan, threads=1)
            with problem.model:
                answer = problem.solve(time_limit=60.)
            assert answer['feasible'] and answer['status'] == 'optimal'
            delta = abs(answer['objective']-reference['value'])
            assert delta < .002
            report['physics_checks'].append(dict(case=name, method=method, boundary_difference_kw=delta))
    assert len(before['queries']) == len(after['queries']) == 42
    for left, right in zip(before['queries'], after['queries']):
        assert all(left[k] == right[k] for k in ('case', 'method', 'query'))
        assert left['status'] == right['status'] == 'optimal'
        assert left['feasible'] and right['feasible'] and right['margin'] >= -1e-8
        error = abs(left['objective']-right['objective'])
        assert error <= (1e-7 if right['query'].startswith('power_') else .002)
        if right['method'] == 'socp':
            assert right['ac_status'] == 1
        report['queries'].append(dict(case=right['case'], method=right['method'], query=right['query'],
                                      before=left['objective'], after=right['objective'], absolute_error=error))
    assert len(before['regions']) == len(after['regions']) == 9
    references = {}
    for left, right in zip(before['regions'], after['regions']):
        assert all(left[k] == right[k] for k in ('case', 'method', 'budget', 'bounds'))
        error = abs(left['max_total']-right['max_total'])
        assert error <= .002
        a, b = np.asarray(left['grid_labels']), np.asarray(right['grid_labels'])
        known = (a != 0) & (b != 0)
        np.testing.assert_array_equal(a[known], b[known])
        method = 'linear' if right['method'] == 'linear' else 'socp'
        key = right['case'], method
        if key not in references:
            equations = PlanningEquations(networks[right['case']], method)
            points = (np.indices((3,)*3).reshape(3, -1).T+.5)*np.asarray(right['bounds'])/3
            labels = []
            for point in points:
                problem = PlanningModel(equations, power=point, budget=right['budget'], threads=1)
                with problem.model:
                    answer = problem.solve(time_limit=60.)
                assert answer is None or (answer['status'] == 'optimal' and answer['feasible'])
                labels.append(-1 if answer is None else 1)
            references[key] = np.asarray(labels)
        np.testing.assert_array_equal(b[b != 0], references[key][b != 0])
        report['regions'].append(dict(case=right['case'], method=right['method'],
            before_status=left['status'], after_status=right['status'], maximum_error_kw=error,
            before_seconds=left['seconds'], after_seconds=right['seconds'],
            known_grid_points=int(np.count_nonzero(b)), unknown_grid_points=int(np.count_nonzero(b == 0)),
            grid_disagreements=0))
    report['direct_grid_queries'] = 27*len(references)
    report['maximum_objective_error'] = max(r['absolute_error'] for r in report['queries'])
    (folder/'audit.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('folder', type=Path)
    args = parser.parse_args()
    with threadpool_limits(limits=1):
        audit(args.folder)
