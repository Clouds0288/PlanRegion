"""Audit saved before/after runs and compare their physical coefficient arrays."""
from argparse import ArgumentParser
from importlib.util import module_from_spec, spec_from_file_location
import json
from pathlib import Path
import sys

import numpy as np
from threadpoolctl import threadpool_limits

from Network.case33bw import Case33
from Network.four_bus_five_corridor import FourBus
from model import PlanningEquations, PlanningModel


def load_source(name, path):
    spec = spec_from_file_location(name, path)
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def audit(folder):
    before, after = [json.loads((folder/name).read_text(encoding='utf-8'))
                     for name in ('before.json', 'after.json')]
    source = folder/'source_before'
    load_source('corridor_before', source/'Network/__init__.py')
    old_four = load_source('corridor_before.four_bus_five_corridor', source/'Network/four_bus_five_corridor.py')
    old_case = load_source('corridor_before.case33bw', source/'Network/case33bw.py')
    old_model = load_source('corridor_model_before', source/'model.py')
    networks = {'fourbus': FourBus(), 'case33_8': Case33(candidate_count=8),
                'case33_16': Case33(candidate_count=16)}
    old_networks = {'fourbus': old_four.FourBus(), 'case33_8': old_case.Case33(candidate_count=8),
                    'case33_16': old_case.Case33(candidate_count=16)}
    fields = ('edges', 'r', 'reactance', 'cost', 'selected', 'senders', 'receivers',
              'c', 'F', 'G', 'upper', 'box_constant', 'box_selection', 'T', 'balance', 'linked', 'relax')
    # 历史快照仍使用旧字段；比较时转换到当前 Ax+By+Cp≼_K b 的约定。
    renamed_fields = {'edges': ('type_corridor', 1), 'selected': ('decision_types', 1),
                      'c': ('b', 1), 'F': ('C', -1), 'G': ('B', -1),
                      'upper': ('y_ub_global', 1), 'box_constant': ('ub_const', 1),
                      'box_selection': ('ub_x', 1), 'relax': ('relax_direction', 1)}
    report = dict(passed=True, matrix_fields=fields, matrix_checks=[], queries=[], regions=[])
    for name, network in networks.items():
        for method in ('linear', 'socp'):
            left, right = old_model.PlanningEquations(old_networks[name], method), PlanningEquations(network, method)
            for field in fields:
                new_field, sign = renamed_fields.get(field, (field, 1))
                np.testing.assert_array_equal(sign*getattr(left, field), getattr(right, new_field),
                                              err_msg=f'{name}/{method}/{field}')
            report['matrix_checks'].append(dict(case=name, method=method, arrays_identical=True))
    assert len(before['queries']) == len(after['queries']) == 42
    for left, right in zip(before['queries'], after['queries']):
        assert all(left[k] == right[k] for k in ('case', 'method', 'query', 'binary_count'))
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
