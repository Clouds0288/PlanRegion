"""case33 四候选的测试侧 16 组合独立核对；不参与正式构域或计时。"""
from argparse import ArgumentParser
from itertools import product
from pathlib import Path
from time import perf_counter
import json

import numpy as np
from threadpoolctl import threadpool_limits

from Network.case33bw import Case33
from plot import sample_region, region_metrics
from region import GEOMETRY_TOL
from model import PlanningEquations, PlanningModel
from vertify import ACPowerFlow
from plot import region_view, save_method_comparison
from tests.reference import _solve, _dispatch_equations, dispatch_support, fixed_topology, upgrade_plan
from plot import BenchmarkResult, METHODS


def independent_feasible(e, method, points):
    points = np.asarray(points).reshape(-1, 3)
    if method == 'linear':
        return np.min(e.linear_c[:, None]+e.linear_F@points.T, axis=0) >= -2e-8
    cones = [('linear', e.linear_count)]+[('soc', size) for size in e.sizes]
    answer = []
    for point in points:
        result = _solve(np.ones(e.G.shape[1]), -e.G, e.c+e.F@point, cones)
        if result is None:
            answer.append(False)
            continue
        # 原始消元方程独立验证内域证书，不用求解器状态替代残量检查。
        residual = e.c+e.F@point+e.G@result.x
        good = residual[:e.linear_count].min() >= -2e-8
        offset = e.linear_count
        for size in e.sizes:
            good &= residual[offset]-np.linalg.norm(residual[offset+1:offset+size]) >= -2e-8
            offset += size
        answer.append(good)
    return np.asarray(answer)


def audit(folder):
    started = perf_counter()
    result = BenchmarkResult.load(folder)
    if result.metadata.get('upgrade_count', result.metadata.get('candidate_count')) != 4:
        raise ValueError('This exhaustive audit is scoped to four upgrades on the initial tree')
    network, bounds = fixed_topology(Case33(upgrade_count=4)), result.bounds
    if result.metadata.get('schema_version', 0) >= 3 and result.metadata.get('network_fingerprint') != network.fingerprint:
        raise ValueError('Upgrade-only enumeration cannot certify a reconfiguration result')
    plans = [upgrade_plan(network, choice) for choice in product((0, 1), repeat=4)]
    designs = [network.tree(network.encode_plan(plan)) for plan in plans]
    equations = {tuple(sorted(plan.items())): _dispatch_equations(d) for plan, d in zip(plans, designs)}
    directions = np.vstack([np.ones(3), np.eye(3), np.array([[1, 6, 1], [3, 1, 4], [1, 1, 5]])])
    references = {method: [[dispatch_support(d, method, normal) for normal in directions] for d in designs]
                  for method in ('linear', 'socp')}
    report = dict(passed=True, designs=16, support_queries=16*2*len(directions), regions=[],
                  inner_vertices_checked=0, ac_grid_points=0)
    for row in result.metadata['continuous']:
        method = 'linear' if row['method'] == 'linear' else 'socp'
        budget = np.inf if row['budget'] is None else row['budget']
        assert row['status'] == 'certified', (row['method'], budget, row['status'])
        assert row['coverage_bound'] is None or row['coverage_bound'] <= GEOMETRY_TOL
        gap = region_metrics(row, bounds)['volume_gap']
        assert -.000001 <= gap <= (.00001 if method == 'linear' else .00601), gap
        indices = [i for i, d in enumerate(designs) if d.cost <= budget]
        support = [references[method][i] for i in indices]
        maximum = max(values[0]['value'] for values in support)
        assert abs(maximum-row['max_total']) <= .003, (maximum, row['max_total'])
        boundary = np.asarray([value['p'] for values in support for value in values])
        assert np.all(sample_region(row, boundary, bounds) != -1), (row['method'], budget, 'boundary outside outer')
        count = 0
        for certificate in region_view(row, bounds)['geometry']:
            points = np.asarray(certificate['inner']['vertices']).reshape(-1, 3)
            assert independent_feasible(equations[tuple(sorted(certificate['choice'].items()))], method, points).all(), (row['method'], budget, certificate['choice'])
            count += len(points)
        report['inner_vertices_checked'] += count
        report['regions'].append(dict(method=row['method'], budget=row['budget'], passed=True,
                                      reference_max_total=maximum, maximum_error_kw=abs(maximum-row['max_total']),
                                      inner_vertices=count, boundary_points=len(boundary)))
        print('Verified', row['method'], row['budget'], 'vertices', count, flush=True)
    # 独立 AC 对全部 16 个方案批量求解，与正式方案搜索/单调传播得到的每个中心交叉核对。
    n = result.metadata['divisions']
    points = (np.indices((n,)*3).reshape(3, -1).T+.5)*bounds/n
    labels = []
    for design in designs:
        oracle = ACPowerFlow(design, threads=1)
        try:
            states = oracle.classify(points)
            for i in np.flatnonzero(states == 0):
                states[i] = oracle.global_status(points[i], None)
        finally:
            oracle.close()
        labels.append(states)
    labels = np.asarray(labels)
    for j, budget in enumerate(result.budgets):
        chosen = labels[[i for i, d in enumerate(designs) if d.cost <= budget]]
        expected = np.where((chosen == 1).any(axis=0), 1, np.where((chosen == -1).all(axis=0), -1, 0))
        np.testing.assert_array_equal(expected, result.states[METHODS.index('ac'), j].reshape(-1))
    report['ac_grid_points'] = len(points)
    report['ac_design_point_checks'] = int(labels.size)
    report['seconds'] = perf_counter()-started
    Path(folder, 'validation.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    result.metadata['audit'] = report
    result.save(folder)
    save_method_comparison(result, folder)
    return report


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('folder', type=Path)
    args = parser.parse_args()
    with threadpool_limits(limits=1):
        print(json.dumps(audit(args.folder), ensure_ascii=False, indent=2))
