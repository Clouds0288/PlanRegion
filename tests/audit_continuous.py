"""case33 四候选的测试侧 16 组合独立核对；不参与正式构域或计时。"""
from argparse import ArgumentParser
from itertools import product
from pathlib import Path
from time import perf_counter
import json

import clarabel
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits

from Network.case33bw import Case33
from plot import sample_region, region_metrics
from region import GEOMETRY_TOL
from model import PlanningEquations, PlanningModel
from vertify import ACPowerFlow
from plot import region_view, save_method_comparison
from tests.reference import _dispatch_equations, dispatch_support
from plot import BenchmarkResult, METHODS


def independent_feasible(e, method, points):
    points = np.asarray(points).reshape(-1, 3)
    if method == 'linear':
        return np.min(e.linear_c[:, None]+e.linear_F@points.T, axis=0) >= -2e-8
    settings = clarabel.DefaultSettings()
    settings.verbose = False
    settings.max_threads = 1
    settings.tol_gap_abs = settings.tol_gap_rel = settings.tol_feas = 1e-12
    n = e.G.shape[1]
    cones = ([clarabel.NonnegativeConeT(e.linear_count)] + [clarabel.SecondOrderConeT(k) for k in e.sizes]
             + [clarabel.NonnegativeConeT(1)])
    matrix = -np.vstack([np.c_[e.G, e.relax], np.r_[np.zeros(n), 1.]])
    def physically_valid(current):
        residual = e.c+e.F@p+e.G@np.asarray(current)
        good = np.isfinite(residual).all() and residual[:e.linear_count].min() >= -2e-8
        offset = e.linear_count
        for size in e.sizes:
            good &= residual[offset]-np.linalg.norm(residual[offset+1:offset+size]) >= -2e-8
            offset += size
        return good
    answer = []
    for p in points:
        solver = clarabel.DefaultSolver(sparse.csc_matrix((n+1, n+1)), np.r_[np.zeros(n), 1.],
                                        sparse.csc_matrix(matrix), np.r_[e.c+e.F@p, 0.], cones, settings)
        result = solver.solve()
        # 只检查未松弛的原始独立方程，不能凭 eta 或求解器状态接受顶点。
        good = physically_valid(result.x[:n])
        if not good:
            # 独立消元模型的边界 phase I 也可能停滞；直接求原约束中的最小电流。
            # 不调用正式 PlanningSP，不调整点或放宽同一个原始余量阈值。
            direct = clarabel.DefaultSolver(sparse.csc_matrix((n, n)), np.ones(n),
                                            sparse.csc_matrix(-e.G), e.c+e.F@p,
                                            cones[:-1], settings).solve()
            good = physically_valid(direct.x)
        answer.append(good)
    return np.asarray(answer)


def audit(folder):
    started = perf_counter()
    result = BenchmarkResult.load(folder)
    if result.metadata['candidate_count'] != 4:
        raise ValueError('This exhaustive audit is scoped to four candidate lines')
    network, bounds = Case33(candidate_count=4), result.bounds
    plans = [network.initial_plan | {corridor.id: corridor.types[k].id
             for corridor, k in zip(network.planning_corridors, choice)}
             for choice in product((0, 1), repeat=4)]
    designs = [network.design(plan) for plan in plans]
    equations = {tuple(sorted(d.x.items())): _dispatch_equations(d) for d in designs}
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
