"""已有网架优先：孔洞、低维剩余和真实 SOCP 调度证书。"""
import unittest

import numpy as np
from threadpoolctl import threadpool_limits

from experiments.fourbus_known_first import remaining_cells, known_candidates, build_known_first_region
from model import GridPhysics, RemainingRegionModel
from monitor import RunMonitor
from Network.four_bus_five_corridor import FourBus
from region import RegionState, halfspaces, contains, GEOMETRY_TOL


SQUARE = np.array([[0., 0.], [1., 0.], [1., 1.], [0., 1.]])


class KnownFirstTests(unittest.TestCase):
    def test_vertices_covered_still_finds_interior_hole(self):
        region = RegionState(np.ones(2), 2., 0.)
        polygons = [SQUARE*[.4, 1.], SQUARE*[.4, 1.]+[.6, 0.],
                    SQUARE*[1., .4], SQUARE*[1., .4]+[0., .6]]
        for index, poly in enumerate(polygons):
            region.add_scheme([index], {}, 0.)
            region.add_point([index], poly)
        self.assertTrue(all(owner is not None for owner in region.covering_schemes(SQUARE)))
        candidates, interior = known_candidates(region)
        self.assertTrue(interior)
        self.assertEqual(len(candidates), 4)
        for _, point in candidates:
            np.testing.assert_allclose(point, [.5, .5], atol=2e-8)
            self.assertIsNone(region.covering_schemes([point])[0])

    def test_union_coverage_needs_no_single_scheme_to_cover_whole_outer(self):
        left = SQUARE*[.6, 1.]
        right = SQUARE*[.6, 1.]+[.4, 0.]
        self.assertEqual(remaining_cells(SQUARE, [halfspaces(left), halfspaces(right)]), [])

    def test_degenerate_segment_gap_is_not_dropped_as_zero_area(self):
        segment = np.array([[0., .5], [1., .5]])
        cells = remaining_cells(segment, [halfspaces(SQUARE*[.4, 1.]),
                                         halfspaces(SQUARE*[.4, 1.]+[.6, 0.])])
        self.assertTrue(cells)
        self.assertTrue(any(contains([[.5, .5]], halfspaces(cell))[0] for cell in cells))

    def test_tau_uses_same_radial_coverage_rule(self):
        region = RegionState(np.ones(2), 2., .005)
        region.add_scheme([0], {}, 0.)
        region.add_point([0], .995*SQUARE)
        candidates, interior = known_candidates(region)
        self.assertTrue(interior)
        self.assertEqual(candidates, [])

    def test_fourbus_global_calls_follow_complete_known_coverage(self):
        network = FourBus(load_nodes=(1, 2))
        bounds = np.full(2, network.power_limit)
        monitor = RunMonitor()
        with threadpool_limits(limits=1):
            result = build_known_first_region(network, 'socp', 20000., bounds, threads=1, progress=monitor)
            self.assertEqual(result['status'], 'certified')
            self.assertLessEqual(result['coverage_bound'], GEOMETRY_TOL)
            calls = 0
            for index, item in enumerate(monitor.history):
                if item['patch']['event'] != 'residual_start':
                    continue
                snapshot = monitor.frame(index)
                equations = [halfspaces(np.asarray(row['inner'])/bounds)
                             for row in snapshot['schemes'].values() if row['inner']]
                cuts = [np.asarray(row['cut']) for row in snapshot.get('cut_history', {}).values()]
                for row in snapshot['schemes'].values():
                    # 独立 MILP 审计：固定已知 x，几何 max-gamma 的全局上界也必须 <= 容差。
                    problem = RemainingRegionModel(GridPhysics(network, 'socp'), 20000., bounds,
                        snapshot['total_bound'], cuts, equations, result['tau'],
                        axis_bounds=snapshot['axis_bounds'], mode='light', threads=1)
                    with problem.model:
                        problem.model.addConstr(problem.problem.x == np.asarray(row['x']))
                        answer = problem.solve(GEOMETRY_TOL)
                    self.assertTrue(answer['complete'])
                calls += 1
            self.assertEqual(calls, result['counts']['global_search'])


if __name__ == '__main__':
    unittest.main()
