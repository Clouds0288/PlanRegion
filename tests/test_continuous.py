"""连续几何、MP1 自由分配、未发现方案覆盖及完整历史回归。"""
from io import StringIO
from itertools import product
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
import base64
import gzip
import json
import re
import unittest

import numpy as np
from threadpoolctl import threadpool_limits

from main import Case33, build_continuous_region
from model import evaluation_bounds, RemainingRegionModel
from tests.planning_checks import margin
from model import PlanningEquations, PlanningModel
from plot import sample_region
from plot import RunMonitor, json_value, region_metrics
from region import halfspaces, contains
from plot import union_volume
from region import polytope_volume, clip_polytope
from tests.reference import dispatch_support, fixed_topology, upgrade_plan


class ContinuousTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.threads = threadpool_limits(limits=1)
        cls.network = fixed_topology(Case33(upgrade_count=4))
        cls.bounds = evaluation_bounds(cls.network, threads=1)


    def test_deadline_preserves_unknown_domain(self):
        result = build_continuous_region(self.network, 'socp', 0., self.bounds, time_limit=0., threads=1)
        self.assertEqual(result['status'], 'time_limit')
        self.assertEqual(result['inner'], [])
        self.assertGreater(region_metrics(result, self.bounds)['outer_volume'], 0.)
        self.assertEqual(result['counts']['sp'], 0)

    def test_complete_query_reuses_verified_state_without_sp(self):
        e = PlanningEquations(self.network, 'socp')
        with patch('model.PlanningSP.solve', side_effect=AssertionError('Redundant SP')):
            problem = PlanningModel(e, budget=0., threads=1)
            with problem.model:
                answer = problem.solve()
        self.assertTrue(answer['feasible'])
        self.assertGreaterEqual(margin(e, answer['x'], answer['p'], answer['state']), -1e-8)

    def test_residual_modes_agree_on_small_certified_union(self):
        points = (np.indices((5,)*3).reshape(3, -1).T+.5)*self.bounds/5
        results = [build_continuous_region(self.network, 'socp', 0., self.bounds,
                                          residual_mode=mode, threads=1) for mode in ('light', 'physical')]
        for result in results:
            self.assertEqual(result['status'], 'certified')
            self.assertGreater(result['timing']['total_seconds'], 0.)
        a, b = [sample_region(r, points, self.bounds) for r in results]
        self.assertFalse(np.any(a*b == -1))
        self.assertLess(abs(region_metrics(results[0], self.bounds)['inner_volume']-region_metrics(results[1], self.bounds)['inner_volume'])/region_metrics(results[0], self.bounds)['inner_volume'], .01)

    @classmethod
    def tearDownClass(cls):
        cls.threads.restore_original_limits()

    def test_union_preserves_overlap_and_nonconvex_gap(self):
        cube = np.array(list(product((0., 1.), repeat=3)))
        self.assertAlmostEqual(union_volume([cube, cube+[.5, .5, 0]]), 1.75, places=9)
        separated = [cube, cube+[2., 0., 0.]]
        self.assertAlmostEqual(union_volume(separated), 2., places=9)
        self.assertFalse(any(contains([[1.5, .5, .5]], halfspaces(p))[0] for p in separated))

    def test_nearly_coplanar_geometry_preserves_small_volume(self):
        # 极薄但仍三维的倾斜盒；不能靠随机抖动或丢弃体积通过。
        basis, _ = np.linalg.qr(np.array([[1., 2., 3.], [3., 1., 4.], [2., 5., 1.]]))
        points = np.array(list(product((0., 1.), repeat=3)))*[1., .5, 2e-9]
        points = points@basis.T + [.3, .4, .5]
        self.assertAlmostEqual(polytope_volume(points)/1e-9, 1., places=6)
        self.assertTrue(contains(points, halfspaces(points)).all())
        clipped = clip_polytope(points, .5+np.dot([.3, .4, .5], basis[:, 0]), -basis[:, 0])
        self.assertAlmostEqual(polytope_volume(clipped)/5e-10, 1., places=6)

    def test_four_bus_thin_overlap_regression(self):
        # 四节点 40000 元 LP 外包络相减产生的真实退化输入。
        points = np.array([
            [1., 0., 2.2616019600602275e-8],
            [1., 1.5831213716380983e-8, -2.8195154956956677e-24],
            [.6162329432184923, 0., .5482386751547624],
            [.8493614574978721, .15063855833412398, 0.],
            [.4655943848852564, .15063855833574227, .5482386751522994],
            [1., 0., 2.3405555953412012e-8],
            [.8493614582499969, .15063855813510713, 0.],
            [.46559438614650345, .15063855813639246, .5482386744267795],
            [1., 1.638388916320674e-8, -2.917946166542045e-24],
            [.6162329442795602, 0., .5482386744298092]])
        self.assertGreater(polytope_volume(points), 0.)
        self.assertLess(polytope_volume(points), 1e-8)
        self.assertTrue(contains(points, halfspaces(points)).all())

    def test_mp1_total_allows_load_redistribution(self):
        equations = PlanningEquations(self.network, 'linear')
        problem = PlanningModel(equations, min_total=3800., threads=1)
        with problem.model:
            answer = problem.solve()
        self.assertTrue(answer['feasible'])
        self.assertEqual(answer['objective'], 1.)
        self.assertGreaterEqual(sum(answer['p']), 3800.-1e-6)
        direct = PlanningModel(equations, min_total=3800., threads=1)
        with direct.model:
            reference = direct.solve()
        self.assertEqual(answer['objective'], reference['objective'])

    def test_continuous_domain_covers_other_plans_and_ignores_grid(self):
        result = build_continuous_region(self.network, 'linear', 1., self.bounds, threads=1)
        self.assertEqual(result['status'], 'certified')
        self.assertGreater(len(result['inner']), 1)
        self.assertLessEqual(result['coverage_bound'], 1e-8)
        before = json.dumps(json_value(result['inner']), sort_keys=True)
        for n in (3, 11):
            grid = (np.indices((n,)*3).reshape(3, -1).T+.5)*self.bounds/n
            sample_region(result, grid, self.bounds)
        self.assertEqual(before, json.dumps(json_value(result['inner']), sort_keys=True))
        points = []
        for choice in product((0, 1), repeat=4):
            plan = upgrade_plan(self.network, choice)
            net = self.network.tree(self.network.encode_plan(plan))
            if net.cost <= 1:
                for normal in (np.eye(3).tolist()+[[1., 1., 1.]]):
                    points.append(dispatch_support(net, 'linear', normal)['p'])
        self.assertTrue(np.all(sample_region(result, points, self.bounds) != -1))

    def test_unfinished_global_search_does_not_certify_domain(self):
        with patch.object(RemainingRegionModel, 'solve', return_value=dict(complete=False, bound=1., x=None, p=None)):
            result = build_continuous_region(self.network, 'linear', 0., self.bounds, threads=1)
        self.assertEqual(result['status'], 'unknown')
        self.assertGreater(region_metrics(result, self.bounds)['outer_volume'], region_metrics(result, self.bounds)['inner_volume'])

    def test_complete_replay_round_trip_beyond_old_limit(self):
        with TemporaryDirectory() as directory:
            with RunMonitor(record=True, output=directory, stream=StringIO()) as monitor:
                geometry = np.array([[0., 1., 2.]])
                monitor('phase_start', method='linear', phase='linear', geometry=geometry)
                geometry[:] = 99
                for i in range(701):
                    monitor('point', point=[float(i), 0., 0.], query=i+1)
                monitor('phase_start', method='hybrid', phase='socp')
                monitor('completed', message='complete')
                monitor.save_snapshot()
                state = {}
                for frame in monitor.history[:352]:
                    state.update(frame['patch'])
                self.assertEqual(state['point'], [350., 0., 0.])
                self.assertEqual(monitor.history[0]['patch']['geometry'], [[0., 1., 2.]])
                html = (Path(directory)/'live_view.html').read_text(encoding='utf-8')
                encoded = re.search(r'window.SAVED_REPLAY_GZIP="([^"]+)";', html)[1]
                saved = json.loads(gzip.decompress(base64.b64decode(encoded)))
                self.assertEqual(len(saved['history']), 704)
                self.assertEqual(saved['history'], monitor.history)
                self.assertEqual(saved['geometry'], [])
                self.assertEqual(len((Path(directory)/'events.jsonl').read_text(encoding='utf-8').splitlines()), 704)


if __name__ == '__main__':
    unittest.main()
