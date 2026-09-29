"""独立双向净功率的数值回归；不改变生产默认模型。"""
import unittest
import numpy as np
from gurobipy import GRB
from threadpoolctl import threadpool_limits

from Network.four_bus_five_corridor import FourBus
from model import GridPhysics, MasterProblem, SubProblem, PLANNING_TOL
from experiments.signed_power_compare import (SignedPhysics, SignedMaster, SignedKkt,
    SignedRegion, Trace, initialize, expanded_faces)
from experiments.fourbus_outer_volume import budget_schemes
from region import contains


class SignedPowerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pool = threadpool_limits(limits=1)
        cls.net = FourBus(load_nodes=(1, 2))
        cls.equations = SignedPhysics(cls.net, np.full(2, -150.), np.full(2, 150.))
        cls.region = initialize(cls.equations, 20000., 4, Trace('test', 120.))

    @classmethod
    def tearDownClass(cls):
        cls.pool.restore_original_limits()

    def test_signed_load_and_original_defaults_are_separate(self):
        problem = SignedMaster(self.equations, budget=20000., power=np.array([-30., 20.]), threads=4)
        with problem.model:
            answer = problem.solve(time_limit=30.)
        self.assertIsNotNone(answer)
        np.testing.assert_allclose(answer['p'], [-30., 20.], atol=1e-8)
        p, q = self.net.loads([-30., 20.])
        self.assertLess(p[0, 0], 0.)
        self.assertLess(q[0, 0], 0.)
        original = MasterProblem(GridPhysics(self.net, 'socp'), cuts_only=True, threads=4)
        with original.model:
            original.model.update()
            np.testing.assert_array_equal(original.power.LB, 0.)

    def test_affine_geometry_and_negative_side_cut(self):
        region = SignedRegion(np.array([-20., -30.]), np.array([40., 50.]), [])
        x = np.array([1, 0])
        region.add_scheme(x, {}, 0.)
        region.apply_cut(np.array([10., 1., 0., 0., 0.]))  # p1 >= -10
        physical = region.physical(region.records[tuple(x)]['outer'])
        self.assertAlmostEqual(physical[:, 0].min(), -10.)
        points = np.array([[-10., -15.], [10., -15.], [10., 5.], [-10., 5.]])
        region.add_point(x, region.normalize(points))
        np.testing.assert_allclose(region.physical(region.normalize(points)), points)
        self.assertTrue(region.covered(region.normalize([[-5., -5.]]), .005)[0])
        self.assertFalse(region.covered(region.normalize([[-18., -5.]]), .005)[0])

    def test_kkt_matches_independent_all_scheme_vertex_sp(self):
        region = SignedRegion(self.region.axis_lower, self.region.axis_bounds, self.region.cuts)
        search = SignedKkt(self.equations, 20000., region, threads=4)
        with search.model:
            answer = search.solve(60., 1e-7)
        maximum = 0.
        for x in budget_schemes(self.equations, 20000., 4):
            problem = SignedMaster(self.equations, budget=20000., cuts_only=True,
                                  fixed_plan=self.net.decode_plan(x), threads=4)
            with problem.model:
                eta = problem.model.addVar(name='violation')
                self.equations.add_operation(problem.model, problem.choices, problem.loads, eta)
                problem.model.setObjective(eta, GRB.MINIMIZE)
                for power in region.physical(region.outer_for(x)):
                    problem.power.LB = problem.power.UB = power
                    problem.model.optimize()
                    self.assertEqual(problem.model.Status, GRB.OPTIMAL)
                    maximum = max(maximum, problem.model.ObjVal)
        self.assertAlmostEqual(answer['objective'], maximum, delta=2e-6)
        self.assertAlmostEqual(answer['bound'], maximum, delta=2e-6)

    def test_kkt_excludes_union_and_cut_preserves_all_signed_plans(self):
        search = SignedKkt(self.equations, 20000., self.region, threads=4)
        with search.model:
            answer = search.solve(60., 1e-7)
        point = self.region.normalize(answer['p'])
        for row in self.region.records.values():
            faces = expanded_faces(row['inner'], 1e-5)
            self.assertGreaterEqual(float(np.max(faces[:, :2]@point+faces[:, 2])), -1e-7)
        checked = SubProblem(self.equations, threads=4).solve(answer['x'], answer['p'], time_limit=30.)
        self.assertFalse(checked['feasible'])
        cut = checked['cut']
        self.assertLess(cut[0]+cut[1:3]@answer['p']+cut[3:]@answer['x'], -1e-9)
        problem = SignedMaster(self.equations, threads=4)  # 不限投资：检验切割对其他合法网架仍有效。
        with problem.model:
            problem.model.setObjective(cut[0]+cut[1:3]@problem.power+cut[3:]@problem.x, GRB.MINIMIZE)
            problem.model.optimize()
            self.assertEqual(problem.model.Status, GRB.OPTIMAL)
            self.assertGreaterEqual(problem.model.ObjBound, -1e-7)


if __name__ == '__main__':
    unittest.main()
