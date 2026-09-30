"""Case33 背景负荷、反送、联合割与独立扫描回归。"""
import unittest

from gurobipy import GRB
import numpy as np
from threadpoolctl import threadpool_limits

from Network.case33bw import Case33
from model import MasterProblem
from model import PortPhysics, PortSubProblem, port_bounds, ray_support
from vertify import ac_scan_line


class Case33SignedTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pool = threadpool_limits(limits=1)

    @classmethod
    def tearDownClass(cls):
        cls.pool.restore_original_limits()

    def test_background_balances_and_source_export(self):
        net = Case33(load_nodes=(18, 25))
        equations = PortPhysics(net, [-1, -1])
        problem = MasterProblem(equations, power=[1000., 4000.], budget=7, threads=1)
        with problem.model:
            self.assertIsNotNone(problem.solve())
            operation = problem.operation
            self.assertAlmostEqual(net.fixed_p.sum(), 3205.)
            self.assertAlmostEqual(net.fixed_q.sum(), 2060.)
            source_p = sum(operation.P[key].X for key in equations.keys
                           if equations.ends[key[0]][0] == net.root)
            source_q = sum(operation.Q[key].X for key in equations.keys
                           if equations.ends[key[0]][0] == net.root)
            loss_p = sum(equations.r[key]*operation.ell[key].X for key in equations.keys)
            loss_q = sum(equations.reactance[key]*operation.ell[key].X for key in equations.keys)
            self.assertLess(source_p, 0.)
            self.assertAlmostEqual((source_p-loss_p)*net.base, 3205.-5000., delta=1e-4)
            self.assertAlmostEqual((source_q-loss_q)*net.base, 2060., delta=1e-4)
            self.assertTrue(np.isfinite(equations.y_ub_global).all())

    def test_background_origin_and_rays_in_two_and_three_dimensions(self):
        for nodes in ((18, 25), (18, 25, 30)):
            net = Case33(load_nodes=nodes)
            equations = PortPhysics(net, [1, -1, -1][:len(nodes)])
            x = net.encode_plan(net.initial_plan)
            answer = PortSubProblem(equations, threads=1).solve(x, np.zeros(len(nodes)), score_only=True)
            self.assertTrue(answer['feasible'])
            for direction in (np.ones(len(nodes)), np.eye(len(nodes))[0]):
                answer = ray_support(equations, 7 if isinstance(net, Case33) else 20000., x, np.zeros(len(nodes)), direction*port_bounds(net), threads=1, time_limit=10.)
                self.assertTrue(answer['feasible'])
                self.assertGreater(answer['ray_fraction'], 0.)

    def test_cut_excludes_point_and_preserves_all_feasible_topologies(self):
        net = Case33(load_nodes=(18, 25))
        equations = PortPhysics(net, [-1, 1])
        x = net.encode_plan(net.initial_plan)
        power = np.array([8000., 12000.])
        answer = PortSubProblem(equations, threads=1).solve(x, power)
        self.assertFalse(answer['feasible'])
        cut = answer['cut']
        self.assertLess(cut[0]+cut[1:3]@power+cut[3:]@x, -1e-7)
        problem = MasterProblem(equations, budget=7, threads=1)
        with problem.model as model:
            model.setObjective(cut[0]+cut[1:3]@problem.power+cut[3:]@problem.x, GRB.MINIMIZE)
            model.optimize()
            self.assertEqual(model.Status, GRB.OPTIMAL)
            self.assertGreaterEqual(model.ObjBound, -1e-7)

    def test_independent_scan_classifies_import_export_and_infeasible_point(self):
        for power, expected in (([100., 100.], 1), ([-1000., -4000.], 1), ([15000., 15000.], -1)):
            power = np.array(power)
            answer = ac_scan_line((np.array([0]), power[None]), network=Case33(), budget=7, schemes=None)
            self.assertEqual(answer['states'][0], expected)

    def test_small_current_boundary_certificate(self):
        net = Case33(load_nodes=(18, 25, 30))
        equations = PortPhysics(net, [1, 1, 1])
        x = np.array([1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1,
                      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0])
        answer = PortSubProblem(equations, threads=1).solve(x,
            [582.1563873981617, .0025, .0008559204862438141], score_only=True)
        self.assertTrue(answer['feasible'])



if __name__ == '__main__':
    unittest.main()
