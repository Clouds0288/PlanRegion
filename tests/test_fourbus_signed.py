"""分区物理、有效割、坐标迁移与独立扫描的数值回归。"""
import unittest

from gurobipy import GRB
import numpy as np
from threadpoolctl import threadpool_limits

from Network.four_bus_five_corridor import FourBus
from model import GridPhysics, MasterProblem, SubProblem
from model import PortSubProblem, PortPhysics, port_bounds, ray_support
from monitor import signed_values, RunMonitor, _cut_polygon
from vertify import scan_line


class SignedPartitionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pool = threadpool_limits(limits=1)

    @classmethod
    def tearDownClass(cls):
        cls.pool.restore_original_limits()

    def test_positive_support_matches_original_model(self):
        for d in (2, 3):
            nodes = tuple(range(1, d+1))
            for direction in [*np.eye(d), np.ones(d)]:
                values = []
                for signed in (False, True):
                    net = FourBus(load_nodes=nodes)
                    equations = PortPhysics(net, np.ones(d)) if signed else GridPhysics(net, 'socp')
                    problem = MasterProblem(equations, budget=20000., direction=direction, threads=1)
                    with problem.model:
                        values.append(problem.solve()['objective'])
                self.assertAlmostEqual(*values, delta=2e-4)

    def test_signed_active_and_fixed_pf_reactive_balances(self):
        for sign in ([1, -1], [-1, 1], [-1, -1]):
            net = FourBus(load_nodes=(1, 2))
            equations = PortPhysics(net, sign)
            problem = MasterProblem(equations, power=[10., 10.], budget=20000., threads=1)
            with problem.model:
                answer = problem.solve()
                self.assertIsNotNone(answer)
                operation = problem.operation
                for flow, impedance, fixed, ratios in ((operation.P, equations.r, net.fixed_p, np.ones(2)),
                    (operation.Q, equations.reactance, net.fixed_q, net.q_ratio)):
                    for j, node in enumerate(net.load_nodes):
                        balance = sum(flow[e, k].X-impedance[e, k]*operation.ell[e, k].X
                            for e in equations.incoming[node] for k in equations.types[e])
                        balance -= sum(flow[e, k].X for e in equations.outgoing[node] for k in equations.types[e])
                        self.assertAlmostEqual(balance*net.base, fixed[net.node_index[node]]+sign[j]*10.*ratios[j], delta=2e-6)

    def test_cut_is_valid_for_all_topologies_in_its_partition(self):
        net = FourBus(load_nodes=(1, 2))
        equations = PortPhysics(net, [-1, 1])
        x = net.encode_plan(net.initial_plan)
        power = .9*port_bounds(net)
        answer = PortSubProblem(equations, threads=1).solve(x, power)
        self.assertFalse(answer['feasible'])
        cut = answer['cut']
        self.assertLess(cut[0]+cut[1:3]@power+cut[3:]@x, -1e-9)
        problem = MasterProblem(equations, threads=1)
        with problem.model as model:
            model.setObjective(cut[0]+cut[1:3]@problem.power+cut[3:]@problem.x, GRB.MINIMIZE)
            model.optimize()
            self.assertEqual(model.Status, GRB.OPTIMAL)
            self.assertGreaterEqual(model.ObjBound, -1e-7)

    def test_coordinate_and_cut_transform_preserves_values(self):
        sign = np.array([-1, 1, -1])
        cut = np.array([2., 3., -4., 5., -6., 7.])
        power, x = np.array([10., 20., 30.]), np.array([1., 0.])
        mapped = np.array(signed_values(cut.tolist(), sign, '-+-:', 'cut'))
        self.assertAlmostEqual(cut[0]+cut[1:4]@power+cut[4:]@x,
                               mapped[0]+mapped[1:4]@(sign*power)+mapped[4:]@x)

    def test_sp_objective_matches_original_solver(self):
        for sign in ([1, 1], [1, -1], [-1, 1]):
            net = FourBus(load_nodes=(1, 2))
            equations = PortPhysics(net, sign)
            x = net.encode_plan(net.initial_plan)
            values = [solver(equations, threads=1).solve(x, [80., 80.], score_only=True)['eta']
                      for solver in (SubProblem, PortSubProblem)]
            self.assertAlmostEqual(*values, delta=2e-7)

    def test_degenerate_three_dimensional_sp_points(self):
        x = np.array([0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0])
        cases = [([1, 1, 1], [54.893577353852464, 40.23711476094518, .0011271656155714817]),
                 ([1, 1, -1], [95.55900970476527, 0., 0.]),
                 ([1, 1, 1], [3.5127387337866605, 33.207723361860296, 8.489479055925155])]
        for sign, power in cases:
            equations = PortPhysics(FourBus(load_nodes=(1, 2, 3)), sign)
            answer = PortSubProblem(equations, threads=1).solve(x, power)
            self.assertFalse(answer['feasible'])
            self.assertLess(answer['cut'][0]+answer['cut'][1:4]@power+answer['cut'][4:]@x, -1e-9)

    def test_three_dimensional_rays_satisfy_original_constraints(self):
        for sign in ([1, 1, 1], [1, -1, 1], [-1, 1, -1], [-1, -1, -1]):
            net = FourBus(load_nodes=(1, 2, 3))
            equations = PortPhysics(net, sign)
            x = net.encode_plan(net.initial_plan)
            for direction in [*np.eye(3), np.ones(3)]:
                power = direction*port_bounds(net)
                answer = ray_support(equations, 20000., x, np.zeros(3), power, threads=1, time_limit=10.)
                np.testing.assert_allclose(answer['p'], answer['ray_fraction']*power, atol=2e-6)
                self.assertTrue(answer['feasible'])
        equations = PortPhysics(FourBus(load_nodes=(1, 2, 3)), [1, 1, 1])
        x = np.array([0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0])
        answer = ray_support(equations, 20000., x,
            np.array([40.56558247100506, 15.775955603107535, 11.992124112103376]),
            np.array([11.604715917462963, 33.208561013363024, .0010351113785352407]), threads=1, time_limit=10.)
        self.assertTrue(answer['feasible'])

    def test_scan_line_matches_independent_point_queries(self):
        lower, upper, n = np.array([-30., -30.]), np.array([30., 30.]), 4
        for index in ((0,), (3,)):
            _, states = scan_line(((1, 2), 1, lower, upper, n, index))
            for j, label in enumerate(states):
                power = lower+(np.array([*index, j])+.5)*(upper-lower)/n
                sign = np.where(power >= 0, 1, -1)
                equations = PortPhysics(FourBus(load_nodes=(1, 2)), sign)
                problem = MasterProblem(equations, power=power*sign, budget=20000., threads=1)
                with problem.model:
                    self.assertEqual(label, 1 if problem.solve() is not None else -1)

    def test_signed_metrics_and_three_dimensional_geometry(self):
        square = np.array([[-2., -2.], [0., -2.], [0., 0.], [-2., 0.]])
        reference = dict(axis_lower=[-2., -2.], bounds=[2., 2.], states=[[1, -1], [-1, -1]])
        result = dict(inner=[dict(vertices=square)], outer=[dict(vertices=square)])
        monitor = RunMonitor()
        monitor.validation(reference, result)
        metrics = monitor.validation_state['validation']['metrics']['inner']
        self.assertEqual(metrics['mr_percent'], 0.)
        self.assertEqual(metrics['fr_percent'], 0.)
        face = _cut_polygon([0., 1., 1., 1., 0.], [1.], np.ones(3), -np.ones(3))
        self.assertEqual(len(face), 6)
        self.assertLess(np.max(np.abs(face.sum(axis=1))), 1e-10)

    def test_scan_numerical_regression_points(self):
        cases = [([-66.284375, 80.325], 1),
                 ([-135.15, 27.73125, -94.0625], 1),
                 ([102.65, -90.61875, 10.7375], -1)]
        for power, expected in cases:
            power = np.asarray(power)
            nodes = tuple(range(1, len(power)+1))
            _, states = scan_line((nodes, 1, power-.5, power+.5, 1, (0,)*(len(power)-1)))
            self.assertEqual(states[0], expected)


if __name__ == '__main__':
    unittest.main()
