"""当前 Case33：七开关、Hamming 预算、径向性与含背景负荷的 KKT 上界。"""
import unittest
from unittest.mock import patch

import numpy as np

from Network.case33bw import Case33
from experiments.fourbus_outer import GlobalViolation
from experiments.fourbus_outer_volume import budget_schemes
from model import GridPhysics, MasterProblem


class Case33SwitchingTests(unittest.TestCase):
    def test_only_requested_switches_and_original_types(self):
        net = Case33()
        self.assertEqual((net.root, net.load_nodes, net.n_types), (1, (18, 25), 37))
        self.assertEqual({frozenset(c.endpoints) for c in net.corridors if c.switchable},
                         {frozenset(edge) for edge in ((21, 8), (7, 8), (22, 12), (11, 12),
                                                       (9, 15), (33, 18), (25, 29))})
        self.assertTrue(all(tuple(t.id for t in c.types) == ('existing',) for c in net.corridors))
        self.assertEqual(net.switch_budget, 7)

    def test_budget_counts_both_opening_and_closing(self):
        net = Case33()
        plan = net.initial_plan | {'7-8': None, '21-8': 'existing'}
        x = net.encode_plan(plan)
        self.assertEqual(net.tree(x).cost, 2.)
        self.assertEqual(net.tree(net.encode_plan(net.initial_plan)).cost, 0.)
        for budget in (1, 2):
            problem = MasterProblem(GridPhysics(net, 'socp'), budget=budget,
                                    fixed_plan=plan, power=np.zeros(2), threads=1)
            with problem.model:
                answer = problem.solve()
            if budget == 1:
                self.assertIsNone(answer)
            else:
                self.assertTrue(answer['feasible'])
                self.assertEqual(answer['objective'], 2.)
                self.assertAlmostEqual(answer['bound'], 2.)

    def test_all_radial_schemes_obey_locked_branches_and_hamming_budget(self):
        net = Case33()
        initial = net.encode_plan(net.initial_plan)
        fixed = np.array([not c.switchable for c in net.corridors])
        for budget, count in ((0, 1), (1, 1), (2, 7), (4, 12), (7, 12)):
            schemes = budget_schemes(GridPhysics(net, 'socp'), budget, threads=1)
            self.assertEqual(len(schemes), count)
            for x in schemes:
                self.assertEqual(net.tree(x).n, 32)
                np.testing.assert_array_equal(x[fixed], initial[fixed])
                self.assertEqual(net.cost_offset+net.cost@x, np.count_nonzero(x != initial))
                self.assertLessEqual(net.tree(x).cost, budget)
                self.assertEqual(x[net.type_keys.index(('25-29', 'existing'))], 0)

    def test_fixed_plan_cannot_override_a_locked_line(self):
        net = Case33()
        plan = net.initial_plan | {'2-3': None, '21-8': 'existing'}
        problem = MasterProblem(GridPhysics(net, 'socp'), fixed_plan=plan,
                                cuts_only=True, threads=1)
        with problem.model:
            self.assertIsNone(problem.solve())

    def test_eta_upper_bound_includes_fixed_active_and_reactive_loads(self):
        net = Case33()
        search = GlobalViolation(GridPhysics(net, 'socp'), 7, [1., 1.], 2., threads=1)
        with search.model:
            search.model.update()
            expected = max(np.max(np.abs(v)) for v in net.loads([1., 1.]))
            self.assertEqual(search.violation.UB, expected)
            self.assertGreater(expected, 1./net.base)

    def test_zero_conservative_lower_bound_does_not_skip_candidate_sp(self):
        from threadpoolctl import threadpool_limits
        from experiments.case33_compare import build_kkt_region
        from experiments.fourbus_outer_kkt import UncoveredKktViolation
        from monitor import RunMonitor
        class CheckedCut(Exception):
            pass
        class Monitor(RunMonitor):
            def updated(self, region, event, **values):
                super().updated(region, event, **values)
                if event == 'cut':
                    raise CheckedCut()
        monitor = Monitor()
        with threadpool_limits(limits=1), \
             patch.object(UncoveredKktViolation, 'candidate_lower', return_value=0.), \
             self.assertRaises(CheckedCut):
            build_kkt_region(Case33(), budget=7, seconds=30., epsilon=.0001,
                             threads=4, monitor=monitor)
        self.assertEqual(monitor.state['sp'], 1)
        self.assertEqual(monitor.state['cuts'], 1)
        self.assertGreater(monitor.state['eta'], 0.)


if __name__ == '__main__':
    unittest.main()
