"""当前 Case33：七开关、Hamming 预算、径向性。"""
import unittest

import numpy as np

from Network.case33bw import Case33
from vertify import budget_schemes
from model import GridPhysics, MasterProblem
from tests.planning_checks import fix_plan


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
            problem = MasterProblem(GridPhysics(net, [1, 1]), budget=budget, power=np.zeros(2), threads=1)
            fix_plan(problem, plan)
            with problem.model:
                answer = problem.solve()
            if budget == 1:
                self.assertIsNone(answer)
            else:
                self.assertEqual(answer['objective'], 2.)
                self.assertAlmostEqual(answer['bound'], 2.)

    def test_all_radial_schemes_obey_locked_branches_and_hamming_budget(self):
        net = Case33()
        initial = net.encode_plan(net.initial_plan)
        fixed = np.array([not c.switchable for c in net.corridors])
        for budget, count in ((0, 1), (1, 1), (2, 7), (4, 12), (7, 12)):
            schemes = budget_schemes(net, budget)
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
        problem = MasterProblem(GridPhysics(net, [1, 1]), threads=1)
        fix_plan(problem, plan)
        with problem.model:
            self.assertIsNone(problem.solve())


if __name__ == '__main__':
    unittest.main()
