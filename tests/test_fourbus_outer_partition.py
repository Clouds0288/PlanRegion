"""Independent numerical checks for the no-continuous-product corner oracle."""
import itertools
import unittest

import numpy as np
from gurobipy import GRB

from experiments.fourbus_outer_partition import CornerViolation, PolyhedralCornerViolation
from Network.four_bus_five_corridor import FourBus
from model import GridPhysics, PLANNING_TOL
from tests.test_fourbus_outer import primal_violation


class CornerViolationTests(unittest.TestCase):
    def setUp(self):
        self.equations = GridPhysics(FourBus(), 'socp')
        self.search = CornerViolation(self.equations, 20000., np.array([96., 41., 35.]), 96., threads=1)

    def tearDown(self):
        self.search.model.dispose()
        self.search.witness_problem.model.dispose()

    def test_corner_maximum_matches_eight_independent_primal_sps(self):
        x = self.equations.network.encode_plan(self.equations.network.initial_plan)
        self.search.problem.x.LB = self.search.problem.x.UB = x
        lower, upper = np.array([15., 6., 5.]), np.array([35., 18., 15.])
        answer = self.search.solve_box(lower, upper, 1., 30., 1e-8)
        reference = [primal_violation(self.equations, x, np.where(bits, upper, lower))
                     for bits in itertools.product((0, 1), repeat=3)]
        self.assertEqual(answer['status'], GRB.OPTIMAL)
        self.assertAlmostEqual(self.search.model.ObjVal, max(reference), delta=2e-7)
        # The residual-corrected certificate is a lower bound, not an exact value.
        self.assertLessEqual(answer['objective'], max(reference)+2e-7)
        self.assertGreaterEqual(answer['bound'], max(reference)-2e-7)
        self.assertLessEqual(answer['bound'], max(reference)+2e-7)
        center = primal_violation(self.equations, x, (lower+upper)/2.)
        self.assertGreaterEqual(answer['bound'], center-2e-7)

    def test_corner_outside_outer_is_only_an_upper_bound(self):
        x = self.equations.network.encode_plan(self.equations.network.initial_plan)
        self.search.problem.x.LB = self.search.problem.x.UB = x
        # All actual outer points have tiny loads, while box corners can be large.
        cut = np.r_[1., -np.ones(3), np.zeros(len(x))]
        self.search.add_cut(cut)
        lower, upper = np.zeros(3), np.array([40., 20., 20.])
        answer = self.search.solve_box(lower, upper, 1., 30., 1e-8)
        actual = self.search.witness(x, lower, upper, 10.)
        self.assertGreater(np.sum(answer['corner']), 1.)
        self.assertLessEqual(np.sum(actual['p']), 1.+PLANNING_TOL)
        eta = primal_violation(self.equations, x, actual['p'])
        self.assertLessEqual(actual['objective'], eta+2e-7)
        self.assertGreaterEqual(answer['bound'], eta-2e-7)

    def test_local_bound_is_not_carried_to_sibling(self):
        x = self.equations.network.encode_plan(self.equations.network.initial_plan)
        self.search.problem.x.LB = self.search.problem.x.UB = x
        first = self.search.solve_box(np.zeros(3), np.ones(3), 1., 30., 1e-8)
        second = self.search.solve_box(np.array([30., 10., 10.]), np.array([40., 20., 20.]),
                                       1., 30., 1e-8)
        self.assertGreater(second['bound'], first['bound']+1e-3)

    def test_no_continuous_bilinear_load_constraint_remains(self):
        model = self.search.model
        model.update()
        self.assertEqual(model.Params.NonConvex, 0)
        self.assertNotIn('dual_objective', [row.QCName for row in model.getQConstrs()])
        self.assertEqual(model.NumQNZs, 0)
        self.assertEqual(model.NumQConstrs, len(self.search.template.cone_slices))
        for row in model.getQConstrs():
            expression = model.getQCRow(row)
            self.assertTrue(all(expression.getVar1(i).sameAs(expression.getVar2(i))
                                for i in range(expression.size())))


class PolyhedralCornerTests(CornerViolationTests):
    def setUp(self):
        self.equations = GridPhysics(FourBus(), 'socp')
        self.search = PolyhedralCornerViolation(self.equations, 20000.,
                                                np.array([96., 41., 35.]), 96., threads=1)

    def test_no_continuous_bilinear_load_constraint_remains(self):
        self.search.model.update()
        self.assertEqual(self.search.model.NumQNZs, 0)
        self.assertEqual(self.search.model.NumQConstrs, 0)


if __name__ == '__main__':
    unittest.main()
