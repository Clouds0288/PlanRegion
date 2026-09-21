"""并集去重、内部空隙、证书归属及预算继承的回归。"""
from itertools import product
from unittest.mock import patch
import unittest

import numpy as np
from threadpoolctl import threadpool_limits

from main import Case33, ContinuousRegion, build_continuous_region, joint_benders
from model import PlanningEquations, PlanningSP, RemainingRegionModel, PLANNING_TOL
from region import RegionState, contains, halfspaces


CUBE = np.array(list(product((0., 1.), repeat=3)))


class UnionGeometryTests(unittest.TestCase):
    def test_other_owner_skips_sp_without_claiming_current_scheme(self):
        solver = ContinuousRegion(Case33(candidate_count=4), 'linear', 1., [400., 4670., 630.], joint_benders)
        a, b = [solver.equations.selection(choice).astype(int) for choice in ([0, 0, 0, 0], [0, 0, 0, 1])]
        solver.region.add_scheme(a, [0, 0, 0, 0], 0.)
        # 合成几何夹具只检验归属规则，不用作物理模型证书。
        solver.region.add_point(a, CUBE)
        with patch.object(solver.oracle, 'solve', side_effect=AssertionError('redundant SP')):
            self.assertTrue(solver.refine(b))
        self.assertGreater(solver.counts['sp_skipped'], 0)
        self.assertEqual(len(solver.region.records[tuple(b)]['inner']), 0)

    def test_vertices_covered_does_not_hide_interior_gap(self):
        state = RegionState(np.ones(3), 3., 0.)
        for x in ([0], [1], [2]):
            state.add_scheme(x, x, 0.)
        state.add_point([0], CUBE*[.4, 1., 1.])
        state.add_point([1], CUBE*[.4, 1., 1.]+[.6, 0., 0.])
        self.assertIsNone(state.next_point([2]))
        self.assertIsNone(state.covering_schemes([[.5, .5, .5]])[0])
        e = PlanningEquations(Case33(candidate_count=4), 'linear')
        problem = RemainingRegionModel(e, 0., np.ones(3), 3., [], state.inner_halfspaces(), 0.)
        with problem.model:
            answer = problem.solve(1e-8)
        self.assertFalse(answer['complete'])
        self.assertGreater(answer['bound'], .09)
        self.assertIsNone(state.covering_schemes([answer['p']])[0])
        support = state.witness_support([2], answer['p'])
        self.assertLessEqual(len(support), 4)
        self.assertTrue(contains([answer['p']], halfspaces(support))[0])

    def test_growth_with_same_vertex_count_invalidates_caches(self):
        state = RegionState(np.ones(3), 3., 0.)
        state.add_scheme([0], [0], 0.)
        state.add_point([0], CUBE*.5)
        old = state.geometry()
        before = state.progress
        self.assertIsNone(state.covering_schemes([[.8, .8, .8]])[0])
        state.add_point([0], CUBE)
        self.assertEqual(len(state.records[(0,)]['inner']), 8)
        self.assertNotEqual(before, state.progress)
        self.assertEqual(state.covering_schemes([[.8, .8, .8]])[0], (0,))
        self.assertNotEqual(old, state.geometry())


class UnionSolverTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.threads = threadpool_limits(limits=1)
        cls.network = Case33(candidate_count=4)
        cls.bounds = np.array([400., 4670., 630.])

    @classmethod
    def tearDownClass(cls):
        cls.threads.restore_original_limits()

    def test_ordinary_vertex_queries_are_outside_existing_union(self):
        solver = ContinuousRegion(self.network, 'linear', 1., self.bounds, joint_benders)
        real_check = solver.check
        ordinary = []
        def check(x, point, *, origin='vertex'):
            if origin == 'vertex':
                self.assertIsNone(solver.region.covering_schemes([point])[0])
                ordinary.append(point)
            return real_check(x, point, origin=origin)
        with patch.object(solver, 'check', side_effect=check):
            result = solver.solve()
        self.assertEqual(result['status'], 'certified')
        self.assertTrue(ordinary)
        self.assertGreater(result['counts']['sp_skipped'], 0)
        self.assertLessEqual(result['coverage_bound'], 1e-8)

    def test_maximum_boundary_stagnation_still_requires_physical_certificate(self):
        equations = PlanningEquations(self.network, 'socp')
        x = equations.selection([1, 0, 0, 0])
        point = [154.27833628730312, 3500.5974303868097, 240.65872464525015]
        answer = PlanningSP(equations).solve(x, point)
        self.assertTrue(answer['feasible'])
        self.assertGreaterEqual(equations.margin(x, point, answer['state']), -PLANNING_TOL)

    def test_budget_reuse_keeps_model_and_budget_provenance(self):
        cache = {}
        first = build_continuous_region(self.network, 'linear', 0., self.bounds, reuse=cache)
        second = build_continuous_region(self.network, 'linear', 1., self.bounds, reuse=cache)
        self.assertEqual(second['status'], 'certified')
        self.assertGreater(second['counts']['reused_vertices'], 0)
        self.assertEqual(second['timing_mode'], 'incremental')
        self.assertTrue(all(c['cost'] <= 1. for c in second['certificates']))
        for method, budget, bounds, prior in [('socp', 1., self.bounds, first),
                                              ('linear', 0., self.bounds, second),
                                              ('linear', 1., self.bounds*2, first)]:
            with self.assertRaisesRegex(ValueError, '同一物理模型'):
                ContinuousRegion(self.network, method, budget, bounds, joint_benders, reuse=prior)


if __name__ == '__main__':
    unittest.main()
