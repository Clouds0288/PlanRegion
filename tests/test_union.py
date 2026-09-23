"""并集去重、内部空隙、证书归属及预算继承的回归。"""
from itertools import product
from unittest.mock import patch
import unittest

import numpy as np
from threadpoolctl import threadpool_limits

from main import Case33, build_continuous_region, solve_region
from region import ContinuousRegion
from model import planning_query
from model import PlanningEquations, PlanningSP, RemainingRegionModel, PLANNING_TOL
from plot import region_geometry
from region import RegionState, contains, halfspaces


CUBE = np.array(list(product((0., 1.), repeat=3)))


class UnionGeometryTests(unittest.TestCase):
    def test_other_owner_skips_sp_without_claiming_current_scheme(self):
        solver = ContinuousRegion(Case33(candidate_count=4), 'linear', 1., [400., 4670., 630.], planning_query, threads=1)
        initial = solver.equations.network.initial_plan
        a, b = [solver.equations.selection(plan).astype(int) for plan in
                (initial, initial | {'27-28': 'parallel'})]
        solver.region.add_scheme(a, initial, 0.)
        # 合成几何夹具只检验归属规则，不用作物理模型证书。
        solver.region.add_point(a, CUBE)
        with patch.object(solver.oracle, 'solve', side_effect=AssertionError('redundant SP')):
            self.assertTrue(solver.refine(b))
        self.assertEqual(solver.oracle.calls, 0)
        self.assertEqual(len(solver.region.records[tuple(b)]['inner']), 0)

    def test_vertices_covered_does_not_hide_interior_gap(self):
        state = RegionState(np.ones(3), 3., 0.)
        for x in ([0], [1], [2]):
            state.add_scheme(x, {'line': str(x[0])}, 0.)
        state.add_point([0], CUBE*[.4, 1., 1.])
        state.add_point([1], CUBE*[.4, 1., 1.]+[.6, 0., 0.])
        self.assertIsNone(state.next_point([2]))
        self.assertIsNone(state.covering_schemes([[.5, .5, .5]])[0])
        e = PlanningEquations(Case33(candidate_count=4), 'linear')
        problem = RemainingRegionModel(e, 0., np.ones(3), 3., [], state.inner_halfspaces(), 0., threads=1)
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
        state.add_scheme([0], {'line': 'type0'}, 0.)
        state.add_point([0], CUBE*.5)
        cache = {}
        old = region_geometry(state.records.values(), state.bounds, cache)
        before = state.progress
        self.assertIsNone(state.covering_schemes([[.8, .8, .8]])[0])
        state.add_point([0], CUBE)
        self.assertEqual(len(state.records[(0,)]['inner']), 8)
        self.assertNotEqual(before, state.progress)
        self.assertEqual(state.covering_schemes([[.8, .8, .8]])[0], (0,))
        self.assertNotEqual(old, region_geometry(state.records.values(), state.bounds, cache))


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
        solver = ContinuousRegion(self.network, 'linear', 1., self.bounds, planning_query, threads=1)
        next_point = solver.region.next_point
        ordinary = []
        def candidate(x):
            point = next_point(x)
            if point is not None:
                self.assertIsNone(solver.region.covering_schemes([point])[0])
                ordinary.append(point)
            return point
        with patch.object(solver.region, 'next_point', side_effect=candidate):
            result = solve_region(solver)
        self.assertEqual(result['status'], 'certified')
        self.assertTrue(ordinary)
        self.assertLessEqual(result['coverage_bound'], 1e-8)

    def test_maximum_boundary_stagnation_still_requires_physical_certificate(self):
        equations = PlanningEquations(self.network, 'socp')
        x = equations.selection(equations.network.initial_plan | {'2-3': 'parallel'})
        point = [154.27833628730312, 3500.5974303868097, 240.65872464525015]
        answer = PlanningSP(equations, threads=1).solve(x, point)
        self.assertEqual(set(answer), {'feasible', 'state', 'cut'})
        self.assertTrue(answer['feasible'])
        self.assertGreaterEqual(equations.margin(x, point, answer['state']), -PLANNING_TOL)



if __name__ == '__main__':
    unittest.main()
