"""并集去重、内部空隙、证书归属及预算继承的回归。"""
from itertools import product
from unittest.mock import patch
import unittest

import numpy as np
from threadpoolctl import threadpool_limits

from main import Case33, FourBus, build_continuous_region
from model import PlanningEquations, PlanningModel, PlanningSP, RemainingRegionModel, PLANNING_TOL, evaluation_bounds
from plot import region_geometry
from region import RegionState, contains, halfspaces
from tests.reference import fixed_topology
from tests.planning_checks import margin


CUBE = np.array(list(product((0., 1.), repeat=3)))


class UnionGeometryTests(unittest.TestCase):
    def test_other_owner_skips_sp_without_claiming_current_scheme(self):
        network = Case33(upgrade_count=4)
        state = RegionState([400., 4670., 630.], network.power_limit, 0.)
        initial = network.initial_plan
        a, b = [network.encode_plan(plan).astype(int) for plan in
                (initial, initial | {'27-28': 'parallel'})]
        state.add_scheme(a, initial, 0.)
        # 合成几何夹具只检验归属规则，不用作物理模型证书。
        state.add_point(a, CUBE)
        state.add_scheme(b, network.decode_plan(b), network.cost@b)
        self.assertIsNone(state.next_point(b))
        self.assertEqual(len(state.records[tuple(b)]['inner']), 0)

    def test_vertices_covered_does_not_hide_interior_gap(self):
        state = RegionState(np.ones(3), 3., 0.)
        for x in ([0], [1], [2]):
            state.add_scheme(x, {'line': str(x[0])}, 0.)
        state.add_point([0], CUBE*[.4, 1., 1.])
        state.add_point([1], CUBE*[.4, 1., 1.]+[.6, 0., 0.])
        self.assertIsNone(state.next_point([2]))
        self.assertIsNone(state.covering_schemes([[.5, .5, .5]])[0])
        e = PlanningEquations(Case33(upgrade_count=4), 'linear')
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
        cls.network = fixed_topology(Case33(upgrade_count=4))
        cls.bounds = np.array([400., 4670., 630.])

    @classmethod
    def tearDownClass(cls):
        cls.threads.restore_original_limits()

    def test_certified_witness_still_expands_same_scheme_support(self):
        net = FourBus()
        state = RegionState(np.ones(3), 3., 0.)
        choices = [net.initial_plan | {'01': kind} for kind in ('L', 'M', 'H')]
        left, right, filling = [net.encode_plan(choice).astype(int) for choice in choices]
        # 两个已知凸域覆盖外域所有顶点，但中间仍有空隙；几何夹具不充当物理证书。
        for x, choice, poly in ((left, choices[0], CUBE*[.4, 1., 1.]),
                                (right, choices[1], CUBE*[.4, 1., 1.]+[.6, 0., 0.])):
            state.add_scheme(x, choice, net.cost@x)
            state.add_point(x, poly)
        witness = np.array([.5, .25, .75])
        self.assertIsNone(state.covering_schemes([witness])[0])
        answer = dict(x=left, p=np.array([.4, 1., 1.]), feasible=True, bound=3., status='optimal')
        residuals = [dict(complete=False, bound=1., x=filling, p=witness, feasible=True),
                     dict(complete=True, bound=0., x=None, p=None)]
        with patch('main.RegionState', return_value=state), \
             patch.object(PlanningModel, 'solve', return_value=answer), \
             patch.object(RemainingRegionModel, 'solve', side_effect=residuals), \
             patch.object(PlanningSP, 'solve', return_value=dict(feasible=True)) as check:
            result = build_continuous_region(net, 'linear', np.inf, np.ones(3), threads=1)
        self.assertEqual(result['status'], 'certified')
        self.assertTrue(any(not np.allclose(call.args[1], witness) for call in check.call_args_list))
        self.assertTrue(contains([witness], state.inner_equations(filling))[0])

    def test_fourbus_physical_search_completes_coverage(self):
        net = FourBus()
        bounds = evaluation_bounds(net, threads=1)
        result = build_continuous_region(net, 'linear', np.inf, bounds,
                                         residual_mode='physical', time_limit=20., threads=1)
        self.assertEqual(result['status'], 'certified')
        self.assertLessEqual(result['coverage_bound'], 1e-8)

    def test_ordinary_vertex_queries_are_outside_existing_union(self):
        next_point = RegionState.next_point
        ordinary = []
        def candidate(state, x):
            point = next_point(state, x)
            if point is not None:
                self.assertIsNone(state.covering_schemes([point])[0])
                ordinary.append(point)
            return point
        with patch.object(RegionState, 'next_point', new=candidate):
            result = build_continuous_region(self.network, 'linear', 1., self.bounds, threads=1)
        self.assertEqual(result['status'], 'certified')
        self.assertTrue(ordinary)
        self.assertLessEqual(result['coverage_bound'], 1e-8)

    def test_maximum_boundary_stagnation_still_requires_physical_certificate(self):
        equations = PlanningEquations(self.network, 'socp')
        x = equations.network.encode_plan(equations.network.initial_plan | {'2-3': 'parallel'})
        point = [154.27833628730312, 3500.5974303868097, 240.65872464525015]
        answer = PlanningSP(equations, threads=1).solve(x, point)
        self.assertEqual(set(answer), {'feasible', 'state', 'cut'})
        self.assertTrue(answer['feasible'])
        self.assertGreaterEqual(margin(equations, x, point, answer['state']), -PLANNING_TOL)



if __name__ == '__main__':
    unittest.main()
