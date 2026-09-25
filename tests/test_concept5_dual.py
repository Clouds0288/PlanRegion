"""Dual rollback: two-dimensional certificates, road masks, and multi-type cuts."""
from dataclasses import asdict, replace
import json
from pathlib import Path
import unittest

import numpy as np
from gurobipy import GRB
from shapely.geometry import box

from Network import Network, TypeParameters
from Network.concept5 import Concept5, ROUTES
from main import build_continuous_region
from model import PlanningEquations, PlanningModel, PlanningSP
from region import clip_polytope, halfspaces, contains, polytope_volume
from survey import increase_bounds


class ConceptFiveDualTests(unittest.TestCase):
    def test_case_and_mask_preserve_indices_and_physics(self):
        whole, restricted = Concept5(), Concept5(['A', 'C', 'E'])
        self.assertEqual(whole.n+1, 5)
        self.assertEqual(whole.n_corridors, 10)
        self.assertEqual(whole.type_keys, restricted.type_keys)
        for name in ('r', 'reactance', 'capacity', 'cost', 'original_p', 'original_q'):
            np.testing.assert_array_equal(getattr(whole, name), getattr(restricted, name))
        self.assertEqual(restricted.road_allowed.tolist(), [True]*4+[True, False, True, False, True, False])
        self.assertNotEqual(whole.fingerprint, restricted.fingerprint)
        with self.assertRaises(ValueError):
            Concept5(['G'])
        initial = whole.initial_plan | {'31': None, 'A': 'new'}
        problem = PlanningModel(PlanningEquations(Concept5([]), 'socp'), fixed_plan=initial, budget=4., threads=1)
        with problem.model:
            self.assertIsNone(problem.solve())

    def test_two_dimensional_geometry_and_unknown_plan_ids(self):
        square = np.array([[0.,0.], [1.,0.], [1.,1.], [0.,1.]])
        triangle = clip_polytope(square, 1., np.array([-1.,-1.]))
        self.assertAlmostEqual(polytope_volume(triangle), .5)
        self.assertTrue(contains([[.2,.2]], halfspaces(triangle))[0])
        self.assertFalse(contains([[.8,.8]], halfspaces(triangle))[0])
        self.assertEqual(halfspaces(np.empty((0,2))).shape, (1,3))
        before = dict(inner=box(0,0,10,10), outer=box(0,0,10,10), plan_ids=None)
        after = dict(inner=box(0,0,20,10), outer=box(0,0,20,10), plan_ids=None)
        self.assertEqual(increase_bounds(before, after), (100.,100.,100.))

    def test_baseline_is_certified_without_enumeration(self):
        result = build_continuous_region(Concept5([]), 'socp', 4., np.array([135.,135.]),
                                         tau=1e-5, threads=1, time_limit=15.)
        self.assertEqual(result['status'], 'certified')
        self.assertGreater(result['counts']['cuts'], 0)
        self.assertEqual(len(result['inner']), 1)
        self.assertAlmostEqual(polytope_volume(result['inner'][0]['vertices']), 3023.2, delta=.15)

    def test_multi_type_cut_is_global_and_keeps_other_types(self):
        base = Concept5()
        corridors = list(base.corridors)
        index = next(i for i,c in enumerate(corridors) if c.id == 'A')
        large = corridors[index].types[0]
        small = TypeParameters('small', large.r*2., large.reactance, .4, 1.)
        corridors[index] = replace(corridors[index], types=(small, large))
        values = {name: getattr(base, name) for name in base.__dataclass_fields__}
        network = Network(**(values | {'corridors': tuple(corridors)}))
        equations = PlanningEquations(network, 'socp')
        oracle = PlanningSP(equations, threads=1)
        plan = network.initial_plan | {'31': None, 'A': 'small'}
        x, power = network.encode_plan(plan), np.array([80., 20.])
        answer = oracle.solve(x, power)
        self.assertFalse(answer['feasible'])
        self.assertIsNotNone(answer['cut'])
        cut = answer['cut']
        self.assertLess(cut[0]+cut[1:3]@power+cut[3:]@x, -1e-8)
        larger = network.encode_plan(plan | {'A': 'new'})
        self.assertTrue(oracle.solve(larger, power)['feasible'])
        problem = PlanningModel(equations, budget=4., threads=1)
        with problem.model:
            problem.model.setObjective(cut[0]+cut[1:3]@problem.power+cut[3:]@problem.x, GRB.MINIMIZE)
            problem.model.Params.TimeLimit = 30.
            problem.model.optimize()
            self.assertGreaterEqual(problem.model.ObjBound, -2e-7)


if __name__ == '__main__':
    unittest.main()
