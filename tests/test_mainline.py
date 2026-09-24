"""新主线数值回归：全局界、非凸缺口、全节点覆盖及 SOCP 物理证书。"""
from itertools import product
import unittest

import gurobipy as gp
from gurobipy import GRB
import numpy as np

from Network import Network, Corridor, TypeParameters
from Network.case33bw import Case33
from main import build_continuous_region
from model import PlanningEquations, PlanningModel, PLANNING_TOL
from region import (Cell, apply_disjunction, cell_gap, clip_polytope, contains,
                    cover_partition, downward_facets, halfspaces, hull_distance,
                    polytope_volume, select_cell, subtract_disjunction)
from tests.planning_checks import margin


def two_area_network():
    """两片区各 200/800 kW，可升级一处；近零线路损耗保留完整 SOCP。"""
    types = (TypeParameters('existing', 1e-6, 0., .2, 0.),
             TypeParameters('upgrade', 1e-6, 0., .8, 1.))
    return Network('two_area', 0, (1, 2),
        tuple(Corridor(str(i), (0, i), 'existing', True, types) for i in (1, 2)),
        1000., 10., np.zeros(2), np.zeros(2), (1, 2), [0., 0.], .9**2, 1.1**2,
        1000., source_pmax=1., source_qmax=1., source_smax=1.)


class MainlineModelTests(unittest.TestCase):
    def test_empty_and_timeout_have_distinct_conclusions(self):
        net = two_area_network()
        empty = build_continuous_region(net, -1., epsilon_kw=1., time_limit=5., threads=1)
        self.assertEqual(empty['status'], 'empty')
        self.assertEqual(empty['max_gap_kw'], 0.)
        self.assertEqual(empty['outer'], [])
        timeout = build_continuous_region(net, 1., epsilon_kw=1., time_limit=0., threads=1)
        self.assertEqual(timeout['status'], 'time_limit')
        self.assertGreater(timeout['max_gap_kw'], 1.)
        self.assertTrue(timeout['outer'])

    def test_global_distance_detects_integer_gap_and_query_switch_is_clean(self):
        equations = PlanningEquations(two_area_network(), 'socp')
        problem = PlanningModel(equations, budget=1., strengthen=True, threads=1)
        with problem.model:
            right = problem.solve(weights=[1., 0.])
            left = problem.solve(weights=[0., 1.], incumbent=right)
            middle = (right['p']+left['p'])/2
            distance = problem.solve(target=np.array([500., 500.]), gap_kw=.001)
            self.assertGreater(distance['bound'], 299.99)
            self.assertLess(distance['objective'], 300.01)
            self.assertTrue(distance['feasible'])
            self.assertGreater(middle.min(), 399.)  # 两个方案的凸组合不能当作内域。
            projected = problem.solve(target=np.array([500., 500.]), cut_normals=np.eye(2))
            self.assertAlmostEqual(projected['objective'], distance['objective'], delta=.002)
            # 跨目标切换必须取消上一轮退让约束/变量。
            again = problem.solve(weights=[1., 0.])
            self.assertAlmostEqual(right['objective'], again['objective'], delta=.002)
            self.assertEqual(problem.distance_kw.UB, 0.)
            self.assertEqual(problem.projected_deficit_kw.UB, 0.)
            self.assertGreaterEqual(again['bound'], again['objective']-1e-6)
            self.assertGreaterEqual(margin(equations, again['x'], again['p'], again['state']), -PLANNING_TOL)

    def test_warm_start_does_not_impose_the_old_cost_ceiling(self):
        equations = PlanningEquations(two_area_network(), 'socp')
        fixed = PlanningModel(equations, power=[0., 0.], budget=1., threads=1)
        with fixed.model:
            cheap = fixed.solve()
        problem = PlanningModel(equations, budget=1., threads=1)
        with problem.model:
            answer = problem.solve(weights=[1., 0.], incumbent=cheap)
        self.assertAlmostEqual(cheap['objective'], 0.)
        self.assertGreater(answer['objective'], 799.9)
        self.assertEqual(equations.network.cost @ answer['x'], 1.)

    def test_two_dimensional_end_to_end_keeps_nonconvex_gap(self):
        result = build_continuous_region(two_area_network(), 1., epsilon_kw=1.,
            time_limit=20., query_time_limit=2., threads=1, max_supports=8)
        self.assertEqual(result['status'], 'certified')
        self.assertLessEqual(result['max_gap_kw'], 1.)
        self.assertTrue(contains([[500., 500.]], halfspaces(result['stage1_vertices']))[0])
        self.assertFalse(any(contains([[500., 500.]], halfspaces(c['vertices']))[0] for c in result['outer']))
        self.assertGreater(len(result['branch_cuts']), 0)
        # 对解析近似域的代表点逐个核对，没有误切两个不同建设方案。
        for p in ([799., 199.], [199., 799.], [0., 799.], [199., 199.]):
            self.assertTrue(any(contains([p], halfspaces(c['vertices']))[0] for c in result['outer']))

    def test_case33_fixed_plan_matches_promoted_experimental_objective(self):
        net = Case33(upgrade_count=8)
        answers = []
        for strengthen in (False, True):
            equations = PlanningEquations(net, 'socp')
            problem = PlanningModel(equations, budget=0., fixed_plan=net.initial_plan,
                                    strengthen=strengthen, threads=1)
            with problem.model:
                answer = problem.solve(target=np.array([900., 3000., 900.]), time_limit=10., gap_kw=.001)
            self.assertTrue(answer['feasible'])
            self.assertGreater(answer['bound'], 0.)
            self.assertLess(answer['objective']-answer['bound'], .002)
            self.assertGreaterEqual(margin(equations, answer['x'], answer['p'], answer['state']), -PLANNING_TOL)
            answers.append(answer['objective'])
        self.assertAlmostEqual(*answers, delta=.01)


class MainlineGeometryTests(unittest.TestCase):
    def test_distance_against_independent_gurobi_lp_in_two_and_three_dimensions(self):
        rng = np.random.default_rng(513)
        for d in (2, 3):
            facets = downward_facets(rng.uniform(0., 30., (9, d)))
            queries = rng.uniform(0., 50., (15, d))
            expected = []
            with gp.Model() as model:
                model.Params.OutputFlag = 0
                model.Params.Threads = 1
                p = model.addMVar(d)
                distance = model.addVar()
                model.addConstr(facets[:, :-1] @ p <= -facets[:, -1])
                rows = model.addConstr(p+distance >= np.zeros(d))
                model.setObjective(distance)
                for q in queries:
                    rows.RHS = q
                    model.optimize()
                    self.assertEqual(model.Status, GRB.OPTIMAL)
                    expected.append(model.ObjVal)
            np.testing.assert_allclose(hull_distance(queries, facets), expected, atol=1e-7)

    def test_nonconvex_coverage_and_global_cut_propagation(self):
        for d in (2, 3):
            masks = np.asarray(list(product((0., 1.), repeat=d)))
            one, two = np.ones(d), np.ones(d)
            one[0], two[1] = 10., 10.
            hulls = [((1,), downward_facets([one])), ((2,), downward_facets([two]))]
            gap, _, individual = cell_gap(np.array([one, two]), hulls)
            self.assertLess(individual.max(), .001)
            self.assertGreater(gap, 8.9)
            threshold = np.arange(3., 3.+d)
            cells = [Cell(masks*10.), Cell(masks*10.+2.)]
            cut = apply_disjunction(cells, threshold, np.eye(d))
            for cell in cut:
                self.assertTrue(np.all(np.any(cell.vertices <= threshold+1e-8, axis=1)))
            pieces = subtract_disjunction(masks*10., threshold, np.eye(d))
            self.assertAlmostEqual(sum(polytope_volume(v) for v in pieces), 10.**d-np.prod(10.-threshold), places=6)

    def test_coverage_split_conserves_area_and_selects_true_maximum(self):
        vertices = np.array([[0., 0.], [10., 0.], [0., 10.], [10., 10.]])
        facets = downward_facets([[3., 8.], [8., 3.]])
        children = cover_partition(vertices, facets, 2.)
        self.assertAlmostEqual(sum(polytope_volume(v) for v in children), 100., places=7)
        self.assertLessEqual(hull_distance(children[-1], facets).max(), 2.+1e-8)
        cells = [Cell(v) for v in children]
        chosen, _ = select_cell(cells, [((1,), facets)])
        exact = [cell_gap(c.vertices, [((1,), facets)])[0] for c in cells]
        self.assertAlmostEqual(exact[chosen], max(exact))


if __name__ == '__main__':
    unittest.main()
