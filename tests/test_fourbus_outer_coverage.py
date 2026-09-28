"""独立验证认证并集排除、完整 KKT、换网架认证及剩余域证书。"""
import contextlib
import io
import unittest

import numpy as np
from gurobipy import GRB
from shapely.geometry import MultiPoint
from shapely.ops import unary_union
from threadpoolctl import threadpool_limits

from experiments.fourbus_outer_kkt import (COVERAGE_PAD, UncoveredKktViolation,
                                          certify_load, initialize, run_coverage)
from experiments.fourbus_outer_scan import candidate_eta
from experiments.fourbus_outer_volume import budget_schemes
from model import GridPhysics, MasterProblem, PLANNING_TOL
from Network.four_bus_five_corridor import FourBus
from region import clip_polytope, contains, halfspaces, initial_polytope


class CoverageKktTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pool = threadpool_limits(limits=1)
        cls.pool.__enter__()
        cls.equations = GridPhysics(FourBus(load_nodes=(1, 2)), 'socp')
        cls.initial = initialize(20000., 4)
        cls.bounds = np.asarray(cls.initial['axis_bounds'])
        cls.schemes = budget_schemes(cls.equations, 20000., 1)
        with contextlib.redirect_stdout(io.StringIO()):
            cls.result = run_coverage(cls.initial, epsilon=.01, seconds=100., threads=4)

    @classmethod
    def tearDownClass(cls):
        cls.pool.__exit__(None, None, None)

    def search(self, inner, cuts=()):
        return UncoveredKktViolation(self.equations, 20000., self.bounds,
                                     self.initial['total_bound'], inner, cuts, threads=4)

    def vertex_reference(self, inner, cuts):
        """GEOS 多边形相减独立于 G 的选择变量/KKT；每个顶点用原始 SP。"""
        offsets = COVERAGE_PAD*np.array([[-1, -1], [-1, 1], [1, -1], [1, 1]])
        expanded = [MultiPoint((np.asarray(p)[:, None, :]+offsets).reshape(-1, 2)).convex_hull
                    for p in inner]
        excluded = unary_union(expanded)
        base = initial_polytope(self.bounds, self.initial['total_bound'], self.bounds)

        def vertices(geometry):
            if geometry.is_empty:
                return []
            if geometry.geom_type == 'Polygon':
                return [p for ring in [geometry.exterior, *geometry.interiors] for p in ring.coords]
            if geometry.geom_type in ('LineString', 'Point'):
                return list(geometry.coords)
            return [p for part in geometry.geoms for p in vertices(part)]

        maximum = 0.
        for x in self.schemes:
            poly = base.copy()
            for values in cuts:
                cut = np.asarray(values)
                poly = clip_polytope(poly, cut[0]+cut[3:]@x, cut[1:3]*self.bounds)
            if not len(poly):
                continue
            remaining = MultiPoint(poly).convex_hull.difference(excluded)
            for point in vertices(remaining):
                maximum = max(maximum, candidate_eta(self.equations, x, np.asarray(point)*self.bounds))
        return maximum

    def test_changed_load_lp_matches_independent_remaining_vertices(self):
        inner = [np.asarray(row['vertices'])/self.bounds for row in self.result['inner']]
        cuts = self.result['cuts']
        reference = self.vertex_reference(inner, cuts)
        search = self.search(inner, cuts)
        with search.model:
            answer = search.solve(30., 0.)
            self.assertEqual(answer['status'], GRB.OPTIMAL)
            self.assertAlmostEqual(search.model.ObjVal, reference, delta=2e-6)
            self.assertGreaterEqual(answer['bound'], reference-2e-7)
            self.assertGreater(len(search.selector_terms), 0)
            self.assertEqual(search.model.NumQConstrs, len(search.template.cone_slices))
            self.assertEqual(search.model.NumQNZs, 0)
            self.assertTrue(all(var.UB < GRB.INFINITY for var in search.outer_dual))
        self.assertTrue(self.result['certified'])
        self.assertLessEqual(reference, self.result['bound']+2e-7)
        self.assertLessEqual(self.result['union_bound'], .01)

    def test_hole_is_not_lost_when_all_outer_vertices_are_covered(self):
        # 人工几何回归：原外域顶点全被覆盖，中央仍有缺口；方块不作物理证书。
        base = initial_polytope(self.bounds, self.initial['total_bound'], self.bounds)
        offsets = .12*np.array([[-1., -1.], [-1., 1.], [1., -1.], [1., 1.]])
        inner = [vertex+offsets for vertex in base]
        self.assertTrue(all(any(contains(vertex, halfspaces(poly))[0] for poly in inner)
                            for vertex in base))
        reference = self.vertex_reference(inner, [])
        search = self.search(inner)
        with search.model:
            answer = search.solve(30., 0.)
            self.assertEqual(answer['status'], GRB.OPTIMAL)
            self.assertAlmostEqual(search.model.ObjVal, reference, delta=2e-6)
            self.assertFalse(any(contains(answer['p']/self.bounds, halfspaces(poly))[0] for poly in inner))

    def test_fixed_load_can_switch_away_from_an_infeasible_scheme(self):
        power = np.array([90., 1.])
        weak = np.asarray(self.initial['initial'][1]['x'])
        self.assertGreater(candidate_eta(self.equations, weak, power), .01)
        problem = MasterProblem(self.equations, budget=20000., cuts_only=True, threads=4)
        with problem.model:
            eta = problem.model.addVar(name='violation')
            self.equations.add_operation(problem.model, problem.choices, problem.loads, eta)
            problem.model.setObjective(eta, GRB.MINIMIZE)
            answer = certify_load(problem, eta, power, 30.)
            self.assertTrue(answer['feasible'])
            self.assertFalse(np.array_equal(answer['x'], weak))
            self.assertLessEqual(candidate_eta(self.equations, answer['x'], power), PLANNING_TOL)
            answer = certify_load(problem, eta, np.array([50., 40.]), 30.)
            self.assertFalse(answer['feasible'])
            self.assertGreater(answer['bound'], .01)

    def test_certified_hulls_and_padding_have_independent_sp_certificates(self):
        for row in self.result['inner']:
            x = np.asarray(row['x'])
            for power in np.asarray(row['vertices']):
                self.assertLessEqual(candidate_eta(self.equations, x, power), PLANNING_TOL)
                nearby = power+COVERAGE_PAD*self.bounds
                self.assertLessEqual(candidate_eta(self.equations, x, nearby),
                                     self.result['padding_bound']+2e-8)

    def test_cuts_remain_valid_for_all_schemes(self):
        for values in self.result['cuts']:
            cut = np.asarray(values)
            direct = MasterProblem(self.equations, threads=1)
            with direct.model:
                direct.model.setObjective(cut[0]+cut[1:3]@direct.power+cut[3:]@direct.x, GRB.MINIMIZE)
                direct.model.optimize()
                self.assertEqual(direct.model.Status, GRB.OPTIMAL)
                self.assertGreaterEqual(direct.model.ObjBound, -1e-7)


if __name__ == '__main__':
    unittest.main()
