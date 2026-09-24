"""Geometric certificates, especially nonconvex gaps and coordinate faces."""
import unittest

import numpy as np
from shapely.geometry import Polygon, box
from shapely.ops import unary_union

from tests.case33_two_stage import (Interval, clip_polygon, interval_gap,
                                   maximal, region_from_intervals, scheme_facets, convex_distance,
                                   add_tree_relaxation)
from tests.audit_case33_two_stage import directed_distance


class TwoStageGeometryTests(unittest.TestCase):
    def test_support_clips_only_the_proved_halfspace(self):
        square = np.array([[0., 0.], [10., 0.], [10., 10.], [0., 10.]])
        polygon = Polygon(clip_polygon(square, [1., 1.], 12.))
        self.assertAlmostEqual(polygon.area, 68.)
        self.assertTrue(polygon.covers(box(0., 0., 5., 5.)))

    def test_nonconvex_notch_survives_convex_supports(self):
        inner = unary_union([box(0, 0, 2, 10), box(0, 0, 10, 2)])
        convex = inner.convex_hull
        self.assertTrue(convex.covers(box(5.9, 5.9, 6, 6)))
        self.assertFalse(inner.covers(box(5.9, 5.9, 6, 6)))
        outer = region_from_intervals(convex, [Interval(0, 2, 10), Interval(2, 10, 2)])
        self.assertAlmostEqual(outer.symmetric_difference(inner).area, 0.)

    def test_bound_for_entire_strip_not_only_polygon_vertices(self):
        points = np.array([[2., 10.], [10., 2.]])
        self.assertEqual(interval_gap(Interval(0., 10., 10.), [], points), 8.)
        self.assertEqual(interval_gap(Interval(0., 2., 10.), [], points), 0.)
        self.assertEqual(interval_gap(Interval(2., 10., 2.), [], points), 0.)

    def test_geometry_handles_axis_and_narrow_strip(self):
        points = np.array([[2., 0.]])
        self.assertAlmostEqual(interval_gap(Interval(2., 2.01, .02), [], points), .02)

    def test_maximal_does_not_merge_different_schemes(self):
        points = maximal([[2, 10], [10, 2], [1, 1], [2, 10]])
        self.assertEqual(len(points), 2)
        self.assertFalse(any(np.all(p >= [6, 6]) for p in points))

    def test_gap_upper_bounds_sampled_entire_node(self):
        points = np.array([[2., 9.], [6., 3.], [8., 1.]])
        node = Interval(3., 7., 7.)
        cuts = [dict(weights=[1., 1.], bound=10.)]
        bound = interval_gap(node, cuts, points)
        for first in np.linspace(3., 7., 23):
            for second in np.linspace(0., min(7., 10.-first), 21):
                distance = np.min(np.max(np.maximum([first, second]-points, 0.), axis=1))
                self.assertLessEqual(distance, bound+1e-10)

    def test_full_polygon_distance_finds_a_mid_edge_nonconvex_gap(self):
        inner = unary_union([box(0, 0, 2, 10), box(0, 0, 10, 2)])
        # Every convex-hull vertex is in the union, but (6,6) is not.
        self.assertAlmostEqual(directed_distance(inner.convex_hull, [(2, 10), (10, 2)]), 4., places=5)

    def test_single_scheme_convex_cover_has_correct_min_max_order(self):
        certificates = [dict(x=[0], p=[2., 10., 0.]), dict(x=[1], p=[10., 2., 0.])]
        facets = scheme_facets(certificates)
        outer = np.array([[0., 0.], [10., 0.], [10., 2.], [2., 10.], [0., 10.]])
        gap = interval_gap(Interval(0., 10., 10.), [], np.array([[2., 10.], [10., 2.]]), facets, outer)
        self.assertGreaterEqual(gap, 4.)
        same = scheme_facets([dict(x=[0], p=[2., 10., 0.]), dict(x=[0], p=[10., 2., 0.])])
        self.assertLessEqual(interval_gap(Interval(0., 10., 10.), [], np.array([[2., 10.], [10., 2.]]), same, outer), 2e-5)

    def test_convex_distance_clips_at_coordinate_faces(self):
        facets = np.array([[1., 1., -2.]])
        np.testing.assert_allclose(convex_distance([[0., 10.], [10., 0.], [4., 4.]], facets), [8., 8., 3.])

    def test_tree_extension_preserves_a_fixed_case33_plan(self):
        from Network.case33bw import Case33
        from model import PlanningEquations, PlanningModel
        from vertify import ACPowerFlow
        network = Case33(upgrade_count=8)
        choice = {c.id: c.existing_type if c.initial_active else None for c in network.corridors}
        for strengthened in (False, True):
            problem = PlanningModel(PlanningEquations(network, 'socp'), power=np.zeros(3),
                                    fixed_plan=choice, budget=2., threads=1)
            if strengthened:
                add_tree_relaxation(problem)
            answer = problem.solve(time_limit=10.)
            self.assertTrue(answer['feasible'])
            self.assertEqual(ACPowerFlow(network.tree(answer['x'])).classify(answer['p'])[0], 1)
            problem.model.dispose()


if __name__ == '__main__':
    unittest.main()
