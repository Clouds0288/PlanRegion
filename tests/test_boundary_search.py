"""Proof-condition tests for the experimental nonconvex box boundary search."""
import io
import unittest
from unittest.mock import patch

import numpy as np

from tests.benchmark_boundary_search import (BOUNDS, TAU, RayOracle, boundary_search,
    clip_boxes, coverage, inner_scales, maximal, json_value)
from region import GEOMETRY_TOL
from Network.case33bw import Case33


class BoxBoundaryTests(unittest.TestCase):
    def test_nonconvex_notch_is_not_filled(self):
        inner = np.array([[1., .3, 1.], [.3, 1., 1.]])
        gap, residual = coverage(np.ones((1, 3)), inner)
        self.assertAlmostEqual(gap[0], .7)
        self.assertGreater(residual[0], GEOMETRY_TOL)
        # The midpoint belongs to the convex hull, but not to either inner box.
        self.assertLess(inner_scales(np.array([[.65, .65, 1.]]), inner)[0], 1.)

    def test_axis_bound_applies_globally_and_keeps_boundary(self):
        outer = clip_boxes(np.ones((1, 3)), np.array([1., 0., 0.]), .4, np.ones(3))
        np.testing.assert_allclose(outer, [[.4+GEOMETRY_TOL, 1., 1.]], rtol=0, atol=1e-15)
        self.assertTrue(np.any(np.all(np.array([.4, 1., 1.]) <= outer, axis=1)))

    def test_plane_cut_excludes_zero_coordinate_from_disjunction(self):
        outer = clip_boxes(np.ones((1, 3)), np.array([1., 1., 0.]), .5, np.ones(3))
        self.assertEqual(len(outer), 2)
        self.assertFalse(np.any(np.all(np.array([.8, .8, 0.]) <= outer, axis=1)))
        self.assertTrue(np.any(np.all(np.array([.5, 1., 1.]) <= outer, axis=1)))

    def test_same_cut_is_applied_to_every_outer_box(self):
        outer = np.array([[1., .8, 1.], [.8, 1., 1.]])
        updated = clip_boxes(outer, np.ones(3), .4, np.ones(3))
        self.assertFalse(np.any(np.all(np.array([.6, .6, .6]) <= updated, axis=1)))
        for point in ([.3, .8, 1.], [.8, .3, 1.], [.8, 1., .3]):
            self.assertTrue(np.any(np.all(point <= updated, axis=1)))

    def test_dominance_does_not_use_tolerance_to_shrink_outer(self):
        points = np.array([[1., 1., 1.], [1.+1e-12, 1.-1e-12, 1.], [1., 1., 1.]])
        self.assertEqual(len(maximal(points)), 2)

    def test_new_children_go_after_existing_fifo_vertices(self):
        outer = np.array([[1., 1., .5], [.3, .8, 1.]])
        updated = clip_boxes(outer, np.array([1., 1., 0.]), .5, np.ones(3))
        np.testing.assert_array_equal(updated[0], outer[1])

    def test_all_vertices_required_and_zero_origin_needs_certificate(self):
        outer = np.array([[1., .2, .2], [.2, 1., 1.]])
        _, residual = coverage(outer, np.array([[1., .2, .2]]))
        self.assertLessEqual(residual[0], GEOMETRY_TOL)
        self.assertGreater(residual[1], GEOMETRY_TOL)
        _, residual = coverage(np.zeros((1, 3)), np.empty((0, 3)))
        self.assertTrue(np.isinf(residual[0]))

    def test_contraction_certificate_uses_all_coordinate_residuals(self):
        outer = np.array([[1., 0., .5], [.2, 1., 1.]])
        gap, residual = coverage(outer, (1-TAU)*outer)
        self.assertTrue(np.all(residual <= 0.))
        np.testing.assert_allclose(gap, TAU)

    def test_unknown_remains_null_in_json(self):
        self.assertEqual(json_value(dict(theta_lower=None, theta_upper=np.inf)),
                         dict(theta_lower=None, theta_upper=None))


class RayStatusTests(unittest.TestCase):
    def test_sp_unknown_never_produces_feasible_lower_bound(self):
        from time import perf_counter
        oracle = RayOracle(Case33(upgrade_count=8), 2., threads=1)
        with patch.object(oracle.sp, 'solve', return_value=dict(feasible=False, cut=None, state=None)):
            answer = oracle.query(BOUNDS, perf_counter()+3.)
        self.assertIsNone(answer['theta_lower'])
        self.assertEqual(oracle.certificates, [])
        self.assertTrue(np.isfinite(answer['theta_upper']))

    def test_known_feasible_lower_is_kept_without_a_new_mp_incumbent(self):
        from time import perf_counter
        oracle = RayOracle(Case33(upgrade_count=8), 2., threads=1)
        answer = oracle.query(BOUNDS, perf_counter()-1., certificates=[dict(p=.2*BOUNDS)])
        self.assertAlmostEqual(answer['theta_lower'], .2)
        self.assertGreater(answer['theta_upper'], .2)
        self.assertEqual(oracle.certificates, [])

    def test_missing_origin_certificate_cannot_certify_boundary(self):
        with patch('tests.benchmark_boundary_search.origin', return_value=None):
            result = boundary_search(Case33(upgrade_count=8), 2., 'C', 1., 1, io.StringIO())
        self.assertEqual(result['status'], 'unknown')
        self.assertEqual(result['certificates'], [])


if __name__ == '__main__':
    unittest.main()
