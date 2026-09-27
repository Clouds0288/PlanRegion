"""Analytical union-volume regressions for the per-cut stopping experiment."""
import unittest

import numpy as np

from experiments.fourbus_outer_volume import ProjectedOuterVolume


class ProjectedOuterVolumeTests(unittest.TestCase):
    def test_overlap_is_counted_once_and_units_are_kw_cubed(self):
        volume = ProjectedOuterVolume(np.array([[0], [1]]), np.array([2., 3., 4.]), 9.)
        self.assertAlmostEqual(volume.outer_volume, 24., places=10)

    def test_zero_projection_change_can_precede_a_large_cut(self):
        volume = ProjectedOuterVolume(np.array([[0], [1]]), np.ones(3), 3.)
        first = volume.add_cut(np.array([1., -1., 0., 0., -.5]))
        self.assertAlmostEqual(first['volume_reduction_ratio'], 0., places=12)
        second = volume.add_cut(np.array([.5, 0., -1., 0., .5]))
        self.assertAlmostEqual(second['outer_volume'], .75, places=10)
        self.assertAlmostEqual(second['volume_reduction_ratio'], .25, places=10)
        self.assertLess(second['outer_volume'], .875)  # Convex hull would fill the missing corner.

    def test_ratio_uses_volume_before_the_cut(self):
        volume = ProjectedOuterVolume(np.array([[0]]), np.ones(3), 3.)
        volume.add_cut(np.array([.5, -1., 0., 0., 0.]))
        measured = volume.add_cut(np.array([.49975, -1., 0., 0., 0.]))
        self.assertAlmostEqual(measured['volume_reduction'], .00025, places=12)
        self.assertAlmostEqual(measured['volume_reduction_ratio'], .0005, places=12)


if __name__ == '__main__':
    unittest.main()
