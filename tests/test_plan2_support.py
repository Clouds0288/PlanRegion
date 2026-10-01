"""方案二对照实验的纯几何部分：center 判据、外扩公式与有效性检查（不建求解器）。"""
import unittest
from unittest.mock import patch

import numpy as np

from experiments.plan2_geometry import (allowed_offsets, boundary_faces, chebyshev_center, choose_criterion,
                                        classify_center, cone_faces, expanded_faces, face_margins, faces_polytope,
                                        grid_metrics, interior_point, validity)
from experiments.test_support_face_certification_fourbus_2d import classify_support
from region import GEOMETRY_TOL, contains, halfspaces

TAU = .005


def same_faces(first, second):
    """两组单位法向面方程作为集合相同（与排列次序无关）。"""
    key = lambda faces: sorted(map(tuple, np.round(faces, 9)))
    return key(first) == key(second)


class CenterCriterionTests(unittest.TestCase):
    def setUp(self):
        # 不含原点的解析凸多边形：方形 [0.4, 0.8] x [0.3, 0.7]
        self.square = np.array(((.4, .3), (.8, .3), (.8, .7), (.4, .7)))
        self.faces = halfspaces(self.square)

    def test_support_equal_to_face_value_fails_origin_and_passes_center(self):
        with patch('gurobipy.Model', side_effect=AssertionError('geometry must not build a solver')):
            # 朝向原点的面 -xi_1 <= -0.4：支撑值恰好等于面值
            f = int(np.flatnonzero(np.all(np.isclose(self.faces[:, :-1], (-1., 0.)), axis=1))[0])
            normal, offset = self.faces[f, :-1], self.faces[f, -1]
            beta = -offset
            self.assertAlmostEqual(beta, -.4)
            point = np.array((.4, .5))
            ub = float(normal@point)
            self.assertEqual(classify_support(ub, ub, point, normal, offset, TAU, GEOMETRY_TOL)[0], 'VIOLATED')
            center, source, radius = interior_point(self.square, self.faces)
            self.assertEqual(source, 'chebyshev')
            np.testing.assert_allclose(center, (.6, .5), atol=1e-9)
            self.assertAlmostEqual(radius, .2)
            status, _, rho_ub = classify_center(ub, ub, point, normal, offset, center, TAU, GEOMETRY_TOL)
            self.assertEqual(status, 'SUPPORT_CERTIFIED')
            self.assertAlmostEqual(rho_ub, -TAU*.2)
            # 同一支撑值下两种裕量：origin 为 (1-tau)*UB+b = 0.4*tau > 0，center 为 -tau*(beta-n@c) < 0
            margins_origin = face_margins(self.faces[[f]], [ub], 'origin', center, TAU)
            margins_center = face_margins(self.faces[[f]], [ub], 'center', center, TAU)
            self.assertAlmostEqual(float(margins_origin[0]), .4*TAU)
            self.assertAlmostEqual(float(margins_center[0]), -.2*TAU)

    def test_auto_uses_center_without_origin_and_origin_when_clear(self):
        boundary = boundary_faces(self.faces)
        self.assertFalse(boundary.any())
        self.assertEqual(choose_criterion(self.faces, boundary, 'auto'), 'center')
        # 原点为顶点、两轴为分区边界面、其余面远离原点：origin
        corner = np.array(((0., 0.), (.7, 0.), (.6, .5), (0., .8)))
        faces = halfspaces(corner)
        boundary = boundary_faces(faces)
        self.assertEqual(int(boundary.sum()), 2)
        self.assertEqual(choose_criterion(faces, boundary, 'auto'), 'origin')
        # 一个非边界面经过原点附近（距离 < 1e-6）：center
        wedge = np.array(((0., 0.), (.7, 1e-5), (.5, .5)))
        faces = halfspaces(wedge)
        self.assertEqual(choose_criterion(faces, boundary_faces(faces), 'auto'), 'center')
        self.assertEqual(choose_criterion(faces, boundary_faces(faces), 'origin'), 'origin')

    def test_chebyshev_center_degenerate_falls_back_to_centroid(self):
        segment = np.array(((.2, .2), (.6, .6)))
        faces = halfspaces(segment)
        self.assertIsNone(chebyshev_center(faces))
        center, source, _ = interior_point(segment, faces)
        self.assertEqual(source, 'centroid')
        np.testing.assert_allclose(center, (.4, .4))


class ExpansionTests(unittest.TestCase):
    def test_center_expansion_matches_affine_image(self):
        polygon = np.array(((.2, .1), (.7, .2), (.8, .6), (.45, .85), (.15, .5)))
        faces = halfspaces(polygon)
        center, _, _ = interior_point(polygon, faces)
        expanded = expanded_faces(faces, 'center', center, TAU)
        image = center+(1.+TAU)*(polygon-center)
        self.assertTrue(same_faces(expanded, halfspaces(image)))
        np.testing.assert_allclose(-expanded[:, -1], allowed_offsets(faces, 'center', center, TAU))
        self.assertTrue(contains(image, expanded, 1e-12).all())

    def test_origin_expansion_matches_scaling(self):
        polygon = np.array(((0., 0.), (.6, 0.), (.7, .3), (0., .5)))
        faces = halfspaces(polygon)
        expanded = expanded_faces(faces, 'origin', None, TAU)
        self.assertTrue(same_faces(expanded, halfspaces(polygon/(1.-TAU))))

    def test_cone_faces_keeps_only_faces_cutting_the_region(self):
        region = np.array(((0., 0.), (.5, 0.), (.5, .5), (0., .5)))   # 锥∩盒的一块
        inside = halfspaces(np.array(((.1, .1), (.4, .1), (.4, .4), (.1, .4))))
        whole, disjoint, cut = cone_faces(region, inside)
        self.assertFalse(whole or disjoint)
        self.assertEqual(len(cut), 4)
        whole, disjoint, cut = cone_faces(region, halfspaces(np.array(((0., 0.), (1., 0.), (1., 1.), (0., 1.)))))
        self.assertTrue(whole)
        self.assertEqual(len(cut), 0)
        whole, disjoint, _ = cone_faces(region, halfspaces(np.array(((.6, .6), (.9, .6), (.9, .9), (.6, .9)))))
        self.assertTrue(disjoint)
        # 与区域部分相交：只保留切过区域的两条边
        corner = halfspaces(np.array(((.3, .3), (.9, .3), (.9, .9), (.3, .9))))
        whole, disjoint, cut = cone_faces(region, corner)
        self.assertFalse(whole or disjoint)
        self.assertTrue(same_faces(cut, corner[np.all(np.isclose(corner[:, :-1], (-1., 0.)), axis=1)
                                               | np.all(np.isclose(corner[:, :-1], (0., -1.)), axis=1)]))

    def test_expanded_polytope_is_clipped_to_the_box(self):
        polygon = np.array(((.0, .0), (1., 0.), (1., 1.), (0., 1.)))
        faces = halfspaces(polygon)
        center, _, _ = interior_point(polygon, faces)
        clipped = faces_polytope(expanded_faces(faces, 'center', center, TAU), 2)
        self.assertTrue(np.all(clipped >= -1e-12) and np.all(clipped <= 1.+1e-12))


class ValidityTests(unittest.TestCase):
    def setUp(self):
        # 人工扫描：8 个格。AC: 1 可行、-1 不可行、0 未决；SOCP 包含 AC。
        self.ac = np.array((1, 1, 1, -1, -1, 0, -1, -1))
        self.socp = np.array((1, 1, 1, 1, -1, 1, -1, 1))

    def test_valid_domains_pass(self):
        inner = np.array((1, 1, 0, 1, 0, 0, 0, 0), bool)   # 第 4 格 AC 不可行但 SOCP 可行：FR 计入、不是有效性错误
        outer = np.array((1, 1, 1, 1, 0, 1, 0, 0), bool)   # 第 8 格 SOCP-only 在外界外：只作诊断
        result = validity(inner, outer, self.ac, self.socp)
        self.assertTrue(result['valid'])
        self.assertEqual(result['outer_missed_socp_cells'], 1)
        self.assertEqual(result['undecided_ac_cells'], 1)
        metrics = grid_metrics(inner, outer, self.ac, self.socp)
        self.assertAlmostEqual(metrics['inner_ac']['fr_percent'], 100./3)   # 未决格不计入 AC 参考
        self.assertEqual(metrics['inner_ac']['computed_cells'], 3)

    def test_inner_with_socp_infeasible_cell_fails(self):
        inner = np.array((1, 0, 0, 0, 1, 0, 0, 0), bool)
        outer = np.ones(8, bool)
        result = validity(inner, outer, self.ac, self.socp)
        self.assertFalse(result['valid'])
        self.assertEqual(result['inner_socp_infeasible_cells'], 1)
        self.assertEqual(result['outer_missed_ac_cells'], 0)

    def test_outer_missing_ac_feasible_cell_fails(self):
        inner = np.zeros(8, bool)
        outer = np.array((1, 1, 0, 1, 1, 1, 1, 1), bool)
        result = validity(inner, outer, self.ac, self.socp)
        self.assertFalse(result['valid'])
        self.assertEqual(result['outer_missed_ac_cells'], 1)

    def test_undecided_ac_cell_is_not_infeasible(self):
        inner = np.array((0, 0, 0, 0, 0, 1, 0, 0), bool)   # 只含未决 AC 格
        outer = np.array((1, 1, 1, 0, 0, 0, 0, 0), bool)   # 外界不含未决格
        result = validity(inner, outer, self.ac, self.socp)
        self.assertTrue(result['valid'])
        self.assertEqual(result['inner_undecided_ac_cells'], 1)


if __name__ == '__main__':
    unittest.main()
