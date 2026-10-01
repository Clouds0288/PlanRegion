"""方案二对照实验：center 判据、外扩公式与有效性检查（纯几何，不建求解器），以及 RCUT/RCUT2 的内域裁剪、
停滞规则与主线割平面循环（后者在 Case33 上建求解器）。"""
from copy import deepcopy
import unittest
from unittest.mock import patch

import numpy as np

from experiments.compare_plan2 import Case, CutState, cell_labels, cut_network, run_label, stagnated
from experiments.plan2_geometry import (allowed_offsets, boundary_faces, box_vertices, chebyshev_center,
                                        choose_criterion, classify_center, cone_faces, expanded_faces, face_margins,
                                        faces_polytope, grid_metrics, h_geometry, h_measures, interior_point, validity)
from experiments.test_radial_sandwich import partition_sign
from experiments.test_support_face_certification_fourbus_2d import classify_support
from model import PortPhysics, PortSubProblem, new_model
from region import GEOMETRY_TOL, contains, halfspaces, polytope_volume

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

    def test_approximate_inner_only_reports_the_inner_rule(self):
        inner = np.array((1, 0, 0, 0, 1, 0, 0, 0), bool)   # 含一个 SOCP 不可行格
        outer = np.ones(8, bool)
        exact, approximate = validity(inner, outer, self.ac, self.socp), validity(inner, outer, self.ac, self.socp, False)
        self.assertFalse(exact['valid'])
        self.assertTrue(approximate['valid'])
        self.assertEqual(approximate['inner_socp_infeasible_cells'], 1)
        missed = np.array((1, 1, 0, 1, 1, 1, 1, 1), bool)   # 外界条件照常检查
        self.assertFalse(validity(np.zeros(8, bool), missed, self.ac, self.socp, False)['valid'])
        self.assertTrue(grid_metrics(inner, outer, self.ac, self.socp, inner_exact=False)['valid'])


class CutRegionTests(unittest.TestCase):
    """RCUT/RCUT2 的内域 I=(I_R ∪ ∪N_x)∩O_R、停滞规则与网架集合的接受条件（纯几何）。"""

    def cone(self, d):
        # 整个卦限为一个锥：内三角顶点 0.5·e_j，外界 c@xi<=1.2，c=2·1（弦上 c@v=1），外域为腿长 0.6 的单纯形
        return [dict(U=np.eye(d), V=.5*np.eye(d), halfspace=(np.full(d, 2.), 1.2))]

    def test_clip_intersects_network_sets_with_the_cone_outer_bound(self):
        for d, square, outer in ((2, .64, .18), (3, .512, .036)):
            cube = np.array(list(np.ndindex(*(2,)*d)), float)*.8   # [0, 0.8]^d 伸出外界
            inner, measured_outer = h_measures(self.cone(d), [cube], None, d)
            self.assertAlmostEqual(inner, square)
            self.assertAlmostEqual(measured_outer, outer)
            inner, _ = h_measures(self.cone(d), [cube], None, d, clip=True)
            self.assertAlmostEqual(inner, outer)   # 立方体包含整个外域单纯形
            small = np.array(list(np.ndindex(*(2,)*d)), float)*.1   # 完全在外界内：裁剪不改变测度
            self.assertAlmostEqual(h_measures(self.cone(d), [small], None, d, clip=True)[0],
                                   h_measures(self.cone(d), [small], None, d)[0])
        inner, outer = h_geometry(self.cone(2), [np.array(((0., 0.), (.8, 0.), (.8, .8), (0., .8)))], None, clip=True)
        self.assertAlmostEqual(inner.area, outer.area, places=6)

    def test_cell_labels_clip_the_inner_region_only_when_flagged(self):
        state = dict(started=True, cones=[dict(U=np.eye(2), V=.5*np.eye(2), halfspace=((2., 2.), 1.2))],
                     networks=[dict(vertices=[[0., 0.], [.8, 0.], [.8, .8], [0., .8]])], cover=None)
        xi = np.array(((.1, .1), (.55, .01), (.7, .7), (.95, .95)))
        inner, outer = cell_labels(state, xi)
        np.testing.assert_array_equal(inner, (True, True, True, False))
        np.testing.assert_array_equal(outer, (True, True, False, False))
        inner, _ = cell_labels(dict(state, inner_clip=True), xi)
        np.testing.assert_array_equal(inner, (True, True, False, False))

    def test_stagnation_needs_patience_consecutive_small_cuts(self):
        self.assertTrue(stagnated([.5, .3, .001, .005, .009], .01, 3))
        self.assertFalse(stagnated([.5, .001, .02, .005, .009], .01, 3))   # 中间一次 >= threshold 重新计数
        self.assertFalse(stagnated([.001, .001], .01, 3))
        self.assertFalse(stagnated([.5, .3, .01, .005, .009], .01, 3))     # 等于 threshold 不算小割（同主线）

    def test_rcut_uses_only_stagnated_or_exact_sets_and_rcut2_the_certified_hull(self):
        square = np.array(((0., 0.), (.5, 0.), (.5, .5), (0., .5)))
        state = CutState((1, 0), '10', 1., 2, np.empty((0, 2)), None)
        state.vertices, state.inner = square, square[:3]
        for status, used in (('stagnated', True), ('exact', True), ('slice', False), ('failed', False),
                             ('point_resolution', False)):
            state.status = status
            self.assertEqual(len(state.region(False)) > 0, used, status)
            np.testing.assert_array_equal(state.region(True), square[:3])   # 认证点凸包任何状态都可用
        state.inner = square[:2]   # 低维认证集不计入内域
        self.assertEqual(len(state.region(True)), 0)

    def test_run_label_marks_non_default_thresholds(self):
        args = type('Args', (), dict(threshold=.005))()
        self.assertEqual(run_label('RCUT', args), 'RCUT-t0.5')
        self.assertEqual(run_label('RCUT2', type('Args', (), dict(threshold=.02))()), 'RCUT2-t2')
        self.assertEqual(run_label('RCUT', type('Args', (), dict(threshold=.01))()), 'RCUT')
        self.assertEqual(run_label('RB', args), 'RB')


class CutNetworkTests(unittest.TestCase):
    """Case33 (18,25) nn 分区、网架 0100010：主线割平面循环必须带 OBBT 盒与包络行，割不切掉认证点。"""

    @classmethod
    def setUpClass(cls):
        case = Case('case33', (18, 25))
        cls.bounds, cls.budget = case.bounds, case.budget
        cls.equations = PortPhysics(deepcopy(case.network), np.asarray(partition_sign('nn', 2)))
        network = cls.equations.network
        x = np.asarray(network.encode_plan(network.initial_plan), int).copy()
        switches = [(c, s) for c, s in zip(network.corridors, network.type_slices) if c.switchable]
        for (corridor, s), bit in zip(switches, '0100010'):
            x[s] = 0
            if bit == '1':
                ids = [k.id for k in corridor.types]
                x[s.start+(ids.index(corridor.existing_type) if corridor.existing_type in ids else 0)] = 1
        cls.x = x
        cls.oracle = PortSubProblem(cls.equations, threads=1)

    def run_cut(self, certify):
        return cut_network(self.equations, self.oracle, self.x, self.bounds, self.budget, box_vertices(2),
                           np.empty((0, 2)), None, 60., .01, 3, certify, 1)

    def rows(self):
        with new_model('rows', 1) as model:
            self.oracle._build(model, self.x, np.zeros(2))
            model.update()
            linear = model.NumConstrs
            model.remove(model.getQConstrs())   # 与 SubProblem._cut 相同：割 LP 只去掉锥约束
            model.update()
            return linear, model.NumConstrs

    def test_cut_loop_carries_the_obbt_rows_and_keeps_certified_points(self):
        scheme = tuple(int(v) for v in self.x)
        self.equations.boxes.pop(scheme, None)
        with self.assertRaises(RuntimeError):
            self.run_cut(False)
        plain = self.rows()
        self.equations.obbt(self.x, self.bounds, 1)
        rcut = self.run_cut(False)
        tight = self.rows()
        self.assertEqual(rcut['obbt_rows'], 2*len(self.equations.boxes[scheme])+2*int(self.x.sum()))
        self.assertEqual(tight[0]-plain[0], rcut['obbt_rows'])   # SP 带全部盒约束行与反向锥包络行
        self.assertEqual(tight[1]-plain[1], rcut['obbt_rows'])   # 割 LP 同样保留
        self.assertIn(rcut['status'], ('stagnated', 'exact'))
        self.assertTrue(rcut['history'] and rcut['counts']['lp'] == len(rcut['cuts']))
        rcut2 = self.run_cut(True)
        self.assertGreater(rcut2['counts']['rays'], 0)
        hull = rcut2['inner']
        self.assertGreater(polytope_volume(hull), 0.)
        # 有效割不会切掉认证点；RCUT2 的认证凸包在两种 N_x 内
        for N in (rcut['vertices'], rcut2['vertices']):
            self.assertTrue(contains(hull, halfspaces(N), 1e-7).all())
        self.assertGreaterEqual(polytope_volume(rcut['vertices']), polytope_volume(hull))


if __name__ == '__main__':
    unittest.main()
