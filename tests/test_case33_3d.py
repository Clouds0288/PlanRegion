"""Independent LP and geometry regressions for the three-dimensional experiment."""
import unittest
import numpy as np
from scipy.optimize import linprog

from tests.case33_3d import (MASKS, downward_facets, hull_distance, cell_gap,
                            subtract_orthant, apply_orthant, select_cell, Cell, PlanningOracle3D,
                            subtract_disjunction,cover_partition)
from region import polytope_volume


class Geometry3DTests(unittest.TestCase):
    def test_coordinate_face_distance(self):
        facets = downward_facets([[0, 0, 3]])
        np.testing.assert_allclose(hull_distance([[1, 2, 10], [0, 0, 3]], facets),
                                   [7.00001, .00001], atol=2e-6)
        np.testing.assert_allclose(hull_distance([[1, 2, 1], [3, 0, 1]], facets), [2., 3.], atol=1e-7)

    def test_distance_against_linear_program(self):
        rng = np.random.default_rng(912)
        points = rng.uniform(0, 30, (13, 3))
        facets = downward_facets(points)
        queries = np.r_[rng.uniform(0, 45, (50, 3)), [[0, 40, 0], [40, 0, 0], [0, 0, 40]]]
        expected = []
        for q in queries:
            constraints = np.r_[np.c_[facets[:, :3], np.zeros(len(facets))], np.c_[-np.eye(3), -np.ones(3)]]
            rhs = np.r_[-facets[:, 3], -q]
            answer = linprog([0, 0, 0, 1], A_ub=constraints, b_ub=rhs, bounds=[(0, None)]*4,
                             method='highs')
            self.assertTrue(answer.success)
            expected.append(answer.fun)
        np.testing.assert_allclose(hull_distance(queries, facets), expected, atol=1e-7)

    def test_different_schemes_do_not_cover_middle(self):
        hulls = [((1,), downward_facets([[10, 1, 1]])),
                 ((2,), downward_facets([[1, 10, 1]]))]
        vertices = np.array([[10, 1, 1], [1, 10, 1]])
        gap, _, individual = cell_gap(vertices, hulls)
        self.assertLess(individual.max(), .001)
        self.assertGreater(gap, 8.9)
        self.assertGreater(min(hull_distance([[5.5, 5.5, 1]], h)[0] for _, h in hulls), 4.4)

    def test_orthant_subtraction_volume(self):
        vertices = MASKS*10
        pieces = subtract_orthant(vertices, [3, 4, 5])
        self.assertAlmostEqual(sum(polytope_volume(v) for v in pieces), 1000-7*6*5, places=7)
        for v in pieces:
            self.assertTrue(np.all(np.any(v <= np.array([3, 4, 5])+1e-8, axis=1)))

    def test_orthant_negative_coordinate_threshold(self):
        pieces = subtract_orthant(MASKS*10, [-1, 4, 5])
        self.assertAlmostEqual(sum(polytope_volume(v) for v in pieces), 1000-10*6*5, places=7)

    def test_oblique_disjunction_volume(self):
        normals = np.array([[.5,.5,0.],[0.,0.,1.]])
        # Removed: x+y>10 and z>5, triangular prism of volume 250.
        pieces = subtract_disjunction(MASKS*10, [5.,5.], normals)
        self.assertAlmostEqual(sum(polytope_volume(v) for v in pieces), 750., places=7)
        self.assertTrue(all(np.all(np.any(v @ normals.T <= 5.+1e-8, axis=1)) for v in pieces))

    def test_coverage_partition_conserves_volume_and_certifies_covered_piece(self):
        facets = downward_facets([[3.,10.,0.],[0.,10.,3.]])
        pieces = cover_partition(MASKS*10,facets,2.)
        self.assertIsNotNone(pieces)
        self.assertAlmostEqual(sum(polytope_volume(v) for v in pieces),1000.,places=6)
        self.assertLessEqual(hull_distance(pieces[-1],facets).max(),2.+1e-8)
        # A single shift of the original facets would incorrectly cover this
        # point near a coordinate face; the coordinate masks prevent that.
        self.assertGreater(hull_distance([[0.,0.,6.]],facets)[0],2.)
        self.assertLessEqual(np.max(pieces[-1][:,2]),5.+1e-7)
        self.assertAlmostEqual(polytope_volume(pieces[-1]),205.,delta=.002)

    def test_orthant_outside_cube_preserves_volume(self):
        pieces = subtract_orthant(MASKS*10, [11, 4, 5])
        self.assertAlmostEqual(sum(polytope_volume(v) for v in pieces), 1000, places=7)

    def test_tangent_coverage_does_not_repeat_unchanged_node(self):
        facets = np.array([[1.,0.,0.,-2.]])
        for vertices in (MASKS*10+np.array([3.,0.,0.]),
                         np.array([[3.,1.,1.],[13.,1.,1.]])):
            self.assertIsNone(cover_partition(vertices,facets,1.))
        # A genuine partition within a coordinate face must still be available.
        vertices = np.unique(MASKS*np.array([4.,10.,0.]),axis=0)
        pieces = cover_partition(vertices,facets,1.)
        self.assertIsNotNone(pieces)
        self.assertAlmostEqual(float(pieces[-1][:,0].max()),3.)
        self.assertTrue(any(np.max(piece[:,0]) == 4. for piece in pieces))

    def test_convex_cell_bound_dominates_interior(self):
        rng = np.random.default_rng(713)
        vertices = MASKS*10
        hulls = [((1,), downward_facets([[9, 4, 7], [3, 10, 8]])),
                 ((2,), downward_facets([[8, 9, 2], [7, 2, 10]]))]
        gap, _, _ = cell_gap(vertices, hulls)
        points = rng.uniform(0, 10, (100, 3))
        sampled = np.min([hull_distance(points, h) for _, h in hulls], axis=0)
        self.assertLessEqual(sampled.max(), gap)

    def test_global_cut_reaches_other_nodes(self):
        left = MASKS*np.array([5.,10.,10.])
        right = left+np.array([5.,0.,0.])
        children = apply_orthant([Cell(left, gap=4.), Cell(right, gap=7.)], np.array([3.,4.,5.]))
        self.assertAlmostEqual(sum(polytope_volume(c.vertices) for c in children), 790., places=6)
        self.assertTrue(all(c.gap in (4.,7.) for c in children))

    def test_lazy_selection_matches_fully_recomputed_maximum(self):
        before = [((1,), downward_facets([[8.,2.,2.]])),
                  ((2,), downward_facets([[2.,8.,2.]]))]
        after = [((1,), downward_facets([[8.,2.,2.],[7.,7.,6.]])),
                 ((2,), downward_facets([[2.,8.,2.],[2.,8.,8.]]))]
        cells = [Cell(MASKS*5+offset) for offset in MASKS*5]
        for cell in cells:
            cell.gap,cell.scheme,_ = cell_gap(cell.vertices,before)
        expected = max(cell_gap(cell.vertices,after)[0] for cell in cells)
        chosen,_ = select_cell(cells,after)
        self.assertAlmostEqual(cells[chosen].gap,expected,places=7)

    def test_surface_cancels_internal_cell_face(self):
        from tests.plot_case33_3d import surface
        left = MASKS*np.array([5.,10.,10.])
        mesh = surface([left, left+np.array([5.,0.,0.])])
        vertices = np.asarray(mesh['vertices'])
        triangles = vertices[np.asarray(mesh['triangles'])]
        area = np.linalg.norm(np.cross(triangles[:,1]-triangles[:,0], triangles[:,2]-triangles[:,0]), axis=1).sum()/2
        self.assertAlmostEqual(area, 600., places=6)

    def test_frame_archive_reconstructs_face_sets(self):
        from tests.plot_case33_3d import surface, GeometryArchive
        archive = GeometryArchive()
        state = dict(t=set(),l=set())
        for cells in ([MASKS*10], subtract_orthant(MASKS*10,[3.,4.,5.])):
            mesh = surface(cells)
            patch = archive.encode(mesh,patch=True)
            for key in ('t','l'):
                state[key].difference_update(patch[key+'_remove'])
                state[key].update(patch[key+'_add'])
            complete = archive.encode(mesh)
            self.assertEqual(state['t'],set(complete['t']))
            self.assertEqual(state['l'],set(complete['l']))

    def test_scan_rejections_require_global_evidence(self):
        from tests.audit_case33_3d import infeasibility_mask
        scan = dict(fixed_value=0.,events=[
            dict(fixed={'0':0.,'1':0.,'2':0.},infeasible=False,bound=None),
            dict(fixed={'0':10.,'2':0.},infeasible=False,bound=15.),
            dict(fixed={'0':20.,'1':0.,'2':0.},infeasible=True,bound=None)])
        proof = infeasibility_mask(scan,np.array([0.,10.,20.]),[20.,20.,20.],50.)
        np.testing.assert_equal(proof,[[False,False,False],[False,False,True],[True,True,True]])

    def test_3d_load_distance_preserves_fixed_plan_optimum(self):
        objectives = []
        for strengthen in (False, True):
            oracle = PlanningOracle3D(0., threads=1, strengthen=strengthen)
            net = oracle.network
            choice = {c.id: c.existing_type if c.initial_active else None for c in net.corridors}
            x = net.encode_plan(choice)
            oracle.problem.x.LB = oracle.problem.x.UB = x
            answer = oracle.solve(target=[900.,3000.,900.], time_limit=10., gap_kw=.01,
                                  warm_count=0, label='fixed_plan_3d_regression')
            self.assertTrue(answer['feasible'])
            self.assertGreater(answer['bound'], 0.)
            self.assertLess(answer['objective']-answer['bound'], .011)
            objectives.append(answer['objective'])
            target = answer['p']+10.
            covered = oracle.solve(target=target, time_limit=10., stop_distance_kw=19.999,
                                   label='warm_coverage_regression')
            self.assertTrue(covered['feasible'])
            self.assertFalse(covered['global_solve'])
            self.assertIsNone(covered['bound'])
            self.assertLess(np.max(target-covered['p']), 20.)
            projected = oracle.solve(target=[900.,3000.,900.], cut_normals=np.eye(3),
                                     time_limit=10., gap_kw=.01, label='identity_projection_regression')
            self.assertEqual(projected['mode'],'projected_distance')
            self.assertTrue(projected['feasible'])
            self.assertAlmostEqual(projected['objective'],answer['objective'],delta=.02)
            for certificate in oracle.certificates:
                self.assertGreaterEqual(np.max(np.array([900.,3000.,900.])-projected['bound']-certificate['p']),-1e-3)
            oracle.close()
        self.assertAlmostEqual(*objectives, delta=.02)


if __name__ == '__main__':
    unittest.main()
