"""用原始 SP 和全部外域顶点独立检查 KKT 全局模型的等价性。"""
import contextlib
import io
import unittest

import numpy as np
from gurobipy import GRB
from threadpoolctl import threadpool_limits

from experiments.fourbus_outer import GlobalViolation
from experiments.fourbus_outer_kkt import KktViolation, initialize, run_trial
from experiments.fourbus_outer_scan import candidate_eta
from experiments.fourbus_outer_volume import budget_schemes
from model import GridPhysics, MasterProblem
from Network.four_bus_five_corridor import FourBus
from region import clip_polytope, initial_polytope


class KktViolationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pool = threadpool_limits(limits=1)
        cls.pool.__enter__()
        cls.equations = GridPhysics(FourBus(load_nodes=(1, 2)), 'socp')
        cls.initial = initialize(20000., 4)
        cls.schemes = budget_schemes(cls.equations, 20000., 1)
        with contextlib.redirect_stdout(io.StringIO()):
            cls.result = run_trial('kkt', cls.initial, epsilon=.01, seconds=100., time_limit=5., threads=4)

    @classmethod
    def tearDownClass(cls):
        cls.pool.__exit__(None, None, None)

    def search(self, kind=KktViolation):
        return kind(self.equations, 20000., self.initial['axis_bounds'],
                    self.initial['total_bound'], threads=4)

    def vertex_reference(self, cuts):
        """固定 x 后 eta*(p) 凸，故多面体上的最大值在某个顶点取得。"""
        bounds = np.asarray(self.initial['axis_bounds'])
        base = initial_polytope(bounds, self.initial['total_bound'], bounds)
        maximum = 0.
        for x in self.schemes:
            vertices = base.copy()
            for cut in cuts:
                vertices = clip_polytope(vertices, cut[0]+cut[3:]@x, cut[1:3]*bounds)
                if not len(vertices):
                    break
            for power in vertices*bounds:
                maximum = max(maximum, candidate_eta(self.equations, x, power))
        return maximum

    def test_initial_global_value_matches_all_17_schemes_and_vertices(self):
        reference = self.vertex_reference([])
        for kind in (GlobalViolation, KktViolation):
            search = self.search(kind)
            with search.model, self.subTest(method=kind.__name__):
                answer = search.solve(30., 0.)
                self.assertEqual(answer['status'], GRB.OPTIMAL)
                self.assertAlmostEqual(search.model.ObjVal, reference, delta=2e-6)
                self.assertLessEqual(answer['objective'], reference+2e-7)
                self.assertGreaterEqual(answer['bound'], reference-2e-7)

    def test_after_joint_cuts_kkt_value_and_certificate_match_primal_vertices(self):
        cuts = [np.asarray(cut) for cut in self.result['cuts']]
        self.assertTrue(cuts)
        self.assertTrue(any(np.any(cut[3:]) for cut in cuts))
        reference = self.vertex_reference(cuts)
        search = self.search()
        with search.model:
            for cut in cuts:
                search.add_cut(cut)
            answer = search.solve(30., 0.)
            self.assertEqual(answer['status'], GRB.OPTIMAL)
            self.assertAlmostEqual(search.model.ObjVal, reference, delta=2e-6)
            self.assertLessEqual(answer['objective'], reference+2e-7)
            self.assertGreaterEqual(answer['bound'], reference-2e-7)
        self.assertTrue(self.result['certified'])
        self.assertLessEqual(reference, self.result['bound']+2e-7)
        self.assertLessEqual(self.result['bound'], .01)

    def test_no_continuous_bilinear_terms_or_relaxed_cones(self):
        search = self.search()
        with search.model:
            for cut in self.result['cuts']:
                search.add_cut(np.asarray(cut))
            model = search.model
            model.update()
            self.assertEqual(model.Params.NonConvex, 0)
            self.assertEqual(model.NumQNZs, 0)
            self.assertEqual(model.NumQConstrs, len(search.template.cone_slices))
            self.assertEqual(len(search.outer_dual), 6+len(self.result['cuts']))
            self.assertTrue(all(0. <= var.UB < GRB.INFINITY for var in search.outer_dual))
            for row in model.getQConstrs():
                expression = model.getQCRow(row)
                self.assertTrue(all(expression.getVar1(i).sameAs(expression.getVar2(i))
                                    for i in range(expression.size())))

    def test_every_generated_cut_is_valid_without_budget_restriction(self):
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
