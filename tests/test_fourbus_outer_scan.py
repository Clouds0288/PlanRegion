"""独立扫描的数值回归；参考解、投影归属与全局界分开验证。"""
import contextlib
import io
import unittest

import numpy as np

from experiments.fourbus_outer_scan import (candidate_eta, outer_problem, point_query,
                                             scan_grid, scan_rays)
from Network.four_bus_five_corridor import FourBus
from model import GridPhysics, MasterProblem, SubProblem


class FourBusScanTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.equations = GridPhysics(FourBus(), 'socp')
        initial = []
        for direction in [*np.eye(3), np.ones(3)]:
            problem = MasterProblem(cls.equations, budget=20000, direction=direction, threads=1)
            with problem.model:
                initial.append(problem.solve())
        x = cls.equations.network.encode_plan(cls.equations.network.initial_plan)
        sp = SubProblem(cls.equations, threads=1)
        cuts = [sp.solve(x, np.array(p, dtype=float))['cut'].tolist()
                for p in ((80, 20, 20), (5, 35, 30), (80, 0, 0))]
        cls.result = dict(budget=20000, axis_bounds=[r['bound'] for r in initial[:3]],
                          total_bound=initial[-1]['bound'], cuts=cuts)

    def test_reference_and_candidate_are_distinct_questions(self):
        e = self.equations
        x = e.network.encode_plan(e.network.initial_plan)
        power = np.array([80., 0., 0.])
        self.assertGreater(candidate_eta(e, x, power), 1e-5)
        direct = MasterProblem(e, budget=20000, threads=1)
        direct.model.setObjective(0.)
        with direct.model:
            self.assertTrue(point_query(direct, power))
            self.assertFalse(point_query(direct, np.array([100., 0., 0.])))

    def test_binary_search_matches_every_cut_prefix_projection(self):
        with contextlib.redirect_stdout(io.StringIO()):
            grid = scan_grid(self.equations, self.result, grid_step=25.)
        for k in range(len(self.result['cuts'])+1):
            outer = outer_problem(self.equations, self.result, k)
            outer.model.setObjective(0.)
            with outer.model:
                for p, excluded in zip(grid['power'], grid['first_exclusion']):
                    self.assertEqual(point_query(outer, p), excluded > k, (k, p))

    def test_direction_scan_contains_socp_and_recovers_axis_mp2(self):
        with contextlib.redirect_stdout(io.StringIO()):
            rays = scan_rays(self.equations, self.result, divisions=2)
        ref = np.array(rays['reference_total'])
        self.assertTrue(np.all(np.array(rays['reference_bound']) >= ref-1e-7))
        self.assertTrue(np.all(np.array(rays['outer_total']) >= ref-1e-6))
        for w, total in zip(rays['weights'], ref):
            if 1. in w:
                self.assertAlmostEqual(total, self.result['axis_bounds'][w.index(1.)], delta=1e-5)


if __name__ == '__main__':
    unittest.main()
