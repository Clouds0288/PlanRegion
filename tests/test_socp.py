"""跨网架的 SOCP 状态、内点、有效割与原线性认证区域回归验证。"""
import json
from pathlib import Path
import re
import unittest
import numpy as np
from Network.four_bus_five_corridor import network
from model import SOCPSP, ACPowerFlow
from region import halfspaces, simplex
from vertify import region_membership
from notebook_flow import workflow


class SOCPTests(unittest.TestCase):
    def test_all_designs_and_dual_cuts_against_ac(self):
        rng = np.random.default_rng(20260920)
        cuts = feasible_queries = 0
        for design in network.designs:
            oracle, reference = SOCPSP(design), ACPowerFlow(design)
            validation = rng.dirichlet(np.ones(4), size=30)[:, :3]*network.power_limit
            valid = validation[reference.classify(validation) == 1]
            for power in (np.array([2., 3., 1.]), validation[0]):
                eta, cut, witness, ell = oracle.solve(power)
                P, Q, v, u = oracle.equations.state(witness, ell)
                for expected, value in zip(reference.state(witness, ell), (P, Q, v, u)):
                    np.testing.assert_allclose(value, expected[0], atol=1e-12, rtol=0)
                self.assertLessEqual(np.max(P*P+Q*Q-u*ell), 1e-9)
                self.assertLessEqual(reference.violation(P[None], Q[None], v[None]).max(), 1e-9)
                self.assertEqual(reference.classify(witness)[0], 1)
                self.assertEqual(eta <= 1e-7, reference.classify(power)[0] == 1)
                if cut is not None:
                    cuts += 1
                    self.assertLess(cut[0]+cut[1:]@power, 0)
                    self.assertTrue(np.all(cut[0]+valid@cut[1:] >= -1e-7))
                else:
                    feasible_queries += 1
            oracle.close()
        self.assertGreater(cuts, 0)
        self.assertGreater(feasible_queries, 0)

    def test_radial_stopping_certificate(self):
        for design in network.designs[::43]:
            result = workflow()['cut_region'](SOCPSP(design), simplex(network), .002)
            inside = halfspaces(result['inner'])
            self.assertLessEqual(np.max((1-.002)*result['outer']@inside[:, :3].T
                                        + inside[:, 3]), 1e-8)
            self.assertTrue(np.all(ACPowerFlow(design).classify(result['inner']) == 1))

    def test_zero_load_subtree_does_not_collapse_witness(self):
        oracle = SOCPSP(network.designs[13])
        _, _, witness, currents = oracle.solve(np.array([32.8, 0., 0.]))
        self.assertGreater(witness[0], 32.)
        np.testing.assert_array_equal(currents[1:], 0.)
        oracle.close()

    def test_extracted_linear_solver_reproduces_saved_region(self):
        plans = [d for d in network.designs if d.cost<=20000.]
        records = workflow()['linear_region'](plans,20000.)
        html = Path('results/four_bus_five_corridor.html').read_text(encoding='utf-8')
        data = json.loads(re.search(r'<script id="experiment" type="application/json">(.*?)</script>', html, re.S)[1])
        scenario = next(s for s in data['scenarios'] if s['budget'] == 20000.)
        groups = {}
        for row in scenario['history']:
            if row['x'] is not None and row['eta'] is not None and row['cut'] is None:
                groups.setdefault(tuple(row['x']), []).append(row['p'])
        points = np.random.default_rng(51).dirichlet(np.ones(4), size=2000)[:, :3]*network.power_limit
        np.testing.assert_array_equal(region_membership(points, [np.array(p) for p in groups.values()]),
                                      region_membership(points, [r['outer'] for r in records if len(r['outer'])]))


if __name__ == '__main__':
    unittest.main()
