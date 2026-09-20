"""同一个模型在 33 节点上的数据、固定背景、独立 AC 与区域停止校验。"""
import unittest
import numpy as np
from Network.case33bw import network
from model import BranchEquations, SOCPSP, ACPowerFlow
from region import halfspaces
from vertify import validate_power_flow
from notebook_flow import workflow


class Case33Tests(unittest.TestCase):
    def test_original_data_and_nodal_power_flow(self):
        self.assertEqual(network.n, 32)
        self.assertEqual(len(network.branches), 32)
        self.assertEqual(network.original_p.sum(), 3715.)
        self.assertEqual(network.original_q.sum(), 2300.)
        validation = validate_power_flow(network)
        self.assertLess(validation['maximum_voltage_difference_pu'], 1e-8)
        self.assertAlmostEqual(validation['baseline']['minimum_voltage_pu'], .91309047936, places=9)

    def test_slice_background_is_fixed(self):
        power = np.array([[100., 200., 300.], [250., 400., 50.]])
        p, q = network.loads(power)
        fixed = np.ones(network.n, dtype=bool)
        fixed[network.selected] = False
        np.testing.assert_allclose(p[:, fixed]*network.base,
                                   np.tile(network.original_p[fixed], (2, 1)), atol=1e-12)
        np.testing.assert_allclose(q[:, fixed]*network.base,
                                   np.tile(network.original_q[fixed], (2, 1)), atol=1e-12)
        np.testing.assert_allclose(q[:, network.selected]*network.base, power*network.q_ratio)

    def test_socp_state_witness_and_cuts_against_ac(self):
        reference, equations, oracle = ACPowerFlow(network), BranchEquations(network), SOCPSP(network)
        samples = np.random.default_rng(20260922).random((120, 3))*[350., 1500., 600.]
        valid = samples[reference.classify(samples) == 1]
        cuts = 0
        for power in samples[:30]:
            eta, cut, witness, ell = oracle.solve(power)
            for actual, expected in zip(equations.state(witness, ell), reference.state(witness, ell)):
                np.testing.assert_allclose(actual, expected[0], atol=1e-12, rtol=0)
            self.assertGreaterEqual(equations.margin(equations.c+equations.F@witness+equations.G@ell), -1e-12)
            self.assertEqual(reference.classify(witness)[0], 1)
            self.assertEqual(eta < 1e-7, reference.classify(power)[0] == 1)
            if cut is not None:
                cuts += 1
                self.assertLess(cut[0]+cut[1:]@power, 0)
                self.assertTrue(np.all(cut[0]+valid@cut[1:] >= -1e-8))
        self.assertGreater(cuts, 0)
        oracle.close()

    def test_regions_meet_the_same_stopping_rule(self):
        for method in ('linear','socp','hybrid'):
            records, _ = workflow()['compute_regions']([network],np.inf,method,.002)
            record = records[0]
            inner, outer = np.array(record['inner']), np.array(record['outer'])
            hull = halfspaces(inner)
            tau = 0. if method == 'linear' else .002
            self.assertLessEqual(((1-tau)*outer@hull[:, :3].T+hull[:, 3]).max(), 1e-8)
            if method != 'linear':
                self.assertTrue(np.all(ACPowerFlow(network).classify(inner) == 1))
            else:
                e = BranchEquations(network)
                self.assertGreaterEqual(np.min(e.linear_c+outer@e.linear_F.T), -1e-9)

    def test_investment_changes_only_selected_impedances(self):
        self.assertEqual(len(network.designs),16)
        for design in network.designs:
            factor = np.ones(network.n)
            for selected, project in zip(design.x,network.projects):
                if selected:
                    factor[network.nodes.index(project.branch[1])] = .5
            np.testing.assert_allclose(design.r,network.r*factor)
            np.testing.assert_allclose(design.reactance,network.reactance*factor)
            np.testing.assert_array_equal(design.fixed_p,network.fixed_p)
            np.testing.assert_array_equal(design.fixed_q,network.fixed_q)
            self.assertEqual(design.cost,sum(x*p.cost for x,p in zip(design.x,network.projects)))

    def test_ac_minimum_cost_equals_all_designs_reference(self):
        points=np.random.default_rng(36).random((100,3))*[380,4500,600]
        plans=network.designs
        expected=np.full(len(points),np.inf)
        for design in plans:
            status=ACPowerFlow(design).classify(points)
            self.assertFalse(np.any(status==0))
            expected[status==1]=np.minimum(expected[status==1],design.cost)
        actual,_,_=workflow()['ac_reference'](plans,points,(0.,1.,2.,np.inf))
        np.testing.assert_array_equal(actual,expected)

    def test_planning_union_equals_separately_constructed_lp_domains(self):
        from model import LinearSP
        from region import simplex
        from vertify import region_membership
        plans=network.designs[:8]
        flow=workflow()
        records=flow['linear_region'](plans,2.)
        independent=[flow['cut_region'](LinearSP(d),simplex(d),0.)['outer'] for d in plans]
        points=np.random.default_rng(38).random((3000,3))*[450,5000,700]
        np.testing.assert_array_equal(region_membership(points,[r['outer'] for r in records]),
                                      region_membership(points,independent))


if __name__ == '__main__':
    unittest.main()
