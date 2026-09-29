"""AC 等式的数值核验：假可行点被排除，真实点可恢复节点交流方程。"""
import unittest
import numpy as np

from experiments.signed_ac_equality import build_equality_problem, solve_trial
from experiments.signed_power_compare import Case33, SignedPhysics


class SignedAcEqualityTests(unittest.TestCase):
    def setUp(self):
        self.net = Case33()
        self.equations = SignedPhysics(self.net, [-10000., -10000.], [10000., 10000.])

    def test_rejects_socp_false_positive_across_all_allowed_topologies(self):
        options = dict(budget=7, seconds=30., threads=4, power=[-10000., 0.])
        relaxed = solve_trial(self.equations, 'socp', **options)
        exact = solve_trial(self.equations, 'ac', **options)
        self.assertEqual(relaxed['status'], 'optimal')
        self.assertGreater(relaxed['current_gap_max'], 1e-3)
        self.assertEqual(exact['status'], 'infeasible')

    def test_negative_injection_recovers_independent_nodal_equations(self):
        power = [-1000., 0.]
        problem = build_equality_problem(self.equations, budget=7, power=power, threads=4)
        with problem.model:
            answer = problem.solve(time_limit=30.)
        self.assertIsNotNone(answer)
        tree = self.net.tree(answer['x'])
        state = answer['state']
        P, Q, ell = (state[part][tree.type_indices] for part in (
            self.equations.P_slice, self.equations.Q_slice, self.equations.ell_slice))
        # 参考方向反转时，树送端功率等于 -(P-r*ell)，无功同理。
        P = np.where(tree.direction > 0, P, -P+tree.r*ell)
        Q = np.where(tree.direction > 0, Q, -Q+tree.reactance*ell)
        voltage = np.ones(tree.n+1, dtype=complex)
        admittance = np.zeros((tree.n+1, tree.n+1), dtype=complex)
        for i in tree.order:
            parent = tree.parent[i] if tree.parent[i] >= 0 else tree.n
            impedance = complex(tree.r[i], tree.reactance[i])
            current = np.conj(complex(P[i], Q[i])/voltage[parent])
            voltage[i] = voltage[parent]-impedance*current
            admittance[i, i] += 1/impedance
            admittance[parent, parent] += 1/impedance
            admittance[i, parent] -= 1/impedance
            admittance[parent, i] -= 1/impedance
        p, q = tree.loads(power)
        injection = voltage*np.conj(admittance@voltage)
        np.testing.assert_allclose(injection[:tree.n], -p[0]-1j*q[0], atol=1e-7)
        np.testing.assert_allclose(abs(voltage[:tree.n])**2,
            state[self.equations.v_slice][tree.node_indices], atol=1e-7)


if __name__ == '__main__':
    unittest.main()
