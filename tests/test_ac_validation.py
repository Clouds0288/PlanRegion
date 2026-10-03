"""带符号 AC 独立证书与 Case33 扫描失败点的数值回归。"""
import numpy as np
from threadpoolctl import threadpool_limits

from vertify import signed_ac_witness, scan_line, ac_interval_possible, export_comparison
from vertify import budget_schemes
from Network.case33bw import Case33
from model import LOAD_PF
from monitor import grid_comparison
from vertify import ACPowerFlow, AC_CACHE_METHOD, ac_identity




def test_signed_witness_satisfies_independent_complex_power_balance():
    network = Case33(load_nodes=(18, 25, 30))
    x = network.encode_plan(network.initial_plan)
    power = np.array([[100., 100., 100.], [-1000., -4000., 100.]])
    answer = signed_ac_witness(network, x, power)
    assert answer['feasible'].all()
    tree = network.tree(x)
    ratio = np.where(power >= 0., np.tan(np.arccos(LOAD_PF)), 0.)
    p = (tree.fixed_p+power@tree.E.T)/tree.base
    q = (tree.fixed_q+(power*ratio)@tree.E.T)/tree.base
    P, Q, v, _ = ACPowerFlow(tree)._state(p, q, answer['ell'])
    voltage, current = np.ones_like(P, dtype=complex), np.zeros_like(P, dtype=complex)
    for i in tree.order:
        parent = 1.+0j if tree.parent[i] < 0 else voltage[:, tree.parent[i]]
        current[:, i] = np.conj((P[:, i]+1j*Q[:, i])/parent)
        voltage[:, i] = parent-(tree.r[i]+1j*tree.reactance[i])*current[:, i]
    demand = voltage*np.conj(current-np.column_stack([current[:, children].sum(axis=1) for children in tree.children]))
    np.testing.assert_allclose(demand, p+1j*q, atol=1e-10, rtol=0.)
    np.testing.assert_allclose(abs(voltage)**2, v, atol=1e-10, rtol=0.)
    assert P[1, tree.roots].sum() < 0.


def test_no_signed_witness_does_not_mean_infeasible():
    network = Case33(load_nodes=(18, 25, 30))
    power = np.array([-1000., -4000., 100.])
    # 不给初步见证，让同一个可行反送点走完整 AC 等式模型。
    with threadpool_limits(limits=1):
        answer = scan_line((np.array([0]), power[None]), network=network, budget=7, ac=True)
    assert answer['states'][0] == 1




def test_ac_necessary_intervals_preserve_physical_witnesses():
    network = Case33(load_nodes=(18, 25, 30))
    schemes = budget_schemes(network, 7)
    rng = np.random.default_rng(72)
    power = rng.uniform([-6000., -7000., -6000.], [1000., 3500., 1500.], (160, 3))
    witnesses = 0
    with threadpool_limits(limits=1):
        for x in schemes:
            feasible = signed_ac_witness(network, x, power)['feasible']
            possible = ac_interval_possible(network, x, power)
            assert np.all(possible[feasible])
            witnesses += int(feasible.sum())
    assert witnesses > 100


def test_ac_interval_exclusion_agrees_with_global_equality():
    network = Case33(load_nodes=(18, 25, 30))
    schemes = budget_schemes(network, 7)
    power = np.array([-9866.61875, -10986.8375, -3058.11875])
    with threadpool_limits(limits=1):
        assert all(not ac_interval_possible(network, x, power)[0] for x in schemes)
        answer = scan_line((np.array([0]), power[None]), network=network, budget=7, ac=True)
    assert answer['states'][0] == -1


def test_comparison_export_preserves_coordinates_and_labels(tmp_path):
    network = Case33(load_nodes=(18, 25, 30))
    lower, upper = -np.ones(3), np.ones(3)
    vertices = np.array(list(np.ndindex(2, 2, 2)))*2.-1.
    result = dict(status='time_limit', certified=False,
                  inner=[dict(vertices=vertices)], outer=[dict(vertices=vertices)])
    ac = dict(axis_lower=lower, bounds=upper, states=np.ones((2, 2, 2), dtype=np.int8),
              socp_states=np.ones((2, 2, 2), dtype=np.int8),
              method=AC_CACHE_METHOD, metadata=dict(identity=ac_identity(network, 7)))
    ac['states'][0, 0, 1] = -1
    summary = export_comparison(network, 7, result, ac, grid_comparison(ac, result), tmp_path/'run.json.gz')
    with np.load(tmp_path/'run_comparison'/'comparison.npz') as saved:   # 导出到记录旁
        np.testing.assert_array_equal(saved['power'][1], [-.5, -.5, .5])
        assert saved['ac_states'][1] == -1 and saved['inner'][1]
        np.testing.assert_array_equal(saved['load_nodes'], [18, 25, 30])
        assert np.all(saved['socp_states'] == 1)
    assert not summary['result_certified']
    assert summary['metrics']['inner']['fr_percent'] == 12.5
