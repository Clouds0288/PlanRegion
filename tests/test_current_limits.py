"""电流单位、松弛与独立 AC 的限额一致性、缓存身份回归。"""
from unittest.mock import patch

import numpy as np
import pytest
from threadpoolctl import threadpool_limits

from Network.case33bw import Case33
from vertify import ac_interval_possible, signed_ac_witness
from model import GridPhysics, PortPhysics, PortSubProblem
from vertify import budget_schemes


def test_mainline_and_scans_keep_all_branch_current_limits():
    import main
    from vertify import ac_network, scan_problem
    with patch('main.run') as run:
        main.main('case33', dimension=2)
    network = run.call_args.args[0]
    amperes = np.sqrt(network.ell_limit)*network.base/(np.sqrt(3)*network.voltage_kv)
    np.testing.assert_allclose(amperes, np.full(37, 200.))
    limits = np.arange(180., 217.)
    varied = Case33(current_limit=limits)
    np.testing.assert_array_equal(ac_network(varied).ell_limit, varied.ell_limit)
    for ac in (False, True):
        problem = scan_problem(varied, 7, [-1, 1], ac=ac, power=[0., 0.])
        with problem.model as model:
            constraints = [c for c in model.getConstrs() if c.ConstrName.startswith('current_max[')]
            assert len(constraints) == 37
            for i, (c, limit) in enumerate(zip(constraints, varied.ell_limit)):
                assert -model.getCoeff(c, problem.x[i].item()) == pytest.approx(limit)
            cones = [c for c in model.getQConstrs() if c.QCName.startswith('current_cone[')]
            assert len(cones) == 37
            assert all(c.QCSense == ('=' if ac else '<') for c in cones)


@pytest.mark.parametrize('power', [[605.086876196641, 4646.37201998241], [366.473099, .0025]])
def test_200a_boundary_sp_preserves_original_certificate_tolerance(power):
    from tests.planning_checks import margin
    network = Case33(current_limit=200.)
    equations = PortPhysics(network, [1, -1])
    x = np.array([1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1,
                  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0])
    answer = PortSubProblem(equations, threads=1).solve(x, np.array(power))
    if answer['feasible']:
        assert margin(equations, x, np.array(power)*[1, -1], answer['state']) >= -1e-8
    else:
        assert answer['cut'][0]+answer['cut'][1:3]@power+answer['cut'][3:]@x < 0.


def test_high_export_socp_scan_preserves_feasible_points_in_both_orders():
    from vertify import scan_line
    network = Case33(current_limit=200.)
    power = np.array([[-6329.025, -307.75], [-6329.025, -158.05], [-6241.125, -1056.25]])
    for points in (power, power[::-1]):
        result = scan_line((np.arange(len(points)), points), network=network, budget=7)
        assert np.all(result['states'] == 1)
        assert np.max(result['residual']) <= 1e-8


def test_socp_scan_certifies_high_export_infeasibility():
    from vertify import scan_line
    power = np.array([[-6241.125, -1255.85], [-6241.125, -1455.45]])
    result = scan_line((np.arange(2), power), network=Case33(current_limit=200.), budget=7)
    np.testing.assert_array_equal(result['states'], [-1, -1])


def test_amperes_convert_to_squared_per_unit_in_models_and_tree():
    for limit in (210., 400.):
        network = Case33(current_limit=limit)
        expected = (limit/(10000/(np.sqrt(3)*12.66)))**2
        np.testing.assert_allclose(network.ell_limit, expected, rtol=1e-14)
        tree = network.tree(network.encode_plan(network.initial_plan))
        np.testing.assert_allclose(tree.ell_limit, expected, rtol=1e-14)
        for equations in (GridPhysics(network, 'socp'), PortPhysics(network, [-1, -1])):
            assert max(equations.ellmax.values()) <= expected*(1+1e-14)
    assert np.isinf(Case33().ell_limit).all()
    for value in (0., -210., np.nan):
        with pytest.raises(ValueError):
            Case33(current_limit=value)


def test_ac_current_limit_rejects_overloaded_physical_witness():
    power = [[-1000., -6000.], [90., 420.], [0., 0.]]
    for limit, expected in ((210., [False, True, True]), (400., [True, True, True])):
        network = Case33(current_limit=limit)
        x = network.encode_plan(network.initial_plan)
        answer = signed_ac_witness(network, x, power)
        np.testing.assert_array_equal(answer['feasible'], expected)
        np.testing.assert_array_equal(ac_interval_possible(network, x, power), expected)
        assert answer['residual'].max() < 1e-12
