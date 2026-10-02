"""固定 AC 参考与实际构域成员标签比较；矩形网格和未知标签回归。"""
from unittest.mock import patch
import numpy as np
import pytest
from Network.case33bw import Case33
from Network.four_bus_five_corridor import FourBus
from monitor import RunMonitor
from vertify import (scan_ac_reference, export_comparison, ac_identity, AC_CACHE_METHOD,
                     ac_scan_line, budget_schemes, signed_ac_witness, reference_box)
from model import GridPhysics


def test_fixed_ac_keeps_physical_point_above_artificial_current_cap(tmp_path):
    """同限流口径下，旧不限流可行点必须被实际限额排除。"""
    power = np.array([-1000., -6000.])
    network = Case33(current_limit=140.)
    answer = scan_ac_reference(network, 7, dict(axis_lower=power-.5, bounds=power+.5, shape=(1, 1)),
                               tmp_path/'ac.npz', workers=1)
    assert answer['states'][0, 0] == answer['socp_states'][0, 0] == -1
    assert np.isfinite(network.ell_limit).all()  # 原始 MISOCP 网架未被修改。


def test_ac_scanner_finds_a_different_topology_than_the_initial_plan(tmp_path):
    network = Case33()
    power = np.array([-5000., -3000.])
    initial = network.encode_plan(network.initial_plan)
    assert not signed_ac_witness(network, initial, power)['feasible'][0]
    schemes = budget_schemes(GridPhysics(network, 'socp'), 7)
    assert len(schemes) == 12
    answer = scan_ac_reference(network, 7, dict(axis_lower=power-.5,bounds=power+.5,shape=(1,1)),
                               tmp_path/'ac.npz',workers=1)
    assert answer['states'][0,0] == 1
    assert not np.array_equal(answer['witness_x'][0,0],initial)
    assert signed_ac_witness(network,answer['witness_x'][0,0],power)['feasible'][0]


def test_rectangular_export_coordinates_and_miss_denominators(tmp_path):
    network = Case33()
    ac = dict(axis_lower=np.array([-1., -1.5]), bounds=np.array([1., 1.5]),
              states=np.array([[1, -1, 1], [1, 1, -1]], np.int8), method=AC_CACHE_METHOD,
              socp_states=np.array([[1, 1, 1], [1, 1, -1]], np.int8),
              metadata=dict(identity=ac_identity(network, 7)))
    vertices = np.array([[-1., -1.5], [0., -1.5], [0., 1.5], [-1., 1.5]])
    result = dict(inner=[dict(vertices=vertices)], outer=[dict(vertices=vertices)],
                  status='time_limit', certified=False)
    summary = export_comparison(network, 7, result, ac, tmp_path)
    metrics = summary['metrics']['inner']
    assert metrics['reference_cells'] == 4 and metrics['computed_cells'] == 3
    assert metrics['missed_cells'] == 2 and metrics['extra_cells'] == 1
    assert metrics['mr_percent'] == 50. and metrics['fr_percent'] == pytest.approx(100/3)
    csv = np.genfromtxt(tmp_path/'comparison.csv', delimiter=',', names=True)
    assert csv.dtype.names == ('p_18_kw', 'p_25_kw', 'ac_state', 'socp_state', 'inner', 'outer')
    comparisons = summary['comparisons']
    assert comparisons['result_ac'] == metrics
    assert comparisons['result_socp']['mr_percent'] == 40.
    assert comparisons['result_socp']['fr_percent'] == 0.
    assert comparisons['socp_ac']['mr_percent'] == 0.
    assert comparisons['socp_ac']['fr_percent'] == 20.
    with np.load(tmp_path/'comparison.npz') as saved:
        np.testing.assert_array_equal(saved['socp_states'], csv['socp_state'])
        np.testing.assert_array_equal(saved['power'], [[-.5,-1.],[-.5,0.],[-.5,1.],[.5,-1.],[.5,0.],[.5,1.]])
        np.testing.assert_array_equal(saved['ac_states'], csv['ac_state'])


def test_unknown_ac_points_are_excluded_and_counted():
    monitor = RunMonitor()
    monitor.validation(dict(bounds=[2., 1.], states=[[0], [1]], socp_states=[[0], [1]]),
                       dict(inner=[dict(vertices=[[1., 0.], [2., 0.], [2., 1.], [1., 1.]])]))
    validation = monitor.state['validation']
    assert validation['undecided_cells'] == validation['socp_undecided_cells'] == 1
    assert validation['mr_percent'] == validation['fr_percent'] == 0.   # 只按已决的一格计算


def test_three_dimensional_serial_parallel_and_full_ac_agree(tmp_path):
    network = FourBus()
    grid = dict(axis_lower=np.zeros(3), bounds=np.array([65.,55.,45.]), shape=(2,3,2))
    serial = scan_ac_reference(network, 20000., grid, tmp_path/'serial.npz', workers=1, mode=0)
    parallel = scan_ac_reference(network, 20000., grid, tmp_path/'parallel.npz', workers=2, mode=0)
    np.testing.assert_array_equal(serial['states'], parallel['states'])
    power = (np.indices(grid['shape']).reshape(3,-1).T+.5)*grid['bounds']/np.array(grid['shape'])
    full = ac_scan_line((np.arange(len(power)),power),network=network,budget=20000.,schemes=None,mode=0)
    assert not full['errors']
    np.testing.assert_array_equal(full['states'],serial['states'].ravel())


def test_ac_reference_bound_is_kw_and_independent_of_artificial_cap(tmp_path):
    """同限流后的包络保留 kW 单位，复用必须匹配电流限额。"""
    lower, upper = reference_box(Case33(current_limit=200.), 7, output=tmp_path)
    assert np.all(lower < -1000.) and np.all(upper > 100.)
    with patch('vertify.MasterProblem', side_effect=AssertionError('must reuse bounds')):
        same = reference_box(Case33(current_limit=200.), 7, output=tmp_path)
        with pytest.raises(AssertionError, match='must reuse bounds'):
            reference_box(Case33(current_limit=140.), 7, output=tmp_path)
    np.testing.assert_array_equal(same, (lower, upper))
