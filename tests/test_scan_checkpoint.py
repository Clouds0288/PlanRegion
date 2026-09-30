"""已完成格点持久化；中断与未决状态不得变为不可行。"""
from unittest.mock import patch
import numpy as np
import pytest
from Network.case33bw import Case33
import vertify


def test_interrupted_scan_resumes_only_missing_rows(tmp_path):
    network = Case33()
    grid = dict(axis_lower=[-10., -10.], bounds=[10., 10.], shape=(17, 17))
    path = tmp_path/'ac.npz'
    solve = vertify.ac_scan_line
    calls = []
    def interrupted(args, **kwargs):
        calls.append(args[0])
        if len(calls) == 2:
            raise RuntimeError('deliberate interruption')
        return solve(args, **kwargs)
    with patch('vertify.ac_scan_line', side_effect=interrupted):
        with pytest.raises(RuntimeError, match='deliberate interruption'):
            vertify.scan_ac_reference(network, 7, grid, path, workers=1)
    saved = vertify.load_scan(next(path.glob('region_*.npz')))
    assert np.count_nonzero(saved['states']) == 128
    seen = []
    def resumed(args, **kwargs):
        seen.extend(args[0])
        return solve(args, **kwargs)
    with patch('vertify.ac_scan_line', side_effect=resumed):
        answer = vertify.scan_ac_reference(network, 7, grid, path, workers=1)
    assert set(seen) == set(range(128, 289))
    assert answer['ac_computed_points'] == 161 and answer['socp_computed_points'] == 289
    assert answer['reused_points'] == 0 and answer['computed_points'] == 289
    assert np.isin(answer['states'], [-1, 1]).all()
    assert np.isin(answer['socp_states'], [-1, 1]).all()


def test_unresolved_global_solve_is_saved_as_unknown(tmp_path):
    network = Case33()
    grid = dict(axis_lower=[-1., -1.], bounds=[1., 1.], shape=(1, 1))
    path = tmp_path/'ac.npz'
    with patch('vertify.signed_ac_witness', return_value=dict(feasible=np.array([False]), residual=np.array([np.inf]))), \
         patch('vertify.ac_interval_possible', return_value=np.array([True])), \
         patch('vertify.MasterProblem.solve', side_effect=TimeoutError('deliberate timeout')):
        with pytest.raises(TimeoutError, match='deliberate timeout'):
            vertify.scan_ac_reference(network, 7, grid, path, workers=1)
    saved = vertify.load_scan(next(path.glob('region_*.npz')))
    assert saved['states'][0, 0] == 0
    assert saved['socp_states'][0, 0] == 0


def test_completed_ac_is_reused_when_only_socp_is_missing(tmp_path):
    network = Case33(current_limit=200.)
    grid = dict(axis_lower=[-1., -1.], bounds=[1., 1.], shape=(2, 2))
    answer = vertify.scan_ac_reference(network, 7, grid, tmp_path, workers=1)
    answer['socp_states'].fill(0)
    answer['socp_residual'].fill(np.nan)
    vertify.save_scan(answer['cache_path'], answer)
    with patch('vertify.ac_scan_line', side_effect=AssertionError('AC already complete')):
        resumed = vertify.scan_ac_reference(network, 7, grid, tmp_path, workers=1)
    assert resumed['ac_computed_points'] == 0 and resumed['socp_computed_points'] == 4
    np.testing.assert_array_equal(resumed['states'], answer['states'])
    assert np.isin(resumed['socp_states'], [-1, 1]).all()
