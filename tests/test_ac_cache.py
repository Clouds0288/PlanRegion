"""AC 缓存身份、真实坐标复用、扩界、加密和历史迁入。"""
from copy import deepcopy
from unittest.mock import patch
import numpy as np
import pytest
import vertify
from Network.case33bw import Case33


def test_covering_cache_reuses_without_starting_solver(tmp_path):
    network = Case33()
    grid = dict(axis_lower=[-10., -10.], bounds=[10., 20.], shape=(2, 3))
    path = vertify.scan_path(network, 7, tmp_path)
    first = vertify.scan_ac_reference(network, 7, grid, path, workers=1)
    requested = dict(axis_lower=[-9., -9.], bounds=[9., 19.], shape=(160, 160))
    with patch('vertify.ac_scan_line', side_effect=AssertionError('must reuse')):
        reused = vertify.scan_ac_reference(network, 7, requested, path, workers=1)
    np.testing.assert_array_equal(reused['states'], first['states'])
    np.testing.assert_array_equal(reused['step'], [10., 10.])
    assert reused['computed_points'] == 0 and reused['reused_points'] == 6


@pytest.mark.parametrize('dimension', [2, 3])
def test_expansion_reuses_overlap_and_only_submits_new_coordinates(tmp_path, dimension):
    network = Case33(load_nodes=(18, 25, 30)[:dimension])
    grid = dict(axis_lower=np.full(dimension, -10.), bounds=np.full(dimension, 10.), shape=(2,)*dimension)
    path = vertify.scan_path(network, 7, tmp_path)
    old = vertify.scan_ac_reference(network, 7, grid, path, workers=1)
    lower, upper = np.array(grid['axis_lower']), np.array(grid['bounds'])
    lower[0] -= 10.
    upper[-1] += 20.
    submitted = []
    solve = vertify.ac_scan_line
    def checked(args, **kwargs):
        submitted.extend(args[1].tolist())
        return solve(args, **kwargs)
    with patch('vertify.ac_scan_line', side_effect=checked):
        expanded = vertify.scan_ac_reference(network, 7,
            dict(axis_lower=lower, bounds=upper, shape=(160,)*dimension), path, workers=1)
    expected_shape = (3, 4) if dimension == 2 else (3, 2, 4)
    assert expanded['states'].shape == expected_shape
    assert len(submitted) == np.prod(expected_shape)-2**dimension
    assert expanded['reused_points'] == 2**dimension
    assert not any(all(abs(value) == 5 for value in point) for point in submitted)
    slices = (slice(1, 3), slice(0, 2)) if dimension == 2 else (slice(1, 3), slice(None), slice(0, 2))
    np.testing.assert_array_equal(expanded['states'][slices], old['states'])
    np.testing.assert_array_equal(vertify.load_scan(expanded['cache_path'])['states'], expanded['states'])
    np.testing.assert_array_equal(expanded['socp_states'][slices], old['socp_states'])


def test_cache_identity_excludes_artificial_current_but_includes_physics(tmp_path):
    """旧不限流口径已退役：现在限额必须参与缓存身份。"""
    network = Case33()
    baseline = vertify.scan_path(network, 7, tmp_path)
    assert baseline != vertify.scan_path(Case33(current_limit=140.), 7, tmp_path)
    assert baseline != vertify.scan_path(Case33(current_limit=400.), 7, tmp_path)
    assert vertify.scan_path(Case33(current_limit=200.), 7, tmp_path) != vertify.scan_path(
        Case33(current_limit=[200.]*36+[190.]), 7, tmp_path)
    variants = []
    for name in ('vmin', 'original_p', 'source_qmax'):
        changed = deepcopy(network)
        setattr(changed, name, getattr(changed, name)*.9)
        variants.append(vertify.scan_path(changed, 7, tmp_path))
    variants += [vertify.scan_path(network, 6, tmp_path), vertify.scan_path(network, 7, tmp_path, mode=0),
                 vertify.scan_path(Case33(load_nodes=(25, 18)), 7, tmp_path)]
    assert all(path != baseline for path in variants)


def test_nested_refinement_preserves_old_power_points(tmp_path):
    network = Case33()
    grid = dict(axis_lower=[-10., -10.], bounds=[10., 10.], shape=(2, 2))
    path = tmp_path/'ac.npz'
    original = vertify.scan_ac_reference(network, 7, grid, path, workers=1)
    refined = vertify.scan_ac_reference(network, 7, dict(grid, refinement=2), path, workers=1)
    assert refined['reused_points'] == 4 and refined['computed_points'] == 21
    np.testing.assert_array_equal(refined['step'], original['step']/2)
    np.testing.assert_array_equal(refined['states'][1::2, 1::2], original['states'])


def test_independent_v2_import_preserves_points_and_rejects_socp(tmp_path):
    """旧 v2 不再迁入；同条件的新区域同时导入 AC、SOCP。"""
    network = Case33(current_limit=200.)
    source = tmp_path/'old_ac.npz'
    grid = dict(axis_lower=np.array([-1., -1.]), bounds=np.array([1., 1.]), states=np.ones((2, 2), np.int8))
    np.savez_compressed(source, **grid, method='ac_independent_v2', use_socp_exclusions=False,
        network=network.name, network_fingerprint=network.fingerprint, load_nodes=network.load_nodes,
        budget=7, ell_limit=network.ell_limit,
        power_factors=[vertify.LOAD_PF, vertify.PV_PF, vertify.PV_Q_SIGN])
    target = tmp_path/'paired'
    with pytest.raises(vertify.ACReferenceMismatch):
        vertify.import_ac_reference(network, 7, source, target)
    paired = vertify.scan_ac_reference(network, 7,
        dict(axis_lower=grid['axis_lower'], bounds=grid['bounds'], shape=(2, 2)), tmp_path/'source', workers=1)
    assert vertify.import_ac_reference(network, 7, paired['cache_path'], target) == 8
    with patch('vertify.ac_scan_line', side_effect=AssertionError('must reuse imported points')):
        answer = vertify.scan_ac_reference(network, 7,
            dict(axis_lower=grid['axis_lower'], bounds=grid['bounds'], shape=(2, 2)), target, workers=1)
    assert answer['computed_points'] == 0 and np.isfinite(answer['residual']).all()
    np.testing.assert_array_equal(answer['socp_states'], paired['socp_states'])
    with pytest.raises(vertify.ACReferenceMismatch):
        vertify.import_ac_reference(Case33(), 7, paired['cache_path'], target)


def test_force_rescan_preserves_backup_and_submits_all_points(tmp_path):
    network = Case33()
    grid = dict(axis_lower=[-1., -1.], bounds=[1., 1.], shape=(2, 2))
    path = tmp_path/'ac.npz'
    vertify.scan_ac_reference(network, 7, grid, path, workers=1)
    answer = vertify.scan_ac_reference(network, 7, grid, path, workers=1, force_rescan=True)
    assert answer['reused_points'] == 0 and answer['computed_points'] == 4
    assert len(list(path.glob('before_*.npz'))) == 1
    assert answer['ac_computed_points'] == answer['socp_computed_points'] == 4


def test_wrong_physical_identity_cannot_overwrite_cache(tmp_path):
    network = Case33()
    grid = dict(axis_lower=[-1., -1.], bounds=[1., 1.], shape=(1, 1))
    path = tmp_path/'ac.npz'
    answer = vertify.scan_ac_reference(network, 7, grid, path, workers=1)
    cached = next(path.glob('region_*.npz'))
    original = cached.read_bytes()
    with pytest.raises(ValueError, match='physical identity'):
        vertify.scan_ac_reference(network, 0, grid, path, workers=1)
    assert cached.read_bytes() == original


def test_imported_block_outside_current_grid_is_reused_after_expansion(tmp_path):
    network = Case33()
    path, block = tmp_path/'current.npz', tmp_path/'block.npz'
    first = dict(axis_lower=[-10., -10.], bounds=[10., 10.], shape=(2, 2))
    later = dict(axis_lower=[10., -10.], bounds=[20., 10.], shape=(1, 2))
    vertify.scan_ac_reference(network, 7, first, path, workers=1)
    other = vertify.scan_ac_reference(network, 7, later, block, workers=1)
    assert vertify.import_ac_reference(network, 7, other['cache_path'], path) == 4
    with patch('vertify.ac_scan_line', side_effect=AssertionError('must reuse both blocks')):
        answer = vertify.scan_ac_reference(network, 7,
            dict(axis_lower=[-10., -10.], bounds=[20., 10.], shape=(160, 160)), path, workers=1)
    assert answer['computed_points'] == 0 and answer['reused_points'] == 6
    assert vertify.load_scan(answer['cache_path'])['states'].shape == (3, 2)


def test_conflicting_certificates_are_not_silently_imported(tmp_path):
    network = Case33()
    path = tmp_path/'ac.npz'
    grid = dict(axis_lower=[-1., -1.], bounds=[1., 1.], shape=(1, 1))
    answer = vertify.scan_ac_reference(network, 7, grid, path, workers=1)
    cached = next(path.glob('region_*.npz'))
    original = cached.read_bytes()
    answer['states'] *= -1
    block = tmp_path/'conflict.npz'
    vertify.save_scan(block, answer)
    with pytest.raises(ValueError, match='Conflicting scan certificates'):
        vertify.import_ac_reference(network, 7, block, path)
    assert cached.read_bytes() == original


def test_repeated_refinement_in_load_mode_preserves_nonnegative_coordinates(tmp_path):
    network = Case33()
    path = tmp_path/'load.npz'
    grid = dict(axis_lower=[0., 0.], bounds=[20., 20.], shape=(2, 2))
    first = vertify.scan_ac_reference(network, 7, grid, path, workers=1, mode=0)
    second = vertify.scan_ac_reference(network, 7, dict(grid, refinement=2), path, workers=1, mode=0)
    submitted = []
    solve = vertify.ac_scan_line
    def checked(args, **kwargs):
        submitted.extend(args[1].tolist())
        return solve(args, **kwargs)
    with patch('vertify.ac_scan_line', side_effect=checked):
        third = vertify.scan_ac_reference(network, 7, dict(grid, refinement=2), path, workers=1, mode=0)
    assert first['states'].size == second['reused_points']
    assert second['states'].size == third['reused_points']
    assert np.min(submitted) >= 0
