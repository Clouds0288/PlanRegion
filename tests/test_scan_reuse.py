"""扫描按物理范围复用：原坐标保持、扩界、强制重扫与本次结果的回放。"""
from unittest.mock import patch

import numpy as np
import pytest

import main
from monitor import RunMonitor


@pytest.mark.parametrize('dimension', [2, 3])
def test_covering_scan_keeps_original_grid_when_requested_divisions_change(tmp_path, dimension):
    network = main.Case33(load_nodes=(18, 25, 30)[:dimension])
    path = main.scan_path(network, 7, tmp_path)
    bounds = np.array([50., 40., 30.])[:dimension]
    states = np.where(np.indices((3,)*dimension).sum(axis=0) < 3, 1, -1)
    main.save_scan(path, dict(bounds=bounds, states=states, scan_seconds=12.5))
    original = path.read_bytes()

    with patch('main.validate_socp_region', side_effect=AssertionError('覆盖范围内不应重扫')):
        reference = main.scan_reference(network, 7, 100, bounds-.1, path,
                                        threads=4, workers=2, progress=lambda *args: None)
    np.testing.assert_array_equal(reference['bounds'], bounds)
    np.testing.assert_array_equal(reference['states'], states)
    assert set(reference) == {'bounds', 'states'}
    assert path.read_bytes() == original


@pytest.mark.parametrize('dimension', [2, 3])
def test_range_expansion_and_force_rescan_preserve_covered_box(tmp_path, dimension):
    network = main.Case33(load_nodes=(18, 25, 30)[:dimension])
    path = main.scan_path(network, 7, tmp_path)
    bounds = np.array([10., 20., 30.])[:dimension]
    main.save_scan(path, dict(bounds=bounds, states=np.ones((3,)*dimension), scan_seconds=1.))

    def solve(network, budget, divisions, bounds, **kwargs):
        return dict(bounds=bounds, states=-np.ones((divisions,)*dimension, dtype=np.int8), scan_seconds=2.)

    with patch('main.validate_socp_region', side_effect=solve) as scan:
        reference = main.scan_reference(network, 7, 4, np.array([10.2, 18., 31.7])[:dimension], path,
                                        threads=4, workers=2, progress=lambda *args: None)
        expected = np.array([11., 20., 32.])[:dimension]
        np.testing.assert_array_equal(scan.call_args.args[3], expected)
        assert reference['states'].shape == (4,)*dimension
        reference = main.scan_reference(network, 7, 6, bounds-1., path, force_rescan=True,
                                        threads=4, workers=2, progress=lambda *args: None)
        assert scan.call_count == 2
        np.testing.assert_array_equal(scan.call_args.args[3], expected)
        assert reference['states'].shape == (6,)*dimension
    with np.load(path) as saved:
        assert set(saved.files) == {'bounds', 'states'}
        np.testing.assert_array_equal(saved['bounds'], expected)
        np.testing.assert_array_equal(saved['states'], reference['states'])


@pytest.mark.parametrize('network,budget', [
    (main.Case33(load_nodes=(18, 25)), 6),
    (main.Case33(load_nodes=(25, 18)), 7),
    (main.Case33(load_nodes=(18, 25, 30)), 7),
    (main.FourBus(load_nodes=(1, 2)), 20000.),
])
def test_different_network_nodes_or_budget_do_not_share_reference(tmp_path, network, budget):
    original = main.scan_path(main.Case33(), 7, tmp_path)
    main.save_scan(original, dict(bounds=[5000., 5000.], states=np.ones((3, 3)), scan_seconds=1.))
    d = len(network.load_nodes)
    reference = dict(bounds=np.full(d, 10.), states=np.ones((2,)*d), scan_seconds=1.)
    with patch('main.validate_socp_region', return_value=reference) as scan:
        main.scan_reference(network, budget, 2, np.full(d, 10.), main.scan_path(network, budget, tmp_path),
                            threads=1, workers=1, progress=lambda *args: None)
    scan.assert_called_once()


@pytest.mark.parametrize('force', [False, True])
def test_main_passes_force_rescan_switch(force):
    with patch('main.FORCE_RESCAN', force), patch('main.run') as run:
        main.main('case33')
    assert run.call_args.kwargs['force_rescan'] is force


def test_real_scan_is_reused_and_metrics_are_recomputed_for_new_result(tmp_path):
    network = main.FourBus(load_nodes=(1, 2))
    options = dict(budget=20000., divisions=4, show_ui=False, threads=1, scan_workers=1,
                   scan_output=tmp_path/'scans')
    # 1. 强制扫描忽略旧参考入口，保存完整物理扫描与回放。
    main.run(network, output=tmp_path/'first.json.gz', force_rescan=True,
             reference=tmp_path/'unused.json.gz', **options)
    first = RunMonitor()
    first.load_recording(tmp_path/'first.json.gz')
    # 2. 改构域精度后再次运行；已有扫描覆盖，禁止调用扫描求解器。
    with patch('main.validate_socp_region', side_effect=AssertionError('应复用已保存扫描')):
        result = main.run(network, output=tmp_path/'second.json.gz', tau=.02, **options)
    second = RunMonitor()
    second.load_recording(tmp_path/'second.json.gz')
    before, after = first.state['validation'], second.state['validation']
    np.testing.assert_array_equal(after['bounds'], before['bounds'])
    np.testing.assert_array_equal(after['states'], before['states'])
    # 3. 误差必须来自本次几何与原扫描坐标，而非沿用上次误差。
    reference = dict(bounds=np.asarray(before['bounds']), states=np.asarray(before['states']))
    expected = RunMonitor()
    expected.validation(reference, result)
    assert after['metrics'] == expected.state['validation']['metrics']
    assert after['fr_percent'] == 0.
    assert second.state['status'] == 'completed'
