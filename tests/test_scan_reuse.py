"""扫描缓存复用，以及主线三步（构域、扫描校验、收敛过程图）的端到端核对。"""
import gzip
import json
from unittest.mock import patch

import numpy as np
import pytest

import main
from monitor import RunMonitor, NativeWindow, grid_comparison
from vertify import convergence


def test_real_scan_is_reused_and_metrics_are_recomputed_for_new_result(tmp_path, monkeypatch):
    for name, value in dict(WORKERS=4, SCAN_WORKERS=1, SCAN_OUTPUT=tmp_path/'scans', SCAN_DIVISIONS={2: 8}).items():
        monkeypatch.setattr(main, name, value)
    main.main('fourbus', 2, show=False, output=tmp_path/'first')
    with gzip.open(tmp_path/'first'/'mode_1'/'fourbus_1_2.json.gz', 'rt', encoding='utf-8') as stream:
        first = json.load(stream)
    monkeypatch.setattr(main, 'REGION_TAU', .01)
    result = main.main('fourbus', 2, show=False, output=tmp_path/'second')
    recording = tmp_path/'second'/'mode_1'/'fourbus_1_2.json.gz'
    with gzip.open(recording, 'rt', encoding='utf-8') as stream:
        second = json.load(stream)
    a = first['validation_state']['validation']
    b = second['validation_state']['validation']
    assert (a['bounds'], a['axis_lower'], a['states']) == (b['bounds'], b['axis_lower'], b['states'])
    assert b['computed_points'] == 0 and b['reused_points'] == np.size(b['states'])   # 第二次运行全部复用缓存
    assert b['metrics'] == grid_comparison(a, result)['metrics']
    assert (tmp_path/'second'/'fourbus_1_2_convergence.png').is_file()
    restored = RunMonitor()
    restored.load_recording(recording)
    curves = convergence(restored)   # 只读记录：末点即校验指标
    assert (curves['mr'][-1], curves['fr'][-1]) == pytest.approx((b['mr_percent'], b['fr_percent']))
    state = restored.state
    assert {tuple(row['sign']) for row in state['result']['outer']} == {(1, 1), (1, -1), (-1, 1), (-1, -1)}
    assert min(row[0] for part in state['result']['inner'] for row in part['vertices']) < 0.
    assert sum('result' in event['patch'] for event in restored.history) == 1
    assert len(state['cut_history']) == sum(p['cuts'] for p in result['partitions'])
    assert all(':' in key for key in (*state['schemes'], *state['cones'], *state['cut_history']))
    window = NativeWindow(restored)
    try:
        window.root.withdraw()
        cut_indices = [i for i, event in enumerate(restored.history) if event['patch'].get('event') == 'cut']
        for index in cut_indices[:3]:
            window.seek(index)
            current = restored.frame(index)
            seen = {key for item in restored.history[:index+1] for key in item['patch'].get('cut_history', {})}
            assert set(current['cut_history']) == seen   # 回放不泄露未来的割
            latest = next(iter(restored.history[index]['patch']['cut_history']))
            row, scheme = current['cut_history'][latest], current['active_scheme']
            from plot import cut_segment
            sign = np.asarray(row['sign'])
            ax = window.scheme_views[scheme][1]
            lower, upper = np.asarray([ax.get_xlim(), ax.get_ylim()]).T
            lower = np.where(sign > 0, np.maximum(lower, 0.), lower)
            upper = np.where(sign < 0, np.minimum(upper, 0.), upper)
            segment = cut_segment(row['cut'], current['schemes'][scheme]['x'], upper, lower)
            if len(segment) == 2:
                drawn = next(line for line in ax.lines if line.get_gid() == f'cut-{latest}-{scheme}')
                np.testing.assert_allclose(np.asarray(drawn.get_data()).T, segment)
        if cut_indices:
            window.seek(cut_indices[0])
            window.seek_cut(1)
            assert window.index == (cut_indices[1] if len(cut_indices) > 1 else cut_indices[0])
        window.seek(len(restored.history)-1)
        window.seek(0)
        assert not restored.frame(0).get('cut_history')
        grid = window.axes['C'].collections[0].get_coordinates()
        np.testing.assert_allclose(grid.min(axis=(0, 1)), a['axis_lower'])
        np.testing.assert_allclose(grid.max(axis=(0, 1)), a['bounds'])
    finally:
        window.close()


def test_time_limit_retains_all_partitions_and_saves_the_record_before_validation(tmp_path, monkeypatch):
    monkeypatch.setattr(main, 'WORKERS', 4)
    output = tmp_path/'mode_1'/'fourbus_1_2.json.gz'

    def saved(monitor):
        """校验开始时构域记录已落盘，且尚无校验内容。"""
        with gzip.open(output, 'rt', encoding='utf-8') as stream:
            data = json.load(stream)
        assert data['version'] == 4 and data['validation_state'] == {} and len(data['history']) > 0

    with patch('main.validate', side_effect=saved) as validate, patch('main.convergence'), patch('main.draw_convergence'):
        result = main.main('fourbus', 2, seconds=0., show=False, output=tmp_path)
    validate.assert_called_once()
    assert result['status'] == 'time_limit' and not result['certified']
    assert len(result['outer']) == 4 and result['inner'] == []
    assert [path for path in tmp_path.rglob('*') if path.is_file()] == [output]
