"""正负范围复用、模式隔离，以及新主线记录与扫描的端到端核对。"""
import gzip
import json
from unittest.mock import patch

import numpy as np
import pytest

import main
import vertify
from monitor import RunMonitor, NativeWindow








@pytest.mark.parametrize('force', [False, True])
def test_main_passes_force_rescan_switch(force):
    with patch('main.FORCE_RESCAN', force), patch('main.run') as run:
        main.main('case33')
    assert run.call_args.kwargs['force_rescan'] is force


@pytest.mark.parametrize('relative_cache', [True, False])
def test_archived_reference_resolves_cache_after_changing_working_directory(tmp_path, monkeypatch, relative_cache):
    project = tmp_path/'project'
    project.mkdir()
    other = tmp_path/'elsewhere'
    other.mkdir()
    monkeypatch.chdir(other)
    monkeypatch.setattr(main, 'ROOT', project)
    cache = project/'results/scan/region.npz'
    recording = project/'recording.json.gz'
    stored = cache.relative_to(project) if relative_cache else cache
    with gzip.open(recording, 'wt', encoding='utf-8') as stream:
        json.dump(dict(version=4, history=[dict(elapsed=0., patch=dict(event='start'))],
            validation_state=dict(validation=dict(method=vertify.AC_CACHE_METHOD, cache_path=str(stored)))), stream)
    lower, upper = np.zeros(2), np.ones(2)
    result = dict(inner=[], outer=[], axis_lower=lower, axis_bounds=upper)
    reference = dict(axis_lower=lower, bounds=upper, states=np.ones((1, 1)))
    with patch('main.build_region', return_value=result), \
         patch('main.reference_box', return_value=(lower, upper)), \
         patch('main.scan_ac_reference', return_value=reference), \
         patch('main.import_ac_reference') as imported:
        main.run(main.FourBus(load_nodes=(1, 2)), reference=recording,
                 show_ui=False, scan_output=tmp_path/'cache')
    assert imported.call_args.args[2] == cache


def test_real_scan_is_reused_and_metrics_are_recomputed_for_new_result(tmp_path):
    network = main.FourBus(load_nodes=(1, 2))
    options = dict(budget=20000., divisions=8, show_ui=False, threads=1, scan_workers=1,
                   scan_output=tmp_path/'scans')
    main.run(network, output=tmp_path/'first.json.gz', force_rescan=True,
             reference=tmp_path/'unused.json.gz', **options)
    with gzip.open(tmp_path/'first.json.gz', 'rt', encoding='utf-8') as stream:
        first = json.load(stream)
    with patch('vertify.ProcessPoolExecutor', side_effect=AssertionError('Must reuse saved scan')):
        result = main.run(network, output=tmp_path/'second.json.gz', tau=.01, **options)
    with gzip.open(tmp_path/'second.json.gz', 'rt', encoding='utf-8') as stream:
        second = json.load(stream)
    a = first['validation_state']['validation']
    b = second['validation_state']['validation']
    assert (a['bounds'], a['axis_lower'], a['states']) == (b['bounds'], b['axis_lower'], b['states'])
    check = RunMonitor()
    check.validation(a, result)
    assert b['metrics'] == check.validation_state['validation']['metrics']
    restored = RunMonitor()
    restored.load_recording(tmp_path/'second.json.gz')
    state = restored.state
    assert state['mode'] == 1 and state['result']['certified']
    assert {tuple(row['sign']) for row in state['result']['inner']} == {(1, 1), (1, -1), (-1, 1), (-1, -1)}
    assert min(row[0] for part in state['result']['inner'] for row in part['vertices']) < 0.
    assert sum('result' in event['patch'] for event in restored.history) == 1
    assert len(state['cut_history']) == result['counts']['cuts']
    assert all(':' in key for key in state['schemes'])
    for counter in ('cuts', 'sp', 'global_search'):
        assert state[counter] == result['counts'][counter]
    window = NativeWindow(restored)
    try:
        window.root.withdraw()
        cut_indices = [i for i, event in enumerate(restored.history) if event['patch'].get('event') == 'cut']
        for partition in ('++', '+-', '-+', '--'):
            index = next(i for i in cut_indices if restored.frame(i)['partition'] == partition)
            window.seek(index)
            current = restored.frame(index)
            scheme = current['active_scheme']
            assert all(int(key) <= current['cuts'] for key in current['cut_history'])
            assert window.axes['A'].get_xlim()[0] < 0 < window.axes['A'].get_xlim()[1]
            latest = str(current['cuts'])
            row = current['cut_history'][latest]
            from monitor import _cut_segment
            sign = np.asarray(row['sign'])
            ax = window.scheme_views[scheme][1]
            lower, upper = np.asarray([ax.get_xlim(), ax.get_ylim()]).T
            lower = np.where(sign > 0, np.maximum(lower, 0.), lower)
            upper = np.where(sign < 0, np.minimum(upper, 0.), upper)
            segment = _cut_segment(row['cut'], current['schemes'][scheme]['x'], upper, lower)
            if len(segment) == 2:
                lines = ax.lines
                drawn = next(line for line in lines if line.get_gid() == f'cut-{latest}-{scheme}')
                np.testing.assert_allclose(np.asarray(drawn.get_data()).T, segment)
            window.seek_cut(1)
            assert window.index == min(i for i in cut_indices if i > index)
            window.seek_cut(-1)
            assert window.index == index
        window.seek(len(restored.history)-1)
        window.seek(0)
        assert not restored.frame(0).get('cut_history')
        grid = window.axes['C'].collections[0].get_coordinates()
        np.testing.assert_allclose(grid.min(axis=(0, 1)), a['axis_lower'])
        np.testing.assert_allclose(grid.max(axis=(0, 1)), a['bounds'])
    finally:
        window.close()


def test_positive_mode_is_nonnegative_and_has_independent_replay(tmp_path):
    with patch('main.run') as run:
        main.main('fourbus', mode=0)
    assert run.call_args.kwargs['mode'] == 0
    assert 'mode_0' in str(run.call_args.kwargs['output'])
    assert main.mode == 1
    output = main.recording_path('fourbus', tmp_path, mode=0)
    result = main.run(main.FourBus(load_nodes=(1, 2)), output=output, mode=0,
                      show_ui=False, scan=False, threads=1)
    assert result['certified']
    assert all(np.min(row['vertices']) >= -1e-8 for row in result['inner'])
    restored = RunMonitor()
    restored.load_recording(output)
    assert restored.state['mode'] == 0
    assert np.min(restored.state['axis_lower']) >= 0.
    assert list(output.parent.iterdir()) == [output]


def test_time_limit_retains_all_partitions_and_replay_without_scan(tmp_path):
    output = tmp_path/'timeout.json.gz'
    result = main.run(main.FourBus(load_nodes=(1, 2)), output=output, time_limit=0.,
                      show_ui=False, scan=False, threads=1)
    assert result['status'] == 'time_limit' and not result['certified']
    assert len(result['outer']) == 4 and result['inner'] == []
    with gzip.open(output, 'rt', encoding='utf-8') as stream:
        data = json.load(stream)
    assert data['version'] == 4 and data['validation_state'] == {}
    assert len(data['history']) > 0 and list(tmp_path.iterdir()) == [output]
