"""原生监视器：不可变增量回放、暂停、取消、指标分母及动态网架。"""
from pathlib import Path
from threading import Thread
from time import perf_counter, sleep
from unittest.mock import patch

import numpy as np
import pytest

from monitor import RunMonitor, NativeWindow, _union
from region import RegionState
from Network.four_bus_five_corridor import FourBus


def make_monitor(output=None):
    network = FourBus(load_nodes=(1, 2))
    region = RegionState([100., 100.], 150., .002)
    x = network.encode_plan(network.initial_plan).astype(int)
    region.add_scheme(x, network.initial_plan, 0.)
    monitor = RunMonitor(output=output)
    monitor.begin(network, 'socp', 20000., region, 30.)
    return monitor, region, x


def test_delta_history_is_immutable_and_round_trips(tmp_path):
    monitor, region, x = make_monitor(tmp_path/'monitor.json.gz')
    first = monitor.frame(0)
    region.add_point(x, [[0., 0.], [.7, 0.], [0., .4]])
    monitor.updated(region, 'feasible', x=x, power=[70., 0.])
    monitor.sp_start(x, [20., 40.], 1)
    assert 'schemes' not in monitor.history[-1]['patch']
    assert monitor.frame(0) == first
    assert first['schemes']['A']['inner'] == []
    monitor.save()
    restored = RunMonitor()
    restored.load_recording(monitor.output)
    assert restored.history == monitor.history
    assert restored.state == monitor.state
    assert [p.name for p in tmp_path.iterdir()] == ['monitor.json.gz']


def test_global_and_sp_points_remain_distinct():
    monitor, region, x = make_monitor()
    monitor.global_start(32)
    monitor.global_end(dict(complete=False, x=x, p=[20., 30.], bound=.1), region)
    monitor.sp_start(x, [45., 5.], 33)
    assert monitor.state['global_point']['p'] == [20., 30.]
    assert monitor.state['sp_point']['p'] == [45., 5.]
    assert monitor.state['sp'] == 33
    monitor.sp_start(x, [10., 5.], 34)
    monitor.updated(region, 'feasible', x=x, power=[45., 5.],
                    checked=dict(eta=0., feasible=True))
    assert monitor.state['sp'] == 34
    assert monitor.state['sp_point']['number'] == 33  # 缓存证书不冒充最新 SP 调用。


def test_pause_next_continue_do_not_spend_solver_deadline():
    monitor, _, _ = make_monitor()
    monitor.control('pause')
    worker = Thread(target=lambda: [monitor.selecting(2, False), monitor.selecting(1, True)])
    worker.start()
    end = perf_counter()+2.
    while len(monitor.history) < 2 and perf_counter() < end:
        sleep(.005)
    assert worker.is_alive()
    sleep(.05)
    monitor.control('next')
    end = perf_counter()+2.
    while len(monitor.history) < 3 and perf_counter() < end:
        sleep(.005)
    assert worker.is_alive()
    monitor.control('continue')
    worker.join(2.)
    assert not worker.is_alive()
    assert monitor.paused_seconds >= .05


def test_failed_run_saves_history_and_reraises(tmp_path):
    monitor, _, _ = make_monitor(tmp_path/'monitor.json.gz')
    def fail():
        raise RuntimeError('solver probe')
    with pytest.raises(RuntimeError, match='solver probe'):
        monitor.execute(fail, show_ui=False)
    restored = RunMonitor()
    restored.load_recording(monitor.output)
    assert restored.state['status'] == 'failed'
    assert 'result' not in restored.state
    assert restored.state['error'] == 'solver probe'


def test_metrics_use_their_declared_denominators():
    monitor = RunMonitor()
    # 2x2 中，真值有 3 个点，算法覆盖左列 2 个点：重叠 1、多余 1、遗漏 2。
    reference = dict(bounds=np.array([2., 2.]), states=np.array([[1, -1], [1, 1]]))
    result = dict(inner=[dict(vertices=np.array([[0., 0.], [1., 0.], [1., 2.], [0., 2.]]))])
    monitor.validation(reference, result)
    metrics = monitor.state['validation']
    assert metrics['mr_percent'] == pytest.approx(200/3)
    assert metrics['fr_percent'] == 50.


def test_empty_metrics_do_not_claim_zero_error():
    monitor = RunMonitor()
    monitor.validation(dict(bounds=np.ones(2), states=-np.ones((2, 2))), dict(inner=[]))
    assert monitor.state['validation']['mr_percent'] is None
    assert monitor.state['validation']['fr_percent'] is None


def test_union_keeps_disconnected_domains_and_holes():
    from shapely.geometry import Point
    union = _union([[[0, 0], [1, 0], [1, 1], [0, 1]], [[2, 0], [3, 0], [3, 1], [2, 1]]])
    assert union.area == 2.
    assert not union.covers(Point(1.5, .5))


def test_native_window_adds_schemes_and_rewinds_without_solver():
    monitor, region, x = make_monitor()
    window = NativeWindow(monitor)
    try:
        window.root.withdraw()
        window.seek(0)
        assert set(window.scheme_views) == {'A'}
        y = x.copy()
        y[0], y[1] = 0, 1
        region.add_scheme(y, {'01': 'M'}, 100.)
        monitor.global_start(32)
        monitor.global_end(dict(complete=False, x=y, p=[20., 30.], bound=.1), region)
        monitor.sp_start(y, [30., 25.], 33)
        window.seek(len(monitor.history)-1)
        assert set(window.scheme_views) == {'A', 'B'}
        assert set(window.axes) == {'A', 'C'}
        assert '20.000' in window.global_text.get()
        assert '30.000' in window.sp_text.get()
        window.seek(0)
        assert not window.scheme_views['B'][0].winfo_manager()
        window.play()
        window.tick()
        assert window.index == 1
        window.go_live()
        window.tick()
        assert window.index == len(monitor.history)-1
    finally:
        window.close()


def test_cut_lines_use_each_scheme_and_replay_never_leaks_future_cuts(tmp_path):
    from monitor import _cut_segment
    monitor, region, x = make_monitor(tmp_path/'monitor.json.gz')
    j = np.flatnonzero(x)[0]
    y = x.copy()
    y[j] = 0
    region.add_scheme(y, {'01': None}, 0.)
    monitor.global_end(dict(complete=False, x=y, p=[100., 50.], bound=.1), region)
    before = monitor.frame(len(monitor.history)-1)
    delta = np.zeros(len(x))
    delta[j] = -40.
    cut = np.r_[130., -1., -1., delta]
    monitor.sp_start(x, [100., 50.], 1)
    region.apply_cut(cut)
    monitor.updated(region, 'cut', x=x, power=[100., 50.], checked=dict(cut=cut, eta=1., feasible=False))
    first_cut = len(monitor.history)-1
    for scheme in ('A', 'B'):
        row = monitor.state['schemes'][scheme]
        segment = _cut_segment(cut, row['x'], [100., 100.])
        assert len(segment) == 2
        np.testing.assert_allclose(cut[0]+segment@cut[1:3]+cut[3:]@row['x'], 0., atol=1e-12)
    assert not np.allclose(_cut_segment(cut, x, [100., 100.]), _cut_segment(cut, y, [100., 100.]))
    assert len(_cut_segment(np.r_[1., 0., 0., delta*0.], x, [100., 100.])) == 0
    assert not before.get('cut_history')
    cut2 = np.r_[70., -1., 0., delta*0.]
    region.apply_cut(cut2)
    monitor.updated(region, 'cut', x=x, power=[100., 50.], checked=dict(cut=cut2, eta=.5, feasible=False))
    assert list(monitor.frame(first_cut)['cut_history']) == ['1']
    assert list(monitor.history[-1]['patch']['cut_history']) == ['2']
    monitor.save()
    restored = RunMonitor()
    restored.load_recording(monitor.output)
    assert restored.frame(first_cut) == monitor.frame(first_cut)
    window = NativeWindow(restored)
    try:
        window.root.withdraw()
        window.seek(first_cut)
        assert window._removed(restored.frame(first_cut), 'A').area > 0.
        lines = window.scheme_views['A'][1].lines
        assert any(line.get_gid() == 'cut-1-A' for line in lines)
        assert not any(line.get_gid() == 'cut-2-A' for line in lines)
        window.seek_cut(1)
        assert window.index == first_cut+1
        window.seek_cut(-1)
        assert window.index == first_cut
    finally:
        window.close()


def test_finish_requires_explicit_coverage_certificate():
    monitor, region, _ = make_monitor()
    before = len(monitor.history)
    with pytest.raises(KeyError, match='certified'):
        monitor.finish(dict(coverage_bound=0., **region.finish(True)), region)
    assert len(monitor.history) == before
    assert 'result' not in monitor.state


def test_scan_updates_outside_timeline_and_final_comparison_survives_rewind(tmp_path):
    monitor, region, x = make_monitor(tmp_path/'monitor.json.gz')
    region.add_point(x, [[0., 0.], [.7, 0.], [0., .4]])
    result = dict(certified=True, coverage_bound=0., **region.finish(True))
    monitor.finish(result, region)
    total = len(monitor.history)
    window = NativeWindow(monitor)
    try:
        window.root.withdraw()
        window.seek(0)
        monitor.scanning(5, 10)
        window.show()  # 回放位置没变，校验进度仍刷新。
        assert any('5 / 10' in text.get_text() for text in window.axes['C'].texts)
        monitor.validation(dict(bounds=np.array([100., 100.]), states=np.ones((2, 2))), result)
        window.show()
        title = window.axes['C'].get_title()
        assert '遗漏' in title and '多余' in title
        window.seek(total-1)
        window.seek(0)
        assert window.axes['C'].get_title() == title
        assert len(monitor.history) == total
        assert 'validation' not in monitor.frame(total-1)
        monitor.save()
        restored = RunMonitor()
        restored.load_recording(monitor.output)
        assert restored.state['validation'] == monitor.state['validation']
        assert restored.validation_state == monitor.validation_state
        assert restored.history == monitor.history
    finally:
        window.close()


def test_old_recording_separates_scan_frames_without_inventing_cuts(tmp_path):
    import gzip
    import json
    monitor, _, _ = make_monitor()
    old_history = monitor.history+[
        dict(elapsed=1., patch=dict(event='scan', phase='SOCP 扫描', scan_progress=[0, 4])),
        dict(elapsed=2., patch=dict(event='completed', phase='完成', status='completed',
                                  validation=dict(mr_percent=0., fr_percent=0.)))]
    path = tmp_path/'old.json.gz'
    with gzip.open(path, 'wt', encoding='utf-8') as stream:
        json.dump(dict(version=3, history=old_history), stream)
    restored = RunMonitor()
    restored.load_recording(path)
    assert restored.history == monitor.history
    assert restored.state['status'] == 'completed'
    assert restored.state['validation']['mr_percent'] == 0.
    assert not restored.state.get('cut_history')


def test_case33_actual_socp_cut_matches_displayed_line():
    from Network.case33bw import Case33
    from model import GridPhysics, SubProblem
    from monitor import _cut_segment
    network = Case33(load_nodes=(18, 25))
    x = network.encode_plan(network.initial_plan)
    oracle = SubProblem(GridPhysics(network, 'socp'), threads=1)
    assert oracle.solve(x, np.array([90., 420.]))['feasible']
    power = np.array([5000., 5000.])
    answer = oracle.solve(x, power)
    assert not answer['feasible']
    cut = answer['cut']
    assert cut[0]+cut[1:3]@power+cut[3:]@x < 0.
    segment = _cut_segment(cut, x, [network.power_limit]*2)
    assert segment.shape == (2, 2)
    np.testing.assert_allclose(cut[0]+segment@cut[1:3]+cut[3:]@x, 0., atol=1e-10)


def test_many_schemes_reuse_four_panels_and_follow_the_selected_point():
    monitor, region, x = make_monitor()
    for value in range(1, 12):
        y = np.array([(value >> j) & 1 for j in range(len(x))])
        region.add_scheme(y, {'01': 'M'}, float(value))
    monitor.global_end(dict(complete=False, x=y, p=[30., 20.], bound=.1), region)
    window = NativeWindow(monitor)
    try:
        window.root.withdraw()
        window.seek(len(monitor.history)-1)
        focus = monitor.state['global_point']['scheme']
        assert focus in window.scheme_views
        assert len(window.scheme_views) <= 4
        window.page_number.set(1)
        window.choose_page()
        assert window.scheme_page == 0
        assert len(window.scheme_views) == 4
        frames = {id(frame) for frame, _, _ in window.scheme_views.values()}
        window.change_page(1)
        assert window.scheme_page == 1
        assert {id(frame) for frame, _, _ in window.scheme_views.values()} == frames
        window.seek(0)
        assert window.scheme_page == 0
        assert sum(bool(frame.winfo_manager()) for frame, _, _ in window.scheme_views.values()) == 1
    finally:
        window.close()


def test_cut_excluding_the_whole_box_does_not_invent_a_boundary_line():
    monitor, region, x = make_monitor()
    cut = np.r_[-1., np.zeros(2+len(x))]
    monitor.sp_start(x, [50., 50.], 1)
    region.apply_cut(cut)
    monitor.updated(region, 'cut', x=x, power=[50., 50.], checked=dict(cut=cut, eta=1., feasible=False))
    window = NativeWindow(monitor)
    try:
        window.root.withdraw()
        window.seek(len(monitor.history)-1)
        ax = window.scheme_views['A'][1]
        assert any('本框全部排除' in text.get_text() for text in ax.texts)
        assert not any(line.get_gid() == 'cut-1-A' for line in ax.lines)
    finally:
        window.close()


@pytest.mark.parametrize('dimension', [2, 3])
def test_native_view_fits_current_domains_cuts_and_rays_when_global_box_is_loose(dimension):
    from itertools import product
    from copy import deepcopy
    cube = np.array(list(product((0., 1.), repeat=dimension)))
    positive = 10.+10.*cube
    monitor = RunMonitor()
    row = dict(x=[], choice={}, cost=0., sign=[1]*dimension,
               inner=(10.+5.*cube).tolist(), outer=positive.tolist())
    monitor._emit('initial_bounds', partition='+'*dimension,
                  load_nodes=list(range(1, dimension+1)), bounds=[1e6]*dimension,
                  axis_lower=[-1e6]*dimension, axis_bounds=[100.]*dimension,
                  total_bound=dimension*1e6, global_outer=[(100.*cube).tolist(), (-1e6*cube).tolist()],
                  schemes={'A': row}, cut_history={})
    first = deepcopy(monitor.frame(0))
    window = NativeWindow(monitor)
    try:
        window.root.withdraw()
        window.seek(0)
        axis = window.scheme_views['A'][1]
        initial_limits = (-6., 106.)
        np.testing.assert_allclose(axis.get_xlim(), initial_limits)
        np.testing.assert_allclose(window.axes['A'].get_xlim(), initial_limits)
        assert window.overview.get_xlim()[0] < -1e6
        assert any(a.get_gid() == 'global-outer' for a in window.axes['A'].get_children())
        assert any(a.get_gid() == 'overview-outer' for a in window.overview.get_children())
        if dimension == 3:
            window.axes['A'].view_init(31., 70.)
            np.testing.assert_allclose(axis.get_zlim(), initial_limits)
        smaller = positive.copy()
        smaller[:, 0] = np.minimum(smaller[:, 0], 15.)
        monitor._emit('cut', schemes={'A': {**row, 'outer': smaller.tolist()}},
                      cut_history={'1': dict(cut=[15., -1.]+[0.]*(dimension-1),
                                              scheme='A', sign=[1]*dimension)})
        window.seek(1)
        axis = window.scheme_views['A'][1]
        np.testing.assert_allclose(axis.get_xlim(), initial_limits)
        artists = axis.collections if dimension == 3 else axis.lines
        assert any(item.get_gid() == 'cut-1-A' for item in artists)
        monitor._emit('selection')
        window.seek(2)
        np.testing.assert_allclose(window.scheme_views['A'][1].get_xlim(), initial_limits)
        monitor._emit('ray_end', phase='射线补点', active_scheme='A', stage=1, ray_fraction=.5,
                      ray=dict(scheme='A', anchor=[10.]*dimension,
                               target=[80.]*dimension, p=[45.]*dimension))
        window.seek(3)
        np.testing.assert_allclose(window.scheme_views['A'][1].get_xlim(), initial_limits)
        assert {'ray-path', 'ray-anchor', 'ray-target', 'ray-point'} <= {
            a.get_gid() for a in window.axes['A'].get_children()}
        window.full_extent.set(True)
        window.refresh_extent()
        assert window.axes['A'].get_xlim()[0] < -1e6
        assert not window.overview.get_visible()
        # 展开全局只切换 A，小图仍保持分区的固定坐标。
        np.testing.assert_allclose(window.scheme_views['A'][1].get_xlim(), initial_limits)
        window.full_extent.set(False)
        window.refresh_extent()
        monitor._emit('initial_bounds', partition='-'*dimension, phase='初始化', active_scheme=None, ray=None,
                      global_outer=[(100.*cube).tolist(), (-200.*cube).tolist()],
                      schemes={'B': {**row, 'sign': [-1]*dimension,
                                     'inner': (-10.*cube).tolist(), 'outer': (-200.*cube).tolist()}})
        window.seek(4)
        np.testing.assert_allclose(window.axes['A'].get_xlim(), [-212., 12.])
        np.testing.assert_allclose(window.scheme_views['A'][1].get_xlim(), initial_limits)
        window.seek(0)
        np.testing.assert_allclose(window.axes['A'].get_xlim(), initial_limits)
        assert window.overview.get_xlim()[0] < -1e6
        if dimension == 3:
            assert window.axes['A'].elev == 31. and window.axes['A'].azim == 70.
        assert monitor.frame(0) == first
        for canvas in window.canvases.values():
            canvas.draw()
    finally:
        window.close()


def test_global_rejections_candidates_and_ray_process_are_distinct_and_rewindable():
    monitor, region, x = make_monitor()
    monitor._emit('residual_end', phase='线性外域候选',
                  global_point=dict(scheme='A', p=[80., 60.]))
    candidate_index = len(monitor.history)-1
    monitor.sp_start(x, [80., 60.], 1)
    monitor.sp_end(dict(eta=.2, feasible=False))
    rejected_index = len(monitor.history)-1
    monitor._emit('certification_end', phase='全网架认证', query_point=[60., 70.],
                  rejected_points=[dict(p=[20., 80.])], unknown_points=[[90., 10.]])
    proof_index = len(monitor.history)-1
    monitor._emit('ray_start', phase='首轮射线认证', active_scheme='A', stage=1,
                  ray=dict(scheme='A', anchor=[10., 10.], target=[80., 60.], p=None))
    ray_index = len(monitor.history)-1
    window = NativeWindow(monitor)
    try:
        window.root.withdraw()
        window.seek(rejected_index)
        ax = window.axes['A']
        assert 'candidate-rejected' in {a.get_gid() for a in ax.collections}
        assert 'global-rejected' not in {a.get_gid() for a in ax.collections}
        limits = ax.get_xlim(), ax.get_ylim()
        window.seek(proof_index)
        assert {'global-rejected', 'candidate-rejected', 'unknown-points', 'query-point'} <= {
            a.get_gid() for a in ax.collections}
        labels = {text.get_text() for text in ax.get_legend().get_texts()}
        assert {'全网架不可行', '全局候选：该网架不可行'} <= labels
        window.seek(ray_index)
        for panel in (ax, window.scheme_views['A'][1]):
            assert {'ray-path', 'ray-anchor', 'ray-target'} <= {a.get_gid() for a in panel.get_children()}
            assert 'ray-point' not in {a.get_gid() for a in panel.get_children()}
        np.testing.assert_allclose((ax.get_xlim(), ax.get_ylim()), limits)
        window.seek(candidate_index)
        assert not window.failed_candidates
        assert 'candidate-rejected' not in {a.get_gid() for a in ax.collections}
        assert 'global-rejected' not in {a.get_gid() for a in ax.collections}
        assert 'ray-path' not in {a.get_gid() for a in ax.lines}
        assert 'global_point' in {a.get_gid() for a in ax.collections}
        window.canvases['A'].draw()
    finally:
        window.close()


def test_validation_view_includes_feasible_cell_edges_outside_computed_region():
    monitor, region, _ = make_monitor()
    result = dict(certified=False, coverage_bound=1., **region.finish(False))
    monitor.finish(result, region)
    # 计算域在 [0, 100]，参考扫描在其左侧还有一格可行域。
    states = -np.ones((20, 20), dtype=int)
    states[4, 5] = 1
    monitor.validation(dict(axis_lower=[-100., -100.], bounds=[100., 100.], states=states), result)
    window = NativeWindow(monitor)
    try:
        window.root.withdraw()
        window.seek(0)
        ax = window.axes['C']
        assert -71. < ax.get_xlim()[0] < -60.
        assert -61. < ax.get_ylim()[0] < -50.
        assert ax.get_xlim()[1] > 100.
        assert ax.get_ylim()[1] > 100.
        window.canvases['C'].draw()
    finally:
        window.close()
