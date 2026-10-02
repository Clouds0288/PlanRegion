"""原生监视器：增量回放、分区转发、子进程通道、暂停与取消、指标分母及锥 / 网架 / 割的显示。"""
import queue
import threading
from threading import Thread
from time import perf_counter, sleep

import numpy as np
import pytest

import monitor as monitor_module
from monitor import NativeWindow, RunMonitor, _cap, _union, signed_values
from region import contains, covered, halfspaces
from tests.planning_checks import recorded_monitor


def test_delta_history_is_immutable_and_round_trips(tmp_path):
    monitor, _, _ = recorded_monitor(2, tmp_path/'monitor.json.gz')
    root = monitor.frame(2)
    assert root['cones']['--:0'] is not None and root['partition'] == '--'
    assert monitor.state['cones']['--:0'] is None   # 被细分的父锥置空，子锥按键增量加入
    assert {'--:1', '--:2'} <= set(monitor.state['cones'])
    assert 'cones' not in monitor.history[3]['patch']
    assert monitor.frame(2) == root
    monitor.save()
    restored = RunMonitor()
    restored.load_recording(monitor.output)
    assert restored.history == monitor.history
    assert restored.state == monitor.state
    assert [p.name for p in tmp_path.iterdir()] == ['monitor.json.gz']


def test_forward_prefixes_partition_and_maps_signs():
    monitor, x, cut = recorded_monitor(2)
    state = monitor.state
    assert np.all(np.asarray(state['cones']['--:1']['inner']) <= 0.)
    assert state['schemes']['--:1']['x'] == x.tolist() and state['schemes']['--:1']['status'] == 'eps_B'
    assert state['cut_history']['--:1']['sign'] == [-1, -1] and state['cut_history']['--:1']['scheme'] == '--:1'
    assert state['schemes']['--:1']['sign'] == [-1, -1]
    mapped = np.asarray(state['cut_history']['--:1']['cut'])
    power = np.array([40., 60.])
    assert cut[0]+cut[1:3]@power == pytest.approx(mapped[0]+mapped[1:3]@(-power))
    assert monitor.frame(2)['global_point'] == dict(scheme='--:1', p=[-30., -30.])
    assert signed_values(dict(inner=[]), np.array([-1, 1]), '-+:') == dict(inner=[])


def test_child_monitor_sends_patches_and_follows_shared_controls():
    channel = queue.Queue(), threading.Event(), threading.Semaphore(0), threading.Event()
    monitor_module._connect(channel)
    try:
        child = RunMonitor(sign=(1, -1))
        cones = {'0': dict(inner=[[0., 0.]], outer=[[0., 0.]], scheme='1', mu=1.)}
        child._emit('cone', cones=cones)
        assert channel[0].get_nowait() == ('+-', dict(cones=cones, event='cone'))
        child._emit('cone', cones=cones)
        assert channel[0].get_nowait() == ('+-', dict(event='cone'))   # 未变的键不重复发送
        channel[1].set()
        worker = Thread(target=lambda: child._emit('point', eta=1.))
        worker.start()
        sleep(.3)
        assert worker.is_alive()
        channel[2].release()
        worker.join(2.)
        assert not worker.is_alive() and child.paused_seconds >= .25
        channel[3].set()
        with pytest.raises(KeyboardInterrupt):
            child._emit('point', eta=2.)
        child.close()
        items = [channel[0].get_nowait() for _ in range(channel[0].qsize())]
        assert items[-1] == ('+-', None)
    finally:
        monitor_module._connect(None)


def test_pause_next_continue_do_not_spend_solver_deadline():
    monitor, _, _ = recorded_monitor(2)
    start = len(monitor.history)
    monitor.control('pause')
    worker = Thread(target=lambda: [monitor._emit('point', eta=1.), monitor._emit('point', eta=2.)])
    worker.start()
    end = perf_counter()+2.
    while len(monitor.history) < start+1 and perf_counter() < end:
        sleep(.005)
    assert worker.is_alive()
    sleep(.05)
    monitor.control('next')
    end = perf_counter()+2.
    while len(monitor.history) < start+2 and perf_counter() < end:
        sleep(.005)
    assert worker.is_alive()
    monitor.control('continue')
    worker.join(2.)
    assert not worker.is_alive()
    assert monitor.paused_seconds >= .05


def test_failed_run_saves_history_and_reraises(tmp_path):
    monitor, _, _ = recorded_monitor(2, tmp_path/'monitor.json.gz')

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


def test_covered_matches_row_by_row_containment():
    rng = np.random.default_rng(3)
    points = rng.uniform(-1., 2., (2000, 3))
    rows = [dict(vertices=rng.uniform(0., 1., (6, 3))+shift) for shift in (0., .8, -.5)]
    expected = np.zeros(len(points), dtype=bool)
    for row in rows:
        expected |= contains(points, halfspaces(row['vertices']))
    np.testing.assert_array_equal(covered(points, rows), expected)


def test_cap_keeps_the_far_face_in_cyclic_order():
    simplex = np.array([[0., 0., 0.], [3., 0., 0.], [0., 3., 0.], [0., 0., 3.]])
    clipped = np.array([[0., 0., 0.], [1., 0., 0.], [1., .5, 0.], [.5, 1., 0.], [0., 1., 0.], [0., 0., 1.]])
    assert len(_cap(simplex)) == 3
    ring = _cap(clipped)[:, :2]
    signed = .5*np.sum(ring[:, 0]*np.roll(ring[:, 1], -1)-ring[:, 1]*np.roll(ring[:, 0], -1))
    assert abs(signed) > 0. and len(ring) == 5


def test_native_window_draws_cones_networks_and_rewinds():
    monitor, _, _ = recorded_monitor(2)
    window = NativeWindow(monitor)
    try:
        window.root.withdraw()
        last = len(monitor.history)-1
        window.seek(last)
        assert set(window.scheme_views) == {'--:1'} and set(window.axes) == {'A', 'C'}
        assert '分区 --' in window.status.get() and '叶锥 2' in window.status.get()
        assert window.axes['A'].patches
        cut_frame = next(i for i, item in enumerate(monitor.history) if item['patch']['event'] == 'cut')
        window.seek(cut_frame)
        scheme_axis = window.scheme_views['--:1'][1]
        assert any(line.get_gid() == 'cut---:1---:1' for line in scheme_axis.lines)
        assert '-100.000, -100.000' in window.sp_text.get()
        window.seek(0)
        assert not window.scheme_views['--:1'][0].winfo_manager()
        window.play()
        window.tick()
        assert window.index == 1
        window.go_live()
        window.tick()
        assert window.index == last
    finally:
        window.close()


def test_three_dimensional_window_draws_cone_caps_and_current_cut():
    monitor, _, _ = recorded_monitor(3)
    window = NativeWindow(monitor)
    try:
        window.root.withdraw()
        cut_frame = next(i for i, item in enumerate(monitor.history) if item['patch']['event'] == 'cut')
        window.seek(cut_frame)
        assert window.axes['A'].name == window.axes['C'].name == '3d'
        assert {'cone-outer', 'cone-inner'} <= {item.get_gid() for item in window.axes['A'].collections}
        scheme_axis = window.scheme_views['---:1'][1]
        assert scheme_axis.name == '3d' and any(item.get_gid() == 'cut----:1----:1' for item in scheme_axis.collections)
        window.axes['A'].view_init(31., 70.)
        window.seek(len(monitor.history)-1)
        assert window.axes['A'].elev == 31. and window.axes['A'].azim == 70.
        for canvas in window.canvases.values():
            canvas.draw()
    finally:
        window.close()
