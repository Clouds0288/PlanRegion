"""原生监视器：增量回放、分区转发、子进程通道、暂停与取消、指标分母及锥 / 网架 / 割的显示。"""
import queue
import threading
from threading import Thread
from time import perf_counter, sleep

import numpy as np
import pytest

import monitor as monitor_module
from geometry import contains, covered, halfspaces, polygon_union
from monitor import NativeWindow, RunMonitor, grid_comparison, signed_values
from plot import cap
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
    assert state['schemes']['--:1']['x'] == x.tolist() and state['schemes']['--:1']['status'] == 'stagnated'
    assert state['cut_history']['--:1']['sign'] == [-1, -1] and state['cut_history']['--:1']['scheme'] == '--:1'
    mapped = np.asarray(state['cut_history']['--:1']['cut'])
    power = np.array([40., 60.])
    assert cut[0]+cut[1:3]@power == pytest.approx(mapped[0]+mapped[1:3]@(-power))
    step = monitor.frame(2)['step']   # 过程点在 step 中：坐标乘 sign，网架带分区前缀
    assert step['p'] == [-30., -30.] and step['scheme'] == '--:1'
    assert signed_values(dict(inner=[]), np.array([-1, 1]), '-+:') == dict(inner=[])


def test_child_monitor_sends_patches_and_follows_shared_controls():
    channel = monitor_module.Channel(queue.Queue(), threading.Event(), threading.Semaphore(0), threading.Event(), None)
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
    # 2x2 中，真值有 3 个点，算法覆盖左列 2 个点：重叠 1、多余 1、遗漏 2。
    reference = dict(axis_lower=np.zeros(2), bounds=np.array([2., 2.]), states=np.array([[1, -1], [1, 1]]),
                     socp_states=np.ones((2, 2)))
    result = dict(inner=[dict(vertices=np.array([[0., 0.], [1., 0.], [1., 2.], [0., 2.]]))], outer=[])
    metrics = grid_comparison(reference, result)['metrics']['inner']
    assert metrics['mr_percent'] == pytest.approx(200/3)
    assert metrics['fr_percent'] == 50.


def test_empty_metrics_do_not_claim_zero_error():
    metrics = grid_comparison(dict(axis_lower=np.zeros(2), bounds=np.ones(2), states=-np.ones((2, 2)),
                                   socp_states=-np.ones((2, 2))), dict(inner=[], outer=[]))['metrics']['inner']
    assert metrics['mr_percent'] is None and metrics['fr_percent'] is None


def test_union_keeps_disconnected_domains_and_holes():
    from shapely.geometry import Point
    union = polygon_union([[[0, 0], [1, 0], [1, 1], [0, 1]], [[2, 0], [3, 0], [3, 1], [2, 1]]])
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
    assert len(cap(simplex)) == 3
    ring = cap(clipped)[:, :2]
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
        assert '-100.0, -100.0' in window.step_text.get()   # 取割顶点 bounds 乘 sign
        assert any(item.get_gid() == 'step-point' for item in scheme_axis.collections)
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


def test_main_view_keeps_axes_and_steps_through_one_partition():
    """主图坐标取最新帧 K^OUT 的范围、逐帧不变，网架面板取其中本分区的卦限；步骤栏、上一步/下一步与播放按所选分区走，
    当前步骤画在主图上（坐标外的点贴边并注明）。"""
    monitor, _, _ = recorded_monitor(2)
    steps = [i for i, item in enumerate(monitor.history) if item['patch'].get('step')]
    window = NativeWindow(monitor)
    try:
        window.root.withdraw()
        limits = set()
        for index in range(len(monitor.history)):
            window.seek(index)
            limits.add((window.axes['A'].get_xlim(), window.axes['A'].get_ylim()))
        assert len(limits) == 1
        (x0, x1), (y0, y1) = limits.pop()   # '--' 的锥外块到 -60 kW，'++' 尚无锥取分区盒到 100 kW，两侧各留 10%
        assert (x0, x1, y0, y1) == pytest.approx((-76., 116., -76., 116.))
        assert window.scheme_views['--:1'][1].get_xlim() == pytest.approx((-76., .04*192))
        window.chosen.set('--')
        window.seek(0)
        window.step(1)
        assert window.index == steps[0] and '分区 --' in window.step_text.get()
        assert {line.get_gid() for line in window.axes['A'].lines} >= {'step-line'}
        window.step(1)
        window.step(1)
        assert window.index == steps[2] and window.monitor.history[steps[2]]['patch']['event'] == 'point'
        assert any('min η（坐标外）' == text.get_text() for text in window.axes['A'].texts)   # 评分顶点 (-100,-100) 在坐标外
        assert window.log.curselection() == (window.log_rows.index(steps[2]),)
        assert any(patch.get_gid() == 'step-network' for patch in window.axes['A'].patches)
        window.seek_cut(1)
        assert window.index == steps[3] and '取割 #1' in window.step_text.get()
        window.seek(steps[-1])
        assert any(patch.get_gid() == 'step-cone' for patch in window.axes['A'].patches)
        window.log.selection_set(0)
        window.choose_step(None)
        assert window.index == window.log_rows[0]
        window.chosen.set('++')
        window.redraw()
        assert window.log.size() == 0 and window.step_text.get() == '步骤：—'
        window.seek(0)
        window.play()
        window.tick()
        assert window.index == 0 and not window.playing   # '++' 没有步骤，播放不动
        window.full_box.set(True)
        window.refit()
        assert window.axes['A'].get_xlim() == pytest.approx((-120., 120.))
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
