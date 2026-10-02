"""双窗口共用播放索引；校验面板和未完成证书不随回放伪造。"""
from monitor import SynchronizedReplay
from tests.planning_checks import recorded_monitor


def test_two_windows_seek_play_and_hold_shorter_recording(tmp_path):
    paths = []
    for index in range(2):
        monitor, _, _ = recorded_monitor(2, tmp_path/f'{index}.json.gz')
        for number in range(index+1):
            monitor.forward('--', dict(event='point', phase='网架支撑', active_scheme='1',
                                       support=dict(certified=1, faces=4, ratio=None, status='active'),
                                       sp_point=dict(scheme='1', p=[number+1., 2.])))
        monitor._emit('region_end', phase='构域停止', result=dict(status='time_limit', certified=False,
                                                              inner=[], outer=[]), partition=None)
        assert monitor.state['result']['certified'] is False
        monitor.save()
        paths.append(monitor.output)
    replay = SynchronizedReplay(paths)
    try:
        for window in replay.windows:
            window.root.withdraw()
        replay.windows[1].seek(1)
        assert [window.index for window in replay.windows] == [1, 1]
        replay.play()
        replay.tick()
        assert [window.index for window in replay.windows] == [2, 2]
        replay.seek(replay.total-1)
        shorter = len(replay.windows[0].monitor.history)
        assert [window.index for window in replay.windows] == [shorter-1, shorter]
        assert all(window.frame_text.get().endswith(str(len(window.monitor.history)))
                   for window in replay.windows)
        replay.windows[0].seek(0)
        assert [window.index for window in replay.windows] == [0, 0]
    finally:
        replay.close()
