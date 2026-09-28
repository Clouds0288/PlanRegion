"""双窗口共用播放索引；校验面板和未完成证书不随回放伪造。"""
from monitor import SynchronizedReplay
from tests.test_monitor import make_monitor


def test_two_windows_seek_play_and_hold_shorter_recording(tmp_path):
    paths = []
    for index in range(2):
        monitor, region, x = make_monitor(tmp_path/f'{index}.json.gz')
        for number in range(index+1):
            monitor.sp_start(x, [number+1., 2.], number+1)
        monitor.stopped('time_limit')
        assert monitor.state['coverage_complete'] is False
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
        assert [window.index for window in replay.windows] == [2, 3]
        assert all(window.frame_text.get().endswith(str(len(window.monitor.history)))
                   for window in replay.windows)
        replay.windows[0].seek(0)
        assert [window.index for window in replay.windows] == [0, 0]
    finally:
        replay.close()
