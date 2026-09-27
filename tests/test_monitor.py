"""真实步骤门、人工等待计时、二维记录与中断恢复。"""
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
from threading import Thread
from unittest.mock import patch
from urllib.request import build_opener, ProxyHandler, Request
import json
import unittest

import numpy as np

from main import FourBus, build_continuous_region
from monitor import RunMonitor
from plot import cut_slice, json_value, surface


class MonitorTests(unittest.TestCase):
    def test_next_releases_one_checkpoint_and_continue_releases_the_rest(self):
        now = [1000.]
        with patch('monitor.perf_counter', side_effect=lambda: now[0]):
            with RunMonitor(step_by_step=True, stream=StringIO()) as monitor:
                def compute():
                    for event in ('point', 'sp_end', 'cut'):
                        monitor(event, checkpoint=True)
                worker = Thread(target=compute, daemon=True)
                worker.start()
                with monitor.condition:
                    self.assertTrue(monitor.condition.wait_for(lambda: monitor.state['waiting'], timeout=3.))
                first = monitor.state['revision']
                now[0] += 600.  # 阅读十分钟仍不消耗计算时钟。
                self.assertEqual(monitor.clock(), 1000.)
                self.assertFalse(monitor.control('next', first-1))
                self.assertTrue(monitor.control('next', first))
                with monitor.condition:
                    self.assertTrue(monitor.condition.wait_for(lambda: monitor.state['revision'] == 2, timeout=3.))
                self.assertEqual(len(monitor.history), 2)
                self.assertTrue(monitor.state['waiting'])
                self.assertFalse(monitor.control('next', first))
                monitor.control('continue')
                worker.join(3.)
                self.assertFalse(worker.is_alive())
                self.assertEqual(len(monitor.history), 3)
                self.assertEqual(monitor.paused_seconds, 600.)
                self.assertEqual(json.loads(monitor.snapshot())['elapsed'], 0.)

    def test_http_control_operates_on_the_live_revision(self):
        opener = build_opener(ProxyHandler({}))
        with RunMonitor(show_ui=True, open_browser=False, step_by_step=True, stream=StringIO()) as monitor:
            worker = Thread(target=lambda: monitor('point', checkpoint=True), daemon=True)
            worker.start()
            with monitor.condition:
                self.assertTrue(monitor.condition.wait_for(lambda: monitor.state['waiting'], timeout=3.))
            with opener.open(monitor.url+'state', timeout=3.) as response:
                state = json.load(response)
            request = Request(monitor.url+'control', data=json.dumps(dict(action='next', revision=state['revision'])).encode(),
                              headers={'Content-Type': 'application/json'}, method='POST')
            with opener.open(request, timeout=3.) as response:
                self.assertEqual(response.status, 204)
            worker.join(3.)
            self.assertFalse(worker.is_alive())

    def test_headless_journal_restores_every_committed_step(self):
        with TemporaryDirectory() as folder:
            with RunMonitor(output=folder, show_ui=False, stream=StringIO()) as monitor:
                for i in range(105):
                    monitor('point', point=np.array([i, i+1.]), bounds=[135., 135.])
                path = Path(folder)/'steps.jsonl'
                frames = [json.loads(line) for line in path.read_text(encoding='utf-8').splitlines()]
                self.assertEqual(len(frames), 105)
                self.assertEqual(frames[100]['patch']['bounds'], [135., 135.])
                with RunMonitor(record=False, stream=StringIO()) as restored:
                    restored.load_recording(path)
                    self.assertEqual(restored.state['point'], [104., 105.])
                    self.assertEqual(len(restored.history), 105)
            with path.open('a', encoding='utf-8') as stream:
                stream.write('{"id":105')
            with RunMonitor(record=False, stream=StringIO()) as restored:
                restored.load_recording(path)
                self.assertEqual(len(restored.history), 105)

    def test_2d_geometry_and_plan_conditioned_cut_keep_kw_units(self):
        polygon = surface([[1., 1.], [0., 0.], [0., 1.], [1., 0.]], [2., 3.])
        self.assertEqual(polygon['faces'], [])
        self.assertEqual({tuple(p) for p in polygon['vertices']}, {(0., 0.), (2., 0.), (2., 3.), (0., 3.)})
        plane = cut_slice([-3., 1., 0., 2.], [1], [2., 3.])
        np.testing.assert_allclose(np.asarray(plane['vertices'])[:, 0], 1.)
        self.assertEqual(plane['normal'], [1., 0.])
        self.assertEqual(plane['constant'], -1.)

    def test_single_step_preserves_numeric_result_and_deadline(self):
        network = FourBus()
        plain = build_continuous_region(network, 'socp', 0., [150.]*3, threads=1, time_limit=60.)
        now, received, errors = [1000.], [], []
        with patch('monitor.perf_counter', side_effect=lambda: now[0]):
            with RunMonitor(step_by_step=True, stream=StringIO()) as monitor:
                def compute():
                    try:
                        received.append(build_continuous_region(network, 'socp', 0., [150.]*3, threads=1,
                                        time_limit=60., progress=monitor, clock=monitor.clock))
                    except BaseException as error:
                        errors.append(error)
                worker = Thread(target=compute, daemon=True)
                worker.start()
                with monitor.condition:
                    self.assertTrue(monitor.condition.wait_for(lambda: monitor.state['waiting'], timeout=3.))
                now[0] += 600.
                monitor.control('next', monitor.state['revision'])
                with monitor.condition:
                    self.assertTrue(monitor.condition.wait_for(lambda: monitor.state['waiting'] and monitor.state['revision'] > 1, timeout=3.))
                monitor.control('continue')
                worker.join(30.)
                self.assertFalse(worker.is_alive())
        self.assertEqual(errors, [])
        self.assertEqual(json_value({k:v for k,v in plain.items() if k != 'timing'}),
                         json_value({k:v for k,v in received[0].items() if k != 'timing'}))


if __name__ == '__main__':
    unittest.main()
