"""迁移后的入口、进度语义和本地监视服务回归。"""
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
from threading import Event
from unittest.mock import patch
from urllib.error import HTTPError
from urllib.request import build_opener, ProxyHandler
import json
import subprocess
import sys
import unittest

import numpy as np
from threadpoolctl import threadpool_limits

import main
from plot import RunMonitor, cut_slice, pack_replay
from vertify import validate_ac_region


class ProgressTests(unittest.TestCase):
    def test_packed_replay_preserves_every_frame_and_numeric_value(self):
        geometry = {'vertices': [[0., -0., 1e-12], [1., 2., 3.]], 'faces': [[0, 1, 0]]}
        original = dict(history=[dict(id=i, patch=dict(geometry=geometry, query=i,
                                                      feasible=i % 2 == 0)) for i in range(1000)],
                        budget=None, recording_version=1)
        packed = pack_replay(original)
        restored = []
        def resolve(value):
            return restored[value['ref']] if isinstance(value, dict) else value
        for kind, node in packed['objects']:
            restored.append({k: resolve(v) for k, v in node.items()} if kind else [resolve(v) for v in node])
        decoded = resolve(packed['root'])
        self.assertEqual(json.dumps(decoded), json.dumps(original))
        self.assertIs(decoded['history'][0]['patch']['geometry'], decoded['history'][-1]['patch']['geometry'])

    def test_import_and_help_do_not_start_computation(self):
        for args in (['-c', 'import main; print("imported")'],):
            process = subprocess.run([sys.executable, '-X', 'utf8', *args], cwd=main.ROOT,
                                     capture_output=True, text=True, encoding='utf-8', timeout=10)
            self.assertEqual(process.returncode, 0, process.stderr)
            self.assertNotIn('license', process.stdout.lower())
            self.assertNotIn('实时监视：', process.stdout)

    def test_plane_is_conditioned_on_selection(self):
        cut = np.array([-3., 1., 0., 0., 2.])
        for x, expected in (([0], 3.), ([1], 1.)):
            sliced = cut_slice(cut, x, [5., 6., 7.])
            vertices = np.asarray(sliced['vertices'])
            self.assertEqual(len(vertices), 4)
            np.testing.assert_allclose(vertices[:, 0], expected)
            np.testing.assert_allclose(vertices@cut[1:4]+sliced['constant'], 0., atol=1e-10)
        self.assertEqual(cut_slice([1., 0., 0., 0., -2.], [1], [5.]*3)['vertices'], [])


    def test_unknown_phase_never_reports_full_classification(self):
        answer = dict(bound=None, feasible=False, status='unknown')
        with patch('vertify.ac_planning_query', return_value=answer):
            states = validate_ac_region(main.FourBus(), [0., np.inf], 2, [150.]*3, threads=1)
            self.assertTrue(np.all(states == 0))


    def test_http_snapshot_history_and_cleanup(self):
        opener = build_opener(ProxyHandler({}))
        with TemporaryDirectory() as folder:
            monitor = RunMonitor(show_ui=True, open_browser=False, output=folder, stream=StringIO())
            with monitor:
                states = np.zeros((2, 2, 2, 2), dtype=np.int8)
                monitor('phase_start', states=states, bounds=[5., 6., 7.], budgets=[0., np.inf], phase='linear')
                states[:] = 1  # 验证快照是独立副本，不泄露可变计算数组。
                monitor('point', point=[4., 2., 3.], query=1)
                monitor('sp_start', choice={'line': 'parallel'})
                monitor('cut', cut=[-3., 1., 0., 0., 2.], selection=[1], point=[4., 2., 3.], pool_size=1)
                with opener.open(monitor.url+'state', timeout=3) as response:
                    state = json.load(response)
                self.assertEqual(state['budgets'], [0., None])
                self.assertTrue(all(v == 0 for v in state['states'][0]))
                self.assertEqual(state['events'][-1]['cut']['constant'], -1.)
                with opener.open(monitor.url, timeout=3) as response:
                    self.assertIn('计算过程', response.read().decode('utf-8'))
                with self.assertRaises(HTTPError):
                    opener.open(monitor.url+'result/../main.py', timeout=3)
                monitor('completed', message='done')
                first = json.loads(monitor.snapshot())['elapsed']
                self.assertEqual(json.loads(monitor.snapshot())['elapsed'], first)
                monitor.save_snapshot()
                html = (Path(folder)/'live_view.html').read_text(encoding='utf-8')
                self.assertIn('window.SAVED_REPLAY_GZIP=', html)
                self.assertNotIn('<script src="/plotly.min.js">', html)
            self.assertFalse(monitor.server_thread.is_alive())
            self.assertFalse(monitor.heartbeat_thread.is_alive())

    def test_heartbeat_and_failure_cleanup(self):
        received = Event()
        class Output(StringIO):
            def write(self, value):
                if '此步骤已等待' in value:
                    received.set()
                return super().write(value)
        monitor = RunMonitor(stream=Output(), heartbeat_seconds=.02)
        with self.assertRaisesRegex(RuntimeError, 'probe'):
            with monitor:
                monitor('mp_start', message='求解主问题')
                self.assertTrue(received.wait(2.))
                raise RuntimeError('probe')
        self.assertEqual(monitor.state['status'], 'failed')
        self.assertFalse(monitor.heartbeat_thread.is_alive())


    def test_final_result_round_trip_has_only_core_records(self):
        with TemporaryDirectory() as folder, redirect_stdout(StringIO()), patch('plot.webbrowser.open') as browser:
            result = main.run(main.Case33(candidate_count=4), budgets=[0.], divisions=2,
                              show_ui=True, output=folder, threads=1)
            self.assertEqual(result.states.shape, (4, 1, 2, 2, 2))
            self.assertTrue((Path(folder)/'region_comparison.html').exists())
            browser.assert_called_once()
            self.assertFalse((Path(folder)/'events.jsonl').exists())
            self.assertFalse((Path(folder)/'replay.json').exists())
            for region in result.metadata['continuous']:
                self.assertEqual(set(region['counts']), {'sp', 'cuts'})
                self.assertEqual(set(region['timing']), {'total_seconds'})
                self.assertNotIn('certificates', region)
                self.assertNotIn('queries', region)
                self.assertNotIn('cuts', region)
                for derived in ('inner_volume', 'outer_volume', 'volume_gap', 'validation', 'tau', 'residual_mode'):
                    self.assertNotIn(derived, region)
                self.assertTrue(all('choice' in p and 'cost' in p for p in region['inner']))
            with patch('main.evaluation_bounds', side_effect=AssertionError('unexpected solve')):
                loaded = main.run(main.FourBus(), recompute=False, output=folder, show_ui=False)
            np.testing.assert_array_equal(loaded.states, result.states)
            code = 'import sys; sys.path.insert(0, '+repr(str(main.ROOT))+'); import main; main.run(main.FourBus(), recompute=False, show_ui=False, output='+repr(folder)+')'
            process = subprocess.run([sys.executable, '-X', 'utf8', '-c', code], cwd=folder,
                                     capture_output=True, text=True, encoding='utf-8', timeout=10)
            self.assertEqual(process.returncode, 0, process.stderr)


if __name__ == '__main__':
    unittest.main()
