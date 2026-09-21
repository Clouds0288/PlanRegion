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
from plot import RunMonitor, cut_slice
from vertify import validate_ac_region


class ProgressTests(unittest.TestCase):
    def test_import_and_help_do_not_start_computation(self):
        for args in (['-c', 'import main; print("imported")'], ['main.py', '--help']):
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

    def test_real_observer_does_not_change_regions(self):
        network = main.Case33(candidate_count=4)
        bounds = [400., 4670., 630.]
        with threadpool_limits(limits=1), RunMonitor(stream=StringIO()) as monitor:
            for method in ('linear', 'socp', 'hybrid'):
                args = (network, method, 0., bounds)
                baseline = main.build_continuous_region(*args)
                observed = main.build_continuous_region(*args, observer=monitor)
                self.assertEqual(observed['status'], 'certified')
                self.assertEqual(observed['inner'], baseline['inner'])
                self.assertEqual(observed['outer'], baseline['outer'])
                self.assertEqual(observed['counts'], baseline['counts'])
            args = (network, network.budgets, 2, bounds)
            np.testing.assert_array_equal(validate_ac_region(*args), validate_ac_region(*args, observer=monitor))
            self.assertIsNone(monitor.server)
            self.assertEqual(len(monitor.events), 0)
            self.assertNotIn('states', monitor.state)

    def test_unknown_phase_never_reports_full_classification(self):
        answer = dict(bound=None, feasible=False, status='unknown')
        with patch('vertify.ac_planning_query', return_value=answer):
            with RunMonitor(stream=StringIO()) as monitor:
                states = validate_ac_region(main.FourBus(), [0., np.inf], 2, [150.]*3, observer=monitor)
                self.assertTrue(np.all(states == 0))
                self.assertEqual(monitor.state['event'], 'phase_end')
                self.assertTrue(all(c['unknown'] == c['total'] for c in monitor.state['counts']))

    def test_hybrid_phase_reset_is_visible(self):
        phases = []
        original = main.joint_benders
        with threadpool_limits(limits=1), RunMonitor(stream=StringIO()) as monitor:
            def query(equations, **kwargs):
                if equations.method == 'socp':
                    self.assertGreater(len(kwargs['cuts']), 0)
                    return dict(bound=None, feasible=False), []
                return original(equations, **kwargs)
            def observe(event, **data):
                monitor(event, **data)
                if event == 'phase_start':
                    phases.append((data['phase'], monitor.state['geometry']))
            with patch('main.joint_benders', side_effect=query):
                result = main.build_continuous_region(main.Case33(candidate_count=4), 'hybrid', 0.,
                                                       [400., 4670., 630.], observer=observe)
        self.assertEqual([p[0] for p in phases], ['linear', 'socp'])
        self.assertEqual(phases[1][1], [])
        self.assertEqual(result['linear_status'], 'certified')
        self.assertEqual(result['status'], 'unknown')
        self.assertEqual(result['inner'], [])

    def test_http_snapshot_history_and_cleanup(self):
        opener = build_opener(ProxyHandler({}))
        with TemporaryDirectory() as folder:
            monitor = RunMonitor(show_ui=True, open_browser=False, output=folder, stream=StringIO())
            with monitor:
                states = np.zeros((2, 2, 2, 2), dtype=np.int8)
                monitor('phase_start', states=states, bounds=[5., 6., 7.], budgets=[0., np.inf], phase='linear')
                states[:] = 1  # 验证快照是独立副本，不泄露可变计算数组。
                monitor('point', point=[4., 2., 3.], query=1)
                monitor('sp_start', choice=[1])
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

    def test_full_ui_run_and_load_from_another_working_directory(self):
        # 全流程在独立目录保存；读取路径不应再次触发优化器。
        with TemporaryDirectory() as folder, redirect_stdout(StringIO()):
            result = main.run(main.Case33(candidate_count=4), budgets=[0.], divisions=2, show_ui=True,
                              open_browser=False, output=folder, plots=False)
            self.assertEqual(result.states.shape, (4, 1, 2, 2, 2))
            self.assertIn('main.py', result.metadata['hashes'])
            self.assertTrue((Path(folder)/'result.npz').exists())
            self.assertTrue((Path(folder)/'live_view.html').exists())
            with patch('main.PlanningModel', side_effect=AssertionError('unexpected solve')):
                loaded = main.run(main.FourBus(), recompute=False, output=folder, plots=False)
            np.testing.assert_array_equal(loaded.states, result.states)
            process = subprocess.run([sys.executable, str(main.ROOT/'main.py'), '--load', '--no-plots',
                                      '--no-ui', '--output', folder], cwd=folder, capture_output=True,
                                     text=True, encoding='utf-8', timeout=10)
            self.assertEqual(process.returncode, 0, process.stderr)


if __name__ == '__main__':
    unittest.main()
