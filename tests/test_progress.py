"""原生监控不会改变数学结果，主线导入没有计算或 GUI 副作用。"""
import json
import subprocess
import sys
import unittest
from unittest.mock import patch

import numpy as np
from threadpoolctl import threadpool_limits

import main
import continuous
from monitor import RunMonitor
from plot import cut_slice, json_value, pack_replay
from vertify import validate_ac_region


class ProgressTests(unittest.TestCase):
    def test_import_and_help_do_not_start_computation(self):
        process = subprocess.run([sys.executable, '-X', 'utf8', '-c', 'import main; print("imported")'],
            cwd=main.ROOT, capture_output=True, text=True, encoding='utf-8', timeout=10)
        self.assertEqual(process.returncode, 0, process.stderr)
        self.assertEqual(process.stdout.strip(), 'imported')

    def test_plane_is_conditioned_on_selection(self):
        cut = np.array([-3., 1., 0., 0., 2.])
        for x, expected in (([0], 3.), ([1], 1.)):
            sliced = cut_slice(cut, x, [5., 6., 7.])
            vertices = np.asarray(sliced['vertices'])
            np.testing.assert_allclose(vertices[:, 0], expected)
            np.testing.assert_allclose(vertices@cut[1:4]+sliced['constant'], 0., atol=1e-10)

    def test_live_cuts_preserve_solver_result_and_show_clipped_geometry(self):
        net, bounds, cuts = main.FourBus(load_nodes=(1, 2)), [142.5, 142.5], []
        def callback(event, **data):
            if event == 'cut':
                cut = data['cut'].copy()
                for row in data['records']:
                    values = cut[0]+row['outer']@np.diag(bounds)@cut[1:3]+cut[3:]@row['x']
                    self.assertTrue(np.all(values >= -1e-8))
                cuts.append(cut)
        monitor = RunMonitor(callback=callback)
        with threadpool_limits(limits=1):
            shown = continuous.build_continuous_region(net, 'socp', 20000., bounds, threads=1, progress=monitor)
            plain = continuous.build_continuous_region(net, 'socp', 20000., bounds, threads=1)
        self.assertEqual(json_value({k: v for k, v in shown.items() if k != 'timing'}),
                         json_value({k: v for k, v in plain.items() if k != 'timing'}))
        self.assertGreater(len(cuts), 0)
        self.assertEqual(len(cuts), shown['counts']['cuts'])
        self.assertEqual(sum(f['patch']['event'] == 'cut' for f in monitor.history), len(cuts))

    def test_unknown_phase_never_reports_full_classification(self):
        with patch('vertify.ac_planning_query', side_effect=RuntimeError('AC not converged')):
            with self.assertRaisesRegex(RuntimeError, 'AC not converged'):
                validate_ac_region(main.FourBus(), [0., np.inf], 2, [150.]*3, threads=1)

    def test_replay_does_not_start_solver_or_browser(self):
        with patch('model.MasterProblem', side_effect=AssertionError('unexpected solve')), \
             patch('webbrowser.open', side_effect=AssertionError('unexpected browser')):
            monitor = RunMonitor()
            monitor._emit('point', point=[1., 2.])
            self.assertEqual(monitor.frame(0)['point'], [1., 2.])


if __name__ == '__main__':
    unittest.main()
