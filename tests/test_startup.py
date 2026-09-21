"""Windows 线程检测失败以及主入口候选数量切换的回归。"""
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import Mock, patch
import json
import unittest

import main


class StartupTests(unittest.TestCase):
    def test_transient_thread_detection_failure_retries_and_restores(self):
        limiter = Mock()
        failures = [OSError('GetModuleFileNameEx failed'),
                    OSError('EnumProcessModulesEx failed'), limiter]
        with patch('main.threadpool_limits', side_effect=failures) as configure:
            with main.numerical_threads() as status:
                self.assertTrue(status['applied'])
                self.assertEqual(status['attempts'], 3)
                limiter.restore_original_limits.assert_not_called()
        self.assertEqual(configure.call_count, 3)
        limiter.restore_original_limits.assert_called_once_with()

    def test_detection_failure_does_not_abort_real_computation(self):
        with TemporaryDirectory() as folder, redirect_stdout(StringIO()):
            with patch('main.threadpool_limits', side_effect=OSError('GetModuleFileNameEx failed')):
                with self.assertWarnsRegex(RuntimeWarning, '线程检测重试 3 次仍失败'):
                    result = main.run(main.FourBus(), budgets=[0.], divisions=2, show_ui=False,
                                      output=folder, plots=False)
            status = result.metadata['thread_control']
            self.assertFalse(status['applied'])
            self.assertEqual(status['attempts'], 3)
            self.assertTrue(all(r['status'] == 'certified' for r in result.metadata['continuous']))
            log = [json.loads(line) for line in (Path(folder)/'events.jsonl').read_text(encoding='utf-8').splitlines()]
            self.assertTrue(any(e['patch'].get('event') == 'runtime_warning' for e in log))
            self.assertEqual(log[-1]['patch']['event'], 'completed')

    def test_unrelated_initialization_errors_are_not_suppressed(self):
        with patch('main.threadpool_limits', side_effect=OSError('unrelated problem')) as configure:
            with self.assertRaisesRegex(OSError, 'unrelated problem'):
                with main.numerical_threads():
                    self.fail('unexpected computation')
        self.assertEqual(configure.call_count, 1)

    def test_computation_errors_propagate_and_restore_threads(self):
        limiter = Mock()
        with patch('main.threadpool_limits', return_value=limiter) as configure:
            with self.assertRaisesRegex(OSError, 'GetModuleFileNameEx failed'):
                with main.numerical_threads():
                    raise OSError('GetModuleFileNameEx failed')
        self.assertEqual(configure.call_count, 1)
        limiter.restore_original_limits.assert_called_once_with()

    def test_candidate_constant_switches_complete_entrypoint(self):
        with TemporaryDirectory() as folder, redirect_stdout(StringIO()):
            for count in (4, 8, 16):
                with self.subTest(candidates=count), patch('main.CANDIDATE_COUNT', count):
                    output = Path(folder)/str(count)
                    code = main.main(['--network', 'case33', '--budgets', '0', '--divisions', '2',
                                      '--recompute', '--no-ui', '--no-plots', '--no-hold',
                                      '--output', str(output)])
                    self.assertEqual(code, 0)
                    result = main.BenchmarkResult.load(output)
                    self.assertEqual(result.metadata['candidate_count'], count)
                    self.assertEqual(result.states.shape, (4, 1, 2, 2, 2))
                    self.assertEqual(len(result.metadata['continuous']), 3)
                    for region in result.metadata['continuous']:
                        self.assertEqual(region['status'], 'certified')
                        self.assertEqual(len(region['max_choice']), count)
                        self.assertLessEqual(region['coverage_bound'], 1e-8)


if __name__ == '__main__':
    unittest.main()
