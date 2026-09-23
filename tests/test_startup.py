"""简洁入口、统一参数和求解线程配置回归。"""
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
import unittest

import main
from model import DEFAULT_SOLVER_THREADS, PlanningEquations, PlanningModel, PlanningSP, RemainingRegionModel
from vertify import ACPowerFlow


class StartupTests(unittest.TestCase):
    def test_default_threads_and_top_level_settings_reach_run(self):
        with patch('main.run') as run:
            main.main()
        self.assertEqual(run.call_args.kwargs['threads'], 20)
        self.assertEqual(main.SOLVER_THREADS, DEFAULT_SOLVER_THREADS)
        self.assertEqual(len(run.call_args.args[0].planning_corridors), main.CANDIDATE_COUNT)

    def test_thread_detection_error_is_not_retried_or_hidden(self):
        with patch('main.threadpool_limits', side_effect=OSError('GetModuleFileNameEx failed')) as configure:
            with self.assertRaisesRegex(OSError, 'GetModuleFileNameEx failed'):
                main.run(main.FourBus(), budgets=[0.], show_ui=False)
        configure.assert_called_once_with(limits=1)

    def test_threads_reach_each_solver(self):
        equations = PlanningEquations(main.FourBus(), 'socp')
        for options, expected in (({}, 20), ({'threads': 3}, 3)):
            with self.subTest(options=options):
                problem = PlanningModel(equations, **options)
                with problem.model:
                    self.assertEqual(problem.model.Params.Threads, expected)
                self.assertEqual(PlanningSP(equations, **options).settings.max_threads, expected)
                residual = RemainingRegionModel(equations, 0., [100.]*3, 300., [], [], .002, **options)
                with residual.model:
                    self.assertEqual(residual.model.Params.Threads, expected)
                network = main.FourBus()
                ac = ACPowerFlow(network.design(network.initial_plan), **options)
                try:
                    ac._build_global(None)
                    self.assertEqual(ac.model[0].Params.Threads, expected)
                finally:
                    ac.close()

    def test_input_boundary_rejects_invalid_parameters_before_solving(self):
        for options in (dict(threads=0), dict(divisions=0), dict(tau=1.),
                        dict(time_limit=-1.), dict(budgets=[1., 0.]), dict(residual_mode='invalid')):
            with self.subTest(options=options), patch('main.evaluation_bounds') as solve:
                with self.assertRaises(ValueError):
                    main.run(main.FourBus(), show_ui=False, **options)
                solve.assert_not_called()

    def test_candidate_constant_switches_complete_entrypoint(self):
        with TemporaryDirectory() as folder, redirect_stdout(StringIO()):
            for count in (4, 8, 16, 32):
                output = Path(folder)/str(count)
                with self.subTest(candidates=count), patch.multiple(main, CANDIDATE_COUNT=count,
                        NETWORK='case33', BUDGETS=[0.], DIVISIONS=2, SOLVER_THREADS=1,
                        RECOMPUTE=True, SHOW_UI=False, OUTPUT=output):
                    main.main()
                    result = main.BenchmarkResult.load(output)
                    self.assertEqual(result.metadata['candidate_count'], count)
                    self.assertEqual(result.states.shape, (4, 1, 2, 2, 2))
                    for region in result.metadata['continuous']:
                        self.assertEqual(region['status'], 'certified')
                        self.assertEqual(region['max_choice'], main.Case33(candidate_count=count).initial_plan)
                        self.assertLessEqual(region['coverage_bound'], 1e-8)
                    self.assertEqual({p.name for p in output.iterdir()}, {'result.npz'})


if __name__ == '__main__':
    unittest.main()
