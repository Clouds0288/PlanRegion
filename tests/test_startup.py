"""配置驱动的维数、线程设置及主入口。"""
from unittest.mock import patch
import unittest
import numpy as np

import main
import model
from model import GridPhysics, MasterProblem, SubProblem


class StartupTests(unittest.TestCase):
    def test_case33_entry_uses_two_loads_and_preserves_background(self):
        with patch('main.run') as run:
            main.main('case33', dimension=2)
        network = run.call_args.args[0]
        self.assertIsInstance(network, main.Case33)
        self.assertEqual(network.load_nodes, (18, 25))
        self.assertEqual(network.n_types, 37)
        self.assertEqual(run.call_args.kwargs['budget'], network.switch_budget)
        self.assertEqual(run.call_args.kwargs['divisions'], main.SCAN_DIVISIONS[2])
        self.assertEqual(run.call_args.kwargs['output'].name, 'case33_18_25.json.gz')
        fixed = ~np.isin(network.nodes, network.load_nodes)
        np.testing.assert_array_equal(network.fixed_p[fixed], network.original_p[fixed])
        np.testing.assert_array_equal(network.fixed_q[fixed], network.original_q[fixed])
        self.assertEqual(network.fixed_p[network.nodes.index(33)], 60.)

    def test_default_threads_and_top_level_settings_reach_run(self):
        with patch('main.run') as run:
            main.main('fourbus')
        network = run.call_args.args[0]
        self.assertIsInstance(network, main.FourBus)
        self.assertEqual(network.load_nodes, (1, 2, 3)[:main.DIMENSION])
        self.assertEqual(run.call_args.kwargs['budget'], 20000.)
        self.assertEqual(run.call_args.kwargs['threads'], main.SOLVER_THREADS)
        self.assertEqual(run.call_args.kwargs['time_limit'], main.CASE_TIME_LIMIT)
        self.assertEqual(network.fixed_p[2], 0.)
        self.assertEqual(network.fixed_q[2], 0.)
        self.assertTrue(network.required.all())

    def test_partition_pool_error_is_not_retried_or_hidden(self):
        with patch('region.ProcessPoolExecutor', side_effect=OSError('pool probe')) as pool:
            with self.assertRaisesRegex(OSError, 'pool probe'):
                main.run(main.FourBus(load_nodes=(1, 2)), show_ui=False, output=None)
        pool.assert_called_once()

    def test_threads_reach_each_solver(self):
        equations = GridPhysics(main.FourBus(), 'socp')
        for options, expected in (({}, 20), ({'threads': 3}, 3)):
            with self.subTest(options=options):
                problem = MasterProblem(equations, **options)
                with problem.model:
                    self.assertEqual(problem.model.Params.Threads, expected)
                oracle = SubProblem(equations, **options)
                with patch('model.new_model', wraps=model.new_model) as create:
                    network = equations.network
                    oracle.solve(network.encode_plan(network.initial_plan), np.zeros(3))
                self.assertEqual(create.call_args.args, ('planning_SP', expected))


if __name__ == '__main__':
    unittest.main()
