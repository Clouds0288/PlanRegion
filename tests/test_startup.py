"""配置驱动的维数、线程设置及主入口。"""
from pathlib import Path
import tempfile
from unittest.mock import patch
import unittest
import numpy as np

import main
import model
from model import GridPhysics, MasterProblem, SubProblem
from tests.planning_checks import mainline_stubs


class StartupTests(unittest.TestCase):
    def test_case33_entry_uses_two_loads_and_preserves_background(self):
        with tempfile.TemporaryDirectory() as folder, mainline_stubs() as calls:
            main.main('case33', 2, show=False, output=Path(folder))
        network, kwargs = calls.build.call_args.args[0], calls.build.call_args.kwargs
        self.assertIsInstance(network, main.Case33)
        self.assertEqual(network.load_nodes, (18, 25))
        self.assertEqual(network.n_types, 37)
        self.assertEqual(kwargs['budget'], network.switch_budget)
        self.assertEqual(kwargs['monitor'].output.name, 'case33_18_25.json.gz')
        fixed = ~np.isin(network.nodes, network.load_nodes)
        np.testing.assert_array_equal(network.fixed_p[fixed], network.original_p[fixed])
        np.testing.assert_array_equal(network.fixed_q[fixed], network.original_q[fixed])
        self.assertEqual(network.fixed_p[network.nodes.index(33)], 60.)

    def test_top_level_settings_reach_construction(self):
        with tempfile.TemporaryDirectory() as folder, mainline_stubs() as calls:
            main.main('fourbus', show=False, output=Path(folder))
        network, kwargs = calls.build.call_args.args[0], calls.build.call_args.kwargs
        self.assertIsInstance(network, main.FourBus)
        self.assertEqual(network.load_nodes, (1, 2, 3)[:main.DIMENSION])
        self.assertEqual(kwargs['budget'], 20000.)
        self.assertEqual(kwargs['seconds'], main.CASE_TIME_LIMIT)
        self.assertEqual(kwargs['settings'].workers, main.WORKERS)
        self.assertEqual((kwargs['settings'].tau, kwargs['settings'].threads), (main.REGION_TAU, main.SOLVER_THREADS))
        self.assertEqual(network.fixed_p[2], 0.)
        self.assertEqual(network.fixed_q[2], 0.)
        self.assertTrue(network.required.all())

    def test_partition_pool_error_is_not_retried_or_hidden(self):
        with tempfile.TemporaryDirectory() as folder, \
             patch('region.ProcessPoolExecutor', side_effect=OSError('pool probe')) as pool:
            with self.assertRaisesRegex(OSError, 'pool probe'):
                main.main('fourbus', 2, show=False, output=Path(folder))
        pool.assert_called_once()

    def test_threads_reach_each_solver(self):
        equations = GridPhysics(main.FourBus(), [1, 1, 1])
        for threads in (1, 3):
            with self.subTest(threads=threads):
                problem = MasterProblem(equations, threads=threads)
                with problem.model:
                    self.assertEqual(problem.model.Params.Threads, threads)
                oracle = SubProblem(equations, threads=threads)
                with patch('model.new_model', wraps=model.new_model) as create:
                    network = equations.network
                    oracle.solve(network.encode_plan(network.initial_plan), np.zeros(3))
                self.assertEqual(create.call_args.args, ('planning_SP', threads))


if __name__ == '__main__':
    unittest.main()
