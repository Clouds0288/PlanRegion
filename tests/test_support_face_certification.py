"""Numerical certificates and failure semantics of the isolated FourBus experiment."""
from contextlib import redirect_stdout
from io import StringIO
from itertools import product
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import gurobipy as gp
from gurobipy import GRB
import numpy as np

from experiments.test_support_face_certification_fourbus_2d import (
    SupportOracle, cached_support, classify_support, compare_global,
    conditional_support_cut, expanded_polygon,
    geometry_margins, normal_key, run_experiment)
from model import GridPhysics, MasterProblem, PLANNING_TOL, RemainingRegionModel, SubProblem
from monitor import RunMonitor
from Network.four_bus_five_corridor import FourBus
from region import GEOMETRY_TOL, contains, halfspaces
from tests.planning_checks import margin


class SupportGeometryTests(unittest.TestCase):
    def test_geometry_uses_all_faces_and_vertices_without_solver(self):
        inner = np.array(((0., 0.), (.6, 0.), (.6, .4), (0., .4)))
        outer = np.array(((0., 0.), (.7, 0.), (.7, .5), (0., .5)))
        faces = halfspaces(inner)
        with patch('gurobipy.Model', side_effect=AssertionError('Geometry must not build a solver')):
            upper, rho = geometry_margins(faces, outer, .005)
        np.testing.assert_allclose(upper, [max(v @ f[:2] for v in outer) for f in faces])
        np.testing.assert_allclose(rho, .995*upper+faces[:, -1])
        self.assertEqual(int(np.count_nonzero(rho <= GEOMETRY_TOL)), 2)

    def test_low_dimensional_hull_is_not_vacuously_certified(self):
        faces = halfspaces(np.array(((0., 0.), (1., 0.))))
        outer = np.array(((0., 0.), (1., 0.), (0., 1.)))
        _, rho = geometry_margins(faces, outer, .005)
        self.assertGreater(rho.max(), GEOMETRY_TOL)
        _, rho = geometry_margins(faces, np.empty((0, 2)), .005)
        self.assertTrue(np.isposinf(rho).all())

    def test_incumbent_below_face_is_not_a_certificate(self):
        result, _, _ = classify_support(.8, 1.2, np.array((.8, 0.)),
                                       np.array((1., 0.)), -1., 0., GEOMETRY_TOL)
        self.assertEqual(result, 'UNRESOLVED')
        result, _, _ = classify_support(-np.inf, .9, None, np.array((1., 0.)),
                                       -1., 0., GEOMETRY_TOL)
        self.assertEqual(result, 'SUPPORT_CERTIFIED')

    def test_actual_audited_point_required_for_violation(self):
        for point, expected in ((None, 'UNRESOLVED'), (np.array((1.1, 0.)), 'VIOLATED')):
            result, _, _ = classify_support(1.1, np.inf, point, np.array((1., 0.)),
                                           -1., 0., GEOMETRY_TOL)
            self.assertEqual(result, expected)

    def test_nearby_direction_cache_adjusts_bound_and_preserves_orientation(self):
        normal = np.array((.6, .8))
        changed = np.array((.6+2e-12, .8-1e-12))
        self.assertEqual(normal_key(normal), normal_key(changed))
        self.assertNotEqual(normal_key(normal), normal_key(-normal))
        # Exact support of [0,1]^2; the perturbed direction must not get a smaller UB.
        entry = dict(normal=normal, lb=1.4, ub=1.4, point=np.ones(2), status=GRB.TIME_LIMIT)
        reused = cached_support(entry, 3*changed, np.ones(2))
        self.assertGreaterEqual(reused['ub']+1e-15, np.maximum(3*changed, 0.).sum())
        self.assertAlmostEqual(reused['lb'], (3*changed).sum())
        # The same support value with a changed offset requires a new face decision.
        self.assertEqual(classify_support(1.4, 1.4, np.ones(2), normal,
                         -1.5, 0., GEOMETRY_TOL)[0], 'SUPPORT_CERTIFIED')
        self.assertEqual(classify_support(1.4, 1.4, np.ones(2), normal,
                         -1.3, 0., GEOMETRY_TOL)[0], 'VIOLATED')

    def test_expanded_polygon_matches_the_requested_halfspace_definition(self):
        bounds = np.array((100., 200.))
        inner = np.array(((0., 0.), (30., 0.), (20., 50.), (0., 80.)))
        tau, epsilon = .1, 1e-8
        polygon = expanded_polygon(inner, bounds, tau, epsilon)
        faces = halfspaces(inner/bounds)
        self.assertTrue(contains((1-tau)*polygon/bounds, faces, epsilon+1e-12).all())
        self.assertTrue(contains(inner/bounds, halfspaces(polygon/bounds)).all())

    def test_conditional_cut_does_not_constrain_other_binary_schemes(self):
        x = np.array((1, 0, 1))
        bounds, axis_bounds = np.array((100., 200.)), np.array((80., 120.))
        normal, ub = np.array((-.6, .8)), .2
        cut = conditional_support_cut(x, normal, ub, bounds, axis_bounds)
        for bits in product((0, 1), repeat=3):
            for corner in product((0, 1), repeat=2):
                power = np.asarray(corner)*axis_bounds
                value = cut[0]+cut[1:3] @ power+cut[3:] @ bits
                if np.array_equal(bits, x):
                    self.assertAlmostEqual(value, ub-normal @ (power/bounds), places=14)
                else:
                    self.assertGreaterEqual(value, -1e-14)


class SupportNumericalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.network = FourBus(load_nodes=(1, 2))
        cls.equations = GridPhysics(cls.network, 'socp')
        cls.bounds = np.full(2, cls.network.power_limit)

    def test_normalized_objective_and_cached_scaling_match_original_master(self):
        net = self.network
        x = net.encode_plan(net.initial_plan)
        normal = np.array((.6, .8))
        oracle = SupportOracle(self.equations, x, self.bounds, 20000., 1, 10.)
        try:
            answer, hit = oracle.get(normal)
            self.assertFalse(hit)
            self.assertIsNotNone(answer['point'])
            self.assertLessEqual(answer['max_violation'], PLANNING_TOL)
            self.assertAlmostEqual(answer['objective'], normal @ answer['point'], places=10)
            self.assertGreaterEqual(margin(self.equations, x, answer['point']*self.bounds,
                                         answer['state']), -PLANNING_TOL)
            again, hit = oracle.get(2*normal)
            self.assertTrue(hit)
            self.assertEqual(oracle.calls, 1)
            self.assertAlmostEqual(again['ub'], 2*answer['ub'], places=12)
            problem = MasterProblem(self.equations, fixed_plan=net.decode_plan(x),
                                    budget=20000., direction=normal/self.bounds, threads=1)
            with problem.model:
                reference = problem.solve()
            self.assertLessEqual(reference['objective'], answer['ub']+1e-8)
            self.assertLessEqual(answer['lb'], reference['bound']+1e-8)
            self.assertLess(answer['ub']-answer['lb'], 1e-7)
        finally:
            oracle.close()

    def test_time_limit_without_bounds_remains_unresolved(self):
        x = self.network.encode_plan(self.network.initial_plan)
        oracle = SupportOracle(self.equations, x, self.bounds, 20000., 1, 0.)
        try:
            answer, _ = oracle.get(np.array((1., 0.)))
            self.assertEqual(answer['status'], GRB.TIME_LIMIT)
            self.assertIsNone(answer['point'])
            self.assertEqual(classify_support(answer['lb'], answer['ub'], answer['point'],
                             np.array((1., 0.)), -1., .005, GEOMETRY_TOL)[0], 'UNRESOLVED')
        finally:
            oracle.close()

    def test_all_schemes_current_faces_and_every_saved_witness_are_certified(self):
        module = 'experiments.test_support_face_certification_fourbus_2d'
        forbidden = AssertionError('Adaptive construction must not enumerate or exclude schemes')
        with (TemporaryDirectory() as directory, redirect_stdout(StringIO()),
              patch(f'{module}.budget_schemes', create=True, side_effect=forbidden),
              patch('vertify.budget_schemes', side_effect=forbidden),
              patch.object(MasterProblem, 'exclude', side_effect=forbidden),
              patch.object(SubProblem, 'solve', side_effect=AssertionError('No candidate SP')),
              patch(f'{module}.RemainingRegionModel', wraps=RemainingRegionModel) as search):
            result = run_experiment(output=Path(directory), scan=False, plots=False, verbose=False)
        self.assertTrue(search.call_args_list)
        self.assertTrue(all(call.args[0].method == 'socp' for call in search.call_args_list))
        self.assertEqual(result['coverage']['mode'], 'physical')
        self.assertEqual(result['totals']['global_milp_calls'], 0)
        self.assertEqual(result['totals']['discovery_sp_calls'], 0)
        self.assertEqual(result['totals']['discovery_cut_lp_calls'], 0)
        self.assertEqual(result['totals']['global_physical_calls'], len(search.call_args_list))
        for step in result['discovery_history']:
            if step.get('candidate_feasible'):
                self.assertLessEqual(step['max_violation'], PLANNING_TOL)
                self.assertGreaterEqual(margin(self.equations, step['x'], step['p'], step['state']),
                                        -PLANNING_TOL)
        self.assertTrue(result['totals']['global_certified'])
        self.assertTrue(result['coverage']['complete'])
        self.assertEqual(result['totals']['enumeration_mip_calls'], 0)
        self.assertEqual(result['totals']['enumeration_seconds'], 0.)
        self.assertGreater(result['totals']['schemes'], 0)
        self.assertLess(result['totals']['schemes'], 17)
        self.assertLess(result['totals']['refinement_support_calls'], result['totals']['faces_seen']/3)
        for row in result['schemes']:
            self.assertTrue(row['certified'])
            self.assertLessEqual(row['max_remaining_margin'], GEOMETRY_TOL)
            self.assertEqual(row['support_calls'], len(row['support_history']))
            self.assertTrue(set(row['final_face_statuses']) <= {'GEOMETRY_CERTIFIED', 'SUPPORT_CERTIFIED'})
            for witness in row['certificates']:
                self.assertGreaterEqual(margin(self.equations, row['x'], witness['p'],
                                             witness['state']), -PLANNING_TOL)
            # Convex combinations stay within the same scheme, including their state.
            a, b = row['certificates'][0], row['certificates'][-1]
            self.assertGreaterEqual(margin(self.equations, row['x'], (a['p']+b['p'])/2,
                                         (a['state']+b['state'])/2), -PLANNING_TOL)
            for step in row['iteration_history']:
                history = [f for f in row['face_history'] if f['iteration'] == step['iteration']]
                if any(f['status'] == 'VIOLATED' for f in history):
                    # No second solve on stale faces after accepting the first violation.
                    self.assertEqual(sum(f['ub'] is not None for f in history), 1)
        # A fresh physical query without generated cuts rechecks the certificate.
        settings = result['settings']
        reference = compare_global(self.equations, result['schemes'], 20000.,
            settings['bounds'], settings['axis_bounds'], settings['total_bound'],
            settings['tau'], settings['epsilon_geom'], 1, 30.)
        self.assertTrue(reference['complete'], reference)

    def test_local_certification_without_global_completion_cannot_stop(self):
        with TemporaryDirectory() as directory, redirect_stdout(StringIO()):
            result = run_experiment(output=directory, max_discovery_iterations=1,
                                    scan=False, plots=False, verbose=False)
        self.assertTrue(all(row['certified'] for row in result['schemes']))
        self.assertFalse(result['totals']['global_certified'])
        self.assertFalse(result['coverage']['complete'])
        self.assertEqual(result['coverage']['status'], 'DISCOVERY_LIMIT')

    def test_global_timeout_keeps_local_results_without_certifying_union(self):
        with (TemporaryDirectory() as directory, redirect_stdout(StringIO()),
              patch.object(RemainingRegionModel, 'solve', side_effect=TimeoutError('test limit'))):
            result = run_experiment(output=directory, scan=False, plots=False, verbose=False)
            self.assertTrue((Path(directory)/'summary.json').is_file())
        self.assertTrue(result['schemes'][0]['certified'])
        self.assertFalse(result['totals']['global_certified'])
        self.assertEqual(result['coverage']['status'], 'TIME_LIMIT')

    def test_all_joint_cuts_are_valid_over_the_complete_physical_model(self):
        with TemporaryDirectory() as directory, redirect_stdout(StringIO()):
            result = run_experiment(output=directory, scan=False, plots=False, verbose=False)
        self.assertTrue(result['totals']['global_certified'])
        # Maximize each cut's violation over the original mixed-integer SOCP,
        # without the generated cuts or any restriction to discovered schemes.
        problem = MasterProblem(self.equations, budget=20000., threads=1)
        with problem.model as model:
            model.Params.TimeLimit = 30.
            for cut in result['cuts']:
                expression = float(cut[0])+gp.quicksum(float(c)*p for c, p in
                    zip(cut[1:3], problem.loads.values()))+gp.quicksum(float(c)*x for c, x in
                    zip(cut[3:], problem.choices.values()))
                model.setObjective(-expression, GRB.MAXIMIZE)
                model.optimize()
                self.assertEqual(model.Status, GRB.OPTIMAL)
                self.assertLessEqual(model.ObjBound, PLANNING_TOL)

    def test_iteration_limit_cannot_become_a_coverage_certificate(self):
        with TemporaryDirectory() as directory, redirect_stdout(StringIO()):
            result = run_experiment(output=directory, budget=0., max_iterations=1,
                                    scan=False, plots=False, verbose=False)
        self.assertFalse(result['totals']['global_certified'])
        self.assertEqual(result['schemes'][0]['status'], 'ITERATION_LIMIT')


class SupportReplayTests(unittest.TestCase):
    def test_native_recording_preserves_algorithm_and_does_not_leak_future(self):
        with TemporaryDirectory() as directory, redirect_stdout(StringIO()):
            baseline = run_experiment(output=Path(directory)/'plain', scan=False,
                                      plots=False, verbose=False)
            recorded = run_experiment(output=Path(directory)/'recorded', scan=False,
                                      plots=False, verbose=False, record=True)
            restored = RunMonitor()
            restored.load_recording(recorded['recording'])
        for key in ('schemes', 'global_certified', 'support_calls', 'global_milp_calls', 'global_physical_calls',
                    'discovery_sp_calls', 'discovery_cut_lp_calls', 'faces_seen'):
            self.assertEqual(baseline['totals'][key], recorded['totals'][key])
        for a, b in zip(baseline['schemes'], recorded['schemes']):
            np.testing.assert_array_equal(a['x'], b['x'])
            np.testing.assert_allclose(a['inner'], b['inner'], atol=1e-9, rtol=0.)
            np.testing.assert_allclose(a['outer'], b['outer'], atol=1e-9, rtol=0.)
        first = restored.frame(0)
        self.assertFalse(first.get('schemes'))
        self.assertFalse(first.get('cut_history'))
        self.assertNotIn('result', first)
        self.assertEqual(len(restored.state['schemes']), recorded['totals']['schemes'])
        self.assertEqual(restored.state['sp'], recorded['totals']['discovery_sp_calls'])
        self.assertEqual(restored.state['global_search'], recorded['totals']['global_physical_calls'])
        self.assertEqual(restored.state['support_calls'], recorded['totals']['support_calls'])
        self.assertEqual(len(restored.state['cut_history']), len(recorded['cuts']))
        self.assertTrue(restored.state['result']['certified'])
        self.assertTrue(restored.state['coverage_complete'])
        kinds = {row['kind'] for row in restored.state['cut_history'].values()}
        self.assertEqual(kinds, {'support'})
        for i, item in enumerate(restored.history):
            frame = restored.frame(i)
            event = item['patch']['event']
            if event == 'scheme_start':
                key = frame['active_scheme']
                self.assertEqual(len(frame['schemes'][key]['inner']), 1)
            if event == 'residual_end' and frame['global_point'] is not None:
                self.assertIn('原约束复核通过', frame['phase'])
                self.assertFalse(frame['coverage_complete'])
        self.assertFalse(any(item['patch']['event'] in ('point', 'sp_end') for item in restored.history))


if __name__ == '__main__':
    unittest.main()
