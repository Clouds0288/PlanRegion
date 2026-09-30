"""3D geometry and signed physical certificates for the isolated experiment."""
from itertools import product
import unittest
from unittest.mock import patch

import numpy as np

from experiments.test_support_face_certification_fourbus_2d import (
    SupportOracle, certify_scheme, discover_schemes, expanded_polygon, initial_bounds, zero_certificate)
from model import MasterProblem, PortPhysics, PLANNING_TOL, SubProblem, port_bounds
from Network.four_bus_five_corridor import FourBus
from region import contains, halfspaces, GEOMETRY_TOL
from tests.planning_checks import margin


class Support3DTests(unittest.TestCase):
    def test_mixed_sign_discovery_never_calls_candidate_sp(self):
        sign = np.array((1, -1, 1))
        network = FourBus(load_nodes=(1, 2, 3))
        equations = PortPhysics(network, sign)
        bounds, axis, total, seed = initial_bounds(equations, 20000., 1, 10., bounds=port_bounds(network))
        with patch.object(SubProblem, 'solve', side_effect=AssertionError('No candidate SP')):
            result = discover_schemes(equations, seed, bounds, axis, total,
                budget=20000., tau=.005, epsilon_geom=GEOMETRY_TOL, threads=1,
                support_seconds=10., global_seconds=10., max_iterations=300,
                max_discovery_iterations=1, outer_cuts=True, verbose=False)
        self.assertIsNone(result['error'])
        self.assertEqual(result['global_milp_calls'], 0)
        self.assertEqual(result['global_physical_calls'], 1)
        self.assertEqual(result['discovery_sp_calls'], 0)
        step = result['history'][0]
        self.assertTrue(step['candidate_feasible'])
        self.assertLessEqual(step['max_violation'], PLANNING_TOL)
        self.assertGreaterEqual(margin(equations, step['x'], step['p']*sign, step['state']), -PLANNING_TOL)

    def test_three_dimensional_expansion_matches_all_faces(self):
        bounds = np.array((100., 200., 150.))
        inner = np.array(((0., 0., 0.), (40., 0., 0.), (0., 75., 0.), (0., 0., 90.), (20., 30., 40.)))
        expanded = expanded_polygon(inner, bounds, .005, GEOMETRY_TOL)
        self.assertEqual(expanded.shape[1], 3)
        self.assertTrue(contains(.995*expanded/bounds, halfspaces(inner/bounds), GEOMETRY_TOL+1e-12).all())
        self.assertTrue(contains(inner/bounds, halfspaces(expanded/bounds)).all())

    def test_each_signed_partition_zero_and_support_are_physical(self):
        for sign in product((1, -1), repeat=3):
            with self.subTest(sign=sign):
                network = FourBus(load_nodes=(1, 2, 3))
                equations = PortPhysics(network, sign)
                bounds = port_bounds(network)
                x = network.encode_plan(network.initial_plan)
                zero = zero_certificate(equations)
                self.assertGreaterEqual(margin(equations, x, zero['p'], zero['state']), -PLANNING_TOL)
                oracle = SupportOracle(equations, x, bounds, 20000., 1, 10.)
                try:
                    answer, _ = oracle.get(np.array((.3, .4, .5)))
                    self.assertIsNotNone(answer['point'])
                    self.assertLessEqual(answer['max_violation'], PLANNING_TOL)
                    self.assertGreaterEqual(margin(equations, x, answer['point']*bounds*sign,
                                                   answer['state']), -PLANNING_TOL)
                finally:
                    oracle.close()

    def test_3d_current_faces_certify_without_scheme_enumeration(self):
        network = FourBus(load_nodes=(1, 2, 3))
        equations = PortPhysics(network, (1, 1, 1))
        bounds = port_bounds(network)
        with patch.object(MasterProblem, 'exclude', side_effect=AssertionError('No enumeration')):
            bounds, axis, total, seed = initial_bounds(equations, 20000., 1, 10., bounds=bounds)
            row = certify_scheme(equations, seed['x'], 1, bounds, axis, total,
                budget=20000., tau=.005, epsilon_geom=GEOMETRY_TOL, threads=1,
                support_seconds=10., max_iterations=100, verbose=False)
        self.assertTrue(row['certified'])
        self.assertEqual(row['final_faces'].shape[1], 4)
        self.assertTrue(all('normal_3' in face for face in row['face_history']))
        self.assertEqual(row['initial_support_calls'], 4)
        for certificate in row['certificates']:
            self.assertGreaterEqual(margin(equations, seed['x'], certificate['p'], certificate['state']),
                                    -PLANNING_TOL)


if __name__ == '__main__':
    unittest.main()
