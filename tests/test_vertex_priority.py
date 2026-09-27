"""Independent SP scores, cache semantics, and a choice where eta reverses load order."""
from unittest.mock import patch
import unittest

import numpy as np

from experiments.fourbus_vertex_priority import ScoredSubProblem, build_vertex_priority_region
from experiments.fourbus_outer_scan import candidate_eta
from model import GridPhysics, MasterProblem, SubProblem, RemainingRegionModel, PLANNING_TOL
from Network.four_bus_five_corridor import FourBus
from region import RegionState
from tests.planning_checks import margin


class VertexPriorityTests(unittest.TestCase):
    def test_score_matches_independent_primal_and_original_certificate(self):
        network = FourBus()
        equations = GridPhysics(network, 'socp')
        x = network.encode_plan(network.initial_plan).astype(int)
        scored = ScoredSubProblem(equations, threads=1)
        for power in (np.zeros(3), np.array([10., 20., 40.])):
            with self.subTest(power=power):
                answer = scored.solve(x, power)
                independent = candidate_eta(equations, x, power)
                original = SubProblem(equations, threads=1).solve(x, power)
                self.assertAlmostEqual(answer['eta'], independent, delta=1e-7)
                self.assertEqual(answer['feasible'], original['feasible'])
                if answer['feasible']:
                    self.assertGreaterEqual(margin(equations, x, power, answer['state']), -PLANNING_TOL)
                else:
                    cut = answer['cut']
                    self.assertLess(cut[0]+cut[1:4]@power+cut[4:]@x, -1e-9)
                    np.testing.assert_allclose(cut, original['cut'], atol=1e-8, rtol=1e-7)

    def test_cache_reuses_only_identical_parameters(self):
        network = FourBus()
        equations = GridPhysics(network, 'socp')
        x = network.encode_plan(network.initial_plan).astype(int)
        oracle = ScoredSubProblem(equations, threads=1)
        power = np.array([10., 20., 40.])
        first = oracle.solve(x, power)
        self.assertIs(oracle.solve(x.copy(), power.copy()), first)
        self.assertEqual((oracle.calls, oracle.cache_hits), (1, 1))
        oracle.solve(x, power+np.array([0., 0., 1e-5]))
        self.assertEqual((oracle.calls, oracle.cache_hits), (2, 1))

    def test_larger_violation_can_select_the_lower_load_vertex(self):
        # 合成几何与评分只检验决策，不冒充真实电气证书。
        network = FourBus()
        x = network.encode_plan(network.initial_plan).astype(int)
        state = RegionState(np.ones(3), 3., 0.)
        state.add_scheme(x, network.initial_plan, network.cost@x)
        state.records[tuple(x)]['outer'] = np.array([[.7, .1, .1], [.2, .6, .6]])
        seed = dict(x=x, p=np.zeros(3), bound=3., feasible=True, status='optimal')
        chosen = np.r_[-.5, np.zeros(3+len(x))]
        other = np.r_[-.2, np.zeros(3+len(x))]

        def score(oracle, choice, power, time_limit=None):
            return dict(eta=.5 if power[0] > .5 else .2, feasible=False, state=None,
                        cut=chosen if power[0] > .5 else other)

        with patch('experiments.fourbus_vertex_priority.RegionState', return_value=state), \
             patch.object(MasterProblem, 'solve', return_value=seed), \
             patch.object(ScoredSubProblem, 'solve', new=score), \
             patch.object(state, 'apply_cut', side_effect=RuntimeError('selected')) as apply:
            with self.assertRaisesRegex(RuntimeError, 'selected'):
                build_vertex_priority_region(network, 20000., np.ones(3), tau=0., threads=1)
        np.testing.assert_array_equal(apply.call_args.args[0], chosen)

    def test_eta_trigger_requires_the_largest_pending_score_to_be_small(self):
        for maximum, should_trigger in ((.5, False), (.008, True)):
            with self.subTest(maximum=maximum):
                network = FourBus()
                x = network.encode_plan(network.initial_plan).astype(int)
                state = RegionState(np.ones(3), 3., 0.)
                state.add_scheme(x, network.initial_plan, network.cost@x)
                state.records[tuple(x)]['outer'] = np.array([[.7, .1, .1], [.2, .6, .6]])
                seed = dict(x=x, p=np.zeros(3), bound=3., feasible=True, status='optimal')
                cut = np.r_[-.5, np.zeros(3+len(x))]
                scores = [dict(eta=maximum, feasible=False, state=None, cut=cut),
                          dict(eta=.003, feasible=False, state=None, cut=cut),
                          RuntimeError('another local batch')]
                complete = dict(complete=True, bound=0., x=None, p=None, feasible=False)
                with patch('experiments.fourbus_vertex_priority.RegionState', return_value=state), \
                     patch.object(MasterProblem, 'solve', return_value=seed), \
                     patch.object(ScoredSubProblem, 'solve', side_effect=scores), \
                     patch.object(state, 'apply_cut'), \
                     patch.object(RemainingRegionModel, 'solve', return_value=complete) as search:
                    if should_trigger:
                        result = build_vertex_priority_region(network, 20000., np.ones(3),
                                                              tau=0., threads=1, eta_trigger=.01)
                        self.assertEqual(result['trace'][-1]['trigger'], ['eta'])
                        self.assertEqual(search.call_count, 1)
                    else:
                        with self.assertRaisesRegex(RuntimeError, 'another local batch'):
                            build_vertex_priority_region(network, 20000., np.ones(3),
                                                         tau=0., threads=1, eta_trigger=.01)
                        search.assert_not_called()


if __name__ == '__main__':
    unittest.main()
