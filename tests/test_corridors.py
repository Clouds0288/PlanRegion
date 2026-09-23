"""Corridor state, type selection and investment accounting."""
from dataclasses import replace
import unittest

import numpy as np

from Network.case33bw import Case33
from Network.four_bus_five_corridor import FourBus
from model import PlanningEquations, PlanningModel


class CorridorTests(unittest.TestCase):
    def test_no_extra_binary_variables(self):
        for network, count in ((FourBus(), 15), (Case33(candidate_count=8), 16),
                               (Case33(candidate_count=16), 32)):
            with self.subTest(network=network.name, count=count):
                problem = PlanningModel(PlanningEquations(network, 'linear'), threads=1)
                with problem.model:
                    problem.model.update()
                    self.assertEqual(problem.model.NumBinVars, count)

    def test_two_indices_and_vector_share_the_same_variables(self):
        fixed = FourBus().design(FourBus().initial_plan)
        for network in (FourBus(), Case33(candidate_count=8), Case33(candidate_count=16), fixed):
            equations = PlanningEquations(network, 'linear')
            problem = PlanningModel(equations, threads=1)
            with problem.model:
                problem.model.update()
                self.assertEqual(list(problem.x), equations.variable_keys)
                self.assertEqual(problem.x_vector.shape, (len(equations.cost),))
                self.assertEqual(problem.model.NumBinVars, len(equations.cost))
                for i, key in enumerate(equations.variable_keys):
                    self.assertTrue(problem.x[key].sameAs(problem.x_vector[i].item()))
                if equations.variable_keys:
                    key = equations.variable_keys[0]
                    problem.x[key].UB = 0.
                    problem.model.update()
                    self.assertEqual(float(problem.x_vector[0].UB), 0.)

    def test_corridor_constraints_use_only_their_own_x_slice(self):
        for network in (FourBus(), Case33(candidate_count=8), Case33(candidate_count=16)):
            equations = PlanningEquations(network, 'linear')
            problem = PlanningModel(equations, threads=1)
            with problem.model:
                problem.model.update()
                for corridor in network.corridors:
                    name = ('required_' if corridor.must_use else 'selectable_') + corridor.id
                    constraint = problem.model.getConstrByName(name)
                    if corridor.must_use and len(corridor.types) == 1:
                        self.assertIsNone(constraint)
                        continue
                    row = problem.model.getRow(constraint)
                    actual = {row.getVar(i).VarName: row.getCoeff(i) for i in range(row.size())}
                    expected = {problem.x[corridor.id, line.id].VarName: 1.
                                for line in corridor.types}
                    self.assertEqual(actual, expected)
                    self.assertEqual(constraint.Sense, '=' if corridor.must_use else '<')
                    self.assertEqual(constraint.RHS, 1.)

    def test_fixed_plan_exposes_state_and_type_expressions(self):
        network = FourBus()
        plan = {'01': None, '12': 'M', '13': 'H', '02': 'H', '23': None}
        equations = PlanningEquations(network, 'linear')
        problem = PlanningModel(equations, fixed_plan=plan, power=[3., 4., 5.], threads=1)
        with problem.model:
            answer = problem.solve()
            self.assertEqual(answer['status'], 'optimal')
            self.assertEqual(equations.choice(answer['x']), plan)
            self.assertAlmostEqual(answer['objective'], network.design(plan).cost)
            for corridor in network.corridors:
                self.assertAlmostEqual(problem.x.sum(corridor.id, '*').getValue(),
                                       float(plan[corridor.id] is not None))
                for line in corridor.types:
                    self.assertAlmostEqual(problem.x[corridor.id, line.id].X,
                                           float(plan[corridor.id] == line.id))
        self.assertTrue(network.corridors[0].initial_active)
        self.assertFalse(network.corridors[3].initial_active)

    def test_required_corridor_can_reverse_and_cannot_open(self):
        network = FourBus()
        network.corridors = tuple(replace(c, must_use=True) if c.id == '12' else c
                                  for c in network.corridors)
        equations = PlanningEquations(network, 'linear')
        plans = (
            {'01': None, '12': 'L', '13': 'L', '02': 'L', '23': None},
            {'01': 'L', '12': None, '13': 'L', '02': 'L', '23': None},
        )
        for plan, feasible in zip(plans, (True, False)):
            problem = PlanningModel(equations, fixed_plan=plan, power=[3., 4., 5.], threads=1)
            with problem.model:
                answer = problem.solve()
                self.assertEqual(answer is not None, feasible)
                if feasible:
                    self.assertEqual(equations.choice(answer['x']), plan)
                    index = equations.type_indices['12']['L']
                    self.assertLess(answer['state'][equations.P_slice][index], 0.)
                    self.assertLess(answer['state'][equations.Q_slice][index], 0.)

    def test_fixed_type_is_constant_and_its_cost_counts(self):
        network = FourBus().design(FourBus().initial_plan)
        network.corridors = tuple(replace(c, types=(replace(c.types[0], investment_cost=7.),))
                                  for c in network.corridors)
        equations = PlanningEquations(network, 'linear')
        self.assertEqual(len(equations.cost), 0)
        self.assertEqual(equations.investment(np.zeros(0)), 21.)
        for budget, feasible in ((20., False), (21., True)):
            problem = PlanningModel(equations, budget=budget, power=[3., 4., 5.], threads=1)
            with problem.model:
                answer = problem.solve()
                self.assertEqual(answer is not None, feasible)
                if feasible:
                    self.assertEqual(answer['objective'], 21.)
                    self.assertEqual(equations.choice(answer['x']), network.initial_plan)

    def test_existing_normally_open_line_is_not_an_empty_corridor(self):
        network = FourBus()
        network.corridors = tuple(replace(c, existing_type='L') if c.id == '02' else c
                                  for c in network.corridors)
        corridor = next(c for c in network.corridors if c.id == '02')
        self.assertEqual(corridor.existing_type, 'L')
        self.assertIsNone(network.initial_plan['02'])
        self.assertFalse(corridor.must_use)


if __name__ == '__main__':
    unittest.main()
