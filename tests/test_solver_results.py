"""直接结果、eta 判定和求解次数的回归；不允许后处理改写负荷或状态。"""
from unittest.mock import patch

import gurobipy as gp
import numpy as np
import pytest

from Network.four_bus_five_corridor import FourBus
from Network.case33bw import Case33
from model import PLANNING_TOL, PlanningEquations, PlanningModel, PlanningSP
from tests.planning_checks import margin
from tests.reference import fixed_topology


@pytest.mark.parametrize('method', ['linear', 'socp'])
def test_mp_and_query_keep_the_solver_power_and_state(method):
    network = FourBus()
    equations = PlanningEquations(network, method)
    solve = PlanningModel.solve
    captured = []

    def record(problem, *args, **kwargs):
        answer = solve(problem, *args, **kwargs)
        assert answer['feasible']
        captured.append((problem.power.X.copy(), problem.state.X.copy(), answer['objective']))
        np.testing.assert_array_equal(answer['state'], problem.state.X)
        return answer

    with patch.object(PlanningModel, 'solve', new=record), \
         patch.object(PlanningSP, 'solve', side_effect=AssertionError('redundant SP')):
        problem = PlanningModel(equations, fixed_plan=network.initial_plan, threads=1)
        with problem.model:
            answer = problem.solve(radial_gap_kw=1.)
    assert len(captured) == 1 and answer['p'].sum() > 1.
    np.testing.assert_array_equal(answer['p'], captured[0][0])
    np.testing.assert_array_equal(answer['state'], captured[0][1])
    assert answer['objective'] == captured[0][2] == answer['p'].sum()
    assert margin(equations, answer['x'], answer['p'], answer['state']) >= -PLANNING_TOL


@pytest.mark.parametrize('feasible', [True, False])
def test_sp_uses_eta_and_keeps_the_existing_dual_cut_path(feasible):
    network = FourBus()
    equations = PlanningEquations(network, 'socp')
    if feasible:
        plan = {'01': 'H', '12': None, '13': None, '02': 'L', '23': 'M'}
        power = np.array([47.47692657884853, 14.400743716323246, 7.840840158490067])
    else:
        plan, power = network.initial_plan, np.array([10., 10., 10.])
    x, original = network.encode_plan(plan), power.copy()
    optimize, add_operation = gp.Model.optimize, equations.add_operation
    operations, solves = [], []

    def operation(*args, **kwargs):
        result = add_operation(*args, **kwargs)
        operations.append(result)
        return result

    def record(model, *args, **kwargs):
        optimize(model, *args, **kwargs)
        solves.append(dict(eta=model.getVarByName('violation').X, quality=model.MaxVio,
                           state=operations[0].state.X.copy(), cones=model.NumQConstrs))

    with patch.object(equations, 'add_operation', side_effect=operation), \
         patch.object(gp.Model, 'optimize', new=record):
        answer = PlanningSP(equations, threads=1).solve(x, power)

    assert answer['feasible'] == feasible
    np.testing.assert_array_equal(power, original)
    assert solves[0]['cones'] > 0  # 第一次就是带 eta 的 SOCP，不预求原问题。
    if feasible:
        assert len(solves) == 1
        assert solves[0]['eta'] + solves[0]['quality'] <= PLANNING_TOL
        assert answer['cut'] is None
        np.testing.assert_array_equal(answer['state'], solves[0]['state'])
        assert margin(equations, x, power, answer['state']) >= -PLANNING_TOL
    else:
        assert solves[0]['eta'] > PLANNING_TOL
        assert len(solves) == 2 and solves[1]['cones'] == 0  # 第二次仅为原有取割 LP。
        assert answer['state'] is None
        cut = answer['cut']
        assert cut is not None and cut[0]+cut[1:4]@power+cut[4:]@x < -1e-9


def test_sp_timeout_without_solution_stays_unknown():
    network = FourBus()
    equations = PlanningEquations(network, 'socp')
    answer = PlanningSP(equations, threads=1).solve(
        network.encode_plan(network.initial_plan), np.zeros(3), time_limit=0.)
    assert answer == dict(feasible=False, state=None, cut=None)


@pytest.mark.parametrize('upgrades,power', [
    ({}, [69.10864684326633, 3000.4352915681134, 93.62828475260056]),
    ({'2-3': 'parallel'}, [154.27833243608373, 3500.5976872770457, 240.6587165097658]),
])
def test_sp_returns_an_accurate_raw_boundary_state_without_repair(upgrades, power):
    network = fixed_topology(Case33(upgrade_count=4))
    equations = PlanningEquations(network, 'socp')
    x = network.encode_plan(network.initial_plan | upgrades)
    # 默认数值设置曾在这些点返回近零 eta 但超限的原始误差，导致联合割停滞。
    optimize = gp.Model.optimize
    calls = []

    def record(model, *args, **kwargs):
        calls.append(model.ModelName)
        return optimize(model, *args, **kwargs)

    with patch.object(gp.Model, 'optimize', new=record):
        answer = PlanningSP(equations, threads=1).solve(x, power)
    assert calls == ['planning_SP']
    assert answer['feasible'] and answer['cut'] is None
    assert margin(equations, x, power, answer['state']) >= -PLANNING_TOL


def test_positive_eta_without_a_valid_cut_stays_unknown():
    network = FourBus()
    equations = PlanningEquations(network, 'socp')
    with patch.object(PlanningSP, '_separating_cut', return_value=None) as separate:
        answer = PlanningSP(equations, threads=1).solve(
            network.encode_plan(network.initial_plan), np.array([10., 10., 10.]))
    separate.assert_called_once()
    assert answer == dict(feasible=False, state=None, cut=None)
