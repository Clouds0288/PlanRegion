"""直接结果、eta 判定和求解次数的回归；不允许后处理改写负荷或状态。"""
from unittest.mock import patch

import gurobipy as gp
import numpy as np
import pytest

from Network.four_bus_five_corridor import FourBus
from model import PLANNING_TOL, GridPhysics, MasterProblem, SubProblem
from tests.planning_checks import margin


@pytest.mark.parametrize('method', ['linear', 'socp'])
def test_mp_and_query_keep_the_solver_power_and_state(method):
    network = FourBus()
    equations = GridPhysics(network, method)
    solve = MasterProblem.solve
    captured = []

    def record(problem, *args, **kwargs):
        answer = solve(problem, *args, **kwargs)
        assert answer['feasible']
        captured.append((problem.power.X.copy(), problem.state.X.copy(), answer['objective']))
        np.testing.assert_array_equal(answer['state'], problem.state.X)
        return answer

    with patch.object(MasterProblem, 'solve', new=record), \
         patch.object(SubProblem, 'solve', side_effect=AssertionError('redundant SP')):
        problem = MasterProblem(equations, fixed_plan=network.initial_plan, threads=1)
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
    equations = GridPhysics(network, 'socp')
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
        answer = SubProblem(equations, threads=1).solve(x, power)

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
    equations = GridPhysics(network, 'socp')
    with pytest.raises(TimeoutError, match='time limit'):
        SubProblem(equations, threads=1).solve(network.encode_plan(network.initial_plan), np.zeros(3), time_limit=0.)


def test_positive_eta_without_a_valid_cut_stays_unknown():
    network = FourBus()
    equations = GridPhysics(network, 'socp')
    with patch.object(SubProblem, '_separating_cut', side_effect=RuntimeError('Invalid separating cut')) as separate:
        with pytest.raises(RuntimeError, match='Invalid separating cut'):
            SubProblem(equations, threads=1).solve(
                network.encode_plan(network.initial_plan), np.array([10., 10., 10.]))
    separate.assert_called_once()


@pytest.mark.parametrize('x,power,feasible', [
    ([0,1,0,1,0,0,0,1,0,0,0,0,0,0,0], [15.053118409066482,19.12279604388464,11.165188097835491], False),
    ([0,1,0,1,0,0,0,1,0,0,0,0,0,0,0], [0.,22.413542764447403,17.558865929406192], True),
])
def test_socp_new_vertex_path_needs_no_numeric_retry(x, power, feasible):
    equations = GridPhysics(FourBus(), 'socp')
    optimize, calls = gp.Model.optimize, []
    def solve(model, *args, **kwargs):
        model.update()
        calls.append(model.NumQConstrs)
        return optimize(model, *args, **kwargs)
    with patch.object(gp.Model, 'optimize', new=solve):
        answer = SubProblem(equations, threads=1).solve(np.array(x), np.array(power))
    assert answer['feasible'] == feasible
    assert len(calls) == (1 if feasible else 2)
    assert calls[0] > 0
    if feasible:
        assert margin(equations, x, power, answer['state']) >= -PLANNING_TOL
    else:
        assert calls[1] == 0  # 正常取割 LP，不是重试 SOCP。
        cut = answer['cut']
        assert cut[0]+cut[1:4]@power+cut[4:]@x < -1e-9
