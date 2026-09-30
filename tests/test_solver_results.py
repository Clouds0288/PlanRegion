"""直接结果、eta 判定和求解次数的回归；不允许后处理改写负荷或状态。"""
from unittest.mock import patch

import gurobipy as gp
import numpy as np
import pytest

from Network.four_bus_five_corridor import FourBus
from tests.legacy_case33 import Case33
from model import PLANNING_TOL, GridPhysics, MasterProblem, SubProblem, RemainingRegionModel
from tests.planning_checks import margin
from tests.reference import fixed_topology


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


@pytest.mark.parametrize('upgrades,power', [
    ({}, [69.10864684326633, 3000.4352915681134, 93.62828475260056]),
    ({'2-3': 'parallel'}, [154.27833243608373, 3500.5976872770457, 240.6587165097658]),
])
def test_sp_returns_an_accurate_raw_boundary_state_without_repair(upgrades, power):
    network = fixed_topology(Case33(upgrade_count=4))
    equations = GridPhysics(network, 'socp')
    x = network.encode_plan(network.initial_plan | upgrades)
    # 默认数值设置曾在这些点返回近零 eta 但超限的原始误差，导致联合割停滞。
    optimize = gp.Model.optimize
    calls = []

    def record(model, *args, **kwargs):
        calls.append(model.ModelName)
        return optimize(model, *args, **kwargs)

    with patch.object(gp.Model, 'optimize', new=record):
        answer = SubProblem(equations, threads=1).solve(x, power)
    assert calls == ['planning_SP']
    assert answer['feasible'] and answer['cut'] is None
    assert margin(equations, x, power, answer['state']) >= -PLANNING_TOL


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


def test_remaining_region_always_contains_full_physics_and_returns_a_physical_witness():
    network = FourBus(load_nodes=(1, 2))
    equations = GridPhysics(network, 'socp')
    problem = RemainingRegionModel(equations, 20000., np.full(2, 300.), 300., [], [], .005, threads=1)
    with problem.model:
        problem.model.update()
        assert problem.model.NumQConstrs > 0
        assert problem.problem.state.shape == (len(equations.y_lb_global),)
        answer = problem.solve(1e-8)
        assert not answer['complete'] and answer['feasible']
        assert margin(equations, answer['x'], answer['p'], problem.problem.state.X) >= -PLANNING_TOL


def test_remaining_region_excludes_nonphysical_points_without_any_joint_cuts():
    network = FourBus(load_nodes=(1, 2))
    equations = GridPhysics(network, 'socp')
    problem = RemainingRegionModel(equations, 20000., np.full(2, 300.), 600., [], [], .005, threads=1)
    with problem.model:
        problem.problem.power.LB = [250., 0.]
        answer = problem.solve(1e-8)
        assert answer == dict(complete=True, bound=None, x=None, p=None, feasible=False)


def test_remaining_region_rejects_linear_equations():
    with pytest.raises(ValueError, match='complete SOCP'):
        RemainingRegionModel(GridPhysics(FourBus(), 'linear'), 20000., np.full(3, 300.),
                             300., [], [], .005, threads=1)
