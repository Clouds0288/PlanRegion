"""Signed reference flows, reverse sending limits and global cut validity."""
from dataclasses import replace

import numpy as np
import pytest
from gurobipy import GRB
from threadpoolctl import threadpool_limits

from Network.four_bus_five_corridor import FourBus
from model import PlanningEquations, PlanningModel
from tests.reference import dispatch_support
from tests.test_corridors import modified


REVERSE_PLAN = {'01': None, '12': 'L', '13': 'H', '02': 'H', '23': None}


@pytest.fixture(autouse=True)
def single_thread():
    with threadpool_limits(limits=1):
        yield


@pytest.mark.parametrize('method', ['linear', 'socp'])
def test_reference_orientation_does_not_change_physics(method):
    network = FourBus()
    flipped = FourBus()
    flipped = modified(flipped, corridors=tuple(replace(c, endpoints=c.endpoints[::-1])
                                                for c in flipped.corridors))
    answers = []
    for net in (network, flipped):
        equations = PlanningEquations(net, method)
        problem = PlanningModel(equations, fixed_plan=REVERSE_PLAN,
                                power=[3., 4., 5.], threads=1)
        with problem.model:
            answer = problem.solve()
        assert answer['status'] == 'optimal'
        assert answer['feasible']
        assert answer['objective'] == network.tree(network.encode_plan(REVERSE_PLAN)).cost
        assert equations.network.decode_plan(answer['x']) == REVERSE_PLAN
        answers.append((equations, answer))
    original, reversed_ = answers
    np.testing.assert_array_equal(original[1]['x'], reversed_[1]['x'])
    for index, (equations, answer) in enumerate(answers):
        branch = equations.network.type_keys.index(('12', 'L'))
        assert (-1 if index == 0 else 1)*answer['state'][equations.P_slice][branch] > 0.


def test_reverse_sending_capacity_includes_losses():
    network = FourBus()
    network.vmin = np.full(network.n, .5)
    equations = PlanningEquations(network, 'socp')
    problem = PlanningModel(equations, fixed_plan=REVERSE_PLAN, threads=1)
    problem.power.UB = [network.power_limit, 0., 0.]
    with problem.model:
        answer = problem.solve()
    reference = dispatch_support(network.tree(network.encode_plan(REVERSE_PLAN)), 'socp', [1., 0., 0.])
    assert answer['status'] == 'optimal'
    assert abs(answer['objective']-reference['value']) < .002
    k = network.type_keys.index(('12', 'L'))
    p = answer['state'][equations.P_slice][k]
    loss = network.r[k]*answer['state'][equations.ell_slice][k]
    capacity = 35./network.base
    assert p < 0. and loss > 1e-4
    assert abs(-p+loss-capacity) < 1e-7
    assert -p < capacity-1e-4


@pytest.mark.parametrize('method', ['linear', 'socp'])
def test_required_single_type_can_carry_reverse_flow(method):
    network = FourBus()
    network = modified(network, corridors=tuple(replace(c, types=c.types[:1])
                         if c.id == '12' else c for c in network.corridors))
    equations = PlanningEquations(network, method)
    assert network.n_types == 13
    problem = PlanningModel(equations, fixed_plan=REVERSE_PLAN, power=[3., 4., 5.], threads=1)
    with problem.model:
        answer = problem.solve()
    assert answer['status'] == 'optimal'
    assert equations.network.decode_plan(answer['x']) == REVERSE_PLAN
    assert answer['state'][equations.P_slice][network.type_keys.index(('12', 'L'))] < 0.


def test_connected_cycle_is_rejected_even_at_zero_load():
    equations = PlanningEquations(FourBus(), 'linear')
    plan = {'01': 'L', '12': 'L', '13': 'L', '02': 'L', '23': None}
    problem = PlanningModel(equations, fixed_plan=plan, power=np.zeros(3),
                            threads=1)
    with problem.model:
        assert problem.solve() is None
