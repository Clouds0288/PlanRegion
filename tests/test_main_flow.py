"""直接模型调用的行为回归：即时记录、方案去重、失败直接传播。"""

import numpy as np
import pytest

import main
from model import GridPhysics, MasterProblem


def test_mp_accepts_cuts_and_incumbent_without_a_query_wrapper():
    network = main.FourBus()
    equations = GridPhysics(network, 'socp')
    cut = np.r_[20., -np.ones(3), np.zeros(network.n_types)]
    problem = MasterProblem(equations, cuts=[cut], threads=1)
    with problem.model:
        answer = problem.solve()
    assert answer['feasible']
    assert answer['p'].sum() == pytest.approx(20., abs=1e-6)

    problem = MasterProblem(equations, min_total=20.-1e-6, cuts=[cut], threads=1)
    with problem.model:
        cheapest = problem.solve(incumbent=answer)
    assert cheapest['feasible'] and cheapest['status'] == 'optimal'
    assert cheapest['p'].sum() >= 20.-2e-6
    assert cheapest['objective'] <= network.cost@answer['x']+1e-7
