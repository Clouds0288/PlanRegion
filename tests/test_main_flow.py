"""直接模型调用的行为回归：即时保存证书、方案去重、显式 SP 补救。"""
from unittest.mock import patch

import numpy as np
import pytest

import main
from model import PlanningEquations, PlanningModel, PlanningSP, RemainingRegionModel
from region import RegionState, contains, halfspaces


def test_mp_accepts_cuts_and_incumbent_without_a_query_wrapper():
    network = main.FourBus()
    equations = PlanningEquations(network, 'linear')
    cut = np.r_[20., -np.ones(3), np.zeros(network.n_types)]
    problem = PlanningModel(equations, cuts=[cut], threads=1)
    with problem.model:
        answer = problem.solve()
    assert answer['feasible']
    assert answer['p'].sum() == pytest.approx(20., abs=1e-6)

    problem = PlanningModel(equations, min_total=20.-1e-6, cuts=[cut], threads=1)
    with problem.model:
        cheapest = problem.solve(incumbent=answer)
    assert cheapest['feasible'] and cheapest['status'] == 'optimal'
    assert cheapest['p'].sum() >= 20.-2e-6
    assert cheapest['objective'] <= network.cost@answer['x']+1e-7


@pytest.mark.parametrize('repair', [False, True])
def test_mp2_certificate_is_saved_before_mp1_can_time_out(repair):
    network = main.FourBus()
    original, calls = PlanningModel.solve, []

    def solve(problem, *args, **kwargs):
        calls.append(problem)
        if len(calls) == 2:
            raise main.RegionTimeout()
        answer = original(problem, *args, **kwargs)
        assert answer['feasible']
        if repair:
            answer.update(feasible=False, status='unknown')
        return answer

    oracle = PlanningSP(PlanningEquations(network, 'linear'), threads=1)
    with patch.object(PlanningModel, 'solve', new=solve), \
         patch('main.PlanningSP', return_value=oracle):
        result = main.build_continuous_region(network, 'linear', np.inf, [100.]*3, threads=1)
    assert len(calls) == 2
    assert result['status'] == 'time_limit' and result['inner']
    assert oracle.calls == int(repair)
    assert any(contains([result['max_point']], halfspaces(row['vertices']))[0] for row in result['inner'])


def test_identical_mp1_mp2_scheme_is_refined_only_once():
    network = main.FourBus()
    original, answers = PlanningModel.solve, []

    def solve(problem, *args, **kwargs):
        if not answers:
            answers.append(original(problem, *args, **kwargs))
        return answers[0].copy()

    with patch.object(PlanningModel, 'solve', new=solve), \
         patch.object(RegionState, 'next_point', return_value=None) as select, \
         patch.object(PlanningSP, 'solve', side_effect=AssertionError('redundant SP')), \
         patch.object(RemainingRegionModel, 'solve', return_value=dict(complete=False, bound=1., x=None)):
        result = main.build_continuous_region(network, 'linear', np.inf, [100.]*3, threads=1)
    assert select.call_count == 1
    assert result['status'] == 'unknown'  # 选点结束不能代替全局覆盖证明。
    assert len(result['inner']) == 1
