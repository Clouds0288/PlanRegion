"""直接模型调用的行为回归：即时记录、方案去重、失败直接传播。"""
from io import StringIO
from unittest.mock import patch

import numpy as np
import pytest

import main
from model import GridPhysics, MasterProblem, SubProblem, RemainingRegionModel
from region import RegionState
from monitor import RunMonitor


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


@pytest.mark.parametrize('failed_call', [1, 2, 4])
def test_axis_certificate_is_saved_before_later_mp2_can_time_out(failed_call):
    network = main.FourBus()
    original, calls = MasterProblem.solve, []

    def solve(problem, *args, **kwargs):
        calls.append(problem)
        if len(calls) == failed_call:
            raise main.RegionTimeout('MP timeout')
        answer = original(problem, *args, **kwargs)
        assert answer['feasible']
        return answer

    oracle = SubProblem(GridPhysics(network, 'socp'), threads=1)
    monitor = RunMonitor(stream=StringIO())
    with pytest.raises(main.RegionTimeout, match='MP timeout'), monitor, \
         patch.object(MasterProblem, 'solve', new=solve), patch('main.SubProblem', return_value=oracle):
        main.build_continuous_region(network, 'socp', 20000., [150.]*3, threads=1, progress=monitor)
    assert len(calls) == failed_call
    assert oracle.calls == 0
    assert monitor.state['status'] == 'failed'
    assert sum(e['event'] == 'feasible' for e in monitor.events) == failed_call-1


def test_directional_seeds_of_same_scheme_share_one_region():
    network = main.FourBus()
    seeds = []
    def progress(event, **data):
        if event == 'feasible' and data.get('point_reason', '').startswith('MP2'):
            seeds.append((len(data['records']), data['point'].copy()))
    result = main.build_continuous_region(network, 'socp', 0., [150.]*3, threads=1, progress=progress)
    assert [count for count, _ in seeds] == [1, 1, 1, 1]
    assert len({tuple(p) for _, p in seeds}) > 1
    assert result['status'] == 'certified' and len(result['inner']) == 1
