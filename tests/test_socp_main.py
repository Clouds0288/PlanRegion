"""SOCP 主线：顶点评分、32 次检查、二维切片及独立扫描。"""
from unittest.mock import patch

import numpy as np
import pytest

import main
import continuous
from model import GridPhysics, MasterProblem, RemainingRegionModel, SubProblem
from monitor import RunMonitor
from region import RegionState, contains, halfspaces
from vertify import validate_socp_region


def test_direction_objective_keeps_other_loads_free():
    problem = MasterProblem(GridPhysics(main.FourBus(), 'socp'), budget=20000.,
                            direction=[1., 0., 0.], threads=1)
    with problem.model:
        problem.model.addConstr(problem.power[1] >= 2.)
        answer = problem.solve()
    assert answer['feasible'] and answer['p'][1] >= 2.-1e-8
    assert answer['objective'] == answer['p'][0]
    assert answer['objective'] < answer['p'].sum()-1.
    assert answer['bound'] >= answer['objective']-1e-7


def test_vertex_priority_uses_kw_not_normalized_sum():
    # 新语义：违反量优先，即使它的总负荷较小；相同 (x,p) 的第二轮不得重复求 SP。
    network = main.FourBus(load_nodes=(1, 2))
    x = network.encode_plan(network.initial_plan).astype(int)
    state = RegionState([1., 1.], 2., 0.)
    state.add_scheme(x, network.initial_plan, 0.)
    state.records[tuple(x)]['outer'] = np.array([[.7, .1], [.2, .9]])
    seed = dict(x=x, p=np.zeros(2), bound=2., feasible=True, status='optimal')
    chosen, other = np.r_[-.5, np.zeros(2+len(x))], np.r_[-.2, np.zeros(2+len(x))]
    def score(oracle, choice, power, time_limit=None):
        oracle.calls += 1
        return dict(eta=.008 if power[0] > .5 else .003, feasible=False, state=None,
                    cut=chosen if power[0] > .5 else other)
    oracle = SubProblem(GridPhysics(network, 'socp'), threads=1)
    with patch('continuous.RegionState', return_value=state), \
         patch('continuous.SubProblem', return_value=oracle), \
         patch.object(MasterProblem, 'solve', return_value=seed), \
         patch.object(SubProblem, 'solve', new=score), \
         patch.object(state, 'apply_cut', side_effect=[None, RuntimeError('second local cut')]) as apply, \
         patch.object(RemainingRegionModel, 'solve', side_effect=AssertionError('0.01 hard trigger')):
        with pytest.raises(RuntimeError, match='second local cut'):
            continuous.build_continuous_region(network, 'socp', 20000., [1., 1.], tau=0., threads=1)
    assert oracle.calls == 2
    np.testing.assert_array_equal(apply.call_args_list[0].args[0], chosen)
    np.testing.assert_array_equal(apply.call_args_list[1].args[0], chosen)


def test_axis_bounds_clip_existing_new_and_final_outer():
    state = RegionState([100., 200.], 300., .002)
    state.add_scheme([0], {}, 0.)
    state.tighten_bounds([40., 60.], 80.)
    state.add_scheme([1], {}, 1.)
    for row in state.records.values():
        p = row['outer']*state.bounds
        assert np.all(p <= [40.+1e-9, 60.+1e-9])
        assert np.all(p.sum(axis=1) <= 80.+1e-9)
    state.add_point([0], [[0., 0.], [.4, 0.], [0., .3]])
    for certified in (False, True):
        for row in state.finish(certified)['outer']:
            assert np.all(row['vertices'] <= [40.+1e-9, 60.+1e-9])


def test_socp_flow_agrees_with_independent_fixed_point_queries():
    network = main.FourBus(load_nodes=(1, 2))
    bounds = np.full(2, network.power_limit)
    monitor = RunMonitor()
    directions, upper = [], []
    solve = MasterProblem.solve
    def mp(problem, *args, **kwargs):
        answer = solve(problem, *args, **kwargs)
        directions.append(problem.direction.copy())
        upper.append(answer['bound'])
        return answer
    with patch.object(MasterProblem, 'solve', new=mp):
        result = continuous.build_continuous_region(network, 'socp', 20000., bounds, threads=1, progress=monitor)
    np.testing.assert_array_equal(directions, np.vstack([np.eye(2), np.ones(2)]))
    np.testing.assert_allclose(result['axis_bounds'], np.minimum(bounds, upper[:2]))
    assert result['status'] == 'certified' and result['coverage_bound'] <= 1e-8
    assert continuous.REFINEMENT_CHECKS == 32
    reference = validate_socp_region(network, 20000., 12, result['axis_bounds'], threads=1)
    points = (np.indices((12, 12)).reshape(2, -1).T+.5)*result['axis_bounds']/12
    inner = np.zeros(len(points), dtype=bool)
    outer = inner.copy()
    for row in result['inner']:
        inner |= contains(points, halfspaces(row['vertices']))
    for row in result['outer']:
        outer |= contains(points, halfspaces(row['vertices']))
    truth = reference['states'].ravel() == 1
    assert not np.any(inner & ~truth)
    assert not np.any(truth & ~outer)
    assert len(monitor.state['schemes']) >= 2  # 两个不同网架已可覆盖当前二维案例。
    assert len({tuple(f['patch']['sp_point']['p']) for f in monitor.history
                if f['patch']['event'] == 'point' and 'sp_point' in f['patch']}) > 5


def test_interval_checks_actual_new_sp_after_batch():
    # 三维保留为核心回归：该预算确实触发 32 间隔，二维默认可能先把顶点查完。
    monitor = RunMonitor()
    result = continuous.build_continuous_region(main.FourBus(), 'socp', 20000., [142.5]*3,
                                          threads=1, progress=monitor)
    first = next(f['patch'] for f in monitor.history if f['patch']['event'] == 'residual_start')
    assert first['sp_since_global'] >= 32
    assert first['sp_since_global'] < 96
    assert result['coverage_bound'] <= 1e-8


def test_default_run_finishes_socp_before_ac(tmp_path):
    import gzip
    import json
    from region import build_sequential_region
    from monitor import RunMonitor
    build, scan = build_sequential_region, main.scan_reference
    calls = []
    def region(*args, **kwargs):
        result = build(*args, **kwargs)
        calls.append(('region', kwargs['mode'], result['status']))
        return result
    def validate(*args, **kwargs):
        assert calls == [('region', 1, 'certified')]*4
        calls.append(('scan',))
        return scan(*args, **kwargs)
    output = tmp_path/'monitor.json.gz'
    with patch('region.build_sequential_region', side_effect=region), \
         patch('main.scan_reference', side_effect=validate):
        result = main.run(main.FourBus(load_nodes=(1, 2)), divisions=4, output=output,
                          show_ui=False, threads=1, scan_workers=1, scan_output=tmp_path/'scans')
    with gzip.open(output, 'rt', encoding='utf-8') as stream:
        data = json.load(stream)
    assert calls[-1] == ('scan',)
    assert data['version'] == 4
    assert data['validation_state']['validation']['metrics']['inner']['fr_percent'] == 0.
    assert result['status'] == 'certified'
    assert sorted(p.name for p in tmp_path.iterdir()) == ['monitor.json.gz', 'scans']
    restored = RunMonitor()
    restored.load_recording(output)
    assert restored.state['result']['certified']
