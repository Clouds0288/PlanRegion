"""SOCP 主线：方向上界、真实负荷选点、全局覆盖及 AC 调用次序。"""
from unittest.mock import patch

import numpy as np
import pytest

import main
from model import GridPhysics, MasterProblem, RemainingRegionModel, SubProblem
from plot import sample_region
from region import RegionState, contains, halfspaces


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
    network, bounds = main.Concept5(), np.array([100., 10.])
    x = network.encode_plan(network.initial_plan).astype(int)
    state = RegionState(bounds, 110., 0.)
    state.add_scheme(x, network.initial_plan, 0.)
    state.records[tuple(x)]['outer'] = np.array([[.7, .1], [.2, .9]])
    seed = dict(x=x, p=np.zeros(2), bound=100., feasible=True, status='optimal')
    with patch('main.RegionState', return_value=state), \
         patch.object(MasterProblem, 'solve', return_value=seed), \
         patch.object(SubProblem, 'solve', side_effect=RuntimeError('selected')) as check:
        with pytest.raises(RuntimeError, match='selected'):
            main.build_continuous_region(network, 'socp', 4., bounds, tau=0., threads=1)
    np.testing.assert_array_equal(check.call_args.args[1], [70., 1.])


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
    network = main.FourBus()
    bounds = np.full(3, network.power_limit)
    events, directions, upper = [], [], []
    solve = MasterProblem.solve

    def mp(problem, *args, **kwargs):
        assert kwargs.get('incumbent') is None
        answer = solve(problem, *args, **kwargs)
        directions.append(problem.direction.copy())
        upper.append(answer['bound'])
        return answer

    def progress(event, **data):
        events.append(event)
        if data.get('answer') is not None:
            assert 'state' not in data['answer'] and 'cut' not in data['answer']
        if event == 'point' and data['point_reason'].startswith('已知网架'):
            point = data['point']/bounds
            assert not any(contains([point], halfspaces(row['inner']))[0] for row in data['records'])

    with patch.object(MasterProblem, 'solve', new=mp):
        result = main.build_continuous_region(network, 'socp', 20000., bounds, threads=1, progress=progress)
    np.testing.assert_array_equal(directions, np.vstack([np.eye(3), np.ones(3)]))
    np.testing.assert_allclose(result['axis_bounds'], np.minimum(bounds, upper[:3]))
    assert result['status'] == 'certified' and result['coverage_bound'] <= 1e-8
    assert events.count('mp_start') == 4 and 'residual_start' in events
    assert events.index('point') > [i for i, e in enumerate(events) if e == 'query_end'][-1]
    # 当前真实负荷点由完整、固定 p 的整数 SOCP 独立判断，不使用割或内域。
    points = np.array([[5., 5., 5.], [70., 10., 5.], [0., 30., 5.], [0., 0., 40.], [80., 20., 20.]])
    labels = sample_region(result, points, bounds)
    assert np.all(labels != 0)
    for point, label in zip(points, labels):
        problem = MasterProblem(GridPhysics(network, 'socp'), power=point, budget=20000., threads=1)
        with problem.model:
            answer = problem.solve()
        assert (label == 1) == (answer is not None)


def test_default_run_finishes_socp_before_ac(tmp_path):
    build, validate = main.build_continuous_region, main.validate_ac_region
    calls = []
    def region(*args, **kwargs):
        result = build(*args, **kwargs)
        calls.append(('region', args[1], result['status']))
        return result
    def ac(*args, **kwargs):
        assert calls == [('region', 'socp', 'certified')]
        calls.append(('ac',))
        return validate(*args, **kwargs)
    with patch('main.build_continuous_region', side_effect=region), patch('main.validate_ac_region', side_effect=ac):
        result = main.run(main.FourBus(), budgets=[20000.], divisions=3, output=tmp_path, show_ui=False, threads=1)
    assert calls[-1] == ('ac',)
    assert set(result.metadata['seconds']) == {'socp', 'ac'}
    assert {r['method'] for r in result.summary} == {'socp', 'ac'}
    assert (tmp_path/'result.npz').exists()
