"""SOCP 主线：逐网架认证、二维切片及独立扫描。"""
from unittest.mock import patch

import numpy as np

import main
from model import GridPhysics, MasterProblem
from monitor import RunMonitor
from region import RegionState, build_sequential_region, contains, halfspaces
from vertify import scan_ac_reference


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


def test_socp_flow_agrees_with_independent_fixed_point_queries(tmp_path):
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
        result = build_sequential_region(network, budget=20000., monitor=monitor, threads=1)
    np.testing.assert_array_equal(directions, np.vstack([np.eye(2), np.ones(2)]))
    np.testing.assert_allclose(result['axis_bounds'], np.minimum(bounds, upper[:2]))
    assert result['status'] == 'certified' and result['coverage_bound'] <= 1e-8
    reference = scan_ac_reference(network, 20000.,
        dict(axis_lower=np.zeros(2), bounds=result['axis_bounds'], shape=(12, 12)),
        tmp_path/'ac.npz', workers=1, mode=0)
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


def test_default_run_finishes_socp_before_ac(tmp_path):
    import gzip
    import json
    from region import build_sequential_region
    from monitor import RunMonitor
    build, scan = build_sequential_region, main.scan_ac_reference
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
         patch('main.scan_ac_reference', side_effect=validate):
        result = main.run(main.FourBus(load_nodes=(1, 2)), divisions=4, output=output,
                          show_ui=False, threads=1, scan_workers=1, scan_output=tmp_path/'scans')
    with gzip.open(output, 'rt', encoding='utf-8') as stream:
        data = json.load(stream)
    assert calls[-1] == ('scan',)
    assert data['version'] == 4
    assert data['validation_state']['validation']['method'] == 'ac_socp_grid_v4'
    assert set(data['validation_state']['validation']['comparisons']) == {'result_ac', 'result_socp', 'socp_ac'}
    assert result['status'] == 'certified'
    assert sorted(p.name for p in tmp_path.iterdir()) == ['monitor.json.gz', 'monitor_comparison', 'scans']
    restored = RunMonitor()
    restored.load_recording(output)
    assert restored.state['result']['certified']
