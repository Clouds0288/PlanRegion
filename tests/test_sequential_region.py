"""逐网架主线：真实数值证书、共享割、单网架调度、超时与原生回放。"""
from unittest.mock import patch

import gurobipy as gp
from gurobipy import GRB
import numpy as np
import pytest
from threadpoolctl import threadpool_limits

from main import (build_sequential_region, stage_candidates, ray_support,
                                          register_power, ray_gain, coverage_halfspaces)
from Network.four_bus_five_corridor import FourBus
from Network.case33bw import Case33
from model import GridPhysics, MasterProblem, SubProblem, PLANNING_TOL
from monitor import RunMonitor, NativeWindow, _union
from region import RegionState, contains, halfspaces
from tests.planning_checks import margin


@pytest.fixture(scope='module', params=('fourbus', 'case33'))
def experiment(request, tmp_path_factory):
    network, budget, focus = ((FourBus(load_nodes=(1, 2)), 20000., 0) if request.param == 'fourbus'
                              else (Case33(), 7, 3))
    monitor = RunMonitor(output=tmp_path_factory.mktemp(request.param)/'monitor.json.gz')
    with threadpool_limits(limits=1):
        result = build_sequential_region(network, budget=budget, monitor=monitor,
                                         seconds=60., threads=4, numeric_focus=focus)
    monitor.save()
    return network, budget, monitor, result


def test_only_active_scheme_gets_sp_and_stopping_needs_global_certificate(experiment):
    _, _, monitor, result = experiment
    assert result['certified'] and result['coverage_bound'] <= 1e-8
    assert sum(row['patch']['event'] == 'scheme_start' for row in monitor.history) >= 2
    checked = set()
    for index, item in enumerate(monitor.history):
        event = item['patch'].get('event')
        if event == 'point':
            state = monitor.frame(index)
            point = state['sp_point']
            assert point['scheme'] == state['active_scheme']
            key = point['scheme'], tuple(point['p'])
            assert key not in checked
            checked.add(key)
        elif event == 'scheme_end':
            state = monitor.frame(index)
            assert not state.get('coverage_complete', False)
            if state['stage_reason'] == 'area_stagnation':
                assert state['small_cuts'] == 3
    assert len(checked) == result['counts']['sp']
    assert result['counts']['global_search'] >= 1


def test_area_ratio_uses_frozen_union_and_preserves_cut_away_candidates(experiment):
    _, _, monitor, _ = experiment
    for index, item in enumerate(monitor.history):
        event = item['patch'].get('event')
        if event == 'cut_measure':
            state, previous = monitor.frame(index), monitor.frame(index-2)
            scheme = state['active_scheme']
            before = _union([previous['schemes'][scheme]['outer']])
            after = _union([state['schemes'][scheme]['outer']])
            inner = _union([row['inner'] for row in previous['schemes'].values()])
            removed = before.area-after.area if state['stage'] == 1 else before.difference(inner).area-after.difference(inner).area
            denominator = before.area if state['stage'] == 1 else inner.area
            if denominator > 0.:
                assert state['area_ratio'] == pytest.approx(removed/denominator, abs=1e-9)
    # 历史割来源留在回放，切掉的顶点不再进入全网架定点认证队列。
    assert not any(item['patch'].get('event') == 'handoff_candidate' for item in monitor.history)


def test_shared_cuts_valid_for_complete_budgeted_physical_model(experiment):
    network, budget, monitor, _ = experiment
    equations = GridPhysics(network, 'socp')
    with threadpool_limits(limits=1):
        for row in monitor.state['cut_history'].values():
            cut = np.asarray(row['cut'])
            problem = MasterProblem(equations, budget=budget, threads=4)
            with problem.model as model:
                expression = cut[0]+gp.quicksum(float(c)*v for c, v in zip(cut[1:3], problem.loads.values()))
                expression += gp.quicksum(float(c)*v for c, v in zip(cut[3:], problem.choices.values()))
                model.setObjective(expression, GRB.MINIMIZE)
                model.Params.TimeLimit = 20.
                model.optimize()
                assert model.Status == GRB.OPTIMAL
                assert model.ObjBound >= -PLANNING_TOL


def test_red_points_have_global_proof_and_replay_does_not_leak_future(experiment):
    _, _, monitor, result = experiment
    assert set(result) == {'status', 'certified', 'axis_bounds', 'coverage_bound',
                           'counts', 'timing', 'inner', 'outer'}
    assert set(result['counts']) == {'initial', 'sp', 'cuts', 'ray', 'global_search'}
    assert set(result['timing']) == {'total_seconds'}
    restored = RunMonitor()
    restored.load_recording(monitor.output)
    assert restored.history == monitor.history
    boundary = next(i for i, item in enumerate(monitor.history) if item['patch'].get('event') == 'boundary_start')
    assert not restored.frame(boundary).get('coverage_complete', False)
    assert restored.frame(len(monitor.history)-1)['coverage_complete']
    window = NativeWindow(restored)
    try:
        window.root.withdraw()
        window.seek(boundary)
        assert not restored.frame(window.index).get('rejected_points')
        ray = next(i for i in range(boundary+1, len(monitor.history)) if monitor.history[i]['patch'].get('event') == 'ray_end')
        window.seek(ray)
        assert 'λ=' in window.status.get()
        window.seek(0)
        assert not restored.frame(window.index).get('rejected_points')
    finally:
        window.close()


def test_ray_returns_same_scheme_physical_state(experiment):
    network, budget, monitor, _ = experiment
    equations = GridPhysics(network, 'socp')
    with threadpool_limits(limits=1):
        for index, item in enumerate(monitor.history):
            if item['patch'].get('event') != 'ray_start':
                continue
            state = monitor.frame(index)
            ray = state['ray']
            x = np.asarray(state['schemes'][ray['scheme']]['x'])
            answer = ray_support(equations, budget, x, np.asarray(ray['anchor']), np.asarray(ray['target']),
                                 threads=4, time_limit=20., numeric_focus=state['numeric_focus'])
            expected = np.asarray(ray['anchor'])+answer['ray_fraction']*(np.asarray(ray['target'])-ray['anchor'])
            np.testing.assert_allclose(answer['p']/network.base, expected/network.base,
                                       atol=PLANNING_TOL, rtol=0.)
            assert margin(equations, x, answer['p'], answer['state']) >= -PLANNING_TOL


def test_difference_uses_true_faces_and_keeps_interior_gaps():
    region = RegionState([1., 1.], 2., .005)
    region.add_scheme([1], {}, 0.)
    region.add_point([1], [[0., 0.], [.51, .23], [1., 0.]])
    candidates = stage_candidates(region, [1])
    np.testing.assert_array_equal(candidates, [[0., 1.], [1., 1.]])
    region.add_scheme([0], {}, 0.)
    region.add_point([0], [[0., 1.], [1., 1.], [.5, .8]])
    # 其他网架覆盖了顶点，也不会改变本网架的切割/补边界目标。
    np.testing.assert_array_equal(stage_candidates(region, [1]), candidates)
    assert region.covering_schemes([[.5, .5]]) == [None]


def test_each_new_scheme_finishes_global_vertex_rays_before_sp_or_cut(experiment):
    from region import initial_polytope
    _, _, monitor, _ = experiment
    completed, targets = set(), []
    sweeping = None
    for index, item in enumerate(monitor.history):
        event = item['patch']['event']
        state = monitor.frame(index)
        key = state.get('active_scheme')
        if event == 'initial_sweep_start':
            assert key not in completed and sweeping is None
            sweeping, targets = key, []
            bounds = np.asarray(state['bounds'])
            vertices = initial_polytope(bounds, state['total_bound'], state['axis_bounds'])*bounds
            initial_inner = halfspaces(np.asarray(state['schemes'][key]['inner'])/bounds)
        elif event == 'ray_start' and sweeping is not None:
            target = np.asarray(state['ray']['target'])
            assert key == sweeping and state['phase'] == '首轮射线认证'
            assert np.min(np.max(abs(vertices-target), axis=1)) < 1e-7
            assert not any(np.max(abs(target-p)) <= 1e-4 for p in targets)
            targets.append(target)
        elif event == 'initial_sweep_end':
            final_inner = halfspaces(np.asarray(state['schemes'][key]['inner'])/bounds)
            for vertex in vertices:
                assert (contains([vertex/bounds], final_inner)[0]
                        or any(np.max(abs(vertex-p)) <= 1e-4 for p in targets))
            assert len(targets) >= np.count_nonzero(~contains(vertices/bounds, initial_inner))
            completed.add(key)
            sweeping = None
        elif event in ('point', 'cut_start', 'ray_skip'):
            assert sweeping is None and key in completed
    assert completed == set(monitor.state['schemes'])


def test_time_limit_retains_safe_outer_and_never_marks_unknown_red():
    monitor = RunMonitor()
    result = build_sequential_region(FourBus(load_nodes=(1, 2)), budget=20000., monitor=monitor, seconds=0.)
    assert result['status'] == 'time_limit' and not result['certified']
    assert result['outer'] and 'rejected_points' not in result
    monitor = RunMonitor()
    with threadpool_limits(limits=1), patch('main.RemainingRegionModel.solve', side_effect=TimeoutError):
        result = build_sequential_region(FourBus(load_nodes=(1, 2)), budget=20000., monitor=monitor, seconds=30.)
    assert result['status'] == 'time_limit' and not result['certified']
    assert 'unknown_points' not in monitor.state and 'rejected_points' not in result
    assert any(item['patch'].get('event') == 'boundary_end' for item in monitor.history)


def test_point_representatives_are_fixed_and_do_not_chain_or_round():
    powers = []
    original = np.array([.000049, 7.])
    assert register_power(original, powers, 1e-4) == 0
    assert register_power(np.array([.000051, 7.]), powers, 1e-4) == 0
    assert register_power(original+[.00009, 0.], powers, 1e-4) == 0
    # 与上一近点很近，但已远离原代表，不能沿链扩大合并范围。
    assert register_power(original+[.00018, 0.], powers, 1e-4) == 1
    original[:] = 0.
    np.testing.assert_array_equal(powers[0], [.000049, 7.])


def test_sp_frontier_and_certification_share_distinct_representatives(experiment):
    _, _, monitor, result = experiment
    local, used_cuts = {}, set()
    for index, item in enumerate(monitor.history):
        event = item['patch'].get('event')
        assert event not in ('certification_start', 'certification_end', 'handoff_candidate')
        if event not in ('point', 'cut'):
            continue
        state = monitor.frame(index)
        if event == 'point':
            point = state['sp_point']
            power = np.asarray(point['p'])
            previous = local.setdefault(point['scheme'], [])
            assert all(np.max(np.abs(power-other)) > state['point_tol'] for other in previous)
            # 包括低维线段和本批刚刚扩大的凸包，均不应再调用 SP。
            inner = np.asarray(state['schemes'][point['scheme']]['inner'])/state['bounds']
            assert not contains([power/state['bounds']], halfspaces(inner))[0]
            previous.append(power)
        else:
            point = state['sp_point']
            key = point['scheme'], tuple(point['p'])
            assert key not in used_cuts
            used_cuts.add(key)
    assert 'certification' not in result['counts']
    assert len(used_cuts) == result['counts']['cuts']


def test_boundary_targets_and_anchors_belong_to_current_scheme(experiment):
    _, _, monitor, result = experiment
    starts, completed, rays = [], [], 0
    for index, item in enumerate(monitor.history):
        event = item['patch'].get('event')
        state = monitor.frame(index)
        if event == 'boundary_start':
            starts.append(state['stage'])
            targets = []
            before = _union([row['inner'] for row in state['schemes'].values()]).area
        elif event == 'ray_start' and state['phase'] == '边界补充':
            ray = state['ray']
            row = state['schemes'][state['active_scheme']]
            bounds = np.asarray(state['bounds'])
            assert ray['scheme'] == state['active_scheme']
            assert contains([np.asarray(ray['anchor'])/bounds], halfspaces(np.asarray(row['inner'])/bounds))[0]
            target = np.asarray(ray['target'])
            assert np.min(np.max(np.abs(np.asarray(row['outer'])-target), axis=1)) < 1e-7
            assert all(np.max(np.abs(target-other)) > state['point_tol'] for other in targets)
            targets.append(target)
            rays += 1
        elif event == 'boundary_end':
            completed.append(state['stage'])
            after = _union([row['inner'] for row in state['schemes'].values()]).area
            assert after >= before-1e-6
        elif event == 'scheme_end':
            assert monitor.history[index+1]['patch']['event'] == 'residual_start'
    assert starts == completed == list(range(1, len(starts)+1))
    assert rays > 0
    assert result['counts']['global_search'] == len(completed)


def test_global_witness_is_feasible_and_outside_previous_coverage(experiment):
    network, _, monitor, _ = experiment
    oracle = SubProblem(GridPhysics(network, 'socp'), threads=4, numeric_focus=monitor.state['numeric_focus'])
    with threadpool_limits(limits=1):
        for index, item in enumerate(monitor.history):
            if item['patch'].get('event') != 'residual_end':
                continue
            state = monitor.frame(index)
            if state['coverage_complete']:
                continue
            previous = monitor.frame(index-1)
            witness = state['global_point']
            power = np.asarray(witness['p'])
            point = (1-state['tau'])*power/state['bounds']
            assert not any(contains([point], halfspaces(np.asarray(row['inner'])/state['bounds']))[0]
                           for row in previous['schemes'].values())
            x = np.asarray(state['schemes'][witness['scheme']]['x'])
            assert oracle.solve(x, power, time_limit=20.)['feasible']


def test_midpoint_is_skipped_when_endpoints_are_certified_in_the_same_batch():
    network = FourBus(load_nodes=(1, 2))
    monitor = RunMonitor()
    calls = 0

    def candidates(region, x):
        nonlocal calls
        calls += 1
        if calls == 1:
            return np.array([[0., 0.], [0., 10.], [0., 5.]])/region.bounds
        return stage_candidates(region, x)

    with threadpool_limits(limits=1), patch('main.stage_candidates', side_effect=candidates):
        result = build_sequential_region(network, budget=20000., monitor=monitor, seconds=60.)
    checked = [monitor.frame(i)['sp_point'] for i, row in enumerate(monitor.history)
               if row['patch'].get('event') == 'point']
    assert result['certified']
    # 原点已在首轮射线中认证，本批只需认证另一端点；中点仍由凸性跳过。
    np.testing.assert_allclose(checked[0]['p'], [0., 10.], atol=1e-12)
    assert not any(point['scheme'] == 'A' and np.allclose(point['p'], [0., 0.]) for point in checked)
    assert not any(point['scheme'] == 'A' and np.allclose(point['p'], [0., 5.]) for point in checked)


def test_convexity_is_used_only_within_one_scheme():
    region = RegionState([1., 1.], 2., .005)
    region.add_scheme([1], {}, 0.)
    region.add_scheme([0], {}, 0.)
    region.add_point([1], [[0., 1.]])
    region.add_point([0], [[1., 0.]])
    assert region.covering_schemes([[.5, .5]]) == [None]
    region.add_point([1], [[1., 0.]])
    assert region.covering_schemes([[.5, .5]]) == [(1,)]
    assert not contains([[.5, .5]], region.inner_equations([0]))[0]
    # 点合并距离不能变成认证凸包的扩边距离。
    assert not contains([[.5, .50005]], region.inner_equations([1]))[0]


def test_initial_sweep_ignores_ray_gain_threshold():
    monitor = RunMonitor()
    with threadpool_limits(limits=1), patch('main.RemainingRegionModel.solve', side_effect=TimeoutError):
        build_sequential_region(FourBus(load_nodes=(1, 2)), budget=20000., monitor=monitor,
                                seconds=10., ray_threshold=1e6)
    sweeping, small_rays = False, 0
    for index, item in enumerate(monitor.history):
        event = item['patch']['event']
        if event == 'initial_sweep_start':
            sweeping = True
        elif event == 'initial_sweep_end':
            break
        elif sweeping:
            assert event not in ('point', 'cut', 'ray_skip')
            if event == 'ray_start':
                state = monitor.frame(index)
                small_rays += state['ray_gain_bound'] < state['ray_gain_threshold']
    assert small_rays > 0


@pytest.mark.parametrize('network,focus', [(FourBus(load_nodes=(1, 2)), 0), (Case33(), 3)])
def test_deferred_cut_uses_cached_cone_directions_without_resolving_socp(network, focus):
    equations = GridPhysics(network, 'socp')
    oracle = SubProblem(equations, threads=4, numeric_focus=focus)
    x, power = network.encode_plan(network.initial_plan), np.full(2, network.power_limit)
    optimize, solves = gp.Model.optimize, []

    def record(model, *args, **kwargs):
        model.update()
        solves.append(model.NumQConstrs)
        return optimize(model, *args, **kwargs)

    with threadpool_limits(limits=1), patch.object(gp.Model, 'optimize', new=record):
        answer = oracle.solve(x, power, score_only=True)
        assert len(solves) == 1 and solves[0] > 0
        assert not answer['feasible'] and answer['cut'] is None and answer['state'] is None
        assert answer['eta'] > PLANNING_TOL and oracle.cut_calls == 0
        cut = oracle.generate_cut(x, power, answer['cone_normals'])
        assert len(solves) == 2 and solves[1] == 0
        assert oracle.calls == oracle.cut_calls == 1
        assert cut[0]+cut[1:3]@power+cut[3:]@x < -1e-9


def test_ray_geometry_bounds_the_gain_and_only_selected_points_generate_cuts(experiment):
    _, _, monitor, result = experiment
    starts, skipped, cuts = 0, 0, 0
    for index, item in enumerate(monitor.history):
        event = item['patch'].get('event')
        state = monitor.frame(index)
        if event == 'ray_start':
            starts += 1
            inner = np.asarray(state['schemes'][state['active_scheme']]['inner'])
            upper = ray_gain(inner, state['ray']['target'])
            assert state['ray_gain_bound'] == pytest.approx(upper, abs=1e-6)
            end = monitor.frame(index+1)
            assert end['event'] == 'ray_end'
            grown = np.asarray(end['schemes'][state['active_scheme']]['inner'])
            assert _union([grown]).area-_union([inner]).area <= upper+1e-6
        elif event == 'ray_skip':
            skipped += 1
            assert state['ray_gain_bound'] <= state['ray_gain_threshold']
            previous = monitor.frame(index-1)
            assert state['schemes'] == previous['schemes']
        elif event == 'cut_start':
            cuts += 1
            assert monitor.frame(index+1)['event'] == 'cut'
    assert starts == result['counts']['ray']
    assert skipped > 0
    assert cuts == result['counts']['cuts']


def test_coverage_pruning_keeps_one_equal_domain_and_preserves_nonconvex_gaps():
    region = RegionState([1., 1.], 2., .005)
    polygons = [
        [[.1, .1], [.2, .1], [.2, .2]],
        [[0., 0.], [.6, 0.], [.6, .6], [0., .6]],
        [[0., 0.], [.6, 0.], [.6, .6], [0., .6]],
        [[.4, .4], [1., .4], [1., 1.], [.4, 1.]],
        [[.1, .1], [.3, .3]],
    ]
    for x, polygon in zip(np.eye(5, dtype=int), polygons):
        region.add_scheme(x, {}, 0.)
        region.add_point(x, polygon)
    kept, equations = coverage_halfspaces(region)
    assert len(kept) == 2 and sum(len(eq) for eq in equations) == 8
    original = _union([row['inner'] for row in region.records.values()])
    retained = _union([region.records[key]['inner'] for key in kept])
    assert original.symmetric_difference(retained).area == 0.

    region = RegionState([1., 1.], 2., .005)
    polygons = [
        [[0., 0.], [.4, 0.], [.4, 1.], [0., 1.]],
        [[.6, 0.], [1., 0.], [1., 1.], [.6, 1.]],
        [[.2, .2], [.8, .2], [.8, .8]],
    ]
    for x, polygon in zip(np.eye(3, dtype=int), polygons):
        region.add_scheme(x, {}, 0.)
        region.add_point(x, polygon)
    # 三角形顶点分别被左右矩形覆盖，但其内部仍跨过缺口，不能删掉。
    assert len(coverage_halfspaces(region)[0]) == 3
