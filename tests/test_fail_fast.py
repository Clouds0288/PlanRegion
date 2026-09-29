"""求解或几何失败必须停止，不重求、不换模型、不发布完成结果。"""
from contextlib import redirect_stdout
from io import StringIO
from types import SimpleNamespace
from unittest.mock import Mock, patch
import json

import numpy as np
import pytest
from gurobipy import GRB
from scipy.spatial import QhullError

import main
import region
import survey
from model import GridPhysics, MasterProblem, RemainingRegionModel
from vertify import ACPowerFlow, ac_planning_query


@pytest.mark.parametrize('status,violation,message', [
    (GRB.SUBOPTIMAL, 0., 'status=13'),
    (GRB.NUMERIC, 0., 'status=12'),
    (GRB.OPTIMAL, 1e-4, 'MaxVio'),
])
def test_mp_rejects_bad_solver_result_after_one_solve(status, violation, message):
    problem = object.__new__(MasterProblem)
    problem.equations = None
    problem.model = SimpleNamespace(Params=SimpleNamespace(), optimize=Mock(), ModelName='MP probe',
                                    Status=status, SolCount=1, MaxVio=violation)
    with pytest.raises(RuntimeError, match=message):
        problem.solve()
    problem.model.optimize.assert_called_once()


def test_residual_timeout_is_not_a_coverage_certificate():
    problem = RemainingRegionModel(GridPhysics(main.FourBus(), 'linear'), 0., np.array([150.]*3),
                                   450., [], [], 0., threads=1)
    with problem.model, pytest.raises(TimeoutError, match='time limit'):
        problem.solve(region.GEOMETRY_TOL, time_limit=0.)


def test_road_domain_failure_does_not_switch_to_physical():
    cache = survey.DomainCache()
    result = dict(status='unknown', inner=[], outer=[])
    with patch('survey.build_continuous_region', return_value=result) as build:
        with pytest.raises(RuntimeError, match='Roads'):
            cache[frozenset()]
    build.assert_called_once()
    assert build.call_args.kwargs['residual_mode'] == 'light'
    assert not cache and not cache.regions and not cache.cuts


def test_qhull_failure_is_not_retried_in_other_coordinates():
    points = np.array([[0.,0.,0.], [1.,0.,0.], [0.,1.,0.], [0.,0.,1.]])
    with patch('region.ConvexHull', side_effect=QhullError('hull probe')) as hull:
        with pytest.raises(QhullError, match='hull probe'):
            region.polytope_volume(points)
    hull.assert_called_once()


def test_singular_support_is_not_silently_skipped():
    state = region.RegionState([1.,1.], 2., 0.)
    state.add_scheme([0], {}, 0.)
    with patch('region.np.linalg.solve', side_effect=np.linalg.LinAlgError('support probe')) as solve:
        with pytest.raises(np.linalg.LinAlgError, match='support probe'):
            state.witness_support([0], [.25,.25])
    solve.assert_called_once()


def test_missing_witness_support_raises():
    state = region.RegionState([1.,1.], 2., 0.)
    state.add_scheme([0], {}, 0.)
    with pytest.raises(RuntimeError, match='No simplex supports'):
        state.witness_support([0], [2.,2.])


def test_ac_nonconvergence_does_not_start_global_solver():
    with patch('vertify.AC_ITERATIONS', 0), patch.object(ACPowerFlow, 'global_status') as global_solver:
        with pytest.raises(RuntimeError, match='did not converge'):
            ac_planning_query(GridPhysics(main.FourBus(), 'linear'), np.zeros(3), threads=1)
    global_solver.assert_not_called()


def test_main_records_failure_and_does_not_publish_result(tmp_path):
    from monitor import RunMonitor
    output = tmp_path/'monitor.json.gz'
    with patch('main.build_sequential_region', side_effect=RuntimeError('solver probe')) as solve:
        with pytest.raises(RuntimeError, match='solver probe'):
            main.run(main.FourBus(load_nodes=(1, 2)), output=output, show_ui=False, threads=1)
    solve.assert_called_once()
    monitor = RunMonitor()
    monitor.load_recording(output)
    assert monitor.state['status'] == 'failed'
    assert 'solver probe' in monitor.state['error']
    assert 'result' not in monitor.state
    assert [p.name for p in tmp_path.iterdir()] == ['monitor.json.gz']
