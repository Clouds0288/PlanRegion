"""求解或几何失败必须停止，不重求、不换模型、不发布完成结果。"""
from types import SimpleNamespace
from unittest.mock import Mock, patch

import numpy as np
import pytest
from gurobipy import GRB
from scipy.spatial import QhullError

import main
import region
from model import GridPhysics, MasterProblem, RemainingRegionModel
from vertify import ACPowerFlow


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
    problem = RemainingRegionModel(GridPhysics(main.FourBus(), 'socp'), 0., np.array([150.]*3),
                                   450., [], [], 0., threads=1)
    with problem.model, pytest.raises(TimeoutError, match='time limit'):
        problem.solve(region.GEOMETRY_TOL, time_limit=0.)


def test_qhull_failure_is_not_retried_in_other_coordinates():
    points = np.array([[0.,0.,0.], [1.,0.,0.], [0.,1.,0.], [0.,0.,1.]])
    with patch('region.ConvexHull', side_effect=QhullError('hull probe')) as hull:
        with pytest.raises(QhullError, match='hull probe'):
            region.polytope_volume(points)
    hull.assert_called_once()


def test_ac_nonconvergence_does_not_start_global_solver():
    with patch('vertify.AC_ITERATIONS', 0), patch.object(ACPowerFlow, 'global_status') as global_solver:
        with pytest.raises(RuntimeError, match='did not converge'):
            ACPowerFlow(main.FourBus().tree(main.FourBus().encode_plan(main.FourBus().initial_plan))).classify(np.zeros(3))
    global_solver.assert_not_called()


def test_main_records_failure_and_does_not_publish_result(tmp_path):
    output = tmp_path/'monitor.json.gz'
    with patch('region.build_sequential_region', side_effect=RuntimeError('solver probe')) as solve:
        with pytest.raises(RuntimeError, match='solver probe'):
            main.run(main.FourBus(load_nodes=(1, 2)), output=output, show_ui=False, threads=1)
    solve.assert_called_once()
    from monitor import RunMonitor
    monitor = RunMonitor()
    monitor.load_recording(output)
    assert monitor.state['status'] == 'failed'
    assert 'result' not in monitor.state
    assert list(tmp_path.iterdir()) == [output]
