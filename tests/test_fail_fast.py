"""求解或几何失败必须停止，不重求、不换模型、不发布完成结果。"""
from types import SimpleNamespace
from unittest.mock import Mock, patch

import numpy as np
import pytest
from gurobipy import GRB
from scipy.spatial import QhullError

import geometry
import main
from model import MasterProblem
from tests.planning_checks import mainline_stubs


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


def test_mp_quality_gate_follows_the_given_tolerance():
    problem = object.__new__(MasterProblem)
    problem.equations = None
    problem.model = SimpleNamespace(Params=SimpleNamespace(), optimize=Mock(), ModelName='MP probe',
                                    Status=GRB.OPTIMAL, SolCount=1, MaxVio=5e-8)
    with pytest.raises(RuntimeError, match='MaxVio=5e-08 > 1e-08'):
        problem.solve()
    with pytest.raises(AttributeError, match='ObjBound'):   # 扫描门槛 1e-6：通过质量检查，继续读取目标界
        problem.solve(tolerance=1e-6)


def test_qhull_failure_is_not_retried_in_other_coordinates():
    points = np.array([[0.,0.,0.], [1.,0.,0.], [0.,1.,0.], [0.,0.,1.]])
    with patch('geometry.ConvexHull', side_effect=QhullError('hull probe')) as hull:
        with pytest.raises(QhullError, match='hull probe'):
            geometry.polytope_volume(points)
    hull.assert_called_once()


def test_main_records_failure_and_does_not_publish_result(tmp_path):
    with mainline_stubs() as calls:
        calls.build.side_effect = RuntimeError('solver probe')
        with pytest.raises(RuntimeError, match='solver probe'):
            main.main('fourbus', 2, show=False, output=tmp_path)
    calls.build.assert_called_once()
    calls.validate.assert_not_called()   # 构域失败时不进入校验，也不画图
    calls.draw.assert_not_called()
    output = tmp_path/'mode_1'/'fourbus_1_2.json.gz'
    from monitor import RunMonitor
    monitor = RunMonitor()
    monitor.load_recording(output)
    assert monitor.state['status'] == 'failed'
    assert 'result' not in monitor.state
    assert [path for path in tmp_path.rglob('*') if path.is_file()] == [output]
