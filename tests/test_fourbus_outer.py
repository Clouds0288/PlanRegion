"""独立外域实验的数学回归：完整锥对偶、跨网架割及全局界。"""
import contextlib
import io
import tempfile
import unittest

import numpy as np
from gurobipy import GRB

from experiments.fourbus_outer import GlobalViolation, run_experiment
from Network.four_bus_five_corridor import FourBus
from model import GridPhysics, MasterProblem, PLANNING_TOL, new_model


def primal_violation(equations, x, power):
    """独立求原始 SOCP 最小 eta，不使用 G 的系数导出或对偶模型。"""
    with new_model('primal_SP_audit', 1) as model:
        model.Params.Aggregate = 0
        model.Params.ScaleFlag = 0
        model.Params.BarQCPConvTol = 1e-9
        eta = model.addVar(name='violation')
        equations.add_operation(model, dict(zip(equations.keys, x)),
                                dict(zip(equations.network.load_nodes, power)), eta)
        model.setObjective(eta)
        model.optimize()
        if model.Status != GRB.OPTIMAL or model.MaxVio > PLANNING_TOL:
            raise RuntimeError(f'primal audit status={model.Status}, MaxVio={model.MaxVio}')
        return eta.X


class FourBusOuterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.equations = GridPhysics(FourBus(), 'socp')
        with tempfile.TemporaryDirectory() as output, contextlib.redirect_stdout(io.StringIO()):
            cls.result = run_experiment(iterations=6, time_limit=5., threads=4, output=output)

    def test_full_conic_dual_matches_primal_sp(self):
        e, result = self.equations, self.result
        points = [(np.array(item['x']), np.array(item['p'])) for item in result['trace']]
        points += [(np.array(item['x']), np.array(item['p'])) for item in result['initial'][:2]]
        points.append((e.network.encode_plan(e.network.initial_plan), np.zeros(3)))
        for x, power in points:
            search = GlobalViolation(e, 20000., result['axis_bounds'], result['total_bound'], threads=1)
            with search.model, self.subTest(x=x.tolist(), p=power.tolist()):
                template = search.template
                search.model.update()
                search.model.remove([row for row in search.model.getQConstrs()
                                     if row.QCName == 'dual_objective'])
                # 将固定数值代入对偶目标，直接求凸 SOCP，消去本已固定的乘积。
                objective = ((template.sp_x@x+template.sp_p@power-template.sp_rhs)@search.row_dual
                             -template.cone_constant@search.cone_dual)
                search.model.setObjective(objective, GRB.MAXIMIZE)
                search.problem.x.LB = search.problem.x.UB = x
                search.problem.power.LB = search.problem.power.UB = power
                search.model.Params.NonConvex = 0
                search.model.Params.BarQCPConvTol = 1e-9
                search.model.optimize()
                self.assertEqual(search.model.Status, GRB.OPTIMAL)
                self.assertLessEqual(search.model.MaxVio, PLANNING_TOL)
                eta = primal_violation(e, x, power)
                self.assertAlmostEqual(search.model.ObjVal, eta, delta=2e-7)
                self.assertGreaterEqual(search.model.ObjBound, eta-2e-7)

    def test_cuts_are_valid_over_all_topologies_and_budgets(self):
        e = self.equations
        for coefficients, record in zip(self.result['cuts'], self.result['trace']):
            cut = np.array(coefficients)
            self.assertLess(record['cut_at_candidate'], -1e-9)
            # 不加预算，也不加待审核割；直接在完整物理可行域最小化割的余量。
            direct = MasterProblem(e, threads=1)
            with direct.model:
                direct.model.setObjective(cut[0]+cut[1:4]@direct.power+cut[4:]@direct.x, GRB.MINIMIZE)
                direct.model.Params.TimeLimit = 30.
                direct.model.optimize()
                self.assertEqual(direct.model.Status, GRB.OPTIMAL)
                self.assertGreaterEqual(direct.model.ObjBound, -1e-7)

    def test_global_search_bounds_and_incomplete_status(self):
        trace = self.result['trace']
        self.assertFalse(self.result['certified'])
        self.assertEqual(self.result['status'], 'iteration_limit')
        self.assertEqual(len(self.result['cuts']), 6)
        for record in trace:
            # 时限内可能尚未证明最大值；候选原始最小违反量必须位于 G 的上下界内。
            eta = primal_violation(self.equations, record['x'], record['p'])
            self.assertLessEqual(record['objective'], eta+2e-7)
            self.assertGreaterEqual(record['bound'], eta-2e-7)


if __name__ == '__main__':
    unittest.main()
