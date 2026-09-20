"""逐线路紧凑模型与联合割；参考来自旧固定方案消元模型。"""
import unittest
from unittest.mock import patch
import numpy as np
from threadpoolctl import threadpool_limits

from Network.case33bw import Case33, network
from model import PlanningModel, PlanningSP, LinearSP, SOCPSP, dispatch_support
from notebook_flow import workflow
from vertify import verify_joint_cuts


class CompactPlanningTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.threads = threadpool_limits(limits=1)

    @classmethod
    def tearDownClass(cls):
        cls.threads.restore_original_limits()

    def test_solvers_never_access_complete_design_table(self):
        def forbidden(_):
            raise AssertionError('The compact solver attempted to enumerate designs')
        with patch.object(Case33,'designs',property(forbidden)):
            case = Case33()
            for method in ('linear','socp'):
                problem = PlanningModel(case,method,power=[100.,800.,150.])
                with problem.model:
                    answer = problem.solve()
                    self.assertEqual(problem.model.NumBinVars,8)
                self.assertEqual(answer['objective'],0.)
                answer,cuts,_ = workflow()['joint_benders'](case,method,budget=2.,direction=[1.,6.,1.])
                self.assertIsNotNone(answer)
                self.assertTrue(cuts)

    def test_each_fixed_assignment_matches_old_eliminated_model(self):
        direction = np.array([1.,6.,1.])
        for method in ('linear','socp'):
            for design in network.designs:
                problem = PlanningModel(network,method,direction=direction)
                x = problem.equations.selection(design.x)
                with problem.model:
                    problem.x.LB = problem.x.UB = x
                    actual = problem.solve()
                reference = dispatch_support(design,method,np.ones(3),direction=direction)
                self.assertLess(abs(actual['objective']-reference['value']),.002)
                np.testing.assert_array_equal(problem.equations.choice(actual['x']),design.x)
                self.assertGreaterEqual(problem.equations.margin(x,actual['p'],actual['state']),-1e-7)

    def test_minimum_cost_matches_all_sixteen_designs(self):
        points = np.array([[100.,800.,150.],[250.,1800.,350.],[300.,800.,100.],[400.,4500.,700.],
                           [204.62908536,306.94362804,102.31454268],[220.45070568]*3])
        for method in ('linear','socp'):
            reference = np.full(len(points),np.inf)
            for design in network.designs:
                oracle = LinearSP(design) if method=='linear' else SOCPSP(design)
                for i,power in enumerate(points):
                    if oracle.solve(power)[0]<=1e-8:
                        reference[i] = min(reference[i],design.cost)
                oracle.close()
            for i,power in enumerate(points):
                problem = PlanningModel(network,method,power=power)
                with problem.model:
                    answer = problem.solve()
                self.assertEqual(np.inf if answer is None else answer['objective'],reference[i])

    def test_joint_cut_is_valid_over_every_continuous_design_region(self):
        for method in ('linear','socp'):
            problem = PlanningModel(network,method,power=[300.,3000.,500.],cuts_only=True)
            with problem.model:
                e = problem.equations
                x,power = e.selection([0,0,0,0]),np.array([300.,3000.,500.])
                cut = PlanningSP(e).solve(x,power)['cut']
            self.assertLess(cut[0]+cut[1:4]@power+cut[4:]@x,-1e-5)
            self.assertGreater(np.max(np.abs(cut[4:])),1e-4)  # 确实含有建设变量系数。
            self.assertGreaterEqual(verify_joint_cuts(network,method,[cut]).min(),-1e-7)

    def test_joint_queries_match_direct_models_and_reuse_cuts(self):
        for method in ('linear','socp'):
            cuts = []
            for query in (dict(power=[250.,1800.,350.]),dict(budget=2.,direction=[1.,6.,1.]),
                          dict(budget=0.,direction=[0.,1.,0.])):
                actual,new,_ = workflow()['joint_benders'](network,method,cuts=cuts,**query)
                cuts.extend(new)
                problem = PlanningModel(network,method,**query)
                with problem.model:
                    reference = problem.solve()
                self.assertEqual(actual is None,reference is None)
                if actual is not None:
                    self.assertLess(abs(actual['objective']-reference['objective']),.002)


if __name__ == '__main__':
    unittest.main()
