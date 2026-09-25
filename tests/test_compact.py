"""紧凑模型、联合割和固定选型的独立消元核对；不生成建设组合表。"""
import unittest  # 直接运行模型回归。
import numpy as np  # 构造固定查询和型号向量。
from threadpoolctl import threadpool_limits  # 所有比较保持单线程。
from Network.case33bw import Case33  # 四、八、十六条候选使用同一网架定义。
from Network.four_bus_five_corridor import FourBus
from model import PlanningEquations, PlanningModel, PlanningSP  # 正式方程、主问题和连续 SP。
from tests.reference import dispatch_support, fixed_topology, upgrade_plan
from tests.planning_checks import margin
from tests.planning_checks import joint_benders  # 冻结的旧算法，与完整模型对照。


class CompactPlanningTests(unittest.TestCase):  # 核对物理映射、最优值和全域割有效性。
    @classmethod
    def setUpClass(cls):  # 固定数值计算线程数。
        cls.threads = threadpool_limits(limits=1)  # 减少求解路径差异。

    @classmethod
    def tearDownClass(cls):  # 还原调用环境。
        cls.threads.restore_original_limits()  # 测试不永久更改线程配置。

    def test_sp_near_boundary_returns_a_physical_certificate_or_valid_cut(self):
        cases = ((16, {}, [69.1086468340237, 3000.43529371335, 93.6282845836038]),
                 (4, {'2-3': 'parallel'}, [154.2787946320288, 3500.597509354904, 240.65860001331868]))
        for count, upgrades, point in cases:
            network = fixed_topology(Case33(upgrade_count=count))
            equations = PlanningEquations(network, 'socp')
            x, power = network.encode_plan(network.initial_plan | upgrades), np.array(point)
            checked = PlanningSP(equations, threads=1).solve(x, power)
            self.assertTrue(checked['feasible'] or checked['cut'] is not None)
            if checked['feasible']:
                self.assertGreaterEqual(margin(equations, x, power, checked['state']), -1e-8)
            else:
                cut = checked['cut']
                self.assertLess(cut[0]+cut[1:4]@power+cut[4:]@x, -1e-9)
                problem = PlanningModel(equations, threads=1)
                with problem.model:
                    problem.model.setObjective(cut[0]+cut[1:4]@problem.power+cut[4:]@problem.x, 1)
                    problem.model.optimize()
                    self.assertGreaterEqual(problem.model.ObjBound, -1e-7)
            np.testing.assert_array_equal(power, point)

    def test_sp_handles_zero_flow_cones_of_unselected_types(self):
        network = FourBus()
        equations = PlanningEquations(network, 'socp')
        cases = (({'01': 'H', '12': None, '13': None, '02': 'L', '23': 'M'},
                  [47.47692657884853, 14.400743716323246, 7.840840158490067], True),
                 ({'01': 'H', '12': None, '13': 'L', '02': 'H', '23': None},
                  [82.83826215331091, 41.51433150869895, 11.864741576952936], False))
        for plan, point, feasible in cases:
            x, power = network.encode_plan(plan), np.array(point)
            checked = PlanningSP(equations, threads=1).solve(x, power)
            self.assertEqual(checked['feasible'], feasible)
            if feasible:
                self.assertGreaterEqual(margin(equations, x, power, checked['state']), -1e-8)
            else:
                cut = checked['cut']
                self.assertIsNotNone(cut)
                self.assertLess(cut[0]+cut[1:4]@power+cut[4:]@x, -1e-9)

    def test_fixed_assignments_match_independent_physics(self):  # 指定代表选型，不展开全部建设组合。
        for count in (4,8,16):  # 规模变化只改变 Network 的候选集。
            network = Case33(upgrade_count=count)  # 唯一物理输入。
            for choice in (np.zeros(count,dtype=int),np.arange(count)%2,np.ones(count,dtype=int)):  # 基础、交错、全升级三个代表选型。
                choice = upgrade_plan(network, choice)
                design = network.tree(network.encode_plan(choice))  # 独立固定网架直接由型号表生成。
                for method in ('linear','socp'):  # LP 和 SOCP 分别对照。
                    equations = PlanningEquations(network,method)  # 按模型内部顺序展开型号。
                    x = equations.network.encode_plan(choice)  # Network 顺序转内部变量顺序。
                    self.assertEqual(equations.network.decode_plan(x),choice)
                    problem = PlanningModel(equations, threads=1)  # 紧凑模型保留同一物理方程。
                    with problem.model:  # 求解结束后释放优化器。
                        problem.x.LB = problem.x.UB = x  # 固定本次指定型号。
                        actual = problem.solve()  # 求固定方案的最大总负荷。
                        self.assertEqual(problem.x.shape, (37+count,))
                    reference = dispatch_support(design,method,np.ones(3))  # 不读取正式方程的独立消元求解。
                    self.assertLess(abs(actual['objective']-reference['value']),.002)  # 物理边界必须在规定 kW 容差内一致。
                    self.assertGreaterEqual(margin(equations,x,actual['p'],actual['state']),-1e-7)  # 同时复核原始约束。

    def test_joint_queries_match_direct_models_and_reuse_cuts(self):  # 审核从 Notebook 移到测试，不进入正式构域。
        for count in (4,8,16):  # 三种规模共用同一算法。
            for method in ('linear','socp'):  # 两套物理假设分别比较。
                equations,cuts = PlanningEquations(fixed_topology(Case33(upgrade_count=count)),method),[]
                queries = [dict(power=p) for p in ([100.,800.,150.],[250.,1800.,350.],[300.,3000.,500.])]  # 含基础、升级和高负荷点。
                budgets = (0.,1.,2.) if count==16 and method=='socp' else (0.,1.,2.,np.inf)
                queries += [dict(budget=b) for b in budgets]  # 覆盖各预算下的自由最大总负荷查询。
                for query in queries:  # 每次使用完全相同的输入。
                    actual,new = joint_benders(equations,cuts=cuts,**query)  # 正式 MP/SP 查询。
                    cuts.extend(new)  # 后续查询复用同一批有效割。
                    direct = PlanningModel(equations,**query)  # 完整 MILP/MISOCP 独立求解路径。
                    with direct.model:  # 不把直接模型加入联合割池。
                        reference = direct.solve()  # 取得可行目标及有效全局界。
                    self.assertEqual(actual is None,reference is None)  # 可行性结论必须相同。
                    if actual is not None:  # 有解时还须核对目标值与证书。
                        self.assertEqual(actual['status'],'optimal')  # 未确定不算通过。
                        self.assertEqual(reference['status'],'optimal')  # 参考也必须完成求解。
                        self.assertLess(abs(actual['objective']-reference['objective']),.002 if 'budget' in query else 1e-7)  # 使用各自物理单位的精度。

    def test_joint_cut_is_valid_on_the_full_compact_domain(self):  # 对全部整数选型的可行域直接优化，检验联合割没有误切。
        for count in (4,8,16):  # 同一检查覆盖全部三档候选集。
            for method in ('linear','socp'):  # 两类对偶割都需通过。
                equations = PlanningEquations(Case33(upgrade_count=count),method)  # 生成被审核方程。
                x,power = equations.network.encode_plan(equations.network.initial_plan),np.array([300.,3000.,500.])
                cut = PlanningSP(equations, threads=1).solve(x,power)['cut']  # 获得真实 SP 分离割。
                self.assertIsNotNone(cut)  # 测试必须实际产生割。
                self.assertLess(cut[0]+cut[1:4]@power+cut[4:]@x,-1e-5)  # 必须排除生成点。
                problem = PlanningModel(equations, threads=1)  # 全部原始物理约束，不能加入被审核的割。
                with problem.model:  # 全局审核与正式求解隔离。
                    problem.model.setObjective(cut[0]+cut[1:4]@problem.power+cut[4:]@problem.x,1)  # 最小化割余量。
                    problem.model.Params.TimeLimit = 20.  # 测试单次求解上限。
                    problem.model.optimize()  # MILP/MISOCP 搜索全部合法选型。
                    self.assertGreaterEqual(problem.model.ObjBound,-1e-7)  # 全局下界非负才证明无误切。


if __name__=='__main__':  # 支持单独执行本组回归。
    unittest.main()  # 不启动 Notebook 的完整实验。
