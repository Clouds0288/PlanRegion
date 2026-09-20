"""逐线路紧凑模型与联合割；参考来自旧固定方案消元模型。"""
import unittest  # 使用标准测试框架验证紧凑规划模型。
from unittest.mock import patch  # 阻止正式求解器偷偷访问完整方案表。
import numpy as np  # 处理选型编码、样本和边界误差。
from threadpoolctl import threadpool_limits  # 统一矩阵运算线程数，使数值回归可复现。

from Network.case33bw import Case33, network  # 读取候选线路设置及网架类型。
from model import PlanningEquations, PlanningModel, PlanningSP  # 使用三层规划接口，不导入旧固定方案求解类。
from tests.reference import dispatch_support, dispatch_feasible  # 独立消元参考不读取正式规划矩阵。
from notebook_flow import workflow  # 直接测试 Notebook 的联合割循环。
from vertify import verify_joint_cuts  # 在完整连续方案域上审核联合割。


class CompactPlanningTests(unittest.TestCase):  # 检验逐线路整数模型及全局联合割。
    @classmethod  # 整组测试共用一次线程配置。
    def setUpClass(cls):  # 开始测试前统一数值计算环境。
        cls.threads = threadpool_limits(limits=1)  # 固定为单线程以减少浮点路径变化。

    @classmethod  # 结束整组测试后恢复调用环境。
    def tearDownClass(cls):  # 避免测试线程配置影响其他实验。
        cls.threads.restore_original_limits()  # 还原原来的线性代数线程数。

    def test_solvers_never_access_complete_design_table(self):  # 正式求解必须不依赖预先枚举完整建设组合。
        def forbidden(_):  # 任何完整方案表读取都立即使测试失败。
            raise AssertionError('The compact solver attempted to enumerate designs')  # 清楚报告意外枚举的位置。
        with patch.object(Case33,'designs',property(forbidden)):  # 临时把方案表替换成禁止访问的属性。
            case = Case33()  # 使用全新网架避免已有缓存掩盖访问。
            for method in ('linear','socp'):  # 分别检查 MILP 和 MISOCP。
                problem = PlanningModel(PlanningEquations(case,method),power=[100.,800.,150.])  # 固定负荷，求最低逐线路建设投资。
                with problem.model:  # 求解完成后释放直接模型。
                    answer = problem.solve()  # 执行正式紧凑求解。
                    self.assertEqual(problem.model.NumBinVars,8)  # 四条候选线路两个型号恰有 8 个二进制变量。
                self.assertEqual(answer['objective'],0.)  # 该低负荷点无需升级。
                answer,cuts,_ = workflow()['joint_benders'](PlanningEquations(case,method),budget=2.,direction=[1.,6.,1.])  # 联合割边界流程也必须遵守禁止枚举条件。
                self.assertIsNotNone(answer)  # 应找到预算内可行边界。
                self.assertTrue(cuts)  # 保证确实发生了联合割迭代。

    def test_each_fixed_assignment_matches_old_eliminated_model(self):  # 每个整数选型须等价于独立固定网架模型。
        direction = np.array([1.,6.,1.])  # 选择穿过多个瓶颈方向的代表射线。
        for method in ('linear','socp'):  # 线性和 SOCP 分别对照。
            for design in network.designs:  # 16 个完整方案仅在测试中逐一固定。
                problem = PlanningModel(PlanningEquations(network,method),direction=direction)  # 直接模型仍使用逐线路变量表示。
                x = problem.equations.selection(design.x)  # 将参考方案编码为紧凑模型的 one-hot 向量。
                with problem.model:  # 该次直接模型求完即释放。
                    problem.x.LB = problem.x.UB = x  # 将每个建设变量上下界锁定到参考选型。
                    actual = problem.solve()  # 求当前固定组合的方向边界。
                reference = dispatch_support(design,method,np.ones(3),direction=direction)  # 用独立消元方程求同一方向参考边界。
                self.assertLess(abs(actual['objective']-reference['value']),.002)  # 绝对半径差须小于 0.002 kW。
                np.testing.assert_array_equal(problem.equations.choice(actual['x']),design.x)  # 返回的逐线路选型须与固定组合一致。
                self.assertGreaterEqual(problem.equations.margin(x,actual['p'],actual['state']),-1e-7)  # 直接模型返回状态还须通过原始物理约束核验。

    def test_minimum_cost_matches_all_sixteen_designs(self):  # 跨全部 16 个方案核对固定负荷最小投资。
        points = np.array([[100.,800.,150.],[250.,1800.,350.],[300.,800.,100.],[400.,4500.,700.],  # 包含基础可行、需升级及不可行负荷。
                           [204.62908536,306.94362804,102.31454268],[220.45070568]*3])  # 补充接近模型边界的典型查询。
        for method in ('linear','socp'):  # 分别对照线性和 SOCP 的相同物理假设。
            reference = np.full(len(points),np.inf)  # 参考最低费用初始化为无可行方案。
            for design in network.designs:  # 逐个参考方案检查所有点。
                feasible = dispatch_feasible(design,method,points)  # 独立消元模型判定当前网架是否承载负荷。
                reference[feasible] = np.minimum(reference[feasible],design.cost)  # 在可行方案中取最小增量投资。
            for i,power in enumerate(points):  # 正式紧凑模型逐点求投资，无需枚举。
                problem = PlanningModel(PlanningEquations(network,method),power=power)  # 固定当前负荷，并保留全部逐线路选型自由度。
                with problem.model:  # 按查询释放优化模型。
                    answer = problem.solve()  # 由混合整数求解器搜索建设组合。
                self.assertEqual(np.inf if answer is None else answer['objective'],reference[i])  # 紧凑结果或不可行结论须与完整参考一致。

    def test_joint_cut_is_valid_over_every_continuous_design_region(self):  # 一条联合割须对其他建设组合也保持有效。
        for method in ('linear','socp'):  # 同时检查线性和锥模型的对偶割。
            e = PlanningEquations(network,method)  # 审核 SP 的割只需要方程，无须建立无关主问题。
            x,power = e.selection([0,0,0,0]),np.array([300.,3000.,500.])  # 选择基础网架无法承载的高负荷点。
            cut = PlanningSP(e).solve(x,power)['cut']  # 固定选型后调用连续 SP 生成联合割。
            self.assertLess(cut[0]+cut[1:4]@power+cut[4:]@x,-1e-5)  # 该割必须严格排除原来的 (p,x)。
            self.assertGreater(np.max(np.abs(cut[4:])),1e-4)  # 确实含有建设变量系数。
            self.assertGreaterEqual(verify_joint_cuts(network,method,[cut]).min(),-1e-7)  # 在每个连续方案域上优化最小余量，不能只检查随机点。

    def test_joint_queries_match_direct_models_and_reuse_cuts(self):  # 跨点及跨预算复用联合割后仍应得到正确规划查询结果。
        for method in ('linear','socp'):  # 两种物理模型分别维护自己的割池。
            cuts = []  # 从空割池开始，后续查询共享新增割。
            for query in (dict(power=[250.,1800.,350.]),dict(budget=2.,direction=[1.,6.,1.]),  # 覆盖固定点投资查询及有预算的方向查询。
                          dict(budget=0.,direction=[0.,1.,0.])):  # 再检查零预算、轴向负荷的边界情形。
                actual,new,_ = workflow()['joint_benders'](PlanningEquations(network,method),cuts=cuts,**query)  # 执行 Notebook 中正式的 MP/SP 循环。
                cuts.extend(new)  # 仅把本次新增割并入下一次查询。
                problem = PlanningModel(PlanningEquations(network,method),**query)  # 以含全部物理方程的直接模型作独立求解路径对照。
                with problem.model:  # 该次直接模型用完即释放。
                    reference = problem.solve()  # 得到直接 MILP/MISOCP 的最优参考结果。
                self.assertEqual(actual is None,reference is None)  # 两条求解路径须对可行性给出相同结论。
                if actual is not None:  # 两者有解时进一步比较目标值。
                    self.assertLess(abs(actual['objective']-reference['objective']),.002)  # 投资或总负荷差须满足同一绝对容差。


if __name__ == '__main__':  # 允许单独执行紧凑模型测试。
    unittest.main()  # 运行该文件的全部回归检查。
