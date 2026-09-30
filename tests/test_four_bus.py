"""四节点五走廊基础回归：原参数、可选拓扑、联合割及同一构域流程。"""
import unittest  # 与 case33 一起执行标准测试。
import numpy as np  # 选型、网格与数值比较。
from gurobipy import GRB  # 只接受明确的全局不可行或最优证书。
from threadpoolctl import threadpool_limits  # 与正式流程使用相同数值线程数。
from Network.four_bus_five_corridor import FourBus  # 原始基础算例。
from model import GridPhysics, MasterProblem, SubProblem
from vertify import ACPowerFlow
from tests.reference import dispatch_support  # 独立固定网架的消元方程。

PLANS = (
    {'01': 'L', '12': 'L', '13': 'L', '02': None, '23': None},
    {'01': None, '12': 'M', '13': 'H', '02': 'H', '23': None},
    {'01': 'H', '12': None, '13': 'M', '02': None, '23': 'H'},
)


class FourBusTests(unittest.TestCase):  # 不生成建设组合表，审核代表拓扑和完整紧凑模型。
    @classmethod
    def setUpClass(cls):  # 所有数值比较单线程。
        cls.threads = threadpool_limits(limits=1)  # 避免矩阵线程影响运行时间。

    @classmethod
    def tearDownClass(cls):  # 不改变其他测试的线程环境。
        cls.threads.restore_original_limits()  # 恢复原设置。

    def test_original_corridors_types_and_costs(self):  # 防止基础算例再次被删减为固定三线路树。
        c = FourBus()  # 新建网架不应预先求解或生成方案表。
        self.assertEqual([(e.endpoints, e.existing_type, e.initial_active) for e in c.corridors],
                         [((0,1),'L',True),((1,2),'L',True),((1,3),'L',True),
                          ((0,2),None,False),((2,3),None,False)])
        self.assertEqual([(t.name,t.r_ohm_km,t.x_ohm_km,t.capacity_kw,t.cable_cny_m) for t in c.lines],  # 原三种设备。
                         [('L',1.15,.08,35.,28.),('M',.62,.08,65.,43.),('H',.32,.08,100.,65.)])
        self.assertEqual(c.budgets,(20000.,40000.,60000.,np.inf))  # 默认预算保持元单位。
        self.assertEqual(c.cost_unit,'元')  # 不沿用 case33 的相对投资单位。
        self.assertFalse(hasattr(c,'designs'))  # 恢复算例不恢复旧组合枚举接口。
        np.testing.assert_array_equal([t.investment_cost for t in c.corridors[0].types],[0.,9460.,14300.])
        np.testing.assert_array_equal([t.investment_cost for t in c.corridors[3].types],[21900.,26400.,33000.])

    def test_selected_trees_match_independent_physics(self):  # 覆盖原树、12 反向和 23 反向三种代表结构。
        c = FourBus()
        for choice in PLANS:
            tree = c.tree(c.encode_plan(choice))  # 仅创建这一棵树。
            for method in ('linear','socp'):  # 两套规划物理假设分别核对。
                e = GridPhysics(c,method)  # 所有候选走廊都进入紧凑模型。
                x = e.network.encode_plan(choice)  # 只编码型号，潮流方向由功率符号决定。
                self.assertEqual(e.network.decode_plan(x),choice)
                direct = MasterProblem(e, threads=1)  # 保留同一套原始方程。
                with direct.model:  # 每次查询及时释放求解器。
                    direct.x.LB = direct.x.UB = x  # 固定指定拓扑和设备。
                    answer = direct.solve()  # 求最大总负荷。
                reference = dispatch_support(tree,method,np.ones(3))  # 独立消元计算。
                self.assertEqual(answer['status'],'optimal')  # 未确定不视为通过。
                self.assertAlmostEqual(answer['objective'],reference['value'],delta=.002)  # 边界误差不超过 0.002 kW。

    def test_topology_rejects_island_cycle_at_zero_load(self):  # 零负荷也不能让孤岛环通过径向约束。
        e = GridPhysics(FourBus(),'linear')  # 规划中必须供到全部非根节点。
        x = e.network.encode_plan({'01': None, '12': 'L', '13': 'L', '02': None, '23': 'L'})
        problem = MasterProblem(e,power=np.zeros(3), threads=1)  # 完整物理模型在零负荷下也必须排除孤岛环。
        with problem.model:  # 测试后释放模型。
            problem.x.LB = problem.x.UB = x  # 强制孤岛环作为候选。
            self.assertIsNone(problem.solve())  # 连通流必须给出明确不可行证明。

    def test_joint_queries_and_cuts_match_full_planning_model(self):  # 联合割必须适用于所有合法树和型号。
        for method in ('linear', 'socp'):
            e = GridPhysics(FourBus(), method)
            x = e.network.encode_plan(e.network.initial_plan)
            cut = SubProblem(e, threads=1).solve(x,np.array([40.,40.,40.]))['cut']  # 含型号及拓扑变量的分离割。
            self.assertIsNotNone(cut)  # 必须实际形成分离证书。
            direct = MasterProblem(e, threads=1)  # 遍历由求解器隐式搜索的全部合法拓扑。
            with direct.model:  # 测试该割在整个紧凑可行域上的最小余量。
                direct.model.setObjective(cut[0]+cut[1:4]@direct.power+cut[4:]@direct.x,GRB.MINIMIZE)  # 不加入被审核割。
                direct.model.Params.TimeLimit = 20.  # 仅限制测试的审核时间。
                direct.model.optimize()  # 全局下界检查，未靠抽样宣称有效。
                self.assertGreaterEqual(direct.model.ObjBound,-1e-7)  # 不能误切任何合法树上的可行点。

    def test_independent_ac_on_reconfigured_trees(self):  # 独立 AC 必须读取选中树，而非一直使用原始拓扑。
        c = FourBus()  # 仅共享物理配置。
        for choice in PLANS:
            oracle = ACPowerFlow(c.tree(c.encode_plan(choice)), threads=1)  # 直接把当前树交给独立模型。
            try:  # 非凸求解器按需创建并及时释放。
                for power in ([3.,4.,5.],[50.,50.,50.]):  # 确保可行与不可行均被核验。
                    self.assertEqual(int(oracle.classify(power)[0]),oracle.global_status(power,None))  # 不动点与显式 AC 等式一致。
            finally:  # 无论测试是否通过都释放资源。
                oracle.close()  # 不影响后续 case33 测试。


if __name__=='__main__':  # 支持单独运行基础网架回归。
    unittest.main()  # 不启动主 Notebook 的完整实验。
