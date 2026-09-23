"""33 节点原始数据、背景负荷与独立 AC 证书校验。"""
import unittest  # 标准回归测试。
import numpy as np  # 物理数组和可复现查询。
from Network.case33bw import Case33, network  # 唯一物理数据输入。
from model import PlanningEquations, PlanningSP
from vertify import ACPowerFlow
from tests.reference import validate_power_flow


class Case33Tests(unittest.TestCase):  # 检查模型共享的数据和物理证书。
    def test_original_data_and_nodal_power_flow(self):  # 核对原始数据及两套独立 AC 实现。
        self.assertEqual(network.n, 32)  # 33 节点径向网有 32 个非根节点。
        self.assertEqual(len(network.branches), 32)  # 固定在运树应有 32 条支路。
        self.assertEqual(network.original_p.sum(), 3715.)  # 核对原始有功负荷总量，kW。
        self.assertEqual(network.original_q.sum(), 2300.)  # 核对原始无功负荷总量，kvar。
        validation = validate_power_flow(network)  # 执行六个工况的独立潮流对照。
        self.assertLess(validation['maximum_voltage_difference_pu'], 1e-8)  # 最大节点电压幅值差须低于给定精度。
        self.assertAlmostEqual(validation['baseline']['minimum_voltage_pu'], .91309047936, places=9)

    def test_slice_background_is_fixed(self):  # 独立坐标变化不能缩放其他节点背景负荷。
        power = np.array([[100., 200., 300.], [250., 400., 50.]])  # 给出两组不同的独立节点有功。
        p, q = network.loads(power)  # 通过正式网架接口组装完整负荷。
        fixed = np.ones(network.n, dtype=bool)  # 先标记全部非根节点。
        fixed[network.selected] = False  # 排除三个独立变化节点。
        np.testing.assert_allclose(p[:, fixed]*network.base,  # 把背景有功从标幺恢复为 kW。
                                   np.tile(network.original_p[fixed], (2, 1)), atol=1e-12)  # 两次查询均须保持原始背景工况。
        np.testing.assert_allclose(q[:, fixed]*network.base,  # 同样检查背景无功。
                                   np.tile(network.original_q[fixed], (2, 1)), atol=1e-12)  # 无功背景必须等于原始 kvar 数据。
        np.testing.assert_allclose(q[:, network.selected]*network.base, power*network.q_ratio)

    def test_investment_changes_only_selected_impedances(self):  # 指定选型只改变相应支路。
        for count in (4,8,16,32):  # 覆盖三档嵌套候选集。
            network = Case33(candidate_count=count)  # 同一网架读取不同候选集。
            choice = np.arange(count)%2  # 交错选型可暴露项目顺序和支路顺序的混淆。
            plan = network.initial_plan | {c.id: c.types[k].id
                for c, k in zip(network.planning_corridors, choice)}
            design = network.design(plan)
            factor = np.ones(network.n)  # 其余支路保持原阻抗。
            for selected,project in zip(choice,network.projects):  # 按 Network 项目顺序解释选型。
                if selected:  # 增设一回相同线路。
                    factor[network.nodes.index(project.branch[1])] = .5  # 并联后等值阻抗减半。
            np.testing.assert_allclose(design.r,network.r*factor)  # 全网电阻逐项核对。
            np.testing.assert_allclose(design.reactance,network.reactance*factor)  # 全网电抗逐项核对。
            np.testing.assert_array_equal(design.fixed_p,network.fixed_p)  # 建设不能改变背景负荷。
            self.assertEqual(design.cost,sum(k*p.cost for k,p in zip(choice,network.projects)))  # 实际投资只由选中项目确定。

    def test_candidate_sets_are_nested_and_match_active_branches(self):  # 扩展候选不能改动原项目或误用常开联络线。
        full = Case33(candidate_count=32)  # 完整的十六条候选配置。
        active = {(full.root if p<0 else full.nodes[p],full.nodes[i]) for i,p in enumerate(full.parent)}  # 根向在运支路。
        self.assertEqual(len({p.branch for p in full.projects}),32)  # 候选支路不得重复。
        self.assertTrue(all(p.branch in active for p in full.projects))  # 两端及方向均与原始树一致。
        for count in (4,8,16,32):  # 小规模候选必须是同一项目表的前缀。
            case = Case33(candidate_count=count)  # 不另维护基准项目副本。
            self.assertEqual(case.projects,full.projects[:count])  # 原候选的顺序和费用保持一致。
            self.assertEqual(len(case.planning_corridors),count)
            self.assertEqual(len(case.corridors),32)

    def test_full_network_costs_and_affordable_design_count(self):
        from tests.benchmark_ac_search import affordable_designs
        case = Case33(candidate_count=32)
        self.assertEqual(sum(p.cost for p in case.projects), 33.)
        self.assertTrue(all(p.cost == 1. for p in case.projects[16:]))
        designs = affordable_designs(case, 2.)
        self.assertEqual(len(designs), 498)
        self.assertEqual(sum(c <= 1. for c, _, _ in designs), 32)
        upgraded = case.design({c.id: c.types[-1].id for c in case.corridors})
        np.testing.assert_array_equal(upgraded.r, case.r/2)
        np.testing.assert_array_equal(upgraded.reactance, case.reactance/2)
        np.testing.assert_array_equal(upgraded.fixed_p, case.fixed_p)

    def test_socp_certificates_against_independent_ac(self):  # 固定基础网架时仍使用同一套 PlanningSP。
        reference = ACPowerFlow(network, threads=1)  # 独立支路 AC 递推。
        fixed_network = network.design(network.initial_plan)
        equations = PlanningEquations(fixed_network,'socp')  # 固定网架没有离散变量。
        oracle = PlanningSP(equations, threads=1)  # 所有查询共享方程。
        samples = np.random.default_rng(20260922).random((120,3))*[350.,1500.,600.]  # 固定背景下的不同负荷。
        valid = samples[reference.classify(samples)==1]  # 用独立 AC 选出已知可行样本。
        cuts = 0  # 保证实际覆盖分离割。
        for power in samples[:30]:  # 同时覆盖域内与域外。
            result = oracle.solve(np.empty(0),power)  # 直接检查原查询，不移动负荷点。
            if result['feasible']:  # 证书必须满足原始物理约束。
                self.assertGreaterEqual(equations.margin(np.empty(0),power,result['state']),-1e-8)  # 不能只检查 eta。
                self.assertEqual(reference.classify(power)[0],1)  # 当前纯负荷样本还应通过完整 AC。
            if result['cut'] is not None:  # 实际有效割进一步交叉核验。
                cut = result['cut']  # 固定网架时只有负荷系数。
                self.assertLess(cut[0]+cut[1:]@power,0.)  # 排除产生它的查询点。
                self.assertTrue(np.all(cut[0]+valid@cut[1:]>=-1e-8))  # 不得排除独立 AC 可行点。
                cuts += 1  # 累计真实割数量。
        self.assertGreater(cuts,0)  # 防止只测到了容易的可行点。


if __name__=='__main__':  # 支持单独运行物理数据检查。
    unittest.main()  # 执行当前文件测试。
