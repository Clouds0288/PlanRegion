"""33 节点原始数据、背景负荷与独立 AC 证书校验。"""
import unittest  # 标准回归测试。
import numpy as np  # 物理数组和可复现查询。
from Network.case33bw import Case33, network  # 唯一物理数据输入。
from model import PlanningEquations, PlanningSP
from vertify import ACPowerFlow
from tests.reference import fixed_topology, upgrade_plan, validate_power_flow
from tests.planning_checks import margin


class Case33Tests(unittest.TestCase):  # 检查模型共享的数据和物理证书。
    def test_original_data_and_nodal_power_flow(self):  # 核对原始数据及两套独立 AC 实现。
        self.assertEqual(network.n, 32)  # 33 节点径向网有 32 个非根节点。
        self.assertEqual(network.n_corridors, 37)
        self.assertEqual(network.original_p.sum(), 3715.)  # 核对原始有功负荷总量，kW。
        self.assertEqual(network.original_q.sum(), 2300.)  # 核对原始无功负荷总量，kvar。
        validation = validate_power_flow(network.tree(network.encode_plan(network.initial_plan)))
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
            network = Case33(upgrade_count=count)  # 同一网架读取不同候选集。
            choice = np.arange(count)%2  # 交错选型可暴露项目顺序和支路顺序的混淆。
            plan = upgrade_plan(network, choice)
            design = network.tree(network.encode_plan(plan))
            factor = np.ones(network.n)  # 其余支路保持原阻抗。
            for selected,project in zip(choice,network.projects):  # 按 Network 项目顺序解释选型。
                if selected:  # 增设一回相同线路。
                    factor[network.nodes.index(project.branch[1])] = .5  # 并联后等值阻抗减半。
            baseline = network.tree(network.encode_plan(network.initial_plan))
            np.testing.assert_allclose(design.r, baseline.r*factor)
            np.testing.assert_allclose(design.reactance, baseline.reactance*factor)
            np.testing.assert_array_equal(design.fixed_p,network.fixed_p)  # 建设不能改变背景负荷。
            self.assertEqual(design.cost,sum(k*p.cost for k,p in zip(choice,network.projects)))  # 实际投资只由选中项目确定。

    def test_candidate_sets_are_nested_and_match_active_branches(self):  # 扩展候选不能改动原项目或误用常开联络线。
        full = Case33(upgrade_count=32)  # 完整的十六条候选配置。
        active = {c.endpoints for c in full.corridors if c.initial_active}
        self.assertEqual(len({p.branch for p in full.projects}),32)  # 候选支路不得重复。
        self.assertTrue(all(p.branch in active for p in full.projects))  # 两端及方向均与原始树一致。
        for count in (4,8,16,32):  # 小规模候选必须是同一项目表的前缀。
            case = Case33(upgrade_count=count)  # 不另维护基准项目副本。
            self.assertEqual(case.projects,full.projects[:count])  # 原候选的顺序和费用保持一致。
            self.assertEqual(sum(len(c.types) > 1 for c in case.corridors), count)
            self.assertEqual((case.n_corridors, case.n_types), (37, 37+count))

    def test_full_network_costs_and_affordable_design_count(self):
        from tests.planning_checks import affordable_designs
        case = Case33(upgrade_count=32)
        self.assertEqual(sum(p.cost for p in case.projects), 33.)
        self.assertTrue(all(p.cost == 1. for p in case.projects[16:]))
        designs = affordable_designs(fixed_topology(case), 2.)
        self.assertEqual(len(designs), 498)
        self.assertEqual(sum(c <= 1. for c, _, _ in designs), 32)
        upgraded = case.tree(case.encode_plan(upgrade_plan(case, np.ones(32))))
        baseline = case.tree(case.encode_plan(case.initial_plan))
        np.testing.assert_array_equal(upgraded.r, baseline.r/2)
        np.testing.assert_array_equal(upgraded.reactance, baseline.reactance/2)
        np.testing.assert_array_equal(upgraded.fixed_p, case.fixed_p)

    def test_socp_certificates_against_independent_ac(self):  # 固定基础网架时仍使用同一套 PlanningSP。
        x = network.encode_plan(network.initial_plan)
        reference = ACPowerFlow(network.tree(x), threads=1)
        equations = PlanningEquations(network, 'socp')
        oracle = PlanningSP(equations, threads=1)  # 所有查询共享方程。
        samples = np.random.default_rng(20260922).random((120,3))*[350.,1500.,600.]  # 固定背景下的不同负荷。
        valid = samples[reference.classify(samples)==1]  # 用独立 AC 选出已知可行样本。
        cuts = 0  # 保证实际覆盖分离割。
        for power in samples[:30]:  # 同时覆盖域内与域外。
            result = oracle.solve(x, power)
            if result['feasible']:  # 证书必须满足原始物理约束。
                self.assertGreaterEqual(margin(equations, x, power, result['state']), -1e-8)
                self.assertEqual(reference.classify(power)[0],1)  # 当前纯负荷样本还应通过完整 AC。
            if result['cut'] is not None:  # 实际有效割进一步交叉核验。
                cut = result['cut']
                alpha = cut[0]+cut[4:]@x
                self.assertLess(alpha+cut[1:4]@power, 0.)
                self.assertTrue(np.all(alpha+valid@cut[1:4] >= -1e-8))
                cuts += 1  # 累计真实割数量。
        self.assertGreater(cuts,0)  # 防止只测到了容易的可行点。

    def test_boundary_refinement_keeps_the_requested_power(self):
        case = fixed_topology(Case33(upgrade_count=8))
        x = case.encode_plan(upgrade_plan(case, np.ones(8)))
        e = PlanningEquations(case, 'socp')
        power = np.array([241.18, 5387.83, 395.69])
        original = power.copy()
        answer = PlanningSP(e, threads=1).solve(x, power)
        self.assertTrue(answer['feasible'])
        self.assertGreaterEqual(margin(e, x, power, answer['state']), -1e-8)
        np.testing.assert_array_equal(power, original)
        outside = np.array([241.45047366, 5389.44569774, 393.81480695])
        cut = PlanningSP(e, threads=1).solve(x, outside)['cut']
        self.assertIsNotNone(cut)
        self.assertLess(cut[0]+cut[1:4]@outside+cut[4:]@x, -1e-9)


if __name__=='__main__':  # 支持单独运行物理数据检查。
    unittest.main()  # 执行当前文件测试。
