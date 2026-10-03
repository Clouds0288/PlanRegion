"""33 节点原始数据、背景负荷与独立 AC 证书校验。"""
import unittest  # 标准回归测试。
import numpy as np  # 物理数组和可复现查询。
from Network.case33bw import network  # 唯一物理数据输入。
from tests.reference import validate_power_flow


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
        power = np.array([[100., 200.], [250., 400.]])  # 给出两组不同的独立节点有功。
        p, q = network.loads(power)  # 通过正式网架接口组装完整负荷。
        fixed = np.ones(network.n, dtype=bool)  # 先标记全部非根节点。
        fixed[network.selected] = False  # 排除两个独立变化节点。
        np.testing.assert_allclose(p[:, fixed]*network.base,  # 把背景有功从标幺恢复为 kW。
                                   np.tile(network.original_p[fixed], (2, 1)), atol=1e-12)  # 两次查询均须保持原始背景工况。
        np.testing.assert_allclose(q[:, fixed]*network.base,  # 同样检查背景无功。
                                   np.tile(network.original_q[fixed], (2, 1)), atol=1e-12)  # 无功背景必须等于原始 kvar 数据。
        np.testing.assert_allclose(q[:, network.selected]*network.base, power*network.q_ratio)

    def test_expansion_case_switches_candidates_and_budget(self):  # 扩展规划算例：开断免费、候选计建设费、原始方案费用 0。
        from types import SimpleNamespace
        from Network.case33bw import Case33Plan, CURRENT_LIMIT
        from vertify import budget_schemes
        plan = Case33Plan(current_limit=CURRENT_LIMIT)
        variable = {c.id: (c.initial_active, float(c.types[0].investment_cost)) for c in plan.corridors if c.switchable}
        self.assertEqual(variable, {'7-8': (True, 0.), '11-12': (True, 0.), '14-15': (True, 0.), '28-29': (True, 0.),
                                    '32-33': (True, 0.), '21-8': (False, 4.), '9-15': (False, 4.), '12-22': (False, 4.),
                                    '18-33': (False, 1.), '25-29': (False, 1.)})
        x0 = plan.encode_plan(plan.initial_plan)
        self.assertEqual(plan.cost_offset+plan.cost@x0, 0.)
        np.testing.assert_allclose(plan.r, network.r)  # 线路参数与原 Case33 相同（候选即原联络线）
        counts = [len(budget_schemes(SimpleNamespace(network=plan), budget)) for budget in (0, 2, plan.plan_budget)]
        self.assertEqual(counts, [1, 11, 87])

    def test_survey_case_limit_mask_and_schemes(self):  # 勘察算例：全部线路 250 A，道路掩码只屏蔽不可用的候选路。
        from types import SimpleNamespace
        from Network.case33bw import Case33S
        from vertify import budget_schemes
        full, pair, none = Case33S(), Case33S(available=(3, 4)), Case33S(available=())
        self.assertEqual(full.roads, ('21-8', '9-15', '12-22', '18-33', '25-29'))
        self.assertEqual([c.id for c, allowed in zip(pair.corridors, pair.road_allowed) if not allowed],
                         ['21-8', '9-15', '12-22'])
        np.testing.assert_allclose(full.ell_limit, (250./(full.base/(np.sqrt(3)*full.voltage_kv)))**2)
        self.assertEqual([len(budget_schemes(SimpleNamespace(network=n), 14.)) for n in (full, pair, none)], [87, 11, 1])

    def test_survey_value_and_dynamic_program(self):  # 勘察：101 个族、Φ 的阶梯式、DP 上下界与最优轨迹期望（合成 κ 表）。
        import survey
        self.assertEqual(len({survey.family(A, b) for A in survey.SUBSETS for b in survey.LEVELS}), 101)
        # κ(A,b)=0.2+0.1·#{A 中建设费不超过 b 的候选}；{C4}：b=0 时 0.2、b>=1 时 0.3，Φ=1·0.1+20·0.7
        kappa = np.array([[.2+.1*sum(survey.COST[g] <= b for g in survey.members(A)) for b in survey.LEVELS]
                          for A in survey.SUBSETS])
        phi = survey.value(kappa)
        self.assertAlmostEqual(phi[0], 16.)
        self.assertAlmostEqual(phi[1 << 3], 14.1)
        plan = survey.Survey(phi, kappa, {O: [O] for O in survey.SUBSETS})
        self.assertAlmostEqual(plan.weights.sum(), 1.)
        for A, B in survey.STATES:
            self.assertLessEqual(phi[A | survey.FULL & ~(A | B)], plan.J[A, B]+1e-12)
            self.assertLessEqual(plan.J[A, B], phi[A]+1e-12)
        self.assertAlmostEqual(plan.evaluate(plan.optimal)['expected_cost'], plan.J[0, 0])
        self.assertGreaterEqual(plan.evaluate(plan.bundle)['expected_cost'], plan.J[0, 0]-1e-12)


if __name__=='__main__':  # 支持单独运行物理数据检查。
    unittest.main()  # 执行当前文件测试。
