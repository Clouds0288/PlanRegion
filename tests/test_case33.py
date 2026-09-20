"""同一个模型在 33 节点上的数据、固定背景、独立 AC 与区域停止校验。"""
import unittest  # 使用标准测试运行器执行断言。
import numpy as np  # 生成可复现样本并核对物理数组。
from Network.case33bw import network  # 读取唯一的 33 节点网架设置。
from model import PlanningEquations, PlanningSP, ACPowerFlow  # 使用统一运行方程、连续 SP 和独立 AC。
from region import halfspaces  # 检查完整几何停止条件。
from vertify import validate_power_flow  # 调用节点导纳矩阵交叉核验。
from notebook_flow import workflow  # 直接加载 Notebook 的正式流程。


class Case33Tests(unittest.TestCase):  # 覆盖网架、固定背景、切割证书和规划并集。
    def test_original_data_and_nodal_power_flow(self):  # 核对原始数据及两套独立 AC 实现。
        self.assertEqual(network.n, 32)  # 33 节点径向网有 32 个非根节点。
        self.assertEqual(len(network.branches), 32)  # 固定在运树应有 32 条支路。
        self.assertEqual(network.original_p.sum(), 3715.)  # 核对原始有功负荷总量，kW。
        self.assertEqual(network.original_q.sum(), 2300.)  # 核对原始无功负荷总量，kvar。
        validation = validate_power_flow(network)  # 执行六个工况的独立潮流对照。
        self.assertLess(validation['maximum_voltage_difference_pu'], 1e-8)  # 最大节点电压幅值差须低于给定精度。
        self.assertAlmostEqual(validation['baseline']['minimum_voltage_pu'], .91309047936, places=9)  # 基础工况最低电压应复现已知值。

    def test_slice_background_is_fixed(self):  # 独立坐标变化不能缩放其他节点背景负荷。
        power = np.array([[100., 200., 300.], [250., 400., 50.]])  # 给出两组不同的独立节点有功。
        p, q = network.loads(power)  # 通过正式网架接口组装完整负荷。
        fixed = np.ones(network.n, dtype=bool)  # 先标记全部非根节点。
        fixed[network.selected] = False  # 排除三个独立变化节点。
        np.testing.assert_allclose(p[:, fixed]*network.base,  # 把背景有功从标幺恢复为 kW。
                                   np.tile(network.original_p[fixed], (2, 1)), atol=1e-12)  # 两次查询均须保持原始背景工况。
        np.testing.assert_allclose(q[:, fixed]*network.base,  # 同样检查背景无功。
                                   np.tile(network.original_q[fixed], (2, 1)), atol=1e-12)  # 无功背景必须等于原始 kvar 数据。
        np.testing.assert_allclose(q[:, network.selected]*network.base, power*network.q_ratio)  # 独立节点无功仅按给定 Q/P 比例变化。

    def test_socp_state_witness_and_cuts_against_ac(self):  # 同时核验 SP 状态、内点与分离割。
        reference = ACPowerFlow(network)  # 独立 AC 不复用 SP 系数矩阵。
        equations = PlanningEquations(network,'socp',planning=False)  # 使用固定网架模式，不引入投资变量。
        oracle = PlanningSP(equations)  # 同一连续 SP 处理所有负荷样本。
        origin = oracle.solve(np.empty(0),np.zeros(3))['state']  # 取得固定背景下的原点运行证书。
        samples = np.random.default_rng(20260922).random((120, 3))*[350., 1500., 600.]  # 固定种子生成边界内外的查询点。
        valid = samples[reference.classify(samples) == 1]  # 独立选取 AC 可行点用于检查割没有误删。
        cuts = 0  # 统计实际产生的分离割。
        for power in samples[:30]:  # 对部分样本逐点检验完整证书。
            result = workflow()['fixed_query'](oracle,power,origin)  # 调用 Notebook 中正式的内点构造流程。
            eta,cut,witness,state = (result[k] for k in ('eta','cut','witness','state'))  # 区分原查询的违反量、割及生成的内点。
            ell = state[equations.ell]  # 取得内点运行证书中的电流平方。
            actual_state = (state[equations.P],state[equations.Q],state[equations.v])  # 提取 SP 给出的有功、无功和电压平方。
            for actual,expected in zip(actual_state,reference.state(witness,ell)[:3]):  # 用同一电流在独立支路递推中重建状态。
                np.testing.assert_allclose(actual, expected[0], atol=1e-12, rtol=0)  # 两套功率平衡和压降实现必须一致。
            self.assertGreaterEqual(equations.margin(np.empty(0),witness,state), -1e-12)  # 内点须满足未松弛的原始约束。
            self.assertEqual(reference.classify(witness)[0], 1)  # 内点还须通过独立 AC 检查。
            self.assertEqual(eta < 1e-7, reference.classify(power)[0] == 1)  # 在本组远离数值歧义的样本上核对查询判定。
            if cut is not None:  # 存在分离割时再核对其方向和有效性。
                cuts += 1  # 累计实际非空割。
                self.assertLess(cut[0]+cut[1:]@power, 0)  # 割必须严格排除产生它的查询点。
                self.assertTrue(np.all(cut[0]+valid@cut[1:] >= -1e-8))  # 已知 AC 可行样本不能被该割排除。
        self.assertGreater(cuts, 0)  # 保证测试实际覆盖了不可行查询和分离割。

    def test_regions_meet_the_same_stopping_rule(self):  # 三种方法均须达到各自声明的几何停止条件。
        for method in ('linear','socp','hybrid'):  # 分别运行纯线性、纯 SOCP 及混合流程。
            records, _ = workflow()['compute_regions']([network],method,.002)  # SOCP 两法使用相同 0.2% 径向精度。
            record = records[0]  # 此测试固定单个网架方案。
            inner, outer = np.array(record['inner']), np.array(record['outer'])  # 读取最终认证内域与候选外域。
            hull = halfspaces(inner)  # 把内域转为支持半空间。
            tau = 0. if method == 'linear' else .002  # 线性域用完整包含，SOCP 域用规定径向内缩。
            self.assertLessEqual(((1-tau)*outer@hull[:, :3].T+hull[:, 3]).max(), 1e-8)  # 检查全部外域顶点的缩放像均在内域中。
            if method != 'linear':  # SOCP 内点进一步交给独立 AC 核验。
                self.assertTrue(np.all(ACPowerFlow(network).classify(inner) == 1))  # 当前纯负荷算例中每个内点均应 AC 可行。
            else:  # 纯线性结果按自身模型核验，不当作 AC 证书。
                e = PlanningEquations(network,'linear',planning=False)  # 建立同一网架的线性物理方程。
                for point in outer:  # 逐一核验最终外域全部顶点。
                    state = e.restore(np.empty(0),point,np.zeros(len(e.upper)))  # 线性模型电流为零，直接重建功率和压降。
                    self.assertGreaterEqual(e.margin(np.empty(0),point,state),-1e-9)  # 原始线性约束须在数值容差内全部成立。

    def test_investment_changes_only_selected_impedances(self):  # 建设只能改变被选支路的阻抗及对应费用。
        self.assertEqual(len(network.designs),16)  # 四个二元升级项目恰有 16 个参考组合。
        for design in network.designs:  # 逐个检查参考方案的数据实例化。
            factor = np.ones(network.n)  # 非升级支路的阻抗系数保持 1。
            for selected, project in zip(design.x,network.projects):  # 将方案选择与唯一投资项目表对应。
                if selected:  # 仅选择建设的项目会改变支路参数。
                    factor[network.nodes.index(project.branch[1])] = .5  # 并联一回相同线路使 R/X 减半。
            np.testing.assert_allclose(design.r,network.r*factor)  # 检查所有支路电阻的改造规则。
            np.testing.assert_allclose(design.reactance,network.reactance*factor)  # 检查所有支路电抗的改造规则。
            np.testing.assert_array_equal(design.fixed_p,network.fixed_p)  # 建设不得改变有功背景工况。
            np.testing.assert_array_equal(design.fixed_q,network.fixed_q)  # 建设不得改变无功背景工况。
            self.assertEqual(design.cost,sum(x*p.cost for x,p in zip(design.x,network.projects)))  # 方案费用必须等于已选项目费用之和。

    def test_ac_minimum_cost_equals_all_designs_reference(self):  # 按投资排序的 AC 扫描须等价于完整方案逐点检查。
        points=np.random.default_rng(36).random((100,3))*[380,4500,600]  # 使用固定随机种子生成独立验证点。
        plans=network.designs  # 这里只为小算例核对使用完整方案表。
        expected=np.full(len(points),np.inf)  # 没有可行方案的点保持无限投资。
        for design in plans:  # 参考算法逐个方案检查全部点。
            status=ACPowerFlow(design).classify(points)  # 使用完整 AC 等式判定当前方案。
            self.assertFalse(np.any(status==0))  # 此测试样本须获得确定结论。
            expected[status==1]=np.minimum(expected[status==1],design.cost)  # 在全部 AC 可行方案中取最低投资。
        actual,_,_=workflow()['ac_reference'](plans,points,(0.,1.,2.,np.inf))  # 对比 Notebook 的按费用扫描及提前跳过流程。
        np.testing.assert_array_equal(actual,expected)  # 每个点的最小投资须逐项一致。

    def test_planning_union_equals_separately_constructed_lp_domains(self):  # 并集覆盖加速不能改变各方案域的真实并集。
        from region import simplex  # 独立单方案构域从同一有效外界开始。
        from vertify import region_membership  # 使用同一网格归属口径比较两种几何流程。
        plans=network.designs[:8]  # 选取八个方案覆盖并集重叠情况。
        flow=workflow()  # 读取主入口中唯一一套流程定义。
        records=flow['linear_region'](plans)  # 运行共享认证并集的线性构域流程。
        independent=[flow['cut_region'](PlanningSP(PlanningEquations(d,'linear',planning=False)),simplex(d),0.)['outer'] for d in plans]  # 每个方案独立完成全部线性顶点认证作为对照。
        points=np.random.default_rng(38).random((3000,3))*[450,5000,700]  # 在同一批 3000 个点上比较并集成员关系。
        np.testing.assert_array_equal(region_membership(points,[r['outer'] for r in records]),  # 取得并集覆盖流程返回的成员标签。
                                      region_membership(points,independent))  # 其结果须与逐方案独立构域完全相同。


if __name__ == '__main__':  # 允许直接运行本测试文件。
    unittest.main()  # 执行该文件内全部测试。
