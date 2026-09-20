"""跨网架的 SOCP 状态、内点、有效割与原线性认证区域回归验证。"""
import json  # 读取已有实验中保存的独立回归记录。
from pathlib import Path  # 定位历史结果页面。
import re  # 提取页面内嵌的原始实验数据。
import unittest  # 使用标准单元测试框架。
import numpy as np  # 生成固定样本并检验物理状态。
from Network.four_bus_five_corridor import network  # 读取四节点原始小算例及其参考方案。
from model import PlanningEquations, PlanningSP, ACPowerFlow  # 使用统一 SP 与独立 AC 实现。
from region import halfspaces, simplex  # 构造初始外界并检查几何包含关系。
from vertify import region_membership  # 比较完整并集的点归属。
from notebook_flow import workflow  # 测试直接执行 Notebook 中的正式构域函数。


class SOCPTests(unittest.TestCase):  # 跨全部小网架检查锥状态和原线性区域回归。
    def test_all_designs_and_dual_cuts_against_ac(self):  # 每个建设方案都检查可行内点及有效割。
        rng = np.random.default_rng(20260920)  # 固定种子以便复现不可行和可行查询。
        cuts = feasible_queries = 0  # 确认测试实际覆盖两类查询结果。
        for design in network.designs:  # 遍历四节点全部 216 个合法方案。
            e = PlanningEquations(design,'socp',planning=False)  # 固定方案没有离散建设变量。
            oracle,reference = PlanningSP(e),ACPowerFlow(design)  # SOCP 与 AC 只共享网架数据。
            origin = oracle.solve(np.empty(0),np.zeros(3))['state']  # 取得本方案原点的可行运行状态。
            validation = rng.dirichlet(np.ones(4), size=30)[:, :3]*network.power_limit  # 在总负荷单纯形内生成验证点。
            valid = validation[reference.classify(validation) == 1]  # 用独立 AC 筛选已知可行样本。
            for power in (np.array([2., 3., 1.]), validation[0]):  # 低负荷和随机查询覆盖不同运行约束。
                result = workflow()['fixed_query'](oracle,power,origin)  # 调用正式的查询与内点构造流程。
                eta,cut,witness,state = (result[k] for k in ('eta','cut','witness','state'))  # 原查询与生成的内点分别取出。
                P,Q,v,ell = (state[s] for s in (e.P,e.Q,e.v,e.ell))  # 根据统一切片取得完整内点状态。
                u = np.r_[1.,v][design.parent+1]  # 从父节点电压恢复真实支路送端电压。
                for expected, value in zip(reference.state(witness, ell), (P, Q, v, u)):  # 使用独立支路递推核对全部状态分量。
                    np.testing.assert_allclose(value, expected[0], atol=1e-12, rtol=0)  # 功率平衡和压降的实现差须在舍入量级。
                self.assertLessEqual(np.max(P*P+Q*Q-u*ell), 1e-9)  # 各支路电流锥必须满足原始不等式。
                self.assertLessEqual(reference.violation(P[None], Q[None], v[None]).max(), 1e-9)  # 内点状态须满足电压和容量限值。
                self.assertEqual(reference.classify(witness)[0], 1)  # 独立 AC 也应认证这些内点。
                self.assertEqual(eta <= 1e-7, reference.classify(power)[0] == 1)  # 本组测试点的辅助违反量与独立 AC 判定一致。
                if cut is not None:  # 存在分离割时检查方向及不误删性质。
                    cuts += 1  # 累计实际产生的有效割。
                    self.assertLess(cut[0]+cut[1:]@power, 0)  # 割须严格排除原查询。
                    self.assertTrue(np.all(cut[0]+valid@cut[1:] >= -1e-7))  # 独立 AC 可行样本不得被误删。
                else:  # 没有割的样本用于覆盖可行查询路径。
                    feasible_queries += 1  # 统计该组无割查询。
        self.assertGreater(cuts, 0)  # 避免测试只覆盖可行点而未实际检验割。
        self.assertGreater(feasible_queries, 0)  # 避免测试只覆盖不可行点而未检验可行查询。

    def test_radial_stopping_certificate(self):  # 停止必须由完整内外域包含证书决定。
        for design in network.designs[::43]:  # 跨费用和拓扑选取多个代表方案。
            result = workflow()['cut_region'](PlanningSP(PlanningEquations(design,'socp',planning=False)),simplex(network),.002)  # 对每个方案以 0.2% 径向精度构造 SOCP 域。
            inside = halfspaces(result['inner'])  # 取得认证内域的全部支持半空间。
            self.assertLessEqual(np.max((1-.002)*result['outer']@inside[:, :3].T  # 把外域按规定比例向原点缩放。
                                        + inside[:, 3]), 1e-8)  # 所有缩放顶点必须属于内域。
            self.assertTrue(np.all(ACPowerFlow(design).classify(result['inner']) == 1))  # 认证内点逐一通过独立 AC 检查。

    def test_zero_load_subtree_does_not_collapse_witness(self):  # 空负荷子树不能因数值伪损耗把内点缩到原点。
        e = PlanningEquations(network.designs[13],'socp',planning=False)  # 选择曾出现零负荷子树退化的固定方案。
        oracle = PlanningSP(e)  # 运行正式统一子问题。
        origin = oracle.solve(np.empty(0),np.zeros(3))['state']  # 取得用于插值的可行原点。
        result = workflow()['fixed_query'](oracle,np.array([32.8,0.,0.]),origin)  # 只增加一个负荷节点，其余节点需求保持零。
        witness,currents = result['witness'],result['state'][e.ell]  # 提取所得内点及其支路电流证书。
        self.assertGreater(witness[0], 32.)  # 内点应保留该方向的实际承载能力。
        np.testing.assert_array_equal(currents[1:], 0.)  # 完全无负荷的下游支路应没有电流。

    def test_extracted_linear_solver_reproduces_saved_region(self):  # 清理后的统一线性 SP 须复现已保存的区域。
        plans = [d for d in network.designs if d.cost<=20000.]  # 使用历史 20000 元预算对应的相同方案集合。
        records = workflow()['linear_region'](plans)  # 以当前 Notebook 流程重新构造线性并集。
        html = Path('results/four_bus_five_corridor.html').read_text(encoding='utf-8')  # 历史原始记录作为独立的回归基准。
        data = json.loads(re.search(r'<script id="experiment" type="application/json">(.*?)</script>', html, re.S)[1])  # 提取 HTML 中保存的实验 JSON。
        scenario = next(s for s in data['scenarios'] if s['budget'] == 20000.)  # 定位同一预算对应的历史场景。
        groups = {}  # 按完整建设方案隔离历史认证点。
        for row in scenario['history']:  # 只读取当时真实执行的 SP 记录。
            if row['x'] is not None and row['eta'] is not None and row['cut'] is None:  # 按旧记录格式识别已认证的原查询。
                groups.setdefault(tuple(row['x']), []).append(row['p'])  # 同一方案的可行点才能取凸包。
        points = np.random.default_rng(51).dirichlet(np.ones(4), size=2000)[:, :3]*network.power_limit  # 用固定的公共负荷样本比较两次完整区域。
        np.testing.assert_array_equal(region_membership(points, [np.array(p) for p in groups.values()]),  # 由历史同方案认证凸包计算并集标签。
                                      region_membership(points, [r['outer'] for r in records if len(r['outer'])]))  # 当前外域并集标签必须逐项一致。


if __name__ == '__main__':  # 允许直接执行跨网架 SOCP 回归。
    unittest.main()  # 运行此文件的全部测试。
