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


if __name__=='__main__':  # 支持单独运行物理数据检查。
    unittest.main()  # 执行当前文件测试。
