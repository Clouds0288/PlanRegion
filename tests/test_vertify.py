"""完整 AC、区域差集和唯一结果文件的回归验证。"""
import tempfile  # 结果往返测试写入临时目录，不污染实验数据。
import unittest  # 使用标准测试框架。
import numpy as np  # 处理相量、体素与实验数组。
import gurobipy as gp  # 用显式非凸模型交叉核验 AC 判定。
from Network.four_bus_five_corridor import network  # 读取小算例全部合法建设方案。
from plot import voxel_surface  # 检查显示表面没有填平真实孔洞。
from model import ACPowerFlow  # 被测 AC 实现独立于线性及 SOCP 方程。
from vertify import BenchmarkResult, disagreement  # 测试误差定义与唯一结果文件的存取。


class ACReferenceTests(unittest.TestCase):  # 覆盖独立 AC 判定、相量重建及送端容量。
    @classmethod  # 整组测试共用方案列表与求解环境。
    def setUpClass(cls):  # 按每个合法小网架建立独立 AC 接口。
        cls.models = [ACPowerFlow(d) for d in network.designs]  # 四节点全部 216 个方案均纳入验证。
        cls.environment = gp.Env(empty=True)  # 创建静默的共享非凸求解环境。
        cls.environment.setParam('OutputFlag', 0)  # 关闭逐点优化日志。
        cls.environment.start()  # 在执行比较前完成环境初始化。

    @classmethod  # 整组测试结束后统一释放资源。
    def tearDownClass(cls):  # 释放所有按需创建的 AC 模型。
        for model in cls.models:  # 逐一处理全部方案接口。
            model.close()  # 未创建全局模型的接口也可以正常关闭。
        cls.environment.dispose()  # 最后释放共享 Gurobi 环境。

    def test_all_designs_against_explicit_global_ac(self):  # 批量不动点判定须与显式非凸 AC 一致。
        rng = np.random.default_rng(20260919)  # 固定随机种子保证验证负荷可复现。
        counts = {-1: 0, 1: 0}  # 分别统计确定不可行与确定可行的覆盖情况。
        for model in self.models:  # 逐个合法径向方案执行交叉验证。
            power = rng.dirichlet(np.ones(4))[:3]*network.power_limit  # 在总负荷初始外界内选取随机点。
            actual = int(model.classify(power)[0])  # 取得独立不动点算法的三态判定。
            expected = model.global_status(power, self.environment)  # 用显式 AC 等式的全局求解器取得参考判定。
            self.assertNotEqual(expected, 0, 'QCQP 未完成判定')  # 测试不接受未确定结论作为参考。
            self.assertEqual(actual, expected, (model.network.choice, power))  # 两个算法须给出相同的可行性结论。
            counts[actual] += 1  # 累计对应类别的实际样本数。
        self.assertEqual(len(self.models), 216)  # 确认完整覆盖了全部 216 个建设方案。
        self.assertGreater(counts[-1], 0)  # 样本中必须实际包含不可行情况。
        self.assertGreater(counts[1], 0)  # 样本中也必须实际包含可行情况。

    def test_ac_witness_satisfies_complex_nodal_power_flow(self):  # 从电流平方证书重建复数相量，核对节点功率。
        power = np.array([5., 7., 4.])  # 使用一个在代表方案中可行的低负荷点。
        for model in self.models[::9]:  # 跨方案选择代表网架进行相量重建。
            c = model.network  # 取得该固定方案的物理阻抗和拓扑。
            status, ell = model.classify(power, return_currents=True)  # 请求完整 AC 电流证书。
            self.assertEqual(status[0], 1)  # 用于相量重建的查询须已认证可行。
            P, Q, v, _ = model.state(power, ell)  # 独立重建送端功率和节点电压平方。
            voltage, current = np.ones(c.n, dtype=complex), np.zeros(c.n, dtype=complex)  # 根电压从单位相量开始，支路电流初始为零。
            for i in c.order:  # 沿根到叶计算支路电流与节点复电压。
                parent = 1+0j if c.parent[i] < 0 else voltage[c.parent[i]]  # 送端电压取根相量或已经计算的父节点相量。
                current[i] = np.conj((P[0, i]+1j*Q[0, i])/parent)  # 由 S=U·I* 恢复复支路电流。
                voltage[i] = parent-(c.r[i]+1j*c.reactance[i])*current[i]  # 由复阻抗压降 U_j=U_i−Z·I 恢复节点电压。
            demand = voltage*np.conj(current-np.array([current[j].sum() for j in c.children]))  # 节点净注入由本支路电流减去子支路电流计算。
            p, q = c.loads(power)  # 读取同一负荷点的完整节点 P/Q。
            np.testing.assert_allclose(demand, (p+1j*q)[0], atol=1e-10, rtol=0)  # 复数节点功率必须等于给定负荷。
            np.testing.assert_allclose(abs(voltage)**2, v[0], atol=1e-10, rtol=0)  # 复电压幅值平方须与支路状态一致。

    def test_transformer_checks_sending_end_power(self):  # 源端容量应包含网损，不能只检查负荷总和。
        power = np.full(3, network.power_limit/3)  # 将总负荷放在无损有功容量边界上。
        for model in self.models:  # 所有方案都有正损耗，源端容量应越限。
            self.assertEqual(model.classify(power)[0], -1)  # 独立 AC 必须判为不可行。


class RegionComparisonTests(unittest.TestCase):  # 检查 FR/MR、表面几何与去重存储。
    def test_symmetric_difference_counts_both_sides(self):  # 误差定义必须同时计入多余和遗漏。
        labels = np.array([0, 0, 1, 2, 2, 3, 3, 3])  # 构造已知计数的四类标签。
        row = disagreement((labels & 1)>0, (labels & 2)>0)  # 从两种区域成员位生成统计结果。
        self.assertAlmostEqual(row['region_error_percent'], 50.)  # 对称差为 3、并集为 6，不一致率应为 50%。
        self.assertAlmostEqual(row['fr_percent'], 25.)  # 多余 1、计算域 4，FR 应为 25%。
        self.assertAlmostEqual(row['mr_percent'], 40.)  # 遗漏 2、AC 域 5，MR 应为 40%。

    def test_empty_region_has_undefined_conditional_rate(self):  # 空域的条件比例不能伪造为零。
        for label, fr, mr in [(0, None, None), (1, 100., None), (2, None, 100.)]:  # 分别覆盖两域皆空、仅计算域非空和仅 AC 域非空。
            labels = np.full(8, label)  # 生成当前情形的统一标签。
            row = disagreement((labels & 1)>0, (labels & 2)>0)  # 按正式误差函数计算结果。
            self.assertEqual(row['fr_percent'], fr)  # FR 分母为空时必须保持无定义。
            self.assertEqual(row['mr_percent'], mr)  # MR 分母为空时必须保持无定义。

    def test_surface_preserves_cavity_and_volume(self):  # 体素表面应保留内部空腔并给出正确有向体积。
        mask = np.ones((3, 3, 3), dtype=bool)  # 先构造一个实心 3×3×3 体素块。
        mask[1, 1, 1] = False  # 挖去中心体素，形成内部空腔。
        h = np.array([.75, 1.5, 2.])  # 三个坐标使用不同步长，检验物理尺度换算。
        mesh = voxel_surface(mask, h)  # 使用正式绘图表面算法生成三角网格。
        points = np.asarray(mesh['vertices'])  # 读取实际表面顶点坐标。
        triangles = points[np.asarray(mesh['triangles'])]  # 按索引恢复每个三角形的三个顶点。
        volume = np.einsum('ij,ij->i', triangles[:, 0],  # 由有向三角面计算封闭表面围成的体积。
                           np.cross(triangles[:, 1], triangles[:, 2])).sum()/6  # 标量三重积之和除以 6 得到带符号体积。
        self.assertAlmostEqual(volume, mask.sum()*h.prod())  # 结果须等于剩余体素数乘单体素体积。

    def test_result_roundtrip_preserves_masks_and_deduplicates_geometry(self):  # 保存再读取后须保持成员标签且不重复记录几何。
        masks = np.random.default_rng(4).random((4, 2, 3, 3, 3)) > .5  # 生成四方法、两预算的模拟布尔标签。
        ac_cost = np.random.default_rng(5).choice([0., 20000., 40000., np.inf], (3,3,3))  # 每个网格点仅保存一个模拟 AC 最小投资。
        masks[2,0] = ac_cost<=20000.  # 有限预算 AC 域由该最小投资生成。
        masks[2,1] = np.isfinite(ac_cost)  # 无限预算仍须排除投资为 inf 的不可行点。
        record = dict(inner=[[0, 0, 0]], outer=[[1, 2, 3]])  # 构造被两个预算共同引用的同一几何记录。
        result = BenchmarkResult(masks, dict(divisions=3, budgets=[20000., None]),  # 建立标准结果容器及预算配置。
                                 {'socp': [[record], [record]]}, ac_cost)  # 两档预算引用相同几何，测试存储去重。
        with tempfile.TemporaryDirectory(prefix='planregion-test-') as folder:  # 文件仅写在自动清理的专用临时目录。
            result.save(folder)  # 按正式格式保存一次实验。
            restored = BenchmarkResult.load(folder)  # 从唯一压缩文件重新加载。
            np.testing.assert_array_equal(restored.masks, masks)  # 所有方法与预算的成员标签必须完整恢复。
            self.assertEqual(restored.regions, result.regions)  # 逐方案几何及预算引用关系必须恢复。
            import json  # 解析文件中唯一的 JSON 几何记录表。
            with np.load(f'{folder}/result.npz') as data:  # 检查实际存储内容，而非只核对内存容器。
                self.assertEqual(set(data.files), {'masks', 'ac_cost', 'record'})  # 只允许三种基础数据字段，不重复保存派生标签。
                np.testing.assert_array_equal(data['ac_cost'], ac_cost)  # AC 最小投资数组应原样保存。
                self.assertEqual(len(json.loads(str(data['record']))['geometry']), 1)  # 重复引用的几何在文件中只能出现一次。


if __name__ == '__main__':  # 允许直接运行验证与存储测试。
    unittest.main()  # 执行本文件的全部检查。
