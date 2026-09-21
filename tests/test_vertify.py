"""独立 AC 等式、FR/MR、三态网格和唯一结果文件的回归验证。"""
import itertools  # 对少量未知标签的所有赋值检查误差区间。
import tempfile  # 测试文件仅写临时目录。
import unittest  # 标准回归框架。
import numpy as np  # 复相量、表面几何和标签数组。
import gurobipy as gp  # 独立非凸 AC 求解环境。
from Network.case33bw import Case33  # 当前正式网架。
from plot import voxel_surface  # 核对颜色区域的真实体素表面。
from vertify import ACPowerFlow
from main import BenchmarkResult
from vertify import disagreement, disagreement_interval, METHODS


class ACReferenceTests(unittest.TestCase):  # 以代表网架核对独立 AC 的数值证书。
    @classmethod
    def setUpClass(cls):  # 不生成完整建设组合表。
        network = Case33(candidate_count=8)  # 当前八候选线路。
        cls.models = [ACPowerFlow(network.design(x)) for x in (np.zeros(8,dtype=int),np.arange(8)%2,np.ones(8,dtype=int))]  # 基础、交错、全升级网架。
        cls.environment = gp.Env(empty=True)  # 静默的共享全局求解环境。
        cls.environment.setParam('OutputFlag',0)  # 不输出逐点求解日志。
        cls.environment.start()  # 在测试计时之外启动环境。

    @classmethod
    def tearDownClass(cls):  # 释放按需创建的非凸模型。
        for model in cls.models:  # 只处理本组测试的实例。
            model.close()  # 正常回收求解资源。
        cls.environment.dispose()  # 最后释放共享环境。

    def test_fixed_point_matches_explicit_global_ac(self):  # 递推证书与完整 AC 等式全局求解对照。
        counts = {-1:0,1:0}  # 确认可行与不可行都被覆盖。
        for model in self.models:  # 三个物理上不同的固定网架。
            for power in ([100.,800.,150.],[300.,3000.,500.],[500.,6000.,950.]):  # 低、中、高负荷查询。
                actual = int(model.classify(power)[0])  # 独立不动点证书。
                expected = model.global_status(power,self.environment)  # 显式完整 AC 电流等式。
                self.assertNotEqual(expected,0)  # 未完成结论不能作为参考真值。
                self.assertEqual(actual,expected,(model.network.x,power))  # 两套求解路径必须一致。
                counts[actual] += 1  # 统计实际覆盖类别。
        self.assertGreater(counts[-1],0)  # 确实检查域外。
        self.assertGreater(counts[1],0)  # 确实检查域内。

    def test_ac_witness_satisfies_complex_nodal_power_flow(self):  # 将支路证书重建为节点复相量。
        power = np.array([100.,800.,150.])  # 代表网架均可承载的低负荷。
        for model in self.models:  # 不同选型均需满足同一复功率关系。
            c = model.network  # 当前固定物理网架。
            status,ell = model.classify(power,return_currents=True)  # 请求完整电流平方证书。
            self.assertEqual(status[0],1)  # 先取得真实可行证书。
            P,Q,v,_ = model.state(power,ell)  # 通过独立支路递推恢复状态。
            voltage,current = np.ones(c.n,dtype=complex),np.zeros(c.n,dtype=complex)  # 根电压为单位相量。
            for i in c.order:  # 由根到叶依次恢复电流和电压。
                parent = 1+0j if c.parent[i]<0 else voltage[c.parent[i]]  # 已知送端复电压。
                current[i] = np.conj((P[0,i]+1j*Q[0,i])/parent)  # S=U*conj(I)。
                voltage[i] = parent-(c.r[i]+1j*c.reactance[i])*current[i]  # 复数欧姆定律。
            demand = voltage*np.conj(current-np.array([current[j].sum() for j in c.children]))  # 节点净功率等于入线电流减去出线电流。
            p,q = c.loads(power)  # 当前背景和独立负荷的完整节点注入。
            np.testing.assert_allclose(demand,(p+1j*q)[0],atol=1e-9,rtol=0)  # 复数功率守恒。
            np.testing.assert_allclose(abs(voltage)**2,v[0],atol=1e-9,rtol=0)  # 相量幅值平方与模型电压一致。


class RegionComparisonTests(unittest.TestCase):  # 三态标签与物理误差分开处理。
    def test_symmetric_difference_counts_both_sides(self):  # 误差定义必须同时计入多余和遗漏。
        labels = np.array([0, 0, 1, 2, 2, 3, 3, 3])  # 构造已知计数的四类标签。
        row = disagreement((labels & 1)>0, (labels & 2)>0)  # 从两种区域成员位生成统计结果。
        self.assertAlmostEqual(row['region_error_percent'], 50.)  # 对称差为 3、并集为 6，不一致率应为 50%。
        self.assertAlmostEqual(row['fr_percent'], 25.)  # 多余 1、计算域 4，FR 应为 25%。
        self.assertAlmostEqual(row['mr_percent'], 40.)

    def test_empty_region_has_undefined_conditional_rate(self):  # 空域的条件比例不能伪造为零。
        for label, fr, mr in [(0, None, None), (1, 100., None), (2, None, 100.)]:  # 分别覆盖两域皆空、仅计算域非空和仅 AC 域非空。
            labels = np.full(8, label)  # 生成当前情形的统一标签。
            row = disagreement((labels & 1)>0, (labels & 2)>0)  # 按正式误差函数计算结果。
            self.assertEqual(row['fr_percent'], fr)  # FR 分母为空时必须保持无定义。
            self.assertEqual(row['mr_percent'], mr)

    def test_surface_preserves_cavity_and_volume(self):  # 体素表面应保留内部空腔并给出正确有向体积。
        mask = np.ones((3, 3, 3), dtype=bool)  # 先构造一个实心 3×3×3 体素块。
        mask[1, 1, 1] = False  # 挖去中心体素，形成内部空腔。
        h = np.array([.75, 1.5, 2.])  # 三个坐标使用不同步长，检验物理尺度换算。
        mesh = voxel_surface(mask, h)  # 使用正式绘图表面算法生成三角网格。
        points = np.asarray(mesh['vertices'])  # 读取实际表面顶点坐标。
        triangles = points[np.asarray(mesh['triangles'])]  # 按索引恢复每个三角形的三个顶点。
        volume = np.einsum('ij,ij->i', triangles[:, 0],  # 由有向三角面计算封闭表面围成的体积。
                           np.cross(triangles[:, 1], triangles[:, 2])).sum()/6  # 标量三重积之和除以 6 得到带符号体积。
        self.assertAlmostEqual(volume, mask.sum()*h.prod())

    def test_unknown_rates_bound_every_label_completion(self):  # 未确定不允许被默认算作不可行。
        left,right = np.array([1,0,-1,0]),np.array([0,-1,1,0])  # 含已知和未知的两个区域。
        interval = disagreement_interval(left,right)  # 正式 FR/MR 保守区间。
        for bits in itertools.product((-1,1),repeat=int((left==0).sum()+(right==0).sum())):  # 逐一检查少量未知标签的可能真值。
            a,b = left.copy(),right.copy()  # 各次赋值相互独立。
            count = int((a==0).sum())  # 前半部分对应计算域。
            a[a==0],b[b==0] = bits[:count],bits[count:]  # 后半部分对应 AC 域。
            rates = disagreement(a==1,b==1)  # 该确定情形的真实比例。
            for key in ('fr','mr'):  # 两个指标都须被区间覆盖。
                value = rates[key+'_percent']  # 分母为空时比例无定义。
                if value is not None:  # 只检查有定义的可能比例。
                    self.assertLessEqual(interval[key+'_interval'][0],value)  # 下界不得偏高。
                    self.assertGreaterEqual(interval[key+'_interval'][1],value)  # 上界不得偏低。

    def test_result_stores_only_states_and_metadata(self):  # 不保存方案表、重复掩码或派生指标。
        states = np.random.default_rng(5).choice([-1,0,1],(4,2,3,3,3)).astype(np.int8)  # 包含未知的四方法、两预算网格。
        metadata = dict(divisions=3,budgets=[0.,None],bounds=[1.,2.,3.],load_nodes=[18,25,33],seconds={m:float(i+1) for i,m in enumerate(METHODS)})  # 配置和总计时各一份。
        result = BenchmarkResult(states,metadata)  # 简化后的唯一结果容器。
        with tempfile.TemporaryDirectory(prefix='planregion-test-') as folder:  # 不污染正式结果。
            result.save(folder)  # 写入唯一文件。
            restored = BenchmarkResult.load(folder)  # 从相同格式恢复。
            with np.load(f'{folder}/result.npz',allow_pickle=False) as data:  # 检查实际文件字段。
                self.assertEqual(set(data.files),{'states','metadata'})  # 文件中只有两种基础数据。
        np.testing.assert_array_equal(restored.states,states)  # 未知与域内外均原样保存。
        self.assertEqual(restored.metadata,metadata)  # 配置及方法总耗时一致。
        np.testing.assert_array_equal(restored.labels,result.labels)  # 差集颜色按需推导且一致。
        self.assertEqual(len(restored.summary),8)  # 四方法和两预算各生成一行。
        self.assertTrue(np.any(restored.labels==4))  # 未确定区域使用灰色。


if __name__=='__main__':  # 支持单独运行指标与 AC 检查。
    unittest.main()  # 不执行主 Notebook 的完整实验。
