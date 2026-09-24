"""独立 AC 等式、FR/MR、三态网格和唯一结果文件的回归验证。"""
import itertools  # 对少量未知标签的所有赋值检查误差区间。
import tempfile  # 测试文件仅写临时目录。
import unittest  # 标准回归框架。
import numpy as np  # 复相量、表面几何和标签数组。
import gurobipy as gp  # 独立非凸 AC 求解环境。
from Network.case33bw import Case33  # 当前正式网架。
from vertify import ACPowerFlow
from tests.reference import upgrade_plan


class ACReferenceTests(unittest.TestCase):  # 以代表网架核对独立 AC 的数值证书。
    @classmethod
    def setUpClass(cls):  # 不生成完整建设组合表。
        network = Case33(upgrade_count=8)  # 当前八候选线路。
        plans = [upgrade_plan(network, x)
                 for x in (np.zeros(8,dtype=int),np.arange(8)%2,np.ones(8,dtype=int))]
        cls.models = [ACPowerFlow(network.tree(network.encode_plan(plan)), threads=1) for plan in plans]
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
                self.assertEqual(actual, expected, (model.network.type_indices, power))
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



if __name__ == '__main__':
    unittest.main()
