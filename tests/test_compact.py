"""未选型号的零流锥回归。"""
import unittest  # 直接运行模型回归。
import numpy as np  # 构造固定查询和型号向量。
from threadpoolctl import threadpool_limits  # 所有比较保持单线程。
from Network.four_bus_five_corridor import FourBus
from model import GridPhysics, SubProblem  # 正式方程、主问题和连续 SP。
from tests.planning_checks import margin


class CompactPlanningTests(unittest.TestCase):  # 核对物理映射、最优值和全域割有效性。
    @classmethod
    def setUpClass(cls):  # 固定数值计算线程数。
        cls.threads = threadpool_limits(limits=1)  # 减少求解路径差异。

    @classmethod
    def tearDownClass(cls):  # 还原调用环境。
        cls.threads.restore_original_limits()  # 测试不永久更改线程配置。


    def test_sp_handles_zero_flow_cones_of_unselected_types(self):
        network = FourBus()
        equations = GridPhysics(network, 'socp')
        cases = (({'01': 'H', '12': None, '13': None, '02': 'L', '23': 'M'},
                  [47.47692657884853, 14.400743716323246, 7.840840158490067], True),
                 ({'01': 'H', '12': None, '13': 'L', '02': 'H', '23': None},
                  [82.83826215331091, 41.51433150869895, 11.864741576952936], False))
        for plan, point, feasible in cases:
            x, power = network.encode_plan(plan), np.array(point)
            checked = SubProblem(equations, threads=1).solve(x, power)
            self.assertEqual(checked['feasible'], feasible)
            if feasible:
                self.assertGreaterEqual(margin(equations, x, power, checked['state']), -1e-8)
            else:
                cut = checked['cut']
                self.assertIsNotNone(cut)
                self.assertLess(cut[0]+cut[1:4]@power+cut[4:]@x, -1e-9)


if __name__=='__main__':  # 支持单独执行本组回归。
    unittest.main()  # 不启动 Notebook 的完整实验。
