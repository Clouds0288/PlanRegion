"""四／八／十六线路共用主线、三态分类以及测试侧硬时限检查。"""
import ast  # 核对正式 Notebook 与测试设施的依赖边界。
import os
from pathlib import Path  # 临时测试夹具路径。
from tempfile import TemporaryDirectory  # 不污染正式实验输出。
import unittest  # 标准测试框架。
from unittest.mock import patch  # 模拟未获证的 SOCP 查询，检查混合阶段不会误认证。
import nbformat  # 读取正式入口和创建短测试夹具。
import numpy as np  # 统一网格和预算标签。
from threadpoolctl import threadpool_limits  # 所有比较单线程。
from Network.case33bw import Case33  # 三档嵌套候选线路配置。
from model import PlanningEquations, PlanningModel  # 完整规划模型作为直接求解对照。
from plot import sample_region
from vertify import ac_planning_query, validate_ac_region
import main
from tests.watchdog import run_guarded  # 只有测试导入看门狗。
from tests.reference import fixed_topology


class ScaleTests(unittest.TestCase):  # 规模扩展只改变网架设置。
    @classmethod
    def setUpClass(cls):  # 固定矩阵线程数。
        cls.threads = threadpool_limits(limits=1)  # 减少数值路径差异。

    @classmethod
    def tearDownClass(cls):  # 还原线程设置。
        cls.threads.restore_original_limits()  # 不影响其他实验。

    @unittest.skipUnless(os.environ.get('PLANREGION_SLOW_TESTS') == '1', '四／八／十六候选全预算连续构域慢测试')
    def test_grid_matches_direct_point_queries(self):  # 连续域采样逐点对照完整 MILP/MISOCP。
        budgets,bounds = np.array([0.,1.,2.,np.inf]),np.array([540.,5950.,970.])  # 同一公共评价箱。
        points = (np.indices((4,)*3).reshape(3,-1).T+.5)*bounds/4  # 小网格的 64 个中心独立检查。
        for count in (4,8,16):  # 同一连续构域流程处理三种规模；这是较慢的全规模回归。
            network = Case33(upgrade_count=count)  # 仅物理配置变化。
            def build(method):  # 看门狗从测试外部约束正式构域函数。
                states = []
                for budget in budgets:
                    run = run_guarded('build_continuous_region',dict(network=network,method=method,budget=budget,bounds=bounds),900.)
                    self.assertEqual(run['termination'],'returned',run['error'])
                    states.append(sample_region(run['value'],points,bounds))
                return np.asarray(states).reshape((len(budgets),4,4,4))
            expected = {}  # 保存两类直接模型的参考标签。
            for method in ('linear','socp'):  # LP/SOCP 分别对照。
                equations,costs = PlanningEquations(network,method),[]  # 共用固定系数。
                for point in points:  # 测试逐点求解，不复用单调分类结论。
                    problem = PlanningModel(equations,power=point, threads=1)  # 完整运行方程与逐线路变量。
                    with problem.model:  # 每次直接查询及时释放。
                        answer = problem.solve()  # 参考必须完成最优性或不可行证明。
                    self.assertTrue(answer is None or answer['status']=='optimal')  # 未确定不算参考真值。
                    costs.append(np.inf if answer is None else answer['objective'])  # 原始最小投资。
                costs = np.asarray(costs)  # 一次投资结果用于所有预算。
                expected[method] = np.array([np.where(np.isfinite(costs)&(costs<=b),1,-1) for b in budgets]).reshape((4,4,4,4))  # 无限预算仍排除物理不可行点。
                actual = build(method)  # 正式构域在测试看门狗下执行。
                known = actual != 0
                self.assertTrue(known.any())
                np.testing.assert_array_equal(actual[known],expected[method][known])
            hybrid = build('hybrid')  # 自己承担 LP 阶段。
            known = hybrid != 0
            self.assertTrue(known.any())
            np.testing.assert_array_equal(hybrid[known],expected['socp'][known])
            ac = validate_ac_region(network,budgets,4,bounds, threads=1)
            costs = []  # 独立逐点执行 AC 方案搜索，核对成块推断。
            equations = PlanningEquations(network,'socp')  # AC 搜索只借助 SOCP 候选。
            for point in points:  # 每个中心独立验证，不继承整块标签。
                answer = ac_planning_query(equations,point, threads=1)  # 独立 AC 等式认证。
                self.assertTrue(answer is None or answer['status']=='optimal')  # 当前小样本须完整获证。
                costs.append(np.inf if answer is None else answer['objective'])  # 记录该点最小 AC 投资。
            costs = np.asarray(costs)  # 按各预算生成独立参考标签。
            reference = np.array([np.where(np.isfinite(costs)&(costs<=b),1,-1) for b in budgets]).reshape(ac.shape)  # 只用于回归对照。
            np.testing.assert_array_equal(ac,reference)  # 单调推断不改变逐点 AC 结果。
            self.assertTrue(np.all((ac!=1)|(expected['socp']==1)))  # AC 必须包含在 SOCP 外松弛内。

    def test_hybrid_rechecks_linear_interior(self):  # SOCP 没有证书时不得继承 LP 内点。
        original = PlanningModel.solve  # LP 阶段使用真实优化器。
        def query(problem, *args, **kwargs):  # 测试只替换 SOCP 的认证结果。
            if problem.equations.method=='socp':  # 模拟没有候选和证书。
                return dict(feasible=False, x=None, p=None, bound=None, status='unknown')
            return original(problem, *args, **kwargs)  # LP 阶段仍真实切割。
        with patch.object(PlanningModel, 'solve', new=query):
            bounds = np.array([540.,5950.,970.])
            domain = main.build_continuous_region(fixed_topology(Case33(upgrade_count=8)), 'hybrid', 0., bounds, threads=1)
            points = (np.indices((4,)*3).reshape(3,-1).T+.5)*bounds/4
            states = sample_region(domain,points,bounds)
        self.assertTrue(np.any(states==0))  # 尚待 SOCP 认证的区域仍为未知。
        self.assertTrue(np.all(states<=0))  # 没有任何 SOCP 证书就不能发布域内点。

    def test_notebook_does_not_depend_on_test_orchestration(self):  # 测试设施与正式流程保持明确边界。
        from pathlib import Path
        source=Path('main.py').read_text(encoding='utf-8')  # 当前正式入口不依赖已删除的 Notebook。
        tree=ast.parse(source)  # 使用语法结构而非导入副作用。
        imports=[n.module or '' for n in ast.walk(tree) if isinstance(n,ast.ImportFrom)]  # 收集所有 from 导入。
        self.assertFalse(any(m=='tests' or m.startswith('tests.') for m in imports))  # 主入口不得导入看门狗或审核实现。
        self.assertNotIn('run_guarded',source)  # 不保留隐式进程包装。
        # 展示回调现位于构域函数内部；依赖边界按导入检查，不按局部函数名猜测。
        self.assertNotIn('joint_benders',source)  # 旧算法只保留在对照测试。

    def test_solver_timeout_remains_unknown(self):  # 优化器自身时限不等于不可行证明。
        problem = PlanningModel(PlanningEquations(Case33(upgrade_count=8),'socp'),power=[300.,3000.,500.], threads=1)  # 正常规划查询。
        with problem.model:  # 测试后释放优化器。
            answer = problem.solve(time_limit=0.)  # 在搜索开始前触发优化器时限。
        self.assertEqual(answer['status'],'unknown')  # 保留未确定。
        self.assertFalse(answer['feasible'])  # 不得产生伪运行证书。

    def test_external_watchdog_handles_return_error_and_timeout(self):  # 被测函数无需 deadline、report 或检查点接口。
        with TemporaryDirectory() as folder:  # 所有夹具仅用于这个测试。
            path = Path(folder)/'probe.ipynb'  # 小型 Notebook 与正式函数使用相同加载方式。
            source = "from pathlib import Path\nfrom time import sleep\ndef echo(value):\n    return value\ndef hang(marker):\n    Path(marker).write_text('started')\n    sleep(30)\ndef fail():\n    raise ValueError('probe failure')"  # 只有普通函数。
            nbformat.write(nbformat.v4.new_notebook(cells=[nbformat.v4.new_code_cell(source,metadata={'tags':['test-flow']})]),path)  # 写入临时夹具。
            returned = run_guarded('echo',dict(value=17),5.,notebook=path)  # 正常返回不依赖报告通道。
            self.assertEqual((returned['termination'],returned['value']),('returned',17))  # 原样保留结果。
            failed = run_guarded('fail',{},5.,notebook=path)  # 异常原文必须返回测试侧。
            self.assertEqual(failed['termination'],'error')  # 不将错误归为超时或不可行。
            self.assertIn('probe failure',failed['error'])  # 异常堆栈包含真实原因。
            marker = Path(folder)/'started.txt'  # 区分函数确实开始与仅导入超时。
            stopped = run_guarded('hang',dict(marker=str(marker)),3.,notebook=path)  # 从测试侧限制总运行时间。
            self.assertTrue(marker.exists())  # 确保实际终止了被卡住的函数。
        self.assertEqual(stopped['termination'],'hard_timeout')  # 硬超时独立于主线代码。
        self.assertIsNone(stopped['value'])  # 不伪造中断函数的结果。
        self.assertLess(stopped['wall_seconds'],5.)  # 不能等待 30 秒自然结束。


if __name__=='__main__':  # 支持直接运行测试。
    unittest.main()  # 正式主线不导入本文件。
