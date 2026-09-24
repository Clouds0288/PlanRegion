# 配电网可行规划域

给定候选线路、设备和建设预算，求出**存在某个合法建设方案能够供电的负荷组合范围**。不同负荷组合允许选择不同方案。

主线采用两阶段算法：全局支持割收紧凸外包络，再用目标空间分支定界保留凸包内的非凸缺口。模型为 SOCP 规划模型；所有优化问题由 Gurobi 求解。

## 运行

安装 requirements.txt 中的依赖，配置有效 Gurobi 许可证，修改 main.py 顶部参数，再运行：

~~~text
python main.py
~~~

默认 Case33bw、8 个升级项目、负荷坐标 (18,25,33)、预算 (0,2,4)，全域认证容差 EPSILON_KW=20 kW。
LOAD_NODES=(18,25) 可直接运行二维；未选节点恢复 Case33 的原始背景负荷。
QUERY_TIME_LIMIT 是单次优化时限，CASE_TIME_LIMIT 是每档预算两个阶段共用的总时限。
其他网络选 fourbus 或 jiangkou 时，设置 BUDGETS=None 使用该网络的预算单位和默认预算。

每档预算只输出两个文件：

- results/planning/budget_*/result.json：内域、两阶段外域、割、认证状态和耗时。
- results/planning/budget_*/region.html：离线交互图，三维可旋转。

蓝色是第一阶段外包络，绿色是第二阶段外近似，浅绿色是同方案可行内域。
certified 表示整个外域到内域并集的无穷范数距离上界不超过 EPSILON_KW；unknown/time_limit 保留未完成结论。绿色外域不称精确 AC 真域。

## 调用关系

~~~text
main()
  按预算循环
    build_continuous_region(network, budget)
      Stage 1 → PlanningModel.solve(weights=omega)
      Stage 2 → PlanningModel.solve(target=q, cut_normals=W)
                region.py 几何裁剪、分支、认证
    save_result(result)
~~~

PlanningModel.solve 只执行一次 Gurobi optimize；没有 planning_query、PlanningOracle 或 solve_region 求解包装。
固定方案的连续热启动和完整整数查询均在 main.py 显式调用，固定方案界不会用于全局割。

| 文件 | 内容 |
|---|---|
| main.py | 配置与编号流程；1.1—1.4 支持割，2.1—2.6 非凸分支，3.1—3.2 最终认证/输出 |
| model.py | 统一变量和约束 (1)—(13)，支持/退让/投影缺额/固定负荷最低投资四种目标 |
| region.py | 二维/三维纯几何；不调用求解器 |
| plot.py | 核心 JSON 和最终交互图；无服务器和回放事件 |
| Network/ | 候选图、节点、类型、阻抗、容量与费用 |
| vertify.py | 独立 AC 校核，不进入构域求解链 |
| tests/ | 数值回归，以及保留的二维/三维实验和独立扫描审核 |

## 数学和接口

[符号契约](docs/notation.md) 是变量名、单位、状态切片和结果字段的唯一约定。
[算法说明](docs/continuous_region.md) 给出模型、割的证明和覆盖证书。

x/p/state 的含义及布局保持不变。几何统一使用 kW，结果升级为 schema_version=4。
旧 SP、剩余域模型、实时回放及专用基准代码已退出；历史登记保存在 [v1 契约](docs/notation-v1.md)。旧结果不自动加载、不覆盖；旧实验可通过 Git 历史及结果源码快照复现。

本算法利用纯负荷域的下闭包性质。发电反送、负电价诱导的非单调限制或新增最低供电约束，需要重新证明下闭包，不能直接套用认证结论。

## 验证

~~~text
python -m unittest tests.test_notation -v
python -m pytest -q
~~~

回归包括：两片区整数建设产生的非凸缺口、全局界方向、目标切换、热启动不截断预算可行集、二维/三维距离与独立 Gurobi LP 对照、节点覆盖的 min/max 顺序、全局析取传播、Case33 固定方案目标、反向潮流和独立 AC。

pytest 只收集 tests/，不收集 results/ 中的历史源码快照。
