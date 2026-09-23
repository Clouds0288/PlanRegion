# 配电网可规划域

在给定网架和升级预算下，构造**所有合法建设方案的连续可行调度域并集**。正式入口是 `main.py`，主线采用完整 MP2 / MP1、可选择的剩余域搜索与必要的 SP 认证。每个方案维护连续内外域，不同方案取并集；独立 AC 校核在构域后进行。

## 最新结果

| 数据集 | 结果入口 | 保留内容 |
|---|---|---|
| Case33，8 / 16 / 32 候选，预算 0 / 1 / 2 / 无限 | [结果与耗时](results/case33bw/latest/report.md) · [HTML 总览](results/case33bw/latest/report.html) | 12 组 auto 与 2 组模式对照、14 份完整回放、原始事件、独立 AC 审核、运行环境及源文件指纹 |
| 江口二维，线路升级、无限预算、1000 kVA | [实验报告](tests/jiangkou_2d/REPORT.md) | 连续域、断点、AC 样本、图表及独立实验源码快照 |

Case33 的 14 组中 **11 组完成全局覆盖认证**，其中 auto 为 9 / 12；三个预算 2 的 auto 结果仍为未确定或达到时限。保存的 **3453 个内域顶点全部通过独立 SOCP 与严格 AC 检查**。页面中未完成组的区域误差描述当前内域，不代表最终收敛精度。全部 **47829 帧**算法事件保留。

江口实验在 3600 秒构域时限后保留部分结果；4047 个内域网格点及 21 个内域顶点通过同口径 AC 检查。该实验另列严格载流阈值检查，其 AC 参考边界对应一个具体建设方案，不能当作全局 AC 并集的精确边界。

结果目录只保留上述最新批次。旧批次与过时报告已移出项目，不再混入当前结果入口。

## 项目结构

| 文件或目录 | 职责 |
|---|---|
| `Network/` | 网架、原始数据、负荷、升级选项与费用 |
| `main.py` | 顶部参数、MP2 / MP1 / SP 协作、AC 校核及结果输出 |
| `model.py` | LP / SOCP 方程、MP1 / MP2 查询、SP、联合割及剩余域模型 |
| `region.py` | 构域操作、下一候选点、裁剪、内域并集覆盖与全局外包络 |
| `vertify.py` | 独立 AC 潮流、方案校核、单调网格分类与未知状态处理 |
| `plot.py` | 结果采样、并集体积、FR/MR、存取、表面网格和最终 HTML；历史回放工具仍可独立使用 |
| `live_view.html` | 历史实时显示和离线回放模板，新主流程不使用 |
| `tests/` | 回归测试、基准实验、独立数值和几何审核 |
| `tests/jiangkou_2d/` | 江口二维实验；`source/` 为冻结的实验副本，正式入口不导入它 |
| `docs/` | 当前连续构域说明、统一模型与固定方案 SOCP 推导 |
| `results/case33bw/latest/` | 最新 Case33 报告、完整回放与审核数据 |

新算法结果中的 `inner` 仅保存各方案的 `choice`、`cost` 和 kW 顶点，`outer` 保存全局外包络。
不再重复保存归一化证书、完整割池、逐查询历史或三角面。保留状态、最大负荷及上界、覆盖界、SP / 割次数和构域总时间。
结果配置保存求解线程数、关键数值容差、预算、时限和 AC 网格精度；AC 时间单列。
`plot.py` 根据原始结果按需计算体积、体积差、比较指标和展示网格，派生值不重复写入结果文件。
算法返回 NumPy 顶点和具名走廊方案（走廊 ID → 型号 ID 或 `None`）；保存时统一转换为 JSON 数据。`region.py`、`vertify.py` 不依赖展示模块。
与 pandapower 的六工况潮流交叉核验位于 `tests.reference.validate_power_flow`，不进入正式 AC 模块。

## 安装与查看

本次记录的环境使用 Python 3.13。安装依赖并准备有效的 Gurobi 许可证：

```text
python -m pip install -r requirements.txt
```

完整原始历史使用 **Git LFS** 保存；报告、图片和可独立播放的 HTML 留在普通 Git 中。获取原始数据时运行：

```text
git lfs install
git lfs pull
```

GitHub 可直接阅读 Markdown 报告。交互 HTML 需下载后用浏览器打开，或在项目根目录启动本地服务：

```text
python -m http.server 8000 --bind 127.0.0.1
```

随后打开 [Case33 HTML 总览](http://127.0.0.1:8000/results/case33bw/latest/report.html)，点击每一行的“回放”。查看已有 HTML 不调用求解器。

回放包含播放 / 暂停、速度切换、前后单步和全历史时间轴，可查看每一步的内域、外域、候选点及割。大回放通过共享重复几何对象减小页面体积，未删帧；原始 `replay.json` 与 `events.jsonl` 保留原格式。32 候选、预算 2 的 16312 帧页面约 8.15 MB。

## 运行主流程

修改 `main.py` 顶部参数，然后运行：

```text
python main.py
```

- `NETWORK='case33'`，`CANDIDATE_COUNT=16`；也支持四节点网架及 Case33 的 4 / 8 / 32 候选。
- `BUDGETS=None` 使用网架默认预算；每档预算独立计算。
- `SOLVER_THREADS` 默认取 `model.DEFAULT_SOLVER_THREADS=20`，控制 Gurobi、SP 支持的并行后端及 AC 后备求解器的线程上限；模块独立调用也使用同一默认值。主流程只有一个 Python 进程；小型矩阵计算保持单线程。更多线程不保证更快，可直接调整该参数。
- `RESIDUAL_MODE='auto'` 对有限预算使用 light，无限预算使用 physical；可以直接设置任一模式。
- `CASE_TIME_LIMIT=300` 是每种方法、每档预算的构域时限，混合方法的线性与 SOCP 阶段共享时限。
- `REGION_TAU=0.002` 决定 SOCP 径向精度；`DIVISIONS=8` 只决定最终 AC 校核网格。
- `RECOMPUTE=False` 从 `OUTPUT` 读取已有 `result.npz`，不重新求解。
- `SHOW_UI=True` 在计算结束后打开 `region_comparison.html`；页面包含连续域和 AC 网格比较，不需要后台服务。设为 False 时只保存数据并显示终端汇总。
- `OUTPUT=None` 使用 `results/<网架名>/planning_<候选数>`，也可指定结果目录。

主入口不再接受命令行覆盖参数，不再生成实时事件或完整回放，不再跨预算继承证书。历史 HTML 仍可直接打开。
程序只在参数入口检查输入；数值库线程配置失败直接报错，不重试或记录额外诊断状态。

## 电力走廊与型号

所有线路统一为 `network.corridors` 中的 `Corridor`：

| 字段 | 含义 |
|---|---|
| `id`、`endpoints` | 走廊编号及真实节点端点 |
| `existing_type` | 原有线路型号；`None` 表示原来没有线路 |
| `initial_active` | 原始投入状态；与规划决策相互独立 |
| `must_use` | `True` 必须使用，`False` 可使用也可断开 |
| `types` | 可选型号及各自的阻抗、容量、增量投资 |

型号参数统一用 `TypeParameters`，不再使用 `line_options`、`optional` 或旧的整数方案数组。
`initial_active` 只描述初始网架，不固定优化结果。已有常开线路可同时满足 `existing_type` 非空、`initial_active=False`；空走廊则为 `existing_type=None`。
Case33 的 32 条在运走廊均必须使用，其中 8/16 条允许升级；FourBus 的 5 条走廊均允许开断，其中初始投入 3 条。

```python
from Network.four_bus_five_corridor import FourBus
from model import PlanningEquations, PlanningModel

network = FourBus()
print(network.corridors[0])
print(network.initial_plan)
# {'01': 'L', '12': 'L', '13': 'L', '02': None, '23': None}

equations = PlanningEquations(network, 'linear')
problem = PlanningModel(equations, budget=20000, threads=1)
with problem.model:
    answer = problem.solve()
    plan = equations.choice(answer['x'])
    print(plan)
    print(problem.x['01', 'L'].X)
    print(plan['01'] is not None)
    print(plan['01'] == 'L')
```

`problem.x[走廊ID, 型号ID]` 直接对应选型二进制变量，不再包含方向。`problem.x.sum(走廊ID, '*')` 对该走廊的全部型号求和：必用走廊等于 1，可断开走廊不超过 1。必用单型号直接取常数 1，不在 `x` 中创建对应键；完整方案由 `equations.choice(answer['x'])` 读取。
`problem.x_vector` 是同一批变量按 `equations.variable_keys` 排列的向量视图，供潮流矩阵、联合割、批量固定选型和热启动使用，不增加变量。结果 `answer['x']` 保存相同顺序的数值向量。
潮流使用固定参考方向，允许 P/Q 为负；两端送出功率均受容量约束。可重构网络用每条走廊一个连续连通流和树的边数约束保证径向结构。FourBus 共 15 个选型二进制变量、5 个连通流；Case33-8/16 分别为 16/32 个选型变量，不需连通流。
需要固定完整方案时传入 `fixed_plan=network.initial_plan`；`network.design(plan)` 生成独立潮流所用的固定树。
对外完整方案使用走廊 ID 到型号 ID 或 `None` 的映射，不依赖数值向量的位置。

## 连续域与质量控制

1. 完整 MP2 在预算内最大化总负荷，MP1 在满足该总负荷要求时寻找最低投资方案。原约束复核通过的主问题运行证书直接复用，仅未获证时补充 SP。
2. 普通候选先检查所有合法方案的认证内域。已覆盖则跳过重复 SP；其余点由 SP 认证，或生成有效割来缩小外域。
3. 每个方案的可行点构成该方案的连续凸内域；跨方案只取并集。全局剩余搜索持续寻找并集尚未覆盖的部分。
4. 顶点分别被不同方案覆盖时，内部仍可能有空隙。算法为全局见证补充必要的同方案支撑证书，只有全局覆盖上界满足容差才标记完成。
5. 独立 AC 校核比较实际保留的内域并集与 AC 参考域。区域误差为对称差体积 / 两域并集体积；多算率以计算域为分母，漏算率以 AC 域为分母。网格估计不是连续体积误差的严格上界，未确定状态单独保留。

计时按“方法 × 预算”记录，总时间包含建模、求解与几何计算，AC 校核时间单列；准备、采样结果整理和页面导出不计入构域时间。

详细说明：[连续构域与覆盖证书](docs/continuous_region.md)、[统一规划模型与联合割](docs/compact_planning.md)、[固定方案 SOCP 数学推导](docs/socp_model.md)。

## 算例与复现

Case33 的 A–D、A–H、A–P、A–AF 为嵌套候选集，32 候选覆盖全部在运线路。节点 18、25、33 的负荷独立变化，其余负荷固定。升级为并联一回、R/X 减半，五条常开联络线保持断开。费用为相对单位：A 为 2，其余 31 项各为 1，全部升级共 33；原有 16 项顺序与费用不变。

四节点五走廊基础算例保留在 `Network/four_bus_five_corridor.py`。江口二维实验的模型口径、运行方式与限制见其独立报告。

重新跑基准时指定新的输出目录，避免覆盖已审核结果：

```text
python -m tests.benchmark_physical_search --counts 8 16 32 --budgets 0 1 2 inf --variants auto --limit 300 --output results/case33bw/rerun
python -m tests.benchmark_physical_search --counts 32 --budgets 1 --variants direct_all --limit 300 --output results/case33bw/rerun
python -m tests.benchmark_physical_search --counts 8 --budgets inf --variants direct_mp --limit 300 --output results/case33bw/rerun
python -m tests.benchmark_ac_search --counts 8 16 32 --divisions 32 --output results/case33bw/rerun
python -m unittest discover -s tests -p "test_*.py" -v
```

新基准直接调用正式主流程，使用单线程独立计时，只保存最终数值及审核结果。`direct_mp`、`direct_all` 分别强制轻量和物理剩余域；`baseline` 是仅供对照的冻结旧算法。两个长时间测试默认跳过。

基准实验（而非主入口）的 `protocol.json` 仍保存运行时源码指纹，用于独立复现实验。已有历史批次的回放与记录保留。`final_checks.json` 及 `ac_validation/` 保存审核明细。新试跑默认不进入版本控制，已审核的最新目录与江口实验结果随仓库保存。
