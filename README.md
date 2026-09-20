# 配电网可规划域

唯一运行入口是 **main.ipynb**。四节点与 33 节点共用模型、几何和验证代码。

当前默认 `COMPACT_REVIEW=True`，审核 **33 节点四条候选线路的逐线路紧凑 MILP/MISOCP 与联合割**。每条线路两个型号，共 8 个选型二进制变量；求解器不读取完整方案表。16 方案枚举仅作为独立对照。模型、等价性证明和割推导见 [审核文档](docs/compact_planning.md)，唯一审核记录为 `results/case33bw/compact_review/review.npz`。

`COMPACT_REVIEW=False` 运行下面保留的小规模区域/AC 基准。该枚举基准没有被标成可扩展求解器；全网选型及非枚举完整构域待本轮审核后推进。

| 文件 | 职责 |
|---|---|
| Network/ | 原始网架、基础负荷、独立负荷节点、逐线路选项与费用；小算例另提供枚举对照 |
| model.py | 仅四个模型类：统一方程、紧凑规划求解、连续 LP/SOCP 认证与联合割、独立 AC |
| region.py | 顶点、半空间、裁剪、同方案认证及并集覆盖等几何运算 |
| vertify.py | 连续方案域上的联合割核验、紧凑模型核对汇总、独立 AC/FR/MR 与结果存取 |
| plot.py | 区域比较和固定建设方案的真实切割回放 |
| main.ipynb | 配置、MP/SP 迭代、区域切割、预算及方法循环、AC 验证、计时与输出 |
| tests/reference.py | 独立消元参考方程，仅用于 16 方案核对与联合割审核，不参与正式求解 |

`planning.py`、`replay.py`、`NodePower.ipynb` 和旧 `vertify.ipynb` 已移到 `archive/2026-09-20/`，只保留为历史档案，不参与当前代码导入。

## 四个模型类

| 类 | 职责 |
|---|---|
| `PlanningEquations(network, method, planning=True)` | 一次构建固定系数、物理有效上界及型号联接；`method` 为 `linear` 或 `socp` |
| `PlanningModel(equations, ...)` | 共用方程求最低投资或预算内最大负荷；`cuts_only=True` 时只建立联合割 MP |
| `PlanningSP(equations)` | 固定 `(x,p)` 求连续 LP/SOCP，检查原约束并生成对所有选型有效的联合割 |
| `ACPowerFlow(network)` | 用独立支路递推及完整 AC 电流等式核验，不读取规划方程或乘子 |

同一套 `equations` 可被多个查询复用。`planning=False` 表示固定网架，选型向量为空，联合割自然退化为负荷割。旧 `BranchEquations`、`LinearSP`、`SOCPSP`、完整方案选择 MP 及条件割接口已移除。

SP 的 `feasible` 字段表示原查询通过了原始约束核验。`cut=None` 本身不代表可行：数值边界可能既没有可行证书，也没有可靠分离割。Notebook 明确处理这种未确定状态，生成的内点与原查询分别保存。所有当前 Python 文件和 Notebook 代码单元逐行附有中文注释；第三方原始数据和历史归档不改写。

## 使用

在安装 requirements.txt 的环境中，使用 main.ipynb，从头执行所有单元。

- `COMPACT_REVIEW=True`：本轮使用 `CASE="case33bw"`、`PLANNING=True` 和预算 `(0.,1.,2.,np.inf)`；直接运行 Notebook 得到模型核对表。
- `COMPACT_REVIEW=False`、`CASE="case33bw"`、`PLANNING=True`：原 33 节点四项投资区域基准；预算 `(0., 1., 2., np.inf)`，费用为相对投资单位。
- 同一算例设置 `PLANNING=False`：固定零投资方案，复现第一步的可调度域截面。
- `CASE="four_bus_five_corridor"`、`PLANNING=True`：原四节点选线/选型规划，预算改为 `(20000., 40000., 60000., np.inf)`，费用单位元。
- `RECOMPUTE=True` 明确重新求解；`False` 读取该输出目录中当前格式的结果。没有自动缓存回退或历史格式兼容层。
- `DIVISIONS` 控制事后网格评价，`RADIAL_TOLERANCE` 控制 SOCP 连续内外域停止精度；二者独立。

输出位于 `results/<CASE>/planning/` 或 `results/<CASE>/fixed/`。`result.npz` 是唯一原始结果，HTML 是可再生的展示文件。

## 原小规模区域基准的模型与流程

独立负荷参数为 p；Network 给定其节点和 Q/P 比例，其余负荷固定。规划域定义为

    D(B) = union_{design.cost <= B} D(design).

小规模参考方案由 Network 生成一次。预算循环先筛选允许方案，再用 `PlanningEquations(..., planning=False)` 和统一 `PlanningSP` 构建每个固定方案的负荷域。完整区域基准仍需要有限方案表，用于核对非凸并集；正式紧凑规划查询只读取逐线路型号表，不访问 `designs`。

线性流程直接以几何残块搜索选择需要认证的顶点，调用对应固定方案的 LP SP；有效割只收紧该方案的外域。已经被认证并集覆盖的块可跳过，但任何凸组合的可行性证明仍须来自同一建设方案。原来只为这一基准服务的费用前沿 MP1/MP2 及重复割管理已删除。紧凑模型的最低投资与预算内最大负荷查询继续由 `PlanningModel` 提供。

纯 SOCP 与线性＋SOCP 复用 Notebook 中同一个 `cut_region` 循环；后者以各方案的线性外域初始化。保持 I 为已认证内点的凸包，U 为有效割外域，直到 `(1-tau) U` 被 I 包含。不同建设方案始终取并集。

独立 AC 位于 model.ACPowerFlow，以另行实现的支路递推保留完整电流等式，不读取 LP/SOCP 方程矩阵和对偶乘子。未确定点由显式非凸 AC 模型核验，不能被默认归为不可行。

## 33 节点投资设置

基础网架来自 MATPOWER case33bw，32 条闭合支路及 5 条常开联络线状态不变。节点 18、25、33 独立变化，其他 P/Q 为原始工况。电压限值 0.9–1.1 p.u.，原始 rateA=0 不补造热限，10 MVA 仅为标幺基准。

| 变量 | 并联改造支路 | 相对费用 |
|---|---|---:|
| A | 2–3 | 2 |
| B | 12–13 | 1 |
| C | 23–24 | 1 |
| D | 27–28 | 1 |

建成后按两回相同线路等值，将该支路 R、X 减半。费用是合成实验设置，不是原始算例的工程造价。共 16 种建设组合，最低费用方案就是原网架。

## 数据与计时

- 三种计算域掩码以位数组保存；AC 每个网格点只保存一个最小可行投资，各预算的 AC 掩码按需生成。
- 每个方案只记录一次建设向量和费用。区域记录包含初始域、最终内外域及 SP 查询；相同记录在多预算出现时去重。
- FR、MR、总时间和差集颜色均读取基础记录计算，不另写 summary.json 或重复标签文件。
- 各预算、各构域方法重新建模和求解。混合方法计入自己的线性阶段。AC 只按费用顺序扫描一次，各预算显示处理完该预算全部方案的累计时间；它是共享参考成本，不加到三种方法时间上。
- 公共网格创建、许可证启动和图像导出不计入方法总时间；模型构建、求解、几何处理及网格归属判定计入。
- 网格 FR/MR 是数值估计；网格上 MR=0 不表示连续区域误差严格为零。

此前实验结果仍保留在原目录，历史表格见 docs/case33_results.md 和 docs/socp_benchmark_results.md。旧格式不通过新入口读取，也不将旧计时标成新代码的运行时间。

## 验证

    python -m unittest discover -s tests -v

测试直接读取 main.ipynb 的函数单元，避免维护第二套流程实现。包含四节点全部 216 个方案的独立 AC/QCQP 交叉核验、节点电压重建、有效锥割与内点证书、原线性区域回归、33 节点规划并集及最小投资核验、唯一结果文件往返检查。Notebook 还需在新内核中完整执行。
