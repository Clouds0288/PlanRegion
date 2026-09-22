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
| `main.py` | 主流程、MP2 / MP1 / SP 协作、停止判断、计时和结果存取 |
| `model.py` | LP / SOCP 方程、整数主问题、SP、联合割及两种剩余域模型 |
| `region.py` | 连续域、下一候选点、裁剪、并集、覆盖半空间与几何计算 |
| `vertify.py` | 独立 AC 潮流、方案校核、区域误差及未确定标签 |
| `plot.py` | 绘图、进度、完整事件记录、本地服务与 HTML 生成 |
| `live_view.html` | 实时显示和完整离线回放的共用模板 |
| `tests/` | 回归测试、基准实验、独立数值和几何审核 |
| `tests/jiangkou_2d/` | 江口二维实验；`source/` 为冻结的实验副本，正式入口不导入它 |
| `docs/` | 当前连续构域说明、统一模型与固定方案 SOCP 推导 |
| `results/case33bw/latest/` | 最新 Case33 报告、完整回放与审核数据 |

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

```text
python main.py --network case33 --candidates 16 --residual-mode auto --budgets 0 1 2 inf --ui --output results/case33bw/run
python main.py --network case33 --candidates 32 --residual-mode physical --budgets inf --no-ui --no-hold --output results/case33bw/run32
python main.py --load --ui --output results/case33bw/run
```

直接运行时读取 `main.py` 顶部设置，命令行可覆盖它们。当前 `CANDIDATE_COUNT=16`，支持 4、8、16、32；预算默认 0、1、2、无限。

- `RESIDUAL_MODE='auto'`：有限预算选择 `light`，无限预算选择 `physical`；可用 `--residual-mode` 强制任一模式。
- 轻量剩余域使用有效联合割搜索未覆盖部分；物理剩余域在同一预算、割及并集排除约束上加入完整物理方程。两者共用正式模型与构域流程。
- `CASE_TIME_LIMIT=300` / `--time-limit 300`：每种方法、每档预算的构域总时限。混合方法的两个阶段共享时限；超时保留当前内域和未确定部分。
- `REGION_TAU=0.002` / `--tau`：SOCP 连续域径向精度。`--divisions` 控制独立 AC 校核网格，不决定构域顶点。
- 默认每档预算独立计算；`--reuse-budgets` 可复用同一模型在较低预算下的证书和有效割。

`--no-ui` 关闭实时服务，仍导出离线页面；同时设置 `--no-plots` 则跳过页面导出，但保留事件。`--no-hold` 在保存后退出。`--load` 读取主流程生成的 `result.npz` 和回放；基准实验保存的 JSON 数据集直接使用其 HTML 查看。

若 Windows 数值库线程检测失败，入口会重试并记录线程控制状态。该状态会影响耗时比较口径，模型与认证条件保持一致。

## 连续域与质量控制

1. 完整 MP2 在预算内最大化总负荷，MP1 在满足该总负荷要求时寻找最低投资方案。原约束复核通过的主问题运行证书直接复用，仅未获证时补充 SP。
2. 普通候选先检查所有合法方案的认证内域。已覆盖则跳过重复 SP；其余点由 SP 认证，或生成有效割来缩小外域。
3. 每个方案的可行点构成该方案的连续凸内域；跨方案只取并集。全局剩余搜索持续寻找并集尚未覆盖的部分。
4. 顶点分别被不同方案覆盖时，内部仍可能有空隙。算法为全局见证补充必要的同方案支撑证书，只有全局覆盖上界满足容差才标记完成。
5. 独立 AC 校核比较实际保留的内域并集与 AC 参考域。区域误差为对称差体积 / 两域并集体积；多算率以计算域为分母，漏算率以 AC 域为分母。网格估计不是连续体积误差的严格上界，未确定状态单独保留。

计时按“方法 × 预算”记录，总时间包含建模、求解、几何更新和事件记录；公共准备、AC 校核及页面导出单列。回放速度不改变原始求解时间。

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

基准直接调用正式主流程。`direct_mp`、`direct_all` 分别强制轻量和物理剩余域；`baseline` 是仅供对照的冻结旧算法。两个长时间测试默认跳过。

基准结果的 `protocol.json` 保存运行时源码指纹；后续回放压缩修复另记 `presentation_update`，原始数值与求解时间不变。`final_checks.json` 及 `ac_validation/` 保存审核明细。新试跑默认不进入版本控制，已审核的最新目录与江口实验结果随仓库保存。
