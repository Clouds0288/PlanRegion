# 职责整理与四候选复核

2026-09-21。按 `main / model / region / vertify / plot` 分工完成整理，删除独立的 `continuous.py` 和 `progress.py`。根目录 `live_view.html` 模板保持不变。

| 模块 | 当前职责与主要接口 |
|---|---|
| `main.py` | `joint_benders`、`ContinuousRegion`、`build_continuous_region` 组织 MP1/MP2/SP、阶段和停止判断；`BenchmarkResult` 存取结果；`run` 计时与调度 |
| `model.py` | `PlanningEquations`、`PlanningModel`、`PlanningSP`；`RemainingRegionModel` 构建全局未覆盖区域 MILP |
| `region.py` | `RegionState` 管理方案内外域、有效割裁剪、下一候选点、覆盖半空间、外包络及并集体积 |
| `vertify.py` | `ACPowerFlow`、`ac_planning_query`、`validate_ac_region` 进行最终独立 AC 校核；派生 FR/MR 和比较标签 |
| `plot.py` | `RunMonitor`、终端摘要、绘图、本地网页服务、事件与 HTML 回放导出 |

`region.py` 不依赖求解器、绘图或主流程；其余模块不反向导入 `main.py`。主流程不再直接添加优化器变量或约束。旧的 LP/SOCP 网格构域分支已删除，网格仅在连续域完成后进行采样和 AC 对照。混合方法将有效 LP 割传给 SOCP，并重新认证 SOCP 内域。

规划方程、SP、独立 AC 实现、回放记录器及 JSON/切面工具经语法树比对保持不变。几何层额外修复了四节点算例薄片导致的 Qhull 精度异常：原凸包计算失败时采用可逆坐标缩放，半空间和体积还原到原坐标；不移动原始顶点、不随机抖动、不把小体积直接丢弃。已用该算例的真实退化输入和已知体积薄盒验证。

## case33 四候选结果

条件与整理前一致：A–D 四候选，节点 18/25/33，预算 0/1/2/无限，评价箱 `[400, 4670, 630]` kW，SOCP `tau=0.002`，归一化几何容差 `1e-8`，AC 网格 `8³`。

**12 组结果全部 certified。** 与保留的整理前结果逐项比较，最大总负荷、最大点和线路选型、全部内外多面体、方案证书、联合割、覆盖上界、体积、调用计数与 AC 标签完全一致；最大总负荷差为 **0 kW**。未知薄层与 FR/MR 区间也保持一致。

| 方法 | 预算 | 最大总负荷 kW | 本次求解 s |
|---|---:|---:|---:|
| LP | 0 | 3415.844896 | 0.117 |
| LP | 1 | 4064.250526 | 0.494 |
| LP | 2 | 4180.834017 | 0.590 |
| LP | 无限 | 5131.106376 | 0.118 |
| SOCP | 0 | 3163.172217 | 0.665 |
| SOCP | 1 | 3776.747688 | 5.574 |
| SOCP | 2 | 3895.534526 | 12.526 |
| SOCP | 无限 | 4796.797666 | 9.217 |
| LP→SOCP | 0 | 3163.172230 | 0.690 |
| LP→SOCP | 1 | 3776.747678 | 7.284 |
| LP→SOCP | 2 | 3895.534491 | 13.934 |
| LP→SOCP | 无限 | 4796.797678 | 14.047 |

LP 合计 **1.319 s**，SOCP **27.982 s**，混合方法 **35.955 s**；独立 AC 校核 **3.081 s**。包含公共准备等开销、导出前墙钟时间 **68.655 s**；首次 HTML 导出 **3.091 s**。独立审计另耗时 **7.152 s**。本次使用 `--no-ui` 计算并完整记录，整理前运行启用了实时页面；这些实测值不能单独作为重构带来性能提升的证据。

## 验证与回放

- 常规回归：39 项中 37 项通过，2 项完整规模慢测试未运行，总耗时 25.048 s。覆盖连续几何、MP1 自由分配、未发现方案搜索、未知保留、混合阶段重认证、AC、四节点基础预算、结果存取、网页服务和回放。
- 独立审计：全部 16 个四候选组合、224 个 LP/SOCP 支持查询、1565 个内域顶点、48 条割的全局有效性、8192 个 AC 方案—负荷点检查均通过。
- 整理前后均有 **14142 条事件**。逐帧比较中，除实测时间、源码指纹、UI 开关和独立审计元数据外，全部事件内容一致。
- 原结果目录 `planning_4` 的文件 SHA-256 在整理后逐个复核未变。回放格式保持兼容，旧历史可直接由 `plot.RunMonitor` 载入。
- 浏览器验证了首帧、单步、连续播放、方法与预算切换及最终连续域，控制台无错误。离线 HTML 内嵌的解压数据与 `replay.json` 完全一致，保留 14142 帧和 12 组结果。
- 四节点全部预算及四／八／十六候选全部预算的连续构域列为慢测试，可设置 `PLANREGION_SLOW_TESTS=1` 后运行。本次完成的全域正式认证范围为 case33 四候选的全部四档预算。

结果目录：[`planning_4_refactored`](../results/case33bw/planning_4_refactored/live_view.html)。其中 `validation.json` 保存独立审计结果，`refactoring.json` 保存逐项对比和整理前文件指纹。`live_view.html` 内保留完整过程与最终求解时间。

```text
python main.py --load --ui --output results/case33bw/planning_4_refactored
python -m unittest tests.test_continuous tests.test_progress tests.test_case33 tests.test_four_bus tests.test_vertify tests.test_scale -v
python -m tests.audit_continuous results/case33bw/planning_4_refactored
```
