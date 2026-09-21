# 配电网可规划域

正式入口是 **main.py**。当前恢复了 **MP1 / MP2 / SP 连续构域**：每个建设方案形成连续内外多面体，不同方案取并集；整数主问题搜索未覆盖部分并给出全域覆盖证书。网格仅在构域完成后用于独立 AC 校核与 FR/MR 统计。

最终交付页面沿用根目录的 **live_view.html** 模板，可离线回放全部算法事件，查看连续域、最终方案和实际求解时间。并集驱动优化的结果在 `results/case33bw/planning_4_union/`，页面包含原流程、优化后独立求解、跨预算复用三组对照；之前的 `planning_4/` 和 `planning_4_refactored/` 结果保留。

| 文件 | 职责 |
|---|---|
| `Network/` | 网架、负荷、逐线路升级选项与费用 |
| `main.py` | 算法主流程、MP1/MP2/SP 协作、方法与预算阶段、停止判断、计时及结果存取 |
| `model.py` | LP/SOCP 规划方程、整数主问题、SP、联合割和全局残余搜索模型 |
| `region.py` | 连续域状态、下一候选点、裁剪、内外域并集、覆盖半空间和体积；AC 网格几何辅助 |
| `vertify.py` | 独立 AC 潮流与非凸校核、AC 方案搜索及采样、FR/MR 与三态标签 |
| `plot.py` | 绘图、终端进度、事件记录、本地服务、HTML 生成与完整回放导出 |
| `live_view.html` | 唯一的实时监视、过程回放和最终结果页面模板 |
| `tests/` | 独立模型、几何、覆盖证书、割、AC 和回放核验 |

## 运行与回放

使用安装了 `requirements.txt` 且拥有有效 Gurobi 许可证的 Python 环境：

```text
python main.py --network case33 --candidates 4 --divisions 8 --ui
python main.py --network case33 --candidates 4 --divisions 8 --no-ui --no-hold
python main.py --load --ui --output results/case33bw/planning_4
python main.py --load --ui --output results/case33bw/planning_4_refactored
python main.py --load --ui --output results/case33bw/planning_4_union
python main.py --reuse-budgets --output results/case33bw/planning_4_union_incremental
```

直接运行时读取 `main.py` 顶部的常用设置；`CANDIDATE_COUNT` 可直接切换为 4、8、16，命令行 `--candidates` 可覆盖它。默认预算为 0、1、2、无限，`REGION_TAU=0.002`。`--tau` 控制 SOCP 连续域径向精度；`--divisions` 控制 AC 校核网格，不决定连续域的顶点或边界。LP 使用零径向收缩，仅保留数值几何容差。

Windows 若出现 `GetModuleFileNameEx failed`，是数值库线程检测失败，与候选数量无关。入口会重试三次；仍失败时使用数值库当前线程数继续，终端和回放记录提示，结果元数据的 `thread_control` 记录限制是否生效。此时模型、精度和求解认证条件不变，但耗时不能视为受控单线程基准；其他初始化错误和实际求解错误仍正常抛出。

网页支持播放/暂停、速度切换、前后单步、全历史时间轴、方法与预算筛选，以及点击事件还原当时的内域、外域、候选点和割。最终结果表可切换所显示的连续域，侧栏给出最大总负荷对应的节点负荷与建设方案。联合割同时约束负荷 p 和规划变量 x；图中切面固定在生成它的方案下，不能当作整个规划域的统一边界。

`--no-ui` 关闭实时服务，仍记录完整过程并导出离线页面。`--no-plots --no-ui` 跳过 HTML 导出，但保留逐事件日志。`--load` 只读取结果和已有回放，不调用优化器。实时网页仅监听本机；页面及绘图资源均可离线使用。`--no-hold` 在保存后退出；默认实时运行完成后按 Enter 或 Ctrl+C 关闭服务。

## 连续构域

1. 三个 LP 轴向查询给出公共评价箱。
2. MP2 在预算内最大化三个独立节点的总负荷，并用 SP 取得可行下界和全局上界。
3. MP1 固定总负荷要求、允许节点间重新分配，寻找最低投资方案。
4. 对发现的方案维护连续外多面体；普通候选先检查所有合法方案的认证内域，已覆盖则跳过 SP 并保留支撑方案。其余候选由 SP 认证或产生全局割，同方案证书形成连续凸内域。
5. 在含全部整数方案的联合割外域上搜索未覆盖区域，使用返回的方案和负荷见证继续探索。如果顶点已分别覆盖但内部仍有空隙，只补齐支撑该见证所需的至多四点同方案单纯形证书。只有全局上界满足容差才停止。
6. 保存内域并集和由覆盖证书推出的全局外包络，随后开展独立 AC 采样校核。

LP、SOCP、LP→SOCP 各自运行，不借用其他对照方法的成果。默认每档预算独立求解；`--reuse-budgets` 可按递增预算继承同一模型的证书与有效割。case33 四候选实测独立预算模式更快，因此作为默认；复用模式仍保留并单独报告。混合方法包含自己的 LP 阶段，LP 内域不能作为 SOCP 内域；SOCP 只继承有效 LP 割和（开启预算复用时）此前低预算 SOCP 的认证内域。没有证书时保留“未确定”，不把求解失败判为不可行。

覆盖证书、数值容差及外包络的定义见 [连续构域说明](docs/continuous_region.md)；本次速度、SP 去重及质量核对见 [并集优化测试结果](docs/union_optimization_results.md)；基础物理模型和联合割推导见 [模型说明](docs/compact_planning.md)。

## 算例配置

case33 的 A–D、A–H、A–P 是四、八、十六条嵌套候选集。本次正式验证范围是 **四候选 A–D**。节点 18、25、33 独立变化，其余负荷固定。32 条在运支路及五条常开联络线状态不变；候选线路可保持原状或并联一回，R/X 减半。费用是合成相对单位：A 为 2，B/C/D 为 1。

四节点五走廊基础网架仍保留在 `Network/four_bus_five_corridor.py`，既有 01、12、13，候选新建 02、23，保留原 L/M/H 设备、线路长度、造价、0.4 kV 电压、150 kVA 配变和 0.95 功率因数。独立负荷为节点 1、2、3。

## 结果文件与计时

| 文件 | 内容 |
|---|---|
| `live_view.html` | 同一页面模板导出的完整离线回放，内嵌压缩数据和 Plotly |
| `replay.json` | 全部历史增量帧、最终结果与元数据，可重新载入页面 |
| `events.jsonl` | 每次事件立即写入的完整日志，不作条数截断 |
| `result.npz` | 连续内外域、方案证书、割、源文件指纹、实测时间和 AC 校核标签 |
| `validation.json` | 单独运行独立核验后的结果 |
| `optimization.json` | 原流程、优化后独立求解和预算复用的对照，以及逐事件去重证书审核 |

每个“方法 × 预算”有独立总耗时，并列出 MP、SP、几何处理时间；总耗时还包括建模、事件记录等开销。混合方法计入自己的 LP 与 SOCP 两阶段。公共准备、AC 校核、HTML 导出和测试侧审核分别计时。回放速度不影响原始时间。

连续内外域之间的薄层在采样时保留为未知；FR/MR 是公共网格上的 AC 对照统计，不等于连续域精度证明。SOCP 精度由径向收缩参数和全局覆盖证书控制。四候选原结果见 [结果说明](docs/case33_continuous_results.md)，职责整理和前后比较见 [整理复核](docs/refactoring_results.md)。旧 `docs/` 实验和 `planning_16` 页面保留作历史记录，不代表此次连续构域运行。

## 验证

```text
python -m unittest tests.test_continuous tests.test_progress tests.test_case33 tests.test_four_bus tests.test_vertify -v
python -m tests.audit_continuous results/case33bw/planning_4_refactored
python -m unittest tests.test_union -v
python -m unittest tests.test_startup tests.test_progress -v
python -m tests.audit_union results/case33bw/planning_4_union_baseline results/case33bw/planning_4_union results/case33bw/planning_4_union_incremental
```

独立审核在测试侧枚举全部 16 个四线路组合，核对 LP/SOCP 边界、全部内域顶点、代表性割的全局有效性及 AC 校核标签。审核通过后，将结果追加到同一个 `live_view.html` 的最终结果区，保留原始求解时间和事件。正式算法不生成完整方案枚举表、不导入 tests。
