# PlanRegion：有限预算下的可规划域

使用统一候选网架、线路型号和潮流模型，通过方向 MP2 初始化、最大负荷候选点和 SP 联合割构造可规划域，由全局剩余域模型检查覆盖。默认只运行 SOCP，完成构域后独立 AC 校核；LP 和混合模型接口仍保留。

当前默认是五节点六候选道路案例：逐步勘察并更新确认域、乐观域和道路边际信息价值。正式计算不依赖历史实验目录或旧结果。

## 安装

Python 环境需要有效的 Gurobi 许可证。

```text
python -m pip install -r requirements.txt
```

测试依赖单独安装：

```text
python -m pip install -r requirements-dev.txt
```

## 运行

修改 `main.py` 顶部参数后运行：

```text
python main.py
```

| NETWORK | 案例 | 入口行为 |
|---|---|---|
| concept5 | 五节点、4 条既有线、6 条候选道路 | 二维规划域和完整逐轮勘察；当前默认 |
| fourbus | 四节点五走廊、多型号线路 | 三维连续构域及 AC 校核 |
| case33 | Case33，含 5 条联络线 | 可配置升级数量、建设预算及重构 |
| jiangkou | 江口候选网架 | 三维连续构域及 AC 校核 |

五节点采用已验证的概念参数：建设预算 4，p₁/p₂ 为节点 1、2 的负荷，独立可用先验 0.8，单次勘察费 1。`SOLVER_THREADS`、`CASE_TIME_LIMIT`、`OUTPUT`、`RECOMPUTE`、`SHOW_UI`、`STEP_BY_STEP` 由主入口设置；三维网架还支持 `BUDGETS`、`REGION_TAU`、`DIVISIONS`。

默认打开本地调试窗口，并在步骤边界等待。点击“执行下一步”推进一次算法步骤，“连续运行”自动执行，“暂停计算”在下一边界停下。时间轴的上一帧、下一帧只回看记录，不改变计算。人工等待不计入构域时限。`STEP_BY_STEP=False` 自动运行；`SHOW_UI=False` 关闭窗口但继续完整记录。

求解超时、数值精度不足、未收敛或无法取得有效证书时直接报错并停止，不自动重求、切换模型或修补答案。已完成步骤及错误仍写入日志与回放；未完成运行不发布成功结果。

要复现单线程五节点计算及其图件：

```text
python survey.py --threads 1
python plot.py --results results/concept5/<本次时间戳>
```

每次计算由 `main.run` 建立独立时间戳目录。五节点的 `results.json` 包含参数、逐轮域、候选评分、最终方案和核查结果；`steps.jsonl` 在计算期间逐步保存选点、求解返回、割和内外域变化；`live_view.html` 是同一记录的离线回放。绘图另存 PDF、SVG 和 PNG。预期勘察顺序为 A+ → B− → E+ → D− → C+，F 停止勘察。更多定义见 [勘察模型与结果](docs/survey.md)。

三维案例的数值结果保存为 `result.npz`，过程使用相同的 `steps.jsonl` 与 `live_view.html`。中断后可用 `RunMonitor.load_recording` 读取已经完整写入的步骤，也可重新导出回放；这不恢复求解器。比较图仍可由数值结果重新生成。`RECOMPUTE=False` 读取 `OUTPUT` 指定的运行；未指定时读取该网架最新时间戳目录。

## 项目结构

| 文件或目录 | 职责 |
|---|---|
| Network/ | 网架、线路型号、负荷、道路参数和原始数据 |
| main.py | 运行配置、对偶构域、统一记录与保存入口 |
| monitor.py | 步骤记录、暂停/继续、时钟与本地查看服务 |
| model.py | 主问题、物理方程、对偶 SP、联合割与剩余域搜索 |
| region.py | 几何、单方案内外域与并集覆盖 |
| vertify.py | 独立 AC 潮流及可行性校核 |
| survey.py | 按需求域、道路信息价值、逐轮勘察与核查 |
| plot.py、live_view.html | 勘察图、三维结果、实时查看和历史回放；勘察绘图统一调用 render_survey |
| tests/ | 核心回归及独立参考模型；数据夹具只用于测试 |
| docs/ | 当前模型、数学符号和迁移记录 |
| results/ | 已审核结果与本地运行输出 |

主循环全部在 `main.build_continuous_region`，注释按 **1 初始化 → 2 选点 → 3 全局搜索 → 4 SP 更新 → 5 汇总** 编号。二维初始化解 3 次 MP2、三维解 4 次；不再默认解 MP1。`region.py` 只处理几何，`model.py` 只完成单次优化。普通候选按实际 kW 总量排序；内部见证的同网架支撑和最终覆盖证明仍保留。

## 验证

本次主线重构、五步调试入口和 SOCP 验证结果见 [主线验证记录](docs/main_flow_validation.md)。

```text
python -m unittest tests.test_notation -v
python -m pytest tests/test_socp_main.py tests/test_main_flow.py tests/test_solver_results.py -q -k "not linear"
python -m pytest tests -q
```

默认测试包含五节点完整重算。较慢的全规模构域测试由环境变量 `PLANREGION_SLOW_TESTS=1` 启用。浏览器几何回归需要 Node.js。

数学符号、单位、数组顺序和结果字段统一遵守 [notation.md](docs/notation.md)。进一步说明见 [连续构域](docs/continuous_region.md)、[联合割模型](docs/compact_planning.md)、[固定方案 SOCP](docs/socp_model.md)。

## 已审核结果与清理记录

- [五节点对偶复现报告](results/concept5_dual/2026-09-25_run2/report.md)：原始复现用时约 30.03 秒，26 个道路集合、166 条共享割，与枚举结果在数值误差内一致；尚未证明速度优于枚举。
- [Case33 历史报告](results/case33bw/latest/report.md)：保留已审核原始数据。大体积回放由 Git LFS 管理，获取时运行 `git lfs pull`。
- [项目清理记录](docs/cleanup.md)：临时实验和历史试跑已归档到项目外。历史符号登记完整保留，当前代码只引用活动接口。
