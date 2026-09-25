# PlanRegion：有限预算下的可规划域

使用统一候选网架、线路型号和潮流模型，通过 MP2 → MP1 → 对偶 SP 联合割构造可规划域，并以剩余域模型检查覆盖。支持 LP、SOCP、混合计算和独立 AC 点校核。

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

五节点采用已验证的概念参数：建设预算 4，p₁/p₂ 为节点 1、2 的负荷，独立可用先验 0.8，单次勘察费 1。`SOLVER_THREADS`、`CASE_TIME_LIMIT`、`OUTPUT` 可从主入口设置；其他三维网架还支持 `BUDGETS`、`REGION_TAU`、`DIVISIONS`、`RECOMPUTE`、`SHOW_UI`。

要复现单线程五节点计算及其图件：

```text
python survey.py --threads 1
python plot.py --results results/concept5/<本次时间戳>
```

每次五节点计算在时间戳目录中只保存一个 `results.json`，包含参数、逐轮域、候选评分、最终方案和核查结果。绘图另存 PDF、SVG 和 PNG。预期勘察顺序为 A+ → B− → E+ → D− → C+，F 停止勘察。更多定义见 [勘察模型与结果](docs/survey.md)。

三维案例的数值结果保存为 `result.npz`；完整过程回放集中在可离线打开的 `live_view.html`，不再重复保存 `events.jsonl` 和 `replay.json`。比较图仍可由数值结果重新生成。

## 项目结构

| 文件或目录 | 职责 |
|---|---|
| Network/ | 网架、线路型号、负荷、道路参数和原始数据 |
| main.py | 运行配置、对偶构域流程和案例选择 |
| model.py | 主问题、物理方程、对偶 SP、联合割与剩余域搜索 |
| region.py | 几何、单方案内外域与并集覆盖 |
| vertify.py | 独立 AC 潮流及可行性校核 |
| survey.py | 按需求域、道路信息价值、逐轮勘察和结果保存 |
| plot.py、live_view.html | 勘察图、三维结果、实时查看和历史回放；勘察绘图统一调用 render_survey |
| tests/ | 核心回归及独立参考模型；数据夹具只用于测试 |
| docs/ | 当前模型、数学符号和迁移记录 |
| results/ | 已审核结果与本地运行输出 |

## 验证

```text
python -m unittest tests.test_notation -v
python -m pytest tests -q
```

默认测试包含五节点完整重算。较慢的全规模构域测试由环境变量 `PLANREGION_SLOW_TESTS=1` 启用。浏览器几何回归需要 Node.js。

数学符号、单位、数组顺序和结果字段统一遵守 [notation.md](docs/notation.md)。进一步说明见 [连续构域](docs/continuous_region.md)、[联合割模型](docs/compact_planning.md)、[固定方案 SOCP](docs/socp_model.md)。

## 已审核结果与清理记录

- [五节点对偶复现报告](results/concept5_dual/2026-09-25_run2/report.md)：原始复现用时约 30.03 秒，26 个道路集合、166 条共享割，与枚举结果在数值误差内一致；尚未证明速度优于枚举。
- [Case33 历史报告](results/case33bw/latest/report.md)：保留已审核原始数据。大体积回放由 Git LFS 管理，获取时运行 `git lfs pull`。
- [项目清理记录](docs/cleanup.md)：临时实验和历史试跑已归档到项目外。历史符号登记完整保留，当前代码只引用活动接口。
