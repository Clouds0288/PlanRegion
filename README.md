# PlanRegion

构建配电网二维／三维功率可行域（正负功率、符号分区），独立 AC/SOCP 扫描校验，原生前端实时显示与回放。本分支 method-RB 为方法 RB：径向锥夹逼（R）做全局搜索与外界，逐网架支撑查询（B）认证各网架的内多面体。

## 环境

Python 依赖见 `requirements.txt`，测试另装 `requirements-dev.txt`；需可用的 Gurobi 许可证。本机在 CMD 中执行 `conda activate methods` 后可用 `python`。若未配置 conda 命令，直接使用：

```cmd
"C:\Users\10856\miniconda3\envs\methods\python.exe" -X utf8 main.py --help
```

以下命令在项目根目录运行，其中 `python` 可替换成上述完整路径。

## 算法

功率空间按符号分成 2^d 个分区（正为负荷、负为光伏），每个分区一个子进程并行，坐标为幅值 u=|p|，归一化为 xi=u/bounds∈[0,1]^d。

1. R：方向用单纯锥剖分。每个锥的内域 T=conv(0, v_1..v_d)，v_i 是内域网架沿生成方向的 OBBT 紧化射线顶点；x 自由的锥 MISOCP 给出上界 μ̄，外界 O=μ̄·T。每次细分体积缺口 (μ̄^d-1)·vol(T) 最大的锥。
2. A 阶段：R 先做到放宽目标 ε_A=0.15（最多用分区时限的 25%），发现网架。
3. B：对出现过的每个网架 x，以径向阶段的认证点（射线端点、零接入点、审计过的现任解）为 V_x，内多面体 P_x=conv(V_x)。P_x 的每个非分区边界面用固定 x 的支撑 SOCP 的可靠上界认证，上界同时裁剪该网架的支撑外界 O_x；审计过的越界解补进 V_x。网架在全部面认证或 vol(O_x)/vol(P_x)-1<=ε/2 时停止。
4. A+：续跑 R 到 ε=d·tau。内域 I=I_R ∪ ∪P_x，vol(O_R∩盒)-vol(I)<=ε·vol(I) 即该分区认证。

P_x 由认证点张成，RB 的内域是认证内域；外界 O_R 由 MISOCP 上界给出，是有效外界。超时保留已有证据，只有全部分区获证才报告 certified。符号与字段约定见 [docs/notation.md](docs/notation.md)。

## 计算

| 案例 | 二维节点 | 三维节点 | 预算 |
|---|---|---|---|
| FourBus | 1、2 | 1、2、3 | 20000 元 |
| Case33 | 18、25 | 18、25、30 | 7 次开合 |

正功率为负荷（PF=0.95），负功率为光伏（PF=1）。Case33 保持原线路参数与背景负荷，仅七条开关可变。主入口、AC 与 SOCP 校验共用 `Network/case33bw.py` 的 `CURRENT_LIMIT=200` A。

```cmd
python -X utf8 main.py --case case33 --dimension 2
python -X utf8 main.py --case case33 --dimension 3 --output results/mainline/new_run
```

构域默认共享 300 秒总时限（`main.py` 的 `CASE_TIME_LIMIT`，`--seconds` 可改），`WORKERS=16` 个分区进程并行；事后扫描另计时。`--no-ui` 不开实时窗口，`--no-scan` 只构域。扫描格数默认二维 160×160、三维 80³（`--divisions`）。记录写入 `<输出目录>/mode_1/<案例>_<节点>.json.gz`，逐格对比写入旁边的 `_comparison/`。

## 前端与回放

计算时原生窗口实时显示，结束后同一窗口可回放；回放不重新求解。

```cmd
python -X utf8 main.py --case case33 --dimension 3 --replay
python -X utf8 monitor.py results/mainline/mode_1/case33_18_25_30.json.gz
python -X utf8 monitor.py A.json.gz --compare B.json.gz
```

窗口分三部分：总图显示外包络 O_R 与内域（二维为并集，三维为各锥远端面片），网架面板显示各网架的 O_x（灰）、P_x（绿）、当前支撑面与支撑点，校验面板显示最终 AC/SOCP 扫描对比。支持暂停、单步、逐割跳转、拖动进度、网架翻页和三维旋转；计算中的暂停与单步同时作用于各分区子进程。

## 已保存结果

结果、参数、来源与文件哈希见 [results/manifest.json](results/manifest.json)。

| 保存内容 | 路径 |
|---|---|
| Case33 二维 RB 主线（含 AC/SOCP 逐格对比） | `results/mainline/mode_1/case33_18_25.json.gz`、`case33_18_25_comparison/` |
| Case33 三维 RB 主线 | `results/mainline/mode_1/case33_18_25_30.json.gz` |
| R、H、RB、RCUT、RCUT2 的方法对照（Case33 二维、三维，300 s） | `results/methods/`，说明见其 README |
| Case33 二维 160²（及其扩界）与三维 80³ 配对 AC/SOCP 参考扫描 | `results/scan/case33bw/` |

主线运行（300 s，16 进程，内域相对 AC 扫描）：

| 案例 | 认证分区 | 最慢获证分区 | 内域 FR / MR | 外域 FR / MR |
|---|---|---|---|---|
| Case33 二维 (18,25) | 4/4 | 25.5 s | 0.809% / 0.147% | 1.348% / 0% |
| Case33 三维 (18,25,30) | 3/8（其余 5 个到时限，间隙 1.75–3.10%） | 146.1 s | 2.580% / 0.172% | 4.202% / 0% |

RB 的内域是认证内域，相对 AC 的 FR 来自 OBBT 紧化模型与 AC 的差（SOCP 松弛残余）。三维的 SOCP 参考在现行严格数值门槛下补不全（80³ 网格中 5 格的 SOCP 求解 MaxVio 超过 1e-8，扫描按设计停止），所以三维记录用 `--no-scan` 构域，校验只对该网格完整的 AC 标签发布，网格外的外域格不计入。

旧主线（逐网架顺序构域 + 完整物理查漏）、支持面认证实验、径向夹逼实验的代码与结果保存在 tag `mainline-sequential-v1` 和 `results-methods-v1` 中，需用对应 tag 的前端回放。

## 方法分支

| 分支 / tag | 内容 |
|---|---|
| `main`、`method-RCUT` | 方法 RCUT 主线：R + 主线割平面，内域为外侧估计 |
| `method-RB` | 方法 RB 主线（本分支）：R + 逐网架支撑查询，内域为认证内域 |
| `method-<方法>-v1` | 方法对照实验中各方法的代码版本 |
| `results-methods-v1` | 方法对照归档 |
| `mainline-sequential-v1` | 旧主线 |

## 扫描复用

`vertify.py` 统一管理 AC 与 SOCP。缓存位置为 `results/scan/<网络>/<节点>/<物理身份>/region_<格架摘要>.npz`；相同配置和坐标复用已有点，扩界后只补缺失点。改变电流限额、功率因数、预算、节点顺序或模式会隔离缓存。主线的扫描框为 SOCP 全局界与结果外界范围的并。

```cmd
python -X utf8 vertify.py results/mainline/mode_1/case33_18_25.json.gz --workers 20
```

`summary.json` 中 result_ac/result_socp/socp_ac 各含遗漏率（MR）和多余率（FR）：遗漏率分母为参考可行点数，多余率分母为计算域点数；未决标签不作为不可行。扫描校验不会改写构域认证状态。

## 代码与检查

| 文件 | 职责 |
|---|---|
| `main.py` | 参数、计算与回放入口 |
| `model.py` | 物理模型、OBBT 紧化、SP 与联合割、锥 MISOCP、固定方案支撑查询 |
| `region.py` | 径向锥夹逼、逐网架支撑查询、分区并行 |
| `monitor.py`、`plot.py` | 原生前端、回放与三维并集测度 |
| `vertify.py` | 独立 AC/SOCP 参考、缓存、结果比较 |
| `Network/` | FourBus、当前 Case33 与原始数据 |
| `tests/` | 物理、构域、扫描与回放回归 |

```cmd
python -m unittest tests.test_notation -v
python -m pytest -q
```
