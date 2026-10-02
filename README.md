# PlanRegion

构建配电网二维／三维功率可行域（正负功率、符号分区），独立 AC/SOCP 扫描校验，原生前端实时显示与回放。主线为方法 RCUT：径向锥夹逼（R）做全局搜索与外界，主线割平面（CUT）逐网架构域。

## 环境

Python 依赖见 `requirements.txt`，测试另装 `requirements-dev.txt`；需可用的 Gurobi 许可证。本机在 CMD 中执行 `conda activate methods` 后可用 `python`。若未配置 conda 命令，直接使用：

```cmd
"C:\Users\10856\miniconda3\envs\methods\python.exe" -X utf8 main.py --help
```

以下命令在项目根目录运行，其中 `python` 可替换成上述完整路径。

## 算法

功率空间按符号分成 2^d 个分区（正为负荷、负为光伏），每个分区一个子进程并行，坐标为幅值 u=|p|，归一化为 xi=u/bounds∈[0,1]^d。

1. R：方向用单纯锥 K_k 剖分。锥内三角 K^IN_k=conv(0, v_{k,1..d})，v_{k,i} 是内域网架沿生成方向的 OBBT 紧化射线顶点（该网架可行集 R^SOCP_x 的边界点）；网架不固定的锥 MISOCP 给出上界 μ_k，锥外块 K^OUT_k=μ_k·K^IN_k∩盒。每次细分体积缺口 (μ_k^d-1)·vol(K^IN_k) 最大的锥。初始网架由中心射线 MISOCP 选出，其懒惰 OBBT 至多 2 轮，现任仍未紧化时在已紧化网架中取中心射线最远者。
2. A 阶段：R 先做到放宽目标 ε_A=0.15（最多用分区时限的 25%），发现网架。
3. CUT：只割可能改变 **I** 的网架（叶锥的内域网架、覆盖网架与锥 MISOCP 解点网架）。N^CUT_x 从 conv(**K**^OUT)∩盒出发，SP 为其顶点评分，违反量最大的顶点取联合割；连续 3 次割掉的体积都小于 1% 即停滞，N^CUT_x 当作该网架的可行域。每割完一个网架看夹逼间隙：已达 ε，或这个网架只让间隙下降不到 0.1ε，就结束这一轮。
4. A+：续跑 R 到 ε=d·tau，每次锥决策后先割新出现的网架。结果内域 **I**=(**K**^IN ∪ **N**^CUT)∩**K**^OUT，vol(**K**^OUT)-vol(**I**)<=ε·vol(**I**) 即该分区认证。

OBBT 的线程数为 `WORKERS` // 仍在计算的分区数：先结束的分区把线程让给仍在计算的分区。

N^CUT_x 是网架可行域的外近似，所以 RCUT 的结果内域 **I** 是外侧估计，不是认证内域；**K**^OUT 由 MISOCP 上界给出，是有效外界；**K**^IN 只对紧化 SOCP 有证书。各区域的定义与示意图见 docs/notation.md 的“区域与符号”。超时保留已有证据，只有全部分区获证才报告 certified。符号与字段约定见 [docs/notation.md](docs/notation.md)。

## 计算

| 案例 | 二维节点 | 三维节点 | 预算 |
|---|---|---|---|
| FourBus | 1、2 | 1、2、3 | 20000 元 |
| Case33 | 18、25 | 18、25、30 | 7 次开合 |
| Case33Plan | 18、25 | 18、25、30 | 相对建设费 14（全部候选可建） |

正功率为负荷（PF=0.95），负功率为光伏（PF=1）。Case33 保持原线路参数与背景负荷，仅七条开关可变。Case33Plan 是扩展规划算例：S1–S5（7-8、11-12、14-15、28-29、32-33）为可开断的既有线路，基态闭合、开断不计费；C1–C5（原联络线 8-21、9-15、12-22、18-33、25-29，阻抗不变）为基态不建的候选，相对建设费 4、4、4、1、1，共 87 个径向方案。主入口、AC 与 SOCP 校验共用 `Network/case33bw.py` 的 `CURRENT_LIMIT=200` A。

```cmd
python -X utf8 main.py --case case33 --dimension 2
python -X utf8 main.py --case case33 --dimension 3 --output results/mainline/new_run
python -X utf8 main.py --case case33plan --dimension 3 --seconds 1000 --output results/plan/run_1
python -X utf8 main.py --convergence results/plan/run_1/mode_1/case33plan_18_25_30.json.gz results/plan/run_2/mode_1/case33plan_18_25_30.json.gz
```

构域默认共享 300 秒总时限（`main.py` 的 `CASE_TIME_LIMIT`，`--seconds` 可改），`WORKERS=16` 个分区进程并行；事后扫描另计时。`--no-ui` 不开实时窗口，`--no-scan` 只构域。扫描格数默认二维 160×160、三维 80³（`--divisions`）。记录写入 `<输出目录>/mode_1/<案例>_<节点>.json.gz`，逐格对比写入旁边的 `_comparison/`。`--convergence` 把同一算例多次运行的收敛过程（各分区夹逼间隙、内域相对 AC 的 MR/FR 随时间）画成一张图，存于各运行目录的公共上级。

## 前端与回放

计算时原生窗口实时显示，结束后同一窗口可回放；回放不重新求解。

```cmd
python -X utf8 main.py --case case33 --dimension 3 --replay
python -X utf8 monitor.py results/mainline/mode_1/case33_18_25_30.json.gz
python -X utf8 monitor.py A.json.gz --compare B.json.gz
```

窗口左侧是全局总图 A：全部分区的外界 **K**^OUT 与结果内域 **I**（二维为并集，三维为各锥远端面片）在同一张图上，坐标取 K^OUT 的范围并保持不变（回放全程固定；实时运行只在内容明显超出或收缩时调整，可勾选取整个分区盒）。各分区在子进程中并行计算，事件交错到达，因此 A 叠加的是“步骤”：每次求解（中心射线 MISOCP max Σξ、OBBT、p=0 的 SP、射线 max t、反向射线、锥 MISOCP max c·ξ → μ_k 及其远端面、细分、顶点 SP min η、取割、夹逼判据）各成一帧，A 上画最近几步的点与线，当前一步加标签，本帧新建或求界的锥与切割中的 N^CUT_x 描边；下方步骤栏列出前后各步的完整说明，点击跳转。“逐步跟踪”选一个分区后，上一步/下一步、上一割/下一割与播放只走该分区的步骤。右侧标签页 B 是各网架的 N^CUT_x、当前割与当前步骤中属于该网架的点（新网架追加在末页，默认不跟随当前网架翻页），C 是最终的逐格对比：参考可行格为浅色（三维画其表面），遗漏格红色、多余格橙色，可在计算域—AC、计算域—SOCP、SOCP—AC 之间切换。支持暂停、单步、拖动进度、网架翻页和三维旋转；计算中的暂停与单步同时作用于各分区子进程。

## 已保存结果

结果、参数、来源与文件哈希见 [results/manifest.json](results/manifest.json)。

| 保存内容 | 路径 |
|---|---|
| Case33 二维 RCUT 主线（含 AC/SOCP 逐格对比） | `results/mainline/mode_1/case33_18_25.json.gz`、`case33_18_25_comparison/` |
| Case33 三维 RCUT 主线（逐格对比只存 NPZ 与摘要） | `results/mainline/mode_1/case33_18_25_30.json.gz`、`case33_18_25_30_comparison/` |
| RCUT 的方法对照运行（停滞阈值 0.5%、1%、2%，Case33 二维、三维，300 s） | `results/methods/`，说明见其 README |
| Case33 二维 160²（及其扩界）与三维 80³ 配对 AC/SOCP 参考扫描 | `results/scan/case33bw/` |

主线运行（300 s，16 进程，内域相对 AC 扫描）：

| 案例 | 认证分区 | 最慢获证分区 | 内域 FR / MR | 外域 FR / MR |
|---|---|---|---|---|
| Case33 二维 (18,25) | 4/4 | 20.2 s | 0.947% / 0.031% | 1.360% / 0% |
| Case33 三维 (18,25,30) | 7/8（--+ 到时限，间隙 1.82%） | 255.7 s | 4.022% / 0.011% | 5.464% / 0% |

二维、三维均为当前主线（吸收五项改动：按间隙结束 CUT、只割叶锥相关网架、从 conv(**K**^OUT) 出发、限制中心射线懒惰 OBBT、OBBT 线程随剩余分区分配）的运行，扫描网格随外域扩为 161×170 与 142×179×116。

两次都是 `main.py` 默认设置的完整流程（构域、扫描、校验）；扫描框为 SOCP 全局界与外包络范围的并，三维为 97×115×94 格，SOCP 参考只有 1 格未决（AC 不可行，不计入）。这次运行时机器后台负载较高，分区耗时比此前同代码的运行长约 20–40%，三维 --+ 分区到时限时的间隙因此较大（此前同算法 1.77%，内域 FR 4.279%）。内域多余格主要在含光伏反送的分区：二维 -- 123、-+ 35 格；三维 --- 2368、-+- 1493、+-- 1157、--+ 992 格，纯负荷分区 +++、+-+ 几乎没有。RCUT 的 **I** 含 **N**^CUT 外近似，FR 是外侧估计误差（二维 158 个多余格中 123 个来自 SOCP 松弛残余、在 **K**^IN 内）：它完整包含 RB 的认证内域，只多出一层 1–2 格厚的边界壳（见 `results/methods/README.md`）。

R、H、RB、RCUT2 的方法对照运行与汇总报告，以及旧主线（逐网架顺序构域 + 完整物理查漏）、支持面认证实验、径向夹逼实验的代码与结果保存在 tag `results-methods-v1` 和 `mainline-sequential-v1` 中，需用对应 tag 的前端回放。

## 方法分支

| 分支 / tag | 内容 |
|---|---|
| `main`、`method-RCUT` | 方法 RCUT 主线（本分支） |
| `method-RB` | 方法 RB 主线：R + 逐网架支撑查询，内域为认证内域 |
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
| `model.py` | 物理模型、OBBT 紧化、SP 与联合割、锥 MISOCP |
| `region.py` | 选点流程：径向锥夹逼、逐网架割平面、分区并行，以及凸多面体几何与并集测度 |
| `monitor.py` | 过程记录与原生窗口（实时 / 回放） |
| `plot.py` | 绘图基元：配色、二维几何、三维凸域、远端面片、体素表面、割的截线 |
| `vertify.py` | 独立 AC/SOCP 参考、缓存、结果比较 |
| `Network/` | FourBus、Case33、Case33Plan 与原始数据 |
| `tests/` | 物理、构域、扫描与回放回归 |

```cmd
python -m unittest tests.test_notation -v
python -m pytest -q
```
