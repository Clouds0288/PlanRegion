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

1. R：方向用单纯锥剖分。每个锥的内域 T=conv(0, v_1..v_d)，v_i 是内域网架沿生成方向的 OBBT 紧化射线顶点；x 自由的锥 MISOCP 给出上界 μ̄，外界 O=μ̄·T。每次细分体积缺口 (μ̄^d-1)·vol(T) 最大的锥。
2. A 阶段：R 先做到放宽目标 ε_A=0.15（最多用分区时限的 25%），发现网架。
3. CUT：对出现过的每个网架 x，从分区盒出发做主线割平面：SP 为 N_x 顶点评分，违反量最大的顶点取联合割；连续 3 次割掉的体积都小于 1% 即停滞，N_x 当作该网架的可行域。
4. A+：续跑 R 到 ε=d·tau，新出现的网架也先割。内域 I=(I_R ∪ ∪N_x)∩O_R，vol(O_R∩盒)-vol(I)<=ε·vol(I) 即该分区认证。

N_x 是网架可行域的外近似，所以 RCUT 的内域是外侧估计，不是认证内域；外界 O_R 由 MISOCP 上界给出，是有效外界。超时保留已有证据，只有全部分区获证才报告 certified。符号与字段约定见 [docs/notation.md](docs/notation.md)。

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

窗口分三部分：总图显示外包络 O_R 与内域（二维为并集，三维为各锥远端面片），网架面板显示各网架的 N_x、当前割与取割的 SP 顶点，校验面板是最终的逐格对比：参考可行格为浅色（三维画其表面），遗漏格红色、多余格橙色，可在计算域—AC、计算域—SOCP、SOCP—AC 之间切换，标题给出遗漏与多余的格数。支持暂停、单步、逐割跳转、拖动进度、网架翻页和三维旋转；计算中的暂停与单步同时作用于各分区子进程。

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
| Case33 二维 (18,25) | 4/4 | 46.2 s | 0.959% / 0.025% | 1.544% / 0% |
| Case33 三维 (18,25,30) | 7/8（--+ 到时限，间隙 7.33%） | 268.3 s | 4.367% / 0% | 6.938% / 0% |

两次都是 `main.py` 默认设置的完整流程（构域、扫描、校验）；扫描框为 SOCP 全局界与外包络范围的并，三维为 97×115×94 格，SOCP 参考只有 1 格未决（AC 不可行，不计入）。这次运行时机器后台负载较高，分区耗时比此前同代码的运行长约 20–40%，三维 --+ 分区到时限时的间隙因此较大（此前同算法 1.77%，内域 FR 4.279%）。内域多余格主要在含光伏反送的分区：二维 -- 123、-+ 35 格；三维 --- 2368、-+- 1493、+-- 1157、--+ 992 格，纯负荷分区 +++、+-+ 几乎没有。RCUT 的内域含 N_x 外近似，FR 是外侧估计误差：它完整包含 RB 的认证内域，只多出一层 1–2 格厚的边界壳（见 `results/methods/README.md`）。

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
| `region.py` | 径向锥夹逼、逐网架割平面、分区并行 |
| `monitor.py`、`plot.py` | 原生前端、回放与三维并集测度 |
| `vertify.py` | 独立 AC/SOCP 参考、缓存、结果比较 |
| `Network/` | FourBus、当前 Case33 与原始数据 |
| `tests/` | 物理、构域、扫描与回放回归 |

```cmd
python -m unittest tests.test_notation -v
python -m pytest -q
```
