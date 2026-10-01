# 方法对照实验的精选结果（2026-10-01）

本目录归档 `experiments/compare_plan2.py` 在 Case33 上的精选运行，统一登记在 `results/manifest.json`
（`method_archive`、`runs`、`files`，每个文件带 sha256 与原始路径）。各检查点的中间快照和分区日志没有归档，
只留在本地的 `results/compare_plan2/` 里。`metrics.json` 是权威结果，不要从本目录重算。

## 方法与代码版本

| 方法 | 内容 | 代码 tag | 产生结果的提交 |
|---|---|---|---|
| R | 径向锥夹逼：x 自由的锥 MISOCP 外界、近端两网架覆盖、体积准则 | `method-R-v1`（4c5b2a4） | 2D 为 b619864，3D 为 4c5b2a4（R 的代码路径相同） |
| H | 方案二：R 的 A 阶段 + 逐网架支撑查询 B + 覆盖证书 C（对照基线） | `method-H-v1`（b619864） | b619864 |
| RB | R + B，去掉覆盖证书；认证内域 I_R ∪ P_x | `method-RB-v1`（4c5b2a4） | 4c5b2a4 |
| RCUT | R + 主线对偶割（带 OBBT 包络），停滞的 N_x 当作网架区域（外侧估计，不是认证内域） | `method-RCUT-v1`（3e8b0a1） | 3e8b0a1 |
| RCUT2 | RCUT + 主线射线内域认证，内域为认证点凸包 N'_x | `method-RCUT2-v1`（3e8b0a1） | 3e8b0a1 |

整个归档的 tag 为 `results-methods-v1`。各运行的生成提交也记在 `summary.json` 的 `commit` 字段，运行时代码树均为干净状态。

## 目录

```
case33_18_25/        2D：R、H 在 workers 1 和 16 下各 3 次；RB、RCUT、RCUT2 在 workers 16 下各 1 次
case33_18_25_30/     3D：R、RB、RCUT、RCUT2 在 workers 16 下各 1 次
  RCUT-t0.5、RCUT-t2、RCUT2-t0.5、RCUT2-t2   停滞阈值 0.5% 与 2% 的敏感性运行（默认阈值 1%）
  comparison.md、comparison.csv、*.png       汇总报告与图
  <方法>/workers_<w>/run_<r>/
    summary.json               运行汇总：各分区结果、阶段耗时、设置、commit
    metrics.json               各检查点与终态的逐格指标（AC 为主参考，SOCP 为诊断）和有效性
    timeline.csv               内外测度时间线
    solves.csv.gz              每次求解一行
    snapshots/final.json.gz    终态几何：叶锥与各网架集合
```

参考扫描（已在版本库中）：`results/scan/case33bw/18_25/.../region_17ceb84375dc2cc1.npz`（2D，160×160），
`results/scan/case33bw/18_25_30/.../region_3837a4b56856fbc1.npz`（3D，80×80×80）。

## 主要结果（workers 16，每次运行墙钟 300 s，内域相对 AC 扫描）

| 方法 | 2D 认证 / t_cert | 2D 内域 FR / MR | 3D 认证 / 终态间隙 | 3D 内域 FR / MR |
|---|---|---|---|---|
| R | 4/4，28.4 s | 0.75% / 0.28% | 3/8，3.69% | 2.42% / 0.77% |
| H | 4/4，50.8 s（3 次中位） | 0.83% / 0.08% | 未运行 | 未运行 |
| RB | 4/4，27.8 s | 0.83% / 0.14% | 3/8，2.35% | 2.58% / 0.16% |
| RCUT（1%） | 4/4，29.4 s | 0.99% / 0.02% | 7/8，1.70% | 4.28% / 0.00% |
| RCUT2（1%） | 4/4，32.4 s | 0.83% / 0.10% | 3/8，3.64% | 2.55% / 0.29% |

RCUT 的内域含割平面外近似，不是认证内域：`metrics.json` 中 `inner_exact=false`，有效性只检查外界。逐格对比显示，
它的内域完整包含 RB 的认证内域，只多出一层 1–2 格厚的边界壳。阈值敏感性见各目录的 `threshold_sensitivity.png`。

## 复现

```
git checkout method-RB-v1
python -X utf8 experiments/compare_plan2.py --case case33 --nodes 18,25,30 --mode 1 --methods RB     --seconds 300 --repeats 1 --workers 16 --threads 1 --output results/compare_plan2/<新目录>
```

RCUT/RCUT2 用 `method-RCUT-v1` 并加 `--methods RCUT,RCUT2`；其他阈值加 `--threshold 0.005` 或 `--threshold 0.02`。
各次求解固定种子，但运行受墙钟时限约束：计时、未认证分区的终态会随机器负载略有不同。
