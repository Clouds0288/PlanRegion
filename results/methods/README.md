# RCUT 方法对照运行（2026-10-01）

本目录保留方法对照实验（`experiments/compare_plan2.py`，代码在 tag `method-RCUT-v1`）中 RCUT 的精选运行：默认停滞阈值 1% 与敏感性运行 0.5%、2%，Case33 二维、三维各一次，每次墙钟 300 s、workers 16。全部方法（R、H、RB、RCUT、RCUT2）的运行与汇总报告、图在 tag `results-methods-v1` 中；RB 主线在分支 `method-RB`。各文件登记在 `results/manifest.json`（`runs`、`files`，带 sha256 与原始路径），`metrics.json` 是权威结果，不要从本目录重算。

## 停滞阈值敏感性（内域相对 AC 扫描）

| 算例 | 阈值 | 认证分区 | t_cert / 终态间隙 | 内域 FR / MR | 外域 FR |
|---|---|---|---|---|---|
| 二维 (18,25) | 0.5% | 4/4 | 29.9 s | 0.935% / 0.043% | 1.360% |
| 二维 (18,25) | 1%（默认） | 4/4 | 29.4 s | 0.989% / 0.025% | 1.431% |
| 二维 (18,25) | 2% | 4/4 | 28.8 s | 1.115% / 0.018% | 1.431% |
| 三维 (18,25,30) | 0.5% | 7/8 | 1.65% | 3.597% / 0% | 4.830% |
| 三维 (18,25,30) | 1%（默认） | 7/8 | 1.70% | 4.282% / 0% | 5.512% |
| 三维 (18,25,30) | 2% | 8/8 | 272.2 s | 4.542% / 0% | 5.546% |

三维指标在 80³ 参考扫描上计算（未决的 SOCP 格不计入）。阈值越小 N_x 割得越紧、FR 越低，但每个网架的割循环更长；阈值 2% 时三维全部分区在 300 s 内获证。RCUT 的内域含 N_x 外近似，是外侧估计：逐格对比显示它完整包含 RB 的认证内域，只多出一层 1–2 格厚的边界壳。

## 目录

```
case33_18_25/      二维：RCUT、RCUT-t0.5、RCUT-t2（workers 16，各 1 次）
case33_18_25_30/   三维：同上
  <方法>/workers_16/run_1/
    summary.json               运行汇总：各分区结果、阶段耗时、设置、commit
    metrics.json               各检查点与终态的逐格指标（AC 为主参考，SOCP 为诊断）和有效性
    timeline.csv               内外测度时间线
    solves.csv.gz              每次求解一行
    snapshots/final.json.gz    终态几何：叶锥与各网架集合
```

## 复现

```
git checkout method-RCUT-v1
python -X utf8 experiments/compare_plan2.py --case case33 --nodes 18,25,30 --mode 1 --methods RCUT --seconds 300 --repeats 1 --workers 16 --threads 1 --output results/compare_plan2/<新目录>
```

其他阈值加 `--threshold 0.005` 或 `--threshold 0.02`。各次求解固定种子，但运行受墙钟时限约束：计时、未认证分区的终态会随机器负载略有不同。
