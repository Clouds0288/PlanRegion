# PlanRegion

构建配电网二维／三维 SOCP 功率可行域，支持正负功率、独立 AC/SOCP 扫描和原生过程回放。当前保留逐网架主线与完整 physical 支持面认证实验。

## 环境

Python 依赖见 `requirements.txt`，测试另装 `requirements-dev.txt`；需可用的 Gurobi 许可证。本机在 CMD 中执行 `conda activate methods` 后可用 `python`。若未配置 conda 命令，直接使用：

```cmd
"C:\Users\10856\miniconda3\envs\methods\python.exe" -X utf8 main.py --help
```

以下命令在项目根目录运行，其中 `python` 可替换成上述完整路径。

## 计算

| 案例 | 二维节点 | 三维节点 | 预算 |
|---|---|---|---|
| FourBus | 1、2 | 1、2、3 | 20000 元 |
| Case33 | 18、25 | 18、25、30 | 7 次开合 |

`mode=1` 表示正负功率：正功率为负荷（PF=0.95），负功率为光伏（PF=1）；`mode=0` 为非负负荷。Case33 保持原线路参数与背景负荷，仅指定七条开关可变。主入口、AC 与 SOCP 校验共用 `Network/case33bw.py` 的 `CURRENT_LIMIT=200` A。

主线：

```cmd
python -X utf8 main.py --case fourbus --dimension 3 --mode 1 --output results/mainline/new_run
python -X utf8 main.py --case case33 --dimension 2 --mode 1 --output results/mainline/new_run
```

新方案 FourBus 三维，以及同条件重跑两种算法：

```cmd
python -X utf8 experiments/compare_support_face_fourbus_3d.py --algorithm support --seconds 300 --output results/support_face/new_run
python -X utf8 experiments/compare_support_face_fourbus_3d.py --algorithm both --seconds 300 --output results/support_face/new_comparison
```

二维实验为非负负荷；与主线比较时使用 mode=0：

```cmd
python -X utf8 experiments/test_support_face_certification_fourbus_2d.py --record --output results/support_face/new_2d
```

上述输出目录独立于归档结果。构域默认共享 300 秒总时限，事后扫描另计时；`--no-scan` 只构域，主线 `--no-ui` 关闭实时窗口。首次二维扫描 160×160、三维 80³，可用 `--divisions` 调整。实验三维过程写入 `<输出目录>/support_1/monitor.json.gz`。

主线以 SP 切割和射线补点构建局部内域；新实验以当前内域各面的支持上界进行局部认证。两者都由完整 physical MISOCP 做全局查漏。新实验按需发现方案，不逐方案枚举；局部认证完成不等于全局覆盖完成。超时保留已有证据，只有全部分区获证才报告 certified。

## 已保存的正式结果与回放

结果、参数、来源与文件哈希见 [results/manifest.json](results/manifest.json)。

| 保存内容 | 路径 |
|---|---|
| FourBus 二维主线 | `results/mainline/mode_1/fourbus_1_2.json.gz` |
| FourBus 三维对照主线 | `results/mainline/mode_1/fourbus_1_2_3.json.gz` |
| Case33 二维 200 A 主线 | `results/mainline/mode_1/case33_18_25.json.gz` |
| FourBus 三维完整 physical 支持面结果、8 个分区证书与对照 | `results/support_face/fourbus_3d_signed_physical/` |
| 全部保留扫描 | `results/scan/` |

```cmd
python -X utf8 monitor.py results/support_face/fourbus_3d_signed_physical/monitor.json.gz
python -X utf8 main.py --case fourbus --dimension 3 --mode 1 --replay
python -X utf8 main.py --case case33 --dimension 2 --mode 1 --replay
```

同一个原生前端支持逐帧、逐割、拖动进度、网架翻页、三维旋转和最终扫描比较；回放无需重求解。

保存的 FourBus 三维对照使用同一份历史 80³ SOCP 扫描：主线 1/8 分区完成，新实验 5/8 完成；内域遗漏率分别 5.66283% 和 0.07468%，多余率均为 0%。这些是原运行的离散格点指标，具体来源版本见各自 protocol，不代表两者均已全局认证，也不是 AC 对照。

Case33 的 160² 扫描是当前同限流的配对 AC/SOCP 参考，保留标签、拓扑见证、残差、逐点 CSV/NPZ 和三组对比图。历史 FourBus SOCP 格点单独置于 `results/scan/fourbus_socp/`，不导入当前 AC/SOCP 缓存。

## 扫描复用

`vertify.py` 统一管理 AC 与 SOCP。缓存位置为 `results/scan/<网络>/<节点>/<物理身份>/region_<格架摘要>.npz`；相同配置和坐标复用已有点，扩界后只补缺失点。改变电流限额、功率因数、预算、节点顺序或模式会隔离缓存。

```cmd
python -X utf8 vertify.py results/mainline/mode_1/case33_18_25.json.gz --workers 20
```

比较默认写入记录旁的 `<记录名>_comparison/`。`summary.json` 中 result_ac/result_socp/socp_ac 各含遗漏率和多余率：遗漏率分母为参考可行点数，多余率分母为计算域点数；未决标签不作为不可行。扫描校验不会改写原构域认证状态。

## 代码与检查

| 文件 | 职责 |
|---|---|
| `main.py` | 参数、计算与回放入口 |
| `model.py`、`region.py` | 物理模型、求解、构域及覆盖认证 |
| `monitor.py`、`plot.py` | 原生回放及三维并集测度 |
| `vertify.py` | 独立 AC/SOCP 参考、缓存、结果比较 |
| `Network/` | FourBus、当前 Case33 与原始数据 |
| `experiments/` | 支持面认证核心/二维入口、三维对照入口 |
| `tests/` | 当前物理、构域、认证、扫描与回放回归 |

```cmd
python -m unittest tests.test_notation -v
python -m pytest -q
```

数学符号和字段约定见 [docs/notation.md](docs/notation.md)。过时算法、诊断产物和重复源码不再放在工作目录；整理前已提交内容可从 Git 提交 `e79f18b` 恢复。
