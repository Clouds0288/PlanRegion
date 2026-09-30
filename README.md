# PlanRegion

按网架构建二维／三维 SOCP 功率域，支持正负功率、独立扫描和离线过程回放。

## 运行

使用已配置 Gurobi 许可证的 Python 环境，依赖见 `requirements.txt`（包含 Clarabel）。本机使用 `conda activate methods`。

`main.py` 顶部 `mode=1` 为默认：正功率是负荷（PF=0.95），负功率是光伏（PF=1），允许向电源反送。设 `mode=0` 计算原非负负荷域（含零边界），沿用原网架功率因数。

```powershell
python -X utf8 main.py --case fourbus --dimension 2
python -X utf8 main.py --case fourbus --dimension 2 --mode 0
python -X utf8 main.py --case case33 --dimension 3 --seconds 300 --no-ui
```

| 案例 | 二维节点 | 三维节点 | 预算 |
|---|---|---|---|
| FourBus | 1、2 | 1、2、3 | 20000 元 |
| Case33 | 18、25 | 18、25、30 | 7 次开合 |

`--load-nodes` 可指定节点序列。Case33 未选节点保留原始负荷；七个可变开关、原线路参数和径向连通约束保持不变。

默认 FourBus 二维，全部符号区共用50秒构域时限，4个求解线程。二维参考扫描每轴160格，三维80格，扫描最多20个进程，每个进程1个求解线程。`--seconds`、`--divisions` 可覆盖；`--no-scan` 只构域，`--no-ui` 不打开原生窗口。

## 算法与证书

每个符号区共用一套逐网架流程：轴向及总量支撑初始化 → 首轮射线 → 最大 SP 违反量切割 → 边界补充 → 完整物理查漏。固定符号与网架内取认证凸包，跨组只取并集；联合割只在相同符号区共享。

SP 和射线固定使用 Clarabel，原约束残差须满足 `1e-8`；联合割仍由 Gurobi LP 生成，完整网架搜索仍用混合整数 SOCP。无数值重试或备用求解器。`tau=0.005` 为覆盖的径向精度，只有全部分区取得覆盖证书才报告 `certified`。超时保存已有认证内域和有效外域。

## 结果与回放

每组只保留一份原始增量记录，例如：

`results/mainline/mode_1/fourbus_1_2.json.gz`

用现有 `monitor.py` 原生前端回放：

```powershell
python -X utf8 main.py --case fourbus --dimension 2 --replay
python -X utf8 main.py --replay results/mainline/mode_1/fourbus_1_2.json.gz
```

窗口支持播放、逐帧、逐割、拖动进度和网架翻页；三维可旋转。显示全局内外域、当前网架、SP 点、射线、割和最终扫描对比。计算与回放共用同一前端、同一份 version=4 记录，也可直接执行 `python monitor.py <记录>`。正负功率均以真实 kW 坐标显示。

扫描按案例、模式、节点、功率因数及预算存放在 `results/scans/<case>_signed/`；已有上下界包含本次范围时复用原网格。`--force-rescan` 强制重扫；修改物理参数或需要更高网格精度时使用。`--reference <新主线记录>` 可导入相同配置的参考扫描。

内外域分别报告：遗漏率 = 未覆盖参考可行格点 / 参考可行格点；多余率 = 覆盖参考不可行格点 / 计算域格点。参考独立求解完整 SOCP，不使用构域割或已知网架限制；同符号、同网架的可行端点可由凸性认证中间格点。这是 SOCP 网格参考，不是非凸 AC 真值。

## 核心文件与验证

- `main.py`：参数和运行入口，默认 `mode=1`。
- `model.py`：数学模型、MP、SP、联合割与连续锥求解。
- `region.py`：全局搜索、符号分区、逐网架构域和几何。
- `monitor.py`：增量记录及 Python 原生前端回放。
- `vertify.py`：AC 校核和独立扫描；`Network/`：算例参数。
- `continuous.py`、`survey.py`、其余 `experiments/`：勘察及独立研究对照。

```powershell
python -m unittest tests.test_notation -v
python -m pytest tests/test_fourbus_signed.py tests/test_case33_signed.py tests/test_sequential_region.py tests/test_scan_reuse.py tests/test_monitor.py -q
```

数学字段与迁移约定见 [docs/notation.md](docs/notation.md)。
