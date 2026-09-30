# PlanRegion

按网架构建二维／三维 SOCP 功率域，支持正负功率、独立 AC / SOCP 扫描和离线过程回放。

## 运行

使用已配置 Gurobi 许可证的 Python 环境，依赖见 `requirements.txt`（包含 Clarabel）。本机使用 `conda activate methods`。

`main.py` 顶部 `mode=1` 为默认：正功率是负荷（PF=0.95），负功率是光伏（PF=1），允许向电源反送。设 `mode=0` 计算原非负负荷域（含零边界），沿用原网架功率因数。

```powershell
python -X utf8 main.py --case fourbus --dimension 2
python -X utf8 main.py --case fourbus --dimension 2 --mode 0
python -X utf8 main.py --case case33 --dimension 2 --no-ui
python -X utf8 main.py --case case33 --dimension 3 --seconds 300 --no-ui
```

| 案例 | 二维节点 | 三维节点 | 预算 |
|---|---|---|---|
| FourBus | 1、2 | 1、2、3 | 20000 元 |
| Case33 | 18、25 | 18、25、30 | 7 次开合 |

`--load-nodes` 可指定节点序列。Case33 未选节点保留原始负荷；七个可变开关、原线路参数和径向连通约束保持不变。

默认 Case33 三维，全部符号区共用300秒构域时限，20个求解线程。首次二维扫描每轴160格，三维80格，AC 与 SOCP 使用相同坐标；扩界后保持步长并增加格数。扫描最多20个进程，每个进程1个求解线程。`--seconds`、`--divisions` 可覆盖；`--no-scan` 只构域，`--no-ui` 不打开原生窗口。

`Network/case33bw.py` 的 `CURRENT_LIMIT = 200.` 是主线及独立校验入口的共同限流配置，单位 A。全部37条候选线路（包括联络线）逐条约束，闭合线路不超过200 A，断开线路电流为零。可将该值改成37项数组，顺序与 `Network/data/case33bw.m` 的 branch 表一致。这是本实验指定的限额，并非原始 IEEE 33 数据提供的额定电流。

## 算法与证书

每个符号区共用一套逐网架流程：轴向及总量支撑初始化 → 首轮射线 → 最大 SP 违反量切割 → 边界补充 → 完整物理查漏。固定符号与网架内取认证凸包，跨组只取并集；联合割只在相同符号区共享。

全局查漏统一采用完整 physical MISOCP。`RemainingRegionModel` 已移除 light/physical 模式选择，`MasterProblem` 已移除 `cuts_only`，所有主线全局查询均包含建设变量、运行变量及原 SOCP 约束。API 调用不再传这两个选项；`main.mode=0/1` 仍表示非负/正负功率模式。

有限逐线路电流上限时，SP 直接复用既有 Gurobi 连续 SOCP；其他 SP 及射线使用 Clarabel。射线按电流界缩放锥，并在求取内域见证时预留既有微小裕量，最终原约束残差仍须满足 `1e-8`。联合割由 Gurobi LP 生成，完整网架搜索用混合整数 SOCP。求解路径由模型配置预先确定，无数值重试或备用求解器。`tau=0.005` 为覆盖的径向精度，只有全部分区取得覆盖证书才报告 `certified`。超时保存已有认证内域和有效外域。

## 结果与回放

每组只保留一份原始增量记录，例如：

`results/mainline/mode_1/fourbus_1_2.json.gz`

用现有 `monitor.py` 原生前端回放：

```powershell
python -X utf8 main.py --case fourbus --dimension 2 --replay
python -X utf8 main.py --replay results/mainline/mode_1/fourbus_1_2.json.gz
```

窗口支持播放、逐帧、逐割、拖动进度和网架翻页；三维可旋转。显示全局内外域、当前网架、SP 点、射线、割和最终扫描对比。计算与回放共用同一前端、同一份 version=4 记录，也可直接执行 `python monitor.py <记录>`。正负功率均以真实 kW 坐标显示。

主图及网架小图按各自分区初始化后的外包络固定坐标，四周保留 6% 留白；SP、加割和射线过程不改变范围。图 A 画当前分区的真实外包络、条件外域和认证内域，旁边的全局总览包含尚未处理分区，并用蓝框标出主图范围。“主图展开全局”可将 A 切换为全局视图；构域结束后 A 显示最终全局内外域。

过程图层包含全局搜索点、SP 不可行点、初始点、未决点、全网架不可行点以及射线锚点/目标/认证点。全局候选在某个网架下被 SP 否决后以叉号保留，不冒充全网架不可行证明；历史回退仅显示当时已出现的证据。扫描对比按最终域和参考可行格的完整边缘取范围，回退历史帧时保持最终对比。

`monitor.py` 中 `RunMonitor` 负责事件录制、控制与扫描指标；`NativeWindow` 先重建当前帧及其分区视口，再分别刷新状态栏、A/B 过程面板和独立的 C 校验面板。二维与三维共用几何入口和过程图层，坐标选择不依赖过程标记。

独立 AC 与 SOCP 扫描统一由 `vertify.py` 负责。两者复用 `MasterProblem` 的约束装配；SOCP 保留 `P²+Q²≤v·ell`，AC 使用等式。每个扫描点允许选择预算内任一合法径向拓扑，不受主线已发现拓扑限制。AC 的固定拓扑潮流证书可直接接受；只有全部拓扑的 AC 必要条件均不满足，或完整非凸 AC 模型证明不可行，才判不可行。单个拓扑越限、潮流迭代不收敛和求解超时不判不可行。两套标签独立求解，不通过凸性、插值或构域结果填充。

SOCP 固定点采用最小违反量 η 判定：功率平衡及压降等式的绝对误差受 η 约束，其他限制保持原式，最小化 η。η=0 就是原可行问题；正的全局下界证明不可行，可行见证须满足 η 加原模型求解残差不超过 `1e-8`。这避免直接对不可行锥问题求零目标引起的数值误判，未降低验收精度。

Case33 主线、SOCP、AC 使用相同逐线路电流限额、电压界、功率因数、背景负荷和开关预算。同条件下 AC 域应包含于 SOCP 域。初次扫描范围采用完整 SOCP 模型的坐标全局上界，并与本次构域范围合并，覆盖两种参考。电流限额、其他物理参数、节点顺序、预算、功率因数或模式变化都会使用不同缓存。

缓存保存到 `results/scan/<案例>/<节点序列>/<物理配置标识>/region_<格架摘要>.npz`。每个功率区域文件包含相同坐标上的 AC、SOCP 标签、可行拓扑和残差；`states/witness_x/residual` 属于 AC，`socp_states/socp_witness_x/socp_residual` 属于 SOCP。首次确定坐标原点和步长，随后扩界保持旧坐标：全覆盖直接复用，部分覆盖优先复用并补齐缺失点，各轴格数可以不同。每批结果原子保存，中断后分别续算两种参考的零标签点；`1/-1/0` 分别代表可行、已证不可行和未决，未决时不发布完整误差率。

旧不限流、不同限额和旧格式缓存保留原位，不迁入本次同限流参考。`--reference <区域缓存.npz>` 可显式导入同版本、同物理配置的参考，也接受包含区域缓存路径的新记录。`--force-rescan` 备份后重扫；更改 `--divisions` 只影响首次建格，需要加密已有格点时用 `vertify.py --refine 2`，旧格点仍可复用。

主要输出六项指标，实验结果使用最终内域 G′：实验结果—AC、实验结果—SOCP、SOCP—AC 三组比较，各给遗漏率和多余率。每组左侧为计算域，右侧为参考域：遗漏率 = 参考可行而计算域未覆盖的点 / 参考可行点；多余率 = 计算域覆盖而参考不可行的点 / 计算域点。空分母记为 null；原内外域对 AC 的 metrics 字段保留。统计基于离散格点，不代表连续区域的严格体积误差。

前端 C 面板同时显示六项指标，并可切换三组对比图。共同可行区域为蓝色，多余为橙色，遗漏为红色，双方不可行留白；逐格绘制，不用凸包填充空洞。

正常执行 `main.py` 会在构域后自动复用或补算两套扫描。相同配置下的已有记录可单独补算，无需重新构域：

```powershell
python -X utf8 vertify.py results/mainline/mode_1/case33_18_25.json.gz --workers 20
```

逐点对比默认放在构域记录旁的 `<记录名>_comparison/`，例如 `results/mainline/mode_1/case33_18_25_comparison/`；`vertify.py --output` 可指定其他目录。`comparison.npz` 和 `comparison.csv` 保存功率坐标、`ac_states/socp_states`（CSV 列名 `ac_state/socp_state`）及原结果的 `inner/outer` 成员标签；`summary.json` 的 `comparisons` 保存六项主要指标及计数。功率按 `load_nodes` 顺序排列，单位 kW。读取示例：

```python
import numpy as np
with np.load('results/mainline/mode_1/case33_18_25_comparison/comparison.npz') as data:
    power = data['power']
    ac = data['ac_states'] == 1
    socp = data['socp_states'] == 1
    inner, outer = data['inner'], data['outer']
    extra = power[inner & ~ac]
    missed = power[ac & ~inner]
```

校验完成不改变原构域的 `result.certified`。若原结果为 `time_limit`，仍表示尚未获得全部符号区的覆盖证书。

## 核心文件与验证

- `main.py`：参数和运行入口，默认 `mode=1`。
- `model.py`：数学模型、MP、SP、联合割与连续锥求解。
- `region.py`：全局搜索、符号分区、逐网架构域和几何。
- `monitor.py`：增量记录及 Python 原生前端回放。
- `vertify.py`：AC / SOCP 全拓扑扫描、区域缓存复用和三组结果对比；`Network/`：算例参数。
- `plot.py`：主线三维并集测度。
- `experiments/test_support_face_certification_fourbus_2d.py`：按需发现方案的支持面认证实验，用线性联合割外域查漏并取得全局覆盖上界证书；最终结果调用 `vertify.py` 与 AC 比较。`budget_schemes` 统一位于 `vertify.py`，用于独立 AC 扫描和拓扑回归，不限制构域模型搜索。

旧 `continuous.py`、勘察与算法对照入口及配套网页已退役。本次合并了 `ac_validation.py`，移除旧 SOCP 审核入口 `experiments/audit_case33_current.py`、`experiments/compare_case33_current.py` 和已迁移拓扑枚举的 `experiments/fourbus_outer_volume.py`。`tests/` 保留主线数值、AC 扫描、缓存、原生回放、符号契约和支持面实验回归；`legacy_case33.py` 仅作为独立物理校核夹具，不参与主入口。

```powershell
python -m unittest tests.test_notation -v
python -m pytest tests -q
```

数学字段与迁移约定见 [docs/notation.md](docs/notation.md)。
