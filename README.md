# PlanRegion：有限预算下的可规划域

`main.py` 采用逐网架生长算法。每次只深入一个网架，以完整 SOCP 认证可行点，以共享联合割收紧条件外域；全局连续查漏负责发现网架并给出停止证书。

## 设置与运行

需要有效的 Gurobi 许可证、Tkinter 和 `requirements.txt` 中的依赖。本机使用 conda `methods` 环境。在 `main.py` 顶部设置：

```python
NETWORK = Case33       # FourBus / Case33
DIMENSION = 2          # 2 / 3
CASE_TIME_LIMIT = 50.  # 构域秒数，不含独立扫描
SOLVER_THREADS = 20    # 构域及单进程扫描中，每个求解器的线程数
DIVISIONS = 80         # FourBus 每轴扫描格数
SCAN_DIVISIONS = {2: 100, 3: 60}  # Case33 每轴扫描格数
SCAN_WORKERS = 20      # 扫描进程数；并行时每个求解器用 1 个线程
FORCE_RESCAN = False   # 已有功率范围覆盖时复用；True 强制重扫
RAY_THRESHOLD = 1e-2   # 射线收益上界阈值：本网架认证面积/体积的 1%
```

| 案例 | 二维动态负荷节点 | 三维动态负荷节点 | 预算 |
|---|---|---|---|
| FourBus | 1、2 | 1、2、3 | 20000 元 |
| Case33 | 18、25 | 18、25、30 | 7 次开合变动 |

Case33 有七个可变开关，无线路升级；其他节点保留原始负荷。FourBus 未选节点负荷为零，所有节点仍须接入。`DIMENSION` 同时控制计算、扫描与回放。显式 `--load-nodes` 使用指定节点及其维数。

```text
python main.py
python main.py --case fourbus --dimension 2
python main.py --case case33 --dimension 3
python main.py --case both --dimension 2 --no-ui
python main.py --case case33 --dimension 3 --no-ui --no-scan
```

默认打开原生窗口；`--no-ui` 无窗口计算，仍保存回放。`--no-scan` 只构域。`--seconds` 修改构域上限，`--divisions` 指定重新扫描时的每轴格数。`--reference path/to/monitor.json.gz` 可把同案例、节点及预算的旧 SOCP 扫描导入独立扫描目录；强制重扫时忽略该导入。参考不参与构域。

这些参数均在 `main.py` 顶部设置。`CUT_THRESHOLD=0.02`、`CUT_PATIENCE=3` 控制局部小割停止；`POINT_TOL=1e-2` 为 kW 坐标下的近点合并距离。`RAY_THRESHOLD=1e-2` 表示：目标点直接加入本网架认证凸包的面积／体积增量，若不足当前认证测度的 1%，则跳过射线优化；首轮 G 顶点射线不筛选。参数只在回放初始化中记录一次。

## 五步流程

1. **初始化**：求 d 个轴向 MP2 和一个总负荷 MP2。上界形成全局外包络 G，总负荷解给出初始网架 A。
2. **首轮射线**：每个新网架遍历 G 的顶点，从同网架认证内域 N′ₓ 的重心朝顶点求最远可行点。已认证点复用凸性，其余首轮射线不受收益阈值筛选。
3. **最大违反量切割**：只对当前 Nₓ 的未认证顶点评分，先收入本批可行点，再仅为最大 eta 点生成一条共享割，并沿该方向补射线。SP 缓存按网架与固定近点代表去重。
4. **边界补充**：连续三次小割或当前顶点处理完后，对最新 Nₓ 剩余顶点补射线。按新增面积／体积上界排序，上界不足本网架认证测度的 RAY_THRESHOLD 时跳过。不同网架始终取并集，不混合取凸包。
5. **全局查漏**：删除被单个其他认证域包含的冗余排除项，在 G 去掉容许扩边的 G′ 后求完整物理模型，开放全部合法 x、p、y。见证直接加入其网架，再进入局部阶段；只有可靠全局上界证明覆盖，才标为 `certified`。

G′ 为各 N′ₓ 的并集。第一阶段小割比例为减少测度 / 加割前 A 外域测度；以后为 Nₓ 尚未被 G′ 覆盖部分的减少测度 / 加割前 G′ 测度。二维用面积，三维用体积，前后冻结同一个 G′。

`REGION_TAU=0.005` 是覆盖的径向精度，不是物理容差或遗漏百分比。小割、射线跳过、近点合并均不能认证剩余外域。50 秒耗尽保存已有内域和有效外包络，状态为 `time_limit`；不接受超时解、不重试、不切换求解器。G 在获得全局证书后才收紧到证书导出的外包络。

## 原生回放与独立扫描

每次只保存一份增量回放：`results/mainline/<case>_<节点序列>.json.gz`。二维、三维互不覆盖；旧回放仍可通过 `monitor.py` 打开。

```text
python main.py --case fourbus --dimension 2 --replay
python main.py --case case33 --dimension 3 --replay
python main.py --case both --dimension 2 --replay
python monitor.py results/mainline/case33_18_25.json.gz --compare results/mainline/case33_18_25_30.json.gz
```

窗口显示 G、G′、各 Nₓ、N′ₓ、SP 点、射线和割面。支持暂停、计算一步、播放、时间轴、上一割／下一割；三维可旋转并保持视角。回退只显示当时的证据，最终扫描比较独立显示。回放不调用优化器。

参考固定每个网格中心 p，自由选择全部合法 x、y，不使用构域割或已知网架。FourBus 默认每轴 80 格；Case33 默认二维每轴 100 格、三维每轴 60 格。扫描不计入 50 秒构域上限。

扫描独立保存在 `results/scans/<算例名>/<负荷节点>/budget_<预算>.npz`，例如 `results/scans/case33bw/18_25/budget_7.npz`。文件只保存各轴功率上限 `bounds` 和网格状态 `states`；各轴下限均为 0 kW，格数由 states.shape 给出。文件名与复用判断均不使用哈希。

默认 `FORCE_RESCAN=False`：同算例、节点和预算下，只要当前功率范围落在已有扫描范围内，就复用原网格，并针对本次构域重新计算遗漏率和多余率。不会缩放旧网格或把已有标签赋给新格点。超出范围时，以新旧范围的包围盒重新扫描，上限向上取整到 kW，保存后继续供下次复用。

设 `FORCE_RESCAN=True` 后按当前 DIVISIONS / SCAN_DIVISIONS 完整重扫并覆盖该扫描文件，功率范围保留新旧范围的包围盒。调整扫描精度或网架物理参数后，需要重算时使用这个开关。`--no-scan` 且未指定 reference 时仍只构域。原生回放继续携带实际使用的参考网格，可独立打开。

`SCAN_WORKERS` 独立控制扫描进程数，不受 `SOLVER_THREADS` 限制。当前设置启动最多 20 个扫描进程，每个求解器固定 1 个线程；设为 1 时，扫描在当前进程执行，求解器使用 `SOLVER_THREADS`。构域与扫描依次进行。

- 遗漏 = 参考可行而计算域未覆盖的格点数 / 参考可行格点数。
- 多余 = 计算域覆盖而参考不可行的格点数 / 计算域格点数。

G′ 和 G 分别统计。这是同一 SOCP 模型的独立网格参考，不是非凸 AC 真值；0% 网格误差不能代替连续覆盖证书。

最终结果只保留停止状态、证书上界、轴向界、内外域、求解次数和总耗时。过程保留几何增量、射线、割、局部比例和查漏证书；不再保存空队列、重复阶段汇总、分项计时或每条射线后的全局体积。

## 文件职责与验证

| 文件 | 职责 |
|---|---|
| main.py | 逐网架算法、案例与维数设置、计算和回放入口 |
| model.py | MP、SP、联合割、全局剩余域模型 |
| region.py / plot.py | 凸域几何、二维／三维并集测度 |
| monitor.py | 原生界面、增量记录与回放 |
| vertify.py | 独立 SOCP 扫描及 AC 校验接口 |
| Network/ | 网架和物理参数 |
| continuous.py | 旧连续构域，供勘察与历史对照使用 |

```text
python -m unittest tests.test_notation -v
python -m pytest tests/test_sequential_region.py tests/test_dimensions.py tests/test_startup.py tests/test_fail_fast.py -q
```

符号与显式迁移见 [docs/notation.md](docs/notation.md)。历史结果保留原算法和记录，不作机械迁移。
