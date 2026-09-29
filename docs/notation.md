# 数学符号与代码变量规范

## Case33 正负接入实验扩展（2026-09-29）

`experiments/case33_signed.py` 复用 FourBus 分区构域、独立扫描和回放，主线不变。[LOAD_NODES](../experiments/case33_signed.py) 为 (18,25)/(18,25,30)，预算7次开合，全部符号区共用60/300秒。未选节点保留原始有功、无功背景负荷；选中节点负荷PF=0.95、光伏PF=1，与 FourBus 实验一致。

Case33 原数据 rateA=0 表示没有线路额定容量数据，不人为新增容量或逆变器圆。[voltage_flow_bounds](../experiments/fourbus_signed.py) 返回由两端电压上界推出的有效 ell/P 界：ell<=(sqrt(vmax_i)+sqrt(vmax_j))²/(r²+χ²)，两端功率幅值均不超过 sqrt(max(vmax_i,vmax_j)*ell)。[port_bounds](../experiments/fourbus_signed.py) 对无限额线路使用此有效功率界，有限额线路仍使用原额定有功容量，节点幅值上界仍为相邻走廊最大界之和。

`PortPhysics._build_variable_bounds` 对 Case33 再由全网无功平衡得到 ell<=(min(source_qmax,source_smax)-sum(fixed_q)/base-minimum_variable_q)/χ，并与电压推导界取小值；支路P/Q界据此收紧。这些是物理可行集的有效界，不添加运行限制。SP的零潮流松弛量上界补入固定背景负荷最大P/Q。新网架先收入实际见证，再实际检查零接入点，只有SP可行才收入原点；不把Case33原点默认认证。

共享入口 `run` / `main` 新增 network_type、load_nodes、budget 参数；`main.time_limits` 选择各维时限，默认仍为 FourBus；`build_partition.budget`、`scan_line.network_type`/`budget`、`scan_reference.network_type`/`budget` 显式传递，不修改变量、状态或记录布局。Case33 结果和扫描分别位于 results/case33_signed 与 results/scans/case33_signed；FourBus 目录保留。回放标题以接入节点标识案例，旧记录可原样回放。

连续求解的等式消元：`solve_conic` 对固定变量消元后仍存在的 Aeq*y=beq 用 SVD 表示 y=offset+basis*z；[solve_conic.offset](../experiments/port_subproblem.py) 为特解，[solve_conic.basis](../experiments/port_subproblem.py) 为零空间正交基。只求自由坐标z，随后还原原状态并核验全部原始约束；消除压降等式中的数值残差，不改变1e-8认证门槛。无等式的SP仍直接求解原自由变量，不重试、不切换求解器。

Case33 电流锥固定等价变换为 norm((20P,20Q,v-100ell))<=v+100ell，以改善小电流附近的数值尺度；原约束、原锥取割方向和原始残差核验不变。连续SP收敛目标为1e-10，物理接受仍要求 eta+原约束残差<=1e-8。两案例统一使用 `PortSubProblem`，没有失败后重试或求解器切换。

Case33 射线的锥head内缩量随该等价变换乘100，保证原电流锥仍有不小于原先的保守余量；SP的margin仍为0，完整物理域与扫描不变。

Case33 构域与扫描的Gurobi模型从第一次装配起固定 ScaleFlag=1、BarQCPConvTol=1e-10；FourBus 扫描保留 ScaleFlag=0、BarQCPConvTol=1e-8。NumericFocus=0 和原1e-8物理残差门槛不变。

Case33 参考扫描固定 Aggregate=1；Aggregate=0 在 p=(-975.05625,-1861.5875,-11588.36875) kW 处出现可复现的错误不可行判定。该点按相同完整物理模型、Aggregate=1 求得可行解，原约束残差3.46e-10；初次80³参考中对应整行已独立重扫并替换。不是用构域的认证标签填充参考，也不在运行时增加重试或求解器切换。

独立扫描沿每行递归检查区间端点：同一符号区内，若两端完整MP返回同一网架及合格可行解，由凸性认证中间格点；其余区间二分继续求解，绝不从不可行端点推断中间不可行。所有不可行格点均由完整MP实际证明。Case33保留求解器热启动，FourBus仍逐次reset；扫描不利用构域中的割、凸包或标签。该扫描加速不改变采样坐标或分类定义。

回放可读取已有的 reference=None 构域记录：先显示真实迭代，明确标记参考扫描尚未完成；参考图和误差在完整扫描完成后生成，不以零误差代替缺失结果。记录格式不变。

## 独立扩展：固定功率因数的负荷 / 光伏分区（2026-09-29）

`experiments/fourbus_signed.py` 只测试 FourBus，不改变主线默认模型。MODE=0 为原正负荷模型；MODE=1 枚举接入符号 sign∈{+1,-1}^d。持久化 p、vertices、割均为实际带符号 kW：p>0 是负荷，p<0 是光伏；每个节点只有一种接入类型，不允许任意负荷与光伏相互抵消。负荷 PF=0.95，光伏 PF=1，因此 q=tan(arccos(0.95))*p（正区间），q=0（负区间）。固定背景负荷仍由 Network 提供。

为复用原非负几何与求解接口，每个固定符号区内部暂用幅值 u=sign*p>=0。`PortPhysics.add_operation` 将原 p 形参中的幅值乘 sign 后送入原运行方程；此局部适配是显式坐标迁移，不新增自由度。SP 返回的 cut 在幅值坐标有效；写入记录时 beta_signed=beta*sign，所有点、射线、内外域顶点同时乘 sign。状态向量 P/Q/ell/v 的参考方向与切片不变。

| 新增量 | 固定代码映射 | 定义 |
|---|---|---|
| 接入模式 | [MODE](../experiments/fourbus_signed.py) | 0 原正负荷；1 正负接入分区 |
| 接入符号 | [PortPhysics.sign](../experiments/fourbus_signed.py) | +1 负荷，-1 光伏，固定参数，不是整数决策变量 |
| 节点初始功率幅值界 | [port_bounds](../experiments/fourbus_signed.py) | 节点相邻走廊各最大有功容量之和，kW；由两端容量界与节点平衡推出，不用电源容量截断内部交换 |
| 构域总时限 | [TIME_LIMITS](../experiments/fourbus_signed.py) | 二维20秒、三维60秒；全部符号区共用，包含建模、几何与原始记录，扫描与回放导出另计 |
| 扫描精度 | [SCAN_DIVISIONS](../experiments/fourbus_signed.py) | 全部正负范围每轴格数，二维160、三维80；网格中心逐点完整 SOCP，无算法割或认证凸包参与 |
| 扫描范围 | `axis_lower`、`bounds` | 下界和上界均为实际 kW；点为 axis_lower+(index+0.5)*(bounds-axis_lower)/divisions |

MODE=1 接根支路允许反向 P/Q，保留原源端视在功率上限，线路两端均受原有功容量限制。由 -P+r*ell<=capacity 与 P<=capacity 推出 ell<=2*capacity/r；无功界由源端能力和区域内最大无功注入推出。此界不是新增额定电流。各分区用幅值总和上界，不把原纯负荷总量界用于限制负荷与光伏的内部交换。

认证域按 (sign,x) 独立保存；只在同一 sign 内共享联合割。首轮射线、最大 eta 切割、连续小割、补充射线和完整物理查漏沿用主线流程。各区动态分配剩余总时间；超时只保存已取得的内域和有效外域。只有全部符号区获得查漏证书才报告 certified。FourBus 零功率可解析认证，每个已知合法网架加入原点，使径向扩边包含原认证域。

SP 评价点在求解前统一朝同网架认证凸包重心内移至多 POINT_TOL/4 kW（当前0.0025 kW），避开多条割交点上的退化；最大 eta 在这些实际评价点上比较。内移点可行不认证原外顶点，原始外包络与连续查漏不内移；内移点不可行的有效割也分离原顶点，因为锚点可行。SP 加入解析有效上界 eta<=max_i(max(1,abs(q_ratio_i))*port_bounds_i/base)：零潮流、单位电压提供该上界内的可行松弛解，故不改变最小 eta，也不影响 eta=0 时的任何物理可行点。它避免无界松弛量使内点法初始步跨越过大尺度。

连续 SP 和射线固定使用 Clarabel，由 `experiments/port_subproblem.py::solve_conic` 将原线性约束、变量界和 SOCP 转换为 Az+s=b、s∈K，解析消去固定 x、固定 p、未选型号的零潮流/电流和接通走廊的零压降松弛，不改变物理域或输出字段。重新计算原始线性、二次约束与变量界残差。`PortSubProblem` 沿用原 eta 目标及 eta+残差<=PLANNING_TOL=1e-8 的可行点认证门槛；非零 eta 只排序，是否能排除候选由原 Gurobi 取割 LP 严格验证。`fourbus_signed.port_ray` 直接代入 p=anchor+lambda*(target-anchor)，最大化射线参数 lambda；省去冗余的 kW 功率等式，认证门槛仍为原物理残差1e-8。没有失败重试或备用求解器。

射线采用固定 `RAY_CONE_MARGIN=1e-6`，将已选型号和源端锥的 head 内缩该标幺量，避免边界舍入把认证点推到锥外；这是保守内点搜索，仅作用于射线，不修改 SP、全局模型、割或参考扫描。连续求解允许 Solved / AlmostSolved 返回候选，但二者一律按原物理残差核验后才可入内域；不把局部求解状态当作覆盖证书。完整 MP 保留1e-10目标收敛精度；参考扫描使用恒零目标、NumericFocus=0、BarQCPConvTol=1e-8，每个点先清除上一个点的解状态。扫描仍要求完整模型最优或不可行终止，可行解的原约束残差不超过1e-8。全部实验 Gurobi 模型首次设置 BarHomogeneous=1、Aggregate=0、ScaleFlag=0。扫描使用自动预处理；Presolve=0 与 NumericFocus=3 的组合在带符号算例上出现过可复现的错误不可行判定，未用于正式扫描。

独立扫描缓存按 MODE / 有序节点 / PF / 预算存放，使用上下界包含关系复用，不使用哈希。FORCE_RESCAN 强制重扫；修改模型物理参数后需显式重扫。MR=遗漏可行格点/全部参考可行格点；FR=多余不可行格点/计算域内格点，均乘100，inner/outer分别报告。参考为同一SOCP的独立扫描，不称为AC真值。

实验记录 `recording.json.gz` 使用独立 version=1，保存 settings、history、result、reference；原生 version=4 格式不变。配套 `replay.html` 读取转换后的真实带符号坐标，逐帧展示 G、G'、N_x、N'_x、SP点、割和射线；不同符号区及网架只取并集，不跨组取凸包。

## 主线迁移与精简（2026-09-29）

本节覆盖下方历史条目中的入口、默认值和记录约定；历史结果不改写。

扫描复用迁移：主入口默认自动复用同算例、同有序负荷节点、同预算的完整 SOCP 扫描；用各轴功率范围包含关系判断，不使用哈希或源码指纹。仅保存一份独立扫描 `results/scans/<network.name>/<节点序列>/budget_<预算>.npz`，回放文件仍按原格式内嵌参考网格。旧回放可通过 `--reference` 导入该目录；启用强制重新扫描时忽略该导入。`--no-scan` 且未指定 reference 时仍只构域。

| 新增量 / 接口 | 固定代码映射 | 定义 |
|---|---|---|
| 强制重新扫描开关 | [FORCE_RESCAN](../main.py)、[run.force_rescan](../main.py)、[scan_reference.force_rescan](../main.py) | 默认 False；True 时完整重算扫描，覆盖该算例、节点和预算的扫描文件 |
| 独立扫描目录 | [SCAN_OUTPUT](../main.py)、[run.scan_output](../main.py) | 默认 results/scans；与构域回放目录 OUTPUT 分开 |
| 扫描文件位置 | [scan_path](../main.py) | 算例名 / 有序负荷节点 / budget_预算.npz；不同维数、节点顺序及预算分开保存 |
| 保存参考网格 | [save_scan](../main.py) | 只写入既有 bounds、states；扫描下限统一为零，单位 kW；不保存旧构域的误差指标 |
| 范围复用与扫描 | [scan_reference](../main.py) | 逐轴请求上限不超过已有 bounds 时复用原网格；否则对新旧范围的包围盒重新扫描，上限向上取整到 kW，避免数值微差重复扫描；强制重算同样保留已覆盖的范围 |

复用保留扫描原有 bounds 与 states 坐标，不缩放、不插值，也不把旧标签移到新格点。divisions 仅在实际重新扫描时生效；已覆盖范围内改变格数不会自动重扫，需 FORCE_RESCAN=True。独立文件不记录耗时，复用时沿用回放已有的 scan_seconds=None 口径，实际新扫描仍记录其耗时；每次仍按本次构域结果重算 inner/outer 的遗漏率和多余率。物理参数修改后由用户显式要求重扫，不作模型哈希认证。

参考扫描的固定数值设置迁移：`_scan_line` 对 FourBus 设置 Presolve=0，Case33 等其他网络保留自动值 -1；NumericFocus=3、BarHomogeneous=1、Aggregate=0、ScaleFlag=0 及物理精度保持不变。FourBus 三维 60³ 扫描的格点 (15,11,38) 在旧设置下可重复出现 NUMERIC，关闭预处理后取得不可行证明，完整 216000 点扫描通过；Case33 的原扫描设置保留，避免不必要的性能损失。设置在每个模型首次求解前确定，没有逐点重试或数值失败后的分支。

- `experiments/sequential_region.py` 的逐网架算法迁入 `main.py`，函数名、物理模型、坐标、割布局和停止证书不变；删除原实验入口，不保留同义包装。旧 `build_continuous_region` 移至 `continuous.py`，仅供勘察和历史对照使用，不再由主入口调用。
- [DIMENSION](../main.py) 为 2 或 3，控制默认动态负荷节点：FourBus 为 (1,2) / (1,2,3)，Case33 为 (18,25) / (18,25,30)。[main.dimension](../main.py) / `--dimension` 覆盖该值；显式 `load_nodes` / `--load-nodes` 则使用完整节点元组及其维数。未选中的 Case33 节点保持原始负荷。
- 默认 `NETWORK=Case33`、`CASE_TIME_LIMIT=50`、`SOLVER_THREADS=4`。新增参数 [CUT_THRESHOLD](../main.py)=0.02、[CUT_PATIENCE](../main.py)=3、[POINT_TOL](../main.py)=1e-2 kW；分别传给既有 threshold、patience、point_tol，不改变数学语义；POINT_TOL 保留当前用户设置。
- 参数位置迁移：SCAN_DIVISIONS、SCAN_WORKERS 从 Network/case33bw.py 移至 main.py 顶部，RAY_THRESHOLD 从 model.py 移至 main.py 顶部；名称、取值、单位与算法含义不变，旧模块不保留别名。SCAN_DIVISIONS={2:100,3:60} 仍用于 Case33，FourBus 仍用 DIVISIONS=80；SCAN_WORKERS=12 控制扫描进程数，独立于 SOLVER_THREADS=4。并行扫描每个求解器用一个线程，单进程扫描使用 SOLVER_THREADS。RAY_THRESHOLD=1e-4，首轮 G 顶点射线仍不受该阈值筛选。
- 主循环注释分五步：初始化、首次 G 顶点射线、最大 eta 切割、边界补充、完整物理查漏。每阶段只深入一个网架；首轮射线不受收益阈值筛选。
- [recording_path](../main.py) 保留按节点区分文件的规则，输出目录改为 `results/mainline`。原生增量回放仍为 version=4；图形、SP 点、射线、割、查漏上界和最终内外域保持原格式。
- 新结果只包含 status、certified、axis_bounds、coverage_bound、counts、timing、inner、outer。counts 仅保留 initial/sp/cuts/ray/global_search，timing 仅保留 total_seconds。参数在回放初始化中只记录一次；最终内外域保留原 kW 坐标和网架信息。
- 退役的结果字段保留历史定义，但新运行不写入：method/budget/tau/threshold/patience/point_tol/ray_threshold/measure_unit/max_total_bound 的重复结果副本，stages、coverage_checks、rejected_points；counts.cut_lp/ray_skipped/certification/stages，分项 timing，unknown_points。cut_lp 的旧含义为实际取割 LP 次数，新流程每条采用割只解一次 LP，完成运行可从 cuts 读出。
- 退役的阶段汇总包括 boundary_rays、local_gap_area、remaining_area、inner_area；每次射线后不再额外计算全局体积。割后仅保存 area_ratio 和 small_cuts，回放几何足以重算 removed_area/area_denominator；初始和最终几何不变。查漏不再另存排除规模或重复目标表，residual_end 仍记录可靠 coverage_bound 和 complete。
- `run_case` 入口退役，由既有 [run](../main.py) 编排构域和独立扫描、[main](../main.py) 选择案例；不存在算法重试、备用求解器或超时解接受。达到 50 秒仅保存当时有效内外域，不能把 time_limit 标为 certified。

## 独立扩展：逐网架生长实验（2026-09-28）

`experiments/sequential_region.py` 每阶段只认证一个固定网架；x、p、state、eta、联合割与物理约束沿用主线。d 个轴向和一个总负荷方向的自由网架 MP2 提供初始外界，取总负荷方向的网架作为 A。认证内域仍为同网架点的凸包；从不把面积或体积停滞后的外域当成可行内域。讨论记号 G 为全局外包络，N_x / N'_x 对应 records[x] 的 outer / inner，G'=并集 N'_x 对应原全局 W；不重命名既有代码字段。

2026-09-29 显式流程迁移：取消历史候选队列 frontier、全网架定点认证缓存 certifications 和 certify_power 求解流程。下表保留退役记录。stage_candidates 改为当前网架自身外包络的未认证顶点；每阶段结束前补一轮固定网架射线，随后立即进行完整物理连续查漏。历史结果不改写；version=4、counts.certification、timing.certification_seconds、rejected_points、unknown_points 保留，新运行分别为 0、0、空列表、空列表。

射线等式统一除以 network.base，与运行约束采用相同标幺尺度；仅作等价行缩放，power / anchor / 返回 p 仍为 kW，lambda 范围和物理认证阈值不变。射线采用与 SP 相同、运行前固定的 numeric_focus，FourBus 为 0、Case33 为 3；不存在失败后的参数切换。

效率版本统一为所有射线预先设置 Aggregate=0、ScaleFlag=0，与 SP 保持一致，避免预处理聚合和缩放放大原始节点平衡残差；仍只求解一次，并按原模型 MaxVio<=1e-8 接受。物理约束和射线目标不变，运行中不按结果切换参数。

射线子问题的 BarQCPConvTol 从原实验的 1e-8 收紧至 1e-9，与 SP 相同，用于局部目标收敛；原 MP 默认值为 1e-10。该局部最优性设置不作为覆盖证书：射线仍须 OPTIMAL 且原模型 MaxVio<=PLANNING_TOL=1e-8 才收入实际可行点；SP 和全局查漏的参数、物理可行性门槛、覆盖上界门槛均保持不变。局部边界未完全补足时由下一次连续查漏返回见证。

| 数学量 / 定义 | 固定代码映射 | 单位与用途 |
|---|---|---|
| 当前活动网架 | [build_sequential_region.x](../main.py) | 完整二元型号向量；每阶段只有一个 |
| 当前网架外包络的未认证顶点 | [stage_candidates](../main.py) | 公共评价箱归一化坐标；仅扣除本网架 N'_x，不再由其他网架认证边界生成候选 |
| 射线参数 λ | [ray_support.ray_fraction](../main.py) | [0,1]；p=anchor+λ(power-anchor)，anchor 必须来自同网架认证凸包；原点不默认可行 |
| 固定 p、开放全部 x/y 的最小违反量 V(p) | 已退役：`experiments/sequential_region.py::certify_power` | 历史定义保留：原始物理模型最小违反量；本实验新流程不再调用 |
| 局部割面积减少比例 | [build_sequential_region.area_ratio](../main.py) | 第一阶段为减少面积/加割前 O_A 面积；后续为未认证部分减少面积/加割前 W 面积；前后固定同一个 W |
| 连续小割次数 | [build_sequential_region.small_cuts](../main.py) | 连续三次 area_ratio<0.02 时换阶段；零面积分母不使用该判据 |
| 近点合并距离 εp | [build_sequential_region.point_tol](../main.py)、CLI `--point-tol` | 默认 1e-4 kW；以原始负荷坐标的 L-infinity 距离比较，不是标幺值或物理可行性容差 |
| 点的固定代表及编号 | [register_power](../main.py)、[build_sequential_region.powers](../main.py) | powers[index] 保存首次登记的原始 kW 坐标；近点复用已有编号，代表不移动、不传递合并，不舍入求解坐标 |
| SP 缓存索引 | [build_sequential_region.cache](../main.py) | (tuple(x),index)；历史 frontier / certifications 同代表索引已退役；输出 p/x 的单位、形状与含义不变 |
| 已采用的局部割来源 | [build_sequential_region.applied](../main.py) | (tuple(x),index) 集合；近点匹配不得让同一缓存割再次被采用并计入连续小割次数 |
| 换网架候选队列 | 已退役：`experiments/sequential_region.py::build_sequential_region.frontier` | 历史定义保留：按固定代表编号缓存来源 x 与 eta；新流程由连续查漏直接返回可行网架和负荷 |
| 全网架已证不可行点 | 已退役：`build_sequential_region.rejected_points`（2026-09-29 主线精简） | p 为 kW；记录 eta 与可靠下界 bound，回放显示红色实心点；不删除邻域 |
| SP 数值精度侧重 | [SubProblem.numeric_focus](../model.py) | 求解器 NumericFocus；默认 0 保持主线；实验 FourBus 固定 0、Case33 固定 3，从首次 SP 起使用；不改变原始 eta 与 1e-8 接受条件 |
| 射线数值精度侧重 | [ray_support.numeric_focus](../main.py) | 与该案例 SP 的 NumericFocus 一致；默认 0；不改变物理约束和精度 |
| 每阶段补充射线次数 | 已退役：`build_sequential_region.boundary_rays`（2026-09-29 主线精简） | 整数；stages[].boundary_rays；counts.ray 仍累计全部射线求解 |
| 网架未认证面积 | 已退役：`build_sequential_region.local_gap_area`（2026-09-29 主线精简） | area(N_x 差 N'_x)，kW²；stages[].local_gap_area，不改变 remaining_area=area(N_x 差 G') |
| 射线新增面积上界 U(q) | [ray_gain](../main.py) | area(conv(N'_x∪{q}))−area(N'_x)，内部用公共箱归一化面积；回放 ray_gain_bound 乘 prod(bounds) 后为 kW² |
| 射线收益筛选比例 | [RAY_THRESHOLD](../main.py)、[build_sequential_region.ray_threshold](../main.py)；历史 CLI `--ray-threshold` 已随实验入口退役 | 默认 1e-4，即 0.01%；U(q)<ray_threshold*μ(N'_x) 时省去该射线；首轮 G 顶点射线不筛选；低维内域先补足维度，不按零测度筛除；不是停止证书 |
| 非冗余认证域半空间 | [coverage_halfspaces](../main.py) | 返回保留的网架键及其半空间；只删被某一个其他认证凸包包含的排除项，不跨网架取凸包；使用已有 GEOMETRY_TOL |
| 全局查漏逐次记录 | 已退役：`build_sequential_region.coverage_checks`（2026-09-29 主线精简） | 每项记录 total_domains/domains、total_faces/faces、objective、bound、complete；objective 为求解器当次可行解的未覆盖余量，bound 仍是全局上界 |
| 求解与筛选计数 | [build_sequential_region.counts](../main.py) | 新增 initial=初始化 MP2 次数、cut_lp=实际割 LP 次数、ray_skipped=按几何上界省去的射线数；sp/cuts/ray/global_search 原义不变，sp_seconds 仍包含评分和取割时间 |

2026-09-29 效率迁移：SP 增加显式 `score_only=True`，只求原始最小 eta 并缓存锥支撑方向；默认 False 保持所有原调用方“可行证书或有效割”的契约。仅获选的最大 eta 候选调用 generate_cut，复用缓存方向建立支撑 LP，不重复求 SOCP。切割后和阶段末射线均先计算 U(q)；补充轮按最新 U(q) 从大到小处理，每收入一个点后重新计算收益。小收益点不获得认证，最终仍由完整物理全局模型查漏。

包含判定只减少全局查漏的排除项，所有物理网架 x 仍开放、所有 N'_x 仍保存在结果和回放。按已有几何容差删除较小排除项可能扩大搜索域，不会缩小搜索域，因此不会把未覆盖的可行点排除在停止证明之外。相同认证域按面积排序的稳定顺序保留一个。新增 ray_skip、cut_start 事件；回放保留几何收益上界及实际全局目标值、上界和排除规模，version=4 不变。

局部选点使用未收缩的条件外域。每批先收入可行点，再选原始 SP eta 最大的不可行候选采用一条共享割，并沿该方向求固定网架的完整 SOCP 射线支撑。局部阶段结束前，遍历当前 N_x 的剩余未认证顶点，按 1e-4 kW 固定代表编号在本轮内去重，以最新 N'_x 认证凸包重心为锚求射线边界点；入内域的是实际求出的 p，不是近邻或外顶点。之后立即调用主线 physical 剩余域搜索，tau=0.005、GEOMETRY_TOL 不变；返回见证直接加入其网架内域，也允许返回已知网架；只有覆盖证书才允许 certified=True。每案例的 1000 秒预算包含初始化、全部求解、几何与轨迹记录；超时保存已有内域和有效初始外界，不重试、不放宽精度、不自动启动扫描。

固定 x 的原 SP 最优 eta 是 p 的凸函数，外包络极点足以定位最大违反量；共线中间点与内部点不会给出更大值。每次实际 SP 前重新检查当前同网架认证凸包，两个端点刚被认证后，线段内部点即可跳过；不同网架的认证点不得合并为一个凸包。近点只复用固定代表的 SP 结果，不把邻近点赋予相同物理证书。若只剩近点或已采用割的代表，本阶段记为 point_resolution，补边界后进行连续查漏，不能凭近点合并宣称全局覆盖。回放 settings 与 result 的 point_tol 保持原义。

回放仍为 version=4，保留 active_scheme、stage、area_ratio、small_cuts、inner_area、ray、rejected_points、unknown_points；新增 boundary_start / boundary_end 事件标记补充轮，phase='边界补充' 的 ray_start / ray_end 展示同网架射线。旧文件语义不变。每案例仅存一份原生 monitor 压缩文件，支持原有 monitor.py、实验入口 --replay 及双窗口回放。参考扫描只能显式指定已有同参数回放，不计入构域，不影响选点。

### 二维 / 三维入口与参考扫描迁移（2026-09-29）

2026-09-29 快速调试迁移：逐网架实验的 build_sequential_region.seconds、run_case.seconds 和 CLI --seconds 默认统一为 50 秒；主入口 CASE_TIME_LIMIT 同为 50 秒，两个案例都读取该设置。完整预算包含初始化、几何、求解和过程记录。旧结果及其 1000 秒设置保留。新增 [build_sequential_region.initialized](../main.py) 为已完成首轮 G 顶点射线认证的网架键集合，每个网架只执行一次；[build_sequential_region.add_ray.initial_sweep](../main.py) 标记首轮，跳过收益筛选但不跳过物理认证。首轮遍历全局初始外包络 G 的顶点，已在同网架认证凸包内的点复用凸性，其余逐个求射线支撑；首轮不产生对偶割，随后恢复 N_x 顶点的最大 eta 选点切割。新增 initial_sweep_start / initial_sweep_end 回放事件，原 ray / sp / cuts 计数含义不变。

三维体积修正：暴露面求差统一使用 1e-11 的二维投影网格精度，与 _clip_face 的面内距离容差一致，避免约 1e-18 的共线偏差让布尔求差漏扣重叠面。仅规范测度计算的中间面片，不移动认证点、不改变物理或覆盖容差、不对体积作非负截断或累计最大值修补。

Case33 的动态负荷节点由文件顶部 LOAD_NODES 配置，默认 (18,25)，设为 (18,25,30) 切换三维；未选中节点始终保持原始负荷。主入口和逐网架实验均读取该配置，也支持 --load-nodes 18,25,30 覆盖。维数 d=len(load_nodes) 同时控制初始化方向数、几何、联合割切片、扫描张量与原生回放投影。

测度显式迁移：二维 μ 为面积，三维 μ 为体积。历史字段 inner_area、removed_area、area_denominator、remaining_area、local_gap_area、ray_gain_bound、ray_gain_threshold 在二维仍为 kW²，在三维为 kW³；记录增加 measure_unit 标明单位。area_ratio 仍是同维测度之比；stages 内同名字段沿用该规则，不增加同义结果接口。ray_gain 改为凸包的 d 维测度增量，三维补点优先级和小收益阈值均按体积计算。跨网架始终取并集，三维并集按暴露表面积分测度，不取跨网架凸包。

| 新增 / 迁移量 | 固定代码映射 | 定义 |
|---|---|---|
| Case33 默认动态负荷节点 | [LOAD_NODES](../Network/case33bw.py) | 节点 ID 元组；二维 (18,25)，三维 (18,25,30) |
| Case33 每轴参考扫描格数 | [SCAN_DIVISIONS](../main.py) | 按维数索引：2→100，3→60；均为均匀网格中心，可由 --divisions 覆盖 |
| d 维凸域并集测度 | [union_measure](../main.py) | 归一化坐标下二维面积或三维体积；二维用多边形布尔并，三维用既有 union_volume |
| 输出测度单位 | 已退役：`build_sequential_region:measure_unit`（2026-09-29 主线精简） | kW² 或 kW³；回放 settings 同步保存 |
| 按节点配置定位回放 | [recording_path](../main.py) | 文件名 case_节点序列.json.gz，例如 case33_18_25_30.json.gz；二维和三维不相互覆盖，旧文件不改写 |
| 扫描张量与耗时 | [validate_socp_region:states](../vertify.py)、[validate_socp_region:scan_seconds](../vertify.py) | states.shape=(divisions,)*d，索引顺序与 load_nodes 一致，1/-1 为完整 SOCP 可行/不可行；扫描秒数独立于构域计时 |
| 扫描统计计数 | [RunMonitor.validation:missed_cells](../monitor.py)、[RunMonitor.validation:extra_cells](../monitor.py)、[RunMonitor.validation:reference_cells](../monitor.py)、[RunMonitor.validation:computed_cells](../monitor.py) | 各自 inner/outer 指标下保存，MR=100*missed/reference，FR=100*extra/computed；原 mr_percent/fr_percent 不变 |
| 三维割面截多边形 | [_cut_polygon](../monitor.py) | 将联合割代入指定 x，与显示盒的十二条棱求交；仅绘图，坐标 kW |
| 三维凸域绘图 | [_draw_3d](../monitor.py) | 分网架绘制真实凸域表面及低维认证集；不跨网架构造凸包 |
| 扫描体素边界 | [_voxel_faces](../monitor.py) | 按真实 states 的六邻接边界生成体素面，坐标 kW，不对参考网格作平滑或插值 |

新实验默认构域后执行扫描，--no-scan 仅计算并保存轨迹；--reference 仍显式复用同节点、同预算的既有参考。每个扫描行复用一个完整 MP 的结构，只修改固定负荷等式右端；每个网格中心都重新优化，所有合法 x/y 始终自由。不使用算法割、已知网架列表或算法认证域过滤参考点。原生三维回放按当时历史绘制 G、G'、N_x、N'_x、切割平面和射线；最终参考与 inner/outer 的两套 MR/FR 保持独立，不随时间轴倒退。

此处“真值对比”统一标记为同一 SOCP 模型的独立网格参考，不是非凸 AC 真值；0% 网格误差不替代连续域覆盖证书。

Case33 数值设置：numeric_focus=3 的 SP 在首次求解前设置 BarHomogeneous=1，用齐次内点算法处理极限负荷处的退化约束；其余仍使用求解器原自动设置 -1。原 eta、MaxVio、二者共用的 1e-8 认证误差预算、默认 BarQCPConvTol=1e-9 和 LP 联合割均保持不变。不存在逐点失败后重试或参数切换。

参考扫描的 [SCAN_WORKERS](../main.py)（main.py 顶部，默认 12）控制独立进程数；validate_socp_region.workers 默认 1 保持原调用口径，_scan_line 对一个固定前 d−1 维索引的扫描行复用完整模型，返回该行的逐点状态。并行时每个模型仅用一个求解线程；单进程扫描使用 SOLVER_THREADS；进度与结果在主进程合并。

三维并集测度迁移：union_volume 用散度定理 μ(U)=Σ_f area(f_exposed)*(n_f·p_f)/3。将每个原始凸域的边界三角面投影到自身二维坐标，扣除被其他凸域覆盖的面片；同向共面重叠表面按网架顺序只计一次，内部面不计。原始凸域和认证半空间不变，避免反复三维相减产生退化薄片。面片裁剪使用原几何容差；没有随机扰动或失败重试。

面内有序多边形由 [_clip_face](../plot.py) 按 Sutherland–Hodgman 半平面裁剪；面上判据与既有裁剪相同，为单位法向距离 1e-11。它保留环上顺序，因此每次切面不需要重新求凸包；仅供并集测度，不生成物理认证点或全局割。

独立扫描从首次求解起固定 NumericFocus=3、BarHomogeneous=1、Aggregate=0、ScaleFlag=0；完整模型的可行性检查仍为 MaxVio<=1e-8。扫描失败直接停止，不将超限残差归为可行、不作逐点参数重试。

## 独立扩展：AC 电流等式试验（2026-09-28）

`experiments/signed_ac_equality.py` 复用带符号实验的 x、p、state、预算和全部参数，只将支路 `current_cone` 替换为 P²+Q²=v_i*ell，使用 NonConvex=2 直接求解完整混合整数非凸模型。`method='ac'` 仅标记这个独立试验；不迁移生产 `GridPhysics.method` 的取值。没有调用 SOCP 的 KKT、对偶割或认证凸包。

方向试验最大化 direction@p；定点试验固定 p、允许全部合法 x/y 调节并最小化既有费用。两种模型使用相同研究盒、预算、功率因数和电压界。`current_gap_max=max(abs(v_i*ell-P²-Q²))` 为标幺电流关系的最大绝对残差，`voltage_min/voltage_max` 为幅值 p.u.，`loss_kw=sum(r*ell)*base` 为有功损耗。方向上界 `bound` 与原契约相同，单位 kW；定点问题目标和界使用原费用单位。计时包含建模与求解；超时单独标记，不能称为不可行或最优。

输出每案例一份 comparison.json，只记录八方向支撑点和定点检验；支撑点之间不填充凸包，不据此宣称已构建完整 AC 域。

`experiments/plot_signed_ac_equality.py` 只读取这些结果：圆点/方点表示各模型实际求出的支撑点，虚线表示八条 `direction@p<=bound` 的交集，是外包络，不认证内部。Case33 绘图坐标除以 1000 显示为 MW，持久化数据仍为 kW；重合支撑点允许重叠显示。定点检验单独标记 AC 可行/不可行。

## 独立扩展：带符号净功率实验（2026-09-28）

`experiments/signed_power_compare.py` 只在实验实例中扩展 p 的允许范围；原主线、Network 默认参数、旧记录均不迁移。p / power / 返回键 p 仍以 kW 表示节点有功需求，负值为净注入；采用原 q_ratio 的双向固定功率因数关系。x、state、联合割布局不变。

新增实验量：`power_lower/power_upper` 为研究盒下/上限；`axis_lower/axis_bounds` 为四个轴向完整 MP2 的可靠下/上界；`bounds=axis_bounds-axis_lower` 为本实验几何跨度；归一化 `point=(power-axis_lower)/bounds`。转换割时常数为 alpha+beta@axis_lower+delta@x，线性项为 beta*bounds。所有持久化 p、inner/outer.vertices 和 cuts 均为原始 kW 坐标，不保存平移后的伪负荷。

初始化求八个方向的完整 MP2，每个方法独立计时；正负两个坐标始终同时自由。方向上界形成八条全局支撑约束，和 SP 割分开计数。研究盒默认各轴 ±min(source_pmax,source_smax)*base，这是有限研究范围，不是新增线路额定容量。首端 P/Q 允许按原进口上限对称反送，原电压范围不变。

运行有效界：总可用有功为首端进口上限加研究盒中全部节点最大净发电；Q 同理。P/Q 支路界由这两个总可用量及原 active-power capacity 得到。总有功损耗上界为首端进口上限减最小总净需求；ell 上界为 min(总损耗/r,2*pmax/r)。所有参考方向包含首端都允许反向，并限制反向送端功率。Case33 rateA=0 仍表示未提供线路额定限额；不得称上述派生界为额定电流。

主线实验保留顶点 SP 最大违反量优先、32 次新增 SP 后全局检查、同网架见证支撑及完整物理全局搜索。仅将原径向 tau 替换为 `padding=0.005` 的归一化 L-infinity 距离扩边：排除 I_x+[-padding,padding]^2，扩边点不认证。全局覆盖成立后这个扩边并集为有效外包络；最终 outer 同时与全部网架联合割投影相交。该精度与旧径向 tau 数值相同但数学含义不同。

KKT 实验使用原完整 SP 锥对偶与负荷 LP 互补模型；下界行 -p_j<=-axis_lower_j、平移后的并集排除行均登记进原约束、驻点、互补、目标。`coverage_pad=1e-5`；eta 零潮流上界在整个正负盒上计算。停止仍为未覆盖配对上界及认证扩边残差界的最大值 <=epsilon，默认 epsilon=0.01，不能称全配对 R_k 已小于阈值。

两方法默认各最多 1000 秒、4 求解线程、数值库 1 线程。同案例共享事后 80×80 网格中心扫描，每点完整 SOCP 自由选择全部合法 x/y，不用实验割、认证点、正象限单调性或外域筛选。扫描为 SOCP 离散参考，不称 AC 真值。MR=集合外可行点/全部可行点；FR=集合内不可行点/集合内点，各乘100；inner/outer 分开报告。精简轨迹、参数、参考和指标仅存每案例一份 comparison.json.gz；绘图、扫描和事后枚举不计入构域时间。

SP 最终直接复用原 `SubProblem`，求解器参数、原始证书 eta+MaxVio<=1e-8 和 `_separating_cut` 均不变。主线候选顶点朝同网架认证凸包重心移动至多 padding/4（空内域时使用完整物理 G 的可行见证为锚）；实际检查与记录的 p 为移动后的点。全局覆盖与跳点仍检查未移动的外域，不能由内移本身宣称外顶点可行。若新点认证，则原顶点距新证书至多 padding/4；若新点不可行，有效割也排除原顶点，因为其与认证锚的凸组合已被分离。该规则替代旧主线向原点收缩，避免无意义的边界数值重复，不是求解失败后的重试。

版本：2.2，2026-09-28。顶点 SP 违反量优先、32 次新增 SP 的全局检查间隔；当前入口为二维 Case33 + SOCP + Python 原生监视窗口，径向精度 tau=0.005。Case33 采用第 39 节的七开关及开合变动预算；联合割布局及公共归一化坐标保持不变。

本文是本项目数学符号、代码名称、单位、数组顺序和结果字段的统一约定。正式代码为 `Network/`、`model.py`、`region.py`、`vertify.py`、`plot.py`、`main.py`、`survey.py`、`monitor.py`；测试、注释和新文档使用同一约定。项目外归档的源码、已有结果和固定方案历史推导保留原口径，其局部记号须通过附录（原第 10 节）换算。

**同一个量的含义、单位和索引没有改变，就保持原代码名。** 性能优化、拆函数、改求解器、整理代码都不是重新命名数学量的理由。本文登记当前名称，不要求为追求字面一致再做一次全库改名。新增量和迁移按第 11 节办理。

## 1. 最常用的对应关系

下表中的实例名 `network`、`equations`、`problem`、`tree` 表示相应类对象，不强制改掉现有局部对象别名。数学量的固定字段、形参、字典键区分大小写。

| 数学符号 | 含义 | 固定代码名称 | 单位 / 形状 |
|---|---|---|---|
| \(x_{e,k}\)；向量 \(x\) | 走廊及型号选取 | 建模 `x[e,k]` / `problem.choices`；数组 `problem.x`、`answer['x']` | 二进制；扁平视图 `(t,)` |
| \(z=Hx\) | 走廊接通状态 | `z`，由 `network.corridor_types @ x` 派生 | `(m,)` |
| \(p\) | 独立可调有功负荷坐标 | 接口 `power`、`problem.power`；公式局部 `p`；返回键 `'p'` | kW，`(d,)` |
| \(y\) | 完整运行向量 | 接口 / 字段 `state`；装配局部 `y`；返回键 `'state'` | `(3*t+n+2*m,)` |
| \(P,Q,\ell,v\) | 支路有功、无功、电流平方、节点电压平方 | `P`、`Q`、`ell`、`v`，由规定切片访问 | p.u. 或 p.u.² |
| \(r,\chi\) | 线路电阻、电抗 | `r`、`reactance` | p.u.，`(t,)` |
| \(c,\mathcal B\) | 预算线性系数、预算 | `network.cost`、`budget` | `network.cost_unit`；Case33 总变动数含常数 cost_offset，见第 39 节 |
| \(S_{\mathrm b}\) | 功率基值 | `network.base` | kVA |
| \(\lambda\) | 必要线性约束的乘子 | `dual` | 支撑平面 LP 的约束行顺序，见第 6 节 |
| \((\alpha,\beta,\delta)\) | 联合割系数 | `cut = [alpha, *beta, *delta]` 的扁平数组 | `(1+d+t,)` |
| \(\xi=p\oslash b^{\rm box}\) | 评价箱归一化坐标 | 几何接口中的 `point`、`points`、`vertices`，见第 7 节的作用域 | 无量纲，二维或三维 |
| \(\tau\) | 径向收缩比例 | `tau`、`REGION_TAU` | 无量纲 |

禁止把电抗称为 `x`、把电流平方写成 `l`、把预算记作标准式矩阵 `B`、把归一化坐标记作走廊状态 `z`。源数据里的 `x_ohm_km` 等字段按第 3 节处理。

## 2. 集合、维数和唯一索引

数学下标从 1 起写，Python 数组从 0 起；真实节点 ID 由数据给出，不能直接当数组位置。根节点在公式中记为 0，内部哨兵为 `-1`，不占 `nodes`、`v`、`a` 的元素。

| 符号 | 定义 / 顺序 | 代码定义 |
|---|---|---|
| \(\mathcal N,n\)；\(i,j\) | 全部候选非根节点及其索引 | [Network.nodes](../Network/__init__.py)、[Network.n](../Network/__init__.py)、[Network.node_index](../Network/__init__.py)、[Network.root](../Network/__init__.py) |
| \(\mathcal E,m\)；\(e\) | 候选走廊，按输入顺序 | [Network.corridors](../Network/__init__.py)、[Network.n_corridors](../Network/__init__.py) |
| \(\mathcal T,t\)；\(\nu\) | 各走廊内型号连续展开；\(\nu\) 是扁平数组位置，总数 t | [Network.type_keys](../Network/__init__.py)、[Network.n_types](../Network/__init__.py) |
| \((e,k)\) | Gurobi 建模键：走廊 ID 与该走廊内型号 ID；k 不要求全图唯一 | [GridPhysics.keys](../model.py)、[GridPhysics.types](../model.py) |
| \(\mathcal T_e\) | 走廊 e 的型号区间；型号到走廊的映射 | [Network.type_slices](../Network/__init__.py)、[Network.type_corridor](../Network/__init__.py) |
| \(\mathcal L,d\)；\(h\) | 独立负荷节点，`d = len(network.load_nodes)` | [Network.load_nodes](../Network/__init__.py)、[Network.selected](../Network/__init__.py) |
| \(\mathcal N_{\rm req}\) | 必须接入节点的布尔掩码，`(n,)` | [Network.required](../Network/__init__.py) |
| \(n_y=3t+n+2m\) | 标准运行向量长度 | [GridPhysics.slack_slice](../model.py) 的 `.stop` |
| \(n_{\rm row}\) | 割生成 LP 的线性约束数，随支撑平面变化 | [SubProblem._separating_cut.rows](../model.py) |
| \(n_{\rm tree}\) | 当前方案实际接入的非根节点数，也是树的支路数 | [OperatingTree.n](../Network/__init__.py) |

`network.type_keys[nu] == (e,k)` 是具名键与扁平位置的唯一桥梁。求解器内 `p[i]` 用真实负荷节点 ID；数组 `power[h]` 用 `load_nodes[h]` 的位置。`n/m/t/d` 在数学文档中固定为上表维数；`m = model`、`e = equations` 等既有局部对象别名不重定义这些数学符号。

模型装配可处理 `d` 个负荷坐标；正式 `region.py`、剩余域几何和实时查看器支持 `d=2` 或 `3`，AC 网格输出仍固定三维。二维勘察由 `survey.py` 计算，`plot.py` 绘图。

## 3. 网架参数、量纲和数据入口

| 符号 / 含义 | 代码定义 | 单位 / 形状 |
|---|---|---|
| \(S_{\mathrm b}\)，功率基值 | [Network.base](../Network/__init__.py) | kVA，标量；不是自动等于电源容量 |
| \(V_{\mathrm b}\)，线电压基值 | [Network.voltage_kv](../Network/__init__.py) | kV，标量 |
| \(r_{e,k}\)，数组 \(r_\nu\) | [TypeParameters.r](../Network/__init__.py) → [Network.r](../Network/__init__.py)；具名 [GridPhysics.r](../model.py) | p.u.，标量 → `(t,)` / 按键字典 |
| \(\chi_{e,k}\)，数组 \(\chi_\nu\) | [TypeParameters.reactance](../Network/__init__.py) → [Network.reactance](../Network/__init__.py)；具名 [GridPhysics.reactance](../model.py) | p.u.，标量 → `(t,)` / 按键字典 |
| \(P_{e,k}^{\max}\)，输入的送端有功上限 | [TypeParameters.capacity](../Network/__init__.py) → [Network.capacity](../Network/__init__.py) | p.u.，标量 → `(t,)`；不是电流或视在功率 |
| \(c_{e,k}\)，预算线性系数 | [TypeParameters.investment_cost](../Network/__init__.py) → [Network.cost](../Network/__init__.py)；具名 [GridPhysics.cost](../model.py) | `cost_unit`，标量 → `(t,)` / 按键字典；通常为增量投资，Case33 的仿射开合次数见第 39 节 |
| 费用单位、默认预算列表 | [Network.cost_unit](../Network/__init__.py)、[Network.budgets](../Network/__init__.py) | FourBus / 江口为元；当前 Case33 为开合变动次数，历史升级测试夹具仍为相对投资单位 |
| \(p^{\rm original},q^{\rm original}\) | [Network.original_p](../Network/__init__.py)、[Network.original_q](../Network/__init__.py) | kW / kvar，`(n,)` |
| \(p^{\rm fixed},q^{\rm fixed}\) | [Network.fixed_p](../Network/__init__.py)、[Network.fixed_q](../Network/__init__.py) | kW / kvar，`(n,)`；独立负荷节点置零 |
| \(\kappa_h=q_h/p_h\) | [Network.q_ratio](../Network/__init__.py) | 无量纲，`(d,)`；并非功率因数本身 |
| \(\underline v,\overline v\) | [Network.vmin](../Network/__init__.py)、[Network.vmax](../Network/__init__.py) | p.u.²，`(n,)`；已平方 |
| \(P_0^{\max},Q_0^{\max},S_0^{\max}\) | [Network.source_pmax](../Network/__init__.py)、[Network.source_qmax](../Network/__init__.py)、[Network.source_smax](../Network/__init__.py) | p.u.，标量；源端含线路损耗 |
| \(p_\Sigma^{\max}\)，独立负荷总量上界 | [Network.power_limit](../Network/__init__.py) | kW，标量；不与源端容量混用 |

统一换算：

\[
Z_{\mathrm b}=V_{\mathrm b}^2/(S_{\mathrm b}/1000),\quad
d^P(p)=(p^{\rm fixed}+Ep)/S_{\mathrm b},\quad
d^Q(p)=(q^{\rm fixed}+E(\kappa\odot p))/S_{\mathrm b}.
\]

[Network.loads](../Network/__init__.py) 输入 `power` 的单位始终为 kW；输出两组节点标幺负荷，形状均为 `(batch, n)`。单个点仍保留批次轴 `(1,n)`。树的 `loads` 取其中接入节点，列数改为 `n_tree`。

允许输入适配层保留原始物理名称：[LineType.r_ohm_km](../Network/four_bus_five_corridor.py)、[LineType.x_ohm_km](../Network/four_bus_five_corridor.py)、[LineType.capacity_kw](../Network/four_bus_five_corridor.py)、[LineType.cable_cny_m](../Network/four_bus_five_corridor.py)。江口 JSON 的 `x_ohm_per_km`、`max_i_ka` 在 [load_jiangkou](../Network/jiangkou.py) 转成统一参数；[Case33](../Network/case33bw.py) 在入口完成 MATPOWER 单位换算。禁止把这些原始字段直接当 `Network` 中的标幺量。

走廊结构固定为 [Corridor.id](../Network/__init__.py)、[Corridor.endpoints](../Network/__init__.py)、[Corridor.existing_type](../Network/__init__.py)、[Corridor.initial_active](../Network/__init__.py)、[Corridor.types](../Network/__init__.py)；型号标识为 [TypeParameters.id](../Network/__init__.py)。初始开合状态不等于规划决策。

## 4. 拓扑、关联矩阵与运行树

| 符号 | 固定映射 | 形状 / 约定 |
|---|---|---|
| \(H\) | [Network.corridor_types](../Network/__init__.py) | `(m,t)`；`z = H @ x` |
| \(R_{\rm recv},R_{\rm send}\) | [Network.receiving](../Network/__init__.py)、[Network.sending](../Network/__init__.py) | `(n,m)`；根行删除 |
| \(G=R_{\rm recv}-R_{\rm send}\) | [Network.incidence](../Network/__init__.py) | `(n,m)`；流入为正 |
| \(E\) | [Network.E](../Network/__init__.py) | `(n,d)`；独立负荷嵌入全部节点 |
| \(T=R_{\rm recv}H\)、\(J=(R_{\rm send}H)^T\) | 数学推导量；当前具名建模不存储 T/J | `(n,t)` / `(t,n)`；受端关联 / 送端电压提取 |
| \(\rho_\nu\) | 数学推导量，对应参考送端是否为根 | `(t,)`；根出线集合由 [GridPhysics.outgoing](../model.py) 给出 |
| 参考端点、节点入边 / 出边 | [GridPhysics.ends](../model.py)、[GridPhysics.incoming](../model.py)、[GridPhysics.outgoing](../model.py) | 按走廊 / 节点 ID 的字典 |
| \(a\) | [MasterProblem.active_nodes](../model.py)，局部 `a` | `(n,)`，二进制；根恒接入 |
| \(f\) | [MasterProblem.__init__.f](../model.py) | `(m,)`；连通虚拟流，不是电功率 |
| \(x\) | [MasterProblem.x](../model.py) | `(t,)`，完整型号顺序；无第二套压缩选型 |

参考端点数组为 [Network.senders](../Network/__init__.py)、[Network.receivers](../Network/__init__.py)，形状 `(m,)`。接根走廊参考方向朝外，其余按输入端点顺序。拓扑约束统一记为

\[
z=Hx,\quad z_e\le1,\quad z_e\le a_i,a_j,\quad
Gf=a,\quad -nz\le f\le nz,\quad \mathbf1^Tz=\mathbf1^Ta.
\]

外部完整方案 `plan` / `choice` 是“走廊 ID → 型号 ID 或 `None`”的字典；`fixed_plan` 也是此格式。通过 [Network.encode_plan](../Network/__init__.py)、[Network.decode_plan](../Network/__init__.py) 与 `x` 转换。`start` 是型号向量，`incumbent` 是含证书的答案字典，三者不能互换。

选定运行树由 [Network.tree](../Network/__init__.py) 构造：

| 树量 | 代码定义 | 解释 |
|---|---|---|
| 全图映射 | [OperatingTree.node_indices](../Network/__init__.py)、[OperatingTree.type_indices](../Network/__init__.py) | `(n_tree,)`；每个树节点与其入边共用局部索引 |
| 父子、遍历 | [OperatingTree.parent](../Network/__init__.py)、[OperatingTree.order](../Network/__init__.py)、[OperatingTree.children](../Network/__init__.py)、[OperatingTree.roots](../Network/__init__.py) | 树局部索引；`roots` 是接根支路索引，不是根节点 ID |
| \(\sigma\) | [OperatingTree.direction](../Network/__init__.py) | `(n_tree,)`，`+1/-1` 表示根向与参考方向同向 / 反向 |
| \(D\) | [OperatingTree.D](../Network/__init__.py) | `(n_tree,n_tree)`；行入边、列下游节点，含受端自身 |

树上的功率记 \(P^{\rm tree},Q^{\rm tree}\)，候选图上的功率为 \(P,Q\)。反向时必须计损耗：

\[
P_\nu=-P^{\rm tree}_j+r_\nu\ell_\nu,\qquad
Q_\nu=-Q^{\rm tree}_j+\chi_\nu\ell_\nu\quad(\sigma_j=-1).
\]

## 5. 具名运行变量、扁平向量和数学标准式

[GridPhysics](../model.py) 固定状态布局为

\[
y=(P,Q,\ell,v,s^+,s^-),\qquad n_y=3t+n+2m.
\]

建模使用 [GridPhysics.add_operation](../model.py) 中的 `P[e,k]`、`Q[e,k]`、`ell[e,k]`、`v[i]`、`plus[e]`、`minus[e]`。它们返回为 `problem.operation` 的同名字段；`state` 是同一批变量的扁平视图，没有第二批运行变量。

| 具名对象 | 固定代码定义 | 与向量的关系 |
|---|---|---|
| 建设选型、负荷 | [MasterProblem.choices](../model.py)、[MasterProblem.loads](../model.py) | `[e,k]` / `[i]`；[MasterProblem.x](../model.py)、[MasterProblem.power](../model.py) 是共享变量的 MVar 视图 |
| 有功、无功、电流平方、电压平方 | [GridPhysics.add_operation.P](../model.py)、[GridPhysics.add_operation.Q](../model.py)、[GridPhysics.add_operation.ell](../model.py)、[GridPhysics.add_operation.v](../model.py) | 根电压 `v[network.root] = 1` 是常数，不进入状态向量 |
| 正、负开断压降余量 | [GridPhysics.add_operation.plus](../model.py)、[GridPhysics.add_operation.minus](../model.py) | 求解器名称 `drop_plus/drop_minus`；数学符号 \(s^+/s^-\) |
| 完整运行对象、状态 | [MasterProblem.operation](../model.py)、[MasterProblem.state](../model.py) | `cuts_only=True` 时均为 `None` |

不要因采用 `tupledict`、`MVar` 或 NumPy 就改数学量名称；转换严格使用 `type_keys`、`nodes`、`load_nodes`、`corridors` 顺序。

| 分块 | 唯一切片 | Python 范围 | 含义 |
|---|---|---|---|
| \(P\) | [GridPhysics.P_slice](../model.py) | `[0:t]` | 参考送端有功，p.u. |
| \(Q\) | [GridPhysics.Q_slice](../model.py) | `[t:2*t]` | 参考送端无功，p.u. |
| \(\ell\) | [GridPhysics.ell_slice](../model.py) | `[2*t:3*t]` | 电流幅值平方，p.u.² |
| \(v\) | [GridPhysics.v_slice](../model.py) | `[3*t:3*t+n]` | 非根电压幅值平方，p.u.² |
| \(s^+,s^-\) | [GridPhysics.slack_slice](../model.py) | `[3*t+n:3*t+n+2*m]` | 前 m 为正、后 m 为负；开断压降余量，p.u.² |

读取使用 `state[equations.P_slice]` 等切片，不在调用方另写一套偏移。接根型号 `P/Q` 下界为零；非根型号允许带符号。未选型号的 `P/Q/ell` 归零；linear 固定全部 `ell=0`。未接入节点的电压是辅助量，不能作为运行电压导出。

节点平衡与送端电压平方固定写成

\[
(T-J^T)P-T(r\odot\ell)=d^P(p),\quad
(T-J^T)Q-T(\chi\odot\ell)=d^Q(p),\quad u=\rho+Jv.
\]

SOCP 为 \((u_\nu+\ell_\nu,2P_\nu,2Q_\nu,u_\nu-\ell_\nu)\in\mathcal Q_4\)；AC 使用 \(P_\nu^2+Q_\nu^2=u_\nu\ell_\nu\)。源端功率为 \(P_0=\rho^TP,Q_0=\rho^TQ\)，具名模型使用 `source_p/source_q` 对根出线求和。

需要整体矩阵推导时，标准式采用以下符号和符号方向：

\[
w=b-Ax-By-Cp\in\mathcal K,\qquad
l_0+Lx\le y\le u_0+Ux.
\]

这些是数学块名，**当前实现没有 `equations.A/B/C/b`、`lb_x/ub_x` 等属性**；不得为了满足公式字面而重建一套平行模型。当前物理约束由 `add_operation` 的具名约束给出，割计算才从求解器导出线性行系数。

| 数学量 | 代码定义 | 形状 / 说明 |
|---|---|---|
| \(l^{\rm glob},u^{\rm glob}\)，全局变量盒 | [GridPhysics.y_lb_global](../model.py)、[GridPhysics.y_ub_global](../model.py) | `(n_y,)`；保留全部型号，与 x 无关 |
| \(\overline P,\overline Q,\overline\ell\) | [GridPhysics.pmax](../model.py)、[GridPhysics.qmax](../model.py)、[GridPhysics.ellmax](../model.py) | 按 `(e,k)` 的派生有效界 |
| \(\underline P,\underline Q\) | [GridPhysics.pmin](../model.py)、[GridPhysics.qmin](../model.py) | 根出线为 0，其余为负的有效上界 |
| \(\underline v,\overline v\) 的具名视图 | [GridPhysics.vmin](../model.py)、[GridPhysics.vmax](../model.py) | 节点 ID 字典，含根节点常数 1 |
| \(M_e\)，开断压降界 | [GridPhysics.drop_max](../model.py) | 走廊 ID 字典；\(0\le s_e^\pm\le M_e(1-z_e)\) |
| 求解结果最大违反量 | Gurobi `model.MaxVio` | 直接读取所建模型的数值质量；MP 接受条件为 `<= PLANNING_TOL` |
| 测试用原约束最小余量 | [margin](../tests/planning_checks.py) | 保留独立代入检查，仅用于测试和审核，不参与生产求解 |

全局盒不是完整可行域：其中电压下界为 0，实际 `v >= vmin` 另作运行约束。各约束有自己的物理尺度，不能给整个 `y` 或 `w` 一概标注 kW。

派生界固定为 \(\overline P=\min(P^{\max},P_0^{\max},S_0^{\max})\)、\(\overline Q=\min(Q_0^{\max},S_0^{\max})\)、\(\overline\ell=\overline P/r\)（linear 为 0）；不能与原始输入 `capacity` 混称。

## 6. SP、对偶乘子和联合割

固定 \((\hat x,\hat p)\) 后，SP 最小化 \(\eta\)、\(\eta\ge0\)，当前 [SubProblem.solve.eta](../model.py) 的求解器名称为 `violation`。仅功率平衡与压降等式允许 ±eta；选型容量、电压界、开断余量、电源限额和二阶锥保持严格。`state` 不包含 eta。此辅助问题替代旧式“锥首分量也松弛”，不直接比较两者目标或乘子数值；eta=0 对应的物理约束、输出字段和数组布局不变。

SP 先求上述带 eta 的问题；求解必须正常结束，且 `MaxVio <= PLANNING_TOL`。`max(0, eta.X) + model.MaxVio <= PLANNING_TOL` 时直接返回原始运行状态；两项共用误差预算。正常求解得到正 eta 超限时，默认继续执行算法规定的支撑平面 LP 取割；显式 score_only=True 时缓存锥方向，推迟到候选获选后取割。超时、异常终止、证书精度不足或无法分离候选点均报错，不重建状态、不重求、不返回未确定结果继续运行。

| 新增 / 迁移量 | 固定代码映射 | 定义 |
|---|---|---|
| 只评分入口 | [SubProblem.solve.score_only](../model.py) | 默认 False；True 的不可行结果为已求得正 eta，cut=None、state=None，并含 cone_normals，不能将其直接当成一条割 |
| 锥的单位支撑方向 | [SubProblem.solve.cone_normals](../model.py)、[SubProblem.solve:cone_normals](../model.py) | 按 operation.cones 顺序的数值向量列表；锥尾为零时方向为零；不缓存求解器对象，不改变物理 state |
| 获选候选的取割入口 | [SubProblem.generate_cut](../model.py) | 固定原 x,p，以缓存 cone_normals 构造必要支撑 LP，返回原布局 cut；不重解 SOCP |
| SP 共同建模与 LP 求解 | [SubProblem._build](../model.py)、[SubProblem._cut](../model.py) | 将原 solve 内建模与取割步骤原样提取；eta 的定义位置显式迁到 [SubProblem._build.eta](../model.py)，变量名 violation、单位、约束与默认求解顺序不变 |
| 实际取割 LP 次数 | [SubProblem.cut_calls](../model.py) | 次；每次 LP optimize 前加一，与原 calls 的 SP 次数分开 |

SP 直接使用 `eta` 目标、`Aggregate=0`、`ScaleFlag=0`、`BarQCPConvTol=1e-9`。前两项避免预处理聚合及额外缩放，最后一项是内点法目标收敛设置。原始可行证书仍须通过 `eta + MaxVio <= PLANNING_TOL=1e-8`，割仍须通过 LP 和候选分离检查。数值设置固定用于所有点，不按失败结果重试或切换。

当前割仍由锥的有效支撑平面 LP 产生。令该 LP 在 eta=0 的线性行系数为 \((M_x,M_p,M_y)\)，右端为 \(b_{\rm row}\)。乘子方向采用 Gurobi 行约定：`>=` 行非负、`<=` 行非正、等式自由；固定 x/p 的等式乘子先置零。

| 数学量 | 代码定义 / 布局 | 说明 |
|---|---|---|
| \(\lambda\) | [SubProblem._separating_cut.dual](../model.py) | `(n_row,)`，对应 `rows`，不是旧锥标准式的行序 |
| \(M^T\lambda\) | [SubProblem._separating_cut.coefficients](../model.py) | 按求解器变量索引，不能按猜测切分 |
| \(h=M_y^T\lambda\) | [SubProblem._separating_cut.h](../model.py) | `(n_y,)`，通过 `state_ids` 提取全部状态列 |
| \(h_+,h_-\) | `maximum(h,0)`、`minimum(h,0)` | 内联派生式，没有公开字段 |
| \(\alpha,\beta,\delta\) | [SubProblem._separating_cut.cut](../model.py) | `cut[0]`、`cut[1:1+d]`、`cut[1+d:]` |

唯一割方向是

\[
\alpha+\beta^Tp+\delta^Tx\ge0,
\]
\[
\alpha=-\lambda^Tb_{\rm row}+h_+^Tu^{\rm glob}+h_-^Tl^{\rm glob},\quad
\beta=M_p^T\lambda,\quad \delta=M_x^T\lambda.
\]

`p` 仍为 kW，`beta` 必须作用于 kW 坐标。实现先作正比例归一化，再给截距加 `1e-10` 的保守补偿；不得把数值割系数直接解释成未经缩放的原始乘子。符号方向必须与上述行约定同时检查，不能搬用旧式 `-B.T @ dual`。`cut` 的布局、正负号、作用范围都属于契约。

[SubProblem.solve](../model.py) 返回键 [SubProblem.solve:feasible](../model.py)、[SubProblem.solve:state](../model.py)、[SubProblem.solve:cut](../model.py)。默认调用中 `feasible=True` 表示已获原始可行证书，否则必须返回有效割。新增 score_only=True 的确定正 eta 结果按上表延迟取割，不代表求解失败或状态未知；旧版未确定返回仍为异常，见第 23 节。

## 7. 连续域、几何坐标和覆盖证书

| 数学量 | 固定映射 | 单位 / 形状 |
|---|---|---|
| \(b^{\rm box}\) | [RegionState.bounds](../region.py)，由总负荷上界构造或由算例给定有效盒界 | kW，`(d,)`，公共正数坐标尺度 |
| \(b^{\rm axis}(\mathcal B)\) | [RegionState.axis_bounds](../region.py)、[RemainingRegionModel.__init__.axis_bounds](../model.py)、[build_continuous_region:axis_bounds](../continuous.py) | kW，`(d,)`；方向 MP2 全局上界，限制本预算的搜索域 |
| \(w\)，方向目标 \(w^Tp\) | [MasterProblem.direction](../model.py)、[MasterProblem.__init__.direction](../model.py) | 无量纲，`(d,)`；默认全 1，轴向初始化取单位向量 |
| \(p_\Sigma^{\rm ub}\) | [RegionState.total_bound](../region.py) | kW；由输入总量上界和 MP2 全局上界收紧 |
| \(\tau\)、\(s_\tau=1-\tau\) | [RegionState.tau](../region.py)、[build_continuous_region.tau](../continuous.py) | 无量纲；linear 为 0 |
| \(I_x,O_x\) | [RegionState.add_scheme:inner](../region.py)、[RegionState.add_scheme:outer](../region.py)，存在 `records[tuple(x)]` | 归一化顶点 `(n_vertices,d)`，分别为认证内域、候选外域 |
| \((F_x,g_x)\) | [RegionState.inner_equations](../region.py) / [halfspaces](../region.py) 返回 `[F_x, g_x]` | `(n_faces,d+1)`；`F_x @ xi + g_x <= 0` |
| 裁剪常数、系数 | [clip_polytope.constant](../region.py)、[clip_polytope.coefficient](../region.py) | `constant + coefficient @ xi >= 0`，与半空间内侧约定相反 |
| \(\gamma\)，未覆盖量 | [RemainingRegionModel.__init__.delta](../model.py)，求解器名称 `uncovered_distance` | 内部乘 [RemainingRegionModel.distance_scale](../model.py)，返回前还原 |
| \(\overline\gamma\)，全局覆盖上界 | [build_continuous_region:coverage_bound](../continuous.py) | 无量纲；不是负荷上界 |

`bounds` 是向量；查询结果 `bound` 是目标的标量全局界；`lb/ub` 是运行变量盒。三者不得互相重命名或复用含义。

坐标跨模块规则：

| 接口 / 数据 | 输入或存储坐标 |
|---|---|
| `SubProblem.solve(..., power)`、`MasterProblem(..., power=...)` | kW |
| `RegionState.add_point(..., point)`、`covering_schemes`、`witness_support` | \(\xi\)，无量纲 |
| `build_continuous_region` 中的 `point` | \(\xi\)，直接调用 SP 前乘 `bounds` |
| `build_continuous_region` 中的 `witness['p']` | kW；见证转成 `point` 时除以 `bounds` 并乘 `1-tau` |
| `RegionState.finish` 导出的 `inner/outer[*]['vertices']` | kW，已乘 `bounds` |
| `plot.sample_region(..., points, bounds)` | 输入 `points` 为 kW，内部归一化 |
| `halfspaces`、`contains`、`clip_polytope` | 不自动换算；输入必须处于同一坐标系，主流程使用 \(\xi\) |

固定方案的割转换为归一化裁剪：\(\alpha+\delta^Tx+(\beta\odot b^{\rm box})^T\xi\ge0\)。跨方案只取并集，不能混合顶点取凸包。

展示层 [region_geometry:x](../plot.py) 原样导出型号向量的独立列表，顺序与 `network.type_keys`、联合割的 `delta` 一致；几何顶点仍为 kW。这是回放格式的可选新增字段，旧记录缺少 `x` 时只显示原生成方案的条件切面，不推测其他方案的截面。最终结果中未保存的单方案外域不能视为不可行。页面的“已证覆盖”仅用一个三维凸内域包含整个候选外域这一充分条件（归一化面余量容差 `1e-10`），不反馈求解器，不替代全局覆盖证书。

覆盖目标统一写作

\[
\max_{x,p}\min_{x'\in\mathcal X_{\rm known}}\max_f
\{F_{x',f}s_\tau(p\oslash b^{\rm box})+g_{x',f}\}.
\]

这里 \(f\) 是带下标的几何面索引，区别于第 4 节的虚拟流向量。只有剩余域的全局上界不超过 `GEOMETRY_TOL`，或明确不可行，才能报告覆盖完成。`delta` 是剩余模型内既有的局部标量名称；不得传播为联合割的 \(\delta\) 向量。

外扩使用几何中心 \(\xi_c\)、内切半径 \(r_{\rm in}\) 和 \(\varepsilon_{\rm geom}\)，写为
\(\xi_c+(1+\varepsilon_{\rm geom}/r_{\rm in})(I_x-\xi_c)\)，再除以 \(s_\tau\)。代码对应 `RegionState.finish` 内的 `center/radius/envelope`；不将电阻 `r` 或费用 `c` 重新解释为这些量。

## 8. 查询答案、AC 校核和最终结果

### 8.1 MP1 / MP2 的答案

[MasterProblem.__init__.power](../model.py)、[MasterProblem.__init__.min_total](../model.py)、[MasterProblem.solve.radial_gap_kw](../model.py) 均用 kW；[MasterProblem.__init__.budget](../model.py) 用投资单位。没有固定负荷和最低总量时是 MP2（最大化 `direction @ p`，默认总负荷）；有任一条件时是 MP1（最低投资）。割作为 [MasterProblem.__init__.cuts](../model.py) 传入；[MasterProblem.solve.incumbent](../model.py) 与 [MasterProblem.solve.start](../model.py) 延续原有初始解语义，前者优先。方向初始化不传 `incumbent`，避免其成本上界限制另一方向的合法网架。

主流程直接创建模型并调用 `solve`，不保留 `planning_query`。MP 的运行证书和目标间隙在 `MasterProblem.solve` 中判断；超时、无候选、异常终止或精度超限直接报错，不再调用 SP 补救或改写答案。明确不可行仍返回 `None`，主动设置的目标停止允许返回已获证候选。`cuts_only=True` 仍不包含运行证书。

| 字段 | 定义位置 | 固定语义 |
|---|---|---|
| `x`、`p`、`state` | [MasterProblem.solve:x](../model.py)、[MasterProblem.solve:p](../model.py)、[MasterProblem.solve:state](../model.py) | 型号向量、kW 负荷、完整运行向量；`cuts_only` 的 `state=None`；未找到候选改为异常 |
| `objective` | [MasterProblem.solve:objective](../model.py) | MP2 为 `direction @ p`，kW；默认全 1 时仍为总负荷；MP1 为投资 |
| `bound` | [MasterProblem.solve:bound](../model.py) | MP2 为全局上界，MP1 为全局下界；单位随目标变化 |
| `feasible`、`status` | [MasterProblem.solve:feasible](../model.py)、[MasterProblem.solve:status](../model.py) | `optimal/feasible`；质量不合格改为异常，旧 `unknown` 记录保留；不等同于独立 AC 认证 |

`None` 答案表示已证不可行。剩余模型的 [RemainingRegionModel.solve:bound](../model.py) 是 \(\overline\gamma\)，[RemainingRegionModel.solve:complete](../model.py) 是覆盖结论，不能按 MP1 / MP2 的目标值解释。

### 8.2 独立 AC 运行量

[ACPowerFlow](../vertify.py) 接收选定的 `OperatingTree`；其属性虽然叫 `network`，索引是树局部索引。`state(power, ell)` 返回 `(P,Q,v,u)`：均为 `(batch,n_tree)`，功率朝根向外；`v/u` 分别为受端 / 送端电压平方。内部 [ACPowerFlow._state.p](../vertify.py)、[ACPowerFlow._state.q](../vertify.py) 表示标幺节点负荷 \(d^P,d^Q\)，不作为 kW 接口传播。

`ell` 是电流平方，AC 残差为 \(P^2+Q^2-u\ell\)。[ACPowerFlow.classify](../vertify.py) 输出 `1/-1`（可行 / 已证不可行），未收敛直接报错，不转用全局模型；[ac_planning_query](../vertify.py) 返回经过 AC 校核的方案答案，并删除 SOCP 的 `'state'`，避免充作 AC 证书。`global_status` 仅供显式交叉核验，未获证同样报错。

### 8.3 连续域结果和网格

| 字段 / 量 | 代码定义 | 单位 / 布局 |
|---|---|---|
| `max_total`、`max_total_bound` | [build_continuous_region:max_total](../continuous.py)、[build_continuous_region:max_total_bound](../continuous.py) | 最大总负荷的可行值、全局上界，kW；未获值可为 `None` |
| `max_point`、`max_choice`、`max_cost`（2.1 退役） | 原 `build_continuous_region` 字段；迁移至 [RunMonitor.seed](../monitor.py)、[RunMonitor._geometry](../monitor.py)，见第 33 节 | 初始化 kW 点、具名方案与费用不再重复写入最终结果 |
| `inner`、`outer`、`vertices` | [RegionState.finish:inner](../region.py)、[RegionState.finish:outer](../region.py)、[RegionState.finish:vertices](../region.py) | 最终顶点为 kW；`inner` 各项另含 `choice/cost`；`outer` 为全局外包络 |
| `status`、`counts`、`timing` | [build_continuous_region:status](../continuous.py)、[build_continuous_region:counts](../continuous.py)、[build_continuous_region:timing](../continuous.py) | 成功返回 `certified`，未完成改为异常；旧 `unknown/time_limit` 记录保留；`sp/cuts` 次数；`total_seconds` 秒 |
| 三态网格 | [BenchmarkResult.states](../plot.py) | `(len(METHODS), n_budgets, N_g, N_g, N_g)`，`-1/0/1` |
| 方法轴顺序 | [METHODS](../plot.py) | `('socp','hybrid','ac','linear')`；不要猜下标 |
| 配置及元数据 | [BenchmarkResult.metadata](../plot.py) | `schema_version`、`network_fingerprint`、`load_nodes`、`bounds`、`budgets` 等 |
| 每轴网格数 \(N_g\)、网格宽度 | [DIVISIONS](../main.py)、[BenchmarkResult.spacing](../plot.py) | 网格中心为 `(index + 0.5) * bounds / divisions`，kW |
| 展示比较标签 | [comparison_labels](../plot.py) | `0` 域外、`1` 多余、`2` 遗漏、`3` 重合、`4` 未确定；不是三态 `states` |

`state` 单数专指物理运行向量；`states` 复数为分类网格；主流程的 `region` 是 `RegionState` 几何对象。`choice` 保持外部具名方案，不能变成压缩整数编码。JSON 中非有限值转为 `null`；预算中的 `null` 读回 `inf`，其余上下界须保留未知语义，不能普遍将 `null` 变为零。

设报告的认证内域并集为 \(\mathcal R\)，AC 参考域为 \(\mathcal A\)，体积测度为 \(\mu\)。固定指标：

| 数学量 | 输出字段 | 定义及单位 |
|---|---|---|
| \(V_{\rm inner},V_{\rm outer}\) | [region_metrics:inner_volume](../plot.py)、[region_metrics:outer_volume](../plot.py) | 各自并集体积，kW³；重叠部分只计一次 |
| \(g_V\) | [region_metrics:volume_gap](../plot.py) | `(outer-inner)/outer`，比例；外域体积为零时实现返回 0 |
| FR | [disagreement:fr_percent](../plot.py) | \(100\mu(\mathcal R\setminus\mathcal A)/\mu(\mathcal R)\)，百分数 |
| MR | [disagreement:mr_percent](../plot.py) | \(100\mu(\mathcal A\setminus\mathcal R)/\mu(\mathcal A)\)，百分数 |
| 区域误差 | [disagreement:region_error_percent](../plot.py) | \(100\mu(\mathcal R\triangle\mathcal A)/\mu(\mathcal R\cup\mathcal A)\)，百分数 |

网格估计不等于连续几何证明；历史 AC 未确定点由 `disagreement_interval` 保留区间，新运行的 AC 未收敛则报错。FR/MR 的空分母返回 `None`，空并集的区域误差返回 0。不要把 `volume_gap` 比例值直接按 `*_percent` 显示。

## 9. 容差和运行参数

| 数学量 / 用途 | 唯一配置名 | 当前值 / 单位 |
|---|---|---|
| \(\varepsilon_{\rm plan}\)，求解结果接受容差 | [PLANNING_TOL](../model.py) | `1e-8`；Gurobi 所建模型的未缩放违反量，SP 另计 eta；尺度迁移见第 12 节 |
| \(\varepsilon_{\rm geom}\)，归一化覆盖 | [GEOMETRY_TOL](../region.py) | `1e-8`，无量纲；与物理容差独立 |
| AC 运行限值容差 | [AC_TOL](../vertify.py) | `1e-9`，标幺尺度 |
| AC 电流等式迭代容差 | [FIXED_POINT_TOL](../vertify.py) | `1e-12`，\(P^2+Q^2-u\ell\) 残差尺度 |
| AC 显式全局参照模型认证容差 | [GLOBAL_AC_TOL](../vertify.py) | `1e-7`，用于残差与运行违反量 |
| 径向收缩 \(\tau\) | [REGION_TAU](../main.py) | `0.005`，`0 <= tau < 1` |
| MP2 目标间隙容差 | [MasterProblem.solve.radial_gap_kw](../model.py) | `1e-3` kW；保留既有参数名，仅用于目标间隙判定，不再修改负荷，不是 `tau` |
| 单次 MP / SP / 剩余域 / AC 时限 | [MP_TIME_LIMIT](../model.py)、[SP_TIME_LIMIT](../model.py)、[RESIDUAL_TIME_LIMIT](../model.py)、[AC_TIME_LIMIT](../vertify.py) | 秒；SP 按方法取值 |
| 整体时限、线程 | [CASE_TIME_LIMIT](../main.py)、[SOLVER_THREADS](../main.py)、[DEFAULT_SOLVER_THREADS](../model.py) | 秒、正整数；混合阶段共享整体时限 |
| 全局搜索间隔、AC 迭代上限 | [REFINEMENT_CHECKS](../continuous.py)、[AC_ITERATIONS](../vertify.py) | 正整数；前者为自上次全局搜索以来的 SP 次数，活动见证处理完后生效 |
| 算例 / 升级数 / 预算 / 剩余域模式 | 原 `NETWORK/UPGRADE_COUNT/BUDGETS` 入口于 2.1 退役；当前 [main](../main.py)、[BUDGET](../main.py)、[RESIDUAL_MODE](../continuous.py)，见第 33 节 | 默认二维 FourBus、单预算；其他算例仍可单独使用其类 |

数值取值以对应代码常量为准；调整时同步本表及结果元数据说明。不得将不同语义的容差合成一个通用 `tol` 后改变认证口径。

## 附录：既有局部记号与历史文档（原第 10 节）

保留现有接口边界上的 `power ↔ p`、`state ↔ y`，不新增 `p_load/demand/load_power` 等同义接口。局部临时别名可以服务实现，但不能变成第二套公开数学名称。

| 既有局部记号 | 适用范围 | 阅读 / 迁移规则 |
|---|---|---|
| `p,q = tree.loads(power)` | 状态恢复、AC 内核 | 是 \(d^P,d^Q\) 标幺值；`power` 才是 kW 坐标 |
| `delta` | `RemainingRegionModel.__init__` | 覆盖目标 \(\gamma\) 的缩放标量；不是联合割的系数向量 |
| `current` | AC 后备解、SP 电流精修局部 | 仍是 \(\ell\)，不是电流幅值；新接口使用 `ell` |
| `e/net/m/c` 等对象别名 | 单函数内部 | 按对象类型解释；不改变公式中走廊、电阻、矩阵、费用的含义 |
| 历史推导 \(B,d,q,f\) | `socp_model.md` 固定方案推导 | 分别对应预算 \(\mathcal B\)、标幺有功 \(d^P\)、标幺无功 \(d^Q\)、功率因数；不是当前矩阵 B、维数 d 或虚拟流 f |
| 历史推导 \(A,H,J\) | `socp_model.md` 的消元式 | 固定树专用矩阵，维数与当前候选图标准式不同；不可直接赋给 `equations.A` 或 `network.corridor_types` |
| 历史 \(\beta_0+\beta^Td\ge0\) | 固定方案割 | 需换算 kW，且仅对该方案有效；不能当成当前全局 `[alpha,beta,delta]` |

历史报告、源码快照和源数据不做机械重命名。新主线文档统一使用本文符号；引用历史推导时明确它的假设、树索引和单位。

## 11. 后续修改必须遵守的规则

1. **先查表。** 修改数学代码前找到对应条目；新增数学量先登记“符号、中文定义、固定代码名、单位、形状、索引、产生和消费位置”。不要先造一套变量名再补文档。
2. **实现变化保留名称。** 不随求解器 API、矩阵存储方式、函数拆分或个人偏好改名，不把同一个量重新包装成另一套公开字段。
3. **语义变化显式迁移。** 单位、正负号、平方关系、索引、切片、字典键、矩阵朝向变化都算契约变化，即使名字没变。记录旧→新映射、理由和调用方；涉及保存格式时更新 schema 或提供显式兼容读取。
4. **改名作为独立事项。** 用户明确要求改名，或存在必须纠正的语义错误时，单独说明范围；按项目 GitNexus 规则做 impact 与 rename，同步公式、全部调用点、测试、输出字段及本文。不要为了让检查通过直接删掉登记项。
5. **代码就近解释。** 公式型函数的说明写清输入和输出单位、索引域与对应公式。求解器工作数组的简称只留在内部；从外部库取出的量马上转换到约定表示。
6. **交付时验证。** 先运行下述无求解器的名称检查；算法变更再运行涉及的物理、割有效性、几何或输出回归。提交前按 AGENTS.md 做图变更分析。

```text
python -m unittest tests.test_notation -v
```

本文中的代码链接同时作为可检查登记项：`Class.field` 检查类字段 / 属性，`Class.method.parameter` 或 `function.parameter` 检查该作用域的定义，`Class.method:key` 检查该函数显式构造或更新的字典键。检查只读取 Python 语法树，不导入模型、不使用 Gurobi 许可证；也会随常规测试发现运行。

**检查边界：** 它能发现已登记名称在指定作用域消失，以及文档指向错误文件；不能证明公式、单位、数组轴顺序或每条返回路径正确，也不能识别任意新变量是否是同义别名。这些继续由本规范和相关数值测试约束。

## 12. 1.1 迁移记录：直接返回求解器结果

- `PlanningEquations.margin(x, power, state)` 从生产接口迁出，测试与审核改用 `tests.planning_checks.margin(equations, x, power, state)`；不再让生产建模和手工检查共同决定每次求解结果。
- `PlanningEquations.restore` 删除。`state` 的含义、单位、切片及返回键不变，数值直接来自 Gurobi；不再沿树重建、裁剪电流或清零无负荷支路。
- `planning_query.radial_gap_kw` 保留名称和 kW 单位，但取消负荷回退用途，仅保留原有 MP2 目标间隙判定。`p/state/objective` 对应同一次求解结果。
- 数值质量改读 Gurobi `MaxVio`，覆盖模型适用的约束、变量界和整数违反量。二次约束采用求解器所建的二次多项式尺度，不再以手工 SOC 范数残差替代；两者不能视为数值恒等。旧范数口径保留在测试中，用于边界和退化锥回归，`PLANNING_TOL` 数值不放宽。
- SP 的 eta 定义、目标缩放、对偶割方向和 `feasible/state/cut` 输出键均不变；删除状态重建及固定 eta=0 后再次最小化损耗的补救分支。独立 AC 校验不受本迁移影响。
- SP 启用 `NumericFocus=2`。回归发现默认设置在近零 eta 的边界点可能返回略超 `1e-8` 的原始残差；提高数值稳健性后直接解通过原有测试检查，无需重建、重复求解或放宽阈值。

## 13. 1.2 迁移记录：主流程直接调用模型

- 删除生产接口 `planning_query`：其 `power/budget/min_total/fixed_plan/cuts/threads` 由当时的 `PlanningModel` 建模接收，`incumbent/start/radial_gap_kw` 由 `PlanningModel.solve` 接收；参数名、单位、初始解优先级和目标间隙容差不变。`deadline` 由构域入口管理，传给各模型的是该次剩余 `time_limit`。
- 删除 `solve_region` 和 `ContinuousRegion` 调度层。`build_continuous_region` 直接调用 MP、SP、剩余域模型；`region.py` 只保留几何状态、选点、并集覆盖及裁剪。`RegionTimeout` 移到 `main.py`，不增加新的求解包装函数。
- `ContinuousRegion.tau` 迁到构域入口的 `tau` 与 `RegionState.tau`；原 `ContinuousRegion.finish` 的所有结果字段迁到 `build_continuous_region`，名称、单位和存储格式不变。`point/pending` 保持归一化坐标，`seed/witness['p']` 保持 kW。
- MP 可行点立即登记，防止 MP1 超时丢失 MP2 证书；按完整 `x` 合并初始方案，保留同方案的不同可行点。剩余域见证仍先补同方案支撑点，避免物理种子造成自覆盖。每批上限仍为 `REFINEMENT_CHECKS`，全局完成仍只由覆盖证书决定。
- hybrid 仍共享整体时限、只继承有效割，SOCP 阶段独立认证内域。冻结的历史源码、已有结果和 SP 的 eta / 对偶割代数保持原口径。


### 10.7 完整路线勘察与二维信息上下域

正式模块 `survey.py` 固定建设预算，并将勘察费用单列。负荷坐标就是真实节点 1、2 的 p₁、p₂，单位 kW。既有线路仅保留原型号、无免费新增方案，故初始获准域与既有调度域相等；这不是一般升级模型的恒等式。

| 符号 | 实验记录 / 映射 | 单位及边界 |
|---|---|---|
| 完整路线 α、θ_α | `routes`、`road_usable_probability`、`survey_observation` | 一次勘察的原子对象。当前五节点的 A–F 各是一条完整候选道路；历史串联路线例见冻结契约 |
| K⁺、K⁻ | `known_usable`、`known_unusable` | 已证可用及不可用的路线集合；未知真实状态只由观测模拟器读取，不进入评分 |
| D⁻_t、D⁺_t | `confirmed`、`optimistic` | 相同预算下分别只准 K⁺、或允许所有未证不可用路线；不同于数值几何的内外近似 |
| S₀、S⁻_t、S⁺_t | `baseline_area`、`confirmed_area`、`optimistic_area` | 连续二维面积，kW²；用户所称 Soptimal 对应会随负观测缩小的 S⁺_t，非已知真实最优域 |
| ΔS⁺_α、ΔS⁻_α | `gain_if_usable`、`removal_if_unusable` | 确认可用时 D⁻ 的新增面积；确认不可用时 D⁺ 的删除面积，均为非负面积 |
| π_α ΔS⁺_α、(1−π_α)ΔS⁻_α | `expected_expansion`、`expected_exclusion` | 预期能力扩展与预期排除，kW²；后者不是新增供电收益 |
| V_α、V_α/c_α^survey | `marginal_information_value`、`information_efficiency` | V 为上述两项之和，即期望信息上下域间隙减少；kW²、kW²/元；不同于仅计扩展的既有 Δ^info/RIE |
| G_t、τ_S | `information_gap`、`stopping_area_tolerance` | G=area(D⁺\D⁻)，kW²；停止须同时检查边际值及整个剩余间隙上界，不能假定贪心全局最优或边际值递减 |
| I_x、O_x | `inner`、`outer` | 通过对偶主线获得的数值内外域；固定方案内求凸包；方案间只作几何并集，不跨方案取凸包 |
| 数值面积区间 | `_lower`、`_upper` | 由固定方案并集的内外界传播；信息不确定性与数值近似误差分别记录 |

本实验采用有限候选全集、准确二元观测、独立合成先验和单一候选型号。原子路线合并不消除不同完整路线之间所有可能互补性。所有面积属于现有 SOCP 模型；有限 AC 核查不能宣称整个多边形已被 AC 认证。


### 10.9 五节点概念算例与逐轮评分历史

`Network/concept5.py` 与 `survey.py` 为独立设计的说明性算例，不标称 IEEE 或实测系统。沿用 10.7 的符号、字段、SOCP 模型及集合语义；建设预算和勘察费仍分开。节点 0 为电源，1、2 为二维可变负荷，3、4 为必接入中间节点。4 条既有线路加 A–F 六条完整候选路线；不是全网总共六条线路。所有阻抗均为正，电压与送端有功容量约束保留。

| 数学量 / 记录 | 实现映射 | 单位与边界 |
|---|---|---|
| 线路参数 | `capacity_kw`、`r_ohm`、`reactance_ohm` | kW、ohm；分别转换为既有 `capacity`、`r`、`reactance` 标幺参数，不改变核心模型含义 |
| 初始占比 | `baseline_share` | 初始确认面积 / 初始乐观面积，无量纲；用于检查升级幅度是否适合概念图，不是算法性能指标 |
| 非凸见证 | `nonconvexity_witness` 中 `p_a`、`p_b`、`p_mid`、`plan_a`、`plan_b` | 点为长度 2 的 kW 向量；中点为两点平均，两端各有固定建设方案；不跨方案取凸包 |
| 中点排除距离 | `midpoint_distance_lower` | 中点到所有预算内已知可用方案数值外域并集的距离，kW；正值给出非凸性的保守几何证据，另用完整整数模型检验不可行 |
| 非凸缺口面积 | `convex_hull_excess_area` | 已确认域凸包相对已确认域的额外面积，kW²；只是描述量，不能替代中点不可行证据 |
| 逐轮评分 | `all_candidates` 的 `step` | `step=t` 为完成 t 次勘察后、决策第 t+1 次之前；保留所有当时未知路线的既有评分字段。已勘察路线在后续轮次缺失，不能画成价值降至 0 |

合成参数为展示 3 次可用、2 次不可用、1 条不勘察而设计，不能据此宣称普遍效果或贪心最优。F 的建设费在预算内，其停止勘察由残余区域与边际值共同判定；不预先赋予 F 可用/不可用真值。停止容差和数值域误差单独登记。



## 17. 1.3 恢复记录：对偶主线、二维几何和五节点候选案例

- 恢复来源为 `d25c070`；恢复前 v2 契约和源码备份已移入项目外清理归档；1.3 登记完整保留在 [历史契约](notation-history.md)。原第 10 节局部记号表改为附录，保留道路研究第 10.1–10.9 节编号。
- `PlanningSP`、`eta`、`cut=[alpha,*beta,*delta]`、`RemainingRegionModel`、`RegionState` 重新成为当时的活动接口。v2 的 `mode/epsilon_kw/branch_cuts/Cell` 等接口退出主线，历史结果不转换，不把旧 `tau` 和 kW 距离容差混用。
- 几何输入 `points/vertices` 为 `(N,d)`，`bounds` 为 `(d,)`，`d=2` 或 `3`。所有割仍为 `(1+d+t,)`，其中负荷系数始终乘 kW；`RegionState` 内部仍使用 `xi=p/bounds`。二维 `polytope_volume` 是面积，三维是体积，乘 `prod(bounds)` 恢复 kW^d。
- 原三维 AC 网格和实时查看器保留；二维五节点由专用勘察比较入口导出结果与图，不把二维数组伪装成三维网格。

| 数学量 / 记录 | 实现映射 | 单位、形状及边界 |
|---|---|---|
| 五节点候选案例 | [Concept5](../Network/concept5.py) | 完全沿用第 10.9 节的 4 条既有线与 A–F 六条候选路线，预算 4，p=(p1,p2)；不增加型号或改变参数 |
| 道路允许掩码 b_e^road | [Network.road_allowed](../Network/__init__.py)、[Concept5.__init__.allowed_roads](../Network/concept5.py) | 布尔 `(m,)`；默认全部允许；MP 加 z_e≤b_e^road，SP 不加道路掩码，因此物理联合割可跨信息状态复用 |
| 初始共享联合割 | [build_continuous_region.cuts](../continuous.py) | 与既有 `cuts` 相同的数组列表；传入后复制，不修改调用方；仅同一物理模型及有效盒界可复用 |
| 需求域缓存 | [DomainCache](../survey.py) | 按请求的道路允许集合求解，不预枚举建设方案；仅本次运行内缓存，记录各次冷启动及共享割耗时 |
| 比较面积误差 | `symmetric_difference_area`、`area_error_upper` | kW²；分别为新旧数值内域对称差、利用双方数值内外界得到的上界 |
| 面积停止精度 | `domain_area_tolerance` | kW²；数值面积区间须达标才能比较信息价值；不覆盖勘察停止容差 5 kW² |
| 参考计算 | `reference_directory`、`reference_seconds` | 历史比较输出字段，已退役；原值随历史结果保留，正式运行不依赖参考结果文件 |

`Network.fingerprint` 包含新增掩码，防止不同道路状态误用缓存。无模型参数变化时联合割可复用；预算/掩码仅进入 MP。勘察真值仍只由观测模拟器读取。数值 `outer` 与信息 `optimistic` 分别保留，不混同。

对偶域的 `plan_ids=null` 表示未枚举方案，不能凭两个 null 推断域相同；已有枚举实验仍保留原编号表及其等域快捷判定。新域插入缓存时，内域与已知子集内域取并，再与已知超集内域取交；外域与已知超集外域取交，再与已知子集外域取并。此前缓存不变，这些保守操作维持道路集合对应域的单调包含关系；处理后的数值面积区间仍须小于 0.5 kW²。


## 18. 1.4 清理迁移：正式勘察模块与历史实验退役

本次只整理模块和退役实验入口，不改变潮流模型、对偶割、道路评分公式或停止容差。原 1.3 契约的全部登记项原样保存在 [历史契约](notation-history.md)，不通过删除登记项掩盖改名。第 14–16 节和原第 10.1–10.6、10.8 节的实验入口及配套数据已归档，所述独立实验接口不再属于活动主线。

| 迁移对象 | 活动定义 / 处理 | 保持或变化 |
|---|---|---|
| 按需对偶构域 | [DomainCache](../survey.py)、[polygon_union](../survey.py) | 从 concept5_dual 提取，名称、字段、预算及割共享规则不变 |
| 面积增量与评分 | [increase_bounds](../survey.py)、[candidate_values](../survey.py) | 从旧勘察实验提取，形参、概率权重和全部结果键保持不变 |
| 信息状态与观测 | [information_state](../survey.py)、[simulate_surveys](../survey.py) | 确认/乐观域、两种观测分支及停止证书不变 |
| 非凸见证 | [nonconvexity_witness](../survey.py) | 形参和 p_a/p_b/p_mid/plan_a/plan_b 保留；solved 改由本次已发现方案组装，方案编号对应本次 schemes.json，不能索引旧枚举表 |
| 正式五节点入口 | [run_survey](../survey.py) | 新增正式编排入口；替代已退役的 run_comparison 实验入口，不读取旧结果、测试目录或参考面积 |
| 计算参数 | [BUDGET](../survey.py)、[STOPPING_AREA_TOLERANCE](../survey.py)、[DOMAIN_AREA_TOLERANCE](../survey.py)、[SURVEY_REALIZATION](../survey.py) | 预算 4、勘察停止 5 kW²、数值面积区间 0.5 kW²；示例观测仍为 A/C/E 可用、B/D 不可用、F 未赋真值 |
| 结果格式 | [run_survey:query_count](../survey.py)、[run_survey:joint_cut_count](../survey.py) | 既有数值含义不变；protocol.schema=survey-v1。reference_directory/reference_seconds 等旧比较字段不进入新结果；比较由测试夹具独立执行 |
| 输出文件 | results.json、domains.json、cuts.json、schemes.json、audit.json、protocol.json、trace.csv、all_candidates.csv | 数值结果及状态字段不改名；独立 AC 点核查和数值区间写入 audit.json。历史 comparison.json 不改写 |
| 图件入口 | survey_plot.py | 移动原绘图函数，不再导入实验目录或用户机器上的技能脚本；绘图只读取已计算结果 |
| 旧基准报告入口 | `plot.save_benchmark_report` | 仅被已退役的两个 benchmark 脚本调用，随实验入口归档；已生成的历史报告保留 |
| 核心测试依赖 | [joint_benders](../tests/planning_checks.py)、[affordable_designs](../tests/planning_checks.py) | 从临时 benchmark 脚本提取，原名称、形参和逻辑保留；只用于测试，正式主线不导入 |

所有历史实验的符号、单位与输出说明仍可在冻结契约追溯。清理范围、归档定位和保留结构见 [清理记录](cleanup.md)。

## 19. 1.5 整理迁移：计算、绘图与单次保存

本次不改变潮流方程、信息价值公式、数值精度、勘察停止证书或历史结果。第 18 节记录的旧输出布局保留用于追溯，新运行采用以下布局：

| 对象 | 旧布局 → 当前布局 | 含义与调用方 |
|---|---|---|
| 勘察结果 | 多个 JSON、逐查询文件、CSV 和报告 → 单个 `results.json`，`protocol.schema=survey-v2` | 既有 `states/trace/all_candidates`、面积、评分、非凸见证和统计字段不改名；参数、核查和方案分别合并至 `protocol/audit/schemes` |
| 道路域缓存 | `DomainCache(output, ...)` → `DomainCache(...)` | 去掉输出目录参数；道路集合、原始构域结果和共享割只在内存中保存，不再维护 `queries`、逐模式记录或单独的域/割文件 |
| 非凸见证方案 | `schemes.json` → `results.json` 的 `schemes` | `plan_a/plan_b` 仍对应本次方案编号，负荷点仍为 kW |
| 三维回放 | `events.jsonl`、`replay.json` 与内嵌相同内容的 HTML → 单个 `live_view.html` | 原始数值结果仍为 `result.npz`；完整回放只嵌入 HTML 一次，`RunMonitor.load_recording` 从该 HTML 恢复 |
| 绘图 | 全局输出目录、导入时设置样式 → 显式输出目录、局部样式 | `concept_figure/history_figure` 返回图对象，`render_survey` 统一保存 PDF/SVG/PNG；独立绘图只读取 `results.json` |

旧目录及其 JSON/CSV 不自动转换或删除。已有科学结果字段与数学登记项保留；此迁移仅减少重复持久化，移除逐查询日志及源文件指纹副本。结果中的网架指纹、模型参数、评分区间、认证状态和独立 AC 点核查仍保留。集合包含关系和面积间隙分解由回归测试检验，生产流程保留未认证域、无有效停止证书和未通过核查时的明确失败条件。

## 20. 1.6 整理迁移：勘察绘图统一入口

`survey_plot.py` 合并至 `plot.py` 并删除。勘察绘图只对外提供 [render_survey](../plot.py)，原 `concept_figure/history_figure` 的图形组装合并到该函数，几何绘制和曲线绘制仅作为函数内部的共用辅助步骤。原脚本命令改为 `python plot.py --results <结果目录>`。

`render_survey(result, output)` 的参数及 PDF/SVG/PNG 文件名保持不变；输入仍为 `survey-v2` 结果，`p_a/p_b/p_mid`、区域坐标、面积、评分、观测顺序和全部数值字段均不变。保留前两节的历史接口说明用于追溯，不保留旧模块转发或第二套绘图接口。

## 21. 1.7 类名迁移：物理规则、主问题与子问题

按职责将 `PlanningEquations` → `GridPhysics`、`PlanningModel` → `MasterProblem`、`PlanningSP` → `SubProblem`。`GridPhysics.add_operation` 仍向传入的 Gurobi 模型添加运行变量及电气约束；`MasterProblem` 仍建立完整的 MP1/MP2 规划问题（包括物理约束）；`SubProblem` 仍固定 `x,p` 检查运行可行性并生成联合割。只改类名和调用点，不改变 `equations` 等局部对象别名、数学量字段、形参、字典键、单位、状态布局或求解目标。历史源码和结果不改写；前述迁移记录保留其当时使用的类名。


## 22. 1.8 迁移：统一运行记录与单步执行

- `RunMonitor` 从 `plot.py` 迁到 `monitor.py`，类名和既有回放字段不改名；调用方直接导入新模块。
- `main.run` 统一管理二维勘察与三维构域的记录、输出目录和保存；`survey.run_survey` 移除 `output`，只返回既有 `survey-v2` 结果。`survey.py` 命令行委托 `main.run`，不再单独写结果。
- 每次运行的 `steps.jsonl` 逐事件保存版本 2 的增量帧，可在结束前读取；每 100 帧保存完整状态。`live_view.html` 是同一记录的离线导出，`results.json` / `result.npz` 保留既有科研结果格式。历史文件不转换、不删除。
- `checkpoint` 表示算法步骤边界；`step_by_step` 控制是否在边界等待。`waiting`、`paused_seconds` 与 `run_id` 是控制和记录字段，不是数学量。`clock` 返回扣除人工等待后的秒数，构域时限、阶段耗时及勘察耗时使用同一时钟；墙钟耗时另存 `wall_seconds`。
- 事件中的 `answer` 保存模型原返回字典，`point` 仍为展示用 kW，几何内部仍为归一化坐标。`point_reason` 只说明选点来源。`query` 和 `attempt` 区分道路集合构域及其 light / physical 尝试，割必须结合当次方案解释。
- 事件中的 `survey` 使用现有信息状态字段；`all_candidates` 为当前轮带 `step` 的评分，`route` / `survey_observation` 为实际决策和观测。信息上下域与数值内外域分别展示，不跨方案取凸包。
- 二维 `surface` 保留按边界排序的顶点，`faces=[]`；三维继续输出三角面。联合割仍按 `[alpha, *beta, *delta]` 切片，二维画条件直线、三维画条件平面。

## 23. 1.9 迁移：失败直接停止，移除自动补救

1. MP 不再接受异常终止、无候选或超容差结果，不再交给 SP 修复 `state/feasible/status`。已证不可行仍返回 `None`；求解器主动目标停止仍可返回合格候选，单位、字段及目标界含义不变。
2. SP 只保留“正常求解 → 可行证书或正 eta 的 LP 联合割”主线。超时、数值异常、无有效割改为 `RuntimeError`；原 `feasible=False, cut=None` 失败返回退出活动接口，历史记录不改写。
3. 剩余域无覆盖证书且无有效见证、割不产生进展及构域总时限耗尽均报错。`RegionTimeout` 向上传递，不再转成 `time_limit` 结果继续下一阶段；失败现场由 `RunMonitor` 写入步骤日志和回放，成功结果格式不变。
4. 道路域只执行一次 light 构域，不再失败后切换 physical；历史回放中的 `attempt` 保留，新记录恒为 1。保留基于已认证子集/超集的数学包含关系收紧，这不是异常补救。
5. 删除 `_convex_hull` 的异常后坐标变换重算，调用方直接使用 `ConvexHull`；支撑点线性方程奇异或找不到包含见证的单纯形时直接报错。AC 不收敛不再调用全局后备模型；`global_status` 仅作显式参照检查。空集、低维集合、可行/不可行证书、有效割和正常算法迭代均保留。

以上迁移改变失败处理，不放宽任何容差、不改负荷、不修改数学约束。测试侧的未确定点二次 SP 分离也已移除。保留此前各节的迁移历史和旧结果字段，旧回放仍可读取。


## 24. 2.0 迁移：方向初始化与单循环 SOCP 主线

1. `MasterProblem.direction` 新增无量纲 `(d,)` 目标向量。MP2 的 `objective/bound` 对应 `direction @ p`，默认全 1 保持原总负荷含义；成本最小化分支不变。
2. 新增 `axis_bounds`，由本预算方向 MP2 的全局上界产生，作为几何和剩余域约束；`bounds` 仍是公共正数归一化箱。`RegionState.tighten_bounds` 同步已有和后续网架；结果新增同名 kW 字段，旧结果缺少该字段仍按既存顶点显示。
3. 删除默认 MP1 初始化、`seeds/queue/pending/seed` 调度和 `RegionState.next_point`。所有选点与求解分支在 `main.build_continuous_region`，普通候选取真实总负荷最大，内部见证取缺少当前网架认证的最大负荷支撑点。没有有效支撑时直接失败，不退回另一种选点规则。
4. `REFINEMENT_CHECKS` 从 `region.py` 移至 `main.py`，从每个方案的批次上限改成全局搜索间隔；活动见证完成后检查。`evaluation_bounds` 退役：公共箱直接来自总负荷上界或算例，本预算紧外域由方向 MP2 产生；旧函数的调用方已迁移，历史记载保留。
5. 主入口默认仅计算 SOCP，AC 在所有预算构域完成后运行。未计算方法的网格维度仍保留为未知；报表只列实际计算的方法，不改变既有 `METHODS` 索引。
6. 删除主入口的通用输入检查和构域中重复的无进展诊断，不增加自动重试、重求或答案修复。保留可行证书、有效割、覆盖上界和求解状态的数学接受条件。
7. 步骤的 `answer` 不再重复存储运行向量 `state` 和联合割 `cut`；割仅在 `cut` 事件保存。结果去掉道路/型号参数副本、容量统计和容差字典，保留网架标识、预算、评价箱、精度、域、认证结论及耗时；已有字段的单位和历史文件不变。
8. SP 目标由 `1000*eta` 改为原始 `eta`，使用 `Aggregate=0`、`ScaleFlag=0`，不再覆写默认 `NumericFocus`；SP 的 `BarQCPConvTol` 从通用 `1e-10` 改为 `1e-9`。这改变求解器聚合、缩放与目标收敛设置，不改变物理接受阈值 `PLANNING_TOL=1e-8`、最优 eta 或有效割判据。联合割仍作同一正比例归一化；不添加失败后分支。

## 25. 独立 FourBus 外域实验：最大 SP 违反量

`experiments/fourbus_outer.py` 不接入 `main.py`。初始化、拓扑预算约束、物理约束及联合割沿用 `MasterProblem`、`GridPhysics`、`SubProblem`；不改变以上接口。实验外域是满足拓扑、预算、方向上界与所有联合割的 **(x,p) 联合集合**，负荷外域取其在 p 上的投影，不跨方案取凸包。

固定 x,p 的 SP 写为 `sp_y @ y + sp_eta * eta + sp_x @ x + sp_p @ p <= sp_rhs`，并保留 `cone_constant + cone_y @ y` 所属的各个 Lorentz 锥。状态变量的有限上下界也列入线性行；eta >= 0 单独处理。所有系数直接导出现有具名模型，不另写潮流方程。

| 符号 / 定义 | 固定代码映射 | 单位、形状、产生和消费位置 |
|---|---|---|
| SP 线性行系数、右端 | [export_sp](../experiments/fourbus_outer.py) 的 `sp_y/sp_eta/sp_x/sp_p/sp_rhs` | 行数为原线性行加 2*len(y)；列序分别为既有 y、eta、x、p；系数沿用原标幺残差与 kW 输入 |
| 锥仿射映射 | `cone_y/cone_constant/cone_slices` | 各锥依次堆叠头及尾；列按既有 y；由 `operation.cones` 导出，交给完整锥对偶 |
| μ：线性 <= 行的非负乘子 | [GlobalViolation.row_dual](../experiments/fourbus_outer.py) | `(n_rows,)`，归一化 `-sp_eta @ row_dual <= 1`；不是已有取割 LP 的带符号 `dual` |
| s：自对偶 Lorentz 锥乘子 | [GlobalViolation.cone_dual](../experiments/fourbus_outer.py) | 堆叠锥维数；头非负，尾自由；满足每个锥约束及 `sp_y.T @ μ = cone_y.T @ s` |
| 候选点的保守对偶下界 | [GlobalViolation.candidate_lower](../experiments/fourbus_outer.py) | eta 尺度；将 μ 取非负并归一化、将 s 的头提升至不小于尾的范数，再用 y 的有限盒补偿驻点残差，保证 `objective` 是下界而非未经检查的数值目标 |
| x_j (sp_x.T μ)_j | [GlobalViolation.choice_term](../experiments/fourbus_outer.py) | `(t,)`；用二元指示约束精确表示，不人为截断无界乘子 |
| R：全局最大最小违反量 | [GlobalViolation.violation](../experiments/fourbus_outer.py) | 标量，沿用 SP 原始 eta 尺度；目标为 max R，R <= μᵀ(sp_x x+sp_p p-sp_rhs)−sᵀcone_constant；不能解释为 kW 或几何距离 |
| G 候选下界与全局上界 | [GlobalViolation.solve:objective](../experiments/fourbus_outer.py)、[GlobalViolation.solve:bound](../experiments/fourbus_outer.py) | 同 R；objective 由 `candidate_lower` 计算，bound 取本轮求解器界、历史有效界和解析界的最小值，并约束后续 G。外域只收缩，故旧上界仍有效；候选下界不能用来判断停止 |
| 总负荷方向上界 | [run_experiment.total_bound](../experiments/fourbus_outer.py) | kW，来自全 1 方向 MP2 的 bound；`axis_bounds` 沿用第 24 节定义 |
| 外域停止阈值 ε | [run_experiment.epsilon](../experiments/fourbus_outer.py) | 原始 eta 尺度，不是 PLANNING_TOL、tau 或 kW；只认证残差精度，不认证面积误差或 AC 可行性 |
| 已有有效联合割 | [run_experiment.cuts](../experiments/fourbus_outer.py) | 沿用既有 `(1+d+t,)` 数组列表；仅复用同一 FourBus 物理模型的割，默认空；`initial_cut_count` 记录继承条数，不继承旧 G 的目标或上界 |

FourBus 中零潮流、v=1、压降松弛为零满足全部 SP 硬约束，故 eta* <= max(p)/base（该算例无固定负荷且 q_ratio < 1）。G 使用这条解析上界，不使用任意乘子大 M。固定方案去掉未选型号的零变量和冗余锥后，松弛等式问题有严格可行点；完整锥对偶用于表达最小 SP 值。实验须数值核对固定候选的原始与对偶值、割的全局有效性；求解器上界及结论受数值容差约束。

候选下界使用 h=sp_y.T μ−cone_y.T s，计算 μᵀ(sp_x x+sp_p p−sp_rhs)−sᵀcone_constant + h⁺ᵀy_lb_global + h⁻ᵀy_ub_global，再与 0 取最大。它是有界 y 域上的拉格朗日下界；驻点等式精确成立时补偿项为零。不能只凭平方锥约束的 `MaxVio` 把数值目标当成严格对偶值。早期试验的 `scaled_cone_dual` / `CONE_DUAL_SCALE=1000` 因放大全局求解的数值残差已撤回，实际 s 的名称与含义不变；登记保留在此作为实验迁移说明。

新增方案来自 G 对全部预算内二元组合的隐式搜索。某对 (x,p) 被切除不等于该 p 在其他方案中被排除。实验仅保存初始化、各轮候选/界/割、最终停止状态与耗时到独立结果目录；另导出 `outer.lp`，其中只有拓扑、预算、方向上界和联合割，可固定 p 后搜索是否存在 x，未展开方案或顶点。不更改主线输出格式。
## 26. FourBus 外域独立扫描与迭代绘图

`experiments/fourbus_outer_scan.py` 只读第 25 节实验结果，调用未加实验割的完整 SOCP 模型获得独立参考；`experiments/fourbus_outer_plot.py` 只读扫描输出。两者不进入主线。参考始终允许全部预算内建设变量 x 与运行变量 y 自由选择，不枚举或限制为已发现方案。此处“参考”指同一 SOCP 模型的有限采样，不是 AC 真值或连续域的完整证书。

| 符号 / 定义 | 固定代码映射 | 单位、形状及含义 |
|---|---|---|
| 三维扫描点及间距 | [scan_grid.power](../experiments/fourbus_outer_scan.py)、[scan_grid.grid_step](../experiments/fourbus_outer_scan.py) | kW，`(N,3)` 与标量；均匀体素中心，网格箱向上取整覆盖轴向上界 |
| 逐点 SOCP 参考标签 | [scan_grid.reference](../experiments/fourbus_outer_scan.py) | `(N,)`；1 表示完整模型可行，-1 表示全局判定不可行；异常或质量不合格直接停止，不计为不可行 |
| 首次被负荷投影外域排除的轮次 | [scan_grid.first_exclusion](../experiments/fourbus_outer_scan.py) | `(N,)`；0 表示初始化已排除，k 表示第 k 条割后首次排除，K+1 表示最终仍保留；固定 p、放开全部 x 查询联合割投影，利用嵌套性二分查找 |
| 负荷比例 w 与径向容量 ρ | [scan_rays.weights](../experiments/fourbus_outer_scan.py)、[scan_rays.radial_total](../experiments/fourbus_outer_scan.py) | w 无量纲、非负且和为 1；ρ 为 kW，p=ρw；每条射线分别最大化完整 SOCP 与联合割外域中的总负荷 |
| 方向扫描参考、外域容量 | [scan_rays.reference_total](../experiments/fourbus_outer_scan.py)、[scan_rays.outer_total](../experiments/fourbus_outer_scan.py) | kW，两个独立优化问题；同时保存各自全局上界。差值是相同负荷比例的总量差，不是欧氏距离 |
| 候选 SP 最优违反量 | [candidate_eta](../experiments/fourbus_outer_scan.py) | 固定历史 x,p 后重新求原始 SOCP 的最小 eta，沿用原始残差尺度；是事后复核，不能作为原运行时已计算的记录 |
| 有效历史 G 上界 | [run_scan.effective_bound](../experiments/fourbus_outer_scan.py) | eta 尺度；跨 61+20 轮取已有有效界的累积最小值。续算程序本身未继承此界，图注明确区别原记录和后处理 |
| 扫描体积估计 | `grid_step**3 * count` | kW³；体素中心分类的数值估计，不是内/外域连续体积证书；跨方案射线边界的连线/曲面也仅用于显示插值 |

完整迭代链须校验续算前缀割与原 61 条逐项一致。图中方案编号 S01… 按初始化及候选首次出现编号，仅表示实际访问的方案；G 候选的方案不等于已通过 SP 认证的方案。不把其他试运行拼接进这条迭代链。

边界射线及历史候选全局可行性复核统一使用 `NumericFocus=3`，历史候选复核另统一设置 `Aggregate=0`，以控制约束残差与边界点的聚合数值问题；不放宽 `PLANNING_TOL=1e-8`，不添加逐点失败后重试或改变参考物理模型。网格参考点沿用标准求解设置。

## 27. FourBus 负荷分块与角点上界实验

`experiments/fourbus_outer_partition.py` 复用第 25 节的完整锥对偶、拓扑预算及 SP 联合割，不修改主线或原双线性实验。固定 x 时 eta*(x,p) 对 p 凸。对负荷块 C=[lower,upper]，以 `X_C={x: 存在 p 属于 C 且 (x,p) 属于当前联合外域}` 筛选网架，计算 `U_C=max_{x in X_C, c in corners(C)} eta*(x,c)`，因此原联合外域在 C 内的最大违反量不超过 U_C。角点 c 可能不在该网架外域内，只能用于上界；实际送入 SP 的 p 必须位于该网架外域与 C 的交集。

| 数学量 / 定义 | 固定代码映射 | 单位、形状与使用位置 |
|---|---|---|
| 负荷块上下界 | [CornerViolation.solve_box.box_lower](../experiments/fourbus_outer_partition.py)、[CornerViolation.solve_box.box_upper](../experiments/fourbus_outer_partition.py) | kW，`(d,)`；约束 `problem.power`，保持 p 为外域交集见证的原语义 |
| theta=sp_p.T mu | [CornerViolation.theta](../experiments/fourbus_outer_partition.py) | `(d,)`，eta/kW；保留全部原始对偶约束和聚合等式；由 eta 归一化推出上下界，不任意截断乘子 |
| 角点选择位 b | [CornerViolation.corner_bits](../experiments/fourbus_outer_partition.py) | 二元 `(d,)`；c=box_lower+(box_upper-box_lower)*b；一次 MISOCP 隐式选择全部 2**d 个角点 |
| t_j=b_j theta_j | [CornerViolation.corner_term](../experiments/fourbus_outer_partition.py) | `(d,)`；有有效 theta 界的精确二元乘积线性化；不含连续双线性项 |
| 实际评估角点 c | [CornerViolation.solve_box:corner](../experiments/fourbus_outer_partition.py) | kW，`(d,)`；与返回键 p 区分，角点本身不要求满足当前割 |
| 角点下界 / 块上界 | [CornerViolation.solve_box:objective](../experiments/fourbus_outer_partition.py)、[CornerViolation.solve_box:bound](../experiments/fourbus_outer_partition.py) | eta 尺度；objective 只是在角点的对偶下界，不能当作外域反例；bound 取父块有效界、解析界及求解器全局界的最小值 |
| 块内真实候选 | [CornerViolation.witness](../experiments/fourbus_outer_partition.py) | 固定角点模型返回的 x 和乘子，在同一块与外域内最大化对偶负荷线性项；返回既有 x/p 字段，调用原 SP 检查 |
| 分块全局上界 | [run_trial:bound](../experiments/fourbus_outer_partition.py) | 所有待处理块及已认证块有效界的最大值；分裂子块继承父界，加割后旧界仍有效；不得跨不相交块传递较小的局部界 |
| 对照运行时间 | [run_trial:seconds](../experiments/fourbus_outer_partition.py) | 秒；包含本方法建模、求解、SP、切割与保存，不包含两方法共享且单列的 MP2 初始化；形参 seconds 是整体时限 |

原始双线性基线使用既有 GlobalViolation，不更改候选或停止策略。两方法共享初始化、预算、epsilon、线程、单次及整体时限；均只在全局上界不超过 epsilon 时报告认证。达到实验时限保留未认证状态，不作为失败补救。新增割来自同一个原 SP；区域仍为 SOCP 联合外域，不宣称 AC 或几何精度证书。

角点模型固定设置 `NonConvex=0`、`MIPGap=0`、`MIPGapAbs=1e-7`。初始原生 MISOCP 试验尝试过 `PreMIQCPForm=2`（锥分解预处理），首个块可解但后续块仍数值失败，因此不作为正式比较方法。正式循环不含失败后重试。

继续分块时，原生 MISOCP 即使采用上述设置仍出现数值终止。因此比较实验预先选用 [PolyhedralCornerViolation](../experiments/fourbus_outer_partition.py)：用 `s_head >= +/-s_tail[i]` 初始化锥外逼近，反复解 MILP 后补充 `s_head >= direction @ s_tail`，其中 direction 的范数不超过 1。全部平面均对完整对偶锥有效，MILP 的全局界仍是保守上界；候选下界继续用原完整锥投影与状态盒残差修正。这个循环是明确的锥外逼近算法，不是数值失败后切换求解器。`oa_iterations` 记录本次角点查询内部的 MILP 求解数，`oa_planes` 记录累计追加锥切面数。所有角点共享这些有效锥切面；块上界仅在同一块及其子块之间继承。

新增锥切面的整行固定乘 1000，以避免很接近的相邻切面在绝对线性约束容差下失去分辨率；不改变乘子变量或目标尺度。原始锥候选下界仍按未缩放的锥范数修正，`PLANNING_TOL` 不变。

角点模型的 `dual_max_violation` 只记录辅助对偶模型的最大数值残差，不作为原始物理可行性的判据。其候选乘子必须经 `candidate_lower` 的非负化、完整锥投影、eta 归一化及有限状态盒驻点误差补偿后才能用于反例判断；此下界即使原始对偶向量存在残差也有效。构造出的 x 在独立外域见证模型中固定，真实 x,p 的外域可行性和 SP 物理/割检查仍使用 `PLANNING_TOL=1e-8`。这明确改变辅助对偶候选的接受方式，不改变 SP 容差或模型。原生 G 基线保持原实现。

`experiments/fourbus_partition_report.py` 读取两方法的实际运行轨迹，并复用第 26 节独立 SOCP 射线参考。`reference_total/outer_total`、`max_radial_difference/mean_radial_difference/min_radial_difference` 保持原名称和 kW 含义；固定负荷比例下的容量差不等于全局 eta 上界。对全部新割在无预算限制的完整 SOCP 可行域上最小化割余量，逐候选复算原始 min eta，分别检查全局割有效性及候选下界/块上界关系。核查和绘图耗时不计入方法运行时间。

[match_boundary_accuracy](../experiments/fourbus_partition_report.py) 利用加入割后外域嵌套的性质，对两种方法的割前缀分别二分查询，在全部相同 861 条参考射线上寻找最大容量差首次不超过基线最终值的割数量（比较容差 `comparison_tolerance_kw=1e-5` kW）；`accuracy_match` 的 `cut_count/seconds/checked_prefixes` 对应分块法，`baseline_cut_count/baseline_seconds/baseline_checked_prefixes` 对应基线，`target_max_radial_difference` 是共同误差标准。时间取原运行中生成该条割的时刻，不把事后扫描耗时计入，也不把基线首次达标之后的时间算作其达标耗时。`speed_ratio` 是两者首次达标时间之比，仅对应这组有限射线，不是全局连续认证速度。

## 28. FourBus 单次割的外域体积停滞实验

仅独立实验 `run_trial` 增加可选 `volume_threshold`，默认 `None` 保持第 27 节行为。启用时，每次 SP 联合割后计算全部预算内网架外域在负荷空间的并集体积；不建立认证内域，不将各方案体积直接相加，不将并集替换成跨方案凸包。

| 数学量 / 定义 | 固定代码映射 | 单位、形状与使用位置 |
|---|---|---|
| 全部预算内合法建设向量 | [budget_schemes](../experiments/fourbus_outer_volume.py)、[ProjectedOuterVolume.schemes](../experiments/fourbus_outer_volume.py) | 二元 `(N,t)`；仅小规模 FourBus 的精确体积测量枚举 x，复用原拓扑预算模型逐方案排除；全局搜索 G 不使用这份列表 |
| 各方案当前外域顶点 | [ProjectedOuterVolume.polytopes](../experiments/fourbus_outer_volume.py) | 每项 `(n_vertices,3)`，坐标为 p/axis_bounds；逐条联合割裁剪，所有合法方案均保留 |
| V_k：第 k 条割后的投影外域并集体积 | [ProjectedOuterVolume.outer_volume](../experiments/fourbus_outer_volume.py)、[ProjectedOuterVolume.add_cut:outer_volume](../experiments/fourbus_outer_volume.py) | kW³；复用 `union_volume` 扣除重叠，按 `prod(axis_bounds)` 恢复单位；浮点连续多面体几何，非体素/随机估计 |
| 单条割减少的体积 | [ProjectedOuterVolume.add_cut:volume_reduction](../experiments/fourbus_outer_volume.py) | V_(k-1)-V_k，kW³ |
| 单条割的相对体积减少 | [ProjectedOuterVolume.add_cut:volume_reduction_ratio](../experiments/fourbus_outer_volume.py) | (V_(k-1)-V_k)/V_(k-1)，无量纲，以加割前的当前体积为分母 |
| 体积停滞阈值 | [run_trial.volume_threshold](../experiments/fourbus_outer_partition.py) | 本次实验为 0.001，即 0.1%；首条满足严格小于的割即停止，不额外增加连续次数或预热轮数 |
| 体积测量时间 | [ProjectedOuterVolume.add_cut:volume_seconds](../experiments/fourbus_outer_volume.py)、[run_trial:volume_seconds](../experiments/fourbus_outer_partition.py) | 秒；逐割量和累计量。累计含本方法的枚举/几何初始化，包含在 run_trial 原有 seconds 中 |

`status='volume_stagnation'` 是用户指定的启发式停止，`certified=False`；不能写成残差认证、真实域体积误差 <0.1% 或连续域无遗漏证书。一条割可能只收紧单个网架的截面，而被其他方案覆盖，导致投影体积没有变化。共享 MP2 初始化时间另列，算法总耗时为 `seconds+initial_seconds`；事后有效割审核与 861 条 SOCP 参考射线核查不计入算法耗时。此精确体积实现的方案枚举只用于小算例比较，不声称可直接扩展到大量建设变量。

## 29. FourBus：残差上界 0.1 停止与相对边界误差

本次只通过既有 `run_trial.epsilon=0.1` 改变残差停止阈值，不启用 `volume_threshold`，不构造认证内域。两方法均从相同方向 MP2 初始化及零条割开始，首次得到可靠全局 `bound <= epsilon` 时停止；候选值不代替全局上界。共享初始化、各自搜索、事后扫描分别计时。

新增 `experiments/fourbus_threshold_report.py` 只做事后比较。全部 861 条射线重新求完整 SOCP 参考边界，x/y 在预算内自由；两个外域使用相同负荷比例。百分比的分母固定为该射线的 `reference_total`，不以外域容量或坐标上界作分母。所有参考容量均须为正。

| 数学量 / 定义 | 固定代码映射 | 单位、形状与边界 |
|---|---|---|
| e(w)=100[ρ_outer(w)−ρ_ref(w)]/ρ_ref(w) | [make_threshold_report.radial_difference_percent](../experiments/fourbus_threshold_report.py) | %，每方法 `(861,)`；微小负数保留，不截断数据 |
| 射线相对误差均值、最大值、95 分位 | [make_threshold_report:mean_radial_difference_percent](../experiments/fourbus_threshold_report.py)、[make_threshold_report:max_radial_difference_percent](../experiments/fourbus_threshold_report.py)、[make_threshold_report:p95_radial_difference_percent](../experiments/fourbus_threshold_report.py) | %；均值为 861 条离散射线等权均值，不声称按球面面积均匀加权 |
| V_mesh=Σ_T abs(det(p_T0,p_T1,p_T2))/6 | [make_threshold_report.mesh_volume](../experiments/fourbus_threshold_report.py) | kW³；以全部射线前沿点和原点组成三角锥的体积和，仅为插值表面的体积估计 |
| 三角网格体积相对高估 | [make_threshold_report:mesh_volume_difference_percent](../experiments/fourbus_threshold_report.py) | 100(V_outer,mesh−V_ref,mesh)/V_ref,mesh，%；不是连续真实域体积的精确误差，也不是认证内域或停止判据 |
| 含共享初始化的总时间 | [make_threshold_report:total_seconds](../experiments/fourbus_threshold_report.py) | `seconds+initial_seconds`，秒；独立扫描、审割和绘图均不计入算法耗时 |

两张图统一坐标范围、视角及百分比色标；a 为三维射线前沿的三角网格插值比较，b 为完整负荷比例三角形上的相对容量误差。跨网架的点只用于明确标注的插值显示，不取凸包、不作为认证区域。原始点、参考求解器界、实际停止轨迹及全部相对误差保存以供复核。

## 30. 主线与角点分块方法的同案例产出比较

`experiments/fourbus_mainline_comparison.py` 直接调用未改动的 `main.build_continuous_region` 与 `run_trial`。预算 20,000 元、4 个求解器线程、数值库 1 线程；主线使用既有 `tau=0.002`、有限预算 `auto -> light`，分块法使用 `epsilon=0.1`。每种方法独立从零初始化，重复 3 次，分别报告含初始化的耗时、停止证书及相同 861 条参考射线的误差。两种停止标准不等价，不能把时间比称为同精度加速比。

主线比较的是 `RegionState.finish` 实际返回的认证内域与全局外包络，不能将已知网架的条件割外域当成全局外包络，也不能仅用全部割重建的较松联合外域代替主线最终输出。分块法仍只报告联合割投影外域。

| 数学量 / 定义 | 固定代码映射 | 单位、形状与用途 |
|---|---|---|
| ρ_region(w)=max{ρ>=0:ρw 属于给定多面体并集} | [region_radial_capacity](../experiments/fourbus_mainline_comparison.py) | kW，逐射线先计算各多面体的可行区间，再取最大上端点；保留不同方案并集，不混合顶点取凸包 |
| 主线认证内域的逐射线总负荷 | [run_mainline_comparison.inner_total](../experiments/fourbus_mainline_comparison.py) | kW，`(861,)`；与已有 `outer_total/reference_total` 的边界总量意义一致 |
| 100(ρ_ref−ρ_inner)/ρ_ref 的最大、平均值 | [run_mainline_comparison:max_inner_difference_percent](../experiments/fourbus_mainline_comparison.py)、[run_mainline_comparison:mean_inner_difference_percent](../experiments/fourbus_mainline_comparison.py) | %；主线认证内域低估，分块法没有相应内域，不编造该指标 |
| 实际含初始化算法时间 | [run_mainline_comparison:total_seconds](../experiments/fourbus_mainline_comparison.py) | 秒；主线直接取 `timing.total_seconds`，分块取 `seconds+initial_seconds`；事后几何射线分析不计入 |

其余射线容量差、相对误差、三角网格体积估计继续采用第 26、29 节名称和定义。参考射线复用同一未变化物理模型的第 29 节独立求解数据；所有重复运行均计算误差，保存逐次结果和中位数/范围。图形使用其中耗时位于中位数的实际运行展示三维域，误差曲线保留全部重复运行，耗时图显示全部观测值。

## 31. 主线顶点 SP 违反量排序实验

`experiments/fourbus_vertex_priority.py` 独立实现主线的 SOCP 循环，复用物理约束、联合割代数、同网架见证支撑、几何更新及 `RemainingRegionModel`。生产 `main.py/model.py/region.py` 不变。普通候选仍为已知网架收缩外域顶点中未被认证内域并集覆盖的 `(x, point)`，其中 `point` 是归一化坐标；评分在实际负荷 `power=point*bounds` 上求解。

每轮固定候选集合，对所有候选求原始 SOCP 最小 eta；先收集所有可行证书，再在仍未被并集覆盖的不可行候选中按 eta 最大、实际总负荷最大、x/point 字典序依次排序，只加入获选割。其他结果以精确 `(x, power)` 元组缓存，后续相同参数复用；不对浮点坐标四舍五入。SOCP 后的支撑平面 LP 目标不作为评分。活动全局见证仍按主线的同方案支撑顺序处理。全局检查间隔按实际新求解的 SP 次数累计，结束必须取得原来的覆盖证书，绝不使用评分阈值代替。

| 数学量 / 定义 | 固定代码映射 | 单位、形状与用途 |
|---|---|---|
| eta*(x,p)，原始 SOCP 最优违反量 | [ScoredSubProblem.solve:eta](../experiments/fourbus_vertex_priority.py) | 原 SP 残差尺度；在替换锥为支撑平面前保存，用于候选排序 |
| 已计算的固定参数 SP 结果 | [ScoredSubProblem.cache](../experiments/fourbus_vertex_priority.py) | 精确 `(tuple(x), tuple(power))` 为键，值含 eta/feasible/state/cut；割布局、运行状态及可行标准不变 |
| 缓存命中次数 | [ScoredSubProblem.cache_hits](../experiments/fourbus_vertex_priority.py) | 整数；不计入实际 SP 求解次数 `calls` |
| 候选评分轮次 | [build_vertex_priority_region.scoring_rounds](../experiments/fourbus_vertex_priority.py) | 一次固定普通候选集的批量评分算一轮 |
| 全局剩余域调用次数 | [build_vertex_priority_region.global_search](../experiments/fourbus_vertex_priority.py) | 次；仍使用 light 模式、原 tau 和 GEOMETRY_TOL |
| 每轮待评分候选数 | [build_vertex_priority_region:candidate_count](../experiments/fourbus_vertex_priority.py) | 次；包括复用缓存的候选，不能当成新 SP 次数 |
| 实验联合割与计算轨迹 | [build_vertex_priority_region:cuts](../experiments/fourbus_vertex_priority.py)、[build_vertex_priority_region:trace](../experiments/fourbus_vertex_priority.py) | 割保持 `[alpha,*beta,*delta]`；轨迹保存选点/评分/全局覆盖检查，不保存运行状态向量 |

对照直接调用原主线；两方法均为 FourBus、预算 20,000 元、tau=0.002、求解器 4 线程、数值库 1 线程，独立初始化、重复 3 次。含初始化总耗时只统计算法，不含事后参考射线核验与绘图。比较实际返回的认证内域与最终外包络；第 30 节的逐射线容量、内域低估和外域高估定义不变。全部 861 条 SOCP 参考射线仅用于事后评价，不参与候选评分或停止。

## 32. 全局搜索间隔与低违反量触发实验

仅为第 31 节独立实验增加调度参数；默认值保持原算法及证书。`experiments/fourbus_global_schedule.py` 对比间隔 96、32、16，以及每个间隔叠加低违反量触发的配置。全部配置仍须取得原来的全局覆盖证书，eta 阈值只决定何时重新搜索，不接受 eta<0.01 为物理可行或构域完成。

| 数学量 / 定义 | 固定代码映射 | 单位、形状与用途 |
|---|---|---|
| 两次全局搜索之间的实际新增 SP 次数阈值 | [build_vertex_priority_region.refinement_checks](../experiments/fourbus_vertex_priority.py) | 默认 `REFINEMENT_CHECKS=96`；实验取 96/32/16；检查发生于评分批次或当前见证支撑处理结束后，不是逐个 SP 的强制中断 |
| 普通候选最大 eta 的全局搜索触发阈值 | [build_vertex_priority_region.eta_trigger](../experiments/fourbus_vertex_priority.py) | 默认 `None` 禁用；启用为 0.01，使用原始 SP 残差尺度，不是 kW 或几何误差 |
| 是否已请求下一轮全局搜索 | [build_vertex_priority_region.global_requested](../experiments/fourbus_vertex_priority.py) | 布尔调度状态；仅普通评分轮选中点的 eta 严格小于 eta_trigger 时设置，全局搜索后清除 |
| 全局搜索本次触发原因 | [build_vertex_priority_region:trigger](../experiments/fourbus_vertex_priority.py) | 字符串列表：`empty` 为普通候选空，`interval` 为新增 SP 达阈值，`eta` 为低违反量请求；允许多个原因同时成立 |
| 首次全局搜索前的实际 SP 数 | [run_global_schedule_comparison:first_global_sp](../experiments/fourbus_global_schedule.py) | 次；来自首条全局轨迹的累计 sp，用于核实批处理后的实际触发时刻 |
| 各原因实际出现次数 | [run_global_schedule_comparison:trigger_counts](../experiments/fourbus_global_schedule.py) | 按 empty/interval/eta 计数；原因非互斥，不能相加替代 global_search |

低违反量判定为：先收入本批全部可行证书，再在仍未覆盖的不可行候选集合上取最大 eta；这个最大值小于阈值时，先应用已经取得的获选有效割，然后下一轮请求全局搜索。不能因某一个非最大候选 eta 很小而跳过其他严重违反者；活动全局见证的支撑 SP 也不单独触发该规则，避免打断同网架覆盖进展。无剩余不可行候选时沿用候选空触发。

迁移范围：`build_vertex_priority_region` 只增加可选参数及结果元数据 `refinement_checks/eta_trigger`、轨迹 `trigger`；原有字段含义、默认选点、数值容差、联合割和最终包络不变。旧实验记录不补写未知触发原因。比较继续采用第 30、31 节的耗时和逐射线误差字段，普通候选最大 eta 不命名为全局 R_k 上界。

## 33. 主线顶点评分与原生窗口迁移（2026-09-27）

本节覆盖前述历史版本的入口、记录与展示约定，历史实验数据不改写。

- `build_continuous_region` 的普通候选按原始 SP 最优违反量 `eta` 降序选割；先收入该批全部可行点，再移除已被任一内域覆盖的候选。内部见证仍补同一网架支撑，不能跨网架取凸包。
- [REFINEMENT_CHECKS](../continuous.py) 改为 32，计数仅含实际新增 SP 求解，缓存命中不计数；评分批次及活动见证完成后检查。生产算法没有 `eta_trigger`，没有 0.01 调度或停止阈值；终止仍须全局覆盖证书。
- [SubProblem.solve:eta](../model.py) 新增为原始 SP（替换锥之前）的最优违反量，无量纲；已有 `feasible/state/cut` 不改名。主线以 `(tuple(x), tuple(power))` 为精确缓存键，不对负荷取整；缓存属于本次构域，不落盘。
- [FourBus.__init__.load_nodes](../Network/four_bus_five_corridor.py) 新增可选负荷节点参数，类默认仍为 `(1,2,3)`，以保留既有三维实验。`main.main` 显式传 `(1,2)`；节点 3 的有功、无功均固定为原始值 0，仍是必须连接的节点。这是二维切片，不是将第三负荷投影消去。
- `main.run` 改为单预算 SOCP 入口；旧 `budgets/recompute/keep_ui/step_by_step`、多网络选择常量与 `BenchmarkResult` 结果入口退役。`build_continuous_region` 保留物理模型入口和历史调用所需的 `method/residual_mode/cuts/progress/clock`，默认主线只传 `socp/light`。旧勘察仍可单独调用 `survey.py`。
- `max_point/max_choice/max_cost` 退役：初始化点、网架与费用已经保存在监视器的几何帧中，不再在最终结果重复保留。`max_total/max_total_bound/axis_bounds/coverage_bound/status/counts/timing/inner/outer` 名称与单位不变。
- [RunMonitor](../monitor.py) 接管计时、精简阶段帧、原生窗口和回放。构域通过 `begin/initializing/seed/selecting/global_start/global_end/sp_start/sp_end/updated/finish` 跟踪。`progress` 函数回调仍由监视器桥接，算法不组装文字报告。
- 新运行只保存一个 `monitor.json.gz`：基础坐标信息和增量 `history` 帧；最终帧含结果和扫描数组。退出或失败由监视器保存同一文件，不生成 HTML、JSONL、报告或单独结果副本。旧 HTML 回放不作为新窗口输入。
- 橙色点附其实际 SP 调用序号 `number`；采用缓存割时仍显示原序号，顶部 `sp` 仅表示累计新增求解次数，两者不能混淆。
- 总图内域是已认证网架内域的几何并集。全局认证前，安全全局外包络仍为 MP2 的轴向及总量界；不能把已知网架外域并集误称为全体网架的外域。每网架面板显示自身联合割条件外域；全局认证后总图外包络采用 `RegionState.finish(True)`。

| 新增量 | 代码映射 | 定义 |
|---|---|---|
| 单次运行预算 | [BUDGET](../main.py)、[run.budget](../main.py) | FourBus 默认 20000 元 |
| 独立 SOCP 扫描 | [validate_socp_region](../vertify.py) | 固定每个二维网格中心的 `p`，完整 MISOCP 自由选择全部合法 `x,y`；不用构域割或已知方案列表 |
| 扫描标签 | [validate_socp_region:states](../vertify.py) | `(divisions, divisions)`，1 可行，-1 不可行；不是 AC 真值 |
| 多余 / 遗漏百分比 | [RunMonitor.validation:fr_percent](../monitor.py)、[RunMonitor.validation:mr_percent](../monitor.py) | 保留原分母：FR = 多余格点/算法内域格点；MR = 遗漏格点/扫描可行格点，均乘 100；空分母为 None |
| 全局 / SP 点 | [RunMonitor.global_end](../monitor.py)、[RunMonitor.sp_start](../monitor.py) | 实际 kW 坐标分别记录，不互相覆盖；红色菱形 / 橙色圆点；SP 支撑点可不同于全局见证 |

扫描精度由 `DIVISIONS` 决定；网格百分比是离散估计，不能把 0% 解读为连续域完全无误差。最终仍保留径向精度 `tau` 与全局覆盖上界。

几何数值迁移：`polytope_vertices/clip_polytope/halfspaces/polytope_volume` 在仿射独立方向上用坐标列范数作一次确定的缩放后调用 Qhull；没有失败后重试。顶点仍返回原始坐标，体积乘回缩放行列式，半空间法向量换回原坐标后重新归一化为单位向量。原 `GEOMETRY_TOL` 的归一化坐标含义、维数判定和联合割符号保持不变。此预条件化保留凸包的仿射结构，不把极薄非零体积抹成零。

## 34. 二维割线回放与独立校验面板（2026-09-27）

- 联合割仍为 `[alpha,*beta,*delta]`，保留侧为 `alpha+beta@p+delta@x >= 0`。二维网架截线为 `beta@p+(alpha+delta@x)=0`，p 的单位仍为 kW；同一条割在各网架中分别代入该网架 x。负荷系数为零时不存在二维直线，不虚构切割线。
- [RunMonitor._geometry:x](../monitor.py) 在网架记录中增加二进制型号向量，顺序、长度与 `Network.type_keys` 相同，仅用于准确代入联合割。
- [RunMonitor.updated:cut_history](../monitor.py) 按割序号增量记录 `cut` 原始系数及来源 `scheme`，每条割只记录一次。初始化传入的割由 `begin` 记录，来源未知记为 None。绘图函数 [_cut_segment](../monitor.py) 只将该条件直线裁到当前显示框，不改变算法几何或容差。
- 当前割在网架子图为紫色实线，历史割为淡紫虚线；加割帧中浅紫色表示该网架前后条件外域之差。总图只叠加当前割在来源网架下的截线，明确注明网架，绝不把该直线用作全局负荷半空间。已知网架条件外域并集与安全全局外包络分别显示。
- 回放格式迁移到 version=4：`history` 只包含构域帧；[RunMonitor.validation_state](../monitor.py) 单独保存当前扫描进度和最终校验结果，仍写在同一个 `monitor.json.gz`。扫描和完成回调不增加时间轴帧，最终对比不随回放位置变化。旧 version=3 回放读取时分离扫描帧；旧文件未记录的割系数不补造。
- [main.case](../main.py) 支持 `fourbus` / `case33`。Case33 使用既有 `Case33(load_nodes=(18,25))`，默认 4 个升级候选、预算 2 相对投资单位；节点 33 及其余非动态负荷保持原始有功/无功。独立扫描默认 20×20，FourBus 仍为 80×80。两入口仅计算二维切片，沿用 SOCP、32 次新增 SP 调度及相同全局覆盖证书。
- 入口迁移（2026-09-28）：恢复单一案例选择项 [NETWORK](../main.py)，取类 `FourBus` 或 `Case33`，默认 `FourBus`。`main(case=None)` 及未指定 `--case` 的脚本启动使用 `NETWORK`；显式 `case` / `--case` 覆盖它，沿用上条的二维配置与结果路径。
- 初始化 MP2 使用本次构域的剩余时限，不再额外受 20 秒的单次 MP 默认上限截断。FourBus 总构域时限仍为 300 秒，Case33 为 900 秒；SP 和全局搜索仍受原单次上限及剩余总时限约束。超时仍直接报告失败，不接收超时解、不重试、不降低求解精度。
- 网架子图按每页 4 个显示，并复用最多 4 个绘图区；改变选中网架时自动定位页面。分页只影响界面，不删减记录中的网架、割或历史帧。未取得构域证书的运行保留完整切割轨迹，校验面板明确显示未完成，不填造百分比。

## 35. 径向精度默认值调整（2026-09-28）

`REGION_TAU` 从 0.002 调到 0.005，含义仍为径向收缩比例，不增加距离参数或修改 SP 的 eta / 物理接受容差。普通候选外域顶点 p 满足 `(1-tau)*p` 落入某个已认证内域时跳过 SP；全局剩余域模型及最终 `RegionState.finish` 外包络继续使用同一 tau。终止仍须全局上界证书，不能仅凭所有已知顶点被跳过就结束。

被跳过的薄层没有获得运行可行性证书，不加入 `inner`；最终 `outer` 按原公式包含径向精度及几何数值容差。结果中的 `tau` 保留实际使用值，旧回放和历史实验中的 0.002 不改写。该调整以降低边界分辨率减少 SP 次数，不宣称同精度加速，也不把 0.5% 解释为面积遗漏率。

## 36. FourBus 已有网架优先实验（2026-09-28）

`experiments/fourbus_known_first.py` 独立实现二维 SOCP 调度，不修改生产主线。令 `I` 为已有网架认证内域的并集，`I_tol` 为每个内域逐面按既有 `GEOMETRY_TOL` 放宽后再取并集。进入全局搜索前要求每个已有网架满足 `(1-tau)O_x ⊆ I_tol`。普通候选仍为未覆盖的收缩外域顶点；顶点全部被覆盖后，用半空间逐次相减，把剩余域分解为凸块，以凸块顶点均值选取内部未覆盖见证。该均值只作几何见证，SP 复用 `RegionState.witness_support` 选择尚未被当前网架自身认证的最大总负荷支撑点，避免逐次认证中心点只能渐进填补孔洞。不能仅以顶点覆盖、剩余面积较小或 SP 次数判定已有网架完成。

局部 SP 固定 x,p，复用原 SOCP、eta 判据及联合割；每批先收入全部可行点，再应用 eta 最大的一条割。普通顶点已被并集覆盖时跳过该候选；内部见证的支撑点按当前网架自身内域判断，不能把其他网架的认证点直接纳入当前网架凸包，因此其失败割也仍须更新条件外域。几何相减只负责选点及局部覆盖判断，不产生物理可行证书。全局搜索固定采用 `physical`，允许全部合法 x,p,y 自由变化；返回的完整可行点立即加入其自身网架内域，再恢复局部探索。终止仍使用 `RemainingRegionModel` 的全局覆盖上界，tau、GEOMETRY_TOL、PLANNING_TOL 和割布局均保持不变。

| 数学量 / 定义 | 固定代码映射 | 单位、形状及用途 |
|---|---|---|
| 一个凸候选域减去内域并集后的凸块 | [remaining_cells](../experiments/fourbus_known_first.py) | 归一化二维顶点数组列表；逐个内域、逐个面分割，保留低维剩余块，不用面积阈值删除 |
| 已有网架下一批待认证候选 | [known_candidates](../experiments/fourbus_known_first.py) | `(x, point)` 列表；point 为归一化坐标，SP 前乘 bounds；返回的 `interior` 区分顶点批次和内部剩余块批次 |
| 已有网架完整剩余检查次数 | [build_known_first_region:known_checks](../experiments/fourbus_known_first.py) | 次；顶点批次为空后执行一次全域几何检查，不是全局优化次数 |
| 内部剩余候选的实际新增 SP 数 | [build_known_first_region:hole_checks](../experiments/fourbus_known_first.py) | 次；包含孔洞或顶点之间的缺口，缓存命中不计数 |
| 构域期间登记的网架数 | [build_known_first_region:schemes](../experiments/fourbus_known_first.py) | 整数；包括所有登记记录，不按非空内域反推 |
| 重复运行的最简比较记录 | [run_comparison.comparison](../experiments/fourbus_known_first.py) | 保存每次运行的 `total_seconds/counts/coverage_bound/fr_percent/mr_percent/volume_gap`；volume_gap 在二维使用面积测度，仍是比例 |

对照直接调用当前 `main.build_continuous_region`：32 次 SP 调度、相同 physical 全局模型、预算 20000 元、负荷节点 (1,2)、tau=0.005、每种方法独立初始化、4 个求解器线程、数值库 1 线程，默认交替顺序重复 3 次。构域时间含初始化和监视记录，不含独立扫描、事后审计或绘图。双方共同使用一次 80×80 完整 SOCP 扫描，不限制参考网架、不使用构域割。每种方法仅保存耗时中位数对应的一份原生 monitor 回放；比较记录写入其 validation_state，不另建报告或逐轮日志。

## 37. 二维 FourBus：双线性与负荷 LP 互补重写对照

`experiments/fourbus_outer_kkt.py` 为独立实验，保持主线及原 `GlobalViolation` 不变。负荷节点为 (1,2)，节点 3 负荷固定为 0；预算 20000 元，原始 SOCP SP、完整锥对偶及全部二元建设变量不变。两方法求同一个 R_k=max_(x,p in O_k) eta*(x,p)，默认仅当有效全局上界 <=0.1 时认证，不用候选值或面积误差停止。

负荷外域局部记为 `outer_p @ p <= outer_rhs + outer_x @ x`（这些记号不替换第 4 节的网络 H/T）。固定 x 与 SP 对偶乘子后，负荷子问题是线性规划。其乘子 lambda>=0 满足 `outer_p.T @ lambda = sp_p.T @ row_dual`，并与每行外域余量互补。最优性给出 `(sp_p.T @ row_dual) @ p = (outer_rhs + outer_x @ x) @ lambda`。余量非负，二元变量为 0 时令 lambda=0、为 1 时令余量=0；用指示约束表达，不设置任意乘子大 M。新增的 x*lambda 同样用二元指示约束精确表达。保留原完整 Lorentz 锥，设置 NonConvex=0，不作网架连续松弛或锥多面体替换。

| 数学量 / 定义 | 固定代码映射 | 单位、形状与用途 |
|---|---|---|
| 负荷外域各行的 lambda | [KktViolation.outer_dual](../experiments/fourbus_outer_kkt.py) | 非负连续变量列表，按加入的上下界、总量界、联合割顺序 |
| 负荷外域的系数行 | [KktViolation.outer_p](../experiments/fourbus_outer_kkt.py)、[KktViolation.outer_x](../experiments/fourbus_outer_kkt.py) | 每行分别为 p、右端 x 的系数；用于计算互补乘子的有效有限界 |
| 各行互补选择 | [KktViolation.outer_active](../experiments/fourbus_outer_kkt.py) | 二元变量列表；1 允许非零乘子并强制该行紧约束，不代表网架 |
| 负荷 LP 对偶平衡行 | [KktViolation.load_stationarity](../experiments/fourbus_outer_kkt.py) | d 条线性等式；每次加割同步新增 lambda 的系数 |
| outer_x.T @ lambda | [KktViolation.outer_choice_dual](../experiments/fourbus_outer_kkt.py) | `(t,)` 自由连续变量，保留完整建设型号索引 |
| x_j*(outer_x.T @ lambda)_j | [KktViolation.outer_choice_term](../experiments/fourbus_outer_kkt.py) | `(t,)`；与原 `choice_term` 分开，原项仍只指 SP 的 sp_x 贡献 |
| 重写后的线性目标约束 | [KktViolation.dual_objective](../experiments/fourbus_outer_kkt.py) | R <= sum(choice_term)+sum(outer_choice_term)+outer_rhs@lambda-sp_rhs@mu-cone_constant@s |
| 单次构域耗时 | [run_trial:seconds](../experiments/fourbus_outer_kkt.py) | 秒，含本方法建模、G、SP 与加割，不含共享 MP2 初始化、事后扫描和绘图 |
| 含初始化总耗时 | [run_trial:total_seconds](../experiments/fourbus_outer_kkt.py) | seconds+initial_seconds；共同初始化时间为两方法各计一次 |
| 外域扫描多余 / 遗漏百分比 | [compare_reference:fr_percent](../experiments/fourbus_outer_kkt.py)、[compare_reference:mr_percent](../experiments/fourbus_outer_kkt.py) | FR=外域中的扫描不可行点/外域格点数；MR=外域外的扫描可行点/扫描可行点数，均乘100；本节算法集合是外域，与主线内域指标须明确区分 |

每次加割必须同时更新负荷可行域、负荷 LP 驻点条件、互补选择和替换目标，遗漏其中任何一项均不等价。G 的候选下界与全局界继续复用原实现及其数学接受条件；超时、无反例或数值异常不能写成构域完成。实验默认交替顺序重复 3 次，4 个求解线程、数值库 1 线程。共用同一组初始方向界且各自从零条割开始。独立二维网格参考固定 p、自由选择全部合法 x/y，不用构域割；只作 SOCP 离散参考，非 AC 真值。事后枚举预算内网架仅用于准确绘制小算例的外域投影并集与检查，不参与 G 或 SP 的选点。结果、轨迹及扫描保留在一份 `comparison.json.gz`，图另存 PNG/SVG/PDF；耗时中位数那次实际运行用于各方法的区域图。

本次数值回归与正式比较均采用 4 个 G 求解线程：用全部 17 个预算内网架及其外域顶点的原始 SP 独立核验初始和加割后的全局最大值，并检查联合割在不限预算的完整物理域上有效。未加有限乘子界的原型曾在继续加割后出现 Gurobi SUBOPTIMAL（状态 13），仍直接报错，不接收该状态、不自动重试或切换算法。最终版本统一施加下述有效界；R_k<=0.01 的四线程比较及额外单线程核查均达标。初始 R_k<=0.1 记录的 `settings.kkt_dual_bounds=False`，新版本为 True；旧记录的时间和割不改写。按用户追加要求，R_k<=0.01、每方法最多 100 秒的正式对照为各一次运行，保存在 `results/fourbus_outer_kkt_001`，不能称为三次中位数。

互补乘子有限界：二维负荷 LP 的对偶可行集为 `{lambda>=0: outer_p.T@lambda=theta}`。非空有界的原负荷 LP 存在最优基本对偶解，其非零乘子可限制在两个线性无关行上。由原 eta 归一化，theta 属于 0 与 `sp_p[i]/(-sp_eta[i])`（sp_eta[i]<0）的凸包；其余行的 sp_p 必须为零。对每个二行基及这些 theta 极点计算对应基本解分量，逐个 lambda 取所有基/极点的非负最大值，即有至少一个最优对偶解满足的有效上界。`_bound_outer_duals` 用实际浮点系数的精确有理数表示计算二阶行列式和商，再向上舍入为浮点上界；加割后重算所有旧、新乘子界，不能沿用较小的旧界。outer_choice_dual/outer_choice_term 的界由这些非负乘子界及 outer_x 的符号作区间传播。该界不依赖枚举网架，不改变原负荷空间或 eta，且只使用两行基组合而非全部负荷顶点。它在建模时统一施加，不是求解失败后的重试策略。

## 38. KKT 未认证区域优先实验

`fourbus_outer_kkt.py --coverage` 比较原 KKT 与认证并集排除版本，主线不变。I_x 为同一网架下完整 SOCP 认证点的凸包，I 为这些凸包的并集，严禁跨网架取凸包。初始化 MP2 点直接复用；每个首次登记网架补充两个轴和总负荷方向的固定网架 SOCP 支撑点，其求解时间计入构域时间，不预先枚举全部网架。FourBus 零固定负荷下，原点由零潮流、单位电压解析认证。

为避免严格不等式与重复边界点，归一化坐标中以 `I_x + [-coverage_pad,coverage_pad]^2` 作搜索排除域，真实 inner 仍只保留认证凸包。G 搜索联合割外域中每个排除多边形之外（含边界）的网架—负荷组合。排除域以 ConvexHull 的完整半空间表示；每个多边形至少选择一个外侧面。每个选择 s 对应 `a@p <= rhs + selector_rhs*s`，M 仅由已知负荷盒精确计算。固定 x 和所有 s 后仍为负荷 LP；全部排除行同时进入原始约束、KKT 驻点、互补和替换目标，s*lambda 用指示约束精确表达，乘子界按全部新行重算。

| 新增量 | 固定代码映射 | 定义与单位 |
|---|---|---|
| 认证域搜索扩边 δ | [COVERAGE_PAD](../experiments/fourbus_outer_kkt.py) | 1e-5；公共评价箱归一化坐标；不是新增可行点 |
| 面选择与其右端系数 | [KktViolation._add_outer_row.selector](../experiments/fourbus_outer_kkt.py)、[KktViolation._add_outer_row.selector_rhs](../experiments/fourbus_outer_kkt.py) | 二元变量与常数；原方法不传此参数 |
| s*lambda 精确辅助项 | [KktViolation.selector_terms](../experiments/fourbus_outer_kkt.py) | `(lambda 行索引, product)` 对列表，乘子有限界同步传播到 product |
| 排除多边形的半空间 | [exclusion_halfspaces](../experiments/fourbus_outer_kkt.py) | 输入为归一化认证凸包，输出 [F,g]；内侧 F@xi+g<=0 |
| 排除并集后的 G | [UncoveredKktViolation](../experiments/fourbus_outer_kkt.py) | bound 为未认证区域最大配对违反量上界；不是原全域 R_k |
| 固定负荷的全网架认证 | [certify_load](../experiments/fourbus_outer_kkt.py) | 固定 p，自由 x,y，最小 eta；不用实验割限制 x；返回 eta/bound/x/feasible/status |
| 认证薄层的违反量上界 | [run_coverage:padding_bound](../experiments/fourbus_outer_kkt.py) | δ max_i(sum_j(abs(sp_p_ij)*axis_bounds_j)/(-sp_eta_i)) + PLANNING_TOL；仅计 sp_eta_i<0 的行 |
| 整个负荷投影的违反量上界 | [run_coverage:union_bound](../experiments/fourbus_outer_kkt.py) | max(未排除部分的 bound, padding_bound)，上界于 max_p min_x eta*(x,p)，不冒充 max_(x,p) eta* |
| 认证凸包结果 | [run_coverage:inner](../experiments/fourbus_outer_kkt.py) | 列表；每项 x 为建设向量、vertices 为 kW 的同方案认证凸包顶点 |

薄层界来自：保持一个认证点的可行运行 y 不变，只增加 eta 即可容纳盒内负荷扰动；所有非松弛行的 sp_p 必须为零。因此排除薄层没有获得精确可行认证，但其 union 违反量有显式界。只有 union_bound<=epsilon 或等价有效全局证书才能停止。新的 bound_scope=`uncovered_pairs`；旧结果 bound 仍按第 37 节解释，未迁移或改写旧数据。

每轮 G 后，固定它的 x,p 求原 SP 并加联合割，再固定该 p 求全网架认证；发现可行方案后扩充其自身内域并重建 G。全网架认证证明 eta>0 只说明该点不可行，不能据此任意删除邻域；后续仍由有效联合割收紧各条件外域。超时未获证保留未完成状态。最终图与 FR/MR 仍以联合割外域在全部预算内网架上的投影并集计算，不扣除认证域，不用扫描指导选点。

FourBus 正式对照保存在 `results/fourbus_outer_coverage_001`：epsilon=0.01、每次构域最多 100 秒、4 线程、交替顺序重复 3 次、共用 80×80 SOCP 扫描。原 KKT / 新版本耗时中位数为 10.7531 / 3.7930 秒，G 为 17 / 13 次，SP 与割为 16 / 12 次；新版本另计入 12 次全网架认证和 6 次固定网架方向支撑，3 次成功认证均切换了候选网架，最终登记网架仍为 2 个。多余网格点为 14 / 15 个，FR 为 0.3764% / 0.4032%，均未发现遗漏。新版本 remaining bound=0.00596354、padding_bound=6.35238e-6；事后枚举得到其未限制认证域的全配对最大值仍为 0.02054188，这不违反其 union 证书，不能改写成原 R_k<=0.01。上述网架枚举、扫描及审计均在计时结束后进行。

## 39. Case33 限定开关与双算法回放

Case33 仅保留原始线路型号，不再提供升级接口。根节点为 1；允许改变状态的支路为 21-8、7-8、22-12、11-12、9-15、33-18、25-29，其他支路固定为原始状态。完整连通和树边数约束保持不变。旧升级算例只作为 tests/legacy_case33.py 的历史回归夹具，不用于当前入口。

| 数学量 / 定义 | 固定代码映射 | 单位与用途 |
|---|---|---|
| 支路是否允许改变原始状态 | [Corridor.switchable](../Network/__init__.py) | bool；默认 True 保持其他网络原行为；False 时 MP 固定初始型号向量 |
| B，最多开合变动次数 | [Case33.switch_budget](../Network/case33bw.py) | 默认 7；Network/case33bw.py 唯一默认入口 |
| c0，预算仿射常数 | [Network.cost_offset](../Network/__init__.py) | 默认 0；Case33 为允许打开的常闭支路数 2 |
| c(x)=c0+c@x | Network.cost、TypeParameters.investment_cost、OperatingTree.cost | 显式迁移：Case33 的系数 c 对常闭可变支路为 -1、常开可变支路为 +1、其余为 0；总成本即相对初始状态的 Hamming 距离。其他算例仍为原增量投资，c0=0 |
| 二维动态负荷 | Case33.load_nodes | 主入口及对照均为 (18,25)，其余节点保持原始负荷 |
| KKT 未认证区域上界 | [build_kkt_region:union_bound](../experiments/case33_compare.py) | 沿用第 38 节语义，不能称为全配对 R_k；epsilon=0.01 |
| 校验使用的算法集合 | [RunMonitor.validation.region_key](../monitor.py) | inner 或 outer；同一套 MR/FR 定义，明确标记集合，不混淆内外域 |
| 原有内域与外域的两套网格指标 | [RunMonitor.validation:metrics](../monitor.py) | 键 inner/outer；各自 MR/FR 保持原分母，主显示集合由 region_key 指定 |
| 固定 p、自由 x/y 的认证点 | [RunMonitor.certification_start.power](../monitor.py) | kW、(2,)；回放独立阶段，不改写固定方案 SP 的点 |
| 开关算例的联合割投影 | [projected_outer](../experiments/case33_compare.py) | 每个预算内网架的条件外域顶点，kW；只作事后绘图和误差统计，不参与 G |

主线沿用 tau=0.005 的几何覆盖证书，新 KKT 使用未认证区域残差上界 <=0.01；两者独立初始化、各计最多 1000 秒，时间包含初始化和记录，不含事后参考扫描。达到时限只保存当时已认证内域和有效外域，不能记成认证完成。Case33 新网架的原点不默认可行，增加完整 SOCP 的负总负荷方向支撑以取得下侧认证点。G 的零潮流解析 eta 上界包含全部固定和动态有功/无功负荷。

回放继续使用 version=4 的单份增量文件；全网架固定负荷认证、固定网架支撑单独标明阶段。双窗口同步按各自帧序号前进一步，较短轨迹停在末帧；这不是两算法数学步骤一一对应。预算、扫描结果及比较指标均存于两份 monitor 文件，不新增文字报告。

数值与超时约定：G 的 `candidate_lower` 是经过残差修正的保守下界，等于 0 不等于该候选 SP 的最优 eta 等于 0；只要 G 返回候选，仍调用原始 SP 判定，不能据此重复求解同一 G。主线全局搜索可用剩余总时限，不再受旧 120 秒单次上限截断。求解器 TIME_LIMIT 显式抛出 TimeoutError，实验入口保存截止状态，不接收超时运行解、不重求 SP、不放宽容差。

当前 1000 秒上限、epsilon=0.01、4 线程、80×80 共同扫描的实际单次结果：主线 0.9364 秒、15 次 SP、4 条割、3 次全局检查；KKT 排除版 1.5847 秒、8 次固定网架支撑、首次 G 上界 0.00363333 即停止，0 次 SP/割。内域遗漏率分别 0.1599% / 0.0640%，外域多余率分别 0.4299% / 2.0367%；两者内域多余和外域遗漏均为 0 个网格点。两种停止标准不同，不是同精度速度比较。结果位于 results/case33_compare 的两份原生回放；扫描是 SOCP 离散参考，不是连续域或 AC 真值证明。
