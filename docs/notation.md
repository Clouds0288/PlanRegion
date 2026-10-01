# 数学符号与代码变量契约

本文只描述当前主线和完整 physical 支持面实验。相同数学量必须保留登记名称、单位、索引、状态布局和结果键；不能通过删登记绕过检查。新增量先登记，语义或输出格式变化必须显式迁移并同步调用方与回归。

## 维数、索引与单位

| 数学量 | 固定名称 | 单位与形状 |
|---|---|---|
| 非根节点数、候选走廊数、型号数、功率维数 | `n / n_corridors / n_types / len(load_nodes)`，公式 `n/m/t/d` | `d=2/3`；真实节点 ID 不等于数组位置 |
| 建设向量 x | `choices[e,k]`、`x`、答案键 `x` | 二元 `(t,)`；唯一顺序为 `network.type_keys` |
| 走廊状态 z=Hx | `network.corridor_types @ x` | `(m,)`；H 为 `(m,t)` |
| 独立功率 p | `power`、`problem.power`、答案键 `p` | kW，`(d,)`；顺序为 `load_nodes` |
| 运行状态 y | `state` | `(3*t+n+2*m,)`；不含 SP 的 eta |
| 功率基值、电压基值 | `base / voltage_kv` | kVA / kV；不是自动等于电源容量 |
| 电阻、电抗 | `r / reactance` | p.u.，`(t,)`；电抗不得改名 x |
| 原始和固定负荷 | `original_p/original_q`、`fixed_p/fixed_q` | kW / kvar，`(n,)`；固定负荷的可调节点置零 |
| 无功比例 | `q_ratio` | q/p，无量纲；不是功率因数 |
| 电压限值 | `vmin / vmax` | 电压平方 p.u.²，`(n,)` |
| 送端有功上限 | `capacity` | p.u.，`(t,)`；不等于电流或视在功率上限 |
| 源端限值 | `source_pmax/source_qmax/source_smax` | p.u.，含线路损耗 |
| 独立负荷总量界 | `power_limit` | kW |
| 预算 | `budget`，费用 `cost_offset + cost @ x` | FourBus 为元；Case33 为开合次数 |

根不占 `nodes` 和 `state` 的电压元素，父节点哨兵为 -1。型号按走廊输入顺序展开，`type_slices/type_corridor` 表示型号和走廊之间的唯一映射。具名 `loads[i]` 使用真实节点 ID，数组使用 `load_nodes` 顺序。`required` 是必须接入节点掩码；`road_allowed` 限制可用走廊。

阻抗基值为 `voltage_kv**2/(base/1000)`。`Network.loads(power)` 返回标幺节点 p/q，形状 `(batch,n)`；单点也保留批次轴。FourBus 输入层 `r_ohm_km/x_ohm_km/capacity_kw/cable_cny_m` 在入口换算，不能直接当标幺量。

Case33 保留原 MATPOWER 线路与背景负荷，仅七条开关可变；相对初始状态的一次闭合或断开各计一次，预算 `switch_budget=7`。常闭开关成本系数为 -1，常开为 +1，`cost_offset` 补回初始常闭数。

入口共同采用 `CURRENT_LIMIT=200` A（实验指定值，不是原数据额定值）。`Case33(current_limit=...)` 可显式传标量或 37 项数组，构造器默认 inf 保留。`ell_limit=(current_limit/I_base)**2`，`I_base=base/(sqrt(3)*voltage_kv)` A；型号及树上的 `ell_limit` 始终是标幺电流平方。主线与 AC/SOCP 参考使用同一限额。

## 拓扑与运行状态

`receiving/sending/incidence` 形状 `(n,m)`，`incidence=receiving-sending`，流入为正。`E` 为 `(n,d)` 的负荷嵌入矩阵。`active_nodes` 是接入二元变量，`MasterProblem.__init__.f` 是连通虚拟流，与电功率无关。拓扑约束为每走廊至多一型号、边两端接入、`Gf=a`、`-n*z<=f<=n*z`、边数等于接入非根节点数。

外部 `plan/choice/fixed_plan` 为“走廊 ID → 型号 ID 或 None”字典，通过 `encode_plan/decode_plan` 转换；`start` 是型号向量，`incumbent` 是带证书的答案字典。运行树只保存实际接入节点，`node_indices/type_indices` 回指全图，`parent/children/order/roots` 使用树局部索引；`D` 汇总下游节点。`direction=±1` 区分根向和参考方向。反向支路必须计损耗：`P=-P_tree+r*ell`，`Q=-Q_tree+reactance*ell`。

| 状态分块 | 唯一切片 | 单位 |
|---|---|---|
| P | `P_slice = [0:t]` | 参考送端有功 p.u. |
| Q | `Q_slice = [t:2*t]` | 参考送端无功 p.u. |
| ell | `ell_slice = [2*t:3*t]` | 电流平方 p.u.² |
| v | `v_slice = [3*t:3*t+n]` | 非根节点电压平方 p.u.² |
| s+、s- | `slack_slice = [3*t+n:3*t+n+2*m]` | 开断压降余量，先 plus 后 minus |

`operation.P/Q/ell/v/plus/minus` 与 `state` 是同一批变量的具名与扁平视图，根电压恒为 1。断线余量满足 `0<=s±<=drop_max*(1-z)`；未选型号功率和电流为零。全局 `y_lb_global/y_ub_global` 是割计算使用的有效盒，不替代完整运行约束。`pmin/pmax/qmin/qmax/ellmax` 是派生界，有限电流时 ellmax 同时受 ell_limit 限制。带符号模型允许源端反送。

当前约束由 `GridPhysics.add_operation` 装配；数学式中的 A/B/C 不对应另一套代码属性。SOCP 电流约束为 `P²+Q²<=u*ell`；独立 AC 改为等式，其余拓扑、预算、运行约束相同。

## 主线、SP 与联合割

主线按符号分区执行：支撑初始化、首轮射线、SP 违反量选点切割、边界补充、完整物理查漏。`build_region` 在所有分区间分配同一总时限，每区不超过 `partition_seconds`。固定网架和符号内可取认证点凸包，跨网架和分区只取并集。`cache/powers/applied/initialized/counts/small_cuts/area_ratio` 是逐网架流程状态；缓存只复用对应方案与点的证据。

MP2 最大化 `direction @ p`，`objective/bound` 是 kW 的可行值/全局上界；MP1 最小化费用，`objective/bound` 是费用值/全局下界。成功答案包含 `x/p/state`；明确不可行返回 None。`radial_gap_kw=1e-3` 只检查目标间隙，不改写功率。原始求解质量不合格或无可靠答案时停止，不修补或重求。

SP 固定 x,p 后最小化非负 `eta`（求解器名 violation），只允许功率平衡与压降等式 ±eta，其余约束保持严格。`max(eta,0)+MaxVio<=PLANNING_TOL` 才接受原始 `state`。正 eta 通过锥的必要支撑平面 LP 生成有效联合割；`score_only=True` 保存 `cone_normals` 并延迟至 `generate_cut` 取割。超时、数值失败或无法分离均不能作为不可行证书。

唯一联合割方向是 `alpha + beta @ p + delta @ x >= 0`，数组布局 `[alpha, *beta, *delta]`，长度 `1+d+t`；beta 作用于 kW。`_separating_cut.dual` 按求解器线性行排列，`coefficients` 按变量索引，`h` 提取全部 state 列。采用 Gurobi 行乘子方向，固定 x/p 的等式乘子先置零：`alpha=-dual@rhs + max(h,0)@y_ub_global + min(h,0)@y_lb_global`。正比例归一化后截距加 `1e-10` 保守补偿。`calls/cut_calls` 分别计 SP 与取割 LP 次数。

方案紧化：方案首次登记时 `GridPhysics.obbt` 在独立功率 `[0,bounds]` 内对选中型号的 P/Q/ell 与节点 v 做 `OBBT_ROUNDS` 轮 OBBT，盒存于 `GridPhysics.boxes[tuple(x)]`，键 `(量名, 型号键或节点)`、标幺；端点外扩 `OBBT_PAD`，未证得最优的端点保留原界（首轮为全局界）。`add_operation(..., scheme=x)` 对已紧化方案追加盒约束与反向锥包络割 `u_L*ell+ell_L*u-u_L*ell_L <= (P_L+P_U)P-P_L*P_U+(Q_L+Q_U)Q-Q_L*Q_U+ENVELOPE_MARGIN`（另一条取 u_U、ell_U，u 为参考送端电压平方），每行右端加 `M*H(x)`：H 为 x 与该方案的汉明距离，M 为该行在全局界上的最大违反量，故联合割对全部 x 仍有效。SP、射线、OBBT 与固定方案查漏传入 scheme；含紧化行的模型和 OBBT 用 `OBBT_CONV_TOL`。`obbt_extremes` 以 `OBBT_WORKERS` 份模型副本并行求解，结果与串行逐位相同。

有限逐线路电流的 Port SP 固定使用既有 Gurobi 连续 SOCP；其他 Port SP 和射线使用 Clarabel。`solve_conic.basis/offset` 表示消元坐标，`cone_scale` 只是等价锥缩放。射线内域见证使用既有 `1e-6` 锥裕量，有限电流界同时扣除相同裕量；原物理残差仍按 `1e-8` 验收。`numeric_focus` 是求解设置，不改变数学可行性。

## 正负功率、几何与覆盖

`main.mode=0` 使用原非负负荷与 q_ratio；`mode=1` 使用固定符号区 `PortPhysics.sign`，内部非负幅值 u，真实 p=sign*u；负荷 PF=0.95、光伏 PF=1，允许反送。输出、回放和扫描始终是真实 kW；割只在所属符号区有效。

`bounds` 为公共正数坐标尺度 b，内部点 `xi=u/bounds` 无量纲；`axis_bounds` 是该分区的方向全局上界，`total_bound` 是总量上界，单位 kW。结果 `bound` 是查询目标的标量界，不能与 bounds 混用。

`RegionState.records[tuple(x)]` 中 inner/outer 和 `add_point/covering_schemes/witness_support` 用归一化坐标；`finish` 导出的 `inner/outer[*].vertices` 已变回 kW。`halfspaces/contains/clip_polytope` 不自动换算。面方程为 `F@xi+g<=0`；裁剪接口为 `constant+coefficient@xi>=0`。固定 x 的联合割转换为 `alpha+delta@x+(beta*bounds)@xi>=0`。

`tau=0.005` 是径向精度，s=1-tau；允许外扩对应归一化点缩至 s*xi。完整 `RemainingRegionModel` 搜索扩张内域并集以外的物理点，覆盖目标为 `max_(x,p) min_known_scheme max_face (F*s*xi+g)`。局部变量 delta 是未覆盖距离，内部乘 `distance_scale`、返回时还原；与联合割的 delta 向量不同。只有可靠全局上界<=COVERAGE_TOL 或已证不可行才 complete。

查漏由 `remaining_search` 按方案分解：已登记方案以 `RemainingRegionModel(..., scheme=x)` 固定 x（含紧化行）；其余方案以 `exclude` 排除已登记及 `covered` 后联合求解（纯 SOCP），所得新方案先 OBBT 再固定 x 复核，无见证则记入 covered。查漏模型给出边界上的候选点与覆盖上界；`register` 从内点锚（原点已认证时取原点，否则取内域重心，新方案无内点时取 `GridPhysics.center`）沿射线取朝候选点的最远紧化可行点登记，并写回答案 p。

该模型始终包含二元建设变量、完整运行状态和 SOCP 约束，无 light/cuts_only 入口。已知方案全部局部认证不能替代全局覆盖证明；超时保留认证内域和有效外包络。几何外扩用 center/radius/envelope 表示中心、内切半径及包络，不更改物理状态。

## 支持面认证实验

二维入口固定 FourBus(1,2)、mode=0；同一支持核心供三维 (1,2,3)、mode=1 使用。由总负荷 MP2 给出初始方案，随后 `discover_schemes` 用完整 physical 全局查询按需发现新方案；不先枚举建设方案，不用 no-good 排除。

当前内域的每个面 `a@xi+b<=0` 需要满足 `(1-tau)*h_x(a)+b<=epsilon_geom`，其中 `h_x(a)=max a@xi` 是固定方案完整 SOCP 支持值。先用当前外域顶点的几何上界检查全部面，再查询必要支持方向。可靠 UB 可认证该面，通过原约束审计的越界点用于补点；其余为 UNRESOLVED。补点后立即重建面，废弃旧面结论。

`SupportOracle.support_cache` 键为完整方案和有朝向的归一化法向，缓存 lb/ub/point/state/status/solve_seconds/normal，不缓存截距或某一面结论。近方向复用需加 `[0,1]^d` 上 `max (a-a_cached)@xi` 的保守补偿，再还原法向长度；正反法向不合并。`SupportOracle.solve` 目标直接为 `sum(a_i*p_i/bounds_i)`，ObjVal/ObjBound 均无量纲。`BOUND_PAD=1e-10` 只放宽可靠支持上界，不放宽物理门槛。

`audit_incumbent` 直接检查原始约束、变量界、整数性及 MaxVio，不修补。`conditional_support_cut` 用 `Delta=sum(x0)+(1-2*x0)@x`，`M=max(0,max_box(a@xi)-ub)`，构造 `ub+M*Delta-a@xi>=0`，仅在 x=x0 时收紧；仍返回原联合割布局。

inner/outer/point 支持 d=2/3；faces 为 `(F,d+1)`，offset 始终最后一列，三维日志增加 normal_3。area 原名保留表示 d 维测度，单位 kW^d；initial_support_calls 为 d+1。`faces_seen` 是累计面访问数，`unique_faces` 是去重数，最终一轮比例与累计几何率分开。`max_iterations/max_discovery_iterations/deadline` 控制局部/全局迭代及时限，触顶不报认证成功。

schema 为 `support-face-fourbus-physical-v3`，coverage.mode=physical。`global_physical_calls` 统计完整查询；原 `enumeration_mip_calls/enumeration_seconds/global_milp_calls/discovery_sp_calls/discovery_cut_lp_calls/discovery_sp_seconds` 名称保留，现均为零。`conditional_support_cuts/support_calls` 保持原义。`SupportReplay` 使用同一个 RunMonitor 格式，`bounds/monitor` 可由分区入口传入；2D 的 `record` 只决定录制，不影响判据。

## 径向夹逼的体积准则与方案二对照实验

径向夹逼实验（experiments/test_radial_sandwich.py、test_radial_sandwich_3d.py）的 θ、rho、va/vb、mu、g、Q、c、α 仍是脚本内局部量。本节登记 2026-10-01 加入的接口；默认行为不变（二维 `criterion='radial'`）。`RadialSandwich.criterion='volume'` 时与三维相同：体积缺口 Δ_k=(μ̄_k^d-1)·vol(T_k)（`interval_delta`，xi^d），每次细分 Δ_k 最大者，ΣΔ_k<=`epsilon`·Σvol(T_k)（`volume_ratio`）即认证，`epsilon` 默认 d·tau，区间 MISOCP 解到 `mip_gap`。`run_volume` 可续跑：已有叶锥时不重建根锥，先求 pending 锥，再从全部 bounded 锥重建堆。`on_solution(problem, record)` 是每轮 MISOCP 求解后、释放模型前的空钩子；`Interval.final` 为区间最后一轮 MISOCP 记录。

方案二对照实验（experiments/compare_plan2.py，纯几何在 experiments/plan2_geometry.py）比较方法 R（径向锥 + 体积准则 + 近端覆盖，即 B 型两网架证书 T ⊆ R_x̂ ∪ R_y）与方法 H（R 的 A 阶段 + 逐网架支撑查询 B + 覆盖证书 C）。坐标为 xi=u/bounds，分区盒 [0,1]^d；面方程沿用 `[n, b]`，|n|=1，beta=-b。

| 数学量 | 代码 | 单位与含义 |
|---|---|---|
| X*、V_x | `SupportPhase.networks`（方案元组 → `NetworkState`）、`NetworkState.points` | A 阶段出现过的全部网架；该网架已认证点（射线远/近端、零接入点、审计过的已紧化现任解、支撑与覆盖见证），xi |
| P_x、O_x | `NetworkState.vertices`、`NetworkState.outer` | conv(V_x) 的顶点；分区盒逐次被可靠支撑上界裁剪的外界，xi |
| c_x | `NetworkState.center`，`interior_point` | P_x 的 Chebyshev 中心（LP），失败时取顶点形心 |
| 面判据 | `NetworkState.criterion`，`choose_criterion` | origin：(1-tau)·UB+b<=epsilon_geom；center：UB-beta<=tau·(beta-n@c_x)+epsilon_geom；auto：原点在 P_x 内且到每个非分区边界面的距离>=`ORIGIN_CLEARANCE` 时取 origin |
| E_x | `expanded_faces`、`SupportPhase.cover` | origin 为 P_x/(1-tau)，center 为 c_x+(1+tau)(P_x-c_x)；覆盖证书完成时冻结所用各 E_x 的面方程 |
| ε、ε_A、share_A、ε_B | `settings:eps / discovery_eps / discovery_share / network_eps` | 体积目标（默认 d·tau）、A 阶段放宽目标与时限占比、网架停止 vol(O_x)/vol(P_x)-1 |
| I_H、O_H | `h_measures` | I_R ∪ (∪P_x)；覆盖前 O_R∩盒，覆盖后 O_R∩盒∩(∪E_x)；返回 xi^d，乘 prod(bounds) 得 kW^d |

`PartitionOracle` 是 `SupportOracle` 加分区盒 0<=u<=bounds 与种子，只接受 OPTIMAL；非 OPTIMAL 不重求，该面记 UNRESOLVED（代码库没有重试流程）。`boundary_faces` 识别分区边界面（xi_j>=0 或 xi_j<=1，容差 `BOUNDARY_TOL`），它们由分区盒本身认证。覆盖证书用 `RemainingRegionModel` 直接接收 E_x 面并令 tau=0（即 s=1），再加全部已紧化方案的提升行，见证方案未紧化时懒惰 OBBT 后重解；上界<=GEOMETRY_TOL 或已证不可行即覆盖完成。H 的分区认证为：覆盖完成、无 UNRESOLVED 面且 vol(O_H)-vol(I_H)<=ε·vol(I_H)；或径向部分自身满足体积准则（继承 R 的证书，`SupportPhase.certified_now`）。

计时与记录：`--seconds` 为一次运行的墙钟总时限（含 OBBT、求解、几何与记录，不含事后逐格评价），t=0 为运行开始。串行时第 j 个分区得到剩余时间的 1/(余下分区数)，并行时各分区同时开始。`Recorder.log/record` 每次求解写一行（`Recorder._row:type` 为 MISOCP/SOCP/LP，`purpose` 为 discovery/cone_outer/support/coverage/obbt/ray/zero_sp）；`Recorder.changed` 在内外测度变化时写时间线（`inner_measure/outer_measure` 为 kW^d，`gap`=outer/inner-1）；`CHECKPOINTS` 各时刻保存变化前的确切几何快照。运行汇总 `write_run:t_cert` 为全部分区认证时刻（未认证为 None），`write_run:t_gap` 为总间隙首次降到 10%/5%/2% 的时刻。

逐格评价 `grid_metrics` 以 AC 为主参考、SOCP 为诊断（未决标签不计入该参考）。`validity` 的必要条件：内域不含 SOCP 已证不可行格（`validity:inner_socp_infeasible_cells`=0，紧化模型 ⊆ SOCP）；AC 可行格都在外界内（`validity:outer_missed_ac_cells`=0）；AC 未决格单列（`validity:undecided_ac_cells`）。`validity:outer_missed_socp_cells` 是 OBBT 去掉的 SOCP-only 区域，只作诊断，不参与 `validity:valid`。

## 独立扫描、结果与回放

`vertify.py` 是独立 AC/SOCP 扫描入口。`ac_network` 复制完整配置，不删除限流。`budget_schemes` 仅供小算例独立 AC 参考的拓扑审计，构域不调用它。AC 单树只提供可行见证；全拓扑必要条件排除或完整 AC 不可行证书才给负标签。迭代失败、超时和未知均不可当作不可行。

`ACPowerFlow` 接收运行树，state(power,ell) 返回 `(P,Q,v,u)`，均为 `(batch,n_tree)`，v/u 为受端/送端电压平方。内部 `_state.p/q` 是节点标幺负荷；AC 残差为 `P²+Q²-u*ell`。`global_status` 供独立交叉检查。

`AC_CACHE_METHOD=ac_socp_grid_v4`，`ac_identity` 包含物理参数、ell_limit、模式、功率因数、有序节点和预算，不含构域 tau、时间和线程。`reference_box` 用同配置完整 SOCP 的方向全局上界覆盖两种参考。`scan_problem/scan_line` 逐点独立求解；SOCP 固定点使用 eta 目标，接受与排除依赖原门槛；AC 不读取 SOCP 标签。

| 缓存字段 | 固定含义 |
|---|---|
| `axis_lower / bounds` | 真实 kW 扫描边界 `(d,)` |
| `origin / step / start` | kW 原点/步长、整数起始索引；点为 `origin+(start+index+0.5)*step` |
| `states / witness_x / residual` | AC 标签、型号见证、残差；形状分别 `(n1,...,nd)`、`(n1,...,nd,t)`、`(n1,...,nd)` |
| `socp_states / socp_witness_x / socp_residual` | 同坐标 SOCP 的独立结果，形状同上 |
| `1 / -1 / 0` | 可行 / 已证不可行 / 未决；未决时不发布完整误差率 |
| `cache_path` | 实际 NPZ 路径；归档记录可使用相对项目根目录的路径 |

`scan_path` 给出 `results/scan/<网络>/<节点>/<身份>/`，`region_path` 由 origin/step/start/shape 生成区域摘要。扩界保留格点，部分覆盖只补算缺点；AC/SOCP 分别补零标签。`import_ac_reference` 仅接受同版本同身份数据。`force_rescan` 先备份；每批由主进程原子落盘，独占锁避免并发覆盖。

比较使用最终 inner 作为计算域，另保留 outer 指标。`comparison_metrics` 的 MR=`missed/reference`、FR=`extra/computed`，百分数，空分母 None；同时保留 missed_cells/extra_cells/reference_cells/computed_cells。`comparisons` 键为 result_ac/result_socp/socp_ac。导出版本 paired_scan_v2，NPZ 的 power 为 `(N,d)` kW，ac_states/socp_states/inner/outer 对应同坐标；CSV 标签列为 ac_state/socp_state。网格误差不是连续体积证明。

回放 version=4：history 保存增量过程，validation_state 保存最终扫描；version=3 兼容读取。`RunMonitor._geometry:x` 保留完整型号顺序，`cut_history` 保存真实割。符号分区由 partition/sign 标记；各图按初始化包络定范围，显示层不改模型坐标。仅当前帧之前的网架、割、点可出现在过程图；最终扫描独立显示。result 包含 status/certified/inner/outer/counts/timing，分区证书在 partition_end；数值未决和时限未完不能改成 certified。JSON 非有限值为 null，仅预算 null 可按 inf 解释。

## 容差与默认设置

| 配置 | 当前值 |
|---|---|
| `PLANNING_TOL / GEOMETRY_TOL` | 各 1e-8；物理残差/归一化几何门槛，语义独立 |
| `COVERAGE_TOL` | 1e-4；全局查漏覆盖门槛与认证外包络的外扩量（归一化距离） |
| `OBBT_ROUNDS / OBBT_PAD / ENVELOPE_MARGIN / OBBT_CONV_TOL` | 2 / 1e-6 标幺 / 1e-4 标幺² / 1e-8；查漏模型的 barrier 收敛容差 1e-6 |
| `PARTITION_TIME_LIMIT / OBBT_WORKERS` | 20 秒 / 8 |
| `AC_TOL / FIXED_POINT_TOL / GLOBAL_AC_TOL` | 1e-9 / 1e-12 / 1e-7 |
| `REGION_TAU` | 0.005；不是体积误差或物理容差 |
| `NETWORK / DIMENSION / mode` | Case33 / 3 / 1 |
| `CASE_TIME_LIMIT / SOLVER_THREADS` | 300 秒 / 20；总时限不含事后扫描 |
| `CUT_THRESHOLD / CUT_PATIENCE` | 0.02 / 3 |
| `RAY_THRESHOLD / POINT_TOL` | 1e-4 几何收益比例 / 0.01 kW |
| `DIVISIONS / SCAN_DIVISIONS / SCAN_WORKERS` | 160 / {2:160,3:80} / 20；每扫描进程一个求解线程 |

单次 MP/SP/剩余查询/AC 的时限仍用 MP_TIME_LIMIT/SP_TIME_LIMIT/RESIDUAL_TIME_LIMIT/AC_TIME_LIMIT；AC_ITERATIONS 为迭代次数。不得合并不同语义的容差。数值配置变化须同步本表和结果协议。

## 本次整理的显式迁移（2026-09-30）

结果统一在 results/mainline、results/support_face、results/scan；映射、原路径、文件哈希和来源见 results/manifest.json。Case33 记录的 cache_path/ac_cache 改为相对项目根目录的 results/scan 路径，main.run 同步解析；求解数据和扫描数组不变。FourBus 已保存参考仍为历史 SOCP 格点，不能冒充当前配对 AC/SOCP 缓存。

退役内容为 Concept5、load_jiangkou/Jiangkou、legacy_case33 升级夹具、fixed_topology/upgrade_plan/affordable_designs/dispatch_state 测试辅助、旧前端测试和一次性数值诊断脚本。其登记随接口退役归档于 Git 提交 e79f18b；现存量未改名或删除登记。旧迁移叙述和源码副本可从同一提交恢复，当前文档不重复叠加互相覆盖的历史规则。

交付运行 `python -m unittest tests.test_notation -v`；算法修改另跑相关数值回归。

## 当前登记接口索引

以下保留整理前全部现存接口登记，以便名称检查覆盖字段、形参和结果键。`作用域.名称` 表示字段或局部量，`函数:键` 表示返回字典键；语义按上文分组约定。

### Network/case33bw.py

- [CURRENT_LIMIT](../Network/case33bw.py)、[Case33](../Network/case33bw.py)、[Case33.__init__.current_limit](../Network/case33bw.py)、[Case33.switch_budget](../Network/case33bw.py)
- [LOAD_NODES](../Network/case33bw.py)

### vertify.py

- [ACPowerFlow](../vertify.py)、[ACPowerFlow._state.p](../vertify.py)、[ACPowerFlow._state.q](../vertify.py)、[ACPowerFlow.classify](../vertify.py)
- [AC_CACHE_METHOD](../vertify.py)、[AC_ITERATIONS](../vertify.py)、[AC_TIME_LIMIT](../vertify.py)、[AC_TOL](../vertify.py)
- [FIXED_POINT_TOL](../vertify.py)、[GLOBAL_AC_TOL](../vertify.py)、[SCAN_FIELDS](../vertify.py)、[ac_identity](../vertify.py)
- [ac_interval_possible](../vertify.py)、[ac_network](../vertify.py)、[ac_scan_line](../vertify.py)、[ac_scan_line:residual](../vertify.py)
- [ac_scan_line:witness_x](../vertify.py)、[budget_schemes](../vertify.py)、[export_comparison](../vertify.py)、[import_ac_reference](../vertify.py)
- [reference_box](../vertify.py)、[region_path](../vertify.py)、[save_scan](../vertify.py)、[scan_ac_reference](../vertify.py)
- [scan_ac_reference.origin](../vertify.py)、[scan_ac_reference.start](../vertify.py)、[scan_ac_reference.states](../vertify.py)、[scan_ac_reference.step](../vertify.py)
- [scan_line](../vertify.py)、[scan_path](../vertify.py)、[scan_problem](../vertify.py)、[scan_problem.eta](../vertify.py)
- [signed_ac_witness](../vertify.py)

### monitor.py

- [COMPARISONS](../monitor.py)、[RunMonitor](../monitor.py)、[RunMonitor._geometry](../monitor.py)、[RunMonitor._geometry:x](../monitor.py)
- [RunMonitor.global_end](../monitor.py)、[RunMonitor.partition_frame](../monitor.py)、[RunMonitor.seed](../monitor.py)、[RunMonitor.sp_start](../monitor.py)
- [RunMonitor.updated:cut_history](../monitor.py)、[RunMonitor.validation](../monitor.py)、[RunMonitor.validation.region_key](../monitor.py)、[RunMonitor.validation:comparisons](../monitor.py)
- [RunMonitor.validation:metrics](../monitor.py)、[RunMonitor.validation_state](../monitor.py)、[_cut_polygon](../monitor.py)、[_cut_polygon.axis_lower](../monitor.py)
- [_cut_segment](../monitor.py)、[_cut_segment.axis_lower](../monitor.py)、[_draw_3d](../monitor.py)、[_voxel_faces](../monitor.py)
- [comparison_metrics](../monitor.py)、[comparison_metrics:computed_cells](../monitor.py)、[comparison_metrics:extra_cells](../monitor.py)、[comparison_metrics:fr_percent](../monitor.py)
- [comparison_metrics:missed_cells](../monitor.py)、[comparison_metrics:mr_percent](../monitor.py)、[comparison_metrics:reference_cells](../monitor.py)

### model.py

- [DEFAULT_SOLVER_THREADS](../model.py)、[GridPhysics](../model.py)、[GridPhysics.P_slice](../model.py)、[GridPhysics.Q_slice](../model.py)
- [GridPhysics.add_operation](../model.py)、[GridPhysics.add_operation.P](../model.py)、[GridPhysics.add_operation.Q](../model.py)、[GridPhysics.add_operation.ell](../model.py)
- [GridPhysics.add_operation.minus](../model.py)、[GridPhysics.add_operation.plus](../model.py)、[GridPhysics.add_operation.v](../model.py)、[GridPhysics.cost](../model.py)
- [GridPhysics.drop_max](../model.py)、[GridPhysics.ell_slice](../model.py)、[GridPhysics.ellmax](../model.py)、[GridPhysics.ends](../model.py)
- [GridPhysics.incoming](../model.py)、[GridPhysics.keys](../model.py)、[GridPhysics.outgoing](../model.py)、[GridPhysics.pmax](../model.py)
- [GridPhysics.pmin](../model.py)、[GridPhysics.qmax](../model.py)、[GridPhysics.qmin](../model.py)、[GridPhysics.r](../model.py)
- [GridPhysics.reactance](../model.py)、[GridPhysics.slack_slice](../model.py)、[GridPhysics.types](../model.py)、[GridPhysics.v_slice](../model.py)
- [GridPhysics.vmax](../model.py)、[GridPhysics.vmin](../model.py)、[GridPhysics.y_lb_global](../model.py)、[GridPhysics.y_ub_global](../model.py)
- [MP_TIME_LIMIT](../model.py)、[MasterProblem.__init__.budget](../model.py)、[MasterProblem.__init__.cuts](../model.py)、[MasterProblem.__init__.direction](../model.py)
- [MasterProblem.__init__.f](../model.py)、[MasterProblem.__init__.min_total](../model.py)、[MasterProblem.__init__.power](../model.py)、[MasterProblem.active_nodes](../model.py)
- [MasterProblem.choices](../model.py)、[MasterProblem.direction](../model.py)、[MasterProblem.loads](../model.py)、[MasterProblem.operation](../model.py)
- [MasterProblem.power](../model.py)、[MasterProblem.solve.incumbent](../model.py)、[MasterProblem.solve.radial_gap_kw](../model.py)、[MasterProblem.solve.start](../model.py)
- [MasterProblem.solve:bound](../model.py)、[MasterProblem.solve:feasible](../model.py)、[MasterProblem.solve:objective](../model.py)、[MasterProblem.solve:p](../model.py)
- [MasterProblem.solve:state](../model.py)、[MasterProblem.solve:status](../model.py)、[MasterProblem.solve:x](../model.py)、[MasterProblem.state](../model.py)
- [MasterProblem.x](../model.py)、[PLANNING_TOL](../model.py)、[PortPhysics.sign](../model.py)、[RESIDUAL_TIME_LIMIT](../model.py)
- [RemainingRegionModel.__init__.axis_bounds](../model.py)、[RemainingRegionModel.__init__.delta](../model.py)、[RemainingRegionModel.distance_scale](../model.py)、[RemainingRegionModel.solve:bound](../model.py)
- [RemainingRegionModel.solve:complete](../model.py)、[SP_TIME_LIMIT](../model.py)、[SubProblem._build](../model.py)、[SubProblem._build.eta](../model.py)
- [SubProblem._cut](../model.py)、[SubProblem._separating_cut.coefficients](../model.py)、[SubProblem._separating_cut.cut](../model.py)、[SubProblem._separating_cut.dual](../model.py)
- [SubProblem._separating_cut.h](../model.py)、[SubProblem._separating_cut.rows](../model.py)、[SubProblem.cut_calls](../model.py)、[SubProblem.generate_cut](../model.py)
- [SubProblem.numeric_focus](../model.py)、[SubProblem.solve](../model.py)、[SubProblem.solve.cone_normals](../model.py)、[SubProblem.solve.eta](../model.py)
- [SubProblem.solve.score_only](../model.py)、[SubProblem.solve:cone_normals](../model.py)、[SubProblem.solve:cut](../model.py)、[SubProblem.solve:eta](../model.py)
- [SubProblem.solve:feasible](../model.py)、[SubProblem.solve:state](../model.py)、[port_bounds](../model.py)、[ray_support](../model.py)
- [ray_support.numeric_focus](../model.py)、[ray_support.ray_fraction](../model.py)、[solve_conic.basis](../model.py)、[solve_conic.cone_scale](../model.py)
- [solve_conic.offset](../model.py)、[voltage_flow_bounds](../model.py)、[OBBT_ROUNDS](../model.py)、[OBBT_PAD](../model.py)
- [ENVELOPE_MARGIN](../model.py)、[OBBT_CONV_TOL](../model.py)、[OBBT_PARAMS](../model.py)、[obbt_pool](../model.py)
- [obbt_extremes](../model.py)、[GridPhysics.boxes](../model.py)、[GridPhysics.obbt](../model.py)、[GridPhysics.center](../model.py)
- [GridPhysics.add_operation.scheme](../model.py)、[RemainingRegionModel.__init__.scheme](../model.py)、[RemainingRegionModel.__init__.exclude](../model.py)

### Network/__init__.py

- [Corridor.endpoints](../Network/__init__.py)、[Corridor.existing_type](../Network/__init__.py)、[Corridor.id](../Network/__init__.py)、[Corridor.initial_active](../Network/__init__.py)
- [Corridor.switchable](../Network/__init__.py)、[Corridor.types](../Network/__init__.py)、[Network.E](../Network/__init__.py)、[Network.base](../Network/__init__.py)
- [Network.budgets](../Network/__init__.py)、[Network.capacity](../Network/__init__.py)、[Network.corridor_types](../Network/__init__.py)、[Network.corridors](../Network/__init__.py)
- [Network.cost](../Network/__init__.py)、[Network.cost_offset](../Network/__init__.py)、[Network.cost_unit](../Network/__init__.py)、[Network.decode_plan](../Network/__init__.py)
- [Network.ell_limit](../Network/__init__.py)、[Network.encode_plan](../Network/__init__.py)、[Network.fixed_p](../Network/__init__.py)、[Network.fixed_q](../Network/__init__.py)
- [Network.incidence](../Network/__init__.py)、[Network.load_nodes](../Network/__init__.py)、[Network.loads](../Network/__init__.py)、[Network.n](../Network/__init__.py)
- [Network.n_corridors](../Network/__init__.py)、[Network.n_types](../Network/__init__.py)、[Network.node_index](../Network/__init__.py)、[Network.nodes](../Network/__init__.py)
- [Network.original_p](../Network/__init__.py)、[Network.original_q](../Network/__init__.py)、[Network.power_limit](../Network/__init__.py)、[Network.q_ratio](../Network/__init__.py)
- [Network.r](../Network/__init__.py)、[Network.reactance](../Network/__init__.py)、[Network.receivers](../Network/__init__.py)、[Network.receiving](../Network/__init__.py)
- [Network.required](../Network/__init__.py)、[Network.road_allowed](../Network/__init__.py)、[Network.root](../Network/__init__.py)、[Network.selected](../Network/__init__.py)
- [Network.senders](../Network/__init__.py)、[Network.sending](../Network/__init__.py)、[Network.source_pmax](../Network/__init__.py)、[Network.source_qmax](../Network/__init__.py)
- [Network.source_smax](../Network/__init__.py)、[Network.tree](../Network/__init__.py)、[Network.type_corridor](../Network/__init__.py)、[Network.type_keys](../Network/__init__.py)
- [Network.type_slices](../Network/__init__.py)、[Network.vmax](../Network/__init__.py)、[Network.vmin](../Network/__init__.py)、[Network.voltage_kv](../Network/__init__.py)
- [OperatingTree.D](../Network/__init__.py)、[OperatingTree.children](../Network/__init__.py)、[OperatingTree.direction](../Network/__init__.py)、[OperatingTree.ell_limit](../Network/__init__.py)
- [OperatingTree.n](../Network/__init__.py)、[OperatingTree.node_indices](../Network/__init__.py)、[OperatingTree.order](../Network/__init__.py)、[OperatingTree.parent](../Network/__init__.py)
- [OperatingTree.roots](../Network/__init__.py)、[OperatingTree.type_indices](../Network/__init__.py)、[TypeParameters.capacity](../Network/__init__.py)、[TypeParameters.ell_limit](../Network/__init__.py)
- [TypeParameters.id](../Network/__init__.py)、[TypeParameters.investment_cost](../Network/__init__.py)、[TypeParameters.r](../Network/__init__.py)、[TypeParameters.reactance](../Network/__init__.py)

### experiments/test_support_face_certification_fourbus_2d.py

- [BOUND_PAD](../experiments/test_support_face_certification_fourbus_2d.py)、[SupportOracle.solve](../experiments/test_support_face_certification_fourbus_2d.py)、[SupportOracle.support_cache](../experiments/test_support_face_certification_fourbus_2d.py)、[SupportReplay](../experiments/test_support_face_certification_fourbus_2d.py)
- [SupportReplay.__init__.bounds](../experiments/test_support_face_certification_fourbus_2d.py)、[SupportReplay.__init__.monitor](../experiments/test_support_face_certification_fourbus_2d.py)、[audit_incumbent](../experiments/test_support_face_certification_fourbus_2d.py)、[cached_support](../experiments/test_support_face_certification_fourbus_2d.py)
- [certify_scheme](../experiments/test_support_face_certification_fourbus_2d.py)、[certify_scheme.deadline](../experiments/test_support_face_certification_fourbus_2d.py)、[classify_support](../experiments/test_support_face_certification_fourbus_2d.py)、[conditional_support_cut](../experiments/test_support_face_certification_fourbus_2d.py)
- [discover_schemes](../experiments/test_support_face_certification_fourbus_2d.py)、[discover_schemes.deadline](../experiments/test_support_face_certification_fourbus_2d.py)、[discover_schemes:bound](../experiments/test_support_face_certification_fourbus_2d.py)、[discover_schemes:global_physical_calls](../experiments/test_support_face_certification_fourbus_2d.py)
- [geometry_margins](../experiments/test_support_face_certification_fourbus_2d.py)、[initial_bounds](../experiments/test_support_face_certification_fourbus_2d.py)、[initial_bounds.bounds](../experiments/test_support_face_certification_fourbus_2d.py)、[initial_bounds.deadline](../experiments/test_support_face_certification_fourbus_2d.py)
- [run_experiment.max_discovery_iterations](../experiments/test_support_face_certification_fourbus_2d.py)、[run_experiment.record](../experiments/test_support_face_certification_fourbus_2d.py)

### main.py

- [BUDGET](../main.py)、[CASE_TIME_LIMIT](../main.py)、[CUT_PATIENCE](../main.py)、[CUT_THRESHOLD](../main.py)
- [DIMENSION](../main.py)、[DIVISIONS](../main.py)、[FORCE_RESCAN](../main.py)、[NETWORK](../main.py)
- [POINT_TOL](../main.py)、[RAY_THRESHOLD](../main.py)、[REGION_TAU](../main.py)、[SCAN_DIVISIONS](../main.py)
- [SCAN_OUTPUT](../main.py)、[SCAN_WORKERS](../main.py)、[SOLVER_THREADS](../main.py)、[main](../main.py)
- [main.case](../main.py)、[main.dimension](../main.py)、[main.load_nodes](../main.py)、[mode](../main.py)
- [recording_path](../main.py)、[run](../main.py)、[run.budget](../main.py)、[run.force_rescan](../main.py)
- [run.scan_output](../main.py)、[OBBT_WORKERS](../main.py)、[PARTITION_TIME_LIMIT](../main.py)、[run.obbt_workers](../main.py)
- [run.partition_seconds](../main.py)

### region.py

- [GEOMETRY_TOL](../region.py)、[RegionState.add_scheme:inner](../region.py)、[RegionState.add_scheme:outer](../region.py)、[RegionState.axis_bounds](../region.py)
- [RegionState.bounds](../region.py)、[RegionState.finish:inner](../region.py)、[RegionState.finish:outer](../region.py)、[RegionState.finish:vertices](../region.py)
- [RegionState.inner_equations](../region.py)、[RegionState.tau](../region.py)、[RegionState.total_bound](../region.py)、[build_region](../region.py)
- [build_sequential_region](../region.py)、[build_sequential_region.add_ray.initial_sweep](../region.py)、[build_sequential_region.applied](../region.py)、[build_sequential_region.area_ratio](../region.py)
- [build_sequential_region.cache](../region.py)、[build_sequential_region.counts](../region.py)、[build_sequential_region.initialized](../region.py)、[build_sequential_region.mode](../region.py)
- [build_sequential_region.point_tol](../region.py)、[build_sequential_region.powers](../region.py)、[build_sequential_region.ray_threshold](../region.py)、[build_sequential_region.sign](../region.py)
- [build_sequential_region.small_cuts](../region.py)、[build_sequential_region.x](../region.py)、[clip_polytope.coefficient](../region.py)、[clip_polytope.constant](../region.py)
- [coverage_halfspaces](../region.py)、[halfspaces](../region.py)、[ray_gain](../region.py)、[register_power](../region.py)
- [stage_candidates](../region.py)、[union_measure](../region.py)、[COVERAGE_TOL](../region.py)、[remaining_search](../region.py)
- [build_sequential_region.covered](../region.py)、[build_sequential_region.obbt_workers](../region.py)、[build_region.partition_seconds](../region.py)、[build_region.obbt_workers](../region.py)

### plot.py

- [_clip_face](../plot.py)

### Network/four_bus_five_corridor.py

- [FourBus.__init__.load_nodes](../Network/four_bus_five_corridor.py)、[LineType.cable_cny_m](../Network/four_bus_five_corridor.py)、[LineType.capacity_kw](../Network/four_bus_five_corridor.py)、[LineType.r_ohm_km](../Network/four_bus_five_corridor.py)
- [LineType.x_ohm_km](../Network/four_bus_five_corridor.py)

### tests/planning_checks.py

- [margin](../tests/planning_checks.py)

### experiments/test_radial_sandwich.py

- [RadialSandwich.criterion](../experiments/test_radial_sandwich.py)、[RadialSandwich.epsilon](../experiments/test_radial_sandwich.py)、[RadialSandwich.run_volume](../experiments/test_radial_sandwich.py)、[RadialSandwich.on_solution](../experiments/test_radial_sandwich.py)
- [RadialSandwich.make_root](../experiments/test_radial_sandwich.py)、[RadialSandwich.solve_interval](../experiments/test_radial_sandwich.py)、[RadialSandwich.interval_delta](../experiments/test_radial_sandwich.py)、[RadialSandwich.volume_ratio](../experiments/test_radial_sandwich.py)
- [Interval.final](../experiments/test_radial_sandwich.py)

### experiments/plan2_geometry.py

- [ORIGIN_CLEARANCE](../experiments/plan2_geometry.py)、[BOUNDARY_TOL](../experiments/plan2_geometry.py)、[boundary_faces](../experiments/plan2_geometry.py)、[chebyshev_center](../experiments/plan2_geometry.py)
- [interior_point](../experiments/plan2_geometry.py)、[choose_criterion](../experiments/plan2_geometry.py)、[allowed_offsets](../experiments/plan2_geometry.py)、[expanded_faces](../experiments/plan2_geometry.py)
- [face_margins](../experiments/plan2_geometry.py)、[classify_center](../experiments/plan2_geometry.py)、[h_measures](../experiments/plan2_geometry.py)、[grid_metrics](../experiments/plan2_geometry.py)
- [validity](../experiments/plan2_geometry.py)、[validity:valid](../experiments/plan2_geometry.py)、[validity:inner_socp_infeasible_cells](../experiments/plan2_geometry.py)、[validity:outer_missed_ac_cells](../experiments/plan2_geometry.py)
- [validity:outer_missed_socp_cells](../experiments/plan2_geometry.py)、[validity:undecided_ac_cells](../experiments/plan2_geometry.py)

### experiments/compare_plan2.py

- [CHECKPOINTS](../experiments/compare_plan2.py)、[SUPPORT_SECONDS](../experiments/compare_plan2.py)、[SUPPORT_PASS_SHARE](../experiments/compare_plan2.py)、[PartitionOracle](../experiments/compare_plan2.py)
- [NetworkState](../experiments/compare_plan2.py)、[NetworkState.points](../experiments/compare_plan2.py)、[NetworkState.vertices](../experiments/compare_plan2.py)、[NetworkState.outer](../experiments/compare_plan2.py)
- [NetworkState.center](../experiments/compare_plan2.py)、[NetworkState.criterion](../experiments/compare_plan2.py)、[SupportPhase](../experiments/compare_plan2.py)、[SupportPhase.networks](../experiments/compare_plan2.py)
- [SupportPhase.cover](../experiments/compare_plan2.py)、[SupportPhase.certified_now](../experiments/compare_plan2.py)、[SupportPhase.coverage](../experiments/compare_plan2.py)、[Recorder](../experiments/compare_plan2.py)
- [Recorder.log](../experiments/compare_plan2.py)、[Recorder.record](../experiments/compare_plan2.py)、[Recorder.changed](../experiments/compare_plan2.py)、[Recorder._row:type](../experiments/compare_plan2.py)
- [Recorder._row:purpose](../experiments/compare_plan2.py)、[Recorder.changed:inner_measure](../experiments/compare_plan2.py)、[Recorder.changed:outer_measure](../experiments/compare_plan2.py)、[Recorder.changed:gap](../experiments/compare_plan2.py)
- [settings:eps](../experiments/compare_plan2.py)、[settings:discovery_eps](../experiments/compare_plan2.py)、[settings:discovery_share](../experiments/compare_plan2.py)、[settings:network_eps](../experiments/compare_plan2.py)
- [write_run:t_cert](../experiments/compare_plan2.py)、[write_run:t_gap](../experiments/compare_plan2.py)、[run_once](../experiments/compare_plan2.py)、[evaluate_run](../experiments/compare_plan2.py)
