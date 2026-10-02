# 数学符号与代码变量契约

本文描述方法 RB 的主线（分支 method-RB）：径向锥夹逼 R 加逐网架支撑查询 B。相同数学量必须保留登记名称、单位、索引、状态布局和结果键；不能通过删登记绕过检查。新增量先登记，语义或输出格式变化必须显式迁移并同步调用方与回归。

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

## MP、SP、联合割与 OBBT 紧化

MP2 最大化 `direction @ p`，`objective/bound` 是 kW 的可行值/全局上界；MP1 最小化费用，`objective/bound` 是费用值/全局下界。成功答案包含 `x/p/state`；明确不可行返回 None。`radial_gap_kw=1e-3` 只检查目标间隙，不改写功率。原始求解质量不合格或无可靠答案时停止，不修补或重求。

SP 固定 x,p 后最小化非负 `eta`（求解器名 violation），只允许功率平衡与压降等式 ±eta，其余约束保持严格。`max(eta,0)+MaxVio<=PLANNING_TOL` 才接受原始 `state`。正 eta 通过锥的必要支撑平面 LP 生成有效联合割；`score_only=True` 保存 `cone_normals` 并延迟至 `generate_cut` 取割。超时、数值失败或无法分离均不能作为不可行证书。

唯一联合割方向是 `alpha + beta @ p + delta @ x >= 0`，数组布局 `[alpha, *beta, *delta]`，长度 `1+d+t`；beta 作用于 kW（分区内为幅值 u）。`_separating_cut.dual` 按求解器线性行排列，`coefficients` 按变量索引，`h` 提取全部 state 列。采用 Gurobi 行乘子方向，固定 x/p 的等式乘子先置零：`alpha=-dual@rhs + max(h,0)@y_ub_global + min(h,0)@y_lb_global`。正比例归一化后截距加 `1e-10` 保守补偿。`calls/cut_calls` 分别计 SP 与取割 LP 次数。

方案紧化：方案首次使用时 `GridPhysics.obbt` 在独立功率 `[0,bounds]` 内对选中型号的 P/Q/ell 与节点 v 做 `OBBT_ROUNDS` 轮 OBBT，盒存于 `GridPhysics.boxes[tuple(x)]`，键 `(量名, 型号键或节点)`、标幺；端点外扩 `OBBT_PAD`，未证得最优的端点保留原界（首轮为全局界）。`add_operation(..., scheme=x)` 对已紧化方案追加盒约束与反向锥包络割 `u_L*ell+ell_L*u-u_L*ell_L <= (P_L+P_U)P-P_L*P_U+(Q_L+Q_U)Q-Q_L*Q_U+ENVELOPE_MARGIN`（另一条取 u_U、ell_U，u 为参考送端电压平方），每行右端加 `M*H(x)`：H 为 x 与该方案的汉明距离，M 为该行在全局界上的最大违反量，故联合割对全部 x 仍有效。SP、射线与支撑查询（固定方案 MP）传入 scheme；锥 MISOCP 同时带全部已紧化方案的提升行。含紧化行的模型和 OBBT 用 `OBBT_CONV_TOL`。`obbt_extremes` 以 `settings.obbt_workers` 份模型副本并行求解，结果与串行逐位相同。

有限逐线路电流的 Port SP 固定使用既有 Gurobi 连续 SOCP；其他 Port SP 和射线使用 Clarabel。`solve_conic.basis/offset` 表示消元坐标，`cone_scale` 只是等价锥缩放。射线内域见证使用既有 `1e-6` 锥裕量，有限电流界同时扣除相同裕量；原物理残差仍按 `1e-8` 验收。`numeric_focus` 是求解设置，不改变数学可行性。

## 符号分区与几何坐标

主线只用符号分区：`PortPhysics.sign` 为 ±1 向量，内部非负幅值 u，真实 p=sign*u；负荷 PF=0.95、光伏 PF=1，允许反送。`main.mode=1` 只传给扫描接口（`vertify` 的缓存身份仍区分 mode）。输出、回放和扫描始终是真实 kW；支撑上界只对所属符号区的该网架有效。

`bounds=port_bounds(network)` 为公共正数坐标尺度 b（kW），内部点 `xi=u/bounds∈[0,1]^d` 无量纲；分区盒即 `[0,1]^d`。结果 `bound` 是查询目标的标量界，不能与 bounds 混用。叶锥、内三角形、P_x、O_x 与测度都在 xi 中计算，`build_partition` 返回前乘回 bounds 得 kW 幅值，`build_region` 再乘 sign 得带符号 kW。

`halfspaces/contains/clip_polytope` 不自动换算。面方程为 `F@xi+g<=0`；裁剪接口为 `constant+coefficient@xi>=0`。`covered` 先按包围盒筛点再逐行 `contains`，结果与逐行判定相同。`clip_box` 裁到 xi<=1（xi>=0 由锥保证），`cone_clip` 裁到锥 `U⁻¹xi>=0`，`box_vertices` 为分区盒顶点。

## 径向锥夹逼（R）

| 数学量 | 代码 | 单位与含义 |
|---|---|---|
| 方向 u_k | `Radial.directions`、`Radial.midpoints` | xi 中的单位向量，初始为坐标轴；棱 (a,b) 的中点方向由相邻锥共用，射线缓存随之共用 |
| 锥 K=cone(u_keys) | `Cone.keys`、`Radial.U` | d 个方向编号；U 的列为 u_i |
| x̂、y | `Cone.x`、`Cone.cover` | 内域网架；x̂ 原点不可行时的近端覆盖网架，可行时为 None |
| v_i、Q | `Cone.verts` | `(d,d)`，第 i 行 v_i=ρ_x̂(u_i)·u_i（xi）；Q=verts.T |
| c | `Cone.c` | c=1ᵀQ⁻¹，锥内 xi=Σα_i v_i 时 Σα_i=c@xi |
| μ̄ | `Cone.mu`、`Cone.inherited` | 外界因子：锥 MISOCP 的 ObjBound；未求解时沿用父锥 (c, μ̄) |
| T、O | `Cone.halfspace`、`cone_outer` | T=conv(0, v_1..v_d)；O={xi∈K: c@xi<=μ̄}=μ̄·T ⊇ R∩K |
| Δ | `Radial.delta` | (μ̄^d-1)·vol(T)，xi^d |
| ΣΔ/Σvol(T) | `Radial.volume_ratio` | 体积缺口比；尚无锥时为 inf |
| 射线与原点证书 | `Radial.rays`、`Radial.nears`、`Radial.origin`、`Radial.incumbents` | 键 (方案, 方向编号)；远端点与近端点为 xi；原点认证为布尔；锥 MISOCP 出现过的现任方案 |
| 审计过的现任解 | `Radial.incumbent_points`、`cone_misocp:audited` | (方案, xi)：已紧化方案的锥 MISOCP 现任解代回原约束（`audit_incumbent`）通过，作为 B 的认证点 |

内域证书 `Radial.inner`：x̂ 原点可行时 T ⊆ R_x̂；否则需原点可行的 y 在每个生成方向的射线半径不小于 x̂ 的近端半径，T ⊆ R_x̂ ∪ R_y。`Radial.vertex` 从原点朝分区盒边界点 `bounds*u/max(u)` 做紧化射线，半径小于 `MIN_RADIUS` 视为没有顶点；`Radial.near` 是从远端朝原点的反向射线。单次射线与原点 SP 的时限不超过 `RAY_SECONDS`。

外界 `cone_misocp`：x 自由的完整 MP（全部选型、拓扑、预算、运行约束），在 xi 上最大化 `objective`；`rows` 给锥约束 `rows@xi>=0` 与分区盒，缺省时为中心射线 xi_1=…=xi_d 并以 no-good 排除 `exclude`。全部已紧化方案的行按汉明距离提升，未紧化方案为纯 SOCP，故任何终止状态下的 ObjBound 都是 R 在该锥上的有效上界；`CONE_CONV_TOL` 是它的 barrier 收敛容差。返回 `status/bound/x/point/audited`（point 为 kW 幅值），尚无界时 bound=inf。`Radial.misocp` 是懒惰 OBBT：现任方案尚未紧化时先 OBBT 再重解。单次时限 `settings.mip_seconds`，相对间隙 `settings.mip_gap`；中心射线求解用 mip_gap=0。

`Radial.run(epsilon, check)` 用体积准则细分：每次取 Δ 最大的锥，ΣΔ<=ε·Σvol(T) 或 `check()` 成立即 certified。`Radial.options` 给候选剖分：解点方向的锥坐标 λ 全部 >= `SPLIT_MARGIN` 时星形剖分；三维恰有一个 λ 偏小时在对边上按 λ 投影二分；最后总有最长棱中点二分。`Radial.split` 的子锥 x̂ 取 {父 x̂, 解的方案, 父覆盖网架} 中证书成立、半径乘积最大者。锥角直径小于 `settings.min_width`、叶锥数达到 `settings.max_cones` 或没有内域证书的锥记 unresolved，不再细分。`Cone.status` 为 pending / bounded / unresolved；`Radial.run` 返回 certified / unresolved / time_limit，可续跑，阶段时限到达抛出 RegionTimeout 后保留已得证书。`Radial.schemes` 是 X*：叶锥的 x̂ 与覆盖网架、现任方案、已紧化与原点认证的方案，按所占锥体积从大到小。

## 逐网架支撑查询（B）

面方程沿用 `halfspaces` 的 `[n, b]`，|n|=1，内侧 n@xi+b<=0，记 β=-b。h_x(n)=max n@xi 为固定方案 x 的紧化可行集上的支撑值，UB 为其可靠上界。

| 数学量 | 代码 | 单位与含义 |
|---|---|---|
| V_x、P_x | `Network.points`、`Network.vertices`、`Network.full` | 网架 x 的认证点（xi）及其凸包顶点；full 表示 P_x 满维 |
| O_x | `Network.outer` | 分区盒逐次被可靠支撑上界 n@xi<=UB 裁剪的外界（xi），O_x ⊇ R_x |
| 面与状态 | `Network.faces`、`Network.boundary`、`Network.face_status`、`CERTIFIED_FACES` | P_x 的面；分区边界面掩码（`boundary_faces`，容差 `BOUNDARY_TOL`）；BOUNDARY / GEOMETRY_CERTIFIED / SUPPORT_CERTIFIED / UNRESOLVED / PENDING |
| 判据与 c_x | `Network.criterion`、`Network.center`、`choose_criterion`、`chebyshev_center` | auto：原点在 P_x 内且到每个非边界面的距离 >= `ORIGIN_CLEARANCE` 用 origin，否则 center；c_x 满维时为 Chebyshev 中心，低维时为顶点形心 |
| β' | `Network.allowed`、`allowed_offsets` | 允许的支撑值：origin 为 β/(1-tau)（E_x=P_x/(1-tau)），center 为 β+tau·(β-n@c_x)（E_x=c_x+(1+tau)(P_x-c_x)）；UB-β'<=GEOMETRY_TOL 即该面认证 |
| 优先级 | `Network.priority`、`Network.upper`、`face_measures` | 面积×max(UB-β, 0)；UB 取 O_x 顶点上的几何上界与缓存支撑上界的较小者；面积用 `FACE_TOL` 判定面上顶点 |
| 支撑缓存 | `Network.cache`、`normal_key`、`Network.lost` | 键为单位法向（不翻转符号）；复用时加键舍入差的保守补偿 Σmax(n-n_cached,0)；补点后 P_x 未变的面记入 lost，状态为 UNRESOLVED |
| ε_B | `settings.network_eps`、`Network.ratio` | 网架停止阈值 vol(O_x)/vol(P_x)-1<=ε_B，默认 ε/2 |
| 权重 | `Network.weight` | 所占叶锥体积加均值的 10%；分配本轮时间片 |

`SupportOracle` 是固定方案 x 的完整 MP（含 x 的 OBBT 盒与包络行）加分区盒 0<=u<=bounds；`SupportOracle.solve(normal, time_limit)` 的目标直接为 Σ a_i·u_i/b_i，ObjVal/ObjBound 即 xi 坐标的支撑值。只接受 OPTIMAL：ub 为 ObjBound 外扩 `BOUND_PAD`，解代回原约束（`audit_incumbent`：拓扑与预算、变量界、整数性、线性与二次行、原 Lorentz 范数，不修补）通过才返回认证点 point（xi）与 lb=n@point；lb>ub 时报错。其他状态返回 ub=inf、point=None，该面保持待查或记 UNRESOLVED。单次时限不超过 `SUPPORT_SECONDS`。

`Network.add` 处理一次支撑结果：缓存；有限 ub 裁剪 O_x；point 越过某个面的 β' 时补进 V_x 并重建 P_x（`Network.rebuild`，`Network.version` 加一，旧面结论作废）；最后 `Network.check` 判断停止：全部面认证（certified）、vol(O_x)/vol(P_x)-1<=ε_B（eps_B）、没有待查面（unresolved）、尚无认证点且种子方向 ±e_j、±1/√d 查完（no_points）。`Network.next_normal` 取优先级最高的待查面（尚无认证点时取未查的种子方向）。

`Support.run` 是一轮 B：X* 的网架以径向阶段的认证点为 V_x（`Support.seeds`：射线远端与近端、零接入点、审计过的已紧化现任解），按权重从大到小依次查询，各网架分得本轮剩余时间（分区剩余时间的 `PASS_SHARE`）按权重的份额，用完记 paused。`Support.query` 为一次查询并记录一帧；`Support.count` 为支撑次数；`Support.sets` 给出满维的 P_x 及其版本键。

内域 I=I_R ∪ ∪P_x（P_x ⊆ R_x 由认证点保证，是认证内域），`sandwich` 返回 (vol(I), vol(O_R∩盒))，单位 xi^d：二维为 shapely 并集；三维按叶锥分解（叶锥内部互不相交），每锥的并集体积按 (锥几何, 相交 P_x 的版本键) 缓存，并集数值失败时只计 T（内域只会低估）。`piece` 是 P_x 与 锥∩O∩盒 之交，零体积时为空。分区认证 `build_partition.certified`：vol(O_R∩盒)-vol(I)<=ε·vol(I)，ε=d·tau；或径向部分自身满足 ΣΔ<=ε·Σvol(T)。

## 分区流程、并行与结果

`build_partition` 依次为：A 阶段 `Radial.run(settings.discovery_eps)`，时限为分区时限的 `settings.discovery_share`；若体积缺口比仍大于 ε，B 对 X* 做一轮支撑查询；A+ 续跑 `Radial.run(ε, check)`，每次锥决策后检查夹逼判据（新出现的网架不再查询支撑）。返回分区摘要：

| 键 | 含义 |
|---|---|
| `build_partition:sign / status / certified` | 分区符号；certified / unresolved / time_limit；是否获证 |
| `build_partition:how` | radial（径向部分自身满足体积准则）或 support（夹逼判据） |
| `build_partition:seconds / gap / volume_ratio` | 分区墙钟秒；vol(O)/vol(I)-1；ΣΔ/Σvol(T) |
| `build_partition:cones / networks / accepted / supports` | 叶锥数；做过 B 的网架数；其中停止于 certified 或 eps_B 的数目；支撑次数 |
| `build_partition:inner / outer` | 内域块（各锥 T 与按锥裁剪的 P_x）与外域块（各锥 O∩盒，尚无锥时为分区盒）的顶点，kW 幅值 |

`build_region` 为每个分区起一个 spawn 子进程（`_partition`），并发数 min(`workers`, 2^d)，`settings.obbt_workers=max(1, workers//2^d)`；后启动的分区得到剩余时间按并发比例的份额，总时限为 `seconds`。子进程的 `RunMonitor(sign=...)` 经 `_connect` 设置的通道发送事件，主进程 `RunMonitor.share` 建立通道、`RunMonitor.forward` 转发、`RunMonitor.close` 发送结束标记；主进程中断时置取消信号并读空队列。结果 `inner/outer` 为带符号 kW 的 `[dict(vertices, sign)]`，`build_region:axis_lower/axis_bounds` 为外域顶点的范围，`build_region:partitions` 为各分区摘要（去掉 inner/outer），另含 status/certified/timing。

## 独立扫描、结果与回放

`vertify.py` 是独立 AC/SOCP 扫描入口。`ac_network` 复制完整配置，不删除限流。`budget_schemes` 仅供小算例独立 AC 参考的拓扑审计，构域不调用它。AC 单树只提供可行见证；全拓扑必要条件排除或完整 AC 不可行证书才给负标签。迭代失败、超时和未知均不可当作不可行。

`ACPowerFlow` 接收运行树，state(power,ell) 返回 `(P,Q,v,u)`，均为 `(batch,n_tree)`，v/u 为受端/送端电压平方。内部 `_state.p/q` 是节点标幺负荷；AC 残差为 `P²+Q²-u*ell`。`global_status` 供独立交叉检查。

`AC_CACHE_METHOD=ac_socp_grid_v4`，`ac_identity` 包含物理参数、ell_limit、模式、功率因数、有序节点和预算，不含构域 tau、时间和线程。`reference_box` 用同配置完整 SOCP 的方向全局上界覆盖两种参考。`scan_problem/scan_line` 逐点独立求解；SOCP 固定点使用 eta 目标，接受与排除依赖原门槛；AC 不读取 SOCP 标签。`main.run` 的扫描框为 reference_box 与结果外域范围的并。

| 缓存字段 | 固定含义 |
|---|---|
| `axis_lower / bounds` | 真实 kW 扫描边界 `(d,)` |
| `origin / step / start` | kW 原点/步长、整数起始索引；点为 `origin+(start+index+0.5)*step` |
| `states / witness_x / residual` | AC 标签、型号见证、残差；形状分别 `(n1,...,nd)`、`(n1,...,nd,t)`、`(n1,...,nd)` |
| `socp_states / socp_witness_x / socp_residual` | 同坐标 SOCP 的独立结果，形状同上 |
| `1 / -1 / 0` | 可行 / 已证不可行 / 未决；未决时不发布完整误差率 |
| `cache_path` | 实际 NPZ 路径；归档记录可使用相对项目根目录的路径 |

`scan_path` 给出 `results/scan/<网络>/<节点>/<身份>/`，`region_path` 由 origin/step/start/shape 生成区域摘要。扩界保留格点，部分覆盖只补算缺点；AC/SOCP 分别补零标签。`import_ac_reference` 仅接受同版本同身份数据。`force_rescan` 先备份；主进程每 `SAVE_SECONDS` 至多原子落盘一次，结束或中断时必落盘，独占锁避免并发覆盖。

比较使用最终 inner 作为计算域，另保留 outer 指标。`comparison_metrics` 的 MR=`missed/reference`、FR=`extra/computed`，百分数，空分母 None；同时保留 missed_cells/extra_cells/reference_cells/computed_cells。`comparisons` 键为 result_ac/result_socp/socp_ac。导出版本 paired_scan_v2，NPZ 的 power 为 `(N,d)` kW，ac_states/socp_states/inner/outer 对应同坐标；CSV 标签列为 ac_state/socp_state。网格误差不是连续体积证明。RB 的内域是认证内域，相对 AC 的 FR 来自 OBBT 紧化模型与 AC 的差（SOCP 松弛残余）。

回放 version=4：history 保存增量过程，validation_state 保存最终扫描；version=3 兼容读取。子进程事件为 phase_start、cone、point、cut、partition_end，主进程另有 start、region_end。`MERGED` 中的 schemes/cones/cut_history 按键增量合并，键带分区前缀 `<分区>:`（如 `+-:3`），坐标在 `signed_values` 中乘 sign：

| 状态 | 行内容 |
|---|---|
| `cones` | `Radial.publish:inner / outer / scheme / mu`：T 与 O∩盒的顶点、x̂ 标签、μ̄；被细分的父锥置 None |
| `schemes` | `Support.row:x / choice / cost / outer / inner / status` 与 sign：outer 为 O_x，inner 为满维的 P_x，否则为空 |
| `cut_history` | 支撑上界写成的条件割 `Support.query:cut`（布局同联合割，beta 已乘 sign；右端加 M·汉明距离，M 为盒上的最大违反量，只对来源网架起作用）、来源网架 `scheme`、`sign` |

`cone` 帧另有 cone_count、volume_ratio 与 global_point（锥 MISOCP 的解点）；支撑查询给出有限上界时为 `cut` 帧，否则为 `point` 帧，两者都有 `Support.query:support`（certified/faces/ratio/status：已认证面数、面数、vol(O_x)/vol(P_x)-1、网架状态）与 sp_point（审计过的支撑点）。仅当前帧之前的锥、网架、割与点出现在过程图；最终扫描独立显示。result 包含 status/certified/inner/outer/partitions/timing，分区证书在 partition_end；数值未决和时限未完不能改成 certified。JSON 非有限值为 null，仅预算 null 可按 inf 解释。

## 方法归档

方法对照实验的精选运行归档在 `results/methods/<算例>/<方法>/workers_<w>/run_<r>/`（summary.json、metrics.json、timeline.csv、solves.csv.gz、snapshots/final.json.gz），登记于 `results/manifest.json` 的 `method_archive`、`runs`、`files`；各方法的代码版本由 tag `method-<方法>-v1` 固定，归档整体为 `results-methods-v1`。主线分支 main 为方法 RCUT（分支 method-RCUT）；本分支 method-RB 是方法 RB 的主线实现。

## 容差与默认设置

| 配置 | 当前值 |
|---|---|
| `PLANNING_TOL / GEOMETRY_TOL` | 各 1e-8；物理残差/归一化几何门槛，语义独立 |
| `OBBT_ROUNDS / OBBT_PAD / ENVELOPE_MARGIN / OBBT_CONV_TOL` | 2 / 1e-6 标幺 / 1e-4 标幺² / 1e-8 |
| `CONE_CONV_TOL` | 1e-6；锥 MISOCP 的 barrier 收敛容差 |
| `AC_TOL / FIXED_POINT_TOL / GLOBAL_AC_TOL` | 1e-9 / 1e-12 / 1e-7 |
| `REGION_TAU` | 0.005；体积目标 ε=d·tau，不是物理容差 |
| `NETWORK / DIMENSION / mode` | Case33 / 3 / 1 |
| `CASE_TIME_LIMIT / WORKERS / SOLVER_THREADS` | 300 秒 / 16 / 1；总时限不含事后扫描 |
| `DISCOVERY_EPS / DISCOVERY_SHARE` | 0.15 / 0.25 |
| `MIP_SECONDS / MIP_GAP` | 60 秒 / 1e-3 |
| `MIN_WIDTH / MAX_CONES` | {2:1e-4, 3:2e-3} rad / {2:256, 3:2048} |
| `SPLIT_MARGIN / MIN_RADIUS / RAY_SECONDS` | 0.05 / 1e-9（xi）/ 30 秒 |
| `SUPPORT_SECONDS / PASS_SHARE / settings.network_eps` | 30 秒 / 0.75 / ε/2 |
| `BOUNDARY_TOL / ORIGIN_CLEARANCE / FACE_TOL / BOUND_PAD` | 1e-7 / 1e-6 / 1e-9（xi）/ 1e-10（相对） |
| `DIVISIONS / SCAN_DIVISIONS / SCAN_WORKERS / SAVE_SECONDS` | 160 / {2:160,3:80} / 20 / 30 秒；每扫描进程一个求解线程 |

单次 MP/SP/AC 的时限仍用 MP_TIME_LIMIT/SP_TIME_LIMIT/AC_TIME_LIMIT；AC_ITERATIONS 为迭代次数。不得合并不同语义的容差。数值配置变化须同步本表和结果协议。

## 显式迁移（2026-10-02：method-RB）

本分支自 method-RCUT 的 77c9043 分出，把 CUT 换成 B；R、分区并行、前端与扫描不变。

- 去掉 RCUT 的割平面构域：`Cutting`、`CutSlice`、`stagnated`、`register_power`、`CUT_SECONDS`、`CUT_ACCEPTED`、`build_partition:cuts`，以及入口的 `CUT_THRESHOLD`、`CUT_PATIENCE`、`POINT_TOL`；它们留在 main / method-RCUT。分区摘要的支撑次数为 `build_partition:supports`，`accepted` 改为停止于 certified 或 eps_B 的网架数。
- B 的接口自方法对照实验（tag `method-RB-v1`，experiments/compare_plan2.py、plan2_geometry.py、test_support_face_certification_fourbus_2d.py）迁入并精简：`BOUND_PAD`、`audit_incumbent` 语义不变；`SupportOracle` 合并了实验的 PartitionOracle（分区盒、只接受 OPTIMAL），`SupportOracle.solve` 改为 (normal, time_limit)，返回 lb/ub/point/status/normal，不再自带缓存（由 `Network.cache` 按同一 `normal_key` 缓存）；`NetworkState` 改为 `Network`（points/vertices/outer/criterion/center 含义不变），`SupportPhase` 的 B 部分改为 `Support`（一轮、顺序查询，无覆盖证书 C 与进程池）。
- 面判据统一为 UB-β'<=GEOMETRY_TOL：origin 判据原为 (1-tau)·UB+b<=GEOMETRY_TOL，两者只差 GEOMETRY_TOL·tau；`classify_support`、`classify_center`、`face_margins` 并入 `Network.evaluate`，近退化几何的兜底（`GEOMETRY_ERRORS`）去掉，数值失败直接报错。
- `cone_misocp` 增加返回键 audited，`Radial.incumbent_points` 收集审计过的已紧化现任解作为 V_x 的种子。前端阶段名为“网架支撑”，网架面板的灰色为 O_x、绿色为 P_x，支撑上界以条件割显示。

## 显式迁移（2026-10-02：RCUT 成为主线）

主线算法由逐网架顺序构域加完整物理查漏，改为 R + CUT。退役接口及其登记保存在 tag `mainline-sequential-v1`（3c061a2），对应关系如下：

- 构域：`build_sequential_region` 及其局部量（cache/powers/applied/initialized/counts/small_cuts/area_ratio/point_tol/ray_threshold 等）、`RegionState`、`stage_candidates`、`ray_gain`、`coverage_halfspaces` 由 `Radial`（全局搜索与外界）和 `Cutting`（逐网架 N_x）代替；割循环的 SP 评分、取割与停滞规则不变，局部状态现属 `Cutting.cut`，area_ratio 的序列为 `Cutting.cut.history`。主线不再做割后补射线与边界补充，`RAY_THRESHOLD` 无对应量。
- 全局覆盖：`remaining_search`、`RemainingRegionModel`、`COVERAGE_TOL`、`RESIDUAL_TIME_LIMIT`、`GridPhysics.center` 由 x 自由的锥 MISOCP `cone_misocp` 给出的外界 O_R 与体积准则代替；`MasterProblem.__init__.cuts` 随之去掉，MP 不再携带联合割。
- 时限与并行：`PARTITION_TIME_LIMIT`、`build_region.partition_seconds`、`run.partition_seconds` 由分区并行共享总时限 `build_region.seconds` 代替；`OBBT_WORKERS`、`run.obbt_workers`、`build_region.obbt_workers` 由 `WORKERS` 推出 `settings.obbt_workers=max(1, WORKERS//2^d)`。
- 默认值：`CUT_THRESHOLD` 语义不变，由 0.02 改为 0.01（方法对照的 RCUT 设置）；`SOLVER_THREADS` 由 20 改为 1（分区已并行）。主入口去掉 mode=0 非负负荷与 `--mode` 参数。
- 回放：version=4 格式不变，事件改为本文所列；`RunMonitor._geometry`、`global_end`、`partition_frame`、`seed`、`sp_start`、`updated` 由子进程事件与 `RunMonitor.forward` 代替。网架行增加 status，键带分区前缀。旧主线记录只能用 `mainline-sequential-v1` 的前端回放。
- 实验：`experiments/` 下支持面认证、径向夹逼与方案二对照（compare_plan2、plan2_geometry）的代码及登记保存在 tag `results-methods-v1`（ed379fb）与各 `method-<方法>-v1`；旧主线、支持面与径向夹逼的结果归档也在这些 tag 中，主线分支只保留方法对照归档 `results/methods` 与 Case33 参考扫描。

## 2026-09-30 整理的显式迁移

结果统一在 results/mainline、results/support_face、results/scan；映射、原路径、文件哈希和来源见 results/manifest.json。Case33 记录的 cache_path/ac_cache 改为相对项目根目录的 results/scan 路径，main.run 同步解析；求解数据和扫描数组不变。FourBus 已保存参考仍为历史 SOCP 格点，不能冒充当前配对 AC/SOCP 缓存。

退役内容为 Concept5、load_jiangkou/Jiangkou、legacy_case33 升级夹具、fixed_topology/upgrade_plan/affordable_designs/dispatch_state 测试辅助、旧前端测试和一次性数值诊断脚本。其登记随接口退役归档于 Git 提交 e79f18b。

交付运行 `python -m unittest tests.test_notation -v`；算法修改另跑相关数值回归。

## 当前登记接口索引

以下登记当前全部现存接口，以便名称检查覆盖字段、形参和结果键。`作用域.名称` 表示字段或局部量，`函数:键` 表示返回字典键；语义按上文分组约定。

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
- [signed_ac_witness](../vertify.py)、[SAVE_SECONDS](../vertify.py)

### monitor.py

- [COMPARISONS](../monitor.py)、[MERGED](../monitor.py)、[RegionTimeout](../monitor.py)、[RunMonitor](../monitor.py)
- [RunMonitor.__init__.sign](../monitor.py)、[RunMonitor.label](../monitor.py)、[RunMonitor.channel](../monitor.py)、[RunMonitor.shared](../monitor.py)
- [RunMonitor.share](../monitor.py)、[RunMonitor.close](../monitor.py)、[RunMonitor.forward](../monitor.py)、[RunMonitor.forward.label](../monitor.py)
- [RunMonitor.validation](../monitor.py)、[RunMonitor.validation.region_key](../monitor.py)、[RunMonitor.validation:comparisons](../monitor.py)、[RunMonitor.validation:metrics](../monitor.py)
- [RunMonitor.validation_state](../monitor.py)、[_connect](../monitor.py)、[signed_values](../monitor.py)、[_cap](../monitor.py)
- [_cut_polygon](../monitor.py)、[_cut_polygon.axis_lower](../monitor.py)、[_cut_segment](../monitor.py)、[_cut_segment.axis_lower](../monitor.py)
- [_draw_3d](../monitor.py)、[_voxel_faces](../monitor.py)、[comparison_metrics](../monitor.py)、[comparison_metrics:computed_cells](../monitor.py)
- [comparison_metrics:extra_cells](../monitor.py)、[comparison_metrics:fr_percent](../monitor.py)、[comparison_metrics:missed_cells](../monitor.py)、[comparison_metrics:mr_percent](../monitor.py)
- [comparison_metrics:reference_cells](../monitor.py)

### model.py

- [DEFAULT_SOLVER_THREADS](../model.py)、[GridPhysics](../model.py)、[GridPhysics.P_slice](../model.py)、[GridPhysics.Q_slice](../model.py)
- [GridPhysics.add_operation](../model.py)、[GridPhysics.add_operation.P](../model.py)、[GridPhysics.add_operation.Q](../model.py)、[GridPhysics.add_operation.ell](../model.py)
- [GridPhysics.add_operation.minus](../model.py)、[GridPhysics.add_operation.plus](../model.py)、[GridPhysics.add_operation.v](../model.py)、[GridPhysics.cost](../model.py)
- [GridPhysics.drop_max](../model.py)、[GridPhysics.ell_slice](../model.py)、[GridPhysics.ellmax](../model.py)、[GridPhysics.ends](../model.py)
- [GridPhysics.incoming](../model.py)、[GridPhysics.keys](../model.py)、[GridPhysics.outgoing](../model.py)、[GridPhysics.pmax](../model.py)
- [GridPhysics.pmin](../model.py)、[GridPhysics.qmax](../model.py)、[GridPhysics.qmin](../model.py)、[GridPhysics.r](../model.py)
- [GridPhysics.reactance](../model.py)、[GridPhysics.slack_slice](../model.py)、[GridPhysics.types](../model.py)、[GridPhysics.v_slice](../model.py)
- [GridPhysics.vmax](../model.py)、[GridPhysics.vmin](../model.py)、[GridPhysics.y_lb_global](../model.py)、[GridPhysics.y_ub_global](../model.py)
- [BOUND_PAD](../model.py)、[audit_incumbent](../model.py)、[SupportOracle](../model.py)、[SupportOracle.solve](../model.py)
- [SupportOracle.solve.time_limit](../model.py)、[SupportOracle.solve:lb](../model.py)、[SupportOracle.solve:ub](../model.py)、[SupportOracle.solve:point](../model.py)
- [SupportOracle.solve:normal](../model.py)、[cone_misocp:audited](../model.py)、[MP_TIME_LIMIT](../model.py)、[MasterProblem.__init__.budget](../model.py)、[MasterProblem.__init__.direction](../model.py)、[MasterProblem.__init__.f](../model.py)
- [MasterProblem.__init__.min_total](../model.py)、[MasterProblem.__init__.power](../model.py)、[MasterProblem.active_nodes](../model.py)、[MasterProblem.choices](../model.py)
- [MasterProblem.direction](../model.py)、[MasterProblem.loads](../model.py)、[MasterProblem.operation](../model.py)、[MasterProblem.power](../model.py)
- [MasterProblem.solve.incumbent](../model.py)、[MasterProblem.solve.radial_gap_kw](../model.py)、[MasterProblem.solve.start](../model.py)、[MasterProblem.solve:bound](../model.py)
- [MasterProblem.solve:feasible](../model.py)、[MasterProblem.solve:objective](../model.py)、[MasterProblem.solve:p](../model.py)、[MasterProblem.solve:state](../model.py)
- [MasterProblem.solve:status](../model.py)、[MasterProblem.solve:x](../model.py)、[MasterProblem.state](../model.py)、[MasterProblem.x](../model.py)
- [PLANNING_TOL](../model.py)、[PortPhysics.sign](../model.py)、[SP_TIME_LIMIT](../model.py)、[SubProblem._build](../model.py)
- [SubProblem._build.eta](../model.py)、[SubProblem._cut](../model.py)、[SubProblem._separating_cut.coefficients](../model.py)、[SubProblem._separating_cut.cut](../model.py)
- [SubProblem._separating_cut.dual](../model.py)、[SubProblem._separating_cut.h](../model.py)、[SubProblem._separating_cut.rows](../model.py)、[SubProblem.cut_calls](../model.py)
- [SubProblem.generate_cut](../model.py)、[SubProblem.numeric_focus](../model.py)、[SubProblem.solve](../model.py)、[SubProblem.solve.cone_normals](../model.py)
- [SubProblem.solve.eta](../model.py)、[SubProblem.solve.score_only](../model.py)、[SubProblem.solve:cone_normals](../model.py)、[SubProblem.solve:cut](../model.py)
- [SubProblem.solve:eta](../model.py)、[SubProblem.solve:feasible](../model.py)、[SubProblem.solve:state](../model.py)、[port_bounds](../model.py)
- [ray_support](../model.py)、[ray_support.numeric_focus](../model.py)、[ray_support.ray_fraction](../model.py)、[solve_conic.basis](../model.py)
- [solve_conic.cone_scale](../model.py)、[solve_conic.offset](../model.py)、[voltage_flow_bounds](../model.py)、[OBBT_ROUNDS](../model.py)
- [OBBT_PAD](../model.py)、[ENVELOPE_MARGIN](../model.py)、[OBBT_CONV_TOL](../model.py)、[OBBT_PARAMS](../model.py)
- [obbt_pool](../model.py)、[obbt_extremes](../model.py)、[GridPhysics.boxes](../model.py)、[GridPhysics.obbt](../model.py)
- [GridPhysics.add_operation.scheme](../model.py)、[CONE_CONV_TOL](../model.py)、[cone_misocp](../model.py)、[cone_misocp.objective](../model.py)
- [cone_misocp.rows](../model.py)、[cone_misocp.exclude](../model.py)、[cone_misocp.mip_gap](../model.py)、[cone_misocp:status](../model.py)
- [cone_misocp:bound](../model.py)、[cone_misocp:x](../model.py)、[cone_misocp:point](../model.py)

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

### main.py

- [BUDGET](../main.py)、[CASE_TIME_LIMIT](../main.py)、[DIMENSION](../main.py)、[DIVISIONS](../main.py)
- [FORCE_RESCAN](../main.py)、[NETWORK](../main.py)、[REGION_TAU](../main.py)、[SCAN_DIVISIONS](../main.py)
- [SCAN_OUTPUT](../main.py)、[SCAN_WORKERS](../main.py)、[SOLVER_THREADS](../main.py)、[WORKERS](../main.py)
- [DISCOVERY_EPS](../main.py)、[DISCOVERY_SHARE](../main.py)、[MIP_SECONDS](../main.py)、[MIP_GAP](../main.py)
- [MIN_WIDTH](../main.py)、[MAX_CONES](../main.py)、[main](../main.py)、[main.case](../main.py)
- [main.dimension](../main.py)、[main.load_nodes](../main.py)、[mode](../main.py)、[recording_path](../main.py)
- [run](../main.py)、[run.budget](../main.py)、[run.force_rescan](../main.py)、[run.scan_output](../main.py)
- [run.workers](../main.py)、[run.threads](../main.py)、[run.settings](../main.py)

### region.py

- [GEOMETRY_TOL](../region.py)、[SPLIT_MARGIN](../region.py)、[MIN_RADIUS](../region.py)、[RAY_SECONDS](../region.py)
- [SUPPORT_SECONDS](../region.py)、[PASS_SHARE](../region.py)、[BOUNDARY_TOL](../region.py)、[ORIGIN_CLEARANCE](../region.py)
- [FACE_TOL](../region.py)、[CERTIFIED_FACES](../region.py)、[halfspaces](../region.py)、[covered](../region.py)
- [clip_polytope.coefficient](../region.py)、[clip_polytope.constant](../region.py)、[union_measure](../region.py)、[box_vertices](../region.py)
- [clip_box](../region.py)、[cone_clip](../region.py)、[cone_outer](../region.py)、[piece](../region.py)
- [sandwich](../region.py)、[boundary_faces](../region.py)、[face_measures](../region.py)、[chebyshev_center](../region.py)
- [choose_criterion](../region.py)、[allowed_offsets](../region.py)、[normal_key](../region.py)、[Cone](../region.py)
- [Cone.keys](../region.py)、[Cone.x](../region.py)、[Cone.verts](../region.py)、[Cone.cover](../region.py)
- [Cone.inherited](../region.py)、[Cone.mu](../region.py)、[Cone.final](../region.py)、[Cone.status](../region.py)
- [Cone.c](../region.py)、[Cone.halfspace](../region.py)、[Radial](../region.py)、[Radial.bounds](../region.py)
- [Radial.directions](../region.py)、[Radial.midpoints](../region.py)、[Radial.rays](../region.py)、[Radial.nears](../region.py)
- [Radial.origin](../region.py)、[Radial.incumbents](../region.py)、[Radial.incumbent_points](../region.py)、[Radial.cones](../region.py)
- [Radial.deadline](../region.py)、[Radial.U](../region.py)、[Radial.vertex](../region.py)、[Radial.near](../region.py)
- [Radial.inner](../region.py)、[Radial.misocp](../region.py)、[Radial.solve](../region.py)、[Radial.options](../region.py)
- [Radial.split](../region.py)、[Radial.delta](../region.py)、[Radial.volume_ratio](../region.py)、[Radial.run](../region.py)
- [Radial.run.epsilon](../region.py)、[Radial.run.check](../region.py)、[Radial.schemes](../region.py)、[Radial.geometry](../region.py)
- [Radial.publish](../region.py)、[Radial.publish:inner](../region.py)、[Radial.publish:outer](../region.py)、[Radial.publish:scheme](../region.py)
- [Radial.publish:mu](../region.py)、[Network](../region.py)、[Network.points](../region.py)、[Network.vertices](../region.py)
- [Network.full](../region.py)、[Network.outer](../region.py)、[Network.faces](../region.py)、[Network.boundary](../region.py)
- [Network.face_status](../region.py)、[Network.criterion](../region.py)、[Network.center](../region.py)、[Network.allowed](../region.py)
- [Network.areas](../region.py)、[Network.upper](../region.py)、[Network.priority](../region.py)、[Network.cache](../region.py)
- [Network.lost](../region.py)、[Network.status](../region.py)、[Network.version](../region.py)、[Network.weight](../region.py)
- [Network.seeds](../region.py)、[Network.rebuild](../region.py)、[Network.evaluate](../region.py)、[Network.next_normal](../region.py)
- [Network.add](../region.py)、[Network.ratio](../region.py)、[Network.check](../region.py)、[Support](../region.py)
- [Support.networks](../region.py)、[Support.oracles](../region.py)、[Support.count](../region.py)、[Support.seeds](../region.py)
- [Support.run](../region.py)、[Support.query](../region.py)、[Support.query:support](../region.py)、[Support.query:certified](../region.py)
- [Support.query:faces](../region.py)、[Support.query:ratio](../region.py)、[Support.query:status](../region.py)、[Support.query:sp_point](../region.py)
- [Support.query:cut](../region.py)、[Support.row](../region.py)、[Support.row:x](../region.py)、[Support.row:choice](../region.py)
- [Support.row:cost](../region.py)、[Support.row:outer](../region.py)、[Support.row:inner](../region.py)、[Support.row:status](../region.py)
- [Support.sets](../region.py)、[build_partition](../region.py)、[build_partition.certified](../region.py)、[build_partition:sign](../region.py)
- [build_partition:status](../region.py)、[build_partition:certified](../region.py)、[build_partition:how](../region.py)、[build_partition:seconds](../region.py)
- [build_partition:gap](../region.py)、[build_partition:volume_ratio](../region.py)、[build_partition:cones](../region.py)、[build_partition:networks](../region.py)
- [build_partition:accepted](../region.py)、[build_partition:supports](../region.py)、[build_partition:inner](../region.py)、[build_partition:outer](../region.py)
- [_partition](../region.py)、[build_region](../region.py)、[build_region.workers](../region.py)、[build_region.seconds](../region.py)
- [build_region.settings](../region.py)、[build_region:axis_lower](../region.py)、[build_region:axis_bounds](../region.py)、[build_region:partitions](../region.py)

### plot.py

- [_clip_face](../plot.py)

### Network/four_bus_five_corridor.py

- [FourBus.__init__.load_nodes](../Network/four_bus_five_corridor.py)、[LineType.cable_cny_m](../Network/four_bus_five_corridor.py)、[LineType.capacity_kw](../Network/four_bus_five_corridor.py)、[LineType.r_ohm_km](../Network/four_bus_five_corridor.py)
- [LineType.x_ohm_km](../Network/four_bus_five_corridor.py)

### tests/planning_checks.py

- [margin](../tests/planning_checks.py)、[recorded_monitor](../tests/planning_checks.py)
