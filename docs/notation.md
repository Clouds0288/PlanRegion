# 数学符号与代码变量规范

版本：2.0，2026-09-24。适用于正式两阶段 SOCP 可行规划域主线。同一数学量语义不变时，保留名称、单位、维数及索引。全部历史登记保留在 [v1 契约](notation-v1.md)，迁移见第 9 节；历史文件不是现行 API。

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
| \(c,\mathcal B\) | 型号投资向量、预算 | `network.cost`、`budget` | `network.cost_unit` |
| \(S_{\mathrm b}\) | 功率基值 | `network.base` | kVA |

禁止把电抗称为 `x`、把电流平方写成 `l`、把预算记作标准式矩阵 `B`、把归一化坐标记作走廊状态 `z`。源数据里的 `x_ohm_km` 等字段按第 3 节处理。

## 2. 集合、维数和唯一索引

数学下标从 1 起写，Python 数组从 0 起；真实节点 ID 由数据给出，不能直接当数组位置。根节点在公式中记为 0，内部哨兵为 `-1`，不占 `nodes`、`v`、`a` 的元素。

| 符号 | 定义 / 顺序 | 代码定义 |
|---|---|---|
| \(\mathcal N,n\)；\(i,j\) | 全部候选非根节点及其索引 | [Network.nodes](../Network/__init__.py)、[Network.n](../Network/__init__.py)、[Network.node_index](../Network/__init__.py)、[Network.root](../Network/__init__.py) |
| \(\mathcal E,m\)；\(e\) | 候选走廊，按输入顺序 | [Network.corridors](../Network/__init__.py)、[Network.n_corridors](../Network/__init__.py) |
| \(\mathcal T,t\)；\(\nu\) | 各走廊内型号连续展开；\(\nu\) 是扁平数组位置，总数 t | [Network.type_keys](../Network/__init__.py)、[Network.n_types](../Network/__init__.py) |
| \((e,k)\) | Gurobi 建模键：走廊 ID 与该走廊内型号 ID；k 不要求全图唯一 | [PlanningEquations.keys](../model.py)、[PlanningEquations.types](../model.py) |
| \(\mathcal T_e\) | 走廊 e 的型号区间；型号到走廊的映射 | [Network.type_slices](../Network/__init__.py)、[Network.type_corridor](../Network/__init__.py) |
| \(\mathcal L,d\)；\(h\) | 独立负荷节点，`d = len(network.load_nodes)` | [Network.load_nodes](../Network/__init__.py)、[Network.selected](../Network/__init__.py) |
| \(\mathcal N_{\rm req}\) | 必须接入节点的布尔掩码，`(n,)` | [Network.required](../Network/__init__.py) |
| \(n_y=3t+n+2m\) | 标准运行向量长度 | [PlanningEquations.slack_slice](../model.py) 的 `.stop` |
| \(n_{\rm tree}\) | 当前方案实际接入的非根节点数，也是树的支路数 | [OperatingTree.n](../Network/__init__.py) |

`network.type_keys[nu] == (e,k)` 是具名键与扁平位置的唯一桥梁。求解器内 `p[i]` 用真实负荷节点 ID；数组 `power[h]` 用 `load_nodes[h]` 的位置。`n/m/t/d` 在数学文档中固定为上表维数；`m = model`、`e = equations` 等既有局部对象别名不重定义这些数学符号。

模型装配处理 d 个负荷坐标；正式几何与展示支持 d=2、3。独立旧 AC 网格工具仍是三维。

## 3. 网架参数、量纲和数据入口

| 符号 / 含义 | 代码定义 | 单位 / 形状 |
|---|---|---|
| \(S_{\mathrm b}\)，功率基值 | [Network.base](../Network/__init__.py) | kVA，标量；不是自动等于电源容量 |
| \(V_{\mathrm b}\)，线电压基值 | [Network.voltage_kv](../Network/__init__.py) | kV，标量 |
| \(r_{e,k}\)，数组 \(r_\nu\) | [TypeParameters.r](../Network/__init__.py) → [Network.r](../Network/__init__.py)；具名 [PlanningEquations.r](../model.py) | p.u.，标量 → `(t,)` / 按键字典 |
| \(\chi_{e,k}\)，数组 \(\chi_\nu\) | [TypeParameters.reactance](../Network/__init__.py) → [Network.reactance](../Network/__init__.py)；具名 [PlanningEquations.reactance](../model.py) | p.u.，标量 → `(t,)` / 按键字典 |
| \(P_{e,k}^{\max}\)，输入的送端有功上限 | [TypeParameters.capacity](../Network/__init__.py) → [Network.capacity](../Network/__init__.py) | p.u.，标量 → `(t,)`；不是电流或视在功率 |
| \(c_{e,k}\)，增量投资 | [TypeParameters.investment_cost](../Network/__init__.py) → [Network.cost](../Network/__init__.py)；具名 [PlanningEquations.cost](../model.py) | `cost_unit`，标量 → `(t,)` / 按键字典 |
| 费用单位、默认预算列表 | [Network.cost_unit](../Network/__init__.py)、[Network.budgets](../Network/__init__.py) | FourBus / 江口为元；Case33 为相对投资单位 |
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
| \(\rho_\nu\) | 数学推导量，对应参考送端是否为根 | `(t,)`；根出线集合由 [PlanningEquations.outgoing](../model.py) 给出 |
| 参考端点、节点入边 / 出边 | [PlanningEquations.ends](../model.py)、[PlanningEquations.incoming](../model.py)、[PlanningEquations.outgoing](../model.py) | 按走廊 / 节点 ID 的字典 |
| \(a\) | [PlanningModel.active_nodes](../model.py)，局部 `a` | `(n,)`，二进制；根恒接入 |
| \(f\) | [PlanningModel.__init__.f](../model.py) | `(m,)`；连通虚拟流，不是电功率 |
| \(x\) | [PlanningModel.x](../model.py) | `(t,)`，完整型号顺序；无第二套压缩选型 |

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

[PlanningEquations](../model.py) 固定状态布局为

\[
y=(P,Q,\ell,v,s^+,s^-),\qquad n_y=3t+n+2m.
\]

建模使用 [PlanningEquations.add_operation](../model.py) 中的 `P[e,k]`、`Q[e,k]`、`ell[e,k]`、`v[i]`、`plus[e]`、`minus[e]`。它们返回为 `problem.operation` 的同名字段；`state` 是同一批变量的扁平视图，没有第二批运行变量。

| 具名对象 | 固定代码定义 | 与向量的关系 |
|---|---|---|
| 建设选型、负荷 | [PlanningModel.choices](../model.py)、[PlanningModel.loads](../model.py) | `[e,k]` / `[i]`；[PlanningModel.x](../model.py)、[PlanningModel.power](../model.py) 是共享变量的 MVar 视图 |
| 有功、无功、电流平方、电压平方 | [PlanningEquations.add_operation.P](../model.py)、[PlanningEquations.add_operation.Q](../model.py)、[PlanningEquations.add_operation.ell](../model.py)、[PlanningEquations.add_operation.v](../model.py) | 根电压 `v[network.root] = 1` 是常数，不进入状态向量 |
| 正、负开断压降余量 | [PlanningEquations.add_operation.plus](../model.py)、[PlanningEquations.add_operation.minus](../model.py) | 求解器名称 `drop_plus/drop_minus`；数学符号 \(s^+/s^-\) |
| 完整运行对象、状态 | [PlanningModel.operation](../model.py)、[PlanningModel.state](../model.py) | 始终存在；不再提供省略物理约束的模式 |

不要因采用 `tupledict`、`MVar` 或 NumPy 就改数学量名称；转换严格使用 `type_keys`、`nodes`、`load_nodes`、`corridors` 顺序。

| 分块 | 唯一切片 | Python 范围 | 含义 |
|---|---|---|---|
| \(P\) | [PlanningEquations.P_slice](../model.py) | `[0:t]` | 参考送端有功，p.u. |
| \(Q\) | [PlanningEquations.Q_slice](../model.py) | `[t:2*t]` | 参考送端无功，p.u. |
| \(\ell\) | [PlanningEquations.ell_slice](../model.py) | `[2*t:3*t]` | 电流幅值平方，p.u.² |
| \(v\) | [PlanningEquations.v_slice](../model.py) | `[3*t:3*t+n]` | 非根电压幅值平方，p.u.² |
| \(s^+,s^-\) | [PlanningEquations.slack_slice](../model.py) | `[3*t+n:3*t+n+2*m]` | 前 m 为正、后 m 为负；开断压降余量，p.u.² |

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

这些是数学块名，**当前实现没有 `equations.A/B/C/b`、`lb_x/ub_x` 等属性**；不得为了满足公式字面而重建一套平行模型。当前物理约束由 add_operation 的具名约束给出。新主线使用完整整数优化的全局目标界，不导出线性行乘子。

| 数学量 | 代码定义 | 形状 / 说明 |
|---|---|---|
| \(l^{\rm glob},u^{\rm glob}\)，全局变量盒 | [PlanningEquations.y_lb_global](../model.py)、[PlanningEquations.y_ub_global](../model.py) | `(n_y,)`；保留全部型号，与 x 无关 |
| \(\overline P,\overline Q,\overline\ell\) | [PlanningEquations.pmax](../model.py)、[PlanningEquations.qmax](../model.py)、[PlanningEquations.ellmax](../model.py) | 按 `(e,k)` 的派生有效界 |
| \(\underline P,\underline Q\) | [PlanningEquations.pmin](../model.py)、[PlanningEquations.qmin](../model.py) | 根出线为 0，其余为负的有效上界 |
| \(\underline v,\overline v\) 的具名视图 | [PlanningEquations.vmin](../model.py)、[PlanningEquations.vmax](../model.py) | 节点 ID 字典，含根节点常数 1 |
| \(M_e\)，开断压降界 | [PlanningEquations.drop_max](../model.py) | 走廊 ID 字典；\(0\le s_e^\pm\le M_e(1-z_e)\) |
| 求解结果最大违反量 | Gurobi `model.MaxVio` | 直接读取所建模型的数值质量；MP 接受条件为 `<= PLANNING_TOL` |
| 测试用原约束最小余量 | [margin](../tests/planning_checks.py) | 保留独立代入检查，仅用于测试和审核，不参与生产求解 |

全局盒不是完整可行域：其中电压下界为 0，实际 `v >= vmin` 另作运行约束。各约束有自己的物理尺度，不能给整个 `y` 或 `w` 一概标注 kW。

派生界固定为 \(\overline P=\min(P^{\max},P_0^{\max},S_0^{\max})\)、\(\overline Q=\min(Q_0^{\max},S_0^{\max})\)、\(\overline\ell=\overline P/r\)（linear 为 0）；不能与原始输入 `capacity` 混称。


## 6. 两阶段规划问题与全局割

\[
D_{\mathcal B}=\{p\ge0:\exists x\in\{0,1\}^t,\exists y,
x\text{满足径向连通},\ c^\top x\le\mathcal B,\ (x,p,y)\in\mathcal F_{\rm SOCP}\}.
\]
不同 p 可以选择不同 x；这不是一套固定建设方案覆盖整个域。

| 符号 | 固定代码定义 | 单位 / 含义 |
|---|---|---|
| \(\omega\) | [PlanningModel.solve.weights](../model.py) | 非负支持方向，(d,) |
| \(q\) | [PlanningModel.solve.target](../model.py) | 待定界外域点，kW，(d,) |
| \(d(q)\) 的变量 \(d\) | [PlanningModel.distance_kw](../model.py) | 负荷退让量，kW；不属于 state，不是 eta |
| \(W\) | [PlanningModel.solve.cut_normals](../model.py) | (r,d)，每行非负、和为 1 |
| \(\rho\) | [PlanningModel.projected_deficit_kw](../model.py) | W 投影下的缺额，kW；不等于坐标距离 |
| \(\bar h,\underline d,\underline\rho\) | [PlanningModel.solve:bound](../model.py) | 单位和上下界方向由 mode 决定 |
| \(\varepsilon\) | [build_continuous_region.epsilon_kw](../main.py) | 全域认证误差，kW |
| 单次目标间隙 | [PlanningModel.solve.gap_kw](../model.py) | kW 目标的绝对间隙，不是 epsilon |
| \(b^{tree}_{ij},F^{(k)}_{ij}\) | [PlanningModel.__init__.parent_arc](../model.py)、[PlanningModel.__init__.commodity_flow](../model.py) | 连续树松弛辅助变量，不改变 x |

Stage 1：
\[
h(\omega)=\max_{x,p,y}\{\omega^\top p:(x,p,y)\in\mathcal F,\ c^\top x\le\mathcal B\},
\qquad \omega^\top p\le\bar h(\omega).
\]
最大化用 ObjBound 上界；可行候选目标不能作为割截距。初始盒每轴上界为 square_kw=network.power_limit，再以总负荷必要界和全局支持割裁剪。

Stage 2：
\[
d(q)=\min_{x,p,y,d\ge0}d,\quad p+d\mathbf1\ge q,\qquad
\rho_W(q)=\min_{x,p,y,\rho\ge0}\rho,\quad Wp+\rho\mathbf1\ge Wq.
\]
两问题保持完整规划约束。全局下界 L>0 给出
\[
\bigvee_k W_kp\le W_kq-L.
\]
这是 OR 析取，不能把各半空间同时加入模型。W=I 时为正交象限排除。
固定 x 的连续求解只生成热启动/可行点，其局部界不能生成跨方案割。
斜向 rho=0 不能证明 q 可行，此时主线改求正交 d。

## 7. 内外域、几何与认证

全部 points/vertices/threshold 使用 **kW**，d=2/3，函数不隐式换算。
下闭包条件：固定负荷和无功比例非负、r>0、reactance≥0，根电压为 1，电压上限允许回升至 1；不包含发电反送、最低负荷或其他非单调运行约束。入口检查这些数学前提。

\[
I_x=\operatorname{down}\operatorname{conv}\{p^j:x^j=x\},\quad
I=\bigcup_x I_x,\quad I\subseteq D_{\mathcal B}\subseteq O_2\subseteq O_1.
\]

| 数学量 / 操作 | 固定代码定义 | 说明 |
|---|---|---|
| \(I_x\) 面表示 | [downward_facets](../region.py)、[scheme_hulls](../region.py) | 每行 [F,g] 表示 Fp+g≤0；只按同一完整 x 分组 |
| \(O_s\) | [Cell.vertices](../region.py) | (n_v,d)，目标空间凸节点 |
| \(\bar\Delta_s\) | [Cell.gap](../region.py)、[cell_gap](../region.py) | min_x max_vertex dist_inf(vertex,I_x)，kW |
| 覆盖见证方案 | [Cell.scheme](../region.py) | 原始完整 x 的 tuple，不是新选型编码 |
| 距离 | [hull_distance](../region.py) | 对非负面用 F·max(q-d,0)+g≤0，处理坐标面 |
| 最大间隙节点 | [select_cell](../region.py) | 内域扩张使旧界仍有效，惰性更新保持全局最大优先 |
| 半空间裁剪 | [clip_polytope.constant](../region.py)、[clip_polytope.coefficient](../region.py) | constant + coefficient @ p ≥ 0 |
| 析取阈值 | [subtract_disjunction.threshold](../region.py) | Wq-L，kW；分区保留 OR，内部不交 |
| 全局传播 | [apply_disjunction](../region.py) | 割作用于所有相交节点 |
| 覆盖分支 | [cover_partition](../region.py) | 沿 \(I_x^\varepsilon\) 分割，不删除点 |

不能交换 min/max：顶点分别落在不同方案内域，不证明节点内部覆盖。
覆盖分支使用所有非空坐标掩码 S 表示距离邻域：
\[
(F_f)_S q\le-g_f+\varepsilon\sum_i(F_f)_{S,i}.
\]
覆盖部分若与父节点仿射维数不同，改试其他方案或最长坐标二分，避免重复原节点。
仅当 \(\max_s\bar\Delta_s\le\varepsilon\) 才返回 certified。
该证书针对 SOCP 投影；独立 AC 留在 vertify.py，端点 AC 通过不代表整个凸内域为 AC 可行。

## 8. 求解与结果契约

每次 [PlanningModel.solve](../model.py) 只调用一次 Gurobi optimize。
None 表示当前模型已证不可行；超时且无充分证据仍为 unknown。
power 固定时为 cost 模式，保留独立 AC 查询的最低投资语义。incumbent 只提供 Start，不增加成本约束。

| 字段 | 定义位置 | 语义 |
|---|---|---|
| x / p / state | [PlanningModel.solve:x](../model.py)、[PlanningModel.solve:p](../model.py)、[PlanningModel.solve:state](../model.py) | 型号向量、kW 负荷、运行向量；无候选时 None |
| mode | [PlanningModel.solve:mode](../model.py) | support / distance / projected_distance / cost |
| objective | [PlanningModel.solve:objective](../model.py) | 前三模式为 kW，cost 为费用单位 |
| bound | [PlanningModel.solve:bound](../model.py) | support 上界，其余下界；未知为 None |
| feasible / status | [PlanningModel.solve:feasible](../model.py)、[PlanningModel.solve:status](../model.py) | MaxVio≤PLANNING_TOL 接受运行证书；optimal/feasible/unknown |
| schema_version | [build_continuous_region:schema_version](../main.py) | 4；不自动加载或原地修改旧结果 |
| stage1_vertices | [build_continuous_region:stage1_vertices](../main.py) | O1，kW |
| inner / outer | [build_continuous_region:inner](../main.py)、[build_continuous_region:outer](../main.py) | I_x 列表 / O2 凸节点列表，vertices 均为 kW |
| inner.x/choice/cost | [build_continuous_region:x](../main.py)、[build_continuous_region:choice](../main.py)、[build_continuous_region:cost](../main.py) | 完整 x、走廊→型号、投资 |
| supports | [build_continuous_region:supports](../main.py) | weights/bound，Stage 1 全局半空间 |
| branch_cuts | [build_continuous_region:branch_cuts](../main.py) | target/cut_normals/threshold/bound，Stage 2 析取 |
| epsilon_kw/max_gap_kw | [build_continuous_region:epsilon_kw](../main.py)、[build_continuous_region:max_gap_kw](../main.py) | 目标误差 / 实际认证上界 |
| status/seconds | [build_continuous_region:status](../main.py)、[build_continuous_region:seconds](../main.py) | certified/empty/unknown/time_limit；构域墙钟秒，不含绘图 |

inner 保存 x/choice/cost/vertices，outer 只保存 vertices；不保存运行状态、事件和派生体积。
JSON 非有限界写 null，未知不能写零。state 是物理向量，旧扫描 states 是 -1/0/1 分类。

| 容差 / 配置 | 固定代码定义 | 单位 |
|---|---|---|
| 物理接受容差 | [PLANNING_TOL](../model.py) | 1e-8，所建模型最大违反量 |
| 全局界补偿 | [PAD_KW](../model.py) | 1e-4 kW；支持界加，距离界减；cost 保留原费用单位的求解器界 |
| 几何容差 | [GEOMETRY_TOL](../region.py) | 1e-8 kW |
| 全域精度 | [EPSILON_KW](../main.py) | kW |
| 单次/整体时限 | [QUERY_TIME_LIMIT](../main.py)、[CASE_TIME_LIMIT](../main.py) | 秒，两阶段共享整体时限 |
| 线程/预算/坐标 | [SOLVER_THREADS](../main.py)、[BUDGETS](../main.py)、[LOAD_NODES](../main.py) | 预算单位见 network，坐标用真实节点 ID |
| 独立 AC 容差 | [AC_TOL](../vertify.py)、[FIXED_POINT_TOL](../vertify.py)、[GLOBAL_AC_TOL](../vertify.py) | 沿用 v1 |

## 9. 显式迁移与退出接口

[v1 完整登记](notation-v1.md) 保留所有旧符号、接口和历史字段；现行绑定由 tests.test_notation 检查。

| 旧契约 | v2 处理及原因 |
|---|---|
| PlanningSP、eta、dual、cut=[alpha,beta,delta] | 退出主线；旧公式和乘子符号见历史第 6 节。新割来自整数优化全局界，不能把 eta 改名为 d。 |
| RemainingRegionModel、RegionState、gamma | 退出，由 Cell/cell_gap 和目标空间分支定界替代；两者不是同一个覆盖目标。 |
| cuts_only、min_total、cuts、add_cut | 退出；新主线始终带物理模型，按 support/distance/projected_distance/cost 区分目标。 |
| use_incumbent 成本截断 | 删除；Start 不改变可行集合，保留预算内更贵的方案。 |
| xi、tau、REGION_TAU、initial_polytope | 几何改用 kW；epsilon_kw 是绝对距离；main.py 直接建立初始盒，旧归一化包装退出。 |
| evaluation_bounds | 移至 tests.planning_checks 供历史实验使用；新主线直接求 SOCP 支持界。 |
| radial_gap_kw | 原 MP2 单目标间隙；新 gap_kw 适用于三个 kW 目标，旧关键字退出。 |
| BenchmarkResult、RunMonitor、live_view、counts/timing、max_total 等 | 旧比较/回放退出；新 schema=4；历史结果和源码快照不改。 |
| 旧 SP、剩余域、实时回放及基准专用测试 | 随退出接口删除；保留物理、潮流方向、拓扑和 AC 测试，增加新全局界/非凸认证回归。 |
| tests.case33_3d 的几何定义 | 同名量提升到 region.py 并推广 d=2/3；独立实验留作历史审核，不被主线导入。 |

新增量先登记。改单位、索引、状态布局、切片或返回字段，须同步调用方及数值测试。旧归一化说明不能再用来调用新函数。
