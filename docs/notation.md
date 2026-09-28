# 数学符号与代码变量规范

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

SP 先求上述带 eta 的问题；求解必须正常结束，且 `MaxVio <= PLANNING_TOL`。`max(0, eta.X) + model.MaxVio <= PLANNING_TOL` 时直接返回原始运行状态；两项共用误差预算。正常求解得到正 eta 超限时，继续执行算法规定的支撑平面 LP 取割。超时、异常终止、证书精度不足或无法分离候选点均报错，不重建状态、不重求、不返回未确定结果继续运行。

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

[SubProblem.solve](../model.py) 返回键 [SubProblem.solve:feasible](../model.py)、[SubProblem.solve:state](../model.py)、[SubProblem.solve:cut](../model.py)。`feasible=True` 表示已获原始可行证书；否则必须返回有效割。旧版 `feasible=False, cut=None` 的未确定返回改为异常，见第 23 节。

## 7. 连续域、几何坐标和覆盖证书

| 数学量 | 固定映射 | 单位 / 形状 |
|---|---|---|
| \(b^{\rm box}\) | [RegionState.bounds](../region.py)，由总负荷上界构造或由算例给定有效盒界 | kW，`(d,)`，公共正数坐标尺度 |
| \(b^{\rm axis}(\mathcal B)\) | [RegionState.axis_bounds](../region.py)、[RemainingRegionModel.__init__.axis_bounds](../model.py)、[build_continuous_region:axis_bounds](../main.py) | kW，`(d,)`；方向 MP2 全局上界，限制本预算的搜索域 |
| \(w\)，方向目标 \(w^Tp\) | [MasterProblem.direction](../model.py)、[MasterProblem.__init__.direction](../model.py) | 无量纲，`(d,)`；默认全 1，轴向初始化取单位向量 |
| \(p_\Sigma^{\rm ub}\) | [RegionState.total_bound](../region.py) | kW；由输入总量上界和 MP2 全局上界收紧 |
| \(\tau\)、\(s_\tau=1-\tau\) | [RegionState.tau](../region.py)、[build_continuous_region.tau](../main.py) | 无量纲；linear 为 0 |
| \(I_x,O_x\) | [RegionState.add_scheme:inner](../region.py)、[RegionState.add_scheme:outer](../region.py)，存在 `records[tuple(x)]` | 归一化顶点 `(n_vertices,d)`，分别为认证内域、候选外域 |
| \((F_x,g_x)\) | [RegionState.inner_equations](../region.py) / [halfspaces](../region.py) 返回 `[F_x, g_x]` | `(n_faces,d+1)`；`F_x @ xi + g_x <= 0` |
| 裁剪常数、系数 | [clip_polytope.constant](../region.py)、[clip_polytope.coefficient](../region.py) | `constant + coefficient @ xi >= 0`，与半空间内侧约定相反 |
| \(\gamma\)，未覆盖量 | [RemainingRegionModel.__init__.delta](../model.py)，求解器名称 `uncovered_distance` | 内部乘 [RemainingRegionModel.distance_scale](../model.py)，返回前还原 |
| \(\overline\gamma\)，全局覆盖上界 | [build_continuous_region:coverage_bound](../main.py) | 无量纲；不是负荷上界 |

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
| `max_total`、`max_total_bound` | [build_continuous_region:max_total](../main.py)、[build_continuous_region:max_total_bound](../main.py) | 最大总负荷的可行值、全局上界，kW；未获值可为 `None` |
| `max_point`、`max_choice`、`max_cost`（2.1 退役） | 原 `build_continuous_region` 字段；迁移至 [RunMonitor.seed](../monitor.py)、[RunMonitor._geometry](../monitor.py)，见第 33 节 | 初始化 kW 点、具名方案与费用不再重复写入最终结果 |
| `inner`、`outer`、`vertices` | [RegionState.finish:inner](../region.py)、[RegionState.finish:outer](../region.py)、[RegionState.finish:vertices](../region.py) | 最终顶点为 kW；`inner` 各项另含 `choice/cost`；`outer` 为全局外包络 |
| `status`、`counts`、`timing` | [build_continuous_region:status](../main.py)、[build_continuous_region:counts](../main.py)、[build_continuous_region:timing](../main.py) | 成功返回 `certified`，未完成改为异常；旧 `unknown/time_limit` 记录保留；`sp/cuts` 次数；`total_seconds` 秒 |
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
| 全局搜索间隔、AC 迭代上限 | [REFINEMENT_CHECKS](../main.py)、[AC_ITERATIONS](../vertify.py) | 正整数；前者为自上次全局搜索以来的 SP 次数，活动见证处理完后生效 |
| 算例 / 升级数 / 预算 / 剩余域模式 | 原 `NETWORK/UPGRADE_COUNT/BUDGETS` 入口于 2.1 退役；当前 [main](../main.py)、[BUDGET](../main.py)、[RESIDUAL_MODE](../main.py)，见第 33 节 | 默认二维 FourBus、单预算；其他算例仍可单独使用其类 |

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
| 初始共享联合割 | [build_continuous_region.cuts](../main.py) | 与既有 `cuts` 相同的数组列表；传入后复制，不修改调用方；仅同一物理模型及有效盒界可复用 |
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
- [REFINEMENT_CHECKS](../main.py) 改为 32，计数仅含实际新增 SP 求解，缓存命中不计数；评分批次及活动见证完成后检查。生产算法没有 `eta_trigger`，没有 0.01 调度或停止阈值；终止仍须全局覆盖证书。
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
