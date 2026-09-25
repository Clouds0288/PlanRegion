# 历史符号契约（冻结）

以下完整保留清理前的 1.3 契约，不作字段或符号改名。原实验入口已于 2026-09-25 退役；文中的旧路径指向项目外清理归档内的原相对路径，不属于活动代码接口。活动实现及明确迁移见 notation.md 第 18 节。

# 数学符号与代码变量规范

版本：1.3，2026-09-25。按用户要求恢复 d25c070 的 MP1 / MP2 / 对偶 SP 主线，并将同一几何流程推广至二维。三维接口、单位及联合割布局保持不变。

本文是本项目数学符号、代码名称、单位、数组顺序和结果字段的统一约定。正式代码为 `Network/`、`model.py`、`region.py`、`vertify.py`、`plot.py`、`main.py`；测试、注释和新文档使用同一约定。冻结的 `tests/jiangkou_2d/source/`、已归档结果和固定方案历史推导保留原口径，其局部记号须通过附录（原第 10 节）换算。

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
| \(c,\mathcal B\) | 型号投资向量、预算 | `network.cost`、`budget` | `network.cost_unit` |
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
| \((e,k)\) | Gurobi 建模键：走廊 ID 与该走廊内型号 ID；k 不要求全图唯一 | [PlanningEquations.keys](../model.py)、[PlanningEquations.types](../model.py) |
| \(\mathcal T_e\) | 走廊 e 的型号区间；型号到走廊的映射 | [Network.type_slices](../Network/__init__.py)、[Network.type_corridor](../Network/__init__.py) |
| \(\mathcal L,d\)；\(h\) | 独立负荷节点，`d = len(network.load_nodes)` | [Network.load_nodes](../Network/__init__.py)、[Network.selected](../Network/__init__.py) |
| \(\mathcal N_{\rm req}\) | 必须接入节点的布尔掩码，`(n,)` | [Network.required](../Network/__init__.py) |
| \(n_y=3t+n+2m\) | 标准运行向量长度 | [PlanningEquations.slack_slice](../model.py) 的 `.stop` |
| \(n_{\rm row}\) | 割生成 LP 的线性约束数，随支撑平面变化 | [PlanningSP._separating_cut.rows](../model.py) |
| \(n_{\rm tree}\) | 当前方案实际接入的非根节点数，也是树的支路数 | [OperatingTree.n](../Network/__init__.py) |

`network.type_keys[nu] == (e,k)` 是具名键与扁平位置的唯一桥梁。求解器内 `p[i]` 用真实负荷节点 ID；数组 `power[h]` 用 `load_nodes[h]` 的位置。`n/m/t/d` 在数学文档中固定为上表维数；`m = model`、`e = equations` 等既有局部对象别名不重定义这些数学符号。

模型装配可处理 `d` 个负荷坐标；正式 `region.py`、剩余域几何和网格输出当前固定 `d=3`。不得仅把文档写成任意维，就声称实现已经支持任意维。

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
| 完整运行对象、状态 | [PlanningModel.operation](../model.py)、[PlanningModel.state](../model.py) | `cuts_only=True` 时均为 `None` |

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

这些是数学块名，**当前实现没有 `equations.A/B/C/b`、`lb_x/ub_x` 等属性**；不得为了满足公式字面而重建一套平行模型。当前物理约束由 `add_operation` 的具名约束给出，割计算才从求解器导出线性行系数。

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

## 6. SP、对偶乘子和联合割

固定 \((\hat x,\hat p)\) 后，SP 最小化 \(1000\eta\)、\(\eta\ge0\)，当前 [PlanningSP.solve.eta](../model.py) 的求解器名称为 `violation`。仅功率平衡与压降等式允许 ±eta；选型容量、电压界、开断余量、电源限额和二阶锥保持严格。目标的正比例缩放不改变可行域或认证容差，`state` 不包含 eta。此辅助问题替代旧式“锥首分量也松弛”，不直接比较两者目标或乘子数值；eta=0 对应的物理约束、输出字段和数组布局不变。

SP 只先求上述带 eta 的问题。存在解且 `max(0, eta.X) + model.MaxVio <= PLANNING_TOL` 时，直接返回求解器的完整运行状态；两项共用误差预算，避免松弛量和求解误差分别达标、相加却超标。不重建状态，也不追加固定 eta=0 的损耗最小化。正 eta 的未完成求解不能单独证明不可行，仍须有效割；无解或无有效证书时保持未确定。

SP 设置 `NumericFocus=2`，提高边界与零潮流锥的数值稳定性；它作用于同一次求解，不改变 eta 的定义、物理约束、接受容差或取割代数。[Gurobi 数值参数说明](https://docs.gurobi.com/projects/optimizer/en/current/concepts/numericguide/numeric_parameters.html)

当前割仍由锥的有效支撑平面 LP 产生。令该 LP 在 eta=0 的线性行系数为 \((M_x,M_p,M_y)\)，右端为 \(b_{\rm row}\)。乘子方向采用 Gurobi 行约定：`>=` 行非负、`<=` 行非正、等式自由；固定 x/p 的等式乘子先置零。

| 数学量 | 代码定义 / 布局 | 说明 |
|---|---|---|
| \(\lambda\) | [PlanningSP._separating_cut.dual](../model.py) | `(n_row,)`，对应 `rows`，不是旧锥标准式的行序 |
| \(M^T\lambda\) | [PlanningSP._separating_cut.coefficients](../model.py) | 按求解器变量索引，不能按猜测切分 |
| \(h=M_y^T\lambda\) | [PlanningSP._separating_cut.h](../model.py) | `(n_y,)`，通过 `state_ids` 提取全部状态列 |
| \(h_+,h_-\) | `maximum(h,0)`、`minimum(h,0)` | 内联派生式，没有公开字段 |
| \(\alpha,\beta,\delta\) | [PlanningSP._separating_cut.cut](../model.py) | `cut[0]`、`cut[1:1+d]`、`cut[1+d:]` |

唯一割方向是

\[
\alpha+\beta^Tp+\delta^Tx\ge0,
\]
\[
\alpha=-\lambda^Tb_{\rm row}+h_+^Tu^{\rm glob}+h_-^Tl^{\rm glob},\quad
\beta=M_p^T\lambda,\quad \delta=M_x^T\lambda.
\]

`p` 仍为 kW，`beta` 必须作用于 kW 坐标。实现先作正比例归一化，再给截距加 `1e-10` 的保守补偿；不得把数值割系数直接解释成未经缩放的原始乘子。符号方向必须与上述行约定同时检查，不能搬用旧式 `-B.T @ dual`。`cut` 的布局、正负号、作用范围都属于契约。

[PlanningSP.solve](../model.py) 返回键 [PlanningSP.solve:feasible](../model.py)、[PlanningSP.solve:state](../model.py)、[PlanningSP.solve:cut](../model.py)。`feasible=True` 表示已获原始可行证书；`feasible=False, cut=None` 表示未确定，不能当作已证不可行。

## 7. 连续域、几何坐标和覆盖证书

| 数学量 | 固定映射 | 单位 / 形状 |
|---|---|---|
| \(b^{\rm box}\) | [RegionState.bounds](../region.py)，由 [evaluation_bounds](../model.py) 产生或由算例给定有效盒界 | kW，`(d,)`，各轴正上界 |
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
| `PlanningSP.solve(..., power)`、`PlanningModel(..., power=...)` | kW |
| `RegionState.add_point(..., point)`、`next_point`、`covering_schemes`、`witness_support` | \(\xi\)，无量纲 |
| `build_continuous_region` 中的 `point/pending` | \(\xi\)，直接调用 SP 前乘 `bounds` |
| `build_continuous_region` 中的 `seed/witness['p']` | kW；见证转成 `pending` 时除以 `bounds` 并乘 `1-tau` |
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

[PlanningModel.__init__.power](../model.py)、[PlanningModel.__init__.min_total](../model.py)、[PlanningModel.solve.radial_gap_kw](../model.py) 均用 kW；[PlanningModel.__init__.budget](../model.py) 用投资单位。没有固定负荷和最低总量时是 MP2（最大总负荷）；有任一条件时是 MP1（最低投资）。割作为 [PlanningModel.__init__.cuts](../model.py) 传入；[PlanningModel.solve.incumbent](../model.py) 与 [PlanningModel.solve.start](../model.py) 延续原有初始解语义，前者优先。

主流程直接创建模型并调用 `solve`，不保留 `planning_query`。MP 的运行证书和目标间隙在 `PlanningModel.solve` 中判断；数值质量不足且仍有候选时，由构域入口显式调用 SP。补救成功只补充 `state/feasible` 并标记 `status='feasible'`，不额外宣称 MP 最优。

| 字段 | 定义位置 | 固定语义 |
|---|---|---|
| `x`、`p`、`state` | [PlanningModel.solve:x](../model.py)、[PlanningModel.solve:p](../model.py)、[PlanningModel.solve:state](../model.py) | 型号向量、kW 负荷、完整运行向量；未找到候选时可为 `None` |
| `objective` | [PlanningModel.solve:objective](../model.py) | MP2 为可行总负荷 kW；MP1 为投资；须结合 `feasible` 解读 |
| `bound` | [PlanningModel.solve:bound](../model.py) | MP2 为全局上界，MP1 为全局下界；单位随目标变化 |
| `feasible`、`status` | [PlanningModel.solve:feasible](../model.py)、[PlanningModel.solve:status](../model.py) | `optimal/feasible/unknown`；由是否有运行解、求解器质量及目标界共同判断，不等同于独立 AC 认证 |

`None` 答案表示已证不可行。剩余模型的 [RemainingRegionModel.solve:bound](../model.py) 是 \(\overline\gamma\)，[RemainingRegionModel.solve:complete](../model.py) 是覆盖结论，不能按 MP1 / MP2 的目标值解释。

### 8.2 独立 AC 运行量

[ACPowerFlow](../vertify.py) 接收选定的 `OperatingTree`；其属性虽然叫 `network`，索引是树局部索引。`state(power, ell)` 返回 `(P,Q,v,u)`：均为 `(batch,n_tree)`，功率朝根向外；`v/u` 分别为受端 / 送端电压平方。内部 [ACPowerFlow._state.p](../vertify.py)、[ACPowerFlow._state.q](../vertify.py) 表示标幺节点负荷 \(d^P,d^Q\)，不作为 kW 接口传播。

`ell` 是电流平方，AC 残差为 \(P^2+Q^2-u\ell\)。[ACPowerFlow.classify](../vertify.py) 输出 `1/-1/0`（可行 / 已证不可行 / 未确定）；[ac_planning_query](../vertify.py) 返回经过 AC 校核的方案答案，并删除 SOCP 的 `'state'`，避免充作 AC 证书。

### 8.3 连续域结果和网格

| 字段 / 量 | 代码定义 | 单位 / 布局 |
|---|---|---|
| `max_total`、`max_total_bound` | [build_continuous_region:max_total](../main.py)、[build_continuous_region:max_total_bound](../main.py) | 最大总负荷的可行值、全局上界，kW；未获值可为 `None` |
| `max_point`、`max_choice`、`max_cost` | [build_continuous_region:max_point](../main.py)、[build_continuous_region:max_choice](../main.py)、[build_continuous_region:max_cost](../main.py) | kW 向量、具名方案、投资 |
| `inner`、`outer`、`vertices` | [RegionState.finish:inner](../region.py)、[RegionState.finish:outer](../region.py)、[RegionState.finish:vertices](../region.py) | 最终顶点为 kW；`inner` 各项另含 `choice/cost`；`outer` 为全局外包络 |
| `status`、`counts`、`timing` | [build_continuous_region:status](../main.py)、[build_continuous_region:counts](../main.py)、[build_continuous_region:timing](../main.py) | `certified/unknown/time_limit`；`sp/cuts` 次数；`total_seconds` 秒 |
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

网格估计不等于连续几何证明；AC 未确定点由 `disagreement_interval` 保留区间。FR/MR 的空分母返回 `None`，空并集的区域误差返回 0。不要把 `volume_gap` 比例值直接按 `*_percent` 显示。

## 9. 容差和运行参数

| 数学量 / 用途 | 唯一配置名 | 当前值 / 单位 |
|---|---|---|
| \(\varepsilon_{\rm plan}\)，求解结果接受容差 | [PLANNING_TOL](../model.py) | `1e-8`；Gurobi 所建模型的未缩放违反量，SP 另计 eta；尺度迁移见第 12 节 |
| \(\varepsilon_{\rm geom}\)，归一化覆盖 | [GEOMETRY_TOL](../region.py) | `1e-8`，无量纲；与物理容差独立 |
| AC 运行限值容差 | [AC_TOL](../vertify.py) | `1e-9`，标幺尺度 |
| AC 电流等式迭代容差 | [FIXED_POINT_TOL](../vertify.py) | `1e-12`，\(P^2+Q^2-u\ell\) 残差尺度 |
| AC 后备全局模型认证容差 | [GLOBAL_AC_TOL](../vertify.py) | `1e-7`，用于残差与运行违反量 |
| 径向收缩 \(\tau\) | [REGION_TAU](../main.py) | `0.002`，`0 <= tau < 1` |
| MP2 目标间隙容差 | [PlanningModel.solve.radial_gap_kw](../model.py) | `1e-3` kW；保留既有参数名，仅用于目标间隙判定，不再修改负荷，不是 `tau` |
| 单次 MP / SP / 剩余域 / AC 时限 | [MP_TIME_LIMIT](../model.py)、[SP_TIME_LIMIT](../model.py)、[RESIDUAL_TIME_LIMIT](../model.py)、[AC_TIME_LIMIT](../vertify.py) | 秒；SP 按方法取值 |
| 整体时限、线程 | [CASE_TIME_LIMIT](../main.py)、[SOLVER_THREADS](../main.py)、[DEFAULT_SOLVER_THREADS](../model.py) | 秒、正整数；混合阶段共享整体时限 |
| 每批 SP 候选数、AC 迭代上限 | [REFINEMENT_CHECKS](../region.py)、[AC_ITERATIONS](../vertify.py) | 正整数 |
| 算例 / 升级数 / 预算 / 剩余域模式 | [NETWORK](../main.py)、[UPGRADE_COUNT](../main.py)、[BUDGETS](../main.py)、[RESIDUAL_MODE](../main.py) | 配置，不另创造数学决策变量 |

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

- 删除生产接口 `planning_query`：其 `power/budget/min_total/fixed_plan/cuts/threads` 由 `PlanningModel` 建模接收，`incumbent/start/radial_gap_kw` 由 `PlanningModel.solve` 接收；参数名、单位、初始解优先级和目标间隙容差不变。`deadline` 由构域入口管理，传给各模型的是该次剩余 `time_limit`。
- 删除 `solve_region` 和 `ContinuousRegion` 调度层。`build_continuous_region` 直接调用 MP、SP、剩余域模型；`region.py` 只保留几何状态、选点、并集覆盖及裁剪。`RegionTimeout` 移到 `main.py`，不增加新的求解包装函数。
- `ContinuousRegion.tau` 迁到构域入口的 `tau` 与 `RegionState.tau`；原 `ContinuousRegion.finish` 的所有结果字段迁到 `build_continuous_region`，名称、单位和存储格式不变。`point/pending` 保持归一化坐标，`seed/witness['p']` 保持 kW。
- MP 可行点立即登记，防止 MP1 超时丢失 MP2 证书；按完整 `x` 合并初始方案，保留同方案的不同可行点。剩余域见证仍先补同方案支撑点，避免物理种子造成自覆盖。每批上限仍为 `REFINEMENT_CHECKS`，全局完成仍只由覆盖证书决定。
- hybrid 仍共享整体时限、只继承有效割，SOCP 阶段独立认证内域。冻结的历史源码、已有结果和 SP 的 eta / 对偶割代数保持原口径。

## 14. Case33 射线边界实验（独立测试入口）

以下新增量仅用于 `tests/benchmark_boundary_search.py`，不改变生产模型的目标和返回字段。

| 数学量 | 固定代码名 | 单位、形状、产生和消费位置 |
|---|---|---|
| \(\theta\) | [RayOracle.query.theta](../tests/benchmark_boundary_search.py) | 无量纲标量；射线 MP 的目标变量，`power = theta * direction` |
| \(\underline\theta,\overline\theta\) | [RayOracle.query.theta_lower](../tests/benchmark_boundary_search.py) / [RayOracle.query.theta_upper](../tests/benchmark_boundary_search.py) | 无量纲标量；射线查询返回，经可行证书 / 全局 ObjBound 更新，由外层扩域 / 裁域消费 |
| \(p^{\rm out}\) | [RayOracle.query.direction](../tests/benchmark_boundary_search.py) | kW `(3,)`；外顶点乘公共评价箱，传入射线查询；零分量不参与逻辑割 |
| \(\theta_{\mathcal I}(v)\) | [inner_scales.inner_scale](../tests/benchmark_boundary_search.py) | 无量纲；内盒沿外顶点射线的覆盖比例，由几何计算并用于选点 |
| \(\Delta(v),\max_v\Delta(v)\) | [coverage.coverage_gap](../tests/benchmark_boundary_search.py) / [boundary_search:max_coverage_gap](../tests/benchmark_boundary_search.py) | 无量纲；`1-inner_scale`，不与旧覆盖目标 gamma 混用 |
| \(V,\{\xi^{\rm feas}\}\) | [boundary_search.outer_vertices](../tests/benchmark_boundary_search.py) / [boundary_search.inner_points](../tests/benchmark_boundary_search.py) | 无量纲 `(n_points,3)`；外盒右上顶点 / 内盒证书点，分别用于全局覆盖和扩域 |

实验 JSON 使用独立 `schema_version`。`theta_lower=None` 表示尚无可行下界（没有原点证书时也不能填 0）；非有限全局界保存为 null，保留未知语义。既有 `x/p/state/cut` 布局和单位不变。射线目标直接读取 Gurobi theta 与 ObjBound，不调用按投资或 kW 解释目标的 `PlanningModel.solve`。

## 15. Case33 二维两阶段实验（独立测试入口）

`scheme_facets` 仅按完全相同的 `x` 分组，将该方案的认证点和向下闭包取凸包；其 `inner_facets` 为二维 kW 坐标中的 `(F_x,g_x)`，每行 `[F1,F2,g]` 满足 `F @ p + g <= 0`。`interval_gap` 对单个方案计算整个凸节点多边形到该凸内域的最大距离，再在方案间取最小值，得到有效覆盖上界；不交换这个 min/max 次序，也不跨方案凸化。树加强与这种内域加强均只在独立实验中启用。

实验保留原三维模型与 `x/p/state` 布局，仅固定未显示的负荷坐标。几何单独使用所选两轴的 kW 坐标，不修改正式三维几何。第二阶段为目标空间分支定界，节点代表负荷区间；节点定界用完整整数 SOCP，覆盖全部允许建设方案，不声称复用求解器内部的建设决策搜索树。

| 数学量 | 固定代码名 | 单位、形状与用途 |
|---|---|---|
| \(\omega\) | [SliceOracle.solve.weights](../tests/case33_two_stage.py) | 无量纲 `(2,)`，非负支持方向；目标为 `weights @ p[axes] / network.base` |
| \(\overline h(\omega)\) | [SliceOracle.solve:bound](../tests/case33_two_stage.py) | kW，完整规划问题的全局目标上界；不可用时保留未知 |
| \(p_1\ge t\) | [SliceOracle.solve.threshold](../tests/case33_two_stage.py) | kW，下界约束；用于二维最大第二坐标查询 |
| \(U\) | [SliceOracle.solve.upper](../tests/case33_two_stage.py) | kW；在当前 `threshold` 下已认证的第二坐标上界，反馈模型定界，不能跨不适用的区间复用 |
| \([a_s,b_s],U_s\) | [Interval.left](../tests/case33_two_stage.py)、[Interval.right](../tests/case33_two_stage.py)、[Interval.upper](../tests/case33_two_stage.py) | kW；区间与对该区间全部网架有效的第二坐标上界 |
| \(\overline\Delta_s\) | [interval_gap](../tests/case33_two_stage.py) | kW；节点外包盒到可行内盒并集的有向无穷范数距离上界，决定优先级与认证停止 |
| \(\varepsilon\) | [two_stage.epsilon_kw](../tests/case33_two_stage.py) | kW，整个外域到认证内域的距离上界容差；不是 SP eta 或物理残差容差 |
| \(\mathcal O_1,\mathcal O_2,\mathcal I\) | `stage1_outer`、`stage2_outer`、`inner_points` | 二维 kW 几何；前两者是认证外近似，后者是可行内盒右上角；只对同一方案允许凸插值 |
| \(h_g\) | [reference_scan.spacing_kw](../tests/case33_two_stage.py) | kW，独立逐点网格扫描间距；三态 `states` 为 `-1/0/1`，未确定不涂为不可行 |

支持割为 `weights @ p[axes] <= bound`。区间定界产生逻辑排除 `p1 >= threshold => p2 <= bound`，仅删除该条件范围中的区域。支持割和区间切除分别保存；候选目标值不能替代全局上界。绿色结果为第二阶段的外近似，须与认证内域及最大剩余距离一起报告；红色 `Actual region` 是注明网格精度的独立扫描参考，不表示解析真边界。

当前纯负荷径向 Case33 截面向下封闭，因此 `max p2, p1 >= threshold` 与 `max p2, p1 = threshold` 最优值相同；实验 `fixed_first=True` 使用后者定界，其他情况不启用该等价变换。所有面积字段单位为 kW²；`max_gap_kw`、`exact_outer_to_inner_gap_kw`、`outer_to_reference_gap_kw` 均为有向无穷范数距离，单位 kW。

二维实验的可选树松弛加强：\(b^{tree}_{ij}\) 对应 [add_tree_relaxation.parent_arc](../tests/case33_two_stage.py)，无量纲连续 `[0,1]` 根向父弧；\(F^{(k)}_{ij}\) 对应 [add_tree_relaxation.commodity_flow](../tests/case33_two_stage.py)，无量纲，每个非根节点的单位连通流。父弧两方向之和等于原走廊状态 `z`，非根节点入父弧之和为 1，每种流受父弧上界约束。这是原全接入径向树的扩展松弛，不新增建设决策或删除整数树；仅用于测试模型。对纯负荷树，正向支路有 `P >= r*ell, Q >= reactance*ell`，反向参考功率非正；用父弧门控的线性不等式加强这些必然关系。运行电压不高于根电压 1；前提为正阻抗、非负负荷、全部节点接入。

## 16. Case33 三维两阶段实验（独立测试入口）

`tests/case33_3d.py` 保留第 15 节的物理模型与树松弛加强，将全部负荷坐标释放。预算为相对投资单位，选型、负荷、运行状态沿用 `x/p/state`。所有几何坐标用 kW；不调用生产构域入口。新增结果 schema_version=3，不迁移历史文件。

| 数学量 | 固定代码名 | 单位、形状与用途 |
|---|---|---|
| \(\omega\) | [PlanningOracle3D.solve.weights](../tests/case33_3d.py) | 非负 `(3,)`，支持目标系数 |
| \(q\) | [PlanningOracle3D.solve.target](../tests/case33_3d.py) | kW `(3,)`，待定界的外域候选点 |
| \(d(q)\) | [PlanningOracle3D.distance_kw](../tests/case33_3d.py) | kW，最小负荷回退变量；约束 `power[i] + distance_kw >= target[i]`；物理运行约束始终不松弛，不是 SP 的 eta |
| \(\underline d(q)\) | [PlanningOracle3D.solve:bound](../tests/case33_3d.py) | distance 模式是全局下界，support 模式是全局上界；由 `mode` 区分，均为 kW |
| \(t=q-\underline d(q)\mathbf1\) | [subtract_orthant.threshold](../tests/case33_3d.py) | kW `(3,)`；保留析取 `p18 <= t1 OR p25 <= t2 OR p33 <= t3`，不能将其写成三个同时成立的线性割 |
| \(O_s\) | [Cell.vertices](../tests/case33_3d.py) | kW `(n_vertices,3)`，目标空间凸节点；节点之间内部不交 |
| \(\bar\Delta_s\) | [cell_gap](../tests/case33_3d.py) | kW，`min_x max_vertex dist_inf(vertex,I_x)`；完整节点距离上界，不能交换 min/max |
| \(\varepsilon\) | [construct.epsilon_kw](../tests/case33_3d.py) | kW，全部节点的最大认证距离停止容差 |
| \(\varepsilon_{warm}\) | [PlanningOracle3D.solve.stop_distance_kw](../tests/case33_3d.py) | kW，固定已知网架的连续求解已给出足够近可行点时跳过全整数求解；构域取 `epsilon_kw - 0.001`，不产生全局下界或排除割 |
| \(h_g,p_{33}^{slice}\) | [scan_slice.spacing_kw](../tests/case33_3d.py)、[scan_slice.fixed_value](../tests/case33_3d.py) | kW，独立二维网格间距与固定的第三坐标 |

`scheme_hulls` 仅在相同完整 x 内取已认证点的向下凸包。其证书属于 SOCP 投影；端点另经 AC 校核，不把端点 AC 通过当作所有凸插值的 AC 证明。距离使用 `max(q-epsilon,0)` 是否属于同方案向下凸内域计算，处理坐标面。阶段一支持割使用全局上界；阶段二回退问题使用全局下界（均作保守数值补偿）。因纯负荷域向下封闭，若所有坐标严格高于 t，则可以得到小于该下界的回退量，矛盾；因此被排除的开正交象限没有可行点，其边界保留。

如果节点各顶点分别被不同内域覆盖而整个节点尚未认证，继续做几何二分；不能据顶点逐个覆盖就判定整个节点覆盖。每次选择全体当前节点中最大 \(\bar\Delta_s\) 者。证书仅随内域扩张变强，已认证节点的子集继承证书。独立截面扫描不导入构域的割、方案或证书；对全部网格点给出 `-1/0/1` 分类并逐点 AC 回放可行点。截面红线不称完整三维 Actual region。

全局正交象限排除同步传播到所有相交节点。`global_solve=False` 表示仅获得固定方案的连续可行证书，此时 `bound=None`、`status=None`，不能按全整数最优或全局定界解读。第一阶段对已超时的支持方向避免法向分量差不超过 0.02 的近重复查询；被跳过方向不宣称达到 10 kW 凸包误差，最终精度独立由第二阶段的全域距离证书决定。

离线体积证书：令 `audit` 重算的最大距离为 \(e\)，则 \((\mathcal O_2-e\mathbf1)\cap\mathbb R_+^3\subseteq\mathcal I\subseteq\mathcal D\subseteq\mathcal O_2\)。同一平移保持节点内部不交，其体积和给出 [audit.inner_volume_lower](../tests/audit_case33_3d.py)（kW³）。所有已认证点向下凸包 \(K\subseteq\operatorname{conv}(\mathcal D)\)，`K` 中落在外域之外的体积给出 [audit.nonconvex_volume_lower](../tests/audit_case33_3d.py)（kW³）；该凸包只用于证明非凸性，不充作可行内域。

覆盖加速探针 `tests/probe_case33_3d_coverage.py` 使用 \(\gamma_{cover}\)（代码 `coverage_slack_kw`，kW）：在完整规划物理约束内，最大化对各同方案扩张内域的最小最大面违反量。面法向非负且和为 1；每个方案用一个析取选择违反面。只有全局上界非正才证明整个可行规划域被扩张内域并集覆盖；超时或可行候选点不能替代这个证书。`expansion_kw` 是同方案认证凸内域沿三个正坐标方向的等距扩张，单位 kW；此探针未通过前不参与正式构域。

三维审计中的 `outer_to_reference_gap_kw` 是固定同一个 \(p_{33}\) 后，整个绿色截面到独立网格红色参考域的有向无穷范数距离。它不等于允许三个坐标一起改变的三维距离证书。

斜向析取探针 `tests/probe_case33_3d_disjunction.py` 使用 \(W\)（`cut_normals`，每行非负、和为 1）和 \(\rho\)（`projected_deficit_kw`，kW），求完整规划模型下 `min rho, W @ p + rho >= W @ q, rho >= 0`。全局下界 \(\underline\rho\) 给出有效析取 `OR_k W[k] @ p <= W[k] @ q - rho_lower`。它是支持方向上的缺额，不等于负荷空间距离；内外域的停止判据仍由原始 `cell_gap` 给出。`cut_normals=I` 时退化为既有正交象限割。方向取自同方案内域仅是选方向的启发式，割的有效性由包含全部建设变量的全局下界保证。

算法修订 5 将该可选加速纳入 `PlanningOracle3D.solve.cut_normals` 和 `PlanningOracle3D.projected_deficit_kw`，事件模式为 `projected_distance`，其 `bound` 是 rho 的全局下界。既有 `distance_kw` 与 `distance` 模式语义不变。斜向割额外保存 `cut_normals`；`threshold` 对应每个面，长度等于法向行数。没有该字段的历史割仍解释为单位矩阵法向，不机械迁移历史快照。新切面 CSV 改为每面一行，新增 `facet,w18,w25,w33,threshold_kw,global_deficit_kw`；同一个第二阶段 index 的行是 OR 关系，JSON 是主要证据。

二元分区支持探针 `tests/probe_case33_3d_partition.py` 选原建设向量一个下标 `split_index`（不改变 x 索引约定），分别固定该二元量为 0、1。每个分区使用原支持目标并求全局上界，完整域包含在两个支持半空间的并集中。分区定界针对该分区内全部剩余建设变量，不能只求一个代表网架。两分区覆盖整个建设集合，未知上界时不能生成该析取证书。

覆盖分支 `cover_partition` 沿同方案内域的距离邻域 \(I_x^{\varepsilon}=\{q\ge0:\operatorname{dist}_\infty(q,I_x)\le\varepsilon\}\) 分割节点，不删除其中任何点。若原内域为 `a @ p <= b, p >= 0`，距离邻域的约束可写为所有坐标掩码 S 的 `a_S @ q <= b + epsilon_kw * sum(a_S)`，因为 `a @ max(q-epsilon_kw,0)` 是这些掩码表达式的最大值。`coverage_split` 仅是分支记录，不是不可行性割；覆盖子节点仍用原 `hull_distance` 重新计算全节点上界。斜向缺额 rho 等于 0 时不代表目标 q 可行，必须检查实际坐标距离或回退至正交查询。

算法修订 7 要求覆盖分支中已覆盖部分与父节点具有相同仿射维数；仅接触一个面、边或点时改试其他方案或几何二分。否则会重复切出零体积边界而留下完全相同的未覆盖节点。此进展条件不删除任何点、不改变割及停止证书，也不丢弃原有低维边界节点。

独立截面按递增 p18 查询，固定同一个 p33。此前已证明的 p25 上界 `scan_slice.upper`，以及当前最后一个未知网格点的 p25 加一个网格间距，都是当前及后续列的必要上界，可通过 `add_support((0,1,0),upper)` 加入截面查询模型。该条件依赖 p18 不减且 p33 固定，只在同一独立扫描内复用，不能当作完整三维规划域的全局支持割。新证据仍只来自该扫描自己的查询及必要界。

待定点列复核 `tests/recheck_case33_3d_columns.py` 的 `index` 沿用扫描的零基 p18 网格列下标；固定该 p18 与截面 p33，再最大化 p25。`upper` 只取同一独立扫描已证明的条件上界。固定网架副本仅生成 Start，最终模型的全部建设变量仍自由；最终支持上界严格小于待定点 p25 才能证其不可行。副本禁用目标界提前停止，避免停止条件干扰热启动求解。输出 `events/certificates/target/upper/seconds` 沿用既有含义；`started_at/finished_at` 为 UTC Unix 秒，仅用于核对并行复核时段，不重复累计重叠的墙钟时长。

`merge_rechecks` 顺序合并已完成复核的全局事件与可行证书，`column_recheck_files` 记录证据文件相对扫描目录的路径。只按重建的全局不可行性证明排除网格点；若仍有未知点，可由既有 `repair` 使用新可行证书继续检查。`seconds` 增加复核活动时段的并集长度，重叠运行时段只计一次。


## 10. 道路信息价值研究扩展（部分已实现）

本节登记 [道路信息价值模型](road_information_value.md) 使用的新量。下列英文名称是研究实现的建议映射；**第 10.1–10.5 节明确列出的量已在独立实验中实现，其余不是现有属性、形参或结果字段**。本节不改变 v2 主线 API、`schema_version=4` 或任何既有绑定。

| 符号 | 研究含义 / 建议映射 | 单位、索引与边界 |
|---|---|---|
| \(\mathcal G^{\rm road},\alpha,\mathcal E^{\rm road}\) | 候选道路图、道路段 ID、道路段集合；`road_graph` | 道路图不使用已登记的关联矩阵符号 G；道路段索引与走廊 e 分开 |
| \(\mathcal R_e\) | 走廊 e 必须经过的道路段集合；`corridor_road_segments` | 走廊 ID → 道路段 ID 集合；多条可选路线应在候选走廊中分别表示 |
| \(\theta_\alpha\) | 道路段真实可敷设状态；`road_usable` | 二进制随机量；勘察揭示它，不把 0 变成 1 |
| \(\mathcal I,\mathcal E_+^{\rm road}(\mathcal I),O_{\mathcal S}\) | 当前信息、已确认可用道路集合、勘察观测；`road_information`、`known_usable_roads`、`survey_observation` | 区分未调查、确认可用、确认不可用；带下标的 O 为观测，不是几何外域 O 或运行向量 y |
| \(\pi_\alpha(\mathcal I)\) | 后验可敷设概率；`road_usable_probability` | [0,1]；不使用已登记的负荷符号 p |
| \(b_e^{\rm road}(\mathcal I)\) | 走廊获准用于规划的掩码；`road_allowed` | 二进制，按 `network.corridors` 顺序 `(m,)`；新增约束 \(z_e\le b_e^{\rm road}\)，不压缩 x |
| \(\mathcal S,c_\alpha^{\rm survey},C_{\rm survey}(\mathcal S)\) | 勘察集合、单段费用、总费用；`survey_set`、`survey_cost`、`survey_total_cost` | 与 `network.cost_unit` 同单位；不覆盖建设投资 c |
| \(s,K,\mathcal B_{\rm total},\mathcal B_{\rm rem}\) | 实际勘察支出、勘察支出上限、总预算、剩余资金；`survey_spend`、`survey_budget`、`total_budget`、`remaining_budget` | 同费用单位；既有 `budget` / \(\mathcal B\) 始终表示输入构域模型的建设预算 |
| \(D_{\mathcal B}(\mathcal I)\) | 给定信息的获准规划域；研究数学量 | 沿用 p 的 kW 坐标；不同 p 仍可选不同 x；未登记第二套 region 输出 |
| \(\mathcal P_{\rm ref},\mu,R(\mathcal I,\mathcal B)\) | 固定负荷参照集、归一化需求测度、需求覆盖率；`demand_reference`、`demand_measure`、`demand_coverage` | \(\mu(\mathcal P_{\rm ref})=1\)，R 无量纲；几何体积仍为 kW\(^d\)，不能直接当货币 |
| \(\overline R(\mathcal S;\mathcal B),\mathrm{RIV}(\mathcal S;\mathcal B),\Delta_\alpha^{\rm info}\) | 预期覆盖率、道路信息价值、条件边际信息价值；`expected_demand_coverage`、`road_information_value`、`marginal_road_information_value` | 均无量纲；期望针对尚未揭示的道路状态；\(\Delta_\alpha^{\rm info}\) 不同于几何误差 \(\bar\Delta_s\) |
| \(\mathrm{RIE}_\alpha,\mathrm{NG}_\alpha\) | 信息投资效率、共用预算下单步净覆盖增益；`road_information_efficiency`、`net_coverage_gain` | 前者为覆盖率/费用，后者无量纲；二者不是货币收益率 |
| \(\mathrm{NG}_{\mathcal S},\mathrm{NGE}_{\mathcal S}\) | 一次选定勘察组合的净覆盖增益、净勘察效率；沿用 `net_coverage_gain`，新增研究映射 `net_survey_efficiency` | \(\mathrm{NG}_{\mathcal S}=\mathbb E R(\mathcal I_{\mathcal S},\mathcal B_{\rm total}-C_{\rm survey})-R(\mathcal I_0,\mathcal B_{\rm total})\)；费用为元时 \(\mathrm{NGE}=100\,\mathrm{NG}/(C_{\rm survey}/10^4)\)，单位为百分点/万元。不同于未扣建设预算机会损失的 RIE；零费用不求比值 |
| \(\mathcal C_{\min}(p,\mathcal I)\) | 固定负荷、固定道路掩码下的最低建设投资；现有 cost 模式 `objective/bound` 对应的最优值 | 费用单位；仅研究记号，不新增同义求解接口。已证不可行才记为 \(+\infty\)，未完成求解保留上下界与未知状态 |
| \(\Lambda,\mathcal J_{\rm survey},\mathcal C_{\rm eq}\) | 覆盖率价值换算尺度、净货币价值、等效增加建设预算；`coverage_value_scale`、`net_information_value`、`equivalent_construction_budget` | 费用单位；\(\Lambda\) 需外部标定，不能把面积/体积直接换成元；不占用现有矩阵 J |
| \(\beta_\alpha^{\rm road},\mathscr V\) | 道路状态枢纽重要度、最优自适应勘察后的覆盖价值；`road_pivotality`、`survey_policy_value` | 无量纲；前者衡量状态改变，不等同于信息价值；后者含未来勘察的延续价值 |

模型中的辅助函数、条件期望和内外界均为数学派生量，未引入现有类或结果字典的同义接口。Case33 的相对投资单位须先校准，才可与以元报价的勘察费用比较。

本轮聚焦的最小验证见 [道路价值实验计划](road_value_experiment_plan.md)：共用总预算，采用 NGE 衡量每万元勘察带来的净需求覆盖增益；R 与 NG 的既有定义保持不变。

### 10.1 最小实验的实现绑定

[main](../experiments/road_value/run_experiment.py) 实现本节的静态集合实验，不改动主线模型或 `schema_version=4`。未列在下表中的研究接口仍未实现。

| 数学量 / 记录 | 实现字段及定义 | 单位、顺序与边界 |
|---|---|---|
| 道路状态与获准掩码 | `known_usable_roads`、`road_allowed` | 本实验一条道路对应一个同名走廊；掩码按 01、13、02、23，长度 4；未勘察道路禁止采用 |
| 勘察集合与概率 | `survey_set`、`road_usable_probability`、`survey_observation` | 集合按用例的未知段顺序；概率独立，失败观察也支付全部已选勘察费 |
| 费用与预算 | `survey_cost`、`survey_total_cost`、`total_budget`、`remaining_budget` | 元；前者为每条道路报价，次者为集合总费；不勘察总费为 0 |
| 固定需求与覆盖率 | `demand_reference`、`demand_coverage`、`expected_demand_coverage` | 固定方格点集，两个坐标对应节点 1、3（kW），点等权；覆盖率均为 [0,1]，不是连续面积 |
| 结果增益与效率 | `net_coverage_gain`、`net_survey_efficiency` | 前者为无量纲覆盖差，后者为百分点/万元；零勘察费效率为 null |
| 最低建设费证据 | 逐点记录的 `p`、`x`、`state`、`objective`、`bound`、`feasible`、`status` | 沿用 cost 模式同名语义、形状和单位；图不连通另用 `proof=graph_disconnected` 证无解；不可行费用 JSON 为 null，不把未知设为无穷大 |
| 指标数值区间 | 上述覆盖率、增益、效率的 `_lower`、`_upper` 后缀字段 | 与对应指标同单位；未知点分别按 0/1 形成上下界；只有界一致时点估计字段非 null |
| 勘察正增益临界价 | `critical_survey_cost`、`critical_survey_total_cost` | 元/段、元/集合；有限需求集上净增益仍为正的报价上确界，不是连续规划域的精确阈值 |
| 需求网格密度 | `grid_intervals` | 每轴均匀分段数；10、20、40、80 分别对应 121、441、1681、6561 点，逐级复用重合需求点 |

实验保留原始求解上下界、数值违反量、运行时间、输入/源码摘要与断点日志；全局最优组合按预期覆盖率选取，NGE 只用于效率排序。

### 10.2 道路图结构筛查

独立入口 [main](../experiments/road_structure/identify_roads.py) 保留完整走廊 ID 与平行边。IEEE33bw 的电力走廊仅作道路代理；江口图含既有线、道路和入户走廊。原数据没有道路勘察状态，以下信息情景为明示的演示输入，不修改原始网架。

| 数学量 | 实现映射 | 单位 / 边界 |
|---|---|---|
| \(\mathcal K,\mathcal U,\mathcal F^{\rm blocked}\) | `known_usable_roads`、`unknown_roads`、`blocked_roads` | 按走廊 ID 的互斥集合；当前演示每条走廊视作一段道路，未推断平行走廊的独立物理风险 |
| \(h_{st}(\mathcal S\mid\mathcal K)\) | [connected](../experiments/road_structure/identify_roads.py) | 是否在 \(\mathcal K\cup\mathcal S\) 上存在 s–t 路径；测试时令 S 中道路均成功可用、其余未知道路不可用，经济评估另枚举观测；0/1；不等于含电气与预算约束的 R |
| \(\widehat{\mathcal G}_{\mathcal K}\) | [screen.quotient_edges](../experiments/road_structure/identify_roads.py) | 将已知图的连通分量收缩为节点，保留未知走廊 ID；分量内未知边为自环，对本轮纯连通判别无效 |
| \(J_{ef}(\mathcal K)\) | [truth_table.interaction](../experiments/road_structure/identify_roads.py) | \(h(ef)-h(e)-h(f)+h(\varnothing)\)；-1 为本目标/信息情景下的替代，+1 为互补，0 不由该指标分类 |
| \(L_e^{\rm cut}\) | `disconnected_load_kw` | 从完整候选图移去该边后，不再与电源连通的原始有功负荷之和（kW）；全图桥也可能为 0，不是勘察收益或预计停电量 |
| 桥与条件道路组 | `bridges`、`alternative_roads`、`complementary_bundles` | 桥为全候选图属性；其余为指定源、目标和信息情景的关系；同一道路在其他情景可有不同作用 |

筛查结果不替代 NGE。已知图已连接全部负荷时，单纯接入指标恒为 1；若分析绕行/备用路径，必须明确目标通道或被排除的参考边。当前长度为 2 的互补筛查是有界候选生成，不声称枚举了所有较长路径组合。

### 10.3 IEEE33bw 负荷加权生成树结构排名

独立入口 [main](../experiments/road_structure/rank_ieee.py) 对全部候选走廊的生成树做精确枚举。树均等权是人为声明的拓扑基准，不是实际运行概率、勘察成功概率或电气可行方案分布；本节不计算既有 R、NG、NGE，也不将生成树数量当作负荷规划域体积。

| 数学量 | 实现映射 | 单位 / 边界 |
|---|---|---|
| \(\mathscr T_{\rm top},N_{\rm top}\) | tree_masks、spanning_tree_count | 全节点无向生成树集合及数量；树掩码按输入走廊顺序，不改既有 x |
| \(n_{e,i}^{\rm path}\) | [enumerate_tree_paths.path_counts](../experiments/road_structure/rank_ieee.py) | 包含源到 i 路径上 e 的生成树数；数组 (m,n+1)，完整节点 ID 顺序显式保存，包含电源且其列恒为 0 |
| \(n_e^{\rm tree},q_e^{\rm tree}\) | tree_inclusion_count、tree_inclusion_fraction | 含 e 的树数及占总树数比例；等于删除 e 后丢失的拓扑比例；不等于随机道路可用概率 |
| \(L_e^{\tau}\) | initial_tree_load_kw 为原始树上的特例；其余由路径计数汇总 | 某棵树上经过 e 的原始负荷之和，未采用 e 时为 0；不是含损耗潮流或额定容量 |
| \(\overline L_e^{\rm top}\) | mean_tree_load_kw | \(\sum_i p_i^{\rm original}n_{e,i}^{\rm path}/N_{\rm top}\)，kW；允许路径方向随树变化 |
| \(I_e^{\rm top}\) | topology_importance_percent | \(100\overline L_e^{\rm top}/\sum_i p_i^{\rm original}\)，0–100 分；可解释为树与负荷加权目标节点联合抽样时该边在路径上的百分比，各边分数不要求和为 100 |
| 给定采用时的负荷 | conditional_tree_load_kw | mean_tree_load_kw / tree_inclusion_fraction，kW；仅作主指标的乘积分解 |
| 原始树的割与绕行 | initial_downstream_nodes、initial_bypass_roads | 原始树去边后的电源外分量、完整图跨该分量的其他走廊；常开联络边无原始下游，节点/道路列表为空 |
| 矩阵树定理校验 | [kirchhoff_tree_count](../experiments/road_structure/rank_ieee.py) | 整数 Bareiss 行列式；独立校验总树数及逐边删除后的树数 |

主排序使用原始负荷，另以各负荷节点等权、24/25 节点负荷翻倍作假设敏感性分析。新增字段只属于本实验，不替代主线规划结果字段或第 10.2 节全图删边负荷。

### 10.4 IEEE33bw 单条候选联络边增益比较

独立入口 [main](../experiments/road_structure/compare_candidate_metrics.py) 以 32 条初始闭合边构成的树为共同基准，每次仅加入一条初始断开的联络边；另外四条不加入。新增表示允许使用一条候选走廊，不表示将原本已建成的联络线认定为尚未建设。结果不替代第 10.3 节全图的结构排名，也不等于实际信息价值或电气规划域增量。

| 数学量 | 实现映射 | 单位 / 边界 |
|---|---|---|
| \(\mathscr T_0,\mathscr G_e=\mathscr T_0+e\) | `tree`、`augmented` | 简单无向图；基准必须为树；完整节点 ID（含电源）作为键 |
| \(r_e^{\Omega},w_e=1/r_e^{\Omega}\) | `r_ohm`、`conductance` | ohm / S；从原始 existing 型号及基值换算，不使用升级型号，也不改既有 p.u. 字段 r |
| \(D^{\rm path}_{\rho}=\sum_i p_i^{\rm original}d_{\rho}(s,i)\) | `path_distance` 中 `baseline`、`after`、`reduction` | 距离版本 `resistance` 为 kW·ohm；`hops` 为 kW·边数；后者不假定真实道路等长 |
| \(I_e^{\rm path,\rho}\) | `path_distance` 中 `gain_fraction` | (baseline-after)/baseline，无量纲；呈现百分数时乘 100；零基准值不定义，拒绝输入 |
| 逐节点路径证据 | `path_distance.node_evidence` | `node`、`original_p`、`before`、`after`、`reduction`、`weighted_reduction`；负荷 kW，距离单位跟随版本 |
| \(Q_f,\mathcal H,a_f\) | `cut_load_kw`；本次 H 为全部原有边，a_f=1 | Q_f 是原始树断边后的无源分量负荷，kW；不改变第 10.2 节全图删除负荷定义 |
| \(I_e^{\rm cut}\) | `cut_coverage_fraction`、`cut_numerator_kw`、`cut_denominator_kw` | 分子为基本回路上原有边的 Q_f 之和；分母为全部原有边 Q_f 之和；可重复计入同一负荷，不是新增可供电 kW |
| \(k_e,\tau(\mathscr G_e),I_e^{\rm tree}\) | `cycle_path_edges`、`tree_count_after`、`tree_gain` | 基准树数 1；k_e 为原树端点路径边数；树数 k_e+1；增益 k_e（新增棵数与相对增益在此数值相同） |
| \(R_{\rm eff}^{\mathscr T_0}(u,v), I_e^{\rm tree,w}\) | `effective_resistance_ohm`、`weighted_tree_gain` | R_eff 单位 ohm，树上等于端点间路径电阻和；加权树增益为 R_eff/r_e，不是实际树数量 |
| \(\lambda_2(\mathsf L)\) | `algebraic_connectivity` 中 `baseline`、`after`、`gain` | 未归一化图拉普拉斯、所有边权 1；无量纲；这里的图拉普拉斯不与既有标准式矩阵 L 混用 |
| 候选边介数 \(B_e\) | `edge_betweenness` 中 `raw`、`normalized` | 在各自 augmented 图上计算；边距离均为 1；无序节点对求和，normalized=raw/C(n,2)；不含负荷权重 |
| 排序与费用临界比 | `rankings`、`cost_ratio_thresholds` | 各指标分别降序，数值相同保留并列；对两个正增益 A/B，临界 c_A/c_B=I_A/I_B；未给报价，不输出货币收益或实际费用排序 |

仅比较基于同一参考图的增益；不把 37 条边先逐条删除后恢复与上述 5 条共同基准加边混排。输入限定简单树，平行边或已存在的候选边须另行建模，不能悄悄覆盖。

### 10.5 江口新建候选走廊与两段连接

独立入口 [main](../experiments/road_structure/jiangkou_candidate_metrics.py) 排除所有已有型号的走廊，只评价 `line_type is None` 的新建走廊。主基准为全部已建走廊（含可重构的常开线）的连通多重图，只保留这些走廊的端点；原始初始闭合树用于算法和基准敏感性。全部原始正负荷节点须在该核心图内。候选路网中未投入的零负荷节点不强制接入，平行走廊 ID 不合并。

| 数学量 / 记录 | 实现映射 | 单位与边界 |
|---|---|---|
| \(\mathscr G_0,V_0\) | `baseline`、`core_nodes`；[AugmentationBaseline](../experiments/road_structure/jiangkou_candidate_metrics.py) | 已建走廊及其节点；默认 343 条边、343 个节点；不将全部 594 个候选节点强制接入 |
| 参考树上的定向路径 | [AugmentationBaseline.path_vector](../experiments/road_structure/jiangkou_candidate_metrics.py) | 按初始闭合树边顺序的 -1/0/1 向量，用于拉普拉斯更新，不是电气潮流 |
| 负荷加权路径指标 | `path_distance` 的 `resistance`、`length`、`hops` | 沿用 baseline/after/reduction/gain_fraction；分别 kW·ohm、kW·m、kW·边数；仅给原始正负荷计权 |
| 广义原有割负荷 | `cut_load_kw`、`cut_numerator_kw`、`cut_denominator_kw`、`cut_coverage_fraction` | 对原有每条走廊 f，Q_f 是删 f 后新增断源的负荷；非桥 Q_f=0。只计候选连接恢复的原有割，避免重复奖励原有常开线已提供的旁路；分母为所有原有 Q_f 之和 |
| 生成树数量及增益 | `baseline_tree_count`、`tree_count_after`、`tree_count_increase`、`tree_gain` | 固定 V0；tree_count_increase 是新增棵数，tree_gain=(after-baseline)/baseline 为相对增益。IEEE 的基准为 1，二者恰同；江口须分开 |
| 有效电阻、加权树 | `effective_resistance_ohm`、`weighted_tree_gain` | 连通基准图上的有效电阻（ohm），加权相对增益为 R_eff/r_ohm；不将电阻距离与最短路径电阻混用 |
| 单条候选边适用性 | `attachment_count`、`evaluation_status` | 端点在 V0 中的个数为 0/1/2；仅 2 可直接评价同节点集的加边树数。0/1 的负荷路径、原有割增益确为 0，树数及两个全图对照为 null/不适用，不解读为道路无组合价值 |
| 统一候选型号 | `planning_line_type`、`r_ohm`、`construction_cost_yuan` | 主表统一 4x70；按模型长度、电阻率、单价与 new_corridor_cost_factor 计算，元；另对全部 3 档新建型号计算敏感性，不是最优型号或实测报价 |
| 单位建设费增益 | `path_efficiency_pp_per_10k`、`cut_efficiency_pp_per_10k`、`tree_efficiency_per_10k` | 前两者 1000000×gain_fraction/cost，百分点/万元；后者 10000×tree_count_increase/cost，棵/万元；零费用不定义，不是勘察收益率 |
| 两段连接 | `bundles`、`corridor_ids`、`intermediate_node` | 完整枚举两条新边经同一个 V0 外节点连接两个不同 V0 节点的组合；成本、电阻、长度相加；该中间节点原始负荷必须为 0 |
| 两段连接的核心图树数 | `core_tree_count_increase`、`core_tree_gain` | 仅将整条两段路径视作一个原子连接后、固定 V0 的树数比较；不是新增中间节点后的全节点生成树数。单条字段 tree_count_increase/tree_gain 置 null，禁止与原始全节点树数混排 |
| 代数连通度与新边介数 | `algebraic_connectivity`、`edge_betweenness` | 沿用第 10.4 节的单位边权和归一化方式，仅对两端在 V0 的单条新边计算；多重边分别保留 |
| 费用排序 / 原始基准敏感性 | `rankings`、`initial_tree_sensitivity` | 各指标独立降序并保留并列；后者只把初始闭合树作为基准，保持其余输入相同；不与主表混排 |

树计数采用初始树加少数弦的矩阵树行列式更新，并用独立枚举 / 拉普拉斯余子式验证。基准路径距离和割负荷固定，新增边只增加选择，不强制全部闭合。两段筛查不声称覆盖三段及更长互补链。没有输入勘察费用或道路可用概率，没有求解实际电气规划域。

### 10.6 江口候选道路交互展示分

[build](../experiments/road_structure/build_jiangkou_explorer.py) 只把第 10.5 节已经计算并验证的指标转换为可调的 0–100 展示分，不重新定义规划域价值，也不写回 `metrics.json`。综合分仅对 `attachment_count=2`、可单独连接既有核心网的候选走廊定义；其余候选保留“需组合评价”状态。

| 数学量 / 记录 | 实现映射 | 单位与边界 |
|---|---|---|
| \(s_e^{\rm path}\) | `path_score` | \(100 I_e^{\rm path}/\max_f I_f^{\rm path}\)，0–100；分母只取可单独评价候选 |
| \(s_e^{\rm cut}\) | `cut_score` | \(100 I_e^{\rm cut}/\max_f I_f^{\rm cut}\)，0–100；不是新增供电负荷比例 |
| \(s_e^{\rm tree}\) | `tree_score` | \(100\log(1+\Delta\tau_e)/\max_f\log(1+\Delta\tau_f)\)，0–100；对大树数作对数压缩 |
| \(s_e^{\rm cost}\) | `cost_score` | 路径、割覆盖、生成树三项单位建设费展示分的等权平均；前两项按各自最大效率归一化，生成树效率先取 \(\log(1+\eta)\) 再归一化；0–100 |
| \(S_e(w)\) | `composite_score` | \(\sum_j w_j s_e^j/\sum_j w_j\)，0–100；默认权重依次为 35、25、20、20，页面可调且自动按权重和归一化 |

`composite_score` 是多指标浏览与排序工具，不是电气可行规划域增量、道路信息价值、勘察收益率或货币净收益。原始子指标、费用与两段互补结果同时展示，以便追溯总分来源。

### 10.7 完整路线勘察与二维信息上下域

独立实验 `experiments/sequential_survey/` 固定建设预算，并将勘察费用单列。负荷坐标就是真实节点 1、2 的 p₁、p₂，单位 kW。既有线路仅保留原型号、无免费新增方案，故初始获准域与既有调度域相等；这不是一般升级模型的恒等式。

| 符号 | 实验记录 / 映射 | 单位及边界 |
|---|---|---|
| 完整路线 α、θ_α | `routes`、`road_usable_probability`、`survey_observation` | 一次勘察的原子对象。D 路线为 0–5–2，两段保留各自电气约束，同时建设或同时不建；中间节点 5 无负荷、无其他接入 |
| K⁺、K⁻ | `known_usable`、`known_unusable` | 已证可用及不可用的路线集合；未知真实状态只由观测模拟器读取，不进入评分 |
| D⁻_t、D⁺_t | `confirmed`、`optimistic` | 相同预算下分别只准 K⁺、或允许所有未证不可用路线；不同于数值几何的内外近似 |
| S₀、S⁻_t、S⁺_t | `baseline_area`、`confirmed_area`、`optimistic_area` | 连续二维面积，kW²；用户所称 Soptimal 对应会随负观测缩小的 S⁺_t，非已知真实最优域 |
| ΔS⁺_α、ΔS⁻_α | `gain_if_usable`、`removal_if_unusable` | 确认可用时 D⁻ 的新增面积；确认不可用时 D⁺ 的删除面积，均为非负面积 |
| π_α ΔS⁺_α、(1−π_α)ΔS⁻_α | `expected_expansion`、`expected_exclusion` | 预期能力扩展与预期排除，kW²；后者不是新增供电收益 |
| V_α、V_α/c_α^survey | `marginal_information_value`、`information_efficiency` | V 为上述两项之和，即期望信息上下域间隙减少；kW²、kW²/元；不同于仅计扩展的既有 Δ^info/RIE |
| G_t、τ_S | `information_gap`、`stopping_area_tolerance` | G=area(D⁺\D⁻)，kW²；停止须同时检查边际值及整个剩余间隙上界，不能假定贪心全局最优或边际值递减 |
| I_x、O_x | `inner`、`outer` | 固定完整建设方案的 SOCP 可行点下闭凸包与全局支持半空间交；方案间只作几何并集，不跨方案取凸包 |
| 数值面积区间 | `_lower`、`_upper` | 由固定方案并集的内外界传播；信息不确定性与数值近似误差分别记录 |

本实验采用有限候选全集、准确二元观测、独立合成先验和单一候选型号。原子路线合并不消除不同完整路线之间所有可能互补性。所有面积属于现有 SOCP 模型；有限 AC 核查不能宣称整个多边形已被 AC 认证。

### 10.8 改进 Das 15 节点勘察算例

`experiments/classic15_survey/` 保留 MATPOWER `case15da` 的 15 个节点、14 条原有支路阻抗、11 kV 电压等级、原始负荷及电压上下限。增加 A、B、C、F 四条完整候选走廊；其阻抗、建设费、勘察费、先验和观测轨迹都是研究设定，不属于原始测试系统。

- 沿用第 10.7 节的结果字段、集合、域、边际信息价值和数值区间语义；无接口改名。
- `load_nodes=(7,13)`，`power[0]`、`power[1]` 分别为节点 7、13 的总有功负荷；其他节点负荷保持原始值。文中二维轴 p₁、p₂ 是这两个坐标的简称，不表示物理母线 1、2；各节点原始 q/p 比保持不变。
- `power` 与全部计算结果继续使用 kW，面积继续使用 kW²。图形显示 MW、MW² 时分别除以 1000、10⁶，并明确标注节点映射；这不是代码单位迁移。
- `BUDGET=4`，`network.cost_unit='相对投资单位'`；勘察费为独立的相对勘察单位，不从建设预算扣除。新增费用不宣称为工程报价。
- 原始 `rateA=0` 表示未提供线路额定限值，转为 `capacity=inf`，不把它误作零容量；源端 P/Q 限值使用源文件值，基准 MVA 不当作额定电源容量。
- 所有原有节点必须供电，允许原有支路开断以保持径向；不额外限制可开断位置。没有新建走廊时只有原始树，因此初始确认域仍等于原网架调度域。
- 按原有 15 节点选择 14 条支路枚举预算内生成树；对所有方案分别求域并取并集。候选 F 的成本大于预算，不能出现在任何可建设方案中。
- `SUPPORT_BOUND_PADDING_KW` 是支持半空间上界的额外数值外扩量（kW），登记为 `support_bound_padding_kw`；只扩大数值外界，不扩大可行内域。固定方案面积精度为 50 kW²，停止容差为 1000 kW²。

### 10.9 五节点概念算例与逐轮评分历史

`experiments/concept5_survey/` 为独立设计的说明性算例，不标称 IEEE 或实测系统。沿用 10.7 的符号、字段、SOCP 模型及集合语义；建设预算和勘察费仍分开。节点 0 为电源，1、2 为二维可变负荷，3、4 为必接入中间节点。4 条既有线路加 A–F 六条完整候选路线；不是全网总共六条线路。所有阻抗均为正，电压与送端有功容量约束保留。

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

- 恢复来源为 `d25c070`；恢复前完整 v2 契约保存在 [v2 契约快照](notation-v2.md)，代码与未提交资料另有 `outputs/rollback_backups/20260925-171215/` 备份。原第 10 节局部记号表改为附录，保留道路研究第 10.1–10.9 节编号。
- `PlanningSP`、`eta`、`cut=[alpha,*beta,*delta]`、`RemainingRegionModel`、`RegionState` 重新成为活动接口。v2 的 `mode/epsilon_kw/branch_cuts/Cell` 等接口退出主线，历史结果不转换，不把旧 `tau` 和 kW 距离容差混用。
- 几何输入 `points/vertices` 为 `(N,d)`，`bounds` 为 `(d,)`，`d=2` 或 `3`。所有割仍为 `(1+d+t,)`，其中负荷系数始终乘 kW；`RegionState` 内部仍使用 `xi=p/bounds`。二维 `polytope_volume` 是面积，三维是体积，乘 `prod(bounds)` 恢复 kW^d。
- 原三维 AC 网格和实时查看器保留；二维五节点由专用勘察比较入口导出结果与图，不把二维数组伪装成三维网格。

| 数学量 / 记录 | 实现映射 | 单位、形状及边界 |
|---|---|---|
| 五节点候选案例 | [Concept5](../Network/concept5.py) | 完全沿用第 10.9 节的 4 条既有线与 A–F 六条候选路线，预算 4，p=(p1,p2)；不增加型号或改变参数 |
| 道路允许掩码 b_e^road | [Network.road_allowed](../Network/__init__.py)、[Concept5.__init__.allowed_roads](../Network/concept5.py) | 布尔 `(m,)`；默认全部允许；MP 加 z_e≤b_e^road，SP 不加道路掩码，因此物理联合割可跨信息状态复用 |
| 初始共享联合割 | [build_continuous_region.cuts](../main.py) | 与既有 `cuts` 相同的数组列表；传入后复制，不修改调用方；仅同一物理模型及有效盒界可复用 |
| 需求域缓存 | [DomainCache](../experiments/concept5_dual/run_experiment.py) | 按请求的道路允许集合求解，不预枚举建设方案；仅本次运行内缓存，记录各次冷启动及共享割耗时 |
| 比较面积误差 | `symmetric_difference_area`、`area_error_upper` | kW²；分别为新旧数值内域对称差、利用双方数值内外界得到的上界 |
| 面积停止精度 | `domain_area_tolerance` | kW²；数值面积区间须达标才能比较信息价值；不覆盖勘察停止容差 5 kW² |
| 参考计算 | `reference_directory`、`reference_seconds` | 指向冻结的 2026-09-25 枚举结果；不重新写入原目录；计时说明冷/热缓存及查询范围 |

`Network.fingerprint` 包含新增掩码，防止不同道路状态误用缓存。无模型参数变化时联合割可复用；预算/掩码仅进入 MP。勘察真值仍只由观测模拟器读取。数值 `outer` 与信息 `optimistic` 分别保留，不混同。

对偶域的 `plan_ids=null` 表示未枚举方案，不能凭两个 null 推断域相同；已有枚举实验仍保留原编号表及其等域快捷判定。新域插入缓存时，内域与已知子集内域取并，再与已知超集内域取交；外域与已知超集外域取交，再与已知子集外域取并。此前缓存不变，这些保守操作维持道路集合对应域的单调包含关系；处理后的数值面积区间仍须小于 0.5 kW²。
