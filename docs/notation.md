# 数学符号与代码变量规范

版本：1.6，2026-09-25。勘察绘图合并至 plot.py，保留单一入口；MP1 / MP2 / 对偶 SP、单位、索引及联合割布局保持不变。

本文是本项目数学符号、代码名称、单位、数组顺序和结果字段的统一约定。正式代码为 `Network/`、`model.py`、`region.py`、`vertify.py`、`plot.py`、`main.py`、`survey.py`；测试、注释和新文档使用同一约定。项目外归档的源码、已有结果和固定方案历史推导保留原口径，其局部记号须通过附录（原第 10 节）换算。

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

模型装配可处理 `d` 个负荷坐标；正式 `region.py` 和剩余域几何支持 `d=2` 或 `3`，原 AC 网格输出和实时查看器仍固定三维。二维勘察由 `survey.py` 计算，`plot.py` 绘图。不得仅把文档写成任意维，就声称几何实现已经支持任意维。

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
- `PlanningSP`、`eta`、`cut=[alpha,*beta,*delta]`、`RemainingRegionModel`、`RegionState` 重新成为活动接口。v2 的 `mode/epsilon_kw/branch_cuts/Cell` 等接口退出主线，历史结果不转换，不把旧 `tau` 和 kW 距离容差混用。
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
