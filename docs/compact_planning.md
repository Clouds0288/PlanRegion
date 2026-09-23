# case33 逐线路紧凑规划模型与联合割

当前实现使用四、八、十六、三十二条嵌套候选线路，模型不建立完整建设组合表。配置来自 Network，主流程直接搜索逐线路型号；区域通过 MP1/MP2/SP 与全局残余搜索构造连续内外多面体，方案间取并集；网格仅用于独立 AC 校核。详见 [连续构域说明](continuous_region.md)。

## 1. 数据、变量和单位

四线路使用 A–D，八线路使用 A–H，十六线路使用 A–P，三十二线路使用 A–AF；支路端点和费用以 `Network/case33bw.py` 的项目表为准。候选走廊有两个型号：`existing` 保持原状，`parallel` 并联一回相同线路，等值 R/X 减半。统一配置为 `network.corridors`，包含全部 32 条在运走廊，而不只是升级候选。

### 统一走廊状态

`Corridor` 记录端点、原有型号 `existing_type`、初始投入状态 `initial_active`、使用要求 `must_use` 和型号参数 `types`。初始状态是数据，不是优化约束；原来有线或无线的走廊均使用同一规则。

对走廊 e、型号 k 定义二进制变量 $x_{e,k}$，并定义投入状态表达式：

$$x_{e,k}\in\{0,1\},\qquad z_e=\sum_k x_{e,k}.$$

代码使用双下标 `problem.x[走廊ID, 型号ID]`，通过 `problem.x.sum(走廊ID, '*')` 得到 $z_e$ 并直接添加约束，不管理分段下标。要求为

$$z_e=1\quad\text{必须使用},\qquad z_e\le1\quad\text{可使用或断开}.$$

所以投入时必须且只能选一个型号，未投入时不选任何型号。$z$ 是表达式，不额外创建状态变量。必须使用的单型号走廊直接取 $z_e=x_{e,k}=1$，其投资计入常数项；这种走廊不在 `problem.x` 中创建对应键。

`problem.x_vector` 将同一批变量按 `equations.variable_keys` 的顺序组织为向量视图，供矩阵方程、联合割和批量属性赋值使用，不创建第二套变量。`answer['x']` 是该顺序的数值向量；完整具名方案仍由 `equations.choice(answer['x'])` 得到。

每条走廊只保留一个参考方向，接根走廊参考方向向外，其余按端点顺序。可重构网络的非根走廊允许 P/Q 为负，不再建立方向二进制变量。FourBus 的五条走廊各有三种型号，共 15 个选型二进制变量；Case33 固定树按根向顺序，仍有 16 / 32 个变量对应 8 / 16 条升级候选。

对含 n 个非根节点的可重构网络，每条走廊只设一个有符号连续连通流 $f_e$。以 $A$ 表示去掉根节点行的关联矩阵（流入为正）：

$$Af=\mathbf 1,\qquad -nz_e\le f_e\le nz_e,\qquad \sum_e z_e=n.$$

任意不含根的连通分量都无法满足净流入等于其节点数，因此全部节点必须连通；再由边数为 n 可知所选图是树。该约束与实际负荷无关，零负荷也不能形成孤岛环。固定拓扑的 Case33 不需要这些辅助流。

对外方案仅为 `{走廊ID: 型号ID或None}`，`None` 表示未投入；方向由所选树自动确定。不使用整套方案枚举或另一套“新建候选”状态。

### Case33 固定拓扑

固定 32 条在运支路，5 条常开联络线保持断开。节点 18、25、33 的有功负荷组成参数 $p\in\mathbb R_+^3$（kW）；其余节点的负荷固定。定义完整节点负荷：

$$
d^P(p)=\frac{p^{\rm fixed}+Ep}{S_{\rm base}},\qquad
d^Q(p)=\frac{q^{\rm fixed}+E\operatorname{diag}(q/p)p}{S_{\rm base}}.
$$

潮流 P/Q、阻抗及电压使用标幺制；$v_j$ 为电压幅值平方，$\ell_{e,k}$ 为电流幅值平方。根节点 $v_0=1$。

对候选支路定义 $x_{e,k}\in\{0,1\}$，表示线路 e 选择型号 k：

$$
\sum_kx_{e,k}=1,\qquad \sum_{e,k}c_{e,k}x_{e,k}\le B.
$$

非候选支路只有一个固定型号，无需二进制变量。若有 $m$ 条候选线路，模型有 $2m$ 个选型二进制变量、$32+m$ 个支路—型号组合，连续物理变量为

$$y=(P_{e,k},Q_{e,k},\ell_{e,k},v_j),$$

共 $3(32+m)+32$ 个；四／八／十六／三十二候选分别为 140／152／176／224 个，选型二进制变量分别为 8／16／32／64 个。另有 3 个自由负荷变量。线性模型中的 $\ell$ 固定为零，由求解器预处理消除。

## 2. 型号与潮流的联接

给每个型号建立自己的 P/Q/电流变量，共用真实的节点电压。对于候选线路：

$$
0\le P_{e,k}\le\overline P_e x_{e,k},\quad
0\le Q_{e,k}\le\overline Q_e x_{e,k},\quad
0\le\ell_{e,k}\le\overline\ell_{e,k}x_{e,k}.
$$

未选型号的潮流严格为零。非候选线路的同类变量直接使用相应上界。

本算例所有负荷非负、线路阻抗为正，故上游源端功率不小于下游支路功率。取

$$
\overline P_e=\min(P_e^{\max},P_s^{\max},S_s^{\max}),\quad
\overline Q_e=\min(Q_s^{\max},S_s^{\max}),\quad
\overline\ell_{e,k}=\overline P_e/r_{e,k}.
$$

最后一个上界来自支路有功平衡 $P_e\ge r_{e,k}\ell_{e,k}$。未启用的限值按无穷处理。case33 的有限源端 P/Q 上界使这些界有限；原始线路 `rateA=0`，没有添加虚构热限。这些联接系数来自物理有效界，不是任取一个大数。

所有变量上界可统一写为

$$0\le y\le\overline y(x)=u_0+Ux.$$

其中候选型号上界乘对应的 x；固定支路和节点电压上界为常数。

### 可重构网络的带符号潮流

对于参考方向 $i\to j$，另一端向线路注入的功率为

$$P_{ji}=-P_{ij}+r\ell,\qquad Q_{ji}=-Q_{ij}+\chi\ell.$$

因此不能只限制 $|P_{ij}|$；必须同时约束两端注入：

$$P_{ij}\le\overline P x,\quad -P_{ij}+r\ell\le\overline P x,$$
$$Q_{ij}\le\overline Q x,\quad -Q_{ij}+\chi\ell\le\overline Q x.$$

反向容量约束已蕴含 $P_{ij}\ge-\overline P x$、$Q_{ij}\ge-\overline Q x$，直接模型不重复添加这些下界约束。SP 的变量盒仍保存它们，统一记作

$$\underline y(x)=\ell_0+Lx\le y\le u_0+Ux=\overline y(x).$$

代码分别用 `lb_const`、`lb_x`、`ub_const`、`ub_x` 表示 $\ell_0,L,u_0,U$；`y_lb_global`、`y_ub_global` 是不依赖选型的全局外包界。

节点平衡在参考受端扣除 $r\ell,\chi\ell$，压降和锥约束均使用参考送端电压。改变参考方向时同时变换上述功率和两端电压，物理可行域不变。未投入型号用正负压降余量断开两端电压联系，余量上界由电压范围推导并乘 $(1-x)$；未选型号的 P/Q/电流为零。

## 3. 紧凑 MILP／MISOCP 的潮流约束

令 e=(i,j) 为节点 j 的入边，$\mathcal C(j)$ 为 j 的出边集合，$\chi_{e,k}$ 表示电抗。对每个非根节点：

$$
\sum_kP_{e,k}=d_j^P(p)+\sum_{f\in\mathcal C(j)}\sum_kP_{f,k}
+\sum_kr_{e,k}\ell_{e,k},
$$

$$
\sum_kQ_{e,k}=d_j^Q(p)+\sum_{f\in\mathcal C(j)}\sum_kQ_{f,k}
+\sum_k\chi_{e,k}\ell_{e,k},
$$

$$
v_j=v_i-2\sum_k(r_{e,k}P_{e,k}+\chi_{e,k}Q_{e,k})
+\sum_k(r_{e,k}^2+\chi_{e,k}^2)\ell_{e,k}.
$$

保留原来的节点电压上下限、源端 P/Q 上限及三个独立负荷的共同总量外界。源端功率由接根支路的送端功率求和，包含损耗。

**MILP：** 令所有 $\ell_{e,k}=0$，不加入电流锥；其余方程均为线性，选型变量保持整数。

**MISOCP：** 保留 $\ell_{e,k}\ge0$，逐型号加入

$$
P_{e,k}^2+Q_{e,k}^2\le v_i\ell_{e,k},
$$

等价于

$$
\left\|(2P_{e,k},2Q_{e,k},v_i-\ell_{e,k})\right\|_2
\le v_i+\ell_{e,k}.
$$

$v_i\ge0$，所以这是凸旋转二阶锥。这里没有把 $r(x)P$、$r(x)\ell$ 等乘积直接交给非凸求解器；各型号 r、χ 均为常数。

整数 x 下只有一个型号的潮流非零，上述方程立即退化为该建设方案的原固定网架方程；反过来，把任意固定方案的可行潮流填入其选中型号，其他型号置零，就得到紧凑模型的可行解。因此两种表示对每个整数建设组合等价，预算内的规划域并集也等价。此结论来自双向构造，不依赖采样密度。

## 4. 两类规划查询

固定负荷点求最小投资，记总投资 $I(x)=C_0+c^Tx$：

$$C_M(p)=\min_{x,y}I(x),\qquad
\mathcal R_B^M=\{p:C_M(p)\le B\}.$$

也可只给总负荷下限 `min_total`，同时优化负荷分配：

$$\min_{x,y,p}I(x)\quad\text{s.t. }\mathbf 1^Tp\ge P_{\min},\ p\ge0,\ \text{模型 M 可行}.$$

未给定负荷点或总负荷下限时，在预算内最大化总负荷：

$$\max_{x,y,p}\mathbf 1^Tp\quad\text{s.t. }I(x)\le B,\ p\ge0,\ \text{模型 M 可行}.$$

三个独立节点的负荷自由分配，不限定方向或比例；固定背景负荷不参与优化。数值求解时，负荷目标用 $\mathbf 1^Tp/S_{\rm base}$ 改善尺度，返回值和界恢复为 kW。

## 5. 含规划变量的对偶割

将物理约束与变量盒统一写成

$$Ax+By+Cp\preceq_{\mathcal K}b,\qquad \ell_0+Lx\le y\le u_0+Ux.$$

定义锥余量

$$s(x,y,p)=b-Ax-By-Cp\in\mathcal K.$$

$\mathcal K$ 是非负锥与二阶锥的直积；线性行的 $\preceq_{\mathcal K}$ 就是逐元素 $\le$。$x$ 是选型变量，$y$ 是运行状态，$p$ 是 kW 负荷；$A,B,C$ 是对应的系数矩阵，$b$ 为右端常数向量。负荷 $p$ 在主问题中参与优化，在 SP 中固定，因此始终保留独立的 $Cp$ 项。

`PlanningEquations.A`、`.B`、`.C`、`.b` 与上述符号逐一对应。相对于早期余量形式 $c_{\rm old}+F_{\rm old}p+H_{\rm old}x+G_{\rm old}y$，有 $A=-H_{\rm old}$、$B=-G_{\rm old}$、$C=-F_{\rm old}$、$b=c_{\rm old}$；三个系数矩阵取负，物理约束和锥余量不变。`P_slice`、`Q_slice`、`ell_slice`、`v_slice`、`switch_slice` 记录 $y$ 的分块索引，`soc_slices` 记录二阶锥的行区间。

构造顺序为 `_index_network()`、`_build_variable_bounds()`、`_assemble_equations()`：先确定型号参数和变量列顺序，再建立变量盒，最后按物理公式逐块填系数并统一拼接。`type_corridor[k]` 表示型号条目 k 所属的走廊；`decision_types[j]` 表示决策变量 $x_j$ 对应的型号条目。

矩阵行依次为原始等式、等式负号副本、线性不等式和二阶锥。`eq_rows` 选取原始等式，`ineq_rows` 选取真实线性不等式；phase I 使用 `ineq_rows.stop` 之前的全部非负锥行。固定网架通过 `PlanningEquations(network.design(plan), method)` 建立，不再使用 `planning` 模式参数。

$A$ 包含反向送端限值的型号系数；固定拓扑下 $A=0,\underline y(x)=0$，退化为非负潮流模型。等式在 phase I 中表示为正负两行；直接 MILP/MISOCP 只建立一条等式。

固定整数选型 $\hat x$ 和负荷 $\hat p$ 后求

$$
\min_{y,\eta\ge0}\eta\quad
\text{s.t. }A\hat x+By+C\hat p-\eta r\preceq_{\mathcal K}b,
\quad \underline y(\hat x)\le y\le\overline y(\hat x).
$$

松弛方向 $r$ 对应 `relax_direction`，与支路电阻 $r_{e,k}$ 区分。线性行取 1；SOC 只在锥余量的首分量加 η。型号上下界保持严格。未选型号及线性模型的电流变量上下界均为零，直接消去；未选型号的零潮流锥也可省去，其对偶乘子回填为零。送入连续优化器的物理块为 $By-\eta r+s=b-A\hat x-C\hat p$，$s\in\mathcal K$。

取物理约束的对偶锥乘子 $\lambda\in\mathcal K^*$。对任意可行 $(x,p,y)$，有

$$\lambda^T(b-Ax-By-Cp)\ge0.$$

记 $h=-B^T\lambda$、$h_+=\max(h,0)$、$h_-=\min(h,0)$。由变量盒可得

$$-\lambda^TBy=h^Ty\le h_+^T\overline y(x)+h_-^T\underline y(x).$$

因此得到对所有合法建设组合有效的联合割：

$$
\boxed{\underbrace{\lambda^Tb+h_+^Tu_0+h_-^T\ell_0}_{\alpha}
+\underbrace{(-C^T\lambda)^T}_{\beta^T}p
+\underbrace{(-A^T\lambda+U^Th_++L^Th_-)^T}_{\delta^T}x\ge0.}
$$

割系数使用 $\alpha,\beta,\delta$，避免与方程的右端向量 $b$ 混用。代码中的 `cut` 仍按 `[alpha, beta..., delta...]` 排列。

这是盒约束支撑函数的精确上界，不要求浮点乘子的驻点残差严格等于零。代码将非负锥乘子截到非负域、SOC 首分量抬至尾部范数，保证对偶锥归属；常数再加 $10^{-12}$ 并作正比例归一化。只有可靠分离当前点的割才进入 MP。

特别注意：SP 消去了某些未选型号变量，**生成割时仍计算完整 B 的全部列，并使用 $h=-B^T\lambda$**。它们对 $h$ 和 $\delta$ 的贡献保留，使割能够在其他线路型号被选中时继续有效。若只使用当前选中型号的列，或遗漏带符号下界的贡献，就会失去这项跨方案保证。

## 6. 联合割循环与数值证书

1. 主问题仅包含逐线路 x、负荷 p、预算/查询约束及历史联合割，不含完整方案索引。
2. 求 MP 的候选及有效全局界，固定 x、p 调用连续 SP；最优性由可行目标与全局界的间隙决定。
3. 固定负荷查询不改变 p；负荷最大化查询检查比 MP 候选总负荷低 0.001 kW 的同方向内侧点。用 SP 电流重建功率平衡和压降等式，并检查原约束。内侧点通过后返回其负荷作为可行下界，同时保留 MP 的全局上界；两者的差就是停止间隙。
4. 否则添加 $\alpha+\beta^Tp+\delta^Tx\ge0$，继续求解；MP 不可行时，返回该查询不可行。

没有用人工迭代次数作为收敛条件。SP 即使报告数值停滞，也只能通过实际原始约束核验或有效分离割核验继续；既无可行证书又无分离割时返回未确定，不能把它当作认证。

原始标幺约束的联合 SP 容差为 $10^{-8}$，独立 AC 的容差不变。η 是辅助违反量，不是区域距离。最大总负荷查询使用 0.001 kW 的认证内缩，并与完整 MILP/MISOCP 在 0.002 kW 容差内核对。SOCP 连续域另采用径向精度 tau=0.002；几何容差在公共评价箱归一化坐标中为 1e-8。三维校核网格仅影响 FR/MR 统计。

## 7. 核对与运行入口

正式入口为 main.py；模型审核只在 tests 中执行，不进入区域求解主线。

- 指定基础、交错和全升级代表选型，用 `tests/reference.py` 的独立消元方程核对 LP/SOCP 物理边界。
- 对四、八、十六条候选线路，同一组负荷、预算和方向分别运行完整模型与联合割，核对最低投资及边界。
- 在包含全部整数选型和运行约束的紧凑模型上最小化联合割余量，以全局下界检查有效性；不通过完整方案表逐一核对。
- 用独立 AC 等式全局求解、节点导纳矩阵潮流及复相量重建检查 AC 参考实现。
- 对小网格逐中心求解，检查成块分类、预算共享和混合阶段重新认证。

切割有效性的理论依据是第 5 节的对偶锥与变量盒支撑函数。测试是实现核对，不用有限样本替代理论证明。

`main.py` 用完整模型查询 `planning_query` 和 `ContinuousRegion` 组织连续构域，`model.RemainingRegionModel` 求解全局残余问题，`region.RegionState` 维护多面体并筛选候选点。`vertify.ac_planning_query` 与 `vertify.validate_ac_region` 只用于最终 AC 校核，旧的多方法网格构域入口已移除。测试侧可使用 `tests/watchdog.py` 保护某次调用，正式流程不导入它。
