# 统一候选网架、规划方程与联合割

符号、代码名、单位和数组顺序以 [数学符号与代码变量规范](notation.md) 为准。本文统一用 e 表示走廊、k 表示该走廊内的型号；具名模型直接使用 ID `(e,k)`，扁平数组位置另记 ν，通过 `network.type_keys` 转换。预算记为 \(\mathcal B\)，\(B\) 专用于数学标准式的运行矩阵。

所有算例使用同一候选图模型。规划选择走廊及其型号，形成覆盖必须供电节点的根树；初始开合状态不限制规划。Case33 包含 37 条走廊，FourBus 包含 5 条，江口默认读取完整候选图。

## 数据流与索引

`Network` 保存物理参数及唯一索引；`PlanningEquations.add_operation` 添加逐节点、逐走廊的具名 Gurobi 物理约束，MP 与 SP 共用该入口；`PlanningModel` 添加拓扑、预算及查询目标。选定的 `x` 直接生成 `OperatingTree`，供状态恢复和独立 AC 校验使用。

| 索引 | 对象 | 数据与变量 |
|---|---|---|
| \(i\) | 非根节点 | 负荷、required、电压限值、\(v_i\)、\(a_i\) |
| \(e\) | 候选走廊 | 端点、初始状态、\(z_e\)、\(f_e\)、\(s_e^+,s_e^-\) |
| \(k\in\mathcal T_e\) | 走廊 e 内的线路型号 | \(r_{e,k},\chi_{e,k},\overline P_{e,k},c_{e,k},x_{e,k},P_{e,k},Q_{e,k},\ell_{e,k}\) |

`network.type_keys[nu]` 保存走廊与型号 ID；`type_slices[e]` 保存走廊的型号区间；`type_corridor[nu]` 指向其走廊（这里 e 是走廊数组位置）。求解器用 `problem.choices[e,k]` 声明具名型号变量，`problem.x` 是同一批变量的扁平视图；全部型号均有对应元素，包括单型号走廊和固定方案中的型号。

方案字典仅用于输入、展示和保存：`network.encode_plan(plan)`、`network.decode_plan(x)`。内部算法直接传递完整 `x`。总投资为 \(c^Tx\)，没有固定型号投资常数或第二套压缩索引。

`network.tree(x)` 检查径向连通性，并返回选定节点、型号在 Network 中的索引，以及根向 parent、order、direction、D。原始候选图不保存另一套基础树参数。

## 选型与拓扑

令 \(\mathcal T_e\) 为走廊 e 的可选型号集合，则

\[
x_{e,k}\in\{0,1\},\qquad z_e=\sum_{k\in\mathcal T_e}x_{e,k},\qquad 0\le z_e\le1.
\]

\(z_e=0\) 表示走廊断开，\(z_e=1\) 时恰选一个型号。所有走廊使用同一规则。需要固定方案时，传入 `fixed_plan`，直接固定完整 \(x\) 的上下界。

以 \(n\) 为非根节点数、\(\delta^-(i),\delta^+(i)\) 为参考流入、流出节点 i 的走廊集合，\(a_i\) 为节点是否接入，根节点 \(a_0=1\)。建立

\[
a_i\in\{0,1\},\quad a_i=1\ (i\in\mathcal N_{\rm required}),
\]
\[
z_e\le a_i,\quad z_e\le a_j\quad(e=(i,j)),
\]
\[
\sum_{e\in\delta^-(i)}f_e-\sum_{e\in\delta^+(i)}f_e=a_i,\qquad
-nz_e\le f_e\le nz_e,\qquad \sum_e z_e=\sum_i a_i.
\]

连通流要求每个接入节点从根获得一单位虚拟流；边数等于接入非根节点数，故选定图为根树。约束不依赖实际负荷，零负荷的必须供电节点也不能形成孤岛。无负荷的可选中间节点可不接入。

Case33 和 FourBus 的非根节点全部必需。江口的建筑节点必需，非负荷中间节点可选；负荷坐标和带背景负荷的节点自动列为必需。

## 物理方程

负荷参数 \(p\) 用 kW，其余物理量用标幺值。\(v_i\) 为电压幅值平方、\(\ell_{e,k}\) 为电流幅值平方，根电压固定为 1。负荷为

\[
d^P(p)=(p^{\rm fixed}+Ep)/S_{\mathrm b},\qquad
d^Q(p)=(q^{\rm fixed}+E\operatorname{diag}(\kappa)p)/S_{\mathrm b}.
\]

其中 \(S_{\mathrm b}=\texttt{network.base}\)，\(\kappa=\texttt{network.q_ratio}\)。每条走廊采用固定参考方向，接根走廊朝外，其余遵循输入端点顺序。具名建模使用 `incoming/outgoing`，不需要 T/J。节点平衡直接写为

\[
\sum_{e\in\delta^-(i)}\sum_{k\in\mathcal T_e}(P_{e,k}-r_{e,k}\ell_{e,k})
-\sum_{e\in\delta^+(i)}\sum_{k\in\mathcal T_e}P_{e,k}=d_i^P(p),
\]
\[
\sum_{e\in\delta^-(i)}\sum_{k\in\mathcal T_e}(Q_{e,k}-\chi_{e,k}\ell_{e,k})
-\sum_{e\in\delta^+(i)}\sum_{k\in\mathcal T_e}Q_{e,k}=d_i^Q(p).
\]

非根走廊允许带符号 P/Q。两端向线路注入的有功分别为 \(P_{e,k}\) 和 \(-P_{e,k}+r_{e,k}\ell_{e,k}\)，因此约束

\[
-\overline P_{e,k}x_{e,k}\le P_{e,k}\le\overline P_{e,k}x_{e,k},\qquad
-P_{e,k}+r_{e,k}\ell_{e,k}\le\overline P_{e,k}x_{e,k},
\]
\[
-\overline Q_{e,k}x_{e,k}\le Q_{e,k}\le\overline Q_{e,k}x_{e,k},\qquad
-Q_{e,k}+\chi_{e,k}\ell_{e,k}\le\overline Q_{e,k}x_{e,k},
\quad 0\le\ell_{e,k}\le\overline\ell_{e,k}x_{e,k}.
\]

接根型号的 P/Q 下界为零。物理有效界取

\[
\overline P_{e,k}=\min(P_{e,k}^{\max},P_s^{\max},S_s^{\max}),\quad
\overline Q_{e,k}=\min(Q_s^{\max},S_s^{\max}),\quad
\overline\ell_{e,k}=\overline P_{e,k}/r_{e,k}.
\]

上述界适用于本项目的非负 P/Q 负荷、正电阻与非负电抗假设。未提供的容量记为无穷，源端必须提供足够的有限功率界。

每条走廊只写一条压降式：

\[
v_j-v_i+2\sum_{k\in\mathcal T_e}(r_{e,k}P_{e,k}+\chi_{e,k}Q_{e,k})
-\sum_{k\in\mathcal T_e}(r_{e,k}^2+\chi_{e,k}^2)\ell_{e,k}-s_e^++s_e^-=0,
\]
\[
0\le s_e^+,s_e^-\le M_e(1-z_e),\qquad
M_e=\max(\overline v_i-\underline v_j,\overline v_j-\underline v_i).
\]

接根时送端电压上下界均取 1。投入走廊的余量为零；断开走廊的型号潮流为零，余量允许其两端电压独立变化。未接入节点的电压仅为限值盒内的辅助量，不导出为运行结果。

线性模型固定 \(\ell=0\)。SOCP 逐型号增加

\[
(P_{e,k}^2+Q_{e,k}^2)\le v_i\ell_{e,k}
\iff (v_i+\ell_{e,k},2P_{e,k},2Q_{e,k},v_i-\ell_{e,k})\in{\rm SOC}.
\]

同时保留节点电压限值、源端有功/无功限值；SOCP 在提供容量时另有源端视在功率锥。源端功率含线路损耗。

整数选型下，每条投入走廊只有一个型号潮流非零。因此，走廊压降汇总与选定树的原物理方程等价；反向线路通过 \(P_{\rm reverse}=-P+r\ell\) 转换。改变分块方式可能改变分数选型下的连续松弛，不能要求新旧矩阵或分支搜索路径相同。

## 数学标准式与当前 SP

运行状态为

\[
y=(P_{e,k},Q_{e,k},\ell_{e,k},v_i,s_e^+,s_e^-).
\]

统一写为

\[
Ax+By+Cp\preceq_{\mathcal K}b,\qquad
l_0+Lx\le y\le u_0+Ux.
\]

上述 A/B/C、L/U 是整体数学表达；当前 `PlanningEquations` 不提供这些同名属性。Gurobi 建模使用 `P[e,k]`、`Q[e,k]`、`ell[e,k]`、`v[i]`、`plus[e]`、`minus[e]`；`P_slice/Q_slice/ell_slice/v_slice/slack_slice` 仅描述同一运行状态的扁平布局。

当前 SP 固定 x/p，以 `eta`（求解器名称 `violation`）度量功率平衡与压降等式的违反量：仅这些等式允许 ±eta，选型容量、电压界、开断余量、电源限额和锥约束保持严格。未选型号保持零潮流，Gurobi 可直接消去其退化锥。取零潮流、限值内电压和足够大的 eta 即可满足辅助问题，eta=0 恢复原物理约束。不能用旧实现“仿射盒严格、SOC 首分量松弛”的行顺序或乘子解释当前 SP。

SP 存在解且 `max(0, eta.X) + model.MaxVio <= PLANNING_TOL` 时直接返回求解器状态，松弛量和求解误差共用接受容差；不再沿树重建状态或固定 eta=0 后追加损耗最小化。其余情况仍用锥支撑平面构造 LP，读取线性行乘子生成联合割；正 eta 的中间解本身不充当不可行证书。eta 目标乘 1000 改善数值尺度，取割共用原 SP 时限；没有运行证书或有效割则返回未确定。Gurobi 二次约束违反量与旧 SOC 范数残差的尺度区别，以及测试检查器的迁移，见 [规范第 12 节](notation.md#12-11-迁移记录直接返回求解器结果)。

## 联合割的有效性

以下先给出整体仿射盒表达下的代数关系，不代表当前代码储存了 A/B/C 或仿射盒矩阵。

取 \(\lambda\in\mathcal K^*\)，令 \(h=-B^T\lambda\)、\(h_+=\max(h,0)\)、\(h_-=\min(h,0)\)。任意原始可行状态满足

\[
0\le\lambda^T(b-Ax-By-Cp)
\le\lambda^T(b-Ax-Cp)+h_+^T(u_0+Ux)+h_-^T(l_0+Lx).
\]

故对所有合法拓扑及型号成立

\[
\boxed{\alpha+\beta^Tp+\delta^Tx\ge0},
\]
\[
\alpha=\lambda^Tb+h_+^Tu_0+h_-^Tl_0,\quad
\beta=-C^T\lambda,\quad
\delta=-A^T\lambda+U^Th_++L^Th_-.
\]

该公式须使用全部状态列和完整仿射变量盒，不能只保留当前投入型号。当前 Gurobi 路径采用等价目的的另一组系数：线性行乘子按 `>=` 非负、`<=` 非正取向，固定参数等式乘子置零，令 \(h=M_y^T\lambda\)，使用全局运行盒消去 y：

\[
\alpha=-\lambda^Tb_{\rm row}+h_+^Tu^{\rm glob}+h_-^Tl^{\rm glob},\quad
\beta=M_p^T\lambda,\quad \delta=M_x^T\lambda.
\]

这里的 M 为支撑平面 LP 的线性行系数；与上述锥标准式的行序和乘子符号不同。两种推导都保持最终约定 `cut = [alpha, *beta, *delta]`、kW 负荷和 `>= 0` 方向。节点激活与连通流属于拓扑约束，不进入电气 SP；物理割对满足拓扑约束的子集仍有效。当前计算仅接受可靠分离当前点的割，具体代码映射见 [规范第 6 节](notation.md#6-sp对偶乘子和联合割)。

## 查询、验证与历史数据

MP2 在预算内最大化总负荷，MP1 在固定负荷或最低总负荷条件下最小化 \(c^Tx\)。主流程直接调用 `PlanningModel.solve`，以建模参数 `power/min_total` 区分目标；连续域流程见 [连续构域说明](continuous_region.md)。

固定原拓扑、交错升级和全升级的 Case33 方案使用独立消元 LP/SOCP 对照。FourBus 的 8 棵树、216 个型号组合可完整枚举核对。反向容量、零负荷孤岛、可选节点、全域割最小余量和独立节点 AC 潮流另有回归测试。

旧 Case33 的升级组合枚举及已保存结果仅对应固定初始拓扑。测试辅助函数 `fixed_topology` 显式截取该子图，仍使用同一规划模型；这些结果不能代表含 5 条联络线的完整重构域。

结果 schema 为 3，保存网架指纹、走廊/型号数量和升级数量。默认输出目录包含版本与指纹；重新计算入口拒绝把旧 schema 或其他网架的缓存用于当前配置。
