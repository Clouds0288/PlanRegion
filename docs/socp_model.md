# 线性域之后的 SOCP 校正模型与对偶割

**定位：固定方案 SOCP 数学背景。** 当前正式主线与实现入口见第 11 节；早期方案讨论不代表当前默认算法。

**记号作用域：** 本文保留固定树推导的局部符号；当前程序统一命名见 [数学符号与代码变量规范](notation.md)。本文 B 是预算、d/q 是标幺负荷、f 是功率因数，A/H/J 是固定树消元矩阵；不能直接当作当前候选图的同名矩阵、维数或虚拟流。对应关系见规范第 10 节。

本文针对当前四节点算例，推导“先用线性模型形成候选域，再用 SOCP 检查并切割”的数学模型。预算与建设方案沿用原问题；电气子问题补回线路损耗、完整平方电压降和配变实际视在功率约束。

需要审核的核心结论是：

1. 对固定建设方案，SOCP 子问题是凸问题，可以用对偶证书生成有效的线性可行性割。
2. 在当前模型的限定条件下，可以证明负荷投影域满足
   \(\mathcal D_{\mathrm{AC}}=\mathcal D_{\mathrm{SOC}}\subseteq\mathcal D_{\mathrm{LP}}\)。
   这不是对一般含分布式电源网络的结论。
3. 一条固定方案的 SOCP 割只约束该方案。不同方案的域取并集，不能直接共用该割或取凸包。

## 1. 建设方案、假设与标幺值

### 1.1 建设层

用 \(x\) 表示原模型的线路建设二进制变量，\(c^\top x\) 为建设费用。预算 \(B\) 下允许的方案集合为

\[
\mathcal X_B=\{x:\ x\text{ 满足原有连通树及线路选型约束},\ c^\top x\le B\}.
\tag{1}
\]

SOCP 子问题中固定 \(x=\bar x\)，因此拓扑、线路参数、容量均为常数。将所选树从根节点 0 向外定向，记为 \(\mathcal T_{\bar x}\)。非根节点数为 \(n\)，当前算例 \(n=3\)。

当前网络共有 216 个树与选型组合。预算 20000、40000、60000、无上限时，分别有 17、105、206、216 个允许方案。预算决定哪些方案参加并集，不改变某个固定方案的电气子问题。

### 1.2 本稿适用的电气条件

- 平衡三相、径向、无并联支路的支路潮流模型。
- 节点只有非负有功负荷，固定滞后功率因数；没有分布式电源、容性无功或反送电。
- 线路电阻 \(r_e>0\)、电抗 \(\chi_e\ge0\)，不含并联电纳、变比调节和相角差约束。
- 根节点电压固定为 1 p.u.，其他节点电压上限也是 1 p.u.。
- 线路约束保持现有数据中的送端有功上限；配变约束为实际送端视在功率上限。

这里的“AC 域”指与独立验证模型一致的完整、平衡径向 AC 支路潮流可行域。它保留损耗和非线性等式，但不代表额外考虑了三相不平衡等未建模因素。

### 1.3 单位与符号

设 \(p_j\) 为图中坐标，即节点有功负荷，单位 kW。采用三相容量、线电压基准：

\[
S_{\mathrm b}=150\ {\rm kVA},\qquad
V_{\mathrm b}=0.4\ {\rm kV},\qquad
Z_{\mathrm b}=\frac{0.4^2}{0.15}=1.0666667\ \Omega.
\tag{2}
\]

\[
d_j=\frac{p_j}{150},\qquad
q_j=\kappa d_j,\qquad
f=0.95,\qquad
\kappa=\tan(\arccos f)=0.3286841.
\tag{3}
\]

后文功率、阻抗、电流、电压均采用相容的标幺值。\(v_i=|V_i|^2\)，\(\ell_{ij}=|I_{ij}|^2\)。因此

\[
v_0=1,\qquad \underline v=0.93^2=0.8649,\qquad
\overline P_{ij}=\frac{\text{所选线路的容量（kW）}}{150}.
\tag{4}
\]

线路标幺阻抗由实际长度乘单位长度阻抗，再除以 \(Z_{\mathrm b}\) 得到。原数据没有独立的载流量参数，本稿不另行增加电流额定值。

## 2. 从完整 AC 支路潮流到 SOCP

固定负荷 \(d\ge0\) 和建设方案 \(\bar x\)。对支路 \(i\to j\)，令 \(P_{ij},Q_{ij}\) 为送端有功、无功，\(\mathcal C(j)\) 为节点 \(j\) 的子节点集合。

完整 AC 支路潮流方程为

\[
\begin{aligned}
P_{ij}&=d_j+\sum_{k\in\mathcal C(j)}P_{jk}+r_{ij}\ell_{ij},\\
Q_{ij}&=\kappa d_j+\sum_{k\in\mathcal C(j)}Q_{jk}+\chi_{ij}\ell_{ij},\\
v_j&=v_i-2(r_{ij}P_{ij}+\chi_{ij}Q_{ij})
             +(r_{ij}^2+\chi_{ij}^2)\ell_{ij},\\
P_{ij}^2+Q_{ij}^2&=v_i\ell_{ij}.
\end{aligned}
\tag{5}
\]

运行约束为

\[
\begin{aligned}
&\ell_{ij}\ge0,\qquad 0\le P_{ij}\le\overline P_{ij},\qquad Q_{ij}\ge0,\\
&\underline v\le v_j\le1,\\
&\left\|
\begin{pmatrix}
\sum_{j\in\mathcal C(0)}P_{0j}\\
\sum_{j\in\mathcal C(0)}Q_{0j}
\end{pmatrix}
\right\|_2\le1.
\end{aligned}
\tag{6}
\]

根节点可能连接多条支路，配变约束必须对所有根支路求和。它包含网络损耗，不能再以 \(\sum_jp_j\le142.5\) kW 替代。

式 (5) 最后一条是非凸等式。将其放松为

\[
P_{ij}^2+Q_{ij}^2\le v_i\ell_{ij}
\quad\Longleftrightarrow\quad
\left\|
\begin{pmatrix}
2P_{ij}\\2Q_{ij}\\v_i-\ell_{ij}
\end{pmatrix}
\right\|_2\le v_i+\ell_{ij}.
\tag{7}
\]

右侧是标准二阶锥约束。保留式 (5) 的前三条与式 (6)，用式 (7) 替换非凸等式，得到固定方案的 SOCP 可行性模型。它保留完整损耗项，不是再次线性化潮流。

定义固定方案的负荷投影域：

\[
\mathcal D_{\mathrm{SOC}}(\bar x)
=\{d\ge0:\ \exists(P,Q,v,\ell)\text{ 满足上述 SOCP}\}.
\tag{8}
\]

该集合是凸集，因为所有约束共同构成凸集，而线性投影保留凸性。预算下的规划域则为

\[
\mathcal D_{\mathrm{SOC}}(B)
=\bigcup_{\bar x\in\mathcal X_B}\mathcal D_{\mathrm{SOC}}(\bar x),
\tag{9}
\]

它一般不凸。AC 域和线性域采用相同的“固定方案投影后取并集”定义。

## 3. 精确消元：只保留电流平方

本节为了得到紧凑的子问题与清楚的对偶表达，消去线性等式中的 \(P,Q,v\)。这是精确代数消元，不是近似。

支路按其下游节点编号。定义下游矩阵

\[
D_{ej}=
\begin{cases}
1,&j\text{ 位于支路 }e\text{ 的下游子树内},\\
0,&\text{否则}.
\end{cases}
\tag{10}
\]

令 \(R=\operatorname{diag}(r)\)、\(X=\operatorname{diag}(\chi)\)，从叶节点向根累加功率平衡，得到

\[
P(d,\ell)=Dd+DR\ell,\qquad
Q(d,\ell)=\kappa Dd+DX\ell.
\tag{11}
\]

沿根到各节点的路径累加电压降，得到

\[
\begin{aligned}
v(d,\ell)&=\mathbf1-Ad-H\ell,\\
A&=2D^\top(R+\kappa X)D,\\
H&=2D^\top(RDR+XDX)-D^\top(R^2+X^2).
\end{aligned}
\tag{12}
\]

其中 \(D^\top\) 正是从支路量到各节点路径累加的矩阵。定义父节点选择矩阵 \(J\)：若支路 \(e\) 的上游节点是非根节点 \(j\)，则 \(J_{ej}=1\)；上游为根节点时该行全为 0。各支路上游电压平方为

\[
u(d,\ell)=\mathbf1-JAd-JH\ell.
\tag{13}
\]

根节点总送端功率可直接写成

\[
P_0=\mathbf1^\top d+r^\top\ell,\qquad
Q_0=\kappa\mathbf1^\top d+\chi^\top\ell.
\tag{14}
\]

将式 (11)–(14) 代入式 (6)–(7)，固定 \(d\) 后只剩 \(n\) 个电流平方变量。当前算例为 3 个；后面的可行性检查另增加一个 \(\eta\)。实际求解器可以保留辅助变量来显式表达锥，变量消元本身不构成运行时间保证。

式 (11)–(12) 与现有独立 AC 验证程序的电气方程一致。审核时应区分：方程一致是比较公平性的要求，SOCP 认证程序与 AC 非线性残差验证仍应分别执行。

## 4. 当前模型中 AC、SOCP 与线性域的关系

### 4.1 系数的单调性

在第 1 节假设下，\(A\) 和 \(H\) 均逐元素非负。\(A\ge0\) 直接成立。对 \(H\)，有

\[
H_{je}
=2\sum_aD_{aj}D_{ae}(r_ar_e+\chi_a\chi_e)
 -D_{ej}(r_e^2+\chi_e^2)\ge0.
\tag{15}
\]

若 \(D_{ej}=0\)，没有减项；若 \(D_{ej}=1\)，求和中 \(a=e\) 的一项为 \(2(r_e^2+\chi_e^2)\)，已足以抵消减项。

因此，固定 \(d\ge0\) 时，增加 \(\ell\) 会使 \(P,Q\) 逐元素不减，使 \(v,u\) 逐元素不增。

### 4.2 SOCP 域包含在线性域内

将损耗取为 0，可得现有线性模型对应的状态

\[
P^{\mathrm L}=Dd,\quad
Q^{\mathrm L}=\kappa Dd,\quad
v^{\mathrm L}=\mathbf1-Ad.
\tag{16}
\]

若 \((d,\bar\ell)\) 在 SOCP 中可行，则

\[
P^{\mathrm L}\le P(d,\bar\ell)\le\overline P,\qquad
v^{\mathrm L}\ge v(d,\bar\ell)\ge\underline v\mathbf1,\qquad
v^{\mathrm L}\le\mathbf1.
\tag{17}
\]

配变约束还给出

\[
1\ge\|(P_0,Q_0)\|_2
\ge\sqrt{1+\kappa^2}\,\mathbf1^\top d
=\frac{\mathbf1^\top d}{f}.
\tag{18}
\]

所以 \(\mathbf1^\top d\le f\)，即总负荷不超过 142.5 kW。固定方案的线性可行域可以直接写成

\[
\mathcal D_{\mathrm{LP}}(\bar x)=
\left\{d\ge0:
\mathbf1^\top d\le f,\quad Dd\le\overline P,\quad
Ad\le(1-\underline v)\mathbf1\right\}.
\tag{19}
\]

由此得到

\[
\mathcal D_{\mathrm{AC}}(\bar x)
\subseteq\mathcal D_{\mathrm{SOC}}(\bar x)
\subseteq\mathcal D_{\mathrm{LP}}(\bar x).
\tag{20}
\]

第一层包含关系来自放松等式；第二层依赖本稿的非负负荷与网络参数条件。这说明当前 LP 的有效割可以保留，线性候选域适合作为 SOCP 切割的初始外域。

### 4.3 特定条件下，SOCP 可行负荷可以恢复 AC 解

下面给出针对当前模型的证明，供审核。它说明的是负荷投影域相等，不是要求求解器返回的每一个 SOCP 状态都满足 AC 等式。

取任意 SOCP 可行状态 \(\bar\ell\)。在盒子 \(0\le\ell\le\bar\ell\) 中定义

\[
\mathcal F_e(\ell)
=\frac{P_e(d,\ell)^2+Q_e(d,\ell)^2}{u_e(d,\ell)}.
\tag{21}
\]

由第 4.1 节，盒子内 \(u_e(\ell)\ge u_e(\bar\ell)\ge\underline v>0\)，且 \(P,Q\ge0\)。因此 \(\mathcal F\) 连续、逐元素单调不减。SOCP 可行性给出 \(\mathcal F(\bar\ell)\le\bar\ell\)。

从零开始构造

\[
\ell^{(0)}=0,\qquad
\ell^{(k+1)}=\mathcal F(\ell^{(k)}).
\tag{22}
\]

利用单调性归纳可得

\[
0\le\ell^{(k)}\le\ell^{(k+1)}\le\bar\ell.
\tag{23}
\]

该序列逐分量单调有界，收敛到 \(\ell^\star\le\bar\ell\)。连续性进一步给出

\[
\mathcal F(\ell^\star)=\ell^\star,
\tag{24}
\]

即全部支路满足 AC 电流—功率等式。恢复状态满足

\[
P^\star\le P^{\mathrm s},\qquad Q^\star\le Q^{\mathrm s},\qquad
v^\star\ge v^{\mathrm s},\qquad v^\star\le\mathbf1.
\tag{25}
\]

这里 \(P^{\mathrm s},Q^{\mathrm s},v^{\mathrm s}\) 是电流平方为 \(\bar\ell\) 时的原 SOCP 可行状态。因此线路有功上限、最低电压、根节点视在功率上限都保持满足。最后一项上电压界由式 (12) 中 \(A,H,d,\ell^\star\ge0\) 保证。

树结构没有回路相角一致性条件。令 \(S_{ij}=P_{ij}+\mathrm iQ_{ij}\)，固定根节点电压相角后，可依次用 \(I_{ij}=(S_{ij}/V_i)^*\)、\(V_j=V_i-(r_{ij}+\mathrm i\chi_{ij})I_{ij}\) 恢复复电压和电流；式 (5) 保证其模长与上述状态一致。

故在第 1 节假设下，

\[
\boxed{\mathcal D_{\mathrm{AC}}(\bar x)
=\mathcal D_{\mathrm{SOC}}(\bar x)
\subseteq\mathcal D_{\mathrm{LP}}(\bar x).}
\tag{26}
\]

对同一组预算允许方案取并集，式 (26) 仍成立。

**适用边界。** 若引入反送电、容性注入、并联电纳，或把非根节点电压上限设得低于根电压，以上单调性或恢复后的上限保证可能失效，必须重新证明。Chen 论文讨论的更一般调度域存在 SOCP 松弛不精确的情形；不能把本节结论直接移植过去，也不能仅因某次 SOCP 解的锥约束不紧就删除对应负荷点。

现有 AC 全局模型中的辅助界 \(Q_e\le1\)、\(\ell_e\le1/\underline v\)，对恢复出的 AC 状态是冗余的：非负支路功率不大于根送端总功率，故 \(P_e^2+Q_e^2\le1\)，再用式 (24) 即可得到电流界。本稿省略它们不会扩大恢复后的 AC 负荷域。

## 5. 用 \(\eta\) 构造始终可解的 SOCP 检查问题

仅求“可行/不可行”不足以直接得到统一形式的割。仿照 Chen 的可行性松弛思路，对固定 \((\bar x,d)\) 定义

\[
\begin{aligned}
\phi_{\bar x}(d)=\min_{\ell,\eta}\quad&\eta\\
\text{s.t.}\quad
&\ell\ge0,\quad\eta\ge0,\\
&P(d,\ell)\le\overline P+\eta\mathbf1,\\
&(\underline v-\eta)\mathbf1\le v(d,\ell)\le(1+\eta)\mathbf1,\\
&\|(P_0(d,\ell),Q_0(d,\ell))\|_2\le1+\eta,\\
&\|(2P_e(d,\ell),\,2Q_e(d,\ell),\,u_e(d,\ell)-\ell_e)\|_2\\
&\hspace{38mm}\le u_e(d,\ell)+\ell_e+\eta,\qquad e=1,\ldots,n.
\end{aligned}
\tag{27}
\]

其中所有功率平衡和电压降等式已由式 (11)–(14) 精确满足。\(P,Q\ge0\) 随 \(d,\ell\ge0\) 自动成立。电流平方非负保持为硬约束；其他运行不等式和锥约束允许松弛。

各量已标幺化，本稿固定使用式 (27) 的松弛尺度。它与原 LP 中的 \(\eta\) 缩放不同，不直接比较两者数值。更换正的松弛权重不改变零值对应的可行域，但会改变得到的割。

\[
\boxed{\phi_{\bar x}(d)=0
\quad\Longleftrightarrow\quad
d\in\mathcal D_{\mathrm{SOC}}(\bar x).}
\tag{28}
\]

当 \(\eta>0\) 时，式 (27) 的状态只是用于寻找不可行证书，不能解释为满足运行条件的潮流。数值实现中应复核原始约束残差，不能把“求解完成”或一个未经检查的近零目标值直接当作认证。

## 6. 标准锥形式与强对偶

### 6.1 将全部约束堆叠

记标准二阶锥为 \(\mathcal Q_m=\{(t,z):\|z\|_2\le t,\ z\in\mathbb R^{m-1}\}\)。构造

\[
g_{\bar x}(d,\ell)=
\begin{pmatrix}
\ell\\
\overline P-P\\
v-\underline v\mathbf1\\
\mathbf1-v\\
(1,P_0,Q_0)\\
\big((u_e+\ell_e,\,2P_e,\,2Q_e,\,u_e-\ell_e)\big)_{e=1}^n
\end{pmatrix}
=c_{\bar x}+F_{\bar x}d+G_{\bar x}\ell.
\tag{29}
\]

其锥为

\[
\mathcal K=\mathbb R_+^{4n}\times\mathcal Q_3
\times\prod_{e=1}^n\mathcal Q_4.
\tag{30}
\]

各个二阶锥与非负正交锥均自对偶，故 \(\mathcal K^*=\mathcal K\)。为避免和建设费用混淆，\(c_{\bar x}\) 是本式的常数向量，不是式 (1) 的费用向量。

松弛方向为

\[
a_\eta=
\begin{pmatrix}
0_n\\\mathbf1_n\\\mathbf1_n\\\mathbf1_n\\
(1,0,0)\\
\big((1,0,0,0)\big)_{e=1}^n
\end{pmatrix}.
\tag{31}
\]

第一块为零，因为 \(\ell\ge0\) 不松弛。式 (27) 等价于

\[
\min_{\ell\in\mathbb R^n,\ \eta\ge0}
\{\eta:\ c_{\bar x}+F_{\bar x}d+G_{\bar x}\ell+\eta a_\eta\in\mathcal K\}.
\tag{32}
\]

这里把 \(\ell\) 视为自由变量，非负性已包含在第一个锥块中。

为使系数可逐项核查，令 \(A_u=JA\)、\(H_u=JH\)，\(e_e^\top\) 为第 \(e\) 个单位行向量。式 (29) 的系数为：

| 锥块 | 常数 \(c_{\bar x}\) | 负荷系数 \(F_{\bar x}\) | 电流系数 \(G_{\bar x}\) |
|---|---|---|---|
| \(\ell\) | \(0\) | \(0\) | \(I\) |
| \(\overline P-P\) | \(\overline P\) | \(-D\) | \(-DR\) |
| \(v-\underline v\mathbf1\) | \((1-\underline v)\mathbf1\) | \(-A\) | \(-H\) |
| \(\mathbf1-v\) | \(0\) | \(A\) | \(H\) |
| 配变锥 | \((1,0,0)\) | 行依次为 \(0,\mathbf1^\top,\kappa\mathbf1^\top\) | 行依次为 \(0,r^\top,\chi^\top\) |
| 支路 \(e\) 锥 | \((1,0,0,1)\) | 行依次为 \(-(A_u)_e,2D_e,2\kappa D_e,-(A_u)_e\) | 行依次为 \(e_e^\top-(H_u)_e,2(DR)_e,2(DX)_e,-e_e^\top-(H_u)_e\) |

### 6.2 强对偶成立的理由

对任何固定有限 \(d\ge0\)，先取 \(\ell=\epsilon\mathbf1>0\)，再取足够大的 \(\eta>0\)，即可使其他线性不等式及二阶锥均严格可行。因此式 (32) 满足 Slater 条件。注意 \(a_\eta\) 本身并不位于整个乘积锥的内部；严格可行性由这个具体构造保证。

目标下界为 0。又因每条线 \(r_e>0\)，在任何有界目标子水平集内，

\[
r_e\ell_e\le P_e(d,\ell)\le\overline P_e+\eta,
\tag{33}
\]

故 \(\ell\) 有界。可行子水平集闭且有界，原问题最优值能取到。因此式 (28) 不存在“下确界为 0 却没有零松弛可行解”的漏洞。

结合 Slater 条件，原、对偶最优值相等，且对偶最优值可取到。这里采用标准锥对偶结论，参见 Boyd 与 Vandenberghe，第 5.2.3、5.9 节。

## 7. 对偶问题及有效割

暂省略矩阵的 \(\bar x\) 下标。令 \(\lambda\in\mathcal K^*=\mathcal K\)，拉格朗日函数为

\[
\begin{aligned}
L(\ell,\eta;\lambda)
&=\eta-\lambda^\top(c+Fd+G\ell+\eta a_\eta)\\
&=-\lambda^\top(c+Fd)-\ell^\top G^\top\lambda
  +\eta(1-a_\eta^\top\lambda).
\end{aligned}
\tag{34}
\]

对自由变量 \(\ell\) 与非负变量 \(\eta\) 取下确界。有限下确界要求

\[
G^\top\lambda=0,\qquad a_\eta^\top\lambda\le1.
\tag{35}
\]

因此对偶问题为

\[
\boxed{
\begin{aligned}
\max_{\lambda}\quad&-\lambda^\top(c+Fd)\\
\text{s.t.}\quad&G^\top\lambda=0,\\
&a_\eta^\top\lambda\le1,\qquad\lambda\in\mathcal K.
\end{aligned}}
\tag{36}
\]

对偶可行域不依赖待检查负荷 \(d\)，所以最优值 \(\phi_{\bar x}(d)\) 是关于 \(d\) 的凸函数。

在检查点 \(\bar d\) 求得对偶最优解 \(\bar\lambda\)。若 \(\phi_{\bar x}(\bar d)>0\)，定义

\[
\beta_0=\bar\lambda^\top c,\qquad
\beta=F^\top\bar\lambda.
\tag{37}
\]

得到对该建设方案有效的线性割

\[
\boxed{\beta_0+\beta^\top d\ge0.}
\tag{38}
\]

**有效性证明。** 任意 SOCP 可行负荷 \(d\) 都有某个 \(\ell\) 使 \(g(d,\ell)\in\mathcal K\)。利用对偶锥性质与式 (35)，

\[
0\le\bar\lambda^\top g(d,\ell)
=\bar\lambda^\top(c+Fd)
=\beta_0+\beta^\top d.
\tag{39}
\]

而在当前不可行点，

\[
\beta_0+\beta^\top\bar d=-\phi_{\bar x}(\bar d)<0.
\tag{40}
\]

因此该割排除当前点，同时保留同一方案的全部 SOCP 可行负荷。依据式 (26)，也保留同一方案的全部 AC 可行负荷。

割有效只需要对偶可行；强对偶保证在不可行点能找到具有正对偶目标的分离证书。求解超时或数值状态不明不能替代这种证书。

图中使用 kW 坐标，所以实际几何割为

\[
\beta_0+\left(\frac{\beta}{150}\right)^\top p\ge0.
\tag{41}
\]

这里不能再额外把 \(\eta^\star\) 加入截距，也不能直接照搬原 LP 读取线性约束乘子的公式。锥对偶变量必须与式 (29) 的符号、排列和缩放一致。

## 8. 与原线性切割流程的衔接

### 8.1 初始域必须确实包含目标域

对每个允许方案，使用式 (19) 或原流程保存的、经有效 LP 割限制的候选外域，构造 \(U_{\bar x}^{(0)}\)。要求

\[
\mathcal D_{\mathrm{SOC}}(\bar x)\subseteq U_{\bar x}^{(0)}
\subseteq\mathcal D_{\mathrm{LP}}(\bar x).
\tag{42}
\]

这就是“用线性结果初始化”的几何含义。原 LP 的潮流状态一般不满足完整功率平衡，不能直接作为 SOCP 可行证书。

还需区分原结果中的候选外域和已认证点的凸包。原算法可能因为某部分被其他方案覆盖而停止扩展某一方案；此时该方案保存的认证凸包未必等于它完整的线性域，不能将其当作该方案的 SOCP 初始外域。式 (19) 提供了不依赖这一问题的明确初始化方式。

### 8.2 同时保存认证内域和候选外域

对每个方案维护

\[
I_{\bar x}^{(k)}
\subseteq\mathcal D_{\mathrm{AC}}(\bar x)
=\mathcal D_{\mathrm{SOC}}(\bar x)
\subseteq U_{\bar x}^{(k)}.
\tag{43}
\]

其中 \(I_{\bar x}^{(k)}\) 是该方案已经认证的负荷点的凸包，可以从零负荷点开始；\(U_{\bar x}^{(k)}\) 是初始线性外域与累计 SOCP 割的交集。

每次检查的逻辑是：

1. 从尚未认证的外域边界选取候选点 \(\bar d\)，固定方案 \(\bar x\)。
2. 若式 (27) 得到零松弛可行证书，将该点加入 \(I_{\bar x}^{(k)}\)。
3. 若得到正对偶目标证书，将式 (38) 加入 \(U_{\bar x}^{(k)}\)。
4. 更新边界并复用已检查点与已有割，继续检查尚未确定的部分。

因此，并非“每算一次 SP 就先认证成功，再切割”：可行点用于扩大认证内域；不可行点用于缩小候选外域。

同一方案的 SOCP 投影域是凸的，所以若某个外多面体的全部顶点均认证可行，则整个多面体都可行。这一结论不能用于混合不同建设方案的顶点。

预算下分别取并集：

\[
I_B^{(k)}=\bigcup_{\bar x\in\mathcal X_B}I_{\bar x}^{(k)},\qquad
U_B^{(k)}=\bigcup_{\bar x\in\mathcal X_B}U_{\bar x}^{(k)}.
\tag{44}
\]

所有预算可以共用同一方案的检查结果和割。这里优先复用的是几何结果与模型结构，不以求解器一定支持连续 SOCP 的数值热启动为前提。

### 8.3 有限计算下如何解释“更精确”

SOCP 域的边界可能是曲面，有限条直线割一般不能精确表示整个域。不能承诺有限次检查后完全恢复连续 AC 域；Chen 的顶点检查算法也区分了终止时的正确性与有限终止保证。

建议保留式 (43) 的两层结果：绘图和误差对比清楚标明使用认证内域还是候选外域，独立 AC 验证仍作为数值参照。

令 \(\mu\) 表示三维体积，\(R\) 为所报告区域，\(A=\mathcal D_{\mathrm{AC}}(B)\)。Chen 的 FR、MR 用采样点比例定义；在均匀采样口径下，对应以下体积比例：

\[
\mathrm{FR}(R)=\frac{\mu(R\setminus A)}{\mu(R)},\qquad
\mathrm{MR}(R)=\frac{\mu(A\setminus R)}{\mu(A)}.
\tag{45}
\]

实际计算仍用同一验证网格估计这些比例，显示为百分数时乘以 100%。FR 对应图中的黄色多余部分；MR 对应红色遗漏部分。原区域不一致率保留为

\[
E(R)=\frac{\mu(R\triangle A)}{\mu(R\cup A)}.
\tag{46}
\]

在精确算术和本稿假设下，输出 \(I_B^{(k)}\) 可保证 FR 为 0，输出 \(U_B^{(k)}\) 可保证 MR 为 0；有限计算时不保证两者同时为 0。

内外域间的体积差可作为内部停止依据，不另增加对外报告指标。若 \(\mu(U_B^{(k)})>0\)，设

\[
\varepsilon_k=
\frac{\mu(U_B^{(k)}\setminus I_B^{(k)})}{\mu(U_B^{(k)})},
\tag{47}
\]

则由集合包含关系，分别有 \(\mathrm{MR}(I_B^{(k)})\le\varepsilon_k\) 和 \(\mathrm{FR}(U_B^{(k)})\le\varepsilon_k\)，前提是相应分母非零。这里必须计算并集体积，不能直接累加重叠方案的体积；若用网格估算体积，该停止判据及误差界也只是数值估计。

结合“减少多余区域”的目标，建议最终将认证内域作为已认证结果，同时保留外域展示尚未确定的边界厚度。模型本身不预设一个未经验证的精度阈值。

## 9. 建设变量的处理：不能把固定方案割当作全局割

式 (29) 的矩阵随拓扑和导线参数变化。上述推导没有证明它们对建设二进制变量呈原 LP 那样的统一仿射参数关系，因此式 (38) 不能无条件加入所有方案。

**建议首版按方案保存割。** 这与当前算例只有 216 个方案的规模相符，也使割的适用范围清楚。减少实际 SOCP 调用依靠候选边界检查和结果复用，不需要对全部验证网格点求 SOCP。

若后续仍要放回统一建设主问题，可以用条件割。设

\[
\Delta(x,\bar x)
=\sum_{\bar x_i=0}x_i+\sum_{\bar x_i=1}(1-x_i).
\tag{48}
\]

对于二进制建设变量，\(\Delta=0\) 当且仅当 \(x=\bar x\)。在共同负荷外界 \(\Omega=\{d\ge0:\mathbf1^\top d\le f\}\) 上，

\[
\min_{d\in\Omega}(\beta_0+\beta^\top d)
=\beta_0+f\min(0,\min_j\beta_j).
\tag{49}
\]

因此取

\[
M_{\bar x}=
\max\{0,-\beta_0-f\min(0,\min_j\beta_j)\},
\tag{50}
\]

可写出有效条件割

\[
\beta_0+\beta^\top d+M_{\bar x}\Delta(x,\bar x)\ge0.
\tag{51}
\]

它只在 \(x=\bar x\) 时施加原割；对其他方案在整个 \(\Omega\) 上自动成立。这是方案条件化，并不意味着得到了一条能同时加强其他方案的电气割。按方案维护几何区域时无需引入 \(M_{\bar x}\)。

## 10. 审核重点与后续实现边界

本稿建议审核以下四项：

1. **物理口径：** 保留式 (5) 的损耗和完整电压降，线路仍按已有 kW 容量约束，配变采用含损耗的 kVA 约束。
2. **认证依据：** 接受式 (26) 在当前假设下的恢复证明；仍保留独立 AC 残差与区域验证，用于检查实现和数值误差。
3. **割的范围：** 首版固定方案求 SOCP、按方案保存对偶割，最终对所有预算允许方案取并集。
4. **输出含义：** 有限计算保留内、外两层区域；认证内域作为认证结果，外域用于展示剩余不确定部分，继续记录 E、FR、MR。

上述保留最初模型的审核边界。当前统一规划模型与模块职责见下节。

已做的代数核对：对当前 216 个建设方案，式 (11)–(14) 与独立 AC 程序的状态计算一致；数值检查中最大差为 \(1.11\times10^{-16}\)，且 \(A,H\) 逐元素非负。该核对辅助排查推导与代码的记号错误，不替代第 4 节证明，也不是 SOCP 性能实验。

## 11. 与当前主线的关系

本文保留固定建设方案下的 SOCP 方程、对偶割与适用条件，作为数学背景。当前主线已使用统一的规划变量与运行约束，不再采用早期 Notebook 中逐方案枚举的实现；统一联合割的依据见 [紧凑规划模型](compact_planning.md)。

- `main.py` 组织完整 MP2 / MP1、必要 SP 与连续覆盖循环。
- `model.py` 实现统一 LP / SOCP 物理方程、主问题、SP 和轻量 / 物理剩余域搜索。
- `region.py` 管理每个方案的认证内域、候选外域以及跨方案并集。
- `vertify.py` 独立进行 AC 校核；`plot.py` 与 `live_view.html` 保存和回放全部过程。

当前区域误差以**实际保留的认证内域并集**为计算域；覆盖状态、径向容差、未知薄层及超时结果的解释见 [连续构域说明](continuous_region.md)。本稿中针对单个方案的恢复证明，不能替代具体算例的约束核对、全局覆盖证书或独立 AC 检查。

固定初始拓扑的历史 Case33 数值见 [归档测试报告](../results/case33bw/latest/report.md)。该记录不代表当前包含联络线的完整重构域；当前模型及数据结构见 [统一规划模型](compact_planning.md)。

## 参考依据

- Yue Chen and Changhong Zhao. *Improved Approximation of Dispatchable Region in Radial Distribution Networks via Dual SOCP*. IEEE Transactions on Power Systems, 38(6), 5585–5597, 2023. [DOI](https://doi.org/10.1109/TPWRS.2022.3226894)。本稿参考其可行性 SOCP、强对偶、外域切割与顶点检查思路，以及第 5594 页的 FR、MR 定义；第 4.3 节是针对本项目限定模型给出的恢复证明，不归为该论文的一般结论。
- Stephen Boyd and Lieven Vandenberghe. *Convex Optimization*. 第 5.2.3 节（Slater 条件），第 5.9 节（广义不等式与锥对偶）。[作者提供的教材](https://web.stanford.edu/~boyd/cvxbook/bv_cvxbook.pdf)。
- 本项目对照文件：[网络参数](../Network/four_bus_five_corridor.py)、[现有线性 SP](../model.py)、[独立 AC 验证](../vertify.py)。
