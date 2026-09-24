# 两阶段可行规划域主线

实现入口为 main.build_continuous_region。统一符号和索引见 [notation.md](notation.md)。

## 1. 求什么域

\[
D_{\mathcal B}=\bigcup_{x\in\mathcal X,\ c^\top x\le\mathcal B}
\{p\ge0:\exists y,\ (x,p,y)\in\mathcal F_{\rm SOCP}\}.
\]

x 为走廊型号二元向量；p 为选定负荷节点的有功，kW；y=state=(P,Q,ell,v,s⁺,s⁻)。
线路和电源功率为标幺值，负荷统一除以 network.base 后进入节点守恒。
模型约束编号直接对应 model.py：

1. 声明 x、p、接入 a、连通虚拟流 f。
2. 每走廊最多一种类型，接通边两端须接入。
3. 连通流 Gf=a，边数等于非根接入节点数，构成径向树。
4. 建设费用 cᵀx≤budget，总负荷 Σp≤power_limit。
5. 可选固定负荷/固定方案，用于独立校核或热启动。
6. 节点有功、无功守恒，流入功率扣除线路损耗。
7. 电压平方上下限。
8. 压降与开断余量，接通时余量为零。
9. 型号门控容量和电流上界，反向送端功率计入损耗。
10. 线路 SOCP：P²+Q²≤v_i·ell。
11. 源端有功、无功、视在功率限额，包含损耗。
12. 全接入纯负荷树的父弧/多商品流松弛加强。
13. 查询辅助退让量 d 或 rho，不放松物理约束。

构域前提为纯负荷、非负无功比例、正阻抗、允许电压回升至根电压，因此 D 为向下封闭集。
可选节点仍由 a 决定；全接入树加强仅在全部节点必接入时启用。

## 2. Stage 1：全局支持割

从每轴长度 power_limit 的正方形/立方体开始，先加 Σp≤power_limit。
对非负方向 omega，直接求完整规划模型：

\[
h(\omega)=\max\{\omega^\top p:p\in D_{\mathcal B}\}.
\]

Gurobi 的全局上界 \(\bar h\ge h\) 给有效半空间 \(\omega^\top p\le\bar h\)。
三条/两条坐标支持先求，再按当前外域与已知点凸包的缺口选方向。
这里跨方案点凸包只用于选方向，不作为可行内域。

目标值是可行下界，不能作为外切面的截距。超时仍可使用有效全局上界；方向查询不必精确最优。
第一阶段即使精确完成也只能获得凸包，无法删除凸包内非凸缺口。

## 3. Stage 2：全局非凸排除

对外域候选 q，完整整数规划模型求

\[
d(q)=\min_{p\in D_{\mathcal B}}\max_i(q_i-p_i)_+.
\]

用辅助量 d≥0 及 p+d·1≥q 建模。取全局下界 L>0，则

\[
D_{\mathcal B}\subseteq\bigcup_i\{p:p_i\le q_i-L\}.
\]

证明：若某个可行 p 在所有坐标上都严格大于 q-L，则它产生的退让量严格小于 L，与全局下界矛盾。
保留边界上的等号；切掉的是所有坐标同时更大的开正交象限。

可选斜向版本：

\[
\rho_W(q)=\min\rho,\qquad Wp+\rho\mathbf1\ge Wq,\quad \rho\ge0,\ p\in D_{\mathcal B}.
\]
\[
D_{\mathcal B}\subseteq\bigcup_k\{p:W_kp\le W_kq-\underline\rho_W(q)\}.
\]

W 各行非负且归一化；方向来自已知方案面仅是启发式，割有效性来自完整整数模型的下界。
rho 为投影缺额，rho=0 不能说明 q 可行；此时回到 W=I 的真实坐标距离。

**固定方案连续模型的对偶/目标界不能直接代替上述全局界。**
固定方案可行点用于扩张内域、提供 Gurobi Start；完整整数分支定界覆盖所有允许建设方案。
主线不显式枚举全部网架，也不从 MIP 伪造统一的线性对偶乘子。

## 4. 覆盖证书与最大间隙优先

每个已知方案独立维护
\[
I_x=\operatorname{down}\operatorname{conv}\{p^j:x^j=x\},\quad I=\bigcup_x I_x.
\]
同一 x 下 SOCP 约束为凸约束，支持凸插值；不同 x 不可混合为一个可行凸包。

把当前外域划为凸节点 O_s，计算
\[
\bar\Delta_s=\min_x\max_{v\in V(O_s)}\operatorname{dist}_\infty(v,I_x).
\]
固定 x 时，到凸集的距离为凸函数，节点上最大距离不超过顶点最大距离；再选最好的一个 x，因此这是整个节点的有效上界。
min/max 不能交换，否则会漏掉由不同方案顶点围成的非凸缺口。

每轮选择最大上界的节点：

1. 若全部节点上界≤epsilon_kw，结束并标为 certified。
2. 有顶点距离大于 epsilon，优先查询该顶点，增加同方案可行点或用正全局下界切割。
3. 若顶点各自被不同方案覆盖，节点整体尚未认证，沿某方案距离邻域分支。
4. 若只接触低维边界而无法推进，沿节点最长坐标二分。

覆盖分支不删除任何点。不可行性析取割传播至所有相交节点；每次 OR 按顺序分成内部不交的多面体。
最终
\[
I\subseteq D_{\mathcal B}\subseteq O_2\subseteq O_1,\qquad
\sup_{q\in O_2}\operatorname{dist}_\infty(q,I)\le\varepsilon.
\]

这是三维/二维 SOCP 规划域证书，不是高精度扫描或精确 AC 投影的声明。
达时限或无充分证据时保留外域、内域和当前间隙，返回 time_limit/unknown。

## 5. 代码调用和最小结果

main 按预算循环调用 build_continuous_region；构域函数直接调用 PlanningModel.solve。
solve 根据 weights/target/cut_normals 切换目标，只执行一次 Gurobi optimize。
region.py 只做几何，plot.py 只保存核心 JSON 和最终 HTML。

结果 schema=4 保存预算、坐标、epsilon、最终间隙、内域、两阶段外域、支持割、析取割、状态和构域秒数。
不保留逐次模型、运行向量、事件流水、检查点或派生体积。
独立高精度截面和 AC 检查保留在测试工具中，不增加主线的再次求解层。
