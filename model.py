"""规划主问题、LP/SOCP 子问题及独立 AC；迭代流程在 main.ipynb。

符号约定：power 为三个独立节点的有功负荷（kW），其余负荷固定；
P/Q 为支路送端功率，v/u 为支路受端/送端电压幅值的平方，ell 为电流平方。
除 power 外，潮流量均为标幺值。支路编号与其下游节点编号一一对应。
固定方案的有效割统一写成 cut[0] + cut[1:] @ power >= 0。
逐线路规划的联合割为 cut[0] + cut[1:4] @ power + cut[4:] @ x >= 0。
"""
import gurobipy as gp
import numpy as np
from gurobipy import GRB
import clarabel
from scipy import sparse

AC_TOL = 1e-9  # AC 运行限值的标幺容差；不是切割区域的几何误差。
FIXED_POINT_TOL = 1e-12  # 完整电流等式 P²+Q²-u*ell 的收敛容差。
PLANNING_TOL = 1e-8  # 联合 SP 的原始标幺约束容差；边界的 kW 差异另与独立枚举核对。


def add_cut(master, selection, power, design, cut):
    """固定方案的割仅在该方案被选中时生效，不能用于切掉其他方案。"""
    expression = cut[0]+gp.LinExpr(cut[1:], list(power.values()))  # a+bᵀp。
    master.addGenConstrIndicator(selection[design], True, expression >= 0.)  # z_d=1 才启用该方案的割。


def _master(designs, history):
    """有限方案的等价离散主问题：选一个合法方案，保留建设决策及费用。"""
    master = gp.Model()  # MP 只选择建设方案和负荷，不在此重复建立潮流方程。
    master.Params.OutputFlag = 0  # 关闭逐次求解日志。
    master.Params.Threads = 1  # 各比较方法统一使用单线程。
    master.Params.MIPGap = master.Params.MIPGapAbs = 0  # MP 目标必须达到全局最优后才用于构域。
    master.Params.FeasibilityTol = master.Params.IntFeasTol = 1e-9  # 约束及二进制变量的数值容差。
    master.Params.DualReductions = 0  # 区分不可行与无界，避免使用含混的求解状态。
    selection = master.addVars(len(designs), vtype=GRB.BINARY, name="design")  # z_d：是否选择完整方案 d。
    power = master.addVars(len(designs[0].load_nodes), lb=0., name="p")  # p≥0，单位 kW。
    master.addConstr(selection.sum() == 1)  # 一次运行只实现一个方案，区域最后才取并集。
    master.addConstr(power.sum() <= max(d.power_limit for d in designs))  # 所有候选方案共用的有效总负荷外界。
    # 方案表已满足选型、连通性和径向性；建设向量 x=Σz_d*x_d 可直接恢复，无需再建一组变量。
    for row in history:
        if row['cut'] is not None:  # 跨 MP 查询复用已经证明有效的条件割。
            add_cut(master, selection, power, row['design'], row['cut'])
    return master, selection, power


def MP1(designs, power, history=(), *, total=False):
    """固定负荷（或总负荷下限），最小化建设费；返回未迭代的主问题。"""
    master, selection, p = _master(designs, history)  # 使用相同方案集合及历史外割。
    master.ModelName = 'MP1'
    if total:
        master.addConstr(p.sum() >= float(power))  # 总量查询：可重新分配三个节点的负荷。
    else:
        master.addConstrs(p[i] == float(value) for i,value in enumerate(power))  # 点查询：三个坐标均固定。
    master.setObjective(gp.quicksum(d.cost*selection[i] for i,d in enumerate(designs)), GRB.MINIMIZE)  # min ΣC_d*z_d。
    return master, selection, p


def MP2(designs, budget=np.inf, history=()):
    """预算内最大化三个独立节点的总负荷，固定背景负荷由 SP 保留。"""
    master, selection, p = _master(designs, history)  # MP2 与 MP1 共用同一可行方案表示。
    master.ModelName = 'MP2'
    if np.isfinite(budget):  # 无限预算表示不建立投资上限约束。
        master.addConstr(gp.quicksum(d.cost*selection[i] for i,d in enumerate(designs)) <= budget)  # ΣC_d*z_d≤B。
    master.setObjective(p.sum(), GRB.MAXIMIZE)  # 最大化独立节点负荷总量；固定背景负荷不计入此目标。
    return master, selection, p


class BranchEquations:
    """消去 P/Q/v/u，保留仿射约束 g=c+F*power+G*ell。

    前 linear_count 行要求 g≥0；后续各段要求 g 位于二阶锥。
    LP 只读取 ell=0 时的运行限值 linear_c+linear_F*power≥0。
    """

    def __init__(self, network):
        self.network = c = network  # c 在本方法中是网架对象；self.c 才是约束常数向量。
        n, D, E = c.n, c.D, c.E  # D(n×n)：下游汇总；E(n×3)：独立负荷到全网节点的映射。
        R, X = np.diag(c.r), np.diag(c.reactance)  # 支路电阻、电抗对角阵，标幺值。
        P0, Q0 = D@c.fixed_p/c.base, D@c.fixed_q/c.base  # 固定背景负荷对应的无损送端功率。
        Pp, Qp = D@E/c.base, (D@E)*c.q_ratio/c.base  # 独立负荷的功率系数；Q 按各节点固定 Q/P 比例变化。
        Pl, Ql = D@R, D@X  # 下游所有支路的 r*ell、x*ell 累加到本支路送端。
        v0 = 1-2*D.T@(R@P0+X@Q0)  # 根电压平方为 1；沿路径累加背景负荷压降。
        vp = -2*D.T@(R@Pp+X@Qp)  # 独立负荷对受端电压平方的线性系数。
        vl = -2*D.T@(R@Pl+X@Ql)+D.T@(R@R+X@X)  # 包含网损压降及 (r²+x²)*ell 修正。
        J = np.zeros((n, n))  # 将每个非根节点电压映射到其子支路的送端。
        for i, parent in enumerate(c.parent):
            if parent >= 0:  # 接根支路的送端电压是常数 1，不引用节点变量。
                J[i, parent] = 1
        u0, up, ul = 1+J@(v0-1), J@vp, J@vl  # u 为送端电压平方的常数项、负荷系数、电流系数。
        self.state_terms = ((P0, Pp, Pl), (Q0, Qp, Ql), (v0, vp, vl), (u0, up, ul))  # 用于重建物理状态。
        roots = c.roots  # 直接连接电源根节点的支路；其送端功率之和就是电源注入。
        # 依次为 P≤Pmax、v≥vmin、v≤vmax、ΣP_root≤Psmax、ΣQ_root≤Qsmax。
        constant = np.r_[c.capacity-P0, v0-c.vmin, c.vmax-v0,  # 支路有功与节点电压约束的常数余量。
                          c.source_pmax-P0[roots].sum(), c.source_qmax-Q0[roots].sum()]  # 电源 P/Q 约束的常数余量。
        F = np.vstack([-Pp, vp, -vp, -Pp[roots].sum(axis=0), -Qp[roots].sum(axis=0)])  # 每条运行限值的 power 系数。
        G = np.vstack([-Pl, vl, -vl, -Pl[roots].sum(axis=0), -Ql[roots].sum(axis=0)])  # 同一顺序下的 ell 系数。
        # 只建立网架中给定的运行限值；inf 表示该约束未启用。
        finite = np.isfinite(constant)  # 去掉未启用的运行限值，不创建人为的容量约束。
        self.linear_c, self.linear_F = constant[finite], F[finite]  # ell=0 的无损 LP 约束。
        self.c = np.r_[np.zeros(n), self.linear_c]  # 在运行限值前加 n 条 ell≥0。
        self.F = np.vstack([np.zeros((n, E.shape[1])), self.linear_F])  # ell≥0 与独立负荷无关。
        identity = np.eye(n)  # 电流平方各分量的单位方向，所有支路共用这一个单位阵。
        self.G = np.vstack([identity, G[finite]])  # 前 n 行直接提取 ell，其后为运行限值。
        self.linear_count = len(self.c)  # 非负锥段长度，后面的二阶锥另行登记。
        self.relax = np.r_[np.zeros(n), np.ones(len(self.linear_c))]  # eta 只松弛运行限值，ell≥0 保持严格。
        self.cones = []  # 每个 slice 标识一个标准二阶锥 (t,w)：t≥||w||₂。
        if np.isfinite(c.source_smax):  # (Smax,ΣP_root,ΣQ_root)∈SOC，保留电源视在容量限制。
            self._cone([c.source_smax, P0[roots].sum(), Q0[roots].sum()],  # 电源锥的常数项。
                       np.vstack([np.zeros(E.shape[1]), Pp[roots].sum(axis=0), Qp[roots].sum(axis=0)]),  # Smax 固定，P/Q 随负荷变化。
                       np.vstack([np.zeros(n), Pl[roots].sum(axis=0), Ql[roots].sum(axis=0)]))  # P/Q 还包含全网损耗。
        # ||(2P,2Q,u-ell)|| <= u+ell 等价于 P²+Q² <= u ell。
        for i in range(n):
            unit = identity[i]  # e_i*ell=ell_i；在锥首尾分别加、减该电流平方。
            self._cone([u0[i], 2*P0[i], 2*Q0[i], u0[i]],  # (u+ell,2P,2Q,u-ell) 的常数项。
                       np.vstack([up[i], 2*Pp[i], 2*Qp[i], up[i]]),  # 同一四维锥向量对独立负荷的系数。
                       np.vstack([ul[i]+unit, 2*Pl[i], 2*Ql[i], ul[i]-unit]))  # 首尾的 ±unit 实现 ±ell_i。

    def _cone(self, c, F, G):
        self.cones.append(slice(len(self.c), len(self.c)+len(c)))  # 记录新锥在拼接向量中的区间。
        self.c, self.F, self.G = np.r_[self.c, c], np.vstack([self.F, F]), np.vstack([self.G, G])  # 同步拼接仿射系数。
        self.relax = np.r_[self.relax, 1., np.zeros(len(c)-1)]  # eta 只加到锥轴 t，使 t+eta≥||w||₂。

    def state(self, power, ell):
        return tuple(c+F@power+G@ell for c, F, G in self.state_terms)  # 依次重建 P、Q、v、u。

    def margin(self, values):
        branches = values[-4*self.network.n:].reshape(-1, 4)  # 每条支路占 (u+ell,2P,2Q,u-ell) 四行。
        # 非负锥取最小分量，二阶锥取 t-||w||₂；所有余量非负才是一个可行证书。
        margin = min(values[:self.linear_count].min(),  # ell≥0 及线性限值中最小的余量。
                     np.min(branches[:, 0]-np.linalg.norm(branches[:, 1:], axis=1)))  # 所有支路锥中最小的余量。
        if np.isfinite(self.network.source_smax):
            s = self.cones[0]  # 若电源视在容量启用，它排在所有支路锥之前。
            margin = min(margin, values[s.start]-np.linalg.norm(values[s.start+1:s.stop]))  # 再纳入电源容量锥余量。
        return margin


class LinearSP:
    """固定方案无损 LP：min eta，c+F*power+eta≥0，eta≥0。

    查询 power 固定；eta=0 表示该点满足所有 LP 约束。
    返回 (违反量, 有效外割或 None, 可行内点, 电流证书或 None)。
    """

    def __init__(self, network):
        e = BranchEquations(network)  # 与 SOCP 共用物理参数消元，LP 只读取 ell=0 的运行限值。
        self.c = np.r_[e.linear_c, network.power_limit/network.base]  # 加入独立节点总有功负荷上界。
        self.F = np.vstack([e.linear_F, -np.ones(len(network.load_nodes))/network.base])  # 将 kW 总负荷换成标幺值。
        self.model = gp.Model('fixed_linear_SP')
        self.model.Params.OutputFlag = 0
        self.model.Params.Threads = 1
        self.model.Params.FeasibilityTol = self.model.Params.OptimalityTol = 1e-9  # 原始可行性和对偶最优性精度。
        self.eta = self.model.addVar(lb=0, obj=1)  # 唯一优化变量：所有运行限值共用的非负违反量。
        self.rows = [self.model.addConstr(self.eta >= -c) for c in self.c]  # 先建约束，查询时只改 RHS。
        self.calls = 0  # SP 实际求解次数，不是事后绘图网格点数。

    def solve(self, power):
        slope = self.F@power  # 当前负荷对全部运行限值的贡献；后续构造内点继续复用。
        self.model.setAttr('RHS', self.rows, -self.c-slope)  # eta≥-(c+F*power)。
        self.model.optimize()  # 求最小违反量及各行对偶乘子。
        self.calls += 1
        if self.model.Status != GRB.OPTIMAL:  # 非最优结果不能用来认证或构造对偶割。
            raise RuntimeError('Linear SP was not solved')
        eta = self.eta.X  # eta 是标幺约束违反量，不是区域距离或 AC 误差。
        if eta <= 1e-9:  # 容差内的零：直接返回查询点作为 LP 可行证书。
            return eta, None, power.copy(), None
        dual = np.array(self.model.getAttr('Pi', self.rows))  # min 问题的 ≥ 行对应非负乘子 λ。
        decreasing = slope < 0  # 沿 power 方向增载时，只有这些约束的余量会变小。
        alpha = min(1., np.min(-self.c[decreasing]/slope[decreasing]))*(1-1e-10)  # c+alpha*F*power≥0 的最大径向比例。
        return eta, np.r_[dual@self.c, self.F.T@dual], alpha*power, None  # λᵀc+(Fᵀλ)ᵀp≥0；内点保持固定背景负荷。

    def close(self):
        self.model.dispose()


class SOCPSP:
    """固定方案 SOCP：min eta，c+F*power+G*ell+relax*eta∈K，eta≥0。

    K 是非负锥与二阶锥的直积。power 固定，优化变量是 (ell, eta)。
    该模型将 AC 电流等式松弛为 P²+Q²≤u*ell；独立 AC 另行验证等式。
    """

    def __init__(self, network):
        self.equations = e = BranchEquations(network)  # 一个方案构建一次仿射约束。
        n = network.n  # ell 的分量数等于支路数。
        # Clarabel 使用 A*x+s=b；取 A=-[G,relax]、b=c+F*power，即 s=g+relax*eta。
        matrix = sparse.csc_matrix(np.vstack([-np.c_[e.G, e.relax], np.r_[np.zeros(n), -1.]]))
        cones = [clarabel.NonnegativeConeT(e.linear_count)]  # ell≥0 及全部线性运行限值。
        cones += [clarabel.SecondOrderConeT(s.stop-s.start) for s in e.cones]  # 电源容量锥及各支路功率锥。
        cones += [clarabel.NonnegativeConeT(1)]  # matrix 最后一行 -eta+s=0 实现 eta≥0。
        settings = clarabel.DefaultSettings()
        settings.verbose = False
        settings.max_threads = 1
        settings.presolve_enable = False  # 保留原始约束排列，便于逐行更新和解释对偶乘子。
        settings.tol_gap_abs = settings.tol_gap_rel = settings.tol_feas = 1e-10  # 控制原始、对偶和间隙精度。
        # Hessian 为零，目标系数 (0,…,0,1)：仅最小化 eta，不最小化网损。
        self.model = clarabel.DefaultSolver(
            sparse.csc_matrix((n+1, n+1)), np.r_[np.zeros(n), 1.],  # 二次目标为零，线性目标为 eta。
            matrix, np.r_[e.c, 0.], cones, settings)  # 初始 RHS 对应独立负荷为零，查询时再更新。
        self.calls = 0  # 包括固定背景的初始证书查询。
        self.background = bool(np.any(network.fixed_p) or np.any(network.fixed_q))  # p=0 时是否仍有固定负荷。
        self.base_currents = np.zeros(n)  # 无固定负荷时，零电流本身就是原点证书。
        self.current_bound = np.minimum(network.capacity, min(network.source_pmax, network.source_smax))/network.r  # P≥r*ell 给出的有效 ell 上界。
        if self.background:
            # p=0 只清空独立坐标；先求同一固定背景下的可行原点证书。
            initial = self._solve(e.c-1e-7*e.relax)  # 轻微收紧运行限值，获得有数值余量的原点证书。
            self.base_currents = np.maximum(initial.x[:-1], 0.)  # 去除求解容差产生的极小负电流平方。
        self.base_slack = e.c+e.G@self.base_currents  # 此处 power=0，但固定背景已包含在 e.c 中。
        if e.margin(self.base_slack) < 0:  # 后续凸插值依赖可行原点；不能用无效原点继续认证。
            raise RuntimeError('The parameter origin has no feasible SOCP certificate')

    def _solve(self, rhs):
        self.model.update(b=np.r_[rhs, 0.])  # 仅更新查询负荷；最后一行 eta≥0 的 RHS 始终为零。
        result = self.model.solve()  # A、目标及锥结构在不同查询间保持不变。
        self.calls += 1
        if result.status not in (clarabel.SolverStatus.Solved, clarabel.SolverStatus.AlmostSolved):  # 内点和外割仍会单独修正、核验。
            raise RuntimeError(f'SOCP phase I status {result.status}')
        return result

    def feasible_scale(self, power, ell):
        """在参数原点与查询状态间内缩；固定背景负荷在整个线段上不变。"""
        e = self.equations
        slack = e.c+e.F@power+e.G@ell  # 去掉辅助违反量 eta，核查原始 SOCP 约束。
        if self.background:
            low, high = 0., 1.  # alpha=0 是已核验原点，alpha=1 是当前查询状态。
            if e.margin(slack) >= 0:
                low = 1.  # 当前状态已满足原始约束，不需要二分。
            else:
                for _ in range(36):
                    alpha = (low+high)/2  # 沿同一固定背景下的凸线段二分。
                    if e.margin(self.base_slack+alpha*(slack-self.base_slack)) >= 0:
                        low = alpha  # 保存已经逐锥核验可行的一侧。
                    else:
                        high = alpha  # 不可行的一侧只能作为搜索上界。
            alpha = low*(1-1e-8)  # 轻微内缩，避免将浮点边界误差当作内点证书。
        else:
            # 无背景负荷时，所有功率同比缩放，可直接计算锥和线性约束的缩放界。
            P, Q, _, u = e.state(power, ell)  # 无背景时 P/Q/ell 同比缩放，u 变为 1+alpha*(u-1)。
            delta = self.base_slack[:e.linear_count]-slack[:e.linear_count]  # 每条线性余量沿线段的下降幅度。
            positive = delta > 0  # 只需约束会下降的余量。
            bounds = [1., *list(self.base_slack[:e.linear_count][positive]/delta[positive])]  # 线性约束给出的 alpha 上界。
            denominator = P*P+Q*Q+(1-u)*ell  # 支路锥整理为 alpha*[P²+Q²+(1-u)ell]≤ell。
            positive = denominator > 0
            bounds.extend(ell[positive]/denominator[positive])  # 各支路锥的径向上界。
            source = np.hypot(P[e.network.roots].sum(), Q[e.network.roots].sum())  # 送端总视在功率。
            if source > 0:
                bounds.append(e.network.source_smax/source)  # 电源视在容量要求 alpha*S≤Smax。
            alpha = max(0., min(bounds))*(1-1e-8)  # 同时满足所有线性、支路锥及电源锥约束。
        return alpha*np.asarray(power), self.base_currents+alpha*(ell-self.base_currents)  # 返回负荷内点及其电流平方证书。

    def solve(self, power):
        e, c = self.equations, self.equations.network
        result = self._solve(e.c+e.F@power)  # 固定查询负荷，求 (ell,eta)。
        ell = np.maximum(result.x[:-1], 0.)  # 最后一项是 eta，前 n 项是电流平方。
        p, q = c.loads(power)  # 完整节点负荷，用于识别没有任何功率需求的子树。
        ell[c.D@(p[0]+q[0]) == 0] = 0.  # 空子树不保留数值伪损耗。
        witness, currents = self.feasible_scale(power, ell)  # 原始约束中的可行内点；不能直接用含 eta 的解认证。
        dual = np.asarray(result.z[:-1])  # 舍去 eta≥0 那一行的乘子，保留物理约束乘子 y。
        dual[:e.linear_count] = np.maximum(dual[:e.linear_count], 0.)  # 将浮点乘子修正到非负对偶锥。
        for s in e.cones:
            dual[s.start] = max(dual[s.start], np.linalg.norm(dual[s.start+1:s.stop]))  # SOC 自对偶：强制 y_t≥||y_w||₂。
        # yᵀ(c+Fp+G*ell)≥0；用 0≤ell≤L 补偿 Gᵀy 的驻点残差，避免无效外割。
        correction = np.maximum(e.G.T@dual, 0)@self.current_bound  # (Gᵀy)_+ᵀL 是所需的常数补偿上界。
        cut = np.r_[dual@e.c+correction+1e-10, e.F.T@dual]  # 修正后 a+bᵀp≥0 包含所有可行负荷。
        if cut[0]+cut[1:]@power >= -1e-9:  # 该割未可靠分离查询点时，只返回内点，不添加近似重复割。
            cut = None
        return max(0., result.x[-1]), cut, witness, currents  # SOCP 认证不等同于独立 AC 等式认证。

    def close(self):
        self.model = None


class PlanningEquations:
    """固定径向拓扑的逐线路选型模型：c+F*p+G*y∈K，0≤y≤U0+Ux*x。

    x 是各候选线路的 one-hot 型号变量；y=(P_k,Q_k,ell_k,v)。
    各型号保留自己的送端功率和电流平方，共用真实节点电压；未选型号潮流为零。
    适用于当前 case33：非负负荷、正阻抗、有限源端 P/Q 上限；不补造线路热限。
    """

    def __init__(self, network, method):
        self.network = c = network
        groups = {o.branch: o for o in c.line_options}  # 型号、电气参数和费用只从 Network 读取。
        edges, r, reactance, cost, selected, self.groups = [], [], [], [], [], []
        for e in range(c.n):
            o = groups.get(e)
            start = len(cost)
            types = zip(o.r, o.reactance, o.cost) if o else [(c.r[e], c.reactance[e], 0.)]
            for resistance, xline, price in types:
                edges.append(e); r.append(resistance); reactance.append(xline)
                if o:
                    selected.append(len(edges)-1)  # 只有候选线路的型号对应二进制变量。
                    cost.append(price)
            if o:
                self.groups.append(np.arange(start, len(cost)))  # 每条候选线路恰选一个型号。
        self.edges, self.r, self.reactance = np.array(edges), np.array(r), np.array(reactance)
        self.cost, self.selected = np.array(cost), np.array(selected)
        t, nx, npower = len(edges), len(cost), len(c.load_nodes)
        self.P, self.Q, self.ell = slice(0,t), slice(t,2*t), slice(2*t,3*t)
        self.v = slice(3*t,3*t+c.n)
        ny = 3*t+c.n
        T = np.eye(c.n)[:, self.edges]  # 按支路汇总各型号的物理量。
        J = np.zeros((c.n,c.n))
        for e, parent in enumerate(c.parent):
            if parent >= 0:
                J[e,parent] = 1.  # J*v 为非根送端电压；根送端另加常数 1。
        root = (c.parent < 0).astype(float)
        self.T = T
        balance = (np.eye(c.n)-J.T)@T  # 本支路送出量减去直接子支路送出量。
        eq = np.zeros((3*c.n, ny))
        eq[:c.n,self.P], eq[:c.n,self.ell] = balance, -T*self.r  # P-ΣP_child-r*ell=p_node。
        eq[c.n:2*c.n,self.Q], eq[c.n:2*c.n,self.ell] = balance, -T*self.reactance
        eq[2*c.n:3*c.n,self.v] = J-np.eye(c.n)
        eq[2*c.n:3*c.n,self.P], eq[2*c.n:3*c.n,self.Q] = -2*T*self.r, -2*T*self.reactance
        eq[2*c.n:3*c.n,self.ell] = T*(self.r**2+self.reactance**2)  # 完整支路压降的损耗项。
        eqc = np.r_[-c.fixed_p/c.base, -c.fixed_q/c.base, root]
        eqf = np.vstack([-c.E/c.base, -c.E*c.q_ratio/c.base, np.zeros((c.n,npower))])
        self.equal_count = len(eqc)
        # 原模型保留等式；phase I 使用正负两行，把等式残差也纳入 eta。
        self.c, self.F = np.r_[eqc,-eqc], np.vstack([eqf,-eqf])
        self.G = np.vstack([eq,-eq])
        pmax = np.minimum(c.capacity[self.edges], min(c.source_pmax,c.source_smax))
        qmax = np.full(t,min(c.source_qmax,c.source_smax))
        lmax = pmax/self.r if method == 'socp' else np.zeros(t)  # P≥r*ell 给出有效电流平方上界；线性模型令 ell=0。
        self.upper = np.r_[pmax,qmax,lmax,c.vmax]  # 所有 y 的物理有效界，用于型号联接和对偶残差补偿。
        self.box_constant = self.upper.copy()
        self.box_selection = np.zeros((ny,nx))
        self.linked = np.concatenate([offset+self.selected for offset in (0,t,2*t)])
        for offset, bound in ((0,pmax),(t,qmax),(2*t,lmax)):
            self.box_constant[offset+self.selected] = 0.
            self.box_selection[offset+self.selected,np.arange(nx)] = bound[self.selected]  # 型号相关上界 U(x)；未选型号的潮流严格为零。
        linear_c, linear_g = [], []
        g = np.zeros((c.n,ny)); g[:,self.v] = np.eye(c.n)
        linear_c.extend(-c.vmin); linear_g.extend(g)
        for block, limit in ((self.P,c.source_pmax),(self.Q,c.source_qmax)):
            if np.isfinite(limit):
                g = np.zeros(ny); g[block] = -root[self.edges]
                linear_c.append(limit); linear_g.append(g)  # 源端功率包含网损。
        self.c = np.r_[self.c,linear_c]
        self.F = np.vstack([self.F,np.zeros((len(linear_c),npower))])
        self.G = np.vstack([self.G,linear_g])
        self.linear_count, self.cones = len(self.c), []
        if method == 'socp':
            for k in range(t):
                g = np.zeros((4,ny))
                e = self.edges[k]
                g[0,self.v] = g[3,self.v] = J[e]
                g[0,self.ell.start+k] = 1.
                g[1,self.P.start+k] = g[2,self.Q.start+k] = 2.
                g[3,self.ell.start+k] = -1.
                self._cone([root[e],0.,0.,root[e]],g)  # (u_e+ell_k,2P_k,2Q_k,u_e-ell_k)∈SOC；u_e 是真实送端电压平方。
            if np.isfinite(c.source_smax):
                g = np.zeros((3,ny)); g[1,self.P], g[2,self.Q] = root[self.edges],root[self.edges]
                self._cone([c.source_smax,0.,0.],g)
        self.relax = np.zeros(len(self.c)); self.relax[:self.linear_count] = 1.
        for s in self.cones:
            self.relax[s.start] = 1.  # 锥只松弛首分量，不更改锥的方向分量。

    def _cone(self, constant, matrix):
        self.cones.append(slice(len(self.c),len(self.c)+len(constant)))
        self.c = np.r_[self.c,constant]; self.G = np.vstack([self.G,matrix])
        self.F = np.vstack([self.F,np.zeros((len(constant),self.F.shape[1]))])

    def selection(self, choice):
        """线路型号编号转 one-hot；只用于实例化方案和枚举对照。"""
        x = np.zeros(len(self.cost))
        for group,k in zip(self.groups,choice):
            x[group[int(k)]] = 1.
        return x

    def choice(self, x):
        return np.array([np.argmax(x[g]) for g in self.groups])

    def margin(self, x, power, state):
        g = self.c+self.F@power+self.G@state
        values = [g[:self.linear_count].min(),np.min(state),
                  np.min(self.box_constant+self.box_selection@x-state)]
        values.extend(g[s.start]-np.linalg.norm(g[s.start+1:s.stop]) for s in self.cones)
        return min(values)  # 不含 eta 的原约束余量；等式正负两行等价于检查绝对残差。

    def restore(self, x, power, state):
        """用选中型号的电流重建等式精确的状态，消除未选型号退化锥中的浮点噪声。"""
        c = self.network
        active = np.ones(len(self.edges)); active[self.selected] = x
        y = np.zeros_like(state)
        ell = np.clip(state[self.ell],0.,self.upper[self.ell])*active
        p,q = c.loads(power)
        P = c.D@(p[0]+self.T@(self.r*ell))
        Q = c.D@(q[0]+self.T@(self.reactance*ell))
        y[self.ell], y[self.P], y[self.Q] = ell, P[self.edges]*active, Q[self.edges]*active
        v = 1.-2*c.D.T@(self.T@(self.r*y[self.P]+self.reactance*y[self.Q]))
        v += c.D.T@(self.T@((self.r**2+self.reactance**2)*ell))
        y[self.v] = v
        return y  # 是否可行仍由 margin 对原始限值与锥逐一核验。


class PlanningModel:
    """逐线路 MILP/MISOCP；cuts_only=True 时为包含 x、p 的联合割主问题。"""

    def __init__(self, network, method, *, power=None, budget=np.inf, direction=None, cuts_only=False):
        self.equations = e = PlanningEquations(network,method)
        self.model = m = gp.Model('compact_'+method)
        m.Params.OutputFlag = 0; m.Params.Threads = 1
        m.Params.MIPGap = m.Params.MIPGapAbs = 0.
        m.Params.FeasibilityTol = m.Params.IntFeasTol = m.Params.OptimalityTol = 1e-9
        m.Params.BarQCPConvTol = 1e-9
        m.Params.DualReductions = 0
        m.Params.NonConvex = 0  # 若模型误写成非凸二次式，立即暴露，不能悄悄改用非凸求解。
        self.x = m.addMVar(len(e.cost),vtype=GRB.BINARY,name='line_type')
        self.power = p = m.addMVar(len(network.load_nodes),lb=0.,name='p_kw')
        for group in e.groups:
            m.addConstr(self.x[group].sum()==1.)  # 每条线路恰好保留或升级为一个型号。
        m.addConstr(p.sum()/network.base<=network.power_limit/network.base)
        if np.isfinite(budget):
            m.addConstr(e.cost@self.x<=budget)
        if power is not None:
            m.addConstr(p==np.asarray(power)); m.setObjective(e.cost@self.x,GRB.MINIMIZE)
            self.objective_scale = 1.
        else:
            if direction is not None:
                scale = m.addMVar(1,lb=0.,name='radial_total_kw')
                m.addConstr(p==(np.asarray(direction)/np.sum(direction))[:,None]@scale)  # 射线比例固定，背景负荷不缩放。
            m.setObjective(p.sum()/network.base,GRB.MAXIMIZE)  # 目标用标幺量改善锥求解精度；返回时恢复 kW。
            self.objective_scale = network.base
        self.state = None
        if not cuts_only:
            self.state = y = m.addMVar(len(e.upper),lb=0.,ub=e.upper,name='flow')
            g = e.c+e.F@p+e.G@y
            m.addConstr(g[:e.equal_count]==0.)  # 原模型等式只建一次，负号副本仅用于 phase I。
            m.addConstr(g[2*e.equal_count:e.linear_count]>=0.)
            m.addConstr(y[e.linked]<=e.box_selection[e.linked]@self.x)  # 0≤y_ek≤U_ek*x_ek，未选型号为零。
            if method == 'socp':
                for k in range(len(e.edges)):
                    P,Q,l = (y[s.start+k].item() for s in (e.P,e.Q,e.ell))
                    parent = network.parent[e.edges[k]]
                    u = 1. if parent < 0 else y[e.v.start+parent].item()
                    m.addQConstr(P*P+Q*Q<=u*l)  # u、ell 非负，这是凸旋转二阶锥。
                if np.isfinite(network.source_smax):
                    roots = np.isin(e.edges,network.roots)
                    ps,qs = y[e.P][roots].sum().item(),y[e.Q][roots].sum().item()
                    m.addQConstr(ps*ps+qs*qs<=network.source_smax**2)

    def add_cut(self, cut):
        n = len(self.equations.network.load_nodes)
        self.model.addConstr(cut[0]+cut[1:1+n]@self.power+cut[1+n:]@self.x>=0.)  # 联合割对全部合法选型有效。

    def solve(self):
        self.model.optimize()
        if self.model.Status == GRB.INFEASIBLE:
            return None
        if self.model.Status != GRB.OPTIMAL:
            raise RuntimeError(f'Compact planning status {self.model.Status}')
        x = np.rint(self.x.X).astype(int)
        objective = self.equations.cost@x if self.model.ModelSense==GRB.MINIMIZE else self.power.X.sum()
        # 投资按实际整数选型计价，避免 1+1e-10 被误判为超出预算 1。
        return dict(x=x,p=self.power.X,
                    objective=float(objective),bound=self.model.ObjBound*self.objective_scale,
                    state=None if self.state is None else self.state.X)


class PlanningSP:
    """固定整数 x,p 的连续 phase I；盒上界随 x 变化，割为 a+bᵀp+dᵀx≥0。"""

    def __init__(self, equations):
        self.equations = equations
        self.settings = settings = clarabel.DefaultSettings()
        settings.verbose = False; settings.max_threads = 1
        settings.tol_gap_abs = settings.tol_gap_rel = settings.tol_feas = 1e-12
        settings.reduced_tol_feas = 1e-11

    def solve(self, x, power):
        e = self.equations
        upper = e.box_constant+e.box_selection@x  # 选型只改变变量盒的右端项，物理方程系数不变。
        columns = np.flatnonzero(upper>0.)  # 未选型号及 LP 的 ell 恒为零，直接消去，不求解退化锥。
        active = np.ones(len(e.edges),dtype=bool); active[e.selected] = np.asarray(x,dtype=bool)
        cones = [s for k,s in enumerate(e.cones) if k>=len(active) or active[k]]
        rows = np.r_[np.arange(e.linear_count),*[np.arange(s.start,s.stop) for s in cones]].astype(int)
        n = len(columns)
        G = e.G[np.ix_(rows,columns)]
        matrix = np.vstack([-np.c_[G,e.relax[rows]],
                            np.c_[-np.eye(n),np.zeros(n)],np.c_[np.eye(n),np.zeros(n)],
                            np.r_[np.zeros(n),-1.]])  # 物理松弛、y≥0、y≤U(x)、eta≥0。
        rhs = np.r_[(e.c+e.F@power)[rows],np.zeros(n),upper[columns],0.]
        kinds = [clarabel.NonnegativeConeT(e.linear_count)]
        kinds += [clarabel.SecondOrderConeT(s.stop-s.start) for s in cones]
        kinds += [clarabel.NonnegativeConeT(2*n+1)]
        solver = clarabel.DefaultSolver(sparse.csc_matrix((n+1,n+1)),np.r_[np.zeros(n),1.],
            sparse.csc_matrix(matrix),rhs,kinds,self.settings)
        result = solver.solve()
        state = np.zeros(len(e.upper)); state[columns] = result.x[:-1]
        state = e.restore(x,power,state)  # 用电流重建功率平衡与压降等式，再检查原始约束。
        eta = max(0.,result.x[-1])
        if e.margin(x,power,state)>=-PLANNING_TOL:
            return dict(eta=eta,cut=None,state=state)
        if result.status not in (clarabel.SolverStatus.Solved,clarabel.SolverStatus.AlmostSolved,
                                 clarabel.SolverStatus.InsufficientProgress):
            raise RuntimeError(f'Planning phase I status {result.status}')
        dual = np.zeros(len(e.c)); dual[rows] = result.z[:len(rows)]  # 被消去的锥乘子补零，属于其对偶锥。
        dual[:e.linear_count] = np.maximum(dual[:e.linear_count],0.)
        for s in e.cones:
            dual[s.start] = max(dual[s.start],np.linalg.norm(dual[s.start+1:s.stop]))
        residual = np.maximum(e.G.T@dual,0.)  # max_{0≤y≤U(x)} (Gᵀλ)ᵀy = (Gᵀλ)_+ᵀU(x)。
        cut = np.r_[dual@e.c+residual@e.box_constant+1e-12,
                    e.F.T@dual,e.box_selection.T@residual]  # 盒支撑函数给出 x 系数，同时补偿全部驻点残差。
        cut /= max(abs(cut[0]),np.max(np.abs(cut[1:1+len(power)]))*e.network.base,
                   np.max(np.abs(cut[1+len(power):])))  # 正比例归一化，不改变半空间。
        if not (cut[0]+cut[1:1+len(power)]@power+cut[1+len(power):]@x < -1e-9):
            raise RuntimeError('Planning phase I has neither a feasible witness nor a separating cut')
        return dict(eta=eta,cut=cut,state=state)


def dispatch_support(network, method, normal, *, direction=None):
    """枚举核对用的独立连续优化：沿用旧 BranchEquations 消元式，不读取紧凑模型矩阵。

    最大化 normalᵀp；direction 给定时限制在该射线上。返回原始值及对偶上界。
    """
    e = BranchEquations(network)
    basis = np.eye(len(network.load_nodes)) if direction is None else (
        np.asarray(direction)/np.sum(direction))[:,None]
    n = basis.shape[1]
    if method == 'linear':
        c, matrix = e.linear_c, e.linear_F@basis
        cones = [clarabel.NonnegativeConeT(len(c))]
    else:
        c, matrix = e.c, np.c_[e.F@basis,e.G]
        cones = [clarabel.NonnegativeConeT(e.linear_count)]
        cones += [clarabel.SecondOrderConeT(s.stop-s.start) for s in e.cones]
    m = matrix.shape[1]
    limits = np.zeros((n+1,m)); limits[:n,:n] = np.eye(n)
    limits[-1,:n] = -basis.sum(axis=0)/network.base  # 原先共同的总负荷外界。
    matrix = sparse.csc_matrix(-np.vstack([matrix,limits]))
    rhs = np.r_[c,np.zeros(n),network.power_limit/network.base]
    cones += [clarabel.NonnegativeConeT(n+1)]
    objective = np.r_[-np.asarray(normal)@basis,np.zeros(m-n)]
    settings = clarabel.DefaultSettings(); settings.verbose = False; settings.max_threads = 1
    settings.tol_gap_abs = settings.tol_gap_rel = settings.tol_feas = 1e-10
    solver = clarabel.DefaultSolver(sparse.csc_matrix((m,m)),objective,matrix,rhs,cones,settings)
    result = solver.solve()
    if result.status not in (clarabel.SolverStatus.Solved,clarabel.SolverStatus.AlmostSolved):
        raise RuntimeError(f'Enumeration support status {result.status}')
    return dict(p=basis@np.asarray(result.x[:n]),value=-result.obj_val,bound=-result.obj_val_dual)


class ACPowerFlow:
    """独立完整 AC；仅读取网架数据，不使用 LP/SOCP 的方程矩阵或乘子。

    适用于非负 P/Q 负荷、正阻抗、根电压固定为 1 p.u. 的径向网络。
    不读取 BranchEquations 的矩阵；以支路递推独立实现 AC 电流等式。
    """

    def __init__(self, network):
        self.network = network  # AC 与 LP/SOCP 只共享网架及运行限值数据。
        self.model = None  # 非凸模型仅在不动点无法判定时创建。

    def state(self, power, ell):
        p, q = self.network.loads(power)  # 补入固定背景负荷，并将 kW/kvar 换成标幺量。
        return self._state(p, q, ell)  # 可一次重建多组负荷状态。

    def _state(self, p, q, ell):
        c = self.network
        ell = np.asarray(ell).reshape(-1, c.n)  # 行是查询样本，列是支路电流平方。
        P, Q = p+ell*c.r, q+ell*c.reactance  # 本节点负荷加本支路有功/无功损耗。
        for i in reversed(c.order):  # 从叶到根，先汇总子树功率。
            if c.parent[i] >= 0:
                P[:, c.parent[i]] += P[:, i]  # P_parent=p_parent+r_parent*ell_parent+ΣP_child。
                Q[:, c.parent[i]] += Q[:, i]  # 无功平衡使用同样的子树汇总。
        v, u = np.ones_like(P), np.ones_like(P)  # 根电压平方为 1。
        for i in c.order:  # 从根到叶，根据上一级电压逐支路计算压降。
            if c.parent[i] >= 0:
                u[:, i] = v[:, c.parent[i]]  # 当前支路送端电压等于父节点电压。
            # 完整支路压降：v=u-2(rP+xQ)+(r²+x²)*ell，保留二次电流项。
            v[:, i] = (u[:, i]-2*(c.r[i]*P[:, i]+c.reactance[i]*Q[:, i])
                       +(c.r[i]**2+c.reactance[i]**2)*ell[:, i])
        return P, Q, v, u

    def violation(self, P, Q, v):
        c = self.network
        ps, qs = P[:, c.roots].sum(axis=1), Q[:, c.roots].sum(axis=1)  # 含网损的电源送出功率。
        # 所有项目均写为“违反量”：≤0 才满足全部限值，未启用上界产生 -inf。
        return np.maximum.reduce([
            np.max(c.vmin-v, axis=1), np.max(v-c.vmax, axis=1),  # 电压平方下限和上限。
            np.max(P-c.capacity, axis=1), ps-c.source_pmax, qs-c.source_qmax,  # 支路及电源 P/Q 限值。
            np.hypot(ps, qs)-c.source_smax])  # 电源视在功率上限。

    def classify(self, power, return_currents=False):
        """1 可行，-1 已证不可行，0 未确定；不把迭代失败当作不可行。"""
        c = self.network
        power = np.asarray(power).reshape(-1, len(c.load_nodes))  # 每个样本包含三个独立负荷坐标。
        p, q = c.loads(power)  # 背景负荷只组装一次，所有不动点迭代复用。
        ell = np.zeros((len(power), c.n))  # 从零电流开始，构造单调递增的电流平方序列。
        status = np.zeros(len(power), dtype=np.int8)  # 初始均为未确定；0 绝不能直接算作不可行。
        active = np.arange(len(power))  # 后续仅更新尚未得到证书的样本。
        for _ in range(160):
            if not len(active):  # 所有样本均已获得可行或不可行证书。
                break
            P, Q, v, u = self._state(p[active], q[active], ell[active])  # 使用当前电流下界重建支路状态。
            # ell 从零单调递增：功率是下界、电压是上界；仅这些越限能提前拒绝。
            ps, qs = P[:, c.roots].sum(axis=1), Q[:, c.roots].sum(axis=1)  # 当前电源功率也是最终功率的下界。
            bad = (np.any(v < c.vmin-AC_TOL, axis=1)  # 电压上界已低于下限，后续损耗只会进一步降压。
                   | np.any(P > c.capacity+AC_TOL, axis=1)  # 支路功率下界已经超过上限。
                   | (ps > c.source_pmax+AC_TOL) | (qs > c.source_qmax+AC_TOL)  # 电源 P/Q 下界已超过运行上限。
                   | (np.hypot(ps, qs) > c.source_smax+AC_TOL) | np.any(u <= 0, axis=1))  # 视在容量越限或正电压条件已不可能满足。
            residual = np.max(np.abs(P*P+Q*Q-u*ell[active]), axis=1)  # AC 等号残差，不能改成 SOCP 不等号。
            good = ~bad & (residual <= FIXED_POINT_TOL) & np.all(v <= c.vmax+AC_TOL, axis=1)  # 收敛且全部运行限值满足。
            status[active[bad]], status[active[good]] = -1, 1  # 分别保存不可行和可行证书。
            keep = ~(bad | good)  # 仅未确定样本继续迭代。
            ell[active[keep]] = (P[keep]**2+Q[keep]**2)/u[keep]  # 完整 AC 电流等式的不动点更新。
            active = active[keep]  # 数值未收敛的样本在迭代结束后仍保留 0 状态。
        return (status, ell) if return_currents else status  # 可选电流证书用于重建相量交叉核验。

    def scan(self, points):
        status = np.empty(len(points), dtype=np.int8)  # 与输入网格点一一对应的三态判定。
        for start in range(0, len(points), 8192):  # 分批限制电流矩阵的内存，判定方程不变。
            status[start:start+8192] = self.classify(points[start:start+8192])
        return status

    def _build_global(self, environment):
        """显式非凸 AC 等式模型，用于未确定点和独立交叉核验。"""
        c = self.network
        m = gp.Model('independent_AC_reference', env=environment)
        m.Params.OutputFlag = 0
        m.Params.Threads = 1
        m.Params.NonConvex = 2  # 用空间分支定界处理 P²+Q²=u*ell 的非凸等式。
        m.Params.FeasibilityTol = m.Params.OptimalityTol = AC_TOL  # 与独立 AC 限值核验保持同级精度。
        m.Params.DualReductions = 0  # 保留明确的不可行状态。
        m.Params.TimeLimit = 30  # 超时无证书时仍返回未确定，不能据此判为不可行。
        P = m.addVars(c.n, lb=0, ub=c.capacity.tolist())  # 纯负荷径向网中送端有功非负，并受已给定支路上限约束。
        Q = m.addVars(c.n, lb=0, ub=min(c.source_qmax, c.source_smax))  # 支路无功不超过电源总无功/视在上界。
        v = m.addVars(c.n, lb=c.vmin.tolist(), ub=c.vmax.tolist())  # 变量是电压平方，网架限值已平方。
        bound = min(c.source_smax**2, c.source_pmax**2+c.source_qmax**2)/c.vmin.min()  # ell=(P²+Q²)/u 的有效全网上界。
        ell = m.addVars(c.n, lb=0, ub=bound)  # 每条支路一个非负电流平方变量。
        bp, bq = [], []  # 保存功率平衡行，后续查询只修改节点负荷 RHS。
        for i in range(c.n):
            up = 1. if c.parent[i] < 0 else v[int(c.parent[i])]  # 根节点电压固定，其余取父节点变量。
            # P_i-ΣP_child-r_i*ell_i=p_i；Q_i-ΣQ_child-x_i*ell_i=q_i。
            bp.append(m.addConstr(P[i]-gp.quicksum(P[int(j)] for j in c.children[i])-c.r[i]*ell[i] == 0))  # 查询时 RHS 替换为 p_i。
            bq.append(m.addConstr(Q[i]-gp.quicksum(Q[int(j)] for j in c.children[i])-c.reactance[i]*ell[i] == 0))  # 查询时 RHS 替换为 q_i。
            m.addConstr(v[i] == up-2*(c.r[i]*P[i]+c.reactance[i]*Q[i])  # 完整 AC 支路压降。
                        +(c.r[i]**2+c.reactance[i]**2)*ell[i])
            m.addQConstr(P[i]*P[i]+Q[i]*Q[i] == up*ell[i])  # AC 必须保留等号。
        ps, qs = (gp.quicksum(values[int(i)] for i in c.roots) for values in (P, Q))  # 电源送出功率包含网损。
        if np.isfinite(c.source_pmax):
            m.addConstr(ps <= c.source_pmax)  # 电源有功上限。
        if np.isfinite(c.source_qmax):
            m.addConstr(qs <= c.source_qmax)  # 电源无功上限。
        if np.isfinite(c.source_smax):
            m.addQConstr(ps*ps+qs*qs <= c.source_smax**2)  # 电源视在容量圆。
        m.setObjective(0.)  # 这里只判可行性，不额外改变运行目标。
        m.update()  # 提交变量及约束，之后可直接更新 RHS。
        self.model = m, ell, bp, bq  # 同一方案的不同查询复用此非凸模型。

    def global_status(self, power, environment):
        if self.model is None:  # 不动点已能判定的方案无需创建全局求解器。
            self._build_global(environment)
        m, ell, bp, bq = self.model
        p, q = self.network.loads(power)  # 固定本次查询对应的全网节点功率。
        m.setAttr('RHS', bp, p[0])  # 更新有功平衡。
        m.setAttr('RHS', bq, q[0])  # 更新无功平衡。
        m.optimize()  # 寻找 AC 可行证书，或证明该查询不可行。
        if m.SolCount:  # 一个满足等式及限值的可行解就足以证明可行，无需等待最优性。
            current = np.array([ell[i].X for i in ell])  # 提取候选电流平方。
            P, Q, v, u = self.state(power, current)  # 通过独立支路递推重新计算，检查求解器数值误差。
            if np.max(np.abs(P*P+Q*Q-u*current)) <= 1e-7 and self.violation(P, Q, v).max() <= 1e-7:
                return 1  # 电流等式残差和所有运行限值都通过才接受证书。
        return -1 if m.Status == GRB.INFEASIBLE else 0  # 只有明确的不可行证明才返回 -1。

    def close(self):
        if self.model is not None:
            self.model[0].dispose()
