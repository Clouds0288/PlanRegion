"""逐线路规划方程、紧凑主问题、联合割子问题及独立 AC 参考。

规划割统一为 a+bᵀp+dᵀx≥0；固定网架时 x 为空，退化为负荷割。
负荷参数 p 用 kW，运行变量 P/Q/v/ell 用标幺值；流程由 main.ipynb 组织。
"""
import gurobipy as gp  # 求解混合整数规划及必要的非凸 AC 核验。
import numpy as np  # 组装网架矩阵与潮流状态。
from gurobipy import GRB  # 使用变量类型、优化方向和明确的求解状态。
import clarabel  # 求解固定整数选型后的连续 LP/SOCP。
from scipy import sparse  # 将连续子问题矩阵转换为稀疏格式。

AC_TOL = 1e-9  # AC 运行限值的标幺容差。
FIXED_POINT_TOL = 1e-12  # AC 电流等式的收敛容差。
PLANNING_TOL = 1e-8  # 规划 SP 的原始标幺约束容差，不是几何误差。


class PlanningEquations:  # 统一生成线性和 SOCP 系数；不承担预算循环或方案枚举。
    """固定径向拓扑的逐线路选型模型：c+F*p+G*y∈K，0≤y≤U0+Ux*x。

    x 是各候选线路的 one-hot 型号变量；y=(P_k,Q_k,ell_k,v)。
    各型号保留自己的送端功率和电流平方，共用真实节点电压；未选型号潮流为零。
    适用于当前 case33：非负负荷、正阻抗、有限源端 P/Q 上限；不补造线路热限。
    """

    def __init__(self, network, method, *, planning=True):  # planning=False 用同一方程检查一个固定方案。
        self.network = c = network  # 保存网架数据引用，局部别名 c 便于对应物理公式。
        self.method = method  # 在建模和证书核验中使用同一个 LP/SOCP 模式。
        groups = {o.branch: o for o in c.line_options} if planning else {}  # 固定方案没有可变型号。
        edges, r, reactance, cost, selected, self.groups = [], [], [], [], [], []  # 依次收集支路索引、型号阻抗、投资及每条线路的选型分组。
        for e in range(c.n):  # 按非根节点的入边顺序遍历在运支路。
            o = groups.get(e)  # 候选支路读取型号选项，其他支路保持已有参数。
            start = len(cost)  # 记住本支路第一个二进制变量的位置。
            types = zip(o.r, o.reactance, o.cost) if o else [(c.r[e], c.reactance[e], 0.)]  # 固定支路仅保留一个无新增费用的物理型号。
            for resistance, xline, price in types:  # 展开本支路的各个型号，阻抗和费用均为常数。
                edges.append(e)  # 记录每个型号对应的真实支路及标幺阻抗。
                r.append(resistance)  # 记录每个型号对应的真实支路及标幺阻抗。
                reactance.append(xline)  # 记录每个型号对应的真实支路及标幺阻抗。
                if o:  # 只有候选支路需要由主问题选择型号。
                    selected.append(len(edges)-1)  # 只有候选线路的型号对应二进制变量。
                    cost.append(price)  # 费用向量与二进制选型变量采用同一顺序。
            if o:  # 固定支路不建立选择等式。
                self.groups.append(np.arange(start, len(cost)))  # 每条候选线路恰选一个型号。
        self.edges, self.r, self.reactance = np.array(edges), np.array(r), np.array(reactance)  # 将逐型号物理参数转为矩阵运算所需的数组。
        self.cost, self.selected = np.array(cost), np.array(selected, dtype=int)  # 投资向量长度就是二进制变量数；空选型索引仍使用整数类型。
        t, nx, npower = len(edges), len(cost), len(c.load_nodes)  # 分别记支路型号数、选型变量数及独立负荷维数。
        self.P, self.Q, self.ell = slice(0,t), slice(t,2*t), slice(2*t,3*t)  # 连续状态 y 的前三块依次是送端有功、无功和电流平方。
        self.v = slice(3*t,3*t+c.n)  # 最后一块是各非根节点的电压平方。
        ny = 3*t+c.n  # 每个型号三个运行变量，每个节点一个电压变量。
        T = np.eye(c.n)[:, self.edges]  # 按支路汇总各型号的物理量。
        J = np.zeros((c.n,c.n))  # 构造父节点电压的索引矩阵。
        for e, parent in enumerate(c.parent):  # 遍历真实支路的送端节点。
            if parent >= 0:  # 根节点电压固定，不作为运行变量的一列。
                J[e,parent] = 1.  # J*v 为非根送端电压；根送端另加常数 1。
        root = (c.parent < 0).astype(float)  # 标出接根支路，用于固定电压常数和电源功率求和。
        self.T = T  # 保留型号到真实支路的汇总矩阵，用于重建潮流。
        balance = (np.eye(c.n)-J.T)@T  # 本支路送出量减去直接子支路送出量。
        eq = np.zeros((3*c.n, ny))  # 三组节点等式分别为有功、无功平衡和压降。
        eq[:c.n,self.P], eq[:c.n,self.ell] = balance, -T*self.r  # P-ΣP_child-r*ell=p_node。
        eq[c.n:2*c.n,self.Q], eq[c.n:2*c.n,self.ell] = balance, -T*self.reactance  # 无功平衡：Q−ΣQ_child−χ·ell=q_node。
        eq[2*c.n:3*c.n,self.v] = J-np.eye(c.n)  # 电压项为送端电压减受端电压。
        eq[2*c.n:3*c.n,self.P], eq[2*c.n:3*c.n,self.Q] = -2*T*self.r, -2*T*self.reactance  # 压降的一次功率项为 −2(rP+χQ)。
        eq[2*c.n:3*c.n,self.ell] = T*(self.r**2+self.reactance**2)  # 完整支路压降的损耗项。
        eqc = np.r_[-c.fixed_p/c.base, -c.fixed_q/c.base, root]  # 固定负荷移到等式左侧，根送端电压贡献常数 1。
        eqf = np.vstack([-c.E/c.base, -c.E*c.q_ratio/c.base, np.zeros((c.n,npower))])  # E 将独立负荷嵌入节点；无功按给定 Q/P 比例变化。
        self.equal_count = len(eqc)  # 记录原始等式数，供直接模型避免重复建模。
        # 原模型保留等式；phase I 使用正负两行，把等式残差也纳入 eta。
        self.c, self.F = np.r_[eqc,-eqc], np.vstack([eqf,-eqf])  # 等式 h=0 在 phase I 中写为 h≥−eta 与 −h≥−eta。
        self.G = np.vstack([eq,-eq])  # 状态系数采用同样的正负复制，保持行顺序一致。
        pmax = np.minimum(c.capacity[self.edges], min(c.source_pmax,c.source_smax))  # 非负纯负荷条件下，支路有功不超过支路及源端有效上界。
        qmax = np.full(t,min(c.source_qmax,c.source_smax))  # 支路无功不超过源端无功或视在容量上界。
        lmax = pmax/self.r if method == 'socp' else np.zeros(t)  # P≥r*ell 给出有效电流平方上界；线性模型令 ell=0。
        self.upper = np.r_[pmax,qmax,lmax,c.vmax]  # 所有 y 的物理有效界，用于型号联接和对偶残差补偿。
        self.box_constant = self.upper.copy()  # 变量盒的常数项先采用物理上界，候选型号随后改为与 x 联接。
        self.box_selection = np.zeros((ny,nx))  # U(x)=U0+Ux·x，矩阵列与各型号的二进制变量对应。
        self.linked = np.concatenate([offset+self.selected for offset in (0,t,2*t)])  # 只对候选型号的 P/Q/ell 建立联接，不复制节点电压。
        for offset, bound in ((0,pmax),(t,qmax),(2*t,lmax)):  # 分别设置三类型号运行变量的上界。
            self.box_constant[offset+self.selected] = 0.  # 候选型号没有与 x 无关的正上界，未选中时必须为零。
            self.box_selection[offset+self.selected,np.arange(nx)] = bound[self.selected]  # 型号相关上界 U(x)；未选型号的潮流严格为零。
        linear_c, linear_g = [], []  # 收集等式之外的线性运行限值。
        g = np.zeros((c.n,ny))  # 选取电压变量，构造 v−vmin≥0。
        g[:,self.v] = np.eye(c.n)  # 选取电压变量，构造 v−vmin≥0。
        linear_c.extend(-c.vmin)  # 将电压下限的常数项和系数加入约束表。
        linear_g.extend(g)  # 将电压下限的常数项和系数加入约束表。
        for block, limit in ((self.P,c.source_pmax),(self.Q,c.source_qmax)):  # 分别处理电源总有功和总无功上限。
            if np.isfinite(limit):  # 未给定的无限上限不需要建立约束行。
                g = np.zeros(ny)  # 只累加接根型号的送端功率，写成 limit−Σpower≥0。
                g[block] = -root[self.edges]  # 只累加接根型号的送端功率，写成 limit−Σpower≥0。
                linear_c.append(limit)  # 源端功率包含网损。
                linear_g.append(g)  # 源端功率包含网损。
        self.c = np.r_[self.c,linear_c]  # 在等式的正负两组之后追加运行限值常数。
        self.F = np.vstack([self.F,np.zeros((len(linear_c),npower))])  # 这些限值直接作用于状态，对独立负荷没有额外显式系数。
        self.G = np.vstack([self.G,linear_g])  # 按同一行顺序追加线性状态系数。
        self.linear_count, self.cones = len(self.c), []  # 非负锥覆盖前 linear_count 行，其后才是二阶锥。
        if method == 'socp':  # 线性模式令电流为零，并省略全部二阶锥。
            for k in range(t):  # SOCP 对每个支路型号建立一个四维锥。
                g = np.zeros((4,ny))  # 四行依次对应 u+ell、2P、2Q、u−ell。
                e = self.edges[k]  # 查到该型号所属的真实支路。
                g[0,self.v] = g[3,self.v] = J[e]  # 首、末分量均含真实送端电压，不能使用型号虚拟电压。
                g[0,self.ell.start+k] = 1.  # 锥首分量含正的电流平方。
                g[1,self.P.start+k] = g[2,self.Q.start+k] = 2.  # 两个方向分量分别为两倍有功和无功。
                g[3,self.ell.start+k] = -1.  # 最后一个方向分量含负的电流平方。
                self._cone([root[e],0.,0.,root[e]],g)  # (u_e+ell_k,2P_k,2Q_k,u_e-ell_k)∈SOC；u_e 是真实送端电压平方。
            if np.isfinite(c.source_smax):  # 仅在网架给定视在上限时建立源端容量圆。
                g = np.zeros((3,ny))  # 三维锥的两个方向分量是源端总 P 和总 Q。
                g[1,self.P], g[2,self.Q] = root[self.edges],root[self.edges]  # 三维锥的两个方向分量是源端总 P 和总 Q。
                self._cone([c.source_smax,0.,0.],g)  # 容量约束写成 (Smax,ΣP,ΣQ)∈SOC。
        self.relax = np.zeros(len(self.c))  # eta 松弛所有线性行，锥方向分量初始保持零。
        self.relax[:self.linear_count] = 1.  # eta 松弛所有线性行，锥方向分量初始保持零。
        for s in self.cones:  # 逐锥指定统一 phase I 的松弛方向。
            self.relax[s.start] = 1.  # 锥只松弛首分量，不更改锥的方向分量。

    def _cone(self, constant, matrix):  # 向统一系数表追加一个二阶锥块。
        self.cones.append(slice(len(self.c),len(self.c)+len(constant)))  # 记录该锥的行切片，供求解与证书检验共同使用。
        self.c = np.r_[self.c,constant]  # 追加锥常数和状态系数。
        self.G = np.vstack([self.G,matrix])  # 追加锥常数和状态系数。
        self.F = np.vstack([self.F,np.zeros((len(constant),self.F.shape[1]))])  # 锥通过运行状态依赖负荷，因此显式负荷系数为零。

    def selection(self, choice):  # 把逐线路型号编号转换为主问题 one-hot 编码。
        """线路型号编号转 one-hot；只用于实例化方案和枚举对照。"""
        x = np.zeros(len(self.cost))  # 每个候选型号对应一个二进制位置。
        for group,k in zip(self.groups,choice):  # 各条线路的型号编号与其变量分组一一对应。
            x[group[int(k)]] = 1.  # 每组只激活被选型号。
        return x  # 返回选型向量，固定网架时长度为零。

    def choice(self, x):  # 将求解后的 one-hot 向量还原成逐线路型号编号。
        return np.array([np.argmax(x[g]) for g in self.groups])  # 每组最大分量的位置就是该线路所选型号。

    def margin(self, x, power, state):  # 检查不含 eta 的原始约束最小余量。
        g = self.c+self.F@power+self.G@state  # 将实际负荷和运行状态代回统一物理方程。
        values = [g[:self.linear_count].min(),np.min(state),  # 同时检查线性行和运行变量非负性。
                  np.min(self.box_constant+self.box_selection@x-state)]  # 检查包含型号联接关系的变量上界。
        values.extend(g[s.start]-np.linalg.norm(g[s.start+1:s.stop]) for s in self.cones)  # 每个 SOC 的余量为首分量减尾部二范数。
        return min(values)  # 不含 eta 的原约束余量；等式正负两行等价于检查绝对残差。

    def restore(self, x, power, state):  # 按物理等式重建候选状态，再交给 margin 检验。
        """用选中型号的电流重建等式精确的状态，消除未选型号退化锥中的浮点噪声。"""
        c = self.network  # 只读取该实例共享的网架配置。
        active = np.ones(len(self.edges))  # 固定支路始终激活，候选型号由整数 x 决定。
        active[self.selected] = x  # 固定支路始终激活，候选型号由整数 x 决定。
        y = np.zeros_like(state)  # 未激活型号从零状态开始，避免残留数值噪声。
        ell = np.clip(state[self.ell],0.,self.upper[self.ell])*active  # 将电流限制在有效盒中，未选型号电流清零。
        p,q = c.loads(power)  # 取得含固定背景的完整节点 P/Q 负荷。
        ell[(c.D@(p[0]+q[0]))[self.edges]==0.] = 0.  # 纯负荷空子树不保留锥松弛产生的无效循环损耗。
        P = c.D@(p[0]+self.T@(self.r*ell))  # 汇总节点负荷和所有下游有功损耗。
        Q = c.D@(q[0]+self.T@(self.reactance*ell))  # 汇总节点无功负荷和所有下游无功损耗。
        y[self.ell], y[self.P], y[self.Q] = ell, P[self.edges]*active, Q[self.edges]*active  # 将真实支路功率填入选中型号，其余型号仍为零。
        v = 1.-2*c.D.T@(self.T@(self.r*y[self.P]+self.reactance*y[self.Q]))  # 沿根到节点路径累加一次压降。
        v += c.D.T@(self.T@((self.r**2+self.reactance**2)*ell))  # 补入完整支路压降中的二次电流项。
        y[self.v] = v  # 将重建的节点电压平方写回状态向量。
        return y  # 是否可行仍由 margin 对原始限值与锥逐一核验。


class PlanningModel:  # 负责单次规划查询；多次查询和切割循环由 Notebook 组织。
    """逐线路 MILP/MISOCP；cuts_only=True 时为包含 x、p 的联合割主问题。"""

    def __init__(self, equations, *, power=None, budget=np.inf, direction=None, cuts_only=False):  # 多个查询复用同一网架方程。
        self.equations = e = equations  # MP 和 SP 共享固定系数，不在每个查询中重复推导。
        network, method = e.network, e.method  # 查询只改变负荷、预算和目标。
        self.model = m = gp.Model('compact_'+method)  # 为本次投资或边界查询创建优化模型。
        m.Params.OutputFlag = 0  # 关闭日志，统一使用单线程以便比较耗时。
        m.Params.Threads = 1  # 关闭日志，统一使用单线程以便比较耗时。
        m.Params.MIPGap = m.Params.MIPGapAbs = 0.  # 不以相对或绝对整数间隙提前接受次优投资。
        m.Params.FeasibilityTol = m.Params.IntFeasTol = m.Params.OptimalityTol = 1e-9  # 分别设置原始约束、整数性和最优性容差。
        m.Params.BarQCPConvTol = 1e-9  # 设置二阶锥连续子问题的内点法收敛精度。
        m.Params.DualReductions = 0  # 区分不可行与无界，避免模糊状态被当作不可行证明。
        m.Params.NonConvex = 0  # 若模型误写成非凸二次式，立即暴露，不能悄悄改用非凸求解。
        self.x = m.addMVar(len(e.cost),vtype=GRB.BINARY,name='line_type')  # 每个候选型号一个二进制变量，固定网架时长度为零。
        self.power = p = m.addMVar(len(network.load_nodes),lb=0.,name='p_kw')  # 独立负荷非负，单位为 kW。
        for group in e.groups:  # 按候选线路分别建立 one-hot 等式。
            m.addConstr(self.x[group].sum()==1.)  # 每条线路恰好保留或升级为一个型号。
        m.addConstr(p.sum()/network.base<=network.power_limit/network.base)  # 把无损总负荷外界转为标幺尺度以改善数值条件。
        if np.isfinite(budget):  # 无限预算无需建立额外投资约束。
            m.addConstr(e.cost@self.x<=budget)  # 所有线路型号的增量投资之和不得超过预算。
        if power is not None:  # 固定负荷查询求能够承载该点的最低投资。
            m.addConstr(p==np.asarray(power))  # 锁定查询负荷，并以型号投资作为最小化目标。
            m.setObjective(e.cost@self.x,GRB.MINIMIZE)  # 锁定查询负荷，并以型号投资作为最小化目标。
            self.objective_scale = 1.  # 投资目标本来就是原单位，返回时无需缩放。
        else:  # 负荷未固定时执行预算内的负荷最大化查询。
            if direction is not None:  # 给定方向时只沿该方向寻找边界。
                scale = m.addMVar(1,lb=0.,name='radial_total_kw')  # 非负射线半径表示三个独立节点的总负荷。
                m.addConstr(p==(np.asarray(direction)/np.sum(direction))[:,None]@scale)  # 射线比例固定，背景负荷不缩放。
            m.setObjective(p.sum()/network.base,GRB.MAXIMIZE)  # 目标用标幺量改善锥求解精度；返回时恢复 kW。
            self.objective_scale = network.base  # 将标幺目标界恢复为 kW。
        self.state = None  # 仅含联合割的 MP 没有运行状态变量。
        if not cuts_only:  # 直接 MILP/MISOCP 模式才显式建立完整运行模型。
            self.state = y = m.addMVar(len(e.upper),lb=0.,ub=e.upper,name='flow')  # 所有运行变量非负，并使用物理推导的有限上界。
            g = e.c+e.F@p+e.G@y  # 代入同一 PlanningEquations 生成 Gurobi 约束表达式。
            m.addConstr(g[:e.equal_count]==0.)  # 原模型等式只建一次，负号副本仅用于 phase I。
            m.addConstr(g[2*e.equal_count:e.linear_count]>=0.)  # 跳过等式负号副本，仅添加后续真实运行不等式。
            m.addConstr(y[e.linked]<=e.box_selection[e.linked]@self.x)  # 0≤y_ek≤U_ek*x_ek，未选型号为零。
            if method == 'socp':  # 线性模式不建立电流锥；ell 已由零上界固定。
                for k in range(len(e.edges)):  # 逐型号建立真实送端电压下的旋转锥。
                    P,Q,l = (y[s.start+k].item() for s in (e.P,e.Q,e.ell))  # 取得当前型号的三个标量变量。
                    parent = network.parent[e.edges[k]]  # 根据真实支路查找送端节点。
                    u = 1. if parent < 0 else y[e.v.start+parent].item()  # 接根支路取固定电压 1，其余读取父节点电压变量。
                    m.addQConstr(P*P+Q*Q<=u*l)  # u、ell 非负，这是凸旋转二阶锥。
                if np.isfinite(network.source_smax):  # 仅检查网架实际给定的源端视在容量。
                    roots = np.isin(e.edges,network.roots)  # 选出所有接根支路及其型号。
                    ps,qs = y[e.P][roots].sum().item(),y[e.Q][roots].sum().item()  # 源端 P/Q 是接根支路送端功率之和。
                    m.addQConstr(ps*ps+qs*qs<=network.source_smax**2)  # 用凸二次不等式表示源端容量圆。

    def add_cut(self, cut):  # 向主问题追加一条同时含 p 和 x 的有效联合割。
        n = len(self.equations.network.load_nodes)  # 用负荷维数确定割中 p 与 x 系数的分界。
        self.model.addConstr(cut[0]+cut[1:1+n]@self.power+cut[1+n:]@self.x>=0.)  # 联合割对全部合法选型有效。

    def solve(self):  # 完成单次全局规划查询，并返回候选与目标界。
        self.model.optimize()  # 调用 Gurobi 求解当前直接模型或含历史割的 MP。
        if self.model.Status == GRB.INFEASIBLE:  # 只有明确不可行状态才能返回无可行方案。
            return None  # 用 None 表示该查询已被求解器证明不可行。
        if self.model.Status != GRB.OPTIMAL:  # 其他终止状态不能作为精确投资或边界查询的答案。
            raise RuntimeError(f'Compact planning status {self.model.Status}')  # 保留明确的求解失败原因，不伪造最优性证书。
        x = np.rint(self.x.X).astype(int)  # 把满足整数容差的 x 还原为实际离散建设方案。
        objective = self.equations.cost@x if self.model.ModelSense==GRB.MINIMIZE else self.power.X.sum()  # 投资按实际型号计算，最大负荷按 kW 求和。
        # 投资按实际整数选型计价，避免 1+1e-10 被误判为超出预算 1。
        return dict(x=x,p=self.power.X,  # 只返回本次查询的选型和负荷，不重复保存网架数据。
                    objective=float(objective),bound=self.model.ObjBound*self.objective_scale,  # 同时返回可行目标与求解器给出的全局目标界。
                    state=None if self.state is None else self.state.X)  # 联合割 MP 没有潮流状态，直接模型返回对应运行解。


class PlanningSP:  # 同一个连续子问题同时支持 LP 和 SOCP 认证及联合割。
    """固定整数 x,p 的连续 phase I；盒上界随 x 变化，割为 a+bᵀp+dᵀx≥0。"""

    def __init__(self, equations):  # 绑定可复用的方程，并设置连续求解精度。
        self.equations = equations  # 不重新生成网架、费用或物理矩阵。
        self.calls = 0  # 只统计实际调用连续求解器的次数。
        self.settings = settings = clarabel.DefaultSettings()  # 使用 Clarabel 的标准锥优化设置。
        settings.verbose = False  # 关闭求解日志，固定为单线程。
        settings.max_threads = 1  # 关闭求解日志，固定为单线程。
        settings.tol_gap_abs = settings.tol_gap_rel = settings.tol_feas = 1e-12  # 控制辅助问题的间隙与原始可行性误差。
        settings.reduced_tol_feas = 1e-11  # 数值停滞时的较宽松终止精度不替代原约束检验。

    def solve(self, x, power):  # 固定一个整数选型和负荷点，返回证书、有效割或未确定结果。
        e = self.equations  # 读取共享的固定系数及变量切片。
        upper = e.box_constant+e.box_selection@x  # 选型只改变变量盒的右端项，物理方程系数不变。
        columns = np.flatnonzero(upper>0.)  # 未选型号及 LP 的 ell 恒为零，直接消去，不求解退化锥。
        active = np.ones(len(e.edges),dtype=bool)  # 标记保留的支路型号，以去掉未选型号的退化锥。
        active[e.selected] = np.asarray(x,dtype=bool)  # 标记保留的支路型号，以去掉未选型号的退化锥。
        cones = [s for k,s in enumerate(e.cones) if k>=len(active) or active[k]]  # 保留活动支路锥，以及可选的源端视在容量锥。
        rows = np.r_[np.arange(e.linear_count),*[np.arange(s.start,s.stop) for s in cones]].astype(int)  # 保留全部线性行，并按原顺序提取活动锥行。
        n = len(columns)  # 消去恒零变量后的连续状态维数。
        G = e.G[np.ix_(rows,columns)]  # 从完整 G 中取本次实际求解所需的行和列。
        matrix = np.vstack([-np.c_[G,e.relax[rows]],  # Clarabel 使用 A·z+s=b，因此物理系数整体取负。
                            np.c_[-np.eye(n),np.zeros(n)],np.c_[np.eye(n),np.zeros(n)],  # 非负锥分别表示运行变量的下界和上界。
                            np.r_[np.zeros(n),-1.]])  # 物理松弛、y≥0、y≤U(x)、eta≥0。
        rhs = np.r_[(e.c+e.F@power)[rows],np.zeros(n),upper[columns],0.]  # 固定 p 和 x 后组装约束右端项。
        kinds = [clarabel.NonnegativeConeT(e.linear_count)]  # 等式正负行及运行限值都属于非负锥。
        kinds += [clarabel.SecondOrderConeT(s.stop-s.start) for s in cones]  # 按同一行次序声明各个二阶锥的维数。
        kinds += [clarabel.NonnegativeConeT(2*n+1)]  # 变量盒和 eta 非负性合并为最后一块非负锥。
        solver = clarabel.DefaultSolver(sparse.csc_matrix((n+1,n+1)),np.r_[np.zeros(n),1.],  # 目标仅最小化最后一个变量 eta，二次目标矩阵为零。
            sparse.csc_matrix(matrix),rhs,kinds,self.settings)  # 将稠密组装矩阵转为求解器所需的压缩列格式。
        result = solver.solve()  # 实际求解一次连续 phase I。
        self.calls += 1  # 不把绘图点、几何检查或内点插值计成 SP 求解。
        if not np.isfinite(np.r_[result.x,result.z]).all():  # 非有限的原始解或乘子不能用于任何证书。
            raise RuntimeError(f'Planning phase I status {result.status}')  # 非有限结果无法形成可审核的物理或对偶证书。
        state = np.zeros(len(e.upper))  # 将消元后的解回填到完整状态，恒零变量仍为零。
        state[columns] = result.x[:-1]  # 将消元后的解回填到完整状态，恒零变量仍为零。
        state = e.restore(x,power,state)  # 用电流重建功率平衡与压降等式，再检查原始约束。
        eta = max(0.,result.x[-1])  # 去掉 eta 的微小负舍入值，保留实际辅助违反量。
        if e.margin(x,power,state)>=-PLANNING_TOL:  # 实际约束余量通过才认证，不能仅凭 eta 很小。
            return dict(eta=eta,cut=None,state=state,feasible=True)  # 原始约束通过才认证原查询。
        # 未收敛时仍可检验候选证书；有效性由原约束或下方的对偶锥与盒补偿保证。
        dual = np.zeros(len(e.c))  # 被消去的锥乘子补零，属于其对偶锥。
        dual[rows] = result.z[:len(rows)]  # 被消去的锥乘子补零，属于其对偶锥。
        dual[:e.linear_count] = np.maximum(dual[:e.linear_count],0.)  # 将线性约束乘子投影到其非负对偶锥。
        for s in e.cones:  # 逐个修正二阶锥乘子，未参与求解的锥仍为零。
            dual[s.start] = max(dual[s.start],np.linalg.norm(dual[s.start+1:s.stop]))  # SOC 自对偶；首分量至少为尾部范数才有全域有效性。
        residual = np.maximum(e.G.T@dual,0.)  # max_{0≤y≤U(x)} (Gᵀλ)ᵀy = (Gᵀλ)_+ᵀU(x)。
        cut = np.r_[dual@e.c+residual@e.box_constant+1e-12,  # 常数项包含对偶物理项、完整变量盒补偿和浮点余量。
                    e.F.T@dual,e.box_selection.T@residual]  # 盒支撑函数给出 x 系数，同时补偿全部驻点残差。
        cut /= max(np.max(np.abs(cut)),np.max(np.abs(cut[1:1+len(power)]))*e.network.base)  # 固定方案 x 为空时仍使用同一个正比例归一化。
        if not (cut[0]+cut[1:1+len(power)]@power+cut[1+len(power):]@x < -1e-9):  # 只在割严格排除当前 (p,x) 时交给外层使用。
            cut = None  # 浮点边界可能尚未获证且无可靠割；外层必须显式处理未确定状态。
        return dict(eta=eta,cut=cut,state=state,feasible=False)  # 失败查询的状态可用于同方案内点插值，但不能认证查询点。


class ACPowerFlow:  # 通过另一套支路递推和完整电流等式建立 AC 参考。
    """独立完整 AC；仅读取网架数据，不使用 LP/SOCP 的方程矩阵或乘子。

    适用于非负 P/Q 负荷、正阻抗、根电压固定为 1 p.u. 的径向网络。
    不读取 PlanningEquations 的矩阵；以支路递推独立实现 AC 电流等式。
    """

    def __init__(self, network):  # 保存网架引用，按需创建非凸核验模型。
        self.network = network  # AC 与 LP/SOCP 只共享网架及运行限值数据。
        self.model = None  # 非凸模型仅在不动点无法判定时创建。

    def state(self, power, ell):  # 为给定负荷及电流重建独立 AC 状态。
        p, q = self.network.loads(power)  # 补入固定背景负荷，并将 kW/kvar 换成标幺量。
        return self._state(p, q, ell)  # 可一次重建多组负荷状态。

    def _state(self, p, q, ell):  # 内部批量递推使用已换算为标幺值的完整节点负荷。
        c = self.network  # 读取独立 AC 所需的网架参数。
        ell = np.asarray(ell).reshape(-1, c.n)  # 行是查询样本，列是支路电流平方。
        P, Q = p+ell*c.r, q+ell*c.reactance  # 本节点负荷加本支路有功/无功损耗。
        for i in reversed(c.order):  # 从叶到根，先汇总子树功率。
            if c.parent[i] >= 0:  # 接根支路不再向其他非根支路汇总。
                P[:, c.parent[i]] += P[:, i]  # P_parent=p_parent+r_parent*ell_parent+ΣP_child。
                Q[:, c.parent[i]] += Q[:, i]  # 无功平衡使用同样的子树汇总。
        v, u = np.ones_like(P), np.ones_like(P)  # 根电压平方为 1。
        for i in c.order:  # 从根到叶，根据上一级电压逐支路计算压降。
            if c.parent[i] >= 0:  # 接根支路继续使用固定根电压。
                u[:, i] = v[:, c.parent[i]]  # 当前支路送端电压等于父节点电压。
            # 完整支路压降：v=u-2(rP+xQ)+(r²+x²)*ell，保留二次电流项。
            v[:, i] = (u[:, i]-2*(c.r[i]*P[:, i]+c.reactance[i]*Q[:, i])  # 从送端电压减去一次支路压降。
                       +(c.r[i]**2+c.reactance[i]**2)*ell[:, i])  # 加入完整 AC 压降中的电流平方修正。
        return P, Q, v, u  # 返回送端功率、受端电压平方及送端电压平方。

    def violation(self, P, Q, v):  # 计算每个样本对全部运行限值的最大违反量。
        c = self.network  # 读取支路和电源限值。
        ps, qs = P[:, c.roots].sum(axis=1), Q[:, c.roots].sum(axis=1)  # 含网损的电源送出功率。
        # 所有项目均写为“违反量”：≤0 才满足全部限值，未启用上界产生 -inf。
        return np.maximum.reduce([  # 每个样本取所有运行约束的最大违反量。
            np.max(c.vmin-v, axis=1), np.max(v-c.vmax, axis=1),  # 电压平方下限和上限。
            np.max(P-c.capacity, axis=1), ps-c.source_pmax, qs-c.source_qmax,  # 支路及电源 P/Q 限值。
            np.hypot(ps, qs)-c.source_smax])  # 电源视在功率上限。

    def classify(self, power, return_currents=False):  # 独立判定 AC 可行性，并可返回电流证书。
        """1 可行，-1 已证不可行，0 未确定；不把迭代失败当作不可行。"""
        c = self.network  # 加载同一物理网架的纯负荷参数。
        power = np.asarray(power).reshape(-1, len(c.load_nodes))  # 每个样本包含三个独立负荷坐标。
        p, q = c.loads(power)  # 背景负荷只组装一次，所有不动点迭代复用。
        ell = np.zeros((len(power), c.n))  # 从零电流开始，构造单调递增的电流平方序列。
        status = np.zeros(len(power), dtype=np.int8)  # 初始均为未确定；0 绝不能直接算作不可行。
        active = np.arange(len(power))  # 后续仅更新尚未得到证书的样本。
        for _ in range(160):  # 达到迭代上限仍未获证的点保持未确定。
            if not len(active):  # 所有样本均已获得可行或不可行证书。
                break  # 没有待判定样本时结束迭代。
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

    def scan(self, points):  # 分批扫描评价网格，限制矩阵内存占用。
        status = np.empty(len(points), dtype=np.int8)  # 与输入网格点一一对应的三态判定。
        for start in range(0, len(points), 8192):  # 分批限制电流矩阵的内存，判定方程不变。
            status[start:start+8192] = self.classify(points[start:start+8192])  # 每批使用完全相同的三态 AC 判据。
        return status  # 按原网格点顺序返回状态。

    def _build_global(self, environment):  # 为不动点未判定的情况建立独立非凸 AC 核验。
        """显式非凸 AC 等式模型，用于未确定点和独立交叉核验。"""
        c = self.network  # 读取网架原始参数，不导入规划方程矩阵。
        m = gp.Model('independent_AC_reference', env=environment)  # 在共享环境中创建本方案的 AC 等式模型。
        m.Params.OutputFlag = 0  # 关闭逐点核验日志。
        m.Params.Threads = 1  # 固定线程数，避免不同方案计时口径变化。
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
        for i in range(c.n):  # 对每条入边及其受端节点建立平衡与压降方程。
            up = 1. if c.parent[i] < 0 else v[int(c.parent[i])]  # 根节点电压固定，其余取父节点变量。
            # P_i-ΣP_child-r_i*ell_i=p_i；Q_i-ΣQ_child-x_i*ell_i=q_i。
            bp.append(m.addConstr(P[i]-gp.quicksum(P[int(j)] for j in c.children[i])-c.r[i]*ell[i] == 0))  # 查询时 RHS 替换为 p_i。
            bq.append(m.addConstr(Q[i]-gp.quicksum(Q[int(j)] for j in c.children[i])-c.reactance[i]*ell[i] == 0))  # 查询时 RHS 替换为 q_i。
            m.addConstr(v[i] == up-2*(c.r[i]*P[i]+c.reactance[i]*Q[i])  # 完整 AC 支路压降。
                        +(c.r[i]**2+c.reactance[i]**2)*ell[i])  # 补入本支路的二次电流项。
            m.addQConstr(P[i]*P[i]+Q[i]*Q[i] == up*ell[i])  # AC 必须保留等号。
        ps, qs = (gp.quicksum(values[int(i)] for i in c.roots) for values in (P, Q))  # 电源送出功率包含网损。
        if np.isfinite(c.source_pmax):  # 只有给定有功上限时才建立此约束。
            m.addConstr(ps <= c.source_pmax)  # 电源有功上限。
        if np.isfinite(c.source_qmax):  # 只有给定无功上限时才建立此约束。
            m.addConstr(qs <= c.source_qmax)  # 电源无功上限。
        if np.isfinite(c.source_smax):  # 未给定变压器容量时不补造视在约束。
            m.addQConstr(ps*ps+qs*qs <= c.source_smax**2)  # 电源视在容量圆。
        m.setObjective(0.)  # 这里只判可行性，不额外改变运行目标。
        m.update()  # 提交变量及约束，之后可直接更新 RHS。
        self.model = m, ell, bp, bq  # 同一方案的不同查询复用此非凸模型。

    def global_status(self, power, environment):  # 通过非凸求解器补判一个 AC 查询点。
        if self.model is None:  # 不动点已能判定的方案无需创建全局求解器。
            self._build_global(environment)  # 首次遇到未确定点时才付出建模成本。
        m, ell, bp, bq = self.model  # 复用同一方案的电流变量和节点平衡行。
        p, q = self.network.loads(power)  # 固定本次查询对应的全网节点功率。
        m.setAttr('RHS', bp, p[0])  # 更新有功平衡。
        m.setAttr('RHS', bq, q[0])  # 更新无功平衡。
        m.optimize()  # 寻找 AC 可行证书，或证明该查询不可行。
        if m.SolCount:  # 一个满足等式及限值的可行解就足以证明可行，无需等待最优性。
            current = np.array([ell[i].X for i in ell])  # 提取候选电流平方。
            P, Q, v, u = self.state(power, current)  # 通过独立支路递推重新计算，检查求解器数值误差。
            if np.max(np.abs(P*P+Q*Q-u*current)) <= 1e-7 and self.violation(P, Q, v).max() <= 1e-7:  # 复核独立重建状态的等式残差及全部运行限值。
                return 1  # 电流等式残差和所有运行限值都通过才接受证书。
        return -1 if m.Status == GRB.INFEASIBLE else 0  # 只有明确的不可行证明才返回 -1。

    def close(self):  # 释放按需创建的非凸优化模型。
        if self.model is not None:  # 仅不动点核验时没有需要释放的求解器。
            self.model[0].dispose()  # 释放 Gurobi 模型及关联资源。
