"""逐线路规划方程、紧凑主问题、联合割子问题及全局残余搜索。

规划割统一为 α+βᵀp+δᵀx≥0；固定网架时 x 为空，退化为负荷割。
负荷参数 p 用 kW，运行变量 P/Q/v/ell 用标幺值；流程由 main.py 组织。
"""
from time import perf_counter

import clarabel
import gurobipy as gp
import numpy as np
from gurobipy import GRB
from scipy import sparse

DEFAULT_SOLVER_THREADS = 20
MP_TIME_LIMIT = 20.
SP_TIME_LIMIT = {'linear': 5., 'socp': 10.}
RESIDUAL_TIME_LIMIT = 120.
PLANNING_TOL = 1e-8  # 原始标幺约束容差，与几何误差分别控制。


class PlanningEquations:
    """固定参考方向：A*x+B*y+C*p≼_K b，lb(x)≤y≤ub(x)。

    锥余量为 b-A*x-B*y-C*p∈K；线性行的 ≼_K 即普通 ≤。
    A/B/C 分别乘选型 x、运行状态 y、负荷 p；b 为右端常数向量。
    lb(x)=lb_const+lb_x*x，ub(x)=ub_const+ub_x*x。
    x 只选择走廊型号；P/Q 的正负表示流向，y=(P_k,Q_k,ell_k,v,switch)。
    switch 仅在可重构网架中包含正负压降开断余量；p 用 kW，物理状态用标幺值。
    各型号保留自己的送端功率和电流平方，共用真实节点电压；未选型号潮流为零。
    非负负荷、正阻抗和有限源端上界；可选走廊另含压降开断余量，仍用同一种联合割。
    """

    def __init__(self, network, method):
        self.network = network
        self.method = method
        self.corridors = network.corridors
        self.reconfigurable = any(not corridor.must_use for corridor in self.corridors)
        self._index_network()
        self._build_variable_bounds()
        self._assemble_equations()

    def _index_network(self):
        """按走廊、型号顺序展开参数，确定 x 与 y 的列顺序。"""
        net = self.network
        node_index = {node: i for i, node in enumerate(net.nodes)}
        node_index[net.root] = -1  # 根节点电压恒为 1，不进入电压变量。
        type_corridor, decision_types = [], []
        senders, receivers = [], []
        r, reactance, cost = [], [], []
        self.type_indices, self.variable_keys = {}, []
        self.constant_cost = 0.
        for e, corridor in enumerate(self.corridors):
            a, b = (node_index[n] for n in corridor.endpoints)
            sender, receiver = (b, a) if b < 0 else (a, b)  # 接根线路统一以根为参考送端。
            indices = {}
            for line_type in corridor.types:
                k = len(type_corridor)
                indices[line_type.id] = k
                type_corridor.append(e)
                senders.append(sender)
                receivers.append(receiver)
                r.append(line_type.r)
                reactance.append(line_type.reactance)
                if not corridor.must_use or len(corridor.types)>1:
                    self.variable_keys.append((corridor.id, line_type.id))
                    decision_types.append(k)
                    cost.append(line_type.investment_cost)
                else:
                    self.constant_cost += line_type.investment_cost  # 固定单型号不占用 x。
            self.type_indices[corridor.id] = indices
        # type_indices['01']['M'] 给出物理型号位置；variable_keys 给出 x 的具名键。
        self.type_corridor = np.array(type_corridor, dtype=int)
        self.decision_types = np.array(decision_types, dtype=int)
        self.senders, self.receivers = np.array(senders), np.array(receivers)
        self.r, self.reactance, self.cost = np.array(r), np.array(reactance), np.array(cost)
        t = len(type_corridor)
        self.type_selection = np.eye(t)[:,self.decision_types]  # z = 1 + S@(x-1)：固定型号恒为 1。

        # y = (P[0:t], Q[0:t], ell[0:t], v[0:n], s_plus[0:t], s_minus[0:t])。
        self.P_slice, self.Q_slice, self.ell_slice = slice(0,t), slice(t,2*t), slice(2*t,3*t)
        self.v_slice = slice(3*t,3*t+net.n)
        self.switch_slice = slice(3*t+net.n,5*t+net.n) if self.reconfigurable else slice(3*t+net.n,3*t+net.n)

    def _build_variable_bounds(self):
        """全局变量界及 lb(x)=lb_const+lb_x@x、ub(x)=ub_const+ub_x@x。"""
        net = self.network
        t, nx, ny = len(self.type_corridor), len(self.variable_keys), self.switch_slice.stop
        capacity = np.array([line.capacity for corridor in self.corridors for line in corridor.types])

        # 物理上限：Pmax=min(线路容量,源端P上限,源端S上限)，Qmax=min(源端Q上限,源端S上限)。
        pmax = np.minimum(capacity,min(net.source_pmax,net.source_smax))
        qmax = np.full(t,min(net.source_qmax,net.source_smax))
        lmax = pmax/self.r if self.method == 'socp' else np.zeros(t)  # r*ell<=Pmax；线性模型令 ell=0。
        self.y_ub_global = np.r_[pmax,qmax,lmax,net.vmax]
        if self.reconfigurable:
            # |u-v|<=M；根节点 u=1，其余参考送端使用节点电压上下限。
            umax = np.where(self.senders<0,1.,net.vmax[self.senders])
            umin = np.where(self.senders<0,1.,net.vmin[self.senders])
            switch_bound = np.maximum(umax-net.vmin[self.receivers],net.vmax[self.receivers]-umin)
            self.y_ub_global = np.r_[self.y_ub_global,switch_bound,switch_bound]

        # 全局外包界：可反向的 P/Q 允许负值；ell、v、s 非负，v>=vmin 另列为运行约束。
        self.y_lb_global = np.zeros(ny)
        reversible = np.flatnonzero((self.senders >= 0) & self.reconfigurable)
        self.y_lb_global[reversible], self.y_lb_global[t+reversible] = -pmax[reversible], -qmax[reversible]
        self.ub_const = self.y_ub_global.copy()
        self.ub_x = np.zeros((ny,nx))
        self.lb_const = self.y_lb_global.copy()
        self.lb_x = np.zeros((ny,nx))

        # 选型联接：Pmin*z<=P<=Pmax*z，Qmin*z<=Q<=Qmax*z，0<=ell<=ell_max*z。
        for block in (self.P_slice,self.Q_slice,self.ell_slice):
            rows = block.start+self.decision_types
            self.ub_const[rows] = self.lb_const[rows] = 0.
            self.ub_x[rows,np.arange(nx)] = self.y_ub_global[rows]
            self.lb_x[rows,np.arange(nx)] = self.y_lb_global[rows]
        if self.reconfigurable:
            # 开断余量：0<=s_plus,s_minus<=M*(1-z)，投入型号的压降余量恒为零。
            for offset in (self.switch_slice.start,self.switch_slice.start+t):
                self.ub_const[offset:offset+t] = 0.
                rows = offset+self.decision_types
                self.ub_const[rows] = switch_bound[self.decision_types]
                self.ub_x[rows,np.arange(nx)] = -switch_bound[self.decision_types]
        self.linked = np.flatnonzero(np.any(self.ub_x!=0.,axis=1)|(self.ub_const!=self.y_ub_global))

    def _assemble_equations(self):
        """逐组构造物理方程，最后按等式、负等式、不等式、锥块拼接 A/B/C/b。"""
        net = self.network
        t, nx = len(self.type_corridor), len(self.variable_keys)
        ny, npower = len(self.y_ub_global), len(net.load_nodes)
        receiving = np.eye(net.n)[:,self.receivers]  # T(n,t)：型号功率汇入参考受端。
        sending = np.zeros((t,net.n))  # J(t,n)：读取参考送端电压；根节点对应零行。
        nonroot = np.flatnonzero(self.senders>=0)
        sending[nonroot,self.senders[nonroot]] = 1.
        root = (self.senders<0).astype(float)
        self.balance = balance = receiving-sending.T

        # 有功平衡：(T-J.T)@P - T@(r*ell) = (fixed_p + E@p)/base。
        B_active = np.zeros((net.n,ny))
        B_active[:,self.P_slice] = -balance
        B_active[:,self.ell_slice] = receiving*self.r
        C_active = net.E/net.base
        b_active = -net.fixed_p/net.base

        # 无功平衡：(T-J.T)@Q - T@(reactance*ell) = (fixed_q + E@(q_ratio*p))/base。
        B_reactive = np.zeros((net.n,ny))
        B_reactive[:,self.Q_slice] = -balance
        B_reactive[:,self.ell_slice] = receiving*self.reactance
        C_reactive = net.E*net.q_ratio/net.base
        b_reactive = -net.fixed_q/net.base

        # 压降：v_j-u_k + 2*(r_k*P_k+X_k*Q_k) - (r_k^2+X_k^2)*ell_k - s_plus+s_minus = 0。
        if self.reconfigurable:
            drop = np.eye(t)  # 可重构网架逐型号列压降；u=root+J@v。
            voltage_difference = receiving.T-sending
            b_drop = root
        else:
            drop = receiving  # 固定树按受端汇总同一走廊的型号项，每条走廊只列一个压降式。
            voltage_difference = np.eye(net.n)
            children = np.flatnonzero(net.parent>=0)
            voltage_difference[children,net.parent[children]] = -1.
            b_drop = (net.parent<0).astype(float)
        B_drop = np.zeros((len(drop),ny))
        B_drop[:,self.v_slice] = voltage_difference
        B_drop[:,self.P_slice] = 2*drop*self.r
        B_drop[:,self.Q_slice] = 2*drop*self.reactance
        B_drop[:,self.ell_slice] = -drop*(self.r**2+self.reactance**2)
        if self.reconfigurable:
            B_drop[:,self.switch_slice] = np.c_[-np.eye(t),np.eye(t)]
        B_eq = np.vstack([B_active,B_reactive,B_drop])
        C_eq = np.vstack([C_active,C_reactive,np.zeros((len(drop),npower))])
        b_eq = np.r_[b_active,b_reactive,b_drop]

        # 电压下限：-v<=-vmin；电压上限已在变量盒中。
        B_block = np.zeros((net.n,ny))
        B_block[:,self.v_slice] = -np.eye(net.n)
        linear_B = [B_block]
        linear_A = [np.zeros((net.n,nx))]
        linear_b = [-net.vmin]

        # 源端功率上限：root@P<=source_pmax，root@Q<=source_qmax，包含支路损耗。
        for block, limit in ((self.P_slice,net.source_pmax),(self.Q_slice,net.source_qmax)):
            if np.isfinite(limit):
                B_block = np.zeros((1,ny))
                B_block[:,block] = root
                linear_b.append(np.array([limit]))
                linear_B.append(B_block)
                linear_A.append(np.zeros((1,nx)))

        # 反向送端容量：-P+r*ell<=Pmax*z，-Q+X*ell<=Qmax*z；z=constant+selection@x。
        reversible = np.flatnonzero((self.senders>=0) & self.reconfigurable)
        pmax, qmax = self.y_ub_global[self.P_slice], self.y_ub_global[self.Q_slice]
        selection = self.type_selection
        constant = np.ones(t)-selection.sum(axis=1)
        for block, impedance, bound in ((self.P_slice,self.r,pmax),(self.Q_slice,self.reactance,qmax)):
            B_block = np.zeros((len(reversible),ny))
            B_block[np.arange(len(reversible)),block.start+reversible] = -1.
            B_block[np.arange(len(reversible)),self.ell_slice.start+reversible] = impedance[reversible]
            linear_b.append(bound[reversible]*constant[reversible])
            linear_B.append(B_block)
            linear_A.append(-bound[reversible,None]*selection[reversible])
        B_ineq, A_ineq, b_ineq = np.vstack(linear_B), np.vstack(linear_A), np.concatenate(linear_b)

        soc_B, soc_b = [], []
        if self.method == 'socp':
            # 支路电流锥：(u+ell,2P,2Q,u-ell)∈SOC，即 P^2+Q^2<=u*ell，u=root+J@v。
            for k in range(t):
                B_block = np.zeros((4,ny))
                B_block[0,self.v_slice] = B_block[3,self.v_slice] = -sending[k]
                B_block[0,self.ell_slice.start+k] = -1.
                B_block[1,self.P_slice.start+k] = B_block[2,self.Q_slice.start+k] = -2.
                B_block[3,self.ell_slice.start+k] = 1.
                soc_B.append(B_block)
                soc_b.append(np.array([root[k],0.,0.,root[k]]))
            # 源端视在功率锥：(Smax,root@P,root@Q)∈SOC，即 P_source^2+Q_source^2<=Smax^2。
            if np.isfinite(net.source_smax):
                B_block = np.zeros((3,ny))
                B_block[1,self.P_slice], B_block[2,self.Q_slice] = -root,-root
                soc_B.append(B_block)
                soc_b.append(np.array([net.source_smax,0.,0.]))

        # 统一拼接：等式 g=0 在 phase I 中表示为 g>=0、-g>=0，直接模型只读取 eq_rows。
        n_eq, n_linear = len(b_eq), 2*len(b_eq)+len(b_ineq)
        self.eq_rows = slice(0,n_eq)
        self.ineq_rows = slice(2*n_eq,n_linear)
        self.b = np.concatenate([b_eq,-b_eq,b_ineq,*soc_b])
        self.B = np.vstack([B_eq,-B_eq,B_ineq,*soc_B])
        self.A = np.zeros((len(self.b),nx))
        self.A[self.ineq_rows] = A_ineq
        self.C = np.zeros((len(self.b),npower))
        self.C[self.eq_rows], self.C[n_eq:2*n_eq] = C_eq,-C_eq
        self.soc_slices = []
        offset = n_linear
        for rhs in soc_b:
            self.soc_slices.append(slice(offset,offset+len(rhs)))
            offset += len(rhs)

        # phase I：b-Ax-By-Cp+eta*d∈K；线性行松弛 eta，SOC 只松弛首分量。
        self.relax_direction = np.zeros(len(self.b))
        self.relax_direction[:n_linear] = 1.
        for s in self.soc_slices:
            self.relax_direction[s.start] = 1.

    def activation(self, x):
        """完整型号向量：可选型号读取 x，固定单型号恒为 1。"""
        return np.ones(len(self.type_corridor))+self.type_selection@(x-1.)

    def selection(self, choice):
        """走廊方案编码为型号变量，与潮流方向无关。"""
        return np.array([choice[corridor_id] == type_id
                         for corridor_id, type_id in self.variable_keys], dtype=int)

    def choice(self, x):
        """唯一对外方案：走廊 ID -> 型号 ID；None 表示未投入。"""
        values = self.activation(np.asarray(x))
        plan = {}
        for corridor in self.corridors:
            plan[corridor.id] = next((type_id for type_id, index in self.type_indices[corridor.id].items()
                                     if values[index] > .5), None)
        return plan

    def investment(self, x):
        return self.constant_cost+self.cost@x

    def margin(self, x, power, state):
        slack = self.b-self.C@power-self.A@x-self.B@state
        values = [slack[:self.ineq_rows.stop].min(),np.min(state-self.lb_const-self.lb_x@x),
                  np.min(self.ub_const+self.ub_x@x-state)]
        values.extend(slack[s.start]-np.linalg.norm(slack[s.start+1:s.stop]) for s in self.soc_slices)
        return min(values)

    def restore(self, x, power, state):
        """用选中型号的电流重建等式精确的状态，消除未选型号退化锥中的浮点噪声。"""
        net = self.network.design(self.choice(x)) if self.reconfigurable else self.network
        active = self.activation(x)
        chosen = np.flatnonzero(active)
        forward = net.parent[self.receivers[chosen]] == self.senders[chosen]
        child = np.where(forward,self.receivers[chosen],self.senders[chosen])
        y = np.zeros_like(state)
        ell = np.clip(state[self.ell_slice][chosen],0.,self.y_ub_global[self.ell_slice][chosen])
        p,q = net.loads(power)
        ell[(net.D@(p[0]+q[0]))[child]==0.] = 0.  # 纯负荷空子树不保留无效损耗。
        r, reactance = self.r[chosen], self.reactance[chosen]
        receiving = np.eye(net.n)[:,child]
        P = net.D@(p[0]+receiving@(r*ell))
        Q = net.D@(q[0]+receiving@(reactance*ell))
        y[self.ell_slice][chosen] = ell
        y[self.P_slice][chosen] = np.where(forward,P[child],-P[child]+r*ell)
        y[self.Q_slice][chosen] = np.where(forward,Q[child],-Q[child]+reactance*ell)
        drop = 2*(r*P[child]+reactance*Q[child])-(r*r+reactance*reactance)*ell
        v = 1.-net.D.T@(receiving@drop)
        y[self.v_slice] = v
        if self.reconfigurable:
            delta = np.where(self.senders<0,1.,v[self.senders])-v[self.receivers]
            y[self.switch_slice] = np.r_[np.maximum(-delta,0.),np.maximum(delta,0.)]*np.tile(1.-active,2)
        return y


class PlanningModel:
    """逐线路 MILP/MISOCP；cuts_only=True 时为包含 x、p 的联合割主问题。"""

    def __init__(self, equations, *, power=None, budget=np.inf, cuts_only=False,
                  min_total=None, fixed_plan=None, threads=DEFAULT_SOLVER_THREADS):
        self.equations = e = equations
        network, method = e.network, e.method
        self.model = m = gp.Model('compact_'+method)
        m.Params.OutputFlag = 0
        m.Params.Threads = threads
        m.Params.MIPGap = m.Params.MIPGapAbs = 0.  # 不以相对或绝对整数间隙提前接受次优投资。
        m.Params.FeasibilityTol = m.Params.IntFeasTol = m.Params.OptimalityTol = 1e-9  # 分别设置原始约束、整数性和最优性容差。
        m.Params.BarQCPConvTol = 1e-9  # 设置二阶锥连续子问题的内点法收敛精度。
        m.Params.DualReductions = 0  # 区分不可行与无界，避免模糊状态被当作不可行证明。
        m.Params.NonConvex = 0  # 若模型误写成非凸二次式，立即暴露，不能悄悄改用非凸求解。

        # 决策变量：先完整声明，再统一建立表达式和约束。
        self.x = m.addVars(e.variable_keys, vtype=GRB.BINARY, name='line_type')
        self.x_vector = gp.MVar.fromlist([self.x[key] for key in e.variable_keys])  # 与方程列顺序一致的矩阵视图。
        self.power = p = m.addMVar(len(network.load_nodes),lb=0.,name='p_kw')
        connectivity = (m.addMVar(len(e.corridors),lb=-network.n,ub=network.n,name='connectivity')
                        if e.reconfigurable else None)
        self.state = (None if cuts_only else
                      m.addMVar(len(e.y_ub_global),lb=e.y_lb_global,ub=e.y_ub_global,name='flow'))

        # 公共表达式：这里只组合变量，不向模型添加约束。
        investment = e.investment(self.x_vector)
        total_power = p.sum()/network.base
        minimize_investment = power is not None or min_total is not None
        if e.reconfigurable:
            type_active = e.activation(self.x_vector)
            corridor_active = np.eye(len(e.corridors))[:,e.type_corridor]@type_active
            # 同一走廊各型号端点相同，连通流只需取每条走廊的一列关联矩阵。
            first = [next(iter(indices.values())) for indices in e.type_indices.values()]
            corridor_balance = e.balance[:,first]
        if self.state is not None:
            y = self.state
            slack = e.b-e.C@p-e.A@self.x_vector-e.B@y  # b-Ax-By-Cp∈K。

        # 线路选型：必用走廊恰选一个型号，可选走廊至多选择一个型号。
        for corridor in e.corridors:
            if corridor.must_use and len(corridor.types) == 1:
                continue  # 该型号恒为 1，已经作为常量写入方程，不占用 x。
            total = self.x.sum(corridor.id, '*')
            if corridor.must_use:
                m.addConstr(total == 1, name='required_' + corridor.id)
            else:
                m.addConstr(total <= 1, name='selectable_' + corridor.id)

        # 固定方案直接收紧二进制变量上下界，不再添加一组重复等式。
        if fixed_plan is not None:
            self.x_vector.LB = self.x_vector.UB = e.selection(fixed_plan)

        # 可重构网络：n 条投入走廊加单商品连通流，刻画覆盖 n+1 个节点的树。
        if e.reconfigurable:
            m.addConstr(corridor_active.sum()==network.n, name='tree_edges')
            m.addConstr(connectivity<=network.n*corridor_active, name='connectivity_ub')
            m.addConstr(connectivity>=-network.n*corridor_active, name='connectivity_lb')
            m.addConstr(corridor_balance@connectivity==np.ones(network.n), name='connectivity_balance')

        # 完整运行模型；cuts_only 主问题只保留 x、p 和后续加入的联合割。
        if self.state is not None:
            m.addConstr(slack[e.eq_rows]==0., name='flow_equalities')
            m.addConstr(slack[e.ineq_rows]>=0., name='operating_limits')
            # U*x 关闭未选型号潮流，M*(1-x) 仅在走廊断开时允许压降余量。
            m.addConstr(y[e.linked]<=e.ub_const[e.linked]+e.ub_x[e.linked]@self.x_vector,
                        name='state_activation')
            if method == 'socp':
                for k in range(len(e.type_corridor)):
                    P,Q,l = (y[s.start+k].item() for s in (e.P_slice,e.Q_slice,e.ell_slice))
                    parent = e.senders[k]
                    u = 1. if parent < 0 else y[e.v_slice.start+parent].item()
                    m.addQConstr(P*P+Q*Q<=u*l, name=f'branch_cone[{k}]')  # u、ell 非负，因此是凸旋转锥。
                if np.isfinite(network.source_smax):
                    roots = e.senders<0
                    ps,qs = y[e.P_slice][roots].sum().item(),y[e.Q_slice][roots].sum().item()
                    m.addQConstr(ps*ps+qs*qs<=network.source_smax**2, name='source_cone')

        # 负荷、预算和查询条件。
        m.addConstr(total_power<=network.power_limit/network.base, name='power_limit')
        if np.isfinite(budget):
            m.addConstr(investment<=budget, name='budget')
        if power is not None:
            m.addConstr(p==np.asarray(power), name='fixed_power')
        elif min_total is not None:
            m.addConstr(total_power>=float(min_total)/network.base, name='minimum_total')

        # 查询目标：自由负荷查询最大化总量，其余查询在给定负荷要求下最小化投资。
        if minimize_investment:
            m.setObjective(investment,GRB.MINIMIZE)
            self.objective_scale = 1.
        else:
            m.setObjective(total_power,GRB.MAXIMIZE)
            self.objective_scale = network.base

    def add_cut(self, cut):
        n = len(self.equations.network.load_nodes)
        self.model.addConstr(cut[0]+cut[1:1+n]@self.power+cut[1+n:]@self.x_vector>=0.)  # 联合割对全部合法选型有效。


    def use_incumbent(self, incumbent):
        """复用同模型、同预算且满足本次负荷要求的已认证解。"""
        cost = float(self.equations.investment(incumbent['x']))
        self.x_vector.Start = incumbent['x']
        self.power.Start = incumbent['p']
        if self.state is not None:
            self.state.Start = incumbent['state']
        self.model.addConstr(self.equations.investment(self.x_vector) <= cost+1e-9)
        return cost

    def exclude(self, x):
        self.model.addConstr((1-2*x)@self.x_vector>=1-x.sum())  # Hamming 距离至少为一，排除且只排除这个离散建设向量。

    def solve(self, time_limit=MP_TIME_LIMIT):
        self.model.Params.TimeLimit = max(0., time_limit)
        self.model.optimize()
        if self.model.Status == GRB.INFEASIBLE:
            return None
        bound = self.model.ObjBound*self.objective_scale
        if not self.model.SolCount:
            return dict(status='unknown', x=None, p=None, objective=None, bound=bound, state=None, feasible=False)
        x = np.rint(self.x_vector.X).astype(int)
        objective = self.equations.investment(x) if self.model.ModelSense==GRB.MINIMIZE else self.power.X.sum()
        # 投资按实际整数选型计价，避免 1+1e-10 被误判为超出预算 1。
        state = None if self.state is None else self.equations.restore(x,self.power.X,self.state.X)
        feasible = state is not None and self.equations.margin(x,self.power.X,state)>=-PLANNING_TOL
        return dict(x=x,p=self.power.X,objective=float(objective),bound=bound,state=state,feasible=feasible,
                    status='optimal' if self.model.Status==GRB.OPTIMAL and (feasible or self.state is None) else ('feasible' if feasible else 'unknown'))


class PlanningSP:
    """固定 x,p 的 phase I；带符号变量盒产生全局有效的联合割。"""

    def __init__(self, equations, *, threads=DEFAULT_SOLVER_THREADS):
        self.equations = equations
        self.settings = settings = clarabel.DefaultSettings()
        self.threads, self.calls = threads, 0
        settings.verbose = False
        settings.max_threads = threads
        settings.tol_gap_abs = settings.tol_gap_rel = settings.tol_feas = 1e-12  # 控制辅助问题的间隙与原始可行性误差。
        settings.reduced_tol_feas = 1e-11  # 数值停滞时的较宽松终止精度不替代原约束检验。

    def solve(self, x, power, time_limit=None):
        self.calls += 1
        self.settings.time_limit = SP_TIME_LIMIT[self.equations.method] if time_limit is None else max(0.,time_limit)
        deadline = perf_counter()+self.settings.time_limit
        e = self.equations
        upper = e.ub_const+e.ub_x@x
        lower = e.lb_const+e.lb_x@x
        columns = np.flatnonzero(upper>lower)
        active = e.activation(x).astype(bool)
        soc_slices = [s for k,s in enumerate(e.soc_slices) if k>=len(active) or active[k]]
        n_linear = e.ineq_rows.stop  # phase I 的非负锥包含正负等式及运行不等式。
        rows = np.r_[np.arange(n_linear),*[np.arange(s.start,s.stop) for s in soc_slices]].astype(int)
        n = len(columns)
        B = e.B[np.ix_(rows,columns)]
        matrix = np.vstack([np.c_[B,-e.relax_direction[rows]],
                            np.c_[-np.eye(n),np.zeros(n)],np.c_[np.eye(n),np.zeros(n)],
                            np.r_[np.zeros(n),-1.]])  # By-eta*r≼_K b-Ax-Cp，以及严格变量盒和 eta≥0。
        rhs = np.r_[(e.b-e.C@power-e.A@x)[rows],-lower[columns],upper[columns],0.]
        kinds = [clarabel.NonnegativeConeT(n_linear)]
        kinds += [clarabel.SecondOrderConeT(s.stop-s.start) for s in soc_slices]
        kinds += [clarabel.NonnegativeConeT(2*n+1)]
        solver = clarabel.DefaultSolver(sparse.csc_matrix((n+1,n+1)),np.r_[np.zeros(n),1.],
            sparse.csc_matrix(matrix),rhs,kinds,self.settings)
        result = solver.solve()
        if not np.isfinite(np.r_[result.x,result.z]).all():
            return dict(cut=None, state=None, feasible=False)
        state = np.zeros(len(e.y_ub_global))
        state[columns] = result.x[:-1]
        state = e.restore(x,power,state)
        if e.margin(x,power,state)>=-PLANNING_TOL:
            return dict(cut=None, state=state, feasible=True)
        # 未收敛时仍可检验候选证书；有效性由原约束或下方的对偶锥与盒补偿保证。
        dual = np.zeros(len(e.b))
        dual[rows] = result.z[:len(rows)]
        dual[:n_linear] = np.maximum(dual[:n_linear],0.)
        for s in e.soc_slices:
            dual[s.start] = max(dual[s.start],np.linalg.norm(dual[s.start+1:s.stop]))  # SOC 自对偶；首分量至少为尾部范数才有全域有效性。
        h = -e.B.T@dual  # λᵀ(b-Ax-By-Cp) 中 y 的系数。
        positive, negative = np.maximum(h,0.),np.minimum(h,0.)
        # max_{lb(x)≤y≤ub(x)} hᵀy = h_+ᵀub(x)+h_-ᵀlb(x)，包含未选型号的全部列。
        cut = np.r_[dual@e.b+positive@e.ub_const+negative@e.lb_const+1e-12,
                    -e.C.T@dual,-e.A.T@dual+e.ub_x.T@positive+e.lb_x.T@negative]
        scale = max(np.max(np.abs(cut)),np.max(np.abs(cut[1:1+len(power)]))*e.network.base)
        cut /= max(scale,1e-30)
        if cut[0]+cut[1:1+len(power)]@power+cut[1+len(power):]@x < -1e-9:
            return dict(cut=cut, state=None, feasible=False)
        state = self._polish(x, power, columns, rows, soc_slices, lower, upper, deadline)
        return dict(cut=None, state=state, feasible=state is not None)

    def _polish(self, x, power, columns, rows, soc_slices, lower, upper, deadline):
        """边界无有效割时精修原等式；两条路径共用 SP 剩余时限。"""
        if perf_counter() >= deadline:
            return None
        e, n = self.equations, len(columns)
        # 零锥表示严格等式，避免 phase I 成对松弛行的边界退化。
        exact_rows = np.r_[np.arange(e.eq_rows.stop), rows[rows >= e.ineq_rows.start]].astype(int)
        identity = np.eye(n)
        matrix = np.vstack([e.B[np.ix_(exact_rows, columns)], -identity, identity])
        rhs = np.r_[(e.b-e.C@power-e.A@x)[exact_rows], -lower[columns], upper[columns]]
        kinds = [clarabel.ZeroConeT(e.eq_rows.stop),
                 clarabel.NonnegativeConeT(e.ineq_rows.stop-e.ineq_rows.start)]
        kinds += [clarabel.SecondOrderConeT(s.stop-s.start) for s in soc_slices]
        kinds += [clarabel.NonnegativeConeT(2*n)]
        self.settings.time_limit = max(0., deadline-perf_counter())
        exact = clarabel.DefaultSolver(sparse.csc_matrix((n, n)), np.zeros(n),
                                      sparse.csc_matrix(matrix), rhs, kinds, self.settings).solve()
        if np.isfinite(exact.x).all():
            candidate = np.zeros(len(e.y_ub_global))
            candidate[columns] = exact.x
            candidate = e.restore(x, power, candidate)
            if e.margin(x, power, candidate) >= -PLANNING_TOL:
                return candidate
        if perf_counter() >= deadline:
            return None
        problem = PlanningModel(e, power=power, threads=self.threads)
        with problem.model:
            problem.x_vector.LB = problem.x_vector.UB = x
            polished = problem.solve(time_limit=max(0., deadline-perf_counter()))
        return polished['state'] if polished is not None and polished['feasible'] else None


class RemainingRegionModel:
    """共享并集排除约束；light 搜索割外域，physical 再加入完整运行约束。"""

    def __init__(self, equations, budget, bounds, total_bound, cuts, inner_halfspaces, tau,
                  *, mode='light', threads=DEFAULT_SOLVER_THREADS):
        self.mode = mode
        self.problem = problem = PlanningModel(equations, budget=budget, cuts_only=mode == 'light', threads=threads)
        m = self.model = problem.model
        problem.power.UB = bounds
        m.addConstr(problem.power.sum() <= total_bound)
        self.distance_scale = 1000.
        for cut in cuts:
            cut = np.asarray(cut)
            scale = max(np.max(np.abs(np.r_[cut[0], cut[1:4]*bounds, cut[4:]])), 1e-20)
            problem.add_cut(cut*self.distance_scale/scale)
        # 等价缩放，避免 1e-8 的覆盖阈值与求解器绝对容差处于相近量级。
        delta = m.addVar(lb=-4.*self.distance_scale, ub=4.*self.distance_scale, name='uncovered_distance')
        for k, eq in enumerate(inner_halfspaces):
            select = m.addVars(len(eq), vtype=GRB.BINARY, name=f'outside_{k}')
            m.addConstr(select.sum() == 1)
            for f, face in enumerate(eq):
                a = face[:3]*(1-tau)
                # z=p/bounds in [0,1]^3；M 由该盒的精确下界推出。
                big_m = 4.-face[3]-np.minimum(a, 0.).sum()
                expression = gp.quicksum(float(a[j]/bounds[j])*problem.power[j].item() for j in range(3))
                m.addConstr(delta <= self.distance_scale*(expression+float(face[3])+float(big_m)*(1-select[f])))
        m.setObjective(delta, GRB.MAXIMIZE)

    def solve(self, tolerance, time_limit=RESIDUAL_TIME_LIMIT):
        m, problem = self.model, self.problem
        m.Params.TimeLimit = time_limit
        # 探索阶段只需一个可靠未覆盖点；停止证明仍必须检查全局上界。
        m.Params.BestObjStop = max(10*tolerance, 1e-6)*self.distance_scale
        m.optimize()
        details = dict(feasible=False)
        if m.Status == GRB.INFEASIBLE:
            return dict(complete=True, bound=None, x=None, p=None, **details)
        bound = float(m.ObjBound)/self.distance_scale
        if bound <= tolerance:
            return dict(complete=True, bound=bound, x=None, p=None, **details)
        if not m.SolCount or m.ObjVal/self.distance_scale <= tolerance:
            return dict(complete=False, bound=bound, x=None, p=None, **details)
        x, point = np.rint(problem.x_vector.X).astype(int), problem.power.X
        if self.mode == 'physical':
            e = problem.equations
            state = e.restore(x, point, problem.state.X)
            margin = e.margin(x, point, state)
            details.update(feasible=margin >= -PLANNING_TOL)
        return dict(complete=False, bound=bound, x=x, p=point, **details)


def planning_query(equations, *, power=None, budget=np.inf,
                   cuts=(), start=None, radial_gap_kw=1e-3, min_total=None,
                   fixed_plan=None, incumbent=None, deadline=np.inf,
                   threads=DEFAULT_SOLVER_THREADS, oracle=None):
    """完整 MP1/MP2；原约束认证通过后直接复用，必要时才调用 SP 修复。"""
    minimizing = power is not None or min_total is not None
    problem = PlanningModel(equations, power=power, budget=budget,
                            min_total=min_total, fixed_plan=fixed_plan, threads=threads)
    with problem.model:
        for cut in cuts:
            problem.add_cut(cut)
        if incumbent is not None:
            problem.use_incumbent(incumbent)
        elif start is not None:
            problem.x_vector.Start = start
        answer = problem.solve(time_limit=min(MP_TIME_LIMIT, max(0., deadline-perf_counter())))
    if answer is None:
        return None
    if answer['x'] is not None and not minimizing and answer['p'].sum() > 0.:
        point = answer['p']*max(0., 1.-radial_gap_kw/answer['p'].sum())
        state = equations.restore(answer['x'], point, answer['state'])
        if equations.margin(answer['x'], point, state) >= -PLANNING_TOL:
            answer.update(p=point, state=state, feasible=True, objective=float(point.sum()))
    if answer['x'] is not None and not answer['feasible'] and perf_counter() < deadline:
        oracle = oracle or PlanningSP(equations, threads=threads)
        checked = oracle.solve(answer['x'], answer['p'], time_limit=min(SP_TIME_LIMIT[equations.method], deadline-perf_counter()))
        if checked['feasible']:
            answer.update(state=checked['state'], feasible=True)
    if answer['feasible']:
        gap = answer['objective']-answer['bound'] if minimizing else answer['bound']-answer['objective']
        answer['status'] = 'optimal' if gap <= (1e-7 if minimizing else radial_gap_kw+1e-5) else 'feasible'
    else:
        answer['status'] = 'unknown'
    return answer


def evaluation_bounds(network, *, threads=DEFAULT_SOLVER_THREADS):
    """三个无预算 LP 轴向全局上界确定公共评价箱。"""
    equations, bounds = PlanningEquations(network, 'linear'), []
    for axis in np.eye(len(network.load_nodes)):
        problem = PlanningModel(equations, threads=threads)
        problem.power.UB = axis*network.power_limit
        with problem.model:
            answer = problem.solve()
        if answer is None or answer['bound'] is None or not np.isfinite(answer['bound']):
            raise RuntimeError('No finite planning bound for the common evaluation box')
        bounds.append(min(network.power_limit, answer['bound']))
    return np.ceil(np.asarray(bounds)/10.)*10.
