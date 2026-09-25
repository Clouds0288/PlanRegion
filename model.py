"""Gurobi 配电网规划：e 为走廊 ID，k 为该走廊的型号 ID，i 为节点 ID。

MP 与 SP 共用逐节点/逐走廊约束；x[e,k] 选型，P/Q/ell[e,k] 为标幺运行量，
p[i] 为 kW 负荷，v[i] 为电压平方。数组仅用于结果、几何算法和割的传递。
固定符号、变量名、单位和索引约定见 docs/notation.md；同义量不随重构改名。
"""
from time import perf_counter
from types import SimpleNamespace

import gurobipy as gp
import numpy as np
from gurobipy import GRB

DEFAULT_SOLVER_THREADS = 20
MP_TIME_LIMIT = 20.
SP_TIME_LIMIT = {'linear': 5., 'socp': 10.}
RESIDUAL_TIME_LIMIT = 120.
PLANNING_TOL = 1e-8


def new_model(name, threads):
    model = gp.Model(name)
    model.Params.OutputFlag, model.Params.Threads = 0, threads
    model.Params.MIPGap = model.Params.MIPGapAbs = 0.
    model.Params.FeasibilityTol = model.Params.IntFeasTol = model.Params.OptimalityTol = 1e-9
    model.Params.BarQCPConvTol = 1e-10
    model.Params.DualReductions, model.Params.NonConvex = 0, 0
    return model


class PlanningEquations:
    """共用物理参数、索引和约束；不再手工装配标准锥矩阵。"""

    def __init__(self, network, method):
        self.network, self.method = network, method
        self.keys = network.type_keys
        self.types = {c.id: tuple(k.id for k in c.types) for c in network.corridors}
        self.ends = {c.id: (network.root if a < 0 else network.nodes[a], network.nodes[b])
                     for c, a, b in zip(network.corridors, network.senders, network.receivers)}
        self.incoming = {i: tuple(e for e, (_, j) in self.ends.items() if j == i)
                         for i in network.nodes}
        self.outgoing = {i: tuple(e for e, (j, _) in self.ends.items() if j == i)
                         for i in (network.root, *network.nodes)}
        self.r = dict(zip(self.keys, network.r))
        self.reactance = dict(zip(self.keys, network.reactance))
        self.cost = dict(zip(self.keys, network.cost))
        t, n, m = network.n_types, network.n, network.n_corridors
        self.P_slice, self.Q_slice, self.ell_slice = slice(0, t), slice(t, 2*t), slice(2*t, 3*t)
        self.v_slice, self.slack_slice = slice(3*t, 3*t+n), slice(3*t+n, 3*t+n+2*m)
        self._build_variable_bounds()

    def _build_variable_bounds(self):
        net = self.network
        pmax = np.minimum(net.capacity, min(net.source_pmax, net.source_smax))
        qmax = np.full(net.n_types, min(net.source_qmax, net.source_smax))
        ellmax = pmax/net.r if self.method == 'socp' else np.zeros(net.n_types)
        reverse = net.senders[net.type_corridor] >= 0
        self.pmax, self.qmax, self.ellmax = (dict(zip(self.keys, a)) for a in (pmax, qmax, ellmax))
        self.pmin, self.qmin = (dict(zip(self.keys, -a*reverse)) for a in (pmax, qmax))
        self.vmin, self.vmax = (dict(zip(net.nodes, a)) for a in (net.vmin, net.vmax))
        self.vmin[net.root] = self.vmax[net.root] = 1.
        # 开断走廊两端电压可独立变化：|v_j-v_i| <= M_e。
        self.drop_max = {e: max(self.vmax[i]-self.vmin[j], self.vmax[j]-self.vmin[i]) for e, (i, j) in self.ends.items()}
        drop = list(self.drop_max.values())
        self.y_lb_global = np.r_[-pmax*reverse, -qmax*reverse, np.zeros(net.n_types+net.n+2*net.n_corridors)]
        self.y_ub_global = np.r_[pmax, qmax, ellmax, net.vmax, drop, drop]

    def add_operation(self, model, x, p, eta=0.):
        """MP 用 eta=0；SP 只放松功率平衡与压降等式，容量、电压界和锥保持原式。"""
        net = self.network
        P = model.addVars(self.keys, lb=self.pmin, ub=self.pmax, name='P')
        Q = model.addVars(self.keys, lb=self.qmin, ub=self.qmax, name='Q')
        ell = model.addVars(self.keys, ub=self.ellmax, name='ell')
        v = model.addVars(net.nodes, ub={i: self.vmax[i] for i in net.nodes}, name='v')
        v[net.root] = 1.
        plus = model.addVars(self.types, ub=self.drop_max, name='drop_plus')
        minus = model.addVars(self.types, ub=self.drop_max, name='drop_minus')
        z = {e: gp.quicksum(x[e, k] for k in self.types[e]) for e in self.types}

        def equality(expression, name):
            if isinstance(eta, gp.Var):
                model.addConstr(expression <= eta, name=name+'_upper')
                model.addConstr(expression >= -eta, name=name+'_lower')
            else:
                model.addConstr(expression == 0., name=name)

        fixed_p, fixed_q = dict(zip(net.nodes, net.fixed_p)), dict(zip(net.nodes, net.fixed_q))
        ratio = dict(zip(net.load_nodes, net.q_ratio))
        for i in net.nodes:
            # 节点有功/无功守恒：流入送端功率扣除线损，再减流出功率，等于负荷。
            incoming_p = gp.quicksum(P[e, k]-self.r[e, k]*ell[e, k]
                                    for e in self.incoming[i] for k in self.types[e])
            incoming_q = gp.quicksum(Q[e, k]-self.reactance[e, k]*ell[e, k]
                                    for e in self.incoming[i] for k in self.types[e])
            outgoing_p = gp.quicksum(P[e, k] for e in self.outgoing[i] for k in self.types[e])
            outgoing_q = gp.quicksum(Q[e, k] for e in self.outgoing[i] for k in self.types[e])
            equality(incoming_p-outgoing_p-(fixed_p[i]+p.get(i, 0.))/net.base, f'active_balance[{i}]')
            equality(incoming_q-outgoing_q-(fixed_q[i]+ratio.get(i, 0.)*p.get(i, 0.))/net.base,
                     f'reactive_balance[{i}]')
            # 节点电压平方下限；上限已在变量声明中给出。
            model.addConstr(v[i] >= self.vmin[i], name=f'voltage_min[{i}]')

        cones = []
        for e, (i, j) in self.ends.items():
            # 支路压降：v_j-v_i + 2(rP+χQ) - (r²+χ²)ell = s⁺-s⁻。
            drop = v[j]-v[i]+gp.quicksum(
                2*(self.r[e, k]*P[e, k]+self.reactance[e, k]*Q[e, k])
                -(self.r[e, k]**2+self.reactance[e, k]**2)*ell[e, k] for k in self.types[e])
            equality(drop-plus[e]+minus[e], f'voltage_drop[{e}]')
            # 接通时压降松弛为零；开断时允许两端电压不同。
            model.addConstr(plus[e] <= self.drop_max[e]*(1-z[e]), name=f'drop_plus_bound[{e}]')
            model.addConstr(minus[e] <= self.drop_max[e]*(1-z[e]), name=f'drop_minus_bound[{e}]')
            for k in self.types[e]:
                # 未选型号 P=Q=ell=0；非根走廊允许参考方向的反向潮流。
                for flow, lower, upper, name in ((P, self.pmin, self.pmax, 'P'),
                                                  (Q, self.qmin, self.qmax, 'Q')):
                    model.addConstr(flow[e, k] >= lower[e, k]*x[e, k], name=f'{name}_min[{e},{k}]')
                    model.addConstr(flow[e, k] <= upper[e, k]*x[e, k], name=f'{name}_max[{e},{k}]')
                model.addConstr(ell[e, k] <= self.ellmax[e, k]*x[e, k], name=f'current_max[{e},{k}]')
                if i != net.root:
                    # 反向送端功率含线路损耗；两端容量均须满足。
                    model.addConstr(-P[e, k]+self.r[e, k]*ell[e, k] <= self.pmax[e, k]*x[e, k],name=f'reverse_P[{e},{k}]')
                    model.addConstr(-Q[e, k]+self.reactance[e, k]*ell[e, k] <= self.qmax[e, k]*x[e, k],name=f'reverse_Q[{e},{k}]')
                if self.method == 'socp':
                    # 支路电流锥：P²+Q² <= v_i*ell，MP/SP 共用同一凸约束。
                    model.addQConstr(P[e, k]**2+Q[e, k]**2 <= v[i]*ell[e, k], name=f'current_cone[{e},{k}]')
                    cones.append((gp.LinExpr(v[i]+ell[e, k]),[gp.LinExpr(2*P[e, k]), gp.LinExpr(2*Q[e, k]), gp.LinExpr(v[i]-ell[e, k])]))

        source_p = gp.quicksum(P[e, k] for e in self.outgoing[net.root] for k in self.types[e])
        source_q = gp.quicksum(Q[e, k] for e in self.outgoing[net.root] for k in self.types[e])
        # 电源注入含全网损耗：有功、无功以及 SOCP 的视在功率限额。
        for flow, limit, name in ((source_p, net.source_pmax, 'source_P'), (source_q, net.source_qmax, 'source_Q')):
            if np.isfinite(limit):
                model.addConstr(flow <= limit, name=name)
        if self.method == 'socp' and np.isfinite(net.source_smax):
            model.addQConstr(source_p**2+source_q**2 <= net.source_smax**2, name='source_capacity')
            cones.append((gp.LinExpr(net.source_smax), [source_p, source_q]))
        state = gp.MVar.fromlist([*P.values(), *Q.values(), *ell.values(),
                                 *(v[i] for i in net.nodes), *plus.values(), *minus.values()])
        return SimpleNamespace(P=P, Q=Q, ell=ell, v=v, plus=plus, minus=minus, state=state, cones=cones)


class PlanningModel:
    """给定 power 或 min_total 时求 MP1 最小投资，否则求 MP2 最大总负荷。"""

    def __init__(self, equations, *, power=None, budget=np.inf, cuts_only=False,
                 min_total=None, fixed_plan=None, cuts=(), threads=DEFAULT_SOLVER_THREADS):
        # 1. 读取网架参数，创建求解模型
        net = equations.network
        model = new_model('planning_'+equations.method, threads)

        # 2. 声明变量：选型 x、负荷 p、节点接入 a、虚拟连通流 f；z 为走廊接通表达式
        x = model.addVars(equations.keys, vtype=GRB.BINARY, name='x')
        p = model.addVars(net.load_nodes, name='p_kw')
        a = model.addVars(net.nodes, lb=dict(zip(net.nodes, net.required.astype(float))), ub=1., vtype=GRB.BINARY, name='active')
        f = model.addVars(equations.types, lb=-net.n, ub=net.n, name='connectivity')
        z = {e: gp.quicksum(x[e, k] for k in equations.types[e]) for e in equations.types}
        for corridor, allowed in zip(net.corridors, net.road_allowed):
            if not allowed:
                model.addConstr(z[corridor.id] == 0., name=f'road_blocked[{corridor.id}]')

        # 3. 每条走廊最多选一种型号，接通走廊的两端必须接入
        for e, (i, j) in equations.ends.items():
            model.addConstr(z[e] <= 1., name=f'one_type[{e}]')
            model.addConstr(z[e] <= a[j], name=f'receiver_active[{e}]')
            if i != net.root:
                model.addConstr(z[e] <= a[i], name=f'sender_active[{e}]')

        # 4. 虚拟流保证接入节点与根连通，边数约束保证径向结构
        for e in equations.types:
            model.addConstr(f[e] <= net.n*z[e], name=f'connectivity_max[{e}]')
            model.addConstr(f[e] >= -net.n*z[e], name=f'connectivity_min[{e}]')  # -n*z_e <= f_e <= n*z_e：开断走廊不通流。
        for i in net.nodes:
            model.addConstr(gp.quicksum(f[e] for e in equations.incoming[i])-gp.quicksum(f[e] for e in equations.outgoing[i]) == a[i], name=f'connectivity_balance[{i}]')  # 每个接入节点消耗一单位虚拟流。
        model.addConstr(gp.quicksum(z.values()) == gp.quicksum(a.values()), name='tree_edges')  # 边数 = 非根接入节点数。

        # 5. 限制预算和总负荷，并设置给定负荷或最低总负荷
        investment = gp.quicksum(equations.cost[e, k]*x[e, k] for e, k in equations.keys)
        total_power = gp.quicksum(p.values())/net.base
        model.addConstr(total_power <= net.power_limit/net.base, name='power_limit')  # 总负荷采用标幺缩放。
        if np.isfinite(budget):
            model.addConstr(investment <= budget, name='budget')  # 投资不超过预算，单位为算例费用单位。
        if power is not None:
            model.addConstrs((p[i] == value for i, value in zip(net.load_nodes, power)), name='fixed_power')  # 固定各节点负荷，单位 kW。
        elif min_total is not None:
            model.addConstr(total_power >= min_total/net.base, name='minimum_total')  # 仅限制总量，允许节点间重新分配。

        # 6. 如给定建设方案，则固定每个走廊的型号选择
        if fixed_plan is not None:
            for (e, k), value in zip(equations.keys, net.encode_plan(fixed_plan)):
                x[e, k].LB = x[e, k].UB = value

        # 7. 添加运行变量及物理约束；cuts_only 模式省略此部分
        operation = None if cuts_only else equations.add_operation(model, x, p)

        # 8. MP1 最小投资，MP2 最大总负荷；返回时将 MP2 目标换回 kW
        minimizing = power is not None or min_total is not None
        model.setObjective(investment if minimizing else total_power, GRB.MINIMIZE if minimizing else GRB.MAXIMIZE)
        self.objective_scale = 1. if minimizing else net.base

        # 9. 保存变量引用；具名变量与扁平视图共享同一批 Gurobi 变量
        self.equations, self.model = equations, model
        self.choices, self.loads = x, p
        self.x = gp.MVar.fromlist(list(x.values()))
        self.power = gp.MVar.fromlist(list(p.values()))
        self.active_nodes = gp.MVar.fromlist(list(a.values()))
        self.operation = operation
        self.state = None if cuts_only else operation.state
        for cut in cuts:
            self.add_cut(cut)

    def add_cut(self, cut):
        """添加联合可行性割：α + Σβ_i p_i + Σδ_ek x_ek >= 0。"""
        d = len(self.equations.network.load_nodes)
        self.model.addConstr(cut[0]+gp.quicksum(c*p for c, p in zip(cut[1:1+d], self.loads.values()))+gp.quicksum(c*x for c, x in zip(cut[1+d:], self.choices.values())) >= 0.)

    def use_incumbent(self, incumbent):
        """提供已有可行解作为起点，并将其投资作为成本上界。"""
        cost = float(self.equations.network.cost@incumbent['x'])
        self.x.Start, self.power.Start = incumbent['x'], incumbent['p']
        if self.state is not None:
            self.state.Start = incumbent['state']
        self.model.addConstr(self.equations.network.cost@self.x <= cost+1e-9)
        return cost

    def exclude(self, x):
        """要求至少一个型号选择不同，排除指定建设方案。"""
        self.model.addConstr((1-2*x)@self.x >= 1-x.sum())

    def solve(self, time_limit=MP_TIME_LIMIT, *, incumbent=None, start=None, radial_gap_kw=1e-3):
        """直接求 MP1/MP2；初始解、目标间隙及运行证书在此处理，不调用 SP。"""
        # 1. 设置初始解、时限并求解
        model, equations = self.model, self.equations
        if incumbent is not None:
            self.use_incumbent(incumbent)
        elif start is not None:
            self.x.Start = start
        model.Params.TimeLimit = max(0., time_limit)
        model.optimize()

        # 2. 判断是否有候选解，提取目标的全局界
        if model.Status == GRB.INFEASIBLE:
            return None
        bound = model.ObjBound*self.objective_scale
        if not model.SolCount:
            return dict(status='unknown', x=None, p=None, objective=None, bound=bound, state=None, feasible=False)

        # 3. 提取选型、负荷和目标值：投资用原费用单位，负荷用 kW
        x = np.rint(self.x.X).astype(int)
        p = self.power.X
        objective = equations.network.cost@x if model.ModelSense == GRB.MINIMIZE else p.sum()

        # 4. 原样读取运行状态；质量指标覆盖所建模型，cuts_only 没有运行证书。
        state = None if self.state is None else self.state.X
        feasible = state is not None and model.MaxVio <= PLANNING_TOL

        # 5. 返回候选解、目标界和认证状态
        minimizing = model.ModelSense == GRB.MINIMIZE
        gap = objective-bound if minimizing else bound-objective
        tolerance = 1e-7 if minimizing else radial_gap_kw+1e-5
        status = ('optimal' if model.MaxVio <= PLANNING_TOL and (
                      model.Status == GRB.OPTIMAL or feasible and gap <= tolerance)
                  else 'feasible' if feasible else 'unknown')
        return dict(x=x, p=p, objective=float(objective), bound=bound, state=state, feasible=feasible, status=status)


class PlanningSP:
    """固定 x,p 的 Gurobi 可行性问题；锥支撑平面的 LP 对偶生成全局联合割。"""

    def __init__(self, equations, *, threads=DEFAULT_SOLVER_THREADS):
        self.equations, self.threads, self.calls = equations, threads, 0

    def solve(self, x, power, time_limit=None):
        self.calls += 1
        equations, net = self.equations, self.equations.network
        limit = SP_TIME_LIMIT[equations.method] if time_limit is None else time_limit
        deadline = perf_counter()+limit
        with new_model('planning_SP', self.threads) as model:
            # 边界及零潮流锥易病态；提高原始解精度，避免依赖事后状态重建。
            model.Params.NumericFocus = 2
            choice = model.addVars(equations.keys, ub=1., name='x')
            p = model.addVars(net.load_nodes, lb=-GRB.INFINITY, name='p_kw')
            # 固定参数用等式表示；求割时去掉这些等式的乘子，保留 x,p 系数。
            fixed_x = model.addConstrs((choice[e, k] == value for (e, k), value in zip(equations.keys, x)), name='fixed_x')
            fixed_p = model.addConstrs((p[i] == value for i, value in zip(net.load_nodes, power)), name='fixed_p')
            eta = model.addVar(name='violation')
            operation = equations.add_operation(model, choice, p, eta)
            # 正比例缩放目标，避免边界的微小标幺违反量淹没在绝对求解误差中。
            model.setObjective(1000.*eta)
            model.Params.TimeLimit = max(0., deadline-perf_counter())
            model.optimize()
            if not model.SolCount:
                return dict(cut=None, state=None, feasible=False)
            # eta=0 对应原约束；松弛量与求解误差共用接受容差，直接返回原始状态。
            if np.maximum(0., eta.X)+model.MaxVio <= PLANNING_TOL:
                return dict(cut=None, state=operation.state.X, feasible=True)
            # ||tail|| <= head 的支撑平面对整个锥有效；正 eta 仍须用对偶取得有效割。
            planes = []
            for head, tail in operation.cones:
                values = np.array([item.getValue() for item in tail])
                length = np.linalg.norm(values)
                if length:
                    planes.append(head-gp.quicksum(float(a/length)*item for a, item in zip(values, tail)))
            # 未取得可行证书则继续分离；用 LP 乘子形成割，不依赖 QCP 对偶恢复。
            model.remove(model.getQConstrs())
            for plane in planes:
                model.addConstr(plane >= 0., name='cone_support')
            model.Params.TimeLimit = max(0., deadline-perf_counter())
            model.optimize()
            if model.Status != GRB.OPTIMAL or model.ObjVal <= 0.:
                return dict(cut=None, state=None, feasible=False)
            cut = self._separating_cut(model, operation, choice, p, [*fixed_x.values(), *fixed_p.values()], x, power)
            return dict(cut=cut, state=None, feasible=False)

    def _separating_cut(self, model, operation, choice, power_vars, fixed, x, power):
        """按行方向组合必要约束，以运行变量全局盒消去 y，得到 α+βᵀp+δᵀx>=0。"""
        rows = model.getConstrs()
        dual = np.array(model.getAttr('Pi', rows))
        senses = np.array(model.getAttr('Sense', rows))
        dual[senses == '>'] = np.maximum(dual[senses == '>'], 0.)
        dual[senses == '<'] = np.minimum(dual[senses == '<'], 0.)
        dual[[row.index for row in fixed]] = 0.
        # Gurobi 导出系数仅用于对偶代数；MP/SP 的物理建模均为具名约束。
        coefficients = np.asarray(model.getA().T@dual).ravel()
        state_ids = [var.index for var in operation.state.tolist()]
        h = coefficients[state_ids]
        equations = self.equations
        constant = (-dual@np.array(model.getAttr('RHS', rows))
                    +np.maximum(h, 0.)@equations.y_ub_global
                    +np.minimum(h, 0.)@equations.y_lb_global)
        cut = np.r_[constant, coefficients[[v.index for v in power_vars.values()]],
                     coefficients[[v.index for v in choice.values()]]]
        scale = max(np.max(np.abs(cut)), np.max(np.abs(cut[1:1+len(power)]))*equations.network.base)
        cut /= scale
        cut[0] += 1e-10
        if cut[0]+cut[1:1+len(power)]@power+cut[1+len(power):]@x < -1e-9:
            return cut
        return None


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
        d = len(bounds)
        for cut in cuts:
            cut = np.asarray(cut)
            scale = max(np.max(np.abs(np.r_[cut[0], cut[1:1+d]*bounds, cut[1+d:]])), 1e-20)
            problem.add_cut(cut*self.distance_scale/scale)
        # 等价缩放，避免 1e-8 的覆盖阈值与求解器绝对容差处于相近量级。
        delta = m.addVar(lb=-4.*self.distance_scale, ub=4.*self.distance_scale, name='uncovered_distance')
        for k, eq in enumerate(inner_halfspaces):
            select = m.addVars(len(eq), vtype=GRB.BINARY, name=f'outside_{k}')
            m.addConstr(select.sum() == 1)
            for f, face in enumerate(eq):
                a = face[:d]*(1-tau)
                # xi=p/bounds in [0,1]^3；M 由该盒的精确下界推出。
                big_m = 4.-face[d]-np.minimum(a, 0.).sum()
                expression = gp.quicksum(float(a[j]/bounds[j])*problem.power[j].item() for j in range(d))
                m.addConstr(delta <= self.distance_scale*(expression+float(face[d])+float(big_m)*(1-select[f])))
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
        x, point = np.rint(problem.x.X).astype(int), problem.power.X
        if self.mode == 'physical':
            details.update(feasible=m.MaxVio <= PLANNING_TOL)
        return dict(complete=False, bound=bound, x=x, p=point, **details)


def evaluation_bounds(network, *, threads=DEFAULT_SOLVER_THREADS):
    """三个无预算 LP 轴向全局上界确定公共评价箱。"""
    equations, bounds = PlanningEquations(network, 'linear'), []
    for axis in np.eye(len(network.load_nodes)):
        problem = PlanningModel(equations, threads=threads)
        for e in network.corridors:
            for k in e.types:
                if any(t.r <= k.r and t.reactance <= k.reactance and t.capacity >= k.capacity
                       and (t.r < k.r or t.reactance < k.reactance or t.capacity > k.capacity) for t in e.types):
                    problem.choices[e.id, k.id].UB = 0.
        problem.power.UB = axis*network.power_limit
        with problem.model:
            answer = problem.solve()
        if answer is None or answer['bound'] is None or not np.isfinite(answer['bound']):
            raise RuntimeError('No finite planning bound for the common evaluation box')
        bounds.append(min(network.power_limit, answer['bound']))
    return np.ceil(np.asarray(bounds)/10.)*10.
