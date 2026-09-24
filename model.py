"""Gurobi 配电网规划：e 为走廊 ID，k 为该走廊的型号 ID，i 为节点 ID。

支持问题与退让问题共用逐节点/逐走廊约束；x[e,k] 选型，P/Q/ell[e,k] 为标幺运行量，
p[i] 为 kW 负荷，v[i] 为电压平方。数组仅用于结果、几何算法和割的传递。
固定符号、变量名、单位和索引约定见 docs/notation.md；同义量不随重构改名。
"""
from types import SimpleNamespace

import gurobipy as gp
import numpy as np
from gurobipy import GRB

DEFAULT_SOLVER_THREADS = 20
MP_TIME_LIMIT = 20.
PLANNING_TOL = 1e-8
PAD_KW = 1e-4  # 全局界向域外补偿 (kW)。


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

    def add_operation(self, model, x, p):
        """运行约束 (6)—(11)：所有目标下保持原物理约束，不引入违反量 eta。"""
        net = self.network
        P = model.addVars(self.keys, lb=self.pmin, ub=self.pmax, name='P')
        Q = model.addVars(self.keys, lb=self.qmin, ub=self.qmax, name='Q')
        ell = model.addVars(self.keys, ub=self.ellmax, name='ell')
        v = model.addVars(net.nodes, ub={i: self.vmax[i] for i in net.nodes}, name='v')
        v[net.root] = 1.
        plus = model.addVars(self.types, ub=self.drop_max, name='drop_plus')
        minus = model.addVars(self.types, ub=self.drop_max, name='drop_minus')
        z = {e: gp.quicksum(x[e, k] for k in self.types[e]) for e in self.types}

        fixed_p, fixed_q = dict(zip(net.nodes, net.fixed_p)), dict(zip(net.nodes, net.fixed_q))
        ratio = dict(zip(net.load_nodes, net.q_ratio))
        for i in net.nodes:
            # (6) 节点有功/无功守恒：流入送端功率扣除线损，再减流出功率，等于负荷。
            incoming_p = gp.quicksum(P[e, k]-self.r[e, k]*ell[e, k]
                                    for e in self.incoming[i] for k in self.types[e])
            incoming_q = gp.quicksum(Q[e, k]-self.reactance[e, k]*ell[e, k]
                                    for e in self.incoming[i] for k in self.types[e])
            outgoing_p = gp.quicksum(P[e, k] for e in self.outgoing[i] for k in self.types[e])
            outgoing_q = gp.quicksum(Q[e, k] for e in self.outgoing[i] for k in self.types[e])
            model.addConstr(incoming_p-outgoing_p == (fixed_p[i]+p.get(i, 0.))/net.base, name=f'active_balance[{i}]')
            model.addConstr(incoming_q-outgoing_q == (fixed_q[i]+ratio.get(i, 0.)*p.get(i, 0.))/net.base,
                            name=f'reactive_balance[{i}]')
            # (7) vmin <= v_i <= vmax：节点电压平方下限；上限已在变量声明中给出。
            model.addConstr(v[i] >= self.vmin[i], name=f'voltage_min[{i}]')

        for e, (i, j) in self.ends.items():
            # (8) 支路压降：v_j-v_i + 2(rP+χQ) - (r²+χ²)ell = s⁺-s⁻。
            drop = v[j]-v[i]+gp.quicksum(
                2*(self.r[e, k]*P[e, k]+self.reactance[e, k]*Q[e, k])
                -(self.r[e, k]**2+self.reactance[e, k]**2)*ell[e, k] for k in self.types[e])
            model.addConstr(drop == plus[e]-minus[e], name=f'voltage_drop[{e}]')
            # 接通时压降松弛为零；开断时允许两端电压不同。
            model.addConstr(plus[e] <= self.drop_max[e]*(1-z[e]), name=f'drop_plus_bound[{e}]')
            model.addConstr(minus[e] <= self.drop_max[e]*(1-z[e]), name=f'drop_minus_bound[{e}]')
            for k in self.types[e]:
                # (9) 未选型号 P=Q=ell=0；非根走廊允许参考方向的反向潮流。
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
                    # (10) 支路电流锥：P²+Q² <= v_i*ell，固定 x 后为凸运行模型。
                    model.addQConstr(P[e, k]**2+Q[e, k]**2 <= v[i]*ell[e, k], name=f'current_cone[{e},{k}]')

        source_p = gp.quicksum(P[e, k] for e in self.outgoing[net.root] for k in self.types[e])
        source_q = gp.quicksum(Q[e, k] for e in self.outgoing[net.root] for k in self.types[e])
        # (11) 电源注入含全网损耗：有功、无功以及 SOCP 的视在功率限额。
        for flow, limit, name in ((source_p, net.source_pmax, 'source_P'), (source_q, net.source_qmax, 'source_Q')):
            if np.isfinite(limit):
                model.addConstr(flow <= limit, name=name)
        if self.method == 'socp' and np.isfinite(net.source_smax):
            model.addQConstr(source_p**2+source_q**2 <= net.source_smax**2, name='source_capacity')
        state = gp.MVar.fromlist([*P.values(), *Q.values(), *ell.values(),
                                 *(v[i] for i in net.nodes), *plus.values(), *minus.values()])
        return SimpleNamespace(P=P, Q=Q, ell=ell, v=v, plus=plus, minus=minus, state=state)


class PlanningModel:
    """统一 Gurobi 模型；solve 只切换目标和查询约束，不嵌套求解。"""

    def __init__(self, equations, *, power=None, budget=np.inf, fixed_plan=None,
                 strengthen=False, threads=DEFAULT_SOLVER_THREADS):
        net = equations.network
        model = new_model('planning_'+equations.method, threads)
        model.Params.NumericFocus = 2

        # 1. x[e,k]∈{0,1}, p[i]≥0 (kW), a[i]∈{0,1}, f[e]∈[-n,n]。
        x = model.addVars(equations.keys, vtype=GRB.BINARY, name='x')
        p = model.addVars(net.load_nodes, name='p_kw')
        a = model.addVars(net.nodes, lb=dict(zip(net.nodes, net.required.astype(float))),
                          ub=1., vtype=GRB.BINARY, name='active')
        f = model.addVars(equations.types, lb=-net.n, ub=net.n, name='connectivity')
        z = {e: gp.quicksum(x[e, k] for k in equations.types[e]) for e in equations.types}

        # 2. 型号互斥及端点接入：z_e=Σ_k x_ek≤1, z_e≤a_i,a_j。
        for e, (i, j) in equations.ends.items():
            model.addConstr(z[e] <= 1., name=f'one_type[{e}]')
            model.addConstr(z[e] <= a[j], name=f'receiver_active[{e}]')
            if i != net.root:
                model.addConstr(z[e] <= a[i], name=f'sender_active[{e}]')

        # 3. 径向连通：Gf=a, -nz≤f≤nz, Σ_e z_e=Σ_i a_i。
        for e in equations.types:
            model.addConstr(f[e] <= net.n*z[e], name=f'connectivity_max[{e}]')
            model.addConstr(f[e] >= -net.n*z[e], name=f'connectivity_min[{e}]')
        for i in net.nodes:
            model.addConstr(gp.quicksum(f[e] for e in equations.incoming[i])-
                            gp.quicksum(f[e] for e in equations.outgoing[i]) == a[i],
                            name=f'connectivity_balance[{i}]')
        model.addConstr(gp.quicksum(z.values()) == gp.quicksum(a.values()), name='tree_edges')

        # 4. cᵀx≤budget, 1ᵀp≤power_limit；S_b=network.base。
        investment = gp.quicksum(equations.cost[e, k]*x[e, k] for e, k in equations.keys)
        model.addConstr(gp.quicksum(p.values())/net.base <= net.power_limit/net.base, name='power_limit')
        if np.isfinite(budget):
            model.addConstr(investment <= budget, name='budget')

        # 5. 独立点校核/固定方案：power 用 kW；fixed_plan 为走廊→型号字典。
        if power is not None:
            model.addConstrs((p[i] == value for i, value in zip(net.load_nodes, power)), name='fixed_power')
        if fixed_plan is not None:
            for key, value in zip(equations.keys, net.encode_plan(fixed_plan)):
                x[key].LB = x[key].UB = value

        # 6—11. 节点守恒、电压、容量、SOCP、电源限额，编号与 add_operation 对应。
        operation = equations.add_operation(model, x, p)

        # 12. 树松弛加强：b_ij+b_ji=z_e, Σ_j b_ji=1, 0≤F^k_ij≤b_ij。
        # 全接入纯负荷树适用；连续父弧和多商品流保留每个可行整数树。
        if strengthen:
            assert np.all(net.required) and np.all(net.fixed_p >= 0.) and np.all(net.fixed_q >= 0.)
            assert np.all(net.q_ratio >= 0.) and np.all(net.r > 0.) and np.all(net.reactance >= 0.)
            arcs = [arc for ends in equations.ends.values() for arc in (ends, ends[::-1])]
            incoming = {i: [arc for arc in arcs if arc[1] == i] for i in (net.root, *net.nodes)}
            outgoing = {i: [arc for arc in arcs if arc[0] == i] for i in (net.root, *net.nodes)}
            parent_arc = model.addVars(arcs, ub=1., name='parent_arc')
            for arc in incoming[net.root]:
                parent_arc[arc].UB = 0.
            for e, (i, j) in equations.ends.items():
                model.addConstr(parent_arc[i, j]+parent_arc[j, i] == z[e])
            for i in net.nodes:
                model.addConstr(gp.quicksum(parent_arc[arc] for arc in incoming[i]) == 1.)
                operation.v[i].UB = min(1., equations.vmax[i])
            commodity_flow = model.addVars([(i, j, k) for i, j in arcs for k in net.nodes],
                                           ub=1., name='commodity_flow')
            for k in net.nodes:
                for i, j in arcs:
                    model.addConstr(commodity_flow[i, j, k] <= parent_arc[i, j])
                for i in net.nodes:
                    model.addConstr(gp.quicksum(commodity_flow[a, b, k] for a, b in incoming[i])-
                                    gp.quicksum(commodity_flow[a, b, k] for a, b in outgoing[i]) == float(i == k))
            # 正向 P≥r·ell，逆向参考 P≤0；Q 同理。v≤1 是纯负荷树的有效界。
            for e, (i, j) in equations.ends.items():
                for k in equations.types[e]:
                    for flow, loss, upper in ((operation.P[e, k], equations.r[e, k], equations.pmax[e, k]),
                                             (operation.Q[e, k], equations.reactance[e, k], equations.qmax[e, k])):
                        model.addConstr(flow <= upper*parent_arc[i, j])
                        model.addConstr(flow >= loss*operation.ell[e, k]-upper*parent_arc[j, i])

        # 13. d / rho 为负荷退让量 (kW)；支持查询时固定为零。
        self.distance_kw = model.addVar(ub=0., name='load_distance_kw')
        self.projected_deficit_kw = model.addVar(ub=0., name='projected_deficit_kw')
        self.query_rows = []
        self.equations, self.model = equations, model
        self.choices, self.loads = x, p
        self.x = gp.MVar.fromlist(list(x.values()))
        self.power = gp.MVar.fromlist(list(p.values()))
        self.active_nodes = gp.MVar.fromlist(list(a.values()))
        self.operation, self.state = operation, operation.state
        self.investment, self.fixed_power = investment, power is not None
        self.objective_scale = 1. if self.fixed_power else net.base
        model.setObjective(investment if self.fixed_power else gp.quicksum(p.values())/net.base,
                           GRB.MINIMIZE if self.fixed_power else GRB.MAXIMIZE)

    def exclude(self, x):
        """独立 AC 校核排除已证 AC 不可行的方案：至少一个 x 分量不同。"""
        self.model.addConstr((1-2*x)@self.x >= 1-x.sum())

    def solve(self, time_limit=MP_TIME_LIMIT, *, weights=None, target=None, cut_normals=None,
              incumbent=None, start=None, gap_kw=1e-3):
        """support: max omegaᵀp；distance: min d；projected_distance: min rho。"""
        model, net = self.model, self.equations.network
        d = len(net.load_nodes)

        # 1. 更新查询；物理模型和已加入支持割保持不变。
        model.remove(self.query_rows)
        self.query_rows = []
        self.distance_kw.UB = self.projected_deficit_kw.UB = 0.
        mode = 'cost' if self.fixed_power else 'support'
        if target is not None:
            mode = 'distance' if cut_normals is None else 'projected_distance'
            normals = np.eye(d) if cut_normals is None else np.asarray(cut_normals)
            deficit = self.distance_kw if cut_normals is None else self.projected_deficit_kw
            deficit.UB = GRB.INFINITY
            # Wp + rho·1 ≥ Wq；W=I 时 d=min max_i(q_i-p_i)_+。
            self.query_rows = [model.addConstr(
                (gp.quicksum(float(w)*p for w, p in zip(normal, self.loads.values()))+deficit)/net.base
                >= float(normal @ target)/net.base) for normal in normals]
            objective, sense, scale = deficit/net.base, GRB.MINIMIZE, net.base
        elif self.fixed_power:
            objective, sense, scale = self.investment, GRB.MINIMIZE, 1.
        else:
            weights = np.ones(d) if weights is None else np.asarray(weights)
            objective = gp.quicksum(float(w)*p for w, p in zip(weights, self.loads.values()))/net.base
            sense, scale = GRB.MAXIMIZE, net.base
        self.objective_scale = scale
        model.setObjective(objective, sense)

        # 2. Start 只加速求解，不限制成本；所有预算内方案仍可参与竞争。
        if incumbent is not None:
            self.x.Start, self.power.Start, self.state.Start = incumbent['x'], incumbent['p'], incumbent['state']
        elif start is not None:
            self.x.Start = start
        model.Params.TimeLimit = max(0., time_limit)
        model.Params.MIPGapAbs = (1e-7 if mode == 'cost' else gap_kw)/scale
        model.optimize()

        # 3. ObjBound 才是全局证据；超时不等于不可行，无 incumbent 也可能有有效界。
        if model.Status == GRB.INFEASIBLE:
            return None
        bound = float(model.ObjBound*scale)
        if mode != 'cost':
            bound += PAD_KW if sense == GRB.MAXIMIZE else -PAD_KW
        bound = bound if np.isfinite(bound) else None
        answer = dict(mode=mode, status='unknown', x=None, p=None, state=None,
                      objective=None, bound=bound, feasible=False)
        if model.SolCount:
            answer.update(x=np.rint(self.x.X).astype(int), p=self.power.X, state=self.state.X,
                          objective=float(model.ObjVal*scale), feasible=model.MaxVio <= PLANNING_TOL)
            if answer['feasible']:
                answer['status'] = 'optimal' if model.Status == GRB.OPTIMAL else 'feasible'
        return answer
