"""Gurobi 配电网规划：e 为走廊 ID，k 为该走廊的型号 ID，i 为节点 ID。

MP 与 SP 共用逐节点/逐走廊约束；x[e,k] 选型，P/Q/ell[e,k] 为标幺运行量，
p[i] 为 kW 负荷，v[i] 为电压平方。数组仅用于结果、几何算法和割的传递。
固定符号、变量名、单位和索引约定见 docs/notation.md；同义量不随重构改名。
"""
from concurrent.futures import ThreadPoolExecutor
from functools import cache
from time import perf_counter
from types import SimpleNamespace

import clarabel
from scipy import sparse

import gurobipy as gp
import numpy as np
from gurobipy import GRB

DEFAULT_SOLVER_THREADS = 20
MP_TIME_LIMIT = 20.
SP_TIME_LIMIT = {'linear': 5., 'socp': 10.}
PLANNING_TOL = 1e-8
CONE_CONV_TOL = 1e-6     # x 自由锥 MISOCP 的 barrier 收敛容差：带全部已紧化方案的提升行时更严的判据使节点松弛数值失败
OBBT_ROUNDS = 2          # 每个方案的 OBBT 轮数
OBBT_PAD = 1e-6          # OBBT 端点外扩（标幺），吸收最优值误差
ENVELOPE_MARGIN = 1e-4   # 反向锥包络割右端裕量（标幺²）：AC 点严格在内，且使锥与包络之间的可行壳不致过薄
OBBT_CONV_TOL = 1e-8     # OBBT 与紧化模型的 barrier 收敛容差：1e-10 时极值点和贴锥可行集上 barrier 偶发失败
OBBT_PARAMS = ('Threads', 'FeasibilityTol', 'OptimalityTol', 'BarQCPConvTol', 'DualReductions',
               'BarHomogeneous', 'Aggregate', 'ScaleFlag')   # Model.copy 不保留参数，副本逐项复制


def new_model(name, threads):
    model = gp.Model(name)
    model.Params.OutputFlag, model.Params.Threads = 0, threads
    model.Params.MIPGap = model.Params.MIPGapAbs = 0.
    model.Params.FeasibilityTol = model.Params.IntFeasTol = model.Params.OptimalityTol = 1e-9
    model.Params.BarQCPConvTol = 1e-10
    model.Params.DualReductions, model.Params.NonConvex = 0, 0
    return model


@cache
def obbt_pool(workers):
    """OBBT 线程池与等量 Gurobi 环境；每个环境只被一个模型副本使用。"""
    return ThreadPoolExecutor(workers), [gp.Env(params={'OutputFlag': 0}) for _ in range(workers)]


def obbt_extremes(model, indices, workers):
    """各变量（按下标）的 (min, max)，未证得最优的端点为 None：复制 workers 份模型交错分块并行，求解时释放 GIL。"""
    pool, envs = obbt_pool(workers)
    copies = [model.copy(env=env) for env in envs]
    for copy in copies:
        for name in OBBT_PARAMS:
            copy.setParam(name, model.getParamInfo(name)[2])
    parts = list(pool.map(_extremes, copies, [indices[w::workers] for w in range(workers)]))
    for copy in copies:
        copy.dispose()
    ends = [None]*len(indices)
    for w, part in enumerate(parts):
        ends[w::workers] = part
    return ends


def _extremes(model, indices):
    variables, ends = model.getVars(), []
    for index in indices:
        values = []
        for sense in (GRB.MINIMIZE, GRB.MAXIMIZE):
            model.setObjective(variables[index], sense)
            model.optimize()
            values.append(model.ObjVal if model.Status == GRB.OPTIMAL else None)
        ends.append(values)
    return ends


class GridPhysics:
    """保存共用电网物理参数，并向 MP/SP 模型添加运行变量与约束。"""

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
        self.boxes = {}   # 方案 x 元组 → OBBT 盒 {(量名, 键): (下界, 上界)}，标幺
        self._build_variable_bounds()

    def _build_variable_bounds(self):
        net = self.network
        pmax = np.minimum(net.capacity, min(net.source_pmax, net.source_smax))
        qmax = np.full(net.n_types, min(net.source_qmax, net.source_smax))
        ellmax = np.minimum(pmax/net.r, net.ell_limit) if self.method == 'socp' else np.zeros(net.n_types)
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

    def obbt(self, x, bounds, workers):
        """固定方案 x 的 OBBT：OBBT_ROUNDS 轮，每轮在已有盒与包络下求各运行量的 min/max，外扩 OBBT_PAD 后
        与上一轮（首轮为全局界）取交；未证得最优的端点保留原界。独立功率在 [0, bounds] kW 内自由。"""
        scheme = tuple(int(v) for v in x)
        if scheme in self.boxes:
            return
        net = self.network
        names = [(name, key) for key, bit in zip(self.keys, scheme) if bit for name in ('P', 'Q', 'ell')]
        names += [('v', i) for i in net.nodes]
        for _ in range(OBBT_ROUNDS):
            with new_model('obbt', 1) as model:
                model.Params.BarQCPConvTol = OBBT_CONV_TOL
                power = model.addVars(net.load_nodes, ub=dict(zip(net.load_nodes, bounds)), name='p_kw')
                operation = self.add_operation(model, dict(zip(self.keys, scheme)), power, scheme=scheme)
                model.update()
                ends = obbt_extremes(model, [getattr(operation, name)[key].index for name, key in names], workers)
            previous, box = self.boxes.get(scheme), {}
            for item, (low, high) in zip(names, ends):
                old = previous[item] if previous else self._limits(*item)
                box[item] = (old[0] if low is None else max(low-OBBT_PAD, old[0]),
                             old[1] if high is None else min(high+OBBT_PAD, old[1]))
            self.boxes[scheme] = box

    def _limits(self, name, key):
        """运行量对任意 x 成立的全局界，用于提升系数 M。"""
        low, high = {'P': (self.pmin, self.pmax), 'Q': (self.qmin, self.qmax),
                     'ell': (None, self.ellmax), 'v': (self.vmin, self.vmax)}[name]
        return (0. if low is None else low[key]), high[key]

    def _tighten(self, add, x, operation, scheme):
        """方案的盒约束与反向锥包络割 u_L*ell+ell_L*u-u_L*ell_L <= (P_L+P_U)P-P_L*P_U+(Q_L+Q_U)Q-Q_L*Q_U+裕量
        （另一条取 u_U、ell_U）；每行右端加 M*H(x)，H 为 x 与该方案的汉明距离，M 为该行在全局界上的最大违反量。
        add 为加约束的方法：建模时 model.addConstr，行生成回调中 model.cbLazy。"""
        box = self.boxes[scheme]
        hamming = gp.quicksum(1-x[key] if bit else x[key] for key, bit in zip(self.keys, scheme))

        def row(terms, rhs):
            worst = sum(max(c*low, c*high) for c, _, (low, high) in terms)-rhs
            add(gp.quicksum(c*w for c, w, _ in terms) <= rhs+max(worst, 0.)*hamming)
        for (name, key), (low, high) in box.items():
            row([(-1., getattr(operation, name)[key], self._limits(name, key))], -low)
            row([(1., getattr(operation, name)[key], self._limits(name, key))], high)
        for key, bit in zip(self.keys, scheme):
            if bit:
                sender = self.ends[key[0]][0]
                u, (ul, uu), u_limits = ((1., (1., 1.), (1., 1.)) if sender == self.network.root else
                                         (operation.v[sender], box['v', sender], self._limits('v', sender)))
                (pl, pu), (ql, qu), (ll, lu) = box['P', key], box['Q', key], box['ell', key]
                for a, b in ((ul, ll), (uu, lu)):
                    row([(a, operation.ell[key], self._limits('ell', key)), (b, u, u_limits),
                         (-(pl+pu), operation.P[key], self._limits('P', key)),
                         (-(ql+qu), operation.Q[key], self._limits('Q', key))], a*b-pl*pu-ql*qu+ENVELOPE_MARGIN)

    def add_operation(self, model, x, p, eta=0., scheme=None):
        """MP 用 eta=0；SP 只放松功率平衡与压降等式，容量、电压界和锥保持原式；已紧化的 scheme 追加其紧化行。"""
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
        operation = SimpleNamespace(P=P, Q=Q, ell=ell, v=v, plus=plus, minus=minus, state=state, cones=cones)
        if scheme is not None and tuple(scheme) in self.boxes:
            model.Params.BarQCPConvTol = OBBT_CONV_TOL
            self._tighten(model.addConstr, x, operation, tuple(scheme))
        return operation


class MasterProblem:
    """完整 MP：给定 power/min_total 时最小投资，否则最大化 direction@p。"""

    def __init__(self, equations, *, power=None, budget=np.inf,
                 min_total=None, fixed_plan=None, direction=None, threads=DEFAULT_SOLVER_THREADS):
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
            if not corridor.switchable:
                for kind in corridor.types:
                    model.addConstr(x[corridor.id, kind.id] == int(
                        corridor.initial_active and kind.id == corridor.existing_type),
                        name=f'fixed_branch[{corridor.id},{kind.id}]')
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
        investment = net.cost_offset+gp.quicksum(equations.cost[e, k]*x[e, k] for e, k in equations.keys)
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

        # 7. 所有 MP 均包含完整运行变量及物理约束；固定方案时含该方案的紧化行。
        operation = equations.add_operation(model, x, p,
                                            scheme=None if fixed_plan is None else net.encode_plan(fixed_plan))

        # 8. MP1 最小投资；MP2 沿给定方向最大化负荷，其他负荷仍自由。
        minimizing = power is not None or min_total is not None
        self.direction = np.ones(len(net.load_nodes)) if direction is None else np.asarray(direction, dtype=float)
        objective = gp.quicksum(w*value for w, value in zip(self.direction, p.values()))/net.base
        model.setObjective(investment if minimizing else objective, GRB.MINIMIZE if minimizing else GRB.MAXIMIZE)
        self.objective_scale = 1. if minimizing else net.base

        # 9. 保存变量引用；具名变量与扁平视图共享同一批 Gurobi 变量
        self.equations, self.model = equations, model
        self.choices, self.loads = x, p
        self.x = gp.MVar.fromlist(list(x.values()))
        self.power = gp.MVar.fromlist(list(p.values()))
        self.active_nodes = gp.MVar.fromlist(list(a.values()))
        self.operation = operation
        self.state = operation.state

    def use_incumbent(self, incumbent):
        """提供已有可行解作为起点，并将其投资作为成本上界。"""
        net = self.equations.network
        cost = float(net.cost_offset+net.cost@incumbent['x'])
        self.x.Start, self.power.Start = incumbent['x'], incumbent['p']
        self.state.Start = incumbent['state']
        self.model.addConstr(net.cost_offset+net.cost@self.x <= cost+1e-9)
        return cost

    def exclude(self, x):
        """要求至少一个型号选择不同，排除指定建设方案。"""
        self.model.addConstr((1-2*x)@self.x >= 1-x.sum())

    def solve(self, time_limit=MP_TIME_LIMIT, *, incumbent=None, start=None, radial_gap_kw=1e-3, tolerance=PLANNING_TOL):
        """直接求 MP1/MP2；初始解、目标间隙及运行证书在此处理，不调用 SP。解的 MaxVio 超过 tolerance 即报错。"""
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
        if model.Status == GRB.TIME_LIMIT:
            raise TimeoutError(f'{model.ModelName}: time limit')
        if model.Status not in (GRB.OPTIMAL, GRB.USER_OBJ_LIMIT) or not model.SolCount:
            raise RuntimeError(f'{model.ModelName}: status={model.Status}, SolCount={model.SolCount}')
        if model.MaxVio > tolerance:
            raise RuntimeError(f'{model.ModelName}: MaxVio={model.MaxVio:g} > {tolerance:g}')
        bound = model.ObjBound*self.objective_scale

        # 3. 提取选型、负荷和目标值：投资用原费用单位，负荷用 kW
        x = np.rint(self.x.X).astype(int)
        p = self.power.X
        objective = equations.network.cost_offset+equations.network.cost@x if model.ModelSense == GRB.MINIMIZE else self.direction@p

        # 4. 原样读取完整运行证书；质量指标覆盖全部物理约束。
        state = self.state.X
        feasible = True

        # 5. 返回候选解、目标界和认证状态
        minimizing = model.ModelSense == GRB.MINIMIZE
        gap = objective-bound if minimizing else bound-objective
        tolerance = 1e-7 if minimizing else radial_gap_kw+1e-5
        status = 'optimal' if model.Status == GRB.OPTIMAL or gap <= tolerance else 'feasible'
        return dict(x=x, p=p, objective=float(objective), bound=bound, state=state, feasible=feasible, status=status)


class SubProblem:
    """固定 x,p 的运行可行性 SP；锥支撑平面的 LP 对偶生成全局联合割。"""

    def __init__(self, equations, *, threads=DEFAULT_SOLVER_THREADS, numeric_focus=0):
        self.equations, self.threads, self.calls = equations, threads, 0
        self.numeric_focus = numeric_focus
        self.cut_calls = 0

    def _build(self, model, x, power):
        """评分 SOCP 与取割 LP 使用同一组变量、等式和物理约束。"""
        equations, net = self.equations, self.equations.network
        model.Params.Aggregate = 0
        model.Params.ScaleFlag = 0
        model.Params.BarQCPConvTol = 1e-9
        model.Params.BarHomogeneous = 1 if self.numeric_focus == 3 else -1
        model.Params.NumericFocus = self.numeric_focus
        choice = model.addVars(equations.keys, ub=1., name='x')
        p = model.addVars(net.load_nodes, lb=-GRB.INFINITY, name='p_kw')
        # 固定参数用等式表示；求割时去掉这些等式的乘子，保留 x,p 系数。
        fixed_x = model.addConstrs((choice[e, k] == value for (e, k), value in zip(equations.keys, x)), name='fixed_x')
        fixed_p = model.addConstrs((p[i] == value for i, value in zip(net.load_nodes, power)), name='fixed_p')
        eta = model.addVar(name='violation')
        operation = equations.add_operation(model, choice, p, eta, scheme=x)
        model.setObjective(eta)
        return choice, p, eta, operation, [*fixed_x.values(), *fixed_p.values()]

    def solve(self, x, power, time_limit=None, *, score_only=False):
        # 1. 固定 x、p，建立等式松弛问题
        self.calls += 1
        equations = self.equations
        limit = SP_TIME_LIMIT[equations.method] if time_limit is None else time_limit
        deadline = perf_counter()+limit
        with new_model('planning_SP', self.threads) as model:
            choice, p, eta, operation, fixed = self._build(model, x, power)
            # 2. 求解并检查终止状态与精度；失败直接报错。
            model.Params.TimeLimit = max(0., deadline-perf_counter())
            model.optimize()
            if model.Status == GRB.TIME_LIMIT:
                raise TimeoutError('SP: time limit')
            if model.Status != GRB.OPTIMAL:
                raise RuntimeError(f'SP {equations.method}: status={model.Status}, x={np.asarray(x).tolist()}, p={np.asarray(power).tolist()}')
            if model.MaxVio > PLANNING_TOL:
                raise RuntimeError(f'SP {equations.method}: MaxVio={model.MaxVio:g} > {PLANNING_TOL:g}, p={np.asarray(power).tolist()}')
            # 3. 保存原始最小违反量供顶点评分；不能用后续切平面 LP 的值替代。
            value = float(eta.X)
            if np.maximum(0., eta.X)+model.MaxVio <= PLANNING_TOL:
                return dict(cut=None, state=operation.state.X, feasible=True, eta=value)
            if eta.X <= PLANNING_TOL:
                raise RuntimeError(f'SP {equations.method}: eta={eta.X:g}, MaxVio={model.MaxVio:g}; certificate exceeds tolerance')
            # 4. 缓存数值支撑方向；只评分时不解 LP，也不保留求解器对象。
            cone_normals = []
            for head, tail in operation.cones:
                values = np.array([item.getValue() for item in tail])
                length = np.linalg.norm(values)
                cone_normals.append(values/length if length else values)
            if score_only:
                return dict(cut=None, state=None, feasible=False, eta=value, cone_normals=cone_normals)
            cut = self._cut(model, operation, choice, p, fixed, x, power, cone_normals, deadline)
            return dict(cut=cut, state=None, feasible=False, eta=value)

    def generate_cut(self, x, power, cone_normals, time_limit=None):
        """仅对获选点解一次支撑 LP；方向来自该 x,p 已完成的 SOCP。"""
        limit = SP_TIME_LIMIT[self.equations.method] if time_limit is None else time_limit
        deadline = perf_counter()+limit
        with new_model('planning_SP_cut', self.threads) as model:
            choice, p, eta, operation, fixed = self._build(model, x, power)
            return self._cut(model, operation, choice, p, fixed, x, power, cone_normals, deadline)

    def _cut(self, model, operation, choice, p, fixed, x, power, cone_normals, deadline):
        model.update()
        model.remove(model.getQConstrs())
        for (head, tail), normal in zip(operation.cones, cone_normals):
            if np.any(normal):
                model.addConstr(head-gp.quicksum(float(a)*item for a, item in zip(normal, tail)) >= 0., name='cone_support')
        model.Params.TimeLimit = max(0., deadline-perf_counter())
        self.cut_calls += 1
        model.optimize()
        if model.Status == GRB.TIME_LIMIT:
            raise TimeoutError('SP cut LP: time limit')
        if model.Status != GRB.OPTIMAL:
            raise RuntimeError(f'SP cut LP: status={model.Status}, p={np.asarray(power).tolist()}')
        if model.ObjVal <= 0. or model.MaxVio > PLANNING_TOL:
            raise RuntimeError(f'SP cut LP: objective={model.ObjVal:g}, MaxVio={model.MaxVio:g}')
        return self._separating_cut(model, operation, choice, p, fixed, x, power)

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
        if not cut[0]+cut[1:1+len(power)]@power+cut[1+len(power):]@x < -1e-9:
            raise RuntimeError(f'SP cut does not separate the candidate: x={np.asarray(x).tolist()}, p={np.asarray(power).tolist()}')
        return cut


def cone_misocp(equations, budget, bounds, objective, *, tighten, rows=None, exclude=(), time_limit, mip_gap=0.,
                threads=DEFAULT_SOLVER_THREADS):
    """x 自由的锥 MISOCP（行生成）：在 xi=p/bounds 上 max objective@xi，一次分支定界返回全局上界与现任解。

    分支定界自己选网架，不枚举方案。模型不预先带紧化行：每找到一个新的现任网架，若其盒约束与反向锥包络行尚不在模型中，
    先调用 tighten(x, 现任功率 kW, 当前上界) 确保该网架已紧化，再把这些行（按汉明距离提升，对全部 x 有效）作为惰性约束
    加入，分支定界继续。没有行的网架为纯 SOCP，故 ObjBound 在任何终止状态下都是 R 在该锥上的有效上界；达到 mip_gap 时
    现任网架都带行。tighten 抛出的异常（分区到时、取消）在求解结束后重新抛出。
    """
    # 1. 完整 MP：全部 0/1 选型、拓扑、预算与运行约束，x 不固定；不预先带紧化行
    problem = MasterProblem(equations, budget=budget, threads=threads)
    with problem.model as model:
        choice = dict(zip(equations.keys, problem.x.tolist()))
        model.Params.BarQCPConvTol = CONE_CONV_TOL
        xi = [p*(1./float(b)) for p, b in zip(problem.loads.values(), bounds)]

        # 3. rows 给锥约束 rows@xi>=0 与分区盒 p<=bounds；否则为中心射线 xi_1=…=xi_d，并排除已知方案
        if rows is None:
            for j in range(1, len(xi)):
                model.addConstr(xi[0]-xi[j] == 0., name=f'radial_direction[{j}]')
            for scheme in exclude:
                problem.exclude(np.asarray(scheme))
        else:
            for j, row in enumerate(rows):
                model.addConstr(gp.quicksum(float(a)*v for a, v in zip(row, xi)) >= 0., name=f'cone[{j}]')
            problem.power.UB = bounds
        model.setObjective(gp.quicksum(float(a)*v for a, v in zip(objective, xi)), GRB.MAXIMIZE)

        # 3. 行生成：新现任网架尚无紧化行时先紧化，再把它的行作为惰性约束加入
        added, errors = set(), []
        variables, power = problem.x.tolist(), problem.power.tolist()

        def callback(solver, where):
            if where != GRB.Callback.MIPSOL or errors:
                return
            x = tuple(int(round(v)) for v in solver.cbGetSolution(variables))
            if x in added:
                return
            try:
                tighten(x, np.asarray(solver.cbGetSolution(power)), solver.cbGet(GRB.Callback.MIPSOL_OBJBND))
                equations._tighten(solver.cbLazy, choice, problem.operation, x)
                added.add(x)
            except BaseException as error:
                errors.append(error)
                solver.terminate()

        # 4. 解到相对间隙 mip_gap 或时限
        model.Params.MIPGap, model.Params.TimeLimit, model.Params.LazyConstraints = mip_gap, time_limit, 1
        model.optimize(callback)
        if errors:
            raise errors[0]
        if model.Status == GRB.INTERRUPTED:
            raise KeyboardInterrupt

        # 5. 上界（尚无界时为 inf）与现任解：方案元组、功率 kW
        try:
            bound = float(model.ObjBound)
        except gp.GurobiError:
            bound = np.inf
        found = model.SolCount > 0
        return dict(status=model.Status, bound=bound if abs(bound) < GRB.INFINITY else np.inf,
                    x=tuple(int(v) for v in np.rint(problem.x.X)) if found else None,
                    point=np.asarray(problem.power.X, float) if found else None)


LOAD_PF = .95
PV_PF = 1.
PV_Q_SIGN = 1.
RAY_CONE_MARGIN = 1e-6


def voltage_flow_bounds(network):
    """由压降等式和电流锥推导有效界；两端电压均有限，反送也有界。"""
    voltage = np.r_[network.vmax, 1.]
    sending = voltage[network.senders[network.type_corridor]]
    receiving = voltage[network.receivers[network.type_corridor]]
    ellmax = (np.sqrt(sending)+np.sqrt(receiving))**2/(network.r**2+network.reactance**2)
    ellmax = np.minimum(ellmax, network.ell_limit)
    return ellmax, np.sqrt(np.maximum(sending, receiving)*ellmax)


def port_bounds(network):
    """|节点净功率| 不超过相邻走廊各最大两端有功界之和。"""
    _, voltage_bound = voltage_flow_bounds(network)
    capacity = np.minimum(network.capacity, voltage_bound)
    return network.base*np.array([sum(capacity[block].max()
        for c, block in zip(network.corridors, network.type_slices) if node in c.endpoints)
        for node in network.load_nodes])


class PortPhysics(GridPhysics):
    """固定 sign 后用非负幅值装配原方程；真实接入功率 p=sign*u。"""

    def __init__(self, network, sign):
        self.sign = np.asarray(sign)
        network.q_ratio = np.where(self.sign > 0, np.tan(np.arccos(LOAD_PF)),
                                   PV_Q_SIGN*np.tan(np.arccos(PV_PF)))
        network.power_limit = float(port_bounds(network).sum())
        super().__init__(network, 'socp')

    def _build_variable_bounds(self):
        super()._build_variable_bounds()
        net = self.network
        lower_q = np.minimum(net.q_ratio*self.sign*port_bounds(net)/net.base, 0.).sum()
        pmax = net.capacity.copy()
        qmax = np.full(net.n_types, min(net.source_qmax, net.source_smax)-lower_q)
        ellmax = 2*pmax/net.r
        if not np.all(np.isfinite(net.capacity)):
            voltage_ell, _ = voltage_flow_bounds(net)
            ellmax = np.minimum(voltage_ell, (qmax-net.fixed_q.sum()/net.base)/net.reactance)
            voltage = np.r_[net.vmax, 1.]
            end_voltage = np.maximum(voltage[net.senders[net.type_corridor]],
                                     voltage[net.receivers[net.type_corridor]])
            pmax = np.minimum(pmax, np.sqrt(end_voltage*ellmax))
            qmax = np.minimum(qmax, pmax)
        ellmax = np.minimum(ellmax, net.ell_limit)
        self.pmax, self.qmax, self.ellmax = (dict(zip(self.keys, a)) for a in (pmax, qmax, ellmax))
        self.pmin, self.qmin = dict(zip(self.keys, -pmax)), dict(zip(self.keys, -qmax))
        self.y_lb_global = np.r_[-pmax, -qmax, np.zeros(net.n_types+net.n+2*net.n_corridors)]
        drop = list(self.drop_max.values())
        self.y_ub_global = np.r_[pmax, qmax, ellmax, net.vmax, drop, drop]

    def add_operation(self, model, x, p, eta=0., scheme=None):
        # 所有分区求解首次建模即采用同一数值设置，不作失败后的参数切换。
        model.Params.BarHomogeneous, model.Params.Aggregate = 1, 0
        model.Params.ScaleFlag = 1 if self.network.name == 'case33bw' else 0
        if isinstance(eta, gp.Var):
            net = self.network
            eta.UB = float(max(np.max(np.abs(net.fixed_p)), np.max(np.abs(net.fixed_q)),
                np.max(np.maximum(1., np.abs(net.q_ratio))*port_bounds(net)))/net.base)
        signed = {i: int(s)*p[i] for i, s in zip(self.network.load_nodes, self.sign)}
        operation = super().add_operation(model, x, signed, eta, scheme)
        for e in self.outgoing[self.network.root]:
            for k in self.types[e]:
                model.addConstr(-operation.P[e, k]+self.r[e, k]*operation.ell[e, k]
                                <= self.pmax[e, k]*x[e, k])
                model.addConstr(-operation.Q[e, k]+self.reactance[e, k]*operation.ell[e, k]
                                <= self.qmax[e, k]*x[e, k])
        return operation


def solve_conic(model, equations, operation, x, fixed, threads, deadline, tolerance, margin):
    """原 SOCP 的固定变量消元、连续求解及原始残差核验。"""
    model.update()
    variables, rows = model.getVars(), model.getConstrs()
    n = len(variables)

    def affine(expression):
        expression = gp.LinExpr(expression)
        coefficients = np.zeros(n)
        for j in range(expression.size()):
            coefficients[expression.getVar(j).index] += expression.getCoeff(j)
        return coefficients, expression.getConstant()

    # 1. 原始等式、不等式和变量界，统一写成 Az+s=b。
    matrix = model.getA().tocsc()
    rhs = np.asarray(model.getAttr('RHS', rows))
    senses = np.asarray(model.getAttr('Sense', rows))
    lb, ub = np.array(model.getAttr('LB', variables)), np.array(model.getAttr('UB', variables))
    for key, selected in zip(equations.keys, x):
        if not selected:
            fixed.update({flow[key].index: 0. for flow in (operation.P, operation.Q, operation.ell)})
    for e in equations.types:
        if sum(x[j] for j, key in enumerate(equations.keys) if key[0] == e):
            fixed.update({operation.plus[e].index: 0., operation.minus[e].index: 0.})
    free = np.array([j for j in range(n) if j not in fixed])
    solution = np.zeros(n)
    solution[list(fixed)] = list(fixed.values())
    reduced = matrix[:, free]
    keep = np.asarray(abs(reduced).sum(axis=1)).ravel() != 0.
    a, b = reduced[keep], (rhs-matrix@solution)[keep]
    equal = senses[keep] == '='
    direction = np.where(senses[keep] == '>', -1., 1.)
    eye = sparse.eye(len(free), format='csc')
    lower, upper = lb[free] > -GRB.INFINITY, ub[free] < GRB.INFINITY
    blocks = [a[equal], sparse.diags(direction[~equal])@a[~equal], -eye[lower], eye[upper]]
    values = [b[equal], direction[~equal]*b[~equal], -lb[free][lower],
              ub[free][upper]-margin*np.isin(free[upper], [operation.ell[key].index
                  for key, limit in zip(equations.keys, equations.network.ell_limit) if np.isfinite(limit)])]
    cones = [clarabel.ZeroConeT(int(equal.sum())),
             clarabel.NonnegativeConeT(int((~equal).sum()+lower.sum()+upper.sum()))]

    # 2. 每个原二阶锥保存为仿射向量，直接交给连续锥求解器。
    cone_rows = []
    for j, (head, tail) in enumerate(operation.cones):
        expressions = [affine(item) for item in (head, *tail)]
        a, b = np.array([item[0] for item in expressions]), np.array([item[1] for item in expressions])
        cone_rows.append((a, b))
        if j < len(x) and not x[j]:
            continue  # 开断型号 P=Q=ell=0，原电流锥恒成立。
        if j < len(x) and equations.network.name == 'case33bw':
            # 按逐线路电流界平衡电压与电流，保持原二阶锥等价。
            key = equations.keys[j]
            cone_scale = (np.clip(1./equations.ellmax[key], 1., 100.)
                          if np.isfinite(equations.network.ell_limit[j]) else 100.)
            v = operation.v[equations.ends[key[0]][0]]
            ell, P, Q = operation.ell[key], operation.P[key], operation.Q[key]
            expressions = [affine(item) for item in
                           (v+cone_scale*ell, 2*np.sqrt(cone_scale)*P, 2*np.sqrt(cone_scale)*Q, v-cone_scale*ell)]
            a, b = np.array([item[0] for item in expressions]), np.array([item[1] for item in expressions])
        blocks.append(sparse.csc_matrix(-a[:, free]))
        constant = b+a@solution
        constant[0] -= margin*(cone_scale if j < len(x) and equations.network.name == 'case33bw' else 1.)
        values.append(constant)
        cones.append(clarabel.SecondOrderConeT(len(b)))
    objective, _ = affine(model.getObjective())
    objective *= model.ModelSense
    # 射线的线性等式先精确消元，避免小阻抗压降行在锥求解中损失精度。
    if equal.any():
        left, singular, right = np.linalg.svd(blocks.pop(0).toarray(), full_matrices=True)
        rank = int((singular > singular[0]*max(len(free), int(equal.sum()))*np.finfo(float).eps).sum())
        offset = right[:rank].T@((left[:, :rank].T@values.pop(0))/singular[:rank])
        basis = right[rank:].T
        cones.pop(0)
        matrix_free = sparse.vstack(blocks, format='csc')
        conic_matrix = sparse.csc_matrix(matrix_free@basis)
        conic_rhs = np.concatenate(values)-matrix_free@offset
        conic_objective = basis.T@objective[free]
    else:
        conic_matrix = sparse.vstack(blocks, format='csc')
        conic_rhs, conic_objective = np.concatenate(values), objective[free]
    settings = clarabel.DefaultSettings()
    settings.verbose = False
    settings.max_threads = threads
    settings.static_regularization_constant = 1e-10
    if np.isfinite(equations.network.ell_limit).all():
        settings.static_regularization_constant = 1e-7   # 1e-8 在紧化射线上偶发 NumericalError
        settings.direct_solve_method = 'faer'
        settings.max_step_fraction = .95
    settings.tol_gap_abs = settings.tol_gap_rel = tolerance
    settings.tol_feas = tolerance
    settings.time_limit = max(0., deadline-perf_counter())
    count = len(conic_objective)
    answer = clarabel.DefaultSolver(sparse.csc_matrix((count, count)), conic_objective,
        conic_matrix, conic_rhs, cones, settings).solve()
    if answer.status == clarabel.SolverStatus.MaxTime:
        raise TimeoutError('连续子问题达到总时限')
    if answer.status not in (clarabel.SolverStatus.Solved, clarabel.SolverStatus.AlmostSolved):
        raise RuntimeError(f'Continuous SOCP: {answer.status}, x={list(x)}')
    solution[free] = offset+basis@answer.x if equal.any() else answer.x

    # 3. 以原始约束重新核验；锥求解器的缩放残差不直接作为物理证书。
    residual = matrix@solution-rhs
    violation = max(0., np.max(np.where(senses == '=', np.abs(residual),
                    np.where(senses == '>', -residual, residual))),
                    np.max(lb-solution), np.max(solution-ub))
    for constraint in model.getQConstrs():
        expression = model.getQCRow(constraint)
        a, b = affine(expression.getLinExpr())
        value = a@solution+b+sum(expression.getCoeff(j)*solution[expression.getVar1(j).index]
            *solution[expression.getVar2(j).index] for j in range(expression.size()))
        violation = max(violation, value-constraint.QCRHS)
    normals = []
    for a, b in cone_rows:
        tail = (a@solution+b)[1:]
        length = np.linalg.norm(tail)
        normals.append(tail/length if length else tail)
    return solution, violation, normals


class PortSubProblem(SubProblem):
    def solve(self, x, power, time_limit=None, *, score_only=False):
        if np.isfinite(self.equations.network.ell_limit).all():
            return super().solve(x, power, time_limit, score_only=score_only)
        self.calls += 1
        deadline = perf_counter()+(SP_TIME_LIMIT[self.equations.method] if time_limit is None else time_limit)
        with new_model('port_SP', self.threads) as model:
            choice, p, eta, operation, _ = self._build(model, x, power)
            model.update()
            fixed = {v.index: float(a) for v, a in zip([*choice.values(), *p.values()], [*x, *power])}
            solution, violation, normals = solve_conic(model, self.equations, operation, x,
                fixed, self.threads, deadline, 1e-10, 0.)
            value = float(solution[eta.index])
            if max(0., value)+violation <= PLANNING_TOL:
                state = solution[[v.index for v in operation.state.tolist()]]
                return dict(cut=None, state=state, feasible=True, eta=value)
            if value <= PLANNING_TOL:
                raise RuntimeError(f'SP: eta={value:g}, MaxVio={violation:g}; certificate exceeds tolerance')
        if score_only:
            return dict(cut=None, state=None, feasible=False, eta=value, cone_normals=normals)
        cut = self.generate_cut(x, power, normals, time_limit=deadline-perf_counter())
        return dict(cut=cut, state=None, feasible=False, eta=value)


def ray_support(equations, budget, x, anchor, power, *, threads, time_limit, numeric_focus=0):
    """固定网架的近边界射线点；原约束残差达标后才收入认证域。"""
    deadline = perf_counter()+time_limit
    with new_model('port_ray', threads) as model:
        choice = dict(zip(equations.keys, x))
        ray_fraction = model.addVar(ub=1., name='ray_fraction')
        p = {i: float(a)+ray_fraction*float(q-a)
             for i, a, q in zip(equations.network.load_nodes, anchor, power)}
        operation = equations.add_operation(model, choice, p, scheme=x)
        model.setObjective(ray_fraction, gp.GRB.MAXIMIZE)
        solution, violation, _ = solve_conic(model, equations, operation, x, {},
            threads, deadline, 1e-9, RAY_CONE_MARGIN)
        if violation > PLANNING_TOL:
            raise RuntimeError(f'Ray: MaxVio={violation:g}, anchor={list(anchor)}, p={list(power)}')
        return dict(x=x.copy(), p=anchor+solution[ray_fraction.index]*(power-anchor),
                    state=solution[[v.index for v in operation.state.tolist()]],
                    ray_fraction=float(solution[ray_fraction.index]), feasible=True)
