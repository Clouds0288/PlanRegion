"""生产运行状态的离线审核；不参与 MP/SP 的可行性分支。"""
from time import perf_counter
from model import PlanningModel, PlanningSP
from vertify import ACPowerFlow
import numpy as np


def margin(equations, x, power, state):
    """按 add_operation 的顺序复核原约束；每项余量非负表示满足，等式余量为负绝对误差。"""
    # 1. 将扁平输入还原成按节点、走廊和线路型号索引的物理量
    net = equations.network
    P = dict(zip(equations.keys, state[equations.P_slice]))
    Q = dict(zip(equations.keys, state[equations.Q_slice]))
    ell = dict(zip(equations.keys, state[equations.ell_slice]))
    v = dict(zip(net.nodes, state[equations.v_slice])) | {net.root: 1.}
    x, p = dict(zip(equations.keys, x)), dict(zip(net.load_nodes, power))
    plus_values, minus_values = np.split(state[equations.slack_slice], 2)
    plus, minus = dict(zip(equations.types, plus_values)), dict(zip(equations.types, minus_values))
    fixed_p, fixed_q = dict(zip(net.nodes, net.fixed_p)), dict(zip(net.nodes, net.fixed_q))
    q_ratio = dict(zip(net.load_nodes, net.q_ratio))

    # 2. 检查完整运行向量的全局上下界
    residual = []
    residual.append(np.min(state-equations.y_lb_global))  # 全部运行变量满足全局下界：y >= y_lb_global。
    residual.append(np.min(equations.y_ub_global-state))  # 全部运行变量满足全局上界：y <= y_ub_global。

    # 3. 逐节点检查有功平衡、无功平衡和电压下限
    for i in net.nodes:
        dp = (fixed_p[i]+p.get(i, 0.))/net.base
        dq = (fixed_q[i]+q_ratio.get(i, 0.)*p.get(i, 0.))/net.base
        incoming_p = sum(P[e, k]-equations.r[e, k]*ell[e, k] for e in equations.incoming[i] for k in equations.types[e])
        outgoing_p = sum(P[e, k] for e in equations.outgoing[i] for k in equations.types[e])
        incoming_q = sum(Q[e, k]-equations.reactance[e, k]*ell[e, k] for e in equations.incoming[i] for k in equations.types[e])
        outgoing_q = sum(Q[e, k] for e in equations.outgoing[i] for k in equations.types[e])
        residual.append(-abs(incoming_p-outgoing_p-dp))  # 节点有功平衡：流入有功扣除线损和流出有功后等于有功负荷。
        residual.append(-abs(incoming_q-outgoing_q-dq))  # 节点无功平衡：流入无功扣除线损和流出无功后等于无功负荷。
        residual.append(v[i]-equations.vmin[i])  # 节点电压平方下限：v_i >= vmin_i；上限已由全局变量盒检查。

    # 4. 逐走廊检查压降等式、开断余量和各线路型号的运行约束
    for e, (i, j) in equations.ends.items():
        z = sum(x[e, k] for k in equations.types[e])
        drop = v[j]-v[i]+sum(2*(equations.r[e, k]*P[e, k]+equations.reactance[e, k]*Q[e, k])-(equations.r[e, k]**2+equations.reactance[e, k]**2)*ell[e, k] for k in equations.types[e])
        residual.append(-abs(drop-plus[e]+minus[e]))  # 走廊压降等式：v_j-v_i+2(rP+χQ)-(r²+χ²)ell-s⁺+s⁻ = 0。
        residual.append(equations.drop_max[e]*(1-z)-plus[e])  # 正向开断余量上界：s⁺_e <= M_e(1-z_e)。
        residual.append(equations.drop_max[e]*(1-z)-minus[e])  # 负向开断余量上界：s⁻_e <= M_e(1-z_e)。
        for k in equations.types[e]:
            residual.append(P[e, k]-equations.pmin[e, k]*x[e, k])  # 型号有功下界：P_ek >= Pmin_ek*x_ek。
            residual.append(equations.pmax[e, k]*x[e, k]-P[e, k])  # 型号有功上界：P_ek <= Pmax_ek*x_ek。
            residual.append(Q[e, k]-equations.qmin[e, k]*x[e, k])  # 型号无功下界：Q_ek >= Qmin_ek*x_ek。
            residual.append(equations.qmax[e, k]*x[e, k]-Q[e, k])  # 型号无功上界：Q_ek <= Qmax_ek*x_ek。
            residual.append(equations.ellmax[e, k]*x[e, k]-ell[e, k])  # 电流平方上界：ell_ek <= ellmax_ek*x_ek。
            if i != net.root:
                residual.append(equations.pmax[e, k]*x[e, k]+P[e, k]-equations.r[e, k]*ell[e, k])  # 反向有功容量：-P_ek+r_ek*ell_ek <= Pmax_ek*x_ek。
                residual.append(equations.qmax[e, k]*x[e, k]+Q[e, k]-equations.reactance[e, k]*ell[e, k])  # 反向无功容量：-Q_ek+χ_ek*ell_ek <= Qmax_ek*x_ek。
            if equations.method == 'socp':
                residual.append(v[i]+ell[e, k]-np.linalg.norm([2*P[e, k], 2*Q[e, k], v[i]-ell[e, k]]))  # 支路电流锥：||(2P,2Q,v_i-ell)||₂ <= v_i+ell。

    # 5. 检查根节点电源容量，返回所有约束中的最小余量
    ps = sum(P[e, k] for e in equations.outgoing[net.root] for k in equations.types[e])
    qs = sum(Q[e, k] for e in equations.outgoing[net.root] for k in equations.types[e])
    residual.append(net.source_pmax-ps)  # 电源有功上限：P_source <= source_pmax。
    residual.append(net.source_qmax-qs)  # 电源无功上限：Q_source <= source_qmax。
    if equations.method == 'socp':
        residual.append(net.source_smax-np.hypot(ps, qs))  # 电源视在功率上限：sqrt(P_source²+Q_source²) <= source_smax。
    return float(min(residual))


def joint_benders(equations, *, power=None, budget=np.inf,
                  cuts=(), start=None, radial_gap_kw=1e-3,
                  min_total=None, fixed_plan=None, incumbent=None, deadline=np.inf, threads=1, oracle=None):
    """供独立对照使用的 MP1/MP2 与 SP 联合割迭代。"""
    problem = PlanningModel(equations, power=power, budget=budget,
                            cuts_only=True, min_total=min_total, fixed_plan=fixed_plan, threads=threads)
    oracle, generated = oracle or PlanningSP(equations, threads=threads), []
    minimizing = power is not None or min_total is not None
    bound = -np.inf if minimizing else np.inf
    with problem.model:
        for cut in cuts:
            problem.add_cut(cut)
        if start is not None:
            problem.x.Start = start
        if incumbent is not None:
            incumbent_cost = problem.use_incumbent(incumbent)
        while True:
            if perf_counter() >= deadline:
                return dict(feasible=False, x=None, p=None, state=None, objective=None,
                            bound=bound, status='unknown', termination='case_time_limit'), generated
            answer = problem.solve(time_limit=min(20., max(0., deadline-perf_counter())))
            if answer is None:
                return None, generated
            bound = max(bound, answer['bound']) if minimizing else min(bound, answer['bound'])
            answer['bound'] = bound
            if incumbent is not None and bound >= incumbent_cost-1e-7:
                answer = dict(incumbent, objective=incumbent_cost, bound=bound, status='optimal')
                return answer, generated
            if answer['x'] is None:
                return answer, generated
            point = answer['p']
            if not minimizing and point.sum() > 0.:
                point = point * max(0., 1. - radial_gap_kw / point.sum())
            checked = oracle.solve(answer['x'], point, time_limit=min(
                5. if equations.method == 'linear' else 10., max(0., deadline-perf_counter())))
            if not checked['feasible'] and checked['cut'] is None and not minimizing:
                # 缩入点数值未决时，在 MP 原点分离；只复用全局有效割，不把未知点判为可行。
                separated = oracle.solve(answer['x'], answer['p'], time_limit=min(
                    5. if equations.method == 'linear' else 10., max(0., deadline-perf_counter())))
                if separated['cut'] is not None:
                    checked = separated
            if checked['feasible']:
                value = answer['objective'] if minimizing else point.sum()
                gap = value-bound if minimizing else bound-value
                tolerance = 1e-7 if minimizing else radial_gap_kw+1e-5
                answer.update(p=point, state=checked['state'], feasible=True, objective=value,
                              status='optimal' if gap <= tolerance else 'feasible')
                return answer, generated
            if checked['cut'] is None:
                answer.update(status='unknown', feasible=False, objective=None)
                return answer, generated
            generated.append(checked['cut'])
            problem.add_cut(checked['cut'])


def affordable_designs(network, budget):
    """固定初始拓扑的升级枚举；不能作为完整重构域的参考。"""
    if network.n_corridors != network.n or not all(c.initial_active for c in network.corridors):
        raise ValueError('Upgrade enumeration requires a fixed initial topology')
    if not np.isfinite(budget):
        raise ValueError('Unlimited budgets require global search, not enumeration')
    result = []
    upgrades = [c for c in network.corridors if len(c.types) > 1]
    def visit(index, cost, chosen):
        if index == len(upgrades):
            choice = chosen
            result.append((cost, choice, ACPowerFlow(network.tree(network.encode_plan(choice)), threads=1)))
            return
        corridor = upgrades[index]
        for line_type in corridor.types:
            price = line_type.investment_cost
            if price < 0.:
                raise ValueError('Budget pruning requires nonnegative investment costs')
            if cost+price <= budget:
                visit(index+1, cost+price, chosen | {corridor.id: line_type.id})
    visit(0, 0., network.initial_plan)
    return sorted(result, key=lambda row: (row[0], tuple(row[1].items())))
