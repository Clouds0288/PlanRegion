"""生产运行状态的离线审核；不参与 MP/SP 的可行性分支。"""
import numpy as np
from model import DEFAULT_SOLVER_THREADS, PlanningEquations, PlanningModel


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
