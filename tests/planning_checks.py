"""生产运行状态的离线审核与无求解器的合成回放；不参与 MP/SP 的可行性分支。"""
from itertools import product

import numpy as np


def recorded_monitor(d, output=None):
    """合成回放（不调用求解器）：分区 '-'*d 的根锥求界；网架 '1' 的 O_x 从分区盒出发，一次支撑查询未给出上界，
    下一次的上界 Σu<=150 裁剪 O_x 后网架按 ε_B 停止，P_x 取 O_x 的一半；根锥沿 e1,e2 的棱中点二分。各帧以子进程的
    幅值 kW 送入 forward。返回 (monitor, x, cut)。"""
    from monitor import RunMonitor
    from Network.four_bus_five_corridor import FourBus
    from region import clip_polytope
    network = FourBus(load_nodes=(1, 2, 3)[:d])
    bounds, label = np.full(d, 100.), '-'*d
    monitor = RunMonitor(output=output, algorithm='RB')
    monitor._emit('start', phase='初始化', network=network.name, load_nodes=network.load_nodes, budget=20000.,
                  algorithm='RB', time_limit=30., cost_unit=network.cost_unit, initial_plan=network.initial_plan,
                  bounds=bounds, partitions=['+'*d, label], cones={}, schemes={}, cut_history={}, status='running')
    x = network.encode_plan(network.initial_plan).astype(int)
    cut = np.r_[150., -np.ones(d), np.zeros(len(x))]
    box = np.array(list(product((0., 1.), repeat=d)))*bounds
    clipped = clip_polytope(box, 150., -np.ones(d))
    row = lambda status, outer=clipped: dict(x=x, choice=network.initial_plan, cost=0., outer=outer, inner=.5*outer,
                                             status=status)
    info = lambda status, ratio: dict(certified=d, faces=2*d+1, ratio=ratio, status=status)
    cone = lambda U: dict(inner=np.vstack([np.zeros(d), .5*U*bounds]), outer=np.vstack([np.zeros(d), .6*U*bounds]),
                          scheme='1', mu=1.2)
    middle = np.zeros(d)
    middle[:2] = 1/np.sqrt(2)
    halves = [np.vstack([np.eye(d)[:k], middle, np.eye(d)[k+1:]]) for k in (1, 0)]
    point = dict(scheme='1', p=bounds)
    for patch in (dict(event='phase_start', phase='径向搜索 A', status='running', cones={}, schemes={}, cut_history={}),
                  dict(event='cone', phase='径向搜索 A', cones={'0': cone(np.eye(d))}, cone_count=1, volume_ratio=.44,
                       global_point=dict(scheme='1', p=np.full(d, 30.))),
                  dict(event='point', phase='网架支撑', active_scheme='1', support=info('active', 2**d-1.),
                       sp_point=None, schemes={'1': row('active', box)}),
                  dict(event='cut', phase='网架支撑', active_scheme='1', schemes={'1': row('eps_B')},
                       cut_history={'1': dict(cut=cut, scheme='1')}, support=info('eps_B', .004), sp_point=point),
                  dict(event='cone', phase='径向续跑 A+', cones={'0': None, '1': cone(halves[0]), '2': cone(halves[1])},
                       cone_count=2, volume_ratio=.44, global_point=None)):
        monitor.forward(label, patch)
    return monitor, x, cut


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
