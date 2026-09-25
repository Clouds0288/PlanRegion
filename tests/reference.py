"""小规模审核的独立消元模型；只读取 Network，不导入正式 model.py。

固定网架后消去 P/Q/v，仅以电流平方表示损耗；用于核对紧凑模型和割。
"""
from types import SimpleNamespace  # 将参考矩阵作为一次建模的只读数据集合传递。
import gurobipy as gp
from gurobipy import GRB
import numpy as np  # 由网架数据独立推导消元矩阵。


def fixed_topology(network):
    """历史升级基准专用：仅保留初始树走廊，仍使用统一规划模型。"""
    from dataclasses import fields
    from Network import Network
    data = {field.name: getattr(network, field.name) for field in fields(Network)}
    selected = np.array([c.initial_active for c in network.corridors])
    data.update(name=network.name+'_initial_tree', corridors=tuple(c for c in network.corridors if c.initial_active),
                road_allowed=network.road_allowed[selected])
    fixed = Network(**data)
    for name in ('projects', 'upgrade_count', 'budgets', 'cost_unit'):
        if hasattr(network, name):
            setattr(fixed, name, getattr(network, name))
    return fixed


def upgrade_plan(network, choices):
    """历史项目顺序只在测试输入边界转换，运行状态始终按 Network 型号顺序。"""
    by_endpoints = {frozenset(c.endpoints): c for c in network.corridors}
    return network.initial_plan | {by_endpoints[frozenset(p.branch)].id: 'parallel' if k else 'existing'
                                   for p, k in zip(network.projects, choices, strict=True)}


def _dispatch_equations(network):  # 构造固定方案的 P/Q/v 消元式。
    c = network  # c 仅在此函数表示网架对象。
    n, D, E = c.n, c.D, c.E  # 下游汇总矩阵与独立负荷节点映射。
    R, X = np.diag(c.r), np.diag(c.reactance)  # 各支路的固定标幺阻抗。
    P0, Q0 = D@c.fixed_p/c.base, D@c.fixed_q/c.base  # 固定背景造成的无损功率。
    Pp, Qp = D@E/c.base, D@E*c.q_ratio/c.base  # 三个独立负荷对支路功率的贡献。
    Pl, Ql = D@R, D@X  # 各支路电流平方对上游损耗的贡献。
    v0 = 1.-2*D.T@(R@P0+X@Q0)  # 参数原点的无损电压平方。
    vp = -2*D.T@(R@Pp+X@Qp)  # 独立负荷造成的路径压降。
    vl = -2*D.T@(R@Pl+X@Ql)+D.T@(R@R+X@X)  # 完整损耗压降修正。
    J = np.zeros((n,n))  # 受端电压到送端电压的父节点映射。
    for i,parent in enumerate(c.parent):  # 根相邻支路的送端电压另用常数 1。
        if parent >= 0:  # 其余支路引用父节点的电压平方。
            J[i,parent] = 1.  # 每行只包含一个父节点。
    u0, up, ul = 1.+J@(v0-1.), J@vp, J@vl  # 送端电压的常数、负荷、电流系数。
    roots = c.roots  # 接电源支路的功率求和得到源端注入。
    constant = np.r_[c.capacity-P0,v0-c.vmin,c.vmax-v0,  # 支路有功及节点电压余量。
                     c.source_pmax-P0[roots].sum(),c.source_qmax-Q0[roots].sum()]  # 源端 P/Q 余量。
    F = np.vstack([-Pp,vp,-vp,-Pp[roots].sum(axis=0),-Qp[roots].sum(axis=0)])  # 运行限值的负荷系数。
    G = np.vstack([-Pl,vl,-vl,-Pl[roots].sum(axis=0),-Ql[roots].sum(axis=0)])  # 同一顺序下的电流系数。
    finite = np.isfinite(constant)  # 数据未给定的无穷上界不建立约束。
    linear_c, linear_F = constant[finite], F[finite]  # ell=0 时的独立线性参考约束。
    constants = [np.r_[np.zeros(n),linear_c]]  # 在运行限值前保留 ell≥0。
    powers = [np.vstack([np.zeros((n,E.shape[1])),linear_F])]  # ell≥0 与负荷无关。
    currents = [np.vstack([np.eye(n),G[finite]])]  # ell≥0 直接使用单位阵。
    linear_count = len(constants[0])  # 记录全部非负锥约束的行数。
    relax = [np.r_[np.zeros(n),np.ones(len(linear_c))]]  # eta 不松弛电流非负约束。
    if np.isfinite(c.source_smax):  # 仅在网架提供视在功率上限时加入源端锥。
        constants.append(np.array([c.source_smax,P0[roots].sum(),Q0[roots].sum()]))  # 源端锥常数项。
        powers.append(np.vstack([np.zeros(E.shape[1]),Pp[roots].sum(axis=0),Qp[roots].sum(axis=0)]))  # 源端锥负荷系数。
        currents.append(np.vstack([np.zeros(n),Pl[roots].sum(axis=0),Ql[roots].sum(axis=0)]))  # 源端锥损耗系数。
        relax.append(np.array([1.,0.,0.]))  # 只松弛标准锥首分量。
    for i,unit in enumerate(np.eye(n)):  # 每条支路加入一个四维电流锥。
        constants.append(np.array([u0[i],2*P0[i],2*Q0[i],u0[i]]))  # (u+ell,2P,2Q,u-ell) 的常数项。
        powers.append(np.vstack([up[i],2*Pp[i],2*Qp[i],up[i]]))  # 电流锥的负荷系数。
        currents.append(np.vstack([ul[i]+unit,2*Pl[i],2*Ql[i],ul[i]-unit]))  # 首尾分量分别加减 ell_i。
        relax.append(np.array([1.,0.,0.,0.]))  # eta 只扩张锥轴。
    sizes = [len(block) for block in constants[1:]]  # 源端锥及支路锥的维数。
    return SimpleNamespace(c=np.concatenate(constants),F=np.vstack(powers),G=np.vstack(currents),  # 拼接固定方案的完整消元矩阵。
                           linear_c=linear_c,linear_F=linear_F,linear_count=linear_count,  # 保留 LP 子集供解析判定。
                           sizes=sizes,relax=np.concatenate(relax))  # 返回锥维数及 phase I 松弛方向。


def _solve(objective, matrix, rhs, cones):
    """独立消元参考：cones 元素为 ('linear' 或 'soc', 行数)。"""
    with gp.Model('independent_dispatch') as model:
        model.Params.OutputFlag, model.Params.Threads = 0, 1
        model.Params.FeasibilityTol = model.Params.OptimalityTol = 1e-9
        model.Params.BarQCPConvTol = 1e-9
        model.Params.NonConvex, model.Params.DualReductions = 0, 0
        variables = model.addMVar(len(objective), lb=-GRB.INFINITY, name='dispatch')
        slack = model.addMVar(len(rhs), lb=-GRB.INFINITY, name='slack')
        model.addConstr(matrix@variables+slack == rhs)
        offset = 0
        for kind, size in cones:
            block = slack[offset:offset+size]
            if kind == 'linear':
                model.addConstr(block >= 0.)
            else:
                # 标准二阶锥：首分量非负，且不小于其余分量的范数。
                head = block[0].item()
                head.LB = 0.
                model.addQConstr(gp.quicksum(item.item()**2 for item in block[1:]) <= head**2)
            offset += size
        model.setObjective(np.asarray(objective)@variables)
        model.optimize()
        if model.Status == GRB.INFEASIBLE:
            return None
        assert model.Status == GRB.OPTIMAL, model.Status
        return SimpleNamespace(x=variables.X, obj_val=model.ObjVal, obj_val_dual=model.ObjVal)


def dispatch_support(network, method, normal):
    """独立消元模型的支撑值，物理方程不读取正式 MP/SP。"""
    equations = _dispatch_equations(network)
    n = len(network.load_nodes)
    if method == 'linear':
        constant, matrix = equations.linear_c, equations.linear_F
        cones = [('linear', len(constant))]
    else:
        constant, matrix = equations.c, np.c_[equations.F, equations.G]
        cones = [('linear', equations.linear_count)]+[('soc', size) for size in equations.sizes]
    dimension = matrix.shape[1]
    limits = np.zeros((n+1, dimension))
    limits[:n, :n] = np.eye(n)
    limits[-1, :n] = -1./network.base
    matrix = -np.vstack([matrix, limits])
    rhs = np.r_[constant, np.zeros(n), network.power_limit/network.base]
    cones.append(('linear', n+1))
    result = _solve(np.r_[-np.asarray(normal), np.zeros(dimension-n)], matrix, rhs, cones)
    return dict(p=result.x[:n], value=-result.obj_val, bound=-result.obj_val_dual)


def dispatch_state(network, power):
    equations = _dispatch_equations(network)
    cones = [('linear', equations.linear_count)]+[('soc', size) for size in equations.sizes]
    return _solve(np.ones(equations.G.shape[1]), -equations.G,
                  equations.c+equations.F@np.asarray(power), cones).x


def nodal_voltages(network, power, start):
    """Gurobi 直角坐标 AC 节点方程，与支路递推独立交叉验证。"""
    net = network
    p, q = net.loads(power)
    admittance = np.zeros((net.n+1, net.n+1), dtype=complex)
    for i, parent in enumerate(net.parent):
        j = net.n if parent < 0 else parent
        value = 1./complex(net.r[i], net.reactance[i])
        admittance[i, i] += value
        admittance[j, j] += value
        admittance[i, j] -= value
        admittance[j, i] -= value
    with gp.Model('nodal_AC') as model:
        model.Params.OutputFlag, model.Params.Threads = 0, 1
        model.Params.NonConvex = 2
        model.Params.FeasibilityTol = model.Params.OptimalityTol = 1e-9
        real = model.addVars(net.n, lb=.3, ub=1.1, name='voltage_real')
        imag = model.addVars(net.n, lb=-.6, ub=.6, name='voltage_imag')
        real[net.n], imag[net.n] = 1., 0.
        for i in range(net.n):
            real[i].Start, imag[i].Start = start[i].real, start[i].imag
            neighbors = np.flatnonzero(admittance[i])
            current_real = gp.quicksum(admittance[i, j].real*real[j]-admittance[i, j].imag*imag[j] for j in neighbors)
            current_imag = gp.quicksum(admittance[i, j].imag*real[j]+admittance[i, j].real*imag[j] for j in neighbors)
            # S_i=V_i*conj(I_i)；负荷是负注入，根节点承担全网功率平衡。
            model.addQConstr(real[i]*current_real+imag[i]*current_imag == -p[0, i])
            model.addQConstr(imag[i]*current_real-real[i]*current_imag == -q[0, i])
        model.setObjective(0.)
        model.optimize()
        assert model.SolCount, model.Status
        return np.array([complex(real[i].X, imag[i].X) for i in range(net.n)])


def validate_power_flow(network):
    """六个工况：支路递推与 Gurobi 节点 AC 方程对照。"""
    from vertify import ACPowerFlow
    reference = ACPowerFlow(network, threads=1)
    maximum_difference = 0.
    for scale in (1., 0., .5, 1.2, 1.5, 2.):
        power = scale*network.network.original_p[network.network.selected]
        ell = np.zeros((1, network.n))
        for _ in range(160):
            P, Q, v, u = reference.state(power, ell)
            new = (P*P+Q*Q)/u
            if np.max(np.abs(new-ell)) < 1e-13:
                break
            ell = new
        voltage = np.ones(network.n, dtype=complex)
        for i in network.order:
            parent = network.parent[i]
            upstream = 1. if parent < 0 else voltage[parent]
            current = np.conj(complex(P[0, i], Q[0, i])/upstream)
            voltage[i] = upstream-complex(network.r[i], network.reactance[i])*current
        nodal = nodal_voltages(network, power, voltage)
        difference = float(np.max(np.abs(np.sqrt(v[0])-np.abs(nodal))))
        assert difference < 1e-8, 'Branch AC and nodal AC disagree'
        maximum_difference = max(maximum_difference, difference)
        if scale == 1.:
            baseline = dict(minimum_voltage_pu=float(np.sqrt(v.min())),
                            minimum_voltage_bus=int(network.nodes[np.argmin(v)]),
                            active_loss_kw=float((ell*network.r).sum()*network.base))
    return dict(baseline=baseline, nodal_cross_checks=6,
                maximum_voltage_difference_pu=maximum_difference, gurobi=gp.gurobi.version())
