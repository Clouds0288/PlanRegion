"""小规模审核的独立消元模型；只读取 Network，不导入正式 model.py。

固定网架后消去 P/Q/v，仅以电流平方表示损耗；用于核对紧凑模型和割。
"""
from types import SimpleNamespace  # 将参考矩阵作为一次建模的只读数据集合传递。
import clarabel  # 独立求解参考 LP/SOCP，不复用正式规划求解器。
import numpy as np  # 由网架数据独立推导消元矩阵。
from scipy import sparse  # 连续优化接口使用稀疏矩阵。


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


def _solve(objective,matrix,rhs,cones):  # 参考优化器只封装共同的数值设置。
    settings = clarabel.DefaultSettings()  # 使用连续凸优化算法。
    settings.verbose = False  # 核对实验不输出每次内点迭代。
    settings.max_threads = 1  # 所有方法使用统一单线程计时。
    settings.tol_gap_abs = settings.tol_gap_rel = settings.tol_feas = 1e-10  # 控制参考解与对偶界精度。
    solver = clarabel.DefaultSolver(sparse.csc_matrix((len(objective),len(objective))),  # 参考目标没有二次项。
        np.asarray(objective),sparse.csc_matrix(matrix),np.asarray(rhs),cones,settings)  # 输入标准锥形式 A*z+s=b。
    result = solver.solve()  # 返回原始解、对偶解及目标上下界。
    if result.status not in (clarabel.SolverStatus.Solved,clarabel.SolverStatus.AlmostSolved):  # 未求解不能发布参考数值。
        raise RuntimeError(f'Reference solver status {result.status}')  # 保留明确的失败状态供调试。
    return result  # 上层按本次查询的含义解释目标值。




def dispatch_support(network,method,normal):  # 独立求固定方案连续域的支撑值。
    e = _dispatch_equations(network)  # 不读取正式 PlanningEquations 的矩阵。
    n = len(network.load_nodes)  # 每个独立负荷对应一个自由变量。
    if method == 'linear':  # 无损模型已经消去全部运行变量。
        c,matrix = e.linear_c,e.linear_F  # 仅保留负荷变量的线性约束。
        cones = [clarabel.NonnegativeConeT(len(c))]  # 全部约束都是非负余量。
    else:  # SOCP 同时优化负荷和各支路电流平方。
        c,matrix = e.c,np.c_[e.F,e.G]  # 变量按负荷、电流排列。
        cones = [clarabel.NonnegativeConeT(e.linear_count)]  # 先放线性约束。
        cones += [clarabel.SecondOrderConeT(size) for size in e.sizes]  # 再放原始电流锥。
    m = matrix.shape[1]  # 参考问题总变量数。
    limits = np.zeros((n+1,m))  # 单独加入负荷非负和总量外界。
    limits[:n,:n] = np.eye(n)  # 各负荷非负。
    limits[-1,:n] = -1./network.base  # 总负荷不超过共同外界。
    matrix = -np.vstack([matrix,limits])  # 将余量形式转换为 Clarabel 的 A*z+s=b。
    rhs = np.r_[c,np.zeros(n),network.power_limit/network.base]  # 总量约束使用标幺余量。
    cones += [clarabel.NonnegativeConeT(n+1)]  # 负荷边界均为线性约束。
    objective = np.r_[-np.asarray(normal),np.zeros(m-n)]  # 最小化负支撑值，电流不进入目标。
    result = _solve(objective,matrix,rhs,cones)  # 获得连续凸域的原始最优值及对偶界。
    return dict(p=np.asarray(result.x[:n]),value=-result.obj_val,bound=-result.obj_val_dual)  # 恢复最大化方向。


def dispatch_state(network,power):  # 独立求固定负荷点的 SOCP 电流状态。
    e = _dispatch_equations(network)
    n = e.G.shape[1]
    cones = [clarabel.NonnegativeConeT(e.linear_count)]
    cones += [clarabel.SecondOrderConeT(size) for size in e.sizes]
    return np.asarray(_solve(np.ones(n),-e.G,e.c+e.F@np.asarray(power),cones).x)


def validate_power_flow(network):  # 用另一套节点导纳矩阵算法交叉核验支路 AC 实现。
    """用六个工况，将独立支路递推与节点导纳矩阵 Newton–Raphson 潮流核对。"""
    from vertify import ACPowerFlow
    import warnings  # 仅在转换原始数据时屏蔽已知的接口弃用提示。
    import pandapower as pp  # 采用独立实现的 Newton–Raphson 潮流。
    from pandapower.converter.pypower.from_ppc import from_ppc  # 将标准 MATPOWER 数据转换为 pandapower 网架。

    reference = ACPowerFlow(network, threads=1)  # 被核验的是独立 AC 递推，不是 SOCP 方程。
    maximum_difference = 0.  # 记录所有工况、所有节点的最大电压幅值差。
    for scale in (1., 0., .5, 1.2, 1.5, 2.):  # 覆盖原始、零独立负荷及多个放大工况。
        power = scale*network.original_p[network.selected]  # 只缩放三个独立节点，其他背景负荷固定。
        with warnings.catch_warnings():  # 警告过滤仅在本次格式转换的作用域内有效。
            warnings.filterwarnings('ignore', category=FutureWarning,  # 只过滤接口的未来弃用警告。
                                    module='pandapower.converter.pypower.from_ppc')  # 限定来源模块，不屏蔽潮流求解异常。
            net = from_ppc(network.ppc(power), f_hz=50)  # 同一物理网架转换成节点导纳矩阵模型。
        pp.runpp(net, algorithm='nr', tolerance_mva=1e-10, numba=False)  # 另一实现的完整 AC 潮流作为交叉核验。
        # 此处只比较潮流方程；即使电压越限也继续收敛，不能复用 classify 的提前不可行判定。
        ell = np.zeros((1, network.n))  # 支路递推从零电流平方开始。
        for _ in range(160):  # 逐次满足完整 AC 电流等式。
            P, Q, v, u = reference.state(power, ell)  # 当前电流下的送端功率及两端电压平方。
            new = (P*P+Q*Q)/u  # 强制完整 AC 电流等式。
            if np.max(np.abs(new-ell)) < 1e-13:  # 电流平方增量达到比运行限值判定更严格的精度。
                break  # 电流更新已经收敛，停止迭代。
            ell = new  # 更新电流，进行下一次完整功率平衡与压降计算。
        difference = float(np.max(np.abs(np.sqrt(v[0])-net.res_bus.loc[list(network.nodes), 'vm_pu'])))  # v 是平方，NR 输出是幅值。
        assert difference < 1e-8, 'Branch AC and nodal AC disagree'  # 不让方程实现不一致的参考模型进入正式评价。
        maximum_difference = max(maximum_difference, difference)  # 保留全部工况中的最大电压幅值误差。
        if scale == 1:  # 保存原始基础工况的物理校验值。
            baseline = dict(minimum_voltage_pu=float(np.sqrt(v.min())),  # 全网最低电压幅值。
                            minimum_voltage_bus=int(network.nodes[np.argmin(v)]),  # 对应真实节点号。
                            active_loss_kw=float((ell*network.r).sum()*network.base))  # Σr*ell，从标幺还原为 kW。
    return dict(baseline=baseline, nodal_cross_checks=6,  # 返回原始工况和交叉核验次数。
                maximum_voltage_difference_pu=maximum_difference, pandapower=pp.__version__)
