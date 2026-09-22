"""最终独立 AC 校核、AC 采样及 FR/MR；未确定点始终保留为 0。"""
import gurobipy as gp
from gurobipy import GRB
import numpy as np

from model import PlanningEquations, PlanningModel
from region import classify_orthant, split_grid_box

AC_TOL = 1e-9
FIXED_POINT_TOL = 1e-12
METHODS = ('socp', 'hybrid', 'ac', 'linear')
METHOD_NAMES = ('纯 SOCP 切割', '线性 + SOCP 精修', 'AC 数值参考', '纯线性切割')


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

    def global_status(self, power, environment, time_limit=10.):  # 独立非凸核验在给定时限内返回证书或未确定。
        if self.model is None:  # 不动点已能判定的方案无需创建全局求解器。
            self._build_global(environment)  # 首次遇到未确定点时才付出建模成本。
        m, ell, bp, bq = self.model  # 复用同一方案的电流变量和节点平衡行。
        m.Params.TimeLimit = max(0.,time_limit)  # 设置本次独立 AC 求解的时限。
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

def ac_planning_query(equations, power, *, observer=None):
    """SOCP 搜索建设方案；完整 AC 等式负责认证。"""
    notify = observer or (lambda event, **data: None)
    notify( 'query_model', message='建立 AC 方案搜索的 SOCP 外松弛')
    problem = PlanningModel(equations, power=power)
    bound = -np.inf
    with problem.model:
        iteration = 0
        while True:
            iteration += 1
            notify( 'mp_start', iteration=iteration, message='搜索候选建设方案')
            answer = problem.solve()
            if answer is None:
                notify( 'query_end', status='infeasible', message='全部候选已证不可行')
                return None
            bound = max(bound, answer['bound'])
            answer['bound'] = bound
            if answer['x'] is None:
                notify( 'query_end', status='unknown', message='无整数候选，保持未确定')
                return answer
            choice = equations.choice(answer['x'])
            notify( 'ac_start', point=power, choice=choice, bound=bound,
                 objective=answer['objective'], message='独立 AC 潮流认证')
            oracle = ACPowerFlow(equations.network.design(choice))
            try:
                status = int(oracle.classify(power)[0])
                if status == 0:
                    notify( 'ac_global', message='不动点未获证，进行非凸 AC 全局验证')
                    status = oracle.global_status(power, None)
            finally:
                oracle.close()
            if status == 1:
                answer.update(feasible=True, status='optimal' if answer['objective']-bound <= 1e-7 else 'feasible')
                notify( 'query_end', status=answer['status'], message='获得独立 AC 可行证书')
                return answer
            if status == 0:
                answer.update(status='unknown', termination='ac_unknown', feasible=False, objective=None)
                notify( 'query_end', status='unknown', message='AC 未获证，保持未确定')
                return answer
            problem.exclude(answer['x'])
            notify( 'ac_exclude', message='当前方案已证 AC 不可行，仅在本负荷点排除')


def validate_ac_region(network, budgets, divisions, bounds, *, observer=None):
    """连续构域完成后的独立 AC 网格校核；网格不参与规划域构造。"""
    assert (np.all(network.r > 0.) and np.all(network.reactance >= 0.)
            and np.all(network.original_p >= 0.) and np.all(network.original_q >= 0.)
            and np.all(network.q_ratio >= 0.) and np.all(network.vmax >= 1.))
    notify = observer or (lambda event, **data: None)
    bounds, budgets = np.asarray(bounds), np.asarray(budgets)
    states = np.zeros((len(budgets),)+(divisions,)*3, dtype=np.int8)
    notify('phase_start', method='ac', phase='ac', phase_number=1, phase_count=1,
           states=states, bounds=bounds, budgets=budgets, load_nodes=network.load_nodes,
           message='开始 ac 阶段')
    equations = PlanningEquations(network, 'socp')
    visited = set()
    pending = [(j, np.zeros(3, dtype=int), np.full(3, divisions-1, dtype=int))
               for j in reversed(range(len(budgets)))]
    while pending:
        j, lower, upper = pending.pop()
        block = (j,)+tuple(slice(int(a), int(b)+1) for a, b in zip(lower, upper))
        if np.all(states[block] != 0):
            continue
        notify('block', lower=lower, upper=upper, budget_index=j,
               pending=len(pending), message='选择待分类区域块')
        for index in (upper, lower):
            key = tuple(map(int, index))
            if states[(j,)+key] != 0 or key in visited:
                continue
            power = (index+.5)*bounds/divisions
            notify('point', point=power, index=index, query=len(visited)+1,
                   message='选择网格中心并开始验证')
            answer = ac_planning_query(equations, power, observer=observer)
            visited.add(key)
            labels = np.zeros(len(budgets), dtype=np.int8)
            if answer is None:
                labels[:] = -1
            else:
                if answer['bound'] is not None:
                    labels[budgets < answer['bound']-1e-7] = -1
                if answer['feasible']:
                    labels[budgets >= answer['objective']] = 1
            for k, status in enumerate(labels):
                classify_orthant(states[k], index, status)
            notify('classified', states=states, labels=labels, pool_size=0,
                   message='传播单调证书并更新区域')
            if np.all(states[block] != 0):
                break
        if np.any(states[block] == 0) and np.any(upper > lower):
            pending.extend((j, a, b) for a, b in split_grid_box(lower, upper))
            notify('split', pending=len(pending), message='未知区域沿最长轴二分')
    notify('phase_end', states=states, pending=0, message='ac 阶段结束')
    return states


def comparison_labels(states):
    """0 域外、1 多余、2 遗漏、3 重合、4 未确定。"""
    inside = states == 1
    ac = METHODS.index('ac')
    labels = inside.astype(np.uint8)+2*inside[ac].astype(np.uint8)
    labels[(states == 0) | (states[ac] == 0)] = 4
    return labels


def validation_summary(states, metadata):
    """比较实际保留内域与 AC 域；构域未知薄层计入未保留部分，AC 未知给区间。"""
    return [dict(method=method, budget=budget,
                 **disagreement_interval(np.where(states[k, j] == 1, 1, -1),
                                         states[METHODS.index('ac'), j]),
                 method_unknown_cells=int(np.count_nonzero(states[k, j] == 0)),
                 metric_scope='retained_inner_union',
                 total_seconds=metadata['seconds'][method])
            for k, method in enumerate(METHODS) for j, budget in enumerate(metadata['budgets'])]


def validate_power_flow(network):  # 用另一套节点导纳矩阵算法交叉核验支路 AC 实现。
    """用六个工况，将独立支路递推与节点导纳矩阵 Newton–Raphson 潮流核对。"""
    import warnings  # 仅在转换原始数据时屏蔽已知的接口弃用提示。
    import pandapower as pp  # 采用独立实现的 Newton–Raphson 潮流。
    from pandapower.converter.pypower.from_ppc import from_ppc  # 将标准 MATPOWER 数据转换为 pandapower 网架。

    reference = ACPowerFlow(network)  # 被核验的是独立 AC 递推，不是 SOCP 方程。
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


def disagreement(approximation, ac):  # 根据同一网格中的多余和遗漏体素计算 FR/MR。
    """同一等体积网格上的 FR/MR；空域的条件比例无定义。"""
    extra = int(np.count_nonzero(approximation & ~ac))  # 黄色：计算域有、AC 参考域没有。
    missed = int(np.count_nonzero(ac & ~approximation))  # 红色：AC 参考域有、计算域遗漏。
    computed, reference = int(approximation.sum()), int(ac.sum())  # 两个域各自的网格单元数量。
    union = computed+missed  # 计算域∪AC 域；遗漏部分与计算域不相交。
    return dict(fr_percent=100*extra/computed if computed else None,  # FR 分母是计算域。
                mr_percent=100*missed/reference if reference else None,  # MR 分母是 AC 参考域。
                region_error_percent=100*(extra+missed)/union if union else 0.,  # 原区域不一致率：对称差/并集。
                computed_cells=computed, ac_cells=reference, extra_cells=extra, missed_cells=missed)


def disagreement_interval(approximation, ac):  # 三态网格中，未确定点只能给误差区间。
    """输入 -1/0/1 分别表示已证域外、未确定、已证域内；区间不含网格离散误差。"""
    def ratio_bounds(left, right):  # 求 |left\right|/|left| 的保守上下界。
        inside, unknown = left==1, left==0  # 已知域内点必须计入分母，未确定点可取域内或域外。
        certain = np.count_nonzero(inside & (right==-1))  # 已确认的多余点不能被未知标签消除。
        lower_denominator = inside.sum()+np.count_nonzero(unknown & (right!=-1))  # 只增加可能重合的点，使比例最小。
        possible = np.count_nonzero((left!=-1) & (right!=1))  # 所有可能成为多余的点。
        upper_denominator = inside.sum()+np.count_nonzero(unknown & (right!=1))  # 为最大比例只增加可能多余的点。
        return [100*certain/lower_denominator if lower_denominator else 0.,  # 可能空域时不给虚假的严格正下界。
                100*possible/upper_denominator if upper_denominator else 100.]  # 未知空分母使用保守上界。
    unknown = int(np.count_nonzero((approximation==0)|(ac==0)))  # 比较双方中任意一方未确定的单元数。
    exact = disagreement(approximation==1,ac==1) if unknown==0 else dict(fr_percent=None,mr_percent=None,region_error_percent=None)
    fr = ratio_bounds(approximation,ac) if unknown else ([exact['fr_percent']]*2 if exact['fr_percent'] is not None else None)  # 完整空域的比例无定义，不能显示成 0–100%。
    mr = ratio_bounds(ac,approximation) if unknown else ([exact['mr_percent']]*2 if exact['mr_percent'] is not None else None)  # MR 交换两个集合使用同一定义。
    certain_difference = np.count_nonzero(approximation*ac == -1)
    possible_overlap = np.count_nonzero((approximation != -1) & (ac != -1))
    certain_overlap = np.count_nonzero((approximation == 1) & (ac == 1))
    possible_difference = np.count_nonzero(~(((approximation == 1) & (ac == 1)) |
                                            ((approximation == -1) & (ac == -1))))
    low_den, high_den = certain_difference+possible_overlap, certain_overlap+possible_difference
    error = [100*certain_difference/low_den if low_den else 0.,
             100*possible_difference/high_den if high_den else 0.]
    return dict(**exact,fr_interval=fr,mr_interval=mr,region_error_interval=error,unknown_cells=unknown)
