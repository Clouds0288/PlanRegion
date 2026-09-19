"""MP1、MP2、SP：建设主问题与线性电气可行性模型。"""
import gurobipy as gp
import numpy as np
from gurobipy import GRB

TOL = 1e-8


def add_cut(master, x, p, cut):
    """全局可行性割 β0 + βx·x + βp·p ≥ 0。"""
    nx = len(x)
    expression = (cut[0] + gp.LinExpr(cut[1:1+nx], list(x.values()))
                  + gp.LinExpr(cut[1+nx:], list(p.values())))
    master.addConstr(expression / max(1., np.abs(cut).max()) >= 0)


def _master(network, cuts):
    """树结构、非负节点负荷、配变上限、累计割；负荷可自由分配。"""
    master = gp.Model()
    master.Params.OutputFlag = 0
    master.Params.Threads = 1
    master.Params.MIPGap = master.Params.MIPGapAbs = 0
    master.Params.FeasibilityTol = master.Params.IntFeasTol = 1e-9
    master.Params.DualReductions = 0
    ne, nk = len(network.corridors), len(network.lines)
    x = master.addVars(ne, nk, vtype=GRB.BINARY, name="x")
    p = master.addVars(len(network.load_nodes), lb=0, name="p")
    master.addConstrs((x.sum(e, "*") <= 1 for e in range(ne)), name="one_type")
    master.addConstr(x.sum() == len(network.nodes)-1, name="tree_edges")
    # 根节点向每个非根节点发送一单位虚拟流；连通且恰有 n-1 条边即为树。
    flow = master.addVars(ne, lb=-GRB.INFINITY, name="connect_flow")
    for e in range(ne):
        capacity = (len(network.nodes)-1)*x.sum(e, "*")
        master.addConstr(flow[e] <= capacity)
        master.addConstr(-flow[e] <= capacity)
    for node in network.load_nodes:
        master.addConstr(gp.quicksum(flow[e] for e, edge in enumerate(network.corridors) if edge.end == node)
                         - gp.quicksum(flow[e] for e, edge in enumerate(network.corridors) if edge.start == node)
                         == 1, name=f"connect_{node}")
    master.addConstr(p.sum() <= network.power_limit, name="transformer")
    for cut in cuts:
        add_cut(master, x, p, cut)
    return master, x, p


def MP1(network, power, cuts=(), *, total=False, fixed_x=None):
    """最低建设费：固定节点负荷；total=True 时只限定总量，允许自由分配。"""
    master, x, p = _master(network, cuts)
    master.ModelName = "MP1"
    if total:
        master.addConstr(p.sum() >= float(power), name="target_total")
    else:
        master.addConstrs((p[i] == float(value) for i, value in zip(p, power)), name="fixed_power")
    if fixed_x is not None:
        for variable, value in zip(x.values(), fixed_x):
            variable.LB = variable.UB = float(value)
    master.setObjective(gp.LinExpr(network.cost/1000, list(x.values())), GRB.MINIMIZE)
    return master, x, p


def MP2(network, budget=np.inf, cuts=()):
    """给定建设预算，最大化节点总负荷 Σp；np.inf 取消预算上限。"""
    master, x, p = _master(network, cuts)
    master.ModelName = "MP2"
    if np.isfinite(budget):
        master.addConstr(gp.LinExpr(network.cost/1000, list(x.values())) <= budget/1000,
                         name="budget")
    master.setObjective(p.sum(), GRB.MAXIMIZE)
    return master, x, p


class SP:
    """min η：Ww−rη≤h+Tx+Dp；固定功率因数、无损平衡、线性平方电压降。"""

    def __init__(self, network):
        ne, nk = len(network.corridors), len(network.lines)
        pref, span = network.power_limit, 1-network.voltage_min_pu**2
        qr = np.tan(np.arccos(network.power_factor))
        lp = gp.Model("SP")
        lp.Params.OutputFlag = 0
        lp.Params.Threads = 1
        lp.Params.FeasibilityTol = lp.Params.OptimalityTol = 1e-9
        lp.Params.Method = 1
        flow = lp.addVars(ne, lb=-GRB.INFINITY, name="P")
        voltage = lp.addVars(network.load_nodes, lb=-GRB.INFINITY, name="v")
        v = {network.root: 1., **voltage}
        eta = lp.addVar(lb=0, name="eta")
        x = lp.addVars(ne, nk, lb=-GRB.INFINITY, name="x")
        p = lp.addVars(len(network.load_nodes), lb=-GRB.INFINITY, name="p")
        lp.setObjective(eta, GRB.MINIMIZE)
        for e, edge in enumerate(network.corridors):
            capacity = gp.quicksum(line.capacity_kw*x[e, k] for k, line in enumerate(network.lines))
            lp.addConstr(flow[e] <= capacity+pref*eta, name=f"capacity_{edge.name}_upper")
            lp.addConstr(-flow[e] <= capacity+pref*eta, name=f"capacity_{edge.name}_lower")
        for row, i in enumerate(network.load_nodes):
            demand = p[row]
            balance = (gp.quicksum(flow[e] for e, edge in enumerate(network.corridors) if edge.end == i)
                       - gp.quicksum(flow[e] for e, edge in enumerate(network.corridors) if edge.start == i))
            lp.addConstr(balance <= demand+pref*eta, name=f"balance_{i}_upper")
            lp.addConstr(-balance <= -demand+pref*eta, name=f"balance_{i}_lower")
            lp.addConstr(v[i] <= 1+span*eta, name=f"voltage_{i}_upper")
            lp.addConstr(-v[i] <= -network.voltage_min_pu**2+span*eta, name=f"voltage_{i}_lower")
        lp.addConstr(p.sum() <= pref+pref*eta, name="transformer")
        for e, edge in enumerate(network.corridors):
            for k, line in enumerate(network.lines):
                a = 2*(line.r_ohm_km+qr*line.x_ohm_km)*(edge.length_m/1000)/(1000*network.voltage_kv**2)
                big_m = span+a*max(t.capacity_kw for t in network.lines)
                drop = v[edge.end]-v[edge.start]+a*flow[e]
                lp.addConstr(drop <= big_m*(1-x[e, k])+span*eta, name=f"drop_{edge.name}_{line.name}_upper")
                lp.addConstr(-drop <= big_m*(1-x[e, k])+span*eta, name=f"drop_{edge.name}_{line.name}_lower")
        lp.update()
        self.rows = lp.getConstrs()
        states = list(flow.values())+list(voltage.values())
        choices, demands = list(x.values()), list(p.values())
        self.W = np.array([[lp.getCoeff(row, w) for w in states] for row in self.rows])
        self.h = np.array(lp.getAttr("RHS", self.rows))
        self.T = -np.array([[lp.getCoeff(row, z) for z in choices] for row in self.rows])
        self.D = -np.array([[lp.getCoeff(row, z) for z in demands] for row in self.rows])
        self.r = -np.array([lp.getCoeff(row, eta) for row in self.rows])
        lp.remove(choices+demands)
        lp.update()
        self.lp = lp

    def solve(self, x, p):
        self.lp.setAttr("RHS", self.rows, (self.h+self.T@x+self.D@p).tolist())
        self.lp.optimize()
        eta = float(self.lp.ObjVal)
        pi = -np.array(self.lp.getAttr("Pi", self.rows))
        cut = np.r_[pi@self.h, pi@self.T, pi@self.D] if eta > TOL else None
        return eta, cut
