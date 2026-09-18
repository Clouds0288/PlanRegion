"""零基准节点功率规划：MP1 最低建设费、MP2 最大总负荷、SP 电气可行性。

p=(p1,p2,p3)≥0，单位 kW；建设 x 保留在整数主问题中。
电气模型为固定功率因数、无损有功平衡、线性平方电压降。
统一形式 Ww≤h+Tx+Dp；SP 的对偶割为 πᵀ(h+Tx+Dp)≥0。
编号：S 为建设结构，E 为电气约束，Q 为查询，CUT 为可行性割。
独立树上 LP 只作验证与绘图基准，不参与 Benders 选解或生成割。
"""
from __future__ import annotations

import itertools
import math
from dataclasses import dataclass, field

import gurobipy as gp
import numpy as np
import pandas as pd
from gurobipy import GRB


@dataclass(frozen=True)
class LineType:
    name: str                 # 线型 k 的名称。
    r_ohm_km: float            # 单位长度电阻 R_k，Ω/km。
    x_ohm_km: float            # 单位长度电抗 X_k，Ω/km。
    capacity_kw: float        # 有功容量 C_k，kW。
    cable_cny_m: float        # 电缆单价，元/m。


@dataclass(frozen=True)
class Corridor:
    name: str                 # 走廊标签；矩阵列按定义顺序排列。
    start: int                # 参考方向起点 s(e)。
    end: int                  # 参考方向终点 t(e)。
    length_m: float           # 走廊长度 ℓ_e，m。
    existing: bool           # 已有走廊保留 L 型时不新增建设费。


@dataclass(frozen=True)
class DemoConfig:
    power_factor: float = 0.95       # 固定 cosφ；Q_e=tanφ·P_e。
    voltage_kv: float = 0.4          # 额定线电压 U_N，kV。
    voltage_min_pu: float = 0.93     # 电压幅值下限；v_i 为其平方。
    transformer_kva: float = 150.0   # 配变视在功率 S_tr，kVA。
    new_corridor_cny_m: float = 45.0 # 新走廊每米附加费用。
    feasibility_tolerance: float = 1e-8  # SP 最优违反量的接受阈值。

    @property
    def power_limit(self) -> float:
        return self.transformer_kva * self.power_factor  # P_tr=142.5 kW。


DEFAULT_LINES = (
    LineType("L", 1.15, 0.08, 35.0, 28.0),
    LineType("M", 0.62, 0.08, 65.0, 43.0),
    LineType("H", 0.32, 0.08, 100.0, 65.0),
)
DEFAULT_CORRIDORS = (
    Corridor("01", 0, 1, 220.0, True),
    Corridor("12", 1, 2, 180.0, True),
    Corridor("13", 1, 3, 240.0, True),
    Corridor("02", 0, 2, 300.0, False),
    Corridor("23", 2, 3, 160.0, False),
)


@dataclass
class ElectricalModel:
    config: DemoConfig
    lines: tuple[LineType, ...]
    corridors: tuple[Corridor, ...]
    cost: np.ndarray           # c，15 维建设费用，元。
    W: np.ndarray              # 53×8；w=(五条潮流,三个非根平方电压)。
    h: np.ndarray              # 53 维常数；已包含 v0=1 的移项。
    T: np.ndarray              # 53×15；建设变量的右端系数。
    D: np.ndarray              # 53×3；节点功率的右端系数。
    scales: np.ndarray         # r>0；独立于零基准负荷的固定松弛尺度。
    subproblem: gp.Model
    electrical_constraints: list[gp.Constr]

    @property
    def nx(self) -> int:
        return len(self.cost)


@dataclass
class DualCut:
    pi: np.ndarray             # π=-Pi≥0，Wᵀπ=0，rᵀπ≤1。
    constant: float            # β0=πᵀh。
    x_coeff: np.ndarray        # βx=Tᵀπ。
    load_coeff: np.ndarray     # βp=Dᵀπ，三个节点各自的功率系数。
    source_violation: float    # 生成点满足 CUT 左端=-η*。


@dataclass
class SPResult:
    violation: float           # η*；只有 η*=0 时 w 是原始电气可行状态。
    cut: DualCut
    w: np.ndarray              # 潮流和平方电压，按 W 的列顺序。
    rhs: np.ndarray            # 原始右端 h+Tx+Dp，不含 rη。
    residual: np.ndarray       # Ww-rhs；正值表示原始电气行被违反。


@dataclass
class SolveResult:
    mode: str
    status: str
    value: float               # MP1 为元，MP2 为 kW。
    x: np.ndarray
    p: np.ndarray              # 三节点同时实现的功率分配，kW。
    cost_cny: float
    history: list[dict] = field(default_factory=list)


@dataclass
class Design:
    x: np.ndarray              # 结构合法的固定建设方案。
    cost_cny: float
    description: str


def build_model(config: DemoConfig = DemoConfig(), lines=DEFAULT_LINES,
                corridors=DEFAULT_CORRIDORS) -> ElectricalModel:
    """写出全部 53 条电气行，提取 W,h,T,D,r；SP 仅保留 w、η。"""
    ne, nk = len(corridors), len(lines)
    pref = config.power_limit                 # 正的功率尺度 P_ref=142.5，不取 Σp_i⁰。
    span = 1 - config.voltage_min_pu**2        # 电压尺度 V_ref=0.1351。
    qr = math.tan(math.acos(config.power_factor))
    costs = np.array([
        0.0 if edge.existing and k == 0 else
        edge.length_m * (line.cable_cny_m + (0 if edge.existing else config.new_corridor_cny_m))
        for edge in corridors for k, line in enumerate(lines)
    ])
    lp = gp.Model("SP_node_power")
    lp.Params.OutputFlag = 0
    lp.Params.FeasibilityTol = 1e-9
    lp.Params.OptimalityTol = 1e-9
    lp.Params.Method = 1                     # 连续 LP 用对偶单纯形，复用 RHS 查询。
    lp.Params.Threads = 1
    flow = lp.addVars(ne, lb=-GRB.INFINITY, name="P")
    voltage = lp.addVars(range(1, 4), lb=-GRB.INFINITY, name="v")
    v = {0: 1.0, **voltage}                  # 根电压固定；其他电压界写成显式电气行。
    eta = lp.addVar(lb=0.0, name="eta")
    x = lp.addVars(ne, nk, lb=-GRB.INFINITY, name="x")  # 参数占位，提取后删除。
    p = lp.addVars(range(1, 4), lb=-GRB.INFINITY, name="p")
    lp.setObjective(eta, GRB.MINIMIZE)         # SP：min η，固定 x、p 后判断电气可行性。
    for e, edge in enumerate(corridors):
        capacity = gp.quicksum(line.capacity_kw * x[e, k] for k, line in enumerate(lines))
        lp.addConstr(flow[e] <= capacity + pref*eta, name=f"capacity_{edge.name}_1")  # E01—E10。
        lp.addConstr(-flow[e] <= capacity + pref*eta, name=f"capacity_{edge.name}_-1")
    for i in range(1, 4):
        balance = (gp.quicksum(flow[e] for e, edge in enumerate(corridors) if edge.end == i)
                   - gp.quicksum(flow[e] for e, edge in enumerate(corridors) if edge.start == i))
        lp.addConstr(balance <= p[i] + pref*eta, name=f"balance_{i}_1")     # E11—E22：b_i(P)=p_i。
        lp.addConstr(-balance <= -p[i] + pref*eta, name=f"balance_{i}_-1")
        lp.addConstr(v[i] <= 1 + span*eta, name=f"voltage_upper_{i}")
        lp.addConstr(-v[i] <= -config.voltage_min_pu**2 + span*eta, name=f"voltage_lower_{i}")
    lp.addConstr(p.sum() <= config.power_limit + pref*eta, name="transformer")  # E23：Σp≤142.5。
    for e, edge in enumerate(corridors):
        for k, line in enumerate(lines):
            a = 2*(line.r_ohm_km + qr*line.x_ohm_km)*(edge.length_m/1000)/(1000*config.voltage_kv**2)
            big_m = span + a*max(t.capacity_kw for t in lines)  # 固定 M_ek，保留所有候选压降行。
            drop = v[edge.end] - v[edge.start] + a*flow[e]
            lp.addConstr(drop <= big_m*(1-x[e, k]) + span*eta, name=f"drop_upper_{edge.name}_{line.name}")
            lp.addConstr(-drop <= big_m*(1-x[e, k]) + span*eta, name=f"drop_lower_{edge.name}_{line.name}")  # E24—E53。
    lp.update()
    rows = lp.getConstrs()
    states, choices, demands = list(flow.values())+list(voltage.values()), list(x.values()), list(p.values())
    W = np.array([[lp.getCoeff(row, w) for w in states] for row in rows])
    h = np.array(lp.getAttr("RHS", rows))
    T = -np.array([[lp.getCoeff(row, z) for z in choices] for row in rows])  # 左端参数移到右端取负号。
    D = -np.array([[lp.getCoeff(row, z) for z in demands] for row in rows])
    r = -np.array([lp.getCoeff(row, eta) for row in rows])
    lp.remove(choices + demands)             # 真正的 SP 只有 8 个运行变量和 η。
    lp.update()
    return ElectricalModel(config, tuple(lines), tuple(corridors), costs, W, h, T, D, r, lp, rows)


def add_dual_cut(master, x, p, cut: DualCut):
    """CUT：β0+βxᵀx+βpᵀp≥0；正数缩放不改变半空间。"""
    scale = max(1.0, abs(cut.constant), np.max(np.abs(cut.x_coeff)), np.max(np.abs(cut.load_coeff)))
    expression = (cut.constant + gp.LinExpr(cut.x_coeff.tolist(), list(x.values()))
                  + gp.LinExpr(cut.load_coeff.tolist(), list(p.values())))
    master.addConstr(expression/scale >= 0, name=f"cut_{master.NumConstrs}")


def build_master(model: ElectricalModel, mode: str, query, cuts=(), *, fixed_x=None, load_limit=None):
    """MP1: min_cost 固定向量；MP2: max_load 固定预算；总量 MP1 用于二维前沿。"""
    master = gp.Model(mode)
    master.Params.OutputFlag = 0
    master.Params.MIPGap = 0
    master.Params.MIPGapAbs = 0
    master.Params.FeasibilityTol = 1e-9
    master.Params.IntFeasTol = 1e-9
    master.Params.Threads = 1
    master.Params.DualReductions = 0
    ne, nk = len(model.corridors), len(model.lines)
    x = master.addVars(ne, nk, vtype=GRB.BINARY, name="x")
    p = master.addVars(range(1, 4), lb=0, name="p")  # p_i⁰=0；各节点功率独立优化。
    master.addConstrs((x.sum(e, "*") <= 1 for e in range(ne)), name="one_type")  # S01—S05。
    master.addConstr(x.sum() == 3, name="tree_edges")                          # S06。
    for count in range(1, 4):
        for subset in itertools.combinations(range(1, 4), count):
            crossing = [e for e, edge in enumerate(model.corridors)
                        if (edge.start in subset) != (edge.end in subset)]
            master.addConstr(gp.quicksum(x[e, k] for e in crossing for k in range(nk)) >= 1,
                             name="connect_"+"".join(map(str, subset)))         # S07—S13。
    master.addConstr(p.sum() <= model.config.power_limit, name="power_limit") # 初始无割时也有界。
    cost = gp.LinExpr((model.cost/1000).tolist(), list(x.values()))             # 主问题费用以千元缩放。
    if mode == "max_load":
        master.addConstr(cost <= query/1000, name="budget")                  # Q2：K(x)≤B。
        master.setObjective(p.sum(), GRB.MAXIMIZE)
    else:
        if mode == "min_cost":
            master.addConstrs((p[i] == float(query[i-1]) for i in p), name="fixed_power")  # Q1：p=p̄。
        else:
            master.addConstr(p.sum() >= query, name="target_total")          # 总量查询：允许重新分配 p。
        master.setObjective(cost, GRB.MINIMIZE)
    if fixed_x is not None:
        for variable, value in zip(x.values(), fixed_x):
            variable.LB = variable.UB = float(value)                         # 固定设计的顶点证书查询。
    if load_limit is not None:
        master.addConstr(p.sum() <= load_limit, name="display_limit")        # 研究区间，例如 Σp≤140。
    for cut in cuts:
        add_dual_cut(master, x, p, cut)
    master.update()
    return master, x, p


def feasibility_oracle(model: ElectricalModel, x, p) -> SPResult:
    """更新 RHS；π=-Pi；由同一连续 LP 读取运行状态和全局有效割。"""
    rhs = model.h + model.T@x + model.D@p
    lp = model.subproblem
    lp.setAttr("RHS", model.electrical_constraints, rhs.tolist())
    lp.optimize()
    pi = -np.array(lp.getAttr("Pi", model.electrical_constraints))
    w = np.array([v.X for v in lp.getVars() if v.VarName != "eta"])
    cut = DualCut(pi, float(pi@model.h), pi@model.T, pi@model.D, float(lp.ObjVal))
    return SPResult(float(lp.ObjVal), cut, w, rhs, model.W@w-rhs)


def design_description(model, x):
    return "; ".join(f"{edge.name}:{model.lines[k].name}" for e, edge in enumerate(model.corridors)
                     for k in range(len(model.lines)) if x[e*len(model.lines)+k] > .5)


class CutSession:
    """一次查询的持续状态；step() 恰好执行一轮 MP→SP→加割/停止。"""

    def __init__(self, model, mode, query, cuts=None, **master_options):
        self.model, self.mode, self.query = model, mode, query
        self.cuts = [] if cuts is None else cuts
        self.master, self.x, self.p = build_master(model, mode, query, self.cuts, **master_options)
        self.history, self.finished = [], False

    def step(self):
        self.master.optimize()
        record = dict(iteration=len(self.history)+1, mode=self.mode, cuts_before=len(self.cuts))
        if self.master.Status == GRB.INFEASIBLE:       # 数学不可行是查询结果，不是求解器兜底。
            record.update(status="infeasible", sp=None, cut_id=None)
            self.finished = True
        else:
            x = np.rint([v.X for v in self.x.values()])
            p = np.array([v.X for v in self.p.values()])
            sp = feasibility_oracle(self.model, x, p)
            self.finished = sp.violation <= self.model.config.feasibility_tolerance
            record.update(x=x, p=p, cost_cny=float(self.model.cost@x), total_kw=float(p.sum()),
                          sp=sp, violation=sp.violation, design=design_description(self.model, x),
                          status="optimal" if self.finished else "cut", cut_id=None)
            if not self.finished:
                self.cuts.append(sp.cut)
                record["cut_id"] = len(self.cuts)
                add_dual_cut(self.master, self.x, self.p, sp.cut)
        self.history.append(record)
        return record

    def result(self):
        last = self.history[-1]
        if last["status"] == "infeasible":
            return SolveResult(self.mode, "infeasible", math.inf if self.mode != "max_load" else -math.inf,
                               np.array([]), np.array([]), math.inf, self.history)
        value = last["total_kw"] if self.mode == "max_load" else last["cost_cny"]
        return SolveResult(self.mode, "optimal", value, last["x"], last["p"], last["cost_cny"], self.history)


def solve_by_cuts(model, mode, query, cuts=None, **master_options):
    """完整求解一次查询；逐轮演示直接使用 CutSession.step。"""
    session = CutSession(model, mode, query, cuts, **master_options)
    while not session.finished:
        session.step()
    result = session.result()
    session.master.dispose()
    return result


def enumerate_designs(model):
    """仅枚举合法树与线型；不查询电气能力，不产生对偶割。"""
    designs, nk = [], len(model.lines)
    for edges in itertools.combinations(range(len(model.corridors)), 3):
        reached = {0}
        for _ in range(3):
            for e in edges:
                edge = model.corridors[e]
                if edge.start in reached or edge.end in reached:
                    reached.update((edge.start, edge.end))
        if len(reached) != 4:
            continue
        for types in itertools.product(range(nk), repeat=3):
            x = np.zeros(model.nx)
            for e, k in zip(edges, types):
                x[e*nk+k] = 1
            designs.append(Design(x, float(model.cost@x), design_description(model, x)))
    return sorted(designs, key=lambda z: (z.cost_cny, z.description))


def design_inequalities(model, x):
    """独立树上模型 Gp≤b：下游负荷、线路容量、路径压降、配变容量。"""
    nk = len(model.lines)
    selected = {e: int(np.argmax(row)) for e, row in enumerate(x.reshape(-1, nk)) if row.sum() > .5}
    adj = {i: [] for i in range(4)}
    for e in selected:
        edge = model.corridors[e]
        adj[edge.start].append((edge.end, e))
        adj[edge.end].append((edge.start, e))
    parents, order = {}, [0]
    for i in order:
        for j, e in adj[i]:
            if j != 0 and j not in parents:
                parents[j] = (i, e)
                order.append(j)
    downstream, flows = np.vstack([np.zeros(3), np.eye(3)]), {}
    for j in reversed(order[1:]):
        i, e = parents[j]
        flows[e] = downstream[j].copy()
        downstream[i] += downstream[j]
    G, b, drops = [np.ones(3)], [model.config.power_limit], {0: np.zeros(3)}
    qr = math.tan(math.acos(model.config.power_factor))
    for j in order[1:]:
        i, e = parents[j]
        line, edge = model.lines[selected[e]], model.corridors[e]
        a = 2*(line.r_ohm_km+qr*line.x_ohm_km)*(edge.length_m/1000)/(1000*model.config.voltage_kv**2)
        drops[j] = drops[i] + a*flows[e]
        G.extend([flows[e], drops[j]])
        b.extend([line.capacity_kw, 1-model.config.voltage_min_pu**2])
    return np.array(G), np.array(b)


def evaluate_design(model, x):
    """独立验证固定树的最大总功率；不使用 Benders 割。"""
    G, b = design_inequalities(model, x)
    lp = gp.Model("independent_tree_capacity")
    lp.Params.OutputFlag = 0
    p = lp.addMVar(3, lb=0)
    lp.addConstr(G@p <= b)
    lp.setObjective(p.sum(), GRB.MAXIMIZE)
    lp.optimize()
    answer = float(lp.ObjVal), p.X.copy()
    lp.dispose()
    return answer


def frontier_table(model, designs):
    """枚举 LP 基准的完整费用—总功率台阶，仅用于独立比较。"""
    rows, best = [], -math.inf
    for design in designs:
        load, p = evaluate_design(model, design.x)
        if load > best+1e-7:
            row = dict(cost_cny=design.cost_cny, total_kw=load, p1=p[0], p2=p[1], p3=p[2], design=design.description)
            if rows and rows[-1]["cost_cny"] == design.cost_cny:
                rows[-1] = row
            else:
                rows.append(row)
            best = load
    return pd.DataFrame(rows)


def discover_frontier_by_cuts(model, cuts=None, load_limit=140.0):
    """MP2→总量 MP1→预算减一费用粒度，精确恢复研究区间内全部台阶。"""
    pool = [] if cuts is None else cuts
    quantum = math.gcd(*np.rint(model.cost[model.cost > 0]).astype(int))
    budget, rows = float(model.cost.sum()), []
    while budget >= 0:
        capacity = solve_by_cuts(model, "max_load", budget, pool, load_limit=load_limit)
        if capacity.status == "infeasible":
            break
        minimum = solve_by_cuts(model, "min_cost_total", capacity.value, pool, load_limit=load_limit)
        rows.append(dict(cost_cny=minimum.cost_cny, total_kw=capacity.value,
                         p1=minimum.p[0], p2=minimum.p[1], p3=minimum.p[2],
                         design=design_description(model, minimum.x)))
        budget = round(minimum.cost_cny)-quantum
    return pd.DataFrame(rows).sort_values("cost_cny").reset_index(drop=True)
