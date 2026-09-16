# 数学阅读入口：README.md §0 变量—符号表 → §1 完整规划模型 → §2 电气 E01—E53 → §3 紧凑式。
# 行尾编号与 README 的公式编号一致；E 为电气约束，S 为树结构约束，D 为变量域，Q 为查询约束。
# O1/O2 为两种规划目标，O3 为最小违反 LP 目标，O4 为凸包辅助目标；CUT 为对偶可行性割。
# 所有数值展开以默认参数为准。P 是 kW 有功潮流；v 是平方电压；alpha 是 λ；x 是建设选择。
# 53 只统计电气行，不包括 η≥0、主问题二元变量界、13 条树约束、预算/固定倍率和动态电气割。
# 注释中的“无对应优化约束”表示程序组织、求解器配置或结果处理，不把这些语句误写成物理约束。
"""四节点单台区：负荷倍率—建设预算可规划域与对偶可行性割。

教学用合成参数；不是江口数据，也不是工程导线选型建议。
连续模型为固定功率因数、无损有功平衡、线性平方电压降。
算法：15 个线型二元变量的 MILP 主问题 + 最小违反 LP 子问题。
穷举只用于独立校验，不参与主问题选解或生成对偶割。
使用 gurobipy 直接建立变量、约束和目标；需要可用的 Gurobi 许可证。
"""
from __future__ import annotations  # 语言设置：延迟解析类型注解；无对应优化约束。

import itertools  # 组合工具：枚举节点子集、三条边组合和线型笛卡尔积。
import json  # 结果序列化：保存核心实验摘要；无对应优化约束。
import math  # 数学运算：计算 tan(arccos(PF))、无穷大和费用最大公约数。
from dataclasses import dataclass  # 数据容器：保存参数、矩阵和求解结果。
from pathlib import Path  # 文件路径：定位并保存实验输出。
from typing import Literal, Sequence  # 类型标注：限定查询模式及序列类型；不限制数学可行域。

import numpy as np  # 数值数组：表示 x、w、W、h、T、d、r 及向量内积。
import pandas as pd  # 表格工具：保存前沿、迭代轨迹及验证结果。
import gurobipy as gp  # Gurobi 建模接口：构建 MILP 主问题和连续 LP 子问题。
from gurobipy import GRB  # Gurobi 常量：变量类型、目标方向、求解状态和无穷界。


@dataclass(frozen=True)  # 不可变参数数据类；不增加决策变量。
class LineType:  # 线型参数 k∈{L,M,H}，符号表见 README §0.2。
    name: str  # 线型标签 k；不是约束名称。
    r_ohm_km: float  # R_k：单位长度电阻，Ω/km。
    x_ohm_km: float  # X_k：单位长度电抗，Ω/km；不要与二元向量 x 混淆。
    capacity_kw: float  # C_k：有功容量上限，kW。
    cable_cny_m: float  # p_k：电缆单价，CNY/m。


@dataclass(frozen=True)  # 不可变走廊参数数据类。
class Corridor:  # 候选走廊 e=(i_e,j_e)，符号表见 README §0.2。
    name: str  # 走廊标签，例如 01；数组索引 e=0 对应 01。
    start: int  # i_e：固定参考方向的起点；实际潮流可反向。
    end: int  # j_e：固定参考方向的终点。
    length_m: float  # ℓ_e：走廊长度，m。
    existing: bool  # b_e：已有走廊标志；只影响费用，不强制保留该边。


@dataclass(frozen=True)  # 不可变算例参数数据类。
class DemoConfig:  # 四节点合成算例参数；当前循环硬编码节点 0、1、2、3。
    loads_kw: tuple[float, ...] = (0.0, 25.0, 20.0, 15.0)  # p_i⁰=(0,25,20,15) kW；实际有功负荷为 λp_i⁰。
    power_factor: float = 0.95  # cosφ=0.95；无功/有功比例 q=tanφ。
    voltage_kv: float = 0.4  # U_N=0.4 kV；用于线性平方电压降系数的标幺换算。
    voltage_min_pu: float = 0.93  # u_min=0.93；平方电压下界为 0.8649。
    transformer_kva: float = 150.0  # S_tr=150 kVA；固定功率因数下有功上限为 142.5 kW。
    new_corridor_cny_m: float = 45.0  # p_new=45 CNY/m；仅新走廊计取附加建设费。
    lambda_search_max: float = 3.0  # λ_search=3；主问题的搜索上界，不是电气极限。
    feasibility_tolerance: float = 1e-8  # ε_feas=10⁻⁸；当最优违反量 η* 不超过此值时接受候选。


DEFAULT_LINES = (  # 按 L、M、H 排列线型，决定 k=0、1、2 及矩阵列顺序。
    LineType("L", 1.15, 0.08, 35.0, 28.0),  # k=0：L 型，(R_L,X_L,C_L,p_L)=(1.15,0.08,35,28)。
    LineType("M", 0.62, 0.08, 65.0, 43.0),  # k=1：M 型，(R_M,X_M,C_M,p_M)=(0.62,0.08,65,43)。
    LineType("H", 0.32, 0.08, 100.0, 65.0),  # k=2：H 型，(R_H,X_H,C_H,p_H)=(0.32,0.08,100,65)。
)  # 结束三类线型参数表。
DEFAULT_CORRIDORS = (  # 按 01、12、13、02、23 排列走廊，决定 e=0、1、2、3、4。
    Corridor("01", 0, 1, 220.0, True),  # e=0：0→1，220 m，已有；保留 L 型新增费用为 0。
    Corridor("12", 1, 2, 180.0, True),  # e=1：1→2，180 m，已有。
    Corridor("13", 1, 3, 240.0, True),  # e=2：1→3，240 m，已有。
    Corridor("02", 0, 2, 300.0, False),  # e=3：0→2，300 m，新建时需另付走廊费用。
    Corridor("23", 2, 3, 160.0, False),  # e=4：2→3，160 m，新建时需另付走廊费用。
)  # 结束五条候选走廊参数表。


@dataclass  # 模型数据容器；保存数学参数及可复用 LP 对象。
class ElectricalModel:  # 统一电气模型 Ww≤h+Tx+dλ+rη，见 README §3。
    config: DemoConfig  # 算例标量参数和基准负荷 p_i⁰。
    lines: tuple[LineType, ...]  # 线型集合及参数 (R_k,X_k,C_k,p_k)。
    corridors: tuple[Corridor, ...]  # 走廊集合及参数 (i_e,j_e,ℓ_e,b_e)。
    cost: np.ndarray  # c∈R¹⁵：各 (e,k) 的建设费用，CNY；不是主问题中以千元计的 cost。
    W: np.ndarray  # W∈R^(53×8)：有功潮流和非根平方电压的系数。
    h: np.ndarray  # h∈R⁵³：常数右端；已包含根电压 v₀=1 移项产生的 ±1。
    T: np.ndarray  # T∈R^(53×15)：二元选择 x 移到右端后的系数。
    d: np.ndarray  # d∈R⁵³：负荷倍率 λ 移到右端后的系数。
    scales: np.ndarray  # r∈R⁵³：统一违反量 η 的行尺度；功率行 60、电压行 0.1351。
    subproblem: gp.Model  # 最小违反 LP：min η；最终只有 5 个 P、3 个 v、1 个 η。
    electrical_constraints: list[gp.Constr]  # 53 个电气约束对象，用于更新 RHS 和读取 Pi。

    @property  # 只读计算属性；不增加模型约束。
    def nx(self) -> int:  # 返回建设二元变量的数量 n_x。
        return self.cost.size  # n_x=|E||K|=5×3=15。

    @property  # 只读计算属性；不增加模型约束。
    def transformer_limit(self) -> float:  # 配变单独给出的负荷倍率极限。
        return self.config.transformer_kva * self.config.power_factor / sum(self.config.loads_kw)  # λ_tr=S_tr cosφ/Σp_i⁰=142.5/60=2.375；对应 E23。


@dataclass  # 一个完整整数设计的独立解析结果。
class Design:  # 设计 z=(x,K(x),Λ(x),瓶颈,描述)。
    x: np.ndarray  # 15 维建设方案 x，按 e 优先、k 次优先展开。
    cost_cny: float  # K(x)=cᵀx，单位 CNY。
    lambda_max: float  # Λ(x)：该固定树的电气最大倍率，不包含 λ_search 截断。
    bottleneck: str  # 达到最小倍率上界的约束标签；并列时按代码顺序保留一个。
    description: str  # 人可读的“走廊:线型”列表；不参与优化。


@dataclass  # 保存一条由连续 LP 对偶解产生的全局电气可行性割。
class DualCut:  # 割形式 β₀+β_xᵀx+β_λλ≥0，见 README §4。
    # pi^T(h + T x + d lambda) >= 0；统一模型下跨预算、跨倍率均有效。
    pi: np.ndarray  # π∈R⁵³：非负对偶向量，满足 Wᵀπ=0、rᵀπ≤1。
    constant: float  # β₀=πᵀh：割的常数项。
    x_coeff: np.ndarray  # β_x=Tᵀπ：15 维割系数；代码 pi@T 存为一维数组。
    lambda_coeff: float  # β_λ=πᵀd：负荷倍率的割系数。
    source_violation: float  # ν=η*：生成该割时候选点的最小违反量。


@dataclass  # 保存一次固定倍率或固定预算查询的结果。
class SolveResult:  # 求解结果容器；不增加数学约束。
    mode: str  # 模式：min_cost 对应 O1，max_lambda 对应 O2。
    query: float  # 查询输入：min_cost 为 λ̄，max_lambda 为 B̄（CNY）。
    status: str  # 求解状态 optimal 或 infeasible。
    value: float  # 最优目标：O1 为最低费用（CNY），O2 为最大倍率。
    x: np.ndarray  # 求得的整数设计 x；不可行返回的零向量仅为占位值。
    cost_cny: float  # 该设计的费用 K(x)，CNY。
    lambda_value: float  # 该次查询采用的 λ，不一定等于设计的 Λ(x)。
    trace: pd.DataFrame  # 收敛轨迹：累计割数、目标上下界和电气违反量。
    new_cuts: int  # 本次查询新增割数；等于结束割池长度减初始长度。


def build_model(config: DemoConfig = DemoConfig(), lines: tuple[LineType, ...] = DEFAULT_LINES,  # 建立参数化电气模型，默认数据的所有电气行见 README §2（E01—E53）。
                corridors: tuple[Corridor, ...] = DEFAULT_CORRIDORS) -> ElectricalModel:  # 可传入自定义参数；节点数仍固定为四，文档数值展开对应默认值。
    """用 Gurobi 写电气约束，并提取 W、h、T、d 供 RHS 更新及对偶割使用。"""
    ne, nk = len(corridors), len(lines)  # n_e=5、n_k=3；矩阵按走廊优先、线型次优先排列。
    total = sum(config.loads_kw)  # D=Σ_i p_i⁰=60 kW；也是所有功率约束的松弛尺度。
    span = 1 - config.voltage_min_pu ** 2  # Δv=1-u_min²=0.1351；也是所有平方电压约束的松弛尺度。
    qratio = math.tan(math.acos(config.power_factor))  # q=tan(arccos(cosφ))=√39/19；固定无功/有功比例。
    costs = np.array([  # 费用向量 c 的逐项构造，完整费用式见 README（COST）。
        0.0 if edge.existing and k == 0 else  # 已有走廊选择 k=0（默认 L 型）时 c_ek=0。
        edge.length_m * (line.cable_cny_m + (0.0 if edge.existing else config.new_corridor_cny_m))  # 其余 c_ek=ℓ_e[p_k+(1-b_e)p_new]，单位 CNY。
        for edge in corridors for k, line in enumerate(lines)  # 依次生成 c_01,L、c_01,M、c_01,H、c_12,L 等 15 项。
    ])  # 形成一维数组 c∈R¹⁵。
    lp = gp.Model("electrical_feasibility")  # 创建最小违反子问题 LP，不是完整规划 MILP。
    lp.Params.OutputFlag = 0  # 关闭求解日志；无数学约束。
    lp.Params.FeasibilityTol = 1e-9  # 原始约束可行性容差 10⁻⁹；区别于外层接受阈值 10⁻⁸。
    lp.Params.OptimalityTol = 1e-9  # LP 对偶可行性容差 10⁻⁹；数值精度设置。
    lp.Params.Method = 1  # 选择对偶单纯形法求 LP；便于相邻 RHS 查询复用基础。
    lp.Params.Threads = 1  # 单线程求解设置；不改变目标和可行域。
    flow = lp.addVars(ne, lb=-GRB.INFINITY, name="P")  # P_e∈R，e=0,…,4；允许参考方向的反向潮流，容量界由 E01—E10 给出。
    voltage = lp.addVars(range(1, 4), lb=-GRB.INFINITY, name="v")  # v₁,v₂,v₃∈R；显式上下界由 E13/E14、E17/E18、E21/E22 给出。
    v = {0: 1.0, **voltage}  # v₀=1 是常数；仅 v₁、v₂、v₃ 是决策变量，v_i=u_i²。
    eta = lp.addVar(lb=0.0, name="eta")  # η≥0：无量纲统一违反量；仅子问题使用，原始规划模型取 η=0。
    x = lp.addVars(ne, nk, lb=-GRB.INFINITY, name="x")  # 临时连续占位 x_ek：仅供提取 T，稍后删除；真正二元变量在 build_master。
    alpha = lp.addVar(lb=-GRB.INFINITY, name="lambda")  # 临时连续占位 λ（Python 名 alpha）：仅供提取 d，稍后删除。
    lp.setObjective(eta, GRB.MINIMIZE)  # （O3）min η；η*=0 表示给定 (x̄,λ̄) 的电气约束可行。

    # 所有电气行写成 <=；潮流和电压的界也显式建约束。
    for e, edge in enumerate(corridors):  # 遍历 5 条边；每条边按正、负方向各添加一行，共 E01—E10。
        capacity = gp.quicksum(line.capacity_kw * x[e, k] for k, line in enumerate(lines))  # C_e(x)=35x_eL+65x_eM+100x_eH；未选边时容量为 0。
        lp.addConstr(flow[e] <= capacity + total * eta, name=f"capacity_{edge.name}_1")  # （E01/E03/E05/E07/E09）P_e≤C_e(x)+60η；η=0 时限制正向输送。
        lp.addConstr(-flow[e] <= capacity + total * eta, name=f"capacity_{edge.name}_-1")  # （E02/E04/E06/E08/E10）−P_e≤C_e(x)+60η；两行合为 |P_e|≤C_e(x)+60η。
    for node in range(1, 4):  # 遍历三个负荷节点；每节点按平衡正行、平衡负行、电压上界、下界添加。
        balance = (gp.quicksum(flow[e] for e, edge in enumerate(corridors) if edge.end == node)  # b_i(P) 的流入项 Σ_{e:j_e=i}P_e；P 按固定参考方向带符号。
                   - gp.quicksum(flow[e] for e, edge in enumerate(corridors) if edge.start == node))  # 减去流出项 Σ_{e:i_e=i}P_e；所以 b₁=P₀₁−P₁₂−P₁₃ 等。
        demand = config.loads_kw[node] * alpha  # p_i(λ)=p_i⁰λ，即节点 1/2/3 的负荷分别为 25λ/20λ/15λ kW。
        lp.addConstr(balance <= demand + total * eta, name=f"balance_{node}_1")  # （E11/E15/E19）b_i(P)≤p_i⁰λ+60η；有功平衡正行。
        lp.addConstr(-balance <= -demand + total * eta, name=f"balance_{node}_-1")  # （E12/E16/E20）−b_i(P)≤−p_i⁰λ+60η；η=0 时两行严格表达等式。
        lp.addConstr(v[node] <= 1.0 + span * eta, name=f"voltage_upper_{node}")  # （E13/E17/E21）v_i≤1+0.1351η；原始模型上界为 1。
        lp.addConstr(-v[node] <= -config.voltage_min_pu**2 + span * eta, name=f"voltage_lower_{node}")  # （E14/E18/E22）−v_i≤−0.8649+0.1351η；原始模型下界为 0.8649。
    lp.addConstr(total * alpha <= config.transformer_kva * config.power_factor + total * eta,  # （E23）60λ≤142.5+60η；η=0 时 λ≤2.375。
                 name="transformer")  # E23 的求解器行名为 transformer；根节点有功注入由负荷总和隐式确定。
    for e, edge in enumerate(corridors):  # 遍历五条走廊；每条有三类线型、每线型两条电压降行。
        for k, line in enumerate(lines):  # 按 L、M、H 顺序展开 k；共 5×3×2=30 行（E24—E53）。
            coefficient = (2 * (line.r_ohm_km + qratio * line.x_ohm_km)  # a_ek=2(R_k+qX_k)(ℓ_e/1000)/(1000U_N²)，单位 1/kW。
                           * (edge.length_m / 1000) / (1000 * config.voltage_kv**2))  # ℓ_e/1000 把 m 换为 km；分母的 1000 配合 kW、kV 完成标幺换算。
            big_m = span + coefficient * max(t.capacity_kw for t in lines)  # M_ek=Δv+a_ek max_k C_k=0.1351+100a_ek；放松未选线型的压降关系。
            drop = v[edge.end] - v[edge.start] + coefficient * flow[e]  # δ_ek=v_j−v_i+a_ek P_e；选中线型且 η=0 时应为 0。
            lp.addConstr(drop <= big_m * (1 - x[e, k]) + span * eta,  # （E24/E26/…/E52）δ_ek≤M_ek(1−x_ek)+0.1351η。
                         name=f"drop_upper_{edge.name}_{line.name}")  # 上行名称 drop_upper_走廊_线型；逐条编号、展开式见 README §2。
            lp.addConstr(-drop <= big_m * (1 - x[e, k]) + span * eta,  # （E25/E27/…/E53）−δ_ek≤M_ek(1−x_ek)+0.1351η。
                         name=f"drop_lower_{edge.name}_{line.name}")  # 下行名称 drop_lower_走廊_线型；x_ek=1、η=0 时与上行合为电压降等式。
    lp.update()  # 提交建模修改，以便读取规范化后的系数和 RHS。
    rows = lp.getConstrs()  # rows=(E01,…,E53)；实际顺序是 Gurobi 添加顺序。
    states = list(flow.values()) + list(voltage.values())  # w=(P₀₁,P₁₂,P₁₃,P₀₂,P₂₃,v₁,v₂,v₃)ᵀ，共 8 个自由状态变量。
    choices = list(x.values())  # 临时 x 列顺序与费用 c 一致；e*nk+k 映射到一维列。
    W = np.array([[lp.getCoeff(row, w) for w in states] for row in rows])  # W_mj=第 m 行中 w_j 的左端系数，形状 (53,8)。
    h = np.array(lp.getAttr("RHS", rows))  # h_m=规范化 RHS；如 E24 为 M_01,L+1，E25 为 M_01,L−1。
    T = -np.array([[lp.getCoeff(row, choice) for choice in choices] for row in rows])  # T_mj=−(临时 x_j 的左端系数)；移项到右端必须取负号。
    d = -np.array([lp.getCoeff(row, alpha) for row in rows])  # d_m=−(临时 λ 的左端系数)；如 E11 为 +25，E12 为 −25，E23 为 −60。
    scales = -np.array([lp.getCoeff(row, eta) for row in rows])  # r_m=−(η 的左端系数)；功率行 60，平方电压行 0.1351。
    # x、lambda 只用于书写参数项；移除后，LP 仅含 P、v、eta，参数通过 RHS 传入。
    lp.remove(choices + [alpha])  # 移除 15 个 x 占位变量和 λ 占位变量；固定参数由 oracle 写入 RHS。
    lp.update()  # 提交删除；最终 LP 是 Ww−rη≤rhs，η≥0，目标 min η。
    return ElectricalModel(config, lines, corridors, costs, W, h, T, d, scales, lp, rows)  # 返回矩阵、参数和 LP；首次实际电气查询须先设置 rhs=h+Tx̄+dλ̄。


def add_dual_cut(master: gp.Model, x: gp.tupledict, alpha: gp.Var, cut: DualCut) -> None:  # （CUT）将 β₀+β_xᵀx+β_λλ≥0 加入当前 MILP 主问题。
    """把 pi^T(h + T x + d lambda) >= 0 加入整数主问题。"""
    scale = max(1.0, np.max(np.abs(cut.x_coeff)), abs(cut.lambda_coeff), abs(cut.constant))  # σ=max(1,||β_x||∞,|β_λ|,|β₀|)>0；仅进行等价数值缩放。
    expression = gp.LinExpr(cut.x_coeff.tolist(), list(x.values())) + cut.lambda_coeff * alpha  # 构造 β_xᵀx+β_λλ；不含常数项 β₀。
    master.addConstr((expression + cut.constant) / scale >= 0)  # （CUT）(β_xᵀx+β_λλ+β₀)/σ≥0；与未缩放割完全等价。


def build_master(model: ElectricalModel, mode: Literal["min_cost", "max_lambda"], query: float,  # 构造只含 x、λ 的整数主问题（MP）；电气状态 w 被可行性割投影消去。
                 cuts: Sequence[DualCut] = ()) -> tuple[gp.Model, gp.tupledict, gp.Var]:  # 输入模式、查询值和已有割池；返回模型及主决策变量。
    """Gurobi 整数主问题：线型选择、生成树、预算及累计电气割。"""
    master = gp.Model(f"planning_{mode}")  # 创建与查询模式对应的 MILP 主问题。
    master.Params.OutputFlag = 0  # 关闭日志；无数学约束。
    master.Params.MIPGap = 0.0  # 相对 MIP 间隙目标为 0；要求证明最优至求解器数值容差。
    master.Params.MIPGapAbs = 0.0  # 绝对 MIP 间隙目标为 0。
    master.Params.FeasibilityTol = 1e-9  # 线性约束容差 10⁻⁹。
    master.Params.IntFeasTol = 1e-9  # 二元变量整数容差 10⁻⁹。
    master.Params.OptimalityTol = 1e-9  # 连续最优性容差 10⁻⁹。
    master.Params.DualReductions = 0  # 关闭对偶化简，便于区分不可行状态；无新增约束。
    master.Params.Threads = 1  # 主问题使用单线程。
    ne, nk = len(model.corridors), len(model.lines)  # 走廊数 n_e 和线型数 n_k，默认 5 和 3。
    x = master.addVars(ne, nk, vtype=GRB.BINARY, name="x")  # （D1）真正的建设决策 x_ek∈{0,1}，共 15 个。
    alpha = master.addVar(lb=0.0, ub=model.config.lambda_search_max, name="lambda")  # （D2）0≤λ≤λ_search=3；Python 变量名 alpha。
    master.addConstrs((x.sum(e, "*") <= 1 for e in range(ne)), name="one_type")  # （S01—S05）每条边 x_eL+x_eM+x_eH≤1；允许该边不选。
    master.addConstr(x.sum() == 3, name="tree_edges")  # （S06）全部 15 个 x 之和等于 3，即四节点运行树恰好选择三条边。
    # 连通割 + 三条边 = 四节点生成树。
    for count in range(1, 4):  # 枚举不含根节点的非空子集大小 1、2、3。
        for subset in itertools.combinations(range(1, 4), count):  # S={1},{2},{3},{1,2},{1,3},{2,3},{1,2,3}，共七个。
            crossing = [e for e, edge in enumerate(model.corridors)  # 构造跨越 S 与其补集的候选边集合 δ(S)。
                        if (edge.start in subset) != (edge.end in subset)]  # 异或判定：边恰有一个端点在 S 内；方向不影响连通性。
            master.addConstr(gp.quicksum(x[e, k] for e in crossing for k in range(nk)) >= 1,  # （S07—S13）Σ_{e∈δ(S),k}x_ek≥1；七条连通割见 README §1.3。
                             name="connect_" + "".join(map(str, subset)))  # 命名 connect_1、connect_2、connect_3、connect_12、connect_13、connect_23、connect_123。
    cost = gp.LinExpr((model.cost / 1000).tolist(), list(x.values()))  # K(x)/1000：主问题费用以千元表示；model.cost 本身仍为 CNY。
    if mode == "min_cost":  # 选择固定倍率、最低建设费模型 O1。
        master.addConstr(alpha == query, name="fixed_lambda")  # （Q1）λ=λ̄=query；与 D2 共同生效。
        master.setObjective(cost, GRB.MINIMIZE)  # （O1）min K(x)/1000；与 min K(x) 的最优设计相同。
    else:  # 选择固定预算、最大负荷倍率模型 O2。
        master.addConstr(cost <= query / 1000, name="budget")  # （Q2）K(x)/1000≤B̄/1000；query 是 CNY。
        master.setObjective(alpha, GRB.MAXIMIZE)  # （O2）max λ；与 O1 是两种独立查询，不是同时优化。
    for cut in cuts:  # 遍历可跨查询复用的电气割 π^t。
        add_dual_cut(master, x, alpha, cut)  # 逐条添加（CUT），得到当前电气可行域的外近似。
    master.update()  # 提交主问题模型；未包含 P、v 或 η。
    return master, x, alpha  # 返回主问题和 x、λ，供 solve_by_cuts 迭代优化。


def evaluate_design(model: ElectricalModel, x: np.ndarray) -> Design:  # （TREE-EVAL）给定合法生成树 x，独立计算 Λ(x)；公式见 README §5。
    """独立解析：固定树后，下游负荷确定潮流，容量和路径压降直接给出倍率上限。"""
    nk, ne = len(model.lines), len(model.corridors)  # 读取线型数与走廊数。
    selected = np.asarray(x).reshape(ne, nk)  # 把一维 x 还原成 5×3 的选择表。
    chosen = {e: int(np.argmax(selected[e])) for e in range(ne) if selected[e].sum() > .5}  # 提取已选走廊及其唯一线型 k(e)；阈值 0.5 用于识别整数设计。
    adj: dict[int, list[tuple[int, int]]] = {i: [] for i in range(4)}  # 构造无向邻接表；运行树的根方向可能与 Corridor.start/end 不同。
    for e in chosen:  # 遍历已选边。
        edge = model.corridors[e]  # 读取边的端点和长度参数。
        adj[edge.start].append((edge.end, e)); adj[edge.end].append((edge.start, e))  # 双向写入邻接表，保证从根节点能沿任意参考方向遍历。
    parents: dict[int, tuple[int, int]] = {}; order = [0]; seen = {0}  # 从根 0 开始；parents 记录父节点/边，order 记录根到叶遍历顺序。
    for i in order:  # 按不断扩展的 order 从根向外遍历。
        for j, e in adj[i]:  # 读取 i 的邻接节点 j 和边 e。
            if j not in seen:  # 只访问尚未访问的节点，建立树的父子关系。
                seen.add(j); parents[j] = (i, e); order.append(j)  # 记录 j 的父节点 i、连接边 e，并加入遍历队列。
    downstream = np.asarray(model.config.loads_kw, dtype=float).copy(); flows = {}  # 初始化下游负荷和支路负荷表；此处计算的是 λ=1 时的功率大小。
    for j in reversed(order[1:]):  # 按叶到根顺序累计负荷。
        i, e = parents[j]; flows[e] = downstream[j]; downstream[i] += downstream[j]  # f_e=子树负荷总和；向父节点累计；flows[e] 是根向外的非负大小，并非带符号 P_e。
    limit = model.transformer_limit; bottleneck = "transformer"  # 先取 Λ=λ_tr=2.375，再逐个检查是否有更紧的容量/电压限制。
    drops = {0: 0.0}; qr = math.tan(math.acos(model.config.power_factor))  # Δ₀=0，q=tanφ；drops 存 λ=1 的累计平方电压降。
    for e, k in chosen.items():  # 遍历选中线型。
        cap = model.lines[k].capacity_kw / flows[e]  # λ_e^cap=C_{k(e)}/f_e；来自 E01—E10 的 |P_e|≤C_{k(e)}。
        if cap < limit - 1e-12:  # 用 10⁻¹² 判断容量上界是否严格收紧当前 Λ。
            limit = cap; bottleneck = f"capacity:{model.corridors[e].name}"  # 更新 Λ 及容量瓶颈标签。
    for j in order[1:]:  # 按根到叶顺序计算各节点的路径压降。
        i, e = parents[j]; edge, line = model.corridors[e], model.lines[chosen[e]]  # 读取父节点、走廊及该走廊的已选线型。
        drops[j] = drops[i] + 2*(line.r_ohm_km+qr*line.x_ohm_km)*(edge.length_m/1000)*flows[e]/(1000*model.config.voltage_kv**2)  # Δ_j=Δ_i+a_{e,k(e)}f_e；于是 v_j(λ)=1−λΔ_j。
        cap = (1-model.config.voltage_min_pu**2) / drops[j]  # λ_j^v=(1−u_min²)/Δ_j；来自平方电压下界 E14/E18/E22。
        if cap < limit - 1e-12:  # 判断电压下限给出的倍率上界是否更紧；cap 表示倍率上界，不是线路容量。
            limit = cap; bottleneck = f"voltage:{j}"  # 更新 Λ 及电压瓶颈标签。
    description = "; ".join(f"{model.corridors[e].name}:{model.lines[k].name}" for e, k in sorted(chosen.items()))  # 生成设计说明，仅用于输出和稳定排序。
    return Design(np.rint(x).astype(float), float(model.cost@x), float(limit), bottleneck, description)  # 返回 x、K(x)=cᵀx、Λ(x)=min(λ_tr,各容量倍率,各电压倍率) 及瓶颈。


def enumerate_designs(model: ElectricalModel) -> list[Design]:  # （ENUM）独立枚举 X_tree；不用于选择主问题候选或生成电气割。
    """枚举全部生成树与线型；完整真值仅用于独立验证。"""
    nk, ne = len(model.lines), len(model.corridors); designs = []  # 初始化线型/走廊数量及设计集合。
    for edge_set in itertools.combinations(range(ne), 3):  # 从五条候选边中选三条，共 C(5,3)=10 个边集。
        connected = {0}  # 初始化根可达节点集合为 {0}。
        # 独立遍历判断连通性，不依赖主问题约束矩阵。
        for _ in range(3):  # 四节点最多经过三轮可达传播覆盖整棵树。
            for e in edge_set:  # 遍历候选三边集合。
                edge = model.corridors[e]  # 读取边参数。
                if edge.start in connected or edge.end in connected:  # 若至少一端可达，则另一端也可达。
                    connected.update((edge.start, edge.end))  # 扩大根连通分量。
        if len(connected) != 4:  # 若根分量未覆盖四节点，则该三边集合不是生成树。
            continue  # 跳过不连通边集；默认留下 8 棵树。
        for choices in itertools.product(range(nk), repeat=3):  # 每棵树枚举 3³=27 种线型组合。
            x = np.zeros(model.nx)  # 初始化长度 15 的零向量 x。
            for e, k in zip(edge_set, choices):  # 为三条选中边配对线型。
                x[e*nk+k] = 1.0  # 令 x_ek=1；一维下标 e*nk+k，与 c、T 的列顺序一致。
            designs.append(evaluate_design(model, x))  # 独立解析当前设计的费用、最大倍率和瓶颈。
    return designs  # 默认返回 8×27=216 个设计。


def frontier_table(designs: Sequence[Design]) -> pd.DataFrame:  # （FRONT）从全部独立设计提取费用—倍率 Pareto 阶梯边界。
    """保留费用递增且倍率严格改进的全部方案，精确表示阶梯边界。"""
    rows = []; best = -math.inf  # rows 存前沿点；best 存此前费用以内的最大倍率。
    for design in sorted(designs, key=lambda z: (z.cost_cny, -z.lambda_max, z.description)):  # 按费用升序、倍率降序、描述排序；同费用先考察容量最大方案。
        if design.lambda_max > best + 1e-10:  # 仅保留严格改善最大倍率的设计，阈值为 10⁻¹⁰。
            rows.append({"cost_cny": design.cost_cny, "lambda_max": design.lambda_max,  # 保存 (K(x),Λ(x)) 前沿拐点。
                         "bottleneck": design.bottleneck, "design": design.description})  # 附带瓶颈和建设方案描述；不改变前沿坐标。
            best = design.lambda_max  # 更新截至当前费用的最好倍率。
    return pd.DataFrame(rows)  # 返回前沿数据表。


def feasibility_oracle(model: ElectricalModel, x: np.ndarray, alpha: float) -> tuple[float, DualCut]:  # （O3/DUAL）固定 (x̄,λ̄)，求电气违反量与可行性割，见 README §4。
    """更新连续 LP 的 RHS，通过 Gurobi 约束 Pi 属性读取对偶乘子。"""
    rhs = model.h + model.T@x + model.d*alpha  # b(x̄,λ̄)=h+Tx̄+dλ̄∈R⁵³；这里只更新参数右端。
    lp = model.subproblem  # 复用 build_model 已建立的 9 变量连续 LP。
    lp.setAttr("RHS", model.electrical_constraints, rhs.tolist())  # 把 53 条电气行统一更新为 Ww−rη≤b(x̄,λ̄)。
    lp.optimize()  # 求 min η；正尺度和 η≥0 保证默认有限参数下子问题有可行最优解。
    # 最小化 LP 的 <= 约束满足 Pi <= 0；理论中的 pi = -Pi >= 0。
    pi = -np.array(lp.getAttr("Pi", model.electrical_constraints))  # （DUAL）π=−Pi≥0；Gurobi 最小化 LP 的 ≤ 行对应 Pi≤0。
    violation = lp.ObjVal  # ν(x̄,λ̄)=η*=LP 目标值；并非建设费用。
    cut = DualCut(pi, float(pi@model.h), pi@model.T, float(pi@model.d), violation)  # 保留割系数、对偶向量及强对偶校验所需的源违反量。
    return violation, cut  # 返回 ν 和割；只有不可行候选的割才由外层加入主问题。


def solve_by_cuts(model: ElectricalModel, mode: Literal["min_cost", "max_lambda"], query: float,  # （BENDERS）保留二元变量的整数主问题 + 连续电气 LP 可行性切割。
                  cuts: list[DualCut] | None = None) -> SolveResult:  # cuts 可由调用者共享；query 的单位由 mode 决定。
    """保留线型二元变量的精确 MILP 主问题，不对完整 MIP 读取影子价。

    min_cost: query 是固定倍率；max_lambda: query 是固定建设预算(CNY)。
    共享 cuts 仅包含统一模型下的全局有效割，因此可跨预算和倍率复用。
    """
    pool = [] if cuts is None else cuts  # 共享割池 Π_t；若未传入则使用新的空列表。
    initial_cut_count = len(pool); nx = model.nx; trace = []  # 记录初始割数、二元维数和迭代日志。
    trace_columns = (  # 保留旧列，并记录每次真实主问题/子问题调用，供 notebook 逐轮回放。
        "iteration", "cut_count_before", "lower_bound", "upper_bound", "violation",
        "candidate_cost_cny", "candidate_lambda", "design_capacity", "design", "candidate_x",
        "master_bound", "absolute_gap", "cut_added", "cut_id", "cut_count_after",
        "cut_margin_at_candidate", "dual_stationarity_error", "dual_normalization",
        "dual_min_multiplier", "strong_duality_error", "action", "master_status",
    )
    best_lower = 0.0; best_upper = math.inf if mode == "min_cost" else model.config.lambda_search_max  # 初始化目标下界 L=0；费用上界 +∞ 或倍率上界 3。
    master, choices, multiplier = build_master(model, mode, query, pool)  # 建立包含已知割的主问题；choices=x，multiplier=λ。
    while True:  # 重复求主问题与电气 LP，直到可行候选或主问题不可行。
        master.optimize()  # 精确求解当前 MILP 外近似。
        if master.Status == GRB.INFEASIBLE:  # 若外近似已不可行，则原规划问题也不可行（割均有效）。
            trace.append({"iteration": len(trace)+1, "cut_count_before": len(pool),
                          "lower_bound": math.inf if mode == "min_cost" else -math.inf,
                          "upper_bound": math.inf if mode == "min_cost" else -math.inf,
                          "cut_added": False, "cut_count_after": len(pool),
                          "action": "master_infeasible", "master_status": "infeasible"})
            master.dispose()  # 释放当前主问题求解器资源。
            return SolveResult(mode, query, "infeasible", math.inf if mode == "min_cost" else -math.inf,  # 不可行时目标用 +∞（费用）或 −∞（倍率）表示。
                               np.zeros(nx), math.inf, math.nan, pd.DataFrame(trace, columns=trace_columns), len(pool)-initial_cut_count)  # 其余数值为占位值；返回轨迹和新增割数，不代表可行零建设方案。
        if master.Status != GRB.OPTIMAL:  # 不把限时、中断等状态误当作已证明最优。
            status = master.Status
            master.dispose()
            raise RuntimeError(f"主问题未证明最优，Gurobi status={status}")
        x = np.rint([variable.X for variable in choices.values()])  # 读取主问题二元解 x̄ 并舍入整数容差。
        alpha = multiplier.X  # 读取主问题候选倍率 λ̄。
        design = evaluate_design(model, x)  # 独立解析该合法生成树的 K(x̄)、Λ(x̄)，用于校验及倍率可行下界。
        violation, cut = feasibility_oracle(model, x, alpha)  # 求 ν(x̄,λ̄) 和当前 LP 对偶割；π 来自此 LP 而非整数主问题。
        if mode == "min_cost":  # 固定倍率最小化费用时，放宽电气约束的主问题给费用下界。
            best_lower = max(best_lower, master.ObjBound * 1000)  # L←max(L,1000×ObjBound)；将千元目标界换回 CNY。
            if violation <= model.config.feasibility_tolerance:  # 若 ν≤ε_feas，则候选满足原始电气模型至接受容差。
                best_upper = min(best_upper, design.cost_cny)  # U←min(U,K(x̄))；可实施方案给费用上界。
        else:  # 固定预算最大化倍率时，主问题外近似给倍率上界。
            best_upper = min(best_upper, master.ObjBound)  # U←min(U,ObjBound)。
            candidate_lower = min(design.lambda_max, model.config.lambda_search_max)  # 该候选树在 min(Λ(x̄),λ_search) 处可实施，且建设费已满足固定预算。
            best_lower = max(best_lower, candidate_lower)  # L←max(L,该树可实现倍率)；即使主问题候选 λ̄ 过大，树仍可给下界。
        accepted = violation <= model.config.feasibility_tolerance
        margin = cut.constant + float(cut.x_coeff@x) + cut.lambda_coeff*alpha
        trace.append({
            "iteration": len(trace)+1, "cut_count_before": len(pool),
            "lower_bound": best_lower, "upper_bound": best_upper, "violation": violation,
            "candidate_cost_cny": design.cost_cny, "candidate_lambda": alpha,
            "design_capacity": min(design.lambda_max, model.config.lambda_search_max),
            "design": design.description, "candidate_x": json.dumps(x.astype(int).tolist()),
            "master_bound": master.ObjBound * (1000 if mode == "min_cost" else 1),
            "absolute_gap": best_upper-best_lower,
            "cut_added": not accepted, "cut_id": len(pool)+1 if not accepted else None,
            "cut_count_after": len(pool)+(not accepted),
            "cut_margin_at_candidate": margin,  # 生成点处应为 -ν；割不是连接前沿点的直线。
            "dual_stationarity_error": float(np.max(np.abs(model.W.T@cut.pi))),
            "dual_normalization": float(model.scales@cut.pi),
            "dual_min_multiplier": float(np.min(cut.pi)),
            "strong_duality_error": abs(margin+violation),
            "action": "accept" if accepted else "add_cut", "master_status": "optimal",
        })
        if violation <= model.config.feasibility_tolerance:  # 接受首个电气可行的主问题最优候选；此时外近似界与可行目标闭合至容差。
            value = design.cost_cny if mode == "min_cost" else alpha  # 返回目标：O1 为费用 CNY，O2 为 λ̄。
            master.dispose()  # 释放主问题资源。
            return SolveResult(mode, query, "optimal", value, x, design.cost_cny, alpha,  # 打包最优目标、设计和倍率。
                               pd.DataFrame(trace, columns=trace_columns), len(pool)-initial_cut_count)  # 附带收敛轨迹及本次新增割数。
        pool.append(cut)  # 不可行时保留新割；其在生成点取值 −ν<0。
        add_dual_cut(master, choices, multiplier, cut)  # 把新割（CUT）加入同一个主问题，然后重新求解。


def validate_cuts(cuts: Sequence[DualCut], designs: Sequence[Design]) -> float:  # （CUT-CHECK）以独立设计检查每条割在全部真实可行区间上不误删点。
    """每份割在全部整数方案的倍率区间两端检查；仿射性保证区间内部同样有效。"""
    worst = math.inf  # 初始化全部割的最小有效性余量。
    for cut in cuts:  # 逐条检查割 β₀+β_xᵀx+β_λλ≥0。
        margin = min(cut.constant+float(cut.x_coeff@design.x)+min(0.0, cut.lambda_coeff*design.lambda_max)  # 对固定 x，λ∈[0,Λ(x)] 上仿射函数的最小值为 β₀+β_xᵀx+min(0,β_λΛ(x))。
                     for design in designs)  # 再对全部 216 个设计取最小值；端点验证覆盖整个区间。
        worst = min(worst, margin)  # 摘要只保留所有割的最小有效性余量。
        assert margin >= -1e-6, f"对偶割误删真实可行点：margin={margin}"  # 断言 μ≥−10⁻⁶；小负数允许浮点舍入。
    return worst if cuts else 0.0  # 返回最差余量；没有割时约定返回 0。


def discover_frontier_by_cuts(model: ElectricalModel, cuts: list[DualCut], *,
                              query_history: list[SolveResult] | None = None,
                              round_history: list[dict] | None = None) -> pd.DataFrame:  # （EPS）预算 ε-constraint 递推。
    """不用方案穷举和预设倍率网格，以预算 ε-constraint 递推恢复本例完整 Pareto 阶梯。

    本合成例子的所有费用为整数 CNY，gcd 给出精确费用粒度。
    每轮先求当前预算的最大倍率，再求实现该倍率的最低费用；
    下一轮预算设为该费用减去一个费用粒度，避免任何人工倍率步长。
    可选 history 列表只记录过程，不参与优化；原有 DataFrame 返回接口不变。
    """
    rounded = np.rint(model.cost).astype(int)  # 费用四舍五入成整数；精确费用粒度论证依赖默认整数 CNY 参数。
    if not np.allclose(model.cost, rounded, atol=1e-8, rtol=0) or np.any(rounded < 0):
        raise ValueError("完整前沿递推要求非负整数费用；小数费用请先统一换成精确的最小货币单位。")
    positive = [int(c) for c in rounded if c > 0]  # 取严格正费用项，排除已有 L 型的零费用。
    quantum = math.gcd(*positive) if positive else 1  # g=20 CNY；全零费用时一步即可完成。
    budget = float(sum(positive))  # B₀=Σ_{c_ek>0}c_ek，是覆盖所有设计费用的保守初始预算。
    rows = []  # 初始化按预算递减发现的前沿记录。
    while budget >= 0:  # 只要预算非负就继续向更低费用区域搜索。
        cuts_before = len(cuts)
        capacity = solve_by_cuts(model, "max_lambda", budget, cuts)  # 先求 Λ*(B_t)=max λ，满足 K(x)≤B_t；与先前查询共享割池。
        if query_history is not None:
            query_history.append(capacity)
        if capacity.status == "infeasible":  # 若该预算已没有可行设计，则更低预算也不可能可行。
            break  # 结束预算递推。
        # 取经过独立电气验算的候选倍率，避免 LP/MILP 边界浮点误差。
        achieved = min(capacity.lambda_value, evaluate_design(model, capacity.x).lambda_max)  # λ_t 取 MILP 返回倍率与候选树解析 Λ(x) 的较小值，避免边界浮点超限。
        minimum = solve_by_cuts(model, "min_cost", achieved, cuts)  # 再求 K_t=K*(λ_t)，找到实现该倍率的最低建设费。
        if query_history is not None:
            query_history.append(minimum)
        if minimum.status != "optimal":
            raise RuntimeError("MP2 给出的可行倍率未通过 MP1，请检查数值容差和模型一致性。")
        design = evaluate_design(model, minimum.x)  # 读取该最低费设计的瓶颈及描述。
        rows.append({"cost_cny": minimum.value, "lambda_max": achieved, "bottleneck": design.bottleneck,  # 保存前沿点 (K_t,λ_t) 和瓶颈。
                     "design": design.description})  # 每个台阶仅保留费用、倍率、瓶颈和建设方案。
        next_budget = float(round(minimum.value)-quantum)  # B_{t+1}=K_t−g；不是 B_t−g。
        if next_budget >= budget:
            raise RuntimeError("预算递推未严格下降，请检查费用粒度和求解精度。")
        if round_history is not None:
            round_history.append({
                "outer_round": len(rows), "budget_cny": budget,
                "mp2_lambda": achieved, "mp2_design_cost_cny": capacity.cost_cny,
                "mp1_cost_cny": minimum.value, "next_budget_cny": next_budget,
                "quantum_cny": quantum, "mp2_iterations": len(capacity.trace),
                "mp1_iterations": len(minimum.trace), "mp2_new_cuts": capacity.new_cuts,
                "mp1_new_cuts": minimum.new_cuts, "cuts_before": cuts_before,
                "cuts_after": len(cuts), "design": design.description,
            })
        budget = next_budget
    return pd.DataFrame(rows).sort_values("cost_cny").reset_index(drop=True)  # 按费用升序输出发现的阶梯边界。


def convex_hull_cost(frontier: pd.DataFrame, lambdas: Sequence[float]) -> np.ndarray:  # （HULL）计算真实预算域凸包下边界的辅助 LP，完整式见 README §5。
    """完整预算域凸包的下边界，不是整数可行域，也不等同于特定 formulation 的 LP 松弛。"""
    lam = frontier.lambda_max.to_numpy(); cost = frontier.cost_cny.to_numpy(); out = []  # 提取前沿点 (Λ_j,K_j)；out 保存各目标倍率的凸包最低费用。
    lp = gp.Model("convex_hull")  # 创建独立连续凸组合 LP。
    lp.Params.OutputFlag = 0  # 关闭日志。
    weights = lp.addVars(len(lam), lb=0.0, name="weight")  # θ_j≥0：各前沿设计的数学混合权重；不是可施工的二元选线变量。
    lp.addConstr(weights.sum() == 1.0, name="convex_combination")  # （H1）Σ_j θ_j=1，凸组合权重归一化。
    capacity = lp.addConstr(gp.quicksum(lam[i] * weights[i] for i in weights) >= 0.0, name="capacity")  # （H2）Σ_j Λ_jθ_j≥λ；先设 RHS=0，查询时替换。
    lp.setObjective(gp.quicksum(cost[i] * weights[i] for i in weights), GRB.MINIMIZE)  # （O4）min Σ_j K_jθ_j，费用单位 CNY。
    for alpha in lambdas:  # 逐个查询目标倍率 λ。
        capacity.RHS = float(alpha)  # 把 H2 的右端更新为当前倍率。
        lp.optimize()  # 求当前凸包最低费 LP。
        out.append(math.inf if lp.Status == GRB.INFEASIBLE else lp.ObjVal)  # 不可行返回 +∞，否则保存凸包费用目标值。
    lp.dispose()  # 释放辅助 LP。
    return np.asarray(out)  # 返回与输入倍率顺序一致的费用数组。


def run_experiment(output_dir: str | Path, config: DemoConfig = DemoConfig(),  # 运行完整验证实验；此函数组织查询和文件输出，不新增物理模型。
                   lines: tuple[LineType, ...] = DEFAULT_LINES,  # 所有线型参数传到 build_model，避免验证与查询使用不同数据。
                   corridors: tuple[Corridor, ...] = DEFAULT_CORRIDORS) -> dict:  # 所有走廊参数同样传递；数值展开和 216 方案结论针对默认参数。
    """完整执行独立校验，仅保存最优前沿、冷启动收敛轨迹和核心摘要。"""
    output = Path(output_dir); output.mkdir(parents=True, exist_ok=True)  # 创建结果目录；仅文件操作。
    model = build_model(config, lines, corridors)  # 建立统一电气模型及可复用子问题。
    designs = enumerate_designs(model); frontier = frontier_table(designs)  # 独立枚举设计并提取精确前沿，作为之后断言的参照。
    # 零费用可行性与每个方案最大倍率均用 LP 交叉验证，并检查稍微超限后违反量为正。
    for design in designs:  # 逐个检查默认 216 个树与线型设计。
        violation, boundary_cut = feasibility_oracle(model, design.x, design.lambda_max)  # 在解析最大倍率 Λ(x) 处求 LP；理论上 ν=0。
        assert violation <= 1e-7  # 断言边界点违反量≤10⁻⁷。
        over_violation, over_cut = feasibility_oracle(model, design.x, design.lambda_max+1e-4)  # 在 Λ(x)+10⁻⁴ 处求 LP；该点应超过至少一个物理上界。
        assert over_violation > 1e-9  # 断言超限点的违反量>10⁻⁹。
        for alpha, cut in ((design.lambda_max, boundary_cut), (design.lambda_max+1e-4, over_cut)):  # 对边界和超限两个 LP 都核对对偶可行性及强对偶。
            rhs = model.h + model.T@design.x + model.d*alpha  # 重建对应查询的 b=h+Tx+dλ。
            assert np.min(cut.pi) >= -1e-8  # （DUAL）断言 π≥0，容差 10⁻⁸。
            assert np.max(np.abs(model.W.T@cut.pi)) < 1e-7  # （DUAL）断言 Wᵀπ=0，容差 10⁻⁷；因为 8 个状态变量都自由。
            assert float(cut.pi@model.scales) <= 1+1e-7  # （DUAL）断言 rᵀπ≤1，容差 10⁻⁷；来自 η≥0 和目标系数 1。
            assert abs(-float(cut.pi@rhs)-cut.source_violation) < 1e-7  # （STRONG）断言 −πᵀb=ν，容差 10⁻⁷。
    cuts: list[DualCut] = []  # 冷启动及后续查询共用割池；不使用枚举结果构造割。
    cold = solve_by_cuts(model, "max_lambda", 33000.0, cuts)  # 执行 B=33000 CNY 的 max_lambda 查询。
    validations = []  # 查询表仅供交互查看，不重复导出枚举真值和求解统计。
    for alpha in (0.4, 0.7, 1.0, 1.3, 1.6, 1.8, 2.1, 2.4):  # 检查八个固定倍率，含超过 λ_tr 的 2.4。
        result = solve_by_cuts(model, "min_cost", alpha, cuts)  # 用切割算法求该倍率下的最低建设费 K*(λ)。
        feasible = [z.cost_cny for z in designs if z.lambda_max >= alpha-1e-8]  # 独立筛出 Λ(x)≥λ 的设计费用，边界筛选容差 10⁻⁸。
        exact = min(feasible) if feasible else math.inf  # 枚举真值=min 可行设计费用；无方案则 +∞。
        assert (math.isinf(exact) and result.status == "infeasible") or abs(result.value-exact) < 1e-4  # 断言切割法与枚举的可行性/最低费用一致。
        validations.append({"mode": "min_cost", "query": alpha,  # 保留查询条件、最优值和可行状态。
                            "value": result.value, "status": result.status})
    for budget in (0.0, 10000.0, 15000.0, 33000.0, 42000.0, 52000.0, 57000.0):  # 检查七个固定预算，包括零预算与高预算。
        result = solve_by_cuts(model, "max_lambda", budget, cuts)  # 用切割算法求 Λ*(B)。
        exact = max(z.lambda_max for z in designs if z.cost_cny <= budget+1e-8)  # 枚举真值=max_{K(x)≤B}Λ(x)；默认 λ_tr<λ_search，无需截断。
        assert abs(result.value-exact) < 2e-6  # 断言最大倍率与枚举之差<2×10⁻⁶。
        validations.append({"mode": "max_lambda", "query": budget,  # 保留查询条件、最优值和可行状态。
                            "value": result.value, "status": result.status})
    discovered = discover_frontier_by_cuts(model, cuts)  # 用预算递推和切割独立发现完整前沿，不传入枚举 designs。
    assert len(discovered) == len(frontier)  # 断言发现的前沿台阶数与独立枚举一致。
    assert np.allclose(discovered.cost_cny, frontier.cost_cny, atol=1e-6, rtol=0)  # 断言各拐点建设费相同，绝对容差 10⁻⁶ CNY。
    assert np.allclose(discovered.lambda_max, frontier.lambda_max, atol=2e-6, rtol=0)  # 断言各拐点倍率相同，绝对容差 2×10⁻⁶。
    worst = validate_cuts(cuts, designs)  # （CUT-CHECK）验证所有累计割在所有设计可行区间上有效。
    target = 1.1  # 选 λ=1.1 作为整数域与凸包的比较点。
    exact_cost = min(z.cost_cny for z in designs if z.lambda_max >= target)  # 求此倍率的整数最低费用 K*(1.1)。
    hull_cost = float(convex_hull_cost(frontier, [target])[0])  # 求同一倍率的凸包最低费用 K_hull(1.1)。
    summary = {"enumerated_designs": len(designs), "frontier_steps": len(frontier),  # 只保留关键结果、校验证据与适用模型。
               "validation_queries": len(validations), "total_cuts": len(cuts), "worst_cut_margin": worst,  # 已完成的查询校验数及全局割有效性。
               "cold_budget_cny": cold.query, "cold_lambda_max": cold.value,  # 代表性预算对应的最大可行倍率。
               "nonconvex_counterexample_lambda": target,  # 整数域与凸包的比较倍率。
               "exact_cost_at_counterexample": exact_cost, "convex_hull_lower_cost_at_counterexample": hull_cost,  # 记录反例处整数最低费和凸包最低费。
               "ac_power_flow_validated": False, "electrical_model": "lossless fixed-PF linear squared-voltage model"}  # 明确本实验未做 AC 潮流复核；电气结论属于固定 PF 无损线性平方电压模型。
    frontier.to_csv(output/"exact_frontier.csv", index=False)  # 两种方法核验一致后，只导出一份最优前沿。
    cold.trace.to_csv(output/"cold_start_cut_trace.csv", index=False)  # 仅导出一次代表性的收敛轨迹。
    (output/"summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")  # 保存严格 JSON 摘要；禁止 NaN 以便其他程序读取。
    return {"model": model, "designs": designs, "frontier": frontier, "cold_result": cold, "cuts": cuts,  # 返回模型、枚举、前沿、冷启动结果和割池，供 Notebook/绘图使用。
            "validation": pd.DataFrame(validations), "discovered_frontier": discovered, "summary": summary}  # 同时返回校验表、算法发现前沿和摘要。


if __name__ == "__main__":  # 仅直接运行该文件时执行实验；import 时只定义模型及函数。
    import argparse  # 命令行解析模块；无对应优化变量。
    parser = argparse.ArgumentParser(description=__doc__)  # 创建命令行解析器，使用模块说明作为帮助文本。
    parser.add_argument("--output", default="results", help="本地结果目录")  # --output 指定结果目录，默认 results。
    args = parser.parse_args()  # 解析输出目录参数。
    result = run_experiment(args.output)  # 运行默认模型的完整实验和断言。
    summary = result["summary"]  # 控制台只显示完成状态、代表性结果和保存位置。
    print(f"校验通过：{summary['enumerated_designs']} 个方案，{summary['validation_queries']} 个查询，{summary['frontier_steps']} 级最优前沿。")
    print(f"预算 {summary['cold_budget_cny']:,.0f} CNY：最大负荷倍率 {summary['cold_lambda_max']:.6f}。")
    print(f"结果目录：{Path(args.output).resolve()}")
    result["model"].subproblem.dispose()  # 释放返回模型中的电气 LP 资源。
