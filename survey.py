"""Case33-S 上的走廊估值与勘察策略（创新点二）：以 RCUT 可行域族为价值函数，比较勘察策略。

候选路 G={C1..C5}（Case33S.candidates 的顺序），子集记为 5 位掩码。信息状态 s=(A,B)：A 确认可用、B 确认阻断，
U=G∖(A∪B) 未知。设计集 X(A,b) 为只用 A 中候选、建设费 C(x)<=b 的径向网架；R(A,b)=∪_{x∈X(A,b)} R_x 取 RCUT 内域，
端口 z=(p18,p25)，只算 ++ 分区。X(A,b) 只依赖可选候选子集族 {S⊆A: C(S)<=b}：32 个 A × 12 档预算共 101 个族。
κ(A,b)=|R(A,b)∩P|/|P|，P=[0,z̄]，z̄ 为 R(G,14) 的坐标最大值；C*_A(z) 为服务 z 的最小建设费，
Φ(A)=E_μ[min(C*_A,B̄)]=Σ_k b_k[κ(A,b_k)-κ(A,b_{k-1})]+B̄[1-κ(A,b_K)]，预算档为候选子集费用的全部不同值。
勘察：共享 Beta(α0,β0) 先验，预测概率 q_s=(α0+|A|)/(α0+β0+|A|+|B|)，勘察费 c^sur=ρ·c^c；J(s) 为最优期望总成本，
本文策略为束指标 σ_j=Δ_j-C_j/Q_j，对照为 DP、单路比值、不勘察与全知；期望对 2^5 种可用性真值精确求和。
记号见 docs/notation.md「勘察（Case33-S）」。正式输出在 results/survey/case33/，过程信息只打印到终端。
"""
import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
from functools import reduce
import json
from math import factorial
import multiprocessing
from operator import or_
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.special import betaln
from shapely import contains_xy
from shapely.geometry import Point, box
from threadpoolctl import threadpool_limits

from main import region_settings
from model import LOAD_PF
from monitor import RunMonitor
from Network.case33bw import Case33S
from plot import draw_survey_layout, draw_survey_policy, draw_survey_storyboard, draw_survey_valuation
from region import Radial, build_partition, polygon_union
from vertify import budget_schemes

OUTPUT = Path(__file__).resolve().parent/'results'/'survey'/'case33'
ROADS = tuple(f'C{k+1}' for k in range(len(Case33S.candidates)))   # 候选路名；子集掩码的第 k 位为 ROADS[k]
COST = np.array(list(Case33S.candidates.values()))                 # 相对建设费 c^c（1 单位为 C4 造价）
FULL = (1 << len(ROADS))-1
SUBSETS = range(FULL+1)
STATES = [(A, B) for A in SUBSETS for B in SUBSETS if not A & B]    # 信息状态 (A,B)，共 3^5=243 个
SIGN = (1, 1)                       # 只算 ++ 分区
PRIOR = (3., 2.)                    # 共享 Beta(α0,β0) 先验，q0=0.6
RHO = .1                            # 勘察费 c^sur=ρ·c^c
FALLBACK = 20.                      # 兜底费 B̄：不可服务点的成本（相对单位，年化系数 1）
GRID_Q0 = (.3, .5, .7, .9)          # 敏感性网格：先验均值 q0（α0+β0=PRIOR_WEIGHT）与勘察费比例 ρ
GRID_RHO = (.05, .1, .25, .5)
PRIOR_WEIGHT = 5.
EXAMPLE_TRUTH = 0b11001             # 示例轨迹的真值：C2、C3 阻断，其余可用
FAMILY_SECONDS = 300.               # 每个族的 RCUT 时限
WORKERS = 16                        # 族与抽查的并行进程数，每进程单线程
SPOT_RAYS = 128                     # 抽查：每个网架的射线数
SPOT_FAMILIES = (0, FULL, 0b11000)  # 抽查的族（b=14）：无候选、全部候选、{C4,C5}
RASTER = 400                        # 故事板 C*(z) 栅格的每轴点数
TIE = 1e-9                          # DP 并列最优的判定容差
POLICIES = ('No survey', 'Single-road ratio', 'Bundle index (ours)', 'DP optimal', 'Clairvoyant')


def members(S):
    """掩码 → 候选路序号。"""
    return [k for k in range(len(ROADS)) if S >> k & 1]


def subset(names):
    """候选路名 → 掩码。"""
    return sum(1 << ROADS.index(name) for name in names)


def label(S):
    return '{'+','.join(ROADS[k] for k in members(S))+'}' if S else '∅'


def cost(S):
    """候选子集的建设费 C(S)。"""
    return float(COST[members(S)].sum())


LEVELS = sorted({cost(S) for S in SUBSETS})   # 预算档 b_1..b_K：候选子集费用的全部不同值


def family(A, b):
    """(A,b) 的可选候选子集族 {S⊆A: C(S)<=b} 的规范代表 (A',b')：A' 为族中子集的并，b' 为族中最大的费用。"""
    chosen = [S for S in SUBSETS if not S & ~A and cost(S) <= b]
    return reduce(or_, chosen, 0), max(map(cost, chosen))


def built(network, x):
    """网架 x 所建候选路的掩码。"""
    active = {key[0] for key, bit in zip(network.type_keys, x) if bit}
    return sum(1 << k for k, road in enumerate(network.roads) if road in active)


# ---- 可行域族 -----------------------------------------------------------------------------------
def run_family(key):
    """一个族 (A',b') 的 RCUT 内域（++ 分区）：Case33S 只开放 A' 中的候选路，预算 b'，单线程。返回状态、终态间隙、
    RCUT 中出现过的网架（X* 与割过的网架，build_partition:schemes）的候选集与内域块（kW）。"""
    A, b = key
    network = Case33S(available=members(A))
    monitor, settings = RunMonitor(), region_settings(len(SIGN))
    monitor.time_limit, settings.workers = FAMILY_SECONDS, 1
    with threadpool_limits(limits=1):
        result = build_partition(network, SIGN, b, monitor, settings)
    return dict(roads=[ROADS[k] for k in members(A)], budget=b, status=result['status'], gap=result['gap'],
                seconds=result['seconds'], cones=result['cones'],
                networks=[[ROADS[k] for k in members(S)] for S in sorted({built(network, x) for x in result['schemes']})],
                inner=[np.round(vertices, 3).tolist() for vertices in result['inner']])


def compute_families(path):
    """全部族并行跑 RCUT（(G,14) 最先提交，它定出 P），结果写入紧凑 json，供画图与复算。"""
    keys = sorted({family(A, b) for A in SUBSETS for b in LEVELS}, key=lambda key: (-key[0].bit_count(), -key[1]))
    rows = []
    with ProcessPoolExecutor(WORKERS, mp_context=multiprocessing.get_context('spawn')) as pool:
        for row in pool.map(run_family, keys):
            print(f"族 {label(subset(row['roads']))} b={row['budget']:g}：{row['status']}，间隙 {row['gap']:.3%}，"
                  f"{row['seconds']:.1f} s，{row['cones']} 个锥，{len(row['networks'])} 种候选集", flush=True)
            rows.append(row)
    path.write_bytes(json.dumps(rows, ensure_ascii=False, separators=(',', ':')).encode('utf-8'))


def load_families(path):
    """族表 {(A',b'): dict(status, gap, seconds, cones, networks=候选集掩码, inner=内域块)}。"""
    return {(subset(row['roads']), row['budget']): dict(row, networks=[subset(names) for names in row['networks']])
            for row in json.loads(path.read_text(encoding='utf-8'))}


def coverage(families):
    """P=[0,z̄1]×[0,z̄2]（z̄ 为 R(G,14) 的坐标最大值）、各族内域的并与 κ(A,b_k) 表（32×12）。"""
    unions = {key: polygon_union(row['inner']) for key, row in families.items()}
    P = box(0., 0., *unions[family(FULL, LEVELS[-1])].bounds[2:])
    kappa = np.array([[unions[family(A, b)].intersection(P).area for b in LEVELS] for A in SUBSETS])/P.area
    return P, unions, kappa


def value(kappa):
    """Φ(A)=Σ_k b_k[κ(A,b_k)-κ(A,b_{k-1})]+B̄[1-κ(A,b_K)]；b_1=0，κ(A,b_0) 取 0 不影响结果。"""
    return np.diff(kappa, axis=1, prepend=0.)@np.array(LEVELS)+FALLBACK*(1.-kappa[:, -1])


# ---- 勘察决策 -----------------------------------------------------------------------------------
class Survey:
    """一组勘察参数下的 DP 与策略。prior=(α0,β0)，勘察费 c^sur=ρ·c^c；phi、kappa 为 Φ(A) 与 κ(A,b_k) 表，
    networks[O] 为乐观族 R(O,14) 的 RCUT 中出现过的网架的候选集（掩码）。"""

    def __init__(self, phi, kappa, networks, prior=PRIOR, rho=RHO):
        self.phi, self.kappa, self.networks, self.prior = phi, kappa, networks, prior
        self.survey = rho*COST
        # 1. 真值 χ 的概率（beta-binomial）：P(χ)=B(α0+k,β0+5-k)/B(α0,β0)，k 为可用条数
        n = len(ROADS)
        self.weights = np.exp([betaln(prior[0]+S.bit_count(), prior[1]+n-S.bit_count())-betaln(*prior) for S in SUBSETS])
        # 2. DP：按已知路数从多到少递推 J(s)=min{Φ(A), min_g [c^sur_g+q_s J(A∪g,B)+(1-q_s) J(A,B∪g)]}
        self.J, self.best = {}, {}
        for A, B in sorted(STATES, key=lambda state: -(state[0] | state[1]).bit_count()):
            q = self.q(A, B)
            values = {None: phi[A], **{g: self.survey[g]+q*self.J[A | 1 << g, B]+(1-q)*self.J[A, B | 1 << g]
                                       for g in members(FULL & ~(A | B))}}
            self.J[A, B] = min(values.values())
            self.best[A, B] = {action for action, v in values.items() if v <= self.J[A, B]+TIE}
        # 3. 单路比值的勘察次序：s0 时按 c^c_e/Δκ_e 升序，Δκ_e=κ({e},14)-κ(∅,14)<=0 的不排
        gain = kappa[[1 << g for g in range(n)], -1]-kappa[0, -1]
        self.order = sorted((g for g in range(n) if gain[g] > 0.), key=lambda g: COST[g]/gain[g])

    def q(self, A, B):
        """预测概率 q_s=(α0+|A|)/(α0+β0+|A|+|B|)。"""
        return (self.prior[0]+A.bit_count())/(sum(self.prior)+(A | B).bit_count())

    def bundles(self, A, B):
        """束 G_j=cand(x)∖A（x 取乐观族 R(A∪U,14) 的 RCUT 中出现过的网架，去重去空），按 σ_j 从大到小：
        dict(roads=束内按勘察费升序的路, delta=Δ_j, Q=Q_j, C=C_j, sigma=σ_j)。"""
        alpha, beta = self.prior[0]+A.bit_count(), self.prior[1]+B.bit_count()
        rows = []
        for G in sorted({S & ~A for S in self.networks[FULL & ~B]}-{0}):
            order = sorted(members(G), key=lambda g: (self.survey[g], g))
            reach = np.cumprod([1.]+[(alpha+i)/(alpha+beta+i) for i in range(len(order))])   # 前 k 条都可用的概率
            delta, C = self.phi[A]-self.phi[A | G], float(self.survey[order]@reach[:-1])
            rows.append(dict(roads=order, delta=delta, Q=reach[-1], C=C, sigma=delta-C/reach[-1]))
        return sorted(rows, key=lambda row: -row['sigma'])

    def certified(self, A, B):
        """停止证书：Φ(A)-Φ(A∪U)<=min_{g∈U} c^sur_g（U 为空时成立）。"""
        unknown = FULL & ~(A | B)
        return not unknown or self.phi[A]-self.phi[A | unknown] <= self.survey[members(unknown)].min()

    def bundle(self, A, B):
        """本文策略：停止证书成立或 max σ_j<=0 时停止，否则勘察 σ 最大的束中勘察费最低的路。"""
        if self.certified(A, B):
            return None
        rows = self.bundles(A, B)
        return rows[0]['roads'][0] if rows and rows[0]['sigma'] > 0. else None

    def ratio(self, A, B):
        """单路比值：按 s0 的 c^c/Δκ 次序依次勘察，下一条的 q_s[Φ(A)-Φ(A∪e)]-c^sur_e<=0 时停止。"""
        rest = [g for g in self.order if not (A | B) >> g & 1]
        return rest[0] if rest and self.q(A, B)*(self.phi[A]-self.phi[A | 1 << rest[0]]) > self.survey[rest[0]] else None

    def optimal(self, A, B):
        """DP 最优动作；并列时停止优先，其次序号小者。"""
        return min(self.best[A, B], key=lambda action: -1 if action is None else action)

    def evaluate(self, policy):
        """策略在 2^5 种真值下的指标；DP 动作一致率按访问概率加权（每次勘察与最后的停止各算一次决策）。"""
        totals, decisions = np.zeros(5), 0.   # C^tot、勘察次数、勘察费、κ(A_final,14)、与 DP 一致的决策数
        for truth, weight in zip(SUBSETS, self.weights):
            A = B = 0
            spent = surveyed = agreed = steps = 0
            while True:
                action = policy(A, B)
                agreed, steps = agreed+(action in self.best[A, B]), steps+1
                if action is None:
                    break
                spent, surveyed = spent+self.survey[action], surveyed+1
                A, B = (A | 1 << action, B) if truth >> action & 1 else (A, B | 1 << action)
            totals += weight*np.array([spent+self.phi[A], surveyed, spent, self.kappa[A, -1], agreed])
            decisions += weight*steps
        return self.metrics(*totals[:4], agreement=totals[4]/decisions)

    def metrics(self, total, surveys, spent, covered, agreement=None):
        """E[C^tot]、gap=(E[C^tot]-J(s0))/J(s0)、G=Φ(∅)-E[C^tot]、η=G/VPI、E[勘察次数]、E[勘察费]、覆盖增量、DP 动作一致率。"""
        J0, base = self.J[0, 0], self.phi[0]
        vpi = base-self.weights@self.phi
        return dict(expected_cost=total, gap=(total-J0)/J0, gain=base-total, capture=(base-total)/vpi,
                    surveys=surveys, survey_cost=spent, coverage_gain=covered-self.kappa[0, -1], dp_agreement=agreement)

    def comparison(self):
        """表 3：五种做法；全知为免费获知全部可用性，E_χ[Φ(A(χ))]。"""
        return dict(zip(POLICIES, (self.evaluate(lambda A, B: None), self.evaluate(self.ratio), self.evaluate(self.bundle),
                                   self.evaluate(self.optimal),
                                   self.metrics(self.weights@self.phi, 0., 0., self.weights@self.kappa[:, -1]))))


def valuation(survey):
    """表 1（s0 的走廊估值，V(S)=Φ(∅)-Φ(S)）与两两交互 I_ij=V({i,j})-V({i})-V({j})。"""
    V, n = survey.phi[0]-survey.phi, len(ROADS)
    first = survey.bundles(0, 0)
    rows = []
    for e in range(n):
        shapley = sum(factorial(S.bit_count())*factorial(n-1-S.bit_count())/factorial(n)*(V[S | 1 << e]-V[S])
                      for S in SUBSETS if not S >> e & 1)
        best = next((row for row in first if e in row['roads']), None)   # first 按 σ 降序
        rows.append(dict(road=ROADS[e], c_c=COST[e], c_sur=survey.survey[e], standalone=V[1 << e],
                         leave_one_out=V[FULL]-V[FULL & ~(1 << e)], shapley=shapley,
                         ratio_rank=survey.order.index(e)+1 if e in survey.order else None,
                         best_bundle=best and label(sum(1 << g for g in best['roads'])), best_sigma=best and best['sigma']))
    interaction = np.array([[V[1 << i | 1 << j]-V[1 << i]-V[1 << j] if i != j else np.nan for j in range(n)]
                            for i in range(n)])
    return rows, interaction


def trajectory(survey, truth):
    """本文策略在真值 truth 下的轨迹：各状态的前三个束、所测道路与结果、两个 κ、Φ 上下界与到达时的累计勘察费。"""
    A = B = 0
    spent, steps = 0., []
    while True:
        unknown, action = FULL & ~(A | B), survey.bundle(A, B)
        steps.append(dict(A=A, B=B, bundles=survey.bundles(A, B)[:3], action=action,
                          available=None if action is None else bool(truth >> action & 1),
                          kappa=survey.kappa[A, -1], kappa_optimistic=survey.kappa[A | unknown, -1],
                          lower=survey.phi[A | unknown], upper=survey.phi[A], spent=spent))
        if action is None:
            return steps
        spent += survey.survey[action]
        A, B = (A | 1 << action, B) if truth >> action & 1 else (A, B | 1 << action)


def sensitivity(phi, kappa, networks):
    """敏感性网格上本文与单路比值的 gap：[策略, q0, ρ]。"""
    gaps = np.zeros((2, len(GRID_Q0), len(GRID_RHO)))
    for i, q0 in enumerate(GRID_Q0):
        for j, rho in enumerate(GRID_RHO):
            survey = Survey(phi, kappa, networks, prior=(q0*PRIOR_WEIGHT, (1.-q0)*PRIOR_WEIGHT), rho=rho)
            gaps[:, i, j] = survey.evaluate(survey.bundle)['gap'], survey.evaluate(survey.ratio)['gap']
    return gaps


# ---- 抽查与栅格 ---------------------------------------------------------------------------------
def network_region(x):
    """固定网架 x 的紧化可行域 R_x（++ 分区，kW）的边界点：SPOT_RAYS 个方向从原点出发的紧化射线远端点，原点不可行时
    另加反向射线的近端点（Radial.vertex / near / origin_feasible，首次使用前 OBBT）；R_x 凸，凸包即内近似。"""
    monitor, settings = RunMonitor(), region_settings(len(SIGN))
    settings.workers = 1
    radial = Radial(Case33S(), SIGN, LEVELS[-1], monitor, settings)
    keys = [radial.direction((np.cos(t), np.sin(t))) for t in np.linspace(0., np.pi/2, SPOT_RAYS)]
    with threadpool_limits(limits=1):
        points = [radial.vertex(x, k) for k in keys]
        points += [np.zeros(2)] if radial.origin_feasible(x) else [radial.near(x, k) for k in keys]
    return [p*radial.bounds for p in points if p is not None]


def spot_check(unions):
    """抽查：87 个网架逐网架求域，按族取并，与 RCUT 内域面积比较；relative_error=RCUT/并-1。"""
    network = Case33S()
    schemes = [tuple(int(v) for v in x) for x in budget_schemes(SimpleNamespace(network=network), LEVELS[-1])]
    with ProcessPoolExecutor(WORKERS, mp_context=multiprocessing.get_context('spawn')) as pool:
        regions = dict(zip(schemes, pool.map(network_region, schemes)))
    rows = []
    for A in SPOT_FAMILIES:
        key = family(A, LEVELS[-1])
        chosen = [x for x in schemes if not built(network, x) & ~key[0] and cost(built(network, x)) <= key[1]]
        area = polygon_union([regions[x] for x in chosen]).area
        rows.append(dict(family=label(key[0]), budget=key[1], networks=len(chosen), rcut_area=unions[key].area,
                         union_area=area, relative_error=unions[key].area/area-1.))
    return rows


def cost_map(unions, A, P):
    """确认集 A 下的最小服务成本 C*_A(z) 的栅格（按预算档离散，不可服务为 nan）。"""
    X, Y = np.meshgrid(np.linspace(0., P.bounds[2], RASTER), np.linspace(0., P.bounds[3], RASTER))
    level = np.full(X.shape, np.nan)
    for b in reversed(LEVELS):
        level[contains_xy(unions[family(A, b)], X, Y)] = b
    return level


# ---- 输出 ---------------------------------------------------------------------------------------
def write_table(path, rows):
    with open(path, 'w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows({k: round(float(v), 6) if isinstance(v, (float, np.floating)) else v for k, v in row.items()}
                         for row in rows)


def bundle_text(row):
    return f"{label(sum(1 << g for g in row['roads']))}: Δ={row['delta']:.4f}, Q={row['Q']:.4f}, C={row['C']:.4f}, σ={row['sigma']:.4f}"


def main(rerun=False):
    OUTPUT.mkdir(parents=True, exist_ok=True)
    # 1. 全部候选、b=14 的径向网架共 87 个
    network = Case33S()
    assert len(budget_schemes(SimpleNamespace(network=network), LEVELS[-1])) == 87
    # 2. 101 个族的 RCUT 内域（缺失或 rerun 时重算）
    path = OUTPUT/'families.json'
    if rerun or not path.exists():
        compute_families(path)
    families = load_families(path)
    # 3. P、κ 与 Φ；现状点属于 R(∅,0)
    P, unions, kappa = coverage(families)
    assert unions[family(0, 0.)].contains(Point(Case33S.status_quo))
    phi = value(kappa)
    networks = {O: families[family(O, LEVELS[-1])]['networks'] for O in SUBSETS}
    print(f"P=[0,{P.bounds[2]:.1f}]×[0,{P.bounds[3]:.1f}] kW；κ 随 b 的最大回落 "
          f"{max(np.max(np.diff(-kappa, axis=1)), 0.):.2e}，随 A 的最大回落 "
          f"{max(kappa[A, k]-kappa[A | 1 << g, k] for A in SUBSETS for g in range(len(ROADS)) for k in range(len(LEVELS))):.2e}")
    # 4. 主设置：DP、五种做法、表 1、示例轨迹；敏感性网格
    survey = Survey(phi, kappa, networks)
    policies = survey.comparison()
    assert np.isclose(policies['DP optimal']['expected_cost'], survey.J[0, 0])
    rows, interaction = valuation(survey)
    assert np.isclose(sum(row['shapley'] for row in rows), phi[0]-phi[FULL])
    steps = trajectory(survey, EXAMPLE_TRUTH)
    gaps = sensitivity(phi, kappa, networks)
    # 5. 抽查：逐网架求域取并，与 RCUT 比较三个族的面积
    spot = spot_check(unions)
    # 6. 表格
    write_table(OUTPUT/'table1.csv', rows)
    write_table(OUTPUT/'table2.csv', [
        dict(step=k, A=label(step['A']), B=label(step['B']),
             **{f'bundle{j+1}': bundle_text(row) if row else '' for j, row in enumerate((step['bundles']+[None]*3)[:3])},
             surveyed='stop' if step['action'] is None else ROADS[step['action']],
             result='' if step['action'] is None else 'available' if step['available'] else 'blocked',
             kappa_confirmed=step['kappa'], kappa_optimistic=step['kappa_optimistic'],
             phi_lower=step['lower'], phi_upper=step['upper'], cumulative_survey_cost=step['spent'])
        for k, step in enumerate(steps)])
    write_table(OUTPUT/'table3.csv', [dict(policy=name, **values) for name, values in policies.items()])
    # 7. 摘要
    summary = dict(
        case='Case33-S', ports=list(network.load_nodes), partition='++', current_limit_a=Case33S.current_limit,
        voltage_pu=[.9, 1.1], port_pf=LOAD_PF, radial_networks=87, status_quo_kw=list(Case33S.status_quo),
        prior=list(PRIOR), rho=RHO, fallback=FALLBACK, levels=LEVELS, box_kw=list(P.bounds[2:]),
        families=[dict(roads=row['roads'], budget=key[1], area_kw2=unions[key].area, gap=row['gap'], status=row['status'],
                       seconds=row['seconds']) for key, row in families.items()],
        phi={label(A): phi[A] for A in SUBSETS}, J0=survey.J[0, 0], vpi=phi[0]-survey.weights@phi, policies=policies,
        sensitivity=dict(q0=GRID_Q0, rho=GRID_RHO, bundle_gap=gaps[0].tolist(), ratio_gap=gaps[1].tolist()),
        spot_check=spot, example=[dict(A=label(s['A']), B=label(s['B']), surveyed=None if s['action'] is None else
                                       ROADS[s['action']], available=s['available'], kappa=s['kappa'],
                                       kappa_optimistic=s['kappa_optimistic'], phi_lower=s['lower'], phi_upper=s['upper'],
                                       cumulative_survey_cost=s['spent']) for s in steps])
    (OUTPUT/'summary.json').write_bytes(json.dumps(summary, ensure_ascii=False, indent=2,
                                                   default=float).encode('utf-8'))
    # 8. 图（文字为英文）
    draw_survey_layout(network, OUTPUT/'fig_layout')
    draw_survey_valuation(LEVELS, {f"A = {'G' if A == FULL else label(A)}": kappa[A]
                                   for A in (0, *(1 << g for g in range(len(ROADS))), FULL)},
                          ROADS, interaction, OUTPUT/'fig_valuation')
    panels = [dict(level=cost_map(unions, s['A'], P),
                   optimistic=unions[family(FULL & ~s['B'], LEVELS[-1])].intersection(P),
                   title=('start' if k == 0 else f"{ROADS[steps[k-1]['action']]} " +
                          ('available' if steps[k-1]['available'] else 'blocked')),
                   kappa=(s['kappa'], s['kappa_optimistic']), phi=(s['lower'], s['upper'])) for k, s in enumerate(steps)]
    draw_survey_storyboard(panels, P.bounds[2:], LEVELS, Case33S.status_quo,
                           {key: [s[key] for s in steps] for key in ('spent', 'kappa', 'kappa_optimistic', 'lower', 'upper')},
                           OUTPUT/'fig_storyboard')
    draw_survey_policy({name: values['expected_cost'] for name, values in policies.items()}, survey.J[0, 0], gaps,
                       GRID_Q0, GRID_RHO, OUTPUT/'fig_policy')
    print(f"J(s0)={survey.J[0, 0]:.4f}；" + '；'.join(f"{name} {values['expected_cost']:.4f}" for name, values in policies.items()))
    print('抽查：'+'；'.join(f"{row['family']} {row['networks']} 个网架，相对误差 {row['relative_error']:+.3%}" for row in spot))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Case33-S 的走廊估值与勘察策略')
    parser.add_argument('--rerun', action='store_true', help='重算全部族的 RCUT（否则复用 families.json）')
    main(parser.parse_args().rerun)
