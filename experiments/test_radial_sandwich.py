"""径向锥剖分 + 内外夹逼：Case33(18,25)、mode=1、d=2 的构域实验，并与配对 AC/SOCP 扫描逐格比较。

在项目根目录运行::

    python -X utf8 experiments/test_radial_sandwich.py --partitions nn --tau 0.02 --output results/radial/smoke
    python -X utf8 experiments/test_radial_sandwich.py --output results/radial/case33_18_25
    python -X utf8 experiments/test_radial_sandwich.py --inner-cert near --output results/radial/case33_18_25_near \
        --compare-with results/radial/case33_18_25

原理
----
可行域为 R = ∪_x R_x。每个 R_x 经 SOCP + OBBT 紧化后是凸集，且含原点（零接入可行），所以 R 对原点星形。
因此 R 由径向函数 rho(θ) = max_x rho_x(θ) 完全决定。每个方向区间只需解一个 x 自由的 MISOCP，分支定界
自己选网架，不枚举网架，也没有“凸包之外”的面选择变量。内域由“单个网架 + 凸性”给出证书，外域由 MISOCP
的上界给出证书，按区间细分，直到径向相对误差 <= tau。

证书
----
内域：区间 k=[θa,θb] 的方案 x̂ 已 OBBT 紧化且原点经 SP 认证可行；va、vb 是从原点出发、带锥裕量的
紧化射线可行点（不用 MISOCP 现任解），故 T_k = conv(0, va, vb) ⊆ R_x̂ ⊆ R。原点不可行的方案只出现在
外域 MISOCP 中，不作任何区间的 x̂。
外域：Q=[va vb]，c=1ᵀQ⁻¹，锥内 xi=αa·va+αb·vb 满足 αa+αb=c·xi。x 自由 MISOCP 在锥内最大化 c·xi，
已紧化方案带其提升行（按汉明距离对全部 x 有效），其余方案为纯 SOCP，故 mu=ObjBound 在任何终止状态下
都是 R∩锥 上的有效上界：O_k={xi∈锥: c·xi<=mu}=mu·T_k ⊇ R∩锥，间隙 g=mu-1；mu<=1+tau 即认证。

近端覆盖证书（--inner-cert near）
----
Case33 闭合 18-33 的方案零接入欠压，原点 ∉ R_x，但其可行段沿射线 [na, va] 不从 0 开始。对这类 x̂：
远端 va、vb 仍由原点出发的射线求得，近端 na、nb 由远端朝原点的反向射线求得。由 R_x̂ 凸，
conv(na, va, vb, nb) ⊆ R_x̂；若某个原点认证可行的方案 y 在 θa、θb 的半径不小于 |na|、|nb|，则
conv(0, na, nb) ⊆ conv(0, ya, yb) ⊆ R_y。二者之并为 T_k = conv(0, va, vb)，故 T_k ⊆ R_x̂ ∪ R_y ⊆ R。
内域数据结构、外域 MISOCP 与逐格评价都不变；原点可行时 na=nb=0，退化为上面的原证书。候选在
{父 x̂, x*} 之外增加父区间的覆盖方案 y（原点可行，总可作 x̂）。默认 --inner-cert origin 保持原证书。

记号（沿用 docs/notation.md）
----
power/p 为 kW，顺序同 load_nodes；bounds=port_bounds(network)；分区 sign∈{±1}²，内部非负幅值 u，
真实 p=sign*u；xi=u/bounds；θ∈[0,π/2] 定义在 xi 平面，e(θ)=(cos θ, sin θ)；x 为型号向量，
方案开关串按走廊顺序列出七个可变开关（1 为闭合）；tau 为径向精度。θ、rho、va/vb、mu、g、Q、c、α
是本实验的局部量，只在本脚本中使用，未写入登记表。结果文件中的 *_kw 坐标是带符号的真实 kW。

实现约定（任务说明未覆盖或需取舍之处）
----
1. 网络与 main.py 相同取 Case33(load_nodes=(18,25), current_limit=CURRENT_LIMIT)：构造器默认 inf 限流
   的指纹与扫描缓存不一致；启动时校验 network.fingerprint 与 ac_identity。
2. 懒惰 OBBT：每次构建 x 自由 MISOCP 都加全部已紧化方案的行，BarQCPConvTol 取 1e-6；解的方案未紧化
   时 OBBT 后重建重解。上界已 <= 1+tau 时证书与紧化无关，不再触发。
3. π/4 径向 MISOCP 只用于选初始方案：胜出方案原点不可行时以 no-good 排除后重解，直到胜出方案原点
   可行；原始胜出方案记为 x̂0 写入结果。外域 MISOCP 从不排除方案。近端覆盖模式下初始区间的 x̂ 取
   {x̂0, 原点可行的胜出方案} 中满足证书且两端半径乘积最大者，后者同时作覆盖方案。
4. θ* 不在 (θa+δ, θb-δ) 内时取中点，δ=SPLIT_MARGIN·(θb-θa)；某候选点上无可用内顶点时也改取中点。
5. 待处理区间按继承的父区间 g 从大到小处理。分区时限到达时，未求解区间沿用父区间外界 c_p·xi<=mu_p
   （在子锥上仍有效），从未得到上界时外域取评价盒；这些区间保持 pending。
6. --criterion volume（默认 radial 不变）：与三维版相同的体积缺口准则，Δ_k=(μ̄_k^2-1)·area(T_k)，
   area(T_k)=|cross(va,vb)|/2；每次细分 Δ_k 最大的区间，子区间立即求外界，ΣΔ_k <= ε·Σarea(T_k) 时停止，
   ε 默认 d·tau；区间 MISOCP 解到相对间隙 --mip-gap，不用 1+tau 提前停止。run_volume 可续跑：已有区间时
   不重建根区间，先求解 pending 区间，再从全部 bounded 区间重建堆（compare_plan2 的方法 H 在 A 阶段后续跑）。

输出（--output 下）
----
intervals.json：各分区最终区间与细分历史；summary.json：分区状态、计数、方案、逐格指标（AC/SOCP 的 FR/MR、
外域有效性、夹逼带、三角形公式面积，另给外域与评价盒之交的面积 outer_box_area_kw2）以及 stuck_ranges
（相邻未认证区间合并后的角度段、MISOCP 状态与耗时、解方案的原点认证）；radial.png：带符号 p 平面图。
intervals.json 的 history 是每次区间决策后的收敛快照（最大径向间隙上界、面积间隙、认证角度占比、逐格
MR/FR 与外域有效性），画成 convergence.png。--compare-with DIR 另写 comparison.json（同网格逐分区对比）、
comparison.png（左 DIR、右本次）；DIR 也带 history 时再写 convergence_comparison.png（横轴为分区耗时）。
本脚本不修改主线文件，结果不加入 results/manifest.json。
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from dataclasses import dataclass
import heapq
import json
from pathlib import Path
import subprocess
import sys
from time import perf_counter
import traceback

import gurobipy as gp
from gurobipy import GRB
import numpy as np
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Network.case33bw import Case33, CURRENT_LIMIT
from model import MasterProblem, PortPhysics, PortSubProblem, port_bounds, ray_support
from monitor import comparison_metrics
from vertify import ac_identity

LOAD_NODES = (18, 25)
MODE = 1
PARTITIONS = ('pp', 'pn', 'np', 'nn')
SCAN_CACHE = (ROOT/'results'/'scan'/'case33bw'/'18_25'
              /'d0121a1ad0ba8cc1e009da6cc204ad9cf53881d3e5f4c488a439cb11d5f6fb69'/'region_17ceb84375dc2cc1.npz')
DEFAULT_OUTPUT = ROOT/'results'/'radial'/'case33_18_25'
COVER_CONV_TOL = 1e-6    # 紧化后 x 自由 MISOCP 的 barrier 收敛容差，同 RemainingRegionModel
SPLIT_MARGIN = .05       # δ=SPLIT_MARGIN*(θb-θa)
RAY_SECONDS = 30.        # 单次射线与原点 SP 的时限上限，同时受分区剩余时间限制
MIN_RADIUS = 1e-9        # 内顶点半径下限（xi）；更小时 Q 奇异，不作内顶点
ALPHA_TOL = 1e-12        # 格点重心坐标 α>=0 的舍入容差
STATUS = {getattr(GRB.Status, name): name for name in dir(GRB.Status) if name.isupper()}


def direction(theta):
    """e(θ)=(cos θ, sin θ)；θ=0、π/2 时精确取坐标轴。"""
    e = np.array([np.cos(theta), np.sin(theta)])
    e[np.abs(e) < 1e-15] = 0.
    return e


def cross(a, b):
    return float(a[0]*b[1]-a[1]*b[0])


def polygon_area(vertices):
    """多边形面积（鞋带公式）；三角形即 |cross(va, vb)|/2。"""
    v = np.asarray(vertices).reshape(-1, 2)
    if len(v) < 3:
        return 0.
    return .5*abs(float(np.sum(v[:, 0]*np.roll(v[:, 1], -1)-v[:, 1]*np.roll(v[:, 0], -1))))


def clip_box(vertices):
    """第一象限多边形与评价盒 [0,1]² 的交（Sutherland–Hodgman，只需裁 xi_j<=1）。"""
    polygon = [np.asarray(v, float) for v in vertices]
    for axis in (0, 1):
        clipped = []
        for k, current in enumerate(polygon):
            previous = polygon[k-1]
            if (current[axis] <= 1.) != (previous[axis] <= 1.):
                clipped.append(previous+(1.-previous[axis])/(current[axis]-previous[axis])*(current-previous))
            if current[axis] <= 1.:
                clipped.append(current)
        polygon = clipped
    return np.array(polygon).reshape(-1, 2)


def partition_sign(code, d=2):
    """代号顺序同 load_nodes：p 为正功率（负荷），n 为负功率（光伏）；d 为功率维数。"""
    if len(code) != d or set(code)-{'p', 'n'}:
        raise argparse.ArgumentTypeError(f'partition code must be {d} letters from p/n: {code!r}')
    return tuple(1 if c == 'p' else -1 for c in code)


def plain(value):
    """JSON 友好：数组转列表，非有限浮点写 null。"""
    if isinstance(value, dict):
        return {str(k): plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [plain(v) for v in value]
    if isinstance(value, np.ndarray):
        return plain(value.tolist())
    if isinstance(value, np.generic):
        return plain(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


class PartitionTimeout(RuntimeError):
    """分区时限已到；已得证书保留，未求解区间保持 pending。"""


@dataclass(eq=False)
class Interval:
    """方向区间 k=[θa,θb]：内域 T_k=conv(0,va,vb)，外域 O_k={xi∈锥: c·xi<=mu}；va/vb 为 xi 坐标。"""
    id: int
    theta_a: float
    theta_b: float
    x: tuple                          # 内域方案 x̂（0/1 型号向量）
    va: np.ndarray                    # rho_x̂(θa)·e(θa)
    vb: np.ndarray                    # rho_x̂(θb)·e(θb)
    parent: int | None = None
    depth: int = 0
    priority: float = np.inf          # 处理顺序：父区间的 g，初始区间为 inf
    inherited: tuple | None = None    # 父区间外界 (c, mu)；本区间尚无有限上界时沿用
    mu: float = np.inf                # 外界因子
    g: float = np.inf                 # 间隙 mu-1
    status: str = 'pending'           # pending / certified / unresolved
    reason: str | None = None
    misocp: dict | None = None
    cover: tuple | None = None        # 近端覆盖方案 y；x̂ 原点可行时为 None
    near_a: np.ndarray | None = None  # x̂ 在 θa 的可行段近端（xi），原点可行时为 0
    near_b: np.ndarray | None = None
    cover_a: float | None = None      # 覆盖方案 y 在 θa、θb 的半径 rho_y
    cover_b: float | None = None
    final: dict | None = None         # 最后一轮 MISOCP 记录（体积准则细分时用解点与方案）

    @property
    def width(self):
        return self.theta_b-self.theta_a

    @property
    def Q(self):
        return np.column_stack([self.va, self.vb])

    @property
    def c(self):
        """c=1ᵀQ⁻¹：锥内 xi=αa·va+αb·vb 时 αa+αb=c·xi。"""
        return np.linalg.solve(self.Q.T, np.ones(2))

    def halfspace(self):
        """有效外界 (c, mu)：自身有有限上界用自身，否则沿用父区间；均无时为 None（外域取评价盒）。"""
        return (self.c, self.mu) if np.isfinite(self.mu) else self.inherited

    def outer_polygon(self):
        """外域多边形（xi，逆时针，首点为原点）；自身上界时即 (0, mu·va, mu·vb)。"""
        ea, eb = direction(self.theta_a), direction(self.theta_b)
        halfspace = self.halfspace()
        if halfspace is None:
            corner = [np.ones(2)] if self.theta_a < np.pi/4 < self.theta_b else []
            return np.array([np.zeros(2), ea/ea.max(), *corner, eb/eb.max()])
        c, mu = halfspace
        return np.array([np.zeros(2), mu/(c@ea)*ea, mu/(c@eb)*eb])


class RadialSandwich:
    """单个符号分区的径向夹逼：x̂ 只取原点认证可行的已紧化方案，外域 MISOCP 对 x 自由。"""
    CENTER_CONTEXT = 'radial pi/4'

    def __init__(self, network, code, bounds, budget, args):
        self.code, self.sign = code, np.asarray(partition_sign(code, len(network.load_nodes)))
        self.dim = len(network.load_nodes)
        self.bounds, self.budget, self.args, self.tau = np.asarray(bounds, float), budget, args, args.tau
        self.started = perf_counter()
        self.deadline = self.started+args.seconds
        self.network = deepcopy(network)
        self.equations = PortPhysics(self.network, self.sign)
        self.oracle = PortSubProblem(self.equations, threads=args.threads)
        self.switches = [s for c, s in zip(self.network.corridors, self.network.type_slices) if c.switchable]
        self.near_mode = args.inner_cert == 'near'
        self.criterion = getattr(args, 'criterion', 'radial')
        epsilon = getattr(args, 'epsilon', None)
        self.epsilon = epsilon if epsilon is not None else self.dim*args.tau   # 体积缺口容差，默认 d·tau
        self.rays, self.origin = {}, {}          # (方案, θ) → 内顶点或 None；方案 → 原点认证
        self.nears = {}                          # (方案, θ) → 可行段近端（xi）或 None，仅原点不可行方案
        self.counts = dict(misocp=0, interval_solves=0, rays=0, ray_cache_hits=0, ray_failures=0,
                           near_rays=0, near_failures=0, obbt=0, lazy_obbt=0, origin_sp=0)
        self.lazy, self.splits, self.errors, self.history = [], [], [], []
        self.intervals, self.selection = [], {}
        self.status, self.seconds, self.next_id = 'pending', 0., 0

    # ---- 基础求解：紧化、原点认证、射线内顶点 -------------------------------------------------
    def label(self, x):
        x = np.asarray(x)
        return ''.join(str(int(x[s].any())) for s in self.switches)

    def remaining(self):
        left = self.deadline-perf_counter()
        if left <= 0.:
            raise PartitionTimeout(f'{self.code}: partition time limit')
        return left

    def tighten(self, x):
        """方案首次使用前 OBBT；返回耗时，已紧化时为 0。"""
        scheme = tuple(int(v) for v in x)
        if scheme in self.equations.boxes:
            return 0.
        self.remaining()
        start = perf_counter()
        self.equations.obbt(np.asarray(scheme), self.bounds, self.args.obbt_workers)
        self.counts['obbt'] += 1
        return perf_counter()-start

    def _limited(self, solve, failure):
        """单次射线/SP：时限取 RAY_SECONDS 与剩余时间的较小者；分区时限内的数值失败按 failure 记录。"""
        try:
            return solve(min(RAY_SECONDS, self.remaining()))
        except TimeoutError as exc:
            if perf_counter() >= self.deadline:
                raise PartitionTimeout(f'{self.code}: partition time limit') from exc
            return failure(exc)
        except RuntimeError as exc:
            if isinstance(exc, PartitionTimeout):
                raise
            return failure(exc)

    def origin_feasible(self, x):
        """零接入（p=0）的紧化 SP 认证，按方案缓存；数值失败视为未认证。"""
        scheme = tuple(int(v) for v in x)
        if scheme not in self.origin:
            self.tighten(scheme)
            self.counts['origin_sp'] += 1

            def solve(limit):
                answer = self.oracle.solve(np.asarray(scheme), np.zeros(self.dim), time_limit=limit, score_only=True)
                return dict(feasible=bool(answer['feasible']), eta=float(answer['eta']))
            self.origin[scheme] = self._limited(solve, lambda exc: dict(feasible=False, eta=None, error=str(exc)))
        return self.origin[scheme]['feasible']

    def unit(self, theta):
        """方向键 → xi 空间单位方向；二维键为角度 θ，三维子类改为方向编号。"""
        return direction(theta)

    def vertex(self, x, theta):
        """rho_x(θ)·e(θ)（xi）：从原点朝评价盒边界点 bounds*e/max(e) 的紧化射线，按 (方案, θ) 缓存。"""
        scheme = tuple(int(v) for v in x)
        key = scheme, float(theta)
        if key in self.rays:
            self.counts['ray_cache_hits'] += 1
            return self.rays[key]
        self.tighten(scheme)
        e = self.unit(theta)
        target = self.bounds*e/e.max()
        self.counts['rays'] += 1

        def solve(limit):
            answer = ray_support(self.equations, self.budget, np.asarray(scheme), np.zeros(self.dim), target,
                                 threads=self.args.threads, time_limit=limit, numeric_focus=0)
            point = answer['p']/self.bounds
            if np.linalg.norm(point) < MIN_RADIUS:
                raise RuntimeError(f'Ray radius {np.linalg.norm(point):g} below {MIN_RADIUS:g}')
            return point

        def failure(exc):
            self.counts['ray_failures'] += 1
            self.errors.append(dict(kind='ray', scheme=self.label(scheme), theta=float(theta), error=str(exc)))
            return None
        self.rays[key] = self._limited(solve, failure)
        return self.rays[key]

    def near(self, x, theta):
        """原点不可行方案在 θ 方向可行段的近端（xi）：从远端 rho_x(θ)·e(θ) 朝原点的反向紧化射线，按 (方案, θ) 缓存。"""
        scheme = tuple(int(v) for v in x)
        key = scheme, float(theta)
        if key in self.nears:
            return self.nears[key]
        far = self.vertex(scheme, theta)
        if far is None:
            self.nears[key] = None
            return None
        self.counts['near_rays'] += 1

        def solve(limit):
            answer = ray_support(self.equations, self.budget, np.asarray(scheme), far*self.bounds, np.zeros(self.dim),
                                 threads=self.args.threads, time_limit=limit, numeric_focus=0)
            return answer['p']/self.bounds

        def failure(exc):
            self.counts['near_failures'] += 1
            self.errors.append(dict(kind='near_ray', scheme=self.label(scheme), theta=float(theta), error=str(exc)))
            return None
        self.nears[key] = self._limited(solve, failure)
        return self.nears[key]

    def inner_vertices(self, x, lo, hi, covers=()):
        """方案 x 在 [lo, hi] 上的内顶点及其证书，不满足时为 None。原点可行时 T=conv(0,va,vb) ⊆ R_x；
        近端覆盖模式下原点不可行的 x 需某个原点可行方案 y（先试 covers，再试已认证方案）在两端的半径不小于近端。"""
        feasible = self.origin_feasible(x)
        if not (feasible or self.near_mode):
            return None
        va, vb = self.vertex(x, lo), self.vertex(x, hi)
        if va is None or vb is None:
            return None
        if feasible:
            return dict(x=tuple(x), va=va, vb=vb, cover=None, near_a=np.zeros(2), near_b=np.zeros(2),
                        cover_a=None, cover_b=None)
        na, nb = self.near(x, lo), self.near(x, hi)
        if na is None or nb is None:
            return None
        ra, rb = float(np.linalg.norm(na)), float(np.linalg.norm(nb))
        order = [*covers, *(s for s, row in self.origin.items() if row['feasible'])]
        for y in dict.fromkeys(tuple(s) for s in order if s is not None):
            if not self.origin_feasible(y):
                continue
            ya, yb = self.vertex(y, lo), self.vertex(y, hi)
            if ya is not None and yb is not None and np.linalg.norm(ya) >= ra and np.linalg.norm(yb) >= rb:
                return dict(x=tuple(x), va=va, vb=vb, cover=y, near_a=na, near_b=nb,
                            cover_a=float(np.linalg.norm(ya)), cover_b=float(np.linalg.norm(yb)))
        return None

    # ---- x 自由 MISOCP 与懒惰 OBBT ---------------------------------------------------------
    def misocp(self, configure, *, stop, context):
        """每次构建都加全部已紧化方案的行；解的方案未紧化则 OBBT 后重建重解。
        退出：不可行、无解、上界 <= stop、解的方案已紧化，或分区时限已到。返回各轮记录与触发次数。"""
        rounds, triggers = [], 0
        while True:
            if rounds and perf_counter() >= self.deadline:
                break
            limit = min(self.args.mip_seconds, self.remaining())
            start = perf_counter()
            problem = MasterProblem(self.equations, budget=self.budget, threads=self.args.threads)
            model = problem.model
            try:
                x_dict = dict(zip(self.equations.keys, problem.x.tolist()))
                for scheme in self.equations.boxes:
                    self.equations._tighten(model, x_dict, problem.operation, scheme)
                model.Params.BarQCPConvTol = COVER_CONV_TOL
                xi = [problem.loads[i]*(1./float(b)) for i, b in zip(self.network.load_nodes, self.bounds)]
                configure(problem, xi)
                model.Params.TimeLimit = limit
                built = perf_counter()
                model.optimize()
                self.counts['misocp'] += 1
                if model.Status == GRB.INTERRUPTED:
                    raise KeyboardInterrupt
                try:
                    bound = float(model.ObjBound)
                except gp.GurobiError:
                    bound = np.inf
                record = dict(status=STATUS.get(model.Status, str(model.Status)), build_seconds=built-start,
                              solve_seconds=perf_counter()-built, tightened_schemes=len(self.equations.boxes),
                              objective=None, bound=bound if abs(bound) < 1e30 else np.inf, x=None, point=None)
                if model.SolCount:
                    record.update(objective=float(model.ObjVal), x=tuple(int(v) for v in np.rint(problem.x.X)),
                                  point=np.asarray(problem.power.X, float))
                self.on_solution(problem, record)
            finally:
                model.dispose()
            rounds.append(record)
            if (record['status'] == 'INFEASIBLE' or record['x'] is None
                    or (stop is not None and record['bound'] <= stop)
                    or record['x'] in self.equations.boxes or perf_counter() >= self.deadline):
                break
            seconds = self.tighten(record['x'])
            triggers += 1
            self.counts['lazy_obbt'] += 1
            record['lazy_obbt_seconds'] = seconds
            self.lazy.append(dict(context=context, scheme=self.label(record['x']), seconds=seconds,
                                  objective=record['objective'], bound=record['bound']))
            print(f'  [{self.code}] lazy OBBT {context}: scheme {self.label(record["x"])} '
                  f'(incumbent {record["objective"]:.6g}, bound {record["bound"]:.6g}) {seconds:.2f}s', flush=True)
        return rounds, triggers

    def on_solution(self, problem, record):
        """每轮 MISOCP 求解后、释放模型前调用；默认无操作。compare_plan2 用它审计现任解并记录求解。"""

    def center_ray(self, problem, xi):
        """初始方案的径向 MISOCP：二维取 θ=π/4，sin θ·xi1 − cos θ·xi2 = 0，max xi1+xi2。"""
        theta = np.pi/4
        problem.model.addConstr(float(np.sin(theta))*xi[0]-float(np.cos(theta))*xi[1] == 0., name='radial_direction')
        problem.model.setObjective(xi[0]+xi[1], GRB.MAXIMIZE)

    def initial_scheme(self):
        """中心方向的径向 MISOCP 选初始方案；胜出方案原点不可行时 no-good 排除后重解。"""
        excluded, history = [], []
        while True:
            def configure(problem, xi):
                self.center_ray(problem, xi)
                for scheme in excluded:
                    problem.exclude(np.asarray(scheme))
            rounds, triggers = self.misocp(configure, stop=None, context=self.CENTER_CONTEXT)
            final = rounds[-1]
            history.append(dict(excluded=[self.label(s) for s in excluded], lazy_obbt=triggers,
                                rounds=[self.round_json(r) for r in rounds]))
            if final['x'] is None:
                raise RuntimeError(f'{self.code}: radial MISOCP returned no scheme ({final["status"]})')
            self.selection.setdefault('x0', final['x'])
            if self.origin_feasible(final['x']):
                self.selection.update(chosen=final['x'], excluded=list(excluded), history=history)
                return final['x']
            print(f'  [{self.code}] radial winner {self.label(final["x"])} has an infeasible origin; excluded',
                  flush=True)
            excluded.append(final['x'])

    def certify(self, interval):
        """外界 MISOCP：锥约束、幅值上界 bounds、目标 max c·xi。径向准则 BestBdStop=BestObjStop=1+tau；
        体积准则解到相对间隙 --mip-gap，不提前停止（同三维版）。"""
        va, vb, c = interval.va, interval.vb, interval.c
        volume = self.criterion == 'volume'

        def configure(problem, xi):
            model = problem.model
            model.addConstr(float(va[0])*xi[1]-float(va[1])*xi[0] >= 0., name='cone_a')   # cross(va, xi) >= 0
            model.addConstr(float(vb[1])*xi[0]-float(vb[0])*xi[1] >= 0., name='cone_b')   # cross(xi, vb) >= 0
            problem.power.UB = self.bounds
            model.setObjective(float(c[0])*xi[0]+float(c[1])*xi[1], GRB.MAXIMIZE)
            if volume:
                model.Params.MIPGap = self.args.mip_gap
            else:
                model.Params.BestBdStop = model.Params.BestObjStop = 1.+self.tau
        self.counts['interval_solves'] += 1
        rounds, triggers = self.misocp(configure, stop=None if volume else 1.+self.tau,
                                       context=f'interval {interval.id}')
        final = rounds[-1]
        interval.mu = 1. if final['status'] == 'INFEASIBLE' else final['bound']
        interval.g = interval.mu-1.
        interval.misocp = dict(status=final['status'], objective=final['objective'], bound=final['bound'],
                               seconds=sum(r['build_seconds']+r['solve_seconds'] for r in rounds),
                               lazy_obbt=triggers, rounds=rounds)
        return final

    # ---- 细分与主循环 ------------------------------------------------------------------------
    def new_interval(self, theta_a, theta_b, inner, parent=None):
        interval = Interval(self.next_id, float(theta_a), float(theta_b), inner['x'], inner['va'], inner['vb'],
                            cover=inner['cover'], near_a=inner['near_a'], near_b=inner['near_b'],
                            cover_a=inner['cover_a'], cover_b=inner['cover_b'])
        self.next_id += 1
        if parent is not None:
            interval.parent, interval.depth, interval.priority = parent.id, parent.depth+1, parent.g
            if np.isfinite(parent.mu):
                interval.inherited = parent.c, parent.mu
        return interval

    def split(self, interval, final):
        """在 θ*（解点方向）或中点细分；子区间 x̂ 取 {父 x̂, 解的方案 x*}（近端覆盖模式再加父区间覆盖方案 y）
        中满足内域证书且两端半径乘积最大者。"""
        a, b = interval.theta_a, interval.theta_b
        options, delta = [], SPLIT_MARGIN*interval.width
        if final['point'] is not None:
            point = final['point']/self.bounds
            theta = float(np.arctan2(point[1], point[0]))
            if a+delta < theta < b-delta:
                options.append((theta, 'incumbent'))
        options.append((.5*(a+b), 'midpoint'))
        candidates = [interval.x]
        if final['x'] is not None and final['x'] != interval.x:
            candidates.append(final['x'])
        if self.near_mode and interval.cover is not None and interval.cover not in candidates:
            candidates.append(interval.cover)
        covers = (interval.cover, interval.x)
        for theta, rule in options:
            children = []
            for lo, hi in ((a, theta), (theta, b)):
                best = None
                for scheme in candidates:
                    inner = self.inner_vertices(scheme, lo, hi, covers)
                    if inner is None:
                        continue
                    score = float(np.linalg.norm(inner['va'])*np.linalg.norm(inner['vb']))
                    if best is None or score > best[0]:
                        best = score, inner
                if best is None:
                    break
                children.append(self.new_interval(lo, hi, best[1], parent=interval))
            if len(children) == 2:
                incumbent = None if final['x'] is None else self.label(final['x'])
                self.splits.append(dict(
                    id=interval.id, theta_a=a, theta_b=b, scheme=self.label(interval.x), mu=interval.mu,
                    g=interval.g, misocp_status=interval.misocp['status'],
                    misocp_seconds=interval.misocp['seconds'], incumbent=interval.misocp['objective'],
                    incumbent_scheme=incumbent,
                    incumbent_origin_feasible=(None if final['x'] is None else
                                               self.origin.get(final['x'], {}).get('feasible')),
                    split_theta=theta, rule=rule, children=[child.id for child in children],
                    child_schemes=[self.label(child.x) for child in children],
                    child_covers=[None if child.cover is None else self.label(child.cover) for child in children]))
                return children
        return None

    def leaf_gap(self, interval):
        """区间当前的径向间隙上界：自身外界为 mu-1；沿用父区间外界时为 max_θ mu_p(c·e)/(c_p·e)-1，
        线性分式在端点取极值；尚无外界为 inf。"""
        if np.isfinite(interval.mu):
            return interval.mu-1.
        if interval.inherited is None:
            return np.inf
        c_p, mu_p = interval.inherited
        return max(mu_p*(interval.c@e)/(c_p@e) for e in map(direction, (interval.theta_a, interval.theta_b)))-1.

    def snapshot(self, interval, decision, scan):
        """每次区间决策后的收敛快照：最大径向间隙上界、面积间隙、认证角度占比、逐格 MR/FR 与外域有效性。"""
        leaves = self.intervals
        inner, _, outer_box = self.areas()
        row = dict(step=self.counts['interval_solves'], seconds=perf_counter()-self.started,
                   interval=None if interval is None else interval.id, decision=decision,
                   mu=None if interval is None else interval.mu, leaves=len(leaves),
                   certified=sum(i.status == 'certified' for i in leaves),
                   certified_angle=sum(i.width for i in leaves if i.status == 'certified')/(np.pi/2),
                   max_gap=max(self.leaf_gap(i) for i in leaves), inner_area_kw2=inner, outer_box_area_kw2=outer_box,
                   misocp=self.counts['misocp'], rays=self.counts['rays']+self.counts['near_rays'])
        if scan is not None:
            labels = self.classify(scan['power'])
            ac = scan['ac'][labels['cells']] == 1
            accuracy = comparison_metrics(labels['inner'], ac)
            row.update(inner_ac_mr=accuracy['mr_percent'], inner_ac_fr=accuracy['fr_percent'],
                       ac_out=int((ac & ~labels['outer']).sum()))
        self.history.append(row)

    def make_root(self, scan=None):
        """径向 MISOCP 选初始方案，建根区间 [0, π/2]。"""
        x0 = self.initial_scheme()
        candidates = [self.selection['x0'], x0] if self.near_mode and self.selection['x0'] != x0 else [x0]
        best = None
        for scheme in candidates:
            inner = self.inner_vertices(scheme, 0., np.pi/2, (x0,))
            if inner is not None:
                score = float(np.linalg.norm(inner['va'])*np.linalg.norm(inner['vb']))
                if best is None or score > best[0]:
                    best = score, inner
        if best is None:
            raise RuntimeError(f'{self.code}: initial scheme {self.label(x0)} has no ray vertex at 0 or pi/2')
        root = self.new_interval(0., np.pi/2, best[1])
        self.intervals = [root]
        self.snapshot(None, 'initial', scan)
        return root

    def solve_interval(self, interval):
        """求区间外界并保存最后一轮记录（体积准则细分用）。"""
        final = self.certify(interval)
        interval.final = final
        info = interval.misocp
        print(f'  [{self.code}] #{interval.id:<3d} [{interval.theta_a:.5f},{interval.theta_b:.5f}] '
              f'xhat={self.label(interval.x)} {info["status"]:<14s} mu={interval.mu:.6g} '
              f'lazy={info["lazy_obbt"]} {info["seconds"]:.2f}s leaves={len(self.intervals)}', flush=True)
        return final

    def interval_delta(self, interval):
        """体积缺口 Δ=(μ̄^d-1)·area(T)（xi²）；μ̄=1+leaf_gap，尚无外界时为 inf。"""
        gap = self.leaf_gap(interval)
        area = .5*abs(cross(interval.va, interval.vb))
        return max(((1.+gap)**self.dim-1.)*area, 0.) if np.isfinite(gap) else np.inf

    def volume_ratio(self):
        """体积缺口比 ΣΔ_k / Σarea(T_k)；任一叶区间尚无外界时为 inf。"""
        total = sum(.5*abs(cross(i.va, i.vb)) for i in self.intervals)
        return sum(self.interval_delta(i) for i in self.intervals)/total if total > 0. else np.inf

    def run_volume(self, scan=None):
        """体积缺口准则（三维版 run_volume 的二维实例）：每次细分 Δ_k 最大的区间，子区间立即求外界；
        ΣΔ_k <= ε·Σarea(T_k) 时停止。已有区间时续跑：先求 pending 区间，再从全部 bounded 区间重建堆。"""
        try:
            if not self.intervals:
                root = self.make_root(scan)
                self.solve_interval(root)
                root.status = 'bounded'
                self.snapshot(root, 'bounded', scan)
            for interval in [i for i in self.intervals if i.status == 'pending']:
                self.solve_interval(interval)
                interval.status = 'bounded'
                self.snapshot(interval, 'bounded', scan)
            heap = [(-self.interval_delta(i), i.id, i) for i in self.intervals if i.status == 'bounded']
            heapq.heapify(heap)
            while True:
                if self.volume_ratio() <= self.epsilon:
                    self.status = 'certified'
                    break
                if not heap:
                    self.status = 'unresolved'
                    break
                self.remaining()
                _, _, interval = heapq.heappop(heap)
                if interval.width < self.args.min_width:
                    interval.status, interval.reason = 'unresolved', 'min_width'
                    continue
                if len(self.intervals) >= self.args.max_intervals:
                    interval.status, interval.reason = 'unresolved', 'max_intervals'
                    continue
                children = self.split(interval, interval.final)
                if children is None:
                    interval.status, interval.reason = 'unresolved', 'no_inner_vertex'
                    continue
                k = next(j for j, item in enumerate(self.intervals) if item is interval)
                self.intervals[k:k+1] = children
                for child in children:
                    self.solve_interval(child)
                    child.status = 'bounded'
                    heapq.heappush(heap, (-self.interval_delta(child), child.id, child))
                self.snapshot(interval, 'split', scan)
        except PartitionTimeout:
            self.status = 'time_limit'
        except Exception as exc:
            self.status = 'error'
            self.errors.append(dict(kind='fatal', error=repr(exc), traceback=traceback.format_exc()))
            print(traceback.format_exc(), flush=True)
        self.seconds = perf_counter()-self.started
        return self

    def run(self, scan=None):
        if self.criterion == 'volume':
            return self.run_volume(scan)
        try:
            root = self.make_root(scan)
            heap = [(-root.priority, -root.width, root.id, root)]
            while heap:
                self.remaining()
                *_, interval = heapq.heappop(heap)
                final = self.certify(interval)
                info = interval.misocp
                print(f'  [{self.code}] #{interval.id:<3d} [{interval.theta_a:.5f},{interval.theta_b:.5f}] '
                      f'xhat={self.label(interval.x)} {info["status"]:<14s} mu={interval.mu:.6g} '
                      f'lazy={info["lazy_obbt"]} {info["seconds"]:.2f}s leaves={len(self.intervals)}', flush=True)
                children = None
                if interval.mu <= 1.+self.tau:
                    interval.status = 'certified'
                elif interval.width < self.args.min_width:
                    interval.status, interval.reason = 'unresolved', 'min_width'
                elif len(self.intervals) >= self.args.max_intervals:
                    interval.status, interval.reason = 'unresolved', 'max_intervals'
                elif (children := self.split(interval, final)) is None:
                    interval.status, interval.reason = 'unresolved', 'no_inner_vertex'
                else:
                    k = next(j for j, item in enumerate(self.intervals) if item is interval)
                    self.intervals[k:k+1] = children
                    for child in children:
                        heapq.heappush(heap, (-child.priority, -child.width, child.id, child))
                self.snapshot(interval, 'split' if children else interval.reason or interval.status, scan)
            self.status = 'certified' if all(i.status == 'certified' for i in self.intervals) else 'unresolved'
        except PartitionTimeout:
            self.status = 'time_limit'
        except Exception as exc:
            self.status = 'error'
            self.errors.append(dict(kind='fatal', error=repr(exc), traceback=traceback.format_exc()))
            print(traceback.format_exc(), flush=True)
        self.seconds = perf_counter()-self.started
        return self

    # ---- 评估与导出 ------------------------------------------------------------------------
    def classify(self, power):
        """格心按符号归属本分区；按 θ 找区间，α=Q⁻¹xi。内域 α>=0 且 Σα<=1，外域 c·xi<=mu（无上界时为评价盒）。"""
        cells = np.all(np.where(power >= 0., 1, -1) == self.sign, axis=1)
        xi = np.abs(power[cells])/self.bounds
        inner, outer = np.zeros(len(xi), bool), np.ones(len(xi), bool)
        theta, index = np.arctan2(xi[:, 1], xi[:, 0]), np.full(len(xi), -1)
        if self.intervals:
            ends = np.array([interval.theta_b for interval in self.intervals])
            index = np.minimum(np.searchsorted(ends, theta, side='left'), len(ends)-1)
            for k, interval in enumerate(self.intervals):
                chosen = index == k
                if not chosen.any():
                    continue
                alpha = xi[chosen]@np.linalg.inv(interval.Q).T
                inner[chosen] = np.all(alpha >= -ALPHA_TOL, axis=1) & (alpha.sum(axis=1) <= 1.)
                halfspace = interval.halfspace()
                if halfspace is not None:
                    outer[chosen] = xi[chosen]@halfspace[0] <= halfspace[1]
                else:
                    outer[chosen] = np.all(xi[chosen] <= 1., axis=1)
        return dict(cells=cells, inner=inner, outer=outer, xi=xi, theta=theta, index=index)

    def areas(self):
        """内/外域面积（kW²，三角形公式）及外域与评价盒之交的面积；xi→u 为对角缩放，面积乘 prod(bounds)。"""
        scale = float(np.prod(self.bounds))
        if not self.intervals:
            return 0., scale, scale
        inner = sum(.5*abs(cross(i.va, i.vb)) for i in self.intervals)*scale
        outer = sum(polygon_area(i.outer_polygon()) for i in self.intervals)*scale
        outer_box = sum(polygon_area(clip_box(i.outer_polygon())) for i in self.intervals)*scale
        return inner, outer, outer_box

    def incumbent(self, interval):
        """区间最后一轮 MISOCP 的解方案及其原点认证（未检查过为 None）。"""
        final = interval.misocp['rounds'][-1] if interval.misocp else None
        if final is None or final['x'] is None:
            return None, None
        return self.label(final['x']), self.origin.get(final['x'], {}).get('feasible')

    def stuck_ranges(self):
        """相邻未认证区间合并为角度段：区间数、最小角宽、MISOCP 状态/耗时、mu 范围、解方案及其原点认证。"""
        groups = []
        for interval in self.intervals:
            if interval.status == 'certified':
                groups.append(None)
            elif groups and groups[-1] is not None:
                groups[-1].append(interval)
            else:
                groups.append([interval])
        rows = []
        for members in filter(None, groups):
            solved = [m for m in members if m.misocp is not None]
            statuses, reasons, schemes = {}, {}, {}
            for m in members:
                key = m.misocp['status'] if m.misocp else 'NOT_SOLVED'
                statuses[key] = statuses.get(key, 0)+1
                reasons[m.reason or m.status] = reasons.get(m.reason or m.status, 0)+1
                label, feasible = self.incumbent(m)
                if label is not None:
                    row = schemes.setdefault(label, dict(intervals=0, origin_feasible=feasible))
                    row['intervals'] += 1
            seconds = [m.misocp['seconds'] for m in solved]
            rows.append(dict(theta_a=members[0].theta_a, theta_b=members[-1].theta_b, intervals=len(members),
                             min_width=min(m.width for m in members), misocp_status=statuses, reasons=reasons,
                             mu_min=min((m.mu for m in solved), default=None),
                             mu_max=max((m.mu for m in solved), default=None),
                             misocp_seconds_total=sum(seconds), misocp_seconds_max=max(seconds, default=None),
                             inner_schemes=sorted({self.label(m.x) for m in members}),
                             incumbent_schemes=schemes))
        return rows

    def kw(self, point):
        return None if point is None else self.sign*np.asarray(point)*self.bounds

    def round_json(self, record):
        point = record['point']
        theta = None if point is None else float(np.arctan2(point[1]/self.bounds[1], point[0]/self.bounds[0]))
        return dict(status=record['status'], build_seconds=record['build_seconds'],
                    solve_seconds=record['solve_seconds'], tightened_schemes=record['tightened_schemes'],
                    objective=record['objective'], bound=record['bound'],
                    scheme=None if record['x'] is None else self.label(record['x']),
                    point_kw=None if point is None else self.sign*point, theta=theta,
                    lazy_obbt_seconds=record.get('lazy_obbt_seconds'))

    def interval_json(self, interval):
        halfspace = interval.halfspace()
        outer = interval.outer_polygon()
        info = interval.misocp
        source = 'box' if halfspace is None else ('own' if np.isfinite(interval.mu) else 'parent')
        incumbent, feasible = self.incumbent(interval)
        return dict(
            id=interval.id, parent=interval.parent, depth=interval.depth,
            theta_a=interval.theta_a, theta_b=interval.theta_b, width=interval.width,
            scheme=self.label(interval.x), status=interval.status, reason=interval.reason,
            va_xi=interval.va, vb_xi=interval.vb, va_kw=self.kw(interval.va), vb_kw=self.kw(interval.vb),
            rho_a=float(np.linalg.norm(interval.va)), rho_b=float(np.linalg.norm(interval.vb)),
            origin_feasible=interval.cover is None,
            cover_scheme=None if interval.cover is None else self.label(interval.cover),
            near_a_xi=interval.near_a, near_b_xi=interval.near_b,
            near_a_kw=self.kw(interval.near_a), near_b_kw=self.kw(interval.near_b),
            cover_rho_a=interval.cover_a, cover_rho_b=interval.cover_b,
            mu=interval.mu, g=interval.g, priority=interval.priority, outer_source=source,
            outer_kw=[self.kw(v) for v in outer[1:]],
            incumbent_scheme=incumbent, incumbent_origin_feasible=feasible,
            misocp=None if info is None else dict(
                status=info['status'], seconds=info['seconds'], objective=info['objective'],
                bound=info['bound'], lazy_obbt=info['lazy_obbt'],
                rounds=[self.round_json(r) for r in info['rounds']]))

    def result(self, evaluation=None):
        solved = [i for i in self.intervals if i.misocp is not None]
        used = {}
        for interval in self.intervals:
            row = used.setdefault(self.label(interval.x), dict(intervals=0, angle=0.))
            row['intervals'] += 1
            row['angle'] += interval.width
        blocking = sorted({s['incumbent_scheme'] for s in self.splits if s['incumbent_origin_feasible'] is False})
        inner_area, outer_area, outer_box_area = self.areas()
        statuses = [i.status for i in self.intervals]
        selection = self.selection
        covered = [i for i in self.intervals if i.cover is not None]
        cover_pairs = {}
        for interval in covered:
            key = f'{self.label(interval.x)}<-{self.label(interval.cover)}'
            row = cover_pairs.setdefault(key, dict(intervals=0, angle=0.))
            row['intervals'] += 1
            row['angle'] += interval.width
        return dict(
            partition=self.code, sign=self.sign, status=self.status, seconds=self.seconds,
            inner_cert=self.args.inner_cert, criterion=self.criterion, epsilon=self.epsilon,
            volume_gap_ratio=self.volume_ratio() if self.intervals else None, covered_intervals=len(covered),
            covered_angle=sum(i.width for i in covered), cover_pairs=cover_pairs,
            max_g=max((i.g for i in solved), default=None),
            intervals=len(self.intervals), certified=statuses.count('certified'), bounded=statuses.count('bounded'),
            unresolved=statuses.count('unresolved'), pending=statuses.count('pending'),
            unresolved_reasons={r: sum(i.reason == r for i in self.intervals)
                                for r in sorted({i.reason for i in self.intervals if i.reason})},
            counts=dict(self.counts), schemes_used=used,
            tightened_schemes=[self.label(s) for s in self.equations.boxes],
            origin={self.label(s): row for s, row in self.origin.items()},
            initial=dict(x0=None if 'x0' not in selection else self.label(selection['x0']),
                         chosen=None if 'chosen' not in selection else self.label(selection['chosen']),
                         excluded_origin_infeasible=[self.label(s) for s in selection.get('excluded', [])],
                         history=selection.get('history', [])),
            blocking_origin_infeasible_schemes=blocking, stuck_ranges=self.stuck_ranges(),
            lazy_obbt=self.lazy, errors=self.errors,
            inner_area_kw2=inner_area, outer_area_kw2=outer_area, outer_box_area_kw2=outer_box_area,
            metrics=None if evaluation is None else evaluation['metrics'])


# ---- 扫描缓存与逐格指标 -----------------------------------------------------------------------
def load_scan(path=SCAN_CACHE):
    with np.load(path, allow_pickle=False) as saved:
        scan = {key: saved[key].copy() for key in ('states', 'socp_states', 'origin', 'step', 'start')}
        scan['metadata'] = json.loads(str(saved['metadata']))
    shape = scan['states'].shape
    index = np.indices(shape).reshape(len(shape), -1).T
    scan.update(shape=shape, index=index, path=path,
                power=scan['origin']+(scan['start']+index+.5)*scan['step'],   # 格心，带符号 kW
                ac=scan['states'].ravel(), socp=scan['socp_states'].ravel())
    return scan


def metrics(inner, outer, ac, socp, areas):
    """FR=|域∧非可行|/|域|，MR=|可行∧非域|/|可行|（百分数，comparison_metrics 同定义）。"""
    ac, socp = ac == 1, socp == 1
    return dict(cells=int(len(inner)), ac_feasible_cells=int(ac.sum()), socp_feasible_cells=int(socp.sum()),
                inner_cells=int(inner.sum()), outer_cells=int(outer.sum()),
                inner_ac=comparison_metrics(inner, ac), inner_socp=comparison_metrics(inner, socp),
                outer_ac=comparison_metrics(outer, ac), outer_socp=comparison_metrics(outer, socp),
                outer_missed_ac_cells=int((ac & ~outer).sum()),
                sandwich_cells=int((outer & ~inner).sum()),
                inner_area_kw2=areas[0], outer_area_kw2=areas[1], outer_box_area_kw2=areas[2])


def evaluate(part, scan):
    labels = part.classify(scan['power'])
    cells = labels['cells']
    result = metrics(labels['inner'], labels['outer'], scan['ac'][cells], scan['socp'][cells], part.areas())
    missed = np.flatnonzero((scan['ac'][cells] == 1) & ~labels['outer'])
    violations = []
    for j in missed:
        k = int(labels['index'][j])
        interval = part.intervals[k] if k >= 0 else None
        halfspace = None if interval is None else interval.halfspace()
        violations.append(dict(p_kw=scan['power'][cells][j], grid_index=scan['index'][cells][j],
                               xi=labels['xi'][j], theta=labels['theta'][j],
                               interval=None if interval is None else interval.id,
                               c_xi=None if halfspace is None else float(labels['xi'][j]@halfspace[0]),
                               mu=None if halfspace is None else halfspace[1]))
    return dict(metrics=result, labels=labels, violations=violations)


def overall_metrics(evaluations, scan, parts):
    if not evaluations:
        return None
    inner = np.concatenate([e['labels']['inner'] for e in evaluations])
    outer = np.concatenate([e['labels']['outer'] for e in evaluations])
    ac = np.concatenate([scan['ac'][e['labels']['cells']] for e in evaluations])
    socp = np.concatenate([scan['socp'][e['labels']['cells']] for e in evaluations])
    return metrics(inner, outer, ac, socp, np.sum([part.areas() for part in parts], axis=0).tolist())


# ---- 图与终端输出 -----------------------------------------------------------------------------
SURFACE, INK, INK_2, MUTED, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#898781', '#e1e0d9'
AC_FILL, SOCP_FILL, INNER, COVERED, OUTER = '#cde2fb', '#f9d5c3', '#008300', '#1baf7a', '#e34948'
ZOOM_KW = (-1000., 800.)


def part_view(part):
    """绘图几何（xi）：取自本次运行的分区。covered 表示 x̂ 原点不可行、由近端覆盖证书认证。"""
    return dict(code=part.code, status=part.status, sign=part.sign, bounds=part.bounds,
                intervals=[dict(theta_a=i.theta_a, theta_b=i.theta_b, va=i.va, vb=i.vb, outer=i.outer_polygon(),
                                covered=i.cover is not None) for i in part.intervals])


def load_views(directory):
    """从既有结果目录的 intervals.json 重建绘图几何；outer_kw 是带符号 kW，换回 xi。"""
    data = json.loads((Path(directory)/'intervals.json').read_text(encoding='utf-8'))
    bounds = np.asarray(data['provenance']['bounds_kw'], float)
    views = []
    for code, row in data['partitions'].items():
        sign = np.asarray(row['sign'], float)
        intervals = [dict(theta_a=r['theta_a'], theta_b=r['theta_b'], va=np.asarray(r['va_xi']),
                          vb=np.asarray(r['vb_xi']), covered=r.get('cover_scheme') is not None,
                          outer=np.array([np.zeros(2), *(np.asarray(v)*sign/bounds for v in r['outer_kw'])]))
                     for r in row['intervals']]
        views.append(dict(code=code, status=row['status'], sign=sign, bounds=bounds, intervals=intervals))
    return views


def _plt():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 9, 'axes.edgecolor': '#c3c2b7',
                         'axes.labelcolor': INK_2, 'xtick.color': INK_2, 'ytick.color': INK_2})
    return plt


def _draw(ax, scan, views, xlim, ylim, title):
    """一个面板：扫描格（AC 可行 / 仅 SOCP 可行）、内三角、外域多边形、区间射线、扫描盒。"""
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Polygon
    shape, step = scan['shape'], scan['step']
    lower = scan['origin']+scan['start']*step
    edges = [lower[j]+np.arange(shape[j]+1)*step[j] for j in range(2)]
    layer = np.where(scan['states'] == 1, 1, np.where(scan['socp_states'] == 1, 2, 0))
    ax.set_facecolor(SURFACE)
    ax.pcolormesh(edges[0], edges[1], layer.T, cmap=ListedColormap([SURFACE, AC_FILL, SOCP_FILL]),
                  vmin=-.5, vmax=2.5, shading='flat', rasterized=True, zorder=0)
    ax.add_patch(Polygon([[edges[0][0], edges[1][0]], [edges[0][-1], edges[1][0]],
                          [edges[0][-1], edges[1][-1]], [edges[0][0], edges[1][-1]]],
                         closed=True, fill=False, edgecolor='#c3c2b7', lw=.6, zorder=1))
    for view in views:
        if not view['intervals']:
            continue
        scale = view['sign']*view['bounds']
        for interval in view['intervals']:
            ax.add_patch(Polygon(np.array([np.zeros(2), interval['va'], interval['vb']])*scale, closed=True,
                                 facecolor=COVERED if interval['covered'] else INNER,
                                 alpha=.45 if interval['covered'] else .30, edgecolor='none', zorder=2))
        outline = [np.zeros(2)]
        for interval in view['intervals']:
            outline.extend(interval['outer'][1:])
        outline = np.array([*outline, np.zeros(2)])*scale
        ax.plot(outline[:, 0], outline[:, 1], color=OUTER, lw=1.1, zorder=4, solid_joinstyle='miter')
        radii = {}
        for interval in view['intervals']:
            for theta, point in ((interval['theta_a'], interval['outer'][1]), (interval['theta_b'], interval['outer'][-1])):
                radii[theta] = max(radii.get(theta, 0.), float(np.linalg.norm(point)))
        for theta, radius in radii.items():
            end = radius*direction(theta)*scale
            ax.plot([0., end[0]], [0., end[1]], ls=(0, (3, 3)), lw=.35, color=MUTED, alpha=.6, zorder=3)
    ax.axhline(0., color='#c3c2b7', lw=.6, zorder=1)
    ax.axvline(0., color='#c3c2b7', lw=.6, zorder=1)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect('equal')
    ax.set_xlabel('p18 (kW, + load / - PV)')
    ax.set_ylabel('p25 (kW, + load / - PV)')
    ax.set_title(title, color=INK, fontsize=10, loc='left')
    ax.grid(color=GRID, lw=.5, zorder=1)
    ax.set_axisbelow(True)


def _legend(figure, covered):
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    handles = [Patch(facecolor=AC_FILL, label='AC-feasible cell'),
               Patch(facecolor=SOCP_FILL, label='SOCP-only feasible cell'),
               Patch(facecolor=INNER, alpha=.30, label='Inner triangle T_k = conv(0, va, vb)'),
               *([Patch(facecolor=COVERED, alpha=.45, label='Inner triangle, origin-infeasible x (near ends covered)')]
                 if covered else []),
               Line2D([], [], color=OUTER, lw=1.1, label='Outer polygon O_k = mu * T_k'),
               Line2D([], [], color=MUTED, lw=.8, ls=(0, (3, 3)), label='Interval rays'),
               Patch(facecolor='none', edgecolor='#c3c2b7', label='Scan box')]
    columns = len(handles) if len(handles) <= 6 else 4
    figure.legend(handles=handles, loc='lower center', ncol=columns, frameon=False, fontsize=9,
                  labelcolor=INK_2, bbox_to_anchor=(.5, .005))
    return -(-len(handles)//columns)   # 图例行数，供版面留白


def _status(views):
    return ', '.join(f'{v["code"]}: {v["status"]} ({len(v["intervals"])} int.)' for v in views)


def plot(path, scan, views, title):
    plt = _plt()
    figure, axes = plt.subplots(1, 3, figsize=(19., 6.9), gridspec_kw=dict(width_ratios=[.9, 1., 1.]))
    figure.patch.set_facecolor(SURFACE)
    lower, upper = scan['origin']+scan['start']*scan['step'], scan['origin']+(scan['start']+scan['shape'])*scan['step']
    box = 1.05*np.max([v['bounds'] for v in views], axis=0) if views else np.abs(lower)
    windows = [((lower[0], upper[0]), (lower[1], upper[1])), ((-box[0], box[0]), (-box[1], box[1])), (ZOOM_KW, ZOOM_KW)]
    titles = ['Paired AC/SOCP scan box (outer clipped)', 'Evaluation box +-bounds (outer extent)',
              'Zoom near the origin']
    for ax, (xlim, ylim), name in zip(axes, windows, titles):
        _draw(ax, scan, views, xlim, ylim, name)
    rows = _legend(figure, any(i['covered'] for v in views for i in v['intervals']))
    figure.suptitle(f'{title}    [{_status(views)}]', color=INK, fontsize=11, x=.01, ha='left')
    figure.tight_layout(rect=(0, .04+.04*rows, 1, .96))
    figure.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(figure)


def plot_comparison(path, scan, before, after, labels, title):
    """同一扫描盒内并排比较两次运行（左：既有结果，右：本次）。"""
    plt = _plt()
    figure, axes = plt.subplots(1, 2, figsize=(15., 7.6))
    figure.patch.set_facecolor(SURFACE)
    lower, upper = scan['origin']+scan['start']*scan['step'], scan['origin']+(scan['start']+scan['shape'])*scan['step']
    for ax, views, label in zip(axes, (before, after), labels):
        status = '\n'.join(_status(views[k:k+2]) for k in range(0, len(views), 2))
        _draw(ax, scan, views, (lower[0], upper[0]), (lower[1], upper[1]), f'{label}\n{status}')
    rows = _legend(figure, any(i['covered'] for v in (*before, *after) for i in v['intervals']))
    figure.suptitle(title, color=INK, fontsize=11, x=.01, ha='left')
    figure.tight_layout(rect=(0, .03+.035*rows, 1, .95))
    figure.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(figure)


PARTITION_COLORS = dict(pp='#2a78d6', pn='#eb6834', np='#1baf7a', nn='#eda100')   # 分类色前四槽，已过配色校验


def _series(rows, key, xkey):
    """收敛快照中的一条曲线；area_gap=外域∩盒/内域-1，certified_angle 换成百分数，非有限值跳过。"""
    points = []
    for row in rows:
        if key == 'area_gap':
            value = row['outer_box_area_kw2']/row['inner_area_kw2']-1. if row['inner_area_kw2'] > 0 else None
        elif key == 'certified_angle':
            value = 100.*row['certified_angle']
        else:
            value = row.get(key)
        if value is not None and np.isfinite(value):
            points.append((row[xkey], value))
    return [p[0] for p in points], [p[1] for p in points]


def _curve(ax, rows, key, xkey, code, log, end_label=False, **style):
    xs, ys = _series(rows, key, xkey)
    if log:
        xs, ys = [a for a, b in zip(xs, ys) if b > 0], [b for b in ys if b > 0]
    if not xs:
        return
    ax.plot(xs, ys, color=PARTITION_COLORS.get(code, INK_2), drawstyle='steps-post', **style)
    if end_label:   # 低对比度色槽的直接标注（文字用墨色）
        ax.annotate(code, (xs[-1], ys[-1]), xytext=(4, 0), textcoords='offset points', color=INK_2,
                    fontsize=8, va='center')


def _frame(ax, title, xlabel, log):
    ax.set_facecolor(SURFACE)
    if log:
        ax.set_yscale('log')
    ax.set_title(title, color=INK, fontsize=10, loc='left')
    ax.set_xlabel(xlabel)
    ax.grid(color=GRID, lw=.5)
    ax.set_axisbelow(True)


def plot_convergence(path, histories, tau, title):
    """单次运行的收敛过程（横轴为区间 MISOCP 次数）。"""
    plt = _plt()
    from matplotlib.lines import Line2D
    figure, axes = plt.subplots(2, 2, figsize=(13., 8.6))
    figure.patch.set_facecolor(SURFACE)
    xlabel = 'Interval MISOCP solves'
    panels = (('max_gap', 'Max radial gap bound  max_k (mu_k - 1)', True),
              ('area_gap', 'Area gap  (outer∩box - inner) / inner', True),
              ('certified_angle', 'Certified share of the angle range [0, pi/2] (%)', False),
              ('inner_ac', 'Inner region vs AC scan (%): MR solid, FR dotted (zeros not shown)', True))
    for ax, (key, name, log) in zip(axes.flat, panels):
        for code, rows in histories.items():
            if key == 'inner_ac':
                _curve(ax, rows, 'inner_ac_mr', 'step', code, True, end_label=True, lw=2)
                _curve(ax, rows, 'inner_ac_fr', 'step', code, True, lw=1.6, ls=(0, (1, 1.5)))
            else:
                _curve(ax, rows, key, 'step', code, log, end_label=True, lw=2)
        if key == 'max_gap':
            ax.axhline(tau, color=MUTED, lw=1)
            ax.annotate(f'tau = {tau:g}', (1., tau), xycoords=('axes fraction', 'data'), xytext=(-4, 4),
                        textcoords='offset points', ha='right', color=INK_2, fontsize=8)
        _frame(ax, name, xlabel, log)
    handles = [Line2D([], [], color=PARTITION_COLORS[code], lw=2, label=code) for code in histories]
    figure.legend(handles=handles, loc='lower center', ncol=len(handles), frameon=False, fontsize=9,
                  labelcolor=INK_2, bbox_to_anchor=(.5, .005))
    figure.suptitle(title, color=INK, fontsize=11, x=.01, ha='left')
    figure.tight_layout(rect=(0, .05, 1, .96))
    figure.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(figure)


def plot_convergence_comparison(path, current, previous, labels, tau, title):
    """两次运行的收敛过程对比（横轴为分区耗时，实线为本次，虚线为 --compare-with）。"""
    plt = _plt()
    from matplotlib.lines import Line2D
    figure, axes = plt.subplots(1, 3, figsize=(17., 5.6))
    figure.patch.set_facecolor(SURFACE)
    panels = (('max_gap', 'Max radial gap bound  max_k (mu_k - 1)', True),
              ('certified_angle', 'Certified share of [0, pi/2] (%)', False),
              ('inner_ac_mr', 'Inner region vs AC scan: MR (%, zeros not shown)', True))
    labelled = {'np', 'nn'}   # 低对比度色槽直接标注；pp/pn 数秒内结束，由图例区分
    for ax, (key, name, log) in zip(axes, panels):
        for code in current:
            if code in previous:
                _curve(ax, previous[code], key, 'seconds', code, log, lw=1.4, ls=(0, (4, 2)), alpha=.6)
            _curve(ax, current[code], key, 'seconds', code, log, end_label=code in labelled, lw=2)
        if key == 'max_gap':
            ax.axhline(tau, color=MUTED, lw=1)
            ax.annotate(f'tau = {tau:g}', (1., tau), xycoords=('axes fraction', 'data'), xytext=(-4, 4),
                        textcoords='offset points', ha='right', color=INK_2, fontsize=8)
        _frame(ax, name, 'Partition wall-clock time (s)', log)
    handles = [*(Line2D([], [], color=PARTITION_COLORS[code], lw=2, label=code) for code in current),
               Line2D([], [], color=INK_2, lw=2, label=labels[1]),
               Line2D([], [], color=INK_2, lw=1.4, ls=(0, (4, 2)), alpha=.6, label=labels[0])]
    figure.legend(handles=handles, loc='lower center', ncol=len(handles), frameon=False, fontsize=9,
                  labelcolor=INK_2, bbox_to_anchor=(.5, .005))
    figure.suptitle(title, color=INK, fontsize=11, x=.01, ha='left')
    figure.tight_layout(rect=(0, .08, 1, .95))
    figure.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(figure)


def fmt(value, digits=3):
    return '-' if value is None else f'{value:.{digits}f}'


def print_table(rows):
    header = (f'{"partition":<10s}{"in-AC FR%":>10s}{"in-AC MR%":>10s}{"in-SOCP FR%":>12s}{"in-SOCP MR%":>12s}'
              f'{"AC-out":>8s}{"band":>7s}{"inner kW2":>12s}{"outer kW2":>12s}{"outer^box":>12s}')
    print(header)
    print('-'*len(header))
    for name, m in rows:
        if m is None:
            continue
        print(f'{name:<10s}{fmt(m["inner_ac"]["fr_percent"]):>10s}{fmt(m["inner_ac"]["mr_percent"]):>10s}'
              f'{fmt(m["inner_socp"]["fr_percent"]):>12s}{fmt(m["inner_socp"]["mr_percent"]):>12s}'
              f'{m["outer_missed_ac_cells"]:>8d}{m["sandwich_cells"]:>7d}'
              f'{m["inner_area_kw2"]:>12.4g}{m["outer_area_kw2"]:>12.4g}{m["outer_box_area_kw2"]:>12.4g}')


def summary_line(result):
    counts = result['counts']
    schemes = ','.join(sorted(result['schemes_used']))
    volume = result.get('criterion') == 'volume'
    gap = f'volume gap={fmt(result["volume_gap_ratio"], 4)} (eps {result["epsilon"]:g}) ' if volume else ''
    leaves = f'bounded {result["bounded"]}' if volume else f'cert {result["certified"]}'
    return (f'{result["partition"]} {result["status"]:<10s} {gap}max g={fmt(result["max_g"], 5)} '
            f'intervals={result["intervals"]} ({leaves}, unres {result["unresolved"]}, '
            f'pend {result["pending"]}) MISOCP={counts["misocp"]} lazy OBBT={counts["lazy_obbt"]} '
            f'rays={counts["rays"]}+{counts["near_rays"]} near OBBT={counts["obbt"]} '
            f'covered={result["covered_intervals"]} schemes={{{schemes}}} time={result["seconds"]:.1f}s')


def git_commit():
    try:
        return subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=ROOT, capture_output=True, text=True,
                              check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def write_outputs(output, args, parts, results, overall, provenance, violations):
    intervals = {part.code: dict(sign=part.sign, status=part.status,
                                 intervals=[part.interval_json(i) for i in part.intervals], splits=part.splits,
                                 history=part.history)
                 for part in parts}
    (output/'intervals.json').write_text(json.dumps(plain(dict(provenance=provenance, partitions=intervals)),
                                                    ensure_ascii=False, indent=1), encoding='utf-8')
    summary = dict(provenance=provenance,
                   settings=dict(tau=args.tau, mip_seconds=args.mip_seconds, seconds=args.seconds,
                                 max_intervals=args.max_intervals, min_width=args.min_width, threads=args.threads,
                                 obbt_workers=args.obbt_workers, partitions=list(args.partitions),
                                 inner_cert=args.inner_cert, criterion=args.criterion,
                                 epsilon=args.epsilon if args.epsilon is not None else 2*args.tau,
                                 mip_gap=args.mip_gap, split_margin=SPLIT_MARGIN,
                                 cover_conv_tol=COVER_CONV_TOL, ray_seconds=RAY_SECONDS),
                   status=('certified' if parts and all(p.status == 'certified' for p in parts) else
                           'error' if any(p.status == 'error' for p in parts) else
                           'time_limit' if any(p.status == 'time_limit' for p in parts) else 'unresolved'),
                   partitions={r['partition']: r for r in results}, overall=overall,
                   outer_validity_violations=violations)
    summary = plain(summary)
    (output/'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=1), encoding='utf-8')
    return summary


def run_stats(summary, code):
    """运行统计；overall 为各分区之和（max g 取最大，状态取总体状态）。"""
    rows = list(summary['partitions'].values()) if code == 'overall' else [summary['partitions'][code]]
    finite = [p['max_g'] for p in rows if p['max_g'] is not None]
    return dict(status=summary['status'] if code == 'overall' else rows[0]['status'],
                max_g=max(finite) if len(finite) == len(rows) else None,
                intervals=sum(p['intervals'] for p in rows), certified=sum(p['certified'] for p in rows),
                misocp=sum(p['counts']['misocp'] for p in rows),
                rays=sum(p['counts']['rays']+p['counts'].get('near_rays', 0) for p in rows),
                lazy_obbt=sum(p['counts']['lazy_obbt'] for p in rows),
                covered=sum(p.get('covered_intervals', 0) for p in rows), seconds=sum(p['seconds'] for p in rows))


def grid_stats(metrics):
    return dict(inner_ac_fr=metrics['inner_ac']['fr_percent'], inner_ac_mr=metrics['inner_ac']['mr_percent'],
                inner_socp_fr=metrics['inner_socp']['fr_percent'], inner_socp_mr=metrics['inner_socp']['mr_percent'],
                ac_out=metrics['outer_missed_ac_cells'], band=metrics['sandwich_cells'],
                inner_area_kw2=metrics['inner_area_kw2'], outer_area_kw2=metrics['outer_area_kw2'],
                outer_box_area_kw2=metrics['outer_box_area_kw2'])


def compare(before, after):
    """同分区、同扫描网格上两次运行的逐项对比（before 为既有结果）。"""
    codes = [code for code in after['partitions'] if code in before['partitions']]
    same = set(before['partitions']) == set(after['partitions'])   # overall 只在分区集合相同时可比
    rows = {}
    for code in [*codes, *(['overall'] if same else [])]:
        old = before['overall'] if code == 'overall' else before['partitions'][code]['metrics']
        new = after['overall'] if code == 'overall' else after['partitions'][code]['metrics']
        rows[code] = dict(before=dict(**run_stats(before, code), **grid_stats(old)),
                          after=dict(**run_stats(after, code), **grid_stats(new)))
    return rows


def print_comparison(rows, labels):
    def pair(row, key, digits=3, kind='f'):
        values = []
        for side in ('before', 'after'):
            value = row[side][key]
            values.append('-' if value is None else f'{value:.{digits}{kind}}' if isinstance(value, float) else str(value))
        return ' -> '.join(values)
    print(f'\ncomparison: {labels[0]} -> {labels[1]}')
    columns = [('status', 0, 'f'), ('max_g', 4, 'f'), ('intervals', 0, 'f'), ('certified', 0, 'f'),
               ('misocp', 0, 'f'), ('rays', 0, 'f'), ('covered', 0, 'f'), ('seconds', 1, 'f'),
               ('inner_ac_fr', 3, 'f'), ('inner_ac_mr', 3, 'f'), ('inner_socp_mr', 3, 'f'), ('ac_out', 0, 'f'),
               ('band', 0, 'f'), ('inner_area_kw2', 4, 'g'), ('outer_box_area_kw2', 4, 'g')]
    for code, row in rows.items():
        print(f'  {code}: ' + '; '.join(f'{key} {pair(row, key, digits, kind)}' for key, digits, kind in columns))


def parse_args():
    parser = argparse.ArgumentParser(description='Radial cone partition with inner/outer sandwich on Case33 (18,25), mode 1')
    parser.add_argument('--partitions', default=','.join(PARTITIONS),
                        type=lambda value: tuple(code for code in value.split(',') if partition_sign(code)),
                        help='comma list of pp,pn,np,nn (order as load_nodes; p positive, n negative power)')
    parser.add_argument('--tau', type=float, default=.005)
    parser.add_argument('--mip-seconds', type=float, default=60.)
    parser.add_argument('--seconds', type=float, default=600., help='time limit per partition')
    parser.add_argument('--max-intervals', type=int, default=256)
    parser.add_argument('--min-width', type=float, default=1e-4, help='minimum angular width (rad)')
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--obbt-workers', type=int, default=8)
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument('--inner-cert', choices=('origin', 'near'), default='origin',
                        help='origin: x-hat must be origin-feasible; near: origin-infeasible x-hat allowed when '
                             'an origin-feasible scheme covers its near ends')
    parser.add_argument('--compare-with', type=Path,
                        help='existing output directory to compare with (writes comparison.json/.png)')
    parser.add_argument('--criterion', choices=('radial', 'volume'), default='radial',
                        help='radial: every interval mu_k <= 1+tau (default); volume: refine max Delta_k, '
                             'stop at sum Delta_k <= eps*sum area(T_k), as in the 3D script')
    parser.add_argument('--epsilon', type=float, help='volume-gap tolerance (default d*tau)')
    parser.add_argument('--mip-gap', type=float, default=1e-3, help='relative MIPGap of interval MISOCPs (volume mode)')
    args = parser.parse_args()
    if len(set(args.partitions)) != len(args.partitions):
        parser.error('duplicate partition codes')
    return args


def main():
    args = parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    started = perf_counter()
    network = Case33(load_nodes=LOAD_NODES, current_limit=CURRENT_LIMIT)
    budget = Case33.switch_budget
    scan = load_scan()
    metadata = scan['metadata']
    if network.fingerprint != metadata['network_fingerprint']:
        raise SystemExit(f'Network fingerprint {network.fingerprint} differs from scan cache '
                         f'{metadata["network_fingerprint"]}: {SCAN_CACHE}')
    if ac_identity(network, budget, MODE) != metadata['identity']:
        raise SystemExit(f'AC identity differs from scan cache {metadata["identity"]}: {SCAN_CACHE}')
    bounds = port_bounds(network)
    provenance = dict(script=Path(__file__).resolve().relative_to(ROOT).as_posix(), commit=git_commit(),
                      argv=sys.argv[1:], network=network.name, load_nodes=list(LOAD_NODES), mode=MODE,
                      budget=budget, current_limit_a=CURRENT_LIMIT, bounds_kw=bounds,
                      network_fingerprint=network.fingerprint, ac_identity=metadata['identity'],
                      scan_cache=SCAN_CACHE.relative_to(ROOT).as_posix(), scan_shape=list(scan['shape']),
                      switch_order=[c.id for c in network.corridors if c.switchable],
                      scheme_note='switch string over switch_order, 1 = closed; *_kw are signed true kW')
    print(f'Case33 {LOAD_NODES} mode={MODE} budget={budget} bounds={np.round(bounds, 3).tolist()} kW, '
          f'tau={args.tau}, partitions={",".join(args.partitions)}, inner certificate={args.inner_cert}', flush=True)
    before = None
    if args.compare_with is not None:
        before = json.loads((args.compare_with/'summary.json').read_text(encoding='utf-8'))
        if before['provenance']['scan_cache'] != provenance['scan_cache']:
            raise SystemExit(f'{args.compare_with} uses another scan cache: {before["provenance"]["scan_cache"]}')
    parts, results, evaluations, violations, summary = [], [], [], [], None
    with threadpool_limits(limits=1):
        for code in args.partitions:
            print(f'== partition {code} sign={partition_sign(code)}', flush=True)
            part = RadialSandwich(network, code, bounds, budget, args).run(scan)
            evaluation = evaluate(part, scan)
            parts.append(part)
            evaluations.append(evaluation)
            results.append(part.result(evaluation))
            print(summary_line(results[-1]), flush=True)
            if evaluation['violations']:
                violations = [dict(partition=code, **row) for row in evaluation['violations']]
                print(f'!! OUTER VALIDITY FAILED in {code}: {len(violations)} AC-feasible cells outside the outer '
                      f'domain', flush=True)
                for row in violations:
                    print(f'   p={np.asarray(row["p_kw"]).tolist()} kW grid={np.asarray(row["grid_index"]).tolist()} '
                          f'interval={row["interval"]} c.xi={row["c_xi"]} mu={row["mu"]}', flush=True)
            overall = overall_metrics(evaluations, scan, parts)
            provenance['seconds'] = perf_counter()-started
            summary = write_outputs(output, args, parts, results, overall, provenance, violations)
            if violations:
                break
    views = [part_view(part) for part in parts]
    plot(output/'radial.png', scan, views,
         f'Case33 (18, 25), mode 1 - radial sandwich, tau = {args.tau:g}, inner certificate: {args.inner_cert}')
    histories = {part.code: part.history for part in parts}
    plot_convergence(output/'convergence.png', histories, args.tau,
                     f'Case33 (18, 25), mode 1 - convergence, tau = {args.tau:g}, inner certificate: {args.inner_cert}')
    print()
    print_table([*((r['partition'], r['metrics']) for r in results), ('overall', overall)])
    if before is not None and summary is not None:
        labels = (f'{args.compare_with.as_posix()} (inner certificate: {before["settings"].get("inner_cert", "origin")})',
                  f'{args.output.as_posix()} (inner certificate: {args.inner_cert})')
        rows = compare(before, summary)
        print_comparison(rows, labels)
        (output/'comparison.json').write_text(json.dumps(plain(dict(before=labels[0], after=labels[1], rows=rows)),
                                                         ensure_ascii=False, indent=1), encoding='utf-8')
        plot_comparison(output/'comparison.png', scan, load_views(args.compare_with), views, labels,
                        f'Case33 (18, 25), mode 1 - radial sandwich, tau = {args.tau:g}: inner certificate comparison')
        previous = {code: row['history'] for code, row in json.loads(
            (args.compare_with/'intervals.json').read_text(encoding='utf-8'))['partitions'].items() if row.get('history')}
        if previous:
            plot_convergence_comparison(
                output/'convergence_comparison.png', histories, previous,
                (f'{before["settings"].get("inner_cert", "origin")} certificate ({args.compare_with.name})',
                 f'{args.inner_cert} certificate ({args.output.name})'), args.tau,
                f'Case33 (18, 25), mode 1 - convergence comparison, tau = {args.tau:g}')
        else:
            print(f'{args.compare_with} has no convergence history; rerun it with this script to compare curves')
    print(f'\noutputs: {output}  total {perf_counter()-started:.1f}s', flush=True)
    return 2 if violations else 1 if any(p.status == 'error' for p in parts) else 0


if __name__ == '__main__':
    raise SystemExit(main())
