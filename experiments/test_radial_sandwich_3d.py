"""径向锥剖分 + 内外夹逼的三维推广：Case33(18,25,30)、mode=1、d=3，与配对 AC/SOCP 扫描逐格比较。

在项目根目录运行::

    python -X utf8 experiments/test_radial_sandwich_3d.py --output results/radial/case33_18_25_30_near

二维版见 experiments/test_radial_sandwich.py。本脚本复用其紧化、原点认证、射线与反向近端射线、x 自由
MISOCP 与懒惰 OBBT（RadialSandwich），只替换方向剖分、几何与评价；原理、证书与实现约定同二维版，
下面只写三维不同之处。

方向剖分
----
xi 空间第一卦限的方向用球面三角形剖分为单纯锥 K = cone(u1, u2, u3)，u_i 为登记在方向表中的单位方向。
内域：v_i = rho_x̂(u_i)·u_i 为方案 x̂ 的紧化射线顶点，T = conv(0, v1, v2, v3)。x̂ 原点可行时 T ⊆ R_x̂；
近端覆盖模式（默认）下 x̂ 原点不可行时，近端 n_i 由远端朝原点的反向射线求得，conv(n1..n3, v1..v3) ⊆ R_x̂，
若原点可行的 y 满足 rho_y(u_i) >= |n_i|（i=1,2,3），则 conv(0, n1, n2, n3) ⊆ conv(0, y1, y2, y3) ⊆ R_y。
平面 n1n2n3 把 T 切成这两块，故 T ⊆ R_x̂ ∪ R_y ⊆ R。
外域：Q=[v1 v2 v3]，c=1ᵀQ⁻¹；x 自由 MISOCP 在 K 内（Q⁻¹xi >= 0，各行归一化）最大化 c·xi，幅值上界
bounds，BestBdStop=BestObjStop=1+tau，mu=ObjBound，O = {xi∈K: c·xi<=mu} = mu·T ⊇ R∩K，间隙 g=mu-1。
细分：解点方向 u* 的锥坐标 λ（λ>=0，Σλ=1）全部不小于 δ=SPLIT_MARGIN 时，在 u* 星形剖分为 3 个子锥；恰有一个
分量小于 δ 时，在对边上按 λ 投影二分（投影过近端点时取棱中点）；否则（无解或靠近顶点）在最长棱中点二分。
只做内部星形剖分时原棱永不细分，棱方向上的间隙无法收敛，故必须有棱二分；棱中点方向由相邻锥共用。相邻锥
可以不协调（悬挂顶点），每个锥的内外证书独立成立，不影响有效性。子锥 x̂ 取 {父 x̂, x*, 父覆盖方案 y} 中
满足证书且三顶点半径乘积（与四面体体积成正比）最大者。锥数达 --max-cones 或角直径小于 --min-width 时记为
unresolved。

停止准则（--criterion）
----
volume（默认）：每个锥的体积缺口 Δ_k=(μ̄_k^d-1)·vol(T_k)，vol(T_k)=|det[v1 v2 v3]|/d!，同时是本锥内
“漏掉的可行域 R∖I”与“多算的区域 O∖R”的体积上界。每次细分 Δ_k 最大的锥，子锥立即求外界；当
ΣΔ_k <= ε·Σvol(T_k) 时停止，ε 默认 d·tau（μ̄^d-1≈d(μ̄-1)，与径向准则对比公平）。Δ 的排序与判停都需要较紧的
μ̄，故锥 MISOCP 解到相对间隙 --mip-gap，不再用 1+tau 提前停止。折角附近的窄锥体积小，μ̄ 偏大也不必一直细分；
这些方向上不再保证径向误差 <= tau，但外界仍有效。最大径向间隙 max(μ̄_k-1) 作为诊断量照常报告。
radial：原规则，每个锥 μ̄_k <= 1+tau 才认证，待处理锥按继承的父 g 从大到小处理。

评价
----
配对 AC/SOCP 扫描取项目标准的 80³ 网格（vertify.scan_ac_reference，reference_box 给出的 SOCP 坐标全局
界），缓存于 results/scan/case33bw/18_25_30/<identity>/。格心按符号归属分区，按方向找所在锥（锥坐标最小
分量最大者），α=Q⁻¹xi：内域 α>=0 且 Σα<=1，外域 c·xi<=mu。指标与二维相同（内域相对 AC/SOCP 的 FR/MR、
AC 可行∧非外域、夹逼带），面积换为体积（kW³，四面体公式；外域另给与评价盒之交）。每次锥决策后记录收敛
快照（逐格标签增量更新）。

输出（--output 下）
----
cones.json（各分区叶锥、细分历史、收敛快照）、summary.json、slices.png（p30 切片上的逐格分类）、
convergence.png。本脚本不修改主线文件，结果不加入 results/manifest.json。
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import heapq
from itertools import product
import json
from pathlib import Path
import sys
from time import perf_counter
import traceback

import gurobipy as gp
from gurobipy import GRB
import numpy as np
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
for _path in (str(ROOT), str(HERE)):
    if _path not in sys.path:
        sys.path.insert(0, _path)

from Network.case33bw import Case33, CURRENT_LIMIT
from model import port_bounds
from monitor import comparison_metrics
from region import clip_polytope, polytope_volume
from test_radial_sandwich import (ALPHA_TOL, COVER_CONV_TOL, INK, INK_2, MUTED, RAY_SECONDS, SPLIT_MARGIN, SURFACE,
                                  PartitionTimeout, RadialSandwich, _frame, _plt, fmt, git_commit, partition_sign,
                                  plain)
from vertify import ac_identity, ac_network, load_scan, scan_path

LOAD_NODES = (18, 25, 30)
MODE = 1
PARTITIONS = tuple(''.join(code) for code in product('pn', repeat=3))
DEFAULT_OUTPUT = ROOT/'results'/'radial'/'case33_18_25_30_near'
OCTANT = np.pi/2   # 第一卦限的立体角（sr）
EDGES = ((0, 1), (0, 2), (1, 2))


@dataclass(eq=False)
class Cone:
    """方向锥 K=cone(u1,u2,u3)：内域 T=conv(0,v1,v2,v3)，外域 {xi∈K: c·xi<=mu}；verts 第 i 行为 v_i（xi）。"""
    id: int
    keys: tuple                         # 三个方向编号
    x: tuple                            # 内域方案 x̂
    verts: np.ndarray                   # (3,3)
    parent: int | None = None
    depth: int = 0
    priority: float = np.inf            # 父锥的 g；根锥为 inf
    inherited: tuple | None = None      # 父锥外界 (c, mu)；本锥尚无有限上界时沿用
    mu: float = np.inf
    g: float = np.inf
    status: str = 'pending'
    reason: str | None = None
    misocp: dict | None = None
    cover: tuple | None = None          # 近端覆盖方案 y；x̂ 原点可行时为 None
    nears: np.ndarray | None = None     # (3,3) 近端（xi），原点可行时为 0
    cover_radii: np.ndarray | None = None
    volumes: tuple | None = None        # 缓存的 (内域, 外域, 外域∩盒) 体积（xi³）
    diam: float = 0.                    # 缓存的角直径（rad）与立体角（sr）
    solid: float = 0.
    gap: float = np.inf                 # 缓存的径向间隙上界（leaf_gap）
    delta: float = np.inf               # 缓存的体积缺口 Δ=(μ̄^d-1)·vol(T)（xi³）
    final: dict | None = None           # 最后一轮 MISOCP 记录（细分用解点与方案）

    @property
    def Q(self):
        return self.verts.T

    @property
    def c(self):
        """c=1ᵀQ⁻¹：锥内 xi=Σα_i v_i 时 Σα_i=c·xi。"""
        return np.linalg.solve(self.Q.T, np.ones(3))

    def halfspace(self):
        return (self.c, self.mu) if np.isfinite(self.mu) else self.inherited


class RadialSandwich3D(RadialSandwich):
    """三维符号分区的径向夹逼；方向键为方向表编号，其余求解器调用沿用二维类。"""
    CENTER_CONTEXT = 'radial center'

    def __init__(self, network, code, bounds, budget, args):
        super().__init__(network, code, bounds, budget, args)
        self.criterion = args.criterion
        self.epsilon = args.epsilon if args.epsilon is not None else self.dim*args.tau
        self.directions = [np.eye(3)[j] for j in range(3)]
        self.midpoints = {}
        self.cones, self.cell_xi = [], None

    # ---- 方向表 --------------------------------------------------------------------------
    def unit(self, key):
        return self.directions[int(key)]

    def new_direction(self, vector):
        u = np.asarray(vector, float)
        u = u/np.linalg.norm(u)
        u[np.abs(u) < 1e-15] = 0.
        self.directions.append(u)
        return len(self.directions)-1

    def midpoint(self, a, b):
        """棱中点方向，相邻锥共用同一编号（射线缓存一致）。"""
        key = min(a, b), max(a, b)
        if key not in self.midpoints:
            self.midpoints[key] = self.new_direction(self.directions[a]+self.directions[b])
        return self.midpoints[key]

    def U(self, keys):
        return np.column_stack([self.directions[k] for k in keys])

    def diameter(self, keys):
        U = self.U(keys)
        return max(float(np.arccos(np.clip(U[:, a]@U[:, b], -1., 1.))) for a, b in EDGES)

    def solid_angle(self, keys):
        """球面三角形立体角（Van Oosterom–Strackee）。"""
        a, b, c = (self.directions[k] for k in keys)
        return float(2.*np.arctan2(abs(a@np.cross(b, c)), 1.+a@b+b@c+c@a))

    # ---- 初始方案与内域证书 ----------------------------------------------------------------
    def center_ray(self, problem, xi):
        """中心方向 (1,1,1)/√3：xi1=xi2=xi3，max Σxi。"""
        for j in range(1, self.dim):
            problem.model.addConstr(xi[0]-xi[j] == 0., name=f'radial_direction_{j}')
        problem.model.setObjective(gp.quicksum(xi), GRB.MAXIMIZE)

    def inner_vertices3(self, x, keys, covers=()):
        """方案 x 在锥 keys 上的内顶点及证书，不满足时为 None（二维 inner_vertices 的三顶点版）。"""
        feasible = self.origin_feasible(x)
        if not (feasible or self.near_mode):
            return None
        verts = [self.vertex(x, k) for k in keys]
        if any(v is None for v in verts):
            return None
        verts = np.array(verts)
        if feasible:
            return dict(x=tuple(x), verts=verts, cover=None, nears=np.zeros((3, 3)), cover_radii=None)
        nears = [self.near(x, k) for k in keys]
        if any(n is None for n in nears):
            return None
        nears = np.array(nears)
        radii = np.linalg.norm(nears, axis=1)
        order = [*covers, *(s for s, row in self.origin.items() if row['feasible'])]
        for y in dict.fromkeys(tuple(s) for s in order if s is not None):
            if not self.origin_feasible(y):
                continue
            reach = [self.vertex(y, k) for k in keys]
            if any(v is None for v in reach):
                continue
            reach = np.linalg.norm(np.array(reach), axis=1)
            if np.all(reach >= radii):
                return dict(x=tuple(x), verts=verts, cover=y, nears=nears, cover_radii=reach)
        return None

    def new_cone(self, keys, inner, parent=None):
        cone = Cone(self.next_id, tuple(keys), inner['x'], inner['verts'], cover=inner['cover'],
                    nears=inner['nears'], cover_radii=inner['cover_radii'])
        self.next_id += 1
        if parent is not None:
            cone.parent, cone.depth, cone.priority = parent.id, parent.depth+1, parent.g
            if np.isfinite(parent.mu):
                cone.inherited = parent.c, parent.mu
        cone.diam, cone.solid = self.diameter(cone.keys), self.solid_angle(cone.keys)
        self.cache_volumes(cone)
        return cone

    # ---- 外域 MISOCP 与细分 -----------------------------------------------------------------
    def certify(self, cone):
        Q, c = cone.Q, cone.c
        rows = np.linalg.inv(Q)
        rows = rows/np.linalg.norm(rows, axis=1, keepdims=True)

        volume = self.criterion == 'volume'

        def configure(problem, xi):
            model = problem.model
            for j, row in enumerate(rows):   # 锥约束 Q⁻¹xi >= 0
                model.addConstr(gp.quicksum(float(row[k])*xi[k] for k in range(3)) >= 0., name=f'cone_{j}')
            problem.power.UB = self.bounds
            model.setObjective(gp.quicksum(float(c[k])*xi[k] for k in range(3)), GRB.MAXIMIZE)
            if volume:   # Δ 排序与停止判据都要较紧的 μ̄：解到相对间隙 mip_gap，不提前停
                model.Params.MIPGap = self.args.mip_gap
            else:
                model.Params.BestBdStop = model.Params.BestObjStop = 1.+self.tau
        self.counts['interval_solves'] += 1
        rounds, triggers = self.misocp(configure, stop=None if volume else 1.+self.tau, context=f'cone {cone.id}')
        final = rounds[-1]
        cone.mu = 1. if final['status'] == 'INFEASIBLE' else final['bound']
        cone.g = cone.mu-1.
        cone.misocp = dict(status=final['status'], objective=final['objective'], bound=final['bound'],
                           seconds=sum(r['build_seconds']+r['solve_seconds'] for r in rounds),
                           lazy_obbt=triggers, rounds=rounds)
        self.cache_volumes(cone)
        return final

    def split_options(self, cone, final):
        """按优先次序给出候选剖分 [(子锥方向组, 规则)]：星形 / 棱投影 / 棱中点，最后总有最长棱中点二分。"""
        keys, U = cone.keys, self.U(cone.keys)
        options = []
        if final['point'] is not None:
            point = final['point']/self.bounds
            lam = np.maximum(np.linalg.solve(U, point), 0.)
            if lam.sum() > 0.:
                lam = lam/lam.sum()
                small = lam < SPLIT_MARGIN
                if not small.any():
                    m = self.new_direction(point)
                    options.append(([(m, keys[1], keys[2]), (keys[0], m, keys[2]), (keys[0], keys[1], m)],
                                    'stellar'))
                elif small.sum() == 1:
                    a = int(np.flatnonzero(small)[0])
                    b, c = (t for t in range(3) if t != a)
                    share = lam[b]/(lam[b]+lam[c])
                    if SPLIT_MARGIN < share < 1.-SPLIT_MARGIN:
                        m, rule = self.new_direction(lam[b]*U[:, b]+lam[c]*U[:, c]), 'edge'
                    else:
                        m, rule = self.midpoint(keys[b], keys[c]), 'edge-mid'
                    options.append(([(keys[a], keys[b], m), (keys[a], m, keys[c])], rule))
        angles = [float(np.arccos(np.clip(U[:, a]@U[:, b], -1., 1.))) for a, b in EDGES]
        a, b = EDGES[int(np.argmax(angles))]
        r = 3-a-b
        m = self.midpoint(keys[a], keys[b])
        options.append(([(keys[r], keys[a], m), (keys[r], m, keys[b])], 'longest-edge'))
        return options

    def split(self, cone, final):
        candidates = [cone.x]
        if final['x'] is not None and final['x'] != cone.x:
            candidates.append(final['x'])
        if self.near_mode and cone.cover is not None and cone.cover not in candidates:
            candidates.append(cone.cover)
        covers = (cone.cover, cone.x)
        for groups, rule in self.split_options(cone, final):
            children = []
            for keys in groups:
                best = None
                for scheme in candidates:
                    inner = self.inner_vertices3(scheme, keys, covers)
                    if inner is None:
                        continue
                    score = float(np.prod(np.linalg.norm(inner['verts'], axis=1)))
                    if best is None or score > best[0]:
                        best = score, inner
                if best is None:
                    break
                children.append(self.new_cone(keys, best[1], parent=cone))
            if len(children) == len(groups):
                incumbent = None if final['x'] is None else self.label(final['x'])
                self.splits.append(dict(
                    id=cone.id, keys=cone.keys, scheme=self.label(cone.x), mu=cone.mu, g=cone.g,
                    misocp_status=cone.misocp['status'], misocp_seconds=cone.misocp['seconds'],
                    incumbent=cone.misocp['objective'], incumbent_scheme=incumbent,
                    incumbent_origin_feasible=(None if final['x'] is None else
                                               self.origin.get(final['x'], {}).get('feasible')),
                    rule=rule, children=[child.id for child in children],
                    child_schemes=[self.label(child.x) for child in children],
                    child_covers=[None if child.cover is None else self.label(child.cover) for child in children]))
                return children
        return None

    # ---- 体积、间隙与逐格标签 ----------------------------------------------------------------
    def outer_polytope(self, cone):
        """外域多面体（xi）：有外界时为 conv(0, t_i u_i)，t_i=mu/(c·u_i)；否则为锥与评价盒之交。"""
        U = self.U(cone.keys)
        halfspace = cone.halfspace()
        if halfspace is None:
            poly = np.vstack([np.zeros(3), (U*10.).T])
            for j in range(3):
                poly = clip_polytope(poly, 1., -np.eye(3)[j])
            return poly
        c, mu = halfspace
        return np.vstack([np.zeros(3), (U*(mu/(c@U))).T])

    def cache_volumes(self, cone):
        """锥的体积与间隙随外界变化（创建、求解）时缓存，快照只做求和。"""
        outer = self.outer_polytope(cone)
        box = outer
        for j in range(3):
            box = clip_polytope(box, 1., -np.eye(3)[j])
        cone.volumes = (abs(float(np.linalg.det(cone.Q)))/6., polytope_volume(outer), polytope_volume(box))
        cone.gap = self.leaf_gap(cone)
        # 体积缺口 Δ=(μ̄^d-1)·vol(T)：μ̄=1+gap 为本锥（或沿用父锥）外界相对 T 的因子
        cone.delta = max(((1.+cone.gap)**self.dim-1.)*cone.volumes[0], 0.) if np.isfinite(cone.gap) else np.inf

    def volumes(self):
        """内域、外域（四面体公式）及外域∩评价盒的体积（kW³）。"""
        scale = float(np.prod(self.bounds))
        if not self.cones:
            return 0., scale, scale
        return tuple(sum(cone.volumes[j] for cone in self.cones)*scale for j in range(3))

    def leaf_gap(self, cone):
        if np.isfinite(cone.mu):
            return cone.mu-1.
        if cone.inherited is None:
            return np.inf
        c_p, mu_p = cone.inherited
        U = self.U(cone.keys)
        return float(np.max(mu_p*(cone.c@U)/(c_p@U)))-1.   # 线性分式在锥的生成方向上取极值

    def attach_cells(self, scan):
        """本分区的扫描格心（xi）及其 AC/SOCP 标签；内外域标签随剖分增量更新。"""
        cells = np.all(np.where(scan['power'] >= 0., 1, -1) == self.sign, axis=1)
        self.cell_index = np.flatnonzero(cells)
        self.cell_xi = np.abs(scan['power'][cells])/self.bounds
        self.cell_ac, self.cell_socp = scan['ac'][cells] == 1, scan['socp'][cells] == 1
        self.cell_inner = np.zeros(len(self.cell_xi), bool)
        self.cell_outer = np.ones(len(self.cell_xi), bool)
        self.cell_owner = np.full(len(self.cell_xi), -1)
        self.cone_cells = {}

    def label_cells(self, cone):
        if self.cell_xi is None:
            return
        idx = self.cone_cells.get(cone.id, np.empty(0, int))
        xi = self.cell_xi[idx]
        alpha = xi@np.linalg.inv(cone.Q).T
        self.cell_inner[idx] = np.all(alpha >= -ALPHA_TOL, axis=1) & (alpha.sum(axis=1) <= 1.)
        halfspace = cone.halfspace()
        self.cell_outer[idx] = (xi@halfspace[0] <= halfspace[1]) if halfspace is not None else np.all(xi <= 1., axis=1)

    def assign_cells(self, idx, cones):
        """把格心分给所在锥：锥坐标最小分量最大者（边界上任取一侧）。"""
        if self.cell_xi is None:
            return
        best, owner = np.full(len(idx), -np.inf), np.zeros(len(idx), int)
        for k, cone in enumerate(cones):
            score = (self.cell_xi[idx]@np.linalg.inv(self.U(cone.keys)).T).min(axis=1)
            better = score > best
            best[better], owner[better] = score[better], k
        for k, cone in enumerate(cones):
            self.cone_cells[cone.id] = idx[owner == k]
            self.cell_owner[self.cone_cells[cone.id]] = cone.id
            self.label_cells(cone)

    def volume_ratio(self):
        """体积缺口比 ΣΔ_k / Σvol(T_k)；任一叶锥尚无外界时为 inf。"""
        total = sum(c.volumes[0] for c in self.cones)
        return sum(c.delta for c in self.cones)/total if total > 0. else np.inf

    def snapshot(self, cone, decision, scan=None):
        inner, _, outer_box = self.volumes()
        row = dict(step=self.counts['interval_solves'], seconds=perf_counter()-self.started,
                   cone=None if cone is None else cone.id, decision=decision,
                   mu=None if cone is None else cone.mu, leaves=len(self.cones),
                   certified=sum(c.status == 'certified' for c in self.cones),
                   certified_angle=sum(c.solid for c in self.cones if c.status == 'certified')/OCTANT,
                   radial_certified_angle=sum(c.solid for c in self.cones if c.gap <= self.tau)/OCTANT,
                   volume_gap_ratio=self.volume_ratio(),
                   max_gap=max(c.gap for c in self.cones), inner_volume_kw3=inner,
                   outer_box_volume_kw3=outer_box, misocp=self.counts['misocp'],
                   rays=self.counts['rays']+self.counts['near_rays'])
        if self.cell_xi is not None:
            accuracy = comparison_metrics(self.cell_inner, self.cell_ac)
            row.update(inner_ac_mr=accuracy['mr_percent'], inner_ac_fr=accuracy['fr_percent'],
                       ac_out=int((self.cell_ac & ~self.cell_outer).sum()))
        self.history.append(row)

    # ---- 主循环 -----------------------------------------------------------------------------
    def make_root(self):
        """中心方向径向 MISOCP 选初始方案，建根锥 cone(e1,e2,e3) 并分配格心。"""
        x0 = self.initial_scheme()
        candidates = [self.selection['x0'], x0] if self.near_mode and self.selection['x0'] != x0 else [x0]
        best = None
        for scheme in candidates:
            inner = self.inner_vertices3(scheme, (0, 1, 2), (x0,))
            if inner is not None:
                score = float(np.prod(np.linalg.norm(inner['verts'], axis=1)))
                if best is None or score > best[0]:
                    best = score, inner
        if best is None:
            raise RuntimeError(f'{self.code}: initial scheme {self.label(x0)} has no ray vertex on an axis')
        root = self.new_cone((0, 1, 2), best[1])
        self.cones = [root]
        if self.cell_xi is not None:
            self.assign_cells(np.arange(len(self.cell_xi)), [root])
        self.snapshot(None, 'initial')
        return root

    def solve_cone(self, cone):
        final = self.certify(cone)
        cone.final = final
        self.label_cells(cone)
        info = cone.misocp
        print(f'  [{self.code}] #{cone.id:<4d} diam={cone.diam:.4f} xhat={self.label(cone.x)}'
              f'{"" if cone.cover is None else "<-"+self.label(cone.cover)} {info["status"]:<14s} '
              f'mu={cone.mu:.6g} lazy={info["lazy_obbt"]} {info["seconds"]:.2f}s leaves={len(self.cones)}', flush=True)
        return final

    def run(self, scan=None):
        if self.criterion == 'volume':
            return self.run_volume()
        try:
            root = self.make_root()
            heap = [(-root.priority, -root.diam, root.id, root)]
            while heap:
                self.remaining()
                *_, cone = heapq.heappop(heap)
                final = self.solve_cone(cone)
                children = None
                if cone.mu <= 1.+self.tau:
                    cone.status = 'certified'
                elif cone.diam < self.args.min_width:
                    cone.status, cone.reason = 'unresolved', 'min_width'
                elif len(self.cones) >= self.args.max_cones:
                    cone.status, cone.reason = 'unresolved', 'max_cones'
                elif (children := self.split(cone, final)) is None:
                    cone.status, cone.reason = 'unresolved', 'no_inner_vertex'
                else:
                    k = next(j for j, item in enumerate(self.cones) if item is cone)
                    self.cones[k:k+1] = children
                    if self.cell_xi is not None:
                        self.assign_cells(self.cone_cells.pop(cone.id), children)
                    for child in children:
                        heapq.heappush(heap, (-child.priority, -child.diam, child.id, child))
                self.snapshot(cone, 'split' if children else cone.reason or cone.status)
            self.status = 'certified' if all(c.status == 'certified' for c in self.cones) else 'unresolved'
        except PartitionTimeout:
            self.status = 'time_limit'
        except Exception as exc:
            self.status = 'error'
            self.errors.append(dict(kind='fatal', error=repr(exc), traceback=traceback.format_exc()))
            print(traceback.format_exc(), flush=True)
        self.seconds = perf_counter()-self.started
        return self

    def run_volume(self):
        """体积缺口准则：每次细分 Δ_k 最大的锥，子锥立即求外界；ΣΔ_k <= ε·Σvol(T_k) 时停止。
        叶锥状态为 bounded（已有外界）/ pending（时限到时尚未求解，沿用父锥外界）/ unresolved（不能再分）。"""
        try:
            root = self.make_root()
            self.solve_cone(root)
            root.status = 'bounded'
            self.snapshot(root, 'bounded')
            heap = [(-root.delta, root.id, root)]
            while True:
                if self.volume_ratio() <= self.epsilon:
                    self.status = 'certified'
                    break
                if not heap:
                    self.status = 'unresolved'
                    break
                self.remaining()
                _, _, cone = heapq.heappop(heap)
                if cone.diam < self.args.min_width:
                    cone.status, cone.reason = 'unresolved', 'min_width'
                    continue
                if len(self.cones) >= self.args.max_cones:
                    cone.status, cone.reason = 'unresolved', 'max_cones'
                    continue
                children = self.split(cone, cone.final)
                if children is None:
                    cone.status, cone.reason = 'unresolved', 'no_inner_vertex'
                    continue
                k = next(j for j, item in enumerate(self.cones) if item is cone)
                self.cones[k:k+1] = children
                if self.cell_xi is not None:
                    self.assign_cells(self.cone_cells.pop(cone.id), children)
                for child in children:
                    self.solve_cone(child)
                    child.status = 'bounded'
                    heapq.heappush(heap, (-child.delta, child.id, child))
                self.snapshot(cone, 'split')
        except PartitionTimeout:
            self.status = 'time_limit'
        except Exception as exc:
            self.status = 'error'
            self.errors.append(dict(kind='fatal', error=repr(exc), traceback=traceback.format_exc()))
            print(traceback.format_exc(), flush=True)
        self.seconds = perf_counter()-self.started
        return self

    # ---- 导出 -------------------------------------------------------------------------------
    def round_json(self, record):
        point = record['point']
        scaled = None if point is None else point/self.bounds
        direction = None if scaled is None or not np.any(scaled) else scaled/np.linalg.norm(scaled)
        return dict(status=record['status'], build_seconds=record['build_seconds'],
                    solve_seconds=record['solve_seconds'], tightened_schemes=record['tightened_schemes'],
                    objective=record['objective'], bound=record['bound'],
                    scheme=None if record['x'] is None else self.label(record['x']),
                    point_kw=None if point is None else self.sign*point, direction=direction,
                    lazy_obbt_seconds=record.get('lazy_obbt_seconds'))

    def cone_json(self, cone):
        halfspace = cone.halfspace()
        incumbent, feasible = self.incumbent(cone)
        info = cone.misocp
        outer = self.outer_polytope(cone)
        return dict(
            id=cone.id, parent=cone.parent, depth=cone.depth, keys=cone.keys, directions=self.U(cone.keys).T,
            diameter=self.diameter(cone.keys), solid_angle=self.solid_angle(cone.keys),
            scheme=self.label(cone.x), status=cone.status, reason=cone.reason,
            verts_xi=cone.verts, verts_kw=[self.kw(v) for v in cone.verts], rho=np.linalg.norm(cone.verts, axis=1),
            origin_feasible=cone.cover is None, cover_scheme=None if cone.cover is None else self.label(cone.cover),
            nears_xi=cone.nears, cover_radii=cone.cover_radii, mu=cone.mu, g=cone.g, priority=cone.priority,
            outer_source='box' if halfspace is None else ('own' if np.isfinite(cone.mu) else 'parent'),
            outer_kw=[self.kw(v) for v in outer[1:]], incumbent_scheme=incumbent, incumbent_origin_feasible=feasible,
            misocp=None if info is None else dict(status=info['status'], seconds=info['seconds'],
                                                  objective=info['objective'], bound=info['bound'],
                                                  lazy_obbt=info['lazy_obbt'],
                                                  rounds=[self.round_json(r) for r in info['rounds']]))

    def result(self, metrics=None):
        solved = [c for c in self.cones if c.misocp is not None]
        statuses = [c.status for c in self.cones]
        used, pairs = {}, {}
        for cone in self.cones:
            row = used.setdefault(self.label(cone.x), dict(cones=0, solid_angle_share=0.))
            row['cones'] += 1
            row['solid_angle_share'] += self.solid_angle(cone.keys)/OCTANT
            if cone.cover is not None:
                row = pairs.setdefault(f'{self.label(cone.x)}<-{self.label(cone.cover)}', dict(cones=0, solid_angle_share=0.))
                row['cones'] += 1
                row['solid_angle_share'] += self.solid_angle(cone.keys)/OCTANT
        unresolved = {}
        for cone in self.cones:
            if cone.status not in ('certified', 'bounded'):
                label, feasible = self.incumbent(cone)
                key = f'{cone.reason or cone.status}:{label}'
                row = unresolved.setdefault(key, dict(cones=0, solid_angle_share=0., incumbent_origin_feasible=feasible))
                row['cones'] += 1
                row['solid_angle_share'] += self.solid_angle(cone.keys)/OCTANT
        inner, outer, outer_box = self.volumes()
        selection = self.selection
        deltas = [c.delta for c in self.cones]
        return dict(
            partition=self.code, sign=self.sign, status=self.status, seconds=self.seconds,
            inner_cert=self.args.inner_cert, criterion=self.criterion, epsilon=self.epsilon,
            volume_gap_ratio=self.volume_ratio() if self.cones else None,
            max_delta_share=(max(deltas)/sum(c.volumes[0] for c in self.cones)) if self.cones else None,
            max_g=max((c.g for c in solved), default=None),
            max_radial_gap=max((c.gap for c in self.cones), default=None),
            radial_certified_solid_angle_share=sum(c.solid for c in self.cones if c.gap <= self.tau)/OCTANT,
            cones=len(self.cones), certified=statuses.count('certified'), bounded=statuses.count('bounded'),
            unresolved=statuses.count('unresolved'), pending=statuses.count('pending'),
            certified_solid_angle_share=sum(self.solid_angle(c.keys) for c in self.cones
                                            if c.status == 'certified')/OCTANT,
            covered_cones=sum(c.cover is not None for c in self.cones), cover_pairs=pairs,
            min_diameter=min((self.diameter(c.keys) for c in self.cones), default=None),
            directions=len(self.directions), counts=dict(self.counts), schemes_used=used,
            unresolved_by_reason_and_incumbent=unresolved,
            split_rules={rule: sum(s['rule'] == rule for s in self.splits)
                         for rule in sorted({s['rule'] for s in self.splits})},
            tightened_schemes=[self.label(s) for s in self.equations.boxes],
            origin={self.label(s): row for s, row in self.origin.items()},
            initial=dict(x0=None if 'x0' not in selection else self.label(selection['x0']),
                         chosen=None if 'chosen' not in selection else self.label(selection['chosen']),
                         excluded_origin_infeasible=[self.label(s) for s in selection.get('excluded', [])],
                         history=selection.get('history', [])),
            lazy_obbt=self.lazy, errors=self.errors,
            inner_volume_kw3=inner, outer_volume_kw3=outer, outer_box_volume_kw3=outer_box, metrics=metrics)


# ---- 扫描缓存与指标 --------------------------------------------------------------------------
def load_scan3(network, budget, path=None):
    """项目标准扫描缓存：校验物理身份，展开格心（带符号 kW）。"""
    if path is None:
        directory = scan_path(ac_network(network), budget, ROOT/'results'/'scan', MODE)
        sources = sorted(directory.glob('region_*.npz'), key=lambda p: p.stat().st_mtime_ns)
        if not sources:
            raise SystemExit(f'No paired AC/SOCP scan cache in {directory}')
        path = sources[-1]
    scan = load_scan(path, ac_identity(ac_network(network), budget, MODE))
    undecided = {key: int(np.count_nonzero(scan[key] == 0)) for key in ('states', 'socp_states')}
    if any(undecided.values()):
        print(f'scan cache has undecided cells {undecided} (reference solver failures in metadata.errors); '
              f'they are excluded from that reference\'s FR/MR', flush=True)
    shape = scan['states'].shape
    index = np.indices(shape).reshape(len(shape), -1).T
    scan.update(path=Path(path), shape=shape, index=index, power=scan['origin']+(scan['start']+index+.5)*scan['step'],
                ac=scan['states'].ravel(), socp=scan['socp_states'].ravel())
    return scan


def metrics3(inner, outer, ac_states, socp_states, volumes):
    """FR=|域∧非可行|/|域|，MR=|可行∧非域|/|可行|（百分数）；体积为 kW³。参考标签为 0（未决）的格不计入该参考的
    FR/MR，单列 undecided_*_cells。"""
    ac, socp = ac_states == 1, socp_states == 1
    known_ac, known_socp = ac_states != 0, socp_states != 0
    pair = lambda domain, feasible, known: comparison_metrics(domain[known], feasible[known])
    return dict(cells=int(len(inner)), ac_feasible_cells=int(ac.sum()), socp_feasible_cells=int(socp.sum()),
                undecided_ac_cells=int((~known_ac).sum()), undecided_socp_cells=int((~known_socp).sum()),
                inner_cells=int(inner.sum()), outer_cells=int(outer.sum()),
                inner_ac=pair(inner, ac, known_ac), inner_socp=pair(inner, socp, known_socp),
                outer_ac=pair(outer, ac, known_ac), outer_socp=pair(outer, socp, known_socp),
                outer_missed_ac_cells=int((ac & ~outer).sum()), sandwich_cells=int((outer & ~inner).sum()),
                inner_volume_kw3=volumes[0], outer_volume_kw3=volumes[1], outer_box_volume_kw3=volumes[2])


def evaluate3(part, scan):
    result = metrics3(part.cell_inner, part.cell_outer, scan['ac'][part.cell_index], scan['socp'][part.cell_index],
                      part.volumes())
    violations = []
    for j in np.flatnonzero(part.cell_ac & ~part.cell_outer):
        cone = next((c for c in part.cones if c.id == part.cell_owner[j]), None)
        halfspace = None if cone is None else cone.halfspace()
        violations.append(dict(p_kw=scan['power'][part.cell_index[j]], grid_index=scan['index'][part.cell_index[j]],
                               xi=part.cell_xi[j], cone=None if cone is None else cone.id,
                               c_xi=None if halfspace is None else float(part.cell_xi[j]@halfspace[0]),
                               mu=None if halfspace is None else halfspace[1]))
    return result, violations


def overall3(labels, volumes, scan):
    """全部分区的逐格指标；labels 为各分区的 cell_index/inner/outer，volumes 为各分区体积三元组。"""
    if not labels:
        return None
    cat = lambda name: np.concatenate([row[name] for row in labels])
    index = cat('cell_index')
    return metrics3(cat('inner'), cat('outer'), scan['ac'][index], scan['socp'][index],
                    np.sum(volumes, axis=0).tolist())


# ---- 图 --------------------------------------------------------------------------------------
PARTITION_COLORS = dict(zip(PARTITIONS, ('#2a78d6', '#eb6834', '#1baf7a', '#eda100',
                                         '#e87ba4', '#008300', '#4a3aa7', '#e34948')))   # 分类色八槽，固定顺序
CELL_CLASSES = (('Inner and AC-feasible', '#9fd3a8'), ('Inner, AC-infeasible (FR)', '#d55181'),
                ('AC-feasible, not inner (MR)', '#2a78d6'), ('SOCP-only feasible, not inner', '#f9d5c3'),
                ('Outer band beyond SOCP', '#d7d6cf'), ('AC-feasible outside outer (must be none)', INK))


def part_labels(part):
    return dict(cell_index=part.cell_index, inner=part.cell_inner, outer=part.cell_outer)


def cell_classes(scan, labels):
    """全网格逐格分类（见 CELL_CLASSES，0 为不可行/未评价）；labels 为各分区的 cell_index/inner/outer。"""
    inner, outer = np.zeros(len(scan['ac']), bool), np.zeros(len(scan['ac']), bool)
    for row in labels:
        inner[row['cell_index']], outer[row['cell_index']] = row['inner'], row['outer']
    ac, socp = scan['ac'] == 1, scan['socp'] == 1
    classes = np.zeros(len(ac), int)
    classes[outer & ~inner & ~socp] = 5
    classes[socp & ~ac & ~inner] = 4
    classes[ac & ~inner] = 3
    classes[inner & ~ac] = 2
    classes[inner & ac] = 1
    classes[ac & ~outer] = 6
    return classes.reshape(scan['shape'])


def plot_slices(path, scan, labels, title, count=6):
    plt = _plt()
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch
    classes = cell_classes(scan, labels)
    shape, step = scan['shape'], scan['step']
    lower = scan['origin']+scan['start']*step
    edges = [lower[j]+np.arange(shape[j]+1)*step[j] for j in range(3)]
    levels = np.flatnonzero((scan['states'] == 1).any(axis=(0, 1)))
    picks = np.unique(np.round(np.linspace(levels.min(), levels.max(), count)).astype(int))
    figure, axes = plt.subplots(2, -(-len(picks)//2), figsize=(6.*(-(-len(picks)//2)), 11.))
    figure.patch.set_facecolor(SURFACE)
    cmap = ListedColormap([SURFACE, *(color for _, color in CELL_CLASSES)])
    for ax, k in zip(axes.flat, picks):
        ax.set_facecolor(SURFACE)
        ax.pcolormesh(edges[0], edges[1], classes[:, :, k].T, cmap=cmap, vmin=-.5, vmax=len(CELL_CLASSES)+.5,
                      shading='flat', rasterized=True)
        ax.axhline(0., color='#c3c2b7', lw=.6)
        ax.axvline(0., color='#c3c2b7', lw=.6)
        ax.set_aspect('equal')
        ax.set_xlabel('p18 (kW, + load / - PV)')
        ax.set_ylabel('p25 (kW, + load / - PV)')
        ax.set_title(f'p30 = {lower[2]+(k+.5)*step[2]:.0f} kW (slice {k})', color=INK, fontsize=10, loc='left')
    for ax in list(axes.flat)[len(picks):]:
        ax.set_visible(False)
    handles = [Patch(facecolor=color, label=name) for name, color in CELL_CLASSES]
    figure.legend(handles=handles, loc='lower center', ncol=3, frameon=False, fontsize=9, labelcolor=INK_2,
                  bbox_to_anchor=(.5, .005))
    figure.suptitle(title, color=INK, fontsize=11, x=.01, ha='left')
    figure.tight_layout(rect=(0, .07, 1, .96))
    figure.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(figure)


def _series3(rows, key):
    points = []
    for row in rows:
        if key == 'volume_gap':
            value = (row['outer_box_volume_kw3']/row['inner_volume_kw3']-1.) if row['inner_volume_kw3'] > 0 else None
        elif key in ('certified_angle', 'radial_certified_angle'):
            value = None if row.get(key) is None else 100.*row[key]
        else:
            value = row.get(key)
        if value is not None and np.isfinite(value) and (key.endswith('certified_angle') or value > 0):
            points.append((row['step'], value))
    return [p[0] for p in points], [p[1] for p in points]


def plot_convergence3(path, histories, tau, title, criterion='radial', epsilon=None):
    plt = _plt()
    from matplotlib.lines import Line2D
    figure, axes = plt.subplots(2, 2, figsize=(13., 8.6))
    figure.patch.set_facecolor(SURFACE)
    if criterion == 'volume':
        panels = (('volume_gap_ratio', 'Volume gap ratio  sum Delta_k / sum vol(T_k)  (stopping criterion)', True),
                  ('max_gap', 'Max radial gap bound  max_k (mu_k - 1)  (diagnostic)', True),
                  ('radial_certified_angle', 'Share of solid angle with mu_k - 1 <= tau (%)', False),
                  ('inner_ac', 'Inner region vs AC scan (%): MR solid, FR dotted (zeros not shown)', True))
    else:
        panels = (('max_gap', 'Max radial gap bound  max_k (mu_k - 1)', True),
                  ('volume_gap', 'Volume gap  (outer∩box - inner) / inner', True),
                  ('certified_angle', 'Certified share of the octant solid angle (%)', False),
                  ('inner_ac', 'Inner region vs AC scan (%): MR solid, FR dotted (zeros not shown)', True))
    thresholds = {'max_gap': ('tau', tau), 'volume_gap_ratio': ('epsilon', epsilon)}
    for ax, (key, name, log) in zip(axes.flat, panels):
        for code, rows in histories.items():
            color = PARTITION_COLORS.get(code, INK_2)
            for series, style in ((('inner_ac_mr', dict(lw=2)), ('inner_ac_fr', dict(lw=1.5, ls=(0, (1, 1.5)))))
                                  if key == 'inner_ac' else ((key, dict(lw=2)),)):
                xs, ys = _series3(rows, series)
                if xs:
                    ax.plot(xs, ys, color=color, drawstyle='steps-post', **style)
                    if style.get('ls') is None:
                        ax.annotate(code, (xs[-1], ys[-1]), xytext=(4, 0), textcoords='offset points',
                                    color=INK_2, fontsize=8, va='center')
        if key in thresholds and thresholds[key][1] is not None:
            label, value = thresholds[key]
            ax.axhline(value, color=MUTED, lw=1)
            ax.annotate(f'{label} = {value:g}', (1., value), xycoords=('axes fraction', 'data'), xytext=(-4, 4),
                        textcoords='offset points', ha='right', color=INK_2, fontsize=8)
        _frame(ax, name, 'Cone MISOCP solves', log)
    handles = [Line2D([], [], color=PARTITION_COLORS[code], lw=2, label=code) for code in histories]
    figure.legend(handles=handles, loc='lower center', ncol=len(handles), frameon=False, fontsize=9,
                  labelcolor=INK_2, bbox_to_anchor=(.5, .005))
    figure.suptitle(title, color=INK, fontsize=11, x=.01, ha='left')
    figure.tight_layout(rect=(0, .05, 1, .96))
    figure.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(figure)


# ---- 可旋转三维页面 ----------------------------------------------------------------------------
PAGE_LAYERS = (('inner', 'Inner region (certified tetrahedra)', '--inner'),
               ('outer', 'Outer bound (mu * T_k)', '--outer'),
               ('ac', 'AC-feasible boundary cells', '--ac'),
               ('mr', 'MR cells: AC-feasible, outside inner', '--mr'),
               ('fr', 'FR cells: inner, AC-infeasible', '--fr'))


def page_data(scan, cone_rows, labels, results, overall, settings):
    """页面数据：各分区内三角/外界三角网格（带符号 kW，取整）、逐格点层（AC 可行边界格、MR 格、FR 格）与指标。"""
    meshes = {'inner': {}, 'outer': {}}
    for code, row in cone_rows.items():
        for kind in meshes:
            x, y, z, faces = [], [], [], []
            for cone in row['cones']:
                vertices = cone['verts_kw'] if kind == 'inner' else cone['outer_kw']
                if vertices is None or len(vertices) != 3:   # 无外界的锥（评价盒多边形）不画外界
                    continue
                base = len(x)
                for v in vertices:
                    x.append(round(v[0]))
                    y.append(round(v[1]))
                    z.append(round(v[2]))
                faces.append((base, base+1, base+2))
            meshes[kind][code] = dict(x=x, y=y, z=z, i=[f[0] for f in faces], j=[f[1] for f in faces],
                                      k=[f[2] for f in faces])
    shape = scan['shape']
    ac_grid = (scan['states'] == 1)
    padded = np.pad(ac_grid, 1, constant_values=False)
    interior = np.ones(shape, bool)
    for axis in range(3):
        for shift in (-1, 1):
            interior &= np.roll(padded, shift, axis=axis)[1:-1, 1:-1, 1:-1]
    boundary = (ac_grid & ~interior).ravel()
    inner = np.zeros(len(scan['ac']), bool)
    evaluated = np.zeros(len(scan['ac']), bool)
    for row in labels.values():
        inner[row['cell_index']] = row['inner']
        evaluated[row['cell_index']] = True
    ac = scan['ac'] == 1
    masks = dict(ac=boundary & evaluated, mr=ac & ~inner & evaluated, fr=inner & (scan['ac'] == -1) & evaluated)
    letters = np.where(scan['power'] >= 0., 'p', 'n')   # 与 attach_cells 的分区归属一致
    codes = np.char.add(np.char.add(letters[:, 0], letters[:, 1]), letters[:, 2])
    points = {}
    for kind, mask in masks.items():
        points[kind] = {}
        for code in cone_rows:
            chosen = mask & (codes == code)
            p = np.rint(scan['power'][chosen]).astype(int)
            points[kind][code] = dict(x=p[:, 0].tolist(), y=p[:, 1].tolist(), z=p[:, 2].tolist())
    table = []
    for r in results:
        m = r['metrics'] or {}
        table.append(dict(partition=r['partition'], status=r['status'], cones=r['cones'],
                          volume_gap=r.get('volume_gap_ratio'), radial_gap=r.get('max_radial_gap'),
                          ac_fr=m.get('inner_ac', {}).get('fr_percent'), ac_mr=m.get('inner_ac', {}).get('mr_percent'),
                          socp_fr=m.get('inner_socp', {}).get('fr_percent'),
                          socp_mr=m.get('inner_socp', {}).get('mr_percent'),
                          ac_out=m.get('outer_missed_ac_cells'), inner_volume=m.get('inner_volume_kw3'),
                          seconds=r['seconds']))
    volumes = [r['metrics']['inner_volume_kw3'] for r in results if r['metrics']]
    gaps = [r.get('volume_gap_ratio') for r in results if r['metrics']]
    total_gap = (sum(g*v for g, v in zip(gaps, volumes))/sum(volumes)
                 if volumes and all(g is not None and np.isfinite(g) for g in gaps) else None)
    return plain(dict(partitions=list(cone_rows), meshes=meshes, points=points, table=table,
                      overall=dict(ac_fr=overall['inner_ac']['fr_percent'], ac_mr=overall['inner_ac']['mr_percent'],
                                   socp_fr=overall['inner_socp']['fr_percent'],
                                   socp_mr=overall['inner_socp']['mr_percent'],
                                   ac_out=overall['outer_missed_ac_cells'], volume_gap=total_gap,
                                   cones=sum(r['cones'] for r in results), cells=overall['cells'],
                                   ac_cells=overall['ac_feasible_cells'],
                                   undecided_socp=overall.get('undecided_socp_cells', 0)),
                      settings=settings))


PAGE = r'''<title>Case33 3D Dispatch Region</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans+Condensed:wght@500;600&family=IBM+Plex+Sans:wght@400;500&display=swap">
<style>
/* Layout: run facts and overall FR/MR on top, the rotatable 3D view beside its layer and partition controls, then the per-partition table. */
:root {
  --bg: #f3f5f4; --panel: #ffffff; --fg: #17201c; --muted: #56625c; --line: #dce2de; --accent: #1b6985;
  --inner: #2b9858; --outer: #d04848; --ac: #3a7bd5; --mr: #1c4c9a; --fr: #bf3b7a;
  --display: "IBM Plex Sans Condensed", "Arial Narrow", system-ui, sans-serif;
  --body: "IBM Plex Sans", system-ui, -apple-system, "Segoe UI", sans-serif;
  --mono: "IBM Plex Mono", ui-monospace, "Cascadia Mono", Consolas, monospace;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --bg: #101513; --panel: #171d1a; --fg: #e4eae6; --muted: #98a49e; --line: #29322d; --accent: #5ab0cf;
    --inner: #4bbb77; --outer: #ec6a6a; --ac: #6da5ec; --mr: #92b8ff; --fr: #e66ca8; color-scheme: dark;
  }
}
:root[data-theme="dark"] {
  --bg: #101513; --panel: #171d1a; --fg: #e4eae6; --muted: #98a49e; --line: #29322d; --accent: #5ab0cf;
  --inner: #4bbb77; --outer: #ec6a6a; --ac: #6da5ec; --mr: #92b8ff; --fr: #e66ca8; color-scheme: dark;
}
body { background: var(--bg); color: var(--fg); font-family: var(--body); font-size: 14px; line-height: 1.5; }
.page { max-width: 1320px; margin: 0 auto; padding-inline: 16px; padding-block: 20px 36px; display: grid; gap: 18px; }
header { display: grid; gap: 6px; }
h1 { font-family: var(--display); font-weight: 600; font-size: clamp(22px, 3vw, 30px); line-height: 1.15; margin: 0;
     letter-spacing: .01em; text-wrap: balance; }
.lede { margin: 0; color: var(--muted); max-width: 72ch; }
.facts { display: flex; flex-wrap: wrap; gap: 4px 18px; font-family: var(--mono); font-size: 12.5px; color: var(--muted); }
.metrics { display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 10px; }
.metric { background: var(--panel); border: 1px solid var(--line); border-radius: 6px; padding: 10px 12px; display: grid; gap: 2px; }
.metric .label { font-size: 11px; letter-spacing: .07em; text-transform: uppercase; color: var(--muted); }
.metric .value { font-family: var(--mono); font-size: 21px; font-variant-numeric: tabular-nums; }
.metric .sub { font-size: 12px; color: var(--muted); }
.viewer { display: grid; grid-template-columns: minmax(0, 1fr) 270px; gap: 14px; align-items: start; }
@media (max-width: 880px) { .viewer { grid-template-columns: minmax(0, 1fr); } }
#plot { background: var(--panel); border: 1px solid var(--line); border-radius: 6px; height: clamp(440px, 74vh, 780px); min-width: 0; }
.panel { background: var(--panel); border: 1px solid var(--line); border-radius: 6px; padding: 12px 14px; display: grid; gap: 16px; }
.panel h2 { font-family: var(--display); font-size: 13px; font-weight: 600; letter-spacing: .08em; text-transform: uppercase;
            margin: 0 0 6px; color: var(--muted); }
.check { display: flex; align-items: flex-start; gap: 8px; padding-block: 3px; cursor: pointer; }
.check input { accent-color: var(--accent); margin-top: 3px; }
.swatch { width: 11px; height: 11px; border-radius: 2px; flex: none; margin-top: 5px; }
.parts { display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 2px 6px; font-family: var(--mono); font-size: 13px; }
.buttons { display: flex; gap: 8px; margin-top: 6px; }
button { font: inherit; font-size: 12.5px; color: var(--fg); background: transparent; border: 1px solid var(--line);
         border-radius: 4px; padding: 3px 10px; cursor: pointer; }
button:hover { border-color: var(--accent); }
label.slider { display: grid; gap: 2px; font-size: 12.5px; color: var(--muted); }
input[type=range] { accent-color: var(--accent); width: 100%; }
:focus-visible { outline: 2px solid var(--accent); outline-offset: 2px; }
.hint { font-size: 12.5px; color: var(--muted); margin: 0; }
.table-wrap { overflow-x: auto; background: var(--panel); border: 1px solid var(--line); border-radius: 6px; }
table { border-collapse: collapse; width: 100%; font-size: 13px; }
caption { text-align: left; padding: 10px 12px 4px; font-family: var(--display); font-weight: 600; font-size: 14px; }
th, td { padding: 7px 12px; text-align: right; border-bottom: 1px solid var(--line); white-space: nowrap; }
th:first-child, td:first-child { text-align: left; }
th { font-weight: 500; color: var(--muted); font-size: 12px; }
td { font-family: var(--mono); font-variant-numeric: tabular-nums; }
tr:last-child td { border-bottom: 0; }
tr.total td { font-weight: 500; border-top: 1px solid var(--muted); }
.pill { font-family: var(--body); font-size: 12px; padding: 1px 9px; border-radius: 999px; border: 1px solid currentColor; }
.pill.certified { color: var(--inner); }
.pill.time_limit, .pill.unresolved, .pill.error { color: var(--outer); }
.notes { color: var(--muted); font-size: 12.5px; max-width: 95ch; margin: 0; display: grid; gap: 4px; }
</style>
<div class="page">
  <header>
    <h1>Case33 3D Dispatch Region</h1>
    <p class="lede">Radial sandwich on nodes 18, 25 and 30: certified inner tetrahedra, MISOCP outer bound and the paired
      AC/SOCP scan, in signed kW (positive = load, negative = PV).</p>
    <div class="facts" id="facts"></div>
  </header>
  <section class="metrics" id="metrics" aria-label="Overall metrics"></section>
  <section class="viewer">
    <div id="plot" role="img" aria-label="Rotatable 3D view of the inner region, outer bound and scan cells"></div>
    <aside class="panel">
      <div><h2>Layers</h2><div id="layers"></div></div>
      <div><h2>Partitions</h2><div class="parts" id="parts"></div>
        <div class="buttons"><button type="button" id="all">All</button><button type="button" id="none">None</button></div></div>
      <div><h2>Opacity</h2>
        <label class="slider" for="inner-opacity">Inner region <input type="range" id="inner-opacity" min="0.1" max="1" step="0.05" value="0.55"></label>
        <label class="slider" for="outer-opacity">Outer bound <input type="range" id="outer-opacity" min="0.05" max="1" step="0.05" value="0.15"></label>
      </div>
      <p class="hint">Drag to rotate, scroll to zoom, right-drag to pan. Partition codes follow nodes 18, 25, 30
        (p = load, n = PV).</p>
    </aside>
  </section>
  <section class="table-wrap"><table id="table"><caption>Per-partition results</caption></table></section>
  <div class="notes" id="notes"></div>
</div>
<script src="https://cdnjs.cloudflare.com/ajax/libs/plotly.js/2.34.0/plotly.min.js"></script>
<script type="application/json" id="data">__DATA__</script>
<script>
(function () {
  const DATA = JSON.parse(document.getElementById('data').textContent);
  const LAYERS = __LAYERS__;
  const state = { layers: {}, parts: new Set(DATA.partitions), inner: 0.55, outer: 0.15 };
  LAYERS.forEach(([id]) => { state.layers[id] = true; });
  const css = name => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
  const pct = v => v === null || v === undefined ? '–' : v.toFixed(2) + '%';
  const num = (v, d) => v === null || v === undefined ? '–' : v.toFixed(d);
  const S = DATA.settings, O = DATA.overall;

  document.getElementById('facts').innerHTML = [
    `criterion ${S.criterion}${S.criterion === 'volume' ? ' (ε = ' + S.epsilon + ')' : ''}`, `τ = ${S.tau}`,
    `inner certificate: ${S.inner_cert}`, `${S.seconds} s per partition`, `${O.cones} cones`,
    `scan grid ${S.scan_shape.join('×')} (${O.cells.toLocaleString()} cells)`].map(t => `<span>${t}</span>`).join('');
  const tiles = [
    ['AC · FR', pct(O.ac_fr), 'inner cells AC-infeasible'], ['AC · MR', pct(O.ac_mr), 'AC-feasible cells missed'],
    ['SOCP · FR', pct(O.socp_fr), 'inner cells SOCP-infeasible'], ['SOCP · MR', pct(O.socp_mr), 'SOCP-feasible cells missed'],
    ['Volume gap', O.volume_gap === null ? '–' : (100 * O.volume_gap).toFixed(2) + '%', 'Σ Δ_k / Σ vol(T_k)'],
    ['AC outside outer', String(O.ac_out), 'must be 0']];
  document.getElementById('metrics').innerHTML = tiles.map(([l, v, s]) =>
    `<div class="metric"><span class="label">${l}</span><span class="value">${v}</span><span class="sub">${s}</span></div>`).join('');

  document.getElementById('layers').innerHTML = LAYERS.map(([id, label, color]) =>
    `<label class="check" for="layer-${id}"><input type="checkbox" id="layer-${id}" checked>` +
    `<span class="swatch" style="background: var(${color})"></span><span>${label}</span></label>`).join('');
  document.getElementById('parts').innerHTML = DATA.partitions.map(code =>
    `<label class="check" for="part-${code}"><input type="checkbox" id="part-${code}" checked><span>${code}</span></label>`).join('');

  const rows = DATA.table.map(r => `<tr><td>${r.partition}</td><td><span class="pill ${r.status}">${r.status.replace('_', ' ')}</span></td>` +
    `<td>${r.cones}</td><td>${r.volume_gap === null ? '–' : (100 * r.volume_gap).toFixed(2) + '%'}</td><td>${num(r.radial_gap, 4)}</td>` +
    `<td>${pct(r.ac_fr)}</td><td>${pct(r.ac_mr)}</td><td>${pct(r.socp_fr)}</td><td>${pct(r.socp_mr)}</td>` +
    `<td>${r.ac_out}</td><td>${r.inner_volume === null ? '–' : r.inner_volume.toExponential(3)}</td><td>${num(r.seconds, 0)}</td></tr>`).join('');
  document.getElementById('table').insertAdjacentHTML('beforeend',
    '<thead><tr><th>Partition</th><th>Status</th><th>Cones</th><th>Volume gap</th><th>Max radial gap</th>' +
    '<th>AC FR</th><th>AC MR</th><th>SOCP FR</th><th>SOCP MR</th><th>AC outside outer</th><th>Inner volume (kW³)</th><th>Time (s)</th></tr></thead>' +
    `<tbody>${rows}<tr class="total"><td>Overall</td><td></td><td>${O.cones}</td>` +
    `<td>${O.volume_gap === null ? '–' : (100 * O.volume_gap).toFixed(2) + '%'}</td><td></td><td>${pct(O.ac_fr)}</td><td>${pct(O.ac_mr)}</td>` +
    `<td>${pct(O.socp_fr)}</td><td>${pct(O.socp_mr)}</td><td>${O.ac_out}</td><td></td><td></td></tr></tbody>`);
  document.getElementById('notes').innerHTML = [
    `FR = inner cells that the reference marks infeasible / inner cells; MR = reference-feasible cells outside the inner region / reference-feasible cells. ` +
    `Cells are the ${S.scan_shape.join('×')} paired scan (${O.ac_cells.toLocaleString()} AC-feasible); ${O.undecided_socp} SOCP cells stayed undecided and are left out of the SOCP rates.`,
    `Volume gap Δ_k = (μ̄_k³ − 1)·vol(T_k) bounds both the missed and the extra volume in cone k. Max radial gap is reported as a diagnostic only.`]
    .map(t => `<p>${t}</p>`).join('');

  function mesh(kind, color, opacity) {
    let x = [], y = [], z = [], i = [], j = [], k = [];
    DATA.partitions.forEach(code => {
      if (!state.parts.has(code)) return;
      const m = DATA.meshes[kind][code];
      if (!m) return;
      const offset = x.length;
      x = x.concat(m.x); y = y.concat(m.y); z = z.concat(m.z);
      i = i.concat(m.i.map(v => v + offset)); j = j.concat(m.j.map(v => v + offset)); k = k.concat(m.k.map(v => v + offset));
    });
    return { type: 'mesh3d', x, y, z, i, j, k, color, opacity, flatshading: true, hoverinfo: 'skip',
             lighting: { ambient: 0.62, diffuse: 0.55, specular: 0.05, roughness: 0.9 } };
  }
  function cloud(kind, color, size, opacity) {
    let x = [], y = [], z = [];
    DATA.partitions.forEach(code => {
      if (!state.parts.has(code)) return;
      const p = DATA.points[kind][code];
      x = x.concat(p.x); y = y.concat(p.y); z = z.concat(p.z);
    });
    return { type: 'scatter3d', mode: 'markers', x, y, z, marker: { size, color, opacity }, hovertemplate:
             'p18 %{x} kW<br>p25 %{y} kW<br>p30 %{z} kW<extra></extra>' };
  }
  function axis(title) {
    return { title: { text: title }, color: css('--muted'), gridcolor: css('--line'), zerolinecolor: css('--muted'),
             showbackground: false, tickfont: { family: css('--mono'), size: 11 } };
  }
  function render() {
    const traces = [];
    if (state.layers.inner) traces.push(mesh('inner', css('--inner'), state.inner));
    if (state.layers.outer) traces.push(mesh('outer', css('--outer'), state.outer));
    if (state.layers.ac) traces.push(cloud('ac', css('--ac'), 1.6, 0.35));
    if (state.layers.mr) traces.push(cloud('mr', css('--mr'), 2.4, 0.85));
    if (state.layers.fr) traces.push(cloud('fr', css('--fr'), 2.6, 0.9));
    const layout = { paper_bgcolor: css('--panel'), plot_bgcolor: css('--panel'), showlegend: false, uirevision: 'keep',
      margin: { l: 0, r: 0, t: 0, b: 0 }, font: { family: css('--body'), color: css('--fg') },
      scene: { aspectmode: 'data', xaxis: axis('p18 (kW)'), yaxis: axis('p25 (kW)'), zaxis: axis('p30 (kW)'),
               camera: { eye: { x: 1.55, y: -1.55, z: 0.95 } } } };
    if (window.Plotly) Plotly.react('plot', traces, layout, { responsive: true, displaylogo: false });
  }
  LAYERS.forEach(([id]) => document.getElementById('layer-' + id).addEventListener('change', e => {
    state.layers[id] = e.target.checked; render(); }));
  DATA.partitions.forEach(code => document.getElementById('part-' + code).addEventListener('change', e => {
    if (e.target.checked) state.parts.add(code); else state.parts.delete(code); render(); }));
  const setAll = on => { DATA.partitions.forEach(code => { document.getElementById('part-' + code).checked = on;
    if (on) state.parts.add(code); else state.parts.delete(code); }); render(); };
  document.getElementById('all').addEventListener('click', () => setAll(true));
  document.getElementById('none').addEventListener('click', () => setAll(false));
  document.getElementById('inner-opacity').addEventListener('input', e => { state.inner = +e.target.value; render(); });
  document.getElementById('outer-opacity').addEventListener('input', e => { state.outer = +e.target.value; render(); });
  window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', render);
  new MutationObserver(render).observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] });
  render();
})();
</script>
'''


def write_region_page(path, data):
    """可旋转三维页面：完整 HTML 文档（本地直接打开），内容部分即 PAGE。"""
    body = (PAGE.replace('__LAYERS__', json.dumps([list(row) for row in PAGE_LAYERS]))
                .replace('__DATA__', json.dumps(data, separators=(',', ':')).replace('</', '<\\/')))
    path.write_text('<!doctype html>\n<html lang="en"><head><meta charset="utf-8">'
                    '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover"></head><body>\n'
                    + body + '\n</body></html>\n', encoding='utf-8')


def render_outputs(output, args, scan, cone_rows, labels, results, overall):
    epsilon = args.epsilon if args.epsilon is not None else 3*args.tau
    plot_convergence3(output/'convergence.png', {code: row['history'] for code, row in cone_rows.items()},
                      args.tau, f'Case33 (18, 25, 30), mode 1 - convergence, criterion: {args.criterion}'
                                f'{f" (epsilon = {epsilon:g})" if args.criterion == "volume" else ""}, '
                                f'tau = {args.tau:g}, inner certificate: {args.inner_cert}',
                      criterion=args.criterion, epsilon=epsilon)
    if scan is None or not labels:
        return
    plot_slices(output/'slices.png', scan, list(labels.values()),
                f'Case33 (18, 25, 30), mode 1 - cell classes on p30 slices, criterion: {args.criterion}, '
                f'tau = {args.tau:g}, inner certificate: {args.inner_cert}')
    settings = dict(criterion=args.criterion, epsilon=epsilon, tau=args.tau, inner_cert=args.inner_cert,
                    seconds=args.seconds, scan_shape=list(scan['shape']))
    write_region_page(output/'region3d.html', page_data(scan, cone_rows, labels, results, overall, settings))


# ---- 终端输出与入口 ----------------------------------------------------------------------------
def print_table3(rows):
    header = (f'{"partition":<10s}{"in-AC FR%":>10s}{"in-AC MR%":>10s}{"in-SOCP FR%":>12s}{"in-SOCP MR%":>12s}'
              f'{"AC-out":>8s}{"band":>8s}{"inner kW3":>12s}{"outer kW3":>12s}{"outer^box":>12s}')
    print(header)
    print('-'*len(header))
    for name, m in rows:
        if m is None:
            continue
        print(f'{name:<10s}{fmt(m["inner_ac"]["fr_percent"]):>10s}{fmt(m["inner_ac"]["mr_percent"]):>10s}'
              f'{fmt(m["inner_socp"]["fr_percent"]):>12s}{fmt(m["inner_socp"]["mr_percent"]):>12s}'
              f'{m["outer_missed_ac_cells"]:>8d}{m["sandwich_cells"]:>8d}'
              f'{m["inner_volume_kw3"]:>12.4g}{m["outer_volume_kw3"]:>12.4g}{m["outer_box_volume_kw3"]:>12.4g}')


def summary_line3(result):
    counts = result['counts']
    if result.get('criterion') == 'volume':
        head = (f'{result["partition"]} {result["status"]:<10s} volume gap={fmt(result["volume_gap_ratio"], 4)} '
                f'(eps {result["epsilon"]:g}) max radial gap={fmt(result["max_radial_gap"], 4)} cones={result["cones"]} '
                f'(bounded {result["bounded"]}, unres {result["unresolved"]}, pend {result["pending"]}, '
                f'radially within tau {100*result["radial_certified_solid_angle_share"]:.1f}%) ')
    else:
        head = (f'{result["partition"]} {result["status"]:<10s} max g={fmt(result["max_g"], 5)} cones={result["cones"]} '
                f'(cert {result["certified"]}, unres {result["unresolved"]}, pend {result["pending"]}, '
                f'certified solid angle {100*result["certified_solid_angle_share"]:.1f}%) ')
    return (head+f'MISOCP={counts["misocp"]} '
            f'lazy OBBT={counts["lazy_obbt"]} rays={counts["rays"]}+{counts["near_rays"]} near '
            f'covered={result["covered_cones"]} schemes={{{",".join(sorted(result["schemes_used"]))}}} '
            f'time={result["seconds"]:.1f}s')


def part_rows(part):
    return dict(sign=part.sign, status=part.status, directions=part.directions,
                cones=[part.cone_json(c) for c in part.cones], splits=part.splits, history=part.history)


def write_outputs3(output, args, cone_rows, results, overall, provenance, violations, labels=None):
    """cones.json、summary.json，以及逐格标签 labels.npz（并行合并与切片图用）。"""
    (output/'cones.json').write_text(json.dumps(plain(dict(provenance=provenance, partitions=cone_rows)),
                                                ensure_ascii=False, indent=1), encoding='utf-8')
    if labels:
        np.savez_compressed(output/'labels.npz', **{f'{code}_{key}': row[key] for code, row in labels.items()
                                                     for key in ('cell_index', 'inner', 'outer')})
    statuses = [r['status'] for r in results]
    status = ('certified' if statuses and all(s == 'certified' for s in statuses) else
              'error' if 'error' in statuses else 'time_limit' if 'time_limit' in statuses else 'unresolved')
    summary = plain(dict(provenance=provenance, status=status,
                         settings=dict(tau=args.tau, mip_seconds=args.mip_seconds, seconds=args.seconds,
                                       max_cones=args.max_cones, min_width=args.min_width, threads=args.threads,
                                       obbt_workers=args.obbt_workers, partitions=list(args.partitions),
                                       inner_cert=args.inner_cert, criterion=args.criterion,
                                       epsilon=args.epsilon if args.epsilon is not None else 3*args.tau,
                                       mip_gap=args.mip_gap, split_margin=SPLIT_MARGIN,
                                       cover_conv_tol=COVER_CONV_TOL, ray_seconds=RAY_SECONDS),
                         partitions={r['partition']: r for r in results}, overall=overall,
                         outer_validity_violations=violations))
    (output/'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=1), encoding='utf-8')
    return summary


def parse_args():
    parser = argparse.ArgumentParser(description='Radial cone partition with inner/outer sandwich on Case33 (18,25,30), mode 1')
    parser.add_argument('--partitions', default=','.join(PARTITIONS),
                        type=lambda value: tuple(code for code in value.split(',') if partition_sign(code, 3)),
                        help='comma list of the 8 sign codes (order as load_nodes; p positive, n negative power)')
    parser.add_argument('--tau', type=float, default=.005)
    parser.add_argument('--mip-seconds', type=float, default=60.)
    parser.add_argument('--seconds', type=float, default=1800., help='time limit per partition')
    parser.add_argument('--max-cones', type=int, default=2048)
    parser.add_argument('--min-width', type=float, default=2e-3, help='minimum angular diameter of a cone (rad)')
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--obbt-workers', type=int, default=8)
    parser.add_argument('--inner-cert', choices=('origin', 'near'), default='near')
    parser.add_argument('--criterion', choices=('volume', 'radial'), default='volume',
                        help='volume: refine max Delta_k, stop at sum Delta_k <= eps*sum vol(T_k); '
                             'radial: every cone mu_k <= 1+tau (previous rule)')
    parser.add_argument('--epsilon', type=float, help='volume-gap tolerance (default d*tau)')
    parser.add_argument('--mip-gap', type=float, default=1e-3, help='relative MIPGap of cone MISOCPs (volume mode)')
    parser.add_argument('--scan', type=Path, help='paired AC/SOCP scan cache (default: newest project cache)')
    parser.add_argument('--no-scan', action='store_true', help='run without grid metrics (development)')
    parser.add_argument('--jobs', type=int, default=1,
                        help='run partitions in this many parallel processes, then merge (each uses --threads)')
    parser.add_argument('--no-plots', action='store_true', help='skip figures (used by parallel workers)')
    parser.add_argument('--render-only', action='store_true',
                        help='re-render figures and the 3D page from the results already in --output')
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT)
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
    scan = None if args.no_scan else load_scan3(network, budget, args.scan)
    if scan is not None and network.fingerprint != scan['metadata']['network_fingerprint']:
        raise SystemExit(f'Network fingerprint differs from scan cache: {scan["path"]}')
    bounds = port_bounds(network)
    provenance = dict(script=Path(__file__).resolve().relative_to(ROOT).as_posix(), commit=git_commit(),
                      argv=sys.argv[1:], network=network.name, load_nodes=list(LOAD_NODES), mode=MODE,
                      budget=budget, current_limit_a=CURRENT_LIMIT, bounds_kw=bounds,
                      network_fingerprint=network.fingerprint,
                      ac_identity=None if scan is None else scan['metadata']['identity'],
                      scan_cache=None if scan is None else scan['path'].resolve().relative_to(ROOT).as_posix(),
                      scan_shape=None if scan is None else list(scan['shape']),
                      switch_order=[c.id for c in network.corridors if c.switchable],
                      scheme_note='switch string over switch_order, 1 = closed; *_kw are signed true kW')
    if args.render_only:   # 只用已保存结果重画图与三维页面，不求解
        summary = json.loads((output/'summary.json').read_text(encoding='utf-8'))
        cone_rows = json.loads((output/'cones.json').read_text(encoding='utf-8'))['partitions']
        with np.load(output/'labels.npz') as saved:
            labels = {code: {key: saved[f'{code}_{key}'] for key in ('cell_index', 'inner', 'outer')}
                      for code in cone_rows if f'{code}_cell_index' in saved}
        args.criterion = summary['settings'].get('criterion', 'radial')
        args.epsilon, args.inner_cert = summary['settings'].get('epsilon'), summary['settings']['inner_cert']
        args.tau, args.seconds = summary['settings']['tau'], summary['settings']['seconds']
        render_outputs(output, args, scan, cone_rows, labels, list(summary['partitions'].values()), summary['overall'])
        print(f'rendered figures and page in {output}', flush=True)
        return 0
    print(f'Case33 {LOAD_NODES} mode={MODE} budget={budget} bounds={np.round(bounds, 3).tolist()} kW, tau={args.tau}, '
          f'partitions={",".join(args.partitions)}, inner certificate={args.inner_cert}, '
          f'scan={None if scan is None else scan["shape"]}, jobs={args.jobs}', flush=True)
    if args.jobs > 1:
        cone_rows, results, labels, violations = run_parallel(args, output)
    else:
        cone_rows, results, labels, violations = {}, [], {}, []
        with threadpool_limits(limits=1):
            for code in args.partitions:
                print(f'== partition {code} sign={partition_sign(code, 3)}', flush=True)
                part = RadialSandwich3D(network, code, bounds, budget, args)
                if scan is not None:
                    part.attach_cells(scan)
                part.run(scan)
                metrics = None
                if scan is not None:
                    metrics, found = evaluate3(part, scan)
                    labels[code] = part_labels(part)
                    if found:
                        violations = [dict(partition=code, **row) for row in found]
                        print(f'!! OUTER VALIDITY FAILED in {code}: {len(violations)} AC-feasible cells outside '
                              f'the outer domain', flush=True)
                        for row in violations[:50]:
                            print(f'   p={np.asarray(row["p_kw"]).tolist()} kW '
                                  f'grid={np.asarray(row["grid_index"]).tolist()} cone={row["cone"]} '
                                  f'c.xi={row["c_xi"]} mu={row["mu"]}', flush=True)
                results.append(part.result(metrics))
                cone_rows[code] = part_rows(part)
                print(summary_line3(results[-1]), flush=True)
                overall = None if scan is None else overall3(list(labels.values()), [
                    (r['inner_volume_kw3'], r['outer_volume_kw3'], r['outer_box_volume_kw3']) for r in results], scan)
                provenance['seconds'] = perf_counter()-started
                write_outputs3(output, args, cone_rows, results, overall, provenance, violations, labels)
                if violations:
                    break
    overall = None if scan is None or not labels else overall3(list(labels.values()), [
        (r['inner_volume_kw3'], r['outer_volume_kw3'], r['outer_box_volume_kw3']) for r in results], scan)
    provenance['seconds'] = perf_counter()-started
    write_outputs3(output, args, cone_rows, results, overall, provenance, violations, labels)
    if not args.no_plots:
        render_outputs(output, args, scan, cone_rows, labels, results, overall)
    if scan is not None:
        print()
        print_table3([*((r['partition'], r['metrics']) for r in results), ('overall', overall)])
    print(f'\noutputs: {output}  total {perf_counter()-started:.1f}s', flush=True)
    return 2 if violations else 1 if any(r['status'] == 'error' for r in results) else 0


def run_parallel(args, output):
    """每个分区一个子进程（同一脚本、单分区、不画图），结束后读回各自的 cones.json/summary.json/labels.npz。"""
    import subprocess
    import time
    workdir = output/'partitions'
    workdir.mkdir(parents=True, exist_ok=True)
    queue, running = list(args.partitions), []
    while queue or running:
        while queue and len(running) < args.jobs:
            code = queue.pop(0)
            command = [sys.executable, '-X', 'utf8', str(Path(__file__).resolve()), '--partitions', code,
                       '--output', str(workdir/code), '--tau', repr(args.tau), '--mip-seconds', repr(args.mip_seconds),
                       '--seconds', repr(args.seconds), '--max-cones', str(args.max_cones),
                       '--min-width', repr(args.min_width), '--threads', str(args.threads),
                       '--obbt-workers', str(args.obbt_workers), '--inner-cert', args.inner_cert, '--no-plots',
                       '--criterion', args.criterion, '--mip-gap', repr(args.mip_gap),
                       *(['--epsilon', repr(args.epsilon)] if args.epsilon is not None else []),
                       *(['--scan', str(args.scan)] if args.scan else []), *(['--no-scan'] if args.no_scan else [])]
            log = (workdir/f'{code}.log').open('w', encoding='utf-8')
            running.append((code, subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, cwd=ROOT), log))
            print(f'started partition {code} (log {workdir/code}.log)', flush=True)
        for item in list(running):
            code, process, log = item
            if process.poll() is not None:
                log.close()
                running.remove(item)
                summary = json.loads((workdir/code/'summary.json').read_text(encoding='utf-8'))
                print(f'finished partition {code}: exit {process.returncode}; '
                      + summary_line3(summary['partitions'][code]), flush=True)
        time.sleep(2.)
    cone_rows, results, labels, violations = {}, [], {}, []
    for code in args.partitions:
        directory = workdir/code
        summary = json.loads((directory/'summary.json').read_text(encoding='utf-8'))
        results.append(summary['partitions'][code])
        violations += summary['outer_validity_violations']
        cone_rows[code] = json.loads((directory/'cones.json').read_text(encoding='utf-8'))['partitions'][code]
        if (directory/'labels.npz').exists():
            with np.load(directory/'labels.npz') as saved:
                labels[code] = {key: saved[f'{code}_{key}'] for key in ('cell_index', 'inner', 'outer')}
    return cone_rows, results, labels, violations


if __name__ == '__main__':
    raise SystemExit(main())
