"""符号分区规划域：径向锥夹逼的全局搜索（R）+ 主线割平面的逐网架构域（CUT），即方法 RCUT。

坐标：符号分区 sign∈{±1}^d 内用非负幅值 u，真实功率 p=sign*u（kW）；几何在 xi=u/bounds∈[0,1]^d 中计算，
bounds 为端口幅值上界。符号（见 docs/notation.md「区域与符号」）：K 与锥相关、N 与网架相关、R 为真实可行域，
上标为物理含义、下标为序号。紧化可行域 R^SOCP=∪_x R^SOCP_x，R^SOCP_x 为网架 x 的 OBBT 紧化可行集（凸），
AC 可行域 R^AC ⊆ R^SOCP。

R：方向用单纯锥 K_k=cone(u_1..u_d) 剖分。锥内三角 K^IN_k=conv(0, v_{k,1..d})，v_{k,i} 为内域网架 x̂ 沿 u_i 的
紧化射线顶点；x̂ 原点不可行时由原点可行的网架 y 覆盖近端，K^IN_k ⊆ R^SOCP_x̂ ∪ R^SOCP_y ⊆ R^SOCP。外界：
Q=[v_1..v_d]，c=1ᵀQ⁻¹，网架不固定的锥 MISOCP 的上界 μ_k 给出 K^OUT_k={xi∈K_k: c@xi<=μ_k}=μ_k·K^IN_k ⊇ R^SOCP∩K_k。
每次细分体积缺口 Δ_k=(μ_k^d-1)·vol(K^IN_k) 最大的锥，ΣΔ_k<=ε·Σvol(K^IN_k) 即认证。初始网架由中心射线 MISOCP
给出。x 自由的 MISOCP 用行生成：一次分支定界中，新现任网架的紧化行以惰性约束加入，尚未紧化的先做 OBBT。
CUT：只割可能改变 I 的网架 X*（叶锥的 x̂、覆盖网架与锥 MISOCP 解点网架）。N^CUT_x 从 conv(K^OUT) ∩ 盒出发
（K^OUT ⊇ R^SOCP ⊇ R^SOCP_x），SP 为顶点评分，违反量 η 最大的顶点由对偶 LP 取联合割（带该网架的 OBBT 盒与包络行）；
连续 patience 次割的体积缩减比例 < threshold 即停滞，N^CUT_x 当作该网架的可行域（外近似，不作内域认证）。一轮 CUT
每割完一个网架求夹逼间隙，已达 ε 或降幅不足 GAP_SHARE·ε 即结束，本轮其余网架不再割。A 阶段（R 的目标放宽为 ε_A）
后割一轮，之后续跑 R（A+），每次锥决策后先割新出现的网架。OBBT 的线程数为 workers // 仍在计算的分区数。
结果内域 I=(K^IN ∪ N^CUT)∩K^OUT；vol(K^OUT)-vol(I)<=ε·vol(I) 即认证。
"""
from concurrent.futures import ProcessPoolExecutor
from copy import deepcopy
from dataclasses import dataclass
import heapq
from itertools import combinations, product
from math import factorial
import multiprocessing
import time

import numpy as np
from gurobipy import GRB
from scipy.spatial import ConvexHull, QhullError
from shapely import set_precision
from shapely.errors import GEOSException
from shapely.geometry import MultiPoint, Polygon
from shapely.ops import unary_union
from threadpoolctl import threadpool_limits

from model import SP_TIME_LIMIT, PortPhysics, PortSubProblem, cone_misocp, port_bounds, ray_support
from monitor import RegionTimeout, RunMonitor, _connect

# 连续几何在公共评价箱归一化后的坐标中计算；与采样网格无关。
GEOMETRY_TOL = 1e-8
SPLIT_MARGIN = .05       # 锥坐标 λ < SPLIT_MARGIN 的生成方向视为“靠近”：决定星形 / 棱 / 最长棱中点剖分
MIN_RADIUS = 1e-9        # 射线顶点半径下限（xi）；更小时 Q 奇异，不作内顶点
RAY_SECONDS = 30.        # 单次射线与原点 SP 的时限上限
CUT_SECONDS = 30.        # 单个网架一次割循环的时限上限
PASS_SHARE = .75         # 一轮 CUT 最多用分区剩余时间的比例，其余留给径向续跑
CUT_ACCEPTED = ('stagnated', 'exact', 'empty')   # 计入 N^CUT 的状态（empty：N^CUT_x=∅，即分区内 R^SOCP_x=∅）
GAP_SHARE = .1           # 一轮 CUT 中某网架使夹逼间隙下降不足 GAP_SHARE·ε 时结束本轮


# ---- 凸多面体几何（xi 或 kW，同一坐标系内运算） -------------------------------------------------
def polytope_vertices(points):
    """任意仿射维数的凸包顶点，保持原始交点坐标。"""
    points = np.asarray(points, dtype=float)
    d = points.shape[1] if points.ndim == 2 else 3
    points = points.reshape(-1, d)
    if not len(points):
        return points
    _, indices = np.unique(np.round(points, 11), axis=0, return_index=True)
    points = points[np.sort(indices)]
    delta = points-points[0]
    rank = np.linalg.matrix_rank(delta, tol=1e-10)
    if rank == 0:
        return points[:1]
    basis = np.linalg.svd(delta, full_matrices=False)[2][:rank]
    coordinates = delta@basis.T
    if rank == 1:
        return points[np.unique([coordinates[:, 0].argmin(), coordinates[:, 0].argmax()])]
    # 一次可逆的仿射缩放，避免极薄多面体的长宽比进入 Qhull；输出仍取原点。
    coordinates /= np.linalg.norm(coordinates, axis=0)
    return points[np.sort(ConvexHull(coordinates).vertices)]


def halfspaces(points):
    """返回 [F,g]，内侧 F@xi+g<=0；低维集以成对不等式表示仿射等式。"""
    points = np.asarray(points)
    d = points.shape[1] if points.ndim == 2 else 3
    if not len(points):
        return np.array([[0.]*d+[1.]])
    center, delta = points[0], points-points[0]
    rank = np.linalg.matrix_rank(delta, tol=1e-10)
    basis = np.linalg.svd(delta, full_matrices=True)[2]
    coordinates = delta@basis[:rank].T
    if rank >= 2:
        scale = np.linalg.norm(coordinates, axis=0)
        hull = ConvexHull(coordinates/scale)
        normal = (hull.equations[:, :rank]/scale)@basis[:rank]
        eq = np.c_[normal, hull.equations[:, rank]-normal@center]
        eq /= np.linalg.norm(normal, axis=1)[:, None]  # 保持 GEOMETRY_TOL 的原单位。
    elif rank == 1:
        normal = np.array([basis[0], -basis[0]])
        eq = np.c_[normal, [-coordinates.max(), coordinates.min()]-normal@center]
    else:
        eq = np.empty((0, d+1))
    normal = basis[rank:]
    eq = np.vstack([eq, np.c_[normal, -normal@center], np.c_[-normal, normal@center]])
    _, indices = np.unique(np.round(eq, 10), axis=0, return_index=True)
    return eq[np.sort(indices)]


def contains(points, equations, tolerance=GEOMETRY_TOL):
    return np.all(np.atleast_2d(points)@equations[:, :-1].T+equations[:, -1] <= tolerance, axis=1)


def covered(points, rows):
    """点是否落在任一凸多面体行（row['vertices']）内：先按包围盒筛点，再逐行判定。"""
    points = np.atleast_2d(points)
    inside = np.zeros(len(points), dtype=bool)
    for row in rows:
        vertices = np.asarray(row['vertices'])
        lower, upper = vertices.min(axis=0)-1e-6, vertices.max(axis=0)+1e-6
        near = np.flatnonzero(~inside & np.all((points >= lower) & (points <= upper), axis=1))
        if len(near):
            inside[near] = contains(points[near], halfspaces(vertices))
    return inside


def clip_polytope(vertices, constant, coefficient):
    """用 constant + coefficient@xi >= 0 裁剪；与输入顶点使用同一坐标系。"""
    d = len(coefficient)
    vertices = np.asarray(vertices).reshape(-1, d)
    if not len(vertices):
        return vertices
    coefficient = np.asarray(coefficient)
    norm = np.linalg.norm(coefficient)
    if norm < 1e-20:
        return vertices if constant >= -1e-12 else np.empty((0, d))
    values = (constant+vertices@coefficient)/norm
    if values.min() >= -1e-11:
        return vertices
    if values.max() < -1e-11:
        return np.empty((0, d))
    points = list(vertices[values >= -1e-11])
    if len(vertices) >= d+1 and np.linalg.matrix_rank(vertices-vertices[0], tol=1e-10) == d:
        coordinates = (vertices-vertices[0])@np.linalg.svd(vertices-vertices[0], full_matrices=False)[2].T
        coordinates /= np.linalg.norm(coordinates, axis=0)
        edges = {tuple(sorted(edge)) for face in ConvexHull(coordinates).simplices
                 for edge in combinations(face, 2)}
    else:
        edges = combinations(range(len(vertices)), 2)
    for a, b in edges:
        if values[a]*values[b] < 0:
            points.append(vertices[a]+(vertices[b]-vertices[a])*values[a]/(values[a]-values[b]))
    return polytope_vertices(points)


def polytope_volume(poly):
    poly = np.asarray(poly)
    d = poly.shape[1] if poly.ndim == 2 else 3
    poly = poly.reshape(-1, d)
    if len(poly) < d+1 or np.linalg.matrix_rank(poly-poly[0], tol=1e-10) < d:
        return 0.
    coordinates = (poly-poly[0])@np.linalg.svd(poly-poly[0], full_matrices=False)[2].T
    scale = np.linalg.norm(coordinates, axis=0)
    return float(ConvexHull(coordinates/scale).volume*np.prod(scale))


def polygon_union(polygons):
    """二维凸多边形（顶点表）的并集：保留不相连部分和孔洞。"""
    return unary_union([MultiPoint(points).convex_hull for points in polygons if len(points)])


def _clip_face(points, constant, coefficient):
    """有序面片裁剪：constant+coefficient@point>=0，不重新求凸包。"""
    values = (constant+points@coefficient)/np.linalg.norm(coefficient)
    values[np.abs(values) <= 1e-11] = 0.
    clipped = []
    for index in range(len(points)):
        previous = index-1
        if values[previous]*values[index] < 0.:
            clipped.append(points[previous]+(points[index]-points[previous])
                           *values[previous]/(values[previous]-values[index]))
        if values[index] >= 0.:
            clipped.append(points[index])
    return np.asarray(clipped).reshape(-1, 2)


def union_volume(polytopes):
    """三维凸多面体并集的体积（散度定理）：只积分并集暴露的边界面，重叠面和内部面不重复计入。"""
    polytopes = [np.asarray(p) for p in polytopes if polytope_volume(p) > 1e-15]
    equations = [halfspaces(p) for p in polytopes]
    volume = 0.
    for i, poly in enumerate(polytopes):
        coordinates = np.linalg.svd(poly-poly[0], full_matrices=False)[0]
        for simplex in ConvexHull(coordinates, qhull_options='Qx').simplices:
            triangle = poly[simplex]
            a, b = triangle[1:]-triangle[0]
            normal = np.cross(a, b)
            area = np.linalg.norm(normal)
            if area == 0.:
                continue
            normal /= area
            if normal@(triangle[0]-poly.mean(axis=0)) < 0.:
                normal = -normal
            offset = -normal@triangle[0]
            basis = np.array([a/np.linalg.norm(a), np.cross(normal, a/np.linalg.norm(a))])
            face = (triangle-triangle[0])@basis.T
            exposed = Polygon(face)
            for j, eq in enumerate(equations):
                if i == j:
                    continue
                # 共面且同向的重叠外表面由索引较小的网架计入一次。
                coincident = ((np.linalg.norm(eq[:, :3]-normal, axis=1) < 1e-10)
                              & (np.abs(eq[:, 3]-offset) < 1e-10))
                if j > i and coincident.any():
                    continue
                values = triangle@eq[:, :3].T+eq[:, 3]
                if np.any(values.min(axis=0) > 1e-10):
                    continue
                if np.all(values <= 1e-10):
                    exposed = Polygon()
                    break
                coefficients = eq[:, :3]@basis.T
                constants = eq[:, :3]@triangle[0]+eq[:, 3]
                parallel = np.linalg.norm(coefficients, axis=1) < 1e-10
                if np.any(constants[parallel] > 1e-10):
                    continue
                overlap = face
                for coefficient, constant in zip(coefficients[~parallel], constants[~parallel]):
                    overlap = _clip_face(overlap, -constant, -coefficient)
                    if not len(overlap):
                        break
                if len(overlap) >= 3:
                    # 共边点统一到面内裁剪精度，避免微小共线偏差使重叠面漏扣。
                    exposed = set_precision(exposed, 1e-11).difference(set_precision(Polygon(overlap), 1e-11))
                    if exposed.is_empty:
                        break
            volume -= offset*exposed.area/3.
    return volume


def union_measure(polytopes, d):
    """二维面积 / 三维体积；重叠部分只计一次。"""
    return polygon_union(polytopes).area if d == 2 else union_volume(polytopes)


def register_power(power, powers, point_tol):
    """以 kW 最大坐标差匹配固定代表；不舍入坐标，也不沿近点链移动代表。"""
    for index, representative in enumerate(powers):
        if np.max(np.abs(power-representative)) <= point_tol:
            return index
    powers.append(np.asarray(power, dtype=float).copy())
    return len(powers)-1


def box_vertices(d):
    return np.array(list(product((0., 1.), repeat=d)))


def clip_box(poly):
    """与评价盒 xi<=1 之交；xi>=0 由锥保证。"""
    for e in np.eye(np.shape(poly)[1]):
        poly = clip_polytope(poly, 1., -e)
    return poly


def cone_clip(poly, U):
    """与锥 cone(U)={U⁻¹xi>=0} 之交。"""
    for row in np.linalg.inv(U):
        poly = clip_polytope(poly, 0., row/np.linalg.norm(row))
    return poly


def cone_outer(U, halfspace):
    """锥的外域（未裁盒）：有外界 (c, mu) 时为 conv(0, mu/(c@u_i)*u_i)，否则为锥∩盒。"""
    d = len(U)
    if halfspace is None:
        return clip_box(np.vstack([np.zeros(d), (U*10.*np.sqrt(d)).T]))
    c, mu = halfspace
    return np.vstack([np.zeros(d), (U*(mu/(c@U))).T])


def stagnated(history, threshold, patience):
    """主线停滞规则：最近 patience 次割的体积缩减比例 area_ratio 都 < threshold。"""
    return len(history) >= patience and all(ratio < threshold for ratio in history[-patience:])


def piece(vertices, U, halfspace):
    """多面体与锥外域 O∩盒 之交（xi）；近退化（零体积）时为空。"""
    d = len(U)
    try:
        poly = cone_clip(vertices, U)
        if halfspace is not None and len(poly):
            poly = clip_polytope(poly, halfspace[1], -halfspace[0])
        poly = clip_box(poly)
        return poly if polytope_volume(poly) > 0. else np.empty((0, d))
    except QhullError:
        return np.empty((0, d))


def sandwich(cones, sets, d, cache):
    """内外测度（xi^d）：I=(K^IN ∪ ∪sets)∩K^OUT 与 K^OUT。cones 为叶锥几何 (U, V, halfspace)，sets 为
    [(版本键, 顶点)]。二维用 shapely 并集；三维按锥分解（叶锥内部互不相交），每锥求并集体积并按
    (锥几何, 相交集合的版本) 缓存；并集数值失败时只计 T（内域只会低估，间隙只会偏大）。"""
    # 1. 各叶锥的外域 O∩盒 与内三角形 T
    outers = [clip_box(cone_outer(U, halfspace)) for U, _, halfspace in cones]
    triangles = [np.vstack([np.zeros(d), V]) for _, V, _ in cones]
    # 2. 二维：内域并集与外域并集取交
    if d == 2:
        inner = polygon_union([*triangles, *(vertices for _, vertices in sets)]).intersection(polygon_union(outers))
        return inner.area, sum(polytope_volume(poly) for poly in outers)
    # 3. 三维：逐锥取与该锥相交的集合，T 与各集合裁到 锥∩O∩盒 后求并集体积
    inner = outer = 0.
    for (U, V, halfspace), triangle, poly in zip(cones, triangles, outers):
        rows = np.linalg.inv(U)
        touching = [(key, vertices) for key, vertices in sets
                    if not np.any(np.all(vertices@rows.T < -1e-12, axis=0))]
        key = (U.tobytes(), V.tobytes(), halfspace and (halfspace[0].tobytes(), halfspace[1]),
               tuple(key for key, _ in touching))
        if key not in cache:
            pieces = [triangle, *(p for p in (piece(vertices, U, halfspace) for _, vertices in touching) if len(p))]
            try:
                cache[key] = union_measure(pieces, d)
            except (QhullError, GEOSException):
                cache[key] = polytope_volume(triangle)
        inner += cache[key]
        outer += polytope_volume(poly)
    return inner, outer


# ---- 径向锥夹逼（R） ----------------------------------------------------------------------------
@dataclass(eq=False)
class Cone:
    """单纯锥 K=cone(u_keys)：内域 T=conv(0, verts)，外界 c@xi<=mu；尚未求解时沿用父锥外界 inherited。"""
    id: int
    keys: tuple                      # d 个方向编号
    x: tuple                         # 内域网架 x̂
    verts: np.ndarray                # (d,d)，第 i 行 v_i=rho_x̂(u_i)·u_i（xi）
    cover: tuple | None = None       # 近端覆盖网架 y；x̂ 原点可行时为 None
    inherited: tuple | None = None   # 父锥外界 (c, mu)
    mu: float = np.inf               # 外界因子 μ_k
    final: dict | None = None        # 最后一轮锥 MISOCP（细分用解点与方案）
    status: str = 'pending'          # pending / bounded / unresolved

    @property
    def c(self):
        """c=1ᵀQ⁻¹，Q=verts.T：锥内 xi=Σα_i v_i 时 Σα_i=c@xi。"""
        return np.linalg.solve(self.verts, np.ones(len(self.keys)))

    def halfspace(self):
        return (self.c, self.mu) if np.isfinite(self.mu) else self.inherited


class Radial:
    """一个符号分区的径向夹逼：方向表、射线与原点认证缓存、x 自由锥 MISOCP 与懒惰 OBBT、体积准则细分。"""

    def __init__(self, network, sign, budget, monitor, settings):
        # 1. 分区物理：非负幅值 u 的端口方程、固定 x,p 的 SP 与公共尺度 bounds（kW）
        self.network, self.budget, self.monitor, self.settings = deepcopy(network), budget, monitor, settings
        self.equations = PortPhysics(self.network, np.asarray(sign))
        self.oracle = PortSubProblem(self.equations, threads=settings.threads)
        self.bounds = port_bounds(self.network)
        self.d = len(self.bounds)
        self.switches = [s for c, s in zip(self.network.corridors, self.network.type_slices) if c.switchable]
        # 2. 方向表：初始为坐标轴；棱中点方向由相邻锥共用，射线缓存随之共用
        self.directions, self.midpoints = list(np.eye(self.d)), {}
        # 3. 证书缓存：射线远端 / 近端（键为 (方案, 方向编号)）与原点认证
        self.rays, self.nears, self.origin = {}, {}, {}
        self.cones, self.next_id = [], 0
        self.deadline, self.phase = np.inf, '径向搜索'

    def label(self, x):
        """网架标签：可变走廊依次写所选型号序号（0 为断开）；Case33 即七个开关的通断串。"""
        x = np.asarray(x)
        return ''.join(str(int(np.argmax(x[s]))+1) if x[s].any() else '0' for s in self.switches)

    def remaining(self, limit):
        """本阶段剩余时间（不超过 limit）；到时抛出 RegionTimeout，已得证书保留。"""
        left = self.deadline-self.monitor.clock()
        if left <= 0.:
            raise RegionTimeout('分区阶段时限已到')
        return min(limit, left)

    def U(self, keys):
        return np.column_stack([self.directions[k] for k in keys])

    def direction(self, vector):
        self.directions.append(np.asarray(vector, float)/np.linalg.norm(vector))
        return len(self.directions)-1

    def midpoint(self, a, b):
        key = min(a, b), max(a, b)
        if key not in self.midpoints:
            self.midpoints[key] = self.direction(self.directions[a]+self.directions[b])
        return self.midpoints[key]

    def step(self, kind, text, **geometry):
        """记录一步求解为 step 帧：kind、说明文字与几何（kW 幅值：点 p、折线 vertices；点所属网架 scheme），回放时
        逐步画出；判据步另带夹逼间隙 gap。"""
        self.monitor._emit('step', phase=self.phase, step=dict(kind=kind, text=text, **geometry))

    # 证书：紧化、原点认证、射线顶点与近端覆盖
    def tighten(self, x):
        """方案首次使用前 OBBT：其后该方案的 SP、射线与锥 MISOCP 都带它的盒与包络行。线程数为
        workers // 仍在计算的分区数，先结束的分区让出线程。"""
        if x not in self.equations.boxes:
            self.remaining(np.inf)
            threads = max(1, self.settings.workers//self.monitor.running())
            self.equations.obbt(np.asarray(x), self.bounds, threads)
            self.step('obbt', f'OBBT 紧化网架 {self.label(x)}（{threads} 线程）：分区盒内求选中型号 P/Q/ell 与节点电压的上下界')

    def solved(self, solve):
        """单次射线 / SP：时限取 RAY_SECONDS 与剩余时间的较小者；不可行或数值失败返回 None，阶段到时抛出。"""
        limit = self.remaining(RAY_SECONDS)
        try:
            return solve(limit)
        except TimeoutError:
            self.remaining(np.inf)
            return None
        except RuntimeError:
            return None

    def origin_feasible(self, x):
        """零接入 p=0 的紧化 SP 认证，按方案缓存。"""
        if x not in self.origin:
            self.tighten(x)
            answer = self.solved(lambda limit: self.oracle.solve(np.asarray(x), np.zeros(self.d), time_limit=limit,
                                                                 score_only=True))
            self.origin[x] = answer is not None and bool(answer['feasible'])
            self.step('origin', f"零接入 SP：固定网架 {self.label(x)}、p=0，min η → "
                      f"{'可行' if self.origin[x] else '不可行'}", p=np.zeros(self.d), scheme=self.label(x))
        return self.origin[x]

    def vertex(self, x, key):
        """v=rho_x(u)·u（xi）：从原点朝分区盒边界点 bounds*u/max(u) 的紧化射线；半径过小视为没有顶点。"""
        if (x, key) not in self.rays:
            self.tighten(x)
            u = self.directions[key]
            answer = self.solved(lambda limit: ray_support(self.equations, self.budget, np.asarray(x), np.zeros(self.d),
                                                           self.bounds*u/u.max(), threads=self.settings.threads,
                                                           time_limit=limit))
            point = None if answer is None else answer['p']/self.bounds
            self.rays[x, key] = v = point if point is not None and np.linalg.norm(point) >= MIN_RADIUS else None
            self.step('ray', f'射线：固定网架 {self.label(x)}，max t  s.t. t·u 对该网架可行（紧化 SOCP）'
                      +(' → 无可行顶点' if v is None else ' → 射线顶点 v'),
                      **({} if v is None else dict(p=v*self.bounds, vertices=np.vstack([np.zeros(self.d), v*self.bounds]),
                                                   scheme=self.label(x))))
        return self.rays[x, key]

    def near(self, x, key):
        """原点不可行方案沿 u 的可行段近端（xi）：从远端朝原点的反向紧化射线。"""
        if (x, key) not in self.nears:
            far = self.vertex(x, key)
            answer = None if far is None else self.solved(
                lambda limit: ray_support(self.equations, self.budget, np.asarray(x), far*self.bounds, np.zeros(self.d),
                                          threads=self.settings.threads, time_limit=limit))
            self.nears[x, key] = near = None if answer is None else answer['p']/self.bounds
            if far is not None:
                self.step('near', f'反向射线：网架 {self.label(x)} 原点不可行，从远端朝原点求可行段的近端'
                          +(' → 无解' if near is None else ''),
                          **({} if near is None else dict(p=near*self.bounds, vertices=np.vstack([far, near])*self.bounds,
                                                          scheme=self.label(x))))
        return self.nears[x, key]

    def inner(self, x, keys, covers):
        """网架 x 在锥 keys 上的内顶点与证书，不满足时为 None。
        原点可行时 T=conv(0, v) ⊆ R_x；否则需原点可行的覆盖网架 y（先试 covers，再试其他原点可行方案）在
        每个生成方向的半径不小于 x 的近端，T ⊆ R_x ∪ R_y。"""
        verts = [self.vertex(x, k) for k in keys]
        if any(v is None for v in verts):
            return None
        if self.origin_feasible(x):
            return dict(x=x, verts=np.array(verts), cover=None)
        radii = [self.near(x, k) for k in keys]
        if any(r is None for r in radii):
            return None
        radii = np.linalg.norm(radii, axis=1)
        for y in dict.fromkeys(y for y in (*covers, *self.origin) if y is not None):
            if self.origin_feasible(y):
                reach = [self.vertex(y, k) for k in keys]
                if all(r is not None for r in reach) and np.all(np.linalg.norm(reach, axis=1) >= radii):
                    return dict(x=x, verts=np.array(verts), cover=y)
        return None

    # 外界：x 自由锥 MISOCP（行生成）
    def misocp(self, objective, rows=None, exclude=()):
        """x 自由 MISOCP：一次分支定界，新现任网架尚未紧化时先 OBBT，其紧化行以惰性约束加入（model.cone_misocp）。"""
        return cone_misocp(self.equations, self.budget, self.bounds, objective, tighten=self.incumbent, rows=rows,
                           exclude=exclude, time_limit=self.remaining(self.settings.mip_seconds),
                           mip_gap=0. if rows is None else self.settings.mip_gap, threads=self.settings.threads)

    def incumbent(self, x, point, bound):
        """行生成回调：记一步，并确保分支定界的新现任网架已紧化（必要时 OBBT）。"""
        self.step('lazy', f'行生成：分支定界的现任网架 {self.label(x)} 尚无紧化行（当前上界 {bound:.4f}）→ '
                          '紧化并以惰性约束加入，分支定界继续', p=point, scheme=self.label(x))
        self.tighten(x)

    def initial(self):
        """中心射线 MISOCP（x 自由，ξ_1=…=ξ_d）选初始网架：胜出方案原点不可行时以 no-good 排除后重解。
        返回 (首个胜出方案, 原点可行的胜出方案)。"""
        excluded = []
        while True:
            answer = self.misocp(np.ones(self.d), exclude=excluded)
            if answer['x'] is None:
                self.remaining(np.inf)
                raise RuntimeError(f'中心射线 MISOCP 没有可行网架（status={answer["status"]}）')
            x, point = answer['x'], answer['point']
            self.step('center', f"中心射线 MISOCP：max Σξ  s.t. ξ_1=…=ξ_d，网架自由 → 网架 {self.label(x)}，"
                      f"上界 {answer['bound']:.3f}"+(f'（已排除 {len(excluded)} 个原点不可行网架）' if excluded else ''),
                      p=point, vertices=np.vstack([np.zeros(self.d), point]), scheme=self.label(x))
            if self.origin_feasible(x):
                return (excluded or [x])[0], x
            excluded.append(x)

    def new_cone(self, keys, inner, parent=None):
        cone = Cone(self.next_id, tuple(keys), inner['x'], inner['verts'], inner['cover'])
        self.next_id += 1
        if parent is not None and np.isfinite(parent.mu):
            cone.inherited = parent.c, parent.mu
        return cone

    def best(self, candidates, keys, covers):
        """候选网架中满足内域证书、且各顶点半径乘积（与 vol(T) 成正比）最大者。"""
        inners = [inner for inner in (self.inner(x, keys, covers) for x in dict.fromkeys(candidates)) if inner]
        return max(inners, key=lambda inner: float(np.prod(np.linalg.norm(inner['verts'], axis=1))), default=None)

    def root(self):
        """根锥 cone(e_1..e_d)：内域网架取首个胜出方案与原点可行胜出方案中证书成立、体积最大者。"""
        first, chosen = self.initial()
        inner = self.best((first, chosen), tuple(range(self.d)), (chosen,))
        if inner is None:
            raise RuntimeError(f'初始网架 {self.label(chosen)} 在坐标轴方向没有射线顶点')
        return self.new_cone(tuple(range(self.d)), inner)

    def far_face(self, cone):
        """锥外块 K^OUT_k 的远端面 c@xi=mu_k（kW 幅值），回放中标出锥 MISOCP 的结果。"""
        return cone.mu*cone.verts*self.bounds

    def solve(self, cone):
        """锥 MISOCP：锥约束 Q⁻¹xi>=0（行归一化）、分区盒，max c@xi；μ_k=ObjBound，不可行时取 1。"""
        rows = np.linalg.inv(cone.verts.T)
        answer = self.misocp(cone.c, rows/np.linalg.norm(rows, axis=1, keepdims=True))
        cone.mu = 1. if answer['status'] == GRB.INFEASIBLE else answer['bound']
        cone.final, cone.status = answer, 'bounded'
        found = ('锥内不可行' if answer['status'] == GRB.INFEASIBLE else '尚无解点' if answer['x'] is None
                 else f"解点网架 {self.label(answer['x'])}")
        geometry = {**({} if answer['point'] is None else dict(p=answer['point'], scheme=self.label(answer['x']))),
                    **(dict(vertices=self.far_face(cone)) if np.isfinite(cone.mu) else {})}
        self.publish([cone], step=dict(kind='cone', text=f'锥 {cone.id} MISOCP：max c·ξ  s.t. ξ∈K_{cone.id}，网架自由 → '
                                                         f'μ={cone.mu:.3f}，{found}', **geometry))

    # 细分
    def options(self, cone):
        """候选剖分（按优先次序）：解点方向 u* 的锥坐标 λ 全部 >= SPLIT_MARGIN 时星形剖分；三维恰有一个分量偏小
        时在对边上按 λ 投影二分（投影过近端点取棱中点）；最后总有最长棱中点二分。"""
        keys, U, d = cone.keys, self.U(cone.keys), self.d
        replace = lambda i, m: keys[:i]+(m,)+keys[i+1:]
        options = []
        point = cone.final['point']
        lam = np.zeros(d) if point is None else np.maximum(np.linalg.solve(U, point/self.bounds), 0.)
        if lam.sum() > 0.:
            lam = lam/lam.sum()
            near = np.flatnonzero(lam >= SPLIT_MARGIN)
            if len(near) == d:
                m = self.direction(point/self.bounds)
                options.append([replace(i, m) for i in range(d)])
            elif len(near) == d-1 >= 2:
                b, c = near
                share = lam[b]/(lam[b]+lam[c])
                m = (self.direction(lam[b]*U[:, b]+lam[c]*U[:, c]) if SPLIT_MARGIN < share < 1.-SPLIT_MARGIN
                     else self.midpoint(keys[b], keys[c]))
                options.append([replace(c, m), replace(b, m)])
        angles = {(a, b): float(np.arccos(np.clip(U[:, a]@U[:, b], -1., 1.))) for a, b in combinations(range(d), 2)}
        a, b = max(angles, key=angles.get)
        m = self.midpoint(keys[a], keys[b])
        options.append([replace(b, m), replace(a, m)])
        return options

    def split(self, cone):
        """子锥 x̂ 取 {父 x̂, 解的方案 x*, 父覆盖网架 y} 中证书成立、体积最大者；第一个全部子锥都成立的剖分生效。"""
        candidates = [cone.x, cone.final['x'], cone.cover]
        for groups in self.options(cone):
            inners = [self.best([x for x in candidates if x is not None], keys, (cone.cover, cone.x)) for keys in groups]
            if all(inners):
                return [self.new_cone(keys, inner, cone) for keys, inner in zip(groups, inners)]
        return None

    def delta(self, cone):
        """体积缺口 Δ_k=(μ_k^d-1)·vol(K^IN_k)（xi^d）：未求解的锥取父锥外界在本锥生成方向上的最大放大倍数，均无为 inf。"""
        halfspace = cone.halfspace()
        if halfspace is None:
            return np.inf
        U = self.U(cone.keys)
        factor = float(np.max(halfspace[1]*(cone.c@U)/(halfspace[0]@U)))
        return max(factor**self.d-1., 0.)*abs(float(np.linalg.det(cone.verts)))/factorial(self.d)

    def volume_ratio(self):
        """体积缺口比 ΣΔ/Σvol(T)；尚无锥时为 inf。"""
        total = sum(abs(float(np.linalg.det(c.verts))) for c in self.cones)/factorial(self.d)
        return sum(self.delta(c) for c in self.cones)/total if total > 0. else np.inf

    def run(self, epsilon, check=None):
        """体积准则细分到 ΣΔ<=ε·Σvol(T)，或 check() 判定认证；可续跑。返回 certified / unresolved / time_limit。"""
        try:
            # 1. 首次运行：中心射线选初始网架，建根锥并求外界
            if not self.cones:
                self.cones = [self.root()]
                self.publish(self.cones, step=dict(kind='root', text=f'根锥 K_0：内域网架 {self.label(self.cones[0].x)}，'
                                                   '坐标轴方向的射线顶点围成 K^IN_0'))
            # 2. 续跑：先求上次时限到达时尚未求解的锥（其间沿用父锥外界）
            for cone in [c for c in self.cones if c.status == 'pending']:
                self.solve(cone)
            heap = [(-self.delta(c), c.id, c) for c in self.cones if c.status == 'bounded']
            heapq.heapify(heap)
            while True:
                # 3. 认证：径向体积准则，或构域器给出的夹逼判据
                if self.volume_ratio() <= epsilon or (check is not None and check()):
                    return 'certified'
                if not heap:
                    return 'unresolved'
                # 4. 细分体积缺口最大的锥；角直径过小、锥数到上限或无内域证书的锥不再细分
                gap, _, cone = heapq.heappop(heap)
                self.remaining(np.inf)
                U = self.U(cone.keys)
                diameter = max(float(np.arccos(np.clip(U[:, a]@U[:, b], -1., 1.)))
                               for a, b in combinations(range(self.d), 2))
                children = (None if diameter < self.settings.min_width or len(self.cones) >= self.settings.max_cones
                            else self.split(cone))
                if children is None:
                    cone.status = 'unresolved'
                    continue
                # 5. 子锥替换父锥并立即求外界
                k = self.cones.index(cone)
                self.cones[k:k+1] = children
                self.publish(children, removed=[cone], step=dict(
                    kind='split', text=f"细分锥 {cone.id}（Δ={-gap:.3g} 为最大）→ 子锥 {'、'.join(str(c.id) for c in children)}，"
                                       f"内域网架 {'、'.join(self.label(c.x) for c in children)}"))
                for child in children:
                    self.solve(child)
                    heapq.heappush(heap, (-self.delta(child), child.id, child))
        except RegionTimeout:
            return 'time_limit'

    # 输出
    def geometry(self):
        """叶锥几何 [(U, V, 外界)]，供测度与结果使用（xi）。"""
        return [(self.U(c.keys), c.verts, c.halfspace()) for c in self.cones]

    def publish(self, cones, removed=(), step=None):
        """记录锥的变化：K^IN_k 与 K^OUT_k∩盒（kW 幅值）、x̂ 与 μ_k，被细分的父锥置空；附 ΣΔ/Σvol(K^IN)。"""
        rows = {str(c.id): None for c in removed}
        for c in cones:
            rows[str(c.id)] = dict(inner=np.vstack([np.zeros(self.d), c.verts])*self.bounds,
                                   outer=clip_box(cone_outer(self.U(c.keys), c.halfspace()))*self.bounds,
                                   scheme=self.label(c.x), mu=c.mu)
        self.monitor._emit('cone', phase=self.phase, cones=rows, volume_ratio=self.volume_ratio(), step=step)

    def schemes(self):
        """X*：可能改变 I 的网架。先取叶锥的 x̂（按所占锥体积从大到小），再接各叶锥的覆盖网架与锥 MISOCP 解点网架
        （按所在锥的 Δ_k 从大到小）；只做过 OBBT 或原点测试的网架不在其中。"""
        volumes = {}
        for cone in self.cones:
            volumes[cone.x] = volumes.get(cone.x, 0.)+abs(float(np.linalg.det(cone.verts)))
        others = [y for cone in sorted(self.cones, key=lambda c: -self.delta(c))
                  for y in (cone.cover, cone.final and cone.final['x']) if y is not None]
        return list(dict.fromkeys([*sorted(volumes, key=lambda x: -volumes[x]), *others]))


# ---- 逐网架构域：主线割平面（CUT） ---------------------------------------------------------------
class CutSlice(Exception):
    """一个网架的割循环用完时间片。"""


class Cutting:
    """逐网架割平面构域；分区内的联合割共享：新网架的 N^CUT_x 先被已有割裁剪，新割也裁剪其他网架的 N^CUT_x
    （割对全部可行点有效，N^CUT_x 仍包含 R^SOCP_x）。"""

    def __init__(self, radial):
        self.radial, self.settings, self.monitor = radial, radial.settings, radial.monitor
        self.epsilon = radial.d*self.settings.tau
        self.cuts, self.networks, self.count = [], {}, 0   # 共享割 α+βᵀu+δᵀx>=0；方案 → 网架状态；割序号
        self.declined, self.measures = set(), {}           # 某轮因间隙降幅不足而不再割的网架；sandwich 逐锥测度缓存

    def clip(self, vertices, x, cuts):
        bounds, d = self.radial.bounds, self.radial.d
        for cut in cuts:
            vertices = clip_polytope(vertices, cut[0]+cut[1+d:]@np.asarray(x), cut[1:1+d]*bounds)
        return vertices

    def row(self, x):
        """网架面板的一行（kW 幅值）：N^CUT_x 与状态；状态在 CUT_ACCEPTED 中时 N^CUT_x 计入 I。"""
        state, radial = self.networks[x], self.radial
        return dict(x=x, choice=radial.network.decode_plan(np.asarray(x)),
                    cost=float(radial.network.cost_offset+radial.network.cost@np.asarray(x)),
                    outer=state['vertices']*radial.bounds, status=state['status'])

    def sets(self):
        """计入 N^CUT 的 N^CUT_x：[(版本键, 顶点)]（xi）。"""
        return [((x, state['version']), state['vertices']) for x, state in self.networks.items()
                if state['status'] in CUT_ACCEPTED and len(state['vertices']) > self.radial.d]

    def gap(self):
        """当前夹逼间隙 vol(K^OUT)/vol(I)-1；I 为空时为 inf。"""
        inner, outer = sandwich(self.radial.geometry(), self.sets(), self.radial.d, self.measures)
        return outer/inner-1. if inner > 0. else np.inf

    def run(self, schemes):
        """一轮 CUT：间隙已达 ε 时不割；否则按 X* 的次序割未割过的网架。本轮最多用剩余时间的 PASS_SHARE，每个网架不超过
        CUT_SECONDS；超过本轮时限仍未开始的网架记 skipped，下一轮重试。"""
        radial, monitor = self.radial, self.monitor
        before = self.gap()
        if before <= self.epsilon:
            return
        end = monitor.clock()+PASS_SHARE*max(radial.deadline-monitor.clock(), 0.)
        queue = [x for x in schemes
                 if x not in self.declined and (x not in self.networks or self.networks[x]['status'] == 'skipped')]
        for k, x in enumerate(queue):
            if monitor.clock() >= end:
                self.networks[x] = dict(vertices=box_vertices(radial.d), status='skipped', version=0)
                continue
            radial.tighten(x)
            self.cut(x, min(CUT_SECONDS, end-monitor.clock()))
            # 每割完一个网架求间隙：已达 ε 或降幅不足 GAP_SHARE·ε 即结束本轮，本轮其余网架之后不再割
            after = self.gap()
            if after <= self.epsilon or before-after < GAP_SHARE*self.epsilon:
                self.declined.update(queue[k+1:])
                reason = '已达 ε' if after <= self.epsilon else f'降幅不足 {GAP_SHARE:g}ε'
                radial.step('check', f'CUT 本轮结束：割完网架 {radial.label(x)} 后间隙 {before:.2%} → {after:.2%}，'
                                     f'{reason}；本轮其余 {len(queue)-k-1} 个网架不再割', gap=after)
                return
            before = after

    def cut(self, x, limit):
        """网架 x 的割平面循环。"""
        radial, settings, monitor = self.radial, self.settings, self.monitor
        d, bounds, label = radial.d, radial.bounds, radial.label(x)
        deadline = monitor.clock()+limit

        def left():
            rest = deadline-monitor.clock()
            if rest <= 0.:
                raise CutSlice()
            return radial.remaining(min(rest, SP_TIME_LIMIT['socp']))

        # 1. N^CUT_x 初值：conv(K^OUT) ∩ 盒（K^OUT ⊇ R^SOCP ⊇ R^SOCP_x）被已有联合割裁剪；先记录初值
        start = polytope_vertices(np.vstack([clip_box(cone_outer(U, halfspace)) for U, _, halfspace in radial.geometry()]))
        state = self.networks[x] = dict(vertices=self.clip(start, x, self.cuts), status='slice', version=0)
        monitor._emit('network', phase='网架切割', active_scheme=label, schemes={label: self.row(x)},
                      step=dict(kind='network', text=f'CUT 网架 {label}：N^CUT_x 从 conv(K^OUT) ∩ 盒出发，'
                                                    f'先被已有 {len(self.cuts)} 刀共享割裁剪'))
        cache, powers, applied, failed, history, fresh = {}, [], set(), set(), [], []

        def score(index):
            """SP 评分（固定 x，带本网架 OBBT 盒与包络行，η 松弛）：可行 / 待割 / 该点数值失败不取割。"""
            if index in applied:
                return 'applied'
            if index in failed:
                return 'failed'
            if index not in cache:
                limit = left()
                try:
                    cache[index] = radial.oracle.solve(np.asarray(x), powers[index], time_limit=limit, score_only=True)
                except RuntimeError:
                    failed.add(index)
                    return 'failed'
                eta, feasible = cache[index]['eta'], cache[index]['feasible']
                monitor._emit('point', phase='网架切割', active_scheme=label,
                              step=dict(kind='sp', p=powers[index], scheme=label, feasible=feasible,
                                        text=f"SP：固定网架 {label} 与 N^CUT_x 的一个顶点，min η → "
                                             f"η={eta:.3g}（{'可行' if feasible else '待割'}）"))
            return 'feasible' if cache[index]['feasible'] else 'pending'

        try:
            while True:
                if not len(state['vertices']):
                    state['status'] = 'empty'
                    break
                # 2. SP 评分 N^CUT_x 的每个顶点（kW 最大坐标差 point_tol 内视为同一点）
                indices = [register_power(point*bounds, powers, settings.point_tol) for point in state['vertices']]
                kinds = {index: score(index) for index in indices}
                pending = [index for index, kind in kinds.items() if kind == 'pending']
                if not pending:
                    values = set(kinds.values())
                    state['status'] = 'failed' if 'failed' in values else 'point_resolution' if 'applied' in values else 'exact'
                    break
                # 3. 违反量最大的顶点取联合割（锥支撑平面 LP 的对偶），裁剪 N^CUT_x
                index = max(pending, key=lambda i: (cache[i]['eta'], float(powers[i].sum()), tuple(powers[i])))
                limit = left()
                try:
                    cut = radial.oracle.generate_cut(np.asarray(x), powers[index], cache[index]['cone_normals'],
                                                     time_limit=limit)
                except RuntimeError:
                    failed.add(index)
                    continue
                before = polytope_volume(state['vertices'])
                state['vertices'] = self.clip(state['vertices'], x, [cut])
                applied.add(index)
                fresh.append(cut)
                history.append(1.-polytope_volume(state['vertices'])/before if before > 0. else 0.)
                self.count += 1
                small = next((k for k, ratio in enumerate(reversed(history)) if ratio >= settings.threshold), len(history))
                monitor._emit('cut', phase='网架切割', active_scheme=label, schemes={label: self.row(x)},
                              cut_history={str(self.count): dict(cut=cut, scheme=label)},
                              step=dict(kind='cut', p=powers[index], scheme=label,
                                        text=f"取割：η 最大的顶点（η={cache[index]['eta']:.3g}）→ 对偶 LP 得联合割 "
                                             f"#{self.count}，割掉 {history[-1]:.1%}，连续小割 {small}/{settings.patience}"))
                # 4. 连续 patience 次割的体积缩减比例 < threshold 即停滞
                if stagnated(history, settings.threshold, settings.patience):
                    state['status'] = 'stagnated'
                    break
        except (CutSlice, TimeoutError):
            radial.remaining(np.inf)
        # 5. 新割裁剪其他网架的 N^CUT_x，记录本网架终态
        self.cuts.extend(fresh)
        rows = {}
        for y, other in self.networks.items():
            if y != x and fresh and len(other['vertices']):
                other['vertices'], other['version'] = self.clip(other['vertices'], y, fresh), other['version']+1
                rows[radial.label(y)] = self.row(y)
        rows[label] = self.row(x)
        monitor._emit('network', phase='网架切割', active_scheme=None, schemes=rows,
                      step=dict(kind='network_end', text=f"网架 {label} 停止：{state['status']}，本网架 {len(fresh)} 刀"))


# ---- 分区与全局 ---------------------------------------------------------------------------------
def build_partition(network, sign, budget, monitor, settings):
    """一个符号分区的 RCUT 构域。"""
    # 1. 径向结构与割平面构域器；记录分区初始帧
    start = monitor.clock()
    end = start+monitor.time_limit
    radial = Radial(network, sign, budget, monitor, settings)
    cutting = Cutting(radial)
    epsilon, d = cutting.epsilon, radial.d
    monitor._emit('phase_start', phase='径向搜索', status='running', cones={}, schemes={}, cut_history={})

    def certified():
        """夹逼判据：vol(K^OUT)-vol(I)<=ε·vol(I)，I=(K^IN ∪ N^CUT)∩K^OUT。"""
        gap = cutting.gap()
        radial.step('check', f"判据：vol(K^OUT)/vol(I)−1 = {gap:.2%}，ε = {epsilon:.1%} → {'认证' if gap <= epsilon else '未满足'}",
                    gap=gap)
        return gap <= epsilon

    # 2. A：径向全局搜索到放宽目标 ε_A，时限为分区时限的 share_A
    radial.phase, radial.deadline = '径向搜索 A', start+settings.discovery_share*(end-start)
    status = radial.run(settings.discovery_eps)
    radial.deadline = end
    try:
        if radial.volume_ratio() > epsilon:
            # 3. CUT：A 中出现的网架割到停滞
            cutting.run(radial.schemes())
            # 4. A+：续跑径向部分到 ε；每次锥决策后先割新网架，再检查夹逼判据
            radial.phase = '径向续跑 A+'
            status = 'certified' if certified() else radial.run(
                epsilon, lambda: cutting.run(radial.schemes()) or certified())
        else:
            status = 'certified'
    except RegionTimeout:
        status = 'time_limit'

    # 5. 汇总：I 为各 K^IN_k 与按锥裁到 K^OUT_k 的 N^CUT_x，外域为 K^OUT；尚无锥时外域为分区盒（kW 幅值）
    gap = cutting.gap()
    rows = dict(inner=[], outer=[] if radial.cones else [box_vertices(d)])
    for U, V, halfspace in radial.geometry():
        rows['inner'] += [np.vstack([np.zeros(d), V]),
                          *(p for p in (piece(vertices, U, halfspace) for _, vertices in cutting.sets()) if len(p))]
        rows['outer'].append(clip_box(cone_outer(U, halfspace)))
    statuses = [state['status'] for state in cutting.networks.values()]
    result = dict(sign=list(sign), status=status, certified=status == 'certified',
                  how='radial' if radial.volume_ratio() <= epsilon else 'cut', seconds=monitor.clock()-start,
                  gap=gap if np.isfinite(gap) else None, volume_ratio=radial.volume_ratio(),
                  cones=len(radial.cones), networks=len(statuses), accepted=sum(s in CUT_ACCEPTED for s in statuses),
                  cuts=cutting.count, schemes=list(dict.fromkeys([*radial.schemes(), *cutting.networks])),
                  inner=[r*radial.bounds for r in rows['inner']],
                  outer=[r*radial.bounds for r in rows['outer']])
    monitor._emit('partition_end', phase='分区完成' if result['certified'] else '分区停止', status=status,
                  active_scheme=None, step=dict(kind='end', gap=result['gap'], text=f"分区结束：{status}，{result['cones']} 个锥、"
                                                                f"{result['cuts']} 刀，间隙 {result['gap']:.2%}"
                                                if result['gap'] is not None else f'分区结束：{status}'))
    return result


def _partition(network, sign, budget, deadline, settings):
    """子进程：一个符号分区的构域；事件经队列送回主进程，结束时发送结束标记。"""
    monitor = RunMonitor(sign=sign)
    monitor.time_limit = deadline-time.time()
    try:
        with threadpool_limits(limits=1):
            return build_partition(network, sign, budget, monitor, settings)
    finally:
        monitor.close()


def build_region(network, *, budget, monitor, seconds, workers, settings):
    """符号分区并行构域，合并为带符号 kW 的内外域。"""
    # 1. 全局初始帧：分区、公共尺度与网络信息
    d = len(network.load_nodes)
    signs = list(product((1, -1), repeat=d))
    labels = [''.join('+' if s > 0 else '-' for s in sign) for sign in signs]
    monitor.time_limit = seconds
    monitor._emit('start', phase='初始化', network=network.name, load_nodes=network.load_nodes, budget=budget,
                  algorithm=monitor.algorithm, time_limit=seconds, cost_unit=network.cost_unit,
                  initial_plan=network.initial_plan, bounds=port_bounds(network), partitions=labels,
                  cones={}, schemes={}, cut_history={}, status='running')
    # 2. 子进程池：每个分区一个任务；并发数不足时后启动的分区得到剩余时间按并发比例的份额
    concurrency = min(workers, len(signs))
    settings.workers = workers
    end = time.time()+seconds
    channel = monitor.share()
    queue, results = list(zip(signs, labels)), {}
    pool = ProcessPoolExecutor(concurrency, mp_context=multiprocessing.get_context('spawn'),
                               initializer=_connect, initargs=(channel,))
    with pool:
        futures = {}

        def submit():
            sign, label = queue.pop(0)
            now = time.time()
            share = min(1., concurrency/(len(queue)+len(futures)+1))
            futures[label] = pool.submit(_partition, network, sign, budget, now+max(end-now, 0.)*share, settings)
            channel[4].value += 1

        while queue and len(futures) < concurrency:
            submit()
        # 3. 转发子进程事件，直到全部分区发来结束标记；仍在计算的分区数供子进程分配 OBBT 线程；
        #    中断时通知子进程停止并读空队列（子进程退出前要送完事件）
        try:
            while futures:
                label, patch = channel[0].get()
                if patch is not None:
                    monitor.forward(label, patch)
                    continue
                results[label] = futures.pop(label).result()
                channel[4].value -= 1
                if queue:
                    submit()
        except BaseException:
            channel[3].set()
            while futures:
                label, patch = channel[0].get()
                if patch is None:
                    futures.pop(label)
            raise
    # 4. 合并：带符号 kW 的内外域与分区摘要
    partitions = [results[label] for label in labels]
    rows = {key: [dict(vertices=np.asarray(v)*np.asarray(p['sign']), sign=p['sign']) for p in partitions for v in p[key]]
            for key in ('inner', 'outer')}
    vertices = np.vstack([row['vertices'] for row in rows['outer']])
    certified = all(p['certified'] for p in partitions)
    result = dict(status='certified' if certified else 'time_limit', certified=certified, timing=monitor.timing(),
                  axis_lower=vertices.min(axis=0), axis_bounds=vertices.max(axis=0),
                  partitions=[{k: v for k, v in p.items() if k not in ('inner', 'outer')} for p in partitions], **rows)
    for label, p in zip(labels, partitions):
        print(f'{network.name} {d}D {label}: {p["status"]} ({p["how"]}), {p["seconds"]:.1f}s, gap {p["gap"]}, '
              f'cones {p["cones"]}, networks {p["accepted"]}/{p["networks"]}', flush=True)
    monitor._emit('region_end', phase='构域完成' if certified else '构域停止', result=result, status=result['status'],
                  active_scheme=None, partition=None)
    return result
