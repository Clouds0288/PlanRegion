"""方案二（支撑查询）对照实验的纯几何部分：不建求解器，供 compare_plan2.py 与 tests/test_plan2_support.py 使用。

坐标均为归一化 xi=u/bounds（分区盒 [0,1]^d）。面方程沿用 region.halfspaces 的 [n, b]：内侧 n@xi+b<=0，
|n|=1；记 beta=-b，即 n@xi<=beta。两种面认证判据（h_x(n) 为固定方案 x 的紧化 SOCP 支撑值，UB 为其可靠上界）：

- origin：(1-tau)*UB+b <= epsilon_geom（原支持面实验判据），外扩集 E_x = P_x/(1-tau)；
- center：UB-beta <= tau*(beta-n@c_x)+epsilon_geom，c_x 为 P_x 内点，E_x = c_x+(1+tau)(P_x-c_x)。

auto：原点在 P_x 内且到每个非分区边界面的距离 >= ORIGIN_CLEARANCE 时取 origin，否则取 center。
记号登记见 docs/notation.md「方案二对照实验」。
"""
from __future__ import annotations

from itertools import product

import numpy as np
from scipy.optimize import linprog
from scipy.spatial import ConvexHull

from monitor import _union, comparison_metrics
from region import GEOMETRY_TOL, clip_polytope, contains, halfspaces, polytope_volume, union_measure

ORIGIN_CLEARANCE = 1e-6   # auto 判据：原点到非分区边界面的最小距离（xi）
BOUNDARY_TOL = 1e-7       # 分区边界面识别容差（法向与截距，xi）
FACE_TOL = 1e-9           # 顶点在面上的判定容差（xi）
ALPHA_TOL = 1e-12         # 锥内重心坐标 α>=0 的舍入容差，同径向脚本


# ---- 面方程与判据 --------------------------------------------------------------------------------
def box_vertices(d):
    return np.asarray(list(product((0., 1.), repeat=d)))


def boundary_faces(faces, tol=BOUNDARY_TOL):
    """分区边界面掩码：xi_j>=0（n=-e_j, b=0）或 xi_j<=1（n=e_j, b=-1）。这些面由分区盒本身认证。"""
    faces = np.asarray(faces, float)
    normal, offset = faces[:, :-1], faces[:, -1]
    mask = np.zeros(len(faces), bool)
    for e in np.eye(normal.shape[1]):
        mask |= (np.abs(normal+e).max(axis=1) <= tol) & (np.abs(offset) <= tol)
        mask |= (np.abs(normal-e).max(axis=1) <= tol) & (np.abs(offset+1.) <= tol)
    return mask


def face_measures(vertices, faces, tol=FACE_TOL):
    """各面的 (d-1) 维测度（xi）：二维为边长，三维为面多边形面积；低维集的成对等式面记 0。"""
    vertices, faces = np.asarray(vertices, float), np.asarray(faces, float)
    d = vertices.shape[1]
    measures = np.zeros(len(faces))
    for f, face in enumerate(faces):
        on = vertices[np.abs(vertices@face[:-1]+face[-1]) <= tol]
        if len(on) < d:
            continue
        basis = np.linalg.svd(np.atleast_2d(face[:-1]), full_matrices=True)[2][1:]   # 面内正交基
        coordinates = on@basis.T
        if d == 2:
            measures[f] = float(np.ptp(coordinates[:, 0]))
        elif len(on) >= 3 and np.linalg.matrix_rank(coordinates-coordinates[0], tol=1e-12) == 2:
            measures[f] = float(ConvexHull(coordinates).volume)
    return measures


def chebyshev_center(faces):
    """P={xi: n_f@xi<=beta_f}（|n_f|=1）的 Chebyshev 中心：max r s.t. n_f@c+r<=beta_f。退化或失败时为 None。"""
    faces = np.asarray(faces, float)
    if not len(faces):
        return None
    d = faces.shape[1]-1
    answer = linprog(np.r_[np.zeros(d), -1.], A_ub=np.c_[faces[:, :-1], np.ones(len(faces))], b_ub=-faces[:, -1],
                     bounds=[(None, None)]*d+[(0., None)], method='highs')
    if answer.status != 0 or not answer.x[-1] > 0.:
        return None
    return answer.x[:d], float(answer.x[-1])


def interior_point(vertices, faces):
    """P_x 的内点 c_x：Chebyshev 中心 LP，失败时取顶点形心。返回 (c_x, 来源, 内切半径)。"""
    center = chebyshev_center(faces)
    if center is not None:
        return center[0], 'chebyshev', center[1]
    return np.asarray(vertices, float).mean(axis=0), 'centroid', 0.


def choose_criterion(faces, boundary, mode):
    """auto：原点在 P_x 内，且到每个非分区边界面的距离 >= ORIGIN_CLEARANCE 时用 origin，否则用 center。"""
    if mode != 'auto':
        return mode
    beta = -np.asarray(faces, float)[:, -1]
    inside = bool(np.all(beta >= -GEOMETRY_TOL))
    clear = bool(np.all(beta[~boundary] >= ORIGIN_CLEARANCE))
    return 'origin' if inside and clear else 'center'


def allowed_offsets(faces, criterion, center, tau):
    """认证允许的支撑值 beta'_f（不含 epsilon_geom）：origin 为 beta/(1-tau)，center 为 beta+tau*(beta-n@c)。"""
    faces = np.asarray(faces, float)
    beta = -faces[:, -1]
    if criterion == 'origin':
        return beta/(1.-tau)
    return beta+tau*(beta-faces[:, :-1]@np.asarray(center, float))


def expanded_faces(faces, criterion, center, tau):
    """E_x 的面方程 [n, -beta']：origin 为 P_x/(1-tau)，center 为 c_x+(1+tau)(P_x-c_x)。"""
    faces = np.asarray(faces, float)
    return np.c_[faces[:, :-1], -allowed_offsets(faces, criterion, center, tau)]


def face_margins(faces, upper, criterion, center, tau):
    """各面认证裕量（<= epsilon_geom 即认证）：origin 为 (1-tau)*UB+b，center 为 UB-beta-tau*(beta-n@c)。"""
    faces, upper = np.asarray(faces, float), np.asarray(upper, float)
    if criterion == 'origin':
        return (1.-tau)*upper+faces[:, -1]
    return upper-allowed_offsets(faces, criterion, center, tau)


def classify_center(lb, ub, point, normal, offset, center, tau, epsilon_geom):
    """center 判据的支撑分类，语义同 classify_support：只有可靠上界能认证，只有审计点能证明越界。"""
    allowed = float(allowed_offsets(np.r_[normal, offset][None], 'center', center, tau)[0])
    rho_lb, rho_ub = lb-allowed, ub-allowed
    violated = point is not None and float(np.asarray(normal)@point)-allowed > epsilon_geom
    if violated and rho_ub <= epsilon_geom:
        raise RuntimeError('Inconsistent support upper bound and certified incumbent')
    if rho_ub <= epsilon_geom:
        return 'SUPPORT_CERTIFIED', rho_lb, rho_ub
    if violated:
        return 'VIOLATED', rho_lb, rho_ub
    return 'UNRESOLVED', rho_lb, rho_ub


def intersect_faces(poly, faces):
    """多面体（顶点）与 {faces 内侧} 之交；空集返回 (0,d)。"""
    for face in np.asarray(faces, float):
        poly = clip_polytope(poly, -face[-1], -face[:-1])
        if not len(poly):
            break
    return poly


def faces_polytope(faces, d):
    """{xi∈[0,1]^d : faces} 的顶点（xi）；空集返回 (0,d)。"""
    return intersect_faces(box_vertices(d), faces)


def ordered_polygon(vertices):
    """二维凸多边形顶点按逆时针排序（绘图用）。"""
    vertices = np.asarray(vertices, float)
    center = vertices.mean(axis=0)
    return vertices[np.argsort(np.arctan2(vertices[:, 1]-center[1], vertices[:, 0]-center[0]))]


# ---- 径向锥（从快照数据重建，二维与三维通用） ------------------------------------------------------
def cone_outer(U, halfspace):
    """锥 cone(U) 的外域多面体（xi，未裁盒）：有外界 (c, mu) 时为 conv(0, mu/(c@u_i)*u_i)，否则为锥∩盒。"""
    U = np.asarray(U, float)
    d = U.shape[0]
    if halfspace is None:
        poly = np.vstack([np.zeros(d), (U*10.*np.sqrt(d)).T])
        for j in range(d):
            poly = clip_polytope(poly, 1., -np.eye(d)[j])
        return poly
    c, mu = np.asarray(halfspace[0], float), float(halfspace[1])
    return np.vstack([np.zeros(d), (U*(mu/(c@U))).T])


def clip_box(poly):
    d = np.asarray(poly).shape[1]
    for j in range(d):
        poly = clip_polytope(poly, 1., -np.eye(d)[j])
    return poly


def cone_clip(poly, U):
    """多面体与锥 cone(U)={U^{-1}xi>=0} 之交。"""
    rows = np.linalg.inv(np.asarray(U, float))
    for row in rows:
        poly = clip_polytope(poly, 0., row/np.linalg.norm(row))
        if not len(poly):
            break
    return poly


def classify_cones(xi, cones):
    """格心（xi，本分区）的径向内域/外域标签。cones 为 dict(U, V, halfspace)：U、V 的列为方向与顶点。
    格心归锥坐标最小分量最大的锥；内域 α=V^{-1}xi>=0 且 Σα<=1，外域 c@xi<=mu（无外界时为评价盒）。"""
    xi = np.asarray(xi, float)
    inner, outer = np.zeros(len(xi), bool), np.ones(len(xi), bool)
    if not cones or not len(xi):
        return inner, outer
    best, owner = np.full(len(xi), -np.inf), np.zeros(len(xi), int)
    for k, cone in enumerate(cones):
        score = (xi@np.linalg.inv(np.asarray(cone['U'], float)).T).min(axis=1)
        better = score > best
        best[better], owner[better] = score[better], k
    for k, cone in enumerate(cones):
        chosen = owner == k
        if not chosen.any():
            continue
        if not cone.get('pseudo'):
            alpha = xi[chosen]@np.linalg.inv(np.asarray(cone['V'], float)).T
            inner[chosen] = np.all(alpha >= -ALPHA_TOL, axis=1) & (alpha.sum(axis=1) <= 1.)
        halfspace = cone['halfspace']
        outer[chosen] = (xi[chosen]@np.asarray(halfspace[0]) <= halfspace[1] if halfspace is not None
                         else np.all(xi[chosen] <= 1., axis=1))
    return inner, outer


# ---- 方法 H 的内外测度（xi^d；乘 prod(bounds) 得 kW^d） ---------------------------------------------
def h_measures(cones, inner_sets, cover_faces, d, keys=None, cache=None):
    """I_H = I_R ∪ (∪_x P_x)，O_H = (O_R∩盒) [∩ (∪_x E_x)]。cones 为径向叶锥 dict(U, V, halfspace)，inner_sets 为
    各 P_x 顶点，cover_faces 为覆盖证书冻结的各 E_x 面方程，覆盖未完成时为 None。二维用 shapely 并集；三维按锥
    分解（叶锥内部互不相交），每锥只对与之相交的多面体做 plot.union_volume，keys（各 P_x 的版本标识）与 cache
    给出时按 (锥几何, 相交的 P_x 版本) 缓存每锥测度。返回 (vol I_H, vol O_H)（xi^d）。"""
    triangles = [None if c.get('pseudo') else np.vstack([np.zeros(d), np.asarray(c['V'], float).T]) for c in cones]
    outers = [clip_box(cone_outer(c['U'], c['halfspace'])) for c in cones]
    chosen = [k for k, p in enumerate(inner_sets) if len(p) > d and polytope_volume(p) > 0.]
    inner_sets = [np.asarray(inner_sets[k], float) for k in chosen]
    keys = [None]*len(inner_sets) if keys is None else [keys[k] for k in chosen]
    if d == 2:
        inner = _union([*(t for t in triangles if t is not None), *inner_sets]).area
        if cover_faces is None:
            return inner, sum(polytope_volume(p) for p in outers)
        covers = [p for p in (faces_polytope(f, d) for f in cover_faces) if len(p) > d]
        return inner, (_union(outers).intersection(_union(covers)).area if covers else 0.)
    covers = None if cover_faces is None else [faces_polytope(f, d) for f in cover_faces]
    inner = outer = 0.
    for cone, triangle, poly in zip(cones, triangles, outers):
        rows = np.linalg.inv(np.asarray(cone['U'], float))
        touching = [k for k, p in enumerate(inner_sets) if meets_cone(p, rows)]
        geometry = (np.asarray(cone['U'], float).tobytes(), np.asarray(cone['V'], float).tobytes(),
                    None if cone['halfspace'] is None else (np.asarray(cone['halfspace'][0], float).tobytes(),
                                                            float(cone['halfspace'][1])))
        key = ('inner', geometry, tuple(keys[k] for k in touching))
        if cache is not None and None not in key[2] and key in cache:
            inner += cache[key]
        else:
            pieces = [cone_clip(inner_sets[k], cone['U']) for k in touching]
            pieces = [*([] if triangle is None else [triangle]), *(p for p in pieces if len(p) > d)]
            value = union_measure(pieces, d) if pieces else 0.
            if cache is not None and None not in key[2]:
                cache[key] = value
            inner += value
        if covers is None:
            outer += polytope_volume(poly)
            continue
        key = ('outer', geometry, len(cover_faces))
        if cache is not None and key in cache:
            outer += cache[key]
            continue
        pieces = [p for p in (intersect_faces(poly, faces) for faces, cover in zip(cover_faces, covers)
                              if len(cover) > d and meets_cone(cover, rows)) if len(p) > d]
        value = union_measure(pieces, d) if pieces else 0.
        if cache is not None:
            cache[key] = value
        outer += value
    return inner, outer


def meets_cone(vertices, rows, tol=1e-12):
    """快速排除：若某条锥约束 row@xi>=0 被全部顶点违反，多面体与锥不相交（否则保守地认为相交）。"""
    return not np.any(np.all(np.asarray(vertices, float)@rows.T < -tol, axis=0))


def h_geometry(cones, inner_sets, cover_faces):
    """二维绘图用的 shapely 几何（xi）：内域 I_R ∪ (∪P_x)，外界 (O_R∩盒) [∩ (∪E_x)]。"""
    triangles = [np.vstack([np.zeros(2), np.asarray(c['V'], float).T]) for c in cones if not c.get('pseudo')]
    outers = [clip_box(cone_outer(c['U'], c['halfspace'])) for c in cones]
    close = lambda geometry: geometry.buffer(1e-9).buffer(-1e-9)   # 相邻锥多边形共边的浮点细缝
    inner = close(_union([*triangles, *(p for p in inner_sets if len(p) > 2)]))
    outer = close(_union([p for p in outers if len(p) > 2]))
    if cover_faces is not None:
        covers = [p for p in (faces_polytope(f, 2) for f in cover_faces) if len(p) > 2]
        outer = outer.intersection(_union(covers))
    return inner, outer


# ---- 逐格指标与有效性 ----------------------------------------------------------------------------
def grid_metrics(inner, outer, ac_states, socp_states):
    """FR=|域∧非可行|/|域|，MR=|可行∧非域|/|可行|（comparison_metrics，百分数）；参考标签 0（未决）不计入该参考。
    AC 为主参考；SOCP 只用于区分误差来自松弛还是模型。"""
    inner, outer = np.asarray(inner, bool), np.asarray(outer, bool)
    ac_states, socp_states = np.asarray(ac_states), np.asarray(socp_states)
    known_ac, known_socp = ac_states != 0, socp_states != 0
    pair = lambda domain, feasible, known: comparison_metrics(domain[known], feasible[known])
    return dict(inner_ac=pair(inner, ac_states == 1, known_ac), inner_socp=pair(inner, socp_states == 1, known_socp),
                outer_ac=pair(outer, ac_states == 1, known_ac), outer_socp=pair(outer, socp_states == 1, known_socp),
                inner_cells=int(inner.sum()), outer_cells=int(outer.sum()),
                **validity(inner, outer, ac_states, socp_states))


def validity(inner, outer, ac_states, socp_states):
    """有效性（逐格的必要条件）：
    1. 内域不含 SOCP 已证不可行格——紧化模型 ⊆ SOCP 松弛，违反即内域证书失效；
    2. AC 可行格都在外界内——真实可行域 ⊆ 紧化模型 ⊆ 外界，违反即外界证书失效；
    3. AC 未决格（0）单列，不计为不可行。SOCP 可行而落在外界外的格是 OBBT 去掉的 SOCP-only 区域，只作诊断。"""
    inner, outer = np.asarray(inner, bool), np.asarray(outer, bool)
    ac_states, socp_states = np.asarray(ac_states), np.asarray(socp_states)
    inner_socp_infeasible = int(np.count_nonzero(inner & (socp_states == -1)))
    outer_missed_ac = int(np.count_nonzero((ac_states == 1) & ~outer))
    return dict(inner_socp_infeasible_cells=inner_socp_infeasible, outer_missed_ac_cells=outer_missed_ac,
                outer_missed_socp_cells=int(np.count_nonzero((socp_states == 1) & ~outer)),
                undecided_ac_cells=int(np.count_nonzero(ac_states == 0)),
                undecided_socp_cells=int(np.count_nonzero(socp_states == 0)),
                inner_undecided_ac_cells=int(np.count_nonzero(inner & (ac_states == 0))),
                valid=inner_socp_infeasible == 0 and outer_missed_ac == 0)


def network_cells(xi, vertices, tolerance=GEOMETRY_TOL):
    """格心是否在多面体 conv(vertices) 内（xi，容差 GEOMETRY_TOL）。"""
    vertices = np.asarray(vertices, float)
    if len(vertices) <= xi.shape[1] or polytope_volume(vertices) <= 0.:
        return np.zeros(len(xi), bool)
    return contains(xi, halfspaces(vertices), tolerance)
