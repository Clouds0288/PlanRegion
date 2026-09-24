"""两阶段构域的纯几何运算；二维/三维，所有负荷坐标为 kW。"""
from dataclasses import dataclass
from itertools import combinations, product
from types import SimpleNamespace

import numpy as np
from scipy.spatial import ConvexHull

# 几何舍入容差 (kW)，不是构域的 epsilon_kw 或物理模型的 PLANNING_TOL。
GEOMETRY_TOL = 1e-8


def _convex_hull(points):
    """在仿射坐标中先作 SVD 缩放，避免狭长小节点的数值病态。"""
    center = points[0]
    _, scales, basis = np.linalg.svd(points-center, full_matrices=False)
    hull = ConvexHull(((points-center)@basis.T)/scales)
    normal = (hull.equations[:, :-1]/scales)@basis
    lengths = np.linalg.norm(normal, axis=1)
    equations = np.c_[normal/lengths[:, None],
                       (hull.equations[:, -1]-normal@center)/lengths]
    return SimpleNamespace(vertices=hull.vertices, simplices=hull.simplices,
                           equations=equations, volume=float(hull.volume*np.prod(scales)))


def polytope_vertices(points):
    """任意仿射维数的凸包顶点，保持原始交点坐标。"""
    points = np.asarray(points, dtype=float)
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
    return points[np.sort(_convex_hull(coordinates).vertices)]


def halfspaces(points):
    """返回 [F,g]，内侧 F@p+g<=0；低维集以成对不等式表示仿射等式。"""
    points = np.asarray(points)
    if not len(points):
        return np.array([[*np.zeros(points.shape[1]), 1.]])
    center, delta = points[0], points-points[0]
    rank = np.linalg.matrix_rank(delta, tol=1e-10)
    basis = np.linalg.svd(delta, full_matrices=True)[2]
    coordinates = delta@basis[:rank].T
    if rank >= 2:
        hull = _convex_hull(coordinates)
        normal = hull.equations[:, :rank]@basis[:rank]
        eq = np.c_[normal, hull.equations[:, rank]-normal@center]
    elif rank == 1:
        normal = np.array([basis[0], -basis[0]])
        eq = np.c_[normal, [-coordinates.max(), coordinates.min()]-normal@center]
    else:
        eq = np.empty((0, points.shape[1]+1))
    normal = basis[rank:]
    eq = np.vstack([eq, np.c_[normal, -normal@center], np.c_[-normal, normal@center]])
    _, indices = np.unique(np.round(eq, 10), axis=0, return_index=True)
    return eq[np.sort(indices)]


def contains(points, equations, tolerance=GEOMETRY_TOL):
    return np.all(np.atleast_2d(points)@equations[:, :-1].T+equations[:, -1] <= tolerance, axis=1)


def clip_polytope(vertices, constant, coefficient):
    """用 constant + coefficient@p >= 0 裁剪；与输入顶点使用同一坐标系。"""
    vertices = np.asarray(vertices)
    dimension = len(coefficient)
    if not len(vertices):
        return vertices
    coefficient = np.asarray(coefficient)
    norm = np.linalg.norm(coefficient)
    if norm < 1e-20:
        return vertices if constant >= -1e-12 else np.empty((0, dimension))
    values = (constant+vertices@coefficient)/norm
    if values.min() >= -1e-11:
        return vertices
    if values.max() < -1e-11:
        return np.empty((0, dimension))
    points = list(vertices[values >= -1e-11])
    if len(vertices) > dimension and np.linalg.matrix_rank(vertices-vertices[0], tol=1e-10) == dimension:
        edges = {tuple(sorted(edge)) for face in _convex_hull(vertices).simplices
                 for edge in combinations(face, 2)}
    else:
        edges = combinations(range(len(vertices)), 2)
    for a, b in edges:
        if values[a]*values[b] < 0:
            points.append(vertices[a]+(vertices[b]-vertices[a])*values[a]/(values[a]-values[b]))
    return polytope_vertices(points)


def polytope_volume(poly):
    poly = np.asarray(poly)
    if len(poly) == 0 or len(poly) <= poly.shape[1] or np.linalg.matrix_rank(poly-poly[0], tol=1e-10) < poly.shape[1]:
        return 0.
    return float(_convex_hull(poly).volume)


def downward_facets(points):
    """Conservative nonnegative-normal facets of one downward convex hull."""
    points = np.maximum(np.asarray(points)-1e-5, 0.)
    dimension = points.shape[1]
    masks = np.asarray(list(product((0., 1.), repeat=dimension)))
    cloud = np.unique((points[:, None, :]*masks).reshape(-1, dimension), axis=0)
    active = np.max(cloud, axis=0) > 1e-8
    rows = []
    if active.sum() >= 2:
        equations = ConvexHull(cloud[:, active]).equations
        for equation in equations:
            if np.min(equation[:-1]) < -1e-9 or equation[:-1].sum() <= 1e-10:
                continue
            weights = np.zeros(dimension)
            weights[active] = np.maximum(equation[:-1], 0.)
            weights /= weights.sum()
            rows.append(np.r_[weights, -max(0., np.max(cloud @ weights)-1e-6)])
    for axis in np.where(~active)[0] if active.sum() >= 2 else range(dimension):
        weights = np.eye(dimension)[axis]
        rows.append(np.r_[weights, -cloud[:, axis].max()])
    rows = np.asarray(rows).reshape(-1, dimension+1)
    _, indices = np.unique(np.round(rows, 9), axis=0, return_index=True)
    return rows[np.sort(indices)]


def scheme_hulls(certificates):
    groups = {}
    for certificate in certificates:
        groups.setdefault(tuple(certificate['x']), []).append(certificate['p'])
    return [(key, downward_facets(points)) for key, points in groups.items()]


def hull_distance(points, facets):
    """Exact L-infinity distance to a downward convex set, including axis faces."""
    points = np.atleast_2d(points)
    weights, bound = facets[:, :-1], -facets[:, -1]
    remaining = np.broadcast_to(weights, (len(points), *weights.shape)).copy()
    distance = np.zeros(remaining.shape[:2])
    for _ in range(points.shape[1]):
        total = remaining.sum(axis=2)
        numerator = np.einsum('pfi,pi->pf', remaining, points)-bound
        distance = np.maximum(0., np.divide(numerator, total, out=distance.copy(), where=total > 0.))
        remaining *= points[:, None, :] > distance[:, :, None]
    return distance.max(axis=1)


@dataclass
class Cell:
    vertices: np.ndarray
    depth: int = 0
    gap: float = np.inf
    scheme: tuple = ()


def cell_gap(vertices, hulls):
    """Δ_s = min_x max_{v∈vertices(O_s)} dist∞(v,I_x)，界住节点内所有点。

    不能交换 min/max：不同顶点属于不同方案时，中间仍可能不可行。
    """
    distances = np.array([hull_distance(vertices, facets) for _, facets in hulls])
    worst = distances.max(axis=1)
    best = int(worst.argmin())
    point_distances = distances.min(axis=0)
    return float(worst[best]+1e-5), hulls[best][0], point_distances


def select_cell(cells, hulls):
    """最大认证间隙优先；内域扩张只会降低旧界，因此可惰性重算最大节点。"""
    while True:
        chosen = int(np.argmax([cell.gap for cell in cells]))
        cell = cells[chosen]
        cell.gap, cell.scheme, point_distances = cell_gap(cell.vertices, hulls)
        if cell.gap >= max(other.gap for other in cells)-1e-9:
            return chosen, point_distances


def subtract_disjunction(vertices, threshold, cut_normals):
    """保留 OR_k W_k p≤threshold_k；逐项分区，节点内部不交，边界保留。"""
    remaining, pieces = np.asarray(vertices), []
    for normal, bound in zip(cut_normals, threshold):
        retained = clip_polytope(remaining, bound, -normal)
        if len(retained):
            pieces.append(retained)
        remaining = clip_polytope(remaining, -bound, normal)
        if not len(remaining):
            break
    return pieces


def apply_disjunction(cells, threshold, cut_normals):
    result = []
    for cell in cells:
        if np.any((cell.vertices @ cut_normals.T).max(axis=0) <= threshold+1e-9):
            result.append(cell)
        else:
            result.extend(Cell(v, cell.depth+1, cell.gap, cell.scheme)
                          for v in subtract_disjunction(cell.vertices, threshold, cut_normals))
    return result


def disjunction_directions(target, hulls):
    candidates = []
    for _, rows in hulls:
        residuals = rows[:, :-1] @ target+rows[:, -1]
        face = int(residuals.argmax())
        candidates.append((residuals[face], rows[face, :-1]))
    directions = []
    for _, normal in sorted(candidates, key=lambda pair: pair[0]):
        if all(np.linalg.norm(normal-old) > .04 for old in directions):
            directions.append(normal)
        if len(directions) == len(target):
            break
    return np.asarray(directions)


def cover_partition(vertices, facets, epsilon_kw):
    """按 I_x^epsilon 分支，不删除点；这是覆盖分支，不是不可行性割。"""
    dimension = vertices.shape[1]
    masks = np.asarray(list(product((0., 1.), repeat=dimension)))
    normals = (facets[:,None,:-1]*masks[None,1:,:]).reshape(-1,dimension)
    scale = normals.sum(axis=1)
    bounds = np.repeat(-facets[:,-1],len(masks)-1)+epsilon_kw*scale
    active = scale > 1e-12
    normals,bounds = normals[active]/scale[active,None],bounds[active]/scale[active]
    residuals = vertices @ normals.T-bounds
    if np.any(residuals.min(axis=0) > 1e-9):
        return None
    indices = np.where(residuals.max(axis=0) > 1e-9)[0]
    indices = indices[np.argsort(-residuals.max(axis=0)[indices])]
    remaining,pieces = vertices,[]
    for index in indices:
        normal,bound = normals[index],bounds[index]
        if np.max(remaining @ normal-bound) <= 1e-9:
            continue
        outside = clip_polytope(remaining,-bound,normal)
        remaining = clip_polytope(remaining,bound,-normal)
        if len(outside):
            pieces.append(outside)
        if not len(remaining):
            return None
    # A tangent face gives the unchanged outside node again; it is not progress.
    if np.linalg.matrix_rank(remaining-remaining[0], tol=1e-7) < np.linalg.matrix_rank(vertices-vertices[0], tol=1e-7):
        return None
    return pieces+[remaining] if pieces else None


def support_direction(vertices, certificates, excluded=()):
    facets = downward_facets([c['p'] for c in certificates])
    if len(excluded):
        facets = facets[[all(np.max(np.abs(row[:-1]-old)) > .02 for old in excluded) for row in facets]]
    if not len(facets):
        return 0., None
    gaps = vertices @ facets[:, :-1].T+facets[:, -1]
    index = int(np.max(gaps, axis=0).argmax())
    return float(gaps[:, index].max()), facets[index, :-1]
