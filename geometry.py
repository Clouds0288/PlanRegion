"""凸多面体与扫描格的几何：凸包顶点、半空间、裁剪、测度与并集，以及扫描格心。

同一坐标系内运算（归一化 xi 或真实 kW），不做单位换算，不依赖求解器；构域（region）、回放（monitor）与扫描比较
（vertify）共用。面方程为 F@xi+g<=0，裁剪接口为 constant+coefficient@xi>=0。
"""
from itertools import combinations, product

import numpy as np
from scipy.spatial import ConvexHull
from shapely import set_precision
from shapely.geometry import MultiPoint, Polygon
from shapely.ops import unary_union

# 连续几何在公共评价箱归一化后的坐标中计算；与采样网格无关。
GEOMETRY_TOL = 1e-8


def polytope_vertices(points):
    """任意仿射维数的凸包顶点，保持原始交点坐标；points 为非空 (n,d)。"""
    points = np.asarray(points, dtype=float)
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
    """非空凸多面体（顶点 (n,d)）的 [F,g]，内侧 F@xi+g<=0；低维集以成对不等式表示仿射等式。"""
    points = np.asarray(points)
    d = points.shape[1]
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
    """用 constant + coefficient@xi >= 0 裁剪；与输入顶点使用同一坐标系，空集保持为空。"""
    d = len(coefficient)
    vertices = np.asarray(vertices).reshape(-1, d)
    if not len(vertices):
        return vertices
    coefficient = np.asarray(coefficient)
    norm = np.linalg.norm(coefficient)
    if norm < 1e-20:   # 系数为零（如只含 x 的联合割）：常数项决定保留或清空
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
    """(n,d) 顶点凸包的 d 维测度；退化或空集为 0。"""
    poly = np.asarray(poly)
    d = poly.shape[1]
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


def box_vertices(d):
    """单位盒 [0,1]^d 的顶点。"""
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


def cell_centers(grid):
    """扫描格（axis_lower、bounds、states 的形状）的格心，真实 kW，按 states 的 C 序展平为 (N,d)。"""
    shape = np.shape(grid['states'])
    lower = np.asarray(grid['axis_lower'], float)
    indices = np.indices(shape).reshape(len(shape), -1).T
    return lower+(indices+.5)*(np.asarray(grid['bounds'], float)-lower)/np.array(shape)
