"""规划域的候选多面体、可行性割、同方案认证与覆盖判断。"""
from heapq import heappop, heappush
from itertools import combinations, count

import numpy as np
from scipy.spatial import ConvexHull

GEOMETRY_TOL = 1e-7  # 单位法向量下的距离容差，单位 kW。


def simplex(network):
    """p ≥ 0、Σp ≤ 配变容量定义三个自由负荷节点的初始候选域。"""
    return np.vstack([np.zeros(3), network.power_limit*np.eye(3)])


def polytope_vertices(points):
    """保留凸包顶点；空集、点、线、面也是合法的切割结果。"""
    points = np.asarray(points).reshape(-1, 3)
    # 舍入只用于去重；保留交点原值，避免扰动共面性、制造额外顶点。
    _, indices = np.unique(np.round(points, 10), axis=0, return_index=True)
    points = points[indices]
    if len(points) < 2:
        return points
    centered = points-points[0]
    rank = np.linalg.matrix_rank(centered, tol=1e-8)
    if rank == 0:
        return points[:1]
    basis = np.linalg.svd(centered, full_matrices=False)[2][:rank]
    coordinates = centered@basis.T
    if rank == 1:
        return points[sorted([coordinates[:, 0].argmin(), coordinates[:, 0].argmax()])]
    return points[np.sort(ConvexHull(coordinates).vertices)]


def solid(points):
    return len(points) >= 4 and np.linalg.matrix_rank(points-points[0], tol=1e-8) == 3


def add_certificate(polytopes, row):
    """只对同一建设方案的 SP 可行点取凸包；history 是唯一原始记录。"""
    if row["eta"] is None or row["cut"] is not None:
        return False
    key = tuple(row["x"])
    previous = polytopes.get(key, np.empty((0, 3)))
    points = polytope_vertices(np.vstack([previous, row["p"]]))
    if np.array_equal(points, previous):
        return False
    polytopes[key] = points
    return True


def covered(polytope, equations):
    """单个认证凸域包含整个候选域的快速充分条件。"""
    return any(contains(polytope, eq).all() for eq in equations)


def halfspaces(points):
    """凸域的单位法向量半空间；同时支持点、线、面。"""
    center = points[0]
    delta = points-center
    rank = np.linalg.matrix_rank(delta, tol=1e-8)
    basis = np.linalg.svd(delta, full_matrices=True)[2]
    coordinates = delta@basis[:rank].T
    if rank >= 2:
        hull = ConvexHull(coordinates)
        normals = hull.equations[:, :rank]@basis[:rank]
        equations = np.c_[normals, hull.equations[:, rank]-normals@center]
    elif rank == 1:
        normals = np.array([basis[0], -basis[0]])
        equations = np.c_[normals, [-coordinates.max(), coordinates.min()]-normals@center]
    else:
        equations = np.empty((0, 4))
    normal = basis[rank:]
    equations = np.vstack([equations, np.c_[normal, -normal@center],
                           np.c_[-normal, normal@center]])
    _, indices = np.unique(np.round(equations, 10), axis=0, return_index=True)
    return equations[indices]


def contains(points, equations):
    return np.all(points@equations[:, :3].T+equations[:, 3] <= GEOMETRY_TOL, axis=1)


def polytope_size(poly):
    """维数、体积（退化时为面积或长度）、最大总负荷；仅用于搜索排序。"""
    delta = poly-poly[0]
    rank = np.linalg.matrix_rank(delta, tol=1e-8)
    if rank >= 2:
        coords = delta@np.linalg.svd(delta, full_matrices=False)[2][:rank].T
        size = ConvexHull(coords).volume
    else:
        size = np.linalg.norm(delta, axis=1).max()
    return rank, size, poly.sum(axis=1).max()


def subtract_polytope(poly, equations):
    """扣除一个认证凸域，返回剩余凸块的闭包，保留交界顶点供同方案认证。"""
    values = poly@equations[:, :3].T+equations[:, 3]
    if np.all(values <= GEOMETRY_TOL):
        return []
    if np.any(np.min(values, axis=0) > GEOMETRY_TOL):
        return [poly]
    pieces = []
    for equation in equations:
        values = poly@equation[:3]+equation[3]
        if values.max() <= GEOMETRY_TOL:
            continue
        if values.min() >= -GEOMETRY_TOL:
            pieces.append(poly)
            return pieces
        pieces.append(clip_polytope(poly, equation[3], equation[:3]))
        poly = clip_polytope(poly, -equation[3], -equation[:3])
    return pieces


def uncovered_cells(poly, certificates):
    """U 减去认证并集；低维认证集不能覆盖高维凸块的相对内部。"""
    rank = polytope_size(poly)[0]
    equations = [eq for cert_rank, eq in certificates if cert_rank >= rank]
    if covered(poly, equations):
        return []
    cells = [poly]
    for eq in equations:
        cells = [part for cell in cells for part in subtract_polytope(cell, eq)]
        if not cells:
            break
    return cells


class ResidualSearch:
    """按需扣除认证并集，优先搜索贡献较大的凸块；块内顶点须在同方案下获证。"""

    def __init__(self, designs, polytopes):
        self.designs = designs
        self.heap, self.sequence = [], count()
        self.cache, self.outer_cache = {}, {}
        self.version, self.active = 0, None
        for i, poly in enumerate(polytopes):
            if len(poly):
                self.push(i, poly, None, -1)

    def push(self, i, poly, outer, version):
        priority = tuple(-value for value in polytope_size(poly))
        heappush(self.heap, (*priority, next(self.sequence), i, poly, outer, version))

    def clip(self, i, poly, outer):
        if i not in self.outer_cache or self.outer_cache[i][0] is not outer:
            self.outer_cache[i] = (outer, halfspaces(outer))
        for eq in self.outer_cache[i][1]:
            poly = clip_polytope(poly, -eq[3], -eq[:3])
            if not len(poly):
                break
        return poly

    def next(self, polytopes, certified):
        for key, poly in certified.items():
            if key not in self.cache or self.cache[key][0] is not poly:
                self.cache[key] = (poly, halfspaces(poly), polytope_size(poly))
                self.version += 1
        ordered = sorted(self.cache.values(), key=lambda value: value[2], reverse=True)
        certificates = [(value[2][0], value[1]) for value in ordered]
        while self.active is not None or self.heap:
            if self.active is not None:
                i, poly, outer = self.active
                if polytopes[i] is not outer:
                    self.active = None
                    if len(polytopes[i]):
                        poly = self.clip(i, poly, polytopes[i])
                        if len(poly):
                            self.push(i, poly, polytopes[i], -1)
                    continue
                own = self.cache.get(tuple(self.designs[i]))
                candidates = poly if own is None else poly[~contains(poly, own[1])]
                if not len(candidates):
                    self.active = None
                    continue
                # 未覆盖顶点优先；交界顶点即使全局可行，也可用于补齐本块的同方案证明。
                known = np.zeros(len(candidates), dtype=bool)
                for _, eq in certificates:
                    known |= contains(candidates, eq)
                index = min(range(len(candidates)), key=lambda j: (known[j], -sum(candidates[j])))
                return self.designs[i], candidates[index]
            _, _, _, _, i, poly, outer, version = heappop(self.heap)
            if not len(polytopes[i]):
                continue
            if polytopes[i] is not outer:
                poly = self.clip(i, poly, polytopes[i])
                if not len(poly):
                    continue
            if version != self.version or polytopes[i] is not outer:
                for cell in uncovered_cells(poly, certificates):
                    self.push(i, cell, polytopes[i], self.version)
                continue
            self.active = i, poly, polytopes[i]
        # 每一块均已被有效割删除或被认证并集覆盖，才能耗尽队列。
        return None


def clip_polytope(vertices, constant, coefficient):
    """constant + coefficient·p ≥ 0；候选顶点由外近似自身产生。"""
    values = constant+vertices@coefficient
    if np.all(values >= -1e-9):
        return vertices
    if np.all(values < 0):
        return np.empty((0, 3))
    points = [p for p, value in zip(vertices, values) if value >= 0]
    for a, b in combinations(range(len(vertices)), 2):
        if values[a]*values[b] < 0:
            points.append(vertices[a]+(vertices[b]-vertices[a])*values[a]/(values[a]-values[b]))
    return polytope_vertices(points)


def cut_polytopes(designs, polytopes, cut):
    nx = designs.shape[1]
    return [clip_polytope(poly, cut[0]+cut[1:1+nx]@x, cut[1+nx:])
            for x, poly in zip(designs, polytopes)]


def update_outer(costs, designs, polytopes, row, query):
    """SP 割及 MP2 已认证总量上界均从原始历史重建，不另存一份事件记录。"""
    if row["cut"] is not None:
        return cut_polytopes(designs, polytopes, row["cut"])
    if query["mode"] == "MP2":
        budget = np.inf if query["target"] is None else query["target"]
        return [poly if cost > budget else
                clip_polytope(poly, sum(row["p"]), -np.ones(3)) if row["p"] is not None else
                np.empty((0, 3)) for cost, poly in zip(costs, polytopes)]
    return polytopes
