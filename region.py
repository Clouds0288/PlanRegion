"""三维负荷空间中的多面体裁剪、同方案认证及并集覆盖。

顶点的坐标单位均为 kW。halfspaces 返回 a·p+b≤0；SP 割使用 a+b·p≥0，
两种符号在调用处显式转换。一个固定方案的 LP/SOCP 投影是凸集，
不同离散建设方案的并集通常不是凸集，不能跨方案对可行点取凸包。
"""
from heapq import heappop, heappush
from itertools import combinations, count

import numpy as np
from scipy.spatial import ConvexHull

GEOMETRY_TOL = 1e-7  # 单位法向量下的距离容差，单位 kW。


def simplex(network):
    """p≥0、Σp≤power_limit 定义初始外域；power_limit 是有功负荷界（kW）。"""
    return np.vstack([np.zeros(3), network.power_limit*np.eye(3)])  # 原点加三个轴截距，构成四面体。


def polytope_vertices(points):
    """保留凸包顶点；空集、点、线、面也是合法的切割结果。"""
    points = np.asarray(points).reshape(-1, 3)  # 每一行是一个三维负荷点。
    # 舍入只用于去重；保留交点原值，避免扰动共面性、制造额外顶点。
    _, indices = np.unique(np.round(points, 10), axis=0, return_index=True)
    points = points[indices]  # 删除数值重复点，同时保留原始坐标。
    if len(points) < 2:  # 空域和单点不需要构造凸包。
        return points
    centered = points-points[0]  # 平移到过原点的仿射子空间。
    rank = np.linalg.matrix_rank(centered, tol=1e-8)  # 确认当前区域是点、线、面还是三维体。
    if rank == 0:  # 点间距离已在秩容差内，保留一个代表点。
        return points[:1]
    basis = np.linalg.svd(centered, full_matrices=False)[2][:rank]  # 非零奇异值对应的正交基。
    coordinates = centered@basis.T  # 在实际维数中求凸包，避免将共面集当成三维体。
    if rank == 1:
        return points[sorted([coordinates[:, 0].argmin(), coordinates[:, 0].argmax()])]  # 线段只需两个端点。
    return points[np.sort(ConvexHull(coordinates).vertices)]  # 去掉面内点及体内点，保留真正的转折顶点。


def add_certificate(polytopes, row):
    """在线性流程中，仅将同一方案已通过 LP SP 的查询点加入认证凸包。"""
    if row["eta"] is None or row["cut"] is not None:  # MP 无解或 SP 返回分离割，都不是可行点证书。
        return
    key = row["design"]  # 认证集按方案编号隔离，禁止跨方案取凸包。
    previous = polytopes.get(key, np.empty((0, 3)))  # 该方案尚未认证时，从空集开始。
    points = polytope_vertices(np.vstack([previous, row["p"]]))  # 凸性保证同方案可行点的凸组合仍可行。
    if not np.array_equal(points, previous):
        polytopes[key] = points  # 仅几何变化时替换数组，供搜索器识别证书版本。


def covered(polytope, equations):
    """单个认证凸域包含整个候选域的快速充分条件。"""
    return any(contains(polytope, eq).all() for eq in equations)  # 顶点全在某个凸域中，则整个候选凸包都在其中。


def halfspaces(points):
    """凸域的单位法向量半空间；同时支持点、线、面。"""
    center = points[0]  # 用一个已知点确定区域所在的仿射子空间。
    delta = points-center
    rank = np.linalg.matrix_rank(delta, tol=1e-8)  # 裁剪后合法地退化为点、线或面时也必须正确表示。
    basis = np.linalg.svd(delta, full_matrices=True)[2]  # 完整正交基同时保留区域方向及其法向方向。
    coordinates = delta@basis[:rank].T  # 投影到区域的实际维数。
    if rank >= 2:
        hull = ConvexHull(coordinates)  # 在二维或三维中求支持半空间 a·q+b≤0。
        normals = hull.equations[:, :rank]@basis[:rank]  # 将法向量映回原三维负荷坐标。
        equations = np.c_[normals, hull.equations[:, rank]-normals@center]  # 补回平移产生的常数项。
    elif rank == 1:
        normals = np.array([basis[0], -basis[0]])  # 线段的两个端点分别给出一个有向界。
        equations = np.c_[normals, [-coordinates.max(), coordinates.min()]-normals@center]  # q≤qmax，q≥qmin。
    else:
        equations = np.empty((0, 4))  # 单点只需后面的三个坐标方向等式。
    normal = basis[rank:]  # 与区域垂直的方向必须满足 n·(p-center)=0。
    equations = np.vstack([equations, np.c_[normal, -normal@center],
                           np.c_[-normal, normal@center]])  # 用正反两个不等式表示每个等式。
    _, indices = np.unique(np.round(equations, 10), axis=0, return_index=True)  # 合并共面的重复支持面。
    return equations[indices]  # 每行前三项为单位法向量，最后一项为常数。


def contains(points, equations):
    return np.all(points@equations[:, :3].T+equations[:, 3] <= GEOMETRY_TOL, axis=1)  # 同时满足全部半空间才在域内。


def polytope_size(poly):
    """维数、体积（退化时为面积或长度）、最大总负荷；仅用于搜索排序。"""
    delta = poly-poly[0]
    rank = np.linalg.matrix_rank(delta, tol=1e-8)  # 高维残余块优先搜索。
    if rank >= 2:
        coords = delta@np.linalg.svd(delta, full_matrices=False)[2][:rank].T  # 使用内在维数计算几何大小。
        size = ConvexHull(coords).volume  # 二维时返回面积，三维时返回体积。
    else:
        size = np.linalg.norm(delta, axis=1).max()  # 低维块的搜索尺度；不作为实验误差指标。
    return rank, size, poly.sum(axis=1).max()  # 同维优先大块，再优先较大总负荷。


def subtract_polytope(poly, equations):
    """扣除一个认证凸域，返回剩余凸块的闭包，保留交界顶点供同方案认证。"""
    values = poly@equations[:, :3].T+equations[:, 3]  # 各候选顶点相对认证凸域各面的带符号距离。
    if np.all(values <= GEOMETRY_TOL):
        return []  # 整块已在认证凸域内，无需继续查询。
    if np.any(np.min(values, axis=0) > GEOMETRY_TOL):
        return [poly]  # 存在一个面将整块隔在域外，两者不相交。
    pieces = []
    for equation in equations:
        values = poly@equation[:3]+equation[3]  # 按认证域的面逐个分割尚未处理的内部残块。
        if values.max() <= GEOMETRY_TOL:
            continue  # 整块满足这一面，本轮没有可留下的域外部分。
        if values.min() >= -GEOMETRY_TOL:
            pieces.append(poly)  # 整块均在本面外侧或交界处，不必再用其他面切分。
            return pieces
        pieces.append(clip_polytope(poly, equation[3], equation[:3]))  # 保留 a·p+b≥0 的域外块。
        poly = clip_polytope(poly, -equation[3], -equation[:3])  # a·p+b≤0 的剩余块继续与下一面相交。
    return pieces


def uncovered_cells(poly, certificates):
    """U 减去认证并集；低维认证集不能覆盖高维凸块的相对内部。"""
    rank = polytope_size(poly)[0]  # 当前候选块的实际维数。
    equations = [eq for cert_rank, eq in certificates if cert_rank >= rank]  # 线、面不能覆盖三维体的内部。
    if covered(poly, equations):
        return []
    cells = [poly]
    for eq in equations:
        cells = [part for cell in cells for part in subtract_polytope(cell, eq)]  # 依次扣除各方案认证凸域，得到并集之外的残块。
        if not cells:
            break
    return cells


class ResidualSearch:
    """按需扣除认证并集，优先搜索贡献较大的凸块；块内顶点须在同方案下获证。"""

    def __init__(self, polytopes):
        self.heap, self.sequence = [], count()  # 优先队列；递增序号用于同优先级时稳定排序。
        self.cache, self.outer_cache = {}, {}  # 分别缓存认证内域与候选外域的半空间表示。
        self.version, self.active = 0, None  # 认证并集版本，以及当前正在补齐顶点证明的凸块。
        for i, poly in enumerate(polytopes):
            if len(poly):
                self.push(i, poly, None, -1)  # i 就是方案编号；新块尚未与当前外域、认证并集同步。

    def push(self, i, poly, outer, version):
        priority = tuple(-value for value in polytope_size(poly))  # 最小堆使用负值，让较大的几何块先出队。
        heappush(self.heap, (*priority, next(self.sequence), i, poly, outer, version))  # 同时记下该块依据的外域及认证版本。

    def clip(self, i, poly, outer):
        if i not in self.outer_cache or self.outer_cache[i][0] is not outer:  # 外域被新割替换后才重建半空间。
            self.outer_cache[i] = (outer, halfspaces(outer))
        for eq in self.outer_cache[i][1]:
            poly = clip_polytope(poly, -eq[3], -eq[:3])  # 将旧残块与最新方案外域求交。
            if not len(poly):
                break
        return poly

    def next(self, polytopes, certified):
        for key, poly in certified.items():
            if key not in self.cache or self.cache[key][0] is not poly:  # 只更新本轮扩大的认证凸域。
                self.cache[key] = (poly, halfspaces(poly), polytope_size(poly))  # 缓存顶点、半空间和搜索优先级。
                self.version += 1  # 所有旧残块需要按新版认证并集重新判断覆盖。
        ordered = sorted(self.cache.values(), key=lambda value: value[2], reverse=True)  # 优先用较大的认证域覆盖候选块。
        certificates = [(value[2][0], value[1]) for value in ordered]  # 扣除并集时需要维数和半空间。
        while self.active is not None or self.heap:
            if self.active is not None:
                i, poly, outer = self.active  # 固定方案 i 内的候选凸块。
                if polytopes[i] is not outer:  # 新割改变了该方案外域，旧块必须先裁剪再入队。
                    self.active = None
                    if len(polytopes[i]):
                        poly = self.clip(i, poly, polytopes[i])
                        if len(poly):
                            self.push(i, poly, polytopes[i], -1)
                    continue
                own = self.cache.get(i)  # 本方案的内点凸包；不能用其他方案证书作凸组合。
                candidates = poly if own is None else poly[~contains(poly, own[1])]  # 只查询本方案尚未覆盖的顶点。
                if not len(candidates):  # 本块全部顶点已由同一方案认证，整个凸块可退出搜索。
                    self.active = None
                    continue
                # 未覆盖顶点优先；交界顶点即使全局可行，也可用于补齐本块的同方案证明。
                known = np.zeros(len(candidates), dtype=bool)  # 记录顶点是否已在任意方案下可行。
                for _, eq in certificates:
                    known |= contains(candidates, eq)  # 这是全局可行性，只影响查询顺序。
                index = min(range(len(candidates)), key=lambda j: (known[j], -sum(candidates[j])))  # 未知点优先，其次较大总负荷。
                return i, candidates[index]  # 交给该方案自己的 SP 判定。
            _, _, _, _, i, poly, outer, version = heappop(self.heap)  # 取优先级最高的残块。
            if not len(polytopes[i]):
                continue
            if polytopes[i] is not outer:
                poly = self.clip(i, poly, polytopes[i])
                if not len(poly):
                    continue
            if version != self.version or polytopes[i] is not outer:  # 外割或内域更新后，之前的覆盖结论已过期。
                for cell in uncovered_cells(poly, certificates):
                    self.push(i, cell, polytopes[i], self.version)  # 只把尚未被认证并集覆盖的凸块放回队列。
                continue
            self.active = i, poly, polytopes[i]  # 本块已同步，可在下一轮逐顶点补齐同方案证明。
        # 每一块均已被有效割删除或被认证并集覆盖，才能耗尽队列。
        return None


def clip_polytope(vertices, constant, coefficient):
    """constant + coefficient·p ≥ 0；候选顶点由外近似自身产生。"""
    values = constant+vertices@coefficient  # 统一使用“余量≥0”的割方向。
    if np.all(values >= -1e-9):  # 全部顶点在保留侧或容差内，整个凸包均保留。
        return vertices
    if np.all(values < 0):  # 整个凸包都位于被删除的一侧。
        return np.empty((0, 3))
    points = [p for p, value in zip(vertices, values) if value >= 0]  # 原有且仍满足割的顶点。
    for a, b in combinations(range(len(vertices)), 2):  # 检查所有点对，包含全部真实边；额外交点后续去掉。
        if values[a]*values[b] < 0:  # 两端异号时，线段与割平面有唯一交点。
            points.append(vertices[a]+(vertices[b]-vertices[a])*values[a]/(values[a]-values[b]))  # 线性插值解 a+b·p=0。
    return polytope_vertices(points)  # 重新取凸包，只留下裁剪后真实顶点。


def update_outer(costs, polytopes, row):
    """方案割只更新对应外域；已认证 MP2 的总量上界适用于其预算内方案。"""
    updated = list(polytopes)  # 只复制容器；未变化方案保留数组身份供搜索缓存识别。
    if row['cut'] is not None:
        i, cut = row['design'], row['cut']
        updated[i] = clip_polytope(updated[i], cut[0], cut[1:])  # 固定方案割只能裁剪方案 i。
    elif row['mode'] == 'MP2':  # 无分离割时，MP2 的最优目标已由 SP 认证。
        for i, cost in enumerate(costs):
            if cost <= row['target']:  # 该最优总量是预算内全部方案的共同上界。
                updated[i] = (clip_polytope(updated[i], sum(row['p']), -np.ones(3))
                              if row['p'] is not None else np.empty((0, 3)))  # MP2 不可行则对应预算内候选域全部为空。
    return updated
