"""三维负荷空间中的多面体裁剪、同方案认证及并集覆盖。

顶点的坐标单位均为 kW。halfspaces 返回 a·p+b≤0；SP 割使用 a+b·p≥0，
两种符号在调用处显式转换。一个固定方案的 LP/SOCP 投影是凸集，
不同离散建设方案的并集通常不是凸集，不能跨方案对可行点取凸包。
"""
from heapq import heappop, heappush  # 用优先队列按剩余区域大小安排查询。
from itertools import combinations, count  # 分别枚举顶点对并生成稳定的队列序号。

import numpy as np  # 处理三维坐标、半空间和线性代数。
from scipy.spatial import ConvexHull  # 计算凸包顶点和支持半空间。

GEOMETRY_TOL = 1e-7  # 单位法向量下的距离容差，单位 kW。


def simplex(network):  # 构造由总负荷上界限定的初始三维候选域。
    """p≥0、Σp≤power_limit 定义初始外域；power_limit 是有功负荷界（kW）。"""
    return np.vstack([np.zeros(3), network.power_limit*np.eye(3)])  # 原点加三个轴截距，构成四面体。


def polytope_vertices(points):  # 从点集中提取实际维数下的凸包顶点。
    """保留凸包顶点；空集、点、线、面也是合法的切割结果。"""
    points = np.asarray(points).reshape(-1, 3)  # 每一行是一个三维负荷点。
    # 舍入只用于去重；保留交点原值，避免扰动共面性、制造额外顶点。
    _, indices = np.unique(np.round(points, 10), axis=0, return_index=True)  # 仅按舍入坐标识别重复项，索引指向未舍入的原值。
    points = points[indices]  # 删除数值重复点，同时保留原始坐标。
    if len(points) < 2:  # 空域和单点不需要构造凸包。
        return points  # 空集或单点本身就是顶点表示。
    centered = points-points[0]  # 平移到过原点的仿射子空间。
    rank = np.linalg.matrix_rank(centered, tol=1e-8)  # 确认当前区域是点、线、面还是三维体。
    if rank == 0:  # 点间距离已在秩容差内，保留一个代表点。
        return points[:1]  # 在几何容差内重合的点只保留第一个。
    basis = np.linalg.svd(centered, full_matrices=False)[2][:rank]  # 非零奇异值对应的正交基。
    coordinates = centered@basis.T  # 在实际维数中求凸包，避免将共面集当成三维体。
    if rank == 1:  # 一维凸包退化为线段。
        return points[sorted([coordinates[:, 0].argmin(), coordinates[:, 0].argmax()])]  # 线段只需两个端点。
    return points[np.sort(ConvexHull(coordinates).vertices)]  # 去掉面内点及体内点，保留真正的转折顶点。


def add_certificate(polytopes, row):  # 将已认证查询纳入其建设方案的内域。
    """在线性流程中，仅将同一方案已通过 LP SP 的查询点加入认证凸包。"""
    if not row['feasible']:  # 只有原约束核验通过才有证书；无分离割不等于可行。
        return  # 未获证查询不能改变认证集合。
    key = row["design"]  # 认证集按方案编号隔离，禁止跨方案取凸包。
    previous = polytopes.get(key, np.empty((0, 3)))  # 该方案尚未认证时，从空集开始。
    points = polytope_vertices(np.vstack([previous, row["p"]]))  # 凸性保证同方案可行点的凸组合仍可行。
    if not np.array_equal(points, previous):  # 数组确有变化时才更新几何版本。
        polytopes[key] = points  # 仅几何变化时替换数组，供搜索器识别证书版本。


def covered(polytope, equations):  # 检查一个候选凸块是否被某一个认证凸域完整包含。
    """单个认证凸域包含整个候选域的快速充分条件。"""
    return any(contains(polytope, eq).all() for eq in equations)  # 顶点全在某个凸域中，则整个候选凸包都在其中。


def halfspaces(points):  # 将任意维数的三维点集转为闭凸域的半空间表示。
    """凸域的单位法向量半空间；同时支持点、线、面。"""
    center = points[0]  # 用一个已知点确定区域所在的仿射子空间。
    delta = points-center  # 平移坐标以确定区域的线性子空间。
    rank = np.linalg.matrix_rank(delta, tol=1e-8)  # 裁剪后合法地退化为点、线或面时也必须正确表示。
    basis = np.linalg.svd(delta, full_matrices=True)[2]  # 完整正交基同时保留区域方向及其法向方向。
    coordinates = delta@basis[:rank].T  # 投影到区域的实际维数。
    if rank >= 2:  # 二维和三维区域在自身子空间内调用凸包算法。
        hull = ConvexHull(coordinates)  # 在二维或三维中求支持半空间 a·q+b≤0。
        normals = hull.equations[:, :rank]@basis[:rank]  # 将法向量映回原三维负荷坐标。
        equations = np.c_[normals, hull.equations[:, rank]-normals@center]  # 补回平移产生的常数项。
    elif rank == 1:  # 一维区域由两个端点界定。
        normals = np.array([basis[0], -basis[0]])  # 线段的两个端点分别给出一个有向界。
        equations = np.c_[normals, [-coordinates.max(), coordinates.min()]-normals@center]  # q≤qmax，q≥qmin。
    else:  # 零维区域随后使用法向等式表示。
        equations = np.empty((0, 4))  # 单点只需后面的三个坐标方向等式。
    normal = basis[rank:]  # 与区域垂直的方向必须满足 n·(p-center)=0。
    equations = np.vstack([equations, np.c_[normal, -normal@center],  # 加入垂直子空间的正向约束。
                           np.c_[-normal, normal@center]])  # 用正反两个不等式表示每个等式。
    _, indices = np.unique(np.round(equations, 10), axis=0, return_index=True)  # 合并共面的重复支持面。
    return equations[indices]  # 每行前三项为单位法向量，最后一项为常数。


def contains(points, equations):  # 逐点判断是否同时满足所有半空间。
    return np.all(points@equations[:, :3].T+equations[:, 3] <= GEOMETRY_TOL, axis=1)  # 同时满足全部半空间才在域内。


def polytope_size(poly):  # 计算残块搜索优先级所需的几何尺度。
    """维数、体积（退化时为面积或长度）、最大总负荷；仅用于搜索排序。"""
    delta = poly-poly[0]  # 平移点集以计算其实际维数。
    rank = np.linalg.matrix_rank(delta, tol=1e-8)  # 高维残余块优先搜索。
    if rank >= 2:  # 面或体使用其内在维度的凸包大小。
        coords = delta@np.linalg.svd(delta, full_matrices=False)[2][:rank].T  # 使用内在维数计算几何大小。
        size = ConvexHull(coords).volume  # 二维时返回面积，三维时返回体积。
    else:  # 线段或单点使用长度作为排序尺度。
        size = np.linalg.norm(delta, axis=1).max()  # 低维块的搜索尺度；不作为实验误差指标。
    return rank, size, poly.sum(axis=1).max()  # 同维优先大块，再优先较大总负荷。


def subtract_polytope(poly, equations):  # 逐面扣除一个凸认证域，返回域外残块。
    """扣除一个认证凸域，返回剩余凸块的闭包，保留交界顶点供同方案认证。"""
    values = poly@equations[:, :3].T+equations[:, 3]  # 各候选顶点相对认证凸域各面的带符号距离。
    if np.all(values <= GEOMETRY_TOL):  # 所有顶点均在认证域内时，由凸性可知整个块已被覆盖。
        return []  # 整块已在认证凸域内，无需继续查询。
    if np.any(np.min(values, axis=0) > GEOMETRY_TOL):  # 某个支持面将全部顶点隔开时，两域没有内部交集。
        return [poly]  # 存在一个面将整块隔在域外，两者不相交。
    pieces = []  # 依次保存各面分割得到的域外部分。
    for equation in equations:  # 尚在认证域各面内侧的残块继续向下处理。
        values = poly@equation[:3]+equation[3]  # 按认证域的面逐个分割尚未处理的内部残块。
        if values.max() <= GEOMETRY_TOL:  # 当前残块完全在这一面的内侧。
            continue  # 整块满足这一面，本轮没有可留下的域外部分。
        if values.min() >= -GEOMETRY_TOL:  # 当前残块已完全落在这一面的外侧或边界。
            pieces.append(poly)  # 整块均在本面外侧或交界处，不必再用其他面切分。
            return pieces  # 整块已分离，无需再检查后续面。
        pieces.append(clip_polytope(poly, equation[3], equation[:3]))  # 保留 a·p+b≥0 的域外块。
        poly = clip_polytope(poly, -equation[3], -equation[:3])  # a·p+b≤0 的剩余块继续与下一面相交。
    return pieces  # 返回所有尚待认证的凸残块。


def uncovered_cells(poly, certificates):  # 计算候选凸块减去各建设方案认证域并集的结果。
    """U 减去认证并集；低维认证集不能覆盖高维凸块的相对内部。"""
    rank = polytope_size(poly)[0]  # 当前候选块的实际维数。
    equations = [eq for cert_rank, eq in certificates if cert_rank >= rank]  # 线、面不能覆盖三维体的内部。
    if covered(poly, equations):  # 先尝试单个凸域完整覆盖，避免不必要的多次分割。
        return []  # 被完整覆盖后没有剩余搜索区域。
    cells = [poly]  # 未被直接覆盖时从整个候选块开始扣除。
    for eq in equations:  # 逐个扣除足够维数的认证域。
        cells = [part for cell in cells for part in subtract_polytope(cell, eq)]  # 依次扣除各方案认证凸域，得到并集之外的残块。
        if not cells:  # 所有残块均被覆盖时停止扣除。
            break  # 无需继续处理后面的认证域。
    return cells  # 返回并集未覆盖部分的凸分块。


class ResidualSearch:  # 维护外域残块队列及按方案隔离的认证缓存。
    """按需扣除认证并集，优先搜索贡献较大的凸块；块内顶点须在同方案下获证。"""

    def __init__(self, polytopes):  # 由各方案初始外域创建待认证搜索队列。
        self.heap, self.sequence = [], count()  # 优先队列；递增序号用于同优先级时稳定排序。
        self.cache, self.outer_cache = {}, {}  # 分别缓存认证内域与候选外域的半空间表示。
        self.version, self.active = 0, None  # 认证并集版本，以及当前正在补齐顶点证明的凸块。
        for i, poly in enumerate(polytopes):  # 各个固定方案对应一个初始凸候选域。
            if len(poly):  # 空方案域没有可搜索的点。
                self.push(i, poly, None, -1)  # i 就是方案编号；新块尚未与当前外域、认证并集同步。

    def push(self, i, poly, outer, version):  # 记录凸块及其几何版本，并按优先级入队。
        priority = tuple(-value for value in polytope_size(poly))  # 最小堆使用负值，让较大的几何块先出队。
        heappush(self.heap, (*priority, next(self.sequence), i, poly, outer, version))  # 同时记下该块依据的外域及认证版本。

    def clip(self, i, poly, outer):  # 把旧搜索块裁剪到该方案的最新候选域中。
        if i not in self.outer_cache or self.outer_cache[i][0] is not outer:  # 外域被新割替换后才重建半空间。
            self.outer_cache[i] = (outer, halfspaces(outer))  # 缓存外域对象及其半空间，避免重复构造凸包。
        for eq in self.outer_cache[i][1]:  # 逐面把旧残块与新外域求交。
            poly = clip_polytope(poly, -eq[3], -eq[:3])  # 将旧残块与最新方案外域求交。
            if not len(poly):  # 交集为空后不再继续裁剪。
                break  # 结束当前残块的逐面求交。
        return poly  # 返回仍可能包含可行点的残块。

    def next(self, polytopes, certified):  # 选择一个需要同方案 SP 认证的候选顶点。
        for key, poly in certified.items():  # 检查哪些方案的认证域发生了扩张。
            if key not in self.cache or self.cache[key][0] is not poly:  # 只更新本轮扩大的认证凸域。
                self.cache[key] = (poly, halfspaces(poly), polytope_size(poly))  # 缓存顶点、半空间和搜索优先级。
                self.version += 1  # 所有旧残块需要按新版认证并集重新判断覆盖。
        ordered = sorted(self.cache.values(), key=lambda value: value[2], reverse=True)  # 优先用较大的认证域覆盖候选块。
        certificates = [(value[2][0], value[1]) for value in ordered]  # 扣除并集时需要维数和半空间。
        while self.active is not None or self.heap:  # 当前块尚未完成或队列未空时继续搜索。
            if self.active is not None:  # 优先补齐正在处理的凸块的顶点证书。
                i, poly, outer = self.active  # 固定方案 i 内的候选凸块。
                if polytopes[i] is not outer:  # 新割改变了该方案外域，旧块必须先裁剪再入队。
                    self.active = None  # 旧活动块失效，重新与外域求交后再排序。
                    if len(polytopes[i]):  # 仅非空的新外域值得继续处理。
                        poly = self.clip(i, poly, polytopes[i])  # 剔除被新增有效割排除的部分。
                        if len(poly):  # 只把仍非空的残块放回队列。
                            self.push(i, poly, polytopes[i], -1)  # 标为待按最新认证并集重新核验的版本。
                    continue  # 重新进入选择流程，不使用已失效的候选点。
                own = self.cache.get(i)  # 本方案的内点凸包；不能用其他方案证书作凸组合。
                candidates = poly if own is None else poly[~contains(poly, own[1])]  # 只查询本方案尚未覆盖的顶点。
                if not len(candidates):  # 本块全部顶点已由同一方案认证，整个凸块可退出搜索。
                    self.active = None  # 本方案已经认证当前凸块的全部顶点。
                    continue  # 结束此块后继续取下一个残块。
                # 未覆盖顶点优先；交界顶点即使全局可行，也可用于补齐本块的同方案证明。
                known = np.zeros(len(candidates), dtype=bool)  # 记录顶点是否已在任意方案下可行。
                for _, eq in certificates:  # 使用已有认证并集标记全局已知可行点。
                    known |= contains(candidates, eq)  # 这是全局可行性，只影响查询顺序。
                index = min(range(len(candidates)), key=lambda j: (known[j], -sum(candidates[j])))  # 未知点优先，其次较大总负荷。
                return i, candidates[index]  # 交给该方案自己的 SP 判定。
            _, _, _, _, i, poly, outer, version = heappop(self.heap)  # 取优先级最高的残块。
            if not len(polytopes[i]):  # 该方案候选域已被割清空。
                continue  # 跳过已不可能贡献可行区域的方案。
            if polytopes[i] is not outer:  # 残块依据的外域对象不是当前版本。
                poly = self.clip(i, poly, polytopes[i])  # 用最新方案外域删除已失效部分。
                if not len(poly):  # 裁剪后为空的块不再参与搜索。
                    continue  # 继续处理队列中的下一块。
            if version != self.version or polytopes[i] is not outer:  # 外割或内域更新后，之前的覆盖结论已过期。
                for cell in uncovered_cells(poly, certificates):  # 按当前全部认证域扣除已覆盖的部分。
                    self.push(i, cell, polytopes[i], self.version)  # 只把尚未被认证并集覆盖的凸块放回队列。
                continue  # 新版残块重新排序后再安排查询。
            self.active = i, poly, polytopes[i]  # 本块已同步，可在下一轮逐顶点补齐同方案证明。
        # 每一块均已被有效割删除或被认证并集覆盖，才能耗尽队列。
        return None  # 队列耗尽意味着候选并集已被认证或排除。


def clip_polytope(vertices, constant, coefficient):  # 以一条有效负荷割裁剪凸多面体。
    """constant + coefficient·p ≥ 0；候选顶点由外近似自身产生。"""
    values = constant+vertices@coefficient  # 统一使用“余量≥0”的割方向。
    if np.all(values >= -1e-9):  # 全部顶点在保留侧或容差内，整个凸包均保留。
        return vertices  # 有效割没有改变当前区域时保留数组身份。
    if np.all(values < 0):  # 整个凸包都位于被删除的一侧。
        return np.empty((0, 3))  # 全部顶点严格在删除侧，返回空域。
    points = [p for p, value in zip(vertices, values) if value >= 0]  # 原有且仍满足割的顶点。
    for a, b in combinations(range(len(vertices)), 2):  # 检查所有点对，包含全部真实边；额外交点后续去掉。
        if values[a]*values[b] < 0:  # 两端异号时，线段与割平面有唯一交点。
            points.append(vertices[a]+(vertices[b]-vertices[a])*values[a]/(values[a]-values[b]))  # 线性插值解 a+b·p=0。
    return polytope_vertices(points)  # 重新取凸包，只留下裁剪后真实顶点。
