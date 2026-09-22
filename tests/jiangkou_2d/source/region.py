"""连续域状态、候选点筛选、裁剪与并集；不跨建设方案取凸包。"""
import numpy as np  # 网格索引和成块分类使用同一坐标顺序。
from itertools import combinations, product
from types import SimpleNamespace
from scipy.spatial import ConvexHull, QhullError

# 连续几何在公共评价箱归一化后的坐标中计算；与采样网格无关。
GEOMETRY_TOL = 1e-8


def _convex_hull(points):
    """近共面输入失败时用可逆坐标变换重算，保留原始点、半空间和体积。"""
    try:
        return ConvexHull(points)
    except QhullError:
        center = points[0]
        _, scales, basis = np.linalg.svd(points-center, full_matrices=False)
        if np.any(scales <= 0.):
            raise
        hull = ConvexHull(((points-center)@basis.T)/scales)
        normal = (hull.equations[:, :-1]/scales)@basis
        lengths = np.linalg.norm(normal, axis=1)
        equations = np.c_[normal/lengths[:, None],
                           (hull.equations[:, -1]-normal@center)/lengths]
        return SimpleNamespace(vertices=hull.vertices, simplices=hull.simplices,
                               equations=equations, volume=float(hull.volume*np.prod(scales)))


def polytope_vertices(points):
    """任意仿射维数的凸包顶点，保持原始交点坐标。"""
    points = np.asarray(points, dtype=float).reshape(-1, 2)
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
    """返回 A z + b <= 0；低维集用成对不等式表示仿射等式。"""
    points = np.asarray(points)
    if not len(points):
        return np.array([[0., 0., 1.]])
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
        eq = np.empty((0, 3))
    normal = basis[rank:]
    eq = np.vstack([eq, np.c_[normal, -normal@center], np.c_[-normal, normal@center]])
    _, indices = np.unique(np.round(eq, 10), axis=0, return_index=True)
    return eq[np.sort(indices)]


def contains(points, equations, tolerance=GEOMETRY_TOL):
    return np.all(np.atleast_2d(points)@equations[:, :2].T+equations[:, 2] <= tolerance, axis=1)


def clip_polytope(vertices, constant, coefficient):
    """用 constant + coefficient*z >= 0 裁剪连续凸域。"""
    vertices = np.asarray(vertices).reshape(-1, 2)
    if not len(vertices):
        return vertices
    coefficient = np.asarray(coefficient)
    norm = np.linalg.norm(coefficient)
    if norm < 1e-20:
        return vertices if constant >= -1e-12 else np.empty((0, 2))
    values = (constant+vertices@coefficient)/norm
    if values.min() >= -1e-11:
        return vertices
    if values.max() < -1e-11:
        return np.empty((0, 2))
    points = list(vertices[values >= -1e-11])
    if polytope_volume(vertices) > 0:
        edges = {tuple(sorted(edge)) for face in _convex_hull(vertices).simplices
                 for edge in combinations(face, 2)}
    else:
        edges = combinations(range(len(vertices)), 2)
    for a, b in edges:
        if values[a]*values[b] < 0:
            points.append(vertices[a]+(vertices[b]-vertices[a])*values[a]/(values[a]-values[b]))
    return polytope_vertices(points)


def initial_polytope(bounds, total):
    return clip_polytope(np.array(list(product((0., 1.), repeat=2))), total, -np.asarray(bounds))


def polytope_volume(poly):
    poly = np.asarray(poly).reshape(-1, 2)
    if len(poly) < 3 or np.linalg.matrix_rank(poly-poly[0], tol=1e-10) < 2:
        return 0.
    return float(_convex_hull(poly).volume)


def subtract_polytope(poly, equations):
    """返回 poly 减去一个凸域的互不重叠内部的凸分块。"""
    if not len(poly):
        return []
    values = poly@equations[:, :2].T+equations[:, 2]
    if np.all(values <= 1e-10):
        return []
    if np.any(values.min(axis=0) >= -1e-10):
        return [poly]
    pieces = []
    for eq in equations:
        values = poly@eq[:2]+eq[2]
        if values.max() <= 1e-10:
            continue
        outside = clip_polytope(poly, eq[2], eq[:2])
        if polytope_volume(outside) > 1e-15:
            pieces.append(outside)
        poly = clip_polytope(poly, -eq[2], -eq[:2])
        if polytope_volume(poly) <= 1e-15:
            break
    return pieces


def union_volume(polytopes):
    """逐域扣除已经计入的部分，避免把重叠方案体积相加。"""
    previous, volume = [], 0.
    for poly in sorted(polytopes, key=polytope_volume, reverse=True):
        if polytope_volume(poly) <= 1e-15:
            continue
        pieces = [poly]
        for eq in previous:
            pieces = [part for p in pieces for part in subtract_polytope(p, eq)]
            if not pieces:
                break
        volume += sum(polytope_volume(p) for p in pieces)
        previous.append(halfspaces(poly))
    return volume


def surface(poly, bounds):
    """用于 HTML 的真实凸多面体表面，退化集也保留坐标。"""
    poly = np.asarray(poly).reshape(-1, 2)
    faces = _convex_hull(poly).simplices.tolist() if polytope_volume(poly) > 0 else []
    return dict(vertices=(poly*np.asarray(bounds)).tolist(), faces=faces)


def classify_orthant(states, index, status):  # 纯负荷单调域中，用一个已证点分类一块网格中心。
    """调用者须验证非负负荷、正阻抗和 vmax≥根电压；不跨方案构造凸包。"""
    if status == 1:  # 同一建设方案在更低负荷下仍可行。
        states[tuple(slice(0,int(i)+1) for i in index)] = 1  # 只覆盖坐标逐项不大于认证点的网格中心。
    elif status == -1:  # 若所有预算内方案在该点不可行，更高负荷也不可能可行。
        states[tuple(slice(int(i),None) for i in index)] = -1


def split_grid_box(lower, upper):  # 把含网格中心的闭索引盒沿最长轴分成两块。
    axis = int(np.argmax(upper-lower))  # 索引尺度对应相同的每轴网格分辨率。
    middle = (lower[axis]+upper[axis])//2  # 整数二分不会丢失或重复中心点。
    left, right = upper.copy(), lower.copy()  # 保留其他两个坐标的完整索引范围。
    left[axis], right[axis] = middle, middle+1  # 左右子盒不重叠，且并集等于父盒。
    return ((lower,left),(right,upper))



class RegionState:
    """固定方案内外多面体、候选点筛选与连续并集；坐标归一化到评价箱。"""

    def __init__(self, bounds, total_bound, tau, cuts=()):
        self.bounds = np.asarray(bounds, dtype=float)
        self.total_bound, self.tau = float(total_bound), float(tau)
        if not np.isfinite(self.tau) or not 0 <= self.tau < 1:
            raise ValueError('tau must be in [0, 1)')
        self.cuts = [np.asarray(c) for c in cuts]
        self.records = {}
        self.revision = 0
        self._surfaces = {}

    def add_scheme(self, x, choice, cost):
        key = tuple(x)
        if key in self.records:
            return False
        outer = initial_polytope(self.bounds, self.total_bound)
        for cut in self.cuts:
            outer = clip_polytope(outer, cut[0]+cut[3:]@x, cut[1:3]*self.bounds)
        self.records[key] = dict(x=np.asarray(x, dtype=int), choice=list(choice),
                                 cost=float(cost), inner=np.empty((0, 2)), outer=outer,
                                 inner_equations=None, inner_box=None)
        return True

    def add_point(self, x, point):
        key = tuple(x)
        row = self.records[key]
        points = np.asarray(point, dtype=float).reshape(-1, 2)
        if not len(points) or (len(row['inner']) and contains(points, self.inner_equations(key)).all()):
            return False
        updated = polytope_vertices(np.vstack([row['inner'], points]))
        if np.array_equal(updated, row['inner']):
            return False
        row['inner'] = updated
        row['inner_equations'] = row['inner_box'] = None
        self._surfaces.pop(key, None)
        self.revision += 1
        return True

    def inner_equations(self, key):
        row = self.records[tuple(key)]
        if row['inner_equations'] is None:
            row['inner_equations'] = halfspaces(row['inner'])
            if len(row['inner']):
                row['inner_box'] = (row['inner'].min(axis=0), row['inner'].max(axis=0))
        return row['inner_equations']

    def covering_schemes(self, points, preferred=None):
        """每点返回一个支撑方案；包围盒只作快速筛选，接受仍须逐面检查。"""
        points = np.asarray(points, dtype=float).reshape(-1, 2)
        owners = [None]*len(points)
        keys = list(self.records)
        if preferred is not None and tuple(preferred) in self.records:
            keys.remove(tuple(preferred))
            keys.insert(0, tuple(preferred))
        for key in keys:
            row = self.records[key]
            if not len(row['inner']):
                continue
            eq = self.inner_equations(key)
            lower, upper = row['inner_box']
            indices = [i for i, owner in enumerate(owners) if owner is None]
            if not indices:
                break
            indices = np.asarray(indices)
            boxed = np.all((points[indices] >= lower-GEOMETRY_TOL)
                           & (points[indices] <= upper+GEOMETRY_TOL), axis=1)
            indices = indices[boxed]
            for i in indices[contains(points[indices], eq)]:
                owners[i] = key
        return owners

    def next_point(self, x, skipped=None):
        row = self.records[tuple(x)]
        targets = (1-self.tau)*row['outer']
        owners = self.covering_schemes(targets, preferred=x)
        if skipped is not None:
            for point, owner in zip(targets, owners):
                if owner is not None and owner != tuple(x):
                    skipped(point, owner)
        indices = [i for i, owner in enumerate(owners) if owner is None]
        if not len(indices):
            return None
        choose = min if not len(row['inner']) else max
        return targets[choose(indices, key=lambda j: float(targets[j].sum()))]

    def witness_support(self, x, point):
        """用同方案收缩外域内的至多三个点支撑未覆盖见证，不跨方案取凸包。

        从凸域重心向各面作单纯形剖分。它只选择待认证点，不赋予可行性。
        """
        poly = (1-self.tau)*self.records[tuple(x)]['outer']
        if not len(poly):
            return []
        center = poly.mean(axis=0)
        rank = np.linalg.matrix_rank(poly-center, tol=1e-10)
        if rank == 0:
            return [poly[0]]
        basis = np.linalg.svd(poly-center, full_matrices=False)[2][:rank]
        coordinates, target = (poly-center)@basis.T, (point-center)@basis.T
        if rank == 1:
            return poly[[coordinates[:, 0].argmin(), coordinates[:, 0].argmax()]]
        for face in _convex_hull(coordinates).simplices:
            try:
                weights = np.linalg.solve(coordinates[face].T, target)
            except np.linalg.LinAlgError:
                continue
            if weights.min() >= -1e-9 and weights.sum() <= 1+1e-9:
                weights = np.r_[1-weights.sum(), weights]
                points = np.vstack([center, poly[face]])
                return points[np.argsort(-weights)[weights[np.argsort(-weights)] > 1e-12]]
        return []

    def apply_cut(self, cut):
        self.cuts.append(np.asarray(cut))
        for key, row in self.records.items():
            updated = clip_polytope(row['outer'], cut[0]+cut[3:]@row['x'], cut[1:3]*self.bounds)
            if not np.array_equal(updated, row['outer']):
                row['outer'] = updated
                self._surfaces.pop(key, None)

    def inner_halfspaces(self):
        return [self.inner_equations(key) for key, r in self.records.items() if len(r['inner'])]

    @property
    def progress(self):
        # 新证书可能替换旧顶点而不改变顶点数；用版本判断实际进展。
        return len(self.cuts), self.revision

    def geometry(self):
        for key, row in self.records.items():
            if key not in self._surfaces:
                self._surfaces[key] = dict(choice=row['choice'], cost=row['cost'],
                                          inner=surface(row['inner'], self.bounds),
                                          outer=surface(row['outer'], self.bounds))
        return [self._surfaces[key] for key in self.records]

    def finish(self, certified):
        """由全局覆盖结论构造外包络；未知时保留整个评价箱内的候选域。"""
        records = list(self.records.values())
        inner = [r['inner'] for r in records if len(r['inner'])]
        outer = []
        if certified:
            for poly in inner:
                eq = halfspaces(poly)
                center = poly.mean(axis=0)
                radius = float(np.min(-eq[:, :2]@center-eq[:, 2]))
                if radius > 1e-10:
                    # 仿射外扩包含每个面外移 GEOMETRY_TOL 的集合，避免近共面裁剪交点。
                    envelope = (center+(1+GEOMETRY_TOL/radius)*(poly-center))/(1-self.tau)
                    for axis in np.eye(2):
                        envelope = clip_polytope(envelope, 0., axis)
                        envelope = clip_polytope(envelope, 1., -axis)
                    envelope = clip_polytope(envelope, self.total_bound, -self.bounds)
                else:
                    envelope = initial_polytope(self.bounds, self.total_bound)
                    for face in eq:
                        envelope = clip_polytope(envelope, (GEOMETRY_TOL-face[2])/(1-self.tau), -face[:2])
                outer.append(envelope)
        else:
            outer = [initial_polytope(self.bounds, self.total_bound)]
        scale = np.prod(self.bounds)
        inner_volume, outer_volume = union_volume(inner)*scale, union_volume(outer)*scale
        return dict(schemes=len(records), inner_volume=inner_volume, outer_volume=outer_volume,
                    volume_gap=(outer_volume-inner_volume)/outer_volume if outer_volume else 0.,
                    geometry=self.geometry(), inner=[surface(p, self.bounds) for p in inner],
                    outer=[surface(p, self.bounds) for p in outer], cuts=[c.tolist() for c in self.cuts],
                    certificates=[dict(x=r['x'].tolist(), choice=r['choice'], cost=r['cost'],
                                       inner=r['inner'].tolist(), outer=r['outer'].tolist()) for r in records])


def sample_region(result, points, bounds):
    """连续域求出后才作采样；薄层保留为未知，采样不影响几何。"""
    points = np.atleast_2d(points)/np.asarray(bounds)
    inside, possible = np.zeros(len(points), bool), np.zeros(len(points), bool)
    for poly in result['inner']:
        inside |= contains(points, halfspaces(np.asarray(poly['vertices'])/bounds))
    for poly in result['outer']:
        possible |= contains(points, halfspaces(np.asarray(poly['vertices'])/bounds))
    return np.where(inside, 1, np.where(possible, 0, -1)).astype(np.int8)
