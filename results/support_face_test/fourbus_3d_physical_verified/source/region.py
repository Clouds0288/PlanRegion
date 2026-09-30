"""连续域几何与全局搜索：逐网架认证、切割、射线补点和符号分区。"""
from copy import deepcopy
from itertools import combinations, product

import numpy as np
from scipy.spatial import ConvexHull
from threadpoolctl import threadpool_limits

from model import GridPhysics, PortPhysics, PortSubProblem, MasterProblem, RemainingRegionModel, port_bounds, ray_support
from monitor import RunMonitor, RegionTimeout, _union
from plot import union_volume

# 连续几何在公共评价箱归一化后的坐标中计算；与采样网格无关。
GEOMETRY_TOL = 1e-8


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


def initial_polytope(bounds, total, axis_bounds):
    vertices = np.array(list(product((0., 1.), repeat=len(bounds))))*axis_bounds/bounds
    return clip_polytope(vertices, total, -np.asarray(bounds))


def polytope_volume(poly):
    poly = np.asarray(poly)
    d = poly.shape[1] if poly.ndim == 2 else 3
    poly = poly.reshape(-1, d)
    if len(poly) < d+1 or np.linalg.matrix_rank(poly-poly[0], tol=1e-10) < d:
        return 0.
    coordinates = (poly-poly[0])@np.linalg.svd(poly-poly[0], full_matrices=False)[2].T
    scale = np.linalg.norm(coordinates, axis=0)
    return float(ConvexHull(coordinates/scale).volume*np.prod(scale))


class RegionState:
    """固定方案内外多面体与连续并集；坐标归一化到公共评价箱。"""

    def __init__(self, bounds, total_bound, tau, cuts=()):
        self.bounds = np.asarray(bounds, dtype=float)
        self.axis_bounds = self.bounds.copy()
        self.total_bound, self.tau = float(total_bound), float(tau)
        self.cuts = [np.asarray(c) for c in cuts]
        self.records = {}

    def add_scheme(self, x, choice, cost):
        key = tuple(x)
        if key in self.records:
            return False
        outer = initial_polytope(self.bounds, self.total_bound, self.axis_bounds)
        d = len(self.bounds)
        for cut in self.cuts:
            outer = clip_polytope(outer, cut[0]+cut[1+d:]@x, cut[1:1+d]*self.bounds)
        self.records[key] = dict(x=np.asarray(x, dtype=int), choice=dict(choice),
                                 cost=float(cost), inner=np.empty((0, d)), outer=outer,
                                 inner_equations=None)
        return True

    def add_point(self, x, point):
        key = tuple(x)
        row = self.records[key]
        points = np.asarray(point, dtype=float).reshape(-1, len(self.bounds))
        if not len(points) or (len(row['inner']) and contains(points, self.inner_equations(key)).all()):
            return False
        updated = polytope_vertices(np.vstack([row['inner'], points]))
        if np.array_equal(updated, row['inner']):
            return False
        row['inner'] = updated
        row['inner_equations'] = None
        return True

    def tighten_bounds(self, axis_bounds, total_bound):
        """MP2 上界以 kW 输入；同步裁剪已登记网架，供新网架复用。"""
        self.axis_bounds = np.minimum(self.axis_bounds, axis_bounds)
        self.total_bound = min(self.total_bound, total_bound)
        for row in self.records.values():
            outer = clip_polytope(row['outer'], self.total_bound, -self.bounds)
            for axis, limit in zip(np.eye(len(self.bounds)), self.axis_bounds/self.bounds):
                outer = clip_polytope(outer, limit, -axis)
            row['outer'] = outer

    def inner_equations(self, key):
        row = self.records[tuple(key)]
        if row['inner_equations'] is None:
            row['inner_equations'] = halfspaces(row['inner'])
        return row['inner_equations']

    def apply_cut(self, cut):
        self.cuts.append(np.asarray(cut))
        d = len(self.bounds)
        for row in self.records.values():
            row['outer'] = clip_polytope(row['outer'], cut[0]+cut[1+d:]@row['x'], cut[1:1+d]*self.bounds)

    def finish(self, certified):
        """由全局覆盖结论构造外包络；未知时保留整个评价箱内的候选域。"""
        records = list(self.records.values())
        inner = [r['inner'] for r in records if len(r['inner'])]
        outer = []
        if certified:
            for poly in inner:
                eq = halfspaces(poly)
                center = poly.mean(axis=0)
                radius = float(np.min(-eq[:, :-1]@center-eq[:, -1]))
                if radius > 1e-10:
                    # 仿射外扩包含每个面外移 GEOMETRY_TOL 的集合，避免近共面裁剪交点。
                    envelope = (center+(1+GEOMETRY_TOL/radius)*(poly-center))/(1-self.tau)
                    for axis, limit in zip(np.eye(len(self.bounds)), self.axis_bounds/self.bounds):
                        envelope = clip_polytope(envelope, 0., axis)
                        envelope = clip_polytope(envelope, limit, -axis)
                    envelope = clip_polytope(envelope, self.total_bound, -self.bounds)
                else:
                    envelope = initial_polytope(self.bounds, self.total_bound, self.axis_bounds)
                    for face in eq:
                        envelope = clip_polytope(envelope, (GEOMETRY_TOL-face[-1])/(1-self.tau), -face[:-1])
                outer.append(envelope)
        else:
            outer = [initial_polytope(self.bounds, self.total_bound, self.axis_bounds)]
        return dict(inner=[dict(choice=r['choice'].copy(), cost=r['cost'],
                                vertices=r['inner']*self.bounds) for r in records if len(r['inner'])],
                    outer=[dict(vertices=p*self.bounds) for p in outer])


def union_measure(polytopes, d):
    """二维面积 / 三维体积；重叠网架只计一次。"""
    return _union(polytopes).area if d == 2 else union_volume(polytopes)


def register_power(power, powers, point_tol):
    """以 kW 最大坐标差匹配固定代表；不舍入坐标，也不沿近点链移动代表。"""
    for index, representative in enumerate(powers):
        if np.max(np.abs(power-representative)) <= point_tol:
            return index
    powers.append(np.asarray(power, dtype=float).copy())
    return len(powers)-1


def stage_candidates(region, x):
    """只检查当前 N_x 的顶点；其他网架的认证边界不参与局部选点。"""
    outer = region.records[tuple(x)]['outer']
    return outer[~contains(outer, region.inner_equations(x))]


def ray_gain(inner, point):
    """射线可行端点在线段上，新增测度不超过把目标点直接收入凸包的增量。"""
    return polytope_volume(np.vstack([inner, point]))-polytope_volume(inner)


def coverage_halfspaces(region):
    """从大到小保留认证域；被单个已保留凸包包含的排除项无需重复建模。"""
    keys = [key for key, row in region.records.items() if len(row['inner'])]
    keys.sort(key=lambda key: polytope_volume(region.records[key]['inner']), reverse=True)
    kept, equations = [], []
    for key in keys:
        if any(contains(region.records[key]['inner'], eq).all() for eq in equations):
            continue
        kept.append(key)
        equations.append(region.inner_equations(key))
    return kept, equations


def build_sequential_region(network, *, budget, monitor, seconds=50., threads=4,
                            threshold=.02, patience=3, tau=.005,
                            numeric_focus=0, point_tol=.01,
                            ray_threshold=1e-4, mode=0, sign=None):
    """新网架先遍历 G 做射线；每阶段小割后补边界，随即由 physical 搜索查漏。"""
    d = len(network.load_nodes)
    sign = np.ones(d, int) if sign is None else np.asarray(sign)
    bounds = port_bounds(network) if mode else np.full(d, network.power_limit)
    region = RegionState(bounds, float(bounds.sum()) if mode else network.power_limit, tau)
    monitor.begin(network, 'socp', budget, region, seconds)
    equations = PortPhysics(network, sign) if mode else GridPhysics(network, 'socp')
    oracle = PortSubProblem(equations, threads=threads, numeric_focus=numeric_focus)
    cache = {}
    powers, applied = [], set()
    initialized = set()
    counts = dict(initial=0, sp=0, cuts=0, ray=0, global_search=0)
    status, certified, coverage = 'time_limit', False, None
    x = None
    stage = 0
    monitor._emit('settings', phase='初始化', threshold=threshold, patience=patience,
                  tau=tau, numeric_focus=numeric_focus, point_tol=point_tol, ray_threshold=ray_threshold,
                  validation_note='独立 SOCP 扫描结果在此显示')

    def register(answer):
        selection = answer['x']
        added = region.add_scheme(selection, network.decode_plan(selection), network.cost_offset+network.cost@selection)
        region.add_point(selection, answer['p']/bounds)
        if added:
            if np.any(network.fixed_p) or np.any(network.fixed_q):
                # Case33 有固定背景负荷：原点也必须实际认证。
                power = np.zeros(d)
                monitor._emit('origin_start', phase='零接入认证',
                              active_scheme=monitor._scheme(selection), **monitor._geometry(region))
                monitor.sp_start(selection, power, oracle.calls+1)
                checked = oracle.solve(selection, power, time_limit=monitor.remaining(seconds), score_only=True)
                monitor.sp_end(checked)
                if checked['feasible']:
                    region.add_point(selection, power)
            else:
                # FourBus 无固定负荷，零潮流和单位电压解析可行。
                region.add_point(selection, np.zeros(d))
        return selection

    def check(point):
        # 每个点都读取最新凸包：本批前面的端点认证后，中间点立即免求解。
        if contains([point], region.inner_equations(x))[0]:
            return None
        power = point*bounds
        if mode:
            # 只认证实际内移点；原外顶点仍交由完整物理模型查漏。
            anchor = region.records[tuple(x)]['inner'].mean(axis=0)*bounds
            power += min(.5, point_tol/(4*np.max(np.abs(anchor-power))))*(anchor-power)
        index = register_power(power, powers, point_tol)
        power = powers[index]
        key = tuple(x), index
        if key in applied or contains([power/bounds], region.inner_equations(x))[0]:
            return None
        if key not in cache:
            monitor.sp_start(x, power, oracle.calls+1)
            cache[key] = oracle.solve(x, power, time_limit=monitor.remaining(seconds), score_only=True)
            monitor.sp_end(cache[key])
        answer = cache[key]
        if answer['feasible']:
            if region.add_point(x, power/bounds):
                monitor.updated(region, 'feasible', x=x, power=power, checked=answer)
            return None
        return index

    def add_ray(power, phase, initial_sweep=False):
        inner = region.records[tuple(x)]['inner']
        gain = ray_gain(inner, power/bounds)
        gain_threshold = ray_threshold*polytope_volume(inner)
        # 低维认证集先补足维度；满维凸包才使用面积/体积收益筛选。
        if not initial_sweep and gain_threshold > 0. and gain < gain_threshold:
            monitor._emit('ray_skip', phase='射线筛选', ray=None, ray_target=power,
                          ray_gain_bound=gain*np.prod(bounds), ray_gain_threshold=gain_threshold*np.prod(bounds))
            return
        anchor = inner.mean(axis=0)*bounds
        monitor._emit('ray_start', phase=phase,
            ray_gain_bound=gain*np.prod(bounds), ray_gain_threshold=gain_threshold*np.prod(bounds),
            ray=dict(scheme=monitor._scheme(x), anchor=anchor, target=power, p=None))
        answer = ray_support(equations, budget, x, anchor, power,
                             threads=threads, time_limit=monitor.remaining(seconds), numeric_focus=numeric_focus)
        counts['ray'] += 1
        region.add_point(x, answer['p']/bounds)
        monitor._emit('ray_end', phase=phase,
            ray=dict(scheme=monitor._scheme(x), anchor=anchor, target=power, p=answer['p']),
            ray_fraction=answer['ray_fraction'],
            **monitor._geometry(region))

    try:
        # 1. d 个轴向 MP2 和一个总负荷 MP2；后者给出初始网架 A。
        for index, direction in enumerate([*np.eye(d), np.ones(d)]):
            monitor.initializing(index)
            problem = MasterProblem(equations, budget=budget, direction=direction, threads=threads)
            with problem.model:
                counts['initial'] += 1
                answer = problem.solve(time_limit=monitor.remaining(seconds))
            if answer is None:
                region.tighten_bounds(np.zeros(d), 0.)
                status, certified = 'empty', True
                break
            axis_bounds = region.axis_bounds.copy()
            if index < d:
                axis_bounds[index] = min(axis_bounds[index], answer['bound'])
            region.tighten_bounds(axis_bounds, answer['bound'] if index == d else region.total_bound)
            monitor._emit('initial_bounds', phase='初始化', **monitor._geometry(region))
            if mode or index == d:
                x = register(answer)
                monitor.seed(answer, region)

        while x is not None:
            # 2. 新网架先遍历 G 顶点做射线认证；每个网架仅执行一次。
            monitor.remaining(seconds)
            stage += 1
            small_cuts = 0
            monitor._emit('scheme_start', phase='单网架构域', stage=stage,
                active_scheme=monitor._scheme(x), small_cuts=0, area_ratio=None,
                ray=None, sp_point=None, global_point=None, **monitor._geometry(region))
            if tuple(x) not in initialized:
                monitor._emit('initial_sweep_start', phase='首轮射线认证')
                for point in initial_polytope(bounds, region.total_bound, region.axis_bounds):
                    if not contains([point], region.inner_equations(x))[0]:
                        add_ray(point*bounds, '首轮射线认证', initial_sweep=True)
                initialized.add(tuple(x))
                monitor._emit('initial_sweep_end', phase='首轮射线认证', ray=None)
            # 3. 只评分当前 N_x 的顶点；收完可行点后，仅为最大 eta 点生成共享割。
            reason = 'covered'
            while True:
                monitor.remaining(seconds)
                points = stage_candidates(region, x)
                if not len(points):
                    break
                monitor.selecting(len(points), False)
                pending = set()
                for point in points:
                    index = check(point)
                    if index is not None:
                        pending.add(index)
                if not len(stage_candidates(region, x)):
                    break
                pending = [index for index in pending
                           if not contains([powers[index]/bounds], region.inner_equations(x))[0]]
                if not pending:
                    reason = 'point_resolution'
                    break
                index = max(pending, key=lambda i: (cache[(tuple(x), i)]['eta'],
                                                    float(powers[i].sum()), tuple(powers[i])))
                power, checked = powers[index], cache[(tuple(x), index)]
                monitor._emit('cut_start', phase='获选点生成割',
                              sp_point=dict(scheme=monitor._scheme(x), p=power), eta=checked['eta'])
                checked['cut'] = oracle.generate_cut(x, power, checked['cone_normals'],
                                                     time_limit=monitor.remaining(seconds))

                # 冻结 G'；μ(N\G')=μ(N∪G')−μ(G')，二维用面积、三维用体积。
                inner = [row['inner'] for row in region.records.values()]
                inner_area = union_measure(inner, d)
                before = region.records[tuple(x)]['outer']
                region.apply_cut(checked['cut'])
                applied.add((tuple(x), index))
                after = region.records[tuple(x)]['outer']
                denominator = polytope_volume(before) if stage == 1 else inner_area
                removed_area = (polytope_volume(before)-polytope_volume(after) if stage == 1 else
                                union_measure([before, *inner], d)-union_measure([after, *inner], d))
                area_ratio = removed_area/denominator if denominator > 0. else None
                small_cuts = small_cuts+1 if area_ratio is not None and area_ratio < threshold else 0
                counts['cuts'] += 1
                monitor.updated(region, 'cut', x=x, power=power, checked=checked)
                monitor._emit('cut_measure', phase='单网架构域', area_ratio=area_ratio,
                    small_cuts=small_cuts)

                # 加割后沿当前方向补可行点；连续小割只结束局部阶段。
                add_ray(power, '射线补点')
                if small_cuts >= patience:
                    reason = 'area_stagnation'
                    break

            # 4. 补一轮最新 N_x 顶点；实际可行射线点加入 N'_x，不能跨网架取凸包。
            monitor._emit('boundary_start', phase='边界补充', ray=None, sp_point=None)
            targets = {}
            for point in stage_candidates(region, x):
                index = register_power(point*bounds, powers, point_tol)
                targets.setdefault(index, point)
            while targets:
                inner = region.records[tuple(x)]['inner']
                index = max(targets, key=lambda i: ray_gain(inner, targets[i]))
                point = targets.pop(index)
                if not contains([point], region.inner_equations(x))[0]:
                    add_ray(point*bounds, '边界补充')
            monitor._emit('boundary_end', phase='边界补充', ray=None)
            monitor._emit('scheme_end', phase='网架阶段结束', stage_reason=reason, sp_point=None)

            # 5. 在 G 去掉容许扩边的 G' 后查漏；x、p、y 自由，见证带完整物理可行解。
            monitor.global_start(0)
            _, inner_halfspaces = coverage_halfspaces(region)
            problem = RemainingRegionModel(equations, budget, bounds, region.total_bound,
                region.cuts, inner_halfspaces, tau, axis_bounds=region.axis_bounds,
                threads=threads)
            with problem.model:
                counts['global_search'] += 1
                answer = problem.solve(GEOMETRY_TOL, time_limit=monitor.remaining(seconds))
            coverage = answer['bound']
            if not answer['complete']:
                x = register(answer)
            monitor.global_end(answer, region)
            if answer['complete']:
                status, certified = 'certified', True
                break
    except (RegionTimeout, TimeoutError):
        status = 'time_limit'

    counts['sp'] = oracle.calls
    result = dict(status=status, certified=certified, axis_bounds=region.axis_bounds,
        coverage_bound=coverage, counts=counts, timing=monitor.timing(), **region.finish(certified))
    monitor.finish(result, region)
    return result


def build_region(network, *, budget, monitor, mode=1, seconds=50., threads=4,
                 tau=.005, threshold=.02, patience=3, point_tol=.01, ray_threshold=1e-4):
    """在各符号分区独立构域；只合并有物理证书的内域与安全外包络。"""
    d = len(network.load_nodes)
    signs = list(product((1, -1), repeat=d)) if mode else [(1,)*d]
    bounds = port_bounds(network) if mode else np.full(d, network.power_limit)
    monitor.envelopes = {''.join('+' if s > 0 else '-' for s in sign):
        [np.asarray(list(product((0., 1.), repeat=d)))*bounds*sign] for sign in signs}
    monitor._emit('start', phase='初始化', mode=mode, network=network.name, load_nodes=network.load_nodes,
        budget=budget, method='socp', algorithm=monitor.algorithm, time_limit=seconds,
        cost_unit=network.cost_unit, initial_plan=network.initial_plan, bounds=bounds,
        axis_bounds=bounds, axis_lower=-bounds if mode else np.zeros(d),
        threshold=threshold, patience=patience, tau=tau, point_tol=point_tol, ray_threshold=ray_threshold,
        schemes={}, cut_history={}, global_outer=[p for polys in monitor.envelopes.values() for p in polys])
    partitions = []
    with threadpool_limits(limits=1):
        for j, sign in enumerate(signs):
            local = RunMonitor(parent=monitor, sign=sign, clock=monitor.clock, algorithm=monitor.algorithm)
            remaining = max(0., (seconds-monitor.timing()['total_seconds'])/(len(signs)-j))
            answer = build_sequential_region(deepcopy(network), budget=budget, monitor=local,
                mode=mode, sign=sign, seconds=remaining, threads=threads, tau=tau,
                threshold=threshold, patience=patience, point_tol=point_tol, ray_threshold=ray_threshold)
            partitions.append((sign, answer))
            print(f'{network.name} {d}D {sign}: {answer["status"]}, '
                  f'{answer["timing"]["total_seconds"]:.2f}s, {answer["counts"]}', flush=True)
    certified = all(answer['certified'] for _, answer in partitions)
    result = dict(status='certified' if certified else 'time_limit', certified=certified,
        timing=monitor.timing(), inner=[], outer=[],
        counts={key: sum(answer['counts'][key] for _, answer in partitions)
                for key in ('initial', 'sp', 'cuts', 'ray', 'global_search')})
    for sign, answer in partitions:
        for key in ('inner', 'outer'):
            result[key].extend({**row, 'sign': sign, 'vertices': np.asarray(row['vertices'])*sign}
                               for row in answer[key])
    vertices = np.vstack([row['vertices'] for row in result['outer']])
    result.update(axis_lower=vertices.min(axis=0) if mode else np.zeros(d), axis_bounds=vertices.max(axis=0))
    monitor._emit('region_end', phase='构域完成' if certified else '构域停止', result=result,
                  coverage_complete=certified, active_scheme=None, ray=None, sp_point=None, global_point=None)
    return result
