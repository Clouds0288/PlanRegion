"""可行规划域：Stage 1 全局支持割；Stage 2 最大认证间隙优先的目标空间分支定界。"""
from itertools import product
from pathlib import Path
from time import perf_counter

import numpy as np
from threadpoolctl import threadpool_limits

from Network.case33bw import Case33
from Network.four_bus_five_corridor import FourBus
from Network.jiangkou import Jiangkou
from model import DEFAULT_SOLVER_THREADS, PlanningEquations, PlanningModel
from region import (Cell, apply_disjunction, cell_gap, clip_polytope, cover_partition,
                    disjunction_directions, hull_distance,
                    scheme_hulls, select_cell, support_direction)
from plot import save_result

# Stage 0 配置：改这里，然后运行 python main.py。
NETWORK = 'case33'
LOAD_NODES = (18, 25, 33)       # Case33；二维可设为 (18,25)，其他节点恢复原始背景负荷。
UPGRADE_COUNT = 8
BUDGETS = (0., 2., 4.)         # Case33 相对投资单位；其他网络用 BUDGETS=None。
EPSILON_KW = 20.               # 全部外域到同方案内域并集的距离上界，不是扫描间距。
QUERY_TIME_LIMIT = 60.
CASE_TIME_LIMIT = 21600.       # 每档预算；两个阶段共享时限。
SOLVER_THREADS = DEFAULT_SOLVER_THREADS
MAX_SUPPORTS = 40
OBLIQUE = True                 # 斜向析取不够强时直接求正交退让量。
OUTPUT = Path(__file__).resolve().parent/'results'/'planning'


def build_continuous_region(network, budget, *, epsilon_kw=EPSILON_KW,
                            time_limit=CASE_TIME_LIMIT, query_time_limit=QUERY_TIME_LIMIT,
                            threads=SOLVER_THREADS, max_supports=MAX_SUPPORTS, oblique=OBLIQUE):
    """构造 I⊆D_budget⊆O；certified 表示 sup_{q∈O} dist∞(q,I)≤epsilon_kw。"""
    # 0.1 下闭包证明的适用条件；不满足时不能用本算法宣称全域认证。
    d = len(network.load_nodes)
    assert d in (2, 3) and epsilon_kw > 2e-4
    assert np.all(network.fixed_p >= 0.) and np.all(network.fixed_q >= 0.)
    assert np.all(network.q_ratio >= 0.) and np.all(network.r > 0.) and np.all(network.reactance >= 0.)
    assert np.all(network.vmin <= 1.) and np.all(network.vmax >= 1.)
    assert np.isfinite(network.power_limit) and network.power_limit > 0.
    start = perf_counter()
    deadline = start+time_limit
    masks = np.asarray(list(product((0., 1.), repeat=d)))
    axes = np.eye(d)
    equations = PlanningEquations(network, 'socp')
    problem = PlanningModel(equations, budget=budget, threads=threads,
                            strengthen=bool(np.all(network.required)))
    certificates, supports, branch_cuts = [], [], []
    status, max_gap_kw = 'unknown', np.inf

    # Stage 1：全局支持割，D⊆O1。每次都允许所有预算内建设方案参与。
    # 1.1 最大正方形/立方体与电源总负荷必要界：p≥0，Σp≤power_limit。
    stage1_vertices = masks*network.power_limit
    weights = np.ones(d)
    supports.append(dict(weights=weights, bound=network.power_limit))
    stage1_vertices = clip_polytope(stage1_vertices, network.power_limit, -weights)

    with problem.model:
        # 1.2 原点可行性：纯负荷下闭包中，原点已证不可行意味着整个域为空。
        problem.power.UB = np.zeros(d)
        answer = problem.solve(min(query_time_limit, max(0., deadline-perf_counter())), weights=np.zeros(d))
        problem.power.UB = np.full(d, network.power_limit)
        if answer is None:
            status, max_gap_kw = 'empty', 0.
            stage1_vertices = np.empty((0, d))
        elif answer['feasible']:
            certificates.append(answer)

        # 1.3 坐标方向先切，再选择当前凸包缺口最大的支持方向。
        # h(omega)=max_{x,p,y} omegaᵀp；只用全局上界 h_bar 切 omegaᵀp≤h_bar。
        excluded = []
        for iteration in range(d+max_supports):
            if status == 'empty' or perf_counter() >= deadline:
                break
            if iteration < d:
                weights = axes[iteration]
                # 下闭包保证轴向最优值可在其他负荷坐标为零时取得。
                problem.power.UB = weights*network.power_limit
            elif certificates:
                gap, weights = support_direction(stage1_vertices, certificates, excluded)
                if weights is None or gap <= epsilon_kw/2:
                    break
            else:
                break

            # 1.4 固定一个已知方案求连续热启动；此界不能生成全局割。
            incumbent = None
            if certificates:
                seed = max(certificates, key=lambda c: float(weights @ c['p']))
                problem.x.LB = problem.x.UB = seed['x']
                warm = problem.solve(min(2., max(0., deadline-perf_counter())), weights=weights, incumbent=seed)
                problem.x.LB, problem.x.UB = np.zeros(network.n_types), np.ones(network.n_types)
                if warm is not None and warm['feasible']:
                    certificates.append(warm)
                    incumbent = warm
            answer = problem.solve(min(query_time_limit, max(0., deadline-perf_counter())),
                                   weights=weights, incumbent=incumbent, gap_kw=epsilon_kw/10)
            problem.power.UB = np.full(d, network.power_limit)
            if answer is None:
                # 查询没有添加负荷下界；若完整模型不可行，则 D 为空。
                status, max_gap_kw = 'empty', 0.
                stage1_vertices = np.empty((0, d))
                break
            if answer['feasible']:
                certificates.append(answer)
            if answer['bound'] is not None:
                bound = answer['bound']
                supports.append(dict(weights=weights.copy(), bound=bound))
                stage1_vertices = clip_polytope(stage1_vertices, bound, -weights)
                problem.model.addConstr(weights @ problem.power <= bound)
            if answer['bound'] is None or answer['objective'] is None or answer['bound']-answer['objective'] > epsilon_kw/2:
                excluded.append(weights.copy())
            print(f'Stage 1 | budget={budget:g} | cut={len(supports)}', flush=True)

        # Stage 2：O2 为凸节点并集，保留非凸结构；内域只在同一 x 内取凸包。
        # 2.1 I_x=down(conv(已知同方案可行点))，I=union_x I_x。
        cells = [Cell(stage1_vertices)] if len(stage1_vertices) else []
        hulls = scheme_hulls(certificates)
        revision = len(certificates)
        while cells and hulls and perf_counter() < deadline:
            if revision != len(certificates):
                hulls = scheme_hulls(certificates)
                revision = len(certificates)

            # 2.2 Δ_s=min_x max_{v∈vertices(O_s)} dist∞(v,I_x)。
            # 选最大 Δ_s；所有 Δ_s≤epsilon 才能宣称全域认证。
            chosen, point_distances = select_cell(cells, hulls)
            cell = cells[chosen]
            max_gap_kw = cell.gap
            if max_gap_kw <= epsilon_kw:
                status = 'certified'
                break

            if point_distances.max() > epsilon_kw:
                # 2.3 查询离当前内域最远的节点顶点 q；先用附近已知方案扩张 I。
                target = cell.vertices[int(point_distances.argmax())]
                seeds = sorted(certificates, key=lambda c: np.maximum(target-c['p'], 0.).max())
                seen, incumbent = set(), None
                for seed in seeds:
                    key = tuple(seed['x'])
                    if key in seen:
                        continue
                    seen.add(key)
                    problem.x.LB = problem.x.UB = seed['x']
                    warm = problem.solve(min(2., max(0., deadline-perf_counter())), target=target,
                                         incumbent=seed, gap_kw=epsilon_kw/5)
                    problem.x.LB, problem.x.UB = np.zeros(network.n_types), np.ones(network.n_types)
                    if warm is not None and warm['feasible']:
                        certificates.append(warm)
                        if incumbent is None or warm['objective'] < incumbent['objective']:
                            incumbent = warm
                    if (incumbent is not None and np.maximum(target-incumbent['p'], 0.).max() <= epsilon_kw-1e-3
                            or len(seen) == 4 or perf_counter() >= deadline):
                        break
                if incumbent is not None and np.maximum(target-incumbent['p'], 0.).max() <= epsilon_kw-1e-3:
                    continue

                # 2.4 完整整数模型：rho*=min rho, Wp+rho·1≥Wq, rho≥0。
                # W=I 时就是 d(q)。正全局下界 L 给 OR_k W_k p≤W_k q-L。
                cut_normals = disjunction_directions(target, hulls) if oblique else axes
                answer = problem.solve(min(query_time_limit, max(0., deadline-perf_counter())),
                                       target=target, cut_normals=cut_normals if oblique else None,
                                       incumbent=incumbent, gap_kw=epsilon_kw/5)
                if answer is not None and answer['feasible']:
                    certificates.append(answer)
                # rho≈0 不说明 q 可行；改求真实坐标退让量，不混淆两个目标。
                if oblique and answer is not None and (answer['bound'] is None or answer['bound'] <= .01):
                    cut_normals = axes
                    answer = problem.solve(min(query_time_limit, max(0., deadline-perf_counter())),
                                           target=target, incumbent=incumbent, gap_kw=epsilon_kw/5)
                    if answer is not None and answer['feasible']:
                        certificates.append(answer)
                if answer is None:
                    status, max_gap_kw, cells = 'empty', 0., []
                    break
                if answer['bound'] is not None and answer['bound'] > .01:
                    threshold = cut_normals @ target-answer['bound']
                    branch_cuts.append(dict(target=target, cut_normals=cut_normals,
                                            threshold=threshold, bound=answer['bound']))
                    # 对所有相交节点传播；OR 析取不能当作 AND 线性约束加入模型。
                    cells = apply_disjunction(cells, threshold, cut_normals)
                elif not answer['feasible'] or np.maximum(target-answer['p'], 0.).max() > epsilon_kw-1e-3:
                    status = 'time_limit' if perf_counter() >= deadline else 'unknown'
                    break
            else:
                # 2.5 顶点分别覆盖不等于节点覆盖。沿 I_x^epsilon 分区，不删除点。
                choices = sorted(hulls, key=lambda pair: -np.sum(
                    hull_distance(cell.vertices, pair[1]) <= epsilon_kw-2e-5))
                children = None
                for key, facets in choices[:3]:
                    children = cover_partition(cell.vertices, facets, epsilon_kw-2e-5)
                    if children is not None:
                        break
                if children is None:
                    # 2.6 仅切到切面而未取得同维覆盖时，沿最长坐标二分。
                    axis = int(np.ptp(cell.vertices, axis=0).argmax())
                    middle = (cell.vertices[:, axis].min()+cell.vertices[:, axis].max())/2
                    children = [clip_polytope(cell.vertices, middle, -axes[axis]),
                                clip_polytope(cell.vertices, -middle, axes[axis])]
                cells[chosen:chosen+1] = [Cell(v, cell.depth+1, cell.gap, cell.scheme)
                                          for v in children if len(v)]
            print(f'Stage 2 | budget={budget:g} | cells={len(cells)} | gap<={max_gap_kw:.3f} kW', flush=True)

    # Stage 3：只保留域、切面及认证结论；不保存逐次求解状态或回放事件。
    # 3.1 用最终内域重新认证；时间耗尽不能用旧的较松界冒充最终证书。
    hulls = scheme_hulls(certificates)
    if cells and hulls:
        max_gap_kw = max(cell_gap(cell.vertices, hulls)[0] for cell in cells)
        if max_gap_kw <= epsilon_kw:
            status = 'certified'
    if status == 'unknown' and perf_counter() >= deadline:
        status = 'time_limit'
    # 3.2 从同方案半空间恢复内域顶点；保留数值收缩后的真实认证内域。
    inner = []
    for key, facets in hulls:
        vertices = masks*network.power_limit
        for row in facets:
            vertices = clip_polytope(vertices, -row[-1], -row[:-1])
        x = np.asarray(key, dtype=int)
        inner.append(dict(x=x, choice=network.decode_plan(x), cost=float(network.cost @ x), vertices=vertices))
    return dict(schema_version=4, method='socp', budget=budget, load_nodes=network.load_nodes,
                status=status, epsilon_kw=epsilon_kw, max_gap_kw=max_gap_kw,
                square_kw=network.power_limit, stage1_vertices=stage1_vertices,
                inner=inner, outer=[dict(vertices=cell.vertices) for cell in cells],
                supports=supports, branch_cuts=branch_cuts, seconds=perf_counter()-start)


def main():
    # 0.2 数据入口；数学模型与构域流程不依赖算例名称。
    network = {'case33': lambda: Case33(load_nodes=LOAD_NODES, upgrade_count=UPGRADE_COUNT),
               'fourbus': FourBus, 'jiangkou': Jiangkou}[NETWORK]()
    with threadpool_limits(limits=1):
        for budget in network.budgets if BUDGETS is None else BUDGETS:
            result = build_continuous_region(network, budget)
            save_result(result, OUTPUT/f'budget_{budget:g}')
            print(f"budget={budget:g} | {result['status']} | gap<={result['max_gap_kw']:.3f} kW | {result['seconds']:.1f}s")


if __name__ == '__main__':
    main()
