"""江口测试（分支 jiangkou-test，主线代码只读）：Jiangkou 网架自检、OBBT 计时基准、二维全网与等值网的构域对比。

二维端口为左下 B000078 × 右上 B000042（Network/jiangkou.py 的 LOAD_NODES），预算单位万元。
  python -X utf8 jiangkou_test.py check     规模、原点可行性、端口单点承载力（未紧化 SOCP 射线：现状 / 路径全换）
  python -X utf8 jiangkou_test.py obbt      OBBT 计时：全网 / 等值网 × 现状 / 全换 × ++ / -- × 线程数，另测中心射线 MISOCP
  python -X utf8 jiangkou_test.py rays      参考：x 自由的射线 MISOCP（每象限 RAYS 个方向），全网与等值网 × 两档预算
  python -X utf8 jiangkou_test.py region    二维 RCUT：全网与等值网 × 两档预算，各存回放
  python -X utf8 jiangkou_test.py compare   对比图、耗时表与 summary.json
输出在 results/jiangkou_test/（不入库）；回放用 python -X utf8 monitor.py results/jiangkou_test/<名>.json.gz。
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
from copy import deepcopy
from functools import cache
from itertools import product
import json
import multiprocessing
from pathlib import Path
import time
from types import SimpleNamespace

import numpy as np
from gurobipy import GRB
from shapely.geometry import LineString
from threadpoolctl import threadpool_limits

import main
import model
from model import GridPhysics, MasterProblem, PortPhysics, PortSubProblem, cone_misocp, port_bounds, ray_support
from monitor import RunMonitor, comparison_metrics
from Network.jiangkou import BUDGET, LOAD_NODES, Jiangkou
from region import build_region, covered, polygon_union
from vertify import (ac_interval_possible, ac_scan_line, budget_schemes, reference_box, scan_ac_reference, scan_path,
                     signed_ac_witness)

OUTPUT = Path(__file__).resolve().parent/'results'/'jiangkou_test'
BUDGETS = (BUDGET, np.inf)          # 万元：二维测试预算与无限预算
SIGNS = list(product((1, -1), repeat=2))
SECONDS = 1800.                     # 每次构域的总时限（四个分区并行）
WORKERS = 16
RAYS = 19                           # 参考射线：每象限 0°–90° 等分
RAY_SECONDS = 300.                  # 单条参考射线 MISOCP 的时限
RAY_GAP = 1e-4                      # 参考射线 MISOCP 的相对间隙
DIVISIONS = main.SCAN_DIVISIONS[2]  # AC 扫描每轴格数（同主线二维）
SCAN = OUTPUT/'scan'                # 扫描缓存（主线格式）
LADDER = OUTPUT/'ladder'            # 预算档：逐方案 AC 真值、各档 RCUT 与面积曲线
LADDER_DIVISIONS = 200              # 预算档公共格架的每轴格数（范围取无限预算的扫描框）
LADDER_STEP = .01                   # 选 RCUT 预算档：真值面积相对上一档累计增长 ≥ 最终面积的 1%
LADDER_SECONDS = 600.               # 每档 RCUT 的时限
CHUNK = 2048                        # 逐方案 AC 判定的每批点数（控制内存）
CATEGORY = {'obbt': 'OBBT', 'center': 'MISOCP', 'cone': 'MISOCP', 'lazy': 'MISOCP',
            'origin': '射线与SP', 'ray': '射线与SP', 'near': '射线与SP',
            'network': '割平面', 'sp': '割平面', 'cut': '割平面', 'network_end': '割平面', 'check': '夹逼判据'}


def region_settings(d):
    """与 main.run 相同的 RCUT 设置（主线默认值）。"""
    return SimpleNamespace(threads=main.SOLVER_THREADS, tau=main.REGION_TAU, discovery_eps=main.DISCOVERY_EPS,
                           discovery_share=main.DISCOVERY_SHARE, mip_seconds=main.MIP_SECONDS, mip_gap=main.MIP_GAP,
                           min_width=main.MIN_WIDTH[d], max_cones=main.MAX_CONES[d], threshold=main.CUT_THRESHOLD,
                           patience=main.CUT_PATIENCE, point_tol=main.POINT_TOL)


@cache
def network(name, ports=LOAD_NODES):
    return Jiangkou(ports, reduced=name == 'reduced')


def plans(net):
    """现状方案 x0 与路径全换方案（可升级走廊取最大型号）。"""
    return (net.encode_plan(net.initial_plan),
            net.encode_plan({c.id: c.types[-1].id for c in net.corridors}))


def tag(name, budget):
    return f"{name}_b{budget:g}" if np.isfinite(budget) else f"{name}_binf"


def plain(value):
    """JSON：数组转列表，非有限数转 null。"""
    if isinstance(value, dict):
        return {str(k): plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [plain(v) for v in value]
    if isinstance(value, (np.floating, float)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    return value


def save(name, data):
    OUTPUT.mkdir(parents=True, exist_ok=True)
    (OUTPUT/name).write_text(json.dumps(plain(data), ensure_ascii=False, indent=1), encoding='utf-8')


def check():
    """1. 规模；2. 原点（端口为零）在现状方案下各分区的 SP 可行性；3. 各端口沿坐标轴的单点承载力（未紧化 SOCP 射线）。"""
    rows = []
    for name in ('full', 'reduced'):
        net = network(name)
        x0, top = plans(net)
        switchable = [c for c in net.corridors if c.switchable]
        print(f"{net.name}: {net.n} 节点 / {net.n_corridors} 走廊 / {net.n_types} 型号变量；可升级 {len(switchable)} 段、"
              f"{sum(len(c.types)-1 for c in switchable)} 个升级型号，全换 {net.cost@top:.2f} 万元")
        origin = [bool(PortSubProblem(PortPhysics(deepcopy(net), np.array(s))).solve(x0, np.zeros(2), score_only=True)['feasible'])
                  for s in SIGNS]
        print(f"  原点（现状方案）SP 可行：{dict(zip(map(str, SIGNS), origin))}")
        for k, port in enumerate(net.load_nodes):
            for s in (1, -1):
                sign = np.ones(2, int)
                sign[k] = s
                equations = PortPhysics(deepcopy(net), sign)
                target = np.zeros(2)
                target[k] = port_bounds(equations.network)[k]
                reach = []
                for x in (x0, top):
                    # 测试记录：主线射线（Clarabel）数值失败时记为 None，构域中同样的失败会使该方向没有顶点
                    try:
                        reach.append(ray_support(equations, BUDGET, x, np.zeros(2), target, threads=4, time_limit=60.)['p'][k])
                    except RuntimeError as error:
                        reach.append(None)
                        print(f"  {port} sign {s}: {str(error)[:60]}")
                rows.append(dict(network=name, port=port, sign=s, existing=reach[0], upgraded=reach[1], box=target[k]))
                print(f"  {port} {'负荷' if s > 0 else '光伏'}：现状 {reach[0]} kW → 路径全换 {reach[1]} kW"
                      f"（分区盒 {target[k]:.1f} kW）")
    save('check.json', rows)


def obbt(threads):
    """OBBT 计时：每个配置新建 PortPhysics，单独计时一次紧化；另测 ++ 分区中心射线 MISOCP（行生成，紧化计时单列）。"""
    counter = dict(solves=0, open=0)
    original = model.obbt_extremes

    def counted(solver, indices, workers):
        ends = original(solver, indices, workers)
        counter['solves'] += 2*len(indices)
        counter['open'] += sum(e is None for pair in ends for e in pair)
        return ends
    model.obbt_extremes = counted
    rows = []
    for name in ('full', 'reduced'):
        net = network(name)
        x0, top = plans(net)
        configurations = [('现状', x0, (1, 1), t) for t in threads]
        configurations += [('现状', x0, (-1, -1), 4), ('全换', top, (1, 1), 4), ('全换', top, (-1, -1), 4)]
        for plan, x, sign, workers in configurations:
            equations = PortPhysics(deepcopy(net), np.array(sign))
            counter.update(solves=0, open=0)
            start = time.perf_counter()
            equations.obbt(x, port_bounds(equations.network), workers)
            seconds = time.perf_counter()-start
            rows.append(dict(kind='obbt', network=name, plan=plan, sign=sign, threads=workers, seconds=seconds,
                             quantities=len(equations.boxes[tuple(x)]), **counter))
            print(f"{net.name} OBBT {plan} {sign} {workers} 线程：{seconds:7.1f} s，{counter['solves']} 次 SOCP"
                  f"（单次 {seconds*workers/counter['solves']*1e3:.1f} ms·线程），未证得最优端点 {counter['open']}")
        # 中心射线 MISOCP：x 自由、ξ1=ξ2，新现任网架先 OBBT（4 线程）再以惰性约束加入
        equations = PortPhysics(deepcopy(net), np.array((1, 1)))
        bounds, spent = port_bounds(equations.network), []

        def tighten(x, point, bound):
            begin = time.perf_counter()
            equations.obbt(np.asarray(x), bounds, 4)
            spent.append(time.perf_counter()-begin)
        start = time.perf_counter()
        answer = cone_misocp(equations, BUDGET, bounds, np.ones(2), tighten=tighten, time_limit=1800., mip_gap=0., threads=4)
        seconds = time.perf_counter()-start
        rows.append(dict(kind='center', network=name, seconds=seconds, obbt_seconds=sum(spent), schemes=len(spent),
                         bound=answer['bound'], point=answer['point'], status=answer['status']))
        print(f"{net.name} 中心射线 MISOCP（++，预算 {BUDGET:g} 万元）：共 {seconds:.1f} s，其中 {len(spent)} 个网架的 OBBT "
              f"{sum(spent):.1f} s；上界 ξ={answer['bound']:.4f}，现任点 {np.round(answer['point'], 1)} kW")
    model.obbt_extremes = original
    save('obbt.json', rows)


def _ray(name, budget, sign, theta):
    """x 自由的射线 MISOCP：max t  s.t. 端口幅值 u=t·(cosθ, sinθ)，预算内存在可行网架（未紧化 SOCP），分区盒内。"""
    equations = PortPhysics(deepcopy(network(name)), np.array(sign))
    direction = np.array([np.cos(theta), np.sin(theta)])
    problem = MasterProblem(equations, budget=budget, threads=1)
    with problem.model as solver:
        reach = solver.addVar(name='reach')
        for value, power in zip(direction, problem.loads.values()):
            solver.addConstr(power == value*reach)
        problem.power.UB = port_bounds(equations.network)
        solver.setObjective(reach, GRB.MAXIMIZE)
        solver.Params.MIPGap, solver.Params.TimeLimit = RAY_GAP, RAY_SECONDS
        start = time.perf_counter()
        solver.optimize()
        return dict(network=name, budget=budget, sign=sign, theta=theta, seconds=time.perf_counter()-start,
                    status=solver.Status, reach=solver.ObjVal if solver.SolCount else None, bound=solver.ObjBound,
                    x=np.rint(problem.x.X).astype(int) if solver.SolCount else None)


def rays():
    """参考射线并行求解（spawn 进程，每进程 1 线程）。"""
    angles = np.linspace(0., np.pi/2, RAYS)
    jobs = [(name, budget, sign, theta) for name in ('full', 'reduced') for budget in BUDGETS
            for sign in SIGNS for theta in angles]
    with ProcessPoolExecutor(WORKERS, mp_context=multiprocessing.get_context('spawn')) as pool:
        rows = list(pool.map(_ray, *zip(*jobs)))
    for name, budget in product(('full', 'reduced'), BUDGETS):
        part = [r for r in rows if r['network'] == name and r['budget'] == budget]
        print(f"{tag(name, budget)}：{len(part)} 条射线，单条 {np.median([r['seconds'] for r in part]):.1f} s（中位）"
              f" / {max(r['seconds'] for r in part):.1f} s（最长），未到最优 {sum(r['status'] != GRB.OPTIMAL for r in part)}")
    save('rays.json', [{k: v for k, v in r.items() if k != 'x'} | dict(cost=None if r['x'] is None else
                        float(network(r['network']).cost@r['x'])) for r in rows])


def breakdown(history):
    """各分区耗时拆分：相邻两帧之间的时间记到后一帧的类别（step 帧按 kind，OBBT/MISOCP/射线与SP/割平面/判据/其他）。"""
    partition, last, seconds, counts = None, {}, {}, {}
    for frame in history:
        patch = frame['patch']
        partition = patch.get('partition', partition)   # 主进程的 start / region_end 不属于任何分区
        if partition is None:
            continue
        if partition in last:
            kind = patch['step']['kind'] if patch['event'] == 'step' and 'step' in patch else patch['event']
            category = CATEGORY.get(kind, '其他')
            seconds.setdefault(partition, {}).setdefault(category, 0.)
            seconds[partition][category] += frame['elapsed']-last[partition]
            counts.setdefault(partition, {}).setdefault(kind, 0)
            counts[partition][kind] += 1
        last[partition] = frame['elapsed']
    return seconds, counts


def region(names, budgets, seconds):
    """二维 RCUT：四个分区并行（WORKERS 进程），各存回放与结果（内外域、分区摘要、耗时拆分）。"""
    for name in names:
        for budget in budgets:
            label = tag(name, budget)
            monitor = RunMonitor(output=OUTPUT/f'{label}.json.gz', algorithm='RCUT · 径向夹逼 + 割平面')
            start = time.perf_counter()
            # 测试记录：某分区抛错时 build_region 取消全部分区，保存出错前的回放与错误信息
            try:
                result = build_region(network(name), budget=budget, monitor=monitor, seconds=seconds, workers=WORKERS,
                                      settings=region_settings(2))
                result = {k: result[k] for k in ('status', 'certified', 'partitions', 'inner', 'outer', 'timing')}
            except Exception as error:
                result = dict(status='error', certified=False, error=repr(error))
            wall = time.perf_counter()-start
            monitor.save()
            split, counts = breakdown(monitor.history)
            save(f'{label}_result.json', dict(network=name, budget=budget, wall=wall, seconds=seconds, breakdown=split,
                                              counts=counts, **result))
            print(f"{label}：{result['status']}，墙钟 {wall:.0f} s；耗时拆分 " +
                  "；".join(f"{p} " + "、".join(f"{k} {v:.0f}s" for k, v in sorted(s.items(), key=lambda kv: -kv[1]))
                           for p, s in split.items()), flush=True)


@cache
def schemes(name, budget):
    """按段枚举的预算内全部网架（主线 budget_schemes，段内走廊同取一个选择）。"""
    return budget_schemes(GridPhysics(deepcopy(network(name)), 'socp'), budget, threads=1)


def _ac_truth(indices, coordinates, budget):
    """全网 AC 真值的一批格点：主线 ac_scan_line（逐网架 AC 潮流见证、必要区间排除、未决点求完整 AC 等式）。"""
    return ac_scan_line((indices, coordinates), network=network('full'), budget=budget, schemes=schemes('full', budget), mode=1)


def level(budget):
    return f'b{budget:g}' if np.isfinite(budget) else 'binf'


def acscan(budgets):
    """与真值对比：1. 格架为等值网 SOCP 全局界与各 RCUT 外域范围之并（同 main.run）；2. 等值网的配对 AC/SOCP 参考
    （主线 scan_ac_reference）；3. 同一格心上全网的 AC 真值；4. RCUT 内域 / 外域与各参考的 MR/FR，存 NPZ 与 JSON。"""
    for budget in budgets:
        runs = {name: json.loads((OUTPUT/f'{tag(name, budget)}_result.json').read_text(encoding='utf-8'))
                for name in ('reduced', 'full') if (OUTPUT/f'{tag(name, budget)}_result.json').exists()}
        runs = {name: run for name, run in runs.items() if run['status'] != 'error'}
        # 1. 格架
        lower, upper = reference_box(network('reduced'), budget, mode=1, output=SCAN)
        for run in runs.values():
            vertices = np.vstack([row['vertices'] for row in run['outer']])
            lower, upper = np.minimum(lower, vertices.min(axis=0)), np.maximum(upper, vertices.max(axis=0))
        # 2. 等值网的配对 AC/SOCP 参考
        start = time.perf_counter()
        paired = scan_ac_reference(network('reduced'), budget, dict(axis_lower=lower, bounds=upper, shape=(DIVISIONS,)*2),
                                   scan_path(network('reduced'), budget, SCAN, 1), workers=20, mode=1)
        paired_seconds = time.perf_counter()-start
        states = np.asarray(paired['states'])
        lower, upper = np.asarray(paired['axis_lower'], float), np.asarray(paired['bounds'], float)
        cells = np.indices(states.shape).reshape(2, -1).T
        points = lower+(cells+.5)*(upper-lower)/np.array(states.shape)
        # 3. 全网 AC 真值（同一格心）
        start = time.perf_counter()
        truth = np.zeros(len(points), np.int8)
        blocks = [np.arange(i, min(i+128, len(points))) for i in range(0, len(points), 128)]
        with ProcessPoolExecutor(20, mp_context=multiprocessing.get_context('spawn')) as pool:
            for line in pool.map(_ac_truth, blocks, [points[b] for b in blocks], [budget]*len(blocks)):
                truth[line['indices']] = line['states']
        truth_seconds = time.perf_counter()-start
        # 4. 对比：计算域（RCUT 内域 / 外域）对各参考，以及参考之间；未决格不计入
        references = dict(ac_full=truth, ac_reduced=states.ravel(), socp_reduced=np.asarray(paired['socp_states']).ravel())
        masks = {f'{name}_{key}': covered(points, run[key]) for name, run in runs.items() for key in ('inner', 'outer')}
        metrics = {}
        for mask, inside in masks.items():
            for reference, labels in references.items():
                known = labels != 0
                metrics[f'{mask}~{reference}'] = comparison_metrics(inside[known], labels[known] == 1)
        for a, b in (('socp_reduced', 'ac_reduced'), ('ac_reduced', 'ac_full'), ('socp_reduced', 'ac_full')):
            known = (references[a] != 0) & (references[b] != 0)
            metrics[f'{a}~{b}'] = comparison_metrics((references[a] == 1)[known], (references[b] == 1)[known])
        undecided = {reference: int((labels == 0).sum()) for reference, labels in references.items()}
        np.savez_compressed(OUTPUT/f'acscan_{level(budget)}.npz', points=points, shape=states.shape, lower=lower,
                            upper=upper, **references, **masks)
        save(f'acscan_{level(budget)}.json', dict(budget=budget, shape=states.shape, lower=lower, upper=upper,
             schemes=len(schemes('full', budget)), paired_seconds=paired_seconds, truth_seconds=truth_seconds,
             undecided=undecided, metrics=metrics))
        print(f"{level(budget)}：格架 {states.shape}，等值网配对扫描 {paired_seconds:.0f} s，全网 AC 真值 {truth_seconds:.0f} s"
              f"（{len(schemes('full', budget))} 个网架），未决 {undecided}", flush=True)
        for key, value in metrics.items():
            print(f"  {key:32s} MR {value['mr_percent']}  FR {value['fr_percent']}  "
                  f"(漏 {value['missed_cells']} / 多 {value['extra_cells']} / 参考 {value['reference_cells']} / 计算 {value['computed_cells']})")
    acplot(budgets)


def acplot(budgets):
    """逐格对比图：每档预算一幅，全网 AC 真值可行格浅色，RCUT 内域（等值网）遗漏格红、多余格橙，内域边界描线。"""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    plt.rcParams['font.sans-serif'], plt.rcParams['axes.unicode_minus'] = ['Microsoft YaHei', 'SimHei'], False
    figure, axes = plt.subplots(1, len(budgets), figsize=(7.5*len(budgets), 7.2), squeeze=False)
    for ax, budget in zip(axes[0], budgets):
        data = np.load(OUTPUT/f'acscan_{level(budget)}.npz')
        report = json.loads((OUTPUT/f'acscan_{level(budget)}.json').read_text(encoding='utf-8'))
        shape, lower, upper = tuple(data['shape']), data['lower'], data['upper']
        truth, inside = data['ac_full'].reshape(shape), data['reduced_inner'].reshape(shape)
        # 0 空白、1 真值可行且在内域、2 遗漏（真值可行、内域外）、3 多余（内域内、真值不可行）、4 未决
        code = np.where(truth == 0, 4, np.where(truth == 1, np.where(inside, 1, 2), np.where(inside, 3, 0)))
        colors = ListedColormap(['#ffffff', '#cfe0f6', '#d93a3a', '#f0a030', '#c9c8bf'])
        ax.imshow(code.T, origin='lower', cmap=colors, vmin=-.5, vmax=4.5, interpolation='nearest', aspect='equal',
                  extent=(lower[0], upper[0], lower[1], upper[1]))
        run = json.loads((OUTPUT/f'{tag("reduced", budget)}_result.json').read_text(encoding='utf-8'))
        shape_inner = polygon_union([np.array(r['vertices']) for r in run['inner']])
        for g in getattr(shape_inner, 'geoms', [shape_inner]):
            ax.plot(*g.exterior.xy, color='#2a78d6', lw=1.2)
        m = report['metrics']['reduced_inner~ac_full']
        ax.set_title(f"预算 {'无限' if not np.isfinite(budget) else f'{budget:g} 万元'}：RCUT 内域（等值网）vs 全网 AC 真值\n"
                     f"MR {m['mr_percent']:.3f}%（遗漏 {m['missed_cells']} 格）  FR {m['fr_percent']:.3f}%（多余 {m['extra_cells']} 格）",
                     fontsize=10.5, loc='left')
        ax.set_xlabel(f'{LOAD_NODES[0]}（左下）注入 kW，正=负荷、负=光伏', fontsize=9, color='#6b6a63')
        ax.set_ylabel(f'{LOAD_NODES[1]}（右上）注入 kW', fontsize=9, color='#6b6a63')
        ax.axhline(0, color='#9a9890', lw=.6)
        ax.axvline(0, color='#9a9890', lw=.6)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in ('#cfe0f6', '#d93a3a', '#f0a030', '#c9c8bf')]
    figure.legend(handles+[plt.Line2D([], [], color='#2a78d6', lw=1.2)],
                  ['真值可行且在内域', '遗漏：真值可行、内域外', '多余：内域内、真值不可行', '真值未决', 'RCUT 内域边界'],
                  loc='lower center', ncol=5, frameon=False, fontsize=9)
    figure.tight_layout(rect=(0, .05, 1, 1))
    figure.savefig(OUTPUT/'acscan.png', dpi=150)
    print(OUTPUT/'acscan.png')


def section_names(net):
    """段的简称：左下端口路径 L、右上端口路径 R，自根向端口编号（L1 为左下首段）。"""
    by_child = {c.endpoints[1]: c for c in net.corridors}
    names = {}
    for letter, port in zip('LR', net.load_nodes):
        path, node = [], port
        while node != net.root:
            path.append(by_child[node].id)
            node = by_child[node].endpoints[0]
        path.reverse()
        for k, section in enumerate(sorted((s for s in net.sections if s[0] in path), key=lambda s: path.index(s[0]))):
            names[section] = f'{letter}{k+1}'
    return names


def describe(net, x):
    """方案的升级段：段名·新型号（70 为 4x70，240 为 2x4x240），无升级为“现状”。"""
    plan, corridors = net.decode_plan(x), {c.id: c for c in net.corridors}
    parts = [f"{name}·{plan[s[0]].split('x')[-1]}" for s, name in sorted(section_names(net).items(), key=lambda kv: kv[1])
             if plan[s[0]] != corridors[s[0]].existing_type]
    return ' '.join(parts) or '现状'


def _scheme_states(index, points):
    """全网一个方案在全部格心上的 AC 标签：潮流见证可行为 1，必要区间排除为 -1，其余未决为 0（主线判据，分批控制内存）。"""
    net, x = network('full'), schemes('full', np.inf)[index]
    states = np.zeros(len(points), np.int8)
    with threadpool_limits(limits=1):
        for start in range(0, len(points), CHUNK):
            block = points[start:start+CHUNK]
            witness = signed_ac_witness(net, x, block, mode=1)['feasible']
            label = witness.astype(np.int8)
            rest = np.flatnonzero(~witness)
            if len(rest):
                label[rest[~ac_interval_possible(net, x, block[rest], mode=1)]] = -1
            states[start:start+len(block)] = label
    return index, states


def ladder(seconds):
    """预算档：1. 公共格架上全网逐方案的 AC 标签（truth.npz，已有则复用）；2. 按费用逐个并入得真值面积阶梯（每个方案
    费用都是跳变点），自 0 起取相对上一档累计增长 ≥ LADDER_STEP·最终面积的跳变点及最高档为 RCUT 预算档；3. 每档在等值网上
    跑 RCUT（ladder.json 中已有的档复用），对该档真值求 MR/FR；4. 画图。"""
    LADDER.mkdir(parents=True, exist_ok=True)
    full = network('full')
    candidates = schemes('full', np.inf)
    costs = full.cost_offset+candidates@full.cost
    # 1. 真值：公共格架（无限预算的扫描框）上逐方案的 AC 标签
    box = json.loads((OUTPUT/'acscan_binf.json').read_text(encoding='utf-8'))
    lower, upper, shape = np.array(box['lower']), np.array(box['upper']), (LADDER_DIVISIONS,)*2
    cells = np.indices(shape).reshape(2, -1).T
    points = lower+(cells+.5)*(upper-lower)/np.array(shape)
    cell_area = float(np.prod((upper-lower)/np.array(shape)))
    if (LADDER/'truth.npz').exists():
        with np.load(LADDER/'truth.npz') as saved:
            states, truth_seconds = saved['states'], float(saved['seconds'])
    else:
        start = time.perf_counter()
        states = np.zeros((len(candidates), len(points)), np.int8)
        with ProcessPoolExecutor(20, mp_context=multiprocessing.get_context('spawn')) as pool:
            for index, row in pool.map(_scheme_states, range(len(candidates)), [points]*len(candidates)):
                states[index] = row
        truth_seconds = time.perf_counter()-start
    print(f"逐方案 AC 真值：{len(candidates)} 个方案 × {len(points)} 格，{truth_seconds:.0f} s；"
          f"单方案未决格合计 {int((states == 0).sum())}", flush=True)
    # 2. 真值阶梯：按费用并入，记面积（总与四个象限）、未决格与新方案的升级段
    quadrants = {label: np.all(np.sign(points) == np.array(sign), axis=1)
                 for label, sign in (('++', (1, 1)), ('+-', (1, -1)), ('-+', (-1, 1)), ('--', (-1, -1)))}
    order = np.argsort(costs)
    feasible, impossible, staircase = np.zeros(len(points), bool), np.ones(len(points), bool), []
    for k in order:
        feasible |= states[k] == 1
        impossible &= states[k] == -1
        staircase.append(dict(cost=float(costs[k]), scheme=describe(full, candidates[k]), area=feasible.sum()*cell_area,
                              quadrants={q: float((feasible & m).sum()*cell_area) for q, m in quadrants.items()},
                              undecided=int((~feasible & ~impossible).sum())))
    final, picks = staircase[-1]['area'], [0]
    for i in range(1, len(staircase)-1):
        if staircase[i]['area']-staircase[picks[-1]]['area'] >= LADDER_STEP*final:
            picks.append(i)
    picks.append(len(staircase)-1)
    print(f"真值阶梯：最终面积 {final:.0f} kW²，{len(picks)} 个 RCUT 预算档：" +
          "、".join(f"{staircase[i]['cost']:.2f}" for i in picks), flush=True)
    np.savez_compressed(LADDER/'truth.npz', states=states, costs=costs, points=points, lower=lower, upper=upper,
                        shape=shape, order=order, seconds=truth_seconds)
    # 3. 各档 RCUT（等值网）与该档真值的逐格对比；已跑过的档复用
    done = ({round(l['budget'], 6): l for l in json.loads((LADDER/'ladder.json').read_text(encoding='utf-8'))['levels']}
            if (LADDER/'ladder.json').exists() else {})
    levels = []
    for i in picks:
        budget = staircase[i]['cost']+1e-6
        if round(budget, 6) in done:
            levels.append(done[round(budget, 6)])
            continue
        labels = np.where(np.any(states[order[:i+1]] == 1, axis=0), 1,
                          np.where(np.all(states[order[:i+1]] == -1, axis=0), -1, 0))
        monitor = RunMonitor(output=LADDER/f'b{budget:.2f}.json.gz', algorithm='RCUT · 径向夹逼 + 割平面')
        begin = time.perf_counter()
        try:
            result = build_region(network('reduced'), budget=budget, monitor=monitor, seconds=seconds, workers=WORKERS,
                                  settings=region_settings(2))
        except Exception as error:
            result = dict(status='error', error=repr(error))
        wall = time.perf_counter()-begin
        monitor.save()
        level = dict(budget=budget, wall=wall, status=result['status'], truth=staircase[i])
        if result['status'] != 'error':
            known = labels != 0
            for key in ('inner', 'outer'):
                inside = covered(points, result[key])
                level[f'{key}_area'] = polygon_union([np.array(r['vertices']) for r in result[key]]).area
                level[f'{key}_metrics'] = comparison_metrics(inside[known], labels[known] == 1)
            level['partitions'] = [{k: p[k] for k in ('sign', 'status', 'seconds', 'gap', 'cones', 'networks')}
                                   for p in result['partitions']]
            level['inner'] = result['inner']
        else:
            level['error'] = result['error']
        levels.append(level)
        print(f"b={budget:6.2f}（{staircase[i]['scheme']}）：{result['status']}，{wall:.0f} s；真值 {staircase[i]['area']:.0f}，"
              f"内域 {level.get('inner_area', float('nan')):.0f}，外界 {level.get('outer_area', float('nan')):.0f} kW²；"
              f"内域 MR/FR {level.get('inner_metrics', {}).get('mr_percent')}/{level.get('inner_metrics', {}).get('fr_percent')}",
              flush=True)
    save('ladder/ladder.json', dict(cell_area=cell_area, lower=lower, upper=upper, shape=shape, truth_seconds=truth_seconds,
                                    schemes=len(candidates), staircase=staircase, levels=levels))
    ladder_plot()


def ladder_plot():
    """1. 面积—预算曲线：真值阶梯、RCUT 内域与外界；四个象限的真值面积；2. 若干预算档的可行域小图。"""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    plt.rcParams['font.sans-serif'], plt.rcParams['axes.unicode_minus'] = ['Microsoft YaHei', 'SimHei'], False
    blue, orange, ink, muted, grid = '#2a78d6', '#eb6834', '#1f1f1e', '#6b6a63', '#e4e3dc'
    data = json.loads((LADDER/'ladder.json').read_text(encoding='utf-8'))
    stairs, levels, scale = data['staircase'], [l for l in data['levels'] if l['status'] != 'error'], 1e4
    cells = lambda l, key: l[f'{key}_metrics']['computed_cells']*data['cell_area']/scale   # 与真值同一格架计面积
    figure, (left, right) = plt.subplots(1, 2, figsize=(17, 6.4), gridspec_kw=dict(width_ratios=[1.35, 1]))
    cost = [s['cost'] for s in stairs]+[max(s['cost'] for s in stairs)*1.02]
    area = [s['area']/scale for s in stairs]
    left.step(cost, area+area[-1:], where='post', color=ink, lw=1.4, label='全网 AC 真值（216 个方案逐个并入）')
    left.scatter([l['budget'] for l in levels], [cells(l, 'inner') for l in levels], s=36, color=blue, zorder=4,
                 label='RCUT 内域（等值网）')
    left.scatter([l['budget'] for l in levels], [cells(l, 'outer') for l in levels], s=44, facecolor='none',
                 edgecolor=orange, lw=1.3, zorder=4, label='RCUT 外界')
    jumps = sorted(range(1, len(stairs)), key=lambda i: stairs[i-1]['area']-stairs[i]['area'])[:6]
    for i in jumps:
        left.annotate(stairs[i]['scheme'], (stairs[i]['cost'], stairs[i]['area']/scale), xytext=(6, -14),
                      textcoords='offset points', fontsize=7.5, color=muted)
    left.set_xlabel('预算 万元（整段换线；L 左下 B000078 路径段、R 右上 B000042 路径段，·70=4x70、·240=2x4x240）',
                    fontsize=9, color=muted)
    left.set_ylabel('二维可行域面积 ×10⁴ kW²（四个分区合计）', fontsize=9, color=muted)
    left.set_title('可行域面积随预算增长', fontsize=11, color=ink, loc='left')
    left.legend(fontsize=8.5, frameon=False, loc='lower right')
    for q, color, label in (('++', blue, '++ 两端负荷'), ('--', orange, '−− 两端光伏'), ('+-', '#1baf7a', '+− 左下负荷、右上光伏'),
                            ('-+', '#e87ba4', '−+ 左下光伏、右上负荷')):
        right.step(cost, [s['quadrants'][q]/scale for s in stairs]+[stairs[-1]['quadrants'][q]/scale], where='post',
                   color=color, lw=1.6, label=label)
    right.set_xlabel('预算 万元', fontsize=9, color=muted)
    right.set_ylabel('真值面积 ×10⁴ kW²', fontsize=9, color=muted)
    right.set_title('按分区拆分（全网 AC 真值）', fontsize=11, color=ink, loc='left')
    right.legend(fontsize=8.5, frameon=False, loc='upper left')
    for ax in (left, right):
        ax.grid(color=grid, lw=.8)
        ax.tick_params(colors=muted, labelsize=8)
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)
    figure.tight_layout()
    figure.savefig(LADDER/'ladder.png', dpi=150)
    # 2. 可行域小图：在 RCUT 档中均匀取 8 档
    truth = np.load(LADDER/'truth.npz')
    states, order, shape = truth['states'], truth['order'], tuple(truth['shape'])
    lower, upper, costs = truth['lower'], truth['upper'], truth['costs']
    chosen = [levels[int(round(k))] for k in np.linspace(0, len(levels)-1, min(8, len(levels)))]
    figure, axes = plt.subplots(2, 4, figsize=(18, 10.6), gridspec_kw=dict(hspace=.42, wspace=.18))
    colors = ListedColormap(['#ffffff', '#cfe0f6', '#d93a3a', '#f0a030'])
    for ax, level in zip(axes.ravel(), chosen):
        affordable = order[costs[order] <= level['budget']]
        labels = np.where(np.any(states[affordable] == 1, axis=0), 1, 0)
        inside = covered(truth['points'], level['inner'])
        code = np.where(labels == 1, np.where(inside, 1, 2), np.where(inside, 3, 0)).reshape(shape)
        ax.imshow(code.T, origin='lower', cmap=colors, vmin=-.5, vmax=3.5, interpolation='nearest', aspect='equal',
                  extent=(lower[0], upper[0], lower[1], upper[1]))
        # 画图用：0.05 kW 的闭合运算合上相邻锥块之间的零宽缝（面积与逐格判定不受影响）
        shape_inner = polygon_union([np.array(r['vertices']) for r in level['inner']]).buffer(.05).buffer(-.05)
        for g in getattr(shape_inner, 'geoms', [shape_inner]):
            ax.plot(*g.exterior.xy, color=blue, lw=1)
        m = level['inner_metrics']
        ax.set_title(f"{level['budget']:.2f} 万元：{level['truth']['scheme']}\n真值 {level['truth']['area']/scale:.2f}、"
                     f"内域 {cells(level, 'inner'):.2f} ×10⁴ kW²；MR {m['mr_percent']:.2f}% FR {m['fr_percent']:.2f}%",
                     fontsize=8.5, loc='left')
        ax.axhline(0, color='#9a9890', lw=.5)
        ax.axvline(0, color='#9a9890', lw=.5)
        ax.tick_params(labelsize=7, colors=muted)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in ('#cfe0f6', '#d93a3a', '#f0a030')]
    figure.legend(handles+[plt.Line2D([], [], color=blue, lw=1)],
                  ['真值可行且在内域', '遗漏：真值可行、内域外', '多余：内域内、真值不可行', 'RCUT 内域边界'],
                  loc='lower center', ncol=4, frameon=False, fontsize=9)
    figure.supxlabel(f'{LOAD_NODES[0]}（左下）注入 kW，正=负荷、负=光伏；纵轴 {LOAD_NODES[1]}（右上）注入 kW',
                     fontsize=9, color=muted, y=.045)
    figure.subplots_adjust(left=.04, right=.99, top=.93, bottom=.11)
    figure.savefig(LADDER/'ladder_regions.png', dpi=140)
    print(LADDER/'ladder.png', LADDER/'ladder_regions.png')


def reach(polygon, theta, sign):
    """星形区域沿 sign·(cosθ, sinθ) 从原点出发的半径。"""
    direction = np.array(sign)*np.array([np.cos(theta), np.sin(theta)])
    part = polygon.intersection(LineString([(0., 0.), tuple(4000.*direction)]))
    points = np.array([c for g in getattr(part, 'geoms', [part]) for c in getattr(g, 'coords', [])])
    return float((points@direction).max()) if len(points) else 0.


def compare():
    """1. 读结果与参考射线；2. 面积、对称差、半径误差；3. 耗时表；4. 对比图。"""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif'], plt.rcParams['axes.unicode_minus'] = ['Microsoft YaHei', 'SimHei'], False
    references = json.loads((OUTPUT/'rays.json').read_text(encoding='utf-8'))
    runs, summary = {}, dict(runs={}, rays={})
    for name, budget in product(('full', 'reduced'), BUDGETS):
        path = OUTPUT/f'{tag(name, budget)}_result.json'
        if path.exists():
            run = json.loads(path.read_text(encoding='utf-8'))
            if run['status'] == 'error':
                summary['runs'][f'{tag(name, budget)}_error'] = dict(error=run['error'], wall=run['wall'], breakdown=run['breakdown'])
                continue
            runs[name, budget] = run
    # 2. 每档预算：全网与等值网的内域面积、对称差；参考射线的等值误差（两边都证得最优的方向，按是否含光伏分开）
    def rays_of(name, budget):
        return [r for r in references if r['network'] == name and (r['budget'] or np.inf) == budget]

    geometry = {}
    for budget in BUDGETS:
        geometry[budget] = {name: {key: polygon_union([np.array(r['vertices']) for r in runs[name, budget][key]])
                                   for key in ('inner', 'outer')} for name in ('full', 'reduced') if (name, budget) in runs}
        if len(geometry[budget]) == 2:
            inner = geometry[budget]['full']['inner'], geometry[budget]['reduced']['inner']
            summary['runs'][f'b{budget:g}'] = dict(
                inner_area=[g.area for g in inner],
                outer_area=[geometry[budget][n]['outer'].area for n in ('full', 'reduced')],
                symmetric_difference=inner[0].symmetric_difference(inner[1]).area/inner[0].area)
        pairs = [(f, g) for f, g in zip(rays_of('full', budget), rays_of('reduced', budget))
                 if f['status'] == GRB.OPTIMAL and g['status'] == GRB.OPTIMAL]
        entry = dict(both_optimal=len(pairs), directions=len(rays_of('full', budget)))
        for label, keep in (('load', lambda r: min(r['sign']) > 0), ('pv', lambda r: min(r['sign']) < 0)):
            error = np.array([g['reach']/f['reach']-1 for f, g in pairs if keep(f)])
            entry[f'{label}_error'] = dict(count=len(error), max=float(np.abs(error).max()), mean=float(np.abs(error).mean()))
        for name in ('full', 'reduced'):
            if (name, budget) in runs:
                rows = [r for r in rays_of(name, budget) if r['status'] == GRB.OPTIMAL]
                error = np.array([reach(geometry[budget][name]['inner'], r['theta'], r['sign'])/r['reach']-1 for r in rows])
                entry[f'{name}_rcut_inner_error'] = dict(max=float(np.abs(error).max()), mean=float(np.abs(error).mean()))
        summary['rays'][f'b{budget:g}'] = entry
    # 3. 耗时表：每次运行每分区的状态、用时、网架数与拆分
    table = []
    for (name, budget), run in runs.items():
        for p in run['partitions']:
            label = ''.join('+' if s > 0 else '-' for s in p['sign'])
            split = run['breakdown'].get(label, {})
            table.append(dict(run=tag(name, budget), partition=label, status=p['status'], seconds=p['seconds'],
                              gap=p['gap'], cones=p['cones'], networks=p['networks'], accepted=p['accepted'], cuts=p['cuts'],
                              obbt_count=run['counts'].get(label, {}).get('obbt', 0),
                              **{f'{k}_s': v for k, v in split.items()}))
    summary['table'] = table
    save('summary.json', summary)
    for row in table:
        print(json.dumps(plain(row), ensure_ascii=False))
    print(json.dumps(plain({k: v for k, v in summary.items() if k != 'table'}), ensure_ascii=False, indent=1))
    # 4. 对比图：两档预算各一幅（参考射线给出的边界：全网实线、等值网虚线，未证得最优的方向空心灰点；有 RCUT 结果时叠加内域），
    #    右侧为计时：一次 OBBT（现状方案 ++）与中心射线 MISOCP、参考射线中位耗时，全网与等值网对照
    figure, axes = plt.subplots(1, 3, figsize=(19, 6.4), gridspec_kw=dict(width_ratios=[1, 1, .9]))
    blue, orange, ink, muted, grey = '#2a78d6', '#eb6834', '#1f1f1e', '#6b6a63', '#b5b3aa'
    for ax, budget in zip(axes[:2], BUDGETS):
        for name, color, style, marker, label in (('full', blue, '-', 'o', '全网'), ('reduced', orange, '--', 'x', '等值网')):
            if (name, budget) in runs:
                shape = geometry[budget][name]['inner']
                for g in getattr(shape, 'geoms', [shape]):
                    ax.fill(*g.exterior.xy, color=color, alpha=.12, lw=0)
            for sign in SIGNS:
                rows = sorted((r for r in rays_of(name, budget) if tuple(r['sign']) == sign and r['reach'] is not None),
                              key=lambda r: r['theta'])
                points = np.array([np.array(sign)*r['reach']*np.array([np.cos(r['theta']), np.sin(r['theta'])]) for r in rows])
                optimal = np.array([r['status'] == GRB.OPTIMAL for r in rows])
                ax.plot(points[optimal, 0], points[optimal, 1], color=color, lw=1.6, ls=style, zorder=3,
                        label=f'{label}（参考射线 MISOCP，连证得最优的方向）' if sign == SIGNS[0] else None)
                ax.scatter(*points[optimal].T, s=18, marker=marker, color=color, zorder=4)
                ax.scatter(*points[~optimal].T, s=36, marker='o', facecolor='none', edgecolor=grey, zorder=5,
                           label='未证得最优（时限或数值）' if sign == SIGNS[0] and name == 'full' else None)
        ax.axhline(0, color='#d9d8d0', lw=1, zorder=0)
        ax.axvline(0, color='#d9d8d0', lw=1, zorder=0)
        ax.set_aspect('equal')
        ax.set_title(f"预算 {'无限' if not np.isfinite(budget) else f'{budget:g} 万元'}：二维规划可行域（未紧化 SOCP）",
                     fontsize=11, color=ink, loc='left')
        ax.set_xlabel(f'{LOAD_NODES[0]}（左下）注入 kW，正=负荷、负=光伏', fontsize=9, color=muted)
        ax.set_ylabel(f'{LOAD_NODES[1]}（右上）注入 kW', fontsize=9, color=muted)
        ax.legend(fontsize=8, frameon=False, loc='upper center', bbox_to_anchor=(.5, -.1), ncol=1)
        ax.tick_params(colors=muted, labelsize=8)
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)
    ax = axes[2]
    timing = json.loads((OUTPUT/'obbt.json').read_text(encoding='utf-8'))

    def seconds(network, **match):
        return next(r['seconds'] for r in timing if r['network'] == network
                    and all(tuple(r.get(k)) == v if isinstance(v, tuple) else r.get(k) == v for k, v in match.items()))
    items = [('一次 OBBT，1 线程', dict(kind='obbt', plan='现状', sign=(1, 1), threads=1)),
             ('一次 OBBT，4 线程', dict(kind='obbt', plan='现状', sign=(1, 1), threads=4)),
             ('一次 OBBT，16 线程', dict(kind='obbt', plan='现状', sign=(1, 1), threads=16)),
             ('中心射线 MISOCP（13 万元）', dict(kind='center'))]
    values = {name: [seconds(name, **match) for _, match in items] for name in ('full', 'reduced')}
    for name in ('full', 'reduced'):
        values[name].append(float(np.median([r['seconds'] for r in rays_of(name, BUDGET)])))
    labels = [label for label, _ in items]+['参考射线 MISOCP 中位（13 万元）']
    position = np.arange(len(labels))
    for offset, name, color, label in ((-.2, 'full', blue, '全网 342 段'), (.2, 'reduced', orange, '等值网 50 段')):
        ax.barh(position+offset, values[name], height=.38, color=color, label=label)
        for y, value in zip(position+offset, values[name]):
            ax.text(value*1.08, y, f'{value:.1f} s', va='center', fontsize=8, color=ink)
    ax.set_xscale('log')
    ax.set_yticks(position, labels, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel('耗时 s（对数轴）', fontsize=9, color=muted)
    ax.set_title('计时：全网 vs 等值网', fontsize=11, color=ink, loc='left')
    ax.legend(fontsize=8, frameon=False, loc='lower right')
    ax.tick_params(colors=muted, labelsize=8)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    figure.tight_layout()
    figure.savefig(OUTPUT/'comparison.png', dpi=150)
    print(OUTPUT/'comparison.png')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='江口测试：网架自检、OBBT 计时、二维全网与等值网对比')
    parser.add_argument('step', choices=('check', 'obbt', 'rays', 'region', 'compare', 'acscan', 'ladder', 'ladder_plot'))
    parser.add_argument('--threads', type=int, nargs='+', default=[1, 4, 16], help='obbt：现状方案 ++ 分区的线程数')
    parser.add_argument('--networks', nargs='+', default=['reduced', 'full'], choices=('full', 'reduced'))
    parser.add_argument('--budgets', type=float, nargs='+', default=list(BUDGETS))
    parser.add_argument('--seconds', type=float, default=SECONDS)
    args = parser.parse_args()
    if args.step == 'check':
        check()
    elif args.step == 'obbt':
        obbt(args.threads)
    elif args.step == 'rays':
        rays()
    elif args.step == 'region':
        region(args.networks, args.budgets, args.seconds)
    elif args.step == 'acscan':
        acscan(args.budgets)
    elif args.step == 'ladder':
        ladder(LADDER_SECONDS)
    elif args.step == 'ladder_plot':
        ladder_plot()
    else:
        compare()
