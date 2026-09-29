"""二维 FourBus：已有网架充分覆盖后，才调用全局 physical 搜索。

运行比较：python experiments/fourbus_known_first.py
回放新法：python experiments/fourbus_known_first.py --replay known_first
主线、物理约束、SP 与监视窗口均复用原实现；本文件只实验新的调度。
"""
import argparse
from pathlib import Path
import sys
from time import perf_counter

import numpy as np
from shapely.geometry import MultiPoint
from shapely.ops import unary_union
from threadpoolctl import threadpool_limits

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from continuous import build_continuous_region, REGION_TAU, REFINEMENT_CHECKS
from model import GridPhysics, MasterProblem, RemainingRegionModel, SubProblem, SP_TIME_LIMIT, RESIDUAL_TIME_LIMIT
from monitor import RunMonitor
from Network.four_bus_five_corridor import FourBus
from region import RegionState, GEOMETRY_TOL, clip_polytope, contains, halfspaces, polytope_volume
from vertify import validate_socp_region


OUTPUT = Path(__file__).resolve().parents[1]/'results'/'fourbus_known_first'


def remaining_cells(vertices, inner_halfspaces):
    """凸候选域减去多个凸内域；逐面分割保留内部孔洞与低维剩余。

    对一个内域，先留下第一个违反面的外侧，再在该面内侧继续检查下一面。
    返回块的内部互不重叠；不能把不同网架的认证点混合取凸包。
    """
    cells = [np.asarray(vertices)] if len(vertices) else []
    for equations in inner_halfspaces:
        remaining = []
        for cell in cells:
            if contains(cell, equations).all():
                continue
            inside = cell
            for face in equations:
                if not len(inside):
                    break
                values = inside@face[:-1]+face[-1]
                if values.max() > GEOMETRY_TOL:
                    outside = clip_polytope(inside, face[-1]-GEOMETRY_TOL, face[:-1])
                    remaining.append(outside)
                inside = clip_polytope(inside, GEOMETRY_TOL-face[-1], -face[:-1])
        cells = remaining
    return cells


def known_candidates(region):
    """优先取未覆盖外域顶点；顶点查完后，从完整差集的凸块中选内部点。"""
    candidates = []
    for row in region.records.values():
        targets = (1-region.tau)*row['outer']
        owners = region.covering_schemes(targets, preferred=row['x'])
        candidates.extend((row['x'], point) for point, owner in zip(targets, owners) if owner is None)
    if candidates:
        return candidates, False

    # 所有顶点都已覆盖也可能留下孔洞；空差集才允许进入全局搜索。
    equations = region.inner_halfspaces()
    for row in region.records.values():
        cells = remaining_cells((1-region.tau)*row['outer'], equations)
        if cells:
            cell = max(cells, key=lambda cell: (polytope_volume(cell), float(cell.mean(axis=0)@region.bounds)))
            point = cell.mean(axis=0)
            assert region.covering_schemes([point])[0] is None, '剩余块选点必须在认证并集外'
            candidates.append((row['x'], point))
    return candidates, True


def build_known_first_region(network, method, budget, bounds, *, tau=REGION_TAU,
                             time_limit=300., threads=4, progress=None, cuts=(), clock=perf_counter):
    """1 初始化 → 2 已有网架选点 → 3 SP 更新 → 4 全局增补 → 5 全局证书。"""
    bounds = np.asarray(bounds, dtype=float)
    assert method == 'socp' and bounds.shape == (2,)
    equations = GridPhysics(network, method)
    region = RegionState(bounds, network.power_limit, tau, cuts)
    oracle = SubProblem(equations, threads=threads)
    monitor = RunMonitor.follow(progress, clock=clock)
    monitor.begin(network, method, budget, region, time_limit)
    cache = {}
    initial_cuts = len(region.cuts)
    maximum_answer = coverage = None
    sp_since_global = global_search = known_checks = hole_checks = 0

    # 1. 初始化与主线一致：两个轴向 MP2 + 总负荷 MP2，均含完整电气约束。
    for index, direction in enumerate([*np.eye(2), np.ones(2)]):
        monitor.initializing(index)
        problem = MasterProblem(equations, budget=budget, direction=direction,
                                cuts=region.cuts, threads=threads)
        with problem.model:
            answer = problem.solve(time_limit=monitor.remaining(time_limit))
        if answer is None:
            region.tighten_bounds(np.zeros(2), 0.)
            break
        axis_bounds = region.axis_bounds.copy()
        if index < 2:
            axis_bounds[index] = min(axis_bounds[index], answer['bound'])
            region.tighten_bounds(axis_bounds, region.total_bound)
        else:
            maximum_answer = answer
            region.tighten_bounds(axis_bounds, answer['bound'])
        x = answer['x']
        region.add_scheme(x, network.decode_plan(x), network.cost_offset+network.cost@x)
        region.add_point(x, answer['p']/bounds)
        monitor.seed(answer, region)

    while maximum_answer is not None:
        # 2. 已有网架优先：先查外域顶点，再检查完整差集；SP 次数不触发换网架。
        candidates, interior = known_candidates(region)
        known_checks += int(interior)
        if candidates:
            if interior:
                # 内部块均值只作未覆盖见证；同网架支撑比反复认证中心点扩大凸包更直接。
                support_candidates = []
                for x, point in candidates:
                    support = [q for q in region.witness_support(x, point)
                               if not contains([q], region.inner_equations(x))[0]]
                    point = max(support, key=lambda q: (float(q@bounds), tuple(q)))
                    support_candidates.append((x, point))
                candidates = support_candidates
            # 3. 固定 x,p 求原 SP；先收入本批全部认证点，再选最大 eta 的有效割。
            monitor.selecting(len(candidates), interior)
            scored = []
            for x, point in candidates:
                power = point*bounds
                key = tuple(x), tuple(power)
                if key not in cache:
                    monitor.sp_start(x, power, oracle.calls+1)
                    cache[key] = oracle.solve(x, power, time_limit=monitor.remaining(SP_TIME_LIMIT[method]))
                    sp_since_global += 1
                    hole_checks += int(interior)
                    monitor.sp_end(cache[key])
                checked = cache[key]
                scored.append((x, point, checked))
                if checked['feasible']:
                    region.add_point(x, point)
                    monitor.updated(region, 'feasible', x=x, power=power, checked=checked)
            if interior:
                pending = [row for row in scored if not row[2]['feasible']]
            else:
                owners = region.covering_schemes([point for _, point, _ in scored])
                pending = [row for row, owner in zip(scored, owners) if not row[2]['feasible'] and owner is None]
            if pending:
                x, point, checked = max(pending, key=lambda row: (row[2]['eta'], float(row[1]@bounds),
                                                                  tuple(row[0]), tuple(row[1])))
                region.apply_cut(checked['cut'])
                monitor.updated(region, 'cut', x=x, power=point*bounds, checked=checked)
            continue

        # 4. 只有所有已知收缩外域都已被并集覆盖，才允许全局 x,p,y 自由搜索。
        monitor.global_start(sp_since_global)
        problem = RemainingRegionModel(equations, budget, bounds, region.total_bound,
            region.cuts, region.inner_halfspaces(), region.tau, axis_bounds=region.axis_bounds,
            mode='physical', threads=threads)
        with problem.model:
            witness = problem.solve(GEOMETRY_TOL, time_limit=monitor.remaining(RESIDUAL_TIME_LIMIT))
        global_search += 1
        coverage = witness['bound']
        if not witness['complete']:
            x = witness['x']
            assert tuple(x) not in region.records, '局部覆盖完成后，全局未覆盖见证不能仍属已知网架'
            region.add_scheme(x, network.decode_plan(x), network.cost_offset+network.cost@x)
            assert witness['feasible']
            region.add_point(x, witness['p']/bounds)
        monitor.global_end(witness, region)
        if witness['complete']:
            break
        sp_since_global = 0

    # 5. 局部完成不等于全局完成；结束仍需原模型的全局覆盖上界证书。
    result = dict(method=method, budget=float(budget), residual_mode='physical', status='certified',
        tau=region.tau, axis_bounds=region.axis_bounds.copy(), coverage_bound=coverage,
        max_total=0. if maximum_answer is None else float(maximum_answer['p'].sum()),
        max_total_bound=region.total_bound,
        counts=dict(sp=oracle.calls, cuts=len(region.cuts)-initial_cuts, global_search=global_search,
                    known_checks=known_checks, hole_checks=hole_checks, schemes=len(region.records)),
        timing=monitor.timing(), **region.finish(True))
    monitor.finish(result, region)
    return result


def run_comparison(output=OUTPUT, *, repeats=3, divisions=80, threads=4):
    """同精度交替重复；独立扫描共用一次，每法只保存一份原生中位耗时回放。"""
    output = Path(output)
    network = FourBus(load_nodes=(1, 2))
    bounds = np.full(2, network.power_limit)
    runs = dict(baseline=[], known_first=[])
    print(f'FourBus SOCP, budget=20000, tau={REGION_TAU}, physical, threads={threads}, '
          f'baseline_interval={REFINEMENT_CHECKS}', flush=True)
    with threadpool_limits(limits=1):
        for repeat in range(repeats):
            order = ('baseline', 'known_first') if repeat % 2 == 0 else ('known_first', 'baseline')
            for name in order:
                monitor = RunMonitor(output=output/name/'monitor.json.gz')
                build = build_continuous_region if name == 'baseline' else build_known_first_region
                options = dict(residual_mode='physical') if name == 'baseline' else {}
                print(f'{name}, run {repeat+1}/{repeats}: start', flush=True)
                result = monitor.execute(lambda: build(network, 'socp', 20000., bounds,
                    tau=REGION_TAU, time_limit=300., threads=threads, progress=monitor, **options), show_ui=False)
                runs[name].append((result, monitor))
                print(f'{name}, run {repeat+1}: {result["timing"]["total_seconds"]:.3f}s, '
                      f'SP={result["counts"]["sp"]}, G={result["counts"]["global_search"]}, '
                      f'schemes={len(monitor.state["schemes"])}', flush=True)

        # 事后校验不参加计时、不限制网架，不将构域割或内域交给参考求解器。
        axis_bounds = runs['baseline'][0][0]['axis_bounds']
        print(f'Independent SOCP scan: {divisions} x {divisions}', flush=True)
        reference = validate_socp_region(network, 20000., divisions, axis_bounds, threads=threads,
            progress=lambda done, total: print(f'scan {done}/{total}', flush=True)
            if done % (10*divisions) == 0 or done == total else None)

    comparison = dict(budget=20000., tau=REGION_TAU, threads=threads, repeats=repeats,
                      divisions=divisions, residual_mode='physical', refinement_checks=REFINEMENT_CHECKS, runs={})
    points = (np.indices(reference['states'].shape).reshape(2, -1).T+.5)*axis_bounds/divisions
    truth = reference['states'].ravel() == 1
    for name, trials in runs.items():
        comparison['runs'][name] = []
        for result, monitor in trials:
            np.testing.assert_allclose(result['axis_bounds'], axis_bounds, atol=1e-8, rtol=0.)
            monitor.validation(reference, result)
            outer_membership = np.zeros(len(points), dtype=bool)
            for row in result['outer']:
                outer_membership |= contains(points, halfspaces(row['vertices']))
            assert not np.any(truth & ~outer_membership), '认证外包络遗漏扫描可行点'
            assert monitor.validation_state['validation']['fr_percent'] == 0., '认证内域出现扫描不可行点'
            inner = unary_union([MultiPoint(row['vertices']).convex_hull for row in result['inner']])
            outer = unary_union([MultiPoint(row['vertices']).convex_hull for row in result['outer']])
            comparison['runs'][name].append(dict(total_seconds=result['timing']['total_seconds'],
                counts=dict(result['counts'], schemes=len(monitor.state['schemes'])),
                coverage_bound=result['coverage_bound'],
                fr_percent=monitor.validation_state['validation']['fr_percent'],
                mr_percent=monitor.validation_state['validation']['mr_percent'],
                volume_gap=(outer.area-inner.area)/outer.area))

    for name, trials in runs.items():
        chosen = sorted(trials, key=lambda trial: trial[0]['timing']['total_seconds'])[len(trials)//2]
        monitor = chosen[1]
        monitor.validation_state['comparison'] = comparison
        monitor.save()
        records = comparison['runs'][name]
        print(f'{name}: seconds={[round(r["total_seconds"], 3) for r in records]}, '
              f'SP={[r["counts"]["sp"] for r in records]}, '
              f'G={[r["counts"]["global_search"] for r in records]}, '
              f'schemes={[r["counts"]["schemes"] for r in records]}, '
              f'MR%={[r["mr_percent"] for r in records]}, '
              f'FR%={[r["fr_percent"] for r in records]}', flush=True)
    return comparison


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--divisions', type=int, default=80)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--output', type=Path, default=OUTPUT)
    parser.add_argument('--replay', choices=('baseline', 'known_first'))
    args = parser.parse_args()
    if args.replay:
        monitor = RunMonitor()
        monitor.load_recording(args.output/args.replay/'monitor.json.gz')
        monitor.replay()
    else:
        run_comparison(args.output, repeats=args.repeats, divisions=args.divisions, threads=args.threads)
