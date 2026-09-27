"""SOCP 规划域：1 初始化 → 2 顶点选点 → 3 全局搜索 → 4 SP 更新 → 5 校验。"""
from pathlib import Path
from time import perf_counter

import numpy as np
from threadpoolctl import threadpool_limits

from Network.four_bus_five_corridor import FourBus
from Network.case33bw import Case33
from model import (DEFAULT_SOLVER_THREADS, SP_TIME_LIMIT, RESIDUAL_TIME_LIMIT,
                   GridPhysics, MasterProblem, SubProblem, RemainingRegionModel)
from region import GEOMETRY_TOL, RegionState, contains
from vertify import validate_socp_region
from monitor import RunMonitor, RegionTimeout

# 案例入口：改为 Case33 即切换二维案例，背景负荷保持原值。
NETWORK = Case33                  # FourBus / Case33  显式 --case 可覆盖
BUDGET = 20000.
SOLVER_THREADS = DEFAULT_SOLVER_THREADS
DIVISIONS = 80                     # 最终独立 SOCP 扫描，每轴网格数
REGION_TAU = .005                 # 径向精度 0.5%；缩回认证内域的外域顶点跳过 SP
REFINEMENT_CHECKS = 32             # 新增 SP 次数；批次 / 活动见证结束后检查
RESIDUAL_MODE = 'physical'
CASE_TIME_LIMIT = 300.
SHOW_UI = True
ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT/'results'/'fourbus_2d'/'monitor.json.gz'


def build_continuous_region(network, method, budget, bounds, *, tau=REGION_TAU,
                            residual_mode=RESIDUAL_MODE, time_limit=CASE_TIME_LIMIT,
                            threads=SOLVER_THREADS, progress=None, cuts=(), clock=perf_counter):
    """普通顶点按 SP 最优违反量选割；只有全局覆盖证书允许终止。"""
    bounds = np.asarray(bounds, dtype=float)
    residual_mode = ('physical' if np.isinf(budget) else 'light') if residual_mode == 'auto' else residual_mode
    equations = GridPhysics(network, method)
    region = RegionState(bounds, network.power_limit, 0. if method == 'linear' else tau, cuts)
    oracle = SubProblem(equations, threads=threads)
    monitor = RunMonitor.follow(progress, clock=clock)
    monitor.begin(network, method, budget, region, time_limit)
    cache = {}
    initial_cuts = len(region.cuts)
    maximum_answer = coverage = witness = None
    sp_since_global = global_search = 0
    d = len(bounds)

    # 1. 初始化：d 个轴向 MP2 + 总负荷 MP2；x、p、y 均自由，点按各自网架认证。
    for index, direction in enumerate([*np.eye(d), np.ones(d)]):
        monitor.initializing(index)
        problem = MasterProblem(equations, budget=budget, direction=direction,
                                cuts=region.cuts, threads=threads)
        with problem.model:
            answer = problem.solve(time_limit=monitor.remaining(time_limit))
        if answer is None:                       # 数学结论：完整 MP2 已证空域。
            region.tighten_bounds(np.zeros(d), 0.)
            break
        axis_bounds = region.axis_bounds.copy()
        if index < d:
            axis_bounds[index] = min(axis_bounds[index], answer['bound'])
            region.tighten_bounds(axis_bounds, region.total_bound)
        else:
            maximum_answer = answer
            region.tighten_bounds(axis_bounds, answer['bound'])
        x = answer['x']
        region.add_scheme(x, network.decode_plan(x), network.cost@x)
        region.add_point(x, answer['p']/bounds)
        monitor.seed(answer, region)

    while maximum_answer is not None:
        candidates = None
        # 2. 选点：先补全局见证的同网架支撑；普通外域顶点缩回 tau 后已覆盖则跳过。
        if witness is not None:
            x = witness['x']
            point = (1-region.tau)*witness['p']/bounds
            if region.covering_schemes([point], preferred=x)[0] is not None:
                if witness['feasible']:
                    region.add_point(x, witness['p']/bounds)
                witness = None
            else:
                support = [q for q in region.witness_support(x, point)
                           if not contains([q], region.inner_equations(x))[0]]
                point = max(support, key=lambda q: (float(q@bounds), tuple(q)))
                candidates = [(x, point)]
        supporting = candidates is not None
        if not supporting:
            candidates = []
            for row in region.records.values():
                targets = (1-region.tau)*row['outer']
                owners = region.covering_schemes(targets, preferred=row['x'])
                candidates.extend((row['x'], point) for point, owner
                                  in zip(targets, owners) if owner is None)

            # 3. 全局搜索：无顶点可查或累计 32 次新增 SP；同一 tau 判覆盖，x、p 自由。
            if not candidates or sp_since_global >= REFINEMENT_CHECKS:
                monitor.global_start(sp_since_global)
                problem = RemainingRegionModel(equations, budget, bounds, region.total_bound,
                    region.cuts, region.inner_halfspaces(), region.tau,
                    axis_bounds=region.axis_bounds, mode=residual_mode, threads=threads)
                with problem.model:
                    witness = problem.solve(GEOMETRY_TOL, time_limit=monitor.remaining(RESIDUAL_TIME_LIMIT))
                global_search += 1
                coverage = witness['bound']
                if not witness['complete']:
                    x = witness['x']
                    region.add_scheme(x, network.decode_plan(x), network.cost@x)
                monitor.global_end(witness, region)
                if witness['complete']:
                    break
                sp_since_global = 0
                continue

        # 4. 固定候选 x、p，让 y 自由求最小 eta；精确缓存避免重复 SP。
        monitor.selecting(len(candidates), supporting)
        scored = []
        for x, point in candidates:
            power = point*bounds
            key = tuple(x), tuple(power)
            if key not in cache:
                monitor.sp_start(x, power, oracle.calls+1)
                cache[key] = oracle.solve(x, power, time_limit=monitor.remaining(SP_TIME_LIMIT[method]))
                sp_since_global += 1
                monitor.sp_end(cache[key])
            checked = cache[key]
            scored.append((x, point, checked))
            if checked['feasible']:
                region.add_point(x, point)
                monitor.updated(region, 'feasible', x=x, power=power, checked=checked)

        # 先收完本批认证点，再在仍未覆盖者中按最大原始 eta 采用一条联合割。
        if supporting:
            pending = [row for row in scored if not row[2]['feasible']]
        else:
            owners = region.covering_schemes([point for _, point, _ in scored])
            pending = [row for row, owner in zip(scored, owners) if not row[2]['feasible'] and owner is None]
        if pending:
            x, point, checked = max(pending, key=lambda row: (row[2]['eta'], float(row[1]@bounds),
                                                            tuple(row[0]), tuple(row[1])))
            cut = checked['cut']
            region.apply_cut(cut)
            monitor.updated(region, 'cut', x=x, power=point*bounds, checked=checked)
            if witness is not None and cut[0]+cut[1:1+d]@witness['p']+cut[1+d:]@witness['x'] < -1e-12:
                witness = None

    # 5. 全局证书给出最终内域并集 / 外包络；独立扫描在整个构域完成后执行。
    result = dict(method=method, budget=float(budget), residual_mode=residual_mode, status='certified',
        tau=region.tau, axis_bounds=region.axis_bounds.copy(), coverage_bound=coverage,
        max_total=0. if maximum_answer is None else float(maximum_answer['p'].sum()),
        max_total_bound=region.total_bound, counts=dict(sp=oracle.calls,
        cuts=len(region.cuts)-initial_cuts, global_search=global_search),
        timing=monitor.timing(), **region.finish(True))
    monitor.finish(result, region)
    return result


def run(network, *, budget=BUDGET, divisions=DIVISIONS, show_ui=SHOW_UI, output=OUTPUT,
        tau=REGION_TAU, time_limit=CASE_TIME_LIMIT, threads=SOLVER_THREADS):
    """窗口与回放由 monitor 承担；主线只编排构域、独立扫描、提交结果。"""
    monitor = RunMonitor(output=output)

    def calculate():
        with threadpool_limits(limits=1):
            bounds = np.full(len(network.load_nodes), network.power_limit)
            result = build_continuous_region(network, 'socp', budget, bounds, tau=tau,
                time_limit=time_limit, threads=threads, progress=monitor)
            reference = validate_socp_region(network, budget, divisions, result['axis_bounds'],
                                             threads=threads, progress=monitor.scanning)
            monitor.validation(reference, result)
            return result

    return monitor.execute(calculate, show_ui=show_ui)


def main(case=None):
    case = {FourBus: 'fourbus', Case33: 'case33'}[NETWORK] if case is None else case
    # 两个二维切片；未选中的节点负荷保持各案例的原始值。
    network_type, load_nodes, budget, divisions, time_limit = {
        'fourbus': (FourBus, (1, 2), BUDGET, DIVISIONS, CASE_TIME_LIMIT),
        'case33': (Case33, (18, 25), 2., 20, 900.),
    }[case]
    return run(network_type(load_nodes=load_nodes), budget=budget, divisions=divisions,
               show_ui=SHOW_UI, output=ROOT/'results'/f'{case}_2d'/'monitor.json.gz', tau=REGION_TAU,
               time_limit=time_limit, threads=SOLVER_THREADS)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='二维 SOCP 规划域与原生割线回放')
    parser.add_argument('--case', choices=('fourbus', 'case33'), help='覆盖顶部 NETWORK 的案例选择')
    main(parser.parse_args().case)
