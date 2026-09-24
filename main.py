"""配电网可规划域：参数 → MP2 / MP1 / SP → AC 校核 → 最终结果。"""
from pathlib import Path
from time import perf_counter

import numpy as np
from threadpoolctl import threadpool_limits

from Network.case33bw import Case33
from Network.four_bus_five_corridor import FourBus
from Network.jiangkou import Jiangkou
from model import (DEFAULT_SOLVER_THREADS, MP_TIME_LIMIT, SP_TIME_LIMIT, RESIDUAL_TIME_LIMIT,
                   PlanningEquations, PlanningModel, PlanningSP, RemainingRegionModel, evaluation_bounds)
from region import GEOMETRY_TOL, REFINEMENT_CHECKS, RegionState, contains
from vertify import validate_ac_region
from plot import BenchmarkResult, RunMonitor, save_method_comparison, show_result

# 常用参数：直接修改这里后运行 main.py。
NETWORK = 'fourbus'                 # 'case33'、'fourbus' 或 'jiangkou'
UPGRADE_COUNT = 8                 # case33 的升级选项数量（0..32），37 条走廊始终可重构
BUDGETS = None                     # None 使用网架默认预算，也可填递增列表
SOLVER_THREADS = DEFAULT_SOLVER_THREADS  # 默认 20；单次求解线程上限
DIVISIONS = 16                      # 仅用于最终 AC 网格校核
REGION_TAU = .002                  # SOCP 连续域径向精度
RESIDUAL_MODE = 'auto'             # 有限预算 light，无限预算 physical
CASE_TIME_LIMIT = 300.             # 每种方法、每档预算的构域时限（秒）
RECOMPUTE = True                   # False 读取已有 result.npz
SHOW_UI = True                     # 实时显示构域过程，结束后可逐步回放
OUTPUT = None                      # None 按网架选择结果目录
ROOT = Path(__file__).resolve().parent


class RegionTimeout(RuntimeError):
    """构域到时；保留已获得的内域证书。"""


def build_continuous_region(network, method, budget, bounds, *, tau=REGION_TAU,
                            residual_mode=RESIDUAL_MODE, time_limit=CASE_TIME_LIMIT, threads=SOLVER_THREADS, progress=None):
    """直接求 MP2 → MP1 → SP / 剩余域；hybrid 只向 SOCP 传递有效割。"""
    started = perf_counter()
    deadline = started+time_limit
    bounds = np.asarray(bounds, dtype=float)
    residual_mode = ('physical' if np.isinf(budget) else 'light') if residual_mode == 'auto' else residual_mode
    phases = {'linear': ('linear',), 'socp': ('socp',), 'hybrid': ('linear', 'socp')}[method]
    cuts, counts = [], dict(sp=0, cuts=0)

    def remaining_time(limit):
        remaining = deadline-perf_counter()
        if remaining <= 0.:
            raise RegionTimeout()
        return min(limit, remaining)

    def report(event, **data):
        if progress is not None:
            progress(event, bounds=bounds, records=region.records.values(),
                     pool_size=len(region.cuts), residual_mode=residual_mode,
                     counts_algorithm=dict(sp=oracle.calls, cuts=len(region.cuts)-initial_cuts), **data)

    for phase in phases:
        phase_started = perf_counter()
        equations = PlanningEquations(network, phase)
        region = RegionState(bounds, network.power_limit, 0. if phase == 'linear' else tau, cuts)
        oracle = PlanningSP(equations, threads=threads)
        initial_cuts = len(region.cuts)
        maximum = maximum_answer = coverage = None
        status, seeds = 'unknown', {}
        report('phase_start', method=method, phase=phase, budget=budget, message=f'{method} / {phase}：开始构域')
        try:
            # 1. MP2 最大总负荷；MP1 保持该总量，允许重分配并最小化投资。
            for mode in ('MP2', 'MP1'):
                remaining_time(np.inf)
                report('mp_start', mode=mode, message=f'{mode}：求解主问题')
                problem = PlanningModel(equations, budget=budget, cuts=region.cuts, threads=threads,
                                        min_total=None if mode == 'MP2' else max(0., maximum-1e-6))
                with problem.model:
                    answer = problem.solve(time_limit=remaining_time(MP_TIME_LIMIT), incumbent=maximum_answer)
                # 数值质量不足时显式调用 SP；正常 MP 解直接使用。
                if answer is not None and answer['x'] is not None and not answer['feasible']:
                    checked = oracle.solve(answer['x'], answer['p'], time_limit=remaining_time(SP_TIME_LIMIT[phase]))
                    if checked['feasible']:
                        answer.update(state=checked['state'], feasible=True, status='feasible')
                report('query_end', message=f'{mode}：已证不可行' if answer is None else f"{mode}：{answer['status']}")
                if answer is None or not answer['feasible']:
                    if mode == 'MP1':
                        continue
                    if answer is None:
                        region.total_bound = maximum = 0.
                        status = 'certified'
                    break
                if mode == 'MP2':
                    maximum_answer = answer
                    maximum = float(answer['p'].sum())
                    region.total_bound = min(region.total_bound, answer['bound'])
                # 证书立即入库；同方案保留所有可行点，但只启动一次初始细化。
                x = answer['x']
                region.add_scheme(x, network.decode_plan(x), network.cost@x)
                region.add_point(x, answer['p']/bounds)
                seeds[tuple(x)] = x
                report('feasible', message=f'{mode}：加入可行点，扩展内域', point=answer['p'])

            # 2. 先细化初始方案，再由剩余域模型提供未覆盖方案与见证。
            queue = list(seeds.values())
            while maximum_answer is not None:
                before = pending = seed = None
                if queue:
                    x = queue.pop()
                else:
                    report('residual_start', mode='剩余域', message='搜索尚未覆盖的负荷区域', point=None, latest_cut=None)
                    problem = RemainingRegionModel(equations, budget, bounds, region.total_bound,
                                                   region.cuts, region.inner_halfspaces(), region.tau,
                                                   mode=residual_mode, threads=threads)
                    with problem.model:
                        witness = problem.solve(GEOMETRY_TOL, time_limit=remaining_time(RESIDUAL_TIME_LIMIT))
                    coverage = witness['bound']
                    if witness['complete']:
                        status = 'certified'
                        break
                    if witness['x'] is None:
                        break
                    x, before = witness['x'], region.progress
                    pending = (1-region.tau)*witness['p']/bounds
                    seed = witness['p'] if witness['feasible'] else None
                region.add_scheme(x, network.decode_plan(x), network.cost@x)
                row = region.records[tuple(x)]
                report('scheme', mode='SP', message='展开当前线路方案', choice=row['choice'], point=None, latest_cut=None)

                # 3. 一批候选：选点 → 直接 SP → 加入内域或用联合割裁剪。
                failed = False
                for _ in range(REFINEMENT_CHECKS):
                    remaining_time(np.inf)
                    if not len(row['outer']):
                        break
                    point = None
                    if pending is not None:
                        if region.covering_schemes([pending], preferred=x)[0] is not None:
                            pending = None
                        else:
                            for candidate in region.witness_support(x, pending):
                                if not contains([candidate], region.inner_equations(x))[0]:
                                    point = candidate
                                    break
                            if point is None:
                                point, pending = pending, None
                    # 见证先补同方案支撑点，避免物理种子立即“自覆盖”。
                    if pending is None and seed is not None:
                        region.add_point(x, seed/bounds)
                        report('feasible', message='加入剩余域可行点，扩展内域', point=seed)
                        seed = None
                    if point is None:
                        point = region.next_point(x)
                    if point is None:
                        break
                    report('point', mode='SP', message='SP：检查候选负荷点', point=point*bounds, choice=row['choice'])
                    checked = oracle.solve(x, point*bounds, time_limit=remaining_time(SP_TIME_LIMIT[phase]))
                    if checked['feasible']:
                        region.add_point(x, point)
                        report('feasible', message='SP：候选点通过认证，扩展内域')
                    elif checked['cut'] is not None:
                        old, cut = row['outer'].copy(), checked['cut']
                        region.apply_cut(cut)
                        report('cut', message='SP：加入新割，裁剪各方案外域', cut=cut, selection=x,
                               point=point*bounds, choice=row['choice'])
                        if pending is not None and cut[0]+cut[1:4]@(pending*bounds/(1-region.tau))+cut[4:]@x < -1e-12:
                            pending = None
                        failed = np.array_equal(old, row['outer'])
                    else:
                        report('unknown', message='SP：尚未取得可行证书或有效割')
                        failed = True
                    if failed:
                        break
                if failed or before is not None and before == region.progress:
                    break
        except RegionTimeout:
            status = 'time_limit'
        if status == 'unknown' and perf_counter() >= deadline:
            status = 'time_limit'

        # 4. 单一出口保存已获证内域、全局界及本阶段统计。
        result = dict(method=phase, budget=float(budget), residual_mode=residual_mode, status=status, tau=region.tau,
                      coverage_bound=coverage, max_total=maximum, max_total_bound=region.total_bound,
                      max_point=None if maximum_answer is None else maximum_answer['p'].copy(),
                      max_choice=None if maximum_answer is None else network.decode_plan(maximum_answer['x']),
                      max_cost=None if maximum_answer is None else float(network.cost@maximum_answer['x']),
                      counts=dict(sp=oracle.calls, cuts=len(region.cuts)-initial_cuts),
                      timing=dict(total_seconds=perf_counter()-phase_started), **region.finish(status == 'certified'))
        report('region_end', message=f'本阶段结束：{status}', region=result, coverage_bound=coverage,
               coverage_complete=status == 'certified', point=None, latest_cut=None)
        cuts = region.cuts
        for key in counts:
            counts[key] += result['counts'][key]
    result.update(method=method, counts=counts, timing=dict(total_seconds=perf_counter()-started))
    return result


def run(network, *, budgets=None, divisions=DIVISIONS, recompute=RECOMPUTE,
        show_ui=SHOW_UI, output=None, tau=REGION_TAU, residual_mode=RESIDUAL_MODE,
        time_limit=CASE_TIME_LIMIT, threads=SOLVER_THREADS, keep_ui=False):
    """依次构域、独立 AC 校核、保存并展示最终结果。"""
    output = Path(output) if output is not None else ROOT/'results'/network.name/('planning_v3_'+network.fingerprint[:12])
    if not recompute:
        result = BenchmarkResult.load(output, network=network)
        show_result(result, output, show_ui)
        return result

    budgets = np.asarray(network.budgets if budgets is None else budgets, dtype=float)
    if (budgets.ndim != 1 or not len(budgets) or np.any(np.isnan(budgets))
            or np.any(budgets < 0) or np.any(budgets[1:] <= budgets[:-1])
            or not isinstance(divisions, (int, np.integer)) or divisions < 1
            or not isinstance(threads, (int, np.integer)) or threads < 1
            or not 0 <= tau < 1 or np.isnan(time_limit) or time_limit < 0
            or residual_mode not in ('auto', 'light', 'physical')):
        raise ValueError('预算须非负递增；网格数、线程数须为正整数；0≤tau<1；时限须非负；剩余域模式须为 auto/light/physical')

    # 小型矩阵运算保持单线程，优化器按 threads 配置并行。
    with threadpool_limits(limits=1), RunMonitor(show_ui=show_ui, output=output) as monitor:
        monitor('preparation', message='准备网架与公共评价箱', network=network.name, load_nodes=network.load_nodes, budgets=budgets)
        bounds = evaluation_bounds(network, threads=threads)
        result = BenchmarkResult.create(network, budgets, divisions, bounds, tau=tau,
                                        residual_mode=residual_mode, time_limit=time_limit, threads=threads)
        for method in ('linear', 'socp', 'hybrid'):
            for index, budget in enumerate(budgets):
                monitor('method_start', method=method, budget=budget, budget_index=index, bounds=bounds,
                        message=f'{method} · 预算 {budget:g} · MP2 → MP1 → SP / 剩余域')
                region = build_continuous_region(network, method, budget, bounds, tau=tau,
                                                 residual_mode=residual_mode, time_limit=time_limit,
                                                 threads=threads, progress=monitor if show_ui else None)
                result.add_region(region, index)
                monitor('method_end', message='本档预算构域结束', seconds=region['timing']['total_seconds'], results=result.metadata['continuous'])
        monitor('method_start', method='ac', phase='AC', mode='AC', message='独立 AC 网格校核', budget_index=0)
        started = perf_counter()
        ac = validate_ac_region(network, budgets, divisions, bounds, threads=threads)
        result.add_validation(ac, perf_counter()-started)

        monitor('saving', message='保存结果与完整过程', states=ac, metadata=result.metadata)
        result.save(output)
        if show_ui:
            save_method_comparison(result, output)
        monitor('completed', message='计算完成，可回放每一步切割', results=result.metadata['continuous'])
        monitor.show_result(result, export=show_ui, keep_ui=keep_ui)
    return result


def main():
    network = {'case33': lambda: Case33(upgrade_count=UPGRADE_COUNT),
               'fourbus': FourBus, 'jiangkou': Jiangkou}[NETWORK]()
    run(network, budgets=BUDGETS, divisions=DIVISIONS, recompute=RECOMPUTE,
        show_ui=SHOW_UI, output=OUTPUT, tau=REGION_TAU, residual_mode=RESIDUAL_MODE,
        time_limit=CASE_TIME_LIMIT, threads=SOLVER_THREADS, keep_ui=SHOW_UI)


if __name__ == '__main__':
    main()
