"""配电网可规划域：参数 → MP2 / MP1 / SP → AC 校核 → 最终结果。"""
from pathlib import Path
from time import perf_counter

import numpy as np
from threadpoolctl import threadpool_limits

from Network.case33bw import Case33
from Network.four_bus_five_corridor import FourBus
from model import DEFAULT_SOLVER_THREADS, evaluation_bounds
from region import ContinuousRegion, RegionTimeout
from vertify import validate_ac_region
from plot import BenchmarkResult, show_result

# 常用参数：直接修改这里后运行 main.py。
NETWORK = 'case33'                 # 'case33' 或 'fourbus'
CANDIDATE_COUNT = 32               # case33 支持 4 / 8 / 16 / 32
BUDGETS = None                     # None 使用网架默认预算，也可填递增列表
SOLVER_THREADS = DEFAULT_SOLVER_THREADS  # 默认 20；单次求解线程上限
DIVISIONS = 8                      # 仅用于最终 AC 网格校核
REGION_TAU = .002                  # SOCP 连续域径向精度
RESIDUAL_MODE = 'auto'             # 有限预算 light，无限预算 physical
CASE_TIME_LIMIT = 300.             # 每种方法、每档预算的构域时限（秒）
RECOMPUTE = True                   # False 读取已有 result.npz
SHOW_UI = True                     # 完成后打开最终结果页面
OUTPUT = None                      # None 按网架选择结果目录
ROOT = Path(__file__).resolve().parent


def solve_region(solver):
    """MP2 最大负荷 → MP1 最小投资 → SP 构域与全局覆盖检查。"""
    maximum = coverage = None
    try:
        # MP2 给出最大总负荷的可行下界和全局上界。
        answer = solver.call_query()
        if answer is None:
            solver.region.total_bound = 0.
            return solver.finish('certified', None, 0.)
        if not answer['feasible']:
            return solver.finish('unknown', None, None)
        solver.maximum = answer
        solver.region.total_bound = min(solver.region.total_bound, answer['bound'])
        maximum = float(answer['p'].sum())

        # MP1 固定总量，允许各节点重新分配负荷，最小化投资。
        cheapest = solver.call_query(min_total=max(0., maximum-1e-6),
                                     start=answer['x'], incumbent=answer)
        seeds = [cheapest, answer] if cheapest is not None and cheapest['feasible'] else [answer]
        for seed in seeds:
            if not solver.refine(seed['x'], seed['p']):
                return solver.finish('unknown', None, maximum)

        # SP 认证候选或生成有效割；剩余域搜索检查全部方案及内部空隙。
        while True:
            witness = solver.residual()
            coverage = witness['bound']
            if witness['complete']:
                return solver.finish('certified', coverage, maximum)
            if witness['x'] is None:
                return solver.finish('unknown', coverage, maximum)
            before = solver.region.progress
            seed = witness['p'] if witness['feasible'] else None
            if (not solver.refine(witness['x'], seed, witness=witness['p'])
                    or before == solver.region.progress):
                return solver.finish('unknown', coverage, maximum)
    except RegionTimeout:
        return solver.finish('time_limit', coverage, maximum)


def build_continuous_region(network, method, budget, bounds, *, tau=REGION_TAU,
                            residual_mode=RESIDUAL_MODE, time_limit=CASE_TIME_LIMIT, threads=SOLVER_THREADS):
    """混合方法先线性后 SOCP，只传递有效割，SOCP 内域重新认证。"""
    started = perf_counter()
    phases = {'linear': ('linear',), 'socp': ('socp',), 'hybrid': ('linear', 'socp')}[method]
    cuts, counts = [], dict(sp=0, cuts=0)
    for phase in phases:
        solver = ContinuousRegion(network, phase, budget, bounds, tau=tau, cuts=cuts,
                                  residual_mode=residual_mode, threads=threads,
                                  time_limit=max(0., time_limit-(perf_counter()-started)))
        result = solve_region(solver)
        cuts = solver.region.cuts
        for key in counts:
            counts[key] += result['counts'][key]
    result.update(method=method, counts=counts, timing=dict(total_seconds=perf_counter()-started))
    return result


def run(network, *, budgets=None, divisions=DIVISIONS, recompute=RECOMPUTE,
        show_ui=SHOW_UI, output=None, tau=REGION_TAU, residual_mode=RESIDUAL_MODE,
        time_limit=CASE_TIME_LIMIT, threads=SOLVER_THREADS):
    """依次构域、独立 AC 校核、保存并展示最终结果。"""
    output = Path(output) if output is not None else ROOT/'results'/network.name/f'planning_{len(network.planning_corridors)}'
    if not recompute:
        result = BenchmarkResult.load(output)
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
    with threadpool_limits(limits=1):
        print('准备网架与公共评价箱', flush=True)
        bounds = evaluation_bounds(network, threads=threads)
        result = BenchmarkResult.create(network, budgets, divisions, bounds, tau=tau,
                                        residual_mode=residual_mode, time_limit=time_limit, threads=threads)
        for method in ('linear', 'socp', 'hybrid'):
            for index, budget in enumerate(budgets):
                print(f'{method} · 预算 {budget:g} · MP2 → MP1 → SP / 剩余域', flush=True)
                region = build_continuous_region(network, method, budget, bounds, tau=tau,
                                                 residual_mode=residual_mode, time_limit=time_limit,
                                                 threads=threads)
                result.add_region(region, index)
        print('独立 AC 校核', flush=True)
        started = perf_counter()
        ac = validate_ac_region(network, budgets, divisions, bounds, threads=threads)
        result.add_validation(ac, perf_counter()-started)

    result.save(output)
    show_result(result, output, show_ui)
    return result


def main():
    network = {'case33': lambda: Case33(candidate_count=CANDIDATE_COUNT), 'fourbus': FourBus}[NETWORK]()
    run(network, budgets=BUDGETS, divisions=DIVISIONS, recompute=RECOMPUTE,
        show_ui=SHOW_UI, output=OUTPUT, tau=REGION_TAU, residual_mode=RESIDUAL_MODE,
        time_limit=CASE_TIME_LIMIT, threads=SOLVER_THREADS)


if __name__ == '__main__':
    main()
