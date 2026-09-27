"""配电网可规划域：方向 MP2 → 选点 / 全局搜索 / SP → 最终 AC 校核。"""
from pathlib import Path
from time import perf_counter
from datetime import datetime
import json

import numpy as np
from threadpoolctl import threadpool_limits

from Network.case33bw import Case33
from Network.concept5 import Concept5
from Network.four_bus_five_corridor import FourBus
from Network.jiangkou import Jiangkou
from model import (DEFAULT_SOLVER_THREADS, MP_TIME_LIMIT, SP_TIME_LIMIT, RESIDUAL_TIME_LIMIT,
                   GridPhysics, MasterProblem, SubProblem, RemainingRegionModel)
from region import GEOMETRY_TOL, RegionState, contains
from vertify import validate_ac_region
from plot import BenchmarkResult, json_value, save_method_comparison, show_result
from monitor import RunMonitor

# 常用参数：直接修改这里后运行 main.py。
NETWORK = 'concept5'                # 'concept5'、'case33'、'fourbus' 或 'jiangkou'
UPGRADE_COUNT = 8                 # case33 的升级选项数量（0..32），37 条走廊始终可重构
BUDGETS = None                     # None 使用网架默认预算，也可填递增列表
SOLVER_THREADS = DEFAULT_SOLVER_THREADS  # 默认 20；单次求解线程上限
DIVISIONS = 16                      # 仅用于最终 AC 网格校核
REGION_TAU = .002                  # SOCP 连续域径向精度
REFINEMENT_CHECKS = 96             # 活动见证处理完后，每累计这些 SP 次数做全局搜索
RESIDUAL_MODE = 'auto'             # 有限预算 light，无限预算 physical
CASE_TIME_LIMIT = 300.             # 每种方法、每档预算的构域时限（秒）
RECOMPUTE = True                   # False 读取已有 result.npz
SHOW_UI = True                     # 实时显示构域过程，结束后可逐步回放
STEP_BY_STEP = True                # 每个算法步骤等待“执行下一步”；连续运行可在窗口切换
OUTPUT = None                      # None 按网架选择结果目录
ROOT = Path(__file__).resolve().parent


class RegionTimeout(RuntimeError):
    """构域超时，终止本次运行；已完成步骤由监视器保存。"""


def build_continuous_region(network, method, budget, bounds, *, tau=REGION_TAU,
                            residual_mode=RESIDUAL_MODE, time_limit=CASE_TIME_LIMIT, threads=SOLVER_THREADS, progress=None,
                            cuts=(), clock=perf_counter):
    """1 初始化 → 2 选点 → 3 全局搜索 → 4 SP 更新 → 5 汇总；AC 在 run 末尾。"""
    started = clock()
    deadline = started+time_limit
    bounds = np.asarray(bounds, dtype=float)
    residual_mode = ('physical' if np.isinf(budget) else 'light') if residual_mode == 'auto' else residual_mode
    phases = {'linear': ('linear',), 'socp': ('socp',), 'hybrid': ('linear', 'socp')}[method]
    cuts, counts = [np.asarray(c).copy() for c in cuts], dict(sp=0, cuts=0)
    d = len(bounds)

    def remaining_time(limit):
        remaining = deadline-clock()
        if remaining <= 0.:
            raise RegionTimeout(f'{method}/{phase}, budget={budget:g}: exceeded {time_limit:g}s')
        return min(limit, remaining)

    def report(event, **data):
        if progress is not None:
            if data.get('answer') is not None:
                data['answer'] = {k: v for k, v in data['answer'].items() if k not in ('state', 'cut')}
            progress(event, checkpoint=True, bounds=bounds, records=region.records.values(),
                     pool_size=len(region.cuts),
                     counts_algorithm=dict(sp=oracle.calls, cuts=len(region.cuts)-initial_cuts), **data)

    for phase in phases:
        # 1. 方向 MP2 初始化：其他负荷自由；上界建外域，可行点按网架入内域。
        phase_started = clock()
        equations = GridPhysics(network, phase)
        region = RegionState(bounds, network.power_limit, 0. if phase == 'linear' else tau, cuts)
        oracle = SubProblem(equations, threads=threads)
        initial_cuts = len(region.cuts)
        maximum_answer = coverage = witness = None
        sp_since_global = 0
        report('phase_start', method=method, phase=phase, budget=budget,
               message=f'1 初始化：{phase}，预算 {budget:g}')
        for index, direction in enumerate([*np.eye(d), np.ones(d)]):
            label = f'p{index+1}' if index < d else '总负荷'
            report('mp_start', mode='MP2', message=f'1 初始化：MP2 最大化{label}')
            problem = MasterProblem(equations, budget=budget, direction=direction,
                                    cuts=region.cuts, threads=threads)
            with problem.model:
                answer = problem.solve(time_limit=remaining_time(MP_TIME_LIMIT))
            report('query_end', answer=answer, status='infeasible' if answer is None else answer['status'],
                   message=f'MP2 {label}：已证不可行' if answer is None else f"MP2 {label}：{answer['status']}")
            if answer is None:
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
            report('feasible', message=f'MP2 {label}：加入认证点', point=answer['p'],
                   choice=network.decode_plan(x), point_reason=f'MP2 {label}的完整可行解')

        while maximum_answer is not None:
            remaining_time(np.inf)
            candidate = None
            # 2. 选点与网架：活动见证补同网架支撑；普通候选取最大实际总负荷。
            if witness is not None:
                x = witness['x']
                point = (1-region.tau)*witness['p']/bounds
                if region.covering_schemes([point], preferred=x)[0] is not None:
                    if witness['feasible']:
                        region.add_point(x, witness['p']/bounds)
                        report('feasible', message='加入全局搜索的物理认证点', point=witness['p'],
                               choice=network.decode_plan(x), point_reason='全局物理模型的可行见证')
                    witness = None
                else:
                    support = [q for q in region.witness_support(x, point)
                               if not contains([q], region.inner_equations(x))[0]]
                    point = max(support, key=lambda q: (float(q@bounds), tuple(q)))
                    candidate = x, point
                    point_reason = '全局见证的同网架支撑点中实际总负荷最大'

            if candidate is None:
                candidates = []
                for row in region.records.values():
                    targets = (1-region.tau)*row['outer']
                    owners = region.covering_schemes(targets, preferred=row['x'])
                    candidates.extend((row['x'], point) for point, owner in zip(targets, owners) if owner is None)

                # 3. 无普通候选或到检查间隔：x、p 都自由，寻找未知网架/内部空隙。
                if not candidates or sp_since_global >= REFINEMENT_CHECKS:
                    report('residual_start', mode='剩余域', point=None, latest_cut=None,
                           message='3 全局搜索：检查尚未覆盖的区域')
                    problem = RemainingRegionModel(equations, budget, bounds, region.total_bound,
                                                   region.cuts, region.inner_halfspaces(), region.tau,
                                                   axis_bounds=region.axis_bounds,
                                                   mode=residual_mode, threads=threads)
                    with problem.model:
                        witness = problem.solve(GEOMETRY_TOL, time_limit=remaining_time(RESIDUAL_TIME_LIMIT))
                    coverage = witness['bound']
                    report('residual_end', answer=witness, coverage_bound=coverage,
                           message='全局覆盖已认证' if witness['complete'] else '发现未覆盖见证')
                    if witness['complete']:
                        break
                    x = witness['x']
                    region.add_scheme(x, network.decode_plan(x), network.cost@x)
                    report('scheme', choice=network.decode_plan(x), message='登记全局搜索返回的网架')
                    sp_since_global = 0
                    continue

                candidate = max(candidates, key=lambda pair: (float(pair[1]@bounds), tuple(pair[0]), tuple(pair[1])))
                point_reason = '已知网架未覆盖顶点中实际总负荷最大'

            # 4. 固定这一对 x、p，调用一次 SP；认证点入该内域，联合割更新各外域。
            x, point = candidate
            report('point', mode='SP', message='4 SP：检查候选点', point=point*bounds,
                   choice=network.decode_plan(x), point_reason=point_reason)
            checked = oracle.solve(x, point*bounds, time_limit=remaining_time(SP_TIME_LIMIT[phase]))
            sp_since_global += 1
            report('sp_end', answer=checked, message='SP 返回可行证书' if checked['feasible'] else 'SP 返回联合有效割')
            if checked['feasible']:
                region.add_point(x, point)
                report('feasible', message='认证点加入当前网架内域')
            else:
                cut = checked['cut']
                region.apply_cut(cut)
                report('cut', message='联合割更新所有已知网架外域', cut=cut, selection=x,
                       point=point*bounds, choice=network.decode_plan(x))
                if witness is not None and cut[0]+cut[1:1+d]@witness['p']+cut[1+d:]@witness['x'] < -1e-12:
                    witness = None

        # 5. 仅由 MP 已证空域或全局覆盖证书结束；下一阶段只继承有效割。
        result = dict(method=phase, budget=float(budget), residual_mode=residual_mode, status='certified', tau=region.tau,
                      axis_bounds=region.axis_bounds.copy(), coverage_bound=coverage,
                      max_total=0. if maximum_answer is None else float(maximum_answer['p'].sum()),
                      max_total_bound=region.total_bound,
                      max_point=None if maximum_answer is None else maximum_answer['p'].copy(),
                      max_choice=None if maximum_answer is None else network.decode_plan(maximum_answer['x']),
                      max_cost=None if maximum_answer is None else float(network.cost@maximum_answer['x']),
                      counts=dict(sp=oracle.calls, cuts=len(region.cuts)-initial_cuts),
                      timing=dict(total_seconds=clock()-phase_started), **region.finish(True))
        report('region_end', message='5 构域完成', region=result, coverage_bound=coverage,
               coverage_complete=True, point=None, latest_cut=None)
        cuts = region.cuts
        for key in counts:
            counts[key] += result['counts'][key]
    result.update(method=method, counts=counts, timing=dict(total_seconds=clock()-started))
    return result


def run(network, *, budgets=None, divisions=DIVISIONS, recompute=RECOMPUTE,
        show_ui=SHOW_UI, output=None, tau=REGION_TAU, residual_mode=RESIDUAL_MODE,
        time_limit=CASE_TIME_LIMIT, threads=SOLVER_THREADS, keep_ui=False, step_by_step=False):
    """依次构域、独立 AC 校核、保存并展示最终结果。"""
    # 1. 创建本次运行目录，或读取指定的历史结果
    root = ROOT/'results'/('concept5' if isinstance(network, Concept5) else network.name)
    output = Path(output) if output is not None else (root/datetime.now().strftime('%Y-%m-%d_%H%M%S_%f')
                                                      if recompute else max(root.glob('20*')))
    if not recompute:
        if isinstance(network, Concept5):
            result = json.loads((output/'results.json').read_text(encoding='utf-8'))
            if show_ui:
                import webbrowser
                webbrowser.open((output/'live_view.html').resolve().as_uri())
            return result
        result = BenchmarkResult.load(output, network=network)
        show_result(result, output, show_ui)
        return result

    budgets = np.asarray(network.budgets if budgets is None else budgets, dtype=float)

    # 2. 启动统一记录和窗口；记录不依赖显示开关
    with threadpool_limits(limits=1), RunMonitor(show_ui=show_ui, output=output, step_by_step=step_by_step) as monitor:
        monitor('preparation', checkpoint=True, message='准备网架与公共评价箱', network=network.name,
                load_nodes=network.load_nodes, budgets=budgets, bounds=[135., 135.] if isinstance(network, Concept5) else None,
                budget_index=0)
        # 3. 二维勘察共用记录器、时钟和保存入口
        if isinstance(network, Concept5):
            from survey import run_survey
            result = run_survey(threads=threads, time_limit=time_limit, progress=monitor, clock=monitor.clock)
            monitor('saving', message='保存勘察结果与完整过程')
            temporary = output/'results.json.tmp'
            temporary.write_text(json.dumps(json_value(result), ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
            temporary.replace(output/'results.json')
            if not result['audit']['passed']:
                raise RuntimeError('Survey validation failed; see results.json audit')
            monitor('completed', message='勘察完成，可回看选点、割和逐轮信息域', survey=result['states'][-1],
                    all_candidates=[row for row in result['all_candidates'] if row['step'] == result['states'][-1]['step']],
                    query_count=result['query_count'], joint_cut_count=result['joint_cut_count'], seconds=result['seconds'],
                    audit=result['audit'], stop_reason=result['stop_reason'],
                    point=None, latest_cut=None)
            monitor.show_result(result, export=show_ui, keep_ui=keep_ui)
            return result
        # 4. 三维构域与独立 AC 校核
        bounds = np.full(len(network.load_nodes), network.power_limit)
        result = BenchmarkResult.create(network, budgets, divisions, bounds, tau=tau,
                                        residual_mode=residual_mode)
        method = 'socp'
        for index, budget in enumerate(budgets):
            monitor('method_start', method=method, budget=budget, budget_index=index, bounds=bounds,
                    message=f'{method} · 预算 {budget:g} · 方向 MP2 → 选点 / 剩余域 / SP')
            region = build_continuous_region(network, method, budget, bounds, tau=tau,
                                             residual_mode=residual_mode, time_limit=time_limit,
                                             threads=threads, progress=monitor, clock=monitor.clock)
            result.add_region(region, index)
            monitor('method_end', message='本档预算构域结束', seconds=region['timing']['total_seconds'], results=result.metadata['continuous'])
        monitor('method_start', method='ac', phase='AC', mode='AC', message='独立 AC 网格校核', budget_index=0)
        started = perf_counter()
        ac = validate_ac_region(network, budgets, divisions, bounds, threads=threads)
        result.add_validation(ac, perf_counter()-started)

        # 5. 保存最终结果，导出同一记录的离线回放
        monitor('saving', message='保存结果与完整过程', states=ac, metadata=result.metadata)
        result.save(output)
        if show_ui:
            save_method_comparison(result, output)
        monitor('completed', message='计算完成，可回放每一步切割', results=result.metadata['continuous'])
        monitor.show_result(result, export=show_ui, keep_ui=keep_ui)
    return result


def main():
    network = {'case33': lambda: Case33(upgrade_count=UPGRADE_COUNT),
               'fourbus': FourBus, 'jiangkou': Jiangkou, 'concept5': Concept5}[NETWORK]()
    return run(network, budgets=BUDGETS, divisions=DIVISIONS, recompute=RECOMPUTE,
        show_ui=SHOW_UI, output=OUTPUT, tau=REGION_TAU, residual_mode=RESIDUAL_MODE,
        time_limit=CASE_TIME_LIMIT, threads=SOLVER_THREADS, keep_ui=SHOW_UI,
        step_by_step=STEP_BY_STEP and SHOW_UI)


if __name__ == '__main__':
    main()
