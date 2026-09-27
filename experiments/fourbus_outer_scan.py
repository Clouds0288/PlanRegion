"""独立扫描 FourBus SOCP 参考域，并重建 61+20 条割的完整迭代历史。

python experiments/fourbus_outer_scan.py
扫描的每个 p 都允许全部合法 x、y 变化；不使用实验割求参考值。
有限扫描不替代连续域证书，也不代表 AC 真值。
"""
import argparse
import csv
import json
from pathlib import Path
import sys
from time import perf_counter

import numpy as np
from gurobipy import GRB

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from Network.four_bus_five_corridor import FourBus
from model import GridPhysics, MasterProblem, PLANNING_TOL, new_model


def point_query(problem, power):
    """固定负荷，保持全部 x/y 自由；异常求解不冒充不可行。"""
    problem.power.LB = problem.power.UB = power
    problem.model.optimize()
    if problem.model.Status == GRB.INFEASIBLE:
        return False
    if problem.model.Status != GRB.OPTIMAL:
        raise RuntimeError(f'Point {power}: status={problem.model.Status}')
    if problem.model.MaxVio > PLANNING_TOL:
        raise RuntimeError(f'Point {power}: status={problem.model.Status}, '
                           f'MaxVio={problem.model.MaxVio}')
    return True


def candidate_eta(equations, x, power):
    """事后复核历史候选：原始 SOCP 的 min eta，仅 y/eta 变化。"""
    with new_model('candidate_SP_audit', 1) as model:
        model.Params.Aggregate = 0
        model.Params.ScaleFlag = 0
        model.Params.BarQCPConvTol = 1e-9
        eta = model.addVar(name='violation')
        equations.add_operation(model, dict(zip(equations.keys, x)),
                                dict(zip(equations.network.load_nodes, power)), eta)
        model.setObjective(eta)
        model.optimize()
        if model.Status != GRB.OPTIMAL or model.MaxVio > PLANNING_TOL:
            raise RuntimeError(f'SP audit status={model.Status}, MaxVio={model.MaxVio}')
        return float(eta.X)


def outer_problem(equations, result, count):
    problem = MasterProblem(equations, budget=result['budget'], cuts_only=True,
                            cuts=np.array(result['cuts'][:count]), threads=1)
    problem.model.addConstr(problem.power <= np.array(result['axis_bounds']))
    problem.model.addConstr(problem.power.sum() <= result['total_bound'])
    return problem


def scan_grid(equations, result, grid_step):
    # 1. 均匀体素中心；参考模型包含完整电气约束，不含任何待检验割。
    axes = [np.arange(grid_step/2, np.ceil(b/grid_step)*grid_step, grid_step)
            for b in result['axis_bounds']]
    power = np.stack(np.meshgrid(*axes, indexing='ij'), axis=-1).reshape(-1, 3)
    reference = np.zeros(len(power), dtype=np.int8)
    x = np.zeros((len(power), len(equations.keys)), dtype=np.int8)
    problem = MasterProblem(equations, budget=result['budget'], threads=1)
    problem.model.setObjective(0.)
    started = perf_counter()
    max_vio = 0.
    with problem.model:
        for index, p in enumerate(power):
            feasible = point_query(problem, p)
            reference[index] = 1 if feasible else -1
            if feasible:
                x[index] = np.rint(problem.x.X).astype(np.int8)
                max_vio = max(max_vio, problem.model.MaxVio)
            if (index+1) % 2000 == 0 or index+1 == len(power):
                print(f'SOCP grid {index+1}/{len(power)}, '
                      f'{perf_counter()-started:.1f} s', flush=True)

    # 2. 每个 p 在全部网架中的外域归属；用前缀割的嵌套性找首次排除轮次。
    cuts = np.array(result['cuts'])
    count = len(cuts)
    first_exclusion = np.full(len(power), count+1, dtype=np.int16)
    outside_initial = ((power > np.array(result['axis_bounds'])+PLANNING_TOL).any(axis=1)
                       | (power.sum(axis=1) > result['total_bound']+PLANNING_TOL))
    first_exclusion[outside_initial] = 0
    margins = cuts[:, 0, None]+cuts[:, 1:4]@power.T+cuts[:, 4:]@x.T
    witness_retained = ((margins >= -PLANNING_TOL).all(axis=0)
                        & (reference == 1) & ~outside_initial)
    problems = {}

    def retained(index, prefix):
        if prefix not in problems:
            problems[prefix] = outer_problem(equations, result, prefix)
            problems[prefix].model.setObjective(0.)
        return point_query(problems[prefix], power[index])

    pending = np.flatnonzero(~outside_initial & ~witness_retained)
    for position, index in enumerate(pending):
        if not retained(index, count):
            left, right = 0, count
            while right-left > 1:
                middle = (left+right)//2
                if retained(index, middle):
                    left = middle
                else:
                    right = middle
            first_exclusion[index] = right
        if (position+1) % 2000 == 0 or position+1 == len(pending):
            print(f'Outer projection grid {position+1}/{len(pending)}', flush=True)
    for problem in problems.values():
        problem.model.dispose()
    return dict(power=power, reference=reference, first_exclusion=first_exclusion,
                shape=np.array([len(a) for a in axes]), grid_step=grid_step,
                reference_max_vio=max_vio, seconds=perf_counter()-started)


def scan_rays(equations, result, divisions):
    # 3. 扫描负荷比例 w；分别在完整 SOCP 与最终联合割外域上求 max ρ。
    weights = np.array([(i, j, divisions-i-j) for i in range(divisions+1)
                        for j in range(divisions+1-i)], dtype=float)/divisions
    problems = [MasterProblem(equations, budget=result['budget'], threads=1),
                outer_problem(equations, result, len(result['cuts']))]
    links = []
    for problem in problems:
        problem.model.Params.NumericFocus = 3
        radial_total = problem.model.addVar(name='radial_total_kw')
        rows = [problem.model.addConstr(p == radial_total/3)
                for p in problem.loads.values()]
        problem.model.update()
        links.append((radial_total, rows))
    reference_total, outer_total, reference_bound, outer_bound = [], [], [], []
    reference_plan, outer_plan = [], []
    max_vio = [0., 0.]
    started = perf_counter()
    for index, w in enumerate(weights):
        for which, (problem, (radial_total, rows)) in enumerate(zip(problems, links)):
            for row, value in zip(rows, w):
                problem.model.chgCoeff(row, radial_total, -float(value))
            problem.model.optimize()
            if problem.model.Status != GRB.OPTIMAL or problem.model.MaxVio > PLANNING_TOL:
                raise RuntimeError(f'Ray {w}, model={which}: status={problem.model.Status}, '
                                   f'MaxVio={problem.model.MaxVio}')
            max_vio[which] = max(max_vio[which], problem.model.MaxVio)
            (reference_total if which == 0 else outer_total).append(float(radial_total.X))
            (reference_bound if which == 0 else outer_bound).append(
                float(problem.model.ObjBound*equations.network.base))
            (reference_plan if which == 0 else outer_plan).append(
                np.rint(problem.x.X).astype(int).tolist())
        if (index+1) % 100 == 0 or index+1 == len(weights):
            print(f'Boundary rays {index+1}/{len(weights)}, '
                  f'{perf_counter()-started:.1f} s', flush=True)
    for problem in problems:
        problem.model.dispose()
    return dict(weights=weights.tolist(), reference_total=reference_total,
                outer_total=outer_total, reference_bound=reference_bound,
                outer_bound=outer_bound, reference_plan=reference_plan,
                outer_plan=outer_plan, max_vio=max_vio, seconds=perf_counter()-started)


def run_scan(seed, final, output, grid_step=2., divisions=40):
    started = perf_counter()
    seed = json.loads(Path(seed).read_text(encoding='utf-8'))
    result = json.loads(Path(final).read_text(encoding='utf-8'))
    inherited = result['initial_cut_count']
    if (seed['cuts'] != result['cuts'][:inherited] or seed['budget'] != result['budget']
            or seed['axis_bounds'] != result['axis_bounds']):
        raise ValueError('The supplied runs are not one nested cut sequence')
    trace = seed['trace']+result['trace']
    if len(trace) != len(result['cuts']):
        raise ValueError('Each plotted iteration must correspond to exactly one cut')
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    equations = GridPhysics(FourBus(), 'socp')

    grid = scan_grid(equations, result, grid_step)
    np.savez_compressed(output/'grid.npz', **grid)
    rays = scan_rays(equations, result, divisions)
    (output/'rays.json').write_text(json.dumps(rays, indent=2), encoding='utf-8')

    # 4. 重放全部 81 轮：固定历史 x,p 审核 min eta，再放开 x 判定 p 是否可行。
    schemes = {}

    def scheme_id(x):
        key = tuple(x)
        if key not in schemes:
            schemes[key] = f'S{len(schemes)+1:02d}'
        return schemes[key]

    initial = [dict(item, scheme=scheme_id(item['x'])) for item in result['initial']]
    effective_bound = np.inf
    iterations = []
    direct = MasterProblem(equations, budget=result['budget'], threads=1)
    direct.model.Params.NumericFocus = 3
    direct.model.Params.Aggregate = 0
    direct.model.setObjective(0.)
    with direct.model:
        for step, item in enumerate(trace, 1):
            effective_bound = min(effective_bound, item['bound'])
            outer = grid['first_exclusion'] > step
            p = np.array(item['p'])
            row = dict(step=step, phase=1 if step <= inherited else 2,
                       scheme=scheme_id(item['x']), p1=p[0], p2=p[1], p3=p[2],
                       eta=candidate_eta(equations, item['x'], p),
                       bound=item['bound'], effective_bound=effective_bound,
                       status=item['status'], seconds=item['seconds'],
                       cut_at_candidate=item['cut_at_candidate'],
                       p_feasible=point_query(direct, p),
                       outer_points=int(outer.sum()),
                       outer_only=int((outer & (grid['reference'] == -1)).sum()),
                       false_exclusion=int((~outer & (grid['reference'] == 1)).sum()))
            iterations.append(row)
    with (output/'iterations.csv').open('w', encoding='utf-8-sig', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(iterations[0]))
        writer.writeheader()
        writer.writerows(iterations)

    # 5. 数值对照与可复核报告；采样指标与未完成的全局证书分别报告。
    difference = np.array(rays['outer_total'])-np.array(rays['reference_total'])
    final_outer = grid['first_exclusion'] > len(trace)
    initial_outer = grid['first_exclusion'] > 0
    summary = dict(budget=result['budget'], epsilon=result['epsilon'],
                   certified=result['certified'], cut_count=len(trace), restart_after=inherited,
                   grid_step=grid_step, grid_shape=grid['shape'].tolist(),
                   grid_points=len(grid['power']), reference_points=int((grid['reference'] == 1).sum()),
                   initial_outer_points=int(initial_outer.sum()), final_outer_points=int(final_outer.sum()),
                   initial_outer_only=int((initial_outer & (grid['reference'] == -1)).sum()),
                   final_outer_only=iterations[-1]['outer_only'],
                   false_exclusion=iterations[-1]['false_exclusion'],
                   ray_count=len(rays['weights']), max_radial_difference=float(difference.max()),
                   mean_radial_difference=float(difference.mean()),
                   min_radial_difference=float(difference.min()),
                   max_relative_radial_difference=float((difference/np.array(rays['reference_total'])).max()),
                   last_recorded_bound=iterations[-1]['bound'], effective_bound=effective_bound,
                   candidate_p_feasible=sum(item['p_feasible'] for item in iterations),
                   visited_schemes=len(schemes), initial=initial,
                   schemes=[dict(id=label, x=list(x), plan=equations.network.decode_plan(np.array(x)))
                            for x, label in schemes.items()],
                   experiment_seconds=seed['seconds']+result['seconds'],
                   scan_seconds=perf_counter()-started,
                   figure_contract=dict(backend='python', width_mm=180, font_pt=8,
                                        exports=['png', 'svg', 'pdf'],
                                        claim='Finite SOCP scan agreement and unresolved global residual bound'))
    (output/'summary.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding='utf-8')
    write_report(output, summary, iterations)
    print(json.dumps({k: v for k, v in summary.items()
                      if k not in ('initial', 'schemes', 'figure_contract')}, indent=2), flush=True)
    return summary


def write_report(output, summary, iterations):
    s = summary
    lines = ['# FourBus 外域：独立 SOCP 扫描与完整迭代记录', '',
             f"预算 {s['budget']:g} 元；{s['cut_count']} 条联合割；"
             f"第 {s['restart_after']} 轮后续算。主线代码未参与扫描编排。", '',
             '## 计算方法', '',
             '1. 初始完整 MP2：分别最大化 p1、p2、p3 和总负荷，x、p、y 自由。',
             '2. 每轮 G 在已有联合割外域内，同时选择 x、p，最大化已对 y 最小化的 eta。',
             '3. 固定这一对 x、p，原始 SP 只调整 y/eta；本次每轮均生成一条联合有效割。',
             '4. 加割后重新求 G，允许换网架；某一对 x、p 不可行不等于这个 p 在所有网架不可行。',
             '5. 独立参考扫描固定 p、放开全部 x/y 求完整 SOCP。射线扫描固定负荷比例，最大化总负荷。', '',
             'eta 沿用原始 SP 残差尺度，不是 kW。表中的 eta 是对历史候选的事后原始 SOCP 复核；'
             '前 61 轮的原始 G 数值目标不用于此列。有效 G 上界取嵌套外域历史界的累积最小值；'
             '续算程序本身未继承第 61 轮的界，原始续算上界另存 CSV。', '',
             '## 扫描对照', '',
             f"- 体素中心间距 {s['grid_step']:g} kW，{s['grid_points']:,} 点；"
             f"参考可行 {s['reference_points']:,} 点。",
             f"- 初始外域多保留 {s['initial_outer_only']:,} 个参考不可行点，最终多保留 "
             f"{s['final_outer_only']:,} 个；最终误排参考可行点 {s['false_exclusion']} 个。",
             f"- {s['ray_count']} 条负荷比例射线：总负荷边界最大超估 "
             f"{s['max_radial_difference']:.6g} kW，平均超估 {s['mean_radial_difference']:.6g} kW；"
             f"最小差 {s['min_radial_difference']:.6g} kW。",
             f"- G 实际访问 {s['visited_schemes']} 个方案（含初始化）；"
             f"{s['candidate_p_feasible']}/{s['cut_count']} 个候选负荷在其他网架下可行。",
             f"- 有效全局上界 {s['effective_bound']:.8g}，停止阈值 {s['epsilon']:g}；"
             f"certified={s['certified']}。有限扫描吻合不能补成连续域认证。", '',
             '这里的参考值是同一 SOCP 模型的数值解，不是 AC 校验。'
             '体素计数乘体素体积仅是体积估计；射线之间绘制的面仅用于显示插值。', '',
             '## 图与图注', '',
             '![外域与 SOCP 参考对照](region_comparison.png)', '',
             '![全部 81 轮迭代](iteration_history.png)', '',
             '![负荷截面的外域收缩](region_evolution.png)', '',
             '[完整图注与数据口径](figure_captions.md)。同名 SVG/PDF 为可编辑矢量图。', '',
             '## 初始化', '', '| 目标 | 网架 | 认证负荷 (kW) |', '|---|---|---|']
    for item in s['initial']:
        lines.append(f"| {item['direction']} | {item['scheme']} | "
                     f"({', '.join(f'{v:.5f}' for v in item['p'])}) |")
    lines += ['', '## 方案编号', '', '| 编号 | 01 | 12 | 13 | 02 | 23 |', '|---|---|---|---|---|---|']
    for item in s['schemes']:
        lines.append('| '+item['id']+' | '+' | '.join(item['plan'][key] or '—'
                                                   for key in ('01', '12', '13', '02', '23'))+' |')
    lines += ['', '## 全部迭代', '',
              '“p 全局可行”来自独立完整 SOCP 检查，不是原算法的选点依据。'
              '最后一列是在这一轮加割之后仍保留的参考不可行网格点数。', '',
              '| 轮次 | 网架 | p (kW) | SP eta | 有效 G 上界 | p 全局可行 | 外域多留点 |',
              '|---|---|---|---|---|---|---|']
    for item in iterations:
        p = ', '.join(f"{item[key]:.4f}" for key in ('p1', 'p2', 'p3'))
        lines.append(f"| {item['step']} | {item['scheme']} | ({p}) | {item['eta']:.6g} | "
                     f"{item['effective_bound']:.6g} | {'是' if item['p_feasible'] else '否'} | "
                     f"{item['outer_only']} |")
    (output/'report.md').write_text('\n'.join(lines)+'\n', encoding='utf-8')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed', default='results/fourbus_outer/run_20000_carried_bound/results.json')
    parser.add_argument('--final', default='results/fourbus_outer/final/results.json')
    parser.add_argument('--output', default='results/fourbus_outer/scan')
    parser.add_argument('--grid-step', type=float, default=2.)
    parser.add_argument('--divisions', type=int, default=40)
    run_scan(**vars(parser.parse_args()))
