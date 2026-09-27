"""Audit and plot the matched FourBus partition/bilinear timing experiment.

Figure contract: a compares global certificates, b compares cut throughput,
c compares paired-ray boundary errors after equal wall time. All observations
are retained. These are single deterministic runs, not statistical replicates.
"""
import argparse
import csv
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from gurobipy import GRB

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from experiments.fourbus_outer_scan import candidate_eta, outer_problem
from model import GridPhysics, MasterProblem, PLANNING_TOL
from Network.four_bus_five_corridor import FourBus


def audit_result(equations, result):
    """Check every cut globally and every candidate against the original primal SP."""
    direct = MasterProblem(equations, threads=1)
    margins = []
    with direct.model:
        direct.model.Params.NumericFocus = 3
        direct.model.Params.Aggregate = 0
        for cut in map(np.array, result['cuts']):
            direct.model.setObjective(cut[0]+cut[1:4]@direct.power+cut[4:]@direct.x, GRB.MINIMIZE)
            direct.model.optimize()
            if direct.model.Status != GRB.OPTIMAL or direct.model.MaxVio > PLANNING_TOL:
                raise RuntimeError(f'Cut audit status={direct.model.Status}, MaxVio={direct.model.MaxVio}')
            margins.append(float(direct.model.ObjBound))
    lower_errors, upper_errors, candidates = [], [], []
    for item in result['trace']:
        if 'candidate_lower' in item:
            eta = candidate_eta(equations, np.array(item['x']), np.array(item['p']))
            lower_errors.append(item['candidate_lower']-eta)
            upper_errors.append(eta-item['local_bound'])
            candidates.append(dict(step=item['step'], eta=eta,
                                   candidate_lower=item['candidate_lower'], local_bound=item['local_bound']))
    assert min(margins) >= -1e-7
    assert max(lower_errors) <= 2e-7
    assert max(upper_errors) <= 2e-7
    return dict(min_cut_margin=min(margins), cut_margins=margins,
                max_lower_error=max(lower_errors), max_upper_error=max(upper_errors),
                candidates=candidates)


def compare_rays(equations, result, reference):
    """Maximize load along every existing reference ray in this method's outer domain."""
    problem = outer_problem(equations, result, len(result['cuts']))
    weights = np.array(reference['weights'])
    outer_total = []
    with problem.model:
        problem.model.Params.NumericFocus = 3
        radial_total = problem.model.addVar(name='radial_total_kw')
        rows = [problem.model.addConstr(p == radial_total/3.) for p in problem.loads.values()]
        problem.model.setObjective(radial_total, GRB.MAXIMIZE)
        for w in weights:
            for row, value in zip(rows, w):
                problem.model.chgCoeff(row, radial_total, -float(value))
            problem.model.optimize()
            if problem.model.Status != GRB.OPTIMAL or problem.model.MaxVio > PLANNING_TOL:
                raise RuntimeError(f'Ray audit status={problem.model.Status}, MaxVio={problem.model.MaxVio}')
            outer_total.append(float(radial_total.X))
    reference_total = np.array(reference['reference_total'])
    difference = np.array(outer_total)-reference_total
    assert difference.min() >= -1e-5
    return dict(outer_total=outer_total, max_radial_difference=float(difference.max()),
                mean_radial_difference=float(difference.mean()),
                min_radial_difference=float(difference.min()))


def match_boundary_accuracy(equations, result, reference, target):
    """First cut prefix meeting the baseline's final error on every reference ray."""
    lower, upper = 0, len(result['cuts'])
    checked = {}
    while lower < upper:
        count = (lower+upper)//2
        candidate = dict(result, cuts=result['cuts'][:count])
        scan = compare_rays(equations, candidate, reference)
        checked[count] = scan['max_radial_difference']
        if checked[count] <= target+1e-5:
            upper = count
        else:
            lower = count+1
    seconds = next(item['seconds'] for item in result['trace'] if item['cuts'] >= lower) if lower else 0.
    return dict(cut_count=lower, seconds=seconds, target_max_radial_difference=target,
                comparison_tolerance_kw=1e-5,
                checked_prefixes=checked)


def render_comparison(results, rays, reference, output):
    """Paper-style 183 mm quantitative triptych; editable text, no invented intervals."""
    plt.rcParams.update({'font.family': 'sans-serif',
                         'font.sans-serif': ['Microsoft YaHei', 'Arial', 'DejaVu Sans'],
                         'font.size': 8, 'axes.labelsize': 8, 'xtick.labelsize': 8,
                         'ytick.labelsize': 8, 'legend.fontsize': 7.5,
                         'axes.spines.right': False, 'axes.spines.top': False,
                         'axes.linewidth': .7, 'lines.linewidth': 1.4,
                         'svg.fonttype': 'none', 'pdf.fonttype': 42})
    fig, axes = plt.subplots(1, 3, figsize=(7.2047244, 2.5984252))  # 183 x 66 mm
    colors = {'bilinear': '#355F7D', 'partition': '#BD743B'}
    names = {'bilinear': '双线性全局搜索', 'partition': '角点分块搜索'}
    styles = {'bilinear': '-', 'partition': '--'}
    for method, result in results.items():
        times = np.r_[0., [item['seconds'] for item in result['trace']]]
        bounds = np.r_[max(result['axis_bounds'])/150., [item['bound'] for item in result['trace']]]
        counts = np.r_[result['initial_cut_count'], [item['cuts'] for item in result['trace']]]
        assert np.all(np.diff(times) >= 0.) and np.all(bounds > 0.)
        axes[0].step(times, bounds, where='post', color=colors[method], linestyle=styles[method], label=names[method])
        axes[1].step(times, counts, where='post', color=colors[method], linestyle=styles[method])
        differences = np.array(rays[method]['outer_total'])-np.array(reference['reference_total'])
        # Tiny negative values from solver tolerance are retained, not clipped.
        axes[2].step(np.sort(differences), 100*np.arange(1, len(differences)+1)/len(differences),
                     where='post', color=colors[method], linestyle=styles[method])
    epsilon = next(iter(results.values()))['epsilon']
    axes[0].axhline(epsilon, color='#767676', linestyle=':', linewidth=.9)
    axes[0].text(.98, .035, 'ε = 0.0001', transform=axes[0].transAxes, ha='right', fontsize=7.5)
    axes[0].set(yscale='log', ylim=(epsilon*.7, 1.), xlabel='运行时间（s）', ylabel='全局上界（η 尺度）')
    axes[1].set(xlabel='运行时间（s）', ylabel='有效割数量')
    axes[2].set(xlabel='径向容量高估（kW）', ylabel='累计射线比例（%）', ylim=(0, 103))
    for ax in axes[:2]:
        ax.set_xlim(0., max(result['seconds'] for result in results.values()))
        ax.set_xticks([0, 100, 200, 300])
    for label, ax in zip('abc', axes):
        ax.text(-.18, 1.06, label, transform=ax.transAxes, fontsize=9, fontweight='bold')
        ax.tick_params(width=.7, length=3)
    fig.legend(*axes[0].get_legend_handles_labels(), loc='upper center', ncol=2,
               frameon=False, bbox_to_anchor=(.5, 1.0))
    fig.subplots_adjust(left=.085, right=.988, bottom=.23, top=.82, wspace=.61)
    fig.savefig(output/'comparison.svg', facecolor='white')
    fig.savefig(output/'comparison.pdf', facecolor='white')
    fig.savefig(output/'comparison.png', dpi=600, facecolor='white')
    plt.close(fig)


def make_report(directory, reference_file):
    output = Path(directory)
    results = {method: json.loads((output/method/'results.json').read_text(encoding='utf-8'))
               for method in ('bilinear', 'partition')}
    left, right = results.values()
    assert left['budget'] == right['budget'] and left['epsilon'] == right['epsilon']
    assert left['initial'] == right['initial'] and left['threads'] == right['threads']
    reference = json.loads(Path(reference_file).read_text(encoding='utf-8'))
    equations = GridPhysics(FourBus(), 'socp')
    audits, rays = {}, {}
    for method, result in results.items():
        audits[method] = audit_result(equations, result)
        rays[method] = compare_rays(equations, result, reference)
        print(method, dict(seconds=result['seconds'], bound=result['bound'], cuts=len(result['cuts']),
                           max_radial_difference=rays[method]['max_radial_difference'],
                           mean_radial_difference=rays[method]['mean_radial_difference']), flush=True)
    summary = {method: dict(seconds=result['seconds'], bound=result['bound'], certified=result['certified'],
                            counts=result['counts'], cuts=len(result['cuts']),
                            first_bound_005_seconds=next((item['seconds'] for item in result['trace']
                                                         if item['bound'] <= .05), None),
                            global_seconds=sum(item['global_seconds'] for item in result['trace']),
                            **{key: value for key, value in rays[method].items() if key != 'outer_total'})
               for method, result in results.items()}
    accuracy_match = match_boundary_accuracy(equations, right, reference,
                                             rays['bilinear']['max_radial_difference'])
    baseline_match = match_boundary_accuracy(equations, left, reference,
                                             rays['bilinear']['max_radial_difference'])
    accuracy_match.update(baseline_cut_count=baseline_match['cut_count'],
                          baseline_seconds=baseline_match['seconds'],
                          baseline_checked_prefixes=baseline_match['checked_prefixes'],
                          speed_ratio=baseline_match['seconds']/accuracy_match['seconds'])
    summary['accuracy_match'] = accuracy_match
    print('accuracy_match', accuracy_match, flush=True)
    (output/'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    (output/'audit.json').write_text(json.dumps(audits, indent=2), encoding='utf-8')
    with (output/'rays_comparison.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream)
        writer.writerow(['w1', 'w2', 'w3', 'reference_total', 'bilinear_outer_total', 'partition_outer_total'])
        for i, w in enumerate(reference['weights']):
            writer.writerow([*w, reference['reference_total'][i],
                             rays['bilinear']['outer_total'][i], rays['partition']['outer_total'][i]])
    with (output/'convergence.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream)
        writer.writerow(['method', 'step', 'seconds', 'bound', 'cuts', 'action'])
        for method, result in results.items():
            for item in result['trace']:
                writer.writerow([method, *[item[key] for key in ('step', 'seconds', 'bound', 'cuts', 'action')]])
    render_comparison(results, rays, reference, output)
    lines = ['# FourBus：双线性与角点分块的同条件比较', '',
             f"预算 {left['budget']:,.0f} 元，SOCP，ε={left['epsilon']}（SP 原始残差尺度），"
             f"{left['threads']} 线程，每方法 {left['wall_limit']:.0f} 秒，单次全局搜索上限 {left['time_limit']:.0f} 秒。",
             f"两者共享同一次 4 个方向 MP2 初始化（{left['initial_seconds']:.3f} 秒），从零条割开始；"
             '比较时间包含各自建模、选点、SP、切割和保存，核查另计。每方法一次实际运行，无重复实验置信区间。', '',
             '| 指标 | 双线性 | 角点分块 |', '|---|---:|---:|']
    for label, key, fmt in [('实际运行秒数', 'seconds', '.3f'), ('全局 η 上界', 'bound', '.8g'),
                             ('有效割数量', 'cuts', 'd'), ('射线最大容量高估 / kW', 'max_radial_difference', '.6f'),
                             ('射线平均容量高估 / kW', 'mean_radial_difference', '.6f')]:
        lines.append(f"| {label} | {summary['bilinear'][key]:{fmt}} | {summary['partition'][key]:{fmt}} |")
    lines += [f"| 达到 ε 停止证书 | {left['certified']} | {right['certified']} |", '',
              f"以基线最终最大射线误差 {accuracy_match['target_max_radial_difference']:.6f} kW 为同一标准，"
              f"分块法第 {accuracy_match['cut_count']} 条割首次达标，对应原运行 {accuracy_match['seconds']:.3f} 秒；"
              f"基线首次达标为第 {accuracy_match['baseline_cut_count']} 条割、"
              f"{accuracy_match['baseline_seconds']:.3f} 秒，首次达到同一精度的用时比约 {accuracy_match['speed_ratio']:.2f}。"
              '前缀比较容差为 0.00001 kW。'
              '这个比值仅针对有限射线扫描，不代表连续域认证速度。', '',
              '![比较曲线](comparison.png)', '',
              'a，所有实际迭代的全局上界，虚线为停止阈值；η 不是 kW 或几何距离。'
              'b，累计有效割数量。c，等时限结束后全部 861 条相同负荷比例射线的径向容量差经验分布；'
              '参考为同一 SOCP 模型下放开全部建设和运行变量求得的独立边界。'
              '数据未筛选；微小负差保留为求解器舍入误差。有限射线不能替代连续域停止证书，也不是 AC 认证。', '',
              '角点分块用 3 个二元变量选择 8 个角点，连续负荷乘积已被消除。'
              '原生 MISOCP 在预试验中出现数值终止，正式分块运行预先选择锥线性外逼近，逐轮补充有效锥切面；'
              '候选乘子使用完整锥和驻点残差的保守下界修正，实际外域见证和 SP 保留原有物理容差。'
              '锥切面与用于切割规划外域的 SP 联合割分别计数。', '',
              '全部新增 SP 联合割均另行在无预算限制的完整 SOCP 可行域上审核；'
              '逐候选重新求原始 min η，检查候选下界和块上界。详见 audit.json。'
              '主线及原双线性实验源码未修改。']
    (output/'report.md').write_text('\n'.join(lines)+'\n', encoding='utf-8')
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', default='results/fourbus_outer_partition/comparison_300s')
    parser.add_argument('--reference', default='results/fourbus_outer/scan/rays.json')
    args = parser.parse_args()
    make_report(args.directory, args.reference)
