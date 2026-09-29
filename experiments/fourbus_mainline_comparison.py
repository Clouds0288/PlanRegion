"""Compare the unmodified mainline with corner partitioning on FourBus.

Each method keeps its own stated stopping certificate: tau=.002 vs epsilon=.1.
This is an output/time comparison, not a matched-accuracy speedup benchmark.
Figure contract: region boundary, paired-ray errors, and every measured runtime.
"""
import argparse
import csv
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator, ScalarFormatter
from matplotlib.tri import Triangulation
import numpy as np
from threadpoolctl import threadpool_limits

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from experiments.fourbus_outer_partition import initialize, run_trial
from experiments.fourbus_partition_report import compare_rays
from continuous import build_continuous_region
from model import GridPhysics
from Network.four_bus_five_corridor import FourBus
from plot import json_value
from region import halfspaces


def region_radial_capacity(polytopes, weights, bounds):
    """Intersect each ray with every convex polytope, then take the union endpoint."""
    outer_total = np.zeros(len(weights))
    for poly in polytopes:
        eq = halfspaces(np.array(poly['vertices'])/bounds)
        coefficients = (weights/bounds) @ eq[:, :-1].T
        rhs = -eq[:, -1]
        positive, negative = coefficients > 1e-12, coefficients < -1e-12
        upper = np.full_like(coefficients, np.inf)
        lower = np.full_like(coefficients, -np.inf)
        np.divide(rhs, coefficients, out=upper, where=positive)
        np.divide(rhs, coefficients, out=lower, where=negative)
        upper, lower = upper.min(axis=1), np.maximum(0., lower.max(axis=1))
        valid = np.all((positive | negative) | (rhs >= -1e-10), axis=1)
        valid &= upper >= lower-1e-8
        assert np.all(np.isfinite(upper[valid]))
        outer_total = np.maximum(outer_total, np.where(valid, upper, 0.))
    return outer_total


def render_mainline_comparison(trials, reference, output):
    """183 x 72 mm quantitative triptych; every repeat and every ray retained."""
    plt.rcParams.update({'font.family': 'sans-serif',
                         'font.sans-serif': ['Microsoft YaHei', 'Arial', 'DejaVu Sans'],
                         'font.size': 8, 'axes.labelsize': 8, 'xtick.labelsize': 7.5,
                         'ytick.labelsize': 7.5, 'legend.fontsize': 7.5,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.linewidth': .7, 'svg.fonttype': 'none', 'pdf.fonttype': 42})
    weights = np.array(reference['weights'])
    reference_total = np.array(reference['reference_total'])
    xy = weights @ np.array([[0., 0.], [1., 0.], [.5, np.sqrt(3.)/2.]])
    triangles = Triangulation(xy[:, 0], xy[:, 1]).triangles
    fig = plt.figure(figsize=(7.2047244, 2.8346457))
    ax = fig.add_axes([.0, .09, .385, .81], projection='3d')
    errors = fig.add_axes([.515, .245, .205, .59])
    times = fig.add_axes([.807, .245, .18, .59])
    colors = {'mainline': '#355F7D', 'partition': '#BD743B'}
    names = {'mainline': '主线', 'partition': '角点分块'}
    points = weights*reference_total[:, None]
    ax.plot_trisurf(*points.T, triangles=triangles, color='#BAC1C5', alpha=.75,
                    edgecolor='none', shade=False)
    for method in ('mainline', 'partition'):
        ordered = sorted(trials[method], key=lambda item: item['total_seconds'])
        representative = ordered[len(ordered)//2]
        points = weights*np.array(representative['outer_total'])[:, None]
        ax.plot_trisurf(*points.T, triangles=triangles, color=colors[method],
                        edgecolor=colors[method], linewidth=.12, alpha=.13, shade=False)
        for coordinate in range(3):
            edge = np.flatnonzero(np.isclose(weights[:, coordinate], 0.))
            edge = edge[np.argsort(weights[edge, (coordinate+1) % 3])]
            ax.plot(*points[edge].T, color=colors[method], linewidth=1.,
                    linestyle='-' if method == 'mainline' else '--')
        for item in trials[method]:
            percent = 100.*(np.array(item['outer_total'])-reference_total)/reference_total
            errors.step(np.sort(percent), 100.*np.arange(1, len(percent)+1)/len(percent),
                        where='post', color=colors[method], linewidth=1.2,
                        alpha=.6, linestyle='-' if method == 'mainline' else '--')
        index = ('mainline', 'partition').index(method)
        values = [item['total_seconds'] for item in trials[method]]
        offsets = np.linspace(-.075, .075, len(values))
        times.scatter(index+offsets, values, color=colors[method], s=16, zorder=3)
        times.plot([index-.19, index+.19], [np.median(values)]*2,
                   color=colors[method], linewidth=1.5)
    bounds = np.max(weights*reference_total[:, None], axis=0)
    ax.set(xlim=(0., bounds[0]), ylim=(0., bounds[1]), zlim=(0., bounds[2]),
           xlabel=r'$p_1$ (kW)', ylabel=r'$p_2$ (kW)', zlabel=r'$p_3$ (kW)')
    ax.set_box_aspect(bounds.copy())
    ax.view_init(elev=24., azim=-55.)
    ax.tick_params(pad=0, labelsize=7)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.set_major_locator(MaxNLocator(nbins=3, integer=True))
        axis.labelpad = 0.
        axis.pane.fill = False
        axis.pane.set_edgecolor('white')
        axis._axinfo['grid'].update(color='#DDDDDD', linewidth=.35)
    # The linear section retains zero and tiny negative solver-roundoff values.
    errors.set_xscale('symlog', linthresh=.2)
    errors.set(xlim=(-.001, 30.), ylim=(0., 103.), xlabel='容量高估 (%)', ylabel='累计射线比例 (%)')
    errors.set_xticks([0., .2, 1., 10.])
    errors.xaxis.set_major_formatter(ScalarFormatter())
    times.set(xlim=(-.4, 1.4), ylim=(0., max(item['total_seconds'] for group in trials.values()
                                          for item in group)*1.18), ylabel='总耗时 (s)',
              xticks=[0, 1], xticklabels=['主线', '角点分块'])
    fig.legend(handles=[Line2D([], [], color='#6B757B', label='SOCP 参考'),
                        Line2D([], [], color=colors['mainline'], label='主线外包络'),
                        Line2D([], [], color=colors['partition'], linestyle='--', label='角点分块外域')],
               loc='upper center', bbox_to_anchor=(.5, 1.), ncol=3, frameon=False)
    for label, x in zip('abc', (.025, .47, .76)):
        fig.text(x, .885, label, fontsize=9, fontweight='bold')
    fig.savefig(output/'mainline_vs_partition.svg', facecolor='white')
    fig.savefig(output/'mainline_vs_partition.pdf', facecolor='white')
    fig.savefig(output/'mainline_vs_partition.png', dpi=600, facecolor='white')
    plt.close(fig)


def run_mainline_comparison(output, reference_file, repeats=3):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    reference = json.loads(Path(reference_file).read_text(encoding='utf-8'))
    network = FourBus()
    bounds = np.full(3, network.power_limit)
    weights = np.array(reference['weights'])
    reference_total = np.array(reference['reference_total'])
    assert reference['budget'] == 20000. and reference['method'] == 'socp'
    trials = {'mainline': [], 'partition': []}
    raw = {'mainline': [], 'partition': []}
    # 1. 交替运行原主线与新方法；两者从零开始，预算、线程和数值库线程相同。
    with threadpool_limits(limits=1):
        for repeat in range(repeats):
            events, cuts = [], []

            def record_mainline(event, **data):
                if event == 'cut':
                    cuts.append(np.array(data['cut']).copy())
                item = dict(event=event, counts=data['counts_algorithm'])
                if event == 'residual_end':
                    item.update(answer=data['answer'])
                events.append(json_value(item))
                if event in ('phase_start', 'residual_end', 'region_end'):
                    print('mainline', repeat+1, event, item, flush=True)

            main_result = build_continuous_region(network, 'socp', 20000., bounds,
                                                   tau=.002, residual_mode='auto',
                                                   time_limit=300., threads=4, progress=record_mainline)
            trial_dir = output/f'trial_{repeat+1}'
            trial_dir.mkdir(exist_ok=True)
            (trial_dir/'mainline.json').write_text(json.dumps(json_value(main_result), indent=2), encoding='utf-8')
            (trial_dir/'mainline_events.json').write_text(json.dumps(events, indent=2), encoding='utf-8')
            (trial_dir/'mainline_cuts.json').write_text(json.dumps(json_value(cuts), indent=2), encoding='utf-8')
            raw['mainline'].append(main_result)
            initialization = initialize(20000., 4)
            partition_result = run_trial('partition', initialization, epsilon=.1, threads=4,
                                         seconds=300., time_limit=5., output=trial_dir/'partition')
            raw['partition'].append(partition_result)
            print('TIMING', repeat+1, main_result['timing']['total_seconds'],
                  partition_result['seconds']+initialization['initial_seconds'], flush=True)

    # 2. 使用同一组独立 SOCP 参考射线，审核主线实际返回的内域/外包络和新方法投影。
    equations = GridPhysics(network, 'socp')
    xy = weights @ np.array([[0., 0.], [1., 0.], [.5, np.sqrt(3.)/2.]])
    triangles = Triangulation(xy[:, 0], xy[:, 1]).triangles
    reference_mesh_volume = float(np.abs(np.linalg.det((weights*reference_total[:, None])[triangles])).sum()/6.)
    for method, group in raw.items():
        for index, result in enumerate(group):
            if method == 'mainline':
                assert result['status'] == 'certified'
                inner_total = region_radial_capacity(result['inner'], weights, bounds)
                outer_total = region_radial_capacity(result['outer'], weights, bounds)
                assert np.max(inner_total-reference_total) <= 1e-4
                item = dict(total_seconds=result['timing']['total_seconds'],
                            counts=result['counts'], coverage_bound=result['coverage_bound'],
                            tau=result['tau'], inner_total=inner_total.tolist(),
                            max_inner_difference_percent=float(np.max(100.*(reference_total-inner_total)/reference_total)),
                            mean_inner_difference_percent=float(np.mean(100.*(reference_total-inner_total)/reference_total)))
            else:
                assert result['certified'] and result['bound'] <= .1
                outer_total = np.array(compare_rays(equations, result, reference)['outer_total'])
                item = dict(total_seconds=result['seconds']+result['initial_seconds'],
                            counts=result['counts'], bound=result['bound'], epsilon=result['epsilon'])
            # 3. 同一分母、同一射线与同一三角网格，记录全部重复运行的误差。
            difference = outer_total-reference_total
            assert np.min(difference) >= -1e-4
            relative = 100.*difference/reference_total
            mesh_volume = float(np.abs(np.linalg.det((weights*outer_total[:, None])[triangles])).sum()/6.)
            item.update(outer_total=outer_total.tolist(),
                        max_radial_difference=float(difference.max()), mean_radial_difference=float(difference.mean()),
                        max_radial_difference_percent=float(relative.max()),
                        mean_radial_difference_percent=float(relative.mean()),
                        mesh_volume_difference_percent=100.*(mesh_volume/reference_mesh_volume-1.))
            trials[method].append(item)
            print('RAYS', method, index+1, {k:v for k,v in item.items() if not isinstance(v, list)}, flush=True)
    (output/'trials.json').write_text(json.dumps(trials, indent=2), encoding='utf-8')
    with (output/'rays_comparison.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream)
        writer.writerow(['method', 'repeat', 'w1', 'w2', 'w3', 'reference_total', 'outer_total', 'inner_total'])
        for method, group in trials.items():
            for repeat, item in enumerate(group):
                for index, w in enumerate(weights):
                    writer.writerow([method, repeat+1, *w, reference_total[index], item['outer_total'][index],
                                     item['inner_total'][index] if method == 'mainline' else ''])
    # 4. 耗时与精度报告中位数和范围；不同停止条件不解释为同精度加速比。
    keys = ('total_seconds', 'max_radial_difference_percent', 'mean_radial_difference_percent',
            'max_radial_difference', 'mesh_volume_difference_percent')
    summary = {method: {key: dict(median=float(np.median([item[key] for item in group])),
                                  minimum=min(item[key] for item in group), maximum=max(item[key] for item in group))
                       for key in keys} for method, group in trials.items()}
    (output/'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    # 5. 科研图展示实际代表运行、全部误差曲线和全部运行时间。
    render_mainline_comparison(trials, reference, output)
    lines = ['# 当前主线与角点分块方法：FourBus 同案例比较', '',
             'SOCP、预算 20,000 元、求解器 4 线程、数值库 1 线程；两者均独立初始化，不使用旧割。'
             f'每方法重复 {repeats} 次，耗时含建模、MP2 初始化、求解和几何处理；事后扫描和绘图另计。', '',
             '主线保持 τ=0.002、light 全局剩余域；新方法保持 ε=0.1。停止证书不同，'
             '这是各自当前设置下的产出比较，不是同精度加速比。参考是同一物理模型的 861 条独立 SOCP 射线。', '',
             '| 指标（中位数；最小—最大） | 主线 | 角点分块 |', '|---|---:|---:|']
    for label, key in [('含初始化总耗时 / s', 'total_seconds'),
                       ('平均容量高估 / %', 'mean_radial_difference_percent'),
                       ('最大容量高估 / %', 'max_radial_difference_percent'),
                       ('最大容量高估 / kW', 'max_radial_difference'),
                       ('三角网格体积高估估计 / %', 'mesh_volume_difference_percent')]:
        cells = [f"{summary[method][key]['median']:.5f} ({summary[method][key]['minimum']:.5f}—"
                 f"{summary[method][key]['maximum']:.5f})" for method in trials]
        lines.append(f'| {label} | '+ ' | '.join(cells)+' |')
    lines += ['', '主线认证内域最大径向低估：'+', '.join(f"{item['max_inner_difference_percent']:.5f}%"
               for item in trials['mainline'])+'。新方法没有认证内域，此项不适用。', '',
              '![当前主线与新方法](mainline_vs_partition.png)', '',
              'a：相同视角下的 SOCP 参考、主线最终外包络及新方法外域；每方法使用耗时位于中位数的实际运行。'
              '边界以全部 861 条射线点插值显示，不作为可行证书。b：每次运行的全部射线高估经验分布；'
              '横轴在 0.2% 内线性、之外对数，保留近零数值误差。c：全部运行时间，横线为中位数，未进行显著性检验。', '',
              '主线实际返回的内域和外包络在全部 861 条射线上分别不超过、不低于完整 SOCP 参考，允许声明的数值容差。'
              '新方法外域也通过相同参考检查。射线均值等权，体积为同一三角网格下的估计；不是 AC 真值或精确连续体积误差。', '',
              '全部主线与新方法源码保持原样；本脚本只调用、记录、比较。逐次输出、全局搜索记录、联合割及逐射线结果均保留。']
    (output/'report.md').write_text('\n'.join(lines)+'\n', encoding='utf-8')
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default='results/fourbus_outer_partition/mainline_comparison')
    parser.add_argument('--reference', default='results/fourbus_outer_partition/epsilon_0p1/reference_rays.json')
    parser.add_argument('--repeats', type=int, default=3)
    args = parser.parse_args()
    run_mainline_comparison(args.output, args.reference, args.repeats)
