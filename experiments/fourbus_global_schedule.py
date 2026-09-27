"""FourBus: 96/32/16 new-SP intervals, with and without a .01 eta trigger.

python experiments/fourbus_global_schedule.py
All configurations retain the same SOCP model, tau and global coverage proof.
"""
import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np
from threadpoolctl import threadpool_limits

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from experiments.fourbus_vertex_priority import build_vertex_priority_region
from experiments.fourbus_mainline_comparison import region_radial_capacity
from Network.four_bus_five_corridor import FourBus
from plot import json_value
from region import GEOMETRY_TOL


def render_global_schedule_comparison(trials, output):
    """183 mm quantitative grid: runtime, global work, and boundary accuracy."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import MaxNLocator
    plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
                         'font.size': 8, 'axes.labelsize': 8, 'xtick.labelsize': 7,
                         'ytick.labelsize': 7, 'legend.fontsize': 7.5,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.linewidth': .7, 'svg.fonttype': 'none', 'pdf.fonttype': 42})
    fig, axes = plt.subplots(1, 3, figsize=(7.2047244, 2.6))
    fig.subplots_adjust(left=.075, right=.99, bottom=.235, top=.78, wspace=.55)
    colors = ['#355F7D', '#BD743B']
    for policy, color in enumerate(colors):
        positions = np.arange(3)+(-.09 if policy == 0 else .09)
        for ax, metric in zip(axes, ('total_seconds', 'global_search', 'mean_radial_difference_percent')):
            medians = []
            for j, interval in enumerate((16, 32, 96)):
                key = f'n{interval}'+('_eta001' if policy else '')
                values = np.array([item[metric] for item in trials[key]])
                offsets = np.linspace(-.035, .035, len(values))
                ax.scatter(positions[j]+offsets, values, color=color,
                           marker='o' if policy == 0 else 's', s=13, zorder=3)
                medians.append(float(np.median(values)))
            ax.plot(positions, medians, color=color, linewidth=.9,
                    linestyle='-' if policy == 0 else '--')
    for ax, label, ylabel in zip(axes, 'abc', ('Runtime (s)', 'Global searches', 'Mean overestimate (%)')):
        ax.set(xlim=(-.35, 2.35), xticks=[0, 1, 2], xticklabels=['16', '32', '96'],
               xlabel='SP-count threshold', ylabel=ylabel)
        ax.set_ylim(bottom=0.)
        ax.text(-.18, 1.10, label, transform=ax.transAxes, fontweight='bold', fontsize=9)
    axes[1].yaxis.set_major_locator(MaxNLocator(integer=True, nbins=5))
    fig.legend(handles=[Line2D([], [], color=colors[0], marker='o', markersize=3, label='Interval only'),
                        Line2D([], [], color=colors[1], marker='s', markersize=3, linestyle='--',
                               label='Interval + low-violation trigger')],
               loc='upper center', bbox_to_anchor=(.5, 1.), ncol=2, frameon=False)
    fig.savefig(output/'comparison.svg', facecolor='white')
    fig.savefig(output/'comparison.pdf', facecolor='white')
    fig.savefig(output/'comparison.png', dpi=600, facecolor='white')
    plt.close(fig)


def run_global_schedule_comparison(output, reference_file, repeats=3):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    network = FourBus()
    bounds = np.full(len(network.load_nodes), network.power_limit)
    reference = json.loads(Path(reference_file).read_text(encoding='utf-8'))
    assert reference['budget'] == 20000. and reference['method'] == 'socp'
    weights = np.array(reference['weights'])
    reference_total = np.array(reference['reference_total'])
    configurations = [(f'n{interval}'+('_eta001' if threshold is not None else ''), interval, threshold)
                      for threshold in (None, .01) for interval in (96, 32, 16)]
    raw = {key: [] for key, _, _ in configurations}

    # 1. 六组只改变调度；每次从零初始化，按轮换顺序运行，避免并发计时。
    with threadpool_limits(limits=1):
        for repeat in range(repeats):
            offset = (2*repeat) % len(configurations)
            for key, interval, threshold in configurations[offset:]+configurations[:offset]:
                result = build_vertex_priority_region(network, 20000., bounds, threads=4,
                                                       refinement_checks=interval, eta_trigger=threshold)
                raw[key].append(result)
                print(key, repeat+1, result['timing'], result['counts'], flush=True)
                (output/'runs.json').write_text(json.dumps(json_value(raw), indent=2), encoding='utf-8')

    # 2. 所有运行均需连续覆盖证书；参考射线只用于事后误差比较。
    trials = {key: [] for key in raw}
    rows = []
    for key, group in raw.items():
        for repeat, result in enumerate(group, 1):
            assert result['status'] == 'certified'
            assert result['coverage_bound'] is None or result['coverage_bound'] <= GEOMETRY_TOL
            inner_total = region_radial_capacity(result['inner'], weights, bounds)
            outer_total = region_radial_capacity(result['outer'], weights, bounds)
            assert np.max(inner_total-reference_total) <= 1e-4
            assert np.min(outer_total-reference_total) >= -1e-4
            relative = 100.*(outer_total-reference_total)/reference_total
            globals_ = [item for item in result['trace'] if item['action'] in ('global', 'complete')]
            trials[key].append(dict(total_seconds=result['timing']['total_seconds'],
                                    **result['counts'], coverage_bound=result['coverage_bound'],
                                    first_global_sp=globals_[0]['sp'],
                                    trigger_counts={cause: sum(cause in item['trigger'] for item in globals_)
                                                    for cause in ('empty', 'interval', 'eta')},
                                    mean_radial_difference_percent=float(relative.mean()),
                                    max_radial_difference_percent=float(relative.max()),
                                    max_radial_difference=float((outer_total-reference_total).max()),
                                    max_inner_difference_percent=float(np.max(100.*(reference_total-inner_total)/reference_total))))
            rows.extend([key, repeat, *w, reference_total[j], inner_total[j], outer_total[j]]
                        for j, w in enumerate(weights))

    # 3. 保存全部射线及每次运行，摘要使用三次中位数和范围。
    with (output/'rays_comparison.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream)
        writer.writerow(['configuration', 'repeat', 'w1', 'w2', 'w3', 'reference_total', 'inner_total', 'outer_total'])
        writer.writerows(rows)
    (output/'trials.json').write_text(json.dumps(trials, indent=2), encoding='utf-8')
    metrics = ('total_seconds', 'sp', 'cuts', 'global_search', 'first_global_sp',
               'mean_radial_difference_percent', 'max_radial_difference_percent',
               'max_radial_difference', 'max_inner_difference_percent')
    summary = {key: {metric: dict(median=float(np.median([item[metric] for item in group])),
                                 minimum=min(item[metric] for item in group),
                                 maximum=max(item[metric] for item in group))
                     for metric in metrics} for key, group in trials.items()}
    (output/'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')

    # 4. 三分面展示全部运行的耗时、全局搜索次数和平均边界高估。
    render_global_schedule_comparison(trials, output)
    # 5. 明确调度阈值不等于停止阈值；原因可同时成立，不累加冒充调用次数。
    lines = ['# FourBus：全局搜索调度对照', '',
             f'固定 SOCP、预算 20,000 元、tau=0.002、求解器 4 线程、数值库 1 线程；每组独立运行 {repeats} 次。'
             '中位总耗时包含初始化、SP 评分、全局搜索和几何更新；事后射线评价与绘图另计。', '',
             '| 间隔 | eta<0.01 触发 | 耗时/s | SP次数 | 割数 | G次数 | 首次G前SP数 | 平均高估/% | 最大高估/% |',
             '|---:|:---:|---:|---:|---:|---:|---:|---:|---:|']
    for key, interval, threshold in configurations:
        s = summary[key]
        lines.append(f"| {interval} | {'开' if threshold is not None else '关'} | {s['total_seconds']['median']:.4f} | "
                     f"{s['sp']['median']:g} | {s['cuts']['median']:g} | {s['global_search']['median']:g} | "
                     f"{s['first_global_sp']['median']:g} | {s['mean_radial_difference_percent']['median']:.5f} | "
                     f"{s['max_radial_difference_percent']['median']:.5f} |")
    lines += ['', '低违反量规则：本批所有候选先求原始 min eta、收入可行证书；仍未覆盖者中最大 eta 小于 0.01 时，'
              '应用本轮获选割后请求下一轮全局搜索。单个非最大候选的小 eta 不触发；活动见证的支撑 SP 不被打断。', '',
              '0.01 只控制搜索切换。所有配置都必须通过原来 1e-8 的全局几何覆盖上界检查；'
              '物理 SP 可行判据仍为 eta+MaxVio<=1e-8，不因调度阈值而放松。', '',
              '![全局搜索调度对照](comparison.png)', '',
              'a：各次含初始化总耗时；b：全局搜索次数；c：全部 861 条独立 SOCP 参考射线上的平均外包络容量高估。'
              '每组均保留全部重复运行，线连接中位数；重合观测使用水平偏移显示。横轴为离散调度配置，不是线性比例尺。'
              '重复运行测量计时波动，不代表独立网络样本，未做显著性检验。', '',
              '参考允许全部合法网架与运行变量变化；有限射线不替代连续覆盖证明，也不表示 AC 真值。'
              'runs.json 保存内外域、有效割和逐步 trigger 原因；trials.json 保存逐次指标；'
              'rays_comparison.csv 保存全部方向。trigger 原因可以同时成立。']
    (output/'report.md').write_text('\n'.join(lines)+'\n', encoding='utf-8')
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default='results/fourbus_global_schedule')
    parser.add_argument('--reference', default='results/fourbus_outer_partition/epsilon_0p1/reference_rays.json')
    parser.add_argument('--repeats', type=int, default=3)
    args = parser.parse_args()
    run_global_schedule_comparison(args.output, args.reference, args.repeats)
