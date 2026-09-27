"""Compare both stopped outer regions with a fresh full-SOCP reference.

Figure contract: two matched quantitative figures, one per method. Panel a
shows every scanned boundary point via a triangular surface; panel b shows
the paired radial percentage error. One run per method; no sampling intervals.
All 861 rays are retained. Surface interpolation is not a feasible certificate.
"""
import argparse
import csv
import json
from pathlib import Path
import sys
from time import perf_counter

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from matplotlib.tri import Triangulation
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from experiments.fourbus_outer_scan import scan_rays
from experiments.fourbus_partition_report import audit_result, compare_rays
from model import GridPhysics
from Network.four_bus_five_corridor import FourBus


def render_threshold_regions(results, reference, rays, summary, output):
    """Matched 183 x 85 mm figures with editable text and all paired observations."""
    plt.rcParams.update({'font.family': 'sans-serif',
                         'font.sans-serif': ['Microsoft YaHei', 'Arial', 'DejaVu Sans'],
                         'font.size': 8, 'axes.labelsize': 8, 'xtick.labelsize': 7,
                         'ytick.labelsize': 7, 'legend.fontsize': 7.5,
                         'axes.linewidth': .7, 'svg.fonttype': 'none', 'pdf.fonttype': 42})
    weights = np.array(reference['weights'])
    reference_total = np.array(reference['reference_total'])
    reference_points = weights*reference_total[:, None]
    height = np.sqrt(3.)/2.
    xy = weights @ np.array([[0., 0.], [1., 0.], [.5, height]])
    triangulation = Triangulation(xy[:, 0], xy[:, 1])
    color_limit = max(1., np.ceil(max(item['max_radial_difference_percent']
                                     for item in summary.values())))
    colors = {'bilinear': '#355F7D', 'partition': '#BD743B'}
    names = {'bilinear': '双线性外域', 'partition': '角点分块外域'}
    axis_bounds = np.array(results['bilinear']['axis_bounds'])
    for method, result in results.items():
        outer_total = np.array(rays[method]['outer_total'])
        outer_points = weights*outer_total[:, None]
        relative = 100.*(outer_total-reference_total)/reference_total
        fig = plt.figure(figsize=(7.2047244, 3.3464567))
        ax = fig.add_axes([.005, .095, .49, .82], projection='3d')
        error_ax = fig.add_axes([.565, .225, .325, .63])
        ax.plot_trisurf(*reference_points.T, triangles=triangulation.triangles,
                        color='#B9C1C6', edgecolor='none', alpha=.76, shade=False)
        ax.plot_trisurf(*outer_points.T, triangles=triangulation.triangles,
                        color=colors[method], edgecolor=colors[method],
                        linewidth=.14, alpha=.15, shade=False)
        # All face-edge rays are emphasized; no observations are omitted.
        for coordinate in range(3):
            edge = np.flatnonzero(np.isclose(weights[:, coordinate], 0.))
            edge = edge[np.argsort(weights[edge, (coordinate+1) % 3])]
            ax.plot(*reference_points[edge].T, color='#505C64', linewidth=1.2)
            ax.plot(*outer_points[edge].T, color=colors[method], linewidth=1.1,
                    linestyle='--')
        ax.set(xlim=(0., axis_bounds[0]), ylim=(0., axis_bounds[1]),
               zlim=(0., axis_bounds[2]), xlabel=r'$p_1$ (kW)', ylabel=r'$p_2$ (kW)',
               zlabel=r'$p_3$ (kW)')
        ax.set_box_aspect(axis_bounds.copy())
        ax.view_init(elev=24., azim=-55.)
        ax.tick_params(pad=0, labelsize=7)
        for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
            axis.set_major_locator(MaxNLocator(nbins=4, integer=True))
            axis.labelpad = 1.
            axis.pane.fill = False
            axis.pane.set_edgecolor('#FFFFFF')
            axis._axinfo['grid'].update(color='#DDDDDD', linewidth=.35)
        colored = error_ax.tripcolor(triangulation, relative, shading='gouraud',
                                     cmap='YlOrBr', norm=Normalize(0., color_limit))
        for fraction in (.2, .4, .6, .8):
            for coordinate in range(3):
                ends = np.zeros((2, 3))
                ends[:, coordinate] = fraction
                ends[0, (coordinate+1) % 3] = 1.-fraction
                ends[1, (coordinate+2) % 3] = 1.-fraction
                segment = ends @ np.array([[0., 0.], [1., 0.], [.5, height]])
                error_ax.plot(*segment.T, color='#7E7465', alpha=.24, linewidth=.4)
        error_ax.plot([0., 1., .5, 0.], [0., 0., height, 0.], color='#666666', linewidth=.7)
        error_ax.text(-.015, -.065, r'$p_1$ 占 100%', ha='left', va='top', fontsize=7.5)
        error_ax.text(1.015, -.065, r'$p_2$ 占 100%', ha='right', va='top', fontsize=7.5)
        error_ax.text(.5, height+.06, r'$p_3$ 占 100%', ha='center', va='bottom', fontsize=7.5)
        error_ax.set(xlim=(-.065, 1.065), ylim=(-.10, height+.115), aspect='equal')
        error_ax.set_axis_off()
        bar_ax = fig.add_axes([.909, .295, .012, .39])
        bar = fig.colorbar(colored, cax=bar_ax)
        bar.set_label('容量高估 (%)', fontsize=8, labelpad=5)
        bar.ax.tick_params(labelsize=7, width=.6, length=2)
        fig.text(.027, .915, 'a', fontsize=9, fontweight='bold')
        fig.text(.545, .915, 'b', fontsize=9, fontweight='bold')
        fig.legend(handles=[Line2D([], [], color='#505C64', linewidth=1.2, label='SOCP 参考'),
                            Line2D([], [], color=colors[method], linewidth=1.2,
                                   linestyle='--', label=names[method])],
                   loc='upper left', bbox_to_anchor=(.055, 1.), ncol=2, frameon=False)
        fig.savefig(output/f'{method}_vs_reference.svg', facecolor='white')
        fig.savefig(output/f'{method}_vs_reference.pdf', facecolor='white')
        fig.savefig(output/f'{method}_vs_reference.png', dpi=600, facecolor='white')
        plt.close(fig)


def make_threshold_report(directory):
    """1 load stopped runs; 2 audit; 3 scan; 4 measure; 5 draw and report."""
    output = Path(directory)
    # 1. 读取实际停止结果；不改写既有实验，不构造认证内域。
    results = {method: json.loads((output/method/'results.json').read_text(encoding='utf-8'))
               for method in ('bilinear', 'partition')}
    left, right = results.values()
    for key in ('budget', 'epsilon', 'threads', 'time_limit', 'initial'):
        assert left[key] == right[key]
    for result in results.values():
        assert result['certified'] and result['bound'] <= result['epsilon']
        assert result['initial_cut_count'] == 0 and 'volume_threshold' not in result
        assert all(row['bound'] > result['epsilon'] for row in result['trace'][:-1])
    equations = GridPhysics(FourBus(), 'socp')
    started = perf_counter()
    # 2. 全部割在完整 SOCP 可行域审查，逐候选重新求解原始 SP。
    audits = {method: audit_result(equations, result) for method, result in results.items()}
    (output/'audit.json').write_text(json.dumps(audits, indent=2), encoding='utf-8')
    # 3. 重新扫描完整模型与双线性外域；新方法使用完全相同的 861 条射线。
    reference = scan_rays(equations, left, 40)
    reference.update(budget=left['budget'], method='socp', divisions=40)
    (output/'reference_rays.json').write_text(json.dumps(reference, indent=2), encoding='utf-8')
    rays = {'bilinear': dict(outer_total=reference['outer_total']),
            'partition': compare_rays(equations, right, reference)}
    weights = np.array(reference['weights'])
    reference_total = np.array(reference['reference_total'])
    assert len(weights) == 861 and np.all(reference_total > 0.)
    assert np.allclose(weights.sum(axis=1), 1.) and np.all(weights >= 0.)
    assert np.max(np.array(reference['reference_bound'])-reference_total) < 1e-3
    xy = weights @ np.array([[0., 0.], [1., 0.], [.5, np.sqrt(3.)/2.]])
    triangles = Triangulation(xy[:, 0], xy[:, 1]).triangles
    reference_points = weights*reference_total[:, None]
    reference_mesh_volume = float(np.abs(np.linalg.det(reference_points[triangles])).sum()/6.)
    summary = {}
    # 4. 容量误差逐射线除以参考值；插值体积单独标注为网格估计。
    for method, result in results.items():
        outer_total = np.array(rays[method]['outer_total'])
        difference = outer_total-reference_total
        radial_difference_percent = 100.*difference/reference_total
        assert difference.min() >= -1e-5
        mesh_volume = float(np.abs(np.linalg.det((weights*outer_total[:, None])[triangles])).sum()/6.)
        summary[method] = dict(seconds=result['seconds'], initial_seconds=result['initial_seconds'],
                               total_seconds=result['seconds']+result['initial_seconds'],
                               bound=result['bound'], epsilon=result['epsilon'],
                               certified=result['certified'], cuts=len(result['cuts']), counts=result['counts'],
                               max_radial_difference=float(difference.max()),
                               mean_radial_difference=float(difference.mean()),
                               min_radial_difference=float(difference.min()),
                               max_radial_difference_percent=float(radial_difference_percent.max()),
                               mean_radial_difference_percent=float(radial_difference_percent.mean()),
                               p95_radial_difference_percent=float(np.percentile(radial_difference_percent, 95.)),
                               mesh_volume=mesh_volume, reference_mesh_volume=reference_mesh_volume,
                               mesh_volume_difference_percent=100.*(mesh_volume/reference_mesh_volume-1.))
        print(method, summary[method], flush=True)
    (output/'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    with (output/'rays_comparison.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream)
        writer.writerow(['w1', 'w2', 'w3', 'reference_total', 'bilinear_outer_total',
                         'partition_outer_total', 'bilinear_difference_percent', 'partition_difference_percent'])
        for index, w in enumerate(weights):
            outer = [rays[method]['outer_total'][index] for method in results]
            relative = [100.*(value-reference_total[index])/reference_total[index] for value in outer]
            writer.writerow([*w, reference_total[index], *outer, *relative])
    with (output/'convergence.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream)
        writer.writerow(['method', 'step', 'seconds', 'bound', 'cuts', 'action'])
        for method, result in results.items():
            for record in result['trace']:
                writer.writerow([method, *[record[key] for key in ('step', 'seconds', 'bound', 'cuts', 'action')]])
    # 5. 分别绘制两个方法；坐标、参考边界、百分比分母和色标完全一致。
    render_threshold_regions(results, reference, rays, summary, output)
    lines = ['# FourBus：全局残差上界 0.1 停止比较', '',
             f"预算 {left['budget']:,.0f} 元，{left['threads']} 线程，单次全局搜索上限 {left['time_limit']:g} 秒。"
             f"两方法均在可靠全局上界首次不超过 {left['epsilon']:g} 时停止。未构造认证内域，也未启用体积停滞停止。", '',
             '| 指标 | 双线性 | 角点分块＋锥外逼近 |', '|---|---:|---:|']
    for label, key, fmt in [('搜索与切割耗时 / s', 'seconds', '.3f'),
                            ('共享初始化 / s', 'initial_seconds', '.3f'),
                            ('含初始化总耗时 / s', 'total_seconds', '.3f'),
                            ('停止时全局上界', 'bound', '.8f'), ('SP 有效割 / 条', 'cuts', 'd'),
                            ('最大径向容量高估 / kW', 'max_radial_difference', '.6f'),
                            ('平均径向容量高估 / kW', 'mean_radial_difference', '.6f'),
                            ('平均径向容量高估 / %', 'mean_radial_difference_percent', '.4f'),
                            ('95 分位径向容量高估 / %', 'p95_radial_difference_percent', '.4f'),
                            ('最大径向容量高估 / %', 'max_radial_difference_percent', '.4f'),
                            ('三角网格体积高估估计 / %', 'mesh_volume_difference_percent', '.4f')]:
        lines.append(f"| {label} | {summary['bilinear'][key]:{fmt}} | {summary['partition'][key]:{fmt}} |")
    lines += ['', '每条射线令 p=ρw、Σw=1，ρ 为总负荷（kW）。相对误差定义为 '
              '100×(ρ_outer−ρ_ref)/ρ_ref；均值为全部 861 条比例网格射线的等权均值。'
              '参考值重新用完整 SOCP 模型独立求解，允许所有预算内 x 和运行 y 变化。', '',
              '![双线性方法与参考值](bilinear_vs_reference.png)', '',
              '![角点分块方法与参考值](partition_vs_reference.png)', '',
              '两图 a：灰色为 SOCP 参考边界，有色网格为停止时外域；采用全部 861 条射线前沿点的三角网格插值。'
              '两图 b：负荷比例三角形上的容量高估百分比，顶点表示单个负荷占 100%；细网格间隔 20%，两图色标一致。'
              '插值面不作为认证可行域，图中未删除或抽稀数据。', '',
              '三角网格体积为各前沿三角形与原点组成的三角锥体积之和；相对参考网格计算体积高估。'
              '这是有限分辨率的体积估计，不是精确连续域体积误差。有限射线也不是 AC 真值或全域边界误差证书。', '',
              '每方法一次实际运行，无重复运行置信区间。算法耗时包含建模、搜索、SP、切割及逐步记录；'
              '共享 MP2 初始化另列，事后审割、扫描与绘图不计入算法耗时。'
              f'本报告事后处理耗时 {perf_counter()-started:.3f} 秒。全部新增割及候选上下界复核通过。', '',
              '复现：', '```powershell',
              'python experiments/fourbus_outer_partition.py --method both --epsilon 0.1 --seconds 300 '
              '--time-limit 5 --threads 4 --output results/fourbus_outer_partition/epsilon_0p1',
              'python experiments/fourbus_threshold_report.py --directory results/fourbus_outer_partition/epsilon_0p1',
              '```']
    (output/'report.md').write_text('\n'.join(lines)+'\n', encoding='utf-8')
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', default='results/fourbus_outer_partition/epsilon_0p1')
    args = parser.parse_args()
    make_threshold_report(args.directory)
