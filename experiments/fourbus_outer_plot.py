"""只读 FourBus 独立扫描数据，绘制正式分面图与完整迭代轨迹。

python experiments/fourbus_outer_plot.py
参考表面来自有限方向扫描的显示插值，不是跨网架凸包或域认证。
"""
import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors, patches, tri
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection
import numpy as np


BLUE = '#3C7B9B'
ORANGE = '#CC8035'
PURPLE = '#795894'
PALE = '#F1F1EF'
REFERENCE = '#8DB5CA'
MISSED = '#AB4666'


def panel_label(ax, label, is_3d=False):
    if is_3d:
        ax.text2D(-0.12, 1.03, label, transform=ax.transAxes, fontweight='bold', fontsize=9)
    else:
        ax.text(-0.16, 1.03, label, transform=ax.transAxes, fontweight='bold', fontsize=9)


def save_figure(fig, output, name):
    # 保持指定物理尺寸；PDF/SVG 文本可编辑，PNG 仅用于预览。
    fig.savefig(output/f'{name}.svg', facecolor='white')
    fig.savefig(output/f'{name}.pdf', facecolor='white')
    fig.savefig(output/f'{name}.png', dpi=600, facecolor='white')
    plt.close(fig)


def slice_panel(ax, grid, step, level):
    shape = tuple(grid['shape'])
    power = grid['power'].reshape(*shape, 3)
    index = int(np.argmin(abs(power[0, 0, :, 2]-level)))
    reference = grid['reference'].reshape(shape)[:, :, index] == 1
    outer = grid['first_exclusion'].reshape(shape)[:, :, index] > step
    # 0=参考不可行且已排除；1=参考可行且保留；2=参考不可行但保留；3=误排。
    category = np.where(reference, np.where(outer, 1, 3), np.where(outer, 2, 0))
    h = float(grid['grid_step'])
    x_edges = np.arange(shape[0]+1)*h
    y_edges = np.arange(shape[1]+1)*h
    cmap = colors.ListedColormap([PALE, REFERENCE, ORANGE, MISSED])
    ax.pcolormesh(x_edges, y_edges, category.T, cmap=cmap, vmin=-0.5, vmax=3.5,
                  shading='flat', rasterized=False)
    ax.set(xlim=(0, x_edges[-1]), ylim=(0, y_edges[-1]), xlabel='p₁ (kW)', ylabel='p₂ (kW)')
    ax.set_xticks([0, 30, 60, 90])
    ax.set_yticks([0, 10, 20, 30, 40])
    ax.set_aspect('equal', adjustable='box')
    return float(power[0, 0, index, 2])


def region_comparison(grid, rays, summary, output):
    # 图 1：整体边界、同负荷比例容量误差、两个独立负荷切片。
    fig = plt.figure(figsize=(7.0866, 4.9213))  # 180 × 125 mm
    layout = fig.add_gridspec(2, 2, left=.10, right=.93, bottom=.15, top=.96,
                            wspace=.39, hspace=.30, height_ratios=[1.45, .65])
    ax = fig.add_subplot(layout[0, 0], projection='3d')
    weights = np.array(rays['weights'])
    reference = weights*np.array(rays['reference_total'])[:, None]
    outer = weights*np.array(rays['outer_total'])[:, None]
    triangles = tri.Triangulation(weights[:, 0], weights[:, 1]).triangles
    ax.add_collection3d(Poly3DCollection(reference[triangles], facecolors=BLUE,
                                        edgecolors='none', alpha=.40))
    edges = np.unique(np.sort(np.vstack([triangles[:, [0, 1]], triangles[:, [1, 2]],
                                        triangles[:, [2, 0]]]), axis=1), axis=0)
    ax.add_collection3d(Line3DCollection(outer[edges], colors=ORANGE, linewidths=.20, alpha=.65))
    ax.set(xlim=(0, 100), ylim=(0, 42), zlim=(0, 36), xlabel='p₁ (kW)',
           ylabel='p₂ (kW)', zlabel='p₃ (kW)')
    ax.set_xticks([0, 50, 100])
    ax.set_yticks([0, 20, 40])
    ax.set_zticks([0, 15, 30])
    ax.tick_params(pad=0)
    ax.xaxis.labelpad = 2
    ax.yaxis.labelpad = 2
    ax.zaxis.labelpad = 1
    ax.set_box_aspect((1.2, 1., .9))
    ax.view_init(elev=24, azim=37)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.fill = False
        axis._axinfo['grid']['color'] = '#DDDDDD'
        axis._axinfo['grid']['linewidth'] = .35
    panel_label(ax, 'a', is_3d=True)
    ax.legend(handles=[patches.Patch(facecolor=BLUE, alpha=.6, label='SOCP 参考边界'),
                       patches.Patch(facecolor='none', edgecolor=ORANGE, label='81 条割外域边界')],
              loc='upper center', bbox_to_anchor=(.52, 1.09), frameon=False, fontsize=7)

    ax = fig.add_subplot(layout[0, 1])
    difference = np.array(rays['outer_total'])-np.array(rays['reference_total'])
    artist = ax.scatter(weights[:, 0], weights[:, 1], c=difference, s=9, marker='s',
                        cmap='Purples', vmin=0, vmax=max(difference.max(), 1e-8), linewidths=0)
    ax.plot([0, 1, 0, 0], [0, 0, 1, 0], color='#555555', lw=.6)
    ax.set(xlabel='负荷比例 w₁', ylabel='负荷比例 w₂', xlim=(-.03, 1.03), ylim=(-.03, 1.03))
    ax.set_xticks([0, .5, 1])
    ax.set_yticks([0, .5, 1])
    ax.set_aspect('equal')
    bar = fig.colorbar(artist, ax=ax, fraction=.045, pad=.035)
    bar.set_label('总负荷边界差 Δρ (kW)', labelpad=5)
    panel_label(ax, 'b')

    for label, slot, level in [('c', layout[1, 0], 9.), ('d', layout[1, 1], 25.)]:
        ax = fig.add_subplot(slot)
        actual = slice_panel(ax, grid, summary['cut_count'], level)
        ax.set_title(f'p₃ = {actual:g} kW', pad=7, fontsize=8)
        panel_label(ax, label)
    handles = [patches.Patch(color=REFERENCE, label='扫描可行且保留'),
               patches.Patch(color=ORANGE, label='扫描不可行但外域保留'),
               patches.Patch(color=PALE, label='扫描不可行且已排除')]
    if summary['false_exclusion']:
        handles.append(patches.Patch(color=MISSED, label='扫描可行却被排除'))
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.5, .015),
               ncol=len(handles), frameon=False, fontsize=7, handlelength=1.)
    save_figure(fig, output, 'region_comparison')


def iteration_history(grid, iterations, summary, output):
    # 图 2：全部轮次均显示，不删去上界停滞或反复访问旧方案的轮次。
    fig, axes = plt.subplots(2, 2, figsize=(7.0866, 5.3543))  # 180 × 136 mm
    fig.subplots_adjust(left=.105, right=.97, bottom=.10, top=.94, wspace=.39, hspace=.43)
    step = np.array([int(r['step']) for r in iterations])
    eta = np.array([float(r['eta']) for r in iterations])
    bound = np.array([float(r['effective_bound']) for r in iterations])
    assert np.all(eta > 0) and np.all(bound > 0), 'Log axis requires positive measured values'
    ax = axes[0, 0]
    ax.semilogy(step, bound, color=PURPLE, lw=1.4, label='历史有效 G 上界')
    ax.semilogy(step, eta, color=BLUE, lw=.6, marker='o', ms=2, label='候选 SP 的最小 η')
    ax.axhline(summary['epsilon'], color='#555555', ls='--', lw=.8, label='停止阈值 ε')
    ax.set(ylabel='SP 残差尺度', xlabel='迭代 k', ylim=(min(eta.min(), summary['epsilon'])/1.8, .5))
    ax.set_yticks([1e-4, 1e-3, 1e-2, 1e-1], ['10⁻⁴', '10⁻³', '10⁻²', '10⁻¹'])
    ax.legend(frameon=False, fontsize=7, loc='upper right')

    ax = axes[0, 1]
    all_steps = np.arange(summary['cut_count']+1)
    extra = np.array([np.sum((grid['first_exclusion'] > k) & (grid['reference'] == -1))
                      for k in all_steps])
    assert np.all(np.diff(extra) <= 0), 'Nested cuts cannot expand grid outer region'
    ax.step(all_steps, extra, where='post', color=ORANGE, lw=1.4)
    ax.set(xlabel='已加入的割数 k', ylabel='外域多保留的不可行扫描点数')
    ax.set_ylim(bottom=0)

    ax = axes[1, 0]
    scheme = np.array([int(r['scheme'][1:]) for r in iterations])
    globally_feasible = np.array([r['p_feasible'] == 'True' for r in iterations])
    ax.scatter(step[globally_feasible], scheme[globally_feasible], c=BLUE, s=10,
               label='负荷 p 在其他网架可行', linewidths=0)
    ax.scatter(step[~globally_feasible], scheme[~globally_feasible], c=ORANGE, s=12,
               marker='x', linewidths=.7, label='负荷 p 在全部网架不可行')
    ax.set(xlabel='迭代 k', ylabel='G 选中的网架')
    labels = np.arange(1, summary['visited_schemes']+1, 2)
    ax.set_yticks(labels, [f'S{i:02d}' for i in labels])
    ax.set_ylim(.3, summary['visited_schemes']+.7)
    ax.legend(frameon=False, fontsize=7, loc='upper left', bbox_to_anchor=(0, 1.20))

    ax = axes[1, 1]
    ax.bar(step, [float(r['seconds']) for r in iterations], color='#A6A6A2', width=.85)
    ax.set(xlabel='迭代 k', ylabel='单轮 G 求解时间 (s)', ylim=(0, 11))
    for index, ax in enumerate(axes.flat):
        ax.axvline(summary['restart_after']+.5, color='#999999', lw=.7, ls=':')
        ax.set_xlim(0, summary['cut_count']+1)
        ax.set_xticks([0, 20, 40, 60, 80])
        panel_label(ax, chr(ord('a')+index))
        ax.spines[['top', 'right']].set_visible(False)
    save_figure(fig, output, 'iteration_history')


def region_evolution(grid, summary, output):
    # 图 3：相同 p3 截面、同一参考数据，显示联合割投影的逐步收缩。
    stages = [0, 1, 5, 10, 20, 40, summary['restart_after'], summary['cut_count']]
    fig, axes = plt.subplots(2, 4, figsize=(7.0866, 2.7))  # 180 × 68.6 mm
    fig.subplots_adjust(left=.075, right=.99, top=.91, bottom=.22, wspace=.31, hspace=.55)
    for index, (ax, step) in enumerate(zip(axes.flat, stages)):
        slice_panel(ax, grid, step, 9.)
        ax.set_title(f'k = {step}', fontsize=8, pad=6)
        ax.set_xticks([0, 40, 80])
        ax.set_yticks([0, 20, 40])
        if index < 4:
            ax.set_xlabel('')
        if index % 4:
            ax.set_ylabel('')
            ax.tick_params(labelleft=False)
        panel_label(ax, chr(ord('a')+index))
    fig.legend(handles=[patches.Patch(color=REFERENCE, label='扫描可行且保留'),
                        patches.Patch(color=ORANGE, label='扫描不可行但外域保留'),
                        patches.Patch(color=PALE, label='扫描不可行且已排除')],
               loc='lower center', bbox_to_anchor=(.5, .025), ncol=3, frameon=False,
               fontsize=7, handlelength=1.)
    save_figure(fig, output, 'region_evolution')


def render_scan(output):
    output = Path(output)
    summary = json.loads((output/'summary.json').read_text(encoding='utf-8'))
    rays = json.loads((output/'rays.json').read_text(encoding='utf-8'))
    with np.load(output/'grid.npz') as source:
        grid = {key: source[key] for key in source.files}
    with (output/'iterations.csv').open(encoding='utf-8-sig', newline='') as stream:
        iterations = list(csv.DictReader(stream))
    with plt.rc_context({'font.family': ['Arial', 'Microsoft YaHei', 'DejaVu Sans'],
                         'font.size': 8, 'axes.labelsize': 8, 'axes.titlesize': 8,
                         'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 7,
                         'axes.linewidth': .6, 'xtick.major.width': .6, 'ytick.major.width': .6,
                         'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
                         'svg.fonttype': 'none', 'pdf.fonttype': 42,
                         'axes.unicode_minus': False, 'savefig.dpi': 300}):
        region_comparison(grid, rays, summary, output)
        iteration_history(grid, iterations, summary, output)
        region_evolution(grid, summary, output)
    captions = '''# 图注与数据来源

**图 1｜最终外域与独立 SOCP 扫描。** a，负荷三维空间中的参考边界与第 81 条割后的外域边界；两者在相同的 861 条负荷比例射线上分别求解，三角面仅为显示插值，不主张射线之间的面经过认证，也不跨方案取凸包。b，相同负荷比例下外域总负荷容量减参考容量，w₃=1−w₁−w₂。c、d，p₃=9、25 kW 截面的体素中心分类，格宽 2 kW。蓝色点通过完整 SOCP 验证，橙色点未通过但仍在外域，浅灰点未通过且已排除；色块表示采样位置，不证明整个色块的状态。参考计算不使用任何实验割，建设变量 x 和运行变量 y 均自由，预算统一为 20,000 元。

**图 2｜全部 81 轮实际迭代。** a，历史有效 G 上界、事后复核的候选最小 SP 违反量与 ε=10⁻⁴；η 是原始 SP 残差尺度，不是 kW。b，加入各轮割后，三维网格中仍被外域保留的参考不可行点数。c，各轮 G 实际选中的网架及候选负荷的全局可行性；每轮该固定网架与负荷组合均被 SP 判为不可行，蓝点表示换一个网架可行。此处负荷的全局可行性是事后独立复核，不是原算法决策步骤。d，单轮 G 实际求解时间。竖虚线表示 61 轮后续算；上界曲线在后处理中保留了重启前的有效界，续算原始上界保存在 CSV。没有省略上界停滞或反复访问旧网架的轮次。

**图 3｜联合割在负荷空间投影后的收缩。** p₃=9 kW 的同一截面，依次显示加入 0、1、5、10、20、40、61、81 条割后的状态。颜色与图 1 相同。每个橙色或浅灰扫描点都在所有合法网架中查询投影归属，未固定为当轮访问的网架。因此一条割排除某个 (x,p) 组合时，负荷点本身可能仍被其他网架保留。

数据源：grid.npz（18,144 点）、rays.json（861 对边界优化）、iterations.csv（81 轮）、summary.json（方案表和统计）。所有数值是确定性求解或有限扫描，非重复随机实验，不使用统计置信区间。有限扫描与 SOCP 松弛均不能替代连续域证书或 AC 校验。
'''
    (output/'figure_captions.md').write_text(captions, encoding='utf-8')
    print(f'Figures saved to {output.resolve()}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default='results/fourbus_outer/scan')
    render_scan(**vars(parser.parse_args()))
