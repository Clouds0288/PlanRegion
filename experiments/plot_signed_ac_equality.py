"""绘制已保存的 AC 等式试验；不重新求解。

python experiments/plot_signed_ac_equality.py

绘图契约：a/b 分别显示 FourBus/Case33 的反向功率极限收缩；采用
quantitative grid。每模型保留全部 8 次确定性方向求解及全部定点检验，
没有抽样、插值或统计置信区间。虚线为 direction@p<=bound 的交集，
不是支撑点凸包，也不表示内部已经全部认证。仅 Case33 显示单位由
kW 转为 MW；原始 JSON 不变。Python 绘图，183 mm 宽，文字至少
8 pt，导出可编辑 SVG/PDF 和 600 dpi PNG。
"""
import argparse
from itertools import combinations
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
COLORS = {'socp': '#BE8547', 'ac': '#24618C'}


def support_envelope(rows, method, scale):
    """1. 相交八条全局方向上界，仅用于绘制外包络。"""
    direction = np.array([row['direction'] for row in rows])
    bound = np.array([row[method]['bound'] for row in rows])/scale
    vertices = []
    for i, j in combinations(range(len(rows)), 2):
        pair = direction[[i, j]]
        if np.linalg.det(pair) == 0:
            continue  # 平行支撑线没有唯一交点。
        p = np.linalg.solve(pair, bound[[i, j]])
        if np.all(direction@p <= bound+1e-8):
            vertices.append(p)
    vertices = np.array(vertices)
    centered = vertices-vertices.mean(axis=0)
    vertices = vertices[np.argsort(np.arctan2(centered[:, 1], centered[:, 0]))]
    return np.vstack((vertices, vertices[0]))


def plot_signed_ac_equality(output=ROOT/'results/signed_ac_equality'):
    """2. 读取结果；3. 支撑点与外包络；4. 定点结果；5. 导出。"""
    output = Path(output)
    plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial'],
        'font.size': 9, 'axes.labelsize': 9, 'xtick.labelsize': 8,
        'ytick.labelsize': 8, 'legend.fontsize': 8, 'svg.fonttype': 'none',
        'pdf.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False,
        'axes.linewidth': .7, 'legend.frameon': False, 'mathtext.fontset': 'dejavusans'})
    fig, axes = plt.subplots(1, 2, figsize=(7.2047244, 4.1732283))  # 183 × 106 mm
    fig.subplots_adjust(left=.105, right=.975, bottom=.25, top=.91, wspace=.29)
    for ax, case, panel in zip(axes, ('fourbus', 'case33'), ('a', 'b')):
        data = json.loads((output/case/'comparison.json').read_text(encoding='utf-8'))
        scale, unit = (1., 'kW') if case == 'fourbus' else (1000., 'MW')
        for method, marker, size in (('socp', 's', 6.5), ('ac', 'o', 4.1)):
            points = np.array([row[method]['p'] for row in data['directions']])/scale
            envelope = support_envelope(data['directions'], method, scale)
            ax.plot(*envelope.T, color=COLORS[method], lw=1.15, ls=(0, (4, 2)), zorder=2)
            ax.plot(*points.T, ls='none', marker=marker, ms=size, color=COLORS[method],
                    mfc='none' if method == 'socp' else COLORS[method], mew=1., zorder=4)
        for row in data['probes']:
            feasible = row['ac']['status'] == 'optimal'
            ax.plot(*(np.array(row['p'])/scale), marker='D' if feasible else 'x',
                    ms=4.8 if feasible else 5.5, ls='none', mew=1.1,
                    color='#548674' if feasible else '#943F48', zorder=5)
        nodes = data['load_nodes']
        ax.set_xlabel(rf'$p_{{{nodes[0]}}}$ ({unit})')
        ax.set_ylabel(rf'$p_{{{nodes[1]}}}$ ({unit})')
        limits = (-165, 145) if case == 'fourbus' else (-11, 4.5)
        ticks = [-150, -100, -50, 0, 50, 100] if case == 'fourbus' else [-10, -5, 0]
        ax.set(xlim=limits, ylim=limits, xticks=ticks, yticks=ticks)
        ax.set_aspect('equal', adjustable='box')
        ax.axhline(0, color='#D8D8D8', lw=.65, zorder=0)
        ax.axvline(0, color='#D8D8D8', lw=.65, zorder=0)
        ax.tick_params(direction='out', width=.7, length=3)
        ax.text(-.18, 1.045, panel, transform=ax.transAxes, weight='bold', fontsize=10)
        ax.text(0, 1.045, 'FourBus' if case == 'fourbus' else 'Case33', transform=ax.transAxes)
    handles = [
        Line2D([], [], color=COLORS['socp'], ls='--', marker='s', mfc='none', ms=6,
               label='SOCP: support / outer bound'),
        Line2D([], [], color=COLORS['ac'], ls='--', marker='o', ms=4,
               label='AC equality: support / outer bound'),
        Line2D([], [], color='#548674', ls='none', marker='D', ms=4.8,
               label='AC-feasible probe'),
        Line2D([], [], color='#943F48', ls='none', marker='x', ms=5.5,
               label='AC-infeasible probe')]
    fig.legend(handles=handles, loc='lower center', ncol=2, bbox_to_anchor=(.51, .045),
               columnspacing=2.2, handlelength=2.5, labelspacing=.8)
    fig.savefig(output/'comparison.svg', facecolor='white')
    fig.savefig(output/'comparison.pdf', facecolor='white')
    fig.savefig(output/'comparison.png', dpi=600, facecolor='white')
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'results/signed_ac_equality')
    plot_signed_ac_equality(**vars(parser.parse_args()))
