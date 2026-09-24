"""Plot measured progress only; time-limited runs are not completed runtimes.

Figure contract: fixed-budget progress and unresolved radial gap answer whether
the screening variants warrant a completion-time comparison. Single-panel
quantitative charts, Python, editable SVG/PDF and PNG preview, all observed runs.
No uncertainty bands for a single seed; raw trajectories are never interpolated.
"""
import argparse
import csv
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SKILL = Path.home()/'.codex/skills/nature-figure'
sys.path.insert(0, str(SKILL/'scripts'))
from audit_panel_alignment import require_matplotlib_panel_alignment


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    args = parser.parse_args()
    directory = Path(args.input)
    data = [json.loads(p.read_text(encoding='utf-8')) for p in sorted(directory.glob('*/result.json'))]
    plt.rcParams.update({'font.family': 'Arial', 'font.size': 9, 'pdf.fonttype': 42,
                         'svg.fonttype': 'none', 'axes.spines.top': False,
                         'axes.spines.right': False, 'legend.frameon': False})
    labels = {'A': 'A: Current method', 'A-L': 'A-L: Linear outer model',
              'B': 'B: FIFO rays', 'C': 'C: Largest-gap rays'}
    colors = {'A': '#555555', 'A-L': '#BB8822', 'B': '#4477AA', 'C': '#009988'}
    styles = {'A': '-', 'A-L': '-.', 'B': '-', 'C': '--'}
    source = []
    for metric, name in [('inner_grid_points', 'inner_progress'), ('max_coverage_gap', 'radial_gap')]:
        fig, ax = plt.subplots(figsize=(6.5, 3.6), layout='constrained')
        for result in data:
            variant = result['variant']
            if metric == 'max_coverage_gap' and variant not in ('B', 'C'):
                continue
            records = result['grid_progress'] if metric == 'inner_grid_points' else result['history']
            times = [r.get('elapsed', r.get('seconds')) for r in records]
            values = [r[metric] for r in records]
            ax.step(times, values, where='post', color=colors[variant], linestyle=styles[variant],
                    linewidth=2 if variant == 'B' else 1.5, label=labels[variant])
            if metric == 'inner_grid_points':
                ax.scatter(times, values, s=18, color=colors[variant], marker='x' if variant == 'C' else 'o')
            source.extend(dict(variant=variant, seed=result['seed'], metric=metric, seconds=t, value=v)
                          for t, v in zip(times, values))
        ax.set_xlabel('Elapsed wall time (s)')
        ax.set_xlim(0, 1.025*max(r['time_limit'] for r in data))
        if metric == 'max_coverage_gap':
            ax.set_ylabel('Maximum radial gap')
            ax.set_yscale('log')
            ax.set_ylim(.001, 1.4)
            ax.axhline(.002, color='#777777', linestyle=':', linewidth=1., label='Target: 0.002')
            ax.legend(loc='center right')
        else:
            ax.set_ylabel('Grid points in the certified inner region')
            ax.set_ylim(bottom=-100)
            ax.legend(loc='best')
        fig.canvas.draw()
        require_matplotlib_panel_alignment(fig, json_out=str(directory/f'{name}.alignment.json'), strict=True)
        fig.savefig(directory/f'{name}.pdf')
        fig.savefig(directory/f'{name}.svg')
        fig.savefig(directory/f'{name}.png', dpi=300)
        plt.close(fig)
    with (directory/'figure_source.csv').open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=['variant', 'seed', 'metric', 'seconds', 'value'])
        writer.writeheader()
        writer.writerows(source)
    (directory/'FIGURE_NOTES.md').write_text(
        '# Figure interpretation\n\n'
        'All screening runs are shown (budget 2, seed 0, one run per variant). '
        'Markers in inner_progress are the recorded 30/60/120-second snapshots; '
        'steps hold the last recorded value and are not estimates of intermediate progress. '
        'The grid has 32³ centers in the common evaluation box. Counts measure certified '
        'inner membership, not coverage of the unknown true domain. '
        'radial_gap shows raw B/C records only: the original method uses a different '
        'coverage objective. The dotted target ignores the additional 1e-8 coordinate '
        'tolerance used by the actual all-vertex certificate. No completed-runtime or '
        'statistical-significance claim is inferred from these screening charts. '
        'Python/matplotlib; single-panel quantitative figures; no omitted runs.\n', encoding='utf-8')


if __name__ == '__main__':
    main()
