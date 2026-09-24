"""Render the measured two-stage experiment, retaining every recorded cut.

Figure contract: show which certified exclusions each stage contributes and
compare the final outer approximation with an independent 5 kW grid reference.
Python quantitative panels; no invented boundary data or curve smoothing.
"""
import argparse
import csv
import io
import json
from pathlib import Path
import sys

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, PathPatch
from matplotlib.path import Path as ArtistPath
import numpy as np
from shapely.geometry import Polygon, LineString, box, shape
from shapely.geometry.polygon import orient
from shapely.ops import unary_union

from tests.case33_two_stage import Interval, region_from_intervals

QA_SCRIPTS = Path.home()/'.codex/skills/nature-figure/scripts'
sys.path.insert(0, str(QA_SCRIPTS))
from audit_panel_alignment import require_matplotlib_panel_alignment

mpl.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 8, 'axes.labelsize': 9, 'axes.titlesize': 10, 'xtick.labelsize': 8,
    'ytick.labelsize': 8, 'legend.fontsize': 8, 'pdf.fonttype': 42, 'svg.fonttype': 'none',
    'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': .7,
    'hatch.linewidth': .5, 'savefig.facecolor': 'white'})
RED, BLUE, GREEN = '#B8323D', '#2363B0', '#138B52'


def polygons(geometry):
    if geometry.is_empty:
        return []
    return [geometry] if geometry.geom_type == 'Polygon' else [p for p in geometry.geoms if p.geom_type == 'Polygon']


def fill(ax, geometry, **kwargs):
    for polygon in polygons(geometry):
        polygon = orient(polygon, sign=1.)
        vertices, codes = [], []
        for ring in [polygon.exterior, *polygon.interiors]:
            points = np.asarray(ring.coords)
            vertices.extend(points)
            codes.extend([ArtistPath.MOVETO]+[ArtistPath.LINETO]*(len(points)-2)+[ArtistPath.CLOSEPOLY])
        ax.add_patch(PathPatch(ArtistPath(vertices, codes), **kwargs))


def boundary(ax, geometry, **kwargs):
    for polygon in polygons(geometry):
        points = np.asarray(polygon.exterior.coords)
        ax.plot(points[:, 0], points[:, 1], **kwargs)


def draw_line(ax, geometry, **kwargs):
    if geometry.is_empty:
        return
    if geometry.geom_type in ('LineString', 'LinearRing'):
        points = np.asarray(geometry.coords)
        ax.plot(points[:, 0], points[:, 1], **kwargs)
    elif hasattr(geometry, 'geoms'):
        for part in geometry.geoms:
            draw_line(ax, part, **kwargs)


def diagonal_fill(ax, geometry, limits):
    """Explicit clipped strokes; PDF pattern tiles otherwise confuse geometry QA."""
    fill(ax, geometry, facecolor='white', edgecolor='none', zorder=1)
    extent = ax.get_window_extent()
    slope = limits[1]/limits[0]*extent.width/extent.height
    step = limits[1]*7./72.*ax.figure.dpi/extent.height
    for intercept in np.arange(-slope*limits[0], limits[1]+step, step):
        line = LineString([(0., intercept), (limits[0], intercept+slope*limits[0])])
        draw_line(ax, line.intersection(geometry), color=RED, linewidth=.5, zorder=1)


def support_line(cut, square_kw):
    wx, wy = cut['weights']
    bound = cut['bound']
    if abs(wy) < 1e-12:
        line = LineString([(bound/wx, 0.), (bound/wx, square_kw)])
    else:
        line = LineString([(0., bound/wy), (square_kw, (bound-wx*square_kw)/wy)])
    return line.intersection(box(0., 0., square_kw, square_kw))


def reference_region(grid):
    coordinates, states = grid['coordinates'], grid['states']
    tops = [(coordinates[i], coordinates[np.where(row == 1)[0][-1]])
            for i, row in enumerate(states) if np.any(row == 1)]
    return unary_union([box(0., 0., first, second) for first, second in tops if first > 0. and second > 0.]), tops


def snapshot_region(snapshot, stage1):
    if snapshot['stage'] < 2:
        return Polygon(snapshot['vertices'])
    return region_from_intervals(stage1, [Interval(**n) for n in snapshot['intervals']])


def render_panel(ax, result, current, actual, *, zoom=False, upto=None, show_reference=True, all_cuts=True, highlight=True):
    metadata = result['metadata']
    square_kw = metadata['square_kw']
    square = box(0., 0., square_kw, square_kw)
    stage1 = Polygon(result['stage1_outer'])
    limits = (max(stage1.bounds[2]*1.07, 1.), max(stage1.bounds[3]*1.06, 1.)) if zoom else (square_kw, square_kw)
    visible = box(0., 0., *limits)
    at_stage2 = current['stage'] == 2
    outer = snapshot_region(current, stage1)
    convex = stage1 if at_stage2 else outer
    fill(ax, square.difference(convex).intersection(visible), facecolor='#F4C7C9', edgecolor='none', zorder=0)
    draw_line(ax, square.boundary.intersection(visible), color='#555555', linewidth=.7, zorder=1)
    if at_stage2:
        diagonal_fill(ax, stage1.difference(outer), limits)
    if current['stage'] > 0:
        shown = result['supports'] if at_stage2 else result['supports'][:current['step']]
        if all_cuts:
            for cut in shown:
                draw_line(ax, support_line(cut, square_kw).intersection(visible), color='#828282', linewidth=.48,
                          alpha=.65, zorder=2)
        if highlight and not at_stage2 and 'cut' in current:
            draw_line(ax, support_line(current['cut'], square_kw).intersection(visible), color='#855FA8', linewidth=1., zorder=3)
        boundary(ax, convex, color=BLUE, linewidth=1.35, linestyle=(0, (6, 2.8)), zorder=4)
    if at_stage2:
        boundary(ax, outer, color=GREEN, linewidth=1.5, linestyle=(2, (4.5, 2.5)), zorder=5)
        if highlight and 'cut' in current:
            cut = current['cut']
            threshold, upper = cut['threshold'], max(0., cut['bound'])
            line = LineString([(threshold, square_kw), (threshold, upper), (square_kw, upper)])
            draw_line(ax, line.intersection(stage1), color='#855FA8', linewidth=.8, zorder=3)
    if show_reference and actual is not None:
        boundary(ax, actual, color=RED, linewidth=1.25, linestyle=(0, (2.2, 2.0)), zorder=6)
    if zoom:
        ax.set_xlim(0., limits[0])
        ax.set_ylim(0., limits[1])
    else:
        ax.set_xlim(0., square_kw)
        ax.set_ylim(0., square_kw)
        ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel('Node 18 load (kW)')
    ax.set_ylabel('Node 25 load (kW)')
    ax.tick_params(direction='out', length=3, width=.6)


def handles(spacing, highlight=False):
    hatch_marker = ArtistPath([(-1., -.6), (-.4, .6), (-.3, -.6), (.3, .6), (.4, -.6), (1., .6)],
                              [ArtistPath.MOVETO, ArtistPath.LINETO]*3)
    entries = [Patch(facecolor='#F4C7C9', edgecolor='none', label='Removed by global supports'),
        Line2D([], [], color=RED, ls='none', marker=hatch_marker, markersize=12.,
               markeredgewidth=.5, label='Removed by stage 2'),
        Line2D([], [], color=BLUE, ls=(0, (6, 2.8)), lw=1.4, label='Stage 1 outer envelope'),
        Line2D([], [], color=GREEN, ls=(2, (4.5, 2.5)), lw=1.5, label='Stage 2 final envelope'),
        Line2D([], [], color=RED, ls=(0, (2.2, 2)), lw=1.3, label=f'Actual region ({spacing:g} kW grid)'),
        Line2D([], [], color='#828282', lw=.6, label='Global support lines')]
    if highlight:
        entries.append(Line2D([], [], color='#855FA8', lw=1., label='Current cut'))
    return entries


def export(fig, output, name):
    fig.canvas.draw()
    require_matplotlib_panel_alignment(fig, json_out=str(output/f'{name}.alignment.json'),
        tolerance_pt=1.5, gutter_tolerance_pt=1.5, strict=True)
    fig.savefig(output/f'{name}.png', dpi=600)
    fig.savefig(output/f'{name}.pdf')
    fig.savefig(output/f'{name}.svg')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', type=Path)
    args = parser.parse_args()
    result = json.loads((args.directory/'run/result.json').read_text(encoding='utf-8'))
    scan = json.loads((args.directory/'reference/scan.json').read_text(encoding='utf-8'))
    grid = np.load(args.directory/'reference/scan_grid.npz')
    output = args.directory/'figures'
    output.mkdir(exist_ok=True)
    actual, tops = reference_region(grid)
    snapshots = result['snapshots']
    final = snapshots[-1]
    spacing = scan['spacing_kw']
    for name, zoom in (('two_stage_overview', False), ('two_stage_boundary', True)):
        fig, ax = plt.subplots(figsize=(7.2, 7.0 if not zoom else 5.8))
        fig.subplots_adjust(left=.13, right=.97, bottom=.235, top=.90)
        render_panel(ax, result, final, actual, zoom=zoom, highlight=False)
        ax.set_title('Case33bw: two stages of certified exclusion', pad=12)
        fig.legend(handles=handles(spacing), loc='lower center', bbox_to_anchor=(.5, .035),
                   ncol=2, frameon=False, columnspacing=1.5, handlelength=2.6, labelspacing=.9)
        fig.text(.5, .015, f'Budget 2 | Node 33 = 0 kW | SOCP L-infinity gap <= {result["max_gap_kw"]:.2f} kW',
                 ha='center', fontsize=8)
        export(fig, output, name)

    first = [s for s in snapshots if s['stage'] == 1]
    second = [s for s in snapshots if s['stage'] == 2]
    chosen = [snapshots[0], first[min(5, len(first)-1)], first[-1],
              second[0], second[len(second)//2], second[-1]]
    titles = ['Initial square', 'Global supporting cuts', 'Stage 1 completed',
              'Largest-gap branching (detail)', 'Nonconvex exclusions (detail)', 'Final envelope and reference (detail)']
    fig, axs = plt.subplots(2, 3, figsize=(11.7, 8.3))
    fig.subplots_adjust(left=.065, right=.985, bottom=.16, top=.935, wspace=.33, hspace=.40)
    for index, (ax, snapshot, title) in enumerate(zip(axs.flat, chosen, titles)):
        render_panel(ax, result, snapshot, actual, show_reference=index == 5, zoom=index >= 3)
        if index >= 3:
            ax.set_box_aspect(1.)
        ax.set_title(title, pad=9, fontsize=9)
        ax.annotate(chr(97+index), xy=(0., 1.), xycoords='axes fraction', xytext=(-17, 9), textcoords='offset points',
                    fontweight='bold', fontsize=10, ha='left', va='bottom', annotation_clip=False)
    fig.legend(handles=handles(spacing, highlight=True), loc='lower center', bbox_to_anchor=(.5, .012),
               ncol=3, frameon=False, handlelength=2.6, columnspacing=2., labelspacing=1.)
    export(fig, output, 'cutting_sequence')

    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    fig.subplots_adjust(left=.12, right=.97, bottom=.18, top=.89)
    history = result['history']
    ax.step([0]+[h['step']+1 for h in history],
            [history[0]['selected_gap_kw']]+[h['max_gap_kw'] for h in history],
            where='post', color=GREEN, lw=1.5)
    ax.axhline(result['metadata']['epsilon_kw'], color='#777777', ls='--', lw=.8)
    ax.set(xlabel='Stage 2 bound query', ylabel='Certified SOCP L-infinity gap (kW)',
           title='Global geometric certificate during branching')
    ax.set_yscale('log')
    export(fig, output, 'certified_gap')

    with (output/'support_lines.csv').open('w', newline='', encoding='utf-8-sig') as stream:
        writer = csv.writer(stream)
        writer.writerow(['Cut', 'Weight_Node18', 'Weight_Node25', 'Global_upper_kW', 'Feasible_objective_kW'])
        for i, cut in enumerate(result['supports'], 1):
            writer.writerow([f'S{i}', *cut['weights'], cut['bound'], cut.get('objective')])
    with (output/'branch_cuts.csv').open('w', newline='', encoding='utf-8-sig') as stream:
        writer = csv.writer(stream)
        writer.writerow(['Cut', 'Node18_at_least_kW', 'Node25_upper_kW', 'Globally_infeasible_strip'])
        for i, cut in enumerate(result['logical_cuts'], 1):
            writer.writerow([f'B{i}', cut['threshold'], cut['bound'], cut['infeasible']])
    with (output/'actual_boundary.csv').open('w', newline='', encoding='utf-8-sig') as stream:
        writer = csv.writer(stream)
        writer.writerow(['Node18_kW', 'Highest_AC_feasible_Node25_kW'])
        writer.writerows(tops)

    # Every step has a vector frame. The HTML embeds them, requiring no server
    # or remote assets and preserving all recorded intermediate exclusions.
    frames = []
    (output/'frames').mkdir(exist_ok=True)
    for frame_index, snapshot in enumerate(snapshots):
        fig, ax = plt.subplots(figsize=(7.2, 5.8))
        fig.subplots_adjust(left=.13, right=.97, bottom=.235, top=.90)
        render_panel(ax, result, snapshot, actual, zoom=snapshot['stage'] == 2,
                     show_reference=snapshot is snapshots[-1])
        title = f'Stage {snapshot["stage"]} | Step {snapshot["step"]}'
        if 'cut' in snapshot:
            cut = snapshot['cut']
            if snapshot['stage'] == 1:
                title += f'\nC{snapshot["step"]}: {cut["weights"][0]:.4f} p18 + {cut["weights"][1]:.4f} p25 <= {cut["bound"]:.2f} kW'
            else:
                title += f'\nB{snapshot["step"]+1}: p18 >= {cut["threshold"]:.2f} kW implies p25 <= {cut["bound"]:.2f} kW'
        ax.set_title(title, pad=10, fontsize=9)
        fig.legend(handles=handles(spacing, highlight=True), loc='lower center', bbox_to_anchor=(.5, .02),
                   ncol=2, frameon=False, handlelength=2.6, labelspacing=.9)
        fig.canvas.draw()
        require_matplotlib_panel_alignment(fig,
            json_out=str(output/'frames'/f'step_{frame_index:03d}.alignment.json'), strict=True)
        fig.savefig(output/'frames'/f'step_{frame_index:03d}.pdf')
        stream = io.StringIO()
        fig.savefig(stream, format='svg')
        frames.append(stream.getvalue()[stream.getvalue().index('<svg'):])
        plt.close(fig)
    html = '''<!doctype html><meta charset="utf-8"><title>Case33 two-stage cuts</title>
<style>body{font:16px system-ui;max-width:1050px;margin:24px auto;color:#222}#plot svg{width:100%;height:75vh}input{width:65%}button{padding:8px 15px;margin-right:10px}</style>
<h2>Case33bw: step-by-step certified cuts</h2>
<p>Stage 1 removes the solid red area. Stage 2 removes the red hatched area. The red dashed reference is an independent grid scan; it appears in the final frame.</p>
<button id="play">Play</button><input id="step" type="range" min="0" value="0"><span id="number"></span><div id="plot"></div>
<script>const frames=FRAME_DATA;const slider=document.getElementById('step');slider.max=frames.length-1;
function show(){document.getElementById('plot').innerHTML=frames[+slider.value];document.getElementById('number').textContent=` ${+slider.value+1} / ${frames.length}`;}slider.oninput=show;show();
let timer=null;document.getElementById('play').onclick=()=>{if(timer){clearInterval(timer);timer=null;return;}timer=setInterval(()=>{slider.value=(+slider.value+1)%frames.length;show();},1100);};</script>'''
    (output/'step_by_step.html').write_text(html.replace('FRAME_DATA', json.dumps(frames)), encoding='utf-8')
    (output/'FIGURE_NOTES.md').write_text(
        '# Figure integrity and contract\n\n'
        'The plots show the measured initial square, global supporting halfspaces, and conditional objective-space branch exclusions. '
        'All recorded cuts are retained in the CSV tables and step viewer; the static sequence selects six representative stages, without deleting measurements. '
        'Green is a certified outer approximation, not an exact feasible set. Red is the downward union of independently AC-certified grid points, with no smoothing. '
        'Grid spacing and the final continuous geometric gap are different quantities. No random-replicate confidence interval applies to these deterministic certificates.\n\n'
        'The overview uses equal kW scales and a square initial domain. The boundary detail and the bottom row of the process sheet use separate axis scales to show the gap. '
        'The six-panel process sheet is a report figure, not a journal-width submission asset. '
        'The interactive frames preserve the same plot geometry and styles; their vector export is for inspection.\n', encoding='utf-8')


if __name__ == '__main__':
    main()
