"""Measured 3-D cut replay and independent 10 kW section figures.

Contract: compare budget-dependent certified outer volumes, expose both cut
stages, and distinguish independently scanned sections from 3-D SOCP evidence.
Python quantitative grid and offline Plotly geometry; no interpolated truth.
"""
import argparse
import base64
import csv
import gzip
import json
from pathlib import Path
from itertools import combinations
import sys

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
import numpy as np
from scipy.spatial import ConvexHull
from shapely import constrained_delaunay_triangles, set_precision, union_all, symmetric_difference
from shapely.geometry import Polygon, LineString, MultiPoint, box
from shapely.ops import unary_union
from plotly.offline import get_plotlyjs
from threadpoolctl import threadpool_limits

from region import clip_polytope, polytope_volume
from tests.case33_3d import MASKS, subtract_disjunction
from tests.case33_two_stage import save
from tests.audit_case33_3d import section, actual_region, final_cells
from tests.plot_case33_two_stage import (fill, boundary, diagonal_fill, draw_line,
                                        support_line, handles, RED, BLUE, GREEN)

sys.path.insert(0, str(Path.home()/'.codex/skills/nature-figure/scripts'))
from audit_panel_alignment import require_matplotlib_panel_alignment

mpl.rcParams.update({'font.family':'sans-serif', 'font.sans-serif':['Arial','DejaVu Sans'],
                     'font.size':8, 'pdf.fonttype':42, 'svg.fonttype':'none'})


def export(fig, output, name):
    fig.canvas.draw()
    require_matplotlib_panel_alignment(fig, json_out=str(output/f'{name}.alignment.json'),
                                       tolerance_pt=1.5, gutter_tolerance_pt=1.5, strict=True)
    fig.savefig(output/f'{name}.png', dpi=600)
    fig.savefig(output/f'{name}.pdf')
    fig.savefig(output/f'{name}.svg')
    plt.close(fig)


def polygon_parts(geometry):
    if geometry.is_empty:
        return []
    if geometry.geom_type == 'Polygon':
        return [geometry]
    return [part for part in getattr(geometry, 'geoms', ()) if part.geom_type == 'Polygon']


def clip_view(vertices):
    for axis,limit in enumerate((1550.,6100.,1500.)):
        vertices = clip_polytope(vertices,limit,-np.eye(3)[axis])
    return vertices


def surface(cells, hatching=False, cache=None):
    """Exterior faces of a disjoint convex-cell union; cancel internal faces."""
    groups = {}
    retained = {}
    for vertices in cells:
        vertices = np.asarray(vertices)
        stamp = vertices.tobytes() if cache is not None else None
        faces = cache.get(stamp) if cache is not None else None
        if faces is None:
            if len(vertices) < 4 or np.linalg.matrix_rank(vertices-vertices[0], tol=1e-7) < 3:
                continue
            hull = ConvexHull(vertices)
            planes, faces = {}, []
            for equation in hull.equations:
                rounded = tuple(np.round(equation, 7))
                planes.setdefault(rounded, equation)
            for equation in planes.values():
                normal, offset = equation[:3], equation[3]
                indices = np.where(np.abs(vertices @ normal+offset) < 1e-4)[0]
                if len(indices) < 3:
                    continue
                axis = int(np.argmax(np.abs(normal)))
                sign = 1 if normal[axis] > 0 else -1
                canonical = equation*sign
                key = tuple(np.round(canonical, 6))
                keep = [i for i in range(3) if i != axis]
                polygon = MultiPoint(np.round(vertices[indices][:, keep], 6)).convex_hull
                if polygon.geom_type == 'Polygon':
                    faces.append((key,canonical,axis,keep,sign,polygon))
        if cache is not None:
            retained[stamp] = faces
        for key,canonical,axis,keep,sign,polygon in faces:
            if key not in groups:
                groups[key] = dict(equation=canonical, axis=axis, keep=keep, positive=[], negative=[])
            groups[key]['positive' if sign > 0 else 'negative'].append(polygon)
    if cache is not None:
        # Keep only the current boundary's cells; old cuts do not accumulate RAM.
        cache.clear()
        cache.update(retained)
    vertices, triangles, lines, hatches = [], [], [], []
    for group in groups.values():
        # Fixed-precision overlays prevent microscopic gaps and invalid slivers.
        positive = union_all(group['positive'], grid_size=1e-5)
        negative = union_all(group['negative'], grid_size=1e-5)
        exposed = symmetric_difference(positive, negative, grid_size=1e-5)
        axis, keep, equation = group['axis'], group['keep'], group['equation']

        def lift(points):
            points = np.asarray(points)
            lifted = np.zeros((len(points), 3))
            lifted[:, keep] = points
            lifted[:, axis] = -(points @ equation[keep]+equation[3])/equation[axis]
            return lifted

        for polygon in polygon_parts(exposed):
            polygon = polygon.simplify(1e-6, preserve_topology=True)
            for triangle in constrained_delaunay_triangles(polygon).geoms:
                pts = lift(np.asarray(triangle.exterior.coords)[:3])
                index = len(vertices)
                vertices.extend(pts)
                triangles.append([index, index+1, index+2])
            for ring in [polygon.exterior, *polygon.interiors]:
                lines.extend([*lift(ring.coords), [None, None, None]])
            left, bottom, right, top = polygon.bounds
            span = max(right-left, top-bottom)
            if hatching and span > 1e-6:
                step = max(span/16, 35.)
                for intercept in np.arange(bottom-right, top-left+step, step):
                    segment = LineString([(left, left+intercept), (right, right+intercept)]).intersection(polygon)
                    pieces = [segment] if segment.geom_type == 'LineString' else getattr(segment, 'geoms', ())
                    for piece in pieces:
                        if piece.geom_type == 'LineString' and piece.length > 1e-6:
                            hatches.extend([*lift(piece.coords), [None, None, None]])
    if vertices:
        unique, inverse = np.unique(np.round(vertices, 5), axis=0, return_inverse=True)
        triangles = inverse[np.asarray(triangles)].tolist()
        vertices = unique.tolist()
    return dict(vertices=vertices, triangles=triangles,
                lines=np.asarray(lines, dtype=object).tolist(), hatches=np.asarray(hatches, dtype=object).tolist())


def support_plane(vertices, cut):
    weights = np.asarray(cut['weights'])
    polygon = clip_polytope(vertices, cut['bound'], -weights)
    polygon = clip_polytope(polygon, -cut['bound'], weights)
    if len(polygon) < 3:
        return dict(vertices=[], triangles=[], lines=[], hatches=[])
    center = polygon.mean(axis=0)
    _, _, basis = np.linalg.svd(polygon-center)
    coordinates = (polygon-center) @ basis[:2].T
    order = ConvexHull(coordinates).vertices
    polygon = polygon[order]
    return dict(vertices=polygon.tolist(), triangles=[[0, i, i+1] for i in range(1, len(polygon)-1)],
                lines=np.r_[polygon, polygon[:1]].tolist(), hatches=[])


class GeometryArchive:
    """Share exact rendered vertices and record only changed boundary faces."""

    def __init__(self):
        self.vertices, self.triangles, self.segments = [], [], []
        self.vertex_ids, self.triangle_ids, self.segment_ids = {}, {}, {}
        self.previous = dict(t=set(), l=set())

    def encode(self, geometry, patch=False):
        def vertex(point):
            key = tuple(np.round(point, 5))
            if key not in self.vertex_ids:
                self.vertex_ids[key] = len(self.vertices)
                self.vertices.append(key)
            return self.vertex_ids[key]

        def segment_lines(rows):
            ids, previous = [], None
            for point in rows:
                if point[0] is None:
                    previous = None
                    continue
                current = vertex(point)
                if previous is not None and previous != current:
                    key = tuple(sorted((previous, current)))
                    if key not in self.segment_ids:
                        self.segment_ids[key] = len(self.segments)
                        self.segments.append(key)
                    ids.append(self.segment_ids[key])
                previous = current
            return sorted(set(ids))

        points = [vertex(point) for point in geometry['vertices']]
        faces = []
        for triangle in geometry['triangles']:
            key = tuple(sorted(points[i] for i in triangle))
            if key not in self.triangle_ids:
                self.triangle_ids[key] = len(self.triangles)
                self.triangles.append(key)
            faces.append(self.triangle_ids[key])
        encoded = dict(t=sorted(set(faces)), l=segment_lines(geometry['lines']),
                       h=segment_lines(geometry['hatches']))
        if not patch:
            return encoded
        delta = {}
        for key in ('t', 'l'):
            current = set(encoded[key])
            delta[key+'_add'] = sorted(current-self.previous[key])
            delta[key+'_remove'] = sorted(self.previous[key]-current)
            self.previous[key] = current
        return delta

    def data(self):
        return dict(vertices=self.vertices, triangles=self.triangles, segments=self.segments)


def make_frames(result, output):
    cube = MASKS*result['square_kw']
    current = cube
    archive = GeometryArchive()
    boundary_cache = {}
    frames = [dict(stage=0, label='初始立方体', formula='0 ≤ p18,p25,p33 ≤ 6855 kW',
                   boundary=archive.encode(surface([cube]), patch=True), removed=None, plane=None)]
    removed_support = []
    for i, cut in enumerate(result['supports']):
        plane = support_plane(current, cut)
        discarded = clip_polytope(current, -cut['bound'], np.asarray(cut['weights']))
        if len(discarded):
            removed_support.append(discarded)
        current = clip_polytope(current, cut['bound'], -np.asarray(cut['weights']))
        weights = cut['weights']
        formula = ' + '.join(f'{w:.5g} p{n}' for w, n in zip(weights, (18,25,33)) if abs(w) > 1e-9)
        frames.append(dict(stage=1, label=f'第一阶段 · 支持割 C{i+1}', formula=f'{formula} ≤ {cut["bound"]:.3f} kW',
                           boundary=archive.encode(surface([current]), patch=True),
                           removed=archive.encode(surface([discarded])),
                           removed_zoom=archive.encode(surface([clip_view(discarded)])), plane=archive.encode(plane)))
    outer = [current]
    removed_branch = []
    for i, cut in enumerate(result['branch_cuts']):
        threshold = np.asarray(cut['threshold'])
        normals = np.asarray(cut.get('cut_normals', np.eye(3)))
        following, discarded = [], []
        for vertices in outer:
            if np.any((vertices @ normals.T).max(axis=0) <= threshold+1e-8):
                following.append(vertices)
                continue
            removed = vertices
            for normal, bound in zip(normals, threshold):
                removed = clip_polytope(removed, -bound, normal)
            if len(removed):
                discarded.append(removed)
            following.extend(subtract_disjunction(vertices, threshold, normals))
        outer = following
        removed_branch.extend(discarded)
        formula = ' OR '.join('('+' + '.join(f'{w:.5g} p{n}' for w,n in zip(normal,(18,25,33)) if w > 1e-9)
                             +f' ≤ {t:.3f})' for normal,t in zip(normals,threshold))
        frames.append(dict(stage=2, label=f'第二阶段 · 析取割 B{i+1}', formula=formula+' kW',
                           boundary=archive.encode(surface(outer, cache=boundary_cache), patch=True),
                           removed=archive.encode(surface(discarded, hatching=True)), plane=None))
        if i % 20 == 0:
            print(f'RENDER budget={result["budget"]:g} cut={i+1}/{len(result["branch_cuts"])}', flush=True)
    certified = result['status'] == 'certified'
    last_stage = 3 if result['branch_cuts'] else 1
    final_label = '最终认证外域' if certified else '当前外域（计算中，尚未完成认证）'
    gap_label = (f'全域 SOCP 有向 L∞ 间隙 ≤ {result["max_gap_kw"]:.2f} kW'
                 if result['status'] != 'stage1' else '第一阶段支持平面仍在收紧')
    frames.append(dict(stage=last_stage, label=final_label, formula=gap_label,
                       boundary=archive.encode(surface(outer, cache=boundary_cache), patch=True), removed=None, plane=None))
    payload = dict(budget=result['budget'], epsilon_kw=result['epsilon_kw'], status=result['status'],
                   square_kw=result['square_kw'], frames=frames, stage1=archive.encode(surface([current])),
                   geometry=archive.data(),
                   outer_volume=sum(polytope_volume(v) for v in outer),
                   max_gap_kw=result['max_gap_kw'] if result['status'] != 'stage1' else None, reference=[])
    payload['removed_final'] = [archive.encode(surface(removed_support)),
                                archive.encode(surface(removed_branch, hatching=True))]
    payload['removed_final_zoom'] = [archive.encode(surface([clip_view(v) for v in removed_support])),
                                    payload['removed_final'][1]]
    payload['reference'] = reference_lines(output.parent,result['budget'])
    save(output/f'budget_{result["budget"]:g}_geometry.json', payload)
    return payload


def reference_lines(directory,budget):
    reference = []
    for fixed_value in (0., 500., 1000.):
        folder = directory/f'budget_{budget:g}'/f'slice_{fixed_value:g}'
        if (folder/'scan_grid.npz').exists():
            geometry = actual_region(np.load(folder/'scan_grid.npz'))
            lines = []
            for polygon in polygon_parts(geometry):
                points = np.asarray(polygon.exterior.coords)
                lines.extend([*np.c_[points, np.full(len(points), fixed_value)], [None,None,None]])
            reference.append(dict(fixed_value=fixed_value, lines=lines))
    return reference


def plot_sections(results, directory, output):
    fig, axes = plt.subplots(3, 3, figsize=(11.2, 9.0), sharex=True, sharey=True)
    fig.subplots_adjust(left=.08, right=.985, top=.92, bottom=.14, wspace=.17, hspace=.27)
    limits = (1550., 6100.)
    source = []
    for column, result in enumerate(results):
        cells = final_cells(result)
        for row, fixed_value in enumerate((0., 500., 1000.)):
            ax = axes[row, column]
            stage1 = section(np.asarray(result['stage1_vertices']), fixed_value)
            final = unary_union([section(v, fixed_value) for v in cells])
            grid = np.load(directory/f'budget_{result["budget"]:g}'/f'slice_{fixed_value:g}/scan_grid.npz')
            actual = actual_region(grid)
            source.append(dict(budget=result['budget'], fixed_value=fixed_value,
                               stage1_area=stage1.area, stage2_area=final.area, actual_area=actual.area))
            # Display geometry only; the audit and exported areas retain raw data.
            stage1, final = (set_precision(g,1e-5) for g in (stage1,final))
            visible = box(0, 0, *limits)
            fill(ax, visible.difference(stage1), facecolor='#F4C7C9', edgecolor='none', zorder=0)
            diagonal_fill(ax, stage1.difference(final), limits)
            for cut in result['supports']:
                adjusted = dict(weights=cut['weights'][:2], bound=cut['bound']-cut['weights'][2]*fixed_value)
                if np.linalg.norm(adjusted['weights']) > 1e-9:
                    line = support_line(adjusted, result['square_kw']).intersection(visible)
                    draw_line(ax, line, color='#888888', lw=.35, alpha=.45, zorder=2)
            boundary(ax, stage1, color=BLUE, lw=1.1, ls=(0,(6,2.8)), zorder=3)
            boundary(ax, final, color=GREEN, lw=1.3, ls=(2,(4.5,2.5)), zorder=4)
            boundary(ax, actual, color=RED, lw=1.15, ls=(0,(2.2,2.)), zorder=5)
            ax.set_xlim(0, limits[0])
            ax.set_ylim(0, limits[1])
            ax.set_xticks([0,500,1000,1500])
            ax.set_yticks([0,2000,4000,6000])
            ax.tick_params(direction='out', length=3, width=.6)
            title = f'Budget {result["budget"]:g} | p33 = {fixed_value:g} kW'
            if actual.is_empty:
                title += '\nNo feasible grid points'
            ax.set_title(title, fontsize=9, pad=7)
            ax.annotate(chr(97+row*3+column), (0.,1.04), xycoords='axes fraction',
                        ha='left', va='bottom', fontsize=10, fontweight='bold')
            if column == 0:
                ax.set_ylabel('Node 25 load (kW)')
            if row == 2:
                ax.set_xlabel('Node 18 load (kW)')
    fig.legend(handles=handles(10), loc='lower center', bbox_to_anchor=(.52,.015),
               ncol=3, columnspacing=2., handlelength=3., frameon=False)
    fig.suptitle('Case33bw | Independent sections of the three-dimensional planning region', fontsize=12, y=.98)
    export(fig, output, 'independent_sections')
    save(output/'section_areas.json', source)


def plot_summary(results, directory, output):
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.7))
    fig.subplots_adjust(left=.08, right=.98, bottom=.2, top=.86, wspace=.32)
    colors = ['#586D82', '#3779AF', '#138B52']
    for result, color in zip(results, colors):
        records = result['history']
        axes[0].plot(np.arange(len(records)), [r['selected_gap_kw'] for r in records],
                     color=color, lw=1.2, label=f'Budget {result["budget"]:g}')
    axes[0].axhline(20, color='#B8323D', ls='--', lw=1)
    axes[0].set_yscale('log')
    axes[0].set_yticks([20,100,1000],['20','100','1000'])
    axes[0].set_xlabel('Objective-space branch / refinement step')
    axes[0].set_ylabel('Largest certified gap (kW)')
    axes[0].legend(fontsize=8, frameon=False)
    for index, (result, color) in enumerate(zip(results, colors)):
        audit = json.loads((directory/f'budget_{result["budget"]:g}/audit.json').read_text())
        low, high = audit['inner_volume_lower']/1e9, audit['stage2_volume']/1e9
        axes[1].plot([index, index], [low, high], color=color, lw=4, solid_capstyle='butt')
        axes[1].plot(index, high, 'o', color=color, ms=5)
        axes[1].plot(index, low, '_', color=color, ms=12)
        axes[1].annotate(f'{low:.3f}–{high:.3f}', (index, high), xytext=(0,8), textcoords='offset points',
                         ha='center', fontsize=8)
    axes[1].set_xticks(range(3), [f'{r["budget"]:g}' for r in results])
    axes[1].set_xlim(-.6, 2.6)
    axes[1].set_ylim(bottom=0)
    axes[1].margins(y=.22)
    axes[1].set_xlabel('Construction budget (relative units)')
    axes[1].set_ylabel('Certified SOCP volume interval (10⁹ kW³)')
    for letter, ax in zip('ab', axes):
        ax.annotate(letter, (0,1.06), xycoords='axes fraction', fontsize=11, fontweight='bold')
    fig.text(.5,.035,'Recorded incremental runs include checkpoint resumes and algorithm revisions.',
             ha='center',fontsize=8,color='#5f6b7a')
    export(fig, output, 'certification_and_budget')


def mesh_on_axes(ax, geometry, color, *, opacity=0., dash='--', width=.65, hatch=False):
    vertices = np.asarray(geometry['vertices'])
    if len(vertices) and opacity:
        collection = Poly3DCollection(vertices[np.asarray(geometry['triangles'])],
                                       facecolors=color, edgecolors='none', alpha=opacity)
        ax.add_collection3d(collection)
    pieces, piece = [], []
    rows = geometry['hatches'] if hatch else geometry['lines']
    for point in rows:
        if point[0] is None:
            if len(piece) >= 2:
                pieces.append(piece)
            piece = []
        else:
            piece.append(point)
    if len(piece) >= 2:
        pieces.append(piece)
    if pieces and width > 0:
        ax.add_collection3d(Line3DCollection(pieces, colors=color, linewidths=width,
                                            linestyles=dash, alpha=.9))


def format_3d_axes(ax, limits=(1550.,6100.,1500.)):
    ax.set_xlim(0, limits[0])
    ax.set_ylim(0, limits[1])
    ax.set_zlim(0, limits[2])
    ax.set_box_aspect((1,1,1))
    ax.view_init(elev=25, azim=-53)
    ax.set_xlabel('p18 (kW)', labelpad=6, fontsize=9)
    ax.set_ylabel('p25 (kW)', labelpad=6, fontsize=9)
    ax.set_zlabel('p33 (kW)', labelpad=6, fontsize=9)
    ax.tick_params(labelsize=8, pad=1)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.set_major_locator(mpl.ticker.MaxNLocator(3))
        axis.pane.set_facecolor((.98,.985,.99,1.))
        axis.pane.set_edgecolor((.9,.92,.94,1.))


def plot_3d_results(results, directory, output, preview=False):
    fig, axes = plt.subplots(1,3, figsize=(12.,4.7), subplot_kw={'projection':'3d'})
    fig.subplots_adjust(left=.025,right=.97,top=.87,bottom=.18,wspace=.06)
    for index, (ax, result) in enumerate(zip(axes, results)):
        outer = final_cells(result) if result['cells'] else [np.asarray(result['stage1_vertices'])]
        mesh_on_axes(ax, surface(outer), GREEN, opacity=.18, width=.55)
        mesh_on_axes(ax, surface([np.asarray(result['stage1_vertices'])]), BLUE, width=.8)
        for fixed_value in (0.,500.,1000.):
            path = directory/f'budget_{result["budget"]:g}'/f'slice_{fixed_value:g}/scan_grid.npz'
            if path.exists():
                reference = actual_region(np.load(path))
                for polygon in polygon_parts(reference):
                    points = np.asarray(polygon.exterior.coords)
                    ax.plot(points[:,0], points[:,1], np.full(len(points),fixed_value),
                            color=RED, lw=1.2, ls=(0,(2.2,2.)))
        format_3d_axes(ax)
        ax.set_title(f'Budget {result["budget"]:g}', fontsize=11, pad=15)
        ax.text2D(0.,1.01,chr(97+index),transform=ax.transAxes,fontsize=11,fontweight='bold')
    legend = [Line2D([],[],color=BLUE,ls='--',lw=1.2,label='Stage 1 outer envelope'),
              Line2D([],[],color=GREEN,ls='--',lw=1.2,label='Current outer envelope' if preview else 'Stage 2 final envelope'),
              Line2D([],[],color=RED,ls=':',lw=1.3,label='Actual region: 10 kW sections')]
    fig.legend(handles=legend,loc='lower center',bbox_to_anchor=(.5,.075),ncol=3,frameon=False)
    fig.text(.5,.035,'Construction preview — certification pending; reference sections shown only where available.' if preview
             else 'Axes use different kW ranges. Red curves are independently scanned sections, not a full 3-D actual surface.',
             ha='center',fontsize=8,color='#5f6b7a')
    export(fig,output,'construction_3d_preview' if preview else 'planning_regions_3d')


def plot_cut_stages(result, output):
    """One representative budget: full initial cube, then enlarged cut geometry."""
    cube = MASKS*result['square_kw']
    current, support_removed = cube, []
    for cut in result['supports']:
        discarded = clip_polytope(current, -cut['bound'], np.asarray(cut['weights']))
        if len(discarded):
            support_removed.append(discarded)
        current = clip_polytope(current, cut['bound'], -np.asarray(cut['weights']))
    outer, branch_removed = [current], []
    for cut in result['branch_cuts']:
        normals = np.asarray(cut.get('cut_normals',np.eye(3)))
        threshold = np.asarray(cut['threshold'])
        following = []
        for vertices in outer:
            if np.any((vertices @ normals.T).max(axis=0) <= threshold+1e-8):
                following.append(vertices)
                continue
            discarded = vertices
            for normal, bound in zip(normals,threshold):
                discarded = clip_polytope(discarded,-bound,normal)
            if len(discarded):
                branch_removed.append(discarded)
            following.extend(subtract_disjunction(vertices,threshold,normals))
        outer = following
    fig, axes = plt.subplots(1,3,figsize=(12.,4.8),subplot_kw={'projection':'3d'})
    fig.subplots_adjust(left=.025,right=.97,top=.86,bottom=.19,wspace=.06)
    mesh_on_axes(axes[0],surface([cube]),'#8A99A9',opacity=.06,width=.9)
    mesh_on_axes(axes[1],surface([clip_view(v) for v in support_removed]),RED,opacity=.10,width=0.)
    mesh_on_axes(axes[1],surface([current]),BLUE,opacity=.12,width=.8)
    removed = surface(branch_removed,hatching=True)
    mesh_on_axes(axes[2],removed,RED,opacity=.08,width=0.)
    mesh_on_axes(axes[2],removed,RED,hatch=True,dash='-',width=.35)
    mesh_on_axes(axes[2],surface(outer),GREEN,opacity=.17,width=.6)
    mesh_on_axes(axes[2],surface([current]),BLUE,width=.8)
    titles = ['Initial cube','Stage 1: global support cuts','Stage 2: nonconvex exclusions']
    for index,(ax,title) in enumerate(zip(axes,titles)):
        format_3d_axes(ax,(6855.,6855.,6855.) if index == 0 else (1550.,6100.,1500.))
        ax.set_title(title,fontsize=10,pad=15)
        ax.text2D(0.,1.01,chr(97+index),transform=ax.transAxes,fontsize=11,fontweight='bold')
    fig.text(.5,.09,f'Budget {result["budget"]:g}. Solid red: stage 1 removed; red hatching: stage 2 removed.',
             ha='center',fontsize=9)
    fig.text(.5,.045,'Panels b–c enlarge the boundary with unequal axis ranges. Interactive replay contains every cut.',
             ha='center',fontsize=8,color='#5f6b7a')
    export(fig,output,'two_stage_3d_cuts')


HTML = r'''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>Case33bw · 三维规划域</title>
<style>*{box-sizing:border-box}body{font-family:Arial,"Microsoft YaHei",sans-serif;margin:0;color:#253346;background:#f4f6f8}header{padding:24px 32px;background:#fff;border-bottom:1px solid #dce2e8}h1{font-size:25px;margin:0 0 10px}header p{margin:0;color:#5f6b7a;line-height:1.6}.layout{display:grid;grid-template-columns:300px 1fr;gap:18px;padding:20px}.card{background:#fff;border:1px solid #e0e5eb;border-radius:12px;padding:20px}h2{font-size:16px;margin:0 0 12px}select,button{padding:9px 12px;border:1px solid #cad3de;border-radius:7px;background:white;font-size:14px}select{width:100%}button{cursor:pointer}button.active{background:#e8f2ec;border-color:#138b52}label{display:block;margin:14px 0;font-size:14px}input[type=range]{width:100%;accent-color:#138b52}#plot{height:690px;width:100%}.small{font-size:12px;line-height:1.65;color:#6b7685}.metric{font-size:22px;color:#138b52;margin:8px 0}.formula{font-family:monospace;font-size:13px;line-height:1.7;overflow-wrap:anywhere;color:#5e3f89;background:#f6f2fa;padding:12px;border-radius:8px}.legend{line-height:2;font-size:13px}.chip{display:inline-block;width:24px;border-top:2px dashed;margin-right:7px;vertical-align:middle}.footer{padding:0 24px 22px;color:#6b7685;font-size:12px}@media(max-width:950px){.layout{grid-template-columns:1fr}#plot{height:560px}}</style>
<header><h1>Case33bw · 三维可行规划域</h1><p>三个建设预算 · 全域 SOCP 距离认证 · 独立 10 kW 二维截面验证<br>拖动图形旋转，滚轮缩放。红色参考边界仅表示已扫描的截面。</p></header>
<div class="layout"><aside class="card"><h2>建设预算</h2><select id="budget"></select><p class="small">相对投资单位；相同的 8 个候选升级项目。</p><h2>逐步切割</h2><input id="step" type="range" min="0" value="0"><div id="stepname" style="margin:10px 0;font-weight:bold"></div><div class="formula" id="formula"></div><p><button id="prev">上一步</button> <button id="play">播放</button> <button id="next">下一步</button></p><p><button id="start">立方体</button> <button id="stage1">阶段一</button> <button id="final">最终结果</button></p><label><input id="zoom" type="checkbox" checked> 放大边界区域</label><label><input id="removed" type="checkbox"> 显示已切除体积</label><label><input id="hatch" type="checkbox" checked> 显示第二阶段红色斜纹</label><label><input id="reference" type="checkbox" checked> 显示独立扫描截面</label><label>固定 p33 截面<select id="slice"><option value="all">全部三个截面</option><option value="0">0 kW</option><option value="500">500 kW</option><option value="1000">1000 kW</option><option value="none">隐藏</option></select></label><h2>最终全域认证间隙</h2><div class="metric" id="gap"></div><p class="small">有向 L∞ 距离，三个负荷坐标允许同时变化；绿色包络是外近似，不能把其中每一点都标为已证可行。</p></aside><main class="card"><div id="plot"></div><div class="legend"><span class="chip" style="color:#2363B0"></span>蓝：第一阶段外包络　<span class="chip" style="color:#138B52"></span>绿：第二阶段外域　<span class="chip" style="color:#B8323D"></span>红：Actual region（10 kW 网格截面）<br>浅红体积：支持割切除；红色斜纹：第二阶段切除；紫色：当前支持切面。</div><p class="small" id="details"></p></main></div><div class="footer">初始立方体为 [0, 6855]³ kW³。阶段二按最大认证间隙选择节点，定界覆盖全部允许网架。预算越大，可用建设集合越大；计算未枚举全部网架。</div>
<script>__PLOTLY__</script><script>
const compressed='__DATA__';let datasets=[],timer=null;
function coords(rows,k){return rows.map(p=>p[k])}function mesh(g,color,opacity,name){return {type:'mesh3d',x:coords(g.vertices,0),y:coords(g.vertices,1),z:coords(g.vertices,2),i:g.triangles.map(t=>t[0]),j:g.triangles.map(t=>t[1]),k:g.triangles.map(t=>t[2]),color,opacity,flatshading:true,hoverinfo:'skip',name,showscale:false,showlegend:false,lighting:{ambient:.72,diffuse:.5,specular:.12}}}
function line(rows,color,dash,width,name){return {type:'scatter3d',mode:'lines',x:coords(rows,0),y:coords(rows,1),z:coords(rows,2),line:{color,dash,width},hoverinfo:'skip',name,showlegend:false}}
const el=id=>document.getElementById(id);function current(){return datasets[+el('budget').value]}
function geometry(d,g){let lookup=new Map(),vertices=[],triangles=[];for(let id of g.t){let face=[];for(let v of d.geometry.triangles[id]){if(!lookup.has(v)){lookup.set(v,vertices.length);vertices.push(d.geometry.vertices[v])}face.push(lookup.get(v))}triangles.push(face)}let lines=[];for(let id of g.l||[]){let s=d.geometry.segments[id];lines.push(d.geometry.vertices[s[0]],d.geometry.vertices[s[1]],[null,null,null])}let hatches=[];for(let id of g.h||[]){let s=d.geometry.segments[id];hatches.push(d.geometry.vertices[s[0]],d.geometry.vertices[s[1]],[null,null,null])}return {vertices,triangles,lines,hatches}}
function boundaryAt(d,index){let state=d.boundaryState;if(!state||state.index>index)state={index:-1,t:new Set(),l:new Set()};for(let n=state.index+1;n<=index;n++){let p=d.frames[n].boundary;for(let k of ['t','l']){for(let id of p[k+'_remove'])state[k].delete(id);for(let id of p[k+'_add'])state[k].add(id)}}state.index=index;d.boundaryState=state;return geometry(d,{t:[...state.t],l:[...state.l],h:[]})}
function draw(){let d=current(),index=+el('step').value,f=d.frames[index],traces=[];el('stepname').textContent=`${index}/${d.frames.length-1} · ${f.label}`;el('formula').textContent=f.formula;el('gap').textContent=d.max_gap_kw===null?'待计算':d.max_gap_kw.toFixed(2)+' kW';let isFinal=f.stage===3;
let shown=boundaryAt(d,index),zoom=el('zoom').checked&&index>=7,removed=[{t:[],l:[],h:[]},{t:[],l:[],h:[]}];if(index===d.frames.length-1){removed=zoom?d.removed_final_zoom:d.removed_final}else{for(let n=1;n<=index;n++){let t=d.frames[n],g=zoom&&t.removed_zoom?t.removed_zoom:t.removed;if(g){let r=removed[t.stage===1?0:1];r.t.push(...g.t);r.h.push(...g.h)}}}
if(el('removed').checked){for(let k=0;k<2;k++)if(removed[k].t.length)traces.push(mesh(geometry(d,removed[k]),k===0?'#df8f94':'#cb6670',.10,'Removed'))}
if(f.stage>=2){let first=geometry(d,d.stage1);traces.push(line(first.lines,'#2363B0','dash',3,'Stage 1'));traces.push(mesh(shown,'#c3e3d1',.24,'Stage 2'));traces.push(line(shown.lines,'#138B52','dash',3,'Stage 2'));if(el('hatch').checked){let hatched=geometry(d,{t:[],l:[],h:removed[1].h});traces.push(line(hatched.hatches,'#B8323D','solid',1.5,'Stage 2 removed'));}}else{traces.push(mesh(shown,f.stage===0?'#d8dee7':'#d0e0f2',.17,'Envelope'));traces.push(line(shown.lines,f.stage===0?'#6f7c8b':'#2363B0','dash',3,'Envelope'));}
if(f.plane&&f.plane.t.length){let plane=geometry(d,f.plane);traces.push(mesh(plane,'#9067b5',.38,'Current cut'));traces.push(line(plane.lines,'#855FA8','solid',4,'Current cut'));}
if(el('reference').checked&&el('slice').value!=='none'){for(let r of d.reference){if(el('slice').value==='all'||+el('slice').value===r.fixed_value)traces.push(line(r.lines,'#B8323D','dot',5,'Actual region'));}}
let ranges=zoom?[[0,1550],[0,6100],[0,1500]]:[[0,6855],[0,6855],[0,6855]];let layout={paper_bgcolor:'white',margin:{l:0,r:0,t:20,b:0},uirevision:'camera',scene:{xaxis:{title:'p18 (kW)',range:ranges[0],backgroundcolor:'#fafbfc',gridcolor:'#e0e5eb'},yaxis:{title:'p25 (kW)',range:ranges[1],backgroundcolor:'#fafbfc',gridcolor:'#e0e5eb'},zaxis:{title:'p33 (kW)',range:ranges[2],backgroundcolor:'#fafbfc',gridcolor:'#e0e5eb'},aspectmode:'cube',camera:{eye:{x:1.5,y:1.6,z:1.1}}}};Plotly.react('plot',traces,layout,{responsive:true,displaylogo:false});el('details').textContent=`预算 ${d.budget}；${d.status}。${zoom?'当前各轴按显示范围缩放，形状比例不等于等 kW 尺度。':'当前三轴采用同一 kW 尺度。'} 当前外边界由真实支持平面和析取割重建；没有平滑非凸凹口。`;}
function choose(){let d=current();el('step').max=d.frames.length-1;el('step').value=d.frames.length-1;draw()}
el('budget').onchange=choose;el('step').oninput=draw;for(let id of ['zoom','removed','hatch','reference','slice'])el(id).onchange=draw;
el('prev').onclick=()=>{el('step').value=Math.max(0,+el('step').value-1);draw()};el('next').onclick=()=>{el('step').value=Math.min(+el('step').max,+el('step').value+1);draw()};el('start').onclick=()=>{el('step').value=0;draw()};el('stage1').onclick=()=>{el('step').value=current().frames.findIndex(f=>f.stage===2)-1;draw()};el('final').onclick=()=>{el('step').value=el('step').max;draw()};el('play').onclick=()=>{if(timer){clearInterval(timer);timer=null;el('play').textContent='播放'}else{if(+el('step').value>=+el('step').max){el('step').value=0;draw()}el('play').textContent='暂停';timer=setInterval(()=>{if(+el('step').value>=+el('step').max){clearInterval(timer);timer=null;el('play').textContent='播放';return}el('next').click()},1000)}};
(async()=>{const bytes=Uint8Array.from(atob(compressed),c=>c.charCodeAt(0));const stream=new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip'));datasets=JSON.parse(await new Response(stream).text());el('budget').innerHTML=datasets.map((d,i)=>`<option value="${i}">预算 ${d.budget}</option>`).join('');choose()})();
</script></html>'''


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', type=Path)
    parser.add_argument('--preview', action='store_true')
    parser.add_argument('--html-only', action='store_true')
    args = parser.parse_args()
    output = args.directory/'figures'
    output.mkdir(parents=True, exist_ok=True)
    results = [json.loads((args.directory/f'budget_{b}/result.json').read_text(encoding='utf-8')) for b in (0,2,4)]
    if not args.preview:
        assert all(r['status'] == 'certified' for r in results)
        audits = [json.loads((args.directory/f'budget_{b}/audit.json').read_text(encoding='utf-8')) for b in (0,2,4)]
        assert all(a['passed'] and len(a['slices']) == 3 for a in audits)
        plot_sections(results, args.directory, output)
        plot_summary(results, args.directory, output)
        plot_3d_results(results, args.directory, output)
        plot_cut_stages(results[1],output)
    payload = ([json.loads((output/f'budget_{b}_geometry.json').read_text(encoding='utf-8')) for b in (0,2,4)]
               if args.html_only else [make_frames(result, output) for result in results])
    for data,result in zip(payload,results):
        if not args.preview:
            assert data['status'] == 'certified'
            assert len(data['frames']) == len(result['supports'])+len(result['branch_cuts'])+2
        data['reference'] = reference_lines(args.directory,result['budget'])
        save(output/f'budget_{result["budget"]:g}_geometry.json',data)
    raw = json.dumps(payload, default=lambda v: v.tolist() if isinstance(v, np.ndarray) else v,
                     separators=(',', ':')).encode()
    encoded = base64.b64encode(gzip.compress(raw)).decode()
    (output/'three_dimensional.html').write_text(HTML.replace('__PLOTLY__', get_plotlyjs()).replace('__DATA__',encoded), encoding='utf-8')
    for result in results:
        with (output/f'budget_{result["budget"]:g}_cuts.csv').open('w', encoding='utf-8-sig', newline='') as handle:
            writer = csv.writer(handle)
            writer.writerow(['stage','index','facet','w18','w25','w33','threshold_kw','global_deficit_kw'])
            for i, cut in enumerate(result['supports']):
                writer.writerow([1,i+1,1,*cut['weights'],cut['bound'],''])
            for i, cut in enumerate(result['branch_cuts']):
                for facet,(normal,threshold) in enumerate(zip(cut.get('cut_normals',np.eye(3)),cut['threshold'])):
                    writer.writerow([2,i+1,facet+1,*normal,threshold,cut['bound']])
    print(f'HTML payload {len(raw)/1e6:.2f} MB, compressed {len(encoded)/1e6:.2f} MB', flush=True)


if __name__ == '__main__':
    with threadpool_limits(limits=1):
        main()
