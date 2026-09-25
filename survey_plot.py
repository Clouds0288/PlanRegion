"""One concept figure and a complete conditional-value history, from solved data."""
from __future__ import annotations
import json
from pathlib import Path
import argparse
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Patch, Rectangle
import numpy as np
from shapely.affinity import affine_transform
from shapely.geometry import LineString, shape

from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch
from shapely.geometry.polygon import orient

from Network.concept5 import ROUTES
from survey import STOPPING_AREA_TOLERANCE

FIGURES = Path('results/concept5/figures')
WIDTH_MM = 183
GREEN, GREEN_LIGHT = '#23845D', '#B5D9C6'
RED, RED_LIGHT = '#BC4D46', '#F2C9C3'
DARK, MUTED, BASE, CONFIRMED, UNKNOWN = '#283C49','#6E808E','#A7B2BA','#ADCADB','#EDF1F5'
COLORS = dict(zip(ROUTES,['#3F6E9E','#A87B42','#7864A5','#B0576A','#32888D','#8A939B']))
MARKERS = dict(zip(ROUTES,['o','s','D','^','v','x']))
mpl.rcParams.update({'font.family':'sans-serif',
    'font.sans-serif':['Microsoft YaHei','Arial','DejaVu Sans'],
    'font.size':7.5,'axes.labelsize':8.,'axes.titlesize':8.,
    'xtick.labelsize':7.,'ytick.labelsize':7.,'legend.fontsize':7.,
    'svg.fonttype':'none','pdf.fonttype':42,'mathtext.fontset':'dejavusans',
    'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.65,
    'legend.frameon':False,'savefig.facecolor':'white','figure.facecolor':'white'})


def polygon_parts(geometry):
    if geometry.is_empty:
        return []
    if geometry.geom_type == 'Polygon':
        return [geometry]
    return [part for piece in geometry.geoms for part in polygon_parts(piece)]


def fill_geometry(ax, geometry, color, *, edgecolor='none', hatch=None, zorder=1):
    for polygon in polygon_parts(geometry):
        polygon = orient(polygon, sign=1.)
        paths = []
        for ring in [polygon.exterior, *polygon.interiors]:
            vertices = np.asarray(ring.coords)
            codes = np.full(len(vertices), MplPath.LINETO, dtype=np.uint8)
            codes[0], codes[-1] = MplPath.MOVETO, MplPath.CLOSEPOLY
            paths.append(MplPath(vertices, codes))
        ax.add_patch(PathPatch(MplPath.make_compound_path(*paths), facecolor=color,
                               edgecolor=edgecolor, linewidth=.35, zorder=zorder))
    if hatch and not geometry.is_empty:
        pending = getattr(ax.figure, '_survey_hatches', [])
        pending.append((ax, geometry, edgecolor, zorder+.1))
        ax.figure._survey_hatches = pending


def boundary(ax, geometry, color, *, linestyle='-', linewidth=1.1, zorder=4):
    for polygon in polygon_parts(geometry):
        coordinates = np.asarray(polygon.exterior.coords)
        ax.plot(coordinates[:, 0], coordinates[:, 1], color=color, linestyle=linestyle,
                linewidth=linewidth, zorder=zorder)


def save_figure(fig, name):
    stem = FIGURES/name
    fig.canvas.draw()
    for ax,geometry,color,zorder in getattr(fig,'_survey_hatches',[]):
        matrix = ax.transData.get_affine().get_matrix()
        transformed = affine_transform(geometry,[matrix[0,0],matrix[0,1],matrix[1,0],
                                                  matrix[1,1],matrix[0,2],matrix[1,2]])
        x0,y0,x1,y1 = transformed.bounds
        spacing = 5.*fig.dpi/72.*np.sqrt(2.)
        segments = []
        for offset in np.arange(y0-x1-spacing,y1-x0+spacing,spacing):
            clipped = transformed.intersection(LineString([(x0-1,x0-1+offset),(x1+1,x1+1+offset)]))
            parts = [clipped] if clipped.geom_type=='LineString' else getattr(clipped,'geoms',[])
            for part in parts:
                if part.geom_type=='LineString' and not part.is_empty:
                    segments.append(ax.transData.inverted().transform(np.asarray(part.coords)))
        ax.add_collection(LineCollection(segments,colors=color,linewidths=.35,zorder=zorder))
    fig.canvas.draw()
    fig.savefig(str(stem)+'.pdf')
    fig.savefig(str(stem)+'.svg')
    fig.savefig(str(stem)+'.png',dpi=300)
    plt.close(fig)


def domain_axes(ax):
    ax.set(xlim=(0,100),ylim=(0,100),xlabel=r'$p_1$ (kW)',ylabel=r'$p_2$ (kW)')
    ax.set_xticks([0,50,100])
    ax.set_yticks([0,50,100])
    ax.tick_params(length=2.4,pad=2)


def topology(ax):
    positions = {0:(0,2.8),3:(-1,1.5),4:(1,1.5),1:(-1,0),2:(1,0)}
    paths = {
        'A':[(0,2.8),(-1.75,2.8),(-1.75,0),(-1,0)],
        'B':[(0,2.8),(1.75,2.8),(1.75,0),(1,0)],
        'C':[(-1,0),(-1,-.48),(1,-.48),(1,0)],
        'D':[(-1,0),(1,1.5)], 'E':[(1,0),(-1,1.5)],
        'F':[(-1,1.5),(1,1.5)]}
    radius = .19
    for a,b in [(0,3),(3,1),(0,4),(4,2)]:
        points = np.asarray([positions[a],positions[b]],float)
        unit = (points[1]-points[0])/np.linalg.norm(points[1]-points[0])
        points[0] += radius*unit
        points[-1] -= radius*unit
        ax.plot(points[:,0],points[:,1],color=DARK,lw=1.1,zorder=1)
    for route,path in paths.items():
        points = np.asarray(path,float)
        points[0] += radius*(points[1]-points[0])/np.linalg.norm(points[1]-points[0])
        points[-1] += radius*(points[-2]-points[-1])/np.linalg.norm(points[-2]-points[-1])
        ax.plot(points[:,0],points[:,1],color=COLORS[route],lw=1.,ls=(0,(3,2)),zorder=2)
    for node,(x,y) in positions.items():
        ax.add_patch(Circle((x,y),radius,facecolor=DARK if node==0 else '#EDF3F7',
                            edgecolor=DARK,lw=.75,zorder=3))
        ax.text(x,y,str(node),ha='center',va='center',fontsize=7.,zorder=4,
                color='white' if node==0 else DARK)
    labels = {'A':(-1.96,1.2),'B':(1.96,1.2),'C':(0,-.74),
              'D':(-.40,.18),'E':(.40,.18),'F':(0,1.77)}
    for road,(x,y) in labels.items():
        ax.text(x,y,road,ha='center',va='center',color=COLORS[road],weight='bold',fontsize=8.)
    ax.text(0,3.19,'电源',ha='center',fontsize=7.)
    ax.text(-1.43,-.34,r'$p_1$',fontsize=8.,ha='center')
    ax.text(1.43,-.34,r'$p_2$',fontsize=8.,ha='center')
    ax.set(xlim=(-2.25,2.25),ylim=(-.92,3.44),aspect='equal')
    ax.axis('off')


def plot_nonconvexity(ax, result):
    final = shape(result['states'][-1]['confirmed'])
    fill_geometry(ax,final.convex_hull,UNKNOWN)
    fill_geometry(ax,final,CONFIRMED)
    boundary(ax,final,DARK,linewidth=1.)
    witness = result['nonconvexity_witness']
    a,b,m = [np.asarray(witness[key]) for key in ('p_a','p_b','p_mid')]
    ax.plot([a[0],b[0]],[a[1],b[1]],color=MUTED,ls='--',lw=.8,zorder=5)
    ax.scatter([a[0],b[0]],[a[1],b[1]],s=15,color=DARK,zorder=6)
    ax.scatter([m[0]],[m[1]],s=22,color=RED,marker='x',linewidth=1.2,zorder=7)
    for p,label,offset in [(a,'U',(-10,-9)),(b,'V',(5,2)),(m,'M',(3,6))]:
        ax.text(*(p+offset),label,color=RED if label=='M' else DARK,fontsize=7.5)
    domain_axes(ax)


def plot_value_history(ax, result, *, detailed=False):
    for route in ROUTES:
        rows = [row for row in result['all_candidates'] if row['route']==route]
        x = [row['step'] for row in rows]
        y = [row['information_efficiency'] for row in rows]
        ax.plot(x,y,marker=MARKERS[route],markersize=4.,lw=1.15,
                color=COLORS[route],label=route,zorder=4 if route!='F' else 2)
        lower = [row['information_efficiency_lower'] for row in rows]
        upper = [row['information_efficiency_upper'] for row in rows]
        ax.fill_between(x,lower,upper,color=COLORS[route],alpha=.18,lw=0)
    for action in result['trace']:
        ax.scatter(action['step']-1,action['information_efficiency'],s=64,facecolors='none',
                   edgecolors=DARK,linewidths=.9,zorder=8)
    ax.axhline(STOPPING_AREA_TOLERANCE,color=MUTED,ls=':',lw=.7,zorder=1)
    ax.set(xlim=(-.16,5.18),ylim=(-90,2000),ylabel='边际信息价值\n(kW² / 勘察单位)')
    ax.set_xticks(range(6),['0','1 (A+)','2 (B−)','3 (E+)','4 (D−)','5 (C+)'])
    ax.set_yticks([0,500,1000,1500,2000])
    ax.set_xlabel('决策时刻：已完成的勘察次数（括号为刚获得的结果）',fontsize=7.)
    ax.grid(axis='y',color='#E6EBEF',lw=.45,zorder=0)
    ax.legend(loc='upper right',ncol=6,columnspacing=.85,handlelength=1.3,
              bbox_to_anchor=(1.005,1.25 if detailed else 1.28),fontsize=7.)
    ax.text(4.80,180,'F ≈ 0；停止',ha='center',va='bottom',fontsize=7.,color=MUTED)
    if detailed:
        ax.text(3.12,1600,'D：E 可用后价值上升',fontsize=7.5,color=COLORS['D'])


def concept_figure(result):
    height_mm = 255
    fig = plt.figure(figsize=(WIDTH_MM/25.4,height_mm/25.4))
    fig.text(.055,.976,'有限预算下，勘察逐步识别可实现的非凸规划域',fontsize=11.,weight='bold',color=DARK)
    fig.text(.055,.952,'五节点合成示例  |  建设预算 4；单次勘察费 1；独立可用先验 0.8',fontsize=7.5,color=MUTED)
    fig.text(.055,.925,'a  五节点网架与六条候选道路',fontsize=8.,weight='bold',color=DARK)
    net = fig.add_axes([.055,.758,.355,.158],label='topology')
    topology(net)
    fig.text(.06,.750,'实线：既有；虚线：候选；交叉处无连接点',fontsize=6.7,color=MUTED)
    fig.text(.485,.925,'b  固定方案域的并集可以非凸',fontsize=8.,weight='bold',color=DARK)
    nonconvex = fig.add_axes([.50,.772,.195,.14],label='nonconvexity')
    plot_nonconvexity(nonconvex,result)
    fig.text(.746,.890,'U、V 分别可行\n中点 M 不可行',fontsize=7.4,color=DARK,linespacing=1.7)
    fig.text(.746,.823,'各点可选择不同方案\n离散建线与径向运行\n使方案并集未必凸',fontsize=6.8,color=MUTED,linespacing=1.7)
    handles = [Patch(facecolor=UNKNOWN,label='尚未排除的可能范围'),
               Patch(facecolor=BASE,label='初始调度域'),
               Patch(facecolor=CONFIRMED,label='此前已确认'),
               Patch(facecolor=GREEN_LIGHT,label='本轮新增 +'),
               Patch(facecolor=RED_LIGHT,edgecolor=RED,label='本轮删除 −')]
    fig.legend(handles=handles,loc='center',bbox_to_anchor=(.5,.725),ncol=5,
               handlelength=1.2,columnspacing=1.05,fontsize=6.5)
    axes = []
    initial_optimistic = shape(result['states'][0]['optimistic'])
    baseline = shape(result['states'][0]['confirmed'])
    width = .218
    height = width*WIDTH_MM/height_mm
    for index,state in enumerate(result['states']):
        row,col = divmod(index,3)
        left,bottom = .07+col*.32,.543-row*.248
        ax = fig.add_axes([left,bottom,width,height],label=f'state_{index}')
        axes.append(ax)
        current,optimistic = shape(state['confirmed']),shape(state['optimistic'])
        fill_geometry(ax,optimistic,UNKNOWN)
        fill_geometry(ax,current,CONFIRMED)
        fill_geometry(ax,baseline,BASE,zorder=2)
        if index:
            before = result['states'][index-1]
            added = current.difference(shape(before['confirmed']))
            removed = shape(before['optimistic']).difference(optimistic)
            fill_geometry(ax,added,GREEN_LIGHT,edgecolor=GREEN,zorder=3)
            fill_geometry(ax,removed,RED_LIGHT,edgecolor=RED,hatch='///',zorder=3)
        boundary(ax,initial_optimistic,'#B9C2CA',linestyle=':',linewidth=.6,zorder=4)
        boundary(ax,optimistic,DARK,linestyle='--',linewidth=.8,zorder=5)
        boundary(ax,current,'#426783',linewidth=1.,zorder=6)
        domain_axes(ax)
        if index==0:
            title = 'c  初始：尚未勘察'
            detail = f'初始调度域占乐观域 {result["baseline_share"]:.1%}'
            color = MUTED
        else:
            action = result['trace'][index-1]
            sign = '+' if action['survey_observation'] else '−'
            status = '可用' if action['survey_observation'] else '不可用'
            title = f'{chr(99+index)}  {index} · {action["route"]} {status}'
            area = action['realized_gain']+action['realized_removal']
            detail = f'{"新增" if action["survey_observation"] else "删除"} {sign}{area:,.1f} kW²'
            color = GREEN if action['survey_observation'] else RED
        fig.text(left-.005,bottom+height+.012,title,fontsize=7.8,weight='bold',color=DARK)
        fig.text(left-.005,bottom-.053,detail,fontsize=6.8,color=color)
        fig.text(left-.005,bottom-.066,f'确认 {state["confirmed_area"]:,.0f} / 乐观 {state["optimistic_area"]:,.0f} kW²',
                 fontsize=6.5,color=MUTED)
    fig.text(.055,.205,'i  动态边际价值：圆环为当轮选择',fontsize=8.,weight='bold',color=DARK)
    value = fig.add_axes([.12,.060,.80,.115],label='value_history')
    plot_value_history(value,result)
    fig.text(.055,.007,f'确认域向外扩、乐观域向内收缩；停止时整体间隙上界 {result["states"][-1]["information_gap_upper"]:.2f} kW² < {STOPPING_AREA_TOLERANCE:g} kW²，F 不再勘察。',
             fontsize=6.8,color=MUTED)
    save_figure(fig,'concept5_overview')


def history_figure(result):
    fig = plt.figure(figsize=(WIDTH_MM/25.4,160/25.4))
    fig.text(.055,.96,'道路价值是条件量：比较同一决策时刻，而非不同道路的最后一次评分',
             fontsize=9.8,weight='bold',color=DARK)
    fig.text(.055,.92,'a  每条路线一条轨迹；圆环 = 当轮选中；勘察后轨迹结束，不能补成零',fontsize=7.8,color=DARK)
    curve = fig.add_axes([.115,.575,.81,.275],label='conditional_value')
    plot_value_history(curve,result,detailed=True)
    fig.text(.055,.456,'b  完整评分矩阵：黑框 = 当轮选择；“—” = 已勘察，不再评分',fontsize=7.8,color=DARK)
    matrix = np.full((6,6),np.nan)
    for row in result['all_candidates']:
        matrix[list(ROUTES).index(row['route']),row['step']] = row['information_efficiency']
    heat = fig.add_axes([.115,.105,.81,.305],label='score_matrix')
    cmap = mpl.colormaps['Blues'].copy()
    cmap.set_bad('#F5F6F7')
    heat.imshow(matrix,aspect='auto',cmap=cmap,vmin=0,vmax=1900,interpolation='none')
    for i in range(6):
        for j in range(6):
            value = matrix[i,j]
            label = '—' if np.isnan(value) else f'{value:,.1f}'
            heat.text(j,i,label,ha='center',va='center',fontsize=7.4,
                      color='white' if value>1100 else MUTED if np.isnan(value) else DARK)
    for row in result['trace']:
        i,j = list(ROUTES).index(row['route']),row['step']-1
        heat.add_patch(Rectangle((j-.465,i-.435),.93,.87,fill=False,edgecolor=DARK,lw=1.1))
    heat.set_yticks(range(6),list(ROUTES))
    heat.set_xticks(range(6),['第 1 次前','第 2 次前','第 3 次前','第 4 次前','第 5 次前','停止检查'])
    heat.tick_params(length=0,pad=6)
    for i,label in enumerate(heat.get_yticklabels()):
        label.set_color(COLORS[list(ROUTES)[i]])
        label.set_weight('bold')
    for spine in heat.spines.values():
        spine.set_visible(False)
    fig.text(.055,.035,'数值单位：kW² / 勘察单位。所有单次勘察费均为 1；区间带为求域数值误差界，小于标记尺寸。',
             fontsize=7.,color=MUTED)
    save_figure(fig,'conditional_road_values')


def render_survey(result, output):
    global FIGURES
    FIGURES = Path(output)
    FIGURES.mkdir(parents=True, exist_ok=True)
    concept_figure(result)
    history_figure(result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, required=True)
    args = parser.parse_args()
    result = json.loads((args.results/'results.json').read_text(encoding='utf-8'))
    render_survey(result, args.results/'figures')


if __name__ == '__main__':
    main()
