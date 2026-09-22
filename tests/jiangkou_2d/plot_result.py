"""One 2D comparison: certified cut polygons, unresolved envelope and AC checks."""
from run_domain import ROOT
import json
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.spatial import ConvexHull
from shapely.geometry import Polygon, LineString, box
from shapely.ops import unary_union
from Network.jiangkou import Jiangkou
from region import initial_polytope, clip_polytope, polytope_volume, contains, halfspaces
from experiment_runtime import save

matplotlib.rcParams.update({'font.family':'sans-serif',
    'font.sans-serif':['Microsoft YaHei','Arial','DejaVu Sans'],
    'font.size':9,'axes.labelsize':10,'axes.titlesize':11,
    'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':8,
    'axes.spines.top':False,'axes.spines.right':False,
    'axes.unicode_minus':False,'pdf.fonttype':42,'svg.fonttype':'none'})


def ordered(poly):
    points=np.asarray(poly,float)
    return points[ConvexHull(points).vertices] if len(points)>=3 else points


def main():
    domain=json.loads((ROOT/'results'/'domain.json').read_text(encoding='utf-8'))
    if domain['status']=='partial_error':
        domain=json.loads((ROOT/'results'/'domain_checkpoint.json').read_text(encoding='utf-8'))
    ac=json.loads((ROOT/'results'/'ac_boundary.json').read_text(encoding='utf-8'))
    c=Jiangkou(load_nodes=('B000025','B000078'))
    bounds=np.repeat(c.power_limit,2)
    # For mandatory one-hot groups, max_x d*x is the sum of per-group maxima.
    # Thus a+max_x(d*x)+b*p>=0 is a valid outer cut for every design.
    projected=initial_polytope(bounds,c.power_limit)
    splits=np.cumsum([len(o.cost) for o in c.line_options])[:-1]
    for cut in np.asarray(domain['cuts']):
        constant=cut[0]+sum(group.max() for group in np.split(cut[3:],splits))
        projected=clip_polytope(projected,constant,cut[1:3]*bounds)
    if domain['status']=='certified':
        envelopes=[np.asarray(p['vertices']) for p in domain['outer']]
    else:
        envelopes=[projected*bounds]
    projected_area=polytope_volume(projected)*np.prod(bounds)
    inner_points=np.array([point for poly in domain['inner'] for point in poly['vertices']])
    assert contains(inner_points/bounds,halfspaces(projected),tolerance=1e-7).all(), 'Outer cuts exclude certified inner points'
    with np.load(ROOT/'results'/'ac_grid.npz') as grid:
        ac_points=grid['points_kw'][grid['ac_verified']==1]
    assert contains(ac_points/bounds,halfspaces(projected),tolerance=1e-7).all(), 'Outer cuts exclude verified AC points'
    save('plot_geometry.json',dict(status=domain['status'],global_outer=envelopes,
        inner=domain['inner'],scheme_geometry=domain['geometry'],
        projected_outer_area_kw2=projected_area,inner_area_kw2=domain['inner_volume'],
        projected_relative_area_gap=1-domain['inner_volume']/projected_area,
        source='Saved globally valid cuts, with exact maximum over independent one-hot line groups'))
    fig,ax=plt.subplots(figsize=(183/25.4,150/25.4))
    fig.subplots_adjust(left=.13,right=.95,bottom=.25,top=.85)
    for poly in envelopes:
        vertices=ordered(poly)
        ax.fill(vertices[:,0],vertices[:,1],facecolor='#eeeeef',edgecolor='#7a7f86',
                linewidth=1,linestyle='--',zorder=1)
    for poly in domain['inner']:
        vertices=ordered(poly['vertices'])
        if len(vertices)>=3:
            ax.fill(vertices[:,0],vertices[:,1],facecolor='#86b6d5',edgecolor='#216b96',
                    linewidth=1.2,zorder=2)
        elif len(vertices)==2:
            ax.plot(vertices[:,0],vertices[:,1],color='#216b96',lw=1.2,zorder=3)
        elif len(vertices)==1:
            ax.scatter(vertices[:,0],vertices[:,1],color='#216b96',s=12,zorder=3)
    scheme_polygons=[Polygon(ordered(row['outer']['vertices'])) for row in domain['geometry']
                     if len(row['inner']['vertices']) and len(row['outer']['vertices'])>=3]
    union=unary_union(scheme_polygons)
    for polygon in ([union] if union.geom_type=='Polygon' else getattr(union,'geoms',[])):
        xy=np.asarray(polygon.exterior.coords)
        ax.plot(xy[:,0],xy[:,1],color='#216b96',lw=.7,ls=':',zorder=3)
    physical=np.asarray(ac['same_contract']['inner_kw'])
    thermal=np.asarray(ac['strict_current']['inner_kw'])
    ax.plot(physical[:,0],physical[:,1],color='#257b53',lw=1.5,zorder=4)
    ax.plot(thermal[:,0],thermal[:,1],color='#d18218',lw=1.1,ls='--',zorder=5)
    nominal=c.original_p[c.selected]
    ax.scatter(*nominal,marker='*',s=85,color='#262b32',zorder=6)
    ax.annotate('原始负荷',nominal,xytext=(9,8),textcoords='offset points',fontsize=8)
    ax.set(xlabel='B000025 有功负荷（kW）',ylabel='B000078 有功负荷（kW）',
           xlim=(0,410),ylim=(0,410),aspect='equal')
    ax.set_xticks(np.arange(0,401,100));ax.set_yticks(np.arange(0,401,100))
    fig.suptitle('江口二维可行规划—调度域与 AC 对比',y=.965,fontsize=13)
    fig.text(.5,.91,'已有线路升级 · 无限预算 · 1000 kVA · 其余负荷固定为 545.45 kW',ha='center',fontsize=8.5)
    handles=[Patch(facecolor='#86b6d5',edgecolor='#216b96',label='切割认证内域'),
        Patch(facecolor='#eeeeef',edgecolor='#7a7f86',linestyle='--',label='全局外界内的未确定部分'),
        Line2D([0],[0],color='#216b96',lw=1,ls=':',label='已有可行证书方案的割外域'),
        Line2D([0],[0],color='#257b53',lw=1.5,label='AC 边界：全最大型号方案'),
        Line2D([0],[0],color='#d18218',lw=1.1,ls='--',label='AC 边界：另加 95% 载流约束')]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.062),
               ncol=2,frameon=False,columnspacing=1.6,handlelength=2.8,labelspacing=.9)
    status={'certified':'覆盖已认证','time_limit':'到达时限，保留部分结果',
            'unknown':'覆盖未完成，保留部分结果','running':'计算进行中'}.get(domain['status'],domain['status'])
    fig.text(.5,.026,f'{status}；有效割 {len(domain["cuts"])} 条。AC 曲线验证一个具体建设方案，不代表全体方案的不可行证明。',
             ha='center',fontsize=7)
    target=ROOT/'figures'/'jiangkou_2d_comparison'
    fig.savefig(target.with_suffix('.pdf'))
    fig.savefig(target.with_suffix('.svg'))
    fig.savefig(target.with_suffix('.png'),dpi=300)
    ax.set_aspect('auto')
    ax.set_xlim(0,125);ax.set_ylim(284,288)
    ax.set_xticks(np.arange(0,126,25));ax.set_yticks(np.arange(284,289))
    # Explicitly clip the detail geometry, so exported vector paths also stay
    # inside the displayed window instead of relying only on a PDF clip mask.
    window=box(0,284,125,288)
    for patch in list(ax.patches):
        clipped=Polygon(patch.get_xy()).intersection(window)
        patch.set_visible(False)
        for part in ([clipped] if clipped.geom_type=='Polygon' else getattr(clipped,'geoms',[])):
            if part.is_empty:
                continue
            xy=np.asarray(part.exterior.coords)
            ax.fill(xy[:,0],xy[:,1],facecolor=patch.get_facecolor(),edgecolor=patch.get_edgecolor(),
                    linewidth=patch.get_linewidth(),linestyle=patch.get_linestyle(),zorder=patch.get_zorder())
    for line in list(ax.lines):
        clipped=LineString(np.c_[line.get_xdata(),line.get_ydata()]).intersection(window)
        line.set_visible(False)
        for part in ([clipped] if clipped.geom_type=='LineString' else getattr(clipped,'geoms',[])):
            if part.is_empty or part.geom_type!='LineString':
                continue
            xy=np.asarray(part.coords)
            ax.plot(xy[:,0],xy[:,1],color=line.get_color(),linewidth=line.get_linewidth(),
                    linestyle=line.get_linestyle(),zorder=line.get_zorder())
    for collection in ax.collections:
        xy=np.asarray(collection.get_offsets())
        inside=(xy[:,0]>=0)&(xy[:,0]<=125)&(xy[:,1]>=284)&(xy[:,1]<=288)
        collection.set_offsets(xy[inside])
    fig.suptitle('江口边界局部：实际载流约束带来的收紧',y=.965,fontsize=13)
    fig.text(.94,.86,'局部放大；两轴比例不同',ha='right',fontsize=7)
    detail=target.with_name(target.name+'_detail')
    fig.savefig(detail.with_suffix('.pdf'))
    fig.savefig(detail.with_suffix('.svg'))
    fig.savefig(detail.with_suffix('.png'),dpi=300)
    plt.close(fig)
    print(str(target.with_suffix('.png')))


if __name__=='__main__':
    main()
