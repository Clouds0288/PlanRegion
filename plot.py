"""区域对比与固定方案切割回放；只消费计算结果。"""
from pathlib import Path
import numpy as np
import plotly.graph_objects as go
from region import clip_polytope, polytope_vertices


def save_replay(result, folder, method='socp', budget_index=-1, design_index=-1):
    record = result.regions[method][budget_index][design_index]
    outer = np.asarray(record['initial'])
    inner = np.zeros((1,3))
    frames = []

    def traces(point):
        mesh = lambda p,c,o,n: go.Mesh3d(x=p[:,0],y=p[:,1],z=p[:,2],alphahull=0,
                                       color=c,opacity=o,name=n,flatshading=True)
        return [mesh(outer,'#AAB7BE',.2,'候选外域'), mesh(inner,'#4286AD',.7,'认证内域'),
                go.Scatter3d(x=[point[0]],y=[point[1]],z=[point[2]],mode='markers',
                             marker=dict(color='#C07B48',size=4),name='本轮查询点')]

    # 先用三条轴向查询建立三维内域；此后每一帧对应真实的一次 SP 查询。
    for i,row in enumerate(record['history']):
        if row['cut'] is not None:
            cut = np.asarray(row['cut'])
            outer = clip_polytope(outer,cut[0],cut[1:])
        inner = polytope_vertices(np.vstack([inner,row['witness']]))
        if i >= 2:
            frames.append(go.Frame(name=str(i+1),data=traces(row['p'])))
    fig = go.Figure(data=frames[0].data,frames=frames)
    axis = lambda k: dict(title=f'节点 {result.load_nodes[k]} 负荷 (kW)',range=[0,result.bounds[k]])
    x = result.metadata['designs'][record['design']]['x']
    fig.update_layout(template='plotly_white',height=670,margin=dict(l=5,r=5,t=45,b=90),
        scene=dict(xaxis=axis(0),yaxis=axis(1),zaxis=axis(2),aspectmode='cube'),
        annotations=[dict(text=f'固定建设方案 x={x}；灰色外域，蓝色认证内域',x=.5,y=1.04,
                          xref='paper',yref='paper',showarrow=False)],
        sliders=[dict(currentvalue=dict(prefix='SP 查询轮次：'),steps=[
            dict(method='animate',label=f.name,args=[[f.name],dict(mode='immediate',
                 frame=dict(duration=0,redraw=True),transition=dict(duration=0))]) for f in frames])],
        updatemenus=[dict(type='buttons',direction='left',x=0,y=-.06,buttons=[
            dict(label='播放',method='animate',args=[None,dict(fromcurrent=True,
                frame=dict(duration=200,redraw=True),transition=dict(duration=0))]),
            dict(label='暂停',method='animate',args=[[None],dict(mode='immediate',
                frame=dict(duration=0,redraw=False))])])])
    path = Path(folder)/'cutting_process.html'
    fig.write_html(path,include_plotlyjs=True,auto_play=False,config=dict(displaylogo=False))
    return path


def voxel_surface(mask, spacing):
    """Exposed cell faces, merged into rectangles without filling holes.

    Geometry and disagreement percentages use the very same labeled cells.
    No global convex hull or smoothing changes the meaning of a colored region.
    """
    mask = np.asarray(mask, dtype=bool)
    vertices, triangles = [], []
    for axis in range(3):
        others = [i for i in range(3) if i != axis]
        oriented = np.moveaxis(mask, axis, 0)
        for side in (-1, 1):
            neighbour = np.zeros_like(oriented)
            if side == 1:
                neighbour[:-1] = oriented[1:]
            else:
                neighbour[1:] = oriented[:-1]
            faces = oriented & ~neighbour
            for plane in range(len(faces)):
                cells = faces[plane].copy()
                for row in range(cells.shape[0]):
                    while cells[row].any():
                        left = int(np.flatnonzero(cells[row])[0])
                        right = left+1
                        while right < cells.shape[1] and cells[row, right]:
                            right += 1
                        bottom = row+1
                        while bottom < cells.shape[0] and cells[bottom, left:right].all():
                            bottom += 1
                        cells[row:bottom, left:right] = False
                        face = np.zeros((4, 3), dtype=float)
                        face[:, axis] = plane+(side == 1)
                        face[:, others[0]] = [row, bottom, bottom, row]
                        face[:, others[1]] = [left, left, right, right]
                        normal = np.cross(face[1]-face[0], face[2]-face[0])
                        if normal[axis]*side < 0:
                            face = face[::-1]
                        first = len(vertices)
                        vertices.extend((face*spacing).tolist())
                        triangles.extend([[first, first+1, first+2], [first, first+2, first+3]])
    return dict(vertices=vertices, triangles=triangles)


def plot_method_comparison(result, budget_index=0):
    """同一预算的四方法对比；曲面与 FR/MR 使用同一份网格标签。"""
    from plotly.subplots import make_subplots
    from vertify import METHOD_NAMES

    fig = make_subplots(rows=2, cols=2, specs=[[{"type": "scene"}]*2]*2,
                        subplot_titles=[f"{letter}  {name}" for letter, name in zip("abcd", METHOD_NAMES)],
                        horizontal_spacing=.02, vertical_spacing=.08)
    categories = [(3, "与 AC 重合", "#4286AD", 1., "common"),
                  (2, "遗漏", "#D43D3D", 1., "missed"),
                  (1, "多余", "#E9B72F", 1., "extra")]
    for panel, labels in enumerate(result.labels[:, budget_index]):
        for code, name, color, opacity, group in categories:
            mesh = voxel_surface(labels == code, result.spacing)
            points = np.asarray(mesh["vertices"]).reshape(-1, 3)
            faces = np.asarray(mesh["triangles"], dtype=int).reshape(-1, 3)
            if len(points):
                trace = go.Mesh3d(x=points[:, 0], y=points[:, 1], z=points[:, 2],
                                 i=faces[:, 0], j=faces[:, 1], k=faces[:, 2],
                                 name="AC 参考域" if panel == 2 else name,
                                 color=color, opacity=opacity, flatshading=True,
                                 legendgroup=group, showlegend=panel == 0,
                                 lighting=dict(ambient=1., diffuse=.25, specular=.05, roughness=.95),
                                 hovertemplate="节点负荷 (%{x:.2f}, %{y:.2f}, %{z:.2f}) kW<extra>%{fullData.name}</extra>")
            else:
                trace = go.Scatter3d(x=[None], y=[None], z=[None], mode="markers", name=name,
                                    marker=dict(color=color, size=7), legendgroup=group,
                                    showlegend=panel == 0, hoverinfo="skip")
            fig.add_trace(trace, row=panel//2+1, col=panel % 2+1)
    axis = lambda node, bound: dict(title=dict(text=f"节点 {node} 负荷 (kW)", font=dict(size=11)),
                             range=[0, bound], nticks=5,
                             tickfont=dict(size=10), backgroundcolor="white", gridcolor="#E2E6E9",
                             zerolinecolor="#B8C1C6", showbackground=True)
    scenes = {"scene" if i == 0 else f"scene{i+1}": dict(
        xaxis=axis(result.load_nodes[0], result.bounds[0]),
        yaxis=axis(result.load_nodes[1], result.bounds[1]),
        zaxis=axis(result.load_nodes[2], result.bounds[2]), aspectmode="cube",
        camera=dict(eye=dict(x=1.5, y=1.6, z=1.2), projection=dict(type="orthographic")))
        for i in range(4)}
    buttons = [dict(label=name, method="update", args=[{
        "visible": [trace.legendgroup in groups for trace in fig.data]}])
        for name, groups in [("全部", {"common", "missed", "extra"}),
                             ("仅差异", {"missed", "extra"}),
                             ("计算域", {"common", "extra"}),
                             ("AC 参考域", {"common", "missed"})]]
    fig.update_layout(**scenes, height=900, template="plotly_white",
                      margin=dict(l=5, r=5, t=72, b=45),
                      font=dict(family="Arial, Microsoft YaHei, sans-serif", size=12),
                      legend=dict(orientation="h", x=.5, xanchor="center", y=-.04,
                                  groupclick="togglegroup"),
                      updatemenus=[dict(type="buttons", direction="right", buttons=buttons,
                                        x=.5, xanchor="center", y=1.09, yanchor="top")],
                      uirevision="method-comparison")
    fig.update_annotations(font=dict(size=13))
    return fig


def save_method_comparison(result, folder):
    import plotly.io as pio
    from vertify import METHODS, METHOD_NAMES

    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    filenames = ['region_comparison.html']+[
        f'methods_{"unlimited" if np.isinf(b) else int(b)}.html' for b in result.budgets[1:]]
    synchronize = r'''
const chart = document.getElementById('{plot_id}');
const scenes = ['scene', 'scene2', 'scene3', 'scene4'];
let synchronizing = false;
chart.on('plotly_relayout', event => {
  if (synchronizing) return;
  const key = Object.keys(event).find(k => /^scene[2-4]?\.camera$/.test(k));
  if (!key) return;
  synchronizing = true;
  const update = Object.fromEntries(scenes.filter(s => s+'.camera' !== key).map(s => [s+'.camera', event[key]]));
  Plotly.relayout(chart, update).finally(() => { synchronizing = false; });
});
'''
    rows = result.summary
    scope = '可规划域' if result.metadata['planning'] else '固定方案可调度域截面'
    spacing = ' × '.join(f'{value:g}' for value in result.spacing)
    cost_unit = result.metadata['cost_unit']
    for index, (budget, filename) in enumerate(zip(result.budgets, filenames)):
        links = ' · '.join(
            f'<a href="{name}" aria-current="{"page" if i == index else "false"}">'
            f'{"无限预算" if np.isinf(b) else f"{b:g} {cost_unit}"}</a>'
            for i, (b, name) in enumerate(zip(result.budgets, filenames)))
        navigation = f'<nav>预算：{links}</nav>' if result.metadata['planning'] else ''
        table = []
        for method, name in zip(METHODS, METHOD_NAMES):
            row = next(row for row in rows if row['method'] == method
                       and row['budget'] == (None if np.isinf(budget) else budget))
            rate = lambda key: '—' if row[key] is None else f'{row[key]:.4f}%'
            table.append(f'<tr><td>{name}</td><td>{rate("fr_percent")}</td>'
                         f'<td>{rate("mr_percent")}</td><td>{row["total_seconds"]:.3f}</td></tr>')
        chart = pio.to_html(plot_method_comparison(result, index), include_plotlyjs=True, full_html=False,
                           div_id='method-comparison', post_script=synchronize,
                           config=dict(responsive=True, displaylogo=False, scrollZoom=True))
        html = f'''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{result.metadata['network']} · {scope}</title><style>
body{{font-family:Arial,"Microsoft YaHei",sans-serif;color:#263641;margin:12px auto;max-width:1200px;padding:0 12px}}
nav{{text-align:center;padding:8px}}a{{color:#426985;text-decoration:none}}a[aria-current=page]{{font-weight:700;text-decoration:underline}}
table{{border-collapse:collapse;margin:10px auto;font-size:14px;min-width:520px}}td,th{{padding:7px 18px;text-align:right;border-bottom:1px solid #dce2e7}}td:first-child,th:first-child{{text-align:left}}
p{{font-size:12px;color:#56616a;line-height:1.7;text-align:center}}
</style></head><body>{navigation}{chart}
<table><thead><tr><th>方法</th><th>FR</th><th>MR</th><th>总计算时间（秒）</th></tr></thead><tbody>{''.join(table)}</tbody></table>
<p>图：{result.metadata['network']} 的{scope}；节点负荷单位为 kW，各面板使用同一组刻度。<br>
网格步长 {spacing} kW；蓝色为重合部分，AC 面板为完整参考域；红色遗漏，黄色多余。<br>
FR = 多余 / 计算域；MR = 遗漏 / AC 域。AC 自比较的零仅表示它是基准；有限网格上的零不表示连续误差严格为零。<br>
总时间含建模、求解、构域及网格判定，混合方法计入线性阶段，AC 扫描成本单列。绘图和导出不计。<br>
SOCP 两法展示切割外域，统一径向精度 {result.metadata['radial_tolerance']:g}。时间来自本次实验记录；AC 各预算为同一独立扫描的累计耗时。</p>
</body></html>'''
        (folder/filename).write_text(html, encoding='utf-8')
    return folder/filenames[0]
