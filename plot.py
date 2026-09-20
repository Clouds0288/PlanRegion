"""区域对比与固定方案切割回放；只消费计算结果。"""
from pathlib import Path  # 定位可再生成的 HTML 输出文件。
import numpy as np  # 处理绘图坐标和体素标签。
import plotly.graph_objects as go  # 创建可旋转的三维图与切割动画。
from region import clip_polytope, polytope_vertices  # 按原始 SP 割和内点重建几何，不额外调用模型。


def save_replay(result, folder, method='socp', budget_index=-1, design_index=-1):  # 回放一个固定方案的真实切割历史。
    record = result.regions[method][budget_index][design_index]  # 选择同一方法、预算及建设方案的唯一记录。
    outer = np.asarray(record['initial'])  # 候选外域从该次实验实际采用的初始域开始。
    inner = np.zeros((1,3))  # 参数原点是该固定方案已认证的初始内点。
    frames = []  # 每个动画帧对应实际 SP 查询。

    def traces(point):  # 用当前内外域和查询点生成一帧的图层。
        mesh = lambda p,c,o,n: go.Mesh3d(x=p[:,0],y=p[:,1],z=p[:,2],alphahull=0,  # 按真实凸包顶点绘制三维表面。
                                       color=c,opacity=o,name=n,flatshading=True)  # 设置颜色、透明度及图例名称。
        return [mesh(outer,'#AAB7BE',.2,'候选外域'), mesh(inner,'#4286AD',.7,'认证内域'),  # 灰色展示候选外域，蓝色展示认证内域。
                go.Scatter3d(x=[point[0]],y=[point[1]],z=[point[2]],mode='markers',  # 单独突出本轮真正接受 SP 查询的负荷点。
                             marker=dict(color='#C07B48',size=4),name='本轮查询点')]  # 查询点使用统一标记，不添加额外表面采样。

    # 先用三条轴向查询建立三维内域；此后每一帧对应真实的一次 SP 查询。
    for i,row in enumerate(record['history']):  # 按求解发生顺序处理历史记录。
        if row['cut'] is not None:  # 只有已产生有效割的轮次才修改外域。
            cut = np.asarray(row['cut'])  # 读取当轮保存的负荷割系数。
            outer = clip_polytope(outer,cut[0],cut[1:])  # 以 a+bᵀp≥0 重建当轮候选域。
        inner = polytope_vertices(np.vstack([inner,row['witness']]))  # 只在这个固定方案内对认证内点取凸包。
        if i >= 2:  # 前三次轴向查询建立三维内域后再开始展示。
            frames.append(go.Frame(name=str(i+1),data=traces(row['p'])))  # 帧编号对应原 SP 历史中的查询轮次。
    fig = go.Figure(data=frames[0].data,frames=frames)  # 以首个完整三维帧初始化动画。
    axis = lambda k: dict(title=f'节点 {result.load_nodes[k]} 负荷 (kW)',range=[0,result.bounds[k]])  # 坐标轴使用实际负荷节点号和公共评价范围。
    x = result.metadata['designs'][record['design']]['x']  # 建设向量从唯一方案表读取。
    fig.update_layout(template='plotly_white',height=670,margin=dict(l=5,r=5,t=45,b=90),  # 采用简洁白底和紧凑留白。
        scene=dict(xaxis=axis(0),yaxis=axis(1),zaxis=axis(2),aspectmode='cube'),  # 三个坐标轴采用统一范围和立方体比例。
        annotations=[dict(text=f'固定建设方案 x={x}；灰色外域，蓝色认证内域',x=.5,y=1.04,  # 标明当前固定建设方案及内外域颜色。
                          xref='paper',yref='paper',showarrow=False)],  # 文字放在图面上方，不添加装饰性箭头。
        sliders=[dict(currentvalue=dict(prefix='SP 查询轮次：'),steps=[  # 滑块按实际 SP 查询轮次切换。
            dict(method='animate',label=f.name,args=[[f.name],dict(mode='immediate',  # 选择轮次后立即跳转到对应帧。
                 frame=dict(duration=0,redraw=True),transition=dict(duration=0))]) for f in frames])],  # 取消插帧，避免把动画过渡误认为求解轨迹。
        updatemenus=[dict(type='buttons',direction='left',x=0,y=-.06,buttons=[  # 添加紧凑的播放和暂停控制。
            dict(label='播放',method='animate',args=[None,dict(fromcurrent=True,  # 播放从当前查看的轮次开始。
                frame=dict(duration=200,redraw=True),transition=dict(duration=0))]),  # 仅按已保存的 SP 帧顺序播放。
            dict(label='暂停',method='animate',args=[[None],dict(mode='immediate',  # 暂停保留当前查询轮次和视角。
                frame=dict(duration=0,redraw=False))])])])  # 停止时不重绘无关帧。
    path = Path(folder)/'cutting_process.html'  # 回放文件是可再生成展示，不是第二份实验数据。
    fig.write_html(path,include_plotlyjs=True,auto_play=False,config=dict(displaylogo=False))  # 内嵌绘图库，页面可离线旋转和回放。
    return path  # 返回 Notebook 可嵌入的本地页面路径。


def voxel_surface(mask, spacing):  # 把体素集合转换为保留孔洞的外表面三角网格。
    """Exposed cell faces, merged into rectangles without filling holes.

    Geometry and disagreement percentages use the very same labeled cells.
    No global convex hull or smoothing changes the meaning of a colored region.
    """
    mask = np.asarray(mask, dtype=bool)  # 统一按体素是否属于目标集合处理。
    vertices, triangles = [], []  # 分别收集表面顶点和三角形索引。
    for axis in range(3):  # 对三个坐标方向分别寻找暴露面。
        others = [i for i in range(3) if i != axis]  # 其余两轴构成当前面的二维坐标。
        oriented = np.moveaxis(mask, axis, 0)  # 将处理方向移到第一个轴，统一后续算法。
        for side in (-1, 1):  # 分别提取负向和正向的外露面。
            neighbour = np.zeros_like(oriented)  # 评价箱外侧视为空体素。
            if side == 1:  # 正方向暴露面由下一层体素决定。
                neighbour[:-1] = oriented[1:]  # 把正向邻居对齐到当前层。
            else:  # 负方向使用上一层体素作为邻居。
                neighbour[1:] = oriented[:-1]  # 把负向邻居对齐到当前层。
            faces = oriented & ~neighbour  # 当前体素存在而邻居为空时才产生表面。
            for plane in range(len(faces)):  # 逐层合并位于同一平面的暴露单元。
                cells = faces[plane].copy()  # 复制当前面掩码，用于标记已被矩形覆盖的单元。
                for row in range(cells.shape[0]):  # 从上到下处理二维单元行。
                    while cells[row].any():  # 本行还有未处理表面单元时继续。
                        left = int(np.flatnonzero(cells[row])[0])  # 找到当前行最左侧剩余单元。
                        right = left+1  # 矩形初始宽度为一个单元。
                        while right < cells.shape[1] and cells[row, right]:  # 向右合并连续且属于同一表面的单元。
                            right += 1  # 扩大当前矩形宽度。
                        bottom = row+1  # 矩形初始高度为一行。
                        while bottom < cells.shape[0] and cells[bottom, left:right].all():  # 后续整行都被表面占据时向下合并。
                            bottom += 1  # 扩大当前矩形高度。
                        cells[row:bottom, left:right] = False  # 标记这块矩形已处理，防止重复画面。
                        face = np.zeros((4, 3), dtype=float)  # 每个合并矩形生成四个三维顶点。
                        face[:, axis] = plane+(side == 1)  # 正向面在体素上边界，负向面在下边界。
                        face[:, others[0]] = [row, bottom, bottom, row]  # 填入矩形沿第一个面内坐标的范围。
                        face[:, others[1]] = [left, left, right, right]  # 填入矩形沿第二个面内坐标的范围。
                        normal = np.cross(face[1]-face[0], face[2]-face[0])  # 计算当前顶点顺序的法向量。
                        if normal[axis]*side < 0:  # 法向量朝向体素内部时需要反转顺序。
                            face = face[::-1]  # 统一使表面法向朝外，保留内部孔洞。
                        first = len(vertices)  # 记住本矩形在总顶点表中的起始索引。
                        vertices.extend((face*spacing).tolist())  # 体素索引乘以网格步长，恢复 kW 坐标。
                        triangles.extend([[first, first+1, first+2], [first, first+2, first+3]])  # 用两个三角形覆盖当前矩形。
    return dict(vertices=vertices, triangles=triangles)  # 返回明确的表面网格，不用跨孔洞凸包替代。


def plot_method_comparison(result, budget_index=0):  # 在同一预算下比较三种计算域与独立 AC 参考。
    """同一预算的四方法对比；曲面与 FR/MR 使用同一份网格标签。"""
    from plotly.subplots import make_subplots  # 建立紧凑的四面板布局。
    from vertify import METHOD_NAMES  # 显示名称与结果数组的方法轴一致。

    fig = make_subplots(rows=2, cols=2, specs=[[{"type": "scene"}]*2]*2,  # 四个面板均为可旋转的三维场景。
                        subplot_titles=[f"{letter}  {name}" for letter, name in zip("abcd", METHOD_NAMES)],  # 使用简洁的 a、b、c、d 面板标记。
                        horizontal_spacing=.02, vertical_spacing=.08)  # 统一控制面板间距。
    categories = [(3, "与 AC 重合", "#4286AD", 1., "common"),  # 重合部分使用克制的蓝色。
                  (2, "遗漏", "#D43D3D", 1., "missed"),  # 遗漏的 AC 可行区域使用红色。
                  (1, "多余", "#E9B72F", 1., "extra")]  # 多余的计算区域使用黄色。
    for panel, labels in enumerate(result.labels[:, budget_index]):  # 读取当前预算下每种方法的差集标签。
        for code, name, color, opacity, group in categories:  # 每个面板分别生成重合、遗漏和多余表面。
            mesh = voxel_surface(labels == code, result.spacing)  # 使用同一套体素表面算法避免不同颜色口径不一致。
            points = np.asarray(mesh["vertices"]).reshape(-1, 3)  # 将输出顶点转为标准三列坐标数组。
            faces = np.asarray(mesh["triangles"], dtype=int).reshape(-1, 3)  # 将三角形转为三列整数索引数组。
            if len(points):  # 类别非空时添加实际表面。
                trace = go.Mesh3d(x=points[:, 0], y=points[:, 1], z=points[:, 2],  # 三维网格坐标使用实际 kW 值。
                                 i=faces[:, 0], j=faces[:, 1], k=faces[:, 2],  # 显式指定三角形，保留非凸区域及孔洞。
                                 name="AC 参考域" if panel == 2 else name,  # AC 自比较面板直接标为参考域。
                                 color=color, opacity=opacity, flatshading=True,  # 各类别使用固定颜色与不透明表面。
                                 legendgroup=group, showlegend=panel == 0,  # 图例只显示一次，按类别控制四面板。
                                 lighting=dict(ambient=1., diffuse=.25, specular=.05, roughness=.95),  # 采用均匀照明以减少颜色识别偏差。
                                 hovertemplate="节点负荷 (%{x:.2f}, %{y:.2f}, %{z:.2f}) kW<extra>%{fullData.name}</extra>")  # 悬停仅显示用户关心的节点负荷和区域类别。
            else:  # 空类别保留图例，避免因某档预算无差异而改变说明。
                trace = go.Scatter3d(x=[None], y=[None], z=[None], mode="markers", name=name,  # 空图层不产生实际数据点。
                                    marker=dict(color=color, size=7), legendgroup=group,  # 沿用该类别的颜色和图例分组。
                                    showlegend=panel == 0, hoverinfo="skip")  # 只在首面板显示图例，关闭空图层悬停。
            fig.add_trace(trace, row=panel//2+1, col=panel % 2+1)  # 按方法顺序填入对应行列。
    axis = lambda node, bound: dict(title=dict(text=f"节点 {node} 负荷 (kW)", font=dict(size=11)),  # 坐标标题使用真实节点号及 kW 单位。
                             range=[0, bound], nticks=5,  # 所有方法采用相同上界和刻度密度。
                             tickfont=dict(size=10), backgroundcolor="white", gridcolor="#E2E6E9",  # 使用白底与浅灰网格保持论文插图风格。
                             zerolinecolor="#B8C1C6", showbackground=True)  # 零线稍作区分，避免装饰性背景。
    scenes = {"scene" if i == 0 else f"scene{i+1}": dict(  # 为四个三维面板建立相同坐标和相机设置。
        xaxis=axis(result.load_nodes[0], result.bounds[0]),  # 第一坐标对应第一个实际独立负荷节点。
        yaxis=axis(result.load_nodes[1], result.bounds[1]),  # 第二坐标对应第二个实际独立负荷节点。
        zaxis=axis(result.load_nodes[2], result.bounds[2]), aspectmode="cube",  # 第三坐标对应第三个节点，统一场景比例。
        camera=dict(eye=dict(x=1.5, y=1.6, z=1.2), projection=dict(type="orthographic")))  # 正交投影便于比较区域形状与尺度。
        for i in range(4)}  # 四个场景共用同样的观察方向。
    buttons = [dict(label=name, method="update", args=[{  # 显示控制只切换已有类别图层。
        "visible": [trace.legendgroup in groups for trace in fig.data]}])  # 按图例分组决定图层可见性。
        for name, groups in [("全部", {"common", "missed", "extra"}),  # 默认同时显示重合、多余及遗漏。
                             ("仅差异", {"missed", "extra"}),  # 可只观察红黄差集。
                             ("计算域", {"common", "extra"}),  # 计算域由重合与多余部分组成。
                             ("AC 参考域", {"common", "missed"})]]  # AC 域由重合与遗漏部分组成。
    fig.update_layout(**scenes, height=900, template="plotly_white",  # 应用四场景设置和统一白底模板。
                      margin=dict(l=5, r=5, t=72, b=45),  # 保持紧凑留白，将长说明放在图注。
                      font=dict(family="Arial, Microsoft YaHei, sans-serif", size=12),  # 采用可编辑的常规中英文字体。
                      legend=dict(orientation="h", x=.5, xanchor="center", y=-.04,  # 图例横向居中放置。
                                  groupclick="togglegroup"),  # 同类别在四个面板中同时切换。
                      updatemenus=[dict(type="buttons", direction="right", buttons=buttons,  # 用简短按钮切换需要比较的集合。
                                        x=.5, xanchor="center", y=1.09, yanchor="top")],  # 控制栏居中，不占用主要绘图区域。
                      uirevision="method-comparison")  # 切换图层时保留用户旋转后的视角。
    fig.update_annotations(font=dict(size=13))  # 统一面板标记字号。
    return fig  # 返回可由 Notebook 展示或导出的图对象。


def save_method_comparison(result, folder):  # 为各预算导出同一原始结果驱动的交互页面。
    import plotly.io as pio  # 把 Plotly 图对象转为 HTML。
    from vertify import METHODS, METHOD_NAMES  # 读取统一的方法顺序及名称。

    folder = Path(folder)  # 使用路径对象组织展示文件。
    folder.mkdir(parents=True, exist_ok=True)  # 创建当前实验输出目录。
    filenames = ['region_comparison.html']+[  # 第一档预算作为入口页面。
        f'methods_{"unlimited" if np.isinf(b) else int(b)}.html' for b in result.budgets[1:]]  # 后续页面按预算命名，无限预算单独标识。
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
    rows = result.summary  # FR/MR 和时间全部由结果容器现场推导。
    scope = '可规划域' if result.metadata['planning'] else '固定方案可调度域截面'  # 根据实验配置区分规划域与固定网架截面。
    spacing = ' × '.join(f'{value:g}' for value in result.spacing)  # 图注标明三个负荷轴的网格步长。
    cost_unit = result.metadata['cost_unit']  # 费用单位来自唯一网架配置。
    for index, (budget, filename) in enumerate(zip(result.budgets, filenames)):  # 逐预算生成相同布局的比较页面。
        links = ' · '.join(  # 构造各预算之间的导航链接。
            f'<a href="{name}" aria-current="{"page" if i == index else "false"}">'  # 为当前页面标记选中状态。
            f'{"无限预算" if np.isinf(b) else f"{b:g} {cost_unit}"}</a>'  # 预算显示使用原始投资单位。
            for i, (b, name) in enumerate(zip(result.budgets, filenames)))  # 导航顺序与实验预算顺序一致。
        navigation = f'<nav>预算：{links}</nav>' if result.metadata['planning'] else ''  # 固定方案实验不需要预算导航。
        table = []  # 收集当前预算四种方法的指标行。
        for method, name in zip(METHODS, METHOD_NAMES):  # 按统一的方法顺序输出比较表。
            row = next(row for row in rows if row['method'] == method  # 定位该方法在当前预算下的派生指标。
                       and row['budget'] == (None if np.isinf(budget) else budget))  # 无限预算的存储形式为 JSON null。
            rate = lambda key: '—' if row[key] is None else f'{row[key]:.4f}%'  # 空域的条件比例没有定义，用横线表示。
            table.append(f'<tr><td>{name}</td><td>{rate("fr_percent")}</td>'  # 输出方法名称和 FR。
                         f'<td>{rate("mr_percent")}</td><td>{row["total_seconds"]:.3f}</td></tr>')  # 输出 MR 和本次实际总计算时间。
        chart = pio.to_html(plot_method_comparison(result, index), include_plotlyjs=True, full_html=False,  # 图形内嵌运行库，离线即可旋转和筛选差集。
                           div_id='method-comparison', post_script=synchronize,  # 同步四个场景的相机，便于同视角比较。
                           config=dict(responsive=True, displaylogo=False, scrollZoom=True))  # 支持窗口尺寸变化与滚轮缩放。
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
        (folder/filename).write_text(html, encoding='utf-8')  # 写入可再生成的展示页面，原始结果仍只存一次。
    return folder/filenames[0]  # 返回默认预算页面作为 Notebook 入口。
