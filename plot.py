"""最终结果输出；不启动服务、不保存事件流，不在展示层进行优化求解。"""
import json
from pathlib import Path

import numpy as np
import plotly.graph_objects as go

from region import _convex_hull


def json_value(value):
    """JSON 只作必要类型转换；非有限界保留为 null，不能写成零。"""
    if isinstance(value, dict):
        return {str(k): json_value(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [json_value(v) for v in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, (int, np.integer)):
        return int(value)
    return value


def surface(vertices, name, color, *, opacity=.15, dashed=False):
    """二维多边形或三维凸节点；不跨方案/节点取凸包。"""
    vertices = np.asarray(vertices)
    if not len(vertices):
        return []
    dimension = vertices.shape[1]
    if dimension == 2:
        center = vertices.mean(axis=0)
        order = np.argsort(np.arctan2(vertices[:, 1]-center[1], vertices[:, 0]-center[0]))
        points = vertices[np.r_[order, order[:1]]]
        return [go.Scatter(x=points[:, 0], y=points[:, 1], name=name, mode='lines',
                           line=dict(color=color, dash='dash' if dashed else 'solid'),
                           fill='toself', opacity=opacity)]
    if np.linalg.matrix_rank(vertices-vertices[0], tol=1e-10) < 3:
        return [go.Scatter3d(x=vertices[:, 0], y=vertices[:, 1], z=vertices[:, 2],
                            name=name, mode='lines+markers', line=dict(color=color))]
    hull = _convex_hull(vertices)
    faces = hull.simplices
    mesh = go.Mesh3d(x=vertices[:, 0], y=vertices[:, 1], z=vertices[:, 2],
                     i=faces[:, 0], j=faces[:, 1], k=faces[:, 2], color=color,
                     opacity=opacity, name=name, showlegend=True, flatshading=True)
    # 同一平面被 Qhull 三角化；这里只画真实棱，避免满屏三角对角线。
    edges = {}
    for face, plane in zip(faces, hull.equations):
        for a, b in zip(face, np.roll(face, -1)):
            edges.setdefault(tuple(sorted((a, b))), []).append(plane[:-1])
    wire = []
    for (a, b), normals in edges.items():
        if len(normals) == 1 or np.linalg.norm(normals[0]-normals[1]) > 1e-6:
            wire.extend([vertices[a], vertices[b], [None]*3])
    line = np.asarray(wire, dtype=object).reshape(-1, 3)
    return [mesh, go.Scatter3d(x=line[:, 0], y=line[:, 1], z=line[:, 2], mode='lines',
                               line=dict(color=color, width=3, dash='dash' if dashed else 'solid'),
                               name=name, showlegend=False)]


def save_result(result, folder):
    """schema=4 的核心数据 + 离线可旋转图；green=O2，非解析真边界。"""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    (folder/'result.json').write_text(json.dumps(json_value(result), ensure_ascii=False, indent=2), encoding='utf-8')
    figure = go.Figure()
    layers = [('Stage 1 outer', '#2463eb', .10, [result['stage1_vertices']]),
              ('Stage 2 outer', '#079455', .18, [item['vertices'] for item in result['outer']]),
              ('Feasible inner', '#50b9a2', .30, [item['vertices'] for item in result['inner']])]
    for name, color, opacity, regions in layers:
        first = True
        for vertices in regions:
            traces = surface(vertices, name, color, opacity=opacity, dashed=True)
            for trace in traces:
                trace.legendgroup = name
                trace.showlegend = first
                figure.add_trace(trace)
                first = False
    labels = [f'p{node} (kW)' for node in result['load_nodes']]
    title = (f"Budget={result['budget']:g} | {result['status']} | "
             f"outer-to-inner gap ≤ {result['max_gap_kw']:.3f} kW | SOCP planning region")
    figure.update_layout(title=title, template='plotly_white', legend=dict(groupclick='togglegroup'))
    if len(labels) == 3:
        figure.update_layout(scene=dict(xaxis_title=labels[0], yaxis_title=labels[1],
                                        zaxis_title=labels[2], aspectmode='data'))
    else:
        figure.update_layout(xaxis_title=labels[0], yaxis_title=labels[1],
                              yaxis=dict(scaleanchor='x', scaleratio=1))
    figure.write_html(folder/'region.html', include_plotlyjs=True)
