"""有符号实验的离线网页回放；不调用求解器，不修改原生回放格式。"""
from copy import deepcopy
from itertools import product
import json
from pathlib import Path

import numpy as np
from scipy.spatial import ConvexHull
from plotly.offline import get_plotlyjs

from monitor import _merge, _voxel_faces


def mesh(points, d):
    """1. 单个凸域的面和棱；从不跨网架或符号区取凸包。"""
    points = np.asarray(points).reshape(-1, d)
    geometry = dict(points=points.tolist(), faces=[], edges=[])
    if len(points) < 2:
        return geometry
    delta = points-points[0]
    rank = np.linalg.matrix_rank(delta, tol=1e-9)
    basis = np.linalg.svd(delta, full_matrices=False)[2][:rank]
    coordinates = delta@basis.T
    if rank < 2:
        geometry['edges'] = [[int(coordinates[:, 0].argmin()), int(coordinates[:, 0].argmax())]]
        return geometry
    hull = ConvexHull(coordinates/np.linalg.norm(coordinates, axis=0))
    if rank == 2:
        order = hull.vertices.tolist()
        geometry['edges'] = [[a, b] for a, b in zip(order, order[1:]+order[:1])]
        geometry['faces'] = [[order[0], order[j], order[j+1]] for j in range(1, len(order)-1)]
        if d == 2:
            geometry['points'] = points[order].tolist()
    else:
        geometry['faces'] = hull.simplices.tolist()
        geometry['edges'] = [np.intersect1d(hull.simplices[i], hull.simplices[j]).tolist()
            for i, neighbors in enumerate(hull.neighbors) for j in neighbors
            if j > i and np.linalg.norm(hull.equations[i]-hull.equations[j]) > 1e-7]
    return geometry


def cut_mesh(cut, x, lower, upper):
    """割仅在所属符号区的盒内绘制。"""
    d = len(lower)
    cut = np.asarray(cut)
    beta, constant = cut[1:1+d], cut[0]+cut[1+d:]@x
    points = []
    for free in range(d):
        if beta[free] == 0.:
            continue
        fixed = [i for i in range(d) if i != free]
        for corner in product((0., 1.), repeat=d-1):
            point = lower.copy()
            point[fixed] += np.asarray(corner)*(upper-lower)[fixed]
            point[free] = 0.
            point[free] = -(constant+beta@point)/beta[free]
            if lower[free]-1e-8 <= point[free] <= upper[free]+1e-8:
                points.append(point)
    return mesh(np.unique(np.asarray(points).reshape(-1, d), axis=0), d)


def export_replay(data, path):
    """2. 为网页预计算单域网格；计算记录保持原始点和增量结构。"""
    path = Path(path)
    payload = deepcopy(data)
    d = len(payload['settings']['load_nodes'])
    state = {}
    for event in payload['history']:
        patch = event['patch']
        _merge(state, deepcopy(patch))
        for row in patch.get('schemes', {}).values():
            row['inner'], row['outer'] = mesh(row['inner'], d), mesh(row['outer'], d)
        if 'global_outer' in patch:
            patch['global_outer'] = [mesh(p, d) for p in patch['global_outer']]
        for row in patch.get('cut_history', {}).values():
            sign = np.asarray(state['sign'])
            lower = np.minimum(sign, 0)*state['bounds']
            upper = np.maximum(sign, 0)*state['bounds']
            row['geometry'] = cut_mesh(row['cut'], np.asarray(state['schemes'][row['scheme']]['x']), lower, upper)
    for key in ('inner', 'outer'):
        payload['result'][key] = [mesh(row['vertices'], d) for row in payload['result'][key]]
    reference = payload['reference']
    if d == 3 and reference is not None:
        lower, upper = np.asarray(reference['axis_lower']), np.asarray(reference['bounds'])
        faces = _voxel_faces(reference['states'], upper-lower)+lower
        reference['surface'] = dict(points=faces.reshape(-1, 3).tolist(),
            faces=np.concatenate([np.array([[0, 1, 2], [0, 2, 3]])+4*j for j in range(len(faces))]).tolist())
        del reference['states']
    library = path.parent/'plotly.min.js'
    library.write_text(get_plotlyjs(), encoding='utf-8')
    template = Path(__file__).with_suffix('.html').read_text(encoding='utf-8')
    path.write_text(template.replace('__DATA__', json.dumps(payload, ensure_ascii=False,
                    allow_nan=False, separators=(',', ':'))), encoding='utf-8')
    return path
