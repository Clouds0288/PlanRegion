"""Notebook 内的离线交互回放：保存真实轮次，按割数复用几何，旋转不改变数据。"""
import base64
import gzip
import json
from pathlib import Path
from html import escape

import numpy as np
import shapely
from IPython.display import IFrame, display
from plotly.offline import get_plotlyjs
from shapely.geometry import Polygon

from node_power_view import region_frame, outer_frontier, _stair


def surface_mesh(faces):
    """在每个真实边界面内部三角化；非凸面也不使用跨面的凸包或三角扇。"""
    vertices, triangles, edges = [], [], []
    for face in faces:
        origin = face[0]
        u = face[1]-origin
        u /= np.linalg.norm(u)
        normal = np.linalg.svd(face-origin)[2][-1]
        basis = np.array([u,np.cross(normal,u)])
        # 外表面布尔运算会留下近重合角点；统一到 1e-8 kW 的显示精度后再三角化。
        polygon = shapely.set_precision(shapely.make_valid(Polygon((face-origin)@basis.T)),1e-8).simplify(1e-8)
        for triangle in shapely.constrained_delaunay_triangles(polygon).geoms:
            points = origin+np.asarray(triangle.exterior.coords)[:3]@basis
            triangles.append(list(range(len(vertices),len(vertices)+3)))
            vertices.extend(points.tolist())
        edges.extend([*face.tolist(),face[0].tolist(),None])
    return dict(vertices=vertices,triangles=triangles,edges=edges)


def replay_data(process):
    """几何快照与原始逐轮数值分开存；可行认证轮次仍保留完整 SP 证据。"""
    model = process.model
    frames = []
    for count in range(len(process.snapshots)):
        panels = region_frame(process,count)
        frontier = _stair(outer_frontier(process.designs,process.snapshots[count]),process.frontier_limit)
        frames.append(dict(frontier=frontier,panels=[dict(mesh=surface_mesh(p['faces']),
                       boundary=surface_mesh(p['boundary']),capacity=p['capacity']) for p in panels]))
    rounds = [dict(iteration=0,cuts=0,phase='初始外近似',status='initial')]
    for r in process.history:
        row = dict(iteration=r['global_iteration'],cuts=r['cuts_after'],phase=r['phase'],status=r['status'])
        if r['sp'] is not None:
            sp = r['sp']
            row.update(x=r['x'],p=r['p'],cost=r['cost_cny'],w=sp.w,eta=sp.violation,
                       pi=sp.cut.pi,rhs=sp.rhs,residual=sp.residual)
        rounds.append(row)
    names = [f'P_{e.name}' for e in model.corridors]+[f'v_{i}' for i in range(1,4)]
    lhs = [' '.join(f'{a:+.6g} {name}' for a,name in zip(row,names) if abs(a)>1e-12).lstrip('+') or '0'
           for row in model.W]
    data = dict(budgets=process.budgets,limit=process.limit,frontier_limit=process.frontier_limit,
                finished=process.finished,target=process.target,frames=frames,rounds=rounds,
                reference=_stair(process.reference,process.frontier_limit),
                cuts=[dict(constant=c.constant,x=c.x_coeff,p=c.load_coeff,eta=c.source_violation) for c in process.cuts],
                model=dict(states=names,rows=[r.ConstrName for r in model.electrical_constraints],lhs=lhs,
                           scales=model.scales,cost=model.cost,edges=[e.name for e in model.corridors],
                           lines=[line.name for line in model.lines]))
    return json.loads(json.dumps(data,default=lambda a:a.tolist()))


def replay_html(datasets, initial='first'):
    """Plotly 与数据一并内嵌；Notebook 和独立 HTML 共用同一个回放界面。"""
    payload = json.dumps(dict(datasets=datasets,initial=initial),ensure_ascii=False,separators=(',',':'),allow_nan=False)
    packed = base64.b64encode(gzip.compress(payload.encode('utf-8'))).decode('ascii')
    template = Path(__file__).with_name('node_power_replay.html').read_text(encoding='utf-8')
    return template.replace('__PLOTLY__',get_plotlyjs()).replace('__PAYLOAD__',packed)


def display_replay(process, saved_replay=None):
    """本次已求出的历史与已保存完整实验可切换；不会把未来割加入当前求解器。"""
    datasets = [dict(name='本次运行',data=replay_data(process))]
    if saved_replay is not None:
        data = json.loads(Path(saved_replay).read_text(encoding='utf-8'))
        datasets.append(dict(name='已保存的完整实验',data=data))
    html = replay_html(datasets,initial='last')
    display(IFrame('about:blank',width='100%',height=1420,
                   extras=['title="节点功率完整回放"',f'srcdoc="{escape(html,quote=True)}"']))
