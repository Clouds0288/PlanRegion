"""Benders 收敛图、费用前沿、三维并集表面及交互绘图。"""
from itertools import combinations

import matplotlib.pyplot as plt
import numpy as np
import shapely
from IPython.display import display
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import ConvexHull
from shapely.geometry import Polygon

from region import solid

BLUE, ORANGE, GRAY = "#48799A", "#C07B48", "#AAB7BE"
plt.rcParams.update({"font.family": "sans-serif", "font.size": 9,
                     "svg.fonttype": "none", "pdf.fonttype": 42,
                     "axes.spines.top": False, "axes.spines.right": False})


def plot_benders(solver):
    records = [r for r in solver.history if r["query"] == solver.query_id and r["x"] is not None]
    rounds = np.arange(1, len(records)+1)
    maximum = solver.queries[solver.query_id]["mode"] == "MP2"
    values = [sum(r["p"]) if maximum else solver.network.cost@r["x"] for r in records]
    fig, axes = plt.subplots(1, 2, figsize=(9, 3), layout="constrained")
    axes[0].plot(rounds, values, "o-", ms=3, color=BLUE,
                 label="Master upper bound" if maximum else "Master lower bound")
    if records and solver.finished and solver.feasible:
        axes[0].scatter(rounds[-1], values[-1], color=ORANGE, s=35, zorder=3, label="SP feasible optimum")
    axes[0].set(xlabel="Benders iteration", ylabel="Total load (kW)" if maximum else "Cost (CNY)")
    axes[0].legend(fontsize=8)
    axes[1].plot(rounds, [r["eta"] for r in records], "o-", ms=3, color=ORANGE)
    axes[1].axhline(0, color=GRAY, lw=.8)
    axes[1].set(xlabel="Benders iteration", ylabel=r"SP violation $\eta$")
    for ax in axes:
        ax.grid(alpha=.2)
    state = "optimal" if solver.finished and solver.feasible else "infeasible" if solver.finished else "in progress"
    fig.suptitle(f'{solver.queries[solver.query_id]["mode"]} + SP | {state}')
    return fig


def show_benders(solver):
    fig = plot_benders(solver)
    display(fig)
    plt.close(fig)


def certified_frame(polytopes):
    points = [p.tolist() for poly in polytopes.values() if not solid(poly) for p in poly]
    return dict(mesh=surface_mesh(union_surface(list(polytopes.values()))), points=points)


def outer_frontier(costs, polytopes):
    """固定设计逐一投影，保留费用—总负荷的整数台阶。"""
    xx, yy, previous = [], [], 0.
    for cost, poly in zip(costs, polytopes):
        if not len(poly):
            continue
        capacity = float(poly.sum(axis=1).max())
        if capacity > previous+1e-7 or not xx:
            xx.extend([previous, capacity])
            yy.extend([float(cost), float(cost)])
            previous = capacity
    return [xx, yy]


def maximal_polytopes(polytopes):
    ranked = sorted(((ConvexHull(p).volume, p) for p in polytopes if solid(p)), key=lambda item: -item[0])
    kept, equations = [], []
    for _, points in ranked:
        if any(np.all(points@eq[:, :3].T+eq[:, 3] <= 2e-7) for eq in equations):
            continue
        kept.append(points)
        equations.append(ConvexHull(points).equations)
    return kept


def union_surface(polytopes):
    """扣除内部覆盖面，绘制固定设计多面体的非凸并集。"""
    polys = maximal_polytopes(polytopes)
    hulls = [ConvexHull(p) for p in polys]
    surfaces = []
    for i, (points, hull) in enumerate(zip(polys, hulls)):
        planes = {}
        for face, equation in zip(hull.simplices, hull.equations):
            planes.setdefault(tuple(np.round(equation, 8)), (equation, set()))[1].update(face)
        for equation, indices in planes.values():
            normal, offset = equation[:3], equation[3]
            u = np.cross(normal, [1., 0., 0.] if abs(normal[0]) < .9 else [0., 1., 0.])
            u /= np.linalg.norm(u)
            basis = np.array([u, np.cross(normal, u)])
            origin = points[next(iter(indices))]
            xy = (points[list(indices)]-origin)@basis.T
            visible = Polygon(xy[ConvexHull(xy).vertices])
            for j, (other, other_hull) in enumerate(zip(polys, hulls)):
                if i == j:
                    continue
                signed = other@normal+offset
                if signed.max() < -1e-7 or signed.min() > 1e-7 or (signed.max() <= 1e-7 and j > i):
                    continue
                section = list(other[np.abs(signed) <= 1e-7])
                edges = {tuple(sorted(pair)) for face in other_hull.simplices
                         for pair in combinations(face, 2)}
                for a, b in edges:
                    if signed[a]*signed[b] < 0:
                        section.append(other[a]+(other[b]-other[a])*signed[a]/(signed[a]-signed[b]))
                if len(section) < 3:
                    continue
                projected = np.unique(np.round((np.array(section)-origin)@basis.T, 8), axis=0)
                if len(projected) < 3 or np.linalg.matrix_rank(projected-projected[0], tol=1e-7) < 2:
                    continue
                visible = visible.difference(Polygon(projected[ConvexHull(projected).vertices]))
                if visible.is_empty:
                    break
            for piece in shapely.get_parts(visible):
                if piece.geom_type != "Polygon" or piece.area < 1e-8:
                    continue
                patches = list(shapely.constrained_delaunay_triangles(piece).geoms) if piece.interiors else [piece]
                for patch in patches:
                    surfaces.append(origin+np.asarray(patch.exterior.coords)[:-1]@basis)
    return surfaces


def surface_mesh(faces):
    vertices, triangles = [], []
    for face in faces:
        origin = face[0]
        basis = np.linalg.svd(face-origin, full_matrices=False)[2][:2]
        polygon = shapely.set_precision(shapely.make_valid(Polygon((face-origin)@basis.T)), 1e-8).simplify(1e-8)
        for triangle in shapely.constrained_delaunay_triangles(polygon).geoms:
            points = origin+np.asarray(triangle.exterior.coords)[:3]@basis
            triangles.append(list(range(len(vertices), len(vertices)+3)))
            vertices.extend(points.tolist())
    return dict(vertices=vertices, triangles=triangles)


def frame_data(costs, designs, polytopes, cut=None):
    faces = union_surface(polytopes)
    boundary = []
    if cut is not None:
        nx = designs.shape[1]
        offsets = cut[0]+designs@cut[1:1+nx]
        boundary = [face for face in faces
                    if np.min(np.max(np.abs(offsets[:, None]+face@cut[1+nx:]), axis=1)) < 2e-7]
    return dict(frontier=outer_frontier(costs, polytopes), mesh=surface_mesh(faces),
                boundary=surface_mesh(boundary))


def plot_frontier(process):
    fig, ax = plt.subplots(figsize=(7, 3), layout="constrained")
    xx, yy = outer_frontier(process.costs, process.polytopes)
    ax.plot(xx, yy, color=GRAY, label="Candidate outer bound")
    points = process.frontier
    if points:
        ax.scatter([p[0] for p in points], [p[1] for p in points], color=ORANGE, s=22,
                   label="Certified frontier points", zorder=3)
    limit = process.network.power_limit
    margin = .02*limit
    ax.set(xlim=(-margin, limit+margin),
           xlabel="Total load (kW)", ylabel="Minimum budget (CNY)")
    ax.grid(alpha=.2)
    ax.legend(fontsize=8)
    return fig


def plot_regions(process, outer=False):
    fig = plt.figure(figsize=(7, 6))
    fig.subplots_adjust(left=.02, right=.88, bottom=.16, top=.94)
    ax = fig.add_subplot(projection="3d")
    ax.add_collection3d(Poly3DCollection(union_surface(list(process.certified.values())),
                       facecolor=BLUE, edgecolor=BLUE, alpha=.3, linewidth=.3), autolim=False)
    if outer:
        ax.add_collection3d(Poly3DCollection(union_surface(process.polytopes),
                           facecolor=GRAY, edgecolor=GRAY, alpha=.1, linewidth=.2), autolim=False)
    for poly in process.certified.values():
        if not solid(poly):
            ax.scatter(*poly.T, color=BLUE, s=10)
    limit = process.network.power_limit
    labels = [f"Node {node} load (kW)" for node in process.network.load_nodes]
    ax.set(xlim=(0, limit), ylim=(0, limit), zlim=(0, limit),
           xlabel=labels[0], ylabel=labels[1], zlabel=labels[2])
    ax.set_box_aspect((1, 1, 1))
    ax.set_title("Certified planning domain" if process.finished else "Certified inner approximation")
    return fig


# 交互绘图也集中在本模块；replay.py 只负责数据打包与播放控制。
PLOT_SCRIPT = r"""
const blue='#48799A', orange='#C07B48', gray='#AAB7BE';
const frontierCeiling=1.05*Math.max(1,...experiment.scenarios.flatMap(s=>s.frames.flatMap(f=>f.frontier[1])));
const meshTrace=(mesh,color,opacity)=>({type:'mesh3d',
 x:mesh.vertices.map(v=>v[0]),y:mesh.vertices.map(v=>v[1]),z:mesh.vertices.map(v=>v[2]),
 i:mesh.triangles.map(t=>t[0]),j:mesh.triangles.map(t=>t[1]),k:mesh.triangles.map(t=>t[2]),
 color,opacity,flatshading:true,hoverinfo:'skip',showscale:false});
function drawFrame(frame, previous, certified, record, points) {
 const frontier=[{x:frame.frontier[0],y:frame.frontier[1],type:'scatter',mode:'lines',
   line:{color:gray,width:2},name:'候选能力上界'},
  {x:points.map(p=>p[0]),y:points.map(p=>p[1]),type:'scatter',mode:'markers',
   marker:{color:orange,size:7},name:'已认证前沿点'}];
 Plotly.react('frontier',frontier,{margin:{l:75,r:25,t:15,b:50},
   xaxis:{title:{text:'总负荷 (kW)'},range:[0,experiment.limit]},
   yaxis:{title:{text:'最低预算 (元)'},range:[0,frontierCeiling]},
   legend:{orientation:'h',y:1.15},uirevision:'frontier'}, {responsive:true,displaylogo:false});
 const traces=[{...meshTrace(certified.mesh,blue,.55),name:'已认证可行区域'}];
 if(certified.points.length)traces.push({type:'scatter3d',mode:'markers',
   x:certified.points.map(p=>p[0]),y:certified.points.map(p=>p[1]),z:certified.points.map(p=>p[2]),
   marker:{color:blue,size:3},name:'已认证可行点'});
 if(document.getElementById('outer').checked){
   traces.unshift(meshTrace(frame.mesh,gray,.13));
   traces.push(meshTrace(frame.boundary,orange,.7));
 }
 if(document.getElementById('before').checked) traces.unshift(meshTrace(previous.mesh,gray,.08));
 if(record && record.p) traces.push({type:'scatter3d',mode:'markers',
   x:[record.p[0]],y:[record.p[1]],z:[record.p[2]],marker:{color:orange,size:4},
   name:record.cut===null?'本轮方案 SP 可行':'本轮方案 SP 不满足',
   hovertemplate:'节点负荷 (%{x:.2f}, %{y:.2f}, %{z:.2f}) kW<extra>%{fullData.name}</extra>'});
 const axis=k=>({title:{text:`节点 ${experiment.load_nodes[k]} 负荷 (kW)`},
   range:[0,experiment.limit],nticks:5});
 Plotly.react('region',traces,{margin:{l:0,r:0,t:5,b:0},showlegend:false,
   scene:{xaxis:axis(0),yaxis:axis(1),zaxis:axis(2),
   aspectmode:'cube'},uirevision:'keep-camera'}, {responsive:true,displaylogo:false});
}
"""
