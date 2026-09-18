"""NodePower 的逐轮教学与图件；主问题和电气方程集中在 planning_domain_demo。

图1：全部合法建设方案的割外近似，投影到总功率—预算平面；独立 LP 给灰色基准。
图2：比较不同预算下的三维并集；各面板保留同一物理上限，不对不同方案取整体凸包。
域证书阶段固定一个合法建设方案，用 MP1/SP 检查其外近似顶点。全部顶点可行
才证明该方案的整个凸多面体可行；所有预算内方案完成后得到完整的并集证书。
"""
from __future__ import annotations

import itertools
import json
import math
from html import escape
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shapely
from IPython.display import HTML, display
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection
from scipy.spatial import ConvexHull
from shapely.geometry import Polygon

import planning_domain_demo as demo


BLUE, ORANGE, GRAY = "#48799A", "#C07B48", "#AAB7BE"
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 8, "axes.labelsize": 9, "legend.fontsize": 8,
    "svg.fonttype": "none", "pdf.fonttype": 42,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": .8, "legend.frameon": False,
})


def simplex(limit):
    """零基准、总功率上限 L 的初始四面体 p≥0，Σp≤L。"""
    return np.vstack([np.zeros(3), float(limit)*np.eye(3)])


def clip_polytope(vertices, constant, coefficient):
    """用 constant+coefficient@p≥0 裁剪凸多面体；只做边与平面的交点。"""
    values = constant + vertices@coefficient
    if np.all(values >= -1e-9):
        return vertices.copy()
    if np.all(values < 0):
        return np.empty((0, 3))
    edges = {tuple(sorted(pair)) for face in ConvexHull(vertices).simplices
             for pair in itertools.combinations(face, 2)}
    points = [p for p, value in zip(vertices, values) if value >= 0]
    for a, b in edges:
        if values[a]*values[b] < 0:
            points.append(vertices[a] + (vertices[b]-vertices[a])*values[a]/(values[a]-values[b]))
    points = np.unique(np.round(points, 10), axis=0)
    return points[ConvexHull(points).vertices]


def cut_polytopes(designs, polytopes, cut):
    return [clip_polytope(poly, cut.constant+cut.x_coeff@design.x, cut.load_coeff)
            for design, poly in zip(designs, polytopes)]


def outer_frontier(designs, polytopes):
    """对每个固定 x 精确取多面体顶点最大总功率，投影仍保留整数设计。"""
    rows, best = [], -1.0
    for design, poly in zip(designs, polytopes):
        capacity = float(np.max(poly.sum(axis=1)))
        if capacity > best+1e-7:
            row = dict(cost_cny=design.cost_cny, total_kw=capacity)
            if rows and rows[-1]["cost_cny"] == design.cost_cny:
                rows[-1] = row
            else:
                rows.append(row)
            best = capacity
    return pd.DataFrame(rows)


def maximal_polytopes(polytopes):
    """绘图时去掉被单个其他多面体完全包含的成员；保持并集不变。"""
    ranked = sorted(((ConvexHull(p).volume, p) for p in polytopes), key=lambda item: -item[0])
    kept, equations = [], []
    for _, points in ranked:
        if any(np.all(points@eq[:, :3].T + eq[:, 3] <= 2e-7) for eq in equations):
            continue
        kept.append(points)
        equations.append(ConvexHull(points).equations)
    return kept


class PowerProcess:
    """保留会话、累计割、各轮原始状态及多面体；重复运行块5不会重置。"""

    def __init__(self, model, target=(40., 35., 25.), budgets=(20000., 30000., 40000., None), frontier_limit=140.):
        self.model, self.target = model, np.asarray(target)
        self.budgets = tuple(budgets)               # None 只取消建设预算；电气约束全部保留。
        self.limit = model.config.power_limit       # 三维域保留完整配变有功上限 142.5 kW。
        self.frontier_limit = float(frontier_limit) # 二维图延续 0—140 kW 的研究区间。
        self.designs = demo.enumerate_designs(model)  # 只用结构与费用选择域证书候选。
        self.budget_designs = [[i for i,z in enumerate(self.designs) if budget is None or z.cost_cny<=budget]
                               for budget in self.budgets]
        self.affordable = sorted(set().union(*map(set,self.budget_designs))) # 各面板所需方案共同认证。
        self.reference = demo.frontier_table(model, self.designs)  # 独立 LP 只供图1基准和最终核验。
        self.polytopes = [simplex(self.limit) for _ in self.designs]
        self.snapshots = [self.polytopes]             # 下标就是割数；无新增割时复用同一几何快照。
        self.region_frames = {}                    # 每个割数只计算一次外表面，供静态图和交互回放共用。
        self.cuts, self.history, self.frontier = [], [], []
        self.checked = {(i, (0., 0., 0.)) for i in self.affordable}  # 原点 P=0、v=1 的解析可行证书。
        self.session = demo.CutSession(model, "min_cost", self.target, self.cuts)
        self.phase, self.finished = "MP1 固定节点功率", False
        self.frontier_budget = float(model.cost.sum())
        self.quantum = math.gcd(*np.rint(model.cost[model.cost > 0]).astype(int))
        self.region_index = None

    def _start(self, phase, mode, query, **options):
        self.phase = phase
        self.session = demo.CutSession(self.model, mode, query, self.cuts, **options)

    def _region_candidate(self):
        """从当前外近似选未认证顶点；不读取独立物理域或其约束。"""
        for i in self.affordable:
            for p in sorted(self.polytopes[i], key=lambda z: -z.sum()):
                key = (i, tuple(np.round(p, 8)))
                if key not in self.checked:
                    self.region_index = i
                    self.region_key = key          # 证书对应原始候选顶点，避免求解器尾数改变顶点身份。
                    self._start("三维域：固定设计顶点 MP1", "min_cost", p,
                                fixed_x=self.designs[i].x, load_limit=self.limit)
                    return
        self.finished = True
        self.phase = "完成：MP1、二维前沿、预算内全部三维域顶点"

    def _advance_phase(self):
        """数学查询的串联：点查询→MP2/总量MP1前沿→逐方案域证书。"""
        result = self.session.result()
        self.session.master.dispose()
        if self.phase == "MP1 固定节点功率":
            self._start("二维前沿 MP2", "max_load", self.frontier_budget, load_limit=self.frontier_limit)
        elif self.phase == "二维前沿 MP2":
            self.capacity = result.value
            self._start("二维前沿 MP1 总量", "min_cost_total", self.capacity, load_limit=self.frontier_limit)
        elif self.phase == "二维前沿 MP1 总量":
            self.frontier.append(dict(cost_cny=result.cost_cny, total_kw=self.capacity,
                                      p1=result.p[0], p2=result.p[1], p3=result.p[2]))
            self.frontier_budget = round(result.cost_cny)-self.quantum
            if self.frontier_budget >= 0:
                self._start("二维前沿 MP2", "max_load", self.frontier_budget, load_limit=self.frontier_limit)
            else:
                self._region_candidate()
        else:
            self._region_candidate()

    def step(self):
        """恰好推进一轮 MP/SP；域阶段被割掉的顶点随后重新选择。"""
        record = self.session.step()
        record.update(global_iteration=len(self.history)+1, phase=self.phase)
        if record["cut_id"] is not None:
            self.polytopes = cut_polytopes(self.designs, self.polytopes, self.cuts[-1])
            self.snapshots.append(self.polytopes)
        record["cuts_after"] = len(self.cuts)
        self.history.append(record)
        if self.region_index is not None:
            if record["status"] == "optimal":
                self.checked.add(self.region_key)
            self.session.master.dispose()
            self._region_candidate()               # 当前顶点被切掉后不再固定该不可行点求解。
        elif self.session.finished:
            self._advance_phase()
        return record


def history_table(records):
    return pd.DataFrame([dict(
        iteration=r["global_iteration"], phase=r["phase"], status=r["status"],
        cost_cny=r.get("cost_cny", np.nan), p1=r.get("p", [np.nan]*3)[0],
        p2=r.get("p", [np.nan]*3)[1], p3=r.get("p", [np.nan]*3)[2],
        total_kw=r.get("total_kw", np.nan), eta=r.get("violation", np.nan),
        cut_id=r["cut_id"], cuts=r["cuts_after"], design=r.get("design", "")) for r in records])


def construction_table(model, x):
    rows = []
    for e, edge in enumerate(model.corridors):
        chosen = x[e*len(model.lines):(e+1)*len(model.lines)]
        k = int(np.argmax(chosen))
        rows.append(dict(corridor=edge.name, selected=bool(chosen.sum()>.5),
                         line=model.lines[k].name if chosen.sum()>.5 else "—",
                         cost_cny=float(model.cost[e*len(model.lines):(e+1)*len(model.lines)]@chosen)))
    return pd.DataFrame(rows)


def state_table(model, sp):
    names = [f"P_{e.name}" for e in model.corridors]+[f"v_{i}" for i in range(1,4)]+["eta"]
    return pd.DataFrame({"variable": names, "value": [*sp.w, sp.violation],
                         "unit": ["kW"]*len(model.corridors)+["p.u.^2"]*3+["dimensionless"]})


def dual_table(model, sp, active_only=False):
    """列出真实电气行、对偶和残差；η>0 的 SP 状态不标为可行潮流。"""
    names = [f"P_{e.name}" for e in model.corridors]+[f"v_{i}" for i in range(1,4)]
    left = [" ".join(f"{a:+.6g} {name}" for a,name in zip(row,names) if abs(a)>1e-12).lstrip("+") or "0"
            for row in model.W]
    table = pd.DataFrame({
        "row": [f"E{i+1:02d}" for i in range(len(model.h))],
        "constraint": [r.ConstrName for r in model.electrical_constraints],
        "SP inequality": [f"{lhs} <= {rhs:.6g} + {r:.6g} eta" for lhs,rhs,r in zip(left,sp.rhs,model.scales)],
        "pi=-Pi": sp.cut.pi, "rhs": sp.rhs, "Ww-rhs": sp.residual,
        "relaxation=r*eta": model.scales*sp.violation,
        "relaxed_slack": model.scales*sp.violation-sp.residual,
    })
    return table.loc[np.abs(sp.cut.pi)>1e-10] if active_only else table


def cut_table(process, x, count=None):
    """累计全部割；g 范围在所示固定 x、p≥0、Σp≤L 的初始单纯形上精确取得。"""
    rows = []
    for s, cut in enumerate(process.cuts[:count], 1):
        offset = cut.constant + cut.x_coeff@x
        rows.append(dict(cut=s, beta0=cut.constant, beta_p1=cut.load_coeff[0],
                         beta_p2=cut.load_coeff[1], beta_p3=cut.load_coeff[2],
                         fixed_x_offset=offset,
                         g_min=offset+process.limit*min(0., min(cut.load_coeff)),
                         g_max=offset+process.limit*max(0., max(cut.load_coeff)),
                         source_g=-cut.source_violation))
    return pd.DataFrame(rows)


def _stair(frontier, limit):
    xx, yy, previous = [], [], 0.
    for row in frontier.sort_values("cost_cny").itertuples():
        end = min(limit, row.total_kw)
        if end > previous+1e-8:
            xx.extend([previous, end]); yy.extend([row.cost_cny/1000]*2)
            previous = end
        if previous >= limit-1e-7:
            break
    return xx, yy


def plot_frontier(process, record=None):
    count = record["cuts_after"] if record else len(process.cuts)
    polygons = process.snapshots[count]
    fig, ax = plt.subplots(figsize=(7.1, 3.25), layout="constrained")
    ref_x, ref_y = _stair(process.reference, process.frontier_limit)
    out_x, out_y = _stair(outer_frontier(process.designs, polygons), process.frontier_limit)
    ax.plot(ref_x, ref_y, color="#555F65", ls="--", lw=1.5, label="Independent LP frontier")
    ax.plot(out_x, out_y, color=BLUE, lw=1.8, label=f"Outer bound: {count} cuts")
    if count:
        old_x, old_y = _stair(outer_frontier(process.designs, process.snapshots[count-1]), process.frontier_limit)
        ax.plot(old_x, old_y, color=ORANGE, lw=.9, alpha=.8, label="Before latest cut")
    ax.set(xlim=(0, process.frontier_limit), ylim=(-1, max(ref_y)*1.12),
           xlabel="Total load (kW)", ylabel="Budget (10³ CNY)")
    ax.legend(loc="upper left")
    ax.grid(axis="y", color="#E4E8EA", lw=.5)
    return fig


def union_surface(polytopes):
    """逐面扣除被其他多面体覆盖的部分，保留非凸并集的真实外表面。"""
    polys = maximal_polytopes(polytopes)
    hulls = [ConvexHull(p) for p in polys]
    surfaces = []
    for i, (points, hull) in enumerate(zip(polys, hulls)):
        planes = {}  # Qhull 三角面按相同支撑平面合并，避免显示人为对角线。
        for face, equation in zip(hull.simplices, hull.equations):
            key = tuple(np.round(equation, 8))
            planes.setdefault(key, (equation, set()))[1].update(face)
        for equation, indices in planes.values():
            normal, offset = equation[:3], equation[3] # 分组可舍入，几何判别必须保留原始平面精度。
            u = np.cross(normal, [1.,0.,0.] if abs(normal[0])<.9 else [0.,1.,0.])
            u /= np.linalg.norm(u)
            basis = np.array([u, np.cross(normal, u)])
            origin = points[next(iter(indices))]
            xy = (points[list(indices)]-origin)@basis.T
            visible = Polygon(xy[ConvexHull(xy).vertices])
            for j, (other, other_hull) in enumerate(zip(polys, hulls)):
                if i == j:
                    continue
                signed = other@normal + offset
                # 同向共面边界只保留编号较小的一份；跨到面外的体积遮住本面。
                if signed.max() < -1e-7 or signed.min() > 1e-7 or (signed.max()<=1e-7 and j>i):
                    continue
                section = list(other[np.abs(signed)<=1e-7])
                edges = {tuple(sorted(pair)) for face in other_hull.simplices
                         for pair in itertools.combinations(face,2)}
                for a,b in edges:
                    if signed[a]*signed[b] < 0:
                        section.append(other[a]+(other[b]-other[a])*signed[a]/(signed[a]-signed[b]))
                if len(section)<3:
                    continue  # 仅点或线相交，没有覆盖面面积。
                projected = np.unique(np.round((np.array(section)-origin)@basis.T,8),axis=0)
                if len(projected)<3 or np.linalg.matrix_rank(projected-projected[0],tol=1e-7)<2:
                    continue
                cover = Polygon(projected[ConvexHull(projected).vertices])
                visible = visible.difference(cover)
                if visible.is_empty:
                    break
            pieces = [visible] if visible.geom_type=='Polygon' else list(visible.geoms)
            for piece in pieces:
                if piece.area < 1e-8:
                    continue
                patches = list(shapely.constrained_delaunay_triangles(piece).geoms) if piece.interiors else [piece]
                for patch in patches:
                    surfaces.append(origin+np.asarray(patch.exterior.coords)[:-1]@basis)
    return surfaces


def _draw_polytopes(ax, surfaces, color, alpha, wire=False):
    if wire:
        edges = [p[[i,(i+1)%len(p)]] for p in surfaces for i in range(len(p))]
        artist = Line3DCollection(edges, colors=color, linewidth=.5, alpha=alpha)
    else:
        artist = Poly3DCollection(surfaces, facecolor=color, edgecolor=color, linewidth=.3, alpha=alpha)
    ax.add_collection3d(artist)


def region_frame(process, count):
    """同一个割池，按预算筛选 x 后取并集；缓存不改变任何几何成员。"""
    if count not in process.region_frames:
        panels = []
        for indices in process.budget_designs:
            polys = [process.snapshots[count][i] for i in indices]
            faces = union_surface(polys)
            boundary = []
            if count:
                cut = process.cuts[count-1]
                offsets = [cut.constant+cut.x_coeff@process.designs[i].x for i in indices]
                boundary = [f for f in faces if any(np.max(np.abs(f@cut.load_coeff+b))<2e-7 for b in offsets)]
            panels.append(dict(faces=faces, boundary=boundary,
                               capacity=max(float(p.sum(axis=1).max()) for p in polys)))
        process.region_frames[count] = panels
    return process.region_frames[count]


def plot_regions(process, record=None):
    """相同坐标范围的预算对照；交互版本使用同一份外表面。"""
    count = record["cuts_after"] if record else len(process.cuts)
    current, previous = region_frame(process,count), region_frame(process,max(0,count-1))
    fig = plt.figure(figsize=(7.1, 7.0))
    fig.subplots_adjust(left=.02, right=.91, bottom=.09, top=.96, wspace=.16, hspace=.23)
    for k, budget in enumerate(process.budgets):
        limit = process.limit
        ax = fig.add_subplot(math.ceil(len(process.budgets)/2), 2, k+1, projection="3d")
        _draw_polytopes(ax, previous[k]['faces'], GRAY, .3, wire=True)
        _draw_polytopes(ax, current[k]['faces'], BLUE, .32)
        if current[k]['boundary']:
            _draw_polytopes(ax, current[k]['boundary'], ORANGE, .6)
        if record and record.get("p") is not None and (budget is None or record["cost_cny"]<=budget):
            ax.scatter(*record["p"], color=ORANGE, s=16, depthshade=False)
        ax.set(xlim=(0,limit), ylim=(0,limit), zlim=(0,limit),
               xlabel=r"$p_1$ (kW)", ylabel=r"$p_2$ (kW)", zlabel=r"$p_3$ (kW)")
        ticks = [0,limit/2,limit]
        ax.set_xticks(ticks, labels=[f'{v:g}' for v in ticks])
        ax.set_yticks(ticks, labels=['',f'{limit/2:g}',f'{limit:g}']) # 共角处省去重复零标记，避免与p1末端重叠。
        ax.set_zticks(ticks, labels=[f'{v:g}' for v in ticks])
        ax.set_box_aspect((1,1,1), zoom=.87)
        ax.tick_params(pad=0, labelsize=7)
        ax.view_init(elev=23, azim=-55)
        label = "No budget limit" if budget is None else f"B = {budget:,.0f} CNY"
        ax.set_title(f"{chr(97+k)}   {label}", loc="left", fontsize=9)
        for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
            axis.pane.fill = False
            axis._axinfo["grid"]["color"] = (.85,.88,.90,.45)
    fig.legend(handles=[Patch(facecolor=BLUE, alpha=.35, label="Union after all cuts"),
                        Line2D([],[],color=GRAY,label="Before latest cut"),
                        Patch(facecolor=ORANGE,alpha=.6,label="Latest cut boundary")],
               loc="lower center", bbox_to_anchor=(.5,.005), ncol=3, fontsize=7)
    return fig


def _details(title, body, opened=False):
    return f'<details {"open" if opened else ""}><summary>{escape(title)}</summary>{body}</details>'


def round_html(process, record):
    if record["sp"] is None:
        return f'<p>第 {record["global_iteration"]} 轮：主问题证明不可行。</p>'
    sp = record["sp"]
    heading = (f'第 {record["global_iteration"]} 轮 · {record["phase"]} · '
               f'费用 {record["cost_cny"]:,.0f} 元 · 总负荷 {record["total_kw"]:.6f} kW · '
               f'η={sp.violation:.3g} · {record["status"]}')
    body = f'<p><b>{escape(heading)}</b></p>'
    body += construction_table(process.model, record["x"]).to_html(index=False)
    body += '<p>SP 运行状态；η&gt;0 时为松弛问题解。</p>'+state_table(process.model, sp).to_html(index=False)
    body += '<p>非零对偶乘子及对应电气约束：</p>'+dual_table(process.model, sp, True).to_html(index=False)
    body += _details('全部 53 条约束、对偶乘子及原始/松弛残差', dual_table(process.model, sp).to_html(index=False))
    return body


def show_process(process, record):
    """刷新画面只替换显示；历史与全部累计割均保存在 process 中。"""
    display(HTML(round_html(process, record)))
    if record.get("x") is not None:
        display(HTML(f'<p>全部累计割：β0+βx·x+βp·p≥0；下表 g 范围固定为本轮 x、p≥0、Σp≤{process.limit:g}。</p>'
                     +cut_table(process, record["x"]).to_html(index=False)))
        coefficients = pd.DataFrame([c.x_coeff for c in process.cuts],
                                    columns=[f'x_{e.name}_{k.name}' for e in process.model.corridors for k in process.model.lines])
        coefficients.index = np.arange(1,len(coefficients)+1)
        display(HTML(_details('全部累计割的 15 个建设系数 βx（行号为割编号）',coefficients.to_html())))
    for plot in (plot_frontier, plot_regions):
        fig = plot(process, record)
        display(fig)
        plt.close(fig)
    display(HTML(f'<p>当前阶段：{escape(process.phase)}；累计 {len(process.history)} 轮、{len(process.cuts)} 条割。'
                 '全部历史在本批计算结束后的交互回放中查看；重跑本块继续，重跑代码块1重新开始。</p>'))


def save_records(process, output):
    """保存原始逐轮证据；JSON 含完整 x、p、w、π、约束残差及所有割系数。"""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    history_table(process.history).to_csv(output/'history.csv', index=False, encoding='utf-8-sig')
    cuts = [dict(cut=i+1, beta0=c.constant, **{f'beta_x{j}':v for j,v in enumerate(c.x_coeff)},
                 **{f'beta_p{j+1}':v for j,v in enumerate(c.load_coeff)}) for i,c in enumerate(process.cuts)]
    pd.DataFrame(cuts).to_csv(output/'cuts.csv', index=False)
    process.reference.to_csv(output/'reference_frontier.csv', index=False)
    pd.DataFrame(process.frontier).to_csv(output/'discovered_frontier.csv', index=False)
    def encode(value):
        return value.tolist() if isinstance(value, np.ndarray) else vars(value)
    data = dict(budgets=process.budgets, power_limit=process.limit, frontier_limit=process.frontier_limit, target=process.target,
                finished=process.finished, cuts=process.cuts, history=process.history,
                constraints=[r.ConstrName for r in process.model.electrical_constraints],
                designs=process.designs, final_polytopes=process.polytopes)
    (output/'process.json').write_text(json.dumps(data, default=encode, ensure_ascii=False, indent=1), encoding='utf-8')


def export_playback(process, output):
    """导出全部原始轮次及可旋转的离线回放；无新增割轮次复用同一份网格。"""
    from node_power_replay import replay_data, replay_html
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    data = replay_data(process)
    (output/'replay.json').write_text(json.dumps(data,ensure_ascii=False,separators=(',',':')),encoding='utf-8')
    (output/'playback.html').write_text(replay_html([dict(name='完整实验',data=data)]),encoding='utf-8')
    for name, plot in [('frontier',plot_frontier), ('regions',plot_regions)]:
        fig = plot(process)
        fig.savefig(output/f'{name}_final.svg')
        fig.savefig(output/f'{name}_final.pdf')
        fig.savefig(output/f'{name}_final.png',dpi=300)
        plt.close(fig)
    print(f'已导出 {len(process.history)} 轮交互回放、{len(process.cuts)} 条割及最终矢量图。',flush=True)
