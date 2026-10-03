"""绘图基元：回放窗口各面板的配色与步骤样式，二维几何、三维凸域、锥远端面片、体素表面与联合割截线的绘制，
以及收敛过程图与 Case33-S 勘察的四张图。只画图，不做构域计算；几何并集等计算在 region.py，勘察的数值在 survey.py。"""
from functools import lru_cache
from itertools import product

import numpy as np

INNER, OUTER, GLOBAL, SP = '#397f85', '#a9b9c4', '#cc3838', '#ef8a23'
REFERENCE = '#c6dbef'   # 对比面板中参考可行的格
CUT, RAY, NETWORK = '#8055a4', '#b8860b', '#555555'
STEP_STYLE = {   # 步骤几何：标记、颜色、主图短标签（完整说明在步骤栏）
    'center': ('D', GLOBAL, 'max Σξ'), 'origin': ('o', NETWORK, 'SP  p=0'), 'ray': ('o', RAY, 'max t'),
    'near': ('s', RAY, '反向 max t'), 'cone': ('D', GLOBAL, 'max c·ξ → μ'), 'lazy': ('P', GLOBAL, '现任网架加入紧化行'),
    'sp': ('x', SP, 'min η'), 'cut': ('X', CUT, '取割')}


def cut_segment(cut, x, bounds, axis_lower=None):
    """alpha + beta @ p + delta @ x = 0 在二维显示框中的截线。"""
    cut, x, bounds = np.asarray(cut), np.asarray(x), np.asarray(bounds)
    lower = np.zeros(2) if axis_lower is None else np.asarray(axis_lower)
    beta, constant = cut[1:3], cut[0]+cut[3:]@x
    points = []
    for fixed in (0, 1):
        free = 1-fixed
        if beta[free] == 0.:
            continue
        for edge in (lower[fixed], bounds[fixed]):
            value = -(constant+beta[fixed]*edge)/beta[free]
            tolerance = 1e-9*(bounds[free]-lower[free])
            if lower[free]-tolerance <= value <= bounds[free]+tolerance:
                point = np.zeros(2)
                point[fixed], point[free] = edge, np.clip(value, lower[free], bounds[free])
                if not any(np.allclose(point, old, rtol=1e-9, atol=1e-9) for old in points):
                    points.append(point)
    return np.asarray(points).reshape(-1, 2)


def cut_polygon(cut, x, bounds, axis_lower=None):
    """三维联合割平面与显示盒十二条棱的交点，保留实际 kW 坐标。"""
    cut, bounds = np.asarray(cut), np.asarray(bounds)
    lower = np.zeros(3) if axis_lower is None else np.asarray(axis_lower)
    beta, constant = cut[1:4], cut[0]+cut[4:]@x
    points = []
    for free in range(3):
        if beta[free] == 0.:
            continue
        fixed = [i for i in range(3) if i != free]
        for corner in product((0., 1.), repeat=2):
            point = np.zeros(3)
            point[fixed] = lower[fixed]+np.asarray(corner)*(bounds-lower)[fixed]
            point[free] = -(constant+beta@point)/beta[free]
            tolerance = 1e-9*(bounds[free]-lower[free])
            if lower[free]-tolerance <= point[free] <= bounds[free]+tolerance:
                points.append(point)
    return np.unique(np.asarray(points).reshape(-1, 3), axis=0)


@lru_cache(maxsize=8192)
def hull_geometry(data, scale):
    """一个凸域的显示几何 (面, 棱, 点)，按顶点字节缓存：回放各帧重复出现的同一凸域只算一次凸包。"""
    from scipy.spatial import ConvexHull
    points = np.frombuffer(data).reshape(-1, 3)
    # 1. 仿射维数：点、线段、平面多边形或三维凸体
    delta = (points-points[0])/np.asarray(scale)
    rank = np.linalg.matrix_rank(delta, tol=1e-10)
    coordinates = delta@np.linalg.svd(delta, full_matrices=False)[2][:rank].T
    if rank == 0:
        return (), (), (points[0],)
    if rank == 1:
        return (), (points[[coordinates[:, 0].argmin(), coordinates[:, 0].argmax()]],), ()
    # 2. 凸包只用于显示：QJ 对近退化的细锥域也给出凸包，扰动远小于识别真实棱的面方程阈值
    hull = ConvexHull(coordinates/np.linalg.norm(coordinates, axis=0), qhull_options='QJ')
    if rank == 2:
        face = points[hull.vertices]
        return (face,), tuple(np.stack([face, np.roll(face, -1, axis=0)], axis=1)), ()
    # 3. 三维凸体：三角面全部填充，只画相邻面方程不同的真实棱
    edges = tuple(points[np.intersect1d(hull.simplices[i], hull.simplices[j])]
                  for i, neighbors in enumerate(hull.neighbors) for j in neighbors
                  if j > i and np.linalg.norm(hull.equations[i]-hull.equations[j]) > 1e-7)
    return tuple(points[hull.simplices]), edges, ()


def draw_3d(ax, polytopes, bounds, *, color, fill=False, alpha=1., linestyle='-', linewidth=.8, gid=None):
    """凸域的面与真实棱：全部凸域合成一个面集合与一个棱集合一次添加（三维每添加一个集合都会重算坐标范围）。"""
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
    faces, edges, dots = [], [], []
    scale = tuple(float(b) for b in np.asarray(bounds, float))
    for points in polytopes:
        points = np.asarray(points, float).reshape(-1, 3)
        if len(points):
            geometry = hull_geometry(points.tobytes(), scale)
            faces.extend(geometry[0])
            edges.extend(geometry[1])
            dots.extend(geometry[2])
    if fill and faces:
        ax.add_collection3d(Poly3DCollection(faces, facecolors=color, edgecolors='none', alpha=alpha, gid=gid))
    if edges:
        ax.add_collection3d(Line3DCollection(edges, colors=color, linewidths=linewidth, linestyles=linestyle,
                                             alpha=min(1., 3*alpha) if fill else alpha, gid=gid))
    if dots:
        ax.scatter(*np.transpose(dots), color=color, s=14, gid=gid)


def voxel_faces(states, bounds):
    """扫描可行体素的六邻接外表面，不插值、不构造跨体素凸包。"""
    truth = np.asarray(states) == 1
    padded = np.pad(truth, 1)
    step = np.asarray(bounds)/np.asarray(truth.shape)
    faces = []
    for axis in range(3):
        fixed = [i for i in range(3) if i != axis]
        for side in (-1, 1):
            neighbors = np.roll(padded, -side, axis=axis)[1:-1, 1:-1, 1:-1]
            starts = np.argwhere(truth & ~neighbors)
            corners = np.zeros((4, 3))
            corners[:, axis] = int(side == 1)
            corners[:, fixed] = [(0, 0), (1, 0), (1, 1), (0, 1)]
            faces.append((starts[:, None, :]+corners)*step)
    return np.concatenate(faces)


def cap(points):
    """锥体远端的面片：去掉原点，其余顶点按绕形心的角度排序（三维一次性绘制）。"""
    return _cap_geometry(np.asarray(points, float).tobytes())


@lru_cache(maxsize=8192)
def _cap_geometry(data):
    points = np.frombuffer(data).reshape(-1, 3)
    points = points[np.linalg.norm(points, axis=1) > 1e-9]
    center = points.mean(axis=0)
    basis = np.linalg.svd(np.eye(3)-np.outer(center, center)/(center@center))[0][:, :2]
    return points[np.argsort(np.arctan2(*((points-center)@basis).T[::-1]))]


def draw_geometry(ax, geometry, *, color, fill=False, alpha=1., linestyle='-', linewidth=1.):
    """二维 shapely 几何（多边形并集、线、点）：保留不相连部分和孔洞。"""
    from matplotlib.path import Path as MplPath
    from matplotlib.patches import PathPatch
    if geometry.is_empty:
        return
    if geometry.geom_type in ('MultiPolygon', 'GeometryCollection', 'MultiLineString', 'MultiPoint'):
        for part in geometry.geoms:
            draw_geometry(ax, part, color=color, fill=fill, alpha=alpha, linestyle=linestyle, linewidth=linewidth)
    elif geometry.geom_type == 'Polygon':
        from shapely.geometry.polygon import orient
        poly = orient(geometry, sign=1.)
        paths = []
        for ring in [poly.exterior, *poly.interiors]:
            points = np.asarray(ring.coords)
            codes = [MplPath.MOVETO]+[MplPath.LINETO]*(len(points)-2)+[MplPath.CLOSEPOLY]
            paths.append(MplPath(points, codes))
        ax.add_patch(PathPatch(MplPath.make_compound_path(*paths), facecolor=color if fill else 'none', edgecolor=color,
                               alpha=alpha, linewidth=linewidth, linestyle=linestyle))
    else:
        points = np.asarray(geometry.coords)
        ax.plot(points[:, 0], points[:, 1], color=color, linewidth=linewidth,
                marker='.' if len(points) == 1 else None, linestyle=linestyle, alpha=alpha)


SERIES = ('#2a78d6', '#eb6834', '#1baf7a', '#4a3aa7')   # 多次运行的配色（分类色前四位，两两通过色觉检验）


def draw_convergence(runs, epsilon, path):
    """收敛过程图。上行：全部分区内域相对 AC 的 MR（对数）与 FR；下面每个分区一格夹逼间隙（对数，虚线为 ε，
    圆点 / 叉为分区获证 / 未获证结束）。runs 为 {图例: vertify.convergence 的结果}，同色为同一次运行。"""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.sans-serif': ['Microsoft YaHei', 'DejaVu Sans'], 'axes.unicode_minus': False,
                         'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.edgecolor': '#c3c2b7', 'axes.labelcolor': '#52514e', 'xtick.color': '#52514e',
                         'ytick.color': '#52514e', 'grid.color': '#e1e0d9', 'grid.linewidth': .6})
    # 1. 版面：上行 MR、FR 各占两列，下面每行四个分区
    labels = list(next(iter(runs.values()))['gaps'])
    rows = 1+(len(labels)+3)//4
    figure = plt.figure(figsize=(14, 3.3*rows+.8))
    grid = figure.add_gridspec(rows, 4)
    mr_axis, fr_axis = figure.add_subplot(grid[0, :2]), figure.add_subplot(grid[0, 2:])
    gap_axes = [figure.add_subplot(grid[1+k//4, k % 4]) for k in range(len(labels))]
    # 2. 各次运行：MR / FR 阶梯线，各分区间隙阶梯线与结束标记
    for color, (name, run) in zip(SERIES, runs.items()):
        mr_axis.step(run['time'], run['mr'], where='post', color=color, lw=1.6, label=name)
        fr_axis.step(run['time'], run['fr'], where='post', color=color, lw=1.6)
        for axis, label in zip(gap_axes, labels):
            time, gap = np.array(run['gaps'][label]).T
            axis.step(time, 100*gap, where='post', color=color, lw=1.6)
            ended, certified = run['ends'][label]
            axis.plot(ended, 100*gap[-1], 'o' if certified else 'x', color=color, ms=5)
    # 3. 坐标与说明：全部子图共用时间轴
    end = 1.03*max(run['time'][-1] for run in runs.values())
    mr_axis.set(yscale='log', xlim=(0., end), title='内域 I 相对 AC 的遗漏率 MR（对数）', xlabel='时间 (s)', ylabel='MR (%)')
    fr_axis.set(xlim=(0., end), title='内域 I 相对 AC 的多余率 FR', xlabel='时间 (s)', ylabel='FR (%)')
    for axis, label in zip(gap_axes, labels):
        axis.axhline(100*epsilon, color='#898781', ls='--', lw=1)
        axis.set(yscale='log', ylim=(.1, 1e3), xlim=(0., end), title=f'分区 {label}', xlabel='时间 (s)')
    gap_axes[0].set_ylabel('夹逼间隙 vol(K^OUT)/vol(I)−1 (%)')
    for axis in (mr_axis, fr_axis, *gap_axes):
        axis.grid(True)
    figure.legend(*mr_axis.get_legend_handles_labels(), loc='lower center', ncol=len(runs), frameon=False, fontsize=8,
                  title=f'间隙虚线 ε={100*epsilon:g}%；圆点为分区获证、叉为未获证的结束时刻', title_fontsize=8)
    figure.tight_layout(rect=(0, .7/(3.3*rows+.8), 1, 1))
    figure.savefig(path, dpi=130)
    plt.close(figure)


# ---- 勘察图（Case33-S；图中文字为英文） -----------------------------------------------------------
CATEGORY = ('#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#4a3aa7', '#e34948')   # 分类色，固定次序
BLUES = ('#cde2fb', '#b7d3f6', '#9ec5f4', '#86b6ef', '#6da7ec', '#5598e7', '#3987e5', '#2a78d6', '#256abf',
         '#1c5cab', '#184f95', '#104281', '#0d366b')                                            # 顺序色，浅 → 深
INK, MUTED, BASELINE, NEUTRAL = '#52514e', '#898781', '#c3c2b7', '#f0efec'
CASE33_LAYOUT = {**{i: (i-1., 0.) for i in range(1, 19)}, **{i: (i-17., -1.2) for i in range(19, 23)},
                 **{i: (i-21., 2.4) for i in range(23, 26)}, **{i: (i-20., 1.2) for i in range(26, 34)}}
# Case33 节点坐标：主馈线 1–18 水平，19–22 在下方，23–25 与 26–33 在上方


def _english():
    """勘察图的样式：英文无衬线、浅色坐标线；返回 pyplot。"""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.sans-serif': ['DejaVu Sans'], 'font.size': 8, 'axes.spines.top': False,
                         'axes.spines.right': False, 'axes.edgecolor': BASELINE, 'axes.labelcolor': INK,
                         'xtick.color': INK, 'ytick.color': INK, 'grid.color': '#e1e0d9', 'grid.linewidth': .6,
                         'pdf.fonttype': 42})
    return plt


def _save(plt, figure, path):
    """同名的 pdf 与 png，裁去空白。"""
    figure.savefig(path.with_suffix('.pdf'), bbox_inches='tight', pad_inches=.05)
    figure.savefig(path.with_suffix('.png'), dpi=200, bbox_inches='tight', pad_inches=.05)
    plt.close(figure)


def draw_survey_layout(network, path):
    """布置图：固定线路实线，可开断线路 S1–S5 加方块，候选路 C1–C5 虚线并标相对建设费，端口画圆环；主馈线上的候选画成弧。"""
    plt = _english()
    from matplotlib.lines import Line2D
    from matplotlib.path import Path as MplPath
    from matplotlib.patches import PathPatch
    figure, ax = plt.subplots(figsize=(7.2, 3.0))
    switches = {frozenset(edge): f'S{k+1}' for k, edge in enumerate(network.switches)}
    roads = {frozenset(edge): f'C{k+1} ({cost:g})' for k, (edge, cost) in enumerate(network.candidates.items())}
    for corridor in network.corridors:
        (x0, y0), (x1, y1) = (np.array(CASE33_LAYOUT[i]) for i in corridor.endpoints)
        key = frozenset(corridor.endpoints)
        if key in roads:
            # 候选路：虚线；两端都在主馈线上时画成向上的弧；标签衬白底放在弧顶或中点
            middle = ((x0+x1)/2, .7 if y0 == y1 == 0. else (y0+y1)/2)
            if y0 == y1 == 0.:
                ax.add_patch(PathPatch(MplPath([(x0, y0), (middle[0], 1.4), (x1, y1)], [1, 3, 3]), fill=False,
                                       edgecolor=CATEGORY[1], linestyle='--', linewidth=1.4))
            else:
                ax.plot([x0, x1], [y0, y1], color=CATEGORY[1], ls='--', lw=1.4)
            ax.text(*middle, roads[key], color=CATEGORY[1], ha='center', va='center', fontsize=7, zorder=6,
                    bbox=dict(facecolor='white', edgecolor='none', pad=.6))
        else:
            ax.plot([x0, x1], [y0, y1], color=INK, lw=1.2, zorder=2)
            if key in switches:
                ax.plot((x0+x1)/2, (y0+y1)/2, marker='s', ms=6, mfc='white', mec=INK, mew=1.1, zorder=3)
                ax.text((x0+x1)/2, (y0+y1)/2+.16, switches[key], color=INK, ha='center', va='bottom', fontsize=7)
    for node, (x, y) in CASE33_LAYOUT.items():
        ax.plot(x, y, 'o', ms=3.2, color=INK, zorder=4)
        ax.text(x, y+(.17 if y > 0. else -.17), str(node), ha='center', va='bottom' if y > 0. else 'top',
                fontsize=5.5, color=MUTED)
    for node in network.load_nodes:
        ax.plot(*CASE33_LAYOUT[node], 'o', ms=11, mfc='none', mec=CATEGORY[0], mew=1.6, zorder=5)
    ax.plot(*CASE33_LAYOUT[1], marker='s', ms=8, color=INK, zorder=5)
    ax.legend(handles=[Line2D([], [], color=INK, lw=1.2, label='Fixed line'),
                       Line2D([], [], color=INK, lw=1.2, marker='s', mfc='white', mec=INK,
                              label='Switchable line S1–S5 (cost 0)'),
                       Line2D([], [], color=CATEGORY[1], lw=1.4, ls='--', label='Candidate road C1–C5 (relative cost)'),
                       Line2D([], [], ls='', marker='o', ms=9, mfc='none', mec=CATEGORY[0], mew=1.6,
                              label='Port: z1 = p18, z2 = p25'),
                       Line2D([], [], ls='', marker='s', ms=7, color=INK, label='Substation (bus 1)')],
              loc='upper center', bbox_to_anchor=(.5, 0.), ncol=3, frameon=False, fontsize=7)
    ax.set_aspect('equal')
    ax.axis('off')
    figure.tight_layout()
    _save(plt, figure, path)


def draw_survey_valuation(levels, curves, roads, interaction, path):
    """估值图：(a) κ(A,b)–b 阶梯前沿（curves 依次为 ∅、各单条候选、全部候选）；(b) 两两交互 I_ij 热图，正为互补、负为替代。"""
    plt = _english()
    from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
    figure, (left, right) = plt.subplots(1, 2, figsize=(7.2, 3.0), gridspec_kw=dict(width_ratios=(1.4, 1.)))
    # 1. 阶梯前沿：∅ 与全部候选用中性色，单条候选用分类色；∅ 的虚线画在最上层，与之重合的曲线仍可辨
    x = np.r_[levels, levels[-1]+1.]
    colors = [MUTED, *CATEGORY[:len(curves)-2], INK]
    for (name, kappa), color in zip(curves.items(), colors):
        left.step(x, np.r_[kappa, kappa[-1]], where='post', color=color, lw=1.4, ls='--' if color == MUTED else '-',
                  zorder=3 if color == MUTED else 2, label=name)
    left.set(xlabel='Budget b (relative build cost)', ylabel='Coverage κ(A, b)', xticks=levels, xlim=(0., x[-1]),
             ylim=(0., 1.05*max(np.max(kappa) for kappa in curves.values())), title='(a) Coverage frontier')
    left.grid(True)
    left.legend(frameon=False, fontsize=6.5, loc='lower right', ncol=2)
    # 2. 交互热图：发散色（红为替代、蓝为互补，零为中性灰），对角线不定义
    top = np.nanmax(np.abs(interaction))
    colormap = LinearSegmentedColormap.from_list('interaction', [CATEGORY[7], NEUTRAL, CATEGORY[0]])
    colormap.set_bad('white')
    image = right.imshow(interaction, cmap=colormap, norm=TwoSlopeNorm(0., -top, top))
    for (i, j), value in np.ndenumerate(interaction):
        right.text(j, i, '—' if np.isnan(value) else f'{value:+.2f}', ha='center', va='center', fontsize=7,
                   color='white' if abs(value) > .6*top else INK)
    right.set(xticks=range(len(roads)), yticks=range(len(roads)), xticklabels=roads, yticklabels=roads,
              title='(b) Pairwise interaction I_ij\n(+ complementary, − substitutive)')
    right.spines[:].set_visible(False)
    figure.colorbar(image, ax=right, fraction=.046, pad=.04, label='I_ij = V({i,j}) − V({i}) − V({j})')
    figure.tight_layout()
    _save(plt, figure, path)


def draw_survey_storyboard(panels, top, levels, status_quo, history, path):
    """故事板。上排每个状态一幅 z1–z2 图：颜色为确认集下的 C*(z)（按预算档离散，白色为不可服务），虚线为乐观域
    R(A∪U,14) 的边界，星号为 z^0；下排为 κ(A_s,14)、κ(A_s∪U_s,14) 与 Φ(A_s)、Φ(A_s∪U_s) 随累计勘察费的阶梯线。"""
    plt = _english()
    from matplotlib.colors import BoundaryNorm, ListedColormap
    n = len(panels)
    figure = plt.figure(figsize=(max(7.2, 2.1*n+1.2), 6.))
    grid = figure.add_gridspec(2, 1, height_ratios=(1.3, 1.), hspace=.55)
    row = grid[0].subgridspec(1, n+1, width_ratios=(*[1.]*n, .08), wspace=.25)
    bottom = grid[1].subgridspec(1, 2, wspace=.3)
    # 1. 各状态的 C*(z) 栅格（按预算档序号着色，每档一格色标）、乐观域边界与现状点
    colormap = ListedColormap(BLUES[1:1+len(levels)])
    colormap.set_bad('white')
    norm = BoundaryNorm(np.arange(len(levels)+1)-.5, colormap.N)
    for k, panel in enumerate(panels):
        ax = figure.add_subplot(row[0, k])
        index = np.where(np.isnan(panel['level']), np.nan, np.searchsorted(levels, panel['level']))
        image = ax.imshow(index, origin='lower', extent=(0., top[0], 0., top[1]), cmap=colormap, norm=norm,
                          aspect='auto', interpolation='nearest')
        draw_geometry(ax, panel['optimistic'], color=INK, linestyle='--', linewidth=1.)
        ax.plot(*status_quo, marker='*', ms=9, color=CATEGORY[1], mec='white', mew=.5)
        ax.set(xlim=(0., top[0]), ylim=(0., top[1]), xlabel='z1 = p18 (kW)', ylabel='z2 = p25 (kW)' if k == 0 else None)
        ax.set_title(f"s{k}: {panel['title']}\nκ conf / opt = {panel['kappa'][0]:.3f} / {panel['kappa'][1]:.3f}\n"
                     f"Φ ∈ [{panel['phi'][0]:.2f}, {panel['phi'][1]:.2f}]", fontsize=7)
    bar = figure.colorbar(image, cax=figure.add_subplot(row[0, n]), ticks=range(len(levels)),
                          label='C*(z): least build cost')
    bar.set_ticklabels([f'{b:g}' for b in levels])
    # 2. 两个 κ 与 Φ 上下界随累计勘察费（到达各状态时）的阶梯线；末段延长以便看清
    spent = np.asarray(history['spent'])
    x = np.r_[spent, spent[-1]+max(.1*spent[-1], .05)]
    for ax, (confirmed, optimistic), ylabel, title in (
            (figure.add_subplot(bottom[0, 0]), ('kappa', 'kappa_optimistic'), 'Coverage κ(·, 14)',
             '(a) Confirmed vs optimistic coverage'),
            (figure.add_subplot(bottom[0, 1]), ('upper', 'lower'), 'Expected cost Φ',
             '(b) Bounds Φ(A_s) ≥ J(s) ≥ Φ(A_s ∪ U_s)')):
        for key, color, style, name in ((confirmed, CATEGORY[0], '-', 'confirmed A_s'),
                                        (optimistic, CATEGORY[1], '--', 'optimistic A_s ∪ U_s')):
            values = np.asarray(history[key])
            ax.step(x, np.r_[values, values[-1]], where='post', color=color, ls=style, lw=1.5, label=name)
            ax.plot(spent, values, 'o', color=color, ms=3.5)
        ax.set(xlabel='Cumulative survey cost', ylabel=ylabel, title=title, xlim=(0., x[-1]))
        ax.grid(True)
        ax.legend(frameon=False, fontsize=7)
    _save(plt, figure, path)


def draw_survey_policy(costs, optimum, gaps, q0, rho, path):
    """策略图：(a) 主设置下五种做法的 E[C^tot]（点线为 J(s0)）；(b) 敏感性网格上本文与单路比值的 gap 热图（%，格内标数值）。"""
    plt = _english()
    from matplotlib.colors import LinearSegmentedColormap
    figure = plt.figure(figsize=(8.4, 3.))
    grid = figure.add_gridspec(1, 4, width_ratios=(1.7, 1., 1., .06), wspace=.4)
    # 1. 柱状图：不勘察、单路比值、本文、DP、全知
    bars = figure.add_subplot(grid[0, 0])
    bars.bar(range(len(costs)), list(costs.values()), width=.62,
             color=[MUTED, CATEGORY[1], CATEGORY[0], CATEGORY[2], BASELINE])
    for i, value in enumerate(costs.values()):
        bars.text(i, value, f'{value:.3f}', ha='center', va='bottom', fontsize=6.5, color=INK)
    bars.axhline(optimum, color=INK, ls=':', lw=.8)
    bars.set_xticks(range(len(costs)), ['\n'.join(name.split(' ', 1)) for name in costs], fontsize=6.5)
    bars.set(ylabel='E[C_tot]', title='(a) Expected total cost (dotted: J(s0))')
    # 2. 敏感性热图：行 q0、列 ρ，两幅共用色标
    colormap = LinearSegmentedColormap.from_list('gap', BLUES)
    top = 100*np.max(gaps)
    for k, (values, title) in enumerate(zip(gaps, ('(b) Bundle index (ours)\ngap to DP (%)',
                                                   '(c) Single-road ratio\ngap to DP (%)'))):
        ax = figure.add_subplot(grid[0, 1+k])
        image = ax.imshow(100*values, cmap=colormap, vmin=0., vmax=top, origin='lower', aspect='auto')
        for (i, j), value in np.ndenumerate(100*values):
            ax.text(j, i, f'{value:.1f}', ha='center', va='center', fontsize=7, color='white' if value > .6*top else INK)
        ax.set(xticks=range(len(rho)), yticks=range(len(q0)), xticklabels=[f'{v:g}' for v in rho],
               yticklabels=[f'{v:g}' for v in q0], xlabel='ρ (survey cost ratio)',
               ylabel='q0 (prior mean)' if k == 0 else None, title=title)
        ax.spines[:].set_visible(False)
    figure.colorbar(image, cax=figure.add_subplot(grid[0, 3]), label='gap (%)')
    _save(plt, figure, path)
