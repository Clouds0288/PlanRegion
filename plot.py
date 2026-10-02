"""绘图基元：回放窗口各面板的配色与步骤样式，二维几何、三维凸域、锥远端面片、体素表面与联合割截线的绘制。
只画图，不做构域计算；几何并集等计算在 region.py。"""
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
