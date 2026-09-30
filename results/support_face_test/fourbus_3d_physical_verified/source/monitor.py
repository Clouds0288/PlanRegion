"""原生监视器：RunMonitor 记录事件，NativeWindow 重建帧并绘图。

录制状态、显示范围和过程图层各自独立；扫描结果不进入构域时间轴。
二维/三维共用外包络、过程标记及回放控制，只在几何绘制处区分维数。
"""
from pathlib import Path
from threading import Condition, Event, Thread
from time import perf_counter
import gzip
import json

import numpy as np


class RegionTimeout(RuntimeError):
    """本次构域未在时限内取得全局证书。"""


def _plain(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [_plain(v) for v in value]
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def _merge(state, patch):
    """每个网架独立更新，历史帧只存发生变化的网架。"""
    for key, value in patch.items():
        if key in ('schemes', 'cut_history'):
            state[key] = {**state.get(key, {}), **value}
        else:
            state[key] = value


def signed_values(value, sign, prefix, key=''):
    """记录时把幅值、割系数、网架标签转成带符号坐标。"""
    if value is None:
        return None
    if key in ('scheme', 'active_scheme'):
        return prefix+value
    if key in ('p', 'anchor', 'target', 'vertices', 'ray_target'):
        return (np.asarray(value)*sign).tolist()
    if key == 'cut':
        cut = np.asarray(value).copy()
        cut[1:1+len(sign)] *= sign
        return cut.tolist()
    if isinstance(value, dict):
        return {k: signed_values(v, sign, prefix, k) for k, v in value.items()}
    if isinstance(value, list):
        return [signed_values(v, sign, prefix) for v in value]
    return value


class RunMonitor:
    """求解线程只提交数值；Tk 与绘图仅在主线程执行。"""

    def __init__(self, *, output=None, clock=perf_counter, algorithm='主线', parent=None, sign=None):
        self.output = None if output is None else Path(output)
        self._clock = clock
        self.algorithm = algorithm
        self.parent, self.sign = parent, None if sign is None else np.asarray(sign)
        self.count_offsets = {key: parent.state.get(key, 0) for key in ('sp', 'cuts', 'global_search')} if parent else {}
        self.condition = Condition()
        self.cancelled = Event()
        self.paused = False
        self.permits = 0
        self.paused_seconds = 0.
        self.started = self._clock()
        self.time_limit = np.inf
        self.state, self.history, self.scheme_ids = {}, [], {}
        self.validation_state = {}
        self.sp_numbers = {}
        self.result = self.error = None
        self.busy = False

    def clock(self):
        return self._clock()-self.paused_seconds

    def remaining(self, limit):
        remaining = self.time_limit-(self.clock()-self.started)
        if remaining <= 0.:
            raise RegionTimeout(f'构域超过 {self.time_limit:g} 秒，未取得全局覆盖证书')
        return min(limit, remaining)

    def timing(self):
        return dict(total_seconds=self.clock()-self.started)

    def _emit(self, event, *, checkpoint=True, **values):
        with self.condition:
            if self.cancelled.is_set():
                raise KeyboardInterrupt()
            values = _plain(values)
            patch = {k: v for k, v in values.items() if self.state.get(k) != v}
            for key in ('schemes', 'cut_history'):
                if key in patch:
                    patch[key] = {k: v for k, v in patch[key].items()
                                  if self.state.get(key, {}).get(k) != v}
                    if not patch[key]:
                        del patch[key]
            patch['event'] = event
            _merge(self.state, patch)
            if self.parent is None:
                self.history.append(dict(elapsed=self.clock()-self.started, patch=patch))
            else:
                self.parent.partition_frame(self, patch)
            before = self._clock()
            while checkpoint and self.paused and not self.permits and not self.cancelled.is_set():
                self.condition.wait()
            self.paused_seconds += self._clock()-before
            if self.permits:
                self.permits -= 1
            if self.cancelled.is_set():
                raise KeyboardInterrupt()

    def partition_frame(self, local, patch):
        """局部分区只保留当前状态；真实坐标下的增量帧统一交给现有时间轴。"""
        from region import initial_polytope
        sign = local.sign
        label = ''.join('+' if s > 0 else '-' for s in sign)
        prefix, d = label+':', len(sign)
        mapped = signed_values(patch, sign, prefix)
        for key in ('bounds', 'axis_bounds', 'total_bound', 'result', 'time_limit'):
            mapped.pop(key, None)
        if patch['event'] == 'phase_start':
            mapped.update(active_scheme=None, ray=None, seed_point=None, coverage_complete=False)
        for key in local.count_offsets:
            if key in mapped:
                mapped[key] += local.count_offsets[key]
        if 'schemes' in patch:
            mapped['schemes'] = {prefix+key: {**row, 'sign': sign,
                'inner': np.asarray(row['inner']).reshape(-1, d)*sign,
                'outer': np.asarray(row['outer']).reshape(-1, d)*sign} for key, row in patch['schemes'].items()}
        if 'cut_history' in mapped:
            mapped['cut_history'] = {str(int(key)+local.count_offsets['cuts']): {**row, 'sign': sign}
                                     for key, row in mapped['cut_history'].items()}
        state = local.state
        if 'result' in patch:
            self.envelopes[label] = [np.asarray(row['vertices'])*sign for row in patch['result']['outer']]
        else:
            self.envelopes[label] = [initial_polytope(state['bounds'], state['total_bound'], state['axis_bounds'])*state['bounds']*sign]
        outer = [p for polys in self.envelopes.values() for p in polys]
        vertices = np.vstack(outer)
        mapped.update(partition=label, global_outer=outer,
                      axis_lower=vertices.min(axis=0) if self.state['mode'] else np.zeros(d),
                      axis_bounds=vertices.max(axis=0))
        event = mapped.pop('event')
        self._emit('partition_end' if event == 'region_end' else event, **mapped)

    def _scheme(self, x):
        key = tuple(x)
        if key not in self.scheme_ids:
            index = len(self.scheme_ids)
            self.scheme_ids[key] = chr(65+index%26)+(str(index//26+1) if index >= 26 else '')
        return self.scheme_ids[key]

    def _geometry(self, region):
        return dict(axis_bounds=region.axis_bounds, total_bound=region.total_bound,
                    cuts=len(region.cuts), schemes={self._scheme(row['x']): dict(
                        x=row['x'], choice=row['choice'], cost=row['cost'], inner=row['inner']*region.bounds,
                        outer=row['outer']*region.bounds) for row in region.records.values()})

    # 1. 初始化：每个完整 MP2 的可行点立即入帧。
    def begin(self, network, method, budget, region, time_limit):
        self.started = self.clock()
        self.time_limit = time_limit
        self._emit('phase_start', phase='初始化', network=network.name, load_nodes=network.load_nodes,
                   bounds=region.bounds, budget=budget, method=method, time_limit=time_limit, status='running',
                   cost_unit=network.cost_unit, initial_plan=network.initial_plan, algorithm=self.algorithm,
                   cut_history={str(i+1): dict(cut=cut, scheme=None) for i, cut in enumerate(region.cuts)},
                   sp=0, global_search=0, global_point=None, sp_point=None,
                   **self._geometry(region))

    def initializing(self, index):
        self._emit('mp_start', phase='初始化', direction=index)

    def seed(self, answer, region):
        scheme = self._scheme(answer['x'])
        self._emit('feasible', phase='初始化', seed_point=dict(scheme=scheme, p=answer['p']),
                   **self._geometry(region))

    # 2. 选点；3. 全局搜索：全局见证与后续 SP 支撑点分别存储。
    def selecting(self, count, supporting):
        self._emit('selection', phase='选点', candidate_count=count, supporting=supporting)

    def global_start(self, sp_since_global):
        self._emit('residual_start', phase='全局搜索', global_search=self.state['global_search']+1,
                   sp_since_global=sp_since_global, global_point=None)

    def global_end(self, answer, region):
        point = None if answer['x'] is None else dict(scheme=self._scheme(answer['x']), p=answer['p'])
        self._emit('residual_end', phase='全局搜索', global_point=point,
                   coverage_bound=answer['bound'], coverage_complete=answer['complete'],
                   **self._geometry(region))

    # 4. SP：橙色点表示固定的 (x,p)；只记录原始 eta 和认证状态，不保存运行向量 y。
    def sp_start(self, x, power, number):
        self.sp_numbers[(tuple(x), tuple(power))] = number
        self._emit('point', phase='SP', sp=number, eta=None, feasible=None,
                   sp_point=dict(scheme=self._scheme(x), p=power, number=number))

    def sp_end(self, answer):
        self._emit('sp_end', phase='SP', eta=answer['eta'], feasible=answer['feasible'])

    def updated(self, region, event, *, x, power, checked=None):
        values = self._geometry(region)
        values.update(sp_point=dict(scheme=self._scheme(x), p=power,
                                   number=self.sp_numbers.get((tuple(x), tuple(power)))), phase='SP')
        if checked is not None:
            values.update(eta=checked['eta'], feasible=checked['feasible'])
        if event == 'cut':
            values.update(cut_history={str(len(region.cuts)): dict(cut=checked['cut'], scheme=self._scheme(x))})
        self._emit(event, **values)

    # 5. 构域完成后独立扫描；最终指标、几何、回放只写同一个压缩文件。
    def finish(self, result, region):
        complete = result['certified']
        self._emit('region_end', phase='构域完成' if complete else '构域停止', result=result, coverage_complete=complete,
                   coverage_bound=result['coverage_bound'], **self._geometry(region))

    def scanning(self, completed, total):
        self._validation_update(phase='SOCP 扫描', scan_progress=[completed, total])

    def _validation_update(self, **values):
        """校验只更新独立面板，不占用构域回放帧，也不等待单步按钮。"""
        with self.condition:
            if self.cancelled.is_set():
                raise KeyboardInterrupt()
            values = _plain(values)
            self.validation_state.update(values)
            self.state.update(values)

    def validation(self, reference, result, *, region_key='inner'):
        from region import contains, halfspaces
        states = np.asarray(reference['states'])
        divisions = states.shape[0]
        indices = np.indices(states.shape).reshape(states.ndim, -1).T
        lower = np.asarray(reference.get('axis_lower', np.zeros(states.ndim)))
        points = lower+(indices+.5)*(np.asarray(reference['bounds'])-lower)/divisions
        truth = states.ravel() == 1
        metrics = {}
        for key in (key for key in ('inner', 'outer') if key in result):
            inside = np.zeros(len(points), dtype=bool)
            for row in result[key]:
                inside |= contains(points, halfspaces(row['vertices']))
            missed, extra = np.count_nonzero(truth & ~inside), np.count_nonzero(inside & ~truth)
            metrics[key] = dict(mr_percent=100.*missed/truth.sum() if truth.any() else None,
                               fr_percent=100.*extra/inside.sum() if inside.any() else None,
                               missed_cells=int(missed), extra_cells=int(extra),
                               reference_cells=int(truth.sum()), computed_cells=int(inside.sum()))
        validation = dict(axis_lower=lower, bounds=reference['bounds'], states=states, region_key=region_key,
                          scan_seconds=reference.get('scan_seconds'),
                          metrics=metrics, **metrics[region_key])
        self._validation_update(phase='完成', status='completed', validation=validation)

    def frame(self, index):
        with self.condition:
            state = {}
            for item in self.history[:index+1]:
                _merge(state, item['patch'])
            return state

    def control(self, action):
        with self.condition:
            if action == 'pause':
                self.paused = True
            elif action == 'continue':
                self.paused, self.permits = False, 0
            elif action == 'next':
                self.paused, self.permits = True, 1
            elif action == 'cancel':
                self.cancelled.set()
            self.condition.notify_all()

    def save(self):
        if self.output is None:
            return
        self.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = self.output.with_suffix('.tmp')
        with gzip.open(temporary, 'wt', encoding='utf-8') as stream:
            json.dump(dict(version=4, history=self.history, validation_state=self.validation_state), stream, ensure_ascii=False,
                      allow_nan=False, separators=(',', ':'))
        temporary.replace(self.output)

    def load_recording(self, path):
        with gzip.open(path, 'rt', encoding='utf-8') as stream:
            data = json.load(stream)
        if data['version'] not in (3, 4):
            raise ValueError('只支持原生窗口的 version=3/4 回放')
        self.history = data['history']
        self.validation_state = data.get('validation_state', {})
        if data['version'] == 3:
            # 旧回放只有几何，没有割系数；仅分离校验，不猜测或重建缺失割。
            self.history = []
            for item in data['history']:
                if item['patch'].get('event') in ('scan', 'completed'):
                    self.validation_state.update({k: v for k, v in item['patch'].items() if k != 'event'})
                else:
                    self.history.append(item)
        self.state = self.frame(len(self.history)-1)
        self.state.update(self.validation_state)

    def _work(self, calculate):
        self.busy = True
        try:
            self.result = calculate()
        except BaseException as error:
            self.error = error
            with self.condition:
                patch = dict(status='failed', phase='中断' if isinstance(error, KeyboardInterrupt) else '失败',
                             error=str(error))
                _merge(self.state, patch)
                if self.validation_state:
                    self.validation_state.update(patch)
                else:
                    self.history.append(dict(elapsed=self.clock()-self.started, patch=patch))
        finally:
            self.save()
            self.busy = False

    def execute(self, calculate, *, show_ui=True):
        if show_ui:
            window = NativeWindow(self)
            worker = Thread(target=self._work, args=(calculate,), daemon=True, name='SOCP')
            worker.start()
            window.root.mainloop()
            if worker.is_alive():
                self.control('cancel')
            worker.join()
        else:
            self._work(calculate)
        if self.error is not None:
            raise self.error
        return self.result

    def replay(self):
        NativeWindow(self).root.mainloop()


# 原生界面与绘图均在此文件中；主线无需管理窗口、计时器或帧播放。
INNER, OUTER, GLOBAL, SP = '#397f85', '#a9b9c4', '#cc3838', '#ef8a23'
CUT = '#8055a4'


def _cut_segment(cut, x, bounds, axis_lower=None):
    """alpha + beta @ p + delta @ x = 0 在二维显示框中的截线；只用于绘图。"""
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


def _cut_polygon(cut, x, bounds, axis_lower=None):
    """三维联合割平面与显示盒十二条棱的交点，保留实际 kW 坐标。"""
    from itertools import product
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


def _draw_3d(ax, polytopes, bounds, *, color, fill=False, alpha=1., linestyle='-', linewidth=.8, gid=None):
    """逐凸域画面和真实棱；点、线、面均保留，不跨网架填充凸包。"""
    from scipy.spatial import ConvexHull
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
    for points in polytopes:
        points = np.asarray(points).reshape(-1, 3)
        if not len(points):
            continue
        delta = (points-points[0])/bounds
        rank = np.linalg.matrix_rank(delta, tol=1e-10)
        basis = np.linalg.svd(delta, full_matrices=False)[2][:rank]
        coordinates = delta@basis.T
        if rank == 0:
            ax.scatter(*points[0], color=color, s=14, gid=gid)
            continue
        if rank == 1:
            ends = points[[coordinates[:, 0].argmin(), coordinates[:, 0].argmax()]]
            ax.plot(*ends.T, color=color, lw=linewidth, ls=linestyle, gid=gid)
            continue
        hull = ConvexHull(coordinates/np.linalg.norm(coordinates, axis=0), qhull_options='Qx')
        if rank == 2:
            face = points[hull.vertices]
            faces = [face]
            edges = np.stack([face, np.roll(face, -1, axis=0)], axis=1)
        else:
            faces = points[hull.simplices]
            edges = []
            for i, neighbors in enumerate(hull.neighbors):
                for j in neighbors:
                    if j > i and np.linalg.norm(hull.equations[i]-hull.equations[j]) > 1e-7:
                        edge = np.intersect1d(hull.simplices[i], hull.simplices[j])
                        edges.append(points[edge])
        if fill:
            ax.add_collection3d(Poly3DCollection(faces, facecolors=color, edgecolors='none', alpha=alpha, gid=gid))
        ax.add_collection3d(Line3DCollection(edges, colors=color, linewidths=linewidth,
                            linestyles=linestyle, alpha=min(1., 3*alpha) if fill else alpha, gid=gid))


def _voxel_faces(states, bounds):
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


def _union(polygons):
    from shapely.geometry import MultiPoint
    from shapely.ops import unary_union
    return unary_union([MultiPoint(points).convex_hull for points in polygons if len(points)])


def _draw(ax, geometry, *, color, fill=False, alpha=1., linestyle='-', linewidth=1.):
    """按几何并集绘制，保留不相连部分和孔洞；不跨方案作凸包。"""
    from matplotlib.path import Path as MplPath
    from matplotlib.patches import PathPatch
    if geometry.is_empty:
        return
    if geometry.geom_type in ('MultiPolygon', 'GeometryCollection', 'MultiLineString', 'MultiPoint'):
        for part in geometry.geoms:
            _draw(ax, part, color=color, fill=fill, alpha=alpha, linestyle=linestyle, linewidth=linewidth)
    elif geometry.geom_type == 'Polygon':
        from shapely.geometry.polygon import orient
        poly = orient(geometry, sign=1.)
        paths = []
        for ring in [poly.exterior, *poly.interiors]:
            points = np.asarray(ring.coords)
            codes = [MplPath.MOVETO]+[MplPath.LINETO]*(len(points)-2)+[MplPath.CLOSEPOLY]
            paths.append(MplPath(points, codes))
        path = MplPath.make_compound_path(*paths)
        ax.add_patch(PathPatch(path, facecolor=color if fill else 'none', edgecolor=color,
                              alpha=alpha, linewidth=linewidth, linestyle=linestyle))
    else:
        points = np.asarray(geometry.coords)
        ax.plot(points[:, 0], points[:, 1], color=color, linewidth=linewidth,
                marker='.' if len(points) == 1 else None, linestyle=linestyle, alpha=alpha)


class NativeWindow:
    """A 总域、B 动态网架、C 独立扫描；实时 / 回放共用同一绘图入口。"""

    def __init__(self, monitor, *, root=None, controller=None):
        import tkinter as tk
        from tkinter import ttk
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        import matplotlib as mpl
        mpl.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Microsoft YaHei', 'DejaVu Sans'],
                             'font.size': 9, 'axes.unicode_minus': False, 'svg.fonttype': 'none'})
        self.monitor = monitor
        self.root = tk.Tk() if root is None else root
        self.controller = controller
        self.root.title('SOCP 规划域 · 过程回放')
        self.root.geometry('1420x900')
        self.root.minsize(1000, 700)
        self.root.protocol('WM_DELETE_WINDOW', self.close)
        self.index, self.last_drawn, self.live, self.playing = 0, -1, True, False
        self.last_validation = None
        self.setting_slider = False
        self.scheme_views = {}
        self.envelope_views, self.failed_candidates = {}, []
        self.overview = None
        self.scheme_page, self.focus_scheme = 0, None
        self.status = tk.StringVar(value='初始化')
        self.global_text, self.sp_text = tk.StringVar(value='全局点：—'), tk.StringVar(value='SP 点：—')
        ttk.Label(self.root, textvariable=self.status, font=('Microsoft YaHei', 11)).pack(anchor='w', padx=12, pady=(8, 2))
        point_bar = ttk.Frame(self.root)
        point_bar.pack(fill='x', padx=12)
        ttk.Label(point_bar, textvariable=self.global_text, foreground=GLOBAL).pack(side='left', padx=(0, 28))
        ttk.Label(point_bar, textvariable=self.sp_text, foreground=SP).pack(side='left')
        self.full_extent = tk.BooleanVar(value=False)
        ttk.Checkbutton(point_bar, text='主图展开全局', variable=self.full_extent,
                        command=self.refresh_extent).pack(side='right')
        ttk.Label(point_bar, text='坐标固定于分区初始外包络').pack(side='right', padx=8)
        body = ttk.Panedwindow(self.root, orient='horizontal')
        body.pack(fill='both', expand=True, padx=8, pady=6)
        left, right = ttk.Frame(body), ttk.LabelFrame(body, text="B  网架 · N'_x 绿色认证 / N_x 灰色外包络 / 紫色割")
        body.add(left, weight=1)
        body.add(right, weight=1)
        self.axes, self.canvases = {}, {}
        for panel, title in (('A', 'A  分区内外域 · 全局总览'), ('C', 'C  SOCP 扫描校验')):
            frame = ttk.LabelFrame(left, text=title)
            frame.pack(fill='both', expand=True, pady=2)
            fig = Figure(figsize=(5.8, 3.3), dpi=100)
            ax = fig.add_subplot(111)
            fig.subplots_adjust(left=.13, right=.97, bottom=.19, top=.86)
            canvas = FigureCanvasTkAgg(fig, master=frame)
            canvas.get_tk_widget().pack(fill='both', expand=True)
            self.axes[panel], self.canvases[panel] = ax, canvas
        pages = ttk.Frame(right)
        pages.pack(fill='x', padx=4, pady=3)
        ttk.Button(pages, text='上一页', command=lambda: self.change_page(-1)).pack(side='left')
        ttk.Button(pages, text='下一页', command=lambda: self.change_page(1)).pack(side='left')
        self.page_number = tk.IntVar(value=1)
        self.page_selector = ttk.Spinbox(pages, from_=1, to=1, width=5, textvariable=self.page_number,
                                         command=self.choose_page)
        self.page_selector.pack(side='left', padx=6)
        self.page_selector.bind('<Return>', lambda event: self.choose_page())
        self.page_text = tk.StringVar(value='0 个网架')
        ttk.Label(pages, textvariable=self.page_text).pack(side='left')
        viewport = ttk.Frame(right)
        viewport.pack(fill='both', expand=True)
        self.scroll = tk.Canvas(viewport, highlightthickness=0)
        scrollbar = ttk.Scrollbar(viewport, orient='vertical', command=self.scroll.yview)
        self.scroll.configure(yscrollcommand=scrollbar.set)
        scrollbar.pack(side='right', fill='y')
        self.scroll.pack(side='left', fill='both', expand=True)
        self.scheme_frame = ttk.Frame(self.scroll)
        self.embedded = self.scroll.create_window((0, 0), window=self.scheme_frame, anchor='nw')
        self.scheme_frame.bind('<Configure>', lambda event: self.scroll.configure(scrollregion=self.scroll.bbox('all')))
        self.scroll.bind('<Configure>', lambda event: self.scroll.itemconfigure(self.embedded, width=event.width))
        self.scheme_frame.columnconfigure((0, 1), weight=1)
        controls = ttk.Frame(self.root)
        controls.pack(fill='x', padx=10, pady=(2, 10))
        actions = (('暂停计算', lambda: monitor.control('pause')),
                             ('计算一步', lambda: monitor.control('next')),
                             ('继续计算', lambda: monitor.control('continue')),
                             ('上一帧', lambda: self.step(-1)),
                             ('下一帧', lambda: self.step(1)),
                             ('上一割', lambda: self.seek_cut(-1)),
                             ('下一割', lambda: self.seek_cut(1)),
                   ('播放 / 暂停', self.play), ('实时', self.go_live))
        for text, action in actions if controller is None else actions[3:]:
            ttk.Button(controls, text=text, command=action).pack(side='left', padx=2)
        self.frame_text = tk.StringVar(value='0 / 0')
        ttk.Label(controls, textvariable=self.frame_text, width=15).pack(side='right')
        self.slider = ttk.Scale(controls, from_=0, to=1, command=self.slide)
        self.slider.pack(side='left', fill='x', expand=True, padx=8)
        if controller is None:
            self.root.after(40, self.tick)

    def close(self):
        if self.controller is not None:
            self.controller.close()
            return
        if self.monitor.busy:
            self.monitor.control('cancel')
        self.root.destroy()

    def seek(self, index):
        if self.controller is not None:
            self.controller.seek(index)
            return
        self.live, self.playing = False, False
        self.index = max(0, min(int(index), len(self.monitor.history)-1))
        self.show()

    def slide(self, value):
        if not self.setting_slider:
            self.seek(float(value))

    def step(self, direction):
        index = self.index if self.controller is None else self.controller.index
        self.seek(index+direction)

    def seek_cut(self, direction):
        indices = [i for i, item in enumerate(self.monitor.history)
                   if item['patch'].get('event') == 'cut' and direction*(i-self.index) > 0]
        if indices:
            self.seek(min(indices) if direction > 0 else max(indices))

    def change_page(self, direction):
        self.page_number.set(self.scheme_page+1+direction)
        self.choose_page()

    def choose_page(self):
        import tkinter as tk
        try:
            self.scheme_page = self.page_number.get()-1
        except (ValueError, tk.TclError):
            return
        self._draw_schemes(self.monitor.frame(self.index))

    def play(self):
        if self.controller is not None:
            self.controller.play()
            return
        self.live = False
        if self.index >= len(self.monitor.history)-1:
            self.index = 0
        self.playing = not self.playing
        self.last_drawn = -1

    def go_live(self):
        if self.controller is not None:
            self.controller.seek(self.controller.total-1)
            return
        self.live, self.playing = True, False

    def tick(self):
        total = len(self.monitor.history)
        if total:
            if self.live:
                self.index = total-1
            elif self.playing:
                self.index = min(self.index+1, total-1)
                if self.index == total-1:
                    self.playing = False
            self.show()
        self.root.after(140 if self.playing else 80, self.tick)

    def refresh_extent(self):
        self.last_drawn, self.last_validation = -1, None
        self.show()

    def _outer_polygons(self, state):
        """真实全局外包络；不能用已知网架外域并集替代未知网架的范围。"""
        from region import initial_polytope
        result = state.get('result')
        if result:
            return [row['vertices'] for row in result['outer']]
        if 'global_outer' in state:
            return state['global_outer']
        return [initial_polytope(state['bounds'], state['total_bound'], state['axis_bounds'])*state['bounds']]

    @staticmethod
    def _in_partition(polygons, partition):
        if not partition:
            return polygons
        sign = np.array([1 if s == '+' else -1 for s in partition])
        return [p for p in polygons if len(p) and np.all(np.asarray(p)*sign >= -1e-8)]

    def _frame_view(self):
        """只读当前帧及之前的事件：重建状态、固定包络和失败候选记录。"""
        state, envelopes, failed = {}, {}, {}
        with self.monitor.condition:
            history = self.monitor.history[:self.index+1]
        for item in history:
            patch = item['patch']
            _merge(state, patch)
            if 'bounds' not in state:
                continue
            partition = state.get('partition', '')
            if partition not in envelopes or patch.get('event') == 'initial_bounds':
                envelopes[partition] = self._in_partition(self._outer_polygons(state), partition)
            point, candidate = state.get('sp_point'), state.get('global_point')
            if (patch.get('feasible') is False and point and candidate
                    and point['scheme'] == candidate['scheme'] and point['p'] == candidate['p']):
                key = (partition, point['scheme'], tuple(point['p']))
                failed[key] = dict(point, partition=partition)
        self.envelope_views, self.failed_candidates = envelopes, list(failed.values())
        return state

    def _view_limits(self, state, *, scheme=None, validation=False):
        """分区初始化完成后固定范围；过程点、割和射线不触发坐标缩放。"""
        d = len(state['bounds'])
        groups = []
        if validation:
            groups = [row['vertices'] for key in ('inner', 'outer')
                      for row in state.get('result', {}).get(key, [])]
            reference = state.get('validation')
            if reference:
                cells = np.argwhere(np.asarray(reference['states']) == 1)
                if len(cells):
                    lower = np.asarray(reference.get('axis_lower', np.zeros(d)))
                    widths = (np.asarray(reference['bounds'])-lower)/np.shape(reference['states'])
                    groups.append(lower+np.array([cells.min(axis=0), cells.max(axis=0)+1])*widths)
        elif self.full_extent.get() and scheme is None:
            groups = self._outer_polygons(state)
        else:
            partition = state.get('partition', '')
            if scheme is not None:
                sign = state['schemes'][scheme].get('sign')
                partition = '' if sign is None else ''.join('+' if s > 0 else '-' for s in sign)
            if state.get('result') and scheme is None:
                groups = self._outer_polygons(state)
            else:
                groups = self.envelope_views.get(partition)
                if groups is None:
                    groups = self._in_partition(self._outer_polygons(state), partition)
        points = [np.asarray(group).reshape(-1, d) for group in groups if len(group)]
        if points:
            points = np.vstack(points)
            lower, upper = points.min(axis=0), points.max(axis=0)
        else:
            lower = np.asarray(state.get('axis_lower', np.zeros(d)))
            upper = np.minimum(state['bounds'], state['axis_bounds'])
        padding = .06*np.maximum(upper-lower, 1.)
        return lower-padding, upper+padding

    def _axes(self, ax, state, *, scheme=None, validation=False):
        d = len(state['load_nodes'])
        if d == 3 and ax.name != '3d':
            figure = ax.figure
            ax.remove()
            ax = figure.add_axes([.03, .12, .80, .70], projection='3d')
            ax.view_init(elev=24, azim=-55)
            ax.set_proj_type('ortho')
        view = (ax.elev, ax.azim, ax.roll) if d == 3 else None
        ax.clear()
        lower, limits = self._view_limits(state, scheme=scheme, validation=validation)
        ax.set(xlim=(lower[0], limits[0]), ylim=(lower[1], limits[1]),
               xlabel=f"$p_{{{state['load_nodes'][0]}}}$ (kW)",
               ylabel=f"$p_{{{state['load_nodes'][1]}}}$ (kW)")
        if d == 3:
            ax.set(zlim=(lower[2], limits[2]), zlabel=f"$p_{{{state['load_nodes'][2]}}}$ (kW)")
            ax.set_box_aspect((1., 1., .8))
            ax.view_init(*view)
            ax.tick_params(pad=0)
            for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
                axis.labelpad = 2
        else:
            ax.set_aspect('auto')  # 刻度仍是物理 kW。
            ax.spines[['top', 'right']].set_visible(False)
        ax.tick_params(labelsize=8)
        ax.grid(alpha=.15, linewidth=.5)
        return ax

    def _markers(self, ax, state, scheme=None, *, partition=None, size_scale=1.):
        """过程图层统一供二维、三维和各网架使用，标记不参与视口计算。"""
        from matplotlib.lines import Line2D
        handles = []
        if scheme is not None:
            sign = state['schemes'][scheme].get('sign')
            partition = None if sign is None else ''.join('+' if s > 0 else '-' for s in sign)

        def scatter(points, label, color, marker, gid, size=32, hollow=False):
            points = np.asarray(points).reshape(-1, len(state['bounds']))
            if partition:
                sign = np.array([1 if s == '+' else -1 for s in partition])
                points = points[np.all(points*sign >= -1e-8, axis=1)]
            if not len(points):
                return
            style = dict(facecolors='none', edgecolors=color) if hollow else dict(c=color)
            ax.scatter(*points.T, marker=marker, s=size*size_scale, zorder=10, gid=gid, **style)
            handles.append(Line2D([], [], color=color, marker=marker, ls='none',
                                  markerfacecolor='none' if hollow else color, label=label))

        scatter([row['p'] for row in state.get('rejected_points', [])],
                '全网架不可行', GLOBAL, 'o', 'global-rejected')
        failed = [row['p'] for row in self.failed_candidates
                  if scheme is None or row['scheme'] == scheme]
        scatter(failed, '全局候选：该网架不可行', GLOBAL, 'x', 'candidate-rejected', 42)
        scatter(state.get('unknown_points', []), '未决点', '#777777', 'x', 'unknown-points')
        if state.get('phase') == '全网架认证' and state.get('query_point') is not None:
            scatter([state['query_point']], '全网架待认证点', GLOBAL, 'o', 'query-point', 65, True)
        for key, color, marker, label in (('seed_point', INNER, '^', '初始可行点'),
                                         ('global_point', GLOBAL, 'D', '全局搜索点'),
                                         ('sp_point', SP, 'o', 'SP / 支撑点')):
            point = state.get(key)
            if point and (scheme is None or point['scheme'] == scheme):
                if key == 'sp_point' and state.get('feasible') is False:
                    marker, label = 'x', 'SP 不可行点'
                scatter([point['p']], label, color, marker, key, 40)
        ray = state.get('ray')
        if (ray and (state.get('event') in ('ray_start', 'ray_end')
                     or state.get('phase') in ('首轮射线认证', '射线补点', '边界补充'))
                and (scheme is None or ray['scheme'] == scheme)):
            ray_partition = self._in_partition([[ray['anchor'], ray['target']]], partition)
            if ray_partition:
                ax.plot(*np.asarray(ray_partition[0]).T, color=INNER, ls=':', lw=1.2, gid='ray-path')
                scatter([ray['anchor']], '射线锚点', INNER, 'o', 'ray-anchor', 24, True)
                scatter([ray['target']], '射线目标', INNER, 'x', 'ray-target', 42)
                if ray['p'] is not None:
                    ax.plot(*np.asarray([ray['anchor'], ray['p']]).T, color=INNER, lw=1.5, gid='ray-certified')
                    scatter([ray['p']], '射线认证点', INNER, 's', 'ray-point', 38)
        return handles

    def show(self):
        if not self.monitor.history:
            return
        self._refresh_validation()
        if self.index == self.last_drawn:
            return
        self.last_drawn = self.index
        state = self._frame_view()
        point = (state.get('global_point') if state.get('event') in ('residual_start', 'residual_end')
                 else state.get('sp_point') or state.get('seed_point'))
        focus = state.get('active_scheme') or (point['scheme'] if point else None)
        if focus and focus != self.focus_scheme:
            self.focus_scheme = focus
            keys = list(state.get('schemes', {}))
            if focus in keys:
                self.scheme_page = keys.index(focus)//4
        self._show_header(state)
        if 'bounds' in state:
            self._draw_union(state)
            self._draw_schemes(state)

    def _refresh_validation(self):
        """校验独立于历史帧；停止拖动时也刷新，回退后仍保留最终对比。"""
        with self.monitor.condition:
            latest = dict(self.monitor.state)
        if 'bounds' not in latest:
            return
        key = (tuple(latest['axis_bounds']), tuple(latest.get('scan_progress', [])),
               id(latest.get('result')), id(latest.get('validation')), latest.get('error'))
        if key != self.last_validation:
            self.last_validation = key
            self._draw_validation(latest)

    def _show_header(self, state):
        d = len(state.get('load_nodes', (0, 1)))
        measure, unit = ('面积', 'kW²') if d == 2 else ('体积', 'kW³')
        self.root.title(f"{state.get('algorithm', '主线')} · {state.get('network', '')} · {d} 维过程回放")
        total = len(self.monitor.history)
        self.frame_text.set(f'{self.index+1} / {total}')
        self.setting_slider = True
        self.slider.configure(to=max(1, total-1 if self.controller is None else self.controller.total-1))
        self.slider.set(self.index if self.controller is None else self.controller.index)
        self.setting_slider = False
        phase = state.get('phase', '初始化')
        detail = ''
        if phase == '初始化' and 'direction' in state:
            index = state['direction']
            detail = f" · MP2 {'p'+str(state['load_nodes'][index]) if index < len(state['load_nodes']) else '总负荷'}"
        elif phase == '选点':
            detail = f" · {state['candidate_count']} 个{'支撑点' if state['supporting'] else '顶点'}"
        elif phase == 'SP' and state.get('eta') is not None:
            detail = f" · η={state['eta']:.3g} · {'认证' if state['feasible'] else '不可行'}"
            if state.get('event') == 'cut':
                detail += ' · 采用割'
        elif phase == '获选点生成割':
            detail = f" · 最大评分 η={state['eta']:.3g} · 只解该点的取割 LP"
        elif phase == '射线筛选':
            detail = (f" · {measure}收益上界 {state['ray_gain_bound']:.4g} < "
                      f"{state['ray_gain_threshold']:.4g} {unit}，免求解")
        elif phase == '全局搜索':
            detail = ' · 在尚未覆盖区域中搜索完整物理可行点'
            if state.get('event') == 'residual_end' and state['coverage_bound'] is not None:
                detail += f" · 全局上界 {state['coverage_bound']:.6g}"
        elif phase == '全网架认证':
            detail = ' · 固定 p=('+', '.join(f'{p:.3f}' for p in state['query_point'])+')，x/y 自由'
            if state.get('certification_eta') is not None:
                detail += f" · η={state['certification_eta']:.3g}"
        elif phase == '网架支撑':
            detail = f" · 固定 {state['support_scheme']}，p/y 自由 · 方向 {state['support_direction']}"
        elif phase in ('首轮射线认证', '射线补点', '边界补充'):
            detail = f" · 固定 {state['active_scheme']}，从认证锚点向候选点求边界"
            if state.get('event') == 'ray_end':
                detail += f" · λ={state['ray_fraction']:.6f}"
        elif phase == '网架阶段结束':
            detail = ' · '+{'covered': '当前网架外顶点已认证，进入连续查漏',
                           'area_stagnation': '边界补充完成，进入连续查漏',
                           'point_resolution': '近点处理完成，进入连续查漏'}[state['stage_reason']]
        if state.get('active_scheme'):
            detail += f" · 阶段 {state.get('stage', 0)} / 网架 {state['active_scheme']}"
            if state.get('area_ratio') is not None:
                detail += f" · 切割{measure} {100*state['area_ratio']:.3f}% · 连续 {state['small_cuts']}/{state['patience']}"
        if state.get('error'):
            detail = ' · '+state['error']
        self.status.set(f"{phase}{detail}    SP 累计 {state.get('sp', 0)} · 全局 #{state.get('global_search', 0)} · 割 {state.get('cuts', 0)}")
        for key, text, label in (('global_point', self.global_text, '全局点'), ('sp_point', self.sp_text, 'SP 点')):
            point = state.get(key)
            if point and key == 'sp_point' and point.get('number') is not None:
                label += f" #{point['number']}"
            coordinates = ', '.join(f'{p:.3f}' for p in point['p']) if point else ''
            text.set(label+'：—' if not point else f"{label}：{point['scheme']}  p=({coordinates}) kW")

    def _draw_cuts(self, ax, state, scheme, *, history=True):
        """同一联合割按当前面板的 x 代入；总图仅展示来源网架的一条截线。"""
        from matplotlib.collections import LineCollection
        cuts = state.get('cut_history', {})
        if not cuts or scheme not in state.get('schemes', {}):
            return
        sign = state['schemes'][scheme].get('sign')
        cuts = {key: row for key, row in cuts.items() if np.array_equal(row.get('sign'), sign)}
        if not cuts:
            return
        latest = next(reversed(cuts))
        intervals = [ax.get_xlim(), ax.get_ylim()]
        if len(state['bounds']) == 3:
            intervals.append(ax.get_zlim())
        lower, limits = np.asarray(intervals).T.copy()
        if sign is not None:
            lower = np.where(np.asarray(sign) > 0, np.maximum(lower, 0.), lower)
            limits = np.where(np.asarray(sign) < 0, np.minimum(limits, 0.), limits)
            if np.any(lower > limits):
                return
        if len(limits) == 3:
            # 历史割已体现在 Nx 的棱面；当前割另画平面，避免几十个透明面遮住认证域。
            polygon = _cut_polygon(cuts[latest]['cut'], state['schemes'][scheme]['x'], limits, lower)
            _draw_3d(ax, [polygon], state['bounds'], color=CUT, fill=True, alpha=.10,
                     linewidth=1., gid=f'cut-{latest}-{scheme}')
            ax.text2D(.02, .97, f'割 #{latest} · 网架 {scheme}', transform=ax.transAxes,
                      color=CUT, fontsize=8, va='top')
            return
        previous = []
        for number, item in cuts.items():
            if not history and number != latest:
                continue
            segment = _cut_segment(item['cut'], state['schemes'][scheme]['x'], limits, lower)
            if len(segment) != 2:
                if number == latest:
                    cut = np.asarray(item['cut'])
                    corners = np.array([lower, [limits[0], lower[1]], [lower[0], limits[1]], limits])
                    slack = cut[0]+corners@cut[1:3]+cut[3:]@state['schemes'][scheme]['x']
                    message = ('本框全部排除' if slack.max() < 0. else
                               '本框均满足此割' if slack.min() >= 0. else '仅与角点相交')
                    ax.text(.02, .94, f'割 #{number}：{message}', transform=ax.transAxes,
                            color=CUT, fontsize=8, va='top')
                continue
            active = number == latest
            if active:
                ax.plot(*segment.T, color=CUT, linewidth=1.6, zorder=5, gid=f'cut-{number}-{scheme}')
                midpoint = segment.mean(axis=0)
                ax.annotate(f'#{number}', midpoint, xytext=(4, 4), textcoords='offset points',
                            color=CUT, fontsize=8, zorder=9)
            else:
                previous.append(segment)
        ax.add_collection(LineCollection(previous, colors=CUT, linestyles='--',
                                         linewidths=.7, alpha=.28, zorder=4))

    def _removed(self, state, scheme=None):
        """只在加割帧显示真实前后差集；无需另存被切掉的多边形。"""
        if state.get('event') != 'cut':
            return _union([])
        before = self.monitor.frame(self.index-1).get('schemes', {})
        after = state['schemes']
        keys = before.keys() if scheme is None else [scheme]
        old = _union([before[key]['outer'] for key in keys if key in before])
        new = _union([after[key]['outer'] for key in keys if key in after])
        return old.difference(new)

    def _draw_polygons(self, ax, polygons, state, *, gid=None, **style):
        """面板共用几何入口；只在此处分派二维并集和三维凸域绘制。"""
        if len(state['bounds']) == 3:
            _draw_3d(ax, polygons, state['bounds'], gid=gid, **style)
        else:
            before = set(ax.get_children())
            _draw(ax, _union(polygons), **style)
            for artist in set(ax.get_children())-before:
                artist.set_gid(gid)

    def _draw_removed(self, ax, state, scheme=None):
        if len(state['bounds']) == 3:
            self._draw_polygons(ax, self._removed_3d(state, scheme), state,
                                color=CUT, fill=True, alpha=.18, gid='removed-region')
        else:
            _draw(ax, self._removed(state, scheme), color=CUT, fill=True, alpha=.22)

    def _draw_union(self, state):
        from matplotlib.lines import Line2D
        ax = self.axes['A'] = self._axes(self.axes['A'], state)
        d = len(state['bounds'])
        partition = None if self.full_extent.get() or state.get('result') else state.get('partition')
        outer = self._in_partition(self._outer_polygons(state), partition)
        rows = list(state.get('schemes', {}).values())
        inner = self._in_partition([row['inner'] for row in rows], partition)
        known = self._in_partition([row['outer'] for row in rows], partition)
        if d == 2:
            self._draw_polygons(ax, outer, state, color=OUTER, fill=True, alpha=.10, gid='global-outer-fill')
        self._draw_polygons(ax, outer, state, color=OUTER, linestyle='--', linewidth=1.3, gid='global-outer')
        if not state.get('result'):
            self._draw_polygons(ax, known, state, color=OUTER, fill=d == 2, alpha=.20, gid='known-outer')
        self._draw_polygons(ax, inner, state, color=INNER, fill=True,
                            alpha=.15 if d == 3 else .38, linewidth=.8, gid='certified-inner')
        self._draw_removed(ax, state)
        cuts = state.get('cut_history', {})
        if cuts:
            self._draw_cuts(ax, state, next(reversed(cuts.values()))['scheme'], history=False)
        handles = [Line2D([], [], color=OUTER, ls='--', label=f'分区 {partition} 外包络' if partition else '全局外包络 G'),
                   Line2D([], [], color=INNER, label='认证内域 G′')]
        if not state.get('result'):
            handles.append(Line2D([], [], color=OUTER, label='网架条件外域'))
        handles.extend(self._markers(ax, state, partition=partition))
        ax.legend(handles=handles, loc='lower left', bbox_to_anchor=(-.08, 1.02),
                  ncol=3, frameon=False, fontsize=7, columnspacing=.8, handlelength=1.5)
        self._draw_overview(state)
        self.canvases['A'].draw_idle()

    def _draw_overview(self, state):
        """总览始终显示真实全局外包络；主图窗口用蓝框标出。"""
        from itertools import product
        d = len(state['bounds'])
        detail = self.axes['A']
        if self.overview is None:
            self.overview = detail.figure.add_axes([.75, .45, .23, .30],
                                                   projection='3d' if d == 3 else None)
        ax = self.overview
        ax.clear()
        ax.set_visible(not self.full_extent.get())
        if self.full_extent.get():
            detail.set_position([.03, .12, .80, .65] if d == 3 else [.13, .19, .84, .59])
            return
        detail.set_position([.02, .12, .68, .65] if d == 3 else [.13, .19, .54, .59])
        outer = self._outer_polygons(state)
        points = [p for p in outer if len(p)]
        points = np.vstack(points) if points else np.zeros((1, d))
        lower, upper = points.min(axis=0), points.max(axis=0)
        padding = .06*np.maximum(upper-lower, 1.)
        ax.set(xlim=(lower[0]-padding[0], upper[0]+padding[0]),
               ylim=(lower[1]-padding[1], upper[1]+padding[1]))
        if d == 3:
            ax.set(zlim=(lower[2]-padding[2], upper[2]+padding[2]))
            ax.set_proj_type('ortho')
            ax.view_init(24, -55)
            ax.set_box_aspect((1., 1., .8))
        self._draw_polygons(ax, outer, state, color=OUTER, linestyle='--', linewidth=.7, gid='overview-outer')
        self._draw_polygons(ax, [row['inner'] for row in state.get('schemes', {}).values()], state,
                            color=INNER, fill=True, alpha=.25, linewidth=.4, gid='overview-inner')
        intervals = [detail.get_xlim(), detail.get_ylim()]
        if d == 3:
            intervals.append(detail.get_zlim())
        box = np.array(list(product(*intervals)))
        self._draw_polygons(ax, [box], state, color='#427cac', linewidth=1.1, gid='overview-viewport')
        self._markers(ax, state, size_scale=.25)
        write = ax.text2D if d == 3 else ax.text
        write(.5, 1.10, '全局外包络 G\n蓝框：主图范围', transform=ax.transAxes,
              ha='center', va='bottom', fontsize=7)
        ax.set_xticks([])
        ax.set_yticks([])
        if d == 3:
            ax.set_zticks([])

    def _draw_schemes(self, state):
        from tkinter import ttk
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        schemes = state.get('schemes', {})
        pages = max(1, (len(schemes)+3)//4)
        self.scheme_page = max(0, min(self.scheme_page, pages-1))
        self.page_number.set(self.scheme_page+1)
        self.page_selector.configure(to=pages)
        self.page_text.set(f'/ {pages} 页 · 共 {len(schemes)} 个网架')
        visible = list(schemes)[self.scheme_page*4:(self.scheme_page+1)*4]
        for key, (frame, ax, canvas) in self.scheme_views.items():
            if key not in visible:
                frame.grid_remove()
        for index, key in enumerate(visible):
            row = schemes[key]
            if key not in self.scheme_views:
                if len(self.scheme_views) == 4:
                    reusable = next(old for old in self.scheme_views if old not in visible)
                    self.scheme_views[key] = self.scheme_views.pop(reusable)
                else:
                    frame = ttk.LabelFrame(self.scheme_frame, text=key)
                    fig = Figure(figsize=(3., 3.0), dpi=100)
                    ax = fig.add_subplot(111)
                    fig.subplots_adjust(left=.20, right=.96, bottom=.18, top=.73)
                    canvas = FigureCanvasTkAgg(fig, master=frame)
                    canvas.get_tk_widget().pack(fill='both', expand=True)
                    self.scheme_views[key] = frame, ax, canvas
            frame, ax, canvas = self.scheme_views[key]
            frame.grid(row=index//2, column=index%2, sticky='nsew', padx=3, pady=3)
            active = ' · 当前认证' if key == state.get('active_scheme') else ''
            frame.configure(text=f"{key} · {row['cost']:g} {state.get('cost_unit', '')}{active}")
            ax = self._axes(ax, state, scheme=key)
            self.scheme_views[key] = frame, ax, canvas
            self._draw_polygons(ax, [row['outer']], state, color=OUTER, fill=True,
                                alpha=.16, linestyle='--', gid='scheme-outer')
            self._draw_polygons(ax, [row['inner']], state, color=INNER, fill=True,
                                alpha=.23 if len(state['bounds']) == 3 else .45, gid='scheme-inner')
            self._draw_removed(ax, state, key)
            self._draw_cuts(ax, state, key)
            self._markers(ax, state, key)
            from textwrap import fill
            initial = state.get('initial_plan', {})
            changes = [f"{edge}:{'断开' if kind is None else '并联' if kind == 'parallel' else '接入' if kind == 'existing' else kind}"
                       for edge, kind in row['choice'].items() if initial.get(edge) != kind]
            label = fill(' · '.join(changes) if changes else '原始网架', width=38)
            ax.set_title(label, fontsize=8, pad=6)
            canvas.draw_idle()

    def _draw_validation(self, state):
        from matplotlib.lines import Line2D
        ax = self.axes['C'] = self._axes(self.axes['C'], state, validation=True)
        d = len(state['bounds'])
        validation = state.get('validation')
        if validation is None:
            progress = state.get('scan_progress')
            message = (state.get('validation_note', '构域完成后独立扫描') if progress is None else
                       f'SOCP 扫描  {progress[0]} / {progress[1]}  ({100*progress[0]/progress[1]:.1f}%)')
            if state.get('error'):
                message = ('校验未完成' if progress is not None else '构域未完成')+'\n'+state['error']
            write = ax.text2D if d == 3 else ax.text
            write(.5, .55, message, transform=ax.transAxes, ha='center', color='#777777')
        elif d == 3:
            self._draw_validation_3d(ax, state)
        else:
            states = np.asarray(validation['states'])
            bounds, n = np.asarray(validation['bounds']), states.shape[0]
            lower = np.asarray(validation.get('axis_lower', np.zeros(2)))
            centers = [a+(np.arange(n)+.5)*(b-a)/n for a, b in zip(lower, bounds)]
            # 参考域按真实扫描格显示；曲线为算法的连续多边形，不平滑扫描结果。
            from matplotlib.colors import ListedColormap
            ax.pcolormesh(np.linspace(lower[0], bounds[0], n+1), np.linspace(lower[1], bounds[1], n+1),
                          (states.T == 1).astype(int), cmap=ListedColormap(['white', '#d5e3ea']),
                          vmin=0, vmax=1, shading='flat', rasterized=True)
            key = validation.get('region_key', 'inner')
            predicted = _union([row['vertices'] for row in state['result'][key]])
            _draw(ax, predicted, color=INNER, linewidth=1.5)
            if key == 'inner' and 'outer' in validation.get('metrics', {}):
                _draw(ax, _union([row['vertices'] for row in state['result']['outer']]),
                      color=OUTER, linestyle='--', linewidth=1.)
            ax.contour(*centers, (states.T == 1).astype(float), levels=[.5], colors=['#526c80'], linewidths=.7)
            fmt = lambda value: '—' if value is None else f'{value:.3f}%'
            title = f"遗漏 {fmt(validation['mr_percent'])}    多余 {fmt(validation['fr_percent'])}    · {n}×{n} 网格"
            if key == 'inner' and 'outer' in validation.get('metrics', {}):
                outer = validation['metrics']['outer']
                title = f"G′：{title}\nG：遗漏 {fmt(outer['mr_percent'])}    多余 {fmt(outer['fr_percent'])}"
            ax.set_title(title, fontsize=9)
            handles = [Line2D([], [], color='#526c80', label='SOCP 扫描参考'),
                       Line2D([], [], color=INNER, label='认证内域 G′' if key == 'inner' else '联合割外域并集')]
            if key == 'inner' and 'outer' in validation.get('metrics', {}):
                handles.append(Line2D([], [], color=OUTER, ls='--', label='外包络 G'))
            ax.legend(handles=handles, loc='upper right', frameon=False, fontsize=8)
        self.canvases['C'].draw_idle()

    def _removed_3d(self, state, scheme=None):
        """当前割在各网架旧 Nx 中实际切掉的部分，仅用于本帧显示。"""
        from region import clip_polytope
        if state.get('event') != 'cut':
            return []
        before = self.monitor.frame(self.index-1)['schemes']
        latest = next(reversed(state['cut_history'].values()))
        cut = np.asarray(latest['cut'])
        bounds = np.asarray(state['bounds'])
        keys = before if scheme is None else [scheme]
        return [clip_polytope(np.asarray(before[key]['outer'])/bounds,
                             -cut[0]-cut[4:]@before[key]['x'], -cut[1:4]*bounds)*bounds
                for key in keys if key in before and np.array_equal(before[key].get('sign'), latest.get('sign'))]

    def _draw_validation_3d(self, ax, state):
        from matplotlib.lines import Line2D
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
        from region import contains, halfspaces
        validation, result = state['validation'], state['result']
        states, bounds = np.asarray(validation['states']), np.asarray(validation['bounds'])
        lower = np.asarray(validation.get('axis_lower', np.zeros(3)))
        ax.set_position([.03, .22, .80, .60])
        ax.add_collection3d(Poly3DCollection(lower+_voxel_faces(states, bounds-lower), facecolors='#7e9bae',
                                            edgecolors='none', alpha=.14, gid='scan-reference'))
        _draw_3d(ax, [row['vertices'] for row in result['outer']], bounds, color=OUTER,
                 linestyle='--', alpha=.75)
        _draw_3d(ax, [row['vertices'] for row in result['inner']], bounds, color=INNER,
                 fill=True, alpha=.06, linewidth=.3)
        points = lower+(np.indices(states.shape).reshape(3, -1).T+.5)*(bounds-lower)/states.shape[0]
        inside = np.zeros(len(points), dtype=bool)
        for row in result['inner']:
            inside |= contains(points, halfspaces(row['vertices']))
        missed = points[(states.ravel() == 1) & ~inside]
        ax.scatter(*missed.T, color=GLOBAL, s=7, marker='.', depthshade=False, gid='missed-cells')
        fmt = lambda value: '—' if value is None else f'{value:.3f}%'
        labels = []
        for key, symbol in (('inner', "G'"), ('outer', 'G')):
            metrics = validation['metrics'][key]
            labels.append(f"{symbol}：遗漏 {fmt(metrics['mr_percent'])} · 多余 {fmt(metrics['fr_percent'])}")
        ax.set_title('\n'.join(labels)+f"  · {states.shape[0]}³ 网格", fontsize=8, pad=8)
        ax.legend(handles=[Line2D([], [], color='#7e9bae', label='SOCP 扫描'),
                           Line2D([], [], color=INNER, label="$G'$"),
                           Line2D([], [], color=OUTER, ls='--', label='$G$'),
                           Line2D([], [], color=GLOBAL, marker='.', ls='none', label='G′ 遗漏点')],
                  loc='lower left', bbox_to_anchor=(-.1, -.27), ncol=4, frameon=False, fontsize=8)


class SynchronizedReplay:
    """两窗口共用帧索引与播放时钟，短轨迹停在末帧。"""

    def __init__(self, recordings):
        import tkinter as tk
        self.root = tk.Tk()
        self.index, self.playing = 0, False
        self.windows = []
        width = max(760, min(1120, self.root.winfo_screenwidth()//2))
        height = min(820, self.root.winfo_screenheight()-80)
        for index, path in enumerate(recordings):
            monitor = RunMonitor()
            monitor.load_recording(path)
            root = self.root if index == 0 else tk.Toplevel(self.root)
            window = NativeWindow(monitor, root=root, controller=self)
            window.live = False
            window.root.minsize(760, 600)
            window.root.geometry(f'{width}x{height}+{index*width}+20')
            self.windows.append(window)
        self.total = max(len(w.monitor.history) for w in self.windows)
        self.seek(0)
        self.root.after(180, self.tick)

    def seek(self, index):
        self.playing = False
        self.index = max(0, min(int(index), self.total-1))
        self.show()

    def show(self):
        for window in self.windows:
            window.index = min(self.index, len(window.monitor.history)-1)
            window.last_drawn = -1
            window.show()

    def play(self):
        if self.index == self.total-1:
            self.index = 0
        self.playing = not self.playing

    def tick(self):
        if self.playing:
            self.index = min(self.index+1, self.total-1)
            self.show()
            if self.index == self.total-1:
                self.playing = False
        self.root.after(180, self.tick)

    def close(self):
        self.root.destroy()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='原生窗口回放，不运行优化器')
    parser.add_argument('recording', nargs='?', type=Path,
                        default=Path(__file__).resolve().parent/'results'/'fourbus_2d'/'monitor.json.gz')
    parser.add_argument('--compare', type=Path, help='第二份轨迹；打开两窗口同步逐帧回放')
    args = parser.parse_args()
    if args.compare is not None:
        SynchronizedReplay([args.recording, args.compare]).root.mainloop()
    else:
        monitor = RunMonitor()
        monitor.load_recording(args.recording)
        monitor.replay()
