"""结果整理、存取与展示：比较标签、表面网格、终端进度、HTML 和回放。"""
from dataclasses import dataclass
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from hashlib import sha256
from pathlib import Path
from threading import Event, RLock, Thread
from time import perf_counter
from urllib.parse import urlsplit, parse_qs
import json
import math
import webbrowser
import base64
import gzip

import numpy as np

ROOT = Path(__file__).resolve().parent
METHODS = ('socp', 'hybrid', 'ac', 'linear')
METHOD_NAMES = ('纯 SOCP 切割', '线性 + SOCP 精修', 'AC 数值参考', '纯线性切割')


def json_value(value):
    """网页采用严格 JSON：无穷预算和非有限求解界使用 null。"""
    if isinstance(value, np.ndarray):
        return json_value(value.tolist())
    if isinstance(value, np.generic):
        return json_value(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {key: json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_value(item) for item in value]
    return value


def subtract_polytope(poly, equations):
    """返回 poly 减去一个凸域的互不重叠内部的凸分块。"""
    from region import clip_polytope, polytope_volume
    if not len(poly):
        return []
    values = poly@equations[:, :3].T+equations[:, 3]
    if np.all(values <= 1e-10):
        return []
    if np.any(values.min(axis=0) >= -1e-10):
        return [poly]
    pieces = []
    for eq in equations:
        values = poly@eq[:3]+eq[3]
        if values.max() <= 1e-10:
            continue
        outside = clip_polytope(poly, eq[3], eq[:3])
        if polytope_volume(outside) > 1e-15:
            pieces.append(outside)
        poly = clip_polytope(poly, -eq[3], -eq[:3])
        if polytope_volume(poly) <= 1e-15:
            break
    return pieces


def union_volume(polytopes):
    """逐域扣除已经计入的部分，避免把重叠方案体积相加。"""
    from region import polytope_volume, halfspaces
    previous, volume = [], 0.
    for poly in sorted(polytopes, key=polytope_volume, reverse=True):
        if polytope_volume(poly) <= 1e-15:
            continue
        pieces = [poly]
        for eq in previous:
            pieces = [part for p in pieces for part in subtract_polytope(p, eq)]
            if not pieces:
                break
        volume += sum(polytope_volume(p) for p in pieces)
        previous.append(halfspaces(poly))
    return volume


def sample_region(result, points, bounds):
    """连续域求出后才作采样；薄层保留为未知，采样不影响几何。"""
    from region import contains, halfspaces
    points = np.atleast_2d(points)/np.asarray(bounds)
    inside, possible = np.zeros(len(points), bool), np.zeros(len(points), bool)
    for poly in result['inner']:
        inside |= contains(points, halfspaces(np.asarray(poly['vertices'])/bounds))
    for poly in result['outer']:
        possible |= contains(points, halfspaces(np.asarray(poly['vertices'])/bounds))
    return np.where(inside, 1, np.where(possible, 0, -1)).astype(np.int8)


def region_metrics(result, bounds):
    """从最终顶点按需计算并集体积，不进入求解循环或持久化记录。"""
    bounds = np.asarray(bounds)
    inner, outer = (union_volume([np.asarray(p['vertices'])/bounds for p in result[kind]])
                    for kind in ('inner', 'outer'))
    scale = np.prod(bounds)
    return dict(inner_volume=inner*scale, outer_volume=outer*scale,
                volume_gap=(outer-inner)/outer if outer else 0.)


def disagreement(approximation, ac):
    """同一等体积网格上的 FR/MR；空域的条件比例无定义。"""
    # 1. 统计多余、遗漏及并集单元数
    extra = int(np.count_nonzero(approximation & ~ac))
    missed = int(np.count_nonzero(ac & ~approximation))
    computed, reference = int(approximation.sum()), int(ac.sum())
    union = computed+missed
    # 2. FR 以计算域为分母，MR 以 AC 域为分母
    return dict(fr_percent=100*extra/computed if computed else None,
                mr_percent=100*missed/reference if reference else None,
                region_error_percent=100*(extra+missed)/union if union else 0.,
                computed_cells=computed, ac_cells=reference, extra_cells=extra, missed_cells=missed)


def disagreement_interval(approximation, ac):
    """输入 -1/0/1 分别表示已证域外、未确定、已证域内；区间不含网格离散误差。"""
    # 1. 根据未知单元的最有利/最不利归属，求 |left\right| / |left| 的界
    def ratio_bounds(left, right):
        inside, unknown = left==1, left==0
        certain = np.count_nonzero(inside & (right==-1))
        lower_denominator = inside.sum()+np.count_nonzero(unknown & (right!=-1))
        possible = np.count_nonzero((left!=-1) & (right!=1))
        upper_denominator = inside.sum()+np.count_nonzero(unknown & (right!=1))
        return [100*certain/lower_denominator if lower_denominator else 0.,
                100*possible/upper_denominator if upper_denominator else 100.]
    # 2. 分别计算 FR、MR；完整空域的条件比例保持无定义
    unknown = int(np.count_nonzero((approximation==0)|(ac==0)))
    exact = disagreement(approximation==1,ac==1) if unknown==0 else dict(fr_percent=None,mr_percent=None,region_error_percent=None)
    fr = ratio_bounds(approximation,ac) if unknown else ([exact['fr_percent']]*2 if exact['fr_percent'] is not None else None)
    mr = ratio_bounds(ac,approximation) if unknown else ([exact['mr_percent']]*2 if exact['mr_percent'] is not None else None)
    # 3. 计算对称差 / 并集的误差区间
    certain_difference = np.count_nonzero(approximation*ac == -1)
    possible_overlap = np.count_nonzero((approximation != -1) & (ac != -1))
    certain_overlap = np.count_nonzero((approximation == 1) & (ac == 1))
    possible_difference = np.count_nonzero(~(((approximation == 1) & (ac == 1)) |
                                            ((approximation == -1) & (ac == -1))))
    low_den, high_den = certain_difference+possible_overlap, certain_overlap+possible_difference
    error = [100*certain_difference/low_den if low_den else 0.,
             100*possible_difference/high_den if high_den else 0.]
    return dict(**exact,fr_interval=fr,mr_interval=mr,region_error_interval=error,unknown_cells=unknown)


def comparison_labels(states):
    """0 域外、1 多余、2 遗漏、3 重合、4 未确定。"""
    inside = states == 1
    ac = METHODS.index('ac')
    labels = inside.astype(np.uint8)+2*inside[ac].astype(np.uint8)
    labels[(states == 0) | (states[ac] == 0)] = 4
    return labels


def validation_summary(states, metadata):
    """由保留内域和 AC 三态网格按需推导比较指标。"""
    times = {(r['method'], r['budget']): r['timing']['total_seconds']
             for r in metadata.get('continuous', [])}
    return [dict(method=method, budget=budget,
                 **disagreement_interval(np.where(states[k, j] == 1, 1, -1),
                                         states[METHODS.index('ac'), j]),
                 method_unknown_cells=int(np.count_nonzero(states[k, j] == 0)),
                 metric_scope='retained_inner_union',
                 total_seconds=times.get((method, budget), metadata['seconds'][method]))
            for k, method in enumerate(METHODS) for j, budget in enumerate(metadata['budgets'])]


@dataclass
class BenchmarkResult:
    """只保存三态网格与元数据；比较标签、指标和图形按需推导。"""

    states: np.ndarray  # (方法, 预算, N, N, N)，-1 域外、0 未确定、1 域内。
    metadata: dict

    @classmethod
    def create(cls, network, budgets, divisions, bounds, **settings):
        from model import PLANNING_TOL
        from region import GEOMETRY_TOL
        from vertify import AC_TOL, FIXED_POINT_TOL, GLOBAL_AC_TOL
        metadata = dict(network=network.name, planning=True, cost_unit=network.cost_unit,
                        load_nodes=list(network.load_nodes), corridor_count=network.n_corridors,
                        type_count=network.n_types, upgrade_count=getattr(network, 'upgrade_count', None),
                        schema_version=3, network_fingerprint=network.fingerprint,
                        required_nodes=[node for node, required in zip(network.nodes, network.required) if required],
                        corridors=[dict(id=c.id, endpoints=c.endpoints,
                            types=[t.id for t in c.types]) for c in network.corridors],
                        budgets=json_value(budgets), divisions=divisions, bounds=bounds.tolist(),
                        seconds={}, continuous=[],
                        tolerances=dict(planning=PLANNING_TOL, geometry=GEOMETRY_TOL,
                                         ac=AC_TOL, fixed_point=FIXED_POINT_TOL,
                                         ac_global=GLOBAL_AC_TOL), **settings)
        return cls(np.zeros((len(METHODS), len(budgets))+(divisions,)*3, dtype=np.int8), metadata)

    def add_region(self, region, index):
        n = self.metadata['divisions']
        points = (np.indices((n,)*3).reshape(3, -1).T+.5)*self.spacing
        method = region['method']
        self.states[METHODS.index(method), index] = sample_region(region, points, self.bounds).reshape((n,)*3)
        record = {key: value for key, value in region.items() if key not in ('tau', 'residual_mode')}
        self.metadata['continuous'].append(json_value(dict(record, budget_index=index)))
        seconds = self.metadata['seconds']
        seconds[method] = seconds.get(method, 0.)+region['timing']['total_seconds']

    def add_validation(self, states, seconds):
        self.states[METHODS.index('ac')] = states
        self.metadata['seconds']['ac'] = seconds

    @property
    def budgets(self):
        return tuple(np.inf if b is None else b for b in self.metadata['budgets'])

    @property
    def bounds(self):
        return np.asarray(self.metadata['bounds'])

    @property
    def load_nodes(self):
        return tuple(self.metadata['load_nodes'])

    @property
    def spacing(self):
        return self.bounds/self.metadata['divisions']

    @property
    def labels(self):
        return comparison_labels(self.states)

    @property
    def summary(self):
        return validation_summary(self.states, self.metadata)

    def save(self, folder):
        """完整写入临时文件后替换，避免中断破坏上一次结果。"""
        folder = Path(folder)
        folder.mkdir(parents=True, exist_ok=True)
        temporary = folder/'result.npz.tmp'
        metadata = json.dumps(json_value(self.metadata), ensure_ascii=False, allow_nan=False)
        with temporary.open('wb') as stream:
            np.savez_compressed(stream, states=self.states, metadata=metadata)
        temporary.replace(folder/'result.npz')

    @classmethod
    def load(cls, folder, network=None):
        with np.load(Path(folder)/'result.npz', allow_pickle=False) as data:
            result = cls(data['states'], json.loads(str(data['metadata'])))
        if network is not None and (result.metadata.get('schema_version') != 3
                or result.metadata.get('network_fingerprint') != network.fingerprint):
            raise ValueError('Saved result does not match the current network and model schema; recompute it')
        return result


def surface(poly, bounds):
    """将归一化凸域转换为 HTML 表面；退化集保留顶点，不伪造三角面。"""
    from region import _convex_hull

    poly = np.asarray(poly, dtype=float).reshape(-1, 3)
    full_dimension = len(poly) >= 4 and np.linalg.matrix_rank(poly-poly[0], tol=1e-10) == 3
    faces = _convex_hull(poly).simplices.tolist() if full_dimension else []
    return dict(vertices=(poly*np.asarray(bounds)).tolist(), faces=faces)


def region_geometry(records, bounds, cache=None):
    """按方案生成展示数据；缓存只保留当前版本，与算法的几何缓存分开。"""
    cache = {} if cache is None else cache
    bounds = np.asarray(bounds, dtype=float)
    current, geometry = {}, []
    for row in records:
        key = tuple(row['x'])
        inner, outer = (np.asarray(row[k], dtype=float).reshape(-1, 3) for k in ('inner', 'outer'))
        signature = (bounds.tobytes(), tuple(row['choice'].items()), row['cost'], inner.tobytes(), outer.tobytes())
        saved = cache.get(key)
        if saved is None or saved[0] != signature:
            saved = (signature, dict(x=list(map(int, key)), choice=dict(row['choice']), cost=row['cost'],
                                     inner=surface(inner, bounds), outer=surface(outer, bounds)))
        current[key] = saved
        geometry.append(saved[1])
    cache.clear()
    cache.update(current)
    return geometry


def region_view(result, bounds):
    """补齐回放需要的表面，兼容已有带三角面的结果，不修改原始证书。"""
    view = json_value(result)
    if 'geometry' not in view:
        if 'certificates' in view:  # 历史结果读取。
            view['geometry'] = region_geometry(view['certificates'], bounds)
        else:
            view['geometry'] = [dict(choice=p['choice'], cost=p['cost'],
                                     inner=surface(np.asarray(p['vertices'])/bounds, bounds),
                                     outer=dict(vertices=[], faces=[])) for p in view['inner']]
    for kind in ('inner', 'outer'):
        view[kind] = [poly if 'faces' in poly else
                      dict(surface(np.asarray(poly['vertices']).reshape(-1, 3)/bounds, bounds),
                           **poly) for poly in view[kind]]
    return view


def view_data(data, bounds, cache=None):
    """统一转换实时事件与离线快照，原始算法数据不参与页面缓存。"""
    view = dict(data)
    records = view.pop('records', None)
    if records is not None:
        view['geometry'] = region_geometry(records, bounds, cache)
    if view.get('region') is not None:
        view['region'] = region_view(view['region'], bounds)
        view.update(geometry=view['region']['geometry'], global_outer=view['region']['outer'])
    if 'results' in view:
        view['results'] = [region_view(region, bounds) for region in view['results']]
    return view


def pack_replay(value):
    """无损共享重复的几何与历史对象，压缩 HTML 内的回放数据。"""
    objects, indices = [], {}
    def encode(item):
        if isinstance(item, dict):
            node = [1, {key: encode(v) for key, v in item.items()}]
        elif isinstance(item, list):
            node = [0, [encode(v) for v in item]]
        else:
            return item
        key = json.dumps(node, ensure_ascii=False, allow_nan=False, separators=(',', ':'))
        if key not in indices:
            indices[key] = len(objects)
            objects.append(node)
        return {'ref': indices[key]}
    root = encode(value)
    return dict(packed_replay_version=1, root=root, objects=objects)


def cut_slice(cut, selection, bounds):
    """联合割在生成方案 x 下的截面；只求平面与评价箱的交多边形。"""
    cut, selection, bounds = map(np.asarray, (cut, selection, bounds))
    normal = cut[1:4]
    constant = float(cut[0] + cut[4:] @ selection)
    # 在单位立方体中相交，避免三个负荷轴的量级差影响容差。
    scaled = normal * bounds
    norm = np.linalg.norm(scaled)
    if norm <= 1e-30:
        return dict(constant=constant, normal=normal.tolist(), vertices=[])
    scaled, offset = scaled/norm, constant/norm
    corners = np.array([[i, j, k] for i in (0., 1.) for j in (0., 1.) for k in (0., 1.)])
    vertices = []
    for i, a in enumerate(corners):
        for b in corners[i+1:]:
            if np.count_nonzero(a != b) != 1:
                continue
            va, vb = offset+scaled@a, offset+scaled@b
            if abs(va) < 1e-10:
                vertices.append(a)
            if abs(vb) < 1e-10:
                vertices.append(b)
            if va*vb < 0:
                vertices.append(a + va/(va-vb)*(b-a))
    if not vertices:
        return dict(constant=constant, normal=normal.tolist(), vertices=[])
    vertices = np.unique(np.round(vertices, 12), axis=0)
    if len(vertices) >= 3:
        center = vertices.mean(axis=0)
        u = vertices[0]-center
        u /= np.linalg.norm(u)
        v = np.cross(scaled, u)
        angles = np.arctan2((vertices-center)@v, (vertices-center)@u)
        vertices = vertices[np.argsort(angles)]
    return dict(constant=constant, normal=normal.tolist(), vertices=(vertices*bounds).tolist())


class RunMonitor:
    """监视数据独立于模型；record 控制完整历史，show_ui 控制实时服务。"""

    def __init__(self, *, show_ui=False, output=None, open_browser=True,
                 heartbeat_seconds=10., stream=None, record=None):
        import sys
        # 1. 显示、记录和终端输出配置
        self.show_ui = show_ui
        self.record = show_ui if record is None else bool(record)
        self.output = Path(output) if output is not None else None
        self.open_browser = open_browser
        self.heartbeat_seconds = heartbeat_seconds
        self.stream = sys.stdout if stream is None else stream
        # 2. 计时与实时服务
        self.started = perf_counter()
        self.finished = None
        self.last_event = self.last_print = self.started
        self.lock = RLock()
        self.stop = Event()
        self.server = self.server_thread = self.heartbeat_thread = None
        self.url = None
        # 3. 当前状态、增量历史与展示几何缓存
        self.events, self.history = [], []
        self.geometry_cache = {}
        self.previous = {}
        self.state = dict(status='running', event='start', message='准备启动',
                          query=0, iteration=0, total_cuts=0, pool_size=0,
                          method_number=0, method_count=4, revision=0, methods=[])

    def __enter__(self):
        try:
            if self.show_ui:
                self.start_server()
            self.heartbeat_thread = Thread(target=self.heartbeat, name='planning-progress', daemon=True)
            self.heartbeat_thread.start()
        except BaseException:
            self.close()
            raise
        return self

    def __exit__(self, kind, error, traceback):
        try:
            if error is not None:
                self('interrupted' if isinstance(error, KeyboardInterrupt) else 'failed',
                     message='用户中断计算' if isinstance(error, KeyboardInterrupt) else f'计算失败：{error}')
                if self.record and self.output is not None:
                    try:
                        self.save_snapshot()
                    except OSError as snapshot_error:
                        print(f'无法保存监视快照：{snapshot_error}', file=self.stream, flush=True)
        finally:
            self.close()

    def __call__(self, event, **data):
        with self.lock:
            # 1. 仅在记录/显示时转换展示几何，终端模式不做绘图计算
            if event in ('method_start', 'phase_start'):
                self.geometry_cache.clear()
            if self.record or self.show_ui:
                data = view_data(data, data.get('bounds', self.state.get('bounds')), self.geometry_cache)
            else:
                data.pop('records', None)
            # 2. 更新计时与求解阶段，清空上一阶段的临时状态
            now = perf_counter()
            step_seconds = now-self.last_event
            self.last_event = now
            self.state.update(event=event, message=data.get('message', event),
                              revision=self.state['revision']+1)
            if event in ('completed', 'failed', 'interrupted'):
                self.state['status'] = event
                self.finished = now
            if event == 'method_start':
                self.state.update(query=0, iteration=0, total_cuts=0, pool_size=0,
                                  point=None, choice=None, latest_cut=None, counts=[], phase='',
                                  states=None, lower=None, upper=None, bound=None, objective=None)
                self.state.update(geometry=[], global_outer=None, coverage_complete=False, region=None,
                                  coverage_bound=None, max_total=None, max_total_bound=None, mode=None,
                                  counts_algorithm={})
            if event in ('phase_start', 'point', 'sp_skip'):
                self.state.update(iteration=0, choice=None,
                                  bound=None, objective=None, eta=None, query_status=None)
                if data.get('mode') != 'SP-gap-support':
                    self.state.update(covered_by=None, covered_cost=None, uncovered_witness=None)
            if event == 'phase_start':
                self.state.update(query=0, point=None, lower=None, upper=None, latest_cut=None)
                self.state.update(geometry=[], global_outer=None, coverage_complete=False, region=None,
                                  states=None, counts=[], counts_algorithm={}, coverage_bound=None,
                                  max_total=None, max_total_bound=None, mode=None)
                if not self.record:
                    self.state.pop('states', None)
            if event == 'query_end':
                self.state['query_status'] = data.get('status')
            # 3. 接收数值结果，按需生成三态计数与割平面
            for key, value in data.items():
                if key not in ('states', 'cut', 'selection', 'status'):
                    self.state[key] = json_value(value)
            if 'states' in data:
                states = np.asarray(data['states'])
                counts = [dict(inside=int(np.count_nonzero(s == 1)),
                               outside=int(np.count_nonzero(s == -1)),
                               unknown=int(np.count_nonzero(s == 0)), total=int(s.size)) for s in states]
                self.state['counts'] = counts
                self.state['divisions'] = states.shape[-1]
                if self.record:
                    self.state['states'] = states.reshape(len(states), -1).tolist()
            if event == 'cut':
                self.state['total_cuts'] += 1
                if self.record and self.state.get('bounds') is not None:
                    sliced = cut_slice(data['cut'], data['selection'], self.state['bounds'])
                    self.state['latest_cut'] = dict(**sliced, point=json_value(data['point']),
                                                    joint_coefficients=json_value(data['cut']),
                                                    selection=json_value(data['selection']),
                                                    choice=self.state.get('choice'),
                                                    iteration=self.state['iteration'],
                                                    number=self.state['total_cuts'])
            if event == 'method_end':
                self.state['methods'] = [*self.state['methods'], dict(method=self.state['method'], seconds=data['seconds'],
                                                  counts=self.state.get('counts', []), cuts=self.state['total_cuts'])]
            # 4. 只在内存中追加变化帧，完成后统一嵌入回放 HTML
            if self.record:
                item = {key: self.state.get(key) for key in
                        ('revision', 'method', 'phase', 'query', 'iteration', 'event', 'message', 'point', 'query_status')}
                item['elapsed'] = now-self.started
                if event == 'cut':
                    item['cut'] = self.state.get('latest_cut')
                self.events.append(item)
                self.state.update(elapsed=now-self.started, step_seconds=step_seconds)
                patch = {k: v for k, v in self.state.items() if k not in self.previous or self.previous[k] != v}
                frame = dict(id=len(self.history), elapsed=now-self.started, patch=patch)
                self.history.append(frame)
                self.previous = dict(self.state)
            # 5. 短查询合并为每秒一次的终端进度
            if event in ('start', 'preparation', 'method_start', 'method_end', 'phase_start', 'phase_end',
                         'saving', 'plotting', 'loaded', 'completed', 'failed', 'interrupted') \
                    or (event == 'query_end' and data.get('status') == 'unknown') \
                    or now-self.last_print >= 1.:
                self.print_status()

    def print_status(self, *, heartbeat=False):
        now = perf_counter()
        state = self.state
        prefix = f"[{now-self.started:8.1f}s]"
        if state.get('method'):
            prefix += f" [{state['method']}/{state.get('phase', '')}]"
        details = []
        if state.get('query'):
            details.append(f"查询 #{state['query']} / MP 轮次 {state['iteration']}")
        if state.get('point') is not None:
            details.append('p=('+', '.join(f'{v:.2f}' for v in state['point'])+') kW')
        counts = state.get('counts', [])
        if counts:
            total = sum(c['total'] for c in counts)
            unknown = sum(c['unknown'] for c in counts)
            details.append(f'已分类 {total-unknown}/{total} ({100*(total-unknown)/total:.1f}%) / 未确定 {unknown}')
        details.append(f"新增割 {state['total_cuts']} / 割池 {state['pool_size']}")
        if state['event'] == 'method_end':
            details.append(f"方法耗时 {state['seconds']:.3f}s")
        if heartbeat:
            details.append(f'此步骤已等待 {now-self.last_event:.1f}s')
        print(f"{prefix} {state['message']} | {' | '.join(details)}", file=self.stream, flush=True)
        self.last_print = now

    def heartbeat(self):
        while not self.stop.wait(self.heartbeat_seconds):
            with self.lock:
                if self.state['status'] == 'running' and perf_counter()-self.last_print >= self.heartbeat_seconds:
                    self.print_status(heartbeat=True)

    def snapshot(self, after=0):
        with self.lock:
            now = perf_counter() if self.finished is None else self.finished
            value = dict(self.state, elapsed=now-self.started,
                         step_seconds=self.state.get('step_seconds', 0.) if self.finished else now-self.last_event,
                         events=self.events[-150:], history=self.history[max(0, after):],
                         history_total=len(self.history), recording_version=1)
            return json.dumps(value, ensure_ascii=False, allow_nan=False).encode('utf-8')

    def load_recording(self, path):
        """从自包含 HTML 恢复完整历史，不调用优化器。"""
        # 1. 读取 HTML 中唯一一份压缩回放
        html = Path(path).read_text(encoding='utf-8')
        encoded = html.split('window.SAVED_REPLAY_GZIP="', 1)[1].split('";', 1)[0]
        value = json.loads(gzip.decompress(base64.b64decode(encoded)))
        # 2. 大记录按对象引用恢复共享几何
        if 'packed_replay_version' in value:
            objects = []
            def resolve(item):
                return objects[item['ref']] if isinstance(item, dict) else item
            for kind, node in value['objects']:
                objects.append({k: resolve(v) for k, v in node.items()} if kind else [resolve(v) for v in node])
            value = resolve(value['root'])
        # 3. 恢复历史、事件列表和计时
        self.history = value.pop('history')
        self.events = value.pop('events', [])
        self.state = value
        self.finished = perf_counter()
        self.started = self.finished-self.state.get('elapsed', 0.)
        self.previous = dict(self.state)

    def start_server(self):
        # 使用现有 Plotly 的本地脚本；页面离线可用，无 CDN 或 Node 服务。
        from plotly.offline import get_plotlyjs
        self.plotly = get_plotlyjs().encode('utf-8')
        self.template = (ROOT/'live_view.html').read_text(encoding='utf-8')
        owner = self

        class Handler(BaseHTTPRequestHandler):
            def do_GET(self):
                path = urlsplit(self.path).path
                if path == '/':
                    payload, mime = owner.template.encode('utf-8'), 'text/html; charset=utf-8'
                elif path == '/state':
                    try:
                        after = int(parse_qs(urlsplit(self.path).query).get('after', ['0'])[0])
                    except ValueError:
                        after = 0
                    payload, mime = owner.snapshot(after), 'application/json; charset=utf-8'
                elif path == '/plotly.min.js':
                    payload, mime = owner.plotly, 'application/javascript; charset=utf-8'
                elif path.startswith('/result/') and owner.output is not None:
                    name = path[len('/result/'):]
                    if '/' in name or '\\' in name or not (name == 'region_comparison.html' or
                            (name.startswith('methods_') and name.endswith('.html'))):
                        self.send_error(404)
                        return
                    try:
                        payload = (owner.output/name).read_bytes()
                    except OSError:
                        self.send_error(404)
                        return
                    mime = 'text/html; charset=utf-8'
                else:
                    self.send_error(404)
                    return
                self.send_response(200)
                self.send_header('Content-Type', mime)
                self.send_header('Content-Length', str(len(payload)))
                self.send_header('Cache-Control', 'no-store')
                self.end_headers()
                try:
                    self.wfile.write(payload)
                except (BrokenPipeError, ConnectionResetError):
                    pass  # 关闭浏览器不影响计算。

            def log_message(self, format, *args):
                pass

        self.server = ThreadingHTTPServer(('127.0.0.1', 0), Handler)
        self.server.daemon_threads = True
        self.url = f'http://127.0.0.1:{self.server.server_port}/'
        self.server_thread = Thread(target=self.server.serve_forever, name='planning-view', daemon=True)
        self.server_thread.start()
        print(f'实时监视：{self.url}', file=self.stream, flush=True)
        if self.open_browser:
            try:
                webbrowser.open(self.url)
            except (OSError, webbrowser.Error) as error:
                print(f'自动打开浏览器失败，可手动访问上方地址：{error}', file=self.stream, flush=True)

    def save_snapshot(self):
        if self.output is None:
            return
        # 1. 从当前数值结果生成展示快照及本地页面资源
        with self.lock:
            self.state = view_data(self.state, self.state.get('bounds'), self.geometry_cache)
        self.output.mkdir(parents=True, exist_ok=True)
        if not hasattr(self, 'template'):
            from plotly.offline import get_plotlyjs
            self.template = (ROOT/'live_view.html').read_text(encoding='utf-8')
            self.plotly = get_plotlyjs().encode('utf-8')
        self.state['viewer_sha256'] = sha256(self.template.encode('utf-8')).hexdigest()
        payload = self.snapshot().decode('utf-8')
        # 2. 大记录先共享跨帧重复几何，再压缩完整回放
        if len(payload) > 10_000_000:
            payload = json.dumps(pack_replay(json.loads(payload)), ensure_ascii=False,
                                 allow_nan=False, separators=(',', ':'))
        data = base64.b64encode(gzip.compress(payload.encode('utf-8'), compresslevel=6, mtime=0)).decode('ascii')
        # 3. 数据与绘图库嵌入同一 HTML，不再另存 JSON/JSONL
        html = self.template.replace('<script src="/plotly.min.js"></script>',
                                     '<script>'+self.plotly.decode('utf-8')+'</script>')
        html = html.replace('/*SNAPSHOT*/', 'window.SAVED_REPLAY_GZIP="'+data+'";')
        path = self.output/'live_view.html'
        path.write_text(html, encoding='utf-8')
        print(f'完整过程回放（{len(self.history)} 条事件）：{path}', file=self.stream, flush=True)

    def close(self):
        self.stop.set()
        if self.server is not None:
            if self.server_thread is not None and self.server_thread.is_alive():
                self.server.shutdown()
                self.server_thread.join(timeout=2.)
            self.server.server_close()
        if self.heartbeat_thread is not None:
            self.heartbeat_thread.join(timeout=2.)


    def show_result(self, result, *, export=True, keep_ui=False):
        print_summary(result)
        if export or self.show_ui:
            exporting = perf_counter()
            self.save_snapshot()
            print(f'回放导出耗时：{perf_counter()-exporting:.3f}s', file=self.stream, flush=True)
        if self.show_ui and keep_ui:
            try:
                input('计算已完成，网页可继续查看；按 Enter 或 Ctrl+C 关闭监视服务。\n')
            except (EOFError, KeyboardInterrupt):
                pass


def print_summary(result):
    """终端显示同一原始结果推导出的 FR/MR；未知时显示区间。"""
    names = dict(zip(METHODS, METHOD_NAMES))
    print('\n预算 | 方法 | 状态 | 最大负荷 (kW) | FR (%) | MR (%) | 耗时 (s)', flush=True)
    continuous = {(r['method'], r['budget']): r for r in result.metadata.get('continuous', [])}
    for row in result.summary:
        rates = []
        for key in ('fr', 'mr'):
            value = row[f'{key}_percent']
            interval = row.get(f'{key}_interval')
            rates.append(f'{value:.5f}' if value is not None else
                         ('—' if interval is None else f'[{interval[0]:.5f}, {interval[1]:.5f}]'))
        budget = '无限' if row['budget'] is None else str(row['budget'])
        region = continuous.get((row['method'], row['budget']))
        seconds = region['timing']['total_seconds'] if region else row['total_seconds']
        status = region['status'] if region else 'AC 采样'
        maximum = f"{region['max_total']:.3f}" if region and region['max_total'] is not None else '—'
        print(f"{budget} | {names[row['method']]} | {status} | {maximum} | {' | '.join(rates)} | {seconds:.3f}", flush=True)


def voxel_surface(mask, spacing):
    """合并同一平面的暴露体素面，保留非凸边界和孔洞；坐标为 kW。"""
    # 1. 沿三个坐标轴提取正负方向的暴露面
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
            # 2. 将每个平面的连续单元合并为互不重叠的矩形
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
                        # 3. 恢复空间坐标并统一朝外法向，用两个三角形表示矩形
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
    import plotly.graph_objects as go

    # 1. 建立四方法面板与统一的区域类别
    fig = make_subplots(rows=2, cols=2, specs=[[{"type": "scene"}]*2]*2,
                        subplot_titles=[f"{letter}  {name}" for letter, name in zip("abcd", METHOD_NAMES)],
                        horizontal_spacing=.02, vertical_spacing=.08)
    categories = [(3, "与 AC 重合", "#4286AD", 1., "common"),
                  (2, "遗漏", "#D43D3D", 1., "missed"),
                  (1, "多余", "#E9B72F", 1., "extra")]
    categories.append((4, "未确定", "#A6A6A6", .22, "unknown"))
    # 2. 同一份比较标签生成各类别表面，空类别只保留图例
    for panel, labels in enumerate(result.labels[:, budget_index]):
        for code, name, color, opacity, group in categories:
            mesh = voxel_surface(labels == code, result.spacing)
            points = np.asarray(mesh["vertices"]).reshape(-1, 3)
            faces = np.asarray(mesh["triangles"], dtype=int).reshape(-1, 3)
            if len(points):
                trace = go.Mesh3d(x=points[:, 0], y=points[:, 1], z=points[:, 2],
                                 i=faces[:, 0], j=faces[:, 1], k=faces[:, 2],
                                 name="AC 参考域" if panel == 2 and code==3 else name,
                                 color=color, opacity=opacity, flatshading=True,
                                 legendgroup=group, showlegend=panel == 0,
                                 lighting=dict(ambient=1., diffuse=.25, specular=.05, roughness=.95),
                                 hovertemplate="节点负荷 (%{x:.2f}, %{y:.2f}, %{z:.2f}) kW<extra>%{fullData.name}</extra>")
            else:
                trace = go.Scatter3d(x=[None], y=[None], z=[None], mode="markers", name=name,
                                    marker=dict(color=color, size=7), legendgroup=group,
                                    showlegend=panel == 0, hoverinfo="skip")
            fig.add_trace(trace, row=panel//2+1, col=panel % 2+1)
    # 3. 统一负荷坐标、观察角度和图层控制
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
        for name, groups in [("全部", {"common", "missed", "extra", "unknown"}),
                             ("仅差异", {"missed", "extra"}),
                             ("计算域", {"common", "extra"}),
                             ("AC 参考域", {"common", "missed"}),
                             ("未确定", {"unknown"})]]
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


def plot_continuous_regions(result, budget_index):
    """直接显示各方案认证内域及全局外包络；不依赖采样网格。"""
    from plotly.subplots import make_subplots
    import plotly.graph_objects as go

    # 1. 选择当前预算的连续域结果
    methods = ('linear', 'socp', 'hybrid')
    names = dict(zip(METHODS, METHOD_NAMES))
    fig = make_subplots(rows=1, cols=3, specs=[[dict(type='scene')]*3],
                        subplot_titles=[f'{letter}  {names[m]}' for letter, m in zip('abc', methods)])
    rows = {r['method']: r for r in result.metadata['continuous'] if r['budget_index'] == budget_index}
    # 2. 绘制方案内域与全局外包络，统一负荷轴及观察角度
    for col, method in enumerate(methods, 1):
        view = region_view(rows[method], result.bounds)
        for kind, color, opacity, label in [('outer', '#7D9CB4', .15, '外包络'),
                                           ('inner', '#338D80', .65, '认证内域')]:
            for poly in view[kind]:
                vertices = np.asarray(poly['vertices']).reshape(-1, 3)
                faces = np.asarray(poly['faces'], dtype=int).reshape(-1, 3)
                if not len(vertices):
                    continue
                name = f"方案 {poly['choice']} · 投资 {poly['cost']:g}" if 'choice' in poly else label
                fig.add_trace(go.Mesh3d(x=vertices[:, 0], y=vertices[:, 1], z=vertices[:, 2],
                              i=faces[:, 0], j=faces[:, 1], k=faces[:, 2], name=name,
                              color=color, opacity=opacity, flatshading=True, showlegend=False), row=1, col=col)
        axes = {axis+'axis': dict(title=f'节点 {node} 负荷 (kW)', range=[0, bound], nticks=4)
                for axis, node, bound in zip('xyz', result.load_nodes, result.bounds)}
        fig.update_scenes(**axes, aspectmode='cube',
                          camera=dict(eye=dict(x=1.5, y=1.6, z=1.2), projection=dict(type='orthographic')),
                          row=1, col=col)
    # 3. 整理版面并返回可复用的图对象
    fig.update_layout(height=460, template='plotly_white', margin=dict(l=0, r=0, t=45, b=0),
                      font=dict(family='Arial, Microsoft YaHei, sans-serif', size=11))
    return fig


def show_result(result, folder, open_browser=True):
    """输出核心指标；需要界面时生成自包含 HTML，无后台服务或事件记录。"""
    print_summary(result)
    if open_browser:
        path = save_method_comparison(result, folder).resolve()
        print(f'最终结果：{path}', flush=True)
        webbrowser.open(path.as_uri())


def save_method_comparison(result, folder):
    """按预算导出交互图与比较表；指标由同一原始结果推导。"""
    import plotly.io as pio

    # 1. 设置各预算的页面名称及联动视角
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    filenames = ['region_comparison.html']+[
        f'methods_{index}.html' for index in range(1, len(result.budgets))]
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
    # 2. 按需计算比较指标与图注
    rows = result.summary
    domains = {(r['method'], r['budget']): r for r in result.metadata.get('continuous', [])}
    scope = '可规划域' if result.metadata['planning'] else '固定方案可调度域截面'
    spacing = ' × '.join(f'{value:g}' for value in result.spacing)
    cost_unit = result.metadata['cost_unit']
    geometry_note = "上排绿色为认证内域、浅蓝为全局外包络，下排为独立 AC 网格比较。灰色为未确定；网格完整不表示连续边界精确。构域耗时按当前预算列出，AC 时间覆盖全部预算。"
    for index, (budget, filename) in enumerate(zip(result.budgets, filenames)):
        # 3. 生成本预算的导航与四方法指标表
        links = ' · '.join(
            f'<a href="{name}" aria-current="{"page" if i == index else "false"}">'
            f'{"无限预算" if np.isinf(b) else f"{b:g} {cost_unit}"}</a>'
            for i, (b, name) in enumerate(zip(result.budgets, filenames)))
        navigation = f'<nav>预算：{links}</nav>' if result.metadata['planning'] else ''
        table = []
        for method, name in zip(METHODS, METHOD_NAMES):
            row = next(row for row in rows if row['method'] == method
                       and row['budget'] == (None if np.isinf(budget) else budget))
            def rate(key):
                if row[key] is not None:
                    return f'{row[key]:.4f}%'
                interval = row.get(key.replace('_percent','_interval'))
                return '—' if interval is None else f'[{interval[0]:.3f}, {interval[1]:.3f}]%'
            domain = domains.get((method, row['budget']))
            status = domain['status'] if domain else 'AC 采样'
            maximum = f"{domain['max_total']:.3f}" if domain and domain['max_total'] is not None else '—'
            counts = f"{domain['counts']['sp']} / {domain['counts']['cuts']}" if domain else '—'
            table.append(f'<tr><td>{name}</td><td>{status}</td><td>{maximum}</td>'
                         f'<td>{rate("region_error_percent")}</td><td>{rate("fr_percent")}</td>'
                         f'<td>{rate("mr_percent")}</td><td>{row["total_seconds"]:.3f}</td>'
                         f'<td>{counts}</td><td>{row.get("unknown_cells",0)}</td></tr>')
        # 4. 嵌入连续域和 AC 网格图，统一保存为离线页面
        chart = pio.to_html(plot_method_comparison(result, index), include_plotlyjs=not bool(result.metadata.get('continuous')), full_html=False,
                           div_id='method-comparison', post_script=synchronize,
                           config=dict(responsive=True, displaylogo=False, scrollZoom=True))
        continuous = ''
        if result.metadata.get('continuous'):
            continuous = pio.to_html(plot_continuous_regions(result, index), include_plotlyjs=False,
                                     full_html=False, config=dict(responsive=True, displaylogo=False))
            # 先加载绘图库，再显示连续域，最后显示 AC 采样比较。
            from plotly.offline import get_plotlyjs
            continuous = '<script>'+get_plotlyjs()+'</script>'+continuous
        html = f'''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{result.metadata['network']} · {scope}</title><style>
body{{font-family:Arial,"Microsoft YaHei",sans-serif;color:#263641;margin:12px auto;max-width:1200px;padding:0 12px}}
nav{{text-align:center;padding:8px}}a{{color:#426985;text-decoration:none}}a[aria-current=page]{{font-weight:700;text-decoration:underline}}
table{{border-collapse:collapse;margin:10px auto;font-size:14px;min-width:520px}}td,th{{padding:7px 18px;text-align:right;border-bottom:1px solid #dce2e7}}td:first-child,th:first-child{{text-align:left}}
p{{font-size:12px;color:#56616a;line-height:1.7;text-align:center}}
</style></head><body>{navigation}{continuous}{chart}
<table><thead><tr><th>方法</th><th>认证状态</th><th>最大负荷 kW</th><th>区域误差</th><th>FR</th><th>MR</th><th>耗时 s</th><th>SP / 割</th><th>AC 未确定</th></tr></thead><tbody>{''.join(table)}</tbody></table>
<p>图：{result.metadata['network']} 的{scope}；节点负荷单位为 kW，各面板使用同一组刻度。<br>
网格步长 {spacing} kW；蓝色为重合部分，AC 面板显示独立参考；红色遗漏，黄色多余。<br>
FR = 多余 / 计算域；MR = 遗漏 / AC 域。AC 自比较的零仅表示它是基准；有限网格上的零不表示连续误差严格为零。<br>
构域时间含建模、求解及几何计算，混合方法计入线性阶段；AC 时间单列。绘图和导出不计。<br>
{geometry_note}</p>
</body></html>'''
        (folder/filename).write_text(html, encoding='utf-8')
    return folder/filenames[0]


def render_survey(result, output):
    """勘察绘图唯一入口：读取结果字典，统一生成总图和条件价值图。"""
    # 1. 本次绘图的依赖、输出目录和样式
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle, Patch, PathPatch, Rectangle
    from matplotlib.path import Path as MplPath
    from matplotlib.ticker import MaxNLocator
    from shapely.geometry import shape
    from shapely.geometry.polygon import orient
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    WIDTH_MM = 183
    GREEN, GREEN_LIGHT = '#23845D', '#B5D9C6'
    RED, RED_LIGHT = '#BC4D46', '#F2C9C3'
    DARK, MUTED, BASE, CONFIRMED, UNKNOWN = '#283C49','#6E808E','#A7B2BA','#ADCADB','#EDF1F5'
    COLORS = dict(zip('ABCDEF',['#3F6E9E','#A87B42','#7864A5','#B0576A','#32888D','#8A939B']))
    MARKERS = dict(zip('ABCDEF',['o','s','D','^','v','x']))
    STYLE = {'font.family':'sans-serif',
        'font.sans-serif':['Microsoft YaHei','Arial','DejaVu Sans'],
        'font.size':7.5,'axes.labelsize':8.,'axes.titlesize':8.,
        'xtick.labelsize':7.,'ytick.labelsize':7.,'legend.fontsize':7.,
        'svg.fonttype':'none','pdf.fonttype':42,'mathtext.fontset':'dejavusans',
        'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.65,
        'legend.frameon':False,'savefig.facecolor':'white','figure.facecolor':'white',
        'hatch.linewidth':.35}

    # 2. 共用的多边形、坐标轴及价值曲线绘制
    def polygon_parts(geometry):
        if geometry.is_empty:
            return []
        if geometry.geom_type == 'Polygon':
            return [geometry]
        return [part for piece in geometry.geoms for part in polygon_parts(piece)]


    def fill_geometry(ax, geometry, color, *, edgecolor='none', hatch=None, zorder=1):
        """按环方向填充多边形，保留孔洞；斜线由 Matplotlib 统一绘制。"""
        for polygon in polygon_parts(geometry):
            polygon = orient(polygon, sign=1.)
            paths = []
            for ring in [polygon.exterior, *polygon.interiors]:
                vertices = np.asarray(ring.coords)
                codes = np.full(len(vertices), MplPath.LINETO, dtype=np.uint8)
                codes[0], codes[-1] = MplPath.MOVETO, MplPath.CLOSEPOLY
                paths.append(MplPath(vertices, codes))
            ax.add_patch(PathPatch(MplPath.make_compound_path(*paths), facecolor=color,
                                   edgecolor=edgecolor, hatch=hatch, linewidth=.35, zorder=zorder))


    def boundary(ax, geometry, color, *, linestyle='-', linewidth=1.1, zorder=4):
        for polygon in polygon_parts(geometry):
            coordinates = np.asarray(polygon.exterior.coords)
            ax.plot(coordinates[:, 0], coordinates[:, 1], color=color, linestyle=linestyle,
                    linewidth=linewidth, zorder=zorder)


    def domain_axes(ax):
        ax.set(xlim=(0,100),ylim=(0,100),aspect='equal',xlabel=r'$p_1$ (kW)',ylabel=r'$p_2$ (kW)')
        ax.set_xticks([0,50,100])
        ax.set_yticks([0,50,100])
        ax.tick_params(length=2.4,pad=2)


    def plot_value_history(ax, result):
        """每条道路只显示尚未勘察时的条件评分；观测顺序来自 trace。"""
        # 1. 绘制各道路评分及数值区间
        for route in result['protocol']['routes']:
            rows = [row for row in result['all_candidates'] if row['route']==route]
            x = [row['step'] for row in rows]
            y = [row['information_efficiency'] for row in rows]
            ax.plot(x,y,marker=MARKERS[route],markersize=4.,lw=1.15,
                    color=COLORS[route],label=route,zorder=4 if route!='F' else 2)
            lower = [row['information_efficiency_lower'] for row in rows]
            upper = [row['information_efficiency_upper'] for row in rows]
            ax.fill_between(x,lower,upper,color=COLORS[route],alpha=.18,lw=0)
        # 2. 圈出每轮实际选择，以观测记录生成横轴标签
        for action in result['trace']:
            ax.scatter(action['step']-1,action['information_efficiency'],s=64,facecolors='none',
                       edgecolors=DARK,linewidths=.9,zorder=8)
        labels = ['0']+[f'{row["step"]} ({row["route"]}{"+" if row["survey_observation"] else "−"})'
                        for row in result['trace']]
        ax.set(xlim=(-.16,len(labels)-.82),ylabel='单位勘察费信息价值\n(kW² / 勘察单位)')
        ax.set_xticks(range(len(labels)), labels)
        ax.yaxis.set_major_locator(MaxNLocator(5))
        ax.set_xlabel('决策时刻：已完成的勘察次数（括号为刚获得的结果）',fontsize=7.)
        ax.grid(axis='y',color='#E6EBEF',lw=.45,zorder=0)
        ax.legend(loc='upper right',ncol=6,columnspacing=.85,handlelength=1.3,
                  bbox_to_anchor=(1.005,1.25),fontsize=7.)


    with mpl.rc_context(STYLE):
        # 3. 总图：网架、非凸见证、逐轮区域及停止证书
        # 按实际状态数安排版面，读取本次预算与停止阈值
        rows = (len(result['states'])+2)//3
        height_mm = 135+65*rows
        protocol = result['protocol']
        fig = plt.figure(figsize=(WIDTH_MM/25.4,height_mm/25.4))
        fig.text(.055,.976,'有限预算下，勘察逐步识别可实现的非凸规划域',fontsize=11.,weight='bold',color=DARK)
        fig.text(.055,.952,f'五节点合成示例  |  建设预算 {protocol["budget"]:g}；按单位勘察费的信息价值选择道路',
                 fontsize=7.5,color=MUTED)
        grid = fig.add_gridspec(rows+2,3,height_ratios=[1.1,*[1.3]*rows,.9],
                                left=.08,right=.96,bottom=.09,top=.91,wspace=.50,hspace=.90)

        # 绘制网架和本次非凸见证
        net = fig.add_subplot(grid[0,0],label='topology')
        ax = net
        # 设置节点及候选道路的示意坐标
        positions = {0:(0,2.8),3:(-1,1.5),4:(1,1.5),1:(-1,0),2:(1,0)}
        paths = {
            'A':[(0,2.8),(-1.75,2.8),(-1.75,0),(-1,0)],
            'B':[(0,2.8),(1.75,2.8),(1.75,0),(1,0)],
            'C':[(-1,0),(-1,-.48),(1,-.48),(1,0)],
            'D':[(-1,0),(1,1.5)], 'E':[(1,0),(-1,1.5)],
            'F':[(-1,1.5),(1,1.5)]}
        radius = .19
        # 绘制既有线路、候选道路和节点
        for a,b in [(0,3),(3,1),(0,4),(4,2)]:
            points = np.asarray([positions[a],positions[b]],float)
            unit = (points[1]-points[0])/np.linalg.norm(points[1]-points[0])
            points[0] += radius*unit
            points[-1] -= radius*unit
            ax.plot(points[:,0],points[:,1],color=DARK,lw=1.1,zorder=1)
        for route,path in paths.items():
            points = np.asarray(path,float)
            points[0] += radius*(points[1]-points[0])/np.linalg.norm(points[1]-points[0])
            points[-1] += radius*(points[-2]-points[-1])/np.linalg.norm(points[-2]-points[-1])
            ax.plot(points[:,0],points[:,1],color=COLORS[route],lw=1.,ls=(0,(3,2)),zorder=2)
        for node,(x,y) in positions.items():
            ax.add_patch(Circle((x,y),radius,facecolor=DARK if node==0 else '#EDF3F7',
                                edgecolor=DARK,lw=.75,zorder=3))
            ax.text(x,y,str(node),ha='center',va='center',fontsize=7.,zorder=4,
                    color='white' if node==0 else DARK)
        # 标注道路、负荷与电源
        labels = {'A':(-1.96,1.2),'B':(1.96,1.2),'C':(0,-.74),
                  'D':(-.40,.18),'E':(.40,.18),'F':(0,1.77)}
        for road,(x,y) in labels.items():
            ax.text(x,y,road,ha='center',va='center',color=COLORS[road],weight='bold',fontsize=8.)
        ax.text(0,3.19,'电源',ha='center',fontsize=7.)
        ax.text(-1.43,-.34,r'$p_1$',fontsize=8.,ha='center')
        ax.text(1.43,-.34,r'$p_2$',fontsize=8.,ha='center')
        ax.set(xlim=(-2.25,2.25),ylim=(-.92,3.44),aspect='equal')
        ax.axis('off')
        net.set_title('a  网架与候选道路',loc='left',weight='bold')
        net.text(.5,-.08,'实线：既有；虚线：候选',ha='center',transform=net.transAxes,fontsize=6.7,color=MUTED)
        nonconvex = fig.add_subplot(grid[0,1],label='nonconvexity')
        ax = nonconvex
        # 叠加最终确认域及其凸包
        final = shape(result['states'][-1]['confirmed'])
        fill_geometry(ax,final.convex_hull,UNKNOWN)
        fill_geometry(ax,final,CONFIRMED)
        boundary(ax,final,DARK,linewidth=1.)
        # 标注两个可行端点与不可行中点
        witness = result['nonconvexity_witness']
        a,b,m = [np.asarray(witness[key]) for key in ('p_a','p_b','p_mid')]
        ax.plot([a[0],b[0]],[a[1],b[1]],color=MUTED,ls='--',lw=.8,zorder=5)
        ax.scatter([a[0],b[0]],[a[1],b[1]],s=15,color=DARK,zorder=6)
        ax.scatter([m[0]],[m[1]],s=22,color=RED,marker='x',linewidth=1.2,zorder=7)
        for p,label,offset in [(a,'U',(-10,-9)),(b,'V',(5,2)),(m,'M',(3,6))]:
            ax.text(*(p+offset),label,color=RED if label=='M' else DARK,fontsize=7.5)
        domain_axes(ax)
        nonconvex.set_title('b  方案域的并集可以非凸',loc='left',weight='bold')
        notes = fig.add_subplot(grid[0,2],label='legend')
        notes.axis('off')
        notes.text(0,1.,'U、V 分别可行，中点 M 不可行\n各负荷点可以选择不同建设方案',
                   va='top',fontsize=7.,color=DARK,linespacing=1.7)
        handles = [Patch(facecolor=UNKNOWN,label='尚未排除的可能范围'),
                   Patch(facecolor=BASE,label='初始调度域'),
                   Patch(facecolor=CONFIRMED,label='此前已确认'),
                   Patch(facecolor=GREEN_LIGHT,label='本轮新增 +'),
                   Patch(facecolor=RED_LIGHT,edgecolor=RED,hatch='///',label='本轮删除 −')]
        notes.legend(handles=handles,loc='lower left',borderaxespad=0,fontsize=6.8)

        # 逐轮叠加确认域、乐观域及本轮增加/删除部分
        initial_optimistic = shape(result['states'][0]['optimistic'])
        baseline = shape(result['states'][0]['confirmed'])
        for index,state in enumerate(result['states']):
            row,col = divmod(index,3)
            ax = fig.add_subplot(grid[row+1,col],label=f'state_{index}')
            current,optimistic = shape(state['confirmed']),shape(state['optimistic'])
            fill_geometry(ax,optimistic,UNKNOWN)
            fill_geometry(ax,current,CONFIRMED)
            fill_geometry(ax,baseline,BASE,zorder=2)
            if index:
                before = result['states'][index-1]
                added = current.difference(shape(before['confirmed']))
                removed = shape(before['optimistic']).difference(optimistic)
                fill_geometry(ax,added,GREEN_LIGHT,edgecolor=GREEN,zorder=3)
                fill_geometry(ax,removed,RED_LIGHT,edgecolor=RED,hatch='///',zorder=3)
            boundary(ax,initial_optimistic,'#B9C2CA',linestyle=':',linewidth=.6,zorder=4)
            boundary(ax,optimistic,DARK,linestyle='--',linewidth=.8,zorder=5)
            boundary(ax,current,'#426783',linewidth=1.,zorder=6)
            domain_axes(ax)
            if index==0:
                title = 'c  初始：尚未勘察'
                detail = f'初始调度域占乐观域 {result["baseline_share"]:.1%}'
                color = MUTED
            else:
                action = result['trace'][index-1]
                sign = '+' if action['survey_observation'] else '−'
                status = '可用' if action['survey_observation'] else '不可用'
                title = f'{chr(99+index)}  {index} · {action["route"]} {status}'
                area = action['realized_gain']+action['realized_removal']
                detail = f'{"新增" if action["survey_observation"] else "删除"} {sign}{area:,.1f} kW²'
                color = GREEN if action['survey_observation'] else RED
            ax.set_title(title,loc='left',fontsize=7.8,weight='bold',color=DARK,pad=8)
            ax.text(0,-.31,detail,transform=ax.transAxes,fontsize=6.8,color=color)
            ax.text(0,-.42,f'确认 {state["confirmed_area"]:,.0f} / 乐观 {state["optimistic_area"]:,.0f} kW²',
                    transform=ax.transAxes,fontsize=6.5,color=MUTED)

        # 展示条件价值与停止证书
        value = fig.add_subplot(grid[-1,:],label='value_history')
        position = value.get_position()
        value.set_position([.12,position.y0,.84,position.height])
        plot_value_history(value,result)
        value.set_title(f'{chr(99+len(result["states"]))}  动态条件价值：圆环为当轮选择',
                        loc='left',fontsize=8.,weight='bold',color=DARK,pad=18)
        uninspected = '、'.join(result['uninspected']) or '无'
        fig.text(.055,.020,f'停止时整体间隙上界 {result["states"][-1]["information_gap_upper"]:.2f} kW²；'
                 f'停止阈值 {protocol["stopping_area_tolerance"]:g} kW²；未勘察道路：{uninspected}。',
                 fontsize=6.8,color=MUTED)
        figures = [(fig, 'concept5_overview')]

        # 4. 条件价值图：完整曲线与逐轮评分矩阵
        # 绘制条件价值曲线
        fig = plt.figure(figsize=(WIDTH_MM/25.4,160/25.4))
        fig.text(.055,.96,'道路价值是条件量：比较同一决策时刻，而非不同道路的最后一次评分',
                 fontsize=9.8,weight='bold',color=DARK)
        fig.text(.055,.92,'a  每条路线一条轨迹；圆环 = 当轮选中；勘察后轨迹结束，不能补成零',fontsize=7.8,color=DARK)
        curve = fig.add_axes([.115,.575,.81,.275],label='conditional_value')
        plot_value_history(curve,result)
        fig.text(.055,.456,'b  完整评分矩阵：黑框 = 当轮选择；“—” = 已勘察，不再评分',fontsize=7.8,color=DARK)
        # 由全部候选评分构造道路 × 决策时刻矩阵
        routes = list(result['protocol']['routes'])
        steps = len(result['states'])
        matrix = np.full((len(routes),steps),np.nan)
        for row in result['all_candidates']:
            matrix[routes.index(row['route']),row['step']] = row['information_efficiency']
        heat = fig.add_axes([.115,.105,.81,.305],label='score_matrix')
        cmap = mpl.colormaps['Blues'].with_extremes(bad='#F5F6F7')
        maximum = np.nanmax(matrix)
        heat.imshow(matrix,aspect='auto',cmap=cmap,vmin=0,vmax=maximum,interpolation='none')
        for i in range(len(routes)):
            for j in range(steps):
                value = matrix[i,j]
                label = '—' if np.isnan(value) else f'{value:,.1f}'
                heat.text(j,i,label,ha='center',va='center',fontsize=7.4,
                          color='white' if value>.6*maximum else MUTED if np.isnan(value) else DARK)
        # 标记实际选择，标注最后一轮停止检查
        for row in result['trace']:
            i,j = routes.index(row['route']),row['step']-1
            heat.add_patch(Rectangle((j-.465,i-.435),.93,.87,fill=False,edgecolor=DARK,lw=1.1))
        heat.set_yticks(range(len(routes)),routes)
        heat.set_xticks(range(steps),[f'第 {i+1} 次前' for i in range(steps-1)]+['停止检查'])
        heat.tick_params(length=0,pad=6)
        for i,label in enumerate(heat.get_yticklabels()):
            label.set_color(COLORS[routes[i]])
            label.set_weight('bold')
        for spine in heat.spines.values():
            spine.set_visible(False)
        fig.text(.055,.035,'数值单位：kW² / 勘察单位。曲线阴影为求域数值误差区间；各轮仅比较当时尚未勘察的道路。',
                 fontsize=7.,color=MUTED)
        figures.append((fig, 'conditional_road_values'))

        # 5. 统一导出 PDF、SVG、PNG
        for fig, name in figures:
            for suffix in ('pdf', 'svg', 'png'):
                fig.savefig(output/f'{name}.{suffix}', dpi=300)
            plt.close(fig)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='从 results.json 绘制勘察结果')
    parser.add_argument('--results', type=Path, required=True)
    args = parser.parse_args()
    result = json.loads((args.results/'results.json').read_text(encoding='utf-8'))
    render_survey(result, args.results/'figures')
