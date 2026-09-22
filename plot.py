"""绘图与只读监视：终端进度、本地 HTML、完整回放和导出。"""
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


def pack_replay(value):
    """无损共享重复的几何与历史对象；只改变 HTML 编码，原始 JSON 不变。"""
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
        self.show_ui = show_ui
        self.record = show_ui if record is None else bool(record)
        self.output = Path(output) if output is not None else None
        self.open_browser = open_browser
        self.heartbeat_seconds = heartbeat_seconds
        self.stream = sys.stdout if stream is None else stream
        self.started = perf_counter()
        self.finished = None
        self.last_event = self.last_print = self.started
        self.lock = RLock()
        self.stop = Event()
        self.server = self.server_thread = self.heartbeat_thread = None
        self.url = None
        self.events, self.history = [], []
        self.previous = {}
        self.journal = None
        self.state = dict(status='running', event='start', message='准备启动',
                          query=0, iteration=0, total_cuts=0, pool_size=0,
                          method_number=0, method_count=4, revision=0, methods=[])

    def __enter__(self):
        try:
            if self.record and self.output is not None:
                self.output.mkdir(parents=True, exist_ok=True)
                self.journal = (self.output/'events.jsonl').open('w', encoding='utf-8')
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
                self.state.update(iteration=0, choice=None, latest_cut=None,
                                  bound=None, objective=None, eta=None, query_status=None)
                if data.get('mode') != 'SP-gap-support':
                    self.state.update(covered_by=None, covered_cost=None, uncovered_witness=None)
            if event == 'phase_start':
                self.state.update(query=0, point=None, lower=None, upper=None)
                self.state.update(geometry=[], global_outer=None, coverage_complete=False, region=None,
                                  states=None, counts=[], counts_algorithm={}, coverage_bound=None,
                                  max_total=None, max_total_bound=None, mode=None)
                if not self.record:
                    self.state.pop('states', None)
            if event == 'query_end':
                self.state['query_status'] = data.get('status')
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
                if self.journal is not None:
                    self.journal.write(json.dumps(frame, ensure_ascii=False, allow_nan=False)+'\n')
                    self.journal.flush()
            # 短查询合并到每秒一次的状态输出，避免快速求解时终端刷屏。
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
        """恢复完整历史，不调用优化器、不生成虚构求解事件。"""
        value = json.loads(Path(path).read_text(encoding='utf-8'))
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
        self.output.mkdir(parents=True, exist_ok=True)
        if not hasattr(self, 'template'):
            from plotly.offline import get_plotlyjs
            self.template = (ROOT/'live_view.html').read_text(encoding='utf-8')
            self.plotly = get_plotlyjs().encode('utf-8')
        self.state['viewer_sha256'] = sha256(self.template.encode('utf-8')).hexdigest()
        payload = self.snapshot().decode('utf-8')
        temporary = self.output/'replay.json.tmp'
        temporary.write_text(payload, encoding='utf-8')
        temporary.replace(self.output/'replay.json')
        # 大记录包含很多跨帧重复几何，先共享对象再压缩，降低浏览器峰值内存。
        if len(payload) > 10_000_000:
            payload = json.dumps(pack_replay(json.loads(payload)), ensure_ascii=False,
                                 allow_nan=False, separators=(',', ':'))
        data = base64.b64encode(gzip.compress(payload.encode('utf-8'), compresslevel=6, mtime=0)).decode('ascii')
        html = self.template.replace('<script src="/plotly.min.js"></script>',
                                     '<script>'+self.plotly.decode('utf-8')+'</script>')
        html = html.replace('/*SNAPSHOT*/', 'window.SAVED_REPLAY_GZIP="'+data+'";')
        path = self.output/'live_view.html'
        path.write_text(html, encoding='utf-8')
        print(f'完整过程回放（{len(self.history)} 条事件）：{path}', file=self.stream, flush=True)

    def close(self):
        self.stop.set()
        if self.journal is not None:
            self.journal.close()
            self.journal = None
        if self.server is not None:
            if self.server_thread is not None and self.server_thread.is_alive():
                self.server.shutdown()
                self.server_thread.join(timeout=2.)
            self.server.server_close()
        if self.heartbeat_thread is not None:
            self.heartbeat_thread.join(timeout=2.)


def print_summary(result):
    """终端显示同一原始结果推导出的 FR/MR；未知时显示区间。"""
    from vertify import METHODS, METHOD_NAMES
    names = dict(zip(METHODS, METHOD_NAMES))
    print('\n预算 | 方法 | FR (%) | MR (%) | 求解耗时 (s)', flush=True)
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
        print(f"{budget} | {names[row['method']]} | {' | '.join(rates)} | {seconds:.3f}", flush=True)


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
    import plotly.graph_objects as go
    from vertify import METHOD_NAMES  # 显示名称与结果数组的方法轴一致。

    fig = make_subplots(rows=2, cols=2, specs=[[{"type": "scene"}]*2]*2,  # 四个面板均为可旋转的三维场景。
                        subplot_titles=[f"{letter}  {name}" for letter, name in zip("abcd", METHOD_NAMES)],  # 使用简洁的 a、b、c、d 面板标记。
                        horizontal_spacing=.02, vertical_spacing=.08)  # 统一控制面板间距。
    categories = [(3, "与 AC 重合", "#4286AD", 1., "common"),  # 重合部分使用克制的蓝色。
                  (2, "遗漏", "#D43D3D", 1., "missed"),  # 遗漏的 AC 可行区域使用红色。
                  (1, "多余", "#E9B72F", 1., "extra")]  # 多余的计算区域使用黄色。
    categories.append((4, "未确定", "#A6A6A6", .22, "unknown"))
    for panel, labels in enumerate(result.labels[:, budget_index]):  # 读取当前预算下每种方法的差集标签。
        for code, name, color, opacity, group in categories:  # 每个面板分别生成重合、遗漏和多余表面。
            mesh = voxel_surface(labels == code, result.spacing)  # 使用同一套体素表面算法避免不同颜色口径不一致。
            points = np.asarray(mesh["vertices"]).reshape(-1, 3)  # 将输出顶点转为标准三列坐标数组。
            faces = np.asarray(mesh["triangles"], dtype=int).reshape(-1, 3)  # 将三角形转为三列整数索引数组。
            if len(points):  # 类别非空时添加实际表面。
                trace = go.Mesh3d(x=points[:, 0], y=points[:, 1], z=points[:, 2],  # 三维网格坐标使用实际 kW 值。
                                 i=faces[:, 0], j=faces[:, 1], k=faces[:, 2],  # 显式指定三角形，保留非凸区域及孔洞。
                                 name="AC 参考域" if panel == 2 and code==3 else name,  # AC 未确定部分仍标灰色。
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
        for name, groups in [("全部", {"common", "missed", "extra", "unknown"}),  # 默认同时显示差集和未确定区域。
                             ("仅差异", {"missed", "extra"}),  # 可只观察红黄差集。
                             ("计算域", {"common", "extra"}),  # 计算域由重合与多余部分组成。
                             ("AC 参考域", {"common", "missed"}),
                             ("未确定", {"unknown"})]]  # 灰色可独立查看或关闭。
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
    geometry_note = "自适应网格逐块认证；灰色为未确定，FR/MR 显示保守区间。网格完整不表示连续边界精确。方法耗时覆盖全部预算，在各预算页重复显示，不应相加。"
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
            def rate(key):
                if row[key] is not None:
                    return f'{row[key]:.4f}%'
                interval = row.get(key.replace('_percent','_interval'))
                return '—' if interval is None else f'[{interval[0]:.3f}, {interval[1]:.3f}]%'
            table.append(f'<tr><td>{name}</td><td>{rate("fr_percent")}</td>'  # 输出方法名称和 FR。
                         f'<td>{rate("mr_percent")}</td><td>{row["total_seconds"]:.3f}</td><td>{row.get("unknown_cells",0)}</td></tr>')  # 未确定点与物理误差分开列出。
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
<table><thead><tr><th>方法</th><th>FR</th><th>MR</th><th>整次方法耗时（秒）</th><th>未确定格点</th></tr></thead><tbody>{''.join(table)}</tbody></table>
<p>图：{result.metadata['network']} 的{scope}；节点负荷单位为 kW，各面板使用同一组刻度。<br>
网格步长 {spacing} kW；蓝色为重合部分，AC 面板显示独立参考；红色遗漏，黄色多余。<br>
FR = 多余 / 计算域；MR = 遗漏 / AC 域。AC 自比较的零仅表示它是基准；有限网格上的零不表示连续误差严格为零。<br>
总时间含建模、求解、构域及网格判定，混合方法计入线性阶段，AC 方法成本单列。绘图和导出不计。<br>
{geometry_note}</p>
</body></html>'''
        (folder/filename).write_text(html, encoding='utf-8')  # 写入可再生成的展示页面，原始结果仍只存一次。
    return folder/filenames[0]  # 返回默认预算页面作为 Notebook 入口。


def save_benchmark_report(output, records, labels, *, ac=None, case_root=None):
    """基准测试与 AC 复核共用报告；原始数据和完整过程仍保存在各案例目录。"""
    from html import escape
    from os.path import relpath
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    case_root = output if case_root is None else Path(case_root)
    checks = {} if ac is None else {(r['candidate_count'], r['budget'], r['variant']): r for r in ac['records']}
    number = lambda v, n=3: '—' if v is None else f'{v:.{n}f}'
    status_names = dict(certified='已完成覆盖认证', unknown='未确定', time_limit='达到时限', timeout='达到时限')
    rows, markdown = [], []
    for record in records:
        r = record['region']
        count, budget, variant = record['candidate_count'], record['budget'], record['variant']
        name = f'n{count}_b{"inf" if budget is None else f"{budget:g}"}_{variant}'
        link = Path(relpath(case_root/name/'live_view.html', output)).as_posix()
        check = checks.get((count, budget, variant), {})
        error = check.get('volume_error', r.get('validation', {}))
        phase = record.get('phases', r.get('phases', {}))
        mode = r.get('residual_mode', {'direct_all': 'physical', 'direct_mp': 'light', 'baseline': 'light'}.get(variant, '—'))
        error_text = number(error.get('region_error_percent'))
        if error.get('region_error_percent') is None and error.get('region_error_interval') is not None:
            error_text = '–'.join(number(v) for v in error['region_error_interval'])
        vertices = check.get('strict_vertices', {})
        quality = f"{vertices.get('feasible', 0)} / {vertices.get('total', 0)}" if vertices else '待 AC 校验'
        values = [count, '无限' if budget is None else f'{budget:g}', labels.get(variant, variant),
                  {'light':'轻量', 'physical':'物理'}.get(mode, mode), status_names.get(r['status'], r['status']),
                  number(record['solve_seconds']), number(phase.get('MP2')), number(phase.get('MP1')),
                  number(phase.get('residual')), r['counts']['sp'], r['counts'].get('sp_skipped', 0),
                  r['counts'].get('reused_certificates', 0), number(r.get('max_total')), error_text,
                  number(error.get('fr_percent')), number(error.get('mr_percent')), quality]
        rows.append('<tr>'+''.join('<td>'+escape(str(v))+'</td>' for v in values)+
                    '<td><a href="'+escape(link, quote=True)+'">完整回放</a></td></tr>')
        markdown.append('|'+ '|'.join(map(str, values))+'|[回放]('+link+')|')
    headers = ['候选线路', '预算', '算法', '剩余域', '完成状态', '构域 s', 'MP2 s', 'MP1 s',
               '剩余搜索 s', 'SP 次数', '覆盖跳过', '证书复用', '最大负荷 kW', '区域误差 %', '多算 %', '漏算 %', 'AC 顶点通过', '过程']
    total = len(records)
    completed = sum(r['status'] == 'certified' for r in records)
    note = f'{total} 组测试，{completed} 组完成全局覆盖认证。有限预算默认轻量搜索，无限预算默认物理搜索。'
    explanation = ('构域时间含几何与过程记录，不含评价箱准备、导出和独立校验。MP2、MP1、剩余搜索时间含各自建模与求解。'
                   '区域误差 = 多算与漏算体积之和 / 两域并集体积；多算率以计算域为分母，漏算率以 AC 域为分母。'
                   '比较对象是当前保留的内域并集；未完成组的误差仅代表当前结果，不能作为最终精度。')
    ac_note = ('尚未附加 AC 复核。' if ac is None else
               f"AC 采用 {ac['protocol']['divisions']}³ 个等体积网格中心独立估计，另逐一复核内域顶点和最大负荷点。"
               '网格误差不是连续体积误差的严格上界；AC 未确定点以误差区间保留。')
    if ac is not None:
        accepted = sum(c['strict_vertices']['feasible'] for c in ac['records'])
        checked = sum(c['strict_vertices']['total'] for c in ac['records'])
        ac_note += f' 严格 AC 顶点通过 {accepted} / {checked}。'
    doc = ('<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">'
           '<title>PlanRegion · 主线框架测试</title><style>body{font:15px/1.7 system-ui,"Microsoft YaHei",sans-serif;margin:32px;background:#f4f7fa;color:#20384b}'
           'h1{font-size:26px}.table{overflow:auto;background:white;border:1px solid #dce5eb;border-radius:12px}'
           'table{border-collapse:collapse;width:100%}th,td{padding:10px 12px;text-align:left;border-bottom:1px solid #e4ebf0;white-space:nowrap;font-size:13px}'
           'th{background:#eaf1f6}a{color:#08758c}p{max-width:1100px}.badge{color:#17795d}</style>'
           '<h1>完整 MP1 / MP2 · 可选择剩余域搜索</h1><p class="badge">'+escape(note)+'</p><p>'+escape(explanation)+'</p><p>'+escape(ac_note)+'</p>'
           '<div class="table"><table><thead><tr>'+''.join('<th>'+h+'</th>' for h in headers)+'</tr></thead><tbody>'+''.join(rows)+
           '</tbody></table></div><p><a href="summary.json">原始汇总数据</a></p></html>')
    (output/'report.html').write_text(doc, encoding='utf-8')
    (output/'report.md').write_text('\n'.join([note, '', explanation, '', ac_note, '',
        '|'+ '|'.join(headers)+'|', '|'+ '|'.join(['---']*len(headers))+'|', *markdown])+'\n', encoding='utf-8')
