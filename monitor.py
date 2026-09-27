"""运行记录、步骤控制与本地查看服务。"""
from pathlib import Path
from threading import Condition, Event, RLock, Thread
from time import perf_counter
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlsplit, parse_qs
from hashlib import sha256
import base64
import gzip
import json
import webbrowser

import numpy as np
from plot import ROOT, json_value, view_data, cut_slice, pack_replay, print_summary


class RunMonitor:
    """同一份记录用于落盘、实时展示和回放；步骤控制只在计算边界等待。"""

    def __init__(self, *, show_ui=False, output=None, open_browser=True,
                  heartbeat_seconds=10., stream=None, record=True, step_by_step=False):
        import sys
        # 1. 显示、记录和终端输出配置
        self.show_ui = show_ui
        self.record = record
        self.output = Path(output) if output is not None else None
        self.open_browser = open_browser
        self.heartbeat_seconds = heartbeat_seconds
        self.stream = sys.stdout if stream is None else stream
        # 2. 计时与实时服务
        self.started = perf_counter()
        self.finished = None
        self.last_event = self.last_print = self.started
        self.lock = RLock()
        self.condition = Condition(self.lock)
        self.step_by_step = step_by_step
        self.paused_seconds = 0.
        self.pause_started = None
        self.journal = None
        self.stop = Event()
        self.server = self.server_thread = self.heartbeat_thread = None
        self.url = None
        # 3. 当前状态、增量历史与展示几何缓存
        self.events, self.history = [], []
        self.geometry_cache = {}
        self.previous = {}
        self.state = dict(status='running', event='start', message='准备启动',
                          query=0, iteration=0, total_cuts=0, pool_size=0,
                          method_number=0, method_count=4, revision=0, methods=[], waiting=False,
                          run_id=self.output.name if self.output else None)

    def __enter__(self):
        try:
            if self.record and self.output is not None:
                self.output.mkdir(parents=True, exist_ok=True)
                self.journal = (self.output/'steps.jsonl').open('w', encoding='utf-8', buffering=1)
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

    def __call__(self, event, *, checkpoint=False, **data):
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
            step_seconds = now-self.paused_seconds-self.last_event
            self.last_event = now-self.paused_seconds
            self.state.update(event=event, message=data.get('message', event),
                              revision=self.state['revision']+1,
                              paused_seconds=self.paused_seconds, wall_seconds=now-self.started,
                              step_by_step=self.step_by_step,
                              waiting=checkpoint and self.step_by_step and not self.stop.is_set())
            self.pause_started = now if self.state['waiting'] else None
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
                                  bound=None, objective=None, eta=None, query_status=None, answer=None)
                if data.get('mode') != 'SP-gap-support':
                    self.state.update(covered_by=None, covered_cost=None, uncovered_witness=None)
            if event == 'phase_start':
                self.state.update(query=0, point=None, lower=None, upper=None, latest_cut=None)
                self.state.update(geometry=[], global_outer=None, coverage_complete=False, region=None,
                                  states=None, counts=[], counts_algorithm={}, coverage_bound=None,
                                  max_total=None, max_total_bound=None, mode=None, domain=None, point_reason=None)
                if not self.record:
                    self.state.pop('states', None)
            if event == 'query_end':
                self.state['query_status'] = data.get('status')
            if event in ('survey_state', 'completed'):
                self.state.update(point_reason=None, answer=None)
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
            # 4. 追加变化帧并立即落盘；每 100 帧保留完整状态
            if self.record:
                item = {key: self.state.get(key) for key in
                        ('revision', 'method', 'phase', 'query', 'iteration', 'event', 'message', 'point', 'query_status')}
                item['elapsed'] = now-self.started-self.paused_seconds
                if event == 'cut':
                    item['cut'] = self.state.get('latest_cut')
                self.events.append(item)
                self.state.update(elapsed=now-self.started-self.paused_seconds, step_seconds=step_seconds)
                patch = {k: v for k, v in self.state.items() if k not in self.previous or self.previous[k] != v}
                frame = dict(id=len(self.history), elapsed=self.state['elapsed'], patch=patch)
                self.history.append(frame)
                self.previous = dict(self.state)
                if self.journal is not None:
                    saved = dict(frame, recording_version=2)
                    if frame['id'] % 100 == 0:
                        saved['patch'] = dict(self.state)
                    self.journal.write(json.dumps(saved, ensure_ascii=False, allow_nan=False)+'\n')
            # 短查询合并为每秒一次的终端进度
            if event in ('start', 'preparation', 'method_start', 'method_end', 'phase_start', 'phase_end',
                         'saving', 'plotting', 'loaded', 'completed', 'failed', 'interrupted') \
                    or (event == 'query_end' and data.get('status') == 'unknown') \
                    or now-self.last_print >= 1.:
                self.print_status()
            self.condition.notify_all()
        # 5. 释放记录锁后等待用户；下一步只放行当前边界
        with self.condition:
            try:
                while self.state['waiting'] and not self.stop.is_set():
                    self.condition.wait()
            finally:
                if self.pause_started is not None:
                    self.paused_seconds += perf_counter()-self.pause_started
                    self.pause_started = None
            if self.stop.is_set() and checkpoint:
                raise KeyboardInterrupt()

    def clock(self):
        """返回扣除人工等待的时钟，供计算时限和耗时共用。"""
        with self.lock:
            return (self.pause_started if self.pause_started is not None else perf_counter())-self.paused_seconds

    def control(self, action, revision=None):
        with self.condition:
            if self.finished is not None:
                return False
            if action == 'next':
                if not self.state['waiting'] or revision != self.state['revision']:
                    return False
                self.state['waiting'] = False
            elif action == 'continue':
                self.step_by_step = False
                self.state['waiting'] = False
            elif action == 'pause':
                self.step_by_step = True
            else:
                raise ValueError('Unknown control action')
            self.condition.notify_all()
            return True

    def print_status(self, *, heartbeat=False):
        now = perf_counter()
        state = self.state
        prefix = f"[{self.clock()-self.started:8.1f}s]"
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
            details.append(f'此步骤已等待 {self.clock()-self.last_event:.1f}s')
        print(f"{prefix} {state['message']} | {' | '.join(details)}", file=self.stream, flush=True)
        self.last_print = now

    def heartbeat(self):
        while not self.stop.wait(self.heartbeat_seconds):
            with self.lock:
                if self.state['status'] == 'running' and not self.state['waiting'] and perf_counter()-self.last_print >= self.heartbeat_seconds:
                    self.print_status(heartbeat=True)

    def snapshot(self, after=0):
        with self.lock:
            now = perf_counter() if self.finished is None else self.finished
            paused = self.paused_seconds+(now-self.pause_started if self.pause_started is not None else 0.)
            value = dict(self.state, elapsed=now-self.started-paused,
                          step_seconds=self.state.get('step_seconds', 0.) if self.finished else now-paused-self.last_event,
                         events=self.events[-150:], history=self.history[max(0, after):],
                          history_total=len(self.history), recording_version=2,
                          paused_seconds=paused,
                          wall_seconds=now-self.started, step_by_step=self.step_by_step)
            return json.dumps(value, ensure_ascii=False, allow_nan=False).encode('utf-8')

    def load_recording(self, path):
        """从逐步日志或自包含 HTML 恢复历史，不调用优化器。"""
        # 1. 读取已提交的日志帧，或 HTML 中的压缩回放
        if Path(path).suffix == '.jsonl':
            self.history, self.state = [], {}
            for line in Path(path).read_text(encoding='utf-8').splitlines(keepends=True):
                if not line.endswith('\n'):
                    break  # 中断写入的最后半行没有提交。
                frame = json.loads(line)
                self.history.append(frame)
                self.state.update(frame['patch'])
            self.previous = dict(self.state)
            self.paused_seconds = self.state.get('paused_seconds', 0.)
            self.finished = perf_counter()
            self.started = self.finished-self.state.get('wall_seconds', self.state.get('elapsed', 0.))
            return
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
        self.paused_seconds = value.get('paused_seconds', 0.)
        self.finished = perf_counter()
        self.started = self.finished-self.state.get('wall_seconds', self.state.get('elapsed', 0.))
        self.previous = dict(self.state)

    def start_server(self):
        # 使用现有 Plotly 的本地脚本；页面离线可用，无 CDN 或 Node 服务。
        from plotly.offline import get_plotlyjs
        self.plotly = get_plotlyjs().encode('utf-8')
        self.template = (ROOT/'live_view.html').read_text(encoding='utf-8')
        owner = self

        class Handler(BaseHTTPRequestHandler):
            def do_POST(self):
                if self.path != '/control':
                    self.send_error(404)
                    return
                command = json.loads(self.rfile.read(int(self.headers['Content-Length'])))
                if command.get('action') not in ('next', 'continue', 'pause'):
                    self.send_error(400)
                    return
                accepted = owner.control(command['action'], command.get('revision'))
                self.send_response(204 if accepted else 409)
                self.end_headers()

            def do_GET(self):
                path = urlsplit(self.path).path
                if path == '/':
                    payload, mime = owner.template.encode('utf-8'), 'text/html; charset=utf-8'
                elif path == '/state':
                    try:
                        after = int(parse_qs(urlsplit(self.path).query).get('after', ['0'])[0])
                    except ValueError:
                        self.send_error(400, 'after must be an integer')
                        return
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
        # 3. 数据与绘图库嵌入同一 HTML，供离线回放
        html = self.template.replace('<script src="/plotly.min.js"></script>',
                                     '<script>'+self.plotly.decode('utf-8')+'</script>')
        html = html.replace('/*SNAPSHOT*/', 'window.SAVED_REPLAY_GZIP="'+data+'";')
        path = self.output/'live_view.html'
        path.write_text(html, encoding='utf-8')
        print(f'完整过程回放（{len(self.history)} 条事件）：{path}', file=self.stream, flush=True)

    def close(self):
        with self.condition:
            self.stop.set()
            self.state['waiting'] = False
            self.condition.notify_all()
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


    def show_result(self, result, *, export=True, keep_ui=False):
        if isinstance(result, dict):
            print(f"勘察完成：{result['seconds']:.3f}s，{result['query_count']} 个道路集合，{result['joint_cut_count']} 条共享割；结果：{self.output}", file=self.stream)
        else:
            print_summary(result)
        if self.record or export or self.show_ui:
            exporting = perf_counter()
            self.save_snapshot()
            print(f'回放导出耗时：{perf_counter()-exporting:.3f}s', file=self.stream, flush=True)
        if self.show_ui and keep_ui:
            try:
                input('计算已完成，网页可继续查看；按 Enter 或 Ctrl+C 关闭监视服务。\n')
            except (EOFError, KeyboardInterrupt):
                pass

