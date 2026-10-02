"""原生监视器：RunMonitor 记录事件，NativeWindow 重建帧并组织各面板（绘图基元在 plot.py）。

录制状态、显示范围和过程图层各自独立；扫描结果不进入构域时间轴。
符号分区在子进程中计算：子进程的 RunMonitor 经队列发送增量帧，主进程加分区前缀、乘符号后并入同一时间轴。
二维/三维共用外包络、过程标记及回放控制，只在几何绘制处区分维数。
"""
from pathlib import Path
from threading import Condition, Event, Thread
from time import perf_counter
import gzip
import json
import multiprocessing

import numpy as np

from plot import (CUT, GLOBAL, INNER, NETWORK, OUTER, REFERENCE, SP, STEP_STYLE, cap, cut_polygon, cut_segment, draw_3d,
                  draw_geometry, voxel_faces)

COMPARISONS = {'result_ac': '实验结果 — AC', 'result_socp': '实验结果 — SOCP', 'socp_ac': 'SOCP — AC'}
MERGED = ('schemes', 'cones', 'cut_history')   # 按键增量合并的状态；锥被细分时其键置 None
FRAME_STRIDE = 64   # 回放状态快照间隔（帧）：任一帧的状态从最近的快照向后合并
_CHANNEL = None   # 子进程：(事件队列, 暂停, 单步, 取消)，由进程池初始化函数 _connect 设置


def _connect(channel):
    """进程池初始化：子进程的 RunMonitor 经此通道发送事件，并响应窗口的暂停、单步与取消。"""
    global _CHANNEL
    _CHANNEL = channel


def comparison_metrics(computed, reference):
    missed, extra = np.count_nonzero(reference & ~computed), np.count_nonzero(computed & ~reference)
    return dict(mr_percent=100.*missed/reference.sum() if reference.any() else None,
                fr_percent=100.*extra/computed.sum() if computed.any() else None,
                missed_cells=int(missed), extra_cells=int(extra),
                reference_cells=int(reference.sum()), computed_cells=int(computed.sum()))


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
    """网架、锥与割按键独立更新，历史帧只存发生变化的键。"""
    for key, value in patch.items():
        if key in MERGED:
            state[key] = {**state.get(key, {}), **value}
        else:
            state[key] = value


def signed_values(value, sign, prefix, key=''):
    """转发时把幅值坐标、割系数、网架标签转成带符号坐标与带分区前缀的标签。"""
    if value is None:
        return None
    if key in ('scheme', 'active_scheme'):
        return prefix+value
    if key == 'p':
        return (np.asarray(value)*sign).tolist()
    if key in ('vertices', 'inner', 'outer'):
        return (np.asarray(value, float).reshape(-1, len(sign))*sign).tolist()
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
    """求解线程只提交数值；Tk 与绘图仅在主线程执行。sign 不为空时是子进程中的分区监视器。"""

    def __init__(self, *, output=None, clock=perf_counter, algorithm='主线', sign=None):
        self.output = None if output is None else Path(output)
        self._clock = clock
        self.algorithm = algorithm
        self.label = None if sign is None else ''.join('+' if s > 0 else '-' for s in sign)
        self.channel = None if sign is None else _CHANNEL
        self.shared = None            # 主进程：与子进程共用的 (队列, 暂停, 单步, 取消)
        self.condition = Condition()
        self.cancelled = Event()
        self.paused = False
        self.permits = 0
        self.paused_seconds = 0.
        self.started = self._clock()
        self.time_limit = np.inf
        self.state, self.history = {}, []
        self.snapshots = [{}]         # snapshots[k]：合并前 k·FRAME_STRIDE 帧后的状态
        self.validation_state = {}
        self.result = self.error = None
        self.busy = False

    def clock(self):
        return self._clock()-self.paused_seconds

    def timing(self):
        return dict(total_seconds=self.clock()-self.started)

    def _emit(self, event, *, checkpoint=True, **values):
        """1. 只保留相对当前状态有变化的键；2. 主进程追加到时间轴，子进程送入队列；
        3. 暂停时在检查点等待，单步放行一个事件，取消时中断计算。"""
        with self.condition:
            if self.cancelled.is_set():
                raise KeyboardInterrupt()
            values = _plain(values)
            patch = {k: v for k, v in values.items() if self.state.get(k) != v}
            for key in MERGED:
                if key in patch:
                    patch[key] = {k: v for k, v in patch[key].items() if self.state.get(key, {}).get(k) != v}
                    if not patch[key]:
                        del patch[key]
            patch['event'] = event
            _merge(self.state, patch)
            before = self._clock()
            if self.channel is None:
                self.history.append(dict(elapsed=self.clock()-self.started, patch=patch))
                while checkpoint and self.paused and not self.permits and not self.cancelled.is_set():
                    self.condition.wait()
                if self.permits:
                    self.permits -= 1
            else:
                queue, pause, step, cancel = self.channel[:4]
                queue.put((self.label, patch))
                while checkpoint and pause.is_set() and not cancel.is_set() and not step.acquire(timeout=.1):
                    pass
                if cancel.is_set():
                    self.cancelled.set()
            self.paused_seconds += self._clock()-before
            if self.cancelled.is_set():
                raise KeyboardInterrupt()

    def share(self):
        """主进程：建立与分区子进程共用的事件队列、暂停 / 单步 / 取消信号与仍在计算的分区数。"""
        context = multiprocessing.get_context('spawn')
        self.shared = context.Queue(), context.Event(), context.Semaphore(0), context.Event(), context.RawValue('i', 0)
        return self.shared

    def running(self):
        """仍在计算的分区数：子进程读主进程维护的共享计数；主进程内直接构域时为 1。"""
        return 1 if self.channel is None else self.channel[4].value

    def close(self):
        """子进程：分区结束标记。"""
        self.channel[0].put((self.label, None))

    def forward(self, label, patch):
        """主进程：子进程分区的增量帧加分区前缀、坐标乘符号后并入总时间轴；网架行与割行附分区符号 sign。"""
        sign = np.array([1 if s == '+' else -1 for s in label])
        prefix = label+':'
        mapped = signed_values(patch, sign, prefix)
        for key in MERGED:
            if key in mapped:
                mapped[key] = {prefix+k: v if v is None or key == 'cones' else {**v, 'sign': sign}
                               for k, v in mapped[key].items()}
        event = mapped.pop('event')
        self._emit(event, partition=label, **mapped)

    # 构域完成后独立扫描；最终指标、几何、回放只写同一个压缩文件。
    def scanning(self, completed, total):
        self._validation_update(phase='AC / SOCP 扫描', scan_progress=[completed, total])

    def _validation_update(self, **values):
        """校验只更新独立面板，不占用构域回放帧，也不等待单步按钮。"""
        with self.condition:
            if self.cancelled.is_set():
                raise KeyboardInterrupt()
            values = _plain(values)
            self.validation_state.update(values)
            self.state.update(values)

    def validation(self, reference, result, *, region_key='inner'):
        """逐格 FR/MR：配对参考的 AC 与 SOCP 标签 1/-1/0 为可行/已证不可行/未决；未决格不计入对应参考，另计格数。"""
        from region import covered
        # 1. 格心与两种参考
        states, socp = np.asarray(reference['states']), np.asarray(reference['socp_states'])
        indices = np.indices(states.shape).reshape(states.ndim, -1).T
        lower = np.asarray(reference.get('axis_lower', np.zeros(states.ndim)))
        points = lower+(indices+.5)*(np.asarray(reference['bounds'])-lower)/np.array(states.shape)
        known, truth = states.ravel() != 0, states.ravel() == 1
        # 2. 内域（主指标）与外包络相对 AC
        metrics, masks = {}, {}
        for key in (key for key in ('inner', 'outer') if key in result):
            inside = covered(points, result[key])
            metrics[key] = comparison_metrics(inside[known], truth[known])
            masks[key] = inside.reshape(states.shape)
        # 3. 三组对比：计算域—AC、计算域—SOCP、SOCP—AC（区分误差来自松弛还是模型）
        both = (socp != 0) & (states != 0)
        comparisons = dict(result_ac=metrics[region_key],
                           result_socp=comparison_metrics(masks[region_key][socp != 0], socp[socp != 0] == 1),
                           socp_ac=comparison_metrics((socp == 1)[both], (states == 1)[both]))
        validation = dict(axis_lower=lower, bounds=reference['bounds'], states=states, socp_states=socp,
                          computed_states=masks[region_key], region_key=region_key,
                          scan_seconds=reference.get('scan_seconds'),
                          method=reference.get('method', 'legacy_reference'),
                          ac_identity=reference.get('metadata', {}).get('identity'),
                          cache_path=reference.get('cache_path'),
                          reused_points=reference.get('reused_points'), computed_points=reference.get('computed_points'),
                          undecided_cells=int(np.count_nonzero(~known)),
                          socp_undecided_cells=int(np.count_nonzero(socp == 0)),
                          metrics=metrics, comparisons=comparisons, **metrics[region_key])
        self._validation_update(phase='完成', status='completed', validation=validation)

    def frame(self, index):
        """第 index 帧的状态：从最近的快照向后合并，途经步长的整数倍时补存快照（_merge 不原地修改，浅拷贝即可）。"""
        with self.condition:
            k = min((index+1)//FRAME_STRIDE, len(self.snapshots)-1)
            state = dict(self.snapshots[k])
            for position in range(k*FRAME_STRIDE, index+1):
                _merge(state, self.history[position]['patch'])
                if (position+1) % FRAME_STRIDE == 0 and (position+1)//FRAME_STRIDE == len(self.snapshots):
                    self.snapshots.append(dict(state))
            return state

    def control(self, action):
        """窗口按钮：暂停 / 继续 / 单步 / 取消；分区子进程经共用信号同步响应。"""
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
        if self.shared is not None:
            _, pause, step, cancel = self.shared[:4]
            if action in ('pause', 'next'):
                pause.set()
            if action == 'next':
                step.release()
            if action == 'continue':
                pause.clear()
            if action == 'cancel':
                cancel.set()

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
        if data['version'] != 4:
            raise ValueError('只支持 version=4 回放')
        self.history, self.validation_state, self.snapshots = data['history'], data['validation_state'], [{}]
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


STEP_TRAIL, STEP_LOG = 6, 8   # 主图保留的最近步骤数；步骤栏行数
FIT_SHRINK = .6   # 内容某一维缩到主图范围的这一比例以下时才收紧坐标
ALL_PARTITIONS = '全部分区'


def accepted(row):
    """网架行的 N^CUT_x 是否计入 I：状态在 region.CUT_ACCEPTED 中且非空。"""
    from region import CUT_ACCEPTED
    return row['status'] in CUT_ACCEPTED and len(row['outer']) > 0


class NativeWindow:
    """A 全局总图与步骤栏、B 动态网架、C 独立扫描；实时 / 回放共用同一绘图入口。"""

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
        self.scheme_views, self.view_signatures = {}, {}   # 网架面板；各面板上次绘制的内容签名
        self.limits = None                                  # 主图坐标 (下界, 上界)，随内容滞回调整
        self.frame_meta, self.partition_at = [], None       # 各帧 (分区, 事件, 是否步骤)；沿历史累积的当前分区
        self.log_rows = []                                  # 步骤栏各行对应的帧号
        self.scheme_page, self.focus_scheme = 0, None
        self.status, self.step_text = tk.StringVar(value='初始化'), tk.StringVar(value='步骤：—')
        ttk.Label(self.root, textvariable=self.status, font=('Microsoft YaHei', 11)).pack(anchor='w', padx=12, pady=(8, 2))
        ttk.Label(self.root, textvariable=self.step_text, font=('Microsoft YaHei', 10, 'bold'), wraplength=1380,
                  justify='left').pack(anchor='w', padx=12)
        options = ttk.Frame(self.root)
        options.pack(fill='x', padx=12)
        self.follow, self.full_box = tk.BooleanVar(value=False), tk.BooleanVar(value=False)
        ttk.Checkbutton(options, text='网架页跟随当前网架', variable=self.follow).pack(side='right')
        ttk.Checkbutton(options, text='坐标取整个分区盒', variable=self.full_box,
                        command=self.refit).pack(side='right', padx=12)
        self.chosen = tk.StringVar(value=ALL_PARTITIONS)
        self.partition_selector = ttk.Combobox(options, textvariable=self.chosen, values=[ALL_PARTITIONS],
                                               state='readonly', width=9)
        self.partition_selector.pack(side='right')
        self.partition_selector.bind('<<ComboboxSelected>>', lambda event: self.redraw())
        ttk.Label(options, text='逐步跟踪').pack(side='right', padx=4)
        body = ttk.Panedwindow(self.root, orient='horizontal')
        body.pack(fill='both', expand=True, padx=8, pady=6)
        left, right = ttk.Frame(body), ttk.Notebook(body)
        body.add(left, weight=3)
        body.add(right, weight=2)
        main = ttk.LabelFrame(left, text='A  全局总图 · 全部分区的外界 K^OUT、结果内域 I 与当前步骤（坐标固定）')
        main.pack(fill='both', expand=True, pady=2)
        steps = ttk.LabelFrame(left, text='步骤 · 所选分区最近的求解（高亮为当前；点击跳转）')
        steps.pack(fill='x', pady=2)
        self.log = tk.Listbox(steps, height=STEP_LOG, activestyle='none', exportselection=False,
                              font=('Microsoft YaHei', 9))
        self.log.pack(fill='x', padx=4, pady=3)
        self.log.bind('<<ListboxSelect>>', self.choose_step)
        networks, comparison = ttk.Frame(right), ttk.Frame(right)
        right.add(networks, text='B  网架 · N^CUT_x 灰色（停滞后计入 I，绿色）/ 紫色割')
        right.add(comparison, text='C  实验 / SOCP / AC 对比')
        self.comparison_mode = tk.StringVar(value=next(iter(COMPARISONS.values())))
        selector = ttk.Combobox(comparison, textvariable=self.comparison_mode,
                                values=list(COMPARISONS.values()), state='readonly', width=24)
        selector.pack(anchor='w', padx=8, pady=3)
        selector.bind('<<ComboboxSelected>>', lambda event: self._refresh_validation())
        self.comparison_text = tk.StringVar()
        ttk.Label(comparison, textvariable=self.comparison_text, justify='left').pack(anchor='w', padx=8)
        self.axes, self.canvases = {}, {}
        for panel, frame, size, margins in (('A', main, (7.6, 5.6), dict(left=.09, right=.98, bottom=.09, top=.92)),
                                            ('C', comparison, (5.8, 4.4), dict(left=.13, right=.97, bottom=.14, top=.90))):
            fig = Figure(figsize=size, dpi=100)
            ax = fig.add_subplot(111)
            fig.subplots_adjust(**margins)
            canvas = FigureCanvasTkAgg(fig, master=frame)
            canvas.get_tk_widget().pack(fill='both', expand=True)
            self.axes[panel], self.canvases[panel] = ax, canvas
        pages = ttk.Frame(networks)
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
        viewport = ttk.Frame(networks)
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
                             ('上一步', lambda: self.step(-1)),
                             ('下一步', lambda: self.step(1)),
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
        """上一步 / 下一步：选定分区时停在该分区相邻的步骤帧，否则逐帧（双窗口同步时总是逐帧）。"""
        if self.controller is not None:
            self.seek(self.controller.index+direction)
        else:
            self.seek(self._neighbor(direction))

    def seek_cut(self, direction):
        index = self._neighbor(direction, 'cut')
        if index != self.index:
            self.seek(index)

    def _neighbor(self, direction, event=None):
        """相邻帧号：未选分区且不限事件时逐帧；否则取所选分区内相邻的步骤帧（或 event 帧），没有则不动。"""
        self._scan_frames()
        chosen = self.chosen.get()
        if event is None and chosen == ALL_PARTITIONS:
            return self.index+direction
        indices = [i for i, (partition, kind, step) in enumerate(self.frame_meta)
                   if (step if event is None else kind == event) and chosen in (ALL_PARTITIONS, partition)
                   and direction*(i-self.index) > 0]
        return (min(indices) if direction > 0 else max(indices)) if indices else self.index

    def _scan_frames(self):
        """逐帧登记 (分区, 事件, 是否步骤)：partition 键只在变化时写入增量帧，沿历史累积得到每帧所属分区。"""
        for item in self.monitor.history[len(self.frame_meta):]:
            patch = item['patch']
            self.partition_at = patch.get('partition', self.partition_at)
            self.frame_meta.append((self.partition_at, patch.get('event'), bool(patch.get('step'))))

    def _steps(self):
        """所选分区（或全部分区）的步骤帧号，升序。"""
        chosen = self.chosen.get()
        return [i for i, (partition, _, step) in enumerate(self.frame_meta)
                if step and chosen in (ALL_PARTITIONS, partition)]

    def choose_step(self, event):
        selection = self.log.curselection()
        if selection:
            self.seek(self.log_rows[selection[0]])

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
            elif self.playing:   # 播放按所选分区逐步前进，到末帧或该分区再无步骤时停
                following = min(self._neighbor(1), total-1)
                self.playing = following not in (self.index, total-1)
                self.index = following
            self.show()
        self.root.after(30 if self.playing else 80, self.tick)   # 播放时只受绘制速度限制

    def redraw(self):
        self.last_drawn = -1
        self.show()

    def refit(self):
        self.limits = None
        self.redraw()

    @staticmethod
    def _boxes(state):
        """各分区盒（带符号 kW）：尚无锥的分区以盒为外界。"""
        from itertools import product
        corners = np.array(list(product((0., 1.), repeat=len(state['bounds']))))*np.asarray(state['bounds'])
        return [corners*np.array([1 if s == '+' else -1 for s in label]) for label in state.get('partitions', [])]

    def _outer_polygons(self, state):
        """外界 K^OUT：各锥外块 K^OUT_k；尚无锥的分区取分区盒。"""
        cones = {key: row for key, row in state.get('cones', {}).items() if row}
        started = {key.split(':')[0] for key in cones}
        return ([row['outer'] for row in cones.values()]
                +[box for label, box in zip(state.get('partitions', []), self._boxes(state)) if label not in started])

    @staticmethod
    def _inner_polygons(state):
        """结果内域 I 的组成：各锥的 K^IN_k 与计入 I 的 N^CUT_x。"""
        return ([row['inner'] for row in state.get('cones', {}).values() if row]
                +[row['outer'] for row in state.get('schemes', {}).values() if accepted(row)])

    def _fit_limits(self):
        """主图坐标：最新帧全部 K^OUT_k（尚无锥的分区取分区盒）的范围，勾选时取全部分区盒；回放时最新帧即最终
        外界，坐标全程不变。实时运行中只在内容超出范围、或某一维缩到范围的 FIT_SHRINK 以下时重设（两侧各留 10%），
        重设后网架面板随之重画。K^OUT ⊇ R^SOCP，射线顶点与 MISOCP 解点都落在最终范围内。"""
        with self.monitor.condition:
            latest = dict(self.monitor.state)
        d = len(latest['bounds'])
        polygons = self._boxes(latest) if self.full_box.get() else self._outer_polygons(latest)
        points = np.vstack([np.asarray(p, float).reshape(-1, d) for p in polygons if len(p)])
        lower, upper = points.min(axis=0), points.max(axis=0)
        if self.limits is not None:
            low, high = self.limits
            if np.all(lower >= low) and np.all(upper <= high) and np.all(upper-lower >= FIT_SHRINK*(high-low)):
                return
        padding = .1*np.maximum(upper-lower, 1.)
        self.limits, self.view_signatures = (lower-padding, upper+padding), {}

    def _view_limits(self, state, *, scheme=None, validation=False):
        """主图取稳定范围 self.limits，网架面板取其中本分区所在卦限的部分（原点一侧留 4%）；校验图取扫描与结果的范围。"""
        d = len(state['bounds'])
        if not validation:
            lower, upper = self.limits
            if scheme is None:
                return lower, upper
            sign, margin = np.array([1 if s == '+' else -1 for s in scheme.split(':')[0]]), .04*(upper-lower)
            return np.where(sign > 0, -margin, lower), np.where(sign < 0, margin, upper)
        groups = ([row['vertices'] for key in ('inner', 'outer') for row in state.get('result', {}).get(key, [])]
                  or self._boxes(state))   # 尚无结果时取全部分区盒
        reference = state.get('validation')
        if reference:
            truth = np.asarray(reference['states']) == 1
            if 'socp_states' in reference:
                truth |= np.asarray(reference['socp_states']) == 1
            cells = np.argwhere(truth)
            if len(cells):
                lower = np.asarray(reference.get('axis_lower', np.zeros(d)))
                widths = (np.asarray(reference['bounds'])-lower)/np.shape(reference['states'])
                groups.append(lower+np.array([cells.min(axis=0), cells.max(axis=0)+1])*widths)
        points = np.vstack([np.asarray(group).reshape(-1, d) for group in groups if len(group)])
        lower, upper = points.min(axis=0), points.max(axis=0)
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

    @staticmethod
    def _markers(ax, state, scheme):
        """网架面板的过程标记：当前步骤的点属于本网架时按步骤样式画出（射线顶点、MISOCP 解点、SP 评分点、取割点等）。"""
        step = state.get('step')
        if step and step.get('scheme') == scheme and 'p' in step:
            marker, color, _ = STEP_STYLE[step['kind']]
            ax.scatter(*np.asarray(step['p']).reshape(-1, len(state['bounds'])).T, color=color, marker=marker, s=40,
                       zorder=10, gid='step-point')

    def show(self):
        if not self.monitor.history:
            return
        self._refresh_validation()
        if self.index == self.last_drawn:
            return
        self.last_drawn = self.index
        self._scan_frames()
        state = self.monitor.frame(self.index)
        focus = state.get('active_scheme')
        if self.follow.get() and focus and focus != self.focus_scheme:   # 默认不跟随：新网架追加在末页，页面不跳
            self.focus_scheme = focus
            keys = list(state.get('schemes', {}))
            if focus in keys:
                self.scheme_page = keys.index(focus)//4
        self._show_header(state)
        if 'bounds' in state:
            self._fit_limits()
            self._draw_union(state)
            self._draw_schemes(state)

    def _refresh_validation(self):
        """校验独立于历史帧；停止拖动时也刷新，回退后仍保留最终对比。"""
        with self.monitor.condition:
            latest = dict(self.monitor.state)
        if 'bounds' not in latest:
            return
        key = (tuple(latest.get('scan_progress', [])),
               id(latest.get('result')), id(latest.get('validation')), latest.get('error'), self.comparison_mode.get())
        if key != self.last_validation:
            self.last_validation = key
            self._draw_validation(latest)

    def _show_header(self, state):
        d = len(state.get('load_nodes', (0, 1)))
        self.root.title(f"{state.get('algorithm', '主线')} · {state.get('network', '')} · {d} 维过程回放")
        total = len(self.monitor.history)
        self.frame_text.set(f'{self.index+1} / {total}')
        self.setting_slider = True
        self.slider.configure(to=max(1, total-1 if self.controller is None else self.controller.total-1))
        self.slider.set(self.index if self.controller is None else self.controller.index)
        self.setting_slider = False
        self.partition_selector.configure(values=[ALL_PARTITIONS, *state.get('partitions', [])])
        self._show_steps()
        phase, detail = state.get('phase', '初始化'), ''
        if phase.startswith('径向') and state.get('volume_ratio') is not None:
            detail = f" · ΣΔ/Σvol(K^IN) = {state['volume_ratio']:.4g}"
        elif phase == '网架切割' and state.get('active_scheme'):
            detail = f" · 网架 {state['active_scheme']}"
        if state.get('error'):
            detail = ' · '+state['error']
        partition = state.get('partition')
        cones = sum(row is not None for row in state.get('cones', {}).values())
        self.status.set(f"{'分区 '+partition+' · ' if partition else ''}{phase}{detail}    "
                        f"叶锥 {cones} · 网架 {len(state.get('schemes', {}))} · 割 {len(state.get('cut_history', {}))}")

    def _show_steps(self):
        """当前步骤（所选分区中不晚于当前帧的最后一步）的说明，步骤栏列出它之前 STEP_LOG-3 步与之后 2 步。"""
        from bisect import bisect_right
        steps = self._steps()
        position = bisect_right(steps, self.index)-1
        self.log_rows = steps[max(0, position-STEP_LOG+3):position+3]
        self.log.delete(0, 'end')
        for index in self.log_rows:
            item = self.monitor.history[index]
            self.log.insert('end', f"{item['elapsed']:8.2f} s   {self.frame_meta[index][0] or '':4s}{item['patch']['step']['text']}")
        if position < 0:
            self.step_text.set('步骤：—')
            return
        current = steps[position]
        row = self.log_rows.index(current)
        self.log.selection_set(row)
        self.log.see(row)
        item = self.monitor.history[current]
        step = item['patch']['step']
        point = '' if 'p' not in step else '    p=('+', '.join(f'{value:.1f}' for value in step['p'])+') kW'
        self.step_text.set(f"步骤 {position+1}/{len(steps)} · 分区 {self.frame_meta[current][0]} · "
                           f"{item['elapsed']:.2f} s · {step['text']}{point}")

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
            polygon = cut_polygon(cuts[latest]['cut'], state['schemes'][scheme]['x'], limits, lower)
            draw_3d(ax, [polygon], state['bounds'], color=CUT, fill=True, alpha=.10,
                     linewidth=1., gid=f'cut-{latest}-{scheme}')
            ax.text2D(.02, .97, f'割 #{latest} · 网架 {scheme}', transform=ax.transAxes,
                      color=CUT, fontsize=8, va='top')
            return
        previous = []
        for number, item in cuts.items():
            if not history and number != latest:
                continue
            segment = cut_segment(item['cut'], state['schemes'][scheme]['x'], limits, lower)
            if len(segment) != 2:
                if number == latest:
                    cut = np.asarray(item['cut'])
                    corners = np.array([lower, [limits[0], lower[1]], [lower[0], limits[1]], limits])
                    slack = cut[0]+corners@cut[1:3]+cut[3:]@state['schemes'][scheme]['x']
                    message = ('显示范围全部排除' if slack.max() < 0. else
                               '割在显示范围外' if slack.min() >= 0. else '仅与角点相交')
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
        from region import polygon_union
        if state.get('event') != 'cut':
            return polygon_union([])
        before = self.monitor.frame(self.index-1).get('schemes', {})
        after = state['schemes']
        keys = before.keys() if scheme is None else [scheme]
        old = polygon_union([before[key]['outer'] for key in keys if key in before])
        new = polygon_union([after[key]['outer'] for key in keys if key in after])
        return old.difference(new)

    def _draw_polygons(self, ax, polygons, state, *, gid=None, **style):
        """面板共用几何入口；只在此处分派二维并集和三维凸域绘制。"""
        from region import polygon_union
        if len(state['bounds']) == 3:
            draw_3d(ax, polygons, state['bounds'], gid=gid, **style)
        else:
            before = set(ax.get_children())
            draw_geometry(ax, polygon_union(polygons), **style)
            for artist in set(ax.get_children())-before:
                artist.set_gid(gid)

    def _draw_removed(self, ax, state, scheme=None):
        if len(state['bounds']) == 3:
            self._draw_polygons(ax, self._removed_3d(state, scheme), state,
                                color=CUT, fill=True, alpha=.18, gid='removed-region')
        else:
            draw_geometry(ax, self._removed(state, scheme), color=CUT, fill=True, alpha=.22)

    def _regions(self, ax, state):
        """全部分区的外界 K^OUT 与结果内域 I：二维为并集（内域再与外包络取交），三维为各锥远端面片与网架内域。"""
        from region import polygon_union
        if len(state['bounds']) == 2:
            hull = polygon_union(self._outer_polygons(state))
            draw_geometry(ax, hull, color=OUTER, fill=True, alpha=.10)
            draw_geometry(ax, hull, color=OUTER, linestyle='--', linewidth=1.3)
            draw_geometry(ax, polygon_union(self._inner_polygons(state)).intersection(hull), color=INNER, fill=True, alpha=.38,
                  linewidth=.8)
            return
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
        cones = [row for row in state.get('cones', {}).values() if row]
        for key, color, alpha in (('outer', OUTER, .12), ('inner', INNER, .30)):
            caps = [cap(row[key]) for row in cones]
            if caps:
                ax.add_collection3d(Poly3DCollection(caps, facecolors=color, edgecolors=color, linewidths=.2,
                                                     alpha=alpha, gid=f'cone-{key}'))
        sets = [row['outer'] for row in state.get('schemes', {}).values() if accepted(row)]
        draw_3d(ax, sets, state['bounds'], color=INNER, fill=True, alpha=.10, linewidth=.4, gid='network-inner')

    def _draw_union(self, state):
        """A 全局总图：全部分区的 K^OUT 与 I、最新割及其切掉的部分、步骤图层；坐标为稳定范围。"""
        from matplotlib.lines import Line2D
        ax = self.axes['A'] = self._axes(self.axes['A'], state)
        if len(state['bounds']) == 3:
            ax.set_position([.01, .03, .96, .86])
        self._regions(ax, state)
        self._draw_removed(ax, state)
        cuts = state.get('cut_history', {})
        if cuts:
            self._draw_cuts(ax, state, next(reversed(cuts.values()))['scheme'], history=False)
        self._draw_steps(ax, state)
        handles = [Line2D([], [], color=OUTER, ls='--', label='外界 K^OUT'), Line2D([], [], color=INNER, label='结果内域 I'),
                   Line2D([], [], color=INNER, lw=2.4, label='本步的锥 K^IN_k / K^OUT_k'),
                   Line2D([], [], color=NETWORK, ls=':', label='切割中的 N^CUT_x')]
        ax.legend(handles=handles, loc='lower left', bbox_to_anchor=(0., 1.01), ncol=4, frameon=False, fontsize=8,
                  columnspacing=.9, handlelength=1.8)
        self.canvases['A'].draw_idle()

    def _draw_steps(self, ax, state):
        """步骤图层：所选分区最近 STEP_TRAIL 个步骤的点与线（越早越淡），当前步骤的点加短标签，落在坐标外的点
        贴边画空心；当前帧新建或求界的锥描出 K^IN_k（粗实线）与 K^OUT_k（虚线），切割中的网架描出 N^CUT_x（点线）。"""
        from bisect import bisect_right
        d, (lower, upper) = len(state['bounds']), self.limits
        # 1. 当前帧的锥与切割中的网架
        patch = self.monitor.history[self.index]['patch']
        for row in (row for row in patch.get('cones', {}).values() if row):
            self._draw_polygons(ax, [row['inner']], state, color=INNER, linewidth=2.4, gid='step-cone')
            self._draw_polygons(ax, [row['outer']], state, color='#5d7686', linestyle='--', linewidth=1.4, gid='step-cone')
        active = state.get('active_scheme')
        if (patch.get('step') or {}).get('kind') in ('network', 'sp', 'cut') and active in state.get('schemes', {}):
            self._draw_polygons(ax, [state['schemes'][active]['outer']], state, color=NETWORK, linestyle=':',
                                linewidth=1.3, gid='step-network')
        # 2. 最近的步骤：射线、反向射线与锥远端面画线，求解点画标记
        steps = self._steps()
        trail = steps[:bisect_right(steps, self.index)][-STEP_TRAIL:]
        for rank, index in enumerate(trail, 1):
            step = self.monitor.history[index]['patch']['step']
            if step['kind'] not in STEP_STYLE:
                continue
            marker, color, label = STEP_STYLE[step['kind']]
            color = INNER if step.get('feasible') else color
            current = rank == len(trail)
            alpha = 1. if current else .2+.5*rank/len(trail)
            if 'vertices' in step:
                line = np.asarray(step['vertices'], float).reshape(-1, d)
                line = np.vstack([line, line[:1]]) if len(line) > 2 else line   # 三维远端面闭合
                ax.plot(*line.T, color=color, alpha=alpha, linewidth=1.8 if current else 1.,
                        linestyle='--' if step['kind'] == 'cone' else '-', gid='step-line')
            if 'p' in step:
                point = np.asarray(step['p'], float)
                shown = np.clip(point, lower, upper)
                outside = not np.allclose(shown, point)
                ax.scatter(*shown.reshape(1, d).T, marker=marker, s=70 if current else 28, alpha=alpha, zorder=11,
                           gid='step-point', **(dict(facecolors='none', edgecolors=color) if outside and marker not in 'xX'
                                                else dict(color=color)))
                if current:
                    text = label+('（坐标外）' if outside else '')
                    if d == 3:
                        ax.text(*shown, '  '+text, color=color, fontsize=9, fontweight='bold', zorder=12)
                    else:
                        ax.annotate(text, shown, xytext=(7, 7), textcoords='offset points', color=color, fontsize=9,
                                    fontweight='bold', zorder=12)

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
                self.view_signatures.pop(key, None)
                if len(self.scheme_views) == 4:
                    reusable = next(old for old in self.scheme_views if old not in visible)
                    self.scheme_views[key] = self.scheme_views.pop(reusable)
                    self.view_signatures.pop(reusable, None)
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
            # 内容签名：网架行、同分区最新的割、本网架的标记与是否当前网架都未变时不重画
            prefix = key.split(':')[0]+':'
            latest = next((k for k in reversed(state.get('cut_history', {})) if k.startswith(prefix)), None)
            step = state.get('step')
            mark = id(step) if step and step.get('scheme') == key else None
            signature = (id(row), latest, state.get('event') == 'cut', key == state.get('active_scheme'), mark)
            if self.view_signatures.get(key) == signature:
                continue
            self.view_signatures[key] = signature
            active = ' · 当前认证' if key == state.get('active_scheme') else ''
            frame.configure(text=f"{key} · {row['cost']:g} {state.get('cost_unit', '')}{active}")
            ax = self._axes(ax, state, scheme=key)
            self.scheme_views[key] = frame, ax, canvas
            self._draw_polygons(ax, [row['outer']], state, color=OUTER, fill=True,
                                alpha=.16, linestyle='--', gid='scheme-outer')
            if accepted(row):
                self._draw_polygons(ax, [row['outer']], state, color=INNER, fill=True,
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
        ax = self.axes['C'] = self._axes(self.axes['C'], state, validation=True)
        validation = state.get('validation')
        if validation is None:
            progress = state.get('scan_progress')
            message = (state.get('validation_note', '构域完成后独立扫描') if progress is None else
                       f'AC / SOCP 扫描  {progress[0]} / {progress[1]}  ({100*progress[0]/progress[1]:.1f}%)')
            if state.get('error'):
                message = ('校验未完成' if progress is not None else '构域未完成')+'\n'+state['error']
            write = ax.text2D if len(state['bounds']) == 3 else ax.text
            write(.5, .55, message, transform=ax.transAxes, ha='center', color='#777777')
        else:
            self._draw_comparison(ax, validation)
        self.canvases['C'].draw_idle()

    def _draw_comparison(self, ax, validation):
        """逐格对比：参考可行格为浅色（三维画其表面），遗漏格红色、多余格橙色；参考未决的格不参与。"""
        from matplotlib.colors import ListedColormap
        from matplotlib.lines import Line2D
        from matplotlib.patches import Patch
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
        # 1. 选定对比的计算域与参考域
        selected = next(key for key, label in COMPARISONS.items() if label == self.comparison_mode.get())
        ac, socp, result = (np.asarray(validation[key]) for key in ('states', 'socp_states', 'computed_states'))
        computed, reference, known = {'result_ac': (result, ac == 1, ac != 0),
                                      'result_socp': (result, socp == 1, socp != 0),
                                      'socp_ac': (socp == 1, ac == 1, (socp != 0) & (ac != 0))}[selected]
        missed, extra = known & reference & ~computed, known & computed & ~reference
        # 2. 二维：格子着色；三维：参考域表面与遗漏/多余格心
        lower, upper = np.asarray(validation['axis_lower']), np.asarray(validation['bounds'])
        step = (upper-lower)/np.array(ac.shape)
        if ac.ndim == 2:
            classes = np.where(missed, 3, np.where(extra, 2, np.where(known & reference, 1, 0)))
            ax.pcolormesh(*(np.linspace(a, b, n+1) for a, b, n in zip(lower, upper, ac.shape)), classes.T,
                          cmap=ListedColormap(['white', REFERENCE, SP, GLOBAL]), vmin=0, vmax=3, shading='flat',
                          rasterized=True)
            handles = [Patch(facecolor=REFERENCE, label='参考可行'), Patch(facecolor=GLOBAL, label='遗漏'),
                       Patch(facecolor=SP, label='多余')]
        else:
            ax.add_collection3d(Poly3DCollection(lower+voxel_faces(known & reference, upper-lower),
                                                 facecolors=REFERENCE, edgecolors='none', alpha=.25, gid='reference'))
            for cells, color, gid in ((missed, GLOBAL, 'missed-cells'), (extra, SP, 'extra-cells')):
                points = lower+(np.argwhere(cells)+.5)*step
                ax.scatter(*points.T, color=color, s=2, alpha=.6, marker='.', depthshade=False, gid=gid)   # 小而半透明：看疏密
            handles = [Patch(facecolor=REFERENCE, label='参考可行（表面）'),
                       Line2D([], [], color=GLOBAL, marker='.', ls='none', label='遗漏'),
                       Line2D([], [], color=SP, marker='.', ls='none', label='多余')]
        # 3. 三组指标；本组遗漏/多余格按符号分区计数（从多到少），看误差主要来自哪些分区
        fmt = lambda value: '—' if value is None else f'{value:.3f}%'

        def by_partition(cells):
            labels = [''.join('+' if v >= 0 else '-' for v in point) for point in lower+(np.argwhere(cells)+.5)*step]
            counts = sorted(((labels.count(label), label) for label in set(labels)), reverse=True)
            return ' · '.join(f'{label} {count}' for count, label in counts) or '—'
        self.comparison_text.set('\n'.join([*(f'{label}：遗漏 {fmt(validation["comparisons"][key]["mr_percent"])}'
            f'    多余 {fmt(validation["comparisons"][key]["fr_percent"])}' for key, label in COMPARISONS.items()),
            f'本组遗漏按分区：{by_partition(missed)}', f'本组多余按分区：{by_partition(extra)}']))
        ax.set_title(f"{COMPARISONS[selected]} · 遗漏 {int(missed.sum())} 格 · 多余 {int(extra.sum())} 格 · "
                     +'×'.join(map(str, ac.shape))+' 网格', fontsize=9)
        ax.legend(handles=handles, loc='upper right', frameon=False, fontsize=8)

    def _removed_3d(self, state, scheme=None):
        """当前割在各网架旧 Nx 中实际切掉的部分，仅用于本帧显示。"""
        from region import clip_polytope
        if state.get('event') != 'cut':
            return []
        before = self.monitor.frame(self.index-1).get('schemes', {})
        latest = next(reversed(state['cut_history'].values()))
        cut = np.asarray(latest['cut'])
        bounds = np.asarray(state['bounds'])
        keys = before if scheme is None else [scheme]
        return [clip_polytope(np.asarray(before[key]['outer'])/bounds,
                             -cut[0]-cut[4:]@before[key]['x'], -cut[1:4]*bounds)*bounds
                for key in keys if key in before and np.array_equal(before[key].get('sign'), latest.get('sign'))]


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
