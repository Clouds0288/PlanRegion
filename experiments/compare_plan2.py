"""方案二（支撑查询）加入后的对照实验：方法 R（径向锥 + 体积准则 + B 型两网架证书）对比方法 H（R + 逐网架支撑
查询 + 覆盖证书），同一物理模型（每个网架带 OBBT 盒与包络行）、同一 tau/ε/总时限/Gurobi 参数与种子。

在项目根目录运行::

    python -X utf8 experiments/compare_plan2.py --smoke
    python -X utf8 experiments/compare_plan2.py --case case33 --nodes 18,25 --mode 1 --methods R,H --seconds 300 \
        --repeats 3 --workers 1,16 --threads 1 --output results/compare_plan2/case33_18_25
    python -X utf8 experiments/compare_plan2.py --report --output results/compare_plan2/case33_18_25

方法 R
----
径向锥夹逼原样调用（test_radial_sandwich.py 的 RadialSandwich、test_radial_sandwich_3d.py 的 RadialSandwich3D，
--criterion volume，近端覆盖即 B 型两网架证书 T ⊆ R_x̂ ∪ R_y）：每次细分体积缺口 Δ_k=(μ̄_k^d-1)·vol(T_k) 最大的锥，
ΣΔ_k <= ε·Σvol(T_k) 时认证。二维的体积准则是三维 run_volume 的二维实例（随本实验加入二维脚本，默认 radial 不变）。

方法 H（同一总时限内 A → B → C，可循环）
----
A  与 R 同一实现、同一种子，目标放宽为 ε_A（--discovery-eps），最多用分区时限的 share_A（--discovery-share），先到
   先停，故 A 结束时 H 的状态即 R 在同一时刻的状态。收集 X*（出现过的全部网架：锥的 x̂、覆盖网架 y、外界 MISOCP
   现任解 x*，以及已 OBBT、做过射线或零接入认证的方案）与种子 V_x（该网架射线的远端与近端、已认证的零接入点、
   通过 audit_incumbent 的已紧化现任解）；OBBT 盒原样复用。
B  逐网架、x 固定只解 SOCP：P_x=conv(V_x)，按 certify_scheme 的逻辑（先用 O_x 顶点的几何上界检查全部面，再查必要
   支撑方向；可靠 UB 认证该面，审计过的越界点补进 V_x 并立即重建面、废弃旧面结论）。支撑 SOCP 为 SupportOracle
   （固定方案的完整 MP，含本网架包络行）加分区盒 0<=u<=bounds，只接受 OPTIMAL。面判据 --criterion origin/center/auto
   见 plan2_geometry。网架按 A 中所占锥体积排队，面按 面积×(UB-β) 从大到小；--workers>1 时用 spawn 进程池并行
   （各进程自建 Gurobi 环境）。网架在 全部非分区边界面认证 / vol(O_x)/vol(P_x)-1<=ε_B / 时间片用完 时停止。
C  x 自由的覆盖证书：RemainingRegionModel 直接接收外扩面 E_x 并令 s=1（tau=0），再按径向脚本的做法加全部已紧化
   方案的提升行（汉明距离 big-M，BarQCPConvTol=1e-6），见证方案未紧化时懒惰 OBBT 后重解。可靠上界 <= GEOMETRY_TOL
   或已证不可行即覆盖完成，冻结当时的 E_x；否则见证点（审计通过，或从内点锚出发的紧化射线认证）补进 X* / V_x，回 B。
输出：I_H = I_R ∪ (∪_x P_x)；覆盖完成前 O_H = O_R，完成后 O_H = O_R ∩ (∪_x E_x)（均与分区盒相交）。
认证：覆盖完成、没有 UNRESOLVED 面且 vol(O_H)-vol(I_H) <= ε·vol(I_H)；或径向部分自身满足 R 的体积准则（此时
O_H ⊆ O_R、I_H ⊇ I_R，H 继承 R 的证书）。覆盖完成而间隙仍大于 ε 时，续跑 A（目标 ε），每次锥决策后复核 H 的间隙。

实现约定（任务说明未覆盖或需取舍之处）
----
1. 时间：--seconds 是每次运行的墙钟总时限（含 OBBT、求解、几何与记录，不含事后逐格评价），t=0 为运行开始。分区
   沿用径向脚本的独立分区时限：并发数 C=min(workers, 分区数)；串行时第 j 个分区得到剩余时间的 1/(余下分区数)（与
   主线 build_region 相同的顺延分配），并行时各分区同时开始、共用总时限。分区内并行度 max(1, workers//分区数)，
   用作 OBBT 线程数与 B 阶段进程池大小；每个求解的 Threads 统一为 --threads。
2. 有效性检查（逐格必要条件，任一不满足即在报告中标红、以退出码 3 结束）：内域不含 SOCP 已证不可行格（紧化模型
   ⊆ SOCP 松弛）；AC 可行格都在外界内（AC ⊆ 紧化模型 ⊆ 外界）；AC 未决格单列。SOCP 可行而在外界外的格是 OBBT
   去掉的 SOCP-only 区域，只作诊断：主指标是相对 AC 的 FR/MR，SOCP 只用来分辨误差来自松弛还是模型。
3. 代码库没有求解重试流程（主线失败即停、SupportOracle 不在失败后改精度），故非 OPTIMAL 的支撑 SOCP 不重求，该面
   记 UNRESOLVED。“未决面”指 UNRESOLVED 面；因 ε_B 或时间片停止而未查的面记 unchecked，不阻止认证（外界的有效性
   由覆盖证书保证）。覆盖见证落在 x∈X* 时“对这个网架重新做 B”：补点后该网架不再按 ε_B 提前停止，查完全部非分区
   边界面才停（否则补点只会让 vol(O_x)/vol(P_x)-1 更小、立即再次停止，未查面附近的见证要靠一轮轮覆盖 MISOCP 逐点补）。
4. 每轮 B 最多用剩余时间的 SUPPORT_PASS_SHARE，其余留给 C；单次支撑 SOCP 不超过 SUPPORT_SECONDS，覆盖 MISOCP 不超过
   --mip-seconds（与锥 MISOCP 相同）。--coverage cone 未实现（任务说明允许先不做）。
5. 检查点快照是该时刻的确切状态：每次几何变化之前处理已到的检查点（锥求解期间到达的检查点按求解前的外界记录）。
   三维 H 的测度按锥分解求并集体积，时间线行最多每 TIMELINE_INTERVAL_3D 秒一行（认证判定与检查点不受影响）。

记号（docs/notation.md「方案二对照实验」）
----
xi=u/bounds；P_x、O_x、E_x、c_x、V_x、X* 同上；ε=--eps（默认 d·tau，即径向脚本的体积目标），ε_A=--discovery-eps，
share_A=--discovery-share，ε_B=--network-eps（默认 ε/2）。测度 *_measure 为 kW^d，*_kw 坐标为带符号真实 kW。

输出
----
<output>/<method>/workers_<w>/run_<r>/：summary.json、solves.csv（每次求解一行）、timeline.csv（内外测度每变化一次
一行；partition=all 为全部分区之和）、snapshots/（各检查点与 final 的几何）、metrics.json（各检查点的逐格指标）。
<output>/：comparison.csv、comparison.md、gap_time.png、mr_fr_time.png、regions.png（二维）、time_breakdown.png。
结果不加入 results/manifest.json。
"""
from __future__ import annotations

import argparse
from argparse import Namespace
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from contextlib import contextmanager
import csv
from copy import deepcopy
from itertools import product
import json
from math import factorial
import multiprocessing
import os
from pathlib import Path
import subprocess
import sys
import time
from time import perf_counter
import traceback

import gurobipy as gp
from gurobipy import GRB
import numpy as np
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
for _path in (str(ROOT), str(HERE)):
    if _path not in sys.path:
        sys.path.insert(0, _path)

from main import BUDGET, REGION_TAU
from model import PLANNING_TOL, PortPhysics, RemainingRegionModel, port_bounds, ray_support
from Network.case33bw import Case33, CURRENT_LIMIT
from Network.four_bus_five_corridor import FourBus
from region import GEOMETRY_TOL, clip_polytope, contains, halfspaces, polytope_vertices, polytope_volume
from vertify import ac_identity, ac_network, load_scan, scan_path
import plan2_geometry as geo
from test_radial_sandwich import (COVER_CONV_TOL, RAY_SECONDS, PartitionTimeout, RadialSandwich, direction,
                                  git_commit, partition_sign, plain)
from test_radial_sandwich_3d import RadialSandwich3D
from test_support_face_certification_fourbus_2d import (SupportOracle, audit_incumbent, cached_support,
                                                         classify_support, normal_key)

CHECKPOINTS = (15., 30., 60., 120., 180., 240., 300.)
SUPPORT_SECONDS = 30.        # 单次支撑 SOCP 的时限上限（同时受分区剩余时间限制）
SUPPORT_PASS_SHARE = .75     # 每轮 B 最多占剩余时间的比例，其余留给覆盖证书 C
TIMELINE_INTERVAL_3D = 10.   # 三维 H 的时间线最短间隔（秒）：每行要按锥求并集体积，检查点不受影响
MAX_ROUNDS = 200             # B→C 循环轮数上限
GAP_LEVELS = (.10, .05, .02)
DEFAULT_OUTPUT = ROOT/'results'/'compare_plan2'
SMOKE = dict(case='case33', nodes=(18, 25), partitions=('nn',), seconds=60., repeats=1, methods=('R', 'H'),
             workers=(1, 16))
CERTIFIED_FACES = ('BOUNDARY', 'GEOMETRY_CERTIFIED', 'SUPPORT_CERTIFIED')
SURFACE, INK, INK_2, MUTED, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#898781', '#e1e0d9'
METHOD_COLORS = dict(R='#2a78d6', H='#eb6834', Hc='#1baf7a')   # 分类色前三槽（径向脚本已校验的配色）
PURPOSES = ('discovery', 'cone_outer', 'support', 'coverage', 'obbt', 'ray', 'zero_sp')
PURPOSE_COLORS = dict(zip(PURPOSES, ('#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#4a3aa7')))
TYPE_COLORS = dict(MISOCP='#2a78d6', SOCP='#eb6834', LP='#1baf7a')
OTHER_COLOR = '#c3c2b7'


# ---- 算例与参考扫描 ------------------------------------------------------------------------------
class Case:
    """算例：网络、预算、公共尺度 bounds（kW）与符号分区代号。"""

    def __init__(self, name, nodes, mode=1):
        if mode != 1:
            raise SystemExit('The radial scripts build signed partitions: only --mode 1 is supported')
        self.name, self.nodes, self.mode, self.d = name, tuple(int(v) for v in nodes), mode, len(nodes)
        if self.d not in (2, 3):
            raise SystemExit('Only d=2 or d=3 is supported')
        if name == 'case33':
            self.network = Case33(load_nodes=self.nodes, current_limit=CURRENT_LIMIT)
            self.budget = float(Case33.switch_budget)
        elif name == 'fourbus':
            self.network = FourBus(load_nodes=self.nodes)
            self.budget = float(BUDGET)
        else:
            raise SystemExit(f'Unknown case {name!r}')
        self.bounds = port_bounds(self.network)
        self.partitions = tuple(''.join(code) for code in product('pn', repeat=self.d))
        self.tag = f'{name}_{"_".join(map(str, self.nodes))}'


def load_reference(case, path=None):
    """项目标准的配对 AC/SOCP 扫描缓存（vertify），校验物理身份与网络指纹；没有缓存时返回 None。"""
    if path is None:
        directory = scan_path(ac_network(case.network), case.budget, ROOT/'results'/'scan', case.mode)
        sources = sorted(directory.glob('region_*.npz'), key=lambda p: p.stat().st_mtime_ns)
        if not sources:
            return None
        path = sources[-1]
    scan = load_scan(path, ac_identity(ac_network(case.network), case.budget, case.mode))
    if case.network.fingerprint != scan['metadata']['network_fingerprint']:
        raise SystemExit(f'Network fingerprint differs from scan cache: {path}')
    shape = scan['states'].shape
    index = np.indices(shape).reshape(len(shape), -1).T
    scan.update(path=Path(path), shape=shape, index=index, power=scan['origin']+(scan['start']+index+.5)*scan['step'],
                ac=scan['states'].ravel(), socp=scan['socp_states'].ravel())
    return scan


# ---- 时钟、求解日志、时间线与检查点 -----------------------------------------------------------------
class Clock:
    """运行时钟：t=0 为本次运行开始，各分区进程共用同一 epoch，计时用 perf_counter。"""

    def __init__(self, epoch):
        self.offset = perf_counter()-(time.time()-epoch)

    def now(self):
        return perf_counter()-self.offset

    def perf(self, t):
        return t+self.offset


def empty_state():
    return dict(started=False, phase=None, cones=[], networks=[], cover=None, certified=False)


def checkpoint_times(seconds):
    return sorted({t for t in CHECKPOINTS if t < seconds} | {float(seconds)})


class Recorder:
    """一个分区的求解日志、时间线与检查点快照；构域期间只记录，不做逐格评价。"""

    def __init__(self, code, method, clock, checkpoints, min_interval=0.):
        self.code, self.method, self.clock, self.min_interval = code, method, clock, min_interval
        self.phase = 'A' if method == 'H' else 'R'
        self.solves, self.timeline, self.snapshots, self.stack = [], [], {}, []
        self.state_fn = lambda override=None: empty_state()
        self.measure = None
        self.last = None
        start = clock.now()
        self.pending = []
        for t in sorted(checkpoints):
            if t <= start:
                self.snapshots[t] = empty_state()   # 分区尚未开始（串行运行中排在后面的分区）
            else:
                self.pending.append(t)

    def _row(self, solver, purpose, info):
        row = dict(partition=self.code, method=self.method, phase=self.phase, type=solver, purpose=purpose,
                   status='ok', process='main')
        row.update(info)
        return row

    @contextmanager
    def log(self, solver, purpose, **info):
        """一次本进程求解：耗时扣除其中嵌套记录的求解（如射线内的 OBBT 单独成行）。"""
        started, row = perf_counter(), self._row(solver, purpose, info)
        row['t_start'] = self.clock.now()
        self.stack.append(0.)
        try:
            yield row
        except PartitionTimeout:
            row['status'] = 'timeout'
            raise
        except Exception as exc:
            row.update(status='error', error=repr(exc)[:300])
            raise
        finally:
            nested, elapsed = self.stack.pop(), perf_counter()-started
            if self.stack:
                self.stack[-1] += elapsed
            row.update(t_end=self.clock.now(), seconds=max(elapsed-nested, 0.), wall_seconds=elapsed)
            self.solves.append(row)

    def record(self, solver, purpose, seconds, *, remote=False, **info):
        """事后记录一次求解（MISOCP 轮次、进程池支撑）；remote 为其他进程的求解，不计入本进程的嵌套耗时。"""
        t = self.clock.now()
        if self.stack and not remote:
            self.stack[-1] += seconds
        row = self._row(solver, purpose, info)
        row.update(t_start=t-seconds, t_end=t, seconds=seconds, wall_seconds=seconds,
                   process='pool' if remote else 'main')
        self.solves.append(row)

    def before_change(self, override=None, empty=False):
        """几何即将变化：已到的检查点取当前（变化前）的状态。"""
        t = self.clock.now()
        while self.pending and self.pending[0] <= t:
            self.snapshots[self.pending.pop(0)] = empty_state() if empty else self.state_fn(override)

    def changed(self, force=False):
        """几何刚变化：时间线追加一行（内外测度 kW^d、间隙、是否认证）；与上一行相同则不记。"""
        t = self.clock.now()
        if not force and self.min_interval and self.last is not None and t-self.last['t'] < self.min_interval:
            return
        inner, outer, certified, extra = self.measure()
        row = dict(t=t, partition=self.code, method=self.method, phase=self.phase, inner_measure=inner,
                   outer_measure=outer, gap=outer/inner-1. if inner > 0 else np.inf, certified=bool(certified), **extra)
        if self.last is not None and all(row[k] == self.last[k] for k in ('inner_measure', 'outer_measure',
                                                                           'certified', 'phase')):
            return
        self.timeline.append(row)
        self.last = row

    def finish(self):
        if self.measure is not None:
            self.changed(force=True)
        state = self.state_fn(None)
        for t in self.pending:
            self.snapshots[t] = state
        self.pending = []
        self.snapshots['final'] = state


# ---- 径向夹逼的观察钩子（方法 R 与 H 的 A 阶段共用，不改变决策） ----------------------------------------
def leaves(radial):
    return radial.cones if isinstance(radial, RadialSandwich3D) else radial.intervals


def cone_rows(radial, override=None):
    """叶锥的方向 U、内顶点 V（列）与有效外界 (c, mu)；override=(锥, 旧 mu) 时该锥按求解前的外界导出。
    根锥建立前以一个无外界的整卦限锥代替（内域为空、外域为分区盒）。"""
    if not leaves(radial):
        return [dict(U=np.eye(radial.dim), V=1e-12*np.eye(radial.dim), halfspace=None, scheme=None, cover=None,
                     pseudo=True)]
    rows = []
    three = isinstance(radial, RadialSandwich3D)
    for cone in leaves(radial):
        U = radial.U(cone.keys) if three else np.column_stack([direction(cone.theta_a), direction(cone.theta_b)])
        halfspace = cone.halfspace()
        if override is not None and cone is override[0]:
            halfspace = (cone.c, override[1]) if np.isfinite(override[1]) else cone.inherited
        rows.append(dict(U=U, V=cone.Q, halfspace=None if halfspace is None else
                         (np.asarray(halfspace[0], float), float(halfspace[1])),
                         scheme=radial.label(cone.x), cover=None if cone.cover is None else radial.label(cone.cover)))
    return rows


def radial_measure(radial):
    """径向内域与外域∩盒的测度（kW^d）。"""
    if isinstance(radial, RadialSandwich3D):
        inner, _, outer = radial.volumes()
    else:
        inner, _, outer = radial.areas()
    return inner, outer


def radial_state(radial, phase, override=None):
    return dict(started=True, phase=phase, cones=cone_rows(radial, override), networks=[], cover=None,
                certified=radial.status == 'certified')


class Instrumented:
    """观察径向夹逼：求解记录、时间线、检查点与（H 的 A 阶段）现任解审计；只读状态，不改变它的任何决策。"""
    audit = False

    def attach(self, recorder, seed):
        self.recorder, self.seed = recorder, seed
        self.incumbent_points, self.incumbent_schemes = [], set()
        self.after_change, self._context = None, None

    def label(self, x):
        """Case33 沿用开关串；其他网络按走廊写所选型号序号（0 为不建）。"""
        if self.network.name == 'case33bw':
            return super().label(x)
        x = np.asarray(x)
        return ''.join(str(int(np.argmax(x[s]))+1) if x[s].any() else '0' for s in self.network.type_slices)

    def tighten(self, x):
        scheme = tuple(int(v) for v in x)
        if scheme in self.equations.boxes:
            return 0.
        with self.recorder.log('SOCP', 'obbt', scheme=self.label(scheme)):
            return super().tighten(x)

    def origin_feasible(self, x):
        scheme = tuple(int(v) for v in x)
        if scheme in self.origin:
            return self.origin[scheme]['feasible']
        self.tighten(scheme)
        with self.recorder.log('SOCP', 'zero_sp', scheme=self.label(scheme)) as row:
            feasible = super().origin_feasible(x)
            row['status'] = 'feasible' if feasible else 'failed' if 'error' in self.origin[scheme] else 'infeasible'
        return feasible

    def vertex(self, x, theta):
        scheme = tuple(int(v) for v in x)
        if (scheme, float(theta)) in self.rays:
            return super().vertex(x, theta)
        self.tighten(scheme)
        with self.recorder.log('SOCP', 'ray', scheme=self.label(scheme), ray='far') as row:
            point = super().vertex(x, theta)
            row['status'] = 'ok' if point is not None else 'failed'
        return point

    def near(self, x, theta):
        scheme = tuple(int(v) for v in x)
        if (scheme, float(theta)) in self.nears or self.vertex(x, theta) is None:
            return super().near(x, theta)
        with self.recorder.log('SOCP', 'ray', scheme=self.label(scheme), ray='near') as row:
            point = super().near(x, theta)
            row['status'] = 'ok' if point is not None else 'failed'
        return point

    def misocp(self, configure, *, stop, context):
        def seeded(problem, xi):
            configure(problem, xi)
            problem.model.Params.Seed = self.seed
        self._context = context
        return super().misocp(seeded, stop=stop, context=context)

    def on_solution(self, problem, record):
        purpose = 'discovery' if self._context == self.CENTER_CONTEXT else 'cone_outer'
        self.recorder.record('MISOCP', purpose, record['build_seconds']+record['solve_seconds'], status=record['status'],
                             scheme=None if record['x'] is None else self.label(record['x']), context=self._context)
        if record['x'] is None:
            return
        self.incumbent_schemes.add(record['x'])
        if self.audit and record['x'] in self.equations.boxes:
            # 只有已紧化方案的现任解才是该网架紧化可行集中的点；原约束逐项代入，不修补。
            residual = audit_incumbent(problem, np.asarray(record['x']))
            record['audit_residual'] = residual
            if residual <= PLANNING_TOL:
                self.incumbent_points.append((record['x'], np.asarray(record['point'], float)/self.bounds))

    def make_root(self, *args):
        self.recorder.before_change()
        return super().make_root(*args)

    def split(self, cone, final):
        children = super().split(cone, final)
        if children is not None:
            self.recorder.before_change()   # run_volume 随即用子锥替换 cone
        return children

    def solve_interval(self, interval):
        return self._observed(interval, super().solve_interval)

    def solve_cone(self, cone):
        return self._observed(cone, super().solve_cone)

    def _observed(self, cone, solve):
        self.recorder.before_change()
        old = cone.mu
        try:
            return solve(cone)
        finally:
            if self.recorder.pending and self.recorder.pending[0] <= self.recorder.clock.now():
                self.recorder.before_change(override=(cone, old))   # 求解期间到达的检查点：本锥仍为旧外界
            self.recorder.changed()

    def snapshot(self, cone, decision, scan=None):
        super().snapshot(cone, decision, scan)
        if decision == 'initial':
            self.recorder.before_change(empty=True)   # 根锥刚建立，此前为空
        self.recorder.changed()
        if self.after_change is not None:
            self.after_change()


class Radial2D(Instrumented, RadialSandwich):
    """二维径向夹逼（体积准则）加观察钩子。"""


class Radial3D(Instrumented, RadialSandwich3D):
    """三维径向夹逼（体积准则）加观察钩子。"""


def make_radial(case, code, args, workers, epsilon, deadline_t, recorder, clock):
    defaults = (1e-4, 256) if case.d == 2 else (2e-3, 2048)
    namespace = Namespace(tau=args.tau, mip_seconds=args.mip_seconds, seconds=max(deadline_t-clock.now(), 0.),
                          threads=args.threads, obbt_workers=workers, inner_cert=args.inner_cert, criterion='volume',
                          epsilon=epsilon, mip_gap=args.mip_gap,
                          min_width=args.min_width if args.min_width is not None else defaults[0],
                          max_intervals=args.max_cones or defaults[1], max_cones=args.max_cones or defaults[1])
    radial = (Radial2D if case.d == 2 else Radial3D)(case.network, code, case.bounds, case.budget, namespace)
    radial.deadline = clock.perf(deadline_t)
    radial.attach(recorder, args.seed)
    return radial


# ---- 方法 R ------------------------------------------------------------------------------------
def run_R(case, code, args, clock, deadline_t, workers):
    recorder = Recorder(code, 'R', clock, checkpoint_times(args.seconds))
    start = clock.now()
    radial = make_radial(case, code, args, workers, args.eps, deadline_t, recorder, clock)
    recorder.state_fn = lambda override=None: radial_state(radial, 'R', override)
    recorder.measure = lambda: (*radial_measure(radial), radial.status == 'certified',
                                dict(volume_ratio=radial.volume_ratio() if leaves(radial) else None,
                                     cones=len(leaves(radial))))
    with threadpool_limits(limits=1):
        radial.run(None)
    end = clock.now()
    recorder.finish()
    certified = radial.status == 'certified'
    inner, outer = radial_measure(radial)
    result = dict(partition=code, method='R', status=radial.status, certified=certified,
                  t_cert=end if certified else None, start=start, end=end, seconds=end-start,
                  inner_measure=inner, outer_measure=outer, gap=outer/inner-1. if inner > 0 else None,
                  volume_ratio=radial.volume_ratio() if leaves(radial) else None, cones=len(leaves(radial)),
                  schemes=len({c.x for c in leaves(radial)}), counts=dict(radial.counts),
                  errors=radial.errors[-5:], phase_seconds=dict(R=end-start))
    return result, recorder


# ---- 方法 H：B 阶段的支撑求解 -----------------------------------------------------------------------
class PartitionOracle(SupportOracle):
    """SupportOracle（固定方案的完整 MP，含本网架包络行）加分区盒 0<=u<=bounds 与种子；只接受 OPTIMAL。"""

    def __init__(self, equations, x, bounds, budget, threads, time_limit, seed=0):
        super().__init__(equations, np.asarray(x, int), bounds, budget, threads, time_limit)
        self.problem.power.UB = self.bounds
        self.problem.model.Params.Seed = seed

    def solve(self, normal):
        started = perf_counter()
        unit = np.asarray(normal, float)/np.linalg.norm(normal)
        try:
            entry = super().solve(normal)
        except (gp.GurobiError, RuntimeError) as exc:
            return dict(lb=-np.inf, ub=np.inf, point=None, state=None, status=-1, solve_seconds=perf_counter()-started,
                        normal=unit, max_violation=np.inf, error=repr(exc)[:300])
        if entry['status'] != GRB.OPTIMAL:
            entry.update(lb=-np.inf, ub=np.inf, point=None, state=None)
        return entry


def slim(entry):
    """进程间传递的支撑记录：不带运行状态向量。"""
    return {k: v for k, v in entry.items() if k not in ('state',)}


_WORKER = {}


def _init_worker(case_name, nodes, code, threads, seed):
    from multiprocessing.util import Finalize
    threadpool_limits(limits=1)
    case = Case(case_name, nodes)
    _WORKER.update(equations=PortPhysics(deepcopy(case.network), np.asarray(partition_sign(code, case.d))),
                   bounds=case.bounds, budget=case.budget, threads=threads, seed=seed, oracles={})
    gp.Model().dispose()   # 提前建立本进程的 Gurobi 环境
    Finalize(None, _close_worker, exitpriority=10)   # 进程退出前释放模型，环境随后正常释放


def _close_worker():
    for oracle in _WORKER.get('oracles', {}).values():
        oracle.close()
    _WORKER['oracles'] = {}


def _warm():
    return os.getpid()


def _support_task(x, box, normal, limit):
    equations = _WORKER['equations']
    x = tuple(int(v) for v in x)
    equations.boxes.setdefault(x, box)
    oracle = _WORKER['oracles'].get(x)
    if oracle is None:
        oracle = _WORKER['oracles'][x] = PartitionOracle(equations, x, _WORKER['bounds'], _WORKER['budget'],
                                                          _WORKER['threads'], limit, _WORKER['seed'])
    oracle.time_limit = limit
    return slim(oracle.solve(np.asarray(normal, float)))


class SupportPool:
    """B 阶段进程池（spawn）；分区开始时即提交预热任务，进程在 A 阶段期间完成导入与 Gurobi 环境。"""

    def __init__(self, case, code, args, workers):
        self.executor = ProcessPoolExecutor(workers, mp_context=multiprocessing.get_context('spawn'),
                                            initializer=_init_worker,
                                            initargs=(case.name, case.nodes, code, args.threads, args.seed))
        self.warm = [self.executor.submit(_warm) for _ in range(workers)]

    def submit(self, x, box, normal, limit):
        return self.executor.submit(_support_task, tuple(int(v) for v in x), box, np.asarray(normal, float), limit)

    def close(self):
        self.executor.shutdown(wait=True, cancel_futures=True)


class NetworkState:
    """固定方案 x：内多面体 P_x=conv(V_x)、支撑外界 O_x（分区盒逐次裁剪）与逐面认证（certify_scheme 的逻辑）。"""

    def __init__(self, x, label, points, weight, d, args, recorder=None):
        self.x, self.label, self.weight, self.d = tuple(int(v) for v in x), label, float(weight), d
        self.recorder = recorder
        self.tau, self.mode, self.epsilon_geom = args.tau, args.criterion, GEOMETRY_TOL
        self.network_eps = args.network_eps
        self.points = np.asarray(points, float).reshape(-1, d)
        self.outer = geo.box_vertices(d)
        self.cache, self.inflight, self.lost = {}, set(), set()
        self.status, self.reason, self.urgent = 'active', None, False
        self.strict = False   # 出过覆盖见证：不再按 ε_B 提前停止，直到全部非分区边界面认证
        self.calls = self.violations = self.rebuilds = self.geometry_failures = 0
        self.seconds = 0.
        self.slice_end = None
        while True:
            try:
                self.rebuild()
                break
            except geo.GEOMETRY_ERRORS:   # 近退化的初始种子：去掉最后一个点重试
                self.points, self.geometry_failures = self.points[:-1], self.geometry_failures+1

    # 几何与面状态
    def rebuild(self):
        self.rebuilds += 1
        d = self.d
        self.vertices = polytope_vertices(self.points) if len(self.points) else np.empty((0, d))
        self.empty = not len(self.vertices)
        if self.empty:   # 尚无认证点：先查 ±e_j 与 ±1/√d 方向（certify_scheme 的初始方向及其反向）
            self.seeds = np.vstack([np.eye(d), -np.eye(d), np.ones(d)/np.sqrt(d), -np.ones(d)/np.sqrt(d)])
            self.faces, self.full = np.empty((0, d+1)), False
            self.face_status, self.priority, self.upper = np.empty(0, object), np.empty(0), np.empty(0)
            self.criterion, self.center, self.center_source = 'center', np.zeros(d), None
            return
        self.faces = halfspaces(self.vertices)
        self.full = len(self.vertices) > d and polytope_volume(self.vertices) > 0.
        self.boundary = geo.boundary_faces(self.faces)
        self.criterion = geo.choose_criterion(self.faces, self.boundary, self.mode) if self.full else 'center'
        if self.recorder is None:
            self.center, self.center_source, self.radius = geo.interior_point(self.vertices, self.faces)
        else:
            with self.recorder.log('LP', 'support', scheme=self.label, ray='chebyshev') as row:
                self.center, self.center_source, self.radius = geo.interior_point(self.vertices, self.faces)
                row['status'] = self.center_source
        self.areas = geo.face_measures(self.vertices, self.faces) if self.full else np.ones(len(self.faces))
        self.evaluate()

    def classify(self, f, entry):
        normal, offset = self.faces[f, :-1], self.faces[f, -1]
        if self.criterion == 'origin':
            return classify_support(entry['lb'], entry['ub'], entry['point'], normal, offset, self.tau,
                                    self.epsilon_geom)[0]
        return geo.classify_center(entry['lb'], entry['ub'], entry['point'], normal, offset, self.center, self.tau,
                                   self.epsilon_geom)[0]

    def evaluate(self):
        """几何先行：O_x 顶点给出各面支撑上界，再用缓存的支撑结果；状态为 BOUNDARY / GEOMETRY_CERTIFIED /
        SUPPORT_CERTIFIED / UNRESOLVED（求解失败或补点在几何精度下丢失）/ PENDING。"""
        faces = self.faces
        upper = (self.outer@faces[:, :-1].T).max(axis=0) if len(self.outer) else np.full(len(faces), np.inf)
        margins = geo.face_margins(faces, upper, self.criterion, self.center, self.tau)
        status = np.where(self.boundary, 'BOUNDARY',
                          np.where(margins <= self.epsilon_geom, 'GEOMETRY_CERTIFIED', 'PENDING')).astype(object)
        for f in np.flatnonzero(status == 'PENDING'):
            key = normal_key(faces[f, :-1])
            if key not in self.cache:
                continue
            entry = cached_support(self.cache[key], faces[f, :-1], np.ones(self.d))
            upper[f] = min(upper[f], entry['ub'])
            label = self.classify(f, entry)
            status[f] = 'UNRESOLVED' if label == 'VIOLATED' or key in self.lost else label
        self.face_status, self.upper = status, upper
        self.priority = self.areas*np.maximum(upper+faces[:, -1], 0.)   # 面积×(UB-β)

    def next_tasks(self, limit):
        if limit <= 0:
            return []
        if self.empty:
            return [n for n in self.seeds if normal_key(n) not in self.cache and normal_key(n) not in self.inflight][:limit]
        order = np.argsort(-self.priority, kind='stable')
        tasks = []
        for f in order:
            if self.face_status[f] != 'PENDING':
                continue
            key = normal_key(self.faces[f, :-1])
            if key in self.inflight or key in self.cache:
                continue
            tasks.append(self.faces[f, :-1].copy())
            if len(tasks) >= limit:
                break
        return tasks

    def add_points(self, points):
        """补点并重建；凸包在近退化输入上数值失败时撤回补点（返回 False，调用方把该面记 UNRESOLVED）。"""
        before, previous = self.vertices, self.points
        self.points = np.vstack([self.points, np.asarray(points, float).reshape(-1, self.d)])
        try:
            self.rebuild()
        except geo.GEOMETRY_ERRORS:
            self.points = previous
            self.rebuild()
            self.geometry_failures += 1
            return False
        return self.empty or len(self.vertices) != len(before) or not np.array_equal(self.vertices, before)

    def on_result(self, normal, entry):
        """缓存支撑结果并裁剪 O_x；审计过的点越过某个当前面的允许外扩时补进 V_x 并重建（废弃旧面结论）。"""
        key = normal_key(normal)
        self.inflight.discard(key)
        self.cache[key] = entry
        self.calls += 1
        self.seconds += float(entry.get('solve_seconds', 0.))
        unit = np.asarray(entry['normal'], float)
        if np.isfinite(entry['ub']):
            try:
                self.outer = clip_polytope(self.outer, entry['ub'], -unit)
            except geo.GEOMETRY_ERRORS:
                self.geometry_failures += 1   # 保留裁剪前的 O_x：仍是有效外界，只是不够紧
        point = entry['point']
        if point is not None:
            point = np.asarray(point, float)
            if self.empty:
                self.add_points([point])
                return
            allowed = geo.allowed_offsets(self.faces, self.criterion, self.center, self.tau)
            if np.any(self.faces[:, :-1]@point-allowed > self.epsilon_geom):
                self.violations += 1
                if not self.add_points([point]):
                    self.lost.add(key)
                    self.evaluate()
                return
        if not self.empty:
            self.evaluate()

    def check_done(self):
        if self.status != 'active':
            return
        if self.empty:
            if not self.next_tasks(1) and not self.inflight:
                self.status, self.reason = 'done', 'no_points'
            return
        status = self.face_status
        if np.isin(status, CERTIFIED_FACES).all():
            self.status, self.reason = 'done', 'certified'
        elif not self.strict and self.full and len(self.outer) > self.d and self.volume_ratio() <= self.network_eps:
            self.status, self.reason = 'done', 'eps_B'
        elif not (status == 'PENDING').any():
            self.status, self.reason = 'done', 'unresolved'

    def volume_ratio(self):
        """vol(O_x)/vol(P_x)-1；数值失败时为 inf（不按 ε_B 停止）。"""
        try:
            return polytope_volume(self.outer)/polytope_volume(self.vertices)-1.
        except geo.GEOMETRY_ERRORS:
            return np.inf

    def expanded(self):
        return geo.expanded_faces(self.faces, self.criterion, self.center, self.tau)

    def counts(self):
        status = list(self.face_status)
        return dict(faces=len(status), unresolved=status.count('UNRESOLVED'), unchecked=status.count('PENDING'),
                    boundary=status.count('BOUNDARY'), geometry=status.count('GEOMETRY_CERTIFIED'),
                    support=status.count('SUPPORT_CERTIFIED'))

    def view(self):
        return dict(scheme=self.label, vertices=self.vertices if self.full else np.empty((0, self.d)),
                    criterion=self.criterion, center=self.center, status=self.status, reason=self.reason)

    def summary(self):
        try:
            volume = polytope_volume(self.vertices) if self.full else 0.
            outer_volume = polytope_volume(self.outer) if len(self.outer) else 0.
        except geo.GEOMETRY_ERRORS:
            volume = outer_volume = None
        return dict(scheme=self.label, weight=self.weight, points=len(self.points), vertices=len(self.vertices),
                    full=self.full, criterion=self.criterion, center_source=self.center_source, strict=self.strict,
                    status=self.status, reason=self.reason, support_calls=self.calls,
                    support_seconds=self.seconds, violations_added=self.violations, rebuilds=self.rebuilds,
                    geometry_failures=self.geometry_failures,
                    inner_volume_xi=volume, outer_volume_xi=outer_volume,
                    **self.counts())


# ---- 方法 H：A → B → C ---------------------------------------------------------------------------
class SupportPhase:
    """方法 H 在 A 阶段之后的支撑查询（B）与覆盖证书（C），以及覆盖完成后续跑的径向部分。"""

    def __init__(self, case, code, radial, args, recorder, clock, deadline_t, pool, workers, coverage='global'):
        self.case, self.code, self.radial, self.args = case, code, radial, args
        self.coverage_mode, self.cone_certificates = coverage, set()
        self.recorder, self.clock, self.deadline, self.pool = recorder, clock, deadline_t, pool
        self.workers = workers if pool is not None else 1
        self.d, self.bounds, self.tau, self.eps = case.d, case.bounds, args.tau, args.eps
        self.networks, self.oracles, self.done = {}, {}, {}
        self.cover, self.coverage_rounds = None, []
        self.status, self.how, self.reason, self.t_cert = 'running', None, None, None
        self.phase_seconds, self.measure_cache = {}, {}

    # 状态与测度
    def state(self, override=None):
        return dict(started=True, phase=self.recorder.phase, cones=cone_rows(self.radial, override),
                    networks=[st.view() for st in self.networks.values() if st.full], cover=self.cover,
                    certified=self.status == 'certified')

    def measure_raw(self):
        nets = [st for st in self.networks.values() if st.full]
        if len(self.measure_cache) > 50000:
            self.measure_cache.clear()
        inner, outer = geo.h_measures(cone_rows(self.radial), [st.vertices for st in nets], self.cover, self.d,
                                      keys=[(st.label, st.rebuilds) for st in nets], cache=self.measure_cache)
        scale = float(np.prod(self.bounds))
        return inner*scale, outer*scale

    def measure(self):
        inner, outer = self.measure_raw()
        return inner, outer, self.status == 'certified', dict(
            networks=len(self.networks), cones=len(leaves(self.radial)), coverage=self.cover is not None,
            volume_ratio=self.radial.volume_ratio() if leaves(self.radial) else None)

    def unresolved(self):
        return sum(st.counts()['unresolved'] for st in self.networks.values() if not st.empty)

    def certified_now(self):
        """H 的分区认证：覆盖完成、无 UNRESOLVED 面且间隙 <= ε；或径向部分自身满足体积准则（继承 R 的证书）。"""
        if leaves(self.radial) and self.radial.volume_ratio() <= self.eps:
            return 'radial'
        if self.cover is None or self.unresolved():
            return None
        inner, outer = self.measure_raw()
        return 'support' if inner > 0 and outer-inner <= self.eps*inner else None

    def certify(self, how):
        self.status, self.how, self.t_cert = 'certified', how, self.clock.now()
        self.recorder.changed(force=True)

    @contextmanager
    def phase(self, name):
        self.recorder.phase = name
        started = self.clock.now()
        try:
            yield
        finally:
            self.phase_seconds[name] = self.phase_seconds.get(name, 0.)+self.clock.now()-started

    # 主流程
    def run(self):
        radial, recorder = self.radial, self.recorder
        if radial.status == 'error':
            self.status, self.reason = 'error', 'phase_A'
            return
        radial.deadline = self.clock.perf(self.deadline)
        recorder.state_fn, recorder.measure = self.state, self.measure
        try:
            # A 在根锥建立前用完时间片时，X* 与 V_x 只来自初始径向 MISOCP（现任解、零接入认证），O_R 为分区盒。
            how = self.certified_now()
            if how:
                self.certify(how)
                return
            with self.phase('B'):
                self.collect()
            for _ in range(MAX_ROUNDS):
                with self.phase('B'):
                    self.support_pass()
                if self.clock.now() >= self.deadline:
                    raise PartitionTimeout('partition time limit')
                with self.phase('C'):
                    answer = self.coverage()
                if answer['complete']:
                    recorder.before_change()
                    self.cover = answer['faces']   # 冻结覆盖证书所用的 E_x：之后 P_x 再变，外界仍取这一组
                    recorder.changed(force=True)
                    how = self.certified_now()
                    if how:
                        self.certify(how)
                        return
                    break
                if answer.get('witness') is None:
                    break   # 覆盖证书超时或失败：外界保持 O_R，续跑径向部分
                with self.phase('B'):
                    self.add_witness(*answer['witness'])
            self.resume_radial()
        except PartitionTimeout:
            if self.status != 'certified':
                self.status, self.reason = 'time_limit', 'deadline'
        except Exception as exc:
            self.status, self.reason = 'error', repr(exc)
            radial.errors.append(dict(kind='fatal', error=repr(exc), traceback=traceback.format_exc()))
            print(traceback.format_exc(), flush=True)

    def resume_radial(self):
        """覆盖完成而间隙仍 > ε，或覆盖证书不可用：续跑 A（目标 ε），每次锥决策后复核 H 的认证条件。"""
        radial = self.radial
        if self.clock.now() >= self.deadline:
            raise PartitionTimeout('partition time limit')
        radial.epsilon = self.eps
        radial.deadline = self.clock.perf(self.deadline)

        def check():
            how = self.certified_now()
            if how:
                self.certify(how)
                raise PartitionTimeout('H certified')
        radial.after_change = check
        with self.phase('A+'):
            radial.run_volume()
        radial.after_change = None
        if self.status == 'certified':
            return
        if radial.status == 'certified':
            self.certify('radial')
        elif radial.status == 'error':
            self.status, self.reason = 'error', 'radial'
        else:
            self.status, self.reason = ('time_limit', 'deadline') if radial.status == 'time_limit' else \
                (radial.status, 'radial')

    # A → B：收集 X* 与 V_x
    def collect(self):
        radial, d = self.radial, self.d
        volumes = {}
        for cone in leaves(radial):
            volumes[cone.x] = volumes.get(cone.x, 0.)+abs(float(np.linalg.det(cone.Q)))/factorial(d)
            if cone.cover is not None:
                volumes.setdefault(cone.cover, 0.)
        for scheme in [*radial.incumbent_schemes, *radial.equations.boxes, *radial.origin,
                       *(s for s, _ in radial.rays), *(s for s, _ in radial.nears)]:
            volumes.setdefault(tuple(scheme), 0.)
        total = sum(volumes.values())
        floor = .1*total/len(volumes) if total > 0. else 1.
        for scheme, volume in volumes.items():
            self.networks[scheme] = NetworkState(scheme, radial.label(scheme), self.known_points(scheme),
                                                 volume+floor, d, self.args, self.recorder)

    def known_points(self, scheme):
        """径向阶段已认证的该网架点：射线远端、近端、零接入点、审计过的已紧化现任解。"""
        radial = self.radial
        points = [p for (s, _), p in radial.rays.items() if s == scheme and p is not None]
        points += [p for (s, _), p in radial.nears.items() if s == scheme and p is not None]
        if radial.origin.get(scheme, {}).get('feasible'):
            points.append(np.zeros(self.d))
        points += [p for s, p in radial.incumbent_points if s == scheme]
        return np.asarray(points, float).reshape(-1, self.d)

    # B：逐网架支撑查询
    def solve_local(self, x, normal, limit):
        oracle = self.oracles.get(x)
        if oracle is None:
            oracle = self.oracles[x] = PartitionOracle(self.radial.equations, x, self.bounds, self.case.budget,
                                                       self.args.threads, limit, self.args.seed)
        oracle.time_limit = limit
        return slim(oracle.solve(normal))

    def submit(self, st, normal):
        self.radial.tighten(st.x)   # 已紧化时为空操作；新方案先 OBBT，支撑模型才带其包络行
        st.inflight.add(normal_key(normal))
        limit = max(min(SUPPORT_SECONDS, self.deadline-self.clock.now()), 0.)
        if self.pool is None:
            token = object()
            with self.recorder.log('SOCP', 'support', scheme=st.label) as row:
                entry = self.solve_local(st.x, normal, limit)
                row['status'] = STATUS_NAMES.get(entry['status'], str(entry['status']))
            self.done[token] = entry
            return token
        return self.pool.submit(st.x, self.radial.equations.boxes[st.x], normal, limit)

    def collect_one(self, inflight):
        for token in inflight:
            if token in self.done:
                return token, self.done.pop(token)
        finished, _ = wait(list(inflight), return_when=FIRST_COMPLETED)
        token = next(iter(finished))
        st, normal = inflight[token]
        try:
            entry = token.result()
        except Exception as exc:
            entry = dict(lb=-np.inf, ub=np.inf, point=None, status=-1, solve_seconds=0.,
                         normal=np.asarray(normal)/np.linalg.norm(normal), error=repr(exc)[:300])
        self.recorder.record('SOCP', 'support', float(entry.get('solve_seconds', 0.)), remote=True, scheme=st.label,
                             status=STATUS_NAMES.get(entry['status'], str(entry['status'])))
        return token, entry

    def support_pass(self):
        nets = [st for st in self.networks.values() if st.status in ('active', 'paused')]
        if not nets:
            return
        for st in nets:
            st.status, st.slice_end = 'active', None
            st.check_done()
        nets.sort(key=lambda st: (not st.urgent, -st.weight))
        start = self.clock.now()
        pass_end = start+SUPPORT_PASS_SHARE*max(self.deadline-start, 0.)
        inflight = {}
        while True:
            now = self.clock.now()
            if now >= pass_end:
                break
            free = self.workers-len(inflight)
            for st in nets:
                if free <= 0:
                    break
                if st.status != 'active':
                    continue
                if st.slice_end is None:
                    waiting = sum(s.weight for s in nets if s.status == 'active' and s.slice_end is None)
                    st.slice_end = now+(pass_end-now)*min(1., self.workers*st.weight/waiting)
                if now >= st.slice_end:
                    if not st.inflight:
                        st.status, st.reason = 'paused', 'slice'
                    continue
                for normal in st.next_tasks(free):
                    inflight[self.submit(st, normal)] = (st, normal)
                    free -= 1
            if not inflight:
                break
            self.apply(inflight)
        while inflight:   # 本轮结束：等待在途求解（各自受时限约束）
            self.apply(inflight)
        for st in nets:
            if st.status == 'active':
                st.status, st.reason = 'paused', 'pass_end'
            st.urgent = False

    def apply(self, inflight):
        token, entry = self.collect_one(inflight)
        st, normal = inflight.pop(token)
        self.recorder.before_change()
        st.on_result(normal, entry)
        st.check_done()
        self.recorder.changed()

    # C：覆盖证书
    def coverage(self):
        """C：x 自由的覆盖证书。global：一个 RemainingRegionModel（全部 E_x）；cone（方法 Hc）：按 A 阶段叶锥分解，
        每锥只保留切过 锥∩盒 的面（按区域顶点逐面判定，精确），某个 E_x 包含整个 锥∩盒 时该锥无需求解；锥证书在
        相关网架未变时复用（紧化只缩小可行集、新网架只扩大覆盖，旧证书仍成立）。"""
        nets = [st for st in self.networks.values() if st.full]
        faces = [st.expanded() for st in nets]
        if not faces:
            return dict(complete=False, witness=None, status='no_inner')
        if self.coverage_mode == 'cone':
            answer = self.coverage_by_cone(nets, faces)
        else:
            answer = self.coverage_solve(faces, None, dict(scope='global'))
        if answer['complete']:
            answer['faces'] = faces
        return answer

    def coverage_by_cone(self, nets, faces):
        cones = cone_rows(self.radial)
        regions = [geo.clip_box(geo.cone_outer(c['U'], None)) for c in cones]
        order = sorted(range(len(cones)),
                       key=lambda k: -polytope_volume(regions[k]) if len(regions[k]) > self.d else 0.)
        bounds = []
        for k in order:
            cone, region = cones[k], regions[k]
            kept, used, covered = [], [], False
            for st, E in zip(nets, faces):
                whole, disjoint, cut = geo.cone_faces(region, E)
                if whole:
                    covered = True
                    break
                if not disjoint:
                    kept.append(cut)
                    used.append((st.label, st.rebuilds))
            if covered:
                continue
            key = (np.asarray(cone['U'], float).tobytes(), frozenset(used))
            if key in self.cone_certificates:
                continue
            rows = np.linalg.inv(np.asarray(cone['U'], float))
            answer = self.coverage_solve(kept, rows/np.linalg.norm(rows, axis=1, keepdims=True),
                                         dict(scope='cone', cone=k, cones=len(cones)))
            if not answer['complete']:
                return answer
            self.cone_certificates.add(key)
            bounds.append(answer['bound'])
        finite = [b for b in bounds if b is not None]
        return dict(complete=True, witness=None, status='complete', bound=max(finite) if finite else None)

    def coverage_solve(self, faces, rows, info):
        """一个覆盖 MISOCP：RemainingRegionModel 直接接收 E_x 面并令 tau=0（s=1），加全部已紧化方案的提升行；
        rows 为锥约束 rows@xi>=0（None 为整个分区）。见证方案未紧化时懒惰 OBBT 后重解。"""
        radial, equations = self.radial, self.radial.equations
        for _ in range(len(equations.keys)+len(self.networks)+8):
            limit = min(self.args.mip_seconds, self.deadline-self.clock.now())
            if limit <= 0.:
                raise PartitionTimeout('partition time limit')
            answer, witness, bound = None, None, None
            row_info = dict(networks=len(faces), faces=int(sum(len(f) for f in faces)),
                            tightened=len(equations.boxes), **info)
            with self.recorder.log('MISOCP', 'coverage', **row_info) as row:
                problem = RemainingRegionModel(equations, self.case.budget, self.bounds, float(self.bounds.sum()), [],
                                               faces, 0., axis_bounds=self.bounds, threads=self.args.threads)
                model = problem.model
                try:
                    choices = dict(zip(equations.keys, problem.problem.x.tolist()))
                    for scheme in equations.boxes:
                        equations._tighten(model, choices, problem.problem.operation, scheme)
                    if rows is not None:
                        power = [problem.problem.power[j].item() for j in range(self.d)]
                        for j, cone_row in enumerate(rows):
                            model.addConstr(gp.quicksum(float(cone_row[i]/self.bounds[i])*power[i] for i in range(self.d))
                                            >= 0., name=f'coverage_cone_{j}')
                    model.Params.BarQCPConvTol = COVER_CONV_TOL
                    model.Params.Seed = self.args.seed
                    try:
                        answer = problem.solve(GEOMETRY_TOL, time_limit=limit)
                        row['status'] = 'complete' if answer['complete'] else 'witness'
                    except TimeoutError:
                        row['status'] = 'TIME_LIMIT'
                    except RuntimeError as exc:
                        row.update(status='ERROR', error=repr(exc)[:300])
                    try:
                        bound = float(model.ObjBound)/problem.distance_scale
                    except (gp.GurobiError, AttributeError):
                        bound = None if answer is None else answer['bound']
                    if answer is not None and not answer['complete']:
                        x = tuple(int(v) for v in answer['x'])
                        residual = audit_incumbent(problem.problem, answer['x'])
                        witness = (x, np.asarray(answer['p'], float)/self.bounds, residual)
                        row.update(scheme=radial.label(x), audit_residual=residual)
                    row['bound'] = bound
                finally:
                    model.dispose()
            self.coverage_rounds.append(dict(t=self.clock.now(), status=row['status'], bound=bound, **row_info,
                                             scheme=None if witness is None else radial.label(witness[0])))
            if answer is None:
                return dict(complete=False, witness=None, status=row['status'], bound=bound)
            if answer['complete']:
                return dict(complete=True, witness=None, status='complete', bound=answer['bound'])
            if witness[0] not in equations.boxes:
                radial.tighten(witness[0])   # 懒惰 OBBT：见证方案加入提升行后重解
                continue
            return dict(complete=False, witness=witness, status='witness', bound=bound)
        return dict(complete=False, witness=None, status='lazy_obbt_limit', bound=None)

    def add_witness(self, x, xi, residual):
        """覆盖见证补进 X*/V_x：审计通过直接用；否则从内点锚出发朝见证点做紧化射线取认证点。"""
        st = self.networks.get(x)
        point = xi if residual <= PLANNING_TOL else self.certify_witness(x, xi, st)
        if st is None:
            points = self.known_points(x)
            if point is not None:
                points = np.vstack([points, point])
            floor = min((s.weight for s in self.networks.values()), default=1.)
            st = self.networks[x] = NetworkState(x, self.radial.label(x), points, floor, self.d, self.args,
                                                 self.recorder)
        elif point is not None:
            self.recorder.before_change()
            st.add_points([point])
            self.recorder.changed()
        # 重新做 B：见证说明按 ε_B 停下的 P_x 的外扩不足以覆盖 R_x，故该网架此后查完全部面才停（时间片仍有效）。
        st.status, st.urgent, st.strict = 'active', True, True

    def certify_witness(self, x, xi, st):
        radial, equations = self.radial, self.radial.equations
        if st is not None and st.full:
            anchor = st.center*self.bounds
        elif radial.origin_feasible(x):
            anchor = np.zeros(self.d)
        else:
            with self.recorder.log('SOCP', 'ray', scheme=radial.label(x), ray='center'):
                anchor = equations.center(np.asarray(x), self.bounds)
        limit = max(min(RAY_SECONDS, self.deadline-self.clock.now()), 0.)
        with self.recorder.log('SOCP', 'ray', scheme=radial.label(x), ray='witness') as row:
            try:
                answer = ray_support(equations, self.case.budget, np.asarray(x), anchor, xi*self.bounds,
                                     threads=self.args.threads, time_limit=limit)
            except (RuntimeError, TimeoutError) as exc:
                row.update(status='failed', error=repr(exc)[:300])
                return None
        return answer['p']/self.bounds

    def close(self):
        for oracle in self.oracles.values():
            oracle.close()
        self.oracles = {}

    def result(self, phase_a, start, end):
        nets = list(self.networks.values())
        counts = [st.counts() for st in nets if not st.empty]
        inner, outer = self.measure_raw()
        return dict(partition=self.code, method='H', status=self.status, certified=self.status == 'certified',
                    how=self.how, reason=self.reason, t_cert=self.t_cert, start=start, end=end, seconds=end-start,
                    inner_measure=inner, outer_measure=outer, gap=outer/inner-1. if inner > 0 else None,
                    volume_ratio=self.radial.volume_ratio() if leaves(self.radial) else None,
                    cones=len(leaves(self.radial)), x_star=len(nets),
                    faces=sum(c['faces'] for c in counts), unresolved_faces=sum(c['unresolved'] for c in counts),
                    unchecked_faces=sum(c['unchecked'] for c in counts),
                    coverage=dict(complete=self.cover is not None, rounds=len(self.coverage_rounds),
                                  history=self.coverage_rounds),
                    phase_a=phase_a, phase_seconds=dict(A=phase_a['seconds'], **self.phase_seconds),
                    networks=[st.summary() for st in nets], counts=dict(self.radial.counts),
                    errors=self.radial.errors[-5:])


STATUS_NAMES = {getattr(GRB.Status, name): name for name in dir(GRB.Status) if name.isupper()}
STATUS_NAMES[-1] = 'ERROR'


def run_H(case, code, args, clock, deadline_t, workers, method='H'):
    """方法 H（coverage=args.coverage）或 Hc（同 H，覆盖证书按锥分解）的一个分区。"""
    recorder = Recorder(code, method, clock, checkpoint_times(args.seconds),
                        TIMELINE_INTERVAL_3D if case.d == 3 else 0.)
    pool = SupportPool(case, code, args, workers) if workers > 1 else None
    start = clock.now()
    deadline_a = min(deadline_t, start+args.discovery_share*max(deadline_t-start, 0.))
    radial = make_radial(case, code, args, workers, args.discovery_eps, deadline_a, recorder, clock)
    radial.audit = True
    recorder.state_fn = lambda override=None: radial_state(radial, 'A', override)
    recorder.measure = lambda: (*radial_measure(radial), False,
                                dict(volume_ratio=radial.volume_ratio() if leaves(radial) else None,
                                     cones=len(leaves(radial))))
    support = SupportPhase(case, code, radial, args, recorder, clock, deadline_t, pool, workers,
                           'cone' if method == 'Hc' else args.coverage)
    phase_a = dict(status='not_started', seconds=0.)
    try:
        with threadpool_limits(limits=1):
            radial.run(None)
            phase_a = dict(status=radial.status, seconds=clock.now()-start, cones=len(leaves(radial)),
                           volume_ratio=radial.volume_ratio() if leaves(radial) else None,
                           incumbents_audited=len(radial.incumbent_points))
            support.run()
    finally:
        support.close()
        if pool is not None:
            pool.close()
    end = clock.now()
    recorder.finish()
    return support.result(phase_a, start, end), recorder


# ---- 一次运行：分区调度、合并、离线评价 -------------------------------------------------------------
def run_partition(method, case, code, args, clock, deadline_t, workers):
    try:
        if method == 'R':
            return run_R(case, code, args, clock, deadline_t, workers)
        return run_H(case, code, args, clock, deadline_t, workers, method)
    except Exception as exc:
        print(traceback.format_exc(), flush=True)
        recorder = Recorder(code, method, clock, checkpoint_times(args.seconds))
        recorder.finish()
        return dict(partition=code, method=method, status='error', certified=False, t_cert=None, start=None,
                    end=clock.now(), seconds=None, error=repr(exc), traceback=traceback.format_exc()), recorder


def partition_payload(result, recorder):
    return plain(dict(result=result, solves=recorder.solves, timeline=recorder.timeline,
                      snapshots={('final' if k == 'final' else f'{float(k):g}'): v for k, v in recorder.snapshots.items()}))


def run_once(args, case, method, workers, repeat, directory):
    """一次运行（方法 × 线程设置 × 重复）：分区串行顺延或并行同时开始，墙钟总时限 args.seconds。"""
    directory.mkdir(parents=True, exist_ok=True)
    codes = args.partitions or case.partitions
    epoch = time.time()
    clock = Clock(epoch)
    concurrency = min(workers, len(codes))
    inner_workers = max(1, workers//len(codes)) if workers >= len(codes) else 1
    print(f'== {case.tag} method {method} workers {workers} run {repeat}: {len(codes)} partitions, '
          f'concurrency {concurrency}, per-partition workers {inner_workers}', flush=True)
    payloads = {}
    if concurrency == 1:
        for j, code in enumerate(codes):
            now = clock.now()
            deadline_t = now+max(args.seconds-now, 0.)/(len(codes)-j)
            print(f'  partition {code}: start {now:.1f}s, deadline {deadline_t:.1f}s', flush=True)
            result, recorder = run_partition(method, case, code, args, clock, deadline_t, inner_workers)
            payloads[code] = partition_payload(result, recorder)
            print(f'  partition {code}: {result["status"]} at {clock.now():.1f}s', flush=True)
    else:
        queue, running = list(codes), []
        work = directory/'partitions'
        work.mkdir(parents=True, exist_ok=True)
        while queue or running:
            while queue and len(running) < concurrency:
                code = queue.pop(0)
                now = clock.now()
                remaining = len(queue)+len(running)+1
                deadline_t = now+max(args.seconds-now, 0.)*min(1., concurrency/remaining)
                spec = dict(case=case.name, nodes=list(case.nodes), code=code, method=method, epoch=epoch,
                            deadline=deadline_t, workers=inner_workers, args=vars(args), output=str(work/code))
                (work/f'{code}.spec.json').write_text(json.dumps(plain(spec), default=str), encoding='utf-8')
                log = (work/f'{code}.log').open('w', encoding='utf-8')
                process = subprocess.Popen([sys.executable, '-X', 'utf8', str(Path(__file__).resolve()), '--worker',
                                            str(work/f'{code}.spec.json')], stdout=log, stderr=subprocess.STDOUT,
                                           cwd=ROOT)
                running.append((code, process, log))
            for item in list(running):
                code, process, log = item
                if process.poll() is not None:
                    log.close()
                    running.remove(item)
                    path = work/code/'partition.json'
                    if path.exists():
                        payloads[code] = json.loads(path.read_text(encoding='utf-8'))
                    else:
                        payloads[code] = dict(result=dict(partition=code, method=method, status='error',
                                                          certified=False, error=f'worker exit {process.returncode}'),
                                              solves=[], timeline=[], snapshots={})
                    print(f'  partition {code}: {payloads[code]["result"]["status"]} (exit {process.returncode}) '
                          f'at {clock.now():.1f}s', flush=True)
            time.sleep(.2)
    wall = clock.now()
    return write_run(directory, args, case, method, workers, repeat, codes, payloads, wall, inner_workers)


def run_timeline(rows, codes, box):
    """全部分区之和的时间线：未开始的分区内域为空、外界为分区盒。"""
    current = {code: (0., box, False) for code in codes}
    merged = []
    for row in sorted(rows, key=lambda r: r['t']):
        current[row['partition']] = (row['inner_measure'], row['outer_measure'], bool(row['certified']))
        inner = sum(v[0] for v in current.values())
        outer = sum(v[1] for v in current.values())
        merged.append(dict(t=row['t'], partition='all', method=row['method'], phase=row['phase'],
                           inner_measure=inner, outer_measure=outer,
                           gap=outer/inner-1. if inner > 0 else None,
                           certified=sum(v[2] for v in current.values())))
    return merged


def write_csv(path, rows):
    columns = list(dict.fromkeys(key for row in rows for key in row))
    with path.open('w', encoding='utf-8-sig', newline='') as stream:
        writer = csv.DictWriter(stream, columns)
        writer.writeheader()
        writer.writerows(rows)


def write_run(directory, args, case, method, workers, repeat, codes, payloads, wall, inner_workers):
    results = {code: payloads[code]['result'] for code in codes}
    solves = [row for code in codes for row in payloads[code]['solves']]
    rows = [row for code in codes for row in payloads[code]['timeline']]
    box = float(np.prod(case.bounds))
    merged = run_timeline(rows, codes, box)
    write_csv(directory/'solves.csv', solves or [dict(partition=None)])
    write_csv(directory/'timeline.csv', [*rows, *merged] or [dict(t=None)])
    snapshots = directory/'snapshots'
    snapshots.mkdir(exist_ok=True)
    keys = sorted({key for code in codes for key in payloads[code]['snapshots']},
                  key=lambda k: np.inf if k == 'final' else float(k))
    for key in keys:
        name = 'final' if key == 'final' else f't{float(key):05.0f}'
        content = {code: payloads[code]['snapshots'].get(key, empty_state()) for code in codes}
        (snapshots/f'{name}.json').write_text(json.dumps(plain(content)), encoding='utf-8')
    certified = all(results[code].get('certified') for code in codes)
    t_cert = max(results[code]['t_cert'] for code in codes) if certified else None
    gaps = {f'{level:g}': next((row['t'] for row in merged if row['gap'] is not None and row['gap'] <= level), None)
            for level in GAP_LEVELS}
    summary = plain(dict(
        case=case.tag, case_name=case.name, network=case.network.name, load_nodes=list(case.nodes), mode=case.mode,
        budget=case.budget,
        bounds_kw=case.bounds, method=method, workers=workers, per_partition_workers=inner_workers, repeat=repeat,
        partitions=list(codes), settings=settings(args, case), commit=git_commit(), dirty=git_dirty(),
        status='certified' if certified else ('error' if any(r.get('status') == 'error' for r in results.values())
                                              else 'time_limit'),
        certified=certified, certified_partitions=sum(bool(r.get('certified')) for r in results.values()),
        t_cert=t_cert, wall_seconds=wall, gap_final=merged[-1]['gap'] if merged else None, t_gap=gaps,
        totals=solve_totals(solves, results), partition_results=results))
    (directory/'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=1), encoding='utf-8')
    return summary


def git_dirty():
    """已跟踪文件是否有未提交的修改（结果对应的代码不止 commit 时为 True）。"""
    try:
        return bool(subprocess.run(['git', 'status', '--porcelain', '--untracked-files=no'], cwd=ROOT,
                                   capture_output=True, text=True, check=True).stdout.strip())
    except (OSError, subprocess.CalledProcessError):
        return None


def settings(args, case):
    return dict(seconds=args.seconds, tau=args.tau, eps=args.eps, discovery_eps=args.discovery_eps,
                discovery_share=args.discovery_share, network_eps=args.network_eps, criterion=args.criterion,
                coverage=args.coverage, threads=args.threads, seed=args.seed, mip_seconds=args.mip_seconds,
                mip_gap=args.mip_gap, inner_cert=args.inner_cert, support_seconds=SUPPORT_SECONDS,
                support_pass_share=SUPPORT_PASS_SHARE, checkpoints=checkpoint_times(args.seconds),
                geometry_tol=GEOMETRY_TOL, planning_tol=PLANNING_TOL, cover_conv_tol=COVER_CONV_TOL,
                gurobi=gp.gurobi.version())


def solve_totals(solves, results):
    """求解次数与耗时：按类型（MISOCP/SOCP/LP）与用途；other=各分区墙钟减去本进程求解耗时。"""
    by_type, by_purpose = {}, {}
    for row in solves:
        for table, key in ((by_type, row['type']), (by_purpose, row['purpose'])):
            entry = table.setdefault(key, dict(count=0, seconds=0.))
            entry['count'] += 1
            entry['seconds'] += float(row['seconds'])
    main_seconds = sum(float(row['seconds']) for row in solves if row.get('process') == 'main')
    wall = sum(float(r.get('seconds') or 0.) for r in results.values())
    phases = {}
    for r in results.values():
        for name, seconds in (r.get('phase_seconds') or {}).items():
            phases[name] = phases.get(name, 0.)+seconds
    return dict(by_type=by_type, by_purpose=by_purpose, partition_wall_seconds=wall,
                main_solve_seconds=main_seconds, other_seconds=max(wall-main_seconds, 0.),
                pool_solve_seconds=sum(float(row['seconds']) for row in solves if row.get('process') == 'pool'),
                phase_seconds=phases,
                x_star=sum(r.get('x_star', 0) or 0 for r in results.values()),
                cones=sum(r.get('cones', 0) or 0 for r in results.values()),
                faces=sum(r.get('faces', 0) or 0 for r in results.values()),
                unresolved_faces=sum(r.get('unresolved_faces', 0) or 0 for r in results.values()))


def worker_main(spec_path):
    """并行运行的分区子进程：读规格，跑一个分区，写 partition.json。"""
    spec = json.loads(Path(spec_path).read_text(encoding='utf-8'))
    args = Namespace(**spec['args'])
    case = Case(spec['case'], spec['nodes'])
    clock = Clock(spec['epoch'])
    output = Path(spec['output'])
    output.mkdir(parents=True, exist_ok=True)
    result, recorder = run_partition(spec['method'], case, spec['code'], args, clock, spec['deadline'], spec['workers'])
    (output/'partition.json').write_text(json.dumps(partition_payload(result, recorder)), encoding='utf-8')
    print(f'partition {spec["code"]} {spec["method"]}: {result["status"]}', flush=True)
    return 0


# ---- 离线逐格评价 ------------------------------------------------------------------------------
def cell_labels(state, xi):
    """一个分区快照的逐格标签：内域 = 径向内域 ∪ 各 P_x；外界 = 径向外界 [∩ 覆盖冻结的 ∪E_x]；未开始为 空/盒。"""
    if not state.get('started') or not state['cones']:
        return np.zeros(len(xi), bool), np.ones(len(xi), bool)
    cones = [dict(U=np.asarray(c['U'], float), V=np.asarray(c['V'], float), pseudo=bool(c.get('pseudo')),
                  halfspace=None if c['halfspace'] is None else (np.asarray(c['halfspace'][0], float),
                                                                  float(c['halfspace'][1])))
             for c in state['cones']]
    inner, outer = geo.classify_cones(xi, cones)
    for network in state['networks']:
        inner |= geo.network_cells(xi, np.asarray(network['vertices'], float).reshape(-1, xi.shape[1]))
    if state['cover'] is not None:
        covered = np.zeros(len(xi), bool)
        for faces in state['cover']:
            covered |= contains(xi, np.asarray(faces, float), GEOMETRY_TOL)
        outer &= covered
    return inner, outer


def partition_cells(case, scan):
    cells = {}
    for code in case.partitions:
        sign = np.asarray(partition_sign(code, case.d))
        index = np.flatnonzero(np.all(np.where(scan['power'] >= 0., 1, -1) == sign, axis=1))
        cells[code] = (index, np.abs(scan['power'][index])/case.bounds)
    return cells


def evaluate_run(directory, case, scan, cells=None):
    """各检查点与 final 的逐格指标（AC 为主参考，SOCP 为诊断）及有效性检查，写 metrics.json。"""
    cells = cells or partition_cells(case, scan)
    summary = json.loads((directory/'summary.json').read_text(encoding='utf-8'))
    codes = summary['partitions']
    report = dict(scan=scan['path'].resolve().relative_to(ROOT).as_posix(), shape=list(scan['shape']), checkpoints={})
    for path in sorted((directory/'snapshots').glob('*.json')):
        snapshot = json.loads(path.read_text(encoding='utf-8'))
        labels, per = [], {}
        for code in codes:
            index, xi = cells[code]
            inner, outer = cell_labels(snapshot.get(code, empty_state()), xi)
            per[code] = geo.grid_metrics(inner, outer, scan['ac'][index], scan['socp'][index])
            labels.append((index, inner, outer))
        index = np.concatenate([row[0] for row in labels])
        inner = np.concatenate([row[1] for row in labels])
        outer = np.concatenate([row[2] for row in labels])
        overall = geo.grid_metrics(inner, outer, scan['ac'][index], scan['socp'][index])
        if path.stem == 'final':
            missed = index[(scan['ac'][index] == 1) & ~outer]
            wrong = index[inner & (scan['socp'][index] == -1)]
            overall['outer_missed_ac_kw'] = scan['power'][missed[:50]]
            overall['inner_socp_infeasible_kw'] = scan['power'][wrong[:50]]
        key = 'final' if path.stem == 'final' else f'{float(path.stem[1:]):g}'
        report['checkpoints'][key] = dict(overall=overall, partitions=per)
    report['valid'] = all(row['overall']['valid'] for row in report['checkpoints'].values())
    (directory/'metrics.json').write_text(json.dumps(plain(report), indent=1), encoding='utf-8')
    return report


# ---- 汇总报告 ----------------------------------------------------------------------------------
def run_rows(output):
    rows = []
    for path in sorted(output.glob('*/workers_*/run_*/summary.json')):
        summary = json.loads(path.read_text(encoding='utf-8'))
        metrics_path = path.parent/'metrics.json'
        metrics = json.loads(metrics_path.read_text(encoding='utf-8')) if metrics_path.exists() else None
        rows.append(dict(directory=path.parent, summary=summary, metrics=metrics))
    rows.sort(key=lambda r: (r['summary']['method'], r['summary']['workers'], r['summary']['repeat']))
    return rows


def flat_row(row):
    s, m = row['summary'], row['metrics']
    final = None if m is None else m['checkpoints'].get('final', {}).get('overall')
    pct = lambda key, kind: None if final is None else final[key][f'{kind}_percent']
    totals = s['totals']
    out = dict(method=s['method'], workers=s['workers'], repeat=s['repeat'], status=s['status'],
               certified_partitions=s['certified_partitions'], partitions=len(s['partitions']), t_cert=s['t_cert'],
               wall_seconds=s['wall_seconds'], gap_final=s['gap_final'],
               **{f't_gap_{k}': v for k, v in s['t_gap'].items()},
               inner_ac_fr=pct('inner_ac', 'fr'), inner_ac_mr=pct('inner_ac', 'mr'),
               inner_socp_fr=pct('inner_socp', 'fr'), inner_socp_mr=pct('inner_socp', 'mr'),
               outer_ac_fr=pct('outer_ac', 'fr'), outer_ac_mr=pct('outer_ac', 'mr'),
               outer_socp_fr=pct('outer_socp', 'fr'), outer_socp_mr=pct('outer_socp', 'mr'),
               valid=None if m is None else m['valid'],
               inner_socp_infeasible_cells=None if final is None else final['inner_socp_infeasible_cells'],
               outer_missed_ac_cells=None if final is None else final['outer_missed_ac_cells'],
               outer_missed_socp_cells=None if final is None else final['outer_missed_socp_cells'],
               undecided_ac_cells=None if final is None else final['undecided_ac_cells'],
               x_star=totals['x_star'], cones=totals['cones'], faces=totals['faces'],
               unresolved_faces=totals['unresolved_faces'], other_seconds=totals['other_seconds'],
               pool_solve_seconds=totals['pool_solve_seconds'])
    for kind in ('MISOCP', 'SOCP', 'LP'):
        entry = totals['by_type'].get(kind, dict(count=0, seconds=0.))
        out[f'{kind.lower()}_count'], out[f'{kind.lower()}_seconds'] = entry['count'], entry['seconds']
    for purpose in PURPOSES:
        entry = totals['by_purpose'].get(purpose, dict(count=0, seconds=0.))
        out[f'{purpose}_count'], out[f'{purpose}_seconds'] = entry['count'], entry['seconds']
    for phase in ('R', 'A', 'B', 'C', 'A+'):
        out[f'phase_{phase}_seconds'] = totals['phase_seconds'].get(phase)
    return out


def median_range(values):
    values = [v for v in values if v is not None]
    if not values:
        return '-'
    low, mid, high = min(values), float(np.median(values)), max(values)
    fmt = (lambda v: f'{v:.3g}') if max(abs(low), abs(high)) < 1000 else (lambda v: f'{v:.0f}')
    return fmt(mid) if low == high else f'{fmt(mid)} [{fmt(low)}, {fmt(high)}]'


def cell(value, digits=3):
    if value is None:
        return '-'
    if isinstance(value, bool):
        return 'yes' if value else 'no'
    if isinstance(value, float):
        return f'{value:.{digits}f}'
    return str(value)


def write_report(output, case_tag):
    rows = run_rows(output)
    if not rows:
        print(f'no runs under {output}', flush=True)
        return 1
    flat = [flat_row(r) for r in rows]
    write_csv(output/'comparison.csv', flat)
    groups = {}
    for row in flat:
        groups.setdefault((row['method'], row['workers']), []).append(row)
    first = rows[0]['summary']
    s = first['settings']
    lines = [f'# Plan 2 comparison: {case_tag}', '',
             f'Commit `{first.get("commit")}`; {first["network"]} load nodes {first["load_nodes"]}, mode {first["mode"]}, '
             f'budget {first["budget"]:g}; partitions {", ".join(first["partitions"])}.', '',
             f'Settings: total {s["seconds"]:g} s wall clock per run; tau {s["tau"]:g}; eps {s["eps"]:g}; '
             f'H: eps_A {s["discovery_eps"]:g}, share_A {s["discovery_share"]:g}, eps_B {s["network_eps"]:g}, '
             f'criterion {s["criterion"]}, coverage {s["coverage"]}; Gurobi threads {s["threads"]}, seed {s["seed"]}, '
             f'MISOCP cap {s["mip_seconds"]:g} s, MIPGap {s["mip_gap"]:g}; Gurobi {".".join(map(str, s["gurobi"]))}.', '',
             'Primary reference: AC scan (FR/MR in %, undecided AC cells excluded). SOCP columns are diagnostics: '
             'the OBBT-tightened model removes SOCP-only area, so outer-vs-SOCP MR > 0 is expected.', '']
    invalid = [r for r in flat if r['valid'] is False]
    if invalid:
        lines += ['**<span style="color:red">VALIDITY FAILED</span>** in: ' +
                  ', '.join(f'{r["method"]}/w{r["workers"]}/run{r["repeat"]}' for r in invalid), '']
    lines += ['## Runs', '', '| method | workers | run | status | certified partitions | t_cert (s) | final gap | '
              't gap<=10% | t gap<=5% | t gap<=2% | valid |', '|---|---|---|---|---|---|---|---|---|---|---|']
    for r in flat:
        mark = 'yes' if r['valid'] else ('**<span style="color:red">NO</span>**' if r['valid'] is False else '-')
        lines.append(f'| {r["method"]} | {r["workers"]} | {r["repeat"]} | {r["status"]} | '
                     f'{r["certified_partitions"]}/{r["partitions"]} | {cell(r["t_cert"], 1)} | {cell(r["gap_final"], 4)} | '
                     f'{cell(r["t_gap_0.1"], 1)} | {cell(r["t_gap_0.05"], 1)} | {cell(r["t_gap_0.02"], 1)} | {mark} |')
    keys = [('t_cert', 't_cert (s)'), ('gap_final', 'final gap'), ('t_gap_0.1', 't gap<=10%'),
            ('t_gap_0.05', 't gap<=5%'), ('t_gap_0.02', 't gap<=2%'),
            ('inner_ac_fr', 'inner-AC FR%'), ('inner_ac_mr', 'inner-AC MR%'),
            ('outer_ac_fr', 'outer-AC FR%'), ('outer_ac_mr', 'outer-AC MR%'),
            ('inner_socp_fr', 'inner-SOCP FR%'), ('inner_socp_mr', 'inner-SOCP MR%'),
            ('outer_socp_fr', 'outer-SOCP FR%'), ('outer_socp_mr', 'outer-SOCP MR%')]
    lines += ['', '## Medians [min, max] over repeats', '',
              '| method | workers | runs certified | ' + ' | '.join(k[1] for k in keys) + ' |',
              '|---|---|---|' + '---|'*len(keys)]
    for (method, workers), items in sorted(groups.items()):
        lines.append(f'| {method} | {workers} | {sum(r["status"] == "certified" for r in items)}/{len(items)} | ' +
                     ' | '.join(median_range([r[k] for r in items]) for k, _ in keys) + ' |')
    lines += ['', '## Validity (final state, cells)', '',
              '| method | workers | run | inner ∧ SOCP-infeasible | AC-feasible ∖ outer | SOCP-feasible ∖ outer '
              '(diagnostic) | undecided AC |', '|---|---|---|---|---|---|---|']
    for r in flat:
        red = lambda v: f'**<span style="color:red">{v}</span>**' if v else str(v)
        lines.append(f'| {r["method"]} | {r["workers"]} | {r["repeat"]} | {red(r["inner_socp_infeasible_cells"])} | '
                     f'{red(r["outer_missed_ac_cells"])} | {r["outer_missed_socp_cells"]} | {r["undecided_ac_cells"]} |')
    time_keys = [('misocp_seconds', 'MISOCP s'), ('socp_seconds', 'SOCP s'), ('lp_seconds', 'LP s'),
                 ('other_seconds', 'other s'), *((f'{p}_seconds', f'{p} s') for p in PURPOSES),
                 ('pool_solve_seconds', 'pool solve s'), ('phase_A_seconds', 'phase A s'), ('phase_B_seconds', 'phase B s'),
                 ('phase_C_seconds', 'phase C s'), ('phase_A+_seconds', 'phase A+ s')]
    lines += ['', '## Time breakdown (medians; solver seconds summed over partitions and processes)', '',
              '| method | workers | ' + ' | '.join(k[1] for k in time_keys) + ' |', '|---|---|' + '---|'*len(time_keys)]
    for (method, workers), items in sorted(groups.items()):
        lines.append(f'| {method} | {workers} | ' + ' | '.join(median_range([r[k] for r in items]) for k, _ in time_keys)
                     + ' |')
    count_keys = [('misocp_count', 'MISOCP'), ('socp_count', 'SOCP'), ('lp_count', 'LP'), ('x_star', '|X*|'),
                  ('cones', 'cones'), ('faces', 'faces'), ('unresolved_faces', 'unresolved faces')]
    lines += ['', '## Counts (medians)', '', '| method | workers | ' + ' | '.join(k[1] for k in count_keys) + ' |',
              '|---|---|' + '---|'*len(count_keys)]
    for (method, workers), items in sorted(groups.items()):
        lines.append(f'| {method} | {workers} | ' + ' | '.join(median_range([r[k] for r in items]) for k, _ in count_keys)
                     + ' |')
    lines += ['', '## Partitions (status, t_cert s)', '']
    for row in rows:
        summary = row['summary']
        parts = '; '.join(f'{code} {p.get("status")}{"" if p.get("t_cert") is None else f" {p["t_cert"]:.1f}"}'
                          + (f' ({p["how"]})' if p.get('how') else '')
                          for code, p in summary['partition_results'].items())
        lines.append(f'- {summary["method"]} w{summary["workers"]} run{summary["repeat"]}: {parts}')
    lines += ['', '## Conclusions', '', *conclusions(groups), '']
    (output/'comparison.md').write_text('\n'.join(lines), encoding='utf-8')
    try:
        plot_report(output, rows)
    except Exception:
        print(traceback.format_exc(), flush=True)
    return 0 if not invalid else 3


def conclusions(groups):
    lines = []
    for workers, method in [(w, m) for w in sorted({w for _, w in groups}) for m in ('H', 'Hc')]:
        r, h = groups.get(('R', workers), []), groups.get((method, workers), [])
        if not r or not h:
            continue
        cert = lambda items: sum(i['status'] == 'certified' for i in items)
        med = lambda items, key: (float(np.median([i[key] for i in items if i[key] is not None]))
                                  if any(i[key] is not None for i in items) else None)
        m = method
        line = (f'- workers {workers}: R certified {cert(r)}/{len(r)} runs, {m} {cert(h)}/{len(h)}; median final gap '
                f'R {cell(med(r, "gap_final"), 4)} vs {m} {cell(med(h, "gap_final"), 4)}; median t_cert '
                f'R {cell(med(r, "t_cert"), 1)} s vs {m} {cell(med(h, "t_cert"), 1)} s; inner-AC MR '
                f'R {cell(med(r, "inner_ac_mr"))}% vs {m} {cell(med(h, "inner_ac_mr"))}%, inner-AC FR '
                f'R {cell(med(r, "inner_ac_fr"))}% vs {m} {cell(med(h, "inner_ac_fr"))}%.')
        phases = {p: med(h, f'phase_{p}_seconds') for p in ('A', 'B', 'C', 'A+')}
        line += (f' {m} phase medians (s): '
                 + ', '.join(f'{p} {cell(v, 1)}' for p, v in phases.items() if v is not None) + '.')
        lines.append(line)
    return lines or ['- (needs both methods for the same workers setting)']


# ---- 图 -----------------------------------------------------------------------------------------
def _plt():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 9, 'axes.edgecolor': '#c3c2b7',
                         'axes.labelcolor': INK_2, 'xtick.color': INK_2, 'ytick.color': INK_2})
    return plt


def _frame(ax, title, xlabel, ylabel=None, log=False):
    ax.set_facecolor(SURFACE)
    if log:
        ax.set_yscale('log')
    ax.set_title(title, color=INK, fontsize=10, loc='left')
    ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    ax.grid(color=GRID, lw=.5)
    ax.set_axisbelow(True)


def read_timeline(directory):
    with (directory/'timeline.csv').open(encoding='utf-8-sig') as stream:
        rows = [r for r in csv.DictReader(stream) if r.get('partition') == 'all']
    return [(float(r['t']), float(r['gap'])) for r in rows if r['gap'] not in ('', 'None', 'inf')]


def plot_report(output, rows):
    plt = _plt()
    from matplotlib.lines import Line2D
    workers = sorted({r['summary']['workers'] for r in rows})
    seconds = rows[0]['summary']['settings']['seconds']
    eps = rows[0]['summary']['settings']['eps']
    # gap(t)
    figure, axes = plt.subplots(1, len(workers), figsize=(6.6*len(workers), 4.8), squeeze=False)
    figure.patch.set_facecolor(SURFACE)
    for ax, w in zip(axes[0], workers):
        for row in rows:
            s = row['summary']
            if s['workers'] != w:
                continue
            series = read_timeline(row['directory'])
            if not series:
                continue
            ts, gaps = zip(*[(t, g) for t, g in series if g > 0])
            color = METHOD_COLORS[s['method']]
            ax.plot(ts, gaps, drawstyle='steps-post', color=color, lw=1.6, alpha=.85)
            if s['t_cert'] is not None:
                ax.plot([s['t_cert']], [gaps[-1]], marker='o', ms=8, color=color, mec=SURFACE, mew=2, zorder=5)
        ax.axhline(eps, color=MUTED, lw=1)
        ax.annotate(f'eps = {eps:g}', (1., eps), xycoords=('axes fraction', 'data'), xytext=(-4, 4),
                    textcoords='offset points', ha='right', color=INK_2, fontsize=8)
        ax.set_xlim(0., seconds*1.02)
        _frame(ax, f'Volume gap (outer∩box - inner) / inner, workers = {w}', 'Run wall-clock time (s)', log=True)
    methods = [m for m in METHOD_COLORS if any(r['summary']['method'] == m for r in rows)]
    handles = [Line2D([], [], color=METHOD_COLORS[m], lw=2, label=f'{m} (each repeat)') for m in methods]
    handles.append(Line2D([], [], color=INK_2, marker='o', ls='none', ms=7, label='certified (all partitions)'))
    figure.legend(handles=handles, loc='lower center', ncol=len(handles), frameon=False, fontsize=9, labelcolor=INK_2)
    figure.tight_layout(rect=(0, .07, 1, 1))
    figure.savefig(output/'gap_time.png', dpi=160, facecolor=SURFACE)
    plt.close(figure)
    # MR/FR at checkpoints
    panels = (('inner', 'mr', 'Inner region MR (%)'), ('outer', 'fr', 'Outer bound FR (%)'))
    figure, axes = plt.subplots(len(workers), 2, figsize=(13., 4.3*len(workers)), squeeze=False)
    figure.patch.set_facecolor(SURFACE)
    for i, w in enumerate(workers):
        for j, (domain, kind, title) in enumerate(panels):
            ax = axes[i][j]
            for row in rows:
                s, m = row['summary'], row['metrics']
                if s['workers'] != w or m is None:
                    continue
                points = sorted((float(k), v['overall']) for k, v in m['checkpoints'].items() if k != 'final')
                for reference, style in (('ac', '-'), ('socp', '--')):
                    xs = [t for t, o in points if o[f'{domain}_{reference}'][f'{kind}_percent'] is not None]
                    ys = [o[f'{domain}_{reference}'][f'{kind}_percent'] for t, o in points
                          if o[f'{domain}_{reference}'][f'{kind}_percent'] is not None]
                    ax.plot(xs, ys, ls=style, marker='o', ms=3.5, lw=1.4, color=METHOD_COLORS[s['method']], alpha=.85)
            _frame(ax, f'{title} vs AC (solid) and SOCP (dashed), workers = {w}', 'Checkpoint time (s)')
            ax.set_xlim(0., seconds*1.02)
    handles = [Line2D([], [], color=METHOD_COLORS[m], lw=2, label=m) for m in methods]
    handles += [Line2D([], [], color=INK_2, lw=1.4, label='vs AC scan'),
                Line2D([], [], color=INK_2, lw=1.4, ls='--', label='vs SOCP scan (diagnostic)')]
    figure.legend(handles=handles, loc='lower center', ncol=len(handles), frameon=False, fontsize=9, labelcolor=INK_2)
    figure.tight_layout(rect=(0, .05, 1, 1))
    figure.savefig(output/'mr_fr_time.png', dpi=160, facecolor=SURFACE)
    plt.close(figure)
    plot_time_breakdown(output, rows)
    if rows[0]['summary'].get('load_nodes') and len(rows[0]['summary']['load_nodes']) == 2:
        plot_regions(output, rows)


def plot_time_breakdown(output, rows):
    plt = _plt()
    from matplotlib.patches import Patch
    groups = {}
    for row in rows:
        s = row['summary']
        groups.setdefault((s['method'], s['workers']), []).append(s['totals'])
    labels = [f'{m}\nworkers {w}' for m, w in sorted(groups)]
    figure, axes = plt.subplots(1, 2, figsize=(13., 4.8))
    figure.patch.set_facecolor(SURFACE)
    for ax, (key, names, colors, title) in zip(axes, (
            ('by_type', ('MISOCP', 'SOCP', 'LP'), TYPE_COLORS, 'Solver time by problem type'),
            ('by_purpose', PURPOSES, PURPOSE_COLORS, 'Solver time by purpose'))):
        bottom = np.zeros(len(groups))
        for name in (*names, 'other'):
            values = []
            for group in sorted(groups):
                items = groups[group]
                if name == 'other':
                    values.append(float(np.median([t['other_seconds'] for t in items])))
                else:
                    values.append(float(np.median([t[key].get(name, {}).get('seconds', 0.) for t in items])))
            values = np.asarray(values)
            ax.bar(range(len(groups)), values, bottom=bottom, width=.6, color=colors.get(name, OTHER_COLOR),
                   edgecolor=SURFACE, lw=2, label=name)
            bottom += values
        ax.set_xticks(range(len(groups)), labels)
        _frame(ax, title, '', 'Seconds (median over repeats, summed over partitions and processes)')
        ax.legend(frameon=False, fontsize=8, labelcolor=INK_2, ncol=2)
    figure.tight_layout()
    figure.savefig(output/'time_breakdown.png', dpi=160, facecolor=SURFACE)
    plt.close(figure)


def plot_regions(output, rows):
    """二维终态：配对扫描为底，叠加 R 与 H（各取 workers 最大、重复序号最小的运行）的内域与外界；H 的 P_x 按网架着色。"""
    plt = _plt()
    from matplotlib.colors import ListedColormap
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from shapely import affinity
    from monitor import _draw
    summary = rows[0]['summary']
    case = Case(summary['case_name'], summary['load_nodes'])
    scan = load_reference(case)
    if scan is None:
        return
    chosen = {}
    for row in sorted(rows, key=lambda r: (-r['summary']['workers'], r['summary']['repeat'])):
        chosen.setdefault(row['summary']['method'], row)
    figure, axes = plt.subplots(1, len(chosen), figsize=(7.4*len(chosen), 7.2), squeeze=False)
    figure.patch.set_facecolor(SURFACE)
    shape, step = scan['shape'], scan['step']
    lower = scan['origin']+scan['start']*step
    edges = [lower[j]+np.arange(shape[j]+1)*step[j] for j in range(2)]
    layer = np.where(scan['states'] == 1, 1, np.where(scan['socp_states'] == 1, 2, 0))
    palette = ('#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#4a3aa7', '#e34948')   # 七个最大的 P_x
    finals = {method: json.loads((row['directory']/'snapshots'/'final.json').read_text(encoding='utf-8'))
              for method, row in chosen.items()}
    areas = {}
    for final in finals.values():
        for state in final.values():
            for network in state.get('networks', []):
                vertices = np.asarray(network['vertices'], float).reshape(-1, 2)
                areas[network['scheme']] = areas.get(network['scheme'], 0.)+polytope_volume(vertices)
    schemes = sorted(areas, key=lambda s: -areas[s])[:len(palette)]
    color_of = lambda scheme: palette[schemes.index(scheme)] if scheme in schemes else MUTED
    for ax, (method, row) in zip(axes[0], sorted(chosen.items())):
        ax.set_facecolor(SURFACE)
        ax.pcolormesh(edges[0], edges[1], layer.T, cmap=ListedColormap([SURFACE, '#cde2fb', '#f9d5c3']), vmin=-.5,
                      vmax=2.5, shading='flat', rasterized=True, zorder=0)
        final = finals[method]
        for code, state in final.items():
            if not state.get('started') or not state['cones']:
                continue
            sign = np.asarray(partition_sign(code, 2), float)
            scale = lambda geometry: affinity.scale(geometry, xfact=sign[0]*case.bounds[0], yfact=sign[1]*case.bounds[1],
                                                    origin=(0., 0.))
            cones = [dict(U=np.asarray(c['U'], float), V=np.asarray(c['V'], float),
                          halfspace=None if c['halfspace'] is None else (np.asarray(c['halfspace'][0], float),
                                                                          float(c['halfspace'][1])))
                     for c in state['cones']]
            sets = [np.asarray(n['vertices'], float).reshape(-1, 2) for n in state['networks']]
            cover = None if state['cover'] is None else [np.asarray(f, float) for f in state['cover']]
            inner, outer = geo.h_geometry(cones, sets, cover)
            _draw(ax, scale(inner), color='#008300', fill=True, alpha=.25, linewidth=0.)
            for network, vertices in zip(state['networks'], sets):
                if len(vertices) >= 3:
                    _draw(ax, scale(geo._union([vertices])), color=color_of(network['scheme']), linewidth=.9)
            _draw(ax, scale(outer), color='#e34948', linewidth=1.3)
        s = row['summary']
        ax.set_xlim(edges[0][0], edges[0][-1])
        ax.set_ylim(edges[1][0], edges[1][-1])
        ax.set_aspect('equal')
        ax.axhline(0., color='#c3c2b7', lw=.6)
        ax.axvline(0., color='#c3c2b7', lw=.6)
        ax.set_xlabel(f'p{case.nodes[0]} (kW, + load / - PV)')
        ax.set_ylabel(f'p{case.nodes[1]} (kW, + load / - PV)')
        ax.set_title(f'Method {method}: {s["status"]}, {s["certified_partitions"]}/{len(s["partitions"])} partitions '
                     f'certified\nworkers {s["workers"]}, run {s["repeat"]}', color=INK, fontsize=10, loc='left')
        ax.grid(color=GRID, lw=.5)
        ax.set_axisbelow(True)
    handles = [Patch(facecolor='#cde2fb', label='AC-feasible cell'),
               Patch(facecolor='#f9d5c3', label='SOCP-only feasible cell'),
               Patch(facecolor='#008300', alpha=.25, label='Inner region (H: radial inner ∪ P_x)'),
               Line2D([], [], color='#e34948', lw=1.3, label='Outer bound (H: O_R ∩ E_x once covered)'),
               *(Line2D([], [], color=palette[k], lw=.9, label=f'P_x, scheme {name}') for k, name in enumerate(schemes)),
               *([Line2D([], [], color=MUTED, lw=.9, label='P_x, other schemes')] if len(areas) > len(schemes) else [])]
    columns = 4
    figure.legend(handles=handles, loc='lower center', ncol=columns, frameon=False, fontsize=8, labelcolor=INK_2)
    figure.tight_layout(rect=(0, .03+.035*(-(-len(handles)//columns)), 1, 1))
    figure.savefig(output/'regions.png', dpi=160, facecolor=SURFACE)
    plt.close(figure)


# ---- 入口 ---------------------------------------------------------------------------------------
def parse_args(argv=None):
    parser = argparse.ArgumentParser(description='Plan 2 (support queries) vs radial sandwich: 300 s comparison')
    parser.add_argument('--case', choices=('case33', 'fourbus'), default='case33')
    parser.add_argument('--nodes', default='18,25', type=lambda v: tuple(int(x) for x in v.split(',')))
    parser.add_argument('--mode', type=int, default=1)
    parser.add_argument('--methods', default='R,H', type=lambda v: tuple(m for m in v.split(',') if m))
    parser.add_argument('--seconds', type=float, default=300., help='wall-clock limit of one run (all partitions)')
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--workers', default='1,16', type=lambda v: tuple(int(x) for x in v.split(',')))
    parser.add_argument('--threads', type=int, default=1, help='Gurobi Threads of every solve')
    parser.add_argument('--partitions', default='all', help='all, or a comma list such as nn,np')
    parser.add_argument('--eps', type=float, help='volume-gap target (default d*tau, the radial volume target)')
    parser.add_argument('--tau', type=float, default=REGION_TAU)
    parser.add_argument('--discovery-eps', type=float, default=.15)
    parser.add_argument('--discovery-share', type=float, default=.25)
    parser.add_argument('--criterion', choices=('origin', 'center', 'auto'), default='auto')
    parser.add_argument('--network-eps', type=float, help='per-network stop vol(O_x)/vol(P_x)-1 (default eps/2)')
    parser.add_argument('--coverage', choices=('global', 'cone'), default='global',
                        help='coverage certificate of method H: one global MISOCP, or one per phase-A cone '
                             '(method Hc always uses cone)')
    parser.add_argument('--seed', type=int, default=0, help='Gurobi Seed of every MISOCP/SOCP model built here')
    parser.add_argument('--mip-seconds', type=float, default=60., help='time cap of one MISOCP (cone or coverage)')
    parser.add_argument('--mip-gap', type=float, default=1e-3)
    parser.add_argument('--inner-cert', choices=('origin', 'near'), default='near')
    parser.add_argument('--min-width', type=float)
    parser.add_argument('--max-cones', type=int)
    parser.add_argument('--scan', type=Path, help='paired AC/SOCP scan cache (default: newest project cache)')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--report', action='store_true', help='only rebuild metrics and the report from existing runs')
    parser.add_argument('--smoke', action='store_true', help='Case33 (18,25), partition nn, 60 s, 1 run, R and H')
    parser.add_argument('--worker', type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    if args.smoke:
        args.case, args.nodes, args.seconds, args.repeats = SMOKE['case'], SMOKE['nodes'], SMOKE['seconds'], SMOKE['repeats']
        args.methods, args.partitions = SMOKE['methods'], ','.join(SMOKE['partitions'])
        if '--workers' not in (argv or sys.argv):
            args.workers = SMOKE['workers']
    d = len(args.nodes)
    args.eps = args.eps if args.eps is not None else d*args.tau
    args.network_eps = args.network_eps if args.network_eps is not None else args.eps/2.
    if set(args.methods)-{'R', 'H', 'Hc'}:
        parser.error('--methods takes R, H and/or Hc (H with the per-cone coverage certificate)')
    tag = f'{args.case}_{"_".join(map(str, args.nodes))}'
    if args.output is None:
        args.output = DEFAULT_OUTPUT/(tag+('_smoke' if args.smoke else ''))
    args.partitions = None if args.partitions == 'all' else tuple(args.partitions.split(','))
    if args.partitions:
        for code in args.partitions:
            partition_sign(code, d)
    return args


def main():
    if '--worker' in sys.argv:
        return worker_main(sys.argv[sys.argv.index('--worker')+1])
    args = parse_args()
    case = Case(args.case, args.nodes, args.mode)
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    scan = load_reference(case, args.scan)
    print(f'{case.tag}: bounds {np.round(case.bounds, 3).tolist()} kW, partitions '
          f'{",".join(args.partitions or case.partitions)}, scan {None if scan is None else scan["path"]}', flush=True)
    cells = None if scan is None else partition_cells(case, scan)
    if not args.report:
        args_stored = argparse.Namespace(**{k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()})
        for repeat in range(1, args.repeats+1):
            for workers in args.workers:
                for method in args.methods:
                    directory = output/method/f'workers_{workers}'/f'run_{repeat}'
                    summary = run_once(args_stored, case, method, workers, repeat, directory)
                    if scan is not None:
                        report = evaluate_run(directory, case, scan, cells)
                        final = report['checkpoints']['final']['overall']
                        print(f'  -> {summary["status"]} t_cert={summary["t_cert"]} gap={summary["gap_final"]} '
                              f'inner-AC FR/MR {final["inner_ac"]["fr_percent"]}/{final["inner_ac"]["mr_percent"]} '
                              f'valid={report["valid"]}', flush=True)
    elif scan is not None:
        for row in run_rows(output):
            evaluate_run(row['directory'], case, scan, cells)
    code = write_report(output, case.tag)
    print(f'report: {output/"comparison.md"}', flush=True)
    errors = [r['summary'] for r in run_rows(output) if any(p.get('status') == 'error'
                                                            for p in r['summary']['partition_results'].values())]
    if errors:
        print('runs with partition errors: ' + ', '.join(f'{s["method"]}/w{s["workers"]}/run{s["repeat"]}'
                                                         for s in errors), flush=True)
    return code or (1 if errors else 0)


if __name__ == '__main__':
    raise SystemExit(main())
