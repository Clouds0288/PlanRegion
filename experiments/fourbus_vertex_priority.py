"""FourBus: replace ordinary mainline vertex ordering with SP violation ordering.

python experiments/fourbus_vertex_priority.py
Both methods retain tau=.002 and the same global coverage certificate.
"""
import argparse
import csv
import json
from pathlib import Path
import sys
from time import perf_counter

import gurobipy as gp
from gurobipy import GRB
import numpy as np
from threadpoolctl import threadpool_limits

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from experiments.fourbus_mainline_comparison import region_radial_capacity
from continuous import build_continuous_region, REFINEMENT_CHECKS, REGION_TAU
from model import (GridPhysics, MasterProblem, SubProblem, RemainingRegionModel,
                   new_model, PLANNING_TOL, MP_TIME_LIMIT, SP_TIME_LIMIT,
                   RESIDUAL_TIME_LIMIT)
from Network.four_bus_five_corridor import FourBus
from plot import json_value
from region import RegionState, contains, GEOMETRY_TOL


class ScoredSubProblem(SubProblem):
    """Original SP formulation/cut algebra, plus raw SOCP eta and exact-key reuse."""

    def __init__(self, equations, *, threads=4):
        super().__init__(equations, threads=threads)
        self.cache, self.cache_hits = {}, 0

    def solve(self, x, power, time_limit=None):
        key = (tuple(x), tuple(power))
        if key in self.cache:
            self.cache_hits += 1
            return self.cache[key]
        self.calls += 1
        equations, net = self.equations, self.equations.network
        limit = SP_TIME_LIMIT['socp'] if time_limit is None else time_limit
        deadline = perf_counter()+limit
        with new_model('vertex_score_SP', self.threads) as model:
            # 1. 与原 SP 相同的变量、固定参数等式及数值设置。
            model.Params.Aggregate = 0
            model.Params.ScaleFlag = 0
            model.Params.BarQCPConvTol = 1e-9
            choice = model.addVars(equations.keys, ub=1., name='x')
            p = model.addVars(net.load_nodes, lb=-GRB.INFINITY, name='p_kw')
            fixed_x = model.addConstrs((choice[e, k] == value for (e, k), value
                                       in zip(equations.keys, x)), name='fixed_x')
            fixed_p = model.addConstrs((p[i] == value for i, value in zip(net.load_nodes, power)),
                                       name='fixed_p')
            eta = model.addVar(name='violation')
            operation = equations.add_operation(model, choice, p, eta)
            model.setObjective(eta)
            model.Params.TimeLimit = max(0., deadline-perf_counter())
            model.optimize()
            if model.Status != GRB.OPTIMAL or model.MaxVio > PLANNING_TOL:
                raise RuntimeError(f'Scoring SP: status={model.Status}, MaxVio={model.MaxVio}')

            # 2. 在替换锥之前保存原始 SOCP eta；可行判据与主线完全相同。
            value = float(eta.X)
            if max(0., value)+model.MaxVio <= PLANNING_TOL:
                answer = dict(eta=value, feasible=True, state=operation.state.X, cut=None)
            else:
                if value <= PLANNING_TOL:
                    raise RuntimeError(f'Scoring certificate: eta={value}, MaxVio={model.MaxVio}')
                # 3. 原 SP 的有效锥支撑平面；正 eta 生成可缓存的联合割。
                planes = []
                for head, tail in operation.cones:
                    values = np.array([item.getValue() for item in tail])
                    length = np.linalg.norm(values)
                    if length:
                        planes.append(head-gp.quicksum(float(a/length)*item
                                                       for a, item in zip(values, tail)))
                model.remove(model.getQConstrs())
                for plane in planes:
                    model.addConstr(plane >= 0., name='cone_support')
                model.Params.TimeLimit = max(0., deadline-perf_counter())
                model.optimize()
                if (model.Status != GRB.OPTIMAL or model.ObjVal <= 0.
                        or model.MaxVio > PLANNING_TOL):
                    raise RuntimeError(f'Scoring cut LP: status={model.Status}, MaxVio={model.MaxVio}')
                # 4. 直接复用主线的对偶割代数及分离检查。
                cut = self._separating_cut(model, operation, choice, p,
                                          [*fixed_x.values(), *fixed_p.values()], x, power)
                answer = dict(eta=value, feasible=False, state=None, cut=cut)
        # 5. 只复用完全相同的 x,p；外域收紧不改变原始 SP 的值或证书。
        self.cache[key] = answer
        return answer


def build_vertex_priority_region(network, budget, bounds, *, tau=REGION_TAU,
                                 time_limit=300., threads=4,
                                 refinement_checks=REFINEMENT_CHECKS, eta_trigger=None):
    """1 初始化 → 2 收集候选 → 3 按需全局搜索 → 4 SP 评分/选割 → 5 汇总。"""
    started = perf_counter()
    deadline = started+time_limit
    bounds = np.asarray(bounds, dtype=float)
    equations = GridPhysics(network, 'socp')
    region = RegionState(bounds, network.power_limit, tau)
    oracle = ScoredSubProblem(equations, threads=threads)
    trace = []
    global_search = scoring_rounds = sp_since_global = 0
    global_requested = False
    witness = None
    d = len(bounds)

    def remaining_time(limit):
        remaining = deadline-perf_counter()
        if remaining <= 0.:
            raise TimeoutError('Vertex priority: no global coverage certificate within time limit')
        return min(limit, remaining)

    # 1. 独立运行同样的 d+1 个完整 MP2；不继承对照组或历史实验的割。
    for index, direction in enumerate([*np.eye(d), np.ones(d)]):
        problem = MasterProblem(equations, budget=budget, direction=direction, threads=threads)
        with problem.model:
            answer = problem.solve(time_limit=remaining_time(MP_TIME_LIMIT))
        if answer is None:
            raise RuntimeError('FourBus initialization is infeasible at the specified budget')
        axis_bounds = region.axis_bounds.copy()
        if index < d:
            axis_bounds[index] = min(axis_bounds[index], answer['bound'])
            region.tighten_bounds(axis_bounds, region.total_bound)
        else:
            maximum_answer = answer
            region.tighten_bounds(axis_bounds, answer['bound'])
        x = answer['x']
        region.add_scheme(x, network.decode_plan(x), network.cost_offset+network.cost@x)
        region.add_point(x, answer['p']/bounds)

    while True:
        remaining_time(np.inf)
        candidate = checked = None
        # 活动见证仍补同网架支撑；不改变主线的内部空隙处理逻辑。
        if witness is not None:
            x = witness['x']
            point = (1-region.tau)*witness['p']/bounds
            if region.covering_schemes([point], preferred=x)[0] is not None:
                witness = None
            else:
                support = [q for q in region.witness_support(x, point)
                           if not contains([q], region.inner_equations(x))[0]]
                point = max(support, key=lambda q: (float(q@bounds), tuple(q)))
                candidate = x, point
                item = dict(action='support', candidate_count=len(support))

        if candidate is None:
            # 2. 固定本轮候选集合：收缩外域顶点，且未被任何已认证内域覆盖。
            candidates = []
            for row in region.records.values():
                targets = (1-region.tau)*row['outer']
                owners = region.covering_schemes(targets, preferred=row['x'])
                candidates.extend((row['x'], point) for point, owner
                                  in zip(targets, owners) if owner is None)

            # 3. 候选空、间隔达到阈值或低违反请求时，调用同一个全局模型。
            trigger = [] if candidates else ['empty']
            if sp_since_global >= refinement_checks:
                trigger.append('interval')
            if global_requested:
                trigger.append('eta')
            if trigger:
                problem = RemainingRegionModel(equations, budget, bounds, region.total_bound,
                                               region.cuts, region.inner_halfspaces(), region.tau,
                                               axis_bounds=region.axis_bounds, mode='light', threads=threads)
                with problem.model:
                    witness = problem.solve(GEOMETRY_TOL,
                                             time_limit=remaining_time(RESIDUAL_TIME_LIMIT))
                global_search += 1
                trace.append(dict(action='complete' if witness['complete'] else 'global',
                                  trigger=trigger,
                                  coverage_bound=witness['bound'], x=witness['x'], p=witness['p'],
                                  sp=oracle.calls, cuts=len(region.cuts), seconds=perf_counter()-started))
                if witness['complete']:
                    coverage = witness['bound']
                    break
                x = witness['x']
                region.add_scheme(x, network.decode_plan(x), network.cost_offset+network.cost@x)
                sp_since_global = 0
                global_requested = False
                continue

            # 4. 全部候选求 SP；先收入所有可行点，再挑仍未覆盖者中 eta 最大的一条割。
            before = oracle.calls
            scored = [(x, point, oracle.solve(x, point*bounds,
                                             time_limit=remaining_time(SP_TIME_LIMIT['socp'])))
                      for x, point in candidates]
            sp_since_global += oracle.calls-before
            scoring_rounds += 1
            for x, point, answer in scored:
                if answer['feasible']:
                    region.add_point(x, point)
            owners = region.covering_schemes([point for _, point, _ in scored])
            pending = [(x, point, answer) for (x, point, answer), owner in zip(scored, owners)
                       if not answer['feasible'] and owner is None]
            item = dict(action='score', candidate_count=len(candidates), new_sp=oracle.calls-before,
                        feasible_points=sum(answer['feasible'] for _, _, answer in scored))
            if not pending:
                trace.append(dict(item, sp=oracle.calls, cuts=len(region.cuts), seconds=perf_counter()-started))
                continue
            x, point, checked = max(pending, key=lambda row: (row[2]['eta'], float(row[1]@bounds),
                                                            tuple(row[0]), tuple(row[1])))
            candidate = x, point

        x, point = candidate
        if checked is None:
            before = oracle.calls
            checked = oracle.solve(x, point*bounds, time_limit=remaining_time(SP_TIME_LIMIT['socp']))
            sp_since_global += oracle.calls-before
            item['new_sp'] = oracle.calls-before
        if checked['feasible']:
            region.add_point(x, point)
        else:
            cut = checked['cut']
            region.apply_cut(cut)
            if (witness is not None
                    and cut[0]+cut[1:1+d]@witness['p']+cut[1+d:]@witness['x'] < -1e-12):
                witness = None
        # 本轮最大原始 eta 已小：用完已取得的割，下一轮转全局；不作为停止判据。
        if eta_trigger is not None and item['action'] == 'score' and checked['eta'] < eta_trigger:
            global_requested = True
        item.update(eta=checked['eta'], feasible=checked['feasible'], x=x, p=point*bounds,
                    sp=oracle.calls, cuts=len(region.cuts), seconds=perf_counter()-started)
        trace.append(item)

    # 5. 只有全局证书才能结束；直接调用原来的内域/外包络输出。
    geometry = region.finish(True)
    return dict(method='socp', budget=float(budget), tau=region.tau, residual_mode='light',
                refinement_checks=refinement_checks, eta_trigger=eta_trigger,
                status='certified', coverage_bound=coverage, axis_bounds=region.axis_bounds.copy(),
                max_total=float(maximum_answer['p'].sum()), max_total_bound=region.total_bound,
                counts=dict(sp=oracle.calls, cuts=len(region.cuts), global_search=global_search,
                            scoring_rounds=scoring_rounds, cache_hits=oracle.cache_hits),
                timing=dict(total_seconds=perf_counter()-started), cuts=region.cuts, trace=trace, **geometry)


def render_vertex_priority_comparison(trials, output):
    """Scientific comparison: all 861 ray errors and all independently timed runs."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
                         'font.size': 8, 'axes.labelsize': 8, 'xtick.labelsize': 7,
                         'ytick.labelsize': 7, 'legend.fontsize': 7,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.linewidth': .7, 'svg.fonttype': 'none', 'pdf.fonttype': 42})
    fig, axes = plt.subplots(1, 2, figsize=(7.2047244, 2.5))
    fig.subplots_adjust(left=.10, right=.98, bottom=.23, top=.90, wspace=.45)
    names = {'mainline': 'Total-load priority', 'vertex_priority': 'SP-violation priority'}
    colors = {'mainline': '#355F7D', 'vertex_priority': '#BD743B'}
    for index, (method, group) in enumerate(trials.items()):
        for repeat, item in enumerate(group):
            values = np.sort(item['difference_percent'])
            axes[0].step(values, 100.*np.arange(1, len(values)+1)/len(values), where='post',
                         color=colors[method], linestyle='-' if index == 0 else '--',
                         linewidth=1.2, alpha=.7, label=names[method] if repeat == 0 else None)
        times = [item['total_seconds'] for item in group]
        axes[1].scatter(index+np.linspace(-.07, .07, len(times)), times,
                        color=colors[method], s=19, zorder=3)
        axes[1].plot([index-.18, index+.18], [np.median(times)]*2, color=colors[method], linewidth=1.5)
    axes[0].set(xlabel='Outer-boundary overestimate (%)', ylabel='Cumulative ray fraction (%)',
                ylim=(0., 103.))
    axes[0].legend(frameon=False, loc='lower right')
    axes[1].set(xlim=(-.4, 1.4), ylim=(0., 1.18*max(item['total_seconds'] for group in trials.values()
                                                  for item in group)), ylabel='Runtime (s)',
                xticks=[0, 1], xticklabels=['Total-load\npriority', 'SP-violation\npriority'])
    for ax, label in zip(axes, 'ab'):
        ax.text(-.17, 1.05, label, transform=ax.transAxes, fontweight='bold', fontsize=9)
    fig.savefig(output/'comparison.svg', facecolor='white')
    fig.savefig(output/'comparison.pdf', facecolor='white')
    fig.savefig(output/'comparison.png', dpi=600, facecolor='white')
    plt.close(fig)


def run_vertex_priority_comparison(output, reference_file, repeats=3):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    reference = json.loads(Path(reference_file).read_text(encoding='utf-8'))
    assert reference['budget'] == 20000. and reference['method'] == 'socp'
    network = FourBus()
    bounds = np.full(len(network.load_nodes), network.power_limit)
    weights = np.array(reference['weights'])
    reference_total = np.array(reference['reference_total'])
    raw = {'mainline': [], 'vertex_priority': []}

    # 1. 同预算、同 tau、同全局证书；计时运行交替先后，全部从零初始化。
    with threadpool_limits(limits=1):
        for repeat in range(repeats):
            order = ('mainline', 'vertex_priority') if repeat % 2 == 0 else ('vertex_priority', 'mainline')
            for method in order:
                if method == 'mainline':
                    trace, cuts = [], []

                    def record(event, **data):
                        if event == 'cut':
                            cuts.append(data['cut'].copy())
                        if event in ('cut', 'residual_end'):
                            trace.append(dict(action=event, **data['counts_algorithm'],
                                              coverage_bound=data.get('coverage_bound')))

                    result = build_continuous_region(network, 'socp', 20000., bounds,
                                                     tau=REGION_TAU, residual_mode='light',
                                                     threads=4, time_limit=300., progress=record)
                    result['counts']['global_search'] = sum(row['action'] == 'residual_end' for row in trace)
                    result.update(trace=trace, cuts=cuts)
                else:
                    result = build_vertex_priority_region(network, 20000., bounds, threads=4)
                raw[method].append(result)
                print(method, repeat+1, result['timing'], result['counts'], flush=True)
                (output/'runs.json').write_text(json.dumps(json_value(raw), indent=2), encoding='utf-8')

    # 2. 事后评价两个实际返回的内域/外包络；参考射线不参与算法。
    trials = {'mainline': [], 'vertex_priority': []}
    for method, group in raw.items():
        for result in group:
            assert result['status'] == 'certified'
            assert result['coverage_bound'] is None or result['coverage_bound'] <= GEOMETRY_TOL
            inner_total = region_radial_capacity(result['inner'], weights, bounds)
            outer_total = region_radial_capacity(result['outer'], weights, bounds)
            assert np.max(inner_total-reference_total) <= 1e-4
            assert np.min(outer_total-reference_total) >= -1e-4
            difference_percent = 100.*(outer_total-reference_total)/reference_total
            trials[method].append(dict(total_seconds=result['timing']['total_seconds'],
                                       counts=result['counts'], coverage_bound=result['coverage_bound'],
                                       inner_total=inner_total, outer_total=outer_total,
                                       difference_percent=difference_percent,
                                       mean_radial_difference_percent=float(difference_percent.mean()),
                                       max_radial_difference_percent=float(difference_percent.max()),
                                       max_radial_difference=float((outer_total-reference_total).max()),
                                       max_inner_difference_percent=float(np.max(100.*(reference_total-inner_total)/reference_total))))

    # 3. 全部射线与关键迭代记录可复核；不导出运行状态 y。
    with (output/'rays_comparison.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream)
        writer.writerow(['method', 'repeat', 'w1', 'w2', 'w3', 'reference_total', 'inner_total', 'outer_total'])
        for method, group in trials.items():
            for repeat, item in enumerate(group, 1):
                writer.writerows([method, repeat, *w, reference_total[j], item['inner_total'][j], item['outer_total'][j]]
                                 for j, w in enumerate(weights))
    keys = ('total_seconds', 'mean_radial_difference_percent', 'max_radial_difference_percent',
            'max_radial_difference', 'max_inner_difference_percent')
    summary = {method: {key: dict(median=float(np.median([item[key] for item in group])),
                                  minimum=min(item[key] for item in group), maximum=max(item[key] for item in group))
                        for key in keys} for method, group in trials.items()}
    for method, group in trials.items():
        summary[method]['counts'] = [item['counts'] for item in group]
        summary[method]['coverage_bounds'] = [item['coverage_bound'] for item in group]
    (output/'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')

    # 4. 科研图保留全部重复运行和射线误差，时间用散点及中位数。
    render_vertex_priority_comparison(trials, output)
    # 5. 报告同一覆盖判据下的耗时、计算量与边界误差。
    lines = ['# FourBus：主线总负荷排序与 SP 违反量排序', '',
             f'SOCP、预算 20,000 元、tau={REGION_TAU}；4 个求解器线程、1 个数值库线程；各独立运行 {repeats} 次。'
             '时间包含初始化与全部评分 SP，不包含事后射线评价和绘图。两方法保留相同全局覆盖证书。', '',
             '| 指标（中位数） | 总负荷排序 | SP 违反量排序 |', '|---|---:|---:|']
    for label, key in [('总耗时 / s', 'total_seconds'), ('平均外域高估 / %', 'mean_radial_difference_percent'),
                       ('最大外域高估 / %', 'max_radial_difference_percent'), ('最大外域高估 / kW', 'max_radial_difference'),
                       ('最大内域低估 / %', 'max_inner_difference_percent')]:
        lines.append(f"| {label} | {summary['mainline'][key]['median']:.6f} | {summary['vertex_priority'][key]['median']:.6f} |")
    for label, key in [('实际 SP 次数', 'sp'), ('采用的联合割数', 'cuts'), ('全局搜索次数', 'global_search')]:
        values = [np.median([row[key] for row in summary[method]['counts']]) for method in trials]
        lines.append(f'| {label} | {values[0]:g} | {values[1]:g} |')
    lines += ['', '评分仅改变普通候选处理：固定本轮全部候选求 SP，收入可行证书，再选仍未覆盖者中 eta 最大的一条割。'
              '其他结果精确缓存；活动全局见证仍使用原来的同网架支撑处理。检查间隔按实际新增 SP 次数累计。', '',
              '![误差与耗时比较](comparison.png)', '',
              'a：各次运行全部 861 条 SOCP 参考射线的外域相对容量高估经验分布；b：各次实际总耗时，横线为中位数。'
              '百分比为 100(外域总负荷−参考总负荷)/参考总负荷。参考允许所有合法网架和运行变量自由选择；'
              '有限射线比较不是 AC 真值，也不代替全局连续域覆盖证明。', '',
              'runs.json 保存各次内外域、采用的割与关键迭代；summary.json 为摘要；rays_comparison.csv 为全部射线数据。']
    (output/'report.md').write_text('\n'.join(lines)+'\n', encoding='utf-8')
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default='results/fourbus_vertex_priority')
    parser.add_argument('--reference', default='results/fourbus_outer_partition/epsilon_0p1/reference_rays.json')
    parser.add_argument('--repeats', type=int, default=3)
    args = parser.parse_args()
    run_vertex_priority_comparison(args.output, args.reference, args.repeats)
