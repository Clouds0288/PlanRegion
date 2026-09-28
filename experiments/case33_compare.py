"""Case33：主线与 KKT 认证并集排除法，共用开关预算和二维 SOCP 参考。

python experiments/case33_compare.py --seconds 1000 --epsilon .01
python experiments/case33_compare.py --method kkt --show
python monitor.py results/case33_compare/mainline.json.gz --compare results/case33_compare/kkt.json.gz
"""
import argparse
from pathlib import Path
import sys

import numpy as np
from gurobipy import GRB
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from Network.case33bw import Case33
from experiments.fourbus_outer_kkt import COVERAGE_PAD, UncoveredKktViolation, certify_load
from experiments.fourbus_outer_volume import budget_schemes
from main import build_continuous_region, REGION_TAU
from model import GridPhysics, MasterProblem, SubProblem, PLANNING_TOL
from monitor import RunMonitor, RegionTimeout, _draw, _union
from region import RegionState, initial_polytope, clip_polytope
from vertify import validate_socp_region


def projected_outer(equations, budget, axis_bounds, total_bound, cuts, threads):
    """仅在构域结束后，枚举小算例合法网架以绘制联合割的精确投影并集。"""
    bounds = np.asarray(axis_bounds)
    result = []
    for x in budget_schemes(equations, budget, threads):
        poly = initial_polytope(bounds, total_bound, bounds)
        for cut in cuts:
            poly = clip_polytope(poly, cut[0]+cut[3:]@x, cut[1:3]*bounds)
        result.append(dict(vertices=poly*bounds))
    return result


def build_kkt_region(network, *, budget, seconds, epsilon, threads, monitor):
    """1 MP2 初始化；2 排除并集的 G；3 固定 x,p 的 SP；4 换网架认证；5 更新。"""
    equations = GridPhysics(network, 'socp')
    bounds = np.full(2, network.power_limit)
    region = RegionState(bounds, network.power_limit, 0.)
    monitor.begin(network, 'socp', budget, region, seconds)
    counts = dict(global_search=0, sp=0, cuts=0, certification=0, support=0)

    def register(x):
        region.add_scheme(x, network.decode_plan(x), network.cost_offset+network.cost@x)

    def support(x):
        # 固定背景负荷不为零：负总负荷方向取得真实下侧证书，不默认认证原点。
        for direction in [-np.ones(2), *np.eye(2), np.ones(2)]:
            monitor.support_start(x, direction)
            problem = MasterProblem(equations, budget=budget, direction=direction,
                                    fixed_plan=network.decode_plan(x), threads=threads)
            with problem.model:
                answer = problem.solve(time_limit=monitor.remaining(seconds))
            counts['support'] += 1
            region.add_point(x, answer['p']/bounds)
            monitor.support_end(answer, region)

    # 1. 三个完整 MP2 独立初始化；x,p,y 全部自由。每个点只加入所属网架。
    for index, direction in enumerate([*np.eye(2), np.ones(2)]):
        monitor.initializing(index)
        problem = MasterProblem(equations, budget=budget, direction=direction, threads=threads)
        with problem.model:
            answer = problem.solve(time_limit=monitor.remaining(seconds))
        axis = region.axis_bounds.copy()
        if index < 2:
            axis[index] = answer['bound']
        region.tighten_bounds(axis, answer['bound'] if index == 2 else region.total_bound)
        register(answer['x'])
        region.add_point(answer['x'], answer['p']/bounds)
        monitor.seed(answer, region)
    supported = set(region.records)
    for key in supported:
        support(np.asarray(key))

    def make_search():
        inner = [row['inner']*bounds/region.axis_bounds for row in region.records.values() if len(row['inner'])]
        return UncoveredKktViolation(equations, budget, region.axis_bounds, region.total_bound,
                                     inner, region.cuts, threads=threads)

    oracle = MasterProblem(equations, budget=budget, cuts_only=True, threads=threads)
    eta = oracle.model.addVar(name='violation')
    equations.add_operation(oracle.model, oracle.choices, oracle.loads, eta)
    oracle.model.setObjective(eta, GRB.MINIMIZE)
    sp = SubProblem(equations, threads=threads)
    search = make_search()
    template = search.template
    relaxed = template.sp_eta < 0.
    padding_bound = float(COVERAGE_PAD*np.max(np.abs(template.sp_p[relaxed])@region.axis_bounds
                                            /(-template.sp_eta[relaxed]))+PLANNING_TOL)
    bound = float(search.violation.UB)
    status, certified = 'time_limit', False
    try:
        with oracle.model:
            while monitor.timing()['total_seconds'] < seconds:
                # 2. G 同时选择 x,p；每个认证内域都被排除。只用全局上界判停。
                monitor.global_start(0)
                answer = search.solve(monitor.remaining(30.), epsilon)
                counts['global_search'] += 1
                bound = answer['bound']
                certified = max(bound, padding_bound) <= epsilon
                if answer['x'] is not None and not certified:
                    register(answer['x'])
                monitor.global_end(dict(answer, complete=certified,
                    x=None if certified else answer['x']), region)
                if certified:
                    status = 'union_residual_certified'
                    break
                if monitor.timing()['total_seconds'] >= seconds:
                    break
                if answer['x'] is None:
                    continue  # G 此次超时尚无候选，继续求解；有候选时交给原始 SP 判定。

                # 3. 固定 G 给出的 x,p，原始 SP 最小化 eta，产生适用于所有 x 的联合割。
                x, power = answer['x'], answer['p']
                monitor.sp_start(x, power, sp.calls+1)
                checked = sp.solve(x, power, time_limit=monitor.remaining(30.))
                counts['sp'] = sp.calls
                monitor.sp_end(checked)
                if checked['feasible']:
                    region.add_point(x, power/bounds)
                    monitor.updated(region, 'feasible', x=x, power=power, checked=checked)
                    answer = dict(feasible=True, x=x)
                else:
                    region.apply_cut(checked['cut'])
                    counts['cuts'] += 1
                    monitor.updated(region, 'cut', x=x, power=power, checked=checked)
                    search.add_cut(checked['cut'])
                    if monitor.timing()['total_seconds'] >= seconds:
                        break

                    # 4. 固定同一个 p，放开全部合法 x,y 认证；失败不直接删除该负荷点。
                    monitor.certification_start(power)
                    answer = certify_load(oracle, eta, power, monitor.remaining(seconds))
                    counts['certification'] += 1
                    if answer['feasible']:
                        register(answer['x'])
                        region.add_point(answer['x'], power/bounds)
                    monitor.certification_end(answer, region)
                    if answer['status'] == 'time_limit':
                        break

                # 5. 认证成功后扩大该网架内域；新网架补方向支撑，再重建并集排除约束。
                if answer['feasible']:
                    if monitor.timing()['total_seconds'] >= seconds:
                        break
                    key = tuple(answer['x'])
                    if key not in supported:
                        support(answer['x'])
                        supported.add(key)
                    search.model.dispose()
                    search = make_search()
                    search.violation.UB = bound
    finally:
        search.model.dispose()
    result = dict(method='socp', budget=budget, status=status, certified=certified,
                  epsilon=epsilon, bound_scope='uncovered_pairs', bound=bound,
                  padding_bound=padding_bound, union_bound=max(bound, padding_bound),
                  coverage_bound=max(bound, padding_bound), axis_bounds=region.axis_bounds,
                  counts=counts, timing=monitor.timing(), **region.finish(False))
    # 事后绘制外域，不让枚举和扫描参与 G 的选点或停止。
    result['outer'] = projected_outer(equations, budget, region.axis_bounds, region.total_bound,
                                      region.cuts, threads)
    monitor.finish(result, region)
    return result


def plot_comparison(monitors, output):
    """紧凑论文分面：同一 SOCP 扫描、各方法实际内域与外域；不跨网架取凸包。"""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.lines import Line2D
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'svg.fonttype': 'none',
                         'pdf.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(1, len(monitors), figsize=(4.1*len(monitors), 3.7), squeeze=False)
    for index, (name, monitor) in enumerate(monitors.items()):
        ax = axes[0, index]
        state = monitor.state
        reference, result = state['validation'], state['result']
        bounds, states = np.asarray(reference['bounds']), np.asarray(reference['states'])
        size = len(states)
        ax.pcolormesh(np.linspace(0, bounds[0], size+1), np.linspace(0, bounds[1], size+1),
                      (states.T == 1).astype(int), cmap=ListedColormap(['white', '#d8e4ed']),
                      vmin=0, vmax=1, rasterized=True)
        _draw(ax, _union([row['vertices'] for row in result['outer']]),
              color='#ba6036', linewidth=1.4, linestyle='--')
        _draw(ax, _union([row['vertices'] for row in result['inner']]),
              color='#268078', linewidth=1.4)
        ax.set(xlim=(0, bounds[0]*1.02), ylim=(0, bounds[1]*1.02),
               xlabel=r'$p_{18}$ (kW)', ylabel=r'$p_{25}$ (kW)')
        ax.text(-.14, 1.03, chr(97+index), weight='bold', transform=ax.transAxes, fontsize=12)
        ax.set_title('Mainline' if name == 'mainline' else 'KKT: certified-union exclusion', fontsize=10)
    handles = [Line2D([], [], color='#a7becf', lw=6, label='SOCP grid reference'),
               Line2D([], [], color='#268078', label='Certified inner union'),
               Line2D([], [], color='#ba6036', ls='--', label='Outer union')]
    fig.legend(handles=handles, loc='lower center', ncol=3, frameon=False, fontsize=8)
    fig.subplots_adjust(left=.09, right=.98, top=.90, bottom=.24, wspace=.29)
    for extension in ('png', 'svg'):
        fig.savefig(output/f'comparison.{extension}', dpi=220)
    plt.close(fig)


def run_comparison(*, method='both', seconds=1000., epsilon=.01, divisions=40, threads=4,
                   output=ROOT/'results/case33_compare', show=False):
    network = Case33(load_nodes=(18, 25))
    budget = network.switch_budget
    output = Path(output)
    monitors = {}
    reference = None
    names = ('mainline', 'kkt') if method == 'both' else (method,)
    with threadpool_limits(limits=1):
        for name in names:
            def progress(event, **values):
                if event == 'residual_end':
                    answer = values['answer']
                    print(f"{name}: G bound={answer['bound']}, complete={answer['complete']}", flush=True)
            monitor = RunMonitor(output=output/f'{name}.json.gz', callback=progress,
                                 algorithm='主线' if name == 'mainline' else 'KKT 认证域排除')
            monitors[name] = monitor

            def calculate():
                nonlocal reference
                try:
                    if name == 'mainline':
                        result = build_continuous_region(network, 'socp', budget,
                            np.full(2, network.power_limit), tau=REGION_TAU, residual_mode='physical',
                            time_limit=seconds, threads=threads, progress=monitor)
                    else:
                        result = build_kkt_region(network, budget=budget, seconds=seconds, epsilon=epsilon,
                                                  threads=threads, monitor=monitor)
                except (RegionTimeout, TimeoutError):
                    result = monitor.stopped('time_limit')
                if show:
                    if reference is None:
                        reference = validate_socp_region(network, budget, divisions, result['axis_bounds'],
                                                         threads=threads, progress=monitor.scanning)
                    monitor.validation(reference, result, region_key='inner' if name == 'mainline' else 'outer')
                return result

            result = monitor.execute(calculate, show_ui=show)
            print(name, result['status'], result['timing'], result['counts'], flush=True)

        # 6. 两算法结束后独立扫描；共享完全相同的格点和全网架物理参考。
        bounds = np.max([m.state['result']['axis_bounds'] for m in monitors.values()], axis=0)
        def scanning(done, total):
            for monitor in monitors.values():
                monitor.scanning(done, total)
            if done % (divisions*5) == 0:
                print(f'SOCP reference {done}/{total}', flush=True)
        if reference is None:
            reference = validate_socp_region(network, budget, divisions, bounds,
                                             threads=threads, progress=scanning)
        for name, monitor in monitors.items():
            result = monitor.state['result']
            monitor.validation(reference, result, region_key='inner' if name == 'mainline' else 'outer')
            monitor.save()
            print(name, monitor.state['validation']['metrics'], flush=True)
    plot_comparison(monitors, output)
    return monitors


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--method', choices=('both', 'mainline', 'kkt'), default='both')
    parser.add_argument('--seconds', type=float, default=1000.)
    parser.add_argument('--epsilon', type=float, default=.01)
    parser.add_argument('--divisions', type=int, default=40)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--output', type=Path, default=ROOT/'results/case33_compare')
    parser.add_argument('--show', action='store_true')
    run_comparison(**vars(parser.parse_args()))
