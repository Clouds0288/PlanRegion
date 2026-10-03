"""规划域主线（方法 RCUT：径向锥夹逼 + 主线割平面）：1 参数；2 构域 → 扫描校验，窗口实时显示；3 收敛过程图。
回放用 monitor.py。"""
import argparse
import json
from functools import partial
from pathlib import Path
from types import SimpleNamespace

from Network.four_bus_five_corridor import FourBus
from Network.case33bw import Case33, Case33Plan, CURRENT_LIMIT
from region import build_region
from monitor import RunMonitor
from plot import draw_convergence
from vertify import convergence, export_comparison, reference_grid

# 1. 参数：算例与维数、输出位置、构域与扫描设置
CASE = 'case33'                   # fourbus / case33 / case33plan
DIMENSION = 3                     # 2 / 3；Case33: (18,25) / (18,25,30)
BUDGET = 20000.                   # FourBus 建设预算；Case33 用其 switch_budget，Case33Plan 用其 plan_budget
ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT/'results'/'mainline'   # 记录与逐格对比在 <OUTPUT>/mode_1/，收敛过程图在 <OUTPUT>/
SHOW_UI = True
CASE_TIME_LIMIT = 300             # 构域总时限（墙钟，各分区并行），含初始化和记录，不含事后扫描
WORKERS = 16                      # 并行分区进程数；OBBT 线程数为 WORKERS//仍在计算的分区数
SOLVER_THREADS = 1                # 每个构域求解器的线程数（分区已并行）
REGION_TAU = .005                 # 径向精度；体积目标 ε=d·tau
DISCOVERY_EPS = .15               # A 阶段（网架发现）的放宽体积目标 ε_A
DISCOVERY_SHARE = .25             # A 阶段最多占分区时限的比例
MIP_SECONDS = 60.                 # 单次锥 MISOCP 的时限上限
MIP_GAP = 1e-3                    # 锥 MISOCP 的相对间隙
MIN_WIDTH = {2: 1e-4, 3: 2e-3}    # 锥角直径下限（rad），更窄的锥不再细分
MAX_CONES = {2: 256, 3: 2048}     # 每个分区的叶锥数上限
CUT_THRESHOLD = .01               # 连续小割的体积缩减比例
CUT_PATIENCE = 3
POINT_TOL = 1e-2                  # 同一几何点的最大坐标差，kW
SCAN_OUTPUT = ROOT/'results'/'scan'   # 扫描缓存：按物理身份存放，覆盖的格点直接复用，只补扫缺失的格点
SCAN_DIVISIONS = {2: 160, 3: 80}  # 某配置首次扫描的每轴格数；之后沿用并扩展缓存的格架
SCAN_WORKERS = 20                 # 扫描进程数；并行时每个求解器用 1 个线程

CASES = {   # 算例键 → (网架构造, 预算, 三维节点)；Case33 系列取主线、AC 与 SOCP 共用的电流限额 CURRENT_LIMIT
    'fourbus': (FourBus, BUDGET, (1, 2, 3)),
    'case33': (partial(Case33, current_limit=CURRENT_LIMIT), Case33.switch_budget, (18, 25, 30)),
    'case33plan': (partial(Case33Plan, current_limit=CURRENT_LIMIT), Case33Plan.plan_budget, (18, 25, 30))}


def validate(monitor):
    """扫描校验一份构域记录：网架、预算与构域结果都取自记录；参考网格复用 SCAN_OUTPUT 下同配置的缓存，只补扫缺失的
    格点；内域（主指标）与外域相对 AC 的逐格对比写回记录，并导出到记录旁。返回导出摘要。"""
    state = monitor.state
    kind = next(kind for kind, _, nodes in CASES.values() if kind(load_nodes=nodes).name == state['network'])
    network, budget, result = kind(load_nodes=tuple(state['load_nodes'])), state['budget'], state['result']
    monitor.validation_state.pop('error', None)   # 上次校验失败留下的信息
    grid = reference_grid(network, budget, result, SCAN_OUTPUT, divisions=SCAN_DIVISIONS[len(network.load_nodes)],
                          workers=SCAN_WORKERS, progress=monitor.scanning)
    return export_comparison(network, budget, result, grid, monitor.validation(grid, result), monitor.output)


def main(case=CASE, dimension=DIMENSION, *, seconds=CASE_TIME_LIMIT, output=OUTPUT, show=SHOW_UI):
    """一次运行：构域并保存记录，再扫描校验（窗口实时显示两步，最后在对比面板呈现）；窗口关闭后画收敛过程图。"""
    kind, budget, nodes = CASES[case]
    network = kind(load_nodes=nodes[:dimension])
    monitor = RunMonitor(output=Path(output)/'mode_1'/f"{case}_{'_'.join(map(str, network.load_nodes))}.json.gz",
                         algorithm='RCUT · 径向夹逼 + 割平面')
    settings = SimpleNamespace(workers=WORKERS, threads=SOLVER_THREADS, tau=REGION_TAU, discovery_eps=DISCOVERY_EPS,
                               discovery_share=DISCOVERY_SHARE, mip_seconds=MIP_SECONDS, mip_gap=MIP_GAP,
                               min_width=MIN_WIDTH[dimension], max_cones=MAX_CONES[dimension],
                               threshold=CUT_THRESHOLD, patience=CUT_PATIENCE, point_tol=POINT_TOL)

    def complete():   # 2. 窗口占主线程，构域与校验在后台线程依次执行
        result = build_region(network, budget=budget, monitor=monitor, seconds=seconds, settings=settings)
        monitor.save()   # 构域记录先落盘，之后的扫描可能很长
        validate(monitor)
        return result

    result = monitor.execute(complete, show_ui=show)
    draw_convergence(convergence(monitor), dimension*REGION_TAU, monitor.output)   # 3. 结果图
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='RCUT 规划域：构域 → 扫描校验 → 收敛过程图；回放用 monitor.py')
    parser.add_argument('--case', choices=tuple(CASES), default=CASE)
    parser.add_argument('--dimension', type=int, choices=(2, 3), default=DIMENSION)
    parser.add_argument('--seconds', type=float, default=CASE_TIME_LIMIT, help='构域总时限')
    parser.add_argument('--output', type=Path, default=OUTPUT, help='记录与对比写入 <output>/mode_1/，收敛过程图写入 <output>/')
    parser.add_argument('--no-ui', action='store_true', help='不打开原生窗口，仍保存记录')
    parser.add_argument('--validate', type=Path, metavar='RECORD',
                        help='不构域：对已有记录扫描校验（复用或补扫缓存），导出对比并重画收敛过程图')
    args = parser.parse_args()
    if args.validate:
        monitor = RunMonitor(output=args.validate)
        monitor.load_recording(args.validate)
        summary = validate(monitor)
        monitor.save()
        draw_convergence(convergence(monitor), len(monitor.state['load_nodes'])*REGION_TAU, monitor.output)
        print(json.dumps(summary, ensure_ascii=False, indent=2))
    else:
        main(args.case, args.dimension, seconds=args.seconds, output=args.output, show=SHOW_UI and not args.no_ui)
