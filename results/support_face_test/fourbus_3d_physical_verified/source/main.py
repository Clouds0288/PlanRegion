"""规划域主线：参数、计算入口与原生回放入口。"""
import argparse
from pathlib import Path

from Network.four_bus_five_corridor import FourBus
from Network.case33bw import Case33
from region import build_region
from monitor import RunMonitor, SynchronizedReplay
from vertify import scan_path, save_scan, scan_reference

# 运行设置：DIMENSION 同时控制计算、扫描及回放的维数。
mode = 1                         # 0: 非负负荷；1: 负光伏至正负荷
NETWORK = Case33                 # FourBus / Case33
DIMENSION = 3                     # 2 / 3；Case33: (18,25) / (18,25,30)
BUDGET = 20000.                   # FourBus 建设预算；Case33 用其 switch_budget
CASE_TIME_LIMIT = 300             # 构域总时限，含初始化和记录，不含事后扫描
SOLVER_THREADS = 20                # 构域求解器线程数；扫描每个进程用1个线程
REGION_TAU = .005                 # 全局覆盖的径向精度；不是可行性容差
CUT_THRESHOLD = .02               # 连续小割的面积/体积比例
CUT_PATIENCE = 3
RAY_THRESHOLD = 1e-4              # 射线收益不足本网架认证测度的 0.01% 时跳过；首轮不筛选
POINT_TOL = 1e-2                  # 同一几何点的最大坐标差，kW
DIVISIONS = 160                   # FourBus 二维每轴扫描格数
SCAN_DIVISIONS = {2: 160, 3: 80}   # Case33 及 FourBus 三维每轴扫描格数
SCAN_WORKERS = 20                 # 扫描进程数；并行时每个求解器用 1 个线程
FORCE_RESCAN = False              # True 强制重扫；False 在已有功率范围内复用
SHOW_UI = True
ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT/'results'/'mainline'
SCAN_OUTPUT = ROOT/'results'/'scans'


def recording_path(case, output=OUTPUT, load_nodes=None, dimension=DIMENSION, mode=mode):
    nodes = ({'fourbus': (1, 2, 3), 'case33': (18, 25, 30)}[case][:dimension]
             if load_nodes is None else load_nodes)
    return Path(output)/f"mode_{mode}"/f"{case}_{'_'.join(map(str, nodes))}.json.gz"


def run(network, *, budget=BUDGET, divisions=DIVISIONS, show_ui=SHOW_UI, output=None,
        tau=REGION_TAU, time_limit=CASE_TIME_LIMIT, threads=SOLVER_THREADS,
        scan=True, reference=None, scan_workers=SCAN_WORKERS,
        force_rescan=FORCE_RESCAN, scan_output=SCAN_OUTPUT, mode=mode):
    monitor = RunMonitor(output=output, algorithm='逐网架主线')

    def calculate():
        result = build_region(network, budget=budget, monitor=monitor, mode=mode, seconds=time_limit,
            threads=threads, tau=tau, threshold=CUT_THRESHOLD, patience=CUT_PATIENCE,
            point_tol=POINT_TOL, ray_threshold=RAY_THRESHOLD)
        monitor.save()
        if scan or reference is not None:
            path = scan_path(network, budget, scan_output, mode)
            if reference is not None and not force_rescan:
                previous = RunMonitor()
                previous.load_recording(reference)
                state = previous.state
                assert (state['network'], state['mode'], state['load_nodes'], state['budget']) == (
                    network.name, mode, list(network.load_nodes), budget)
                save_scan(path, state['validation'])
            reference_grid = scan_reference(network, budget, divisions, result['axis_bounds'], path,
                axis_lower=result['axis_lower'], mode=mode, threads=threads, workers=scan_workers,
                progress=monitor.scanning, force_rescan=force_rescan)
            monitor.validation(reference_grid, result)
        return result

    return monitor.execute(calculate, show_ui=show_ui)


def main(case=None, load_nodes=None, divisions=None, *, dimension=None, seconds=CASE_TIME_LIMIT,
         show=SHOW_UI, scan=True, reference=None, output=OUTPUT, mode=mode, force_rescan=None):
    case = {FourBus: 'fourbus', Case33: 'case33'}[NETWORK] if case is None else case
    dimension = DIMENSION if dimension is None else dimension
    load_nodes = ({'fourbus': (1, 2, 3), 'case33': (18, 25, 30)}[case][:dimension]
                  if load_nodes is None else tuple(load_nodes))
    network_type, budget = {'fourbus': (FourBus, BUDGET), 'case33': (Case33, Case33.switch_budget)}[case]
    divisions = (DIVISIONS if case == 'fourbus' and len(load_nodes) == 2 else SCAN_DIVISIONS[len(load_nodes)]) if divisions is None else divisions
    return run(network_type(load_nodes=load_nodes), budget=budget, divisions=divisions,
        show_ui=show, output=recording_path(case, output, load_nodes, mode=mode), tau=REGION_TAU,
        time_limit=seconds, threads=SOLVER_THREADS, scan=scan, reference=reference, scan_workers=SCAN_WORKERS,
        force_rescan=FORCE_RESCAN if force_rescan is None else force_rescan, scan_output=SCAN_OUTPUT, mode=mode)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='正负功率逐网架规划域与离线回放')
    parser.add_argument('--case', choices=('fourbus', 'case33', 'both'),
                        default={FourBus: 'fourbus', Case33: 'case33'}[NETWORK])
    parser.add_argument('--mode', type=int, choices=(0, 1), default=mode)
    parser.add_argument('--dimension', type=int, choices=(2, 3), default=DIMENSION)
    parser.add_argument('--load-nodes', type=lambda value: tuple(map(int, value.split(','))))
    parser.add_argument('--seconds', type=float, default=CASE_TIME_LIMIT)
    parser.add_argument('--divisions', type=int)
    parser.add_argument('--no-ui', action='store_true', help='不打开原生窗口，仍保存回放')
    parser.add_argument('--no-scan', action='store_true', help='只构域，不启动独立扫描')
    parser.add_argument('--force-rescan', action='store_true', default=FORCE_RESCAN)
    parser.add_argument('--reference', type=Path, help='导入相同配置的新主线回放中的扫描')
    parser.add_argument('--replay', nargs='?', const=True, type=Path, help='只回放；可指定记录路径')
    args = parser.parse_args()
    cases = ('fourbus', 'case33') if args.case == 'both' else (args.case,)
    if args.replay:
        paths = ([args.replay] if args.replay is not True else
                 [recording_path(case, load_nodes=args.load_nodes, dimension=args.dimension, mode=args.mode) for case in cases])
        if len(paths) > 1:
            SynchronizedReplay(paths).root.mainloop()
        else:
            monitor = RunMonitor()
            monitor.load_recording(paths[0])
            monitor.replay()
    else:
        for case in cases:
            main(case, args.load_nodes, args.divisions, dimension=args.dimension, seconds=args.seconds,
                 show=SHOW_UI and not args.no_ui, scan=not args.no_scan, reference=args.reference,
                 mode=args.mode, force_rescan=args.force_rescan)
