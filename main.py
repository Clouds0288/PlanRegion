"""规划域主线（方法 RB：径向锥夹逼 + 逐网架支撑查询）：参数、计算入口与原生回放入口。"""
import argparse
from pathlib import Path
from types import SimpleNamespace
import numpy as np

from Network.four_bus_five_corridor import FourBus
from Network.case33bw import Case33, CURRENT_LIMIT
from region import build_region
from monitor import RunMonitor, SynchronizedReplay
from vertify import (ac_network, scan_path, scan_ac_reference, reference_box,
                     import_ac_reference, export_comparison, AC_CACHE_METHOD)

# 运行设置：DIMENSION 同时控制计算、扫描及回放的维数。
mode = 1                          # 符号分区：负光伏至正负荷
NETWORK = Case33                  # FourBus / Case33
DIMENSION = 3                     # 2 / 3；Case33: (18,25) / (18,25,30)
BUDGET = 20000.                   # FourBus 建设预算；Case33 用其 switch_budget
CASE_TIME_LIMIT = 300             # 构域总时限（墙钟，各分区并行），含初始化和记录，不含事后扫描
WORKERS = 16                      # 并行分区进程数；分区内 OBBT 线程数为 WORKERS//分区数
SOLVER_THREADS = 1                # 每个构域求解器的线程数（分区已并行）
REGION_TAU = .005                 # 径向精度；体积目标 ε=d·tau
DISCOVERY_EPS = .15               # A 阶段（网架发现）的放宽体积目标 ε_A
DISCOVERY_SHARE = .25             # A 阶段最多占分区时限的比例
MIP_SECONDS = 60.                 # 单次锥 MISOCP 的时限上限
MIP_GAP = 1e-3                    # 锥 MISOCP 的相对间隙
MIN_WIDTH = {2: 1e-4, 3: 2e-3}    # 锥角直径下限（rad），更窄的锥不再细分
MAX_CONES = {2: 256, 3: 2048}     # 每个分区的叶锥数上限
DIVISIONS = 160                   # FourBus 二维每轴扫描格数
SCAN_DIVISIONS = {2: 160, 3: 80}  # Case33 及 FourBus 三维每轴扫描格数
SCAN_WORKERS = 20                 # 扫描进程数；并行时每个求解器用 1 个线程
FORCE_RESCAN = False              # True 备份后重扫；False 分别复用 AC/SOCP 点并补齐缺失部分
SHOW_UI = True
ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT/'results'/'mainline'
SCAN_OUTPUT = ROOT/'results'/'scan'


def recording_path(case, output=OUTPUT, load_nodes=None, dimension=DIMENSION, mode=mode):
    nodes = ({'fourbus': (1, 2, 3), 'case33': (18, 25, 30)}[case][:dimension]
             if load_nodes is None else load_nodes)
    return Path(output)/f"mode_{mode}"/f"{case}_{'_'.join(map(str, nodes))}.json.gz"


def run(network, *, budget=BUDGET, divisions=DIVISIONS, show_ui=SHOW_UI, output=None,
        tau=REGION_TAU, time_limit=CASE_TIME_LIMIT, workers=WORKERS, threads=SOLVER_THREADS,
        scan=True, reference=None, scan_workers=SCAN_WORKERS,
        force_rescan=FORCE_RESCAN, scan_output=SCAN_OUTPUT):
    monitor = RunMonitor(output=output, algorithm='RB · 径向夹逼 + 网架支撑')
    ac = ac_network(network)
    d = len(network.load_nodes)
    settings = SimpleNamespace(threads=threads, tau=tau, discovery_eps=DISCOVERY_EPS, discovery_share=DISCOVERY_SHARE,
                               mip_seconds=MIP_SECONDS, mip_gap=MIP_GAP, min_width=MIN_WIDTH[d],
                               max_cones=MAX_CONES[d], network_eps=d*tau/2.)   # ε_B=ε/2：网架的 vol(O_x)/vol(P_x)-1

    def calculate():
        # 1. 符号分区并行构域，保存过程回放
        result = build_region(network, budget=budget, monitor=monitor, seconds=time_limit, workers=workers,
                              settings=settings)
        monitor.save()
        # 2. 独立 AC/SOCP 扫描：复用或导入同配置缓存，扫描框覆盖 SOCP 全局界与外包络
        if scan or reference is not None:
            path = scan_path(ac, budget, scan_output, mode)
            if reference is not None and not force_rescan:
                source = Path(reference)
                if source.suffix != '.npz':
                    previous = RunMonitor()
                    previous.load_recording(source)
                    validation = previous.state.get('validation', {})
                    if validation.get('method') != AC_CACHE_METHOD or not validation.get('cache_path'):
                        raise ValueError('The recording does not identify an independent AC cache')
                    source = Path(validation['cache_path'])
                    if not source.is_absolute():
                        source = ROOT/source
                import_ac_reference(ac, budget, source, path, mode=mode)
            lower, upper = reference_box(ac, budget, mode=mode, output=scan_output)
            reference_grid = scan_ac_reference(ac, budget,
                dict(axis_lower=np.minimum(lower, result['axis_lower']),
                     bounds=np.maximum(upper, result['axis_bounds']), shape=(divisions,)*len(lower)),
                path, mode=mode, workers=scan_workers,
                progress=monitor.scanning, force_rescan=force_rescan)
            # 3. 内域（主指标）与外包络相对 AC 的逐格 FR/MR，导出同坐标对比
            monitor.validation(reference_grid, result)
            if output is not None:
                destination = Path(output).parent/(Path(output).name.removesuffix('.json.gz')+'_comparison')
                export_comparison(ac, budget, result, reference_grid, destination)
        return result

    return monitor.execute(calculate, show_ui=show_ui)


def main(case=None, load_nodes=None, divisions=None, *, dimension=None, seconds=CASE_TIME_LIMIT,
         show=SHOW_UI, scan=True, reference=None, output=OUTPUT, force_rescan=None):
    case = {FourBus: 'fourbus', Case33: 'case33'}[NETWORK] if case is None else case
    dimension = DIMENSION if dimension is None else dimension
    load_nodes = ({'fourbus': (1, 2, 3), 'case33': (18, 25, 30)}[case][:dimension]
                  if load_nodes is None else tuple(load_nodes))
    network_type, budget = {'fourbus': (FourBus, BUDGET), 'case33': (Case33, Case33.switch_budget)}[case]
    divisions = (DIVISIONS if case == 'fourbus' and len(load_nodes) == 2 else SCAN_DIVISIONS[len(load_nodes)]) if divisions is None else divisions
    network = network_type(load_nodes=load_nodes, **({'current_limit': CURRENT_LIMIT} if case == 'case33' else {}))
    return run(network, budget=budget, divisions=divisions,
        show_ui=show, output=recording_path(case, output, load_nodes), tau=REGION_TAU,
        time_limit=seconds, workers=WORKERS, threads=SOLVER_THREADS, scan=scan, reference=reference, scan_workers=SCAN_WORKERS,
        force_rescan=FORCE_RESCAN if force_rescan is None else force_rescan, scan_output=SCAN_OUTPUT)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='符号分区规划域（RB）与离线回放')
    parser.add_argument('--case', choices=('fourbus', 'case33', 'both'),
                        default={FourBus: 'fourbus', Case33: 'case33'}[NETWORK])
    parser.add_argument('--dimension', type=int, choices=(2, 3), default=DIMENSION)
    parser.add_argument('--load-nodes', type=lambda value: tuple(map(int, value.split(','))))
    parser.add_argument('--seconds', type=float, default=CASE_TIME_LIMIT)
    parser.add_argument('--divisions', type=int)
    parser.add_argument('--output', type=Path, default=OUTPUT, help='求解结果及默认回放路径的根目录')
    parser.add_argument('--no-ui', action='store_true', help='不打开原生窗口，仍保存回放')
    parser.add_argument('--no-scan', action='store_true', help='只构域，不启动独立扫描')
    parser.add_argument('--force-rescan', action='store_true', default=FORCE_RESCAN)
    parser.add_argument('--reference', type=Path, help='导入同配置的 AC/SOCP 区域缓存或包含该缓存路径的新记录')
    parser.add_argument('--replay', nargs='?', const=True, type=Path, help='只回放；可指定记录路径')
    args = parser.parse_args()
    cases = ('fourbus', 'case33') if args.case == 'both' else (args.case,)
    if args.replay:
        paths = ([args.replay] if args.replay is not True else
                 [recording_path(case, args.output, load_nodes=args.load_nodes, dimension=args.dimension)
                  for case in cases])
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
                 force_rescan=args.force_rescan, output=args.output)
