"""FourBus 3D signed support-face construction and paired mainline comparison.

python experiments/compare_support_face_fourbus_3d.py --threads 20 --seconds 300
python experiments/compare_support_face_fourbus_3d.py --algorithm support --no-scan

The support-face method discovers schemes with complete physical MISOCP queries.
Eight power-sign partitions are not construction-scheme enumeration. Both methods
keep the original physical constraints, budget and certification tolerances.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import hashlib
from itertools import product
import json
from pathlib import Path
import sys
from time import perf_counter

import numpy as np
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import main as baseline
from Network.four_bus_five_corridor import FourBus
from model import PortPhysics, port_bounds
from monitor import RunMonitor
from region import GEOMETRY_TOL, initial_polytope
from vertify import scan_path, scan_reference
from experiments.test_support_face_certification_fourbus_2d import (
    SupportReplay, discover_schemes, initial_bounds, json_value)

DEFAULT_OUTPUT = ROOT/'results'/'support_face_test'/'fourbus_3d_signed_physical'


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(json_value(value), ensure_ascii=False, indent=2,
                               allow_nan=False), encoding='utf-8')


def source_hashes():
    paths = ('main.py', 'model.py', 'region.py', 'monitor.py', 'vertify.py',
             'Network/four_bus_five_corridor.py',
             'experiments/test_support_face_certification_fourbus_2d.py',
             'experiments/compare_support_face_fourbus_3d.py')
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths}


def run_support(*, output, seconds=300., threads=20, budget=20000., tau=baseline.REGION_TAU,
                max_iterations=300, max_discovery_iterations=500, verbose=False):
    network = FourBus(load_nodes=(1, 2, 3))
    bounds = port_bounds(network)
    signs = list(product((1, -1), repeat=3))
    monitor = RunMonitor(output=output/'monitor.json.gz', algorithm='支持面认证 · 按需发现')
    monitor.envelopes = {''.join('+' if s > 0 else '-' for s in sign):
        [np.asarray(list(product((0., 1.), repeat=3)))*bounds*sign] for sign in signs}
    monitor._emit('start', phase='初始化', mode=1, network=network.name, load_nodes=network.load_nodes,
        budget=budget, method='socp', algorithm=monitor.algorithm, time_limit=seconds,
        cost_unit=network.cost_unit, initial_plan=network.initial_plan, bounds=bounds,
        axis_bounds=bounds, axis_lower=-bounds, tau=tau, schemes={}, cut_history={},
        global_outer=[p for polys in monitor.envelopes.values() for p in polys])
    partitions = []
    with threadpool_limits(limits=1):
        for j, sign in enumerate(signs):
            started = perf_counter()
            allowance = max(0., (seconds-monitor.timing()['total_seconds'])/(len(signs)-j))
            deadline = started+allowance
            local = RunMonitor(parent=monitor, sign=sign, clock=monitor.clock, algorithm=monitor.algorithm)
            equations = PortPhysics(deepcopy(network), sign)
            replay = SupportReplay(equations.network, budget, tau, GEOMETRY_TOL, output/'unused',
                                   bounds=bounds, monitor=local)
            try:
                _, axis_bounds, total_bound, seed = initial_bounds(equations, budget, threads,
                    allowance, bounds=bounds, replay=replay, deadline=deadline)
                discovery = discover_schemes(equations, seed, bounds, axis_bounds, total_bound,
                    budget=budget, tau=tau, epsilon_geom=GEOMETRY_TOL, threads=threads,
                    support_seconds=allowance, global_seconds=allowance,
                    max_iterations=max_iterations, max_discovery_iterations=max_discovery_iterations,
                    outer_cuts=True, verbose=verbose, replay=replay, deadline=deadline)
            except (RuntimeError, TimeoutError) as exc:
                discovery = dict(rows=[], cuts=[], history=[], complete=False, bound=None,
                    status='TIME_LIMIT' if isinstance(exc, TimeoutError) else 'UNRESOLVED', error=str(exc),
                    global_milp_calls=0, global_physical_calls=0, discovery_sp_calls=0, discovery_cut_lp_calls=0,
                    conditional_support_cuts=0)
            rows = discovery['rows']
            totals = dict(global_certified=discovery['complete'] and all(r['certified'] for r in rows),
                status=discovery['status'], schemes=len(rows), enumeration_mip_calls=0,
                support_calls=sum(r['support_calls'] for r in rows),
                initial_global_mp2_calls=4, global_milp_calls=discovery['global_milp_calls'],
                global_physical_calls=discovery['global_physical_calls'],
                discovery_sp_calls=discovery['discovery_sp_calls'],
                discovery_cut_lp_calls=discovery['discovery_cut_lp_calls'],
                conditional_support_cuts=discovery['conditional_support_cuts'])
            summary = dict(schema='support-face-fourbus-physical-v3', sign=sign, bounds=bounds, tau=tau, measure_unit='kW^3',
                totals=totals, schemes=rows, cuts=discovery['cuts'],
                discovery_history=discovery['history'],
                coverage=dict(mode='physical', **{k: discovery[k] for k in ('complete', 'bound', 'status', 'error')}))
            replay.finish(summary)
            result = local.state['result']
            totals['total_algorithm_seconds'] = perf_counter()-started
            partitions.append(dict(sign=sign, summary=summary, result=result))
            label = ''.join('+' if s > 0 else '-' for s in sign)
            write_json(output/f'partition_{label}.json', summary)
            monitor.save()
            print(f'support 3D {label}: {totals["status"]}, '
                  f'{totals["total_algorithm_seconds"]:.3f}s, schemes={len(rows)}, '
                  f'support={totals["support_calls"]}, physical={totals["global_physical_calls"]}, '
                  f'SP={totals["discovery_sp_calls"]}, error={discovery["error"]}', flush=True)
    certified = all(p['result']['certified'] for p in partitions)
    result = dict(status='certified' if certified else 'incomplete', certified=certified,
        timing=monitor.timing(), inner=[], outer=[],
        counts={key: sum(p['summary']['totals'][key] for p in partitions) for key in
                ('schemes', 'support_calls', 'initial_global_mp2_calls', 'global_milp_calls', 'global_physical_calls',
                 'discovery_sp_calls', 'discovery_cut_lp_calls', 'conditional_support_cuts',
                 'enumeration_mip_calls')})
    for p in partitions:
        for key in ('inner', 'outer'):
            result[key].extend({**row, 'sign': p['sign'],
                               'vertices': np.asarray(row['vertices'])*p['sign']}
                              for row in p['result'][key])
    vertices = np.vstack([row['vertices'] for row in result['outer']])
    result.update(axis_lower=vertices.min(axis=0), axis_bounds=vertices.max(axis=0))
    monitor._emit('region_end', phase='构域完成' if certified else '构域停止', result=result,
                  coverage_complete=certified, active_scheme=None, ray=None, sp_point=None, global_point=None)
    monitor.save()
    write_json(output/'result.json', result)
    return result


def run_comparison(*, output=DEFAULT_OUTPUT, seconds=300., threads=20, repetitions=1,
                   divisions=80, workers=4, scan=True, algorithm='both'):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    before = source_hashes()
    protocol = dict(network='fourbus', load_nodes=[1, 2, 3], mode=1, budget=20000.,
        tau=baseline.REGION_TAU, epsilon_geom=GEOMETRY_TOL, threads=threads, seconds=seconds,
        repetitions=repetitions, divisions_requested=divisions, source_hashes=before,
        timing_scope='end-to-end construction with native recording, excludes scan and final save',
        baseline_settings={k: getattr(baseline, k) for k in
                           ('CUT_THRESHOLD', 'CUT_PATIENCE', 'POINT_TOL', 'RAY_THRESHOLD')})
    write_json(output/'protocol.json', protocol)
    runs = []
    for repeat in range(repetitions):
        order = ('support', 'mainline') if repeat % 2 == 0 else ('mainline', 'support')
        for name in order:
            if algorithm not in ('both', name):
                continue
            destination = output/f'{name}_{repeat+1}'
            destination.mkdir(parents=True, exist_ok=True)
            print(f'RUN {name} {repeat+1}/{repetitions}', flush=True)
            started = perf_counter()
            if name == 'support':
                result = run_support(output=destination, seconds=seconds, threads=threads)
            else:
                result = baseline.run(FourBus(load_nodes=(1, 2, 3)), budget=20000.,
                    show_ui=False, output=destination/'monitor.json.gz', mode=1, scan=False,
                    threads=threads, tau=baseline.REGION_TAU, time_limit=seconds)
                write_json(destination/'result.json', result)
            runs.append(dict(algorithm=name, repeat=repeat+1, output=str(destination),
                wall_seconds=perf_counter()-started, result=result))
            write_json(output/'comparison.json', dict(protocol=protocol, runs=runs))
    after = source_hashes()
    if before != after:
        raise RuntimeError('Source changed during comparison; results are not a paired benchmark')
    if scan:
        network = FourBus(load_nodes=(1, 2, 3))
        lower = np.min([r['result']['axis_lower'] for r in runs], axis=0)
        upper = np.max([r['result']['axis_bounds'] for r in runs], axis=0)
        path = scan_path(network, 20000., ROOT/'results'/'scans', 1)
        reference = scan_reference(network, 20000., divisions, upper, path,
            axis_lower=lower, mode=1, threads=threads, workers=workers, progress=lambda a, b: None)
        protocol.update(scan_path=str(path), scan_shape=reference['states'].shape,
                        scan_lower=reference['axis_lower'], scan_upper=reference['bounds'])
        with threadpool_limits(limits=1):
            for row in runs:
                monitor = RunMonitor(output=Path(row['output'])/'monitor.json.gz')
                monitor.load_recording(monitor.output)
                monitor.validation(reference, row['result'])
                monitor.save()
                row['metrics'] = monitor.validation_state['validation']['metrics']
    write_json(output/'comparison.json', dict(protocol=protocol, runs=runs))
    for row in runs:
        print(json.dumps(json_value({k: v for k, v in row.items() if k != 'result'})), flush=True)
        print(row['algorithm'], row['result']['status'], row['result']['timing'], flush=True)
    return runs


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument('--threads', type=int, default=20)
    parser.add_argument('--seconds', type=float, default=300.)
    parser.add_argument('--repetitions', type=int, default=1)
    parser.add_argument('--divisions', type=int, default=80)
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--algorithm', choices=('both', 'support', 'mainline'), default='both')
    parser.add_argument('--no-scan', action='store_true')
    args = parser.parse_args()
    run_comparison(output=args.output, threads=args.threads, seconds=args.seconds,
        repetitions=args.repetitions, divisions=args.divisions, workers=args.workers,
        algorithm=args.algorithm, scan=not args.no_scan)
