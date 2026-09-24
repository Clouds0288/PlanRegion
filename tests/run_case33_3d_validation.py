"""Independent section jobs, plus provenance-preserving reuse of the 5 kW slice."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from time import perf_counter, sleep

import numpy as np
from tests.case33_two_stage import fingerprints, save


def reuse_origin_slice(directory):
    start = perf_counter()
    source = Path('results/case33_two_stage/20260924/reference')
    scan = json.loads((source/'scan.json').read_text(encoding='utf-8'))
    assert scan['source_hashes'] == fingerprints() and scan['independent']
    grid = np.load(source/'scan_grid.npz')
    indices = np.where(np.isclose(np.remainder(grid['coordinates'], 10.), 0.) |
                       (grid['coordinates'] == grid['coordinates'][-1]))[0]
    coordinates = grid['coordinates'][indices]
    states = grid['states'][np.ix_(indices, indices)]
    old_columns = {c['index']: c for c in scan['columns']}
    columns = []
    for i, old_index in enumerate(indices):
        if np.any(states[i] == 1):
            column = old_columns[int(old_index)]
            columns.append(dict(index=i, threshold=float(coordinates[i]),
                lower=float(coordinates[np.where(states[i] == 1)[0][-1]]),
                upper=column['upper'], witness_x=column['witness_x']))
    output = directory/'budget_2/slice_0'
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output/'scan_grid.npz', coordinates=coordinates, states=states)
    save(output/'scan.json', dict(schema_version=3, budget=2., fixed_value=0., spacing_kw=10.,
        independent=True, source_hashes=fingerprints(), columns=columns, events=scan['events'],
        certificates=scan['certificates'], status='classified_grid' if np.all(states != 0) else 'unknown_grid_points',
        counts={str(k): int(np.sum(states == k)) for k in (-1,0,1)},
        seconds=perf_counter()-start, reused=True, source_original_seconds=scan['seconds'],
        original_model_checks=scan.get('original_model_checks',[]),
        source_original_check_seconds=scan.get('original_model_check_seconds',0.),
        provenance=dict(path=str(source), grid_sha256=hashlib.sha256((source/'scan_grid.npz').read_bytes()).hexdigest(),
                        source_spacing_kw=5., rule='take every 10 kW coordinate and retain 6855 kW endpoint')))
    print('Reused independent budget=2, p33=0 slice on the 10 kW subgrid.', flush=True)


def run_job(directory, budget, fixed_value, threads, time_limit):
    folder = directory/f'budget_{budget:g}'/f'slice_{fixed_value:g}'
    previous = json.loads((folder/'scan.json').read_text(encoding='utf-8')) if (folder/'scan.json').exists() else None
    if previous and previous['status'] == 'classified_grid':
        return dict(budget=budget, fixed_value=fixed_value, reused_completed=True)
    command = [sys.executable, '-X', 'utf8', '-m', 'tests.case33_3d', '--mode', 'scan',
               '--budget', str(budget), '--fixed-value', str(fixed_value), '--threads', str(threads),
               '--time-limit', str(time_limit), '--output', str(folder)]
    if previous:
        command.append('--resume')
    log = directory/f'slice_budget_{budget:g}_p33_{fixed_value:g}.log'
    with log.open('a', encoding='utf-8') as handle:
        result = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT)
    print(f'SLICE JOB budget={budget:g} p33={fixed_value:g} exit={result.returncode}', flush=True)
    return dict(budget=budget, fixed_value=fixed_value, returncode=result.returncode)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', type=Path)
    parser.add_argument('--budget', type=float, nargs='+', default=[0.,2.,4.])
    parser.add_argument('--workers', type=int, default=2)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--time-limit', type=float, default=120.)
    parser.add_argument('--wait-budget', type=float)
    parser.add_argument('--reuse-only', action='store_true')
    args = parser.parse_args()
    if 2. in args.budget:
        reuse_origin_slice(args.directory)
    if not args.reuse_only:
        if args.wait_budget is not None:
            path = args.directory/f'budget_{args.wait_budget:g}/result.json'
            while True:
                try:
                    status = json.loads(path.read_text(encoding='utf-8'))['status']
                except (FileNotFoundError, json.JSONDecodeError):
                    status = 'running'
                if status == 'certified':
                    break
                if status in ('unresolved_query', 'time_limit'):
                    raise RuntimeError(f'Construction requires attention: {status}')
                sleep(5)
        jobs = [(budget, fixed) for budget in args.budget for fixed in (0.,500.,1000.)]
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            futures = [executor.submit(run_job, args.directory, budget, fixed, args.threads, args.time_limit)
                       for budget, fixed in jobs]
            results = [future.result() for future in as_completed(futures)]
        save(args.directory/('validation_jobs_'+'_'.join(f'{b:g}' for b in args.budget)+'.json'), results)
