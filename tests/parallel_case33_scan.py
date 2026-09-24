"""Continue an independent grid scan in disjoint column blocks, then merge.

No two-stage cuts or construction certificates are imported. Each worker uses
only the independent scan's own completed prefix and global bound. Six workers
times four Gurobi threads use at most 24 solver threads.
"""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
from time import perf_counter

import numpy as np

from tests.case33_two_stage import ROOT, save


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', type=Path)
    parser.add_argument('--workers', type=int, default=6)
    args = parser.parse_args()
    directory = args.directory.resolve()
    base = json.loads((directory/'scan.json').read_text(encoding='utf-8'))
    grid = np.load(directory/'scan_grid.npz')
    coordinates, initial = grid['coordinates'], grid['states']
    first = base['columns'][-1]['index']+1
    remaining = np.where(np.any(initial[first:] == 0, axis=1))[0]
    if not len(remaining):
        return
    stop = first+int(remaining[-1])+1
    boundaries = np.linspace(first, stop, args.workers+1, dtype=int)
    destination = directory/'parallel_chunks'
    destination.mkdir(exist_ok=False)
    shutil.copyfile(directory/'scan.json', directory/'serial_prefix.json')
    shutil.copyfile(directory/'scan_grid.npz', directory/'serial_prefix_grid.npz')
    jobs = []
    start = perf_counter()
    for worker, (left, right) in enumerate(zip(boundaries[:-1], boundaries[1:])):
        if right <= left:
            continue
        output = destination/f'chunk_{worker:02d}'
        output.mkdir()
        shutil.copyfile(directory/'scan.json', output/'scan.json')
        shutil.copyfile(directory/'scan_grid.npz', output/'scan_grid.npz')
        log = (output/'console.log').open('w', encoding='utf-8')
        command = [sys.executable, '-X', 'utf8', '-m', 'tests.case33_two_stage', '--mode', 'scan',
            '--resume', '--threads', '4', '--time-limit', '120', '--spacing-kw', str(base['spacing_kw']),
            '--start-index', str(left), '--stop-index', str(right), '--output', str(output)]
        process = subprocess.Popen(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT,
                                   creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0))
        jobs.append((process, output, log))
        print(f'Worker {worker}: columns {left}..{right-1}', flush=True)
    save(directory/'parallel_manifest.json', dict(workers=len(jobs), threads_per_worker=4,
        start_index=first, stop_index=stop, boundaries=boundaries,
        prefix_seconds=base['seconds'], prefix_columns=len(base['columns']),
        prefix_events=len(base['events']), prefix_certificates=len(base['certificates'])))
    for process, output, log in jobs:
        code = process.wait()
        log.close()
        if code != 0:
            raise RuntimeError(f'{output.name} exited with {code}; its log and completed points are preserved.')
    parallel_seconds = perf_counter()-start
    merged = initial.copy()
    columns = {c['index']: c for c in base['columns']}
    events, certificates = list(base['events']), list(base['certificates'])
    for _, output, _ in jobs:
        result = json.loads((output/'scan.json').read_text(encoding='utf-8'))
        local = np.load(output/'scan_grid.npz')['states']
        conflict = (merged != 0) & (local != 0) & (merged != local)
        assert not np.any(conflict), f'Contradictory grid certificates from {output.name}'
        merged[local != 0] = local[local != 0]
        for column in result['columns'][len(base['columns']):]:
            columns[column['index']] = column
        events.extend(result['events'][len(base['events']):])
        certificates.extend(result['certificates'][len(base['certificates']):])
    base.update(columns=[columns[k] for k in sorted(columns)], events=events, certificates=certificates,
        seconds=base['seconds']+parallel_seconds, parallel_seconds=parallel_seconds,
        parallel_workers=len(jobs), parallel_threads_per_worker=4,
        counts={str(k): int(np.count_nonzero(merged == k)) for k in (-1, 0, 1)},
        status='classified_grid' if np.all(merged != 0) else 'unknown_grid_points')
    np.savez_compressed(directory/'scan_grid.npz', coordinates=coordinates, states=merged)
    save(directory/'scan.json', base)
    print(json.dumps(dict(status=base['status'], counts=base['counts'], parallel_seconds=parallel_seconds)))


if __name__ == '__main__':
    main()
