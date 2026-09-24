"""Resolve pending grid points using the ORIGINAL unstrengthened integer model.

Use --collect during the scan (only a separate check file is written), or run
without it after the scan has finished to merge the independently proved labels.
"""
import argparse
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from threadpoolctl import threadpool_limits

from tests.case33_two_stage import SliceOracle, save, fingerprints
from vertify import ACPowerFlow


def check(directory, collect=False, time_limit=600.):
    start = perf_counter()
    scan = json.loads((directory/'scan.json').read_text(encoding='utf-8'))
    if not collect:
        assert scan.get('status') in ('classified_grid', 'unknown_grid_points'), 'Finish the scan before merging.'
    grid = np.load(directory/'scan_grid.npz')
    coordinates, states = grid['coordinates'], grid['states'].copy()
    completed = scan['columns'][-1]['index']+1
    pending = np.argwhere(states[:completed] == 0)
    path = directory/'original_model_checks.json'
    checks = json.loads(path.read_text(encoding='utf-8')) if path.exists() else []
    saved = {tuple(item['index']): item for item in checks}
    oracle = SliceOracle(threads=20, seed=29, strengthen=False)
    oracle.certificates = [{**c, **{k: np.asarray(c[k]) for k in ('x', 'p', 'state')}}
                           for c in scan['certificates']]
    for i, j in pending:
        if states[i, j] != 0:
            continue
        existing = saved.get((int(i), int(j)))
        if existing is not None and existing['classification'] != 0:
            record = existing
        else:
            power = np.array([coordinates[i], coordinates[j], 0.])
            answer = oracle.solve(weights=(0., 0.), power=power, time_limit=time_limit,
                                  label=f'original_grid_{i}_{j}')
            ac_feasible = (answer['x'] is not None and
                ACPowerFlow(oracle.network.tree(answer['x'])).classify(power)[0] == 1)
            classification = -1 if answer['infeasible'] else 1 if ac_feasible else 0
            record = dict(index=[int(i), int(j)], p=power, classification=classification,
                x=answer['x'], answer=answer, strengthened=False, source_hashes=fingerprints())
            checks.append(record)
            save(path, checks)
        if record['classification'] == -1:
            assert not np.any(states[i:, j:] == 1), 'Infeasibility proof contradicts a feasible witness.'
            states[i:, j:] = -1
        elif record['classification'] == 1:
            key = np.asarray(record['x'])
            second = coordinates[:j+1]
            power = np.zeros((len(second), 3))
            power[:, 0], power[:, 1] = coordinates[i], second
            passed = ACPowerFlow(oracle.network.tree(key)).classify(power)
            assert np.all(passed == 1)
            states[i, :j+1] = 1
            column = next(c for c in scan['columns'] if c['index'] == i)
            column.update(lower=float(coordinates[j]), witness_x=key)
    if not collect:
        scan['original_model_checks'] = checks
        scan['counts'] = {str(k): int(np.count_nonzero(states == k)) for k in (-1, 0, 1)}
        scan['status'] = 'classified_grid' if np.all(states != 0) else 'unknown_grid_points'
        scan['original_model_check_seconds'] = sum(c['answer']['seconds'] for c in checks)
        np.savez_compressed(directory/'scan_grid.npz', coordinates=coordinates, states=states)
        save(directory/'scan.json', scan)
    oracle.close()
    print(json.dumps(dict(collect=collect, checked=len(checks), elapsed=perf_counter()-start,
                         pending=int(np.count_nonzero(states[:completed] == 0)))))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', type=Path)
    parser.add_argument('--collect', action='store_true')
    parser.add_argument('--time-limit', type=float, default=600.)
    args = parser.parse_args()
    with threadpool_limits(limits=1):
        check(args.directory, args.collect, args.time_limit)
