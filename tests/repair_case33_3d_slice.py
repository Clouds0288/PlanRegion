"""Resolve only unknown scan points with the unstrengthened original MISOCP."""
import argparse
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from threadpoolctl import threadpool_limits

from tests.case33_3d import PlanningOracle3D
from tests.case33_two_stage import save,fingerprints
from vertify import ACPowerFlow


def repair(directory,threads=4,time_limit=300.):
    scan = json.loads((directory/'scan.json').read_text(encoding='utf-8'))
    grid = np.load(directory/'scan_grid.npz')
    coordinates,states = grid['coordinates'],grid['states'].copy()
    assert scan['source_hashes'] == fingerprints() and scan['independent']
    if not np.any(states == 0):
        print('No unknown grid points.',flush=True)
        return True
    start = perf_counter()
    oracle = PlanningOracle3D(scan['budget'],threads=threads,seed=97,strengthen=False)
    oracle.certificates = [{**c,**{k:np.asarray(c[k]) for k in ('x','p','state')}} for c in scan['certificates']]
    for c in oracle.certificates:
        key = tuple(c['x'])
        if key not in oracle.trees:
            oracle.trees[key] = ACPowerFlow(oracle.network.tree(c['x']))
    columns = {column['index']:column for column in scan['columns']}
    checks = scan.setdefault('original_model_checks',[])
    baseline_seconds = scan['seconds']
    baseline_events = list(scan['events'])
    for i,j in np.argwhere(states == 0):
        if states[i,j] != 0:
            continue
        p = np.array([coordinates[i],coordinates[j],scan['fixed_value']])
        witness = next((key for key,tree in oracle.trees.items() if tree.classify(p)[0] == 1),None)
        answer = None
        if witness is None:
            answer = oracle.solve(fixed=dict(enumerate(p)),time_limit=time_limit,
                                  label=f'original_point_{i}_{j}',warm_count=8)
            if answer['infeasible']:
                assert not np.any(states[i:,j:] == 1)
                states[i:,j:] = -1
            elif answer['feasible']:
                witness = tuple(answer['x'])
        if witness is not None:
            indices = np.where((states[i] == 0) & (coordinates <= coordinates[j]))[0]
            batch = np.column_stack((np.full(len(indices),coordinates[i]),coordinates[indices],
                                     np.full(len(indices),scan['fixed_value'])))
            accepted = oracle.trees[witness].classify(batch) == 1
            states[i,indices[accepted]] = 1
            old_top = columns[i].get('lower')
            new_top = float(coordinates[np.where(states[i] == 1)[0][-1]])
            if old_top is None or new_top > old_top:
                columns[i].update(lower=new_top,witness_x=witness)
        checks.append(dict(index=[int(i),int(j)],p=p,classification=int(states[i,j]),
                           answer=None if answer is None else {k:v for k,v in answer.items() if k not in ('state','x','p')},
                           strengthened=False,source_hashes=fingerprints()))
        scan.update(certificates=oracle.certificates,events=baseline_events+oracle.events,
            seconds=baseline_seconds+perf_counter()-start,columns=list(columns.values()),
            counts={str(k):int(np.sum(states == k)) for k in (-1,0,1)},
            status='classified_grid' if np.all(states != 0) else 'unknown_grid_points')
        np.savez_compressed(directory/'scan_grid.npz',coordinates=coordinates,states=states)
        save(directory/'scan.json',scan)
        print(f'ORIGINAL CHECK p={p.tolist()} classification={states[i,j]} remaining={np.sum(states == 0)}',flush=True)
    oracle.close()
    return bool(np.all(states != 0))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('directory',type=Path)
    parser.add_argument('--threads',type=int,default=4)
    parser.add_argument('--time-limit',type=float,default=300.)
    args = parser.parse_args()
    with threadpool_limits(limits=1):
        complete = repair(args.directory,args.threads,args.time_limit)
    raise SystemExit(0 if complete else 1)
