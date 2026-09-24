"""Certify unresolved grid points by an independent global column upper bound."""
import argparse
import json
import os
from pathlib import Path
from time import perf_counter, time

import gurobipy as gp
import numpy as np
from threadpoolctl import threadpool_limits

from tests.case33_3d import PlanningOracle3D, PAD_KW
from tests.audit_case33_3d import infeasibility_mask
from tests.case33_two_stage import save, fingerprints
from Network.case33bw import Case33
from model import evaluation_bounds


def recheck(directory, index, threads=2, time_limit=1200.):
    started_at, start = time(), perf_counter()
    scan = json.loads((directory/'scan.json').read_text(encoding='utf-8'))
    grid = np.load(directory/'scan_grid.npz')
    coordinates, states = grid['coordinates'], grid['states']
    pending = np.where(states[index] == 0)[0]
    assert len(pending) and scan['independent'] and scan['source_hashes'] == fingerprints()
    target = np.array([coordinates[index], coordinates[pending[0]], scan['fixed_value']])
    oracle = PlanningOracle3D(scan['budget'], threads=threads, seed=71)
    proved = infeasibility_mask(scan, coordinates, oracle.bounds, oracle.square_kw)
    upper = float(coordinates[np.where((coordinates > target[1]) & proved[index])[0][0]])
    oracle.add_support((0., 1., 0.), upper)
    answer = oracle.solve(weights=(0, 1, 0), fixed={0: target[0], 2: target[2]},
                          time_limit=.1, gap_kw=.01, warm_count=0, label='column_setup')

    # This copy supplies a feasible start only. Its bound is never exported.
    witness = next(c['witness_x'] for c in scan['columns'] if c['index'] == index)
    assert witness is not None
    model = oracle.problem.model
    warm = model.copy()
    variables = warm.getVars()
    for k, value in zip(oracle.x_indices, witness):
        variables[k].LB = variables[k].UB = value
    warm.Params.BestBdStop = gp.GRB.INFINITY
    warm.Params.TimeLimit, warm.Params.Threads = 15., 1
    warm.optimize()
    if warm.SolCount:
        values = np.asarray(warm.getAttr('X', variables))
        accepted, certificate = oracle.accept(values, float(warm.MaxVio), 'column_fixed_plan_start')
        print(f'WARM p25={values[oracle.p_indices[1]]:.6f} accepted={accepted}', flush=True)
        if accepted:
            model.setAttr('Start', model.getVars(), values)
    warm.dispose()
    assert np.all(oracle.problem.x.LB == 0) and np.all(oracle.problem.x.UB == 1)
    model.Params.BestBdStop = (target[1]-.01)/oracle.network.base
    model.Params.TimeLimit = time_limit
    last_report = [-60.]

    def progress(problem, where):
        if where == gp.GRB.Callback.MIP:
            runtime = problem.cbGet(gp.GRB.Callback.RUNTIME)
            if runtime-last_report[0] >= 60.:
                last_report[0] = runtime
                bound = problem.cbGet(gp.GRB.Callback.MIP_OBJBND)*oracle.network.base
                print(f'COLUMN index={index} elapsed={runtime:.0f}s upper={bound:.6f} kW', flush=True)

    model.optimize(progress)
    answer.update(label=f'column_recheck_{index}', status=int(model.Status),
                  bound=None, objective=None, feasible=False,
                  infeasible=model.Status == gp.GRB.INFEASIBLE, global_solve=True,
                  seconds=perf_counter()-start, nodes=float(model.NodeCount))
    if not answer['infeasible']:
        if np.isfinite(model.ObjBound):
            answer['bound'] = float(model.ObjBound*oracle.network.base+PAD_KW)
        if model.SolCount:
            values = np.asarray(model.getAttr('X', model.getVars()))
            accepted, certificate = oracle.accept(values, float(model.MaxVio), answer['label'])
            answer.update(certificate, feasible=accepted, objective=float(model.ObjVal*oracle.network.base))
    event = {k: v for k, v in answer.items() if k != 'state'}
    destination = directory/f'column_recheck_{index}.json'
    save(destination, dict(budget=scan['budget'], target=target, upper=upper, events=[event], certificates=oracle.certificates,
                           answer=answer, independent=True, source_hashes=fingerprints(),
                           seconds=perf_counter()-start, started_at=started_at, finished_at=time()))
    print({k: answer[k] for k in ('label','status','bound','objective','feasible','seconds')}, flush=True)
    oracle.close()


def merge_rechecks(directory, paths):
    scan_path = directory/'scan.json'
    scan = json.loads(scan_path.read_text(encoding='utf-8'))
    grid = np.load(directory/'scan_grid.npz')
    coordinates, states = grid['coordinates'], grid['states'].copy()
    assert scan['independent'] and scan['source_hashes'] == fingerprints()
    base_path = directory/'support_recheck_timing_base.json'
    if not base_path.exists():
        save(base_path,dict(seconds=scan['seconds'],finished_at=scan_path.stat().st_mtime))
    base = json.loads(base_path.read_text(encoding='utf-8'))
    labels = {e['label'] for e in scan['events']}
    stamps = {tuple(c['x'])+tuple(np.round(c['p'],6)) for c in scan['certificates']}
    intervals = []
    for path in paths:
        proof = json.loads(path.read_text(encoding='utf-8'))
        assert proof['independent'] and proof['source_hashes'] == scan['source_hashes']
        assert proof['budget'] == scan['budget'] and proof['target'][2] == scan['fixed_value']
        for event in proof['events']:
            assert event['global_solve'] and event['mode'] == 'support'
            assert event['weights'] == [0.,1.,0.] and '1' not in event['fixed']
            assert event['fixed']['2'] == scan['fixed_value'] and event['label'] not in labels
            scan['events'].append(event)
            labels.add(event['label'])
            if event['bound'] is not None:
                for column in scan['columns']:
                    if coordinates[column['index']] >= event['fixed']['0']:
                        column['upper'] = min(column['upper'],event['bound'])
        for certificate in proof['certificates']:
            stamp = tuple(certificate['x'])+tuple(np.round(certificate['p'],6))
            if stamp not in stamps:
                scan['certificates'].append(certificate)
                stamps.add(stamp)
        intervals.append((max(base['finished_at'],proof['started_at']),proof['finished_at']))
    net = Case33(upgrade_count=8)
    proved = infeasibility_mask(scan,coordinates,evaluation_bounds(net,threads=1),net.power_limit)
    assert not np.any(proved & (states == 1))
    states[(states == 0) & proved] = -1
    assert not np.any((states == -1) & ~proved)
    elapsed, previous_end = 0., base['finished_at']
    for begin,end in sorted(intervals):
        elapsed += max(0.,end-max(begin,previous_end))
        previous_end = max(previous_end,end)
    scan.update(seconds=base['seconds']+elapsed,
                counts={str(k):int(np.sum(states == k)) for k in (-1,0,1)},
                status='classified_grid' if np.all(states != 0) else 'unknown_grid_points',
                column_recheck_files=[os.path.relpath(p,directory) for p in paths])
    np.savez_compressed(directory/'scan_grid.npz',coordinates=coordinates,states=states)
    save(scan_path,scan)
    print(scan['status'],scan['counts'],flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', type=Path)
    operation = parser.add_mutually_exclusive_group(required=True)
    operation.add_argument('--index', type=int)
    operation.add_argument('--merge', type=Path, nargs='+')
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--time-limit', type=float, default=1200.)
    args = parser.parse_args()
    with threadpool_limits(limits=1):
        if args.merge:
            merge_rechecks(args.directory,args.merge)
        else:
            recheck(args.directory, args.index, args.threads, args.time_limit)
