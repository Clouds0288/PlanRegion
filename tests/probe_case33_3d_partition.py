"""A two-family global disjunctive support-cut probe, without plan enumeration."""
import argparse
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from threadpoolctl import threadpool_limits

from tests.case33_3d import PlanningOracle3D, scheme_hulls, Cell, select_cell, hull_distance
from tests.case33_two_stage import save
from region import clip_polytope, polytope_volume


def probe(directory, threads=2, time_limit=90.):
    data = json.loads((directory/'result.json').read_text(encoding='utf-8'))
    hulls = scheme_hulls(data['certificates'])
    cells = [Cell(np.asarray(c['vertices']),c['depth'],c['gap'],tuple(c['scheme'])) for c in data['cells']]
    chosen, distances = select_cell(cells,hulls)
    target = cells[chosen].vertices[distances.argmax()]
    plans = np.asarray([key for key,_ in hulls])
    normals = np.asarray([rows[np.argmax(rows[:,:3] @ target+rows[:,3]),:3] for _,rows in hulls])
    gaps = np.asarray([hull_distance([target],rows)[0] for _,rows in hulls])
    importance = 1./(1.+(gaps-gaps.min())/20.)**2
    best_score, split_index = -1.,None
    for index in range(plans.shape[1]):
        groups = [np.where(plans[:,index] == bit)[0] for bit in (0,1)]
        if any(not len(group) for group in groups):
            continue
        weights = [importance[group].sum() for group in groups]
        means = [np.average(normals[group],axis=0,weights=importance[group]) for group in groups]
        score = weights[0]*weights[1]/sum(weights)*np.linalg.norm(means[0]-means[1])**2
        if score > best_score:
            best_score,split_index = score,index
    assert split_index is not None
    oracle = PlanningOracle3D(data['budget'],threads=threads)
    for cut in data['supports']:
        oracle.add_support(cut['weights'],cut['bound'])
    certificates = [{**c,**{k:np.asarray(c[k]) for k in ('x','p','state')}} for c in data['certificates']]
    cut_normals, thresholds, events = [],[],[]
    start = perf_counter()
    for bit in (0,1):
        indices = np.where(plans[:,split_index] == bit)[0]
        normal = normals[indices[np.argmin(gaps[indices])]]
        oracle.problem.x[split_index].LB = oracle.problem.x[split_index].UB = bit
        oracle.certificates = [c for c in certificates if c['x'][split_index] == bit]
        answer = oracle.solve(weights=normal,time_limit=time_limit,gap_kw=.5,label=f'partition_{split_index}_{bit}')
        events.append({k:v for k,v in answer.items() if k not in ('state','x','p')})
        if answer['infeasible']:
            continue
        if answer['bound'] is None:
            raise RuntimeError('No global upper bound for one partition.')
        cut_normals.append(normal)
        thresholds.append(answer['bound'])
    cut_normals,thresholds = np.asarray(cut_normals),np.asarray(thresholds)
    removed_volume = 0.
    for cell in cells:
        removed = cell.vertices
        for normal,bound in zip(cut_normals,thresholds):
            removed = clip_polytope(removed,-bound,normal)
            if not len(removed):
                break
        removed_volume += polytope_volume(removed)
    result = dict(budget=data['budget'],split_index=split_index,target=target,
        cut_normals=cut_normals,threshold=thresholds,events=events,seconds=perf_counter()-start,
        removes_target=bool(np.all(cut_normals @ target > thresholds)),
        removed_volume=removed_volume,outer_volume=data['stage2_volume'])
    save(directory.parent/f'partition_probe_{data["budget"]:g}.json',result)
    print(json.dumps({k:v for k,v in result.items() if k not in ('cut_normals','threshold','target','events')}),flush=True)
    oracle.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('directory',type=Path)
    args = parser.parse_args()
    with threadpool_limits(limits=1):
        probe(args.directory)
