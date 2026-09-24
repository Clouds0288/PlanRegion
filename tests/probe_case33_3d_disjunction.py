"""Measure globally valid oblique disjunctions at the worst objective cell."""
import argparse
import json
from pathlib import Path
from time import perf_counter

import gurobipy as gp
from gurobipy import GRB
import numpy as np
from threadpoolctl import threadpool_limits

from tests.case33_3d import PlanningOracle3D, scheme_hulls, Cell, select_cell
from tests.case33_two_stage import save
from region import clip_polytope, polytope_volume


def probe(directory, threads=2, time_limit=90.):
    data = json.loads((directory/'result.json').read_text(encoding='utf-8'))
    hulls = scheme_hulls(data['certificates'])
    cells = [Cell(np.asarray(c['vertices']), c['depth'], c['gap'], tuple(c['scheme'])) for c in data['cells']]
    chosen, distances = select_cell(cells, hulls)
    cell = cells[chosen]
    target = cell.vertices[distances.argmax()]
    candidates = []
    for _, rows in hulls:
        residuals = rows[:, :3] @ target+rows[:, 3]
        face = int(residuals.argmax())
        candidates.append((residuals[face], rows[face, :3]))
    candidates.sort(key=lambda v: v[0])
    normals = []
    for residual, normal in candidates:
        if all(np.linalg.norm(normal-old) > .04 for old in normals):
            normals.append(normal)
        if len(normals) == 3:
            break
    cut_normals = np.asarray(normals)
    oracle = PlanningOracle3D(data['budget'], threads=threads)
    for cut in data['supports']:
        oracle.add_support(cut['weights'], cut['bound'])
    model, power = oracle.problem.model, oracle.problem.power
    power.UB = oracle.bounds
    projected_deficit_kw = model.addVar(ub=oracle.square_kw, name='projected_deficit_kw')
    for normal in cut_normals:
        model.addConstr((gp.quicksum(float(normal[i])*power[i].item() for i in range(3))+
                         projected_deficit_kw)/oracle.network.base >= float(normal @ target)/oracle.network.base)
    model.setObjective(projected_deficit_kw/oracle.network.base, GRB.MINIMIZE)
    model.Params.TimeLimit = time_limit
    model.Params.MIPGapAbs = 1./oracle.network.base
    model.Params.OutputFlag = 1
    start = perf_counter()
    model.optimize()
    lower = max(0., float(model.ObjBound*oracle.network.base)-1e-4)
    thresholds = cut_normals @ target-lower
    removed_volume = 0.
    for node in cells:
        removed = node.vertices
        for normal, threshold in zip(cut_normals, thresholds):
            removed = clip_polytope(removed, -threshold, normal)
            if not len(removed):
                break
        removed_volume += polytope_volume(removed)
    result = dict(budget=data['budget'], status=int(model.Status), seconds=perf_counter()-start,
        selected_gap_kw=cell.gap, target=target, cut_normals=cut_normals,
        bound=lower, objective=float(model.ObjVal*oracle.network.base) if model.SolCount else None,
        removed_volume=removed_volume, outer_volume=data['stage2_volume'])
    save(directory.parent/f'disjunction_probe_{data["budget"]:g}.json', result)
    print(json.dumps({k:v for k,v in result.items() if k not in ('cut_normals','target')}), flush=True)
    oracle.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', type=Path)
    args = parser.parse_args()
    with threadpool_limits(limits=1):
        probe(args.directory)
