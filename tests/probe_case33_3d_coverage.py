"""Bounded global coverage probe; does not modify construction checkpoints."""
import argparse
import json
from pathlib import Path
from time import perf_counter

import gurobipy as gp
from gurobipy import GRB
import numpy as np
from threadpoolctl import threadpool_limits

from tests.case33_3d import PlanningOracle3D, downward_facets
from tests.case33_two_stage import save


def probe(directory, threads=2, time_limit=90., expansion_kw=19.9):
    data = json.loads((directory/'result.json').read_text(encoding='utf-8'))
    groups = {}
    for certificate in data['certificates']:
        groups.setdefault(tuple(certificate['x']), []).append(certificate['p'])
    facets = [downward_facets(np.asarray(points)+expansion_kw) for points in groups.values()]
    oracle = PlanningOracle3D(data['budget'], threads=threads)
    for cut in data['supports']:
        oracle.add_support(cut['weights'], cut['bound'])
    model, power = oracle.problem.model, oracle.problem.power
    power.UB = oracle.bounds
    coverage_slack_kw = model.addVar(lb=-oracle.square_kw, ub=100., name='coverage_slack_kw')
    for index, rows in enumerate(facets):
        select = model.addVars(len(rows), vtype=GRB.BINARY, name=f'outside_{index}')
        model.addConstr(select.sum() == 1)
        for face, row in enumerate(rows):
            # A normalized nonnegative face has minimum residual row[3] at 0.
            big_m = 100.-float(row[3])
            residual = gp.quicksum(float(row[i])*power[i].item() for i in range(3))+float(row[3])
            model.addConstr(coverage_slack_kw <= residual+big_m*(1-select[face]))
    model.setObjective(coverage_slack_kw, GRB.MAXIMIZE)
    model.Params.TimeLimit = time_limit
    model.Params.MIPGapAbs = .01
    model.Params.BestObjStop = 1.
    model.Params.OutputFlag = 1
    start = perf_counter()
    model.optimize()
    result = dict(budget=data['budget'], expansion_kw=expansion_kw, schemes=len(facets),
        face_count=sum(map(len, facets)), status=int(model.Status), seconds=perf_counter()-start,
        global_upper_kw=float(model.ObjBound) if np.isfinite(model.ObjBound) else None,
        candidate_slack_kw=float(model.ObjVal) if model.SolCount else None,
        candidate_p=power.X if model.SolCount else None,
        certified=model.Status == GRB.INFEASIBLE or model.ObjBound <= -1e-4)
    save(directory.parent/f'coverage_probe_{data["budget"]:g}.json', result)
    print(json.dumps({k:v for k,v in result.items() if k != 'candidate_p'}), flush=True)
    oracle.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', type=Path)
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--time-limit', type=float, default=90.)
    args = parser.parse_args()
    with threadpool_limits(limits=1):
        probe(args.directory, args.threads, args.time_limit)
