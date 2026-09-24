"""Bounded local-support probe; no production or construction result is changed."""
import json
from pathlib import Path
import gurobipy as gp
import numpy as np
from threadpoolctl import threadpool_limits
from region import halfspaces, clip_polytope, polytope_volume
from tests.case33_3d import PlanningOracle3D, scheme_hulls, cell_gap, hull_distance
from tests.case33_two_stage import save


def probe(directory, budget):
    result = json.loads((directory/f'budget_{budget}/result.json').read_text())
    certificates = [{**c, **{k:np.asarray(c[k]) for k in ('x','p','state')}} for c in result['certificates']]
    hulls = scheme_hulls(certificates)
    cell = max(result['cells'],key=lambda c:c['gap'])
    vertices = np.asarray(cell['vertices'])
    before, key, distances = cell_gap(vertices,hulls)
    target = vertices[int(distances.argmax())]
    facets = dict(hulls)[key]
    row = max(facets,key=lambda row:hull_distance([target],row[None,:])[0])
    weights = row[:3]
    oracle = PlanningOracle3D(budget,threads=4)
    oracle.certificates = certificates
    for equation in halfspaces(vertices):
        oracle.problem.model.addConstr((gp.quicksum(float(w)*oracle.problem.power[i].item()
                                        for i,w in enumerate(equation[:3]))+float(equation[3]))/oracle.network.base <= 0.)
    answer = oracle.solve(weights=weights,time_limit=40.,gap_kw=2.,warm_count=6,label=f'local_support_budget_{budget}')
    clipped = np.empty((0,3)) if answer['infeasible'] else (
        clip_polytope(vertices,answer['bound'],-weights) if answer['bound'] is not None else vertices)
    after = cell_gap(clipped,scheme_hulls(oracle.certificates))[0] if len(clipped) else 0.
    data = dict(budget=budget,before_gap=before,after_gap=after,before_volume=polytope_volume(vertices),
                after_volume=polytope_volume(clipped),answer=answer,weights=weights,
                vertices=vertices,certificates=oracle.certificates)
    save(directory/f'local_bound_probe_{budget}.json',data)
    print(json.dumps({k:data[k] for k in ('budget','before_gap','after_gap','before_volume','after_volume')}),flush=True)
    oracle.close()


if __name__ == '__main__':
    with threadpool_limits(limits=1):
        for budget in (0,2,4):
            probe(Path('results/case33_3d/20260924'),budget)
