"""Cross-algorithm cut audit and independent checks at disagreement points."""
import argparse
import json
from itertools import product
from pathlib import Path
from time import perf_counter

import numpy as np
from threadpoolctl import threadpool_limits

import model
from Network.case33bw import Case33
from region import GEOMETRY_TOL, contains, halfspaces
from tests.benchmark_boundary_search import BOUNDS, write_json


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--validation', required=True)
    parser.add_argument('--points', type=int, default=6)
    args = parser.parse_args()
    directory = Path(args.input)
    data = [(p.parent.name, json.loads(p.read_text(encoding='utf-8'))) for p in sorted(directory.glob('*/result.json'))]
    validation = json.loads(Path(args.validation).read_text(encoding='utf-8'))
    network = Case33(upgrade_count=8)
    certificates = [row for _, result in data for row in result['certificates']]
    certificates.extend(row for ray in validation for row in ray['certificates'])
    certificates.extend(ray['direct']['certificate'] for ray in validation if ray['direct']['certificate'] is not None)
    values = np.array([np.r_[1., row['p'], row['x']] for row in certificates])
    costs = values[:, 4:]@network.cost
    results = []
    with threadpool_limits(limits=1):
        for name, result in data:
            budget = np.inf if result['budget'] is None else result['budget']
            cut = np.asarray(result['joint_cuts'])
            minimum = min((float(np.min(values[j:j+128]@cut.T)) for j in range(0, len(values), 128)), default=None) if len(cut) else None
            eligible = values[costs <= budget+1e-8, 1:4]/BOUNDS
            failures = 0
            for logical in result.get('logical_cuts', []):
                direction = np.asarray(logical['direction'])/BOUNDS
                active = direction > 0.
                failures += int(np.sum(np.all(eligible[:, active] >
                    logical['theta_upper']*direction[active]+GEOMETRY_TOL, axis=1)))
            outside = 0
            if 'outer_vertices' in result:
                outer = np.asarray(result['outer_vertices'])
                outside = sum(not np.any(np.all(point <= outer+GEOMETRY_TOL, axis=1)) for point in eligible)
            results.append(dict(case=name, known_certificates=len(values), budget_eligible=len(eligible),
                                minimum_joint_cut=minimum, logical_failures=failures, outer_failures=outside))
        report = dict(cross_audit=results, disagreement=[])
        write_json(directory/'cross_audit.json', report)
        assert all((r['minimum_joint_cut'] is None or r['minimum_joint_cut'] >= -model.PLANNING_TOL)
                   and r['logical_failures'] == r['outer_failures'] == 0 for r in results), 'Cross-algorithm cut audit failed'
        # Stratified positions in the common evaluation grid, not hand-picked successful points.
        grid = (np.array(list(product(range(32), repeat=3)))+.5)/32
        for budget_key in sorted({str(r['budget']) for _, r in data}):
            rows = [(name, r) for name, r in data if str(r['budget']) == budget_key]
            budget = np.inf if rows[0][1]['budget'] is None else rows[0][1]['budget']
            masks = []
            for _, result in rows:
                mask = np.zeros(len(grid), dtype=bool)
                if 'inner_points' in result:
                    for point in result['inner_points']:
                        mask |= np.all(grid <= point, axis=1)
                else:
                    for region in result['result']['inner']:
                        mask |= contains(grid, halfspaces(np.asarray(region['vertices'])/BOUNDS))
                masks.append(mask)
            masks = np.asarray(masks)
            indices = np.flatnonzero(masks.any(axis=0) & ~masks.all(axis=0))
            selected = indices[np.unique(np.linspace(0, len(indices)-1, min(args.points, len(indices))).astype(int))] if len(indices) else []
            for index in selected:
                power = grid[index]*BOUNDS
                tick = perf_counter()
                problem = model.PlanningModel(model.PlanningEquations(network, 'socp'), power=power, budget=budget, threads=20)
                with problem.model:
                    answer = problem.solve(time_limit=max(0., 30.-(perf_counter()-tick)))
                report['disagreement'].append(dict(budget=budget, p=power, methods_containing=[name for (name, _), mask in zip(rows, masks) if mask[index]],
                    status='infeasible' if answer is None else 'feasible' if answer['feasible'] else 'unknown',
                    seconds=perf_counter()-tick, answer=answer))
                write_json(directory/'cross_audit.json', report)
                print(json.dumps(dict(event='disagreement', budget=budget_key, point=int(index), status=report['disagreement'][-1]['status'])), flush=True)
                assert answer is not None, 'Independent full model contradicts a certified inner point'
    print(json.dumps(dict(event='audit_done', cases=len(results), certificates=len(certificates))), flush=True)


if __name__ == '__main__':
    main()
