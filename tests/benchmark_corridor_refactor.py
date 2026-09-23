"""Before/after audit of corridor representation on the current working tree."""
from argparse import ArgumentParser
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from threadpoolctl import threadpool_limits

from Network.case33bw import Case33
from Network.four_bus_five_corridor import FourBus
from main import build_continuous_region
from model import PlanningEquations, PlanningModel, evaluation_bounds
from plot import json_value, sample_region
from vertify import ACPowerFlow


ROOT = Path(__file__).resolve().parents[1]


def run(output, region_seconds):
    sources = ['model.py', 'region.py', 'main.py', 'Network/__init__.py',
               'Network/case33bw.py', 'Network/four_bus_five_corridor.py']
    report = dict(source_hashes={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
                  threads=1, query_limit_seconds=60, region_limit_seconds=region_seconds,
                  queries=[], regions=[])
    output.parent.mkdir(parents=True, exist_ok=True)

    def save():
        output.write_text(json.dumps(json_value(report), indent=2, allow_nan=False), encoding='utf-8')

    with threadpool_limits(limits=1):
        for name, net in [('fourbus', FourBus()), ('case33_8', Case33(candidate_count=8)),
                          ('case33_16', Case33(candidate_count=16))]:
            powers = ([10., 10., 10.], [25., 15., 20.], [30., 30., 30.]) if name == 'fourbus' else (
                [100., 800., 150.], [250., 1800., 350.], [300., 3000., 500.])
            for method in ('linear', 'socp'):
                equations = PlanningEquations(net, method)
                queries = [(f'budget_{b}', dict(budget=b)) for b in net.budgets]
                queries += [(f'power_{i}', dict(power=p)) for i, p in enumerate(powers)]
                for key, options in queries:
                    started = perf_counter()
                    problem = PlanningModel(equations, threads=1, **options)
                    with problem.model:
                        answer = problem.solve(time_limit=60.)
                        binary_count = problem.model.NumBinVars
                    row = dict(case=name, method=method, query=key, seconds=perf_counter()-started,
                               binary_count=binary_count, status='infeasible' if answer is None else answer['status'])
                    if answer is not None:
                        row.update({k: answer[k] for k in ('objective', 'bound', 'feasible', 'p')})
                        if answer['x'] is not None:
                            row['margin'] = equations.margin(answer['x'], answer['p'], answer['state'])
                            row['cost'] = float(equations.investment(answer['x']))
                            choice = equations.choice(answer['x'])
                            row['choice'] = choice
                            oracle = ACPowerFlow(net.design(choice), threads=1)
                            try:
                                row['ac_status'] = int(oracle.classify(answer['p'])[0])
                            finally:
                                oracle.close()
                    report['queries'].append(row)
                    save()
                    print(json.dumps(json_value(row)), flush=True)
            if region_seconds:
                bounds = evaluation_bounds(net, threads=1)
                budget = 20000. if name == 'fourbus' else 2.
                points = (np.indices((3,)*3).reshape(3, -1).T+.5)*bounds/3
                for method in ('linear', 'socp', 'hybrid'):
                    region = build_continuous_region(net, method, budget, bounds, threads=1,
                                                     time_limit=region_seconds)
                    row = dict(case=name, method=method, budget=budget, bounds=bounds,
                               status=region['status'], max_total=region['max_total'],
                               seconds=region['timing']['total_seconds'], counts=region['counts'],
                               grid_labels=sample_region(region, points, bounds))
                    report['regions'].append(row)
                    save()
                    print(json.dumps(json_value(row)), flush=True)
    return report


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--region-seconds', type=float, default=30.)
    args = parser.parse_args()
    run(args.output, args.region_seconds)
