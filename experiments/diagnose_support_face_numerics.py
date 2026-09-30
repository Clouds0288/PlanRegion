"""Reproduce saved candidate-SP failures without changing production tolerances."""
import argparse
import json
from pathlib import Path
import sys
from time import perf_counter

import numpy as np
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from model import SubProblem, PortSubProblem, PortPhysics, new_model
from Network.four_bus_five_corridor import FourBus


def diagnose(source, output):
    output.mkdir(parents=True, exist_ok=True)
    rows = []
    with threadpool_limits(limits=1):
        for path in sorted(source.glob('partition_*.json')):
            saved = json.loads(path.read_text(encoding='utf-8'))
            if 'status=13' not in str(saved['coverage']['error']):
                continue
            step = saved['discovery_history'][-1]
            sign = np.asarray(saved['sign'])
            x, power = np.asarray(step['x']), np.asarray(step['p'])
            equations = PortPhysics(FourBus(load_nodes=(1, 2, 3)), sign)
            label = ''.join('+' if s > 0 else '-' for s in sign)
            row = dict(sign=sign.tolist(), x=x.tolist(), power=power.tolist(),
                       real_power=(power*sign).tolist(), inactive_types=int((x == 0).sum()))
            for reduced in (False, True):
                name = 'original' if not reduced else 'fixed_variables'
                oracle = SubProblem(equations, threads=20, numeric_focus=3)
                with new_model('reproduce_candidate_SP', 20) as model:
                    choice, p, eta, operation, fixed = oracle._build(model, x, power)
                    if reduced:
                        for variable, value in zip([*choice.values(), *p.values()], [*x, *power]):
                            variable.LB = variable.UB = float(value)
                        model.remove(fixed)
                        for key, value in zip(equations.keys, x):
                            if not value:
                                for flow in (operation.P, operation.Q, operation.ell):
                                    flow[key].LB = flow[key].UB = 0.
                    model.Params.OutputFlag = 1
                    model.Params.LogToConsole = 0
                    model.Params.LogFile = str(output/f'{label}_{name}.log')
                    model.Params.TimeLimit = 30.
                    model.optimize()
                    item = dict(status=int(model.Status), solution_count=int(model.SolCount),
                        iterations=int(model.BarIterCount), seconds=float(model.Runtime),
                        tolerance=float(model.Params.BarQCPConvTol),
                        matrix_min=float(model.MinCoeff), matrix_max=float(model.MaxCoeff))
                    if model.SolCount:
                        item.update(eta=float(eta.X), max_violation=float(model.MaxVio),
                            cone_violation=max(0., max(float(np.linalg.norm([a.getValue() for a in tail])
                                -head.getValue()) for head, tail in operation.cones)))
                    row[name] = item
            started = perf_counter()
            try:
                answer = PortSubProblem(equations, threads=20, numeric_focus=3).solve(x, power, 30.)
                row['mainline_continuous_form'] = dict(feasible=answer['feasible'], eta=answer['eta'],
                    cut_returned=answer['cut'] is not None, seconds=perf_counter()-started)
            except (RuntimeError, TimeoutError) as exc:
                row['mainline_continuous_form'] = dict(error=str(exc), seconds=perf_counter()-started)
            rows.append(row)
            print(json.dumps(row), flush=True)
    (output/'diagnosis.json').write_text(json.dumps(rows, indent=2), encoding='utf-8')
    return rows


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path,
                        default=ROOT/'results/support_face_test/fourbus_3d_signed_trial/support_1')
    parser.add_argument('--output', type=Path,
                        default=ROOT/'results/support_face_test/physical_diagnosis')
    args = parser.parse_args()
    diagnose(args.source, args.output)
