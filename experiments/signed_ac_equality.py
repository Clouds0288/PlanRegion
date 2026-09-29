"""独立比较 SOCP 与 AC 电流等式；不修改主线，不生成非凸模型的对偶割。

python experiments/signed_ac_equality.py --case both --seconds 60
"""
import argparse
import json
from pathlib import Path
import sys
from time import perf_counter

import numpy as np
from gurobipy import GRB
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from experiments.signed_power_compare import (
    Case33, FourBus, SignedPhysics, SignedMaster, DIRECTIONS)
from model import PLANNING_TOL


def build_equality_problem(equations, **options):
    """1. 保留全部网架及运行约束，仅恢复支路电流等式。"""
    problem = SignedMaster(equations, **options)
    model, operation = problem.model, problem.operation
    model.update()
    model.remove([row for row in model.getQConstrs()
                  if row.QCName.startswith('current_cone[')])
    for e, k in equations.keys:
        i = equations.ends[e][0]
        model.addQConstr(operation.P[e, k]**2+operation.Q[e, k]**2
                         == operation.v[i]*operation.ell[e, k],
                         name=f'current_equality[{e},{k}]')
    model.ModelName = 'planning_signed_ac'
    model.Params.NonConvex = 2
    return problem


def solve_trial(equations, method, *, budget, seconds, threads, **options):
    """2. 同一方向或定点下比较两种物理模型；超时不作不可行判断。"""
    started = perf_counter()
    builder = build_equality_problem if method == 'ac' else SignedMaster
    problem = builder(equations, budget=budget, threads=threads, **options)
    with problem.model as model:
        model.Params.TimeLimit = seconds
        model.Params.MIPGapAbs = 1e-7
        model.optimize()
        statuses = {GRB.OPTIMAL: 'optimal', GRB.INFEASIBLE: 'infeasible',
                    GRB.TIME_LIMIT: 'time_limit'}
        if model.Status not in statuses:
            raise RuntimeError(f'{method}: solver status {model.Status}')
        answer = dict(status=statuses[model.Status], seconds=perf_counter()-started)
        if model.Status == GRB.INFEASIBLE:
            return answer
        answer['bound'] = float(model.ObjBound*problem.objective_scale)
        if model.SolCount:
            if model.MaxVio > PLANNING_TOL:
                raise RuntimeError(f'{method}: constraint violation {model.MaxVio}')
            net, state = equations.network, problem.state.X
            P, Q, ell, v = (state[part] for part in (
                equations.P_slice, equations.Q_slice, equations.ell_slice, equations.v_slice))
            upstream = np.array([1. if equations.ends[e][0] == net.root
                else v[net.node_index[equations.ends[e][0]]] for e, k in equations.keys])
            answer.update(p=problem.power.X.tolist(), x=np.rint(problem.x.X).astype(int).tolist(),
                objective=float(model.ObjVal*problem.objective_scale),
                current_gap_max=float(np.max(np.abs(upstream*ell-P**2-Q**2))),
                voltage_min=float(np.sqrt(v.min())), voltage_max=float(np.sqrt(v.max())),
                loss_kw=float(net.r@ell*net.base))
        return answer


def run_case(case, *, seconds=60., threads=4, output=ROOT/'results/signed_ac_equality'):
    """3. 八方向完整规划；4. 同点检验；5. 保存精简结果。"""
    network = FourBus(load_nodes=(1, 2)) if case == 'fourbus' else Case33(load_nodes=(18, 25))
    budget = 20000. if case == 'fourbus' else network.switch_budget
    study_limit = min(network.source_pmax, network.source_smax)*network.base
    equations = SignedPhysics(network, np.full(2, -study_limit), np.full(2, study_limit))
    data = dict(case=case, load_nodes=network.load_nodes, budget=budget,
                power_lower=equations.power_lower.tolist(), power_upper=equations.power_upper.tolist(),
                seconds_limit=seconds, threads=threads, directions=[], probes=[])
    probes = ([-30., 20.], [-150., 0.], [0., -150.]) if case == 'fourbus' else (
        [90., 420.], [-1000., 0.], [-10000., 0.], [0., -10000.], [-10000., -10000.])
    with threadpool_limits(limits=1):
        for direction in DIRECTIONS:
            row = dict(direction=direction.tolist())
            for method in ('socp', 'ac'):
                row[method] = solve_trial(equations, method, budget=budget, seconds=seconds,
                                         threads=threads, direction=direction)
                print(case, direction.tolist(), method, row[method], flush=True)
            data['directions'].append(row)
        for power in probes:
            row = dict(p=power)
            for method in ('socp', 'ac'):
                row[method] = solve_trial(equations, method, budget=budget, seconds=seconds,
                                         threads=threads, power=power)
                print(case, power, method, row[method], flush=True)
            data['probes'].append(row)
    target = Path(output)/case/'comparison.json'
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(data, indent=2, allow_nan=False), encoding='utf-8')
    return data


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=('fourbus', 'case33', 'both'), default='both')
    parser.add_argument('--seconds', type=float, default=60.)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--output', type=Path, default=ROOT/'results/signed_ac_equality')
    args = vars(parser.parse_args())
    case = args.pop('case')
    for name in ('fourbus', 'case33') if case == 'both' else (case,):
        run_case(name, **args)
