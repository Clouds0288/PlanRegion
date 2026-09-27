"""FourBus: compare the original bilinear oracle with certified corner partitioning.

python experiments/fourbus_outer_partition.py --method both --seconds 300
The production main.py and the original experimental oracle remain unchanged.
"""
import argparse
import heapq
import json
from pathlib import Path
import sys
from time import perf_counter

import gurobipy as gp
from gurobipy import GRB
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from experiments.fourbus_outer import GlobalViolation
from Network.four_bus_five_corridor import FourBus
from model import GridPhysics, MasterProblem, PLANNING_TOL, SubProblem


class CornerViolation(GlobalViolation):
    """max eta*(x,c), over box corners c and schemes intersecting the current box.

problem.power is an outer-domain witness, not the load at which eta* is evaluated.
Only corner_bits * theta is multiplied, and corner_bits is binary.
"""

    def __init__(self, equations, budget, axis_bounds, total_bound, *, threads=4):
        super().__init__(equations, budget, axis_bounds, total_bound, threads=threads)
        model, template = self.model, self.template
        model.update()
        model.remove([row for row in model.getQConstrs() if row.QCName == 'dual_objective'])
        model.ModelName = 'G_box_corners'
        model.Params.NonConvex = 0
        model.Params.MIPGap = 0.
        model.Params.MIPGapAbs = 1e-7

        # 1. Aggregate load multipliers; normalization supplies exact finite bounds.
        relaxed = template.sp_eta < 0.
        assert not np.any(template.sp_p[~relaxed])
        coefficients = template.sp_p[relaxed]/(-template.sp_eta[relaxed, None])
        lower = np.minimum(0., coefficients.min(axis=0))
        upper = np.maximum(0., coefficients.max(axis=0))
        dimension = len(axis_bounds)
        self.theta = model.addMVar(dimension, lb=lower, ub=upper, name='theta')
        model.addConstr(self.theta == template.sp_p.T @ self.row_dual, name='load_dual')

        # 2. One binary per coordinate selects all box corners without enumeration.
        self.corner_bits = model.addMVar(dimension, vtype=GRB.BINARY, name='corner_bits')
        self.corner_term = model.addMVar(dimension, lb=lower, ub=upper, name='corner_term')
        model.addConstr(self.corner_term >= lower*self.corner_bits)
        model.addConstr(self.corner_term <= upper*self.corner_bits)
        model.addConstr(self.corner_term >= self.theta-upper*(1.-self.corner_bits))
        model.addConstr(self.corner_term <= self.theta-lower*(1.-self.corner_bits))

        # 3. Box bounds only change coefficients in this linear dual objective.
        dual_value = (self.choice_term.sum()-template.sp_rhs @ self.row_dual
                      -template.cone_constant @ self.cone_dual).item()
        self.dual_objective = model.addConstr(self.violation <= dual_value, name='corner_objective')
        self.witness_problem = MasterProblem(equations, budget=budget, cuts_only=True, threads=threads)
        self.witness_problem.power.UB = axis_bounds
        self.witness_problem.model.addConstr(self.witness_problem.power.sum() <= total_bound)
        model.update()

    def add_cut(self, cut):
        self.problem.add_cut(cut)
        self.witness_problem.add_cut(cut)

    def solve_box(self, box_lower, box_upper, inherited_bound, time_limit, epsilon):
        """Return an upper bound for this box; a corner can lie outside its x-slice."""
        model = self.model
        self.problem.power.LB = box_lower
        self.problem.power.UB = box_upper
        for j in range(len(box_lower)):
            model.chgCoeff(self.dual_objective, self.theta[j].item(), -float(box_lower[j]))
            model.chgCoeff(self.dual_objective, self.corner_term[j].item(),
                           -float(box_upper[j]-box_lower[j]))
        # This is a local bound: never reuse the previous sibling's smaller bound.
        self.violation.UB = min(inherited_bound, float(max(box_upper)/self.equations.network.base))
        model.Params.TimeLimit = time_limit
        model.Params.BestBdStop = epsilon
        model.optimize()
        if model.Status == GRB.INFEASIBLE:
            return dict(objective=None, bound=0., status=int(model.Status),
                        seconds=float(model.Runtime), x=None, p=None, corner=None)
        if model.Status not in (GRB.OPTIMAL, GRB.TIME_LIMIT, GRB.USER_OBJ_LIMIT):
            raise RuntimeError(f'Corner G status={model.Status}')
        if model.SolCount and model.IntVio > PLANNING_TOL:
            raise RuntimeError(f'Corner G IntVio={model.IntVio:g} > {PLANNING_TOL:g}')
        # Dual residuals are compensated by candidate_lower, never used as a
        # primal feasibility certificate. witness() independently checks x,p.
        x = np.rint(self.problem.x.X).astype(int) if model.SolCount else None
        power = self.problem.power.X if model.SolCount else None
        corner = (box_lower+(box_upper-box_lower)*np.rint(self.corner_bits.X)
                  if model.SolCount else None)
        return dict(objective=self.candidate_lower(x, corner) if model.SolCount else None,
                    bound=min(float(model.ObjBound), self.violation.UB),
                    status=int(model.Status), seconds=float(model.Runtime),
                    x=x, p=power, corner=corner,
                    dual_max_violation=float(model.MaxVio) if model.SolCount else None)

    def witness(self, x, box_lower, box_upper, time_limit):
        """Fix x and duals; maximize their affine certificate at an actual outer point."""
        problem, model = self.witness_problem, self.witness_problem.model
        problem.x.LB = problem.x.UB = x
        problem.power.LB, problem.power.UB = box_lower, box_upper
        model.setObjective(self.theta.X @ problem.power, GRB.MAXIMIZE)
        model.Params.TimeLimit = time_limit
        model.optimize()
        if model.Status != GRB.OPTIMAL or model.MaxVio > PLANNING_TOL:
            raise RuntimeError(f'Box witness status={model.Status}, MaxVio={model.MaxVio:g}')
        power = problem.power.X
        return dict(x=x, p=power, objective=self.candidate_lower(x, power))


class PolyhedralCornerViolation(CornerViolation):
    """Globally valid linear outer approximation of the same dual SOCs.

MILP objective bounds are upper bounds even before cone refinement converges.
The inherited candidate_lower always checks the original full cones.
"""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.model.remove(self.model.getQConstrs())
        for block in self.template.cone_slices:
            head = self.cone_dual[block.start].item()
            for item in self.cone_dual[block.start+1:block.stop].tolist():
                self.model.addConstr(head >= item)
                self.model.addConstr(head >= -item)
        self.oa_planes = 0
        self.model.update()

    def solve_box(self, box_lower, box_upper, inherited_bound, time_limit, epsilon):
        deadline = perf_counter()+time_limit
        started = perf_counter()
        bound = inherited_bound
        iterations = 0
        while True:
            answer = super().solve_box(box_lower, box_upper, bound,
                                       max(0., deadline-perf_counter()), epsilon)
            iterations += 1
            bound = answer['bound']
            if (answer['x'] is None or bound <= epsilon
                    or self.model.ObjVal-answer['objective'] <= min(1e-6, epsilon/10.)
                    or perf_counter() >= deadline):
                break
            values = self.cone_dual.X
            planes = 0
            for block in self.template.cone_slices:
                tail = values[block.start+1:block.stop]
                length = np.linalg.norm(tail)
                if length > values[block.start]+1e-9:
                    direction = tail/(length*(1.+1e-12))
                    self.model.addConstr(1000.*self.cone_dual[block.start] >=
                                         1000.*direction @ self.cone_dual[block.start+1:block.stop])
                    planes += 1
            self.oa_planes += planes
            if not planes:
                break
        answer.update(bound=bound, seconds=perf_counter()-started,
                      oa_iterations=iterations, oa_planes=self.oa_planes)
        return answer


def initialize(budget, threads):
    """Shared full-SOCP MP2 initialization; all x,p,y are free."""
    started = perf_counter()
    equations = GridPhysics(FourBus(), 'socp')
    dimension = len(equations.network.load_nodes)
    initial = []
    for direction in [*np.eye(dimension), np.ones(dimension)]:
        problem = MasterProblem(equations, budget=budget, direction=direction, threads=threads)
        with problem.model:
            answer = problem.solve(time_limit=60.)
        initial.append(dict(direction=direction.tolist(), x=answer['x'].tolist(),
                            p=answer['p'].tolist(), bound=float(answer['bound'])))
    return dict(budget=budget, initial=initial, initial_seconds=perf_counter()-started,
                axis_bounds=[item['bound'] for item in initial[:-1]], total_bound=initial[-1]['bound'])


def run_trial(method, initialization, *, epsilon=1e-4, seconds=300., time_limit=5.,
              threads=4, cuts=(), output='results/fourbus_outer_partition', volume_threshold=None):
    """1 initialize; 2 global search; 3 certify; 4 SP cut; 5 split/repeat."""
    started = perf_counter()
    deadline = started+seconds
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    equations = GridPhysics(FourBus(), 'socp')
    axis_bounds = np.array(initialization['axis_bounds'])
    total_bound, budget = initialization['total_bound'], initialization['budget']
    result = dict(initialization, method=method, epsilon=epsilon, threads=threads,
                  time_limit=time_limit, wall_limit=seconds, trace=[],
                  cuts=[np.asarray(cut).tolist() for cut in cuts], initial_cut_count=len(cuts),
                  status='running', certified=False, bound=float(max(axis_bounds)/equations.network.base),
                  counts=dict(global_search=0, sp=0, splits=0, pruned=0), seconds=0.)

    # 1. 初始化：两方法共用方向 MP2 外域、完整 SOCP 物理模型与 SP 取割。
    search_type = PolyhedralCornerViolation if method == 'partition' else GlobalViolation
    search = search_type(equations, budget, axis_bounds, total_bound, threads=threads)
    sp = SubProblem(equations, threads=threads)
    volume_monitor = None
    if volume_threshold is not None:
        from experiments.fourbus_outer_volume import ProjectedOuterVolume, budget_schemes
        volume_started = perf_counter()
        volume_monitor = ProjectedOuterVolume(
            budget_schemes(equations, budget, threads), axis_bounds, total_bound, cuts=cuts)
        result.update(volume_threshold=volume_threshold,
                      initial_outer_volume=volume_monitor.outer_volume,
                      outer_volume=volume_monitor.outer_volume,
                      volume_seconds=perf_counter()-volume_started,
                      scheme_count=len(volume_monitor.schemes))
    for cut in cuts:
        if method == 'partition':
            search.add_cut(np.asarray(cut))
        else:
            search.problem.add_cut(np.asarray(cut))
    queue = [(-result['bound'], 0, np.zeros(len(axis_bounds)), axis_bounds.copy())]
    serial, certified_bound = 0, 0.
    step = 0
    with search.model:
        while perf_counter() < deadline:
            step += 1
            record = dict(step=step)
            # 2. 全局选点与网架：双线性 G 搜 x,p；分块 G 搜与当前块相交的 x 及最坏角点。
            if method == 'partition':
                negative_bound, _, box_lower, box_upper = heapq.heappop(queue)
                answer = search.solve_box(box_lower, box_upper, -negative_bound,
                                          min(time_limit, deadline-perf_counter()), epsilon)
                record.update(box_lower=box_lower.tolist(), box_upper=box_upper.tolist(),
                              corner=None if answer['corner'] is None else answer['corner'].tolist())
            else:
                answer = search.solve(min(time_limit, deadline-perf_counter()), epsilon)
            result['counts']['global_search'] += 1
            record.update(objective=answer['objective'], local_bound=answer['bound'],
                          status=answer['status'], global_seconds=answer['seconds'])
            if method == 'partition':
                record.update(oa_iterations=answer['oa_iterations'], oa_planes=answer['oa_planes'],
                              dual_max_violation=answer.get('dual_max_violation'))
            candidate = None

            # 3. 认证：仅上界 <= epsilon 才排清当前块；所有块排清后才全局停止。
            if answer['bound'] <= epsilon:
                record['action'] = 'certify'
                if method == 'partition':
                    certified_bound = max(certified_bound, answer['bound'])
                    result['counts']['pruned'] += 1
                else:
                    result.update(certified=True, status='residual_certified')
            elif perf_counter() < deadline and answer['x'] is not None:
                candidate = (search.witness(answer['x'], box_lower, box_upper,
                                            min(time_limit, deadline-perf_counter()))
                             if method == 'partition' else answer)
                record.update(x=candidate['x'].tolist(), p=candidate['p'].tolist(),
                              candidate_lower=candidate['objective'])

            # 4. 固定真实外域候选 x,p 调用原 SP；仅 y/eta 自由，失败返回含 x,p 的有效割。
            if (candidate is not None and candidate['objective'] > PLANNING_TOL
                    and perf_counter() < deadline):
                checked = sp.solve(candidate['x'], candidate['p'],
                                   time_limit=min(30., deadline-perf_counter()))
                result['counts']['sp'] += 1
                if checked['feasible']:
                    raise RuntimeError('Positive conservative lower bound disagrees with SP')
                cut = checked['cut']
                dimension = len(axis_bounds)
                record['cut_at_candidate'] = float(cut[0]+cut[1:1+dimension]@candidate['p']
                                                    +cut[1+dimension:]@candidate['x'])
                result['cuts'].append(cut.tolist())
                record['action'] = 'cut'
                if method == 'partition':
                    search.add_cut(cut)
                else:
                    search.problem.add_cut(cut)
                if volume_monitor is not None:
                    measurement = volume_monitor.add_cut(cut)
                    record.update(measurement)
                    result['outer_volume'] = measurement['outer_volume']
                    result['volume_seconds'] += measurement['volume_seconds']
                    if measurement['volume_reduction_ratio'] < volume_threshold:
                        result['status'] = 'volume_stagnation'

            # 5. 更新与下一轮：加割后重查同一块；无可切反例则二分最长负荷区间。
            if method == 'partition' and answer['bound'] > epsilon:
                if record.get('action') == 'cut' or perf_counter() >= deadline:
                    serial += 1
                    heapq.heappush(queue, (-answer['bound'], serial, box_lower, box_upper))
                    record.setdefault('action', 'time_limit')
                else:
                    dimension = int(np.argmax(box_upper-box_lower))
                    midpoint = (box_lower[dimension]+box_upper[dimension])/2.
                    left_upper, right_lower = box_upper.copy(), box_lower.copy()
                    left_upper[dimension] = right_lower[dimension] = midpoint
                    for lower, upper in ((box_lower, left_upper), (right_lower, box_upper)):
                        serial += 1
                        heapq.heappush(queue, (-answer['bound'], serial, lower, upper))
                    result['counts']['splits'] += 1
                    record.update(action='split', split_dimension=dimension)
            elif method == 'bilinear' and 'action' not in record:
                record['action'] = 'unresolved'

            if method == 'partition':
                result['bound'] = max(certified_bound, -queue[0][0] if queue else 0.)
                result['active_boxes'] = len(queue)
                if not queue:
                    result.update(certified=True, status='residual_certified')
            else:
                result['bound'] = answer['bound']
            result['seconds'] = perf_counter()-started
            record.update(seconds=result['seconds'], bound=result['bound'], cuts=len(result['cuts']))
            result['trace'].append(record)
            (output/'results.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
            print(f"{method} {step}: {record['action']} cuts={len(result['cuts'])} "
                  f"upper={result['bound']:.7g} boxes={len(queue) if method == 'partition' else 1} "
                  f"t={result['seconds']:.2f}s", flush=True)
            if result['certified'] or result['status'] == 'volume_stagnation':
                break
        if result['status'] == 'running':
            result['status'] = 'time_limit'
        result['seconds'] = perf_counter()-started
    if method == 'partition':
        search.witness_problem.model.dispose()
    (output/'results.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--method', choices=('bilinear', 'partition', 'both'), default='both')
    parser.add_argument('--budget', type=float, default=20000.)
    parser.add_argument('--epsilon', type=float, default=1e-4)
    parser.add_argument('--seconds', type=float, default=300.)
    parser.add_argument('--time-limit', type=float, default=5.)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--output', default='results/fourbus_outer_partition')
    parser.add_argument('--cuts-from')
    parser.add_argument('--volume-threshold', type=float)
    args = parser.parse_args()
    initialization = initialize(args.budget, args.threads)
    cuts = (() if args.cuts_from is None
            else json.loads(Path(args.cuts_from).read_text(encoding='utf-8'))['cuts'])
    methods = ('bilinear', 'partition') if args.method == 'both' else (args.method,)
    for method in methods:
        run_trial(method, initialization, epsilon=args.epsilon, seconds=args.seconds,
                  time_limit=args.time_limit, threads=args.threads, cuts=cuts,
                  output=str(Path(args.output)/method), volume_threshold=args.volume_threshold)
