"""FourBus: stop at the first SP cut removing less than 0.1% of projected volume.

python experiments/fourbus_outer_volume.py --output results/fourbus_outer_volume
The volume calculation includes every budget-feasible scheme and counts overlap once.
"""
import argparse
import json
from pathlib import Path
import sys
from time import perf_counter

from gurobipy import GRB
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from model import MasterProblem, PLANNING_TOL
from plot import union_volume
from region import clip_polytope, initial_polytope


def budget_schemes(equations, budget, threads=4):
    """Enumerate x only for exact small-case volume measurement, not for G search."""
    problem = MasterProblem(equations, budget=budget, cuts_only=True, threads=threads)
    schemes = []
    with problem.model:
        problem.power.UB = 0.
        problem.model.setObjective(0.)
        while True:
            problem.model.optimize()
            if problem.model.Status == GRB.INFEASIBLE:
                break
            if problem.model.Status != GRB.OPTIMAL or problem.model.MaxVio > PLANNING_TOL:
                raise RuntimeError(f'Scheme enumeration status={problem.model.Status}')
            x = np.rint(problem.x.X).astype(int)
            schemes.append(x)
            problem.exclude(x)
    return np.array(schemes)


class ProjectedOuterVolume:
    """Continuous 3-D polytope union volume in kW^3; no certified inner region."""

    def __init__(self, schemes, axis_bounds, total_bound, *, cuts=()):
        self.schemes = np.asarray(schemes)
        self.axis_bounds = np.asarray(axis_bounds)
        poly = initial_polytope(self.axis_bounds, total_bound, self.axis_bounds)
        self.polytopes = [poly.copy() for _ in self.schemes]
        self.outer_volume = union_volume(self.polytopes)*np.prod(self.axis_bounds)
        for cut in cuts:
            self.add_cut(np.asarray(cut))

    def add_cut(self, cut):
        """Clip every scheme, then measure (V_before-V_after)/V_before."""
        started = perf_counter()
        dimension = len(self.axis_bounds)
        previous = self.outer_volume
        coefficient = cut[1:1+dimension]*self.axis_bounds
        constants = cut[0]+self.schemes @ cut[1+dimension:]
        self.polytopes = [clip_polytope(poly, constant, coefficient)
                          for poly, constant in zip(self.polytopes, constants)]
        self.outer_volume = union_volume(self.polytopes)*np.prod(self.axis_bounds)
        assert self.outer_volume <= previous*(1.+1e-8)
        return dict(outer_volume=float(self.outer_volume),
                    volume_reduction=float(previous-self.outer_volume),
                    volume_reduction_ratio=float((previous-self.outer_volume)/previous),
                    volume_seconds=perf_counter()-started)


def run_volume_comparison(output, *, threshold=.001, seconds=300., threads=4, time_limit=5.):
    """Run both existing algorithms sequentially with identical volume stopping."""
    from experiments.fourbus_outer_partition import initialize, run_trial
    # 1. 共用完整 SOCP 方向 MP2 初始化。
    initialization = initialize(20000., threads)
    output = Path(output)
    results = {}
    # 2. 按同一条件依次运行，避免两方法争用求解线程。
    for method in ('bilinear', 'partition'):
        # 3. 原全局搜索和 SP 取割不变；每条割之后计算所有方案外域的并集体积。
        result = run_trial(method, initialization, seconds=seconds, threads=threads,
                           time_limit=time_limit, volume_threshold=threshold, output=output/method)
        # 4. 首次减少比例 <0.1% 即停；体积计算已包含在 seconds 中。
        results[method] = result
        print(json.dumps(dict(method=method, status=result['status'], cuts=len(result['cuts']),
                              seconds=result['seconds'], volume_seconds=result['volume_seconds'],
                              schemes=result['scheme_count'], outer_volume=result['outer_volume'])), flush=True)
    # 5. 保存独立结果；后续 SOCP 扫描仅作误差评价，不作为停止条件。
    (output/'runs.json').write_text(json.dumps(results, indent=2), encoding='utf-8')
    return results


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default='results/fourbus_outer_volume')
    parser.add_argument('--threshold', type=float, default=.001)
    parser.add_argument('--seconds', type=float, default=300.)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--time-limit', type=float, default=5.)
    args = parser.parse_args()
    run_volume_comparison(args.output, threshold=args.threshold, seconds=args.seconds,
                          threads=args.threads, time_limit=args.time_limit)
