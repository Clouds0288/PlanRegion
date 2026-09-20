"""对已保存的黄色负荷点做 AC 约束消融；不修改原模型与验证标签。

从根向外的正负荷网络中，零电流起步的单调迭代给出最小损耗解。
先求该完整潮流，再分别检查电压、线路有功、配变容量；避免把首次
触发的提前拒绝条件误当作唯一原因。放宽限制后仍遍历全部允许方案。

运行：python results/vertify/diagnose_causes.py
"""
from hashlib import sha256
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from Network.four_bus_five_corridor import network
from vertify import AC_TOL, FIXED_POINT_TOL, ACPowerFlow, BenchmarkResult


def power_flow(design, points, max_iterations=1200):
    """求完整 AC 潮流；1=等式收敛，-1=不存在正电压解，0=未确定。

不在电压下限、线路容量或配变容量处提前停止，以保留全部越限信息。
若迭代电压上界已不为正，则任何正电压 AC 解都不可能存在。
"""
    ell = np.zeros((len(points), design.network.n))
    status = np.zeros(len(points), dtype=np.int8)
    active = np.arange(len(points))
    for _ in range(max_iterations):
        if not len(active):
            break
        P, Q, v, u = design.state(points[active], ell[active])
        bad = np.min(v, axis=1) <= 0
        residual = np.max(np.abs(P*P+Q*Q-u*ell[active]), axis=1)
        good = (~bad) & (residual <= FIXED_POINT_TOL)
        status[active[bad]] = -1
        status[active[good]] = 1
        keep = ~(bad | good)
        active = active[keep]
        ell[active] = (P[keep]**2+Q[keep]**2)/u[keep]
    return status, ell


def violations(design, points, ell):
    """位标记：1=电压下限，2=线路送端有功，4=配变送端视在功率。"""
    P, Q, v, _ = design.state(points, ell)
    c = design.network
    low_voltage = np.any(v < c.vmin-AC_TOL, axis=1)
    line_overload = np.max(P-c.capacity, axis=1) > AC_TOL
    source_squared = (P[:, c.roots].sum(axis=1)**2
                      + Q[:, c.roots].sum(axis=1)**2)
    transformer_overload = source_squared > c.source_smax**2+AC_TOL
    return (low_voltage.astype(np.uint8) + 2*line_overload.astype(np.uint8)
            + 4*transformer_overload.astype(np.uint8))


def diagnose():
    """记录放宽各组约束后每个点最便宜的可行建设方案，再按预算取并集。"""
    folder = ROOT/"results/vertify"
    source = ROOT/"results"/network.name/"result.npz"
    reference = BenchmarkResult.load(source.parent)
    divisions = reference.metadata["divisions"]
    linear_labels = reference.labels[3]
    selected = np.flatnonzero(np.any(linear_labels == 1, axis=0).ravel())
    grid_indices = np.column_stack(np.unravel_index(selected, linear_labels.shape[1:]))
    points = (grid_indices+.5)*reference.spacing
    original = linear_labels.reshape(len(reference.budgets), -1)[:, selected]
    designs = [ACPowerFlow(d) for d in network.designs]
    minimum_cost = np.full((8, len(points)), np.inf)
    unknown_cost = np.full_like(minimum_cost, np.inf)
    witness = np.full(minimum_cost.shape, -1, dtype=np.int16)
    start = time.perf_counter()
    unresolved_checks = 0
    for number, design in enumerate(designs):
        active = np.flatnonzero(np.any(np.isinf(minimum_cost), axis=0))
        if not len(active):
            break
        status, ell = power_flow(design, points[active])
        flags = violations(design, points[active], ell)
        unresolved_checks += np.count_nonzero(status == 0)
        for removed in range(8):
            # 保留的约束全部满足；removed=0 为原 AC，removed=7 只要求潮流有解。
            compatible = (flags & (7 ^ removed)) == 0
            feasible = active[(status == 1) & compatible]
            fresh = feasible[np.isinf(minimum_cost[removed, feasible])]
            minimum_cost[removed, fresh] = design.network.cost
            witness[removed, fresh] = number
            unknown = active[(status == 0) & compatible]
            unknown_cost[removed, unknown] = np.minimum(unknown_cost[removed, unknown], design.network.cost)
        if (number+1) % 24 == 0 or number == len(designs)-1:
            print(f"Designs {number+1}/{len(designs)}; elapsed {time.perf_counter()-start:.1f}s",
                  flush=True)

    rows = []
    for row, budget in enumerate(reference.budgets):
        feasible = np.isfinite(minimum_cost) & (minimum_cost <= budget)
        unresolved = (~feasible & np.isfinite(unknown_cost) & (unknown_cost <= budget))
        yellow = original[row] == 1
        assert np.array_equal(feasible[0], (original[row] & 2) > 0), "Original AC labels disagree"
        assert not np.any(unresolved[:, yellow]), "Unresolved ablation labels; do not publish percentages"
        assert not np.any(feasible[0, yellow])
        count = int(yellow.sum())
        singles = feasible[[1, 2, 4]][:, yellow]
        single_pattern = singles[0].astype(int)+2*singles[1]+4*singles[2]
        patterns = np.bincount(single_pattern, minlength=8)
        restore = [int(np.count_nonzero(feasible[m, yellow])) for m in range(8)]
        rows.append(dict(
            budget=None if np.isinf(budget) else float(budget),
            yellow_points=count,
            restored_points_by_removed_mask=restore,
            restored_percent_by_removed_mask=[100*n/count for n in restore],
            single_relaxation_pattern_counts=patterns.tolist(),
            requires_combined_relaxation=int(np.count_nonzero(
                feasible[7, yellow] & ~singles.any(axis=0))),
            no_power_flow_witness=int(np.count_nonzero(~feasible[7, yellow])),
        ))

    report = dict(
        method="Constraint-removal counterfactuals over all affordable construction plans",
        spacing_kw=reference.spacing.tolist(),
        source_grid=source.name,
        source_sha256=sha256(source.read_bytes()).hexdigest(),
        script_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        unique_yellow_points=len(points),
        design_count=len(designs),
        removed_mask={"1": "voltage lower bound", "2": "line active-power limits",
                      "4": "transformer apparent-power limit", "7": "all three groups"},
        percentage_denominator="Yellow points in the corresponding budget; categories may overlap",
        nonlinear_equations="Full AC DistFlow, unchanged in every counterfactual",
        no_independent_ampacity_limit=True,
        unresolved_point_design_checks=int(unresolved_checks),
        unresolved_yellow_counterfactuals=0,
        original_ac_labels_reproduced=True,
        seconds=time.perf_counter()-start,
        results=rows,
    )
    output = folder/f"causes_{divisions}"
    output.with_suffix(".json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    np.savez_compressed(output.with_suffix(".npz"), points=points, original=original,
                        minimum_cost=minimum_cost, witness=witness, unknown_cost=unknown_cost)
    print(json.dumps(report, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    diagnose()
