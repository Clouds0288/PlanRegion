"""四节点网络与可规划域科研插图；所有优化计算使用 Gurobi。"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Rectangle
from old_learning.planning_domain_legacy import ElectricalModel, build_master, convex_hull_cost

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 8,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.8,
    "legend.frameon": False,
})


def save_figure(fig, output: Path, name: str):
    fig.savefig(output / f"{name}.png", dpi=300)
    fig.savefig(output / f"{name}.svg")
    fig.savefig(output / f"{name}.pdf")


def plot_network(model: ElectricalModel, output_dir: str | Path):
    """现状三条 L 型线路及两条候选走廊；坐标仅表达拓扑，不代表地理位置。"""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    positions = {0: (0.0, 1.5), 1: (2.0, 1.5), 2: (2.0, -0.1), 3: (4.0, 1.5)}
    existing_color, candidate_color = "#405767", "#A16F43"
    fig, ax = plt.subplots(figsize=(7.2, 3.6), layout="constrained")
    for edge in model.corridors:
        start, end = np.array(positions[edge.start]), np.array(positions[edge.end])
        color = existing_color if edge.existing else candidate_color
        ax.plot([start[0], end[0]], [start[1], end[1]], color=color,
                linewidth=1.6, linestyle="-" if edge.existing else (0, (5, 3)), zorder=1)
        delta = end - start
        normal = np.array([-delta[1], delta[0]]) / np.linalg.norm(delta)
        label = (start + end) / 2 + (0.16 if edge.existing else -0.16) * normal
        angle = np.degrees(np.arctan2(delta[1], delta[0]))
        ax.text(*label, f"{edge.length_m:g} m", color=color, ha="center", va="center",
                rotation=angle, rotation_mode="anchor", fontsize=8)
    for node, (px, py) in positions.items():
        symbol = (Rectangle((px-0.12, py-0.12), 0.24, 0.24) if node == 0
                  else Circle((px, py), 0.12))
        symbol.set(facecolor="white", edgecolor=existing_color, linewidth=1.3, zorder=3)
        ax.add_patch(symbol)
        ax.text(px, py, str(node), ha="center", va="center", fontsize=8, zorder=4)
        if node == 0:
            label = f"Root · {model.config.transformer_kva:g} kVA\n1 p.u."
        else:
            label = f"{model.config.loads_kw[node]:g}λ kW"
        ax.text(px, py + (-0.30 if node == 2 else 0.28), label,
                ha="center", va="top" if node == 2 else "bottom", fontsize=8)
    ax.legend(handles=[
        Line2D([], [], color=existing_color, linewidth=1.6, label="Existing L-type line"),
        Line2D([], [], color=candidate_color, linewidth=1.6, linestyle=(0, (5, 3)),
               label="Candidate new corridor"),
    ], loc="lower center", ncols=2, handlelength=3.0, columnspacing=2.5, fontsize=8)
    ax.set(xlim=(-0.55, 4.55), ylim=(-1.0, 2.35), aspect="equal")
    ax.axis("off")
    save_figure(fig, output, "00_network_topology")
    return fig


def plot_frontier(result: dict, output_dir: str | Path):
    output = Path(output_dir); output.mkdir(parents=True, exist_ok=True)
    frontier = result["frontier"]
    xp = np.r_[0.0, frontier.lambda_max.to_numpy()]
    yp = np.r_[0.0, frontier.cost_cny.to_numpy()/1000]
    grid = np.linspace(0, float(frontier.lambda_max.max()), 160)
    hull = convex_hull_cost(frontier, grid)/1000
    fig, ax = plt.subplots(figsize=(7.2, 4.5), layout="constrained")
    ax.step(xp, yp, where="pre", color="#3C728F", linewidth=2, label="Exact integer planning boundary")
    ax.fill_between(xp, yp, 68, step="pre", color="#3C728F", alpha=0.09, label="Feasible planning region")
    ax.plot(grid, hull, "--", color="#A16F43", linewidth=1.5, label="Convex-hull lower envelope")
    recovered = result["discovered_frontier"]
    ax.scatter(recovered.lambda_max, recovered.cost_cny/1000, marker="x", color="#3C728F", s=45,
               label="Recovered with integer master + dual cuts", zorder=5)
    example = result["summary"]
    a=example["nonconvex_counterexample_lambda"]; c=example["convex_hull_lower_cost_at_counterexample"]/1000
    ax.scatter([a], [c], marker="s", color="#A16F43", s=30, zorder=6)
    ax.set(xlabel="Uniform load multiplier $\\lambda$", ylabel="Construction budget (thousand CNY)",
           xlim=(0,2.48), ylim=(-1,68))
    ax.grid(True, alpha=.25); ax.legend(loc="upper left", fontsize=8.5)
    save_figure(fig, output, "01_planning_domain")
    return fig


def plot_convergence(result: dict, output_dir: str | Path):
    output = Path(output_dir); output.mkdir(parents=True, exist_ok=True)
    trace = result["cold_result"].trace
    fig, ax=plt.subplots(figsize=(7.2,3.9), layout="constrained")
    ax.plot(trace.cut_count_before, trace.upper_bound, marker="o", color="#A16F43",
            label="MILP master upper bound on maximum multiplier")
    ax.plot(trace.cut_count_before, trace.lower_bound, marker="s", color="#3C728F",
            label="Electrically verified feasible multiplier")
    ax.set(xlabel="Number of accumulated dual feasibility cuts", ylabel="Maximum load multiplier")
    ax.grid(True,alpha=.25);ax.legend(fontsize=8.7)
    save_figure(fig, output, "02_cut_convergence")
    return fig


def relaxed_master_capacity(model, cuts, budget):
    """只含累计电气割的外近似；线型变量仍为二元变量。"""
    master, _, _ = build_master(model, "max_lambda", budget, cuts)
    master.optimize()
    bound = master.ObjBound
    master.dispose()
    return bound


def plot_cut_snapshots(result: dict, output_dir: str | Path):
    output=Path(output_dir);output.mkdir(parents=True,exist_ok=True)
    model=result["model"];f=result["frontier"];pool=result["cuts"]
    # 这里只对比主问题外近似；前沿精确重建使用的是 ε-constraint 递推，不是此绘图采样。
    budgets=np.unique(np.r_[np.arange(0,60001,5000),f.cost_cny.to_numpy(),np.maximum(0,f.cost_cny.to_numpy()-20)])
    fig,ax=plt.subplots(figsize=(7.2,4.3),layout="constrained")
    rows=[]
    stages = (0, result["cold_result"].new_cuts, len(pool))
    for count, color in zip(stages, ("#969DA2", "#A16F43", "#3C728F")):
        values=np.array([relaxed_master_capacity(model,pool[:count],float(b)) for b in budgets])
        ax.step(budgets/1000,values,where="post",linestyle="--",color=color,zorder=4,
                label=f"Projected integer master with {count} cuts")
        rows += [{"cuts":count,"budget_cny":float(b),"master_lambda_upper":float(v)} for b,v in zip(budgets,values)]
    ax.step(np.r_[f.cost_cny.to_numpy(),60000]/1000,np.r_[f.lambda_max.to_numpy(),f.lambda_max.iloc[-1]],
            where="post",color="#263B47",linewidth=2,zorder=3,label="Exact feasible capacity")
    ax.set(xlabel="Construction budget (thousand CNY)",ylabel="Maximum load multiplier",ylim=(.4,3.15))
    ax.grid(True,alpha=.25);ax.legend(fontsize=8.5)
    pd.DataFrame(rows).to_csv(output/"cut_snapshot_values.csv",index=False)
    save_figure(fig, output, "03_cut_snapshots")
    return fig
