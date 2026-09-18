"""独立数学验证；--export 同时保存完整过程、科研图件及逐轮 HTML 回放。"""
import argparse
import json
from pathlib import Path

import numpy as np

import planning_domain_demo as demo
from node_power_view import PowerProcess, clip_polytope, simplex, save_records, export_playback, maximal_polytopes, union_surface
from node_power_replay import surface_mesh
from scipy.spatial import ConvexHull
import itertools


def verify(output, export=False):
    model = demo.build_model()
    designs = demo.enumerate_designs(model)
    assert len(designs) == 216
    assert model.W.shape == (53, 8) and model.D.shape == (53, 3)
    assert np.all(model.scales > 0)  # 零基准不能让功率松弛尺度消失。
    physical = [demo.design_inequalities(model, z.x) for z in designs]
    capacities = [demo.evaluate_design(model, z.x)[0] for z in designs]
    for target in ([0.,0.,0.], [40.,35.,25.], [70.,10.,5.], [0.,0.,140.]):
        exact = min((z.cost_cny for z,(G,b) in zip(designs,physical)
                     if np.all(G@target <= b+1e-9)), default=float('inf'))
        solved = demo.solve_by_cuts(model, 'min_cost', target)
        assert solved.value == exact
    for budget in (0.,9460.,14300.,20000.,30000.,32999.,33000.,40000.,40700.):
        solved = demo.solve_by_cuts(model,'max_load',budget)
        exact = max(v for z,v in zip(designs,capacities) if z.cost_cny<=budget)
        assert abs(solved.value-exact)<2e-6
        sp = demo.feasibility_oracle(model,solved.x,solved.p)
        assert sp.violation<=model.config.feasibility_tolerance
    model.subproblem.dispose()

    # 新建模型，验证默认冷启动全过程以及重复 step 的持续状态。
    model = demo.build_model()
    process = PowerProcess(model)
    for _ in range(10):
        process.step()
    assert len(process.history)==10 and not process.finished
    while not process.finished:
        process.step()
    print(f'完整过程：{len(process.history)} 轮，{len(process.cuts)} 条割。',flush=True)

    worst_gap, worst_stationarity = 0., 0.
    for record in process.history:
        assert record['status'] != 'infeasible'
        sp = record['sp']; pi = sp.cut.pi
        gap = abs(pi@sp.rhs+sp.violation)
        stationarity = np.max(np.abs(model.W.T@pi))
        worst_gap, worst_stationarity = max(worst_gap,gap), max(worst_stationarity,stationarity)
        assert pi.min()>=-1e-9 and pi@model.scales<=1+1e-8
        assert stationarity<1e-8 and gap<1e-8
        assert np.max(sp.residual-model.scales*sp.violation)<1e-7
        if record['cut_id']:
            cut = process.cuts[record['cut_id']-1]
            assert abs(cut.constant+cut.x_coeff@record['x']+cut.load_coeff@record['p']+sp.violation)<1e-8

    # 在每个设计的整个真实功率多面体上检验所有割，而非仅检查一条倍率射线。
    worst_margin, worst_domain = float('inf'), 0.
    P = np.array([c.load_coeff for c in process.cuts])
    for i, design in enumerate(process.designs):
        G,b = demo.design_inequalities(model,design.x)
        true_poly = simplex(model.config.power_limit)
        for row,bound in zip(G,b):
            true_poly = clip_polytope(true_poly,bound,-row)
        offsets = np.array([c.constant+c.x_coeff@design.x for c in process.cuts])
        margin = np.min(offsets[:,None]+P@true_poly.T)
        worst_margin = min(worst_margin,float(margin))
        assert margin>=-1e-7
        if i in process.affordable:
            domain_error = np.max(G@process.polytopes[i].T-b[:,None])
            worst_domain = max(worst_domain,float(domain_error))
            assert domain_error<3e-6  # 覆盖整个外近似：仿射行的最坏点一定在顶点。
            assert all((i,tuple(np.round(p,8))) in process.checked for p in process.polytopes[i])

    recovered = sorted(process.frontier,key=lambda row:row['cost_cny'])
    reference = process.reference.copy()
    reference['total_kw'] = np.minimum(reference.total_kw,process.frontier_limit)
    reference = reference.drop_duplicates('total_kw',keep='first')
    assert np.allclose([r['cost_cny'] for r in recovered],reference.cost_cny,atol=1e-6)
    assert np.allclose([r['total_kw'] for r in recovered],reference.total_kw,atol=2e-6)
    save_records(process,output)
    panels = []
    for budget,indices in zip(process.budgets,process.budget_designs):
        outer_capacity = max(process.polytopes[i].sum(axis=1).max() for i in indices)
        exact_capacity = max(capacities[i] for i in indices)
        assert abs(outer_capacity-exact_capacity)<2e-6
        polys = maximal_polytopes([process.polytopes[i] for i in indices])
        # 独立容斥求并集体积，检查交互网格没有跨过非凸边界填面。
        volume = 0.
        for size in range(1,len(polys)+1):
            for group in itertools.combinations(polys,size):
                intersection = group[0]
                for poly in group[1:]:
                    for row in ConvexHull(poly).equations:
                        intersection = clip_polytope(intersection,-row[3],-row[:3])
                volume += (-1)**(size+1)*ConvexHull(intersection).volume
        mesh = surface_mesh(union_surface(polys))
        vertices = np.array(mesh['vertices'])
        rendered = sum(abs(np.linalg.det(vertices[t]))/6 for t in mesh['triangles'])
        print(f'预算 {budget}：并集体积 {volume:.9f}，网格体积 {rendered:.9f}，相对误差 {abs(rendered-volume)/volume:.3g}',flush=True)
        assert abs(rendered-volume)/volume<1e-7
        panels.append(dict(budget_cny=budget,designs=len(indices),maximum_kw=exact_capacity,
                           union_volume=volume,mesh_relative_error=abs(rendered-volume)/volume))
    for small,large in zip(process.budget_designs,process.budget_designs[1:]):
        assert set(small)<=set(large) # 预算增加时域嵌套；所有面板使用同一组逐方案多面体。
    report = dict(iterations=len(process.history), cuts=len(process.cuts), designs=len(process.designs),
                  budgets=process.budgets, panels=panels, certified_designs=len(process.affordable),
                  verified_final_vertices=sum(len(process.polytopes[i]) for i in process.affordable),
                  frontier_steps=len(recovered), worst_duality_gap=worst_gap,
                  worst_stationarity=worst_stationarity, worst_cut_margin=worst_margin,
                  worst_domain_constraint_error=worst_domain,
                  complete=True, ac_power_flow_validated=False)
    output = Path(output)
    (output/'verification.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2),flush=True)
    if export:
        export_playback(process,output)
    model.subproblem.dispose()
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output',default='notebook_results/node_power/complete')
    parser.add_argument('--export',action='store_true')
    args = parser.parse_args()
    verify(args.output,args.export)
