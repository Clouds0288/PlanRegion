"""两区域负荷倍率的完整三维规划域及真实 SP 对偶割回放。

显示坐标为 (节点1倍率 lambda_A, 节点2/3倍率 lambda_B, 建设成本上限 B)。
几何数据内部仍以 (B, lambda_A, lambda_B) 存储，显示时统一置换。
为完整覆盖两个连续负荷方向，本教学扩展枚举有限的建设方案，
逐一分离当前外近似多边形的顶点；不是原 MP1/MP2 的单次查询轨迹。
候选选择不使用电气真值；所有割来自原电气 LP 的对偶乘子。
原单倍率模型、主问题与 notebook 前六个代码块保持原有行为。
"""
from __future__ import annotations

import itertools
import json
import math
from pathlib import Path

import gurobipy as gp
import numpy as np
import shapely
from gurobipy import GRB
from shapely.geometry import Polygon

import old_learning.planning_domain_legacy as demo


OUTPUT = Path(__file__).resolve().parent / "notebook_results" / "three_dimensional"
GEOM_TOL = 2e-8


def grouped_demands(model):
    """D 的两列之和严格还原原来的 d，功率/电压松弛尺度保持不变。"""
    D = np.zeros((len(model.h), 2))
    demands = np.array([[0, 0], [model.config.loads_kw[1], 0],
                        [0, model.config.loads_kw[2]], [0, model.config.loads_kw[3]]])
    for r, row in enumerate(model.electrical_constraints):
        name = row.ConstrName
        if name.startswith("balance_"):
            _, node, sign = name.split("_")
            D[r] = int(sign) * demands[int(node)]
        elif name == "transformer":
            D[r] = -demands.sum(axis=0)
    np.testing.assert_allclose(D.sum(axis=1), model.d, atol=1e-12)
    return D, demands


def tree_choices(model):
    """只枚举结构合法的 x，不计算或使用其真实承载能力。"""
    nk, ne = len(model.lines), len(model.corridors)
    choices = []
    for edges in itertools.combinations(range(ne), 3):
        reached = {0}
        for _ in range(3):
            for e in edges:
                edge = model.corridors[e]
                if edge.start in reached or edge.end in reached:
                    reached.update((edge.start, edge.end))
        if len(reached) != 4:
            continue
        for types in itertools.product(range(nk), repeat=3):
            x = np.zeros(model.nx)
            for e, k in zip(edges, types):
                x[e * nk + k] = 1
            label = "; ".join(f"{model.corridors[e].name}:{model.lines[k].name}"
                              for e, k in zip(edges, types))
            choices.append((float(model.cost @ x), label, x))
    choices.sort(key=lambda item: (item[0], item[1]))
    return np.array([v[2] for v in choices]), np.array([v[0] for v in choices]), [v[1] for v in choices]


def clip_polygon(vertices, constant, coeff):
    """精确线段求交，保留 constant + coeff @ lambda >= 0 的半平面。"""
    if not len(vertices):
        return np.empty((0, 2))
    values = constant + vertices @ coeff
    if np.all(values >= -1e-12):
        return vertices.copy()
    points = []
    for i, end in enumerate(vertices):
        start = vertices[i - 1]
        va, vb = values[i - 1], values[i]
        inside_a, inside_b = va >= 0, vb >= 0
        if inside_a != inside_b:
            points.append(start + (end - start) * va / (va - vb))
        if inside_b:
            points.append(end)
    cleaned = []
    for p in points:
        if not cleaned or np.linalg.norm(p - cleaned[-1]) > 1e-10:
            cleaned.append(p)
    if len(cleaned) > 1 and np.linalg.norm(cleaned[0] - cleaned[-1]) < 1e-10:
        cleaned.pop()
    return np.array(cleaned).reshape(-1, 2)


def separate_vertex(model, D, x, point):
    rhs = model.h + model.T @ x + D @ point
    model.subproblem.setAttr("RHS", model.electrical_constraints, rhs.tolist())
    model.subproblem.optimize()
    if model.subproblem.Status != GRB.OPTIMAL:
        raise RuntimeError(f"SP 未证明最优: {model.subproblem.Status}")
    pi = -np.array(model.subproblem.getAttr("Pi", model.electrical_constraints))
    violation = float(model.subproblem.ObjVal)
    cut = dict(constant=float(pi @ model.h), x_coeff=(pi @ model.T).tolist(),
               load_coeff=(pi @ D).tolist(), pi=pi.tolist(), violation=violation)
    if violation > model.config.feasibility_tolerance:
        assert pi.min() >= -1e-9
        assert np.max(np.abs(model.W.T @ pi)) < 1e-8
        assert pi @ model.scales <= 1 + 1e-8
        assert abs(pi @ rhs + violation) < 1e-8
    return violation, cut


def physical_polygon(model, x, demands):
    """独立物理校验：树上潮流、线路容量、节点电压、配变容量。"""
    nk = len(model.lines)
    chosen = {e: int(np.argmax(row)) for e, row in enumerate(x.reshape(-1, nk)) if row.sum() > .5}
    adj = {i: [] for i in range(4)}
    for e in chosen:
        edge = model.corridors[e]
        adj[edge.start].append((edge.end, e))
        adj[edge.end].append((edge.start, e))
    parents, order = {}, [0]
    for i in order:
        for j, e in adj[i]:
            if j != 0 and j not in parents:
                parents[j] = (i, e)
                order.append(j)
    downstream, flows = demands.copy(), {}
    for j in reversed(order[1:]):
        i, e = parents[j]
        flows[e] = downstream[j].copy()
        downstream[i] += downstream[j]
    constraints = [(model.config.transformer_kva * model.config.power_factor, -demands.sum(axis=0))]
    drops = {0: np.zeros(2)}
    qr = math.tan(math.acos(model.config.power_factor))
    for j in order[1:]:
        i, e = parents[j]
        line, edge = model.lines[chosen[e]], model.corridors[e]
        coefficient = 2 * (line.r_ohm_km + qr * line.x_ohm_km) * (edge.length_m / 1000) / (1000 * model.config.voltage_kv**2)
        drops[j] = drops[i] + coefficient * flows[e]
        constraints.extend([(line.capacity_kw, -flows[e]), (1 - model.config.voltage_min_pu**2, -drops[j])])
    bound = model.config.lambda_search_max
    poly = np.array([[0, 0], [bound, 0], [bound, bound], [0, bound]], dtype=float)
    for constant, coeff in constraints:
        poly = clip_polygon(poly, constant, coeff)
    return poly, constraints


def polygon_geometry(vertices):
    return Polygon(vertices) if len(vertices) >= 3 else Polygon()


def budget_sections(polygons, costs, budget_max):
    """在所有费用断点上求多边形的并集；绝不以凸包填平凹陷。"""
    slabs = []
    current = Polygon()
    for cost in np.unique(costs):
        additions = [polygon_geometry(polygons[i]) for i in np.flatnonzero(costs == cost)]
        grown = shapely.union_all([current, *additions])
        if not slabs or grown.symmetric_difference(current).area > 1e-10:
            if slabs:
                slabs[-1][1] = float(cost)
            slabs.append([float(cost), float(budget_max), grown])
        current = grown
    return slabs


def region_volume(slabs):
    return sum((hi - lo) * poly.area for lo, hi, poly in slabs)


def run_complete_process(output=OUTPUT, progress=True, max_cuts=2000):
    """穷尽每种设计的外近似顶点，形成有限建设集合上的完整分离证书。"""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    model = demo.build_model()
    try:
        D, demands = grouped_demands(model)
        X, costs, labels = tree_choices(model)
        bound = model.config.lambda_search_max
        budget_max = math.ceil(float(max(costs)) / 1000) * 1000
        square = np.array([[0, 0], [bound, 0], [bound, bound], [0, bound]], dtype=float)
        polygons = [square.copy() for _ in X]
        snapshots, cuts, certified = [[p.copy() for p in polygons]], [], set()
        oracle_calls = 0
        for design_id, x in enumerate(X):
            while True:
                # 不依赖独立物理真值：从当前候选域中选尚未核验的最外侧顶点。
                candidates = sorted(polygons[design_id], key=lambda p: (-sum(p), -p[0]))
                found = False
                for point in candidates:
                    key = (design_id, *np.round(point, 10))
                    if key in certified:
                        continue
                    violation, cut = separate_vertex(model, D, x, point)
                    oracle_calls += 1
                    if violation <= model.config.feasibility_tolerance:
                        certified.add(key)
                        continue
                    if len(cuts) >= max_cuts:
                        raise RuntimeError("达到割数上限，尚未完成全域证明；不保存为完整结果。")
                    cut.update(source_design=design_id, source_load=point.tolist(), source_cost=float(costs[design_id]))
                    cuts.append(cut)
                    offsets = cut["constant"] + X @ np.array(cut["x_coeff"])
                    coeff = np.array(cut["load_coeff"])
                    polygons = [clip_polygon(poly, offset, coeff) for poly, offset in zip(polygons, offsets)]
                    assert all(len(p) >= 3 for p in polygons), "零负荷附近应保持可行"
                    snapshots.append([p.copy() for p in polygons])
                    found = True
                    if progress and (len(cuts) <= 5 or len(cuts) % 25 == 0):
                        print(f"cut={len(cuts)}, design={design_id + 1}/{len(X)}, SP={oracle_calls}", flush=True)
                    break
                if not found:
                    break

        # 后续割可能生成新的顶点：早期已证明整个多边形可行，所以子集也可行。
        # 再逐一验证最终顶点，并由独立树物理公式交叉校验所有 216 个多边形。
        max_eta, max_area_error, max_diagonal_error, max_physical_error = 0., 0., 0., 0.
        exact_polygons = []
        final_vertex_checks = 0
        for i, x in enumerate(X):
            exact, constraints = physical_polygon(model, x, demands)
            exact_polygons.append(exact)
            max_area_error = max(max_area_error, polygon_geometry(polygons[i]).symmetric_difference(polygon_geometry(exact)).area)
            for point in polygons[i]:
                eta, _ = separate_vertex(model, D, x, point)
                max_eta = max(max_eta, eta)
                final_vertex_checks += 1
                max_physical_error = max(max_physical_error, *[-(c + a @ point) for c, a in constraints])
            diagonal = min([bound, *[c / -sum(a) for c, a in constraints if sum(a) < -1e-12]])
            original = min(bound, demo.evaluate_design(model, x).lambda_max)
            max_diagonal_error = max(max_diagonal_error, abs(diagonal - original))
        assert max_eta <= model.config.feasibility_tolerance
        assert max_area_error < GEOM_TOL
        assert max_diagonal_error < 1e-9
        assert max_physical_error < 1e-6
        # 每条割对独立真值的所有顶点有效，故不会误删任一真实可行设计。
        min_valid_margin = 0.
        for cut in cuts:
            offsets = cut["constant"] + X @ np.array(cut["x_coeff"])
            for offset, exact in zip(offsets, exact_polygons):
                min_valid_margin = min(min_valid_margin, float(np.min(offset + exact @ cut["load_coeff"])))
        assert min_valid_margin >= -1e-8
        all_sections = [budget_sections(polys, costs, budget_max) for polys in snapshots]
        volumes = np.array([region_volume(s) for s in all_sections])
        assert np.all(np.diff(volumes) <= 1e-5)
        proof = dict(complete=True, cut_count=len(cuts), design_count=len(X),
                     oracle_calls=oracle_calls, final_vertex_checks=final_vertex_checks,
                     maximum_final_eta=max_eta, maximum_polygon_area_error=max_area_error,
                     maximum_physical_violation=max_physical_error, maximum_diagonal_error=max_diagonal_error,
                     minimum_true_region_cut_margin=min_valid_margin,
                     nested_regions=True, exact_polygon_unions=True,
                     algorithm="Finite design enumeration + SP vertex separation; not a single MP1/MP2 run",
                     unchanged_projection_cuts=int(np.sum(np.abs(np.diff(volumes)) < 1e-5)))
        result = dict(version=1, axes=["Budget B (CNY)", "lambda_A (node 1)", "lambda_B (nodes 2,3)"],
                      loads_kw=demands.tolist(), lambda_max=bound, budget_max_cny=budget_max,
                      costs=costs.tolist(), designs=X.astype(int).tolist(), design_labels=labels,
                      cuts=cuts, volumes=volumes.tolist(), verification=proof)
        (output / "cuts_and_verification.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
        np.savez_compressed(output / "polygon_history.npz", **{f"s{s}_d{i}": p for s, polys in enumerate(snapshots) for i, p in enumerate(polys)})
        (output / "verification.json").write_text(json.dumps(proof, ensure_ascii=False, indent=2), encoding="utf-8")
        if progress:
            print(json.dumps(proof, ensure_ascii=False, indent=2), flush=True)
        return result, snapshots, all_sections
    finally:
        model.subproblem.dispose()


def load_process(output=OUTPUT):
    output = Path(output)
    result = json.loads((output / "cuts_and_verification.json").read_text(encoding="utf-8"))
    with np.load(output / "polygon_history.npz") as history:
        snapshots = [[history[f"s{s}_d{i}"] for i in range(len(result["designs"]))]
                     for s in range(len(result["cuts"]) + 1)]
    sections = [budget_sections(p, np.array(result["costs"]), result["budget_max_cny"]) for p in snapshots]
    return result, snapshots, sections


def mesh_from_sections(sections, end_faces=True):
    """按精确多边形构造三角面，不在预算台阶间插值、不使用三维凸包。"""
    vertices, triangles, lookup = [], [], {}

    def vertex(x, yz):
        p = (round(float(x) / 1000, 8), round(float(yz[0]), 9), round(float(yz[1]), 9))
        if p not in lookup:
            lookup[p] = len(vertices)
            vertices.append(p)
        return lookup[p]

    def cap(x, poly):
        if poly.is_empty:
            return
        for tri in shapely.constrained_delaunay_triangles(poly).geoms:
            triangles.append([vertex(x, p) for p in list(tri.exterior.coords)[:3]])

    previous = Polygon()
    for lo, hi, poly in sections:
        if hi <= lo:
            continue
        # GEOS 并集可能保留共线碎点；仅消除远小于校验容差的数值碎片。
        poly = poly.simplify(1e-10, preserve_topology=True)
        # 预算断面只显示新增/消失的面积，避免重复内部面遮挡。
        if end_faces:
            cap(lo, poly.symmetric_difference(previous))
        components = [poly] if poly.geom_type == "Polygon" else list(poly.geoms)
        for part in components:
            if part.is_empty or part.geom_type != "Polygon":
                continue
            for ring in [part.exterior, *part.interiors]:
                coords = list(ring.coords)
                for a, b in zip(coords, coords[1:]):
                    if np.linalg.norm(np.subtract(a, b)) < 1e-9:
                        continue
                    ids = [vertex(lo, a), vertex(hi, a), vertex(hi, b), vertex(lo, b)]
                    triangles.extend([[ids[0], ids[1], ids[2]], [ids[0], ids[2], ids[3]]])
        previous = poly
    if end_faces and sections:
        cap(sections[-1][1], previous)
    return dict(v=vertices, f=triangles)


def removed_sections(previous, current):
    boundaries = sorted({b for slab in previous + current for b in slab[:2]})
    removed = []
    i = j = 0
    for lo, hi in zip(boundaries, boundaries[1:]):
        while i + 1 < len(previous) and previous[i][1] <= lo:
            i += 1
        while j + 1 < len(current) and current[j][1] <= lo:
            j += 1
        diff = previous[i][2].difference(current[j][2])
        if diff.area < 1e-11:
            diff = Polygon()
        removed.append([lo, hi, diff])
    return removed


def export_display_data(result, sections, output=OUTPUT):
    """完整轨迹逐步导出；重复投影共享网格，所有割仍保留各自的步号。"""
    meshes, mesh_ids, steps = [], {}, []

    def register(mesh):
        key = json.dumps(mesh, separators=(",", ":"))
        if key not in mesh_ids:
            mesh_ids[key] = len(meshes)
            meshes.append(mesh)
        return mesh_ids[key]

    for k, current in enumerate(sections):
        unchanged = k > 0 and abs(result["volumes"][k] - result["volumes"][k - 1]) < 1e-5
        allowed = steps[-1]["allowed"] if unchanged else register(mesh_from_sections(current))
        removed = register(mesh_from_sections(removed_sections(sections[k - 1], current))) if k and not unchanged else register(dict(v=[], f=[]))
        cut = result["cuts"][k - 1] if k else None
        steps.append(dict(allowed=allowed, removed=removed, volume=round(result["volumes"][k] / 1000, 8),
                          source=[round(cut["source_cost"] / 1000, 6), *cut["source_load"]] if cut else None,
                          design=cut["source_design"] if cut else None))
    payload = dict(meshes=meshes, steps=steps, budget_max=result["budget_max_cny"] / 1000,
                   lambda_max=result["lambda_max"], design_count=len(result["designs"]),
                   cut_count=len(result["cuts"]), final_verified=result["verification"]["complete"])
    path = Path(output) / "display_meshes.json"
    path.write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")
    return path


def pack_display_meshes(path):
    """共用顶点和三角形并压缩嵌入；坐标显示精度 10^-6，无网格采样。"""
    import base64
    import gzip

    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    vertices, vertex_lookup, faces, face_lookup, packed_meshes = [], {}, [], {}, []
    for mesh in payload.pop("meshes"):
        local_ids = []
        for v in mesh["v"]:
            key = tuple(int(round(c * 1e6)) for c in v)
            if key not in vertex_lookup:
                vertex_lookup[key] = len(vertices)
                vertices.append(key)
            local_ids.append(vertex_lookup[key])
        ids = []
        for face in mesh["f"]:
            key = tuple(sorted(local_ids[i] for i in face))
            if len(set(key)) < 3:
                continue
            if key not in face_lookup:
                face_lookup[key] = len(faces)
                faces.append(key)
            ids.append(face_lookup[key])
        packed_meshes.append(sorted(set(ids)))
    payload.update(vertices=vertices, faces=faces, meshes=packed_meshes, coordinate_scale=1e6)
    encoded = json.dumps(payload, separators=(",", ":")).encode("utf-8")
    packed = base64.b64encode(gzip.compress(encoded, compresslevel=9)).decode("ascii")
    destination = Path(path).with_name("display_meshes.base64.txt")
    destination.write_text(packed, encoding="ascii")
    assert len(packed) < 980000, "嵌入数据超过 1 MB，需要进一步共享几何数据"
    return destination


def plot_region_step(k=None, *, process=None, output=OUTPUT, elev=24, azim=-135, show_removed=True):
    """固定三轴的静态快照；同一 process 可用于逐步回放。"""
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    result, _, sections = load_process(output) if process is None else process
    count = len(result["cuts"])
    k = count if k is None else int(k)
    if not 0 <= k <= count:
        raise ValueError(f"k 必须在 0 到 {count} 之间")
    fig = plt.figure(figsize=(10.2, 7.4))
    ax = fig.add_subplot(111, projection="3d")
    layers = [(mesh_from_sections(sections[k]), "#268a86", .70)]
    if k and show_removed:
        layers.append((mesh_from_sections(removed_sections(sections[k - 1], sections[k])), "#d97842", .22))
    for layer_id, (mesh, color, alpha) in enumerate(layers):
        if mesh["v"]:
            coords = np.array(mesh["v"])
            faces = np.array(mesh["f"])
            if layer_id == 0:
                # 打开人为搜索盒子的顶盖和零负荷侧壁，露出最低成本台阶。
                blocks = coords[faces]
                keep = ~((np.max(np.abs(blocks[:, :, 0] - result["budget_max_cny"] / 1000), axis=1) < 1e-8)
                         | (np.max(np.abs(blocks[:, :, 1]), axis=1) < 1e-8)
                         | (np.max(np.abs(blocks[:, :, 2]), axis=1) < 1e-8))
                faces = faces[keep]
            display_coords = coords[:, [1, 2, 0]]
            ax.add_collection3d(Poly3DCollection(display_coords[faces], facecolors=color, edgecolors="none", alpha=alpha))
    ax.set(xlim=(0, result["lambda_max"]), ylim=(0, result["lambda_max"]), zlim=(0, result["budget_max_cny"] / 1000),
           xlabel="区域 A 负荷倍率 λA", ylabel="区域 B 负荷倍率 λB", zlabel="建设成本上限 B（千元）")
    ax.set_box_aspect((1, 1, 1.1))
    ax.view_init(elev=elev, azim=azim)
    ax.set_title(f"两区域负荷规划域 · 第 {k}/{count} 条割" + (" · 完整域已验证" if k == count else " · 当前外近似"))
    fig.text(.10, .035, "台阶下边界：当前最低成本界；其上为保留域。橙色：本步新排除体积。", fontsize=11)
    fig.subplots_adjust(left=0, right=.93, bottom=.1, top=.91)
    return fig


def play_region_cuts(*, pause=.7, output=OUTPUT):
    """Jupyter 中回放从初始盒子到完整三维域的每一条真实割。"""
    import time
    import matplotlib.pyplot as plt
    from IPython.display import clear_output, display
    process = load_process(output)
    for k in range(len(process[0]["cuts"]) + 1):
        clear_output(wait=True)
        fig = plot_region_step(k, process=process)
        display(fig)
        plt.close(fig)
        if pause > 0:
            time.sleep(pause)


if __name__ == "__main__":
    process = run_complete_process()
    print(export_display_data(process[0], process[2]))
