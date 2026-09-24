"""Independent checks and a Chinese report for the two-stage 2-D experiment."""
import argparse
import json
from pathlib import Path

import numpy as np
from threadpoolctl import threadpool_limits
from shapely.geometry import Polygon, box, shape, MultiPoint
from shapely.ops import unary_union

from Network.case33bw import Case33
from model import PlanningEquations, PLANNING_TOL
from tests.case33_two_stage import Interval, fingerprints, interval_gap, save, scheme_facets, serial
from tests.planning_checks import margin
from vertify import ACPowerFlow


def directed_distance(outer, points, inner_polygons=None):
    """Whole-polygon infinity distance to a downward box union, in kW.

    This covers every edge and interior point, not just outer vertices.
    """
    lower, upper = 0., 10000.
    for _ in range(35):
        middle = (lower+upper)/2.
        if inner_polygons is None:
            expanded = unary_union([box(0., 0., p[0]+middle, p[1]+middle) for p in points])
        else:
            offsets = middle*np.array([[-1., -1.], [-1., 1.], [1., -1.], [1., 1.]])
            expanded = unary_union([MultiPoint((vertices[:, None, :]+offsets).reshape(-1, 2)).convex_hull
                                    for vertices in inner_polygons])
        if expanded.covers(outer):
            upper = middle
        else:
            lower = middle
    return upper


def audit(directory):
    result = json.loads((directory/'run/result.json').read_text(encoding='utf-8'))
    scan = json.loads((directory/'reference/scan.json').read_text(encoding='utf-8'))
    grid = np.load(directory/'reference/scan_grid.npz')
    coordinates, states = grid['coordinates'], grid['states']
    net = Case33(upgrade_count=8)
    equations = PlanningEquations(net, 'socp')
    stage1 = Polygon(result['stage1_outer'])
    stage2 = shape(result['stage2_outer'])
    inner_points = np.asarray(result['inner_points'])
    groups = {}
    for c in result['certificates']:
        groups.setdefault(tuple(c['x']), []).append(c['p'][:2])
    inner_polygons = []
    inner_shapes = []
    for group in groups.values():
        points = np.asarray(group)
        cloud = np.vstack((points, points*[1., 0.], points*[0., 1.], [0., 0.]))
        hull = MultiPoint(cloud).convex_hull
        inner_shapes.append(hull)
        inner_polygons.append(np.asarray(hull.exterior.coords if hull.geom_type == 'Polygon' else hull.coords))
    inner = unary_union(inner_shapes)
    certificates = result['certificates']+scan['certificates']
    min_margin, min_support_slack, min_logical_slack = 0., np.inf, np.inf
    ac_pass, trees = 0, {}
    for certificate in certificates:
        x, p, state = (np.asarray(certificate[k]) for k in ('x', 'p', 'state'))
        assert net.cost @ x <= 2.+1e-8
        min_margin = min(min_margin, margin(equations, x, p, state))
        key = tuple(x.astype(int))
        if key not in trees:
            trees[key] = ACPowerFlow(net.tree(x))
        ac_pass += int(trees[key].classify(p)[0] == 1)
        for cut in result['supports']:
            min_support_slack = min(min_support_slack, cut['bound']-np.dot(cut['weights'], p[:2]))
        for cut in result['logical_cuts']:
            if p[0] >= cut['threshold']+1e-6:
                min_logical_slack = min(min_logical_slack, cut['bound']-p[1])
    grid_tops = np.array([(coordinates[i], coordinates[np.where(row == 1)[0][-1]])
                         for i, row in enumerate(states) if np.any(row == 1)])
    actual = unary_union([box(0., 0., *p) for p in grid_tops if np.min(p) > 0.])
    # Independently replay every feasible grid column with its saved AC witness.
    grid_ac_count, grid_ac_fail = 0, 0
    for column in scan['columns']:
        if column.get('witness_x') is None:
            continue
        i = column['index']
        second = coordinates[states[i] == 1]
        power = np.zeros((len(second), 3))
        power[:, 0], power[:, 1] = coordinates[i], second
        key = tuple(column['witness_x'])
        if key not in trees:
            trees[key] = ACPowerFlow(net.tree(np.asarray(key)))
        checked = trees[key].classify(power)
        grid_ac_count += len(checked)
        grid_ac_fail += int(np.count_nonzero(checked != 1))
    inner_facets = scheme_facets(result['certificates'])
    geometric_bound = max(interval_gap(Interval(**node), result['supports'], inner_points,
                          inner_facets, np.asarray(result['stage1_outer']))
                          for node in result['intervals'])
    square_area = result['metadata']['square_kw']**2
    hull_inner = MultiPoint(np.vstack((inner_points, inner_points*[1., 0.],
                                       inner_points*[0., 1.], [0., 0.]))).convex_hull
    metrics = dict(status=result['status'], certificates=len(certificates), ac_pass=ac_pass,
        min_raw_margin=min_margin, min_support_slack_kw=min_support_slack,
        min_logical_slack_kw=min_logical_slack, grid_counts=scan['counts'],
        grid_ac_replayed=grid_ac_count, grid_ac_fail=grid_ac_fail,
        square_area_kw2=square_area, stage1_area_kw2=stage1.area,
        stage2_area_kw2=stage2.area, inner_area_kw2=inner.area,
        reference_area_kw2=actual.area, stage1_removed_area_kw2=square_area-stage1.area,
        stage2_removed_area_kw2=stage1.area-stage2.area,
        certified_nonconvex_area_kw2=hull_inner.difference(stage2).area,
        max_gap_kw=geometric_bound, exact_outer_to_inner_gap_kw=directed_distance(stage2, inner_points, inner_polygons),
        outer_to_reference_gap_kw=directed_distance(stage2, grid_tops),
        reference_outside_final_area_kw2=actual.difference(stage2.buffer(.001)).area,
        reference_outside_stage1_area_kw2=actual.difference(stage1.buffer(.001)).area,
        stage2_extra_vs_reference_percent=100*stage2.difference(actual).area/actual.area,
        run_seconds=result['seconds'], reference_seconds=scan['seconds'],
        original_model_check_seconds=scan.get('original_model_check_seconds', 0.),
        support_count=len(result['supports']), branch_queries=len(result['logical_cuts']),
        reference_queries=len(scan['events']),
        reference_cached_columns=sum(c.get('cached_AC_witness', False) for c in scan['columns']),
        run_schemes=len({tuple(c['x']) for c in result['certificates']}),
        reference_schemes=len({tuple(c['x']) for c in scan['certificates']}),
        source_unchanged=fingerprints() == result['metadata']['source_hashes'] == scan['source_hashes'])
    checks = dict(complete=result['status'] == 'certified', physical=min_margin >= -PLANNING_TOL,
        ac=ac_pass == len(certificates), support_valid=min_support_slack >= -.001,
        conditional_cuts_valid=min_logical_slack >= -.001, grid_complete=scan['counts']['0'] == 0,
        grid_ac=grid_ac_fail == 0 and grid_ac_count == scan['counts']['1'],
        geometry=geometric_bound <= result['metadata']['epsilon_kw']+1e-6,
        full_polygon_geometry=metrics['exact_outer_to_inner_gap_kw'] <= result['metadata']['epsilon_kw']+1e-6,
        largest_gap_priority=all(abs(current['selected_gap_kw']-previous['max_gap_kw']) < 1e-7
            for previous, current in zip(result['history'], result['history'][1:])),
        reference_contained=metrics['reference_outside_final_area_kw2'] < 1e-5,
        production_unchanged=metrics['source_unchanged'])
    save(directory/'audit.json', dict(checks=checks, metrics=metrics))
    rows = '\n'.join(f'| {name} | {"通过" if passed else "未通过"} |' for name, passed in checks.items())
    report = f'''# Case33bw 二维两阶段测试结果

## 结果

- 截面：节点 18、25；节点 33 固定 0 kW；37 条可重构线路、8 个升级候选、建设预算 2。
- 初始正方形：`[0,6855] × [0,6855]` kW。6855 来自根端有功上限扣除固定负荷，不是变压器铭牌容量。
- 算法状态：**{result['status']}**。最终认证几何间隙 **{geometric_bound:.3f} kW**，目标为 {result['metadata']['epsilon_kw']:g} kW。
- 独立参考：**{scan['spacing_kw']:g} kW 网格**，共 {states.size:,} 个点；可行 {scan['counts']['1']:,}、不可行 {scan['counts']['-1']:,}、未定 {scan['counts']['0']:,}。
- 主线文件哈希核对：{'未改变' if metrics['source_unchanged'] else '有差异'}。

| 指标 | 数值 |
|---|---:|
| 第一阶段支持割（含 3 条初始有效界） | {len(result['supports'])} |
| 第二阶段区间定界查询 | {len(result['logical_cuts'])} |
| 第一阶段切除面积 | {metrics['stage1_removed_area_kw2']:,.1f} kW² |
| 第二阶段再切除面积 | {metrics['stage2_removed_area_kw2']:,.1f} kW² |
| 第二阶段切除占第一阶段外域的比例 | {100*metrics['stage2_removed_area_kw2']/stage1.area:.3f}% |
| 已严格证出的非凸缺口面积下界 | {metrics['certified_nonconvex_area_kw2']:,.1f} kW² |
| 最终外域面积 | {stage2.area:,.1f} kW² |
| 5 kW 扫描参考内域面积 | {actual.area:,.1f} kW² |
| 最终外域比网格参考多出的面积 | {metrics['stage2_extra_vs_reference_percent']:.3f}% |
| 参考可行域落在最终外域之外 | {metrics['reference_outside_final_area_kw2']:.6g} kW² |
| 完整外多边形到算法认证内域的距离 | {metrics['exact_outer_to_inner_gap_kw']:.3f} kW |
| 完整外多边形到网格参考内域的距离 | {metrics['outer_to_reference_gap_kw']:.3f} kW |
| 两阶段总耗时 | {result['seconds']:.1f} s |
| 独立扫描总耗时 | {scan['seconds']:.1f} s |
| 另行原模型未决点核查求解耗时 | {scan.get('original_model_check_seconds', 0.):.1f} s |
| 独立扫描的全网架查询 | {len(scan['events'])} |
| 扫描中复用 AC 证书即可判定的列 | {metrics['reference_cached_columns']} |
| 两阶段保留的不同可行网架 | {metrics['run_schemes']} |

## 图与切割过程

![两阶段总览](figures/two_stage_overview.png)

![边界细节](figures/two_stage_boundary.png)

![分阶段切割](figures/cutting_sequence.png)

[逐步回放](figures/step_by_step.html) · [认证间隙](figures/certified_gap.png) · [支持割系数](figures/support_lines.csv) · [区间逻辑割](figures/branch_cuts.csv) · [参考边界点](figures/actual_boundary.csv)

颜色严格区分阶段：浅红实色是第一阶段切除，红色斜线是第二阶段进一步切除；蓝、绿、红虚线分别是第一阶段外包络、第二阶段最终外近似、独立扫描的 Actual region。斜线区域包含非凸缺口及第一阶段剩余的外近似余量，不能将其全部面积都解释为真实非凸缺口。

绿色是有误差保证的外近似，红色是标明精度的网格参考内域，均不是解析真边界。网格单调性为连续域参考提供 {scan['spacing_kw']:g} kW 的坐标分辨率；算法的 {result['metadata']['epsilon_kw']:g} kW 认证阈值另行计算。这里的距离是外域任一点到内域最近点的有向无穷范数距离：两个负荷坐标允许同时改变。它不是固定横坐标时纵向容量误差不超过 {result['metadata']['epsilon_kw']:g} kW。

算法内域按完全相同的网架分组，利用固定整数方案下 SOCP 可行集的凸性取向下凸包，然后对不同方案取并集；没有跨网架凸化。连续几何认证首先针对本次 SOCP 规划模型；所有原始可行点另经 AC 核验，完整绿色外域到独立 AC 网格内域的距离也单独列出。

非凸性证据采用 `conv(已认证可行点的下闭包) \\ 最终认证外域`。这个差集既在真实规划域的凸包内，又已经被全局界证明在规划域之外，因此其正面积可以证明存在非凸缺口；不依赖将网格折线当成解析边界。

## 方法与调用关系

`two_stage → SliceOracle.solve → 原 PlanningModel 的 Gurobi 模型`。支持目标直接读取全局 ObjBound；之后每个节点代表一个负荷区间，用 `max p25, p18 >= t` 对全部允许建设方案定界，再选最大认证间隙的节点分支。优先级是整个节点的几何上界，不是某一固定网架的 SP 违反量。

第一阶段在预算约束和所有允许的整数网架上计算支持函数的全局上界，再加入 `omega @ p <= bound`。第二阶段对同一完整规划模型计算 `U(t) >= max p25 : p18 >= t`，所以 `p18 >= t 且 p25 > U(t)` 中不存在任何满足预算的网架。这是区间条件割能跨网架排除区域的原因；固定一个网架得到的 SP 割没有被直接用于裁掉整个规划域。

此处实现的是**目标空间分支定界**，每个标量定界问题内部使用整数分支定界。它没有实现跨方向共享同一棵建设决策搜索树的完整 MOMIX；模型对象和可行起点复用不等于搜索树复用。运行未显式枚举全部建设方案。

测试模型加入根向父弧、单位多类连通流、纯负荷方向功率及电压上界等有效加强，不改变允许的整数树。前提是当前算例的非负负荷、正线路阻抗、全部节点接入。固定已知网架的连续 SOCP 只用于提供完整可行初值；全局定界始终允许改变全部建设决策。93 份加强前的原始证书对新增功率方向和电压不等式的最大违反量为 0；另有原模型与加强模型的固定树数值回归。

参考扫描另建模型、使用种子 17，不导入两阶段的支持割、网架或边界。逐列全局上界证明域外网格点，每个域内网格点用独立 AC 潮流检查；单调性及已找到网架的 AC 证据减少重复整数求解。不是对 N² 个点各运行一次完整 MIP，也不是仅扫描两阶段已知网架。

## 审核

| 检查 | 结果 |
|---|---|
{rows}

全部 {len(certificates)} 份原始求解证书记录（含重复网架和热启动记录）进行了原方程代入与 AC 复核；另按各列存储的网架重算 {grid_ac_count:,} 个参考可行网格点。外多边形距离审核覆盖完整边和内部，不使用“顶点都可行所以非凸并集覆盖”的错误判据。

逐次核对本轮所选节点的认证间隙等于上轮结束时全局最大节点间隙，验证最大认证间隙优先的执行顺序。

GitNexus 绑定 PlanRegion（当前工作区；索引提交 `7c367e6b8609d04c8fbc2f2d6fc5d630d58db1e8`）。新增实验文件的 impact 为 UNKNOWN，已用文本检查确认引用限于新测试和文档；符号文档的 impact 为 MEDIUM，直接引用为 6 份项目文档，未触及生产调用流程。未提交变更。

## 复现

```text
python -m tests.case33_two_stage --mode run --time-limit 120 --epsilon-kw 20 --output results/case33_two_stage/20260924/run
python -m tests.case33_two_stage --mode scan --threads 4 --time-limit 120 --spacing-kw 5 --output results/case33_two_stage/20260924/reference
python -m tests.check_case33_reference results/case33_two_stage/20260924/reference
python -m tests.audit_case33_two_stage results/case33_two_stage/20260924
python -m tests.plot_case33_two_stage results/case33_two_stage/20260924
python -m unittest tests.test_notation tests.test_two_stage tests.test_boundary_search -v
```

以上是顺序复现命令，应选择新输出目录，已有结果不会被静默覆盖。本次运行环境为 Intel Core Ultra 9 285K（24 核）、Gurobi 13.0.2、Python 3.13。两阶段构域使用 20 线程；与之重叠的独立扫描使用 4 线程。后续独立扫描可用 `python -m tests.parallel_case33_scan <reference目录>` 将剩余列分为 6 个 4 线程进程，前提是原扫描已停止；各进程只继承扫描自身的已完成证据。原模型未决点核查使用 20 线程，另列求解时间。BLAS 均限制为 1 线程。扫描表中的耗时为已记录串行耗时加并行阶段墙钟耗时，不是各进程 CPU 时间之和。

最终数据来自定型后从初始正方形重新运行的 `run`；`pilot` 是开发期间试验，不纳入最终构域耗时。完整元数据、上下界、求解状态、逐步多边形和证书保存于 `run/result.json` 与 `reference/scan.json`，并保存对应源码快照。正式代码哈希保存在两份结果内。当前测试是一次算例验证，不据此声称普遍加速或最坏情况复杂度降低。
'''
    (directory/'REPORT.md').write_text(report, encoding='utf-8')
    print(json.dumps(serial(dict(checks=checks, metrics=metrics)), ensure_ascii=False, indent=2))
    return all(checks.values())


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', type=Path)
    with threadpool_limits(limits=1):
        raise SystemExit(0 if audit(parser.parse_args().directory) else 1)
