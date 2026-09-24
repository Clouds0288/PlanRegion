"""Publish measured three-budget results only after all certificate audits pass."""
import argparse
import hashlib
import json
from importlib.metadata import version
from pathlib import Path
import platform
import shutil
import sys

import gurobipy as gp
import numpy as np
import scipy

from Network.case33bw import Case33
from tests.case33_two_stage import save, fingerprints, SOURCE_FILES


def duration(seconds):
    seconds = int(round(seconds))
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return f'{hours} h {minutes:02d} min {seconds:02d} s' if hours else f'{minutes} min {seconds:02d} s'


def report(directory):
    network = Case33(upgrade_count=8)
    projects = [dict(corridor=c.id,cost=c.types[1].investment_cost,
                     r_existing=c.types[0].r,reactance_existing=c.types[0].reactance,
                     r_parallel=c.types[1].r,reactance_parallel=c.types[1].reactance,
                     capacity_is_unbounded=bool(np.isinf(c.types[0].capacity)))
                for c in network.corridors if len(c.types)>1]
    results = [json.loads((directory/f'budget_{b}/result.json').read_text(encoding='utf-8')) for b in (0,2,4)]
    audits = [json.loads((directory/f'budget_{b}/audit.json').read_text(encoding='utf-8')) for b in (0,2,4)]
    assert all(r['status'] == 'certified' for r in results)
    assert all(a['passed'] and len(a['slices']) == 3 for a in audits)
    assert all(r['source_hashes'] == fingerprints() for r in results)
    lines = ['# Case33bw：三维可行规划域与三档建设预算', '',
        '坐标为 $(p_{18},p_{25},p_{33})$，单位 kW。建设预算为 0、2、4 个相对投资单位，三个场景使用相同的 8 个升级候选。预算 0 仍允许免费改变开关连接；它不是固定初始网架。', '',
        '并联升级将对应线路的电阻和电抗减半。2–3 线路费用为 2，其余七条候选各为 1，预算不能直接解读为升级条数。候选为 2–3、5–6、12–13、13–14、23–24、24–25、27–28、28–29。', '',
        '本次沿用原始 Case33bw 参数：线路热额定容量未给定，代码将其视为无限额，因此这些结果主要反映电压、源端功率、径向拓扑和预算约束的共同作用。其他节点负荷固定；完整运行方程见[原 SOCP 模型文档](../../../docs/socp_model.md)。', '',
        '**三个场景均通过完整三维 SOCP 构域审计；独立截面采用 10 kW 网格。** 绿色域是认证外近似；红色 Actual region 是独立扫描的网格参考域。', '',
        '## 结果', '',
        '| 预算 | 支持割 / 第二阶段割 | 最大认证间隙 (kW) | 可行体积认证区间 (10⁹ kW³) | 证书网架数 | 累计构域时间 |',
        '|---:|---:|---:|---:|---:|---:|']
    for r,a in zip(results,audits):
        lines.append(f'| {r["budget"]:g} | {len(r["supports"])} / {len(r["branch_cuts"])} | {a["max_gap_kw"]:.4f} | '
                     f'[{a["inner_volume_lower"]/1e9:.6f}, {a["stage2_volume"]/1e9:.6f}] | {a["scheme_count"]} | {duration(r["seconds"])} |')
    lines += ['', '体积区间的下端来自认证内域的保守下界，上端是最终外域体积。证书网架数只统计保存了可行证书的不同方案，不代表求解器内部访问或发现的全部方案；每次全局查询仍允许全部满足预算的建设方案。预算集合在数学上嵌套，有限误差下的独立外近似边界不要求逐点完全嵌套。', '',
        '| 预算 | 第一阶段外域 (10⁹ kW³) | 阶段一凸包误差上界 (kW) | 第二阶段切除 (10⁹ kW³) | 已证凸包内部非凸缺口下界 (10⁹ kW³) |',
        '|---:|---:|---:|---:|---:|']
    for a in audits:
        lines.append(f'| {a["budget"]:g} | {a["stage1_volume"]/1e9:.6f} | {a["convex_support_gap_upper_kw"]:.3f} | {a["removed_stage2_volume"]/1e9:.6f} | {a["nonconvex_volume_lower"]/1e9:.6f} |')
    lines += ['', '红色斜线表示第二阶段实际切除的全部区域，其中也可能包含第一阶段剩余的凸外包络误差。表中最后一列只统计已证明位于真实 SOCP 可行域凸包内、而位于可行域外的部分，避免把所有第二阶段切除都冒称为非凸缺口。', '',
        '## 独立截面', '',
        '| 预算 | 固定 p33 (kW) | 可行网格点 | 未确定点 | 网格内域面积 (10⁶ kW²) | 绿色截面到红色网格域距离 (kW) |',
        '|---:|---:|---:|---:|---:|---:|']
    scanned, feasible = 0, 0
    for a in audits:
        for s in a['slices']:
            gap = '空域' if s['outer_to_reference_gap_kw'] is None else f'{s["outer_to_reference_gap_kw"]:.3f}'
            scanned += sum(s['grid_counts'].values())
            feasible += s['grid_counts']['1']
            lines.append(f'| {a["budget"]:g} | {s["fixed_value"]:g} | {s["grid_counts"]["1"]:,} | '
                         f'{s["grid_counts"]["0"]} | {s["actual_area"]/1e6:.6f} | {gap} |')
    lines += ['', f'共分类 {scanned:,} 个网格点，{feasible:,} 个可行点全部通过独立 AC 潮流回放，未确定点为 0。扫描利用单调性和全局定界批量分类网格点，不对每个点重复启动整数求解。', '',
        '截面扫描从原物理模型重新开始，不使用三维构域的方案、切面或可行证书。中预算 p33=0 使用已有同源独立 5 kW 扫描的 10 kW 子网格，其余截面本次独立计算。', '',
        '三维 20 kW 证书允许三个坐标同时改变；固定 p33 后的二维误差另行计算，不能直接继承 20 kW。红色阶梯线表示网格内近似，不是解析精确边界。', '',
        '## 两阶段数学模型', '',
        r'记预算 $\mathcal B$ 下的完整 MISOCP 可行集为 $\mathcal F_{\mathcal B}$，投影规划域为 $\mathcal D_{\mathcal B}=\{p:\exists x,y,(x,p,y)\in\mathcal F_{\mathcal B}\}$。沿用符号契约：x 为选型、p 为负荷、y 对应代码运行向量 state。初始立方体为 $[0,6855]^3$。6855 kW 是本实验源端有功限额扣除固定负荷后的必要上界，不是从 MATPOWER 自动推定的变压器铭牌容量。', '',
        r'第一阶段：对非负方向 $\omega$ 求 $h_{\mathcal B}(\omega)=\max_{\mathcal F_{\mathcal B}}\omega^Tp$，使用全局上界 $\overline h_{\mathcal B}$ 添加 $\omega^Tp\le\overline h_{\mathcal B}$。所得蓝色多面体包含所有允许建设方案的可行点。', '',
        r'第二阶段：节点是目标空间的凸多面体 $O_s$。同一完整建设方案的已认证点可取向下凸包 $I_x$；不同方案只取并集。节点认证上界为 $\bar\Delta_s=\min_x\max_{v\in\operatorname{vert}(O_s)}\operatorname{dist}_{\infty}(v,I_x)$。始终处理最大上界节点；顶点分别被不同方案覆盖而节点整体未被认证时，沿已知方案的 20 kW 距离邻域边界分支，必要时按最长坐标边二分。分支保留父节点全部体积，不冒称为不可行性割。', '',
        r'全局定界采用 $\min\rho$，满足 $(x,p,y)\in\mathcal F_{\mathcal B}$、$Wp+\rho\mathbf1\ge Wq$、$\rho\ge0$，其中 W 各行非负且和为 1。全局下界 $\underline\rho>0$ 产生析取割 $\bigvee_k\{W_kp\le W_kq-\underline\rho\}$。单位矩阵 W 对应正交象限割；斜向 W 对应倾斜切面。切割依据完整规划问题的全局界，不把固定网架的局部割直接用于其他网架。', '',
        r'最终 $\mathcal I\subseteq\mathcal D_{\mathcal B}\subseteq\mathcal O_2$，且 $\sup_{p\in\mathcal O_2}\operatorname{dist}_{\infty}(p,\mathcal I)\le20$ kW。认证在原 SOCP 模型和求解器数值容差下成立；端点通过 AC 不等于证明全部凸插值都满足原始非凸 AC 方程，独立截面用于检验这部分差异。', '',
        '## 图与切面数据', '',
        '- [可旋转三维、全部切割步骤与截面开关](figures/three_dimensional.html)',
        '- [三预算三维边界](figures/planning_regions_3d.png)',
        '- [从立方体到两阶段结果](figures/two_stage_3d_cuts.png)',
        '- [九张独立截面](figures/independent_sections.png)',
        '- [认证收敛与预算体积区间](figures/certification_and_budget.png)',
        '- 全部切面系数：[预算 0](figures/budget_0_cuts.csv)、[预算 2](figures/budget_2_cuts.csv)、[预算 4](figures/budget_4_cuts.csv)。', '',
        '同名 PDF/SVG 提供矢量版本。budget_*_cuts.csv 每面一行；同一第二阶段 index 的各面以 OR 连接。交互回放保存每次有效切割，压缩仅去除重复图形数据。浏览器安全策略阻止工具打开本地 HTML，因此实际旋转和按钮交互未能在浏览器中验收；已检查 JavaScript 语法、几何回放数据和静态图。', '',
        '## 计时与复现', '',
        '| 预算 | 阶段一 | 阶段二 | 本次新增截面累计 | 复用截面的原始耗时 |',
        '|---:|---:|---:|---:|---:|']
    for r,a in zip(results,audits):
        fresh, historical = 0.,0.
        for s in a['slices']:
            scan = json.loads((directory/f'budget_{r["budget"]:g}'/f'slice_{s["fixed_value"]:g}'/'scan.json').read_text(encoding='utf-8'))
            if scan.get('reused'):
                historical += scan['source_original_seconds']+scan.get('source_original_check_seconds',0.)
            else:
                fresh += scan['seconds']
        lines.append(f'| {r["budget"]:g} | {duration(r["stage1_seconds"])} | {duration(r["seconds"]-r["stage1_seconds"])} | {duration(fresh)} | {duration(historical) if historical else "—"} |')
    lines += ['', '三个预算及截面任务有并行重叠，上表是各任务检查点记录的累计时间，不能相加当作总墙钟时间。本次包含增量算法试验、线程资源调整及断点恢复，已中断但未写入检查点的求解和独立探针开销不计入该表；它不是一次干净启动的速度基准，也不宜直接比较不同预算的计算效率。', '',
        '生产 main.py、model.py、region.py、vertify.py 和网络数据散列均保持不变。最终实验源码、环境版本、生产散列与审计结果保存在 source_bundle/final 和 summary.json。', '',
        '```powershell',
        'python -m tests.case33_3d --mode run --budget 2 --threads 8 --time-limit 90 --oblique --output results/new_3d/budget_2',
        'python -m tests.case33_3d --mode scan --budget 2 --fixed-value 500 --spacing-kw 10 --threads 4 --time-limit 120 --output results/new_3d/budget_2/slice_500',
        'python -m tests.audit_case33_3d results/new_3d/budget_2',
        'python -m unittest tests.test_notation tests.test_case33_3d tests.test_two_stage tests.test_boundary_search -v',
        '```', '',
        '分别运行预算 0、2、4 和固定 p33=0、500、1000 后，再执行 plot_case33_3d 和 report_case33_3d。独立重新运行可以全部从零扫描，无需复用历史截面。']
    (directory/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    environment = dict(python=sys.version, platform=platform.platform(), gurobi=gp.gurobi.version(),
                       numpy=np.__version__, scipy=scipy.__version__,
                       figure_packages={name:version(name) for name in ('matplotlib','plotly','shapely','threadpoolctl')})
    sources = ['tests/case33_3d.py','tests/audit_case33_3d.py','tests/plot_case33_3d.py','tests/recheck_case33_3d_columns.py',
               'tests/report_case33_3d.py','tests/run_case33_3d_validation.py','tests/repair_case33_3d_slice.py','tests/test_case33_3d.py',
               'tests/case33_two_stage.py','tests/plot_case33_two_stage.py','tests/planning_checks.py',
               'tests/audit_case33_two_stage.py','docs/notation.md','docs/case33_3d_test_plan.md','docs/socp_model.md',*SOURCE_FILES]
    frozen = directory/'source_bundle/final'
    manifest = {}
    for source in sources:
        target = frozen/source
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(source,target)
        manifest[source] = hashlib.sha256(target.read_bytes()).hexdigest()
    save(frozen/'manifest.json',dict(environment=environment,sha256=manifest))
    save(directory/'summary.json',dict(environment=environment,audits=audits,
                                       case_config=dict(load_nodes=network.load_nodes,base_kva=network.base,projects=projects),
                                       production_hashes=fingerprints(),grid_points=scanned,ac_feasible_points=feasible))
    print(directory/'report.md')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('directory',type=Path)
    report(parser.parse_args().directory)
