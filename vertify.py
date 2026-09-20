"""独立潮流交叉核验、等体积网格 FR/MR 和单文件结果存取。

这里不构建切割模型。区域成员标签必须使用同一评价箱、同一网格中心，
这样点数之比才等于网格体积之比；有限网格上的零误差不代表连续域精确相同。
"""
from dataclasses import dataclass
import json
from pathlib import Path
import numpy as np
from model import ACPowerFlow, PlanningEquations, dispatch_support
from region import GEOMETRY_TOL, contains, halfspaces, polytope_size

METHODS = ('socp', 'hybrid', 'ac', 'linear')  # 结果数组的方法轴顺序，AC 固定在索引 2。
METHOD_NAMES = ('纯 SOCP 切割', '线性 + SOCP 精修', 'AC 数值参考', '纯线性切割')


def validate_power_flow(network):
    """用六个工况，将独立支路递推与节点导纳矩阵 Newton–Raphson 潮流核对。"""
    import warnings
    import pandapower as pp
    from pandapower.converter.pypower.from_ppc import from_ppc

    reference = ACPowerFlow(network)  # 被核验的是独立 AC 递推，不是 SOCP 方程。
    maximum_difference = 0.  # 记录所有工况、所有节点的最大电压幅值差。
    for scale in (1., 0., .5, 1.2, 1.5, 2.):
        power = scale*network.original_p[network.selected]  # 只缩放三个独立节点，其他背景负荷固定。
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=FutureWarning,
                                    module='pandapower.converter.pypower.from_ppc')
            net = from_ppc(network.ppc(power), f_hz=50)  # 同一物理网架转换成节点导纳矩阵模型。
        pp.runpp(net, algorithm='nr', tolerance_mva=1e-10, numba=False)  # 另一实现的完整 AC 潮流作为交叉核验。
        # 此处只比较潮流方程；即使电压越限也继续收敛，不能复用 classify 的提前不可行判定。
        ell = np.zeros((1, network.n))  # 支路递推从零电流平方开始。
        for _ in range(160):
            P, Q, v, u = reference.state(power, ell)  # 当前电流下的送端功率及两端电压平方。
            new = (P*P+Q*Q)/u  # 强制完整 AC 电流等式。
            if np.max(np.abs(new-ell)) < 1e-13:  # 电流平方增量达到比运行限值判定更严格的精度。
                break
            ell = new  # 更新电流，进行下一次完整功率平衡与压降计算。
        difference = float(np.max(np.abs(np.sqrt(v[0])-net.res_bus.loc[list(network.nodes), 'vm_pu'])))  # v 是平方，NR 输出是幅值。
        assert difference < 1e-8, 'Branch AC and nodal AC disagree'  # 不让方程实现不一致的参考模型进入正式评价。
        maximum_difference = max(maximum_difference, difference)
        if scale == 1:  # 保存原始基础工况的物理校验值。
            baseline = dict(minimum_voltage_pu=float(np.sqrt(v.min())),  # 全网最低电压幅值。
                            minimum_voltage_bus=int(network.nodes[np.argmin(v)]),  # 对应真实节点号。
                            active_loss_kw=float((ell*network.r).sum()*network.base))  # Σr*ell，从标幺还原为 kW。
    return dict(baseline=baseline, nodal_cross_checks=6,
                maximum_voltage_difference_pu=maximum_difference, pandapower=pp.__version__)


def region_membership(points, polytopes):
    """相同几何口径的并集判定；先用包围盒和总负荷上界排除无关网格点。"""
    mask = np.zeros(len(points), dtype=bool)  # 初始所有点均未落入任何方案域。
    totals = points.sum(axis=1)  # 各点的独立节点总有功负荷。
    for poly in sorted(polytopes, key=polytope_size, reverse=True):  # 先判断大域，减少后续重复计算。
        candidate = ~mask & (totals <= poly.sum(axis=1).max()+GEOMETRY_TOL)  # 已在并集中的点无需再检查。
        candidate &= np.all(points >= poly.min(axis=0)-GEOMETRY_TOL, axis=1)  # 先用坐标下界排除域外点。
        candidate &= np.all(points <= poly.max(axis=0)+GEOMETRY_TOL, axis=1)  # 再用坐标上界排除域外点。
        indices = np.flatnonzero(candidate)  # 只有可能落入此多面体的点才进行逐面判断。
        equations = halfspaces(poly)  # 统一转换为 a·p+b≤0 的单位法向量约束。
        for start in range(0, len(indices), 32768):  # 分批控制“点×面”矩阵大小，不改变判定口径。
            take = indices[start:start+32768]
            mask[take] = contains(points[take], equations)  # 复用几何模块的同一包含判据和容差。
    return mask  # 每个 True 只要求至少存在一个可行建设方案。


def disagreement(approximation, ac):
    """同一等体积网格上的 FR/MR；空域的条件比例无定义。"""
    extra = int(np.count_nonzero(approximation & ~ac))  # 黄色：计算域有、AC 参考域没有。
    missed = int(np.count_nonzero(ac & ~approximation))  # 红色：AC 参考域有、计算域遗漏。
    computed, reference = int(approximation.sum()), int(ac.sum())  # 两个域各自的网格单元数量。
    union = computed+missed  # 计算域∪AC 域；遗漏部分与计算域不相交。
    return dict(fr_percent=100*extra/computed if computed else None,  # FR 分母是计算域。
                mr_percent=100*missed/reference if reference else None,  # MR 分母是 AC 参考域。
                region_error_percent=100*(extra+missed)/union if union else 0.,  # 原区域不一致率：对称差/并集。
                computed_cells=computed, ac_cells=reference, extra_cells=extra, missed_cells=missed)


def verify_joint_cuts(network, method, cuts):
    """逐条割在全部 16 个固定方案域上最小化，而非只检查有限样本点。"""
    equations = PlanningEquations(network,method)
    minimum = np.empty((len(cuts),len(network.designs)))
    for j,design in enumerate(network.designs):
        x = equations.selection(design.x)
        for i,cut in enumerate(cuts):
            # min(a+bᵀp+dᵀx)=a+dᵀx-max(-bᵀp)；固定域由旧消元方程独立建立。
            support = dispatch_support(design,method,-np.asarray(cut[1:4]))
            minimum[i,j] = cut[0]+np.asarray(cut[4:])@x-support['bound']
    return minimum  # 使用最大化问题的对偶上界，得到割最小余量的数值下界。


def compact_review_summary(record):
    """从唯一审核记录计算核对表，不另存重复的汇总数据。"""
    costs = np.array(record['design_costs'])
    budgets = [np.inf if b is None else b for b in record['budgets']]
    rows = []
    for method,data in record['models'].items():
        reference = np.array(data['enumerated_cost'])
        actual = np.array(data['compact_cost'])
        reference_radii = np.array([np.max(np.array(data['enumerated_radii'])[costs<=b],axis=0) for b in budgets])
        membership_difference = sum(np.count_nonzero((np.isfinite(reference)&(reference<=b))
                                    !=(np.isfinite(actual)&(actual<=b))) for b in budgets)
        joint_cost_errors, joint_radius_errors = 0, []
        for query in data['joint_queries']:
            if query['kind']=='point':
                joint_cost_errors += int(query['objective'] != reference[query['index']])
            else:
                joint_radius_errors.append(abs(query['objective']-reference_radii[query['budget'],query['index']]))
        rows.append(dict(method=method,point_count=len(reference),
            minimum_cost_mismatches=int(np.count_nonzero(reference!=actual)),
            budget_membership_mismatches=int(membership_difference),
            radial_queries=int(reference_radii.size),
            max_radius_difference_kw=float(np.max(np.abs(np.array(data['compact_radii'])-reference_radii))),
            joint_queries=len(data['joint_queries']),joint_cost_mismatches=joint_cost_errors,
            joint_max_radius_difference_kw=float(max(joint_radius_errors)),
            joint_max_gap_kw=float(max(q['bound']-q['objective'] for q in data['joint_queries'] if q['kind']=='ray')),
            generated_cuts=len(data['cuts']),audited_cuts=len(data['cut_minimum']),
            minimum_cut_margin=float(np.min(data['cut_minimum'])),**data['seconds']))
    return rows


@dataclass
class BenchmarkResult:
    """一个实验一个文件：三法掩码、AC 最小投资、计时和去重后的方案记录。"""
    masks: np.ndarray  # 内存形状：(4 种方法, 预算数, N, N, N)，存网格中心的区域归属。
    metadata: dict  # 配置、方案表、实测耗时和源码指纹。
    regions: dict  # 每种方法、每档预算的逐方案几何及真实 SP 历史。
    ac_cost: np.ndarray  # 每个网格点的最小 AC 可行投资；inf 表示无已允许方案可行。

    @property
    def budgets(self):
        return tuple(np.inf if b is None else b for b in self.metadata['budgets'])  # JSON 中用 null 表示无限预算。

    @property
    def bounds(self):
        return np.array(self.metadata['bounds'])  # 三个独立负荷轴的公共评价上界，单位 kW。

    @property
    def load_nodes(self):
        return tuple(self.metadata['load_nodes'])  # 数组的三个坐标与真实节点号对应。

    @property
    def spacing(self):
        return self.bounds/self.metadata['divisions']  # 每个体素沿三条轴的步长。

    @property
    def labels(self):
        return self.masks.astype(np.uint8)+2*self.masks[2].astype(np.uint8)  # 0 域外，1 多余，2 遗漏，3 重合。

    @property
    def summary(self):
        rows = []
        for k, method in enumerate(METHODS):
            for j, budget in enumerate(self.metadata['budgets']):
                timing = self.metadata['timings'][method][j]  # 对应方法和预算的本次实测耗时。
                rows.append(dict(method=method, budget=budget,
                    **disagreement(self.masks[k, j], self.masks[2, j]), **timing,  # 每次从原始标签计算 FR/MR，不重复保存指标。
                    total_seconds=timing['region_seconds']+timing['evaluation_seconds']))  # 构域加网格归属时间，绘图不计。
        return rows

    def save(self, folder):
        folder = Path(folder)
        folder.mkdir(parents=True, exist_ok=True)
        # 相同方案记录在多档预算中只存一次；AC 不再逐预算重复保存布尔掩码。
        unique = {}  # 规范化记录 → 唯一编号。
        encode = lambda r: json.dumps(r, sort_keys=True, default=lambda x: x.tolist())  # 数组转列表，统一字典键顺序以便去重。
        groups = {m: [[unique.setdefault(encode(r), len(unique)) for r in batch]
                      for batch in batches] for m, batches in self.regions.items()}  # 各预算只保存对唯一记录的编号引用。
        record = dict(metadata=self.metadata, geometry=[json.loads(r) for r in unique], groups=groups)  # 唯一几何表及其引用关系。
        np.savez_compressed(folder/'result.npz', masks=np.packbits(self.masks[[0, 1, 3]]),  # 仅保存三种构域方法，布尔值按位压缩。
                            ac_cost=self.ac_cost, record=json.dumps(record, ensure_ascii=False))  # AC 各预算由最小投资推导。

    @classmethod
    def load(cls, folder):
        with np.load(Path(folder)/'result.npz', allow_pickle=False) as data:
            record = json.loads(str(data['record']))  # 读取当前格式的唯一记录，不猜测或兼容历史格式。
            meta, costs = record['metadata'], data['ac_cost']
            shape = (4, len(meta['budgets']))+(meta['divisions'],)*3  # 恢复方法、预算和三个空间轴。
            masks = np.zeros(shape, dtype=bool)
            masks[[0, 1, 3]] = np.unpackbits(data['masks'], count=3*int(np.prod(shape[1:]))).reshape((3,)+shape[1:])  # count 排除末字节补齐位。
            for j, budget in enumerate(meta['budgets']):
                masks[2, j] = np.isfinite(costs) & (costs <= (np.inf if budget is None else budget))  # 有限最小投资≤预算才属于 AC 域。
        regions = {m: [[record['geometry'][i] for i in batch] for batch in batches]
                   for m, batches in record['groups'].items()}  # 按编号恢复每种方法、预算的方案记录列表。
        return cls(masks, meta, regions, costs)


