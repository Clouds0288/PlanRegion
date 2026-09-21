"""独立 AC 交叉核验、同网格 FR/MR 和单文件结果；未确定点始终保留为 0。"""
from dataclasses import dataclass  # 一个结果容器集中保存基础数据。
from pathlib import Path  # 组织输出文件路径。
from io import BytesIO  # 解码前释放结果文件句柄。
import json  # 配置与方法计时仅序列化一次。
import numpy as np  # 处理三态网格及派生指标。
from model import ACPowerFlow  # 独立 AC 不读取规划方程矩阵。

METHODS = ('socp','hybrid','ac','linear')  # 方法轴顺序；独立 AC 位于索引 2。
METHOD_NAMES = ('纯 SOCP 切割','线性 + SOCP 精修','AC 数值参考','纯线性切割')  # 展示名称保持相同顺序。


def validate_power_flow(network):  # 用另一套节点导纳矩阵算法交叉核验支路 AC 实现。
    """用六个工况，将独立支路递推与节点导纳矩阵 Newton–Raphson 潮流核对。"""
    import warnings  # 仅在转换原始数据时屏蔽已知的接口弃用提示。
    import pandapower as pp  # 采用独立实现的 Newton–Raphson 潮流。
    from pandapower.converter.pypower.from_ppc import from_ppc  # 将标准 MATPOWER 数据转换为 pandapower 网架。

    reference = ACPowerFlow(network)  # 被核验的是独立 AC 递推，不是 SOCP 方程。
    maximum_difference = 0.  # 记录所有工况、所有节点的最大电压幅值差。
    for scale in (1., 0., .5, 1.2, 1.5, 2.):  # 覆盖原始、零独立负荷及多个放大工况。
        power = scale*network.original_p[network.selected]  # 只缩放三个独立节点，其他背景负荷固定。
        with warnings.catch_warnings():  # 警告过滤仅在本次格式转换的作用域内有效。
            warnings.filterwarnings('ignore', category=FutureWarning,  # 只过滤接口的未来弃用警告。
                                    module='pandapower.converter.pypower.from_ppc')  # 限定来源模块，不屏蔽潮流求解异常。
            net = from_ppc(network.ppc(power), f_hz=50)  # 同一物理网架转换成节点导纳矩阵模型。
        pp.runpp(net, algorithm='nr', tolerance_mva=1e-10, numba=False)  # 另一实现的完整 AC 潮流作为交叉核验。
        # 此处只比较潮流方程；即使电压越限也继续收敛，不能复用 classify 的提前不可行判定。
        ell = np.zeros((1, network.n))  # 支路递推从零电流平方开始。
        for _ in range(160):  # 逐次满足完整 AC 电流等式。
            P, Q, v, u = reference.state(power, ell)  # 当前电流下的送端功率及两端电压平方。
            new = (P*P+Q*Q)/u  # 强制完整 AC 电流等式。
            if np.max(np.abs(new-ell)) < 1e-13:  # 电流平方增量达到比运行限值判定更严格的精度。
                break  # 电流更新已经收敛，停止迭代。
            ell = new  # 更新电流，进行下一次完整功率平衡与压降计算。
        difference = float(np.max(np.abs(np.sqrt(v[0])-net.res_bus.loc[list(network.nodes), 'vm_pu'])))  # v 是平方，NR 输出是幅值。
        assert difference < 1e-8, 'Branch AC and nodal AC disagree'  # 不让方程实现不一致的参考模型进入正式评价。
        maximum_difference = max(maximum_difference, difference)  # 保留全部工况中的最大电压幅值误差。
        if scale == 1:  # 保存原始基础工况的物理校验值。
            baseline = dict(minimum_voltage_pu=float(np.sqrt(v.min())),  # 全网最低电压幅值。
                            minimum_voltage_bus=int(network.nodes[np.argmin(v)]),  # 对应真实节点号。
                            active_loss_kw=float((ell*network.r).sum()*network.base))  # Σr*ell，从标幺还原为 kW。
    return dict(baseline=baseline, nodal_cross_checks=6,  # 返回原始工况和交叉核验次数。
                maximum_voltage_difference_pu=maximum_difference, pandapower=pp.__version__)


def disagreement(approximation, ac):  # 根据同一网格中的多余和遗漏体素计算 FR/MR。
    """同一等体积网格上的 FR/MR；空域的条件比例无定义。"""
    extra = int(np.count_nonzero(approximation & ~ac))  # 黄色：计算域有、AC 参考域没有。
    missed = int(np.count_nonzero(ac & ~approximation))  # 红色：AC 参考域有、计算域遗漏。
    computed, reference = int(approximation.sum()), int(ac.sum())  # 两个域各自的网格单元数量。
    union = computed+missed  # 计算域∪AC 域；遗漏部分与计算域不相交。
    return dict(fr_percent=100*extra/computed if computed else None,  # FR 分母是计算域。
                mr_percent=100*missed/reference if reference else None,  # MR 分母是 AC 参考域。
                region_error_percent=100*(extra+missed)/union if union else 0.,  # 原区域不一致率：对称差/并集。
                computed_cells=computed, ac_cells=reference, extra_cells=extra, missed_cells=missed)


def disagreement_interval(approximation, ac):  # 三态网格中，未确定点只能给误差区间。
    """输入 -1/0/1 分别表示已证域外、未确定、已证域内；区间不含网格离散误差。"""
    def ratio_bounds(left, right):  # 求 |left\right|/|left| 的保守上下界。
        inside, unknown = left==1, left==0  # 已知域内点必须计入分母，未确定点可取域内或域外。
        certain = np.count_nonzero(inside & (right==-1))  # 已确认的多余点不能被未知标签消除。
        lower_denominator = inside.sum()+np.count_nonzero(unknown & (right!=-1))  # 只增加可能重合的点，使比例最小。
        possible = np.count_nonzero((left!=-1) & (right!=1))  # 所有可能成为多余的点。
        upper_denominator = inside.sum()+np.count_nonzero(unknown & (right!=1))  # 为最大比例只增加可能多余的点。
        return [100*certain/lower_denominator if lower_denominator else 0.,  # 可能空域时不给虚假的严格正下界。
                100*possible/upper_denominator if upper_denominator else 100.]  # 未知空分母使用保守上界。
    unknown = int(np.count_nonzero((approximation==0)|(ac==0)))  # 比较双方中任意一方未确定的单元数。
    exact = disagreement(approximation==1,ac==1) if unknown==0 else dict(fr_percent=None,mr_percent=None)  # 只有标签齐全才报告单个 FR/MR 值。
    fr = ratio_bounds(approximation,ac) if unknown else ([exact['fr_percent']]*2 if exact['fr_percent'] is not None else None)  # 完整空域的比例无定义，不能显示成 0–100%。
    mr = ratio_bounds(ac,approximation) if unknown else ([exact['mr_percent']]*2 if exact['mr_percent'] is not None else None)  # MR 交换两个集合使用同一定义。
    return dict(**exact,fr_interval=fr,mr_interval=mr,unknown_cells=unknown)


@dataclass  # 仅保存三态网格与复现配置，不维护方案表或逐查询日志。
class BenchmarkResult:  # 指标与颜色由唯一基础数据推导。
    states: np.ndarray  # (方法,预算,N,N,N)，-1 域外、0 未确定、1 域内。
    metadata: dict  # 网架、预算、评价箱、源码指纹和各方法总耗时。

    @property  # 将 JSON 中的无限预算恢复为数值。
    def budgets(self):  # 返回计算使用的预算序列。
        return tuple(np.inf if b is None else b for b in self.metadata['budgets'])  # null 表示无限预算。

    @property  # 坐标上界只在配置中保存一次。
    def bounds(self):  # 返回公共评价箱，单位 kW。
        return np.asarray(self.metadata['bounds'])  # 与三维坐标顺序一致。

    @property  # 负荷节点直接读取配置。
    def load_nodes(self):  # 返回三条坐标轴对应的实际节点。
        return tuple(self.metadata['load_nodes'])  # 不另维护节点副本。

    @property  # 步长由评价箱和分辨率推导。
    def spacing(self):  # 返回三个方向的单元宽度。
        return self.bounds/self.metadata['divisions']  # 等体积网格的长度尺度。

    @property  # 绘图颜色与 FR/MR 使用相同三态网格。
    def labels(self):  # 0 域外、1 多余、2 遗漏、3 重合、4 未确定。
        inside = self.states==1  # 只有认证域内点才进入可行集合。
        labels = inside.astype(np.uint8)+2*inside[2].astype(np.uint8)  # 计算域和 AC 域的二位组合。
        labels[(self.states==0)|(self.states[2]==0)] = 4  # 任一方未知时不绘成红黄误差。
        return labels  # 颜色数组仅按需生成，不重复落盘。

    @property  # 汇总表从原始标签现场计算。
    def summary(self):  # 每个预算分别与相同预算的 AC 域比较。
        return [dict(method=method,budget=budget,  # 保留方法与预算以便筛选展示。
                     **disagreement_interval(self.states[k,j],self.states[2,j]),  # 未知点只能给出 FR/MR 区间。
                     total_seconds=self.metadata['seconds'][method])  # 一个方法的全部预算共享这次总耗时。
                for k,method in enumerate(METHODS) for j,budget in enumerate(self.metadata['budgets'])]  # 不保存重复的派生表。

    def save(self, folder):  # 一次实验只写一个原始结果文件。
        folder = Path(folder)  # 统一路径表示。
        folder.mkdir(parents=True,exist_ok=True)  # 创建本次实验目录。
        temporary = folder/'result.npz.tmp'  # 完整写入后再替换正式文件。
        with temporary.open('wb') as stream:  # 中断写入不破坏上一次完整结果。
            np.savez_compressed(stream,states=self.states,metadata=json.dumps(self.metadata,ensure_ascii=False))  # 仅存三态标签与唯一配置。
        temporary.replace(folder/'result.npz')  # 原子发布本次结果。

    @classmethod  # 无需建立空实例即可读取结果。
    def load(cls, folder):  # 只读取当前简化格式，不保留旧方案表兼容分支。
        with np.load(BytesIO((Path(folder)/'result.npz').read_bytes()),allow_pickle=False) as data:  # 释放文件句柄后解码。
            return cls(data['states'],json.loads(str(data['metadata'])))  # 指标和图形均重新从基础数据生成。
