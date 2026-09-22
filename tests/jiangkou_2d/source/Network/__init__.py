"""唯一网架数据源：基础负荷、径向结构和逐线路型号；功率用 kW/kvar。"""
from dataclasses import dataclass  # 声明物理数据字段。
from functools import cached_property  # 固定网架的派生矩阵只计算一次。
import numpy as np  # 统一节点、支路和型号参数的数组顺序。


@dataclass(frozen=True)  # 投资项目的配置保持不可变。
class Project:  # 固定拓扑下的单条支路改造项目。
    """同走廊增设一回相同线路；建成后的等值阻抗为原值的一半。"""
    name: str  # 改造项目名称。
    branch: tuple[int, int]  # 原支路的两个真实端点。
    cost: float  # 项目的增量投资，单位由母网定义。


@dataclass(frozen=True)  # 每条候选支路的型号表不可变。
class LineOptions:  # 紧凑模型直接读取的逐线路型号配置。
    """一条走廊的可选型号；branch 是 planning_branches 中的索引。

    r/reactance 为标幺阻抗，cost 为各型号的增量投资；第 0 型保持原状。
    """
    branch: int  # 支路受端在非根节点数组中的索引。
    r: tuple  # 各型号的标幺电阻。
    reactance: tuple  # 各型号的标幺电抗。
    cost: tuple  # 各型号相对基础网架的新增费用。
    capacity: tuple = ()  # 可选的逐型号送端有功上限；为空时沿用基础支路限值。
    optional: bool = False  # True 允许这条走廊不投入，由主问题决定径向连接。


@dataclass  # 由字段声明生成固定径向网架的初始化函数。
class RadialNetwork:  # 所有运行模型读取的统一径向网架数据。
    """固定方案；阻抗和限值用 p.u.，负荷用 kW/kvar。

    nodes 只列非根节点，load_nodes 是独立变化的负荷节点。
    capacity 是支路送端有功上限；未给定的限值用 inf，不补造额定值。
    """

    name: str  # 所属算例名称。
    root: int  # 固定电压的电源根节点。
    nodes: tuple  # 只含非根节点，决定所有节点数组的顺序。
    branches: object  # 原始支路端点表，初始化时定向为根到叶。
    r: object  # 原始支路顺序下的标幺电阻。
    reactance: object  # 原始支路顺序下的标幺电抗。
    base: float  # 数值为 kVA；P/Q 用相同基准转换成标幺值。
    voltage_kv: float  # 基准线电压，kV。
    original_p: np.ndarray  # 非根节点的原始有功负荷，kW。
    original_q: np.ndarray  # 非根节点的原始无功负荷，kvar。
    load_nodes: tuple  # 用于区域坐标的独立变化节点。
    q_ratio: np.ndarray  # 独立节点的固定 Q/P 比例。
    vmin: object  # 节点电压幅值下限的平方（p.u.²）。
    vmax: object  # 节点电压幅值上限的平方（p.u.²）。
    power_limit: float  # 三个独立节点的总负荷外界（kW），不直接等同于源端容量。
    capacity: object = np.inf  # 支路送端有功上限，p.u.；未给定则不限制。
    source_pmax: float = np.inf  # 源端送出有功上限，p.u.。
    source_qmax: float = np.inf  # 源端送出无功上限，p.u.。
    source_smax: float = np.inf  # 源端视在容量上限，p.u.。
    x: object = ()  # 本固定方案的离散建设记录。
    cost: float = 0.  # 本方案增量投资，基础固定方案默认零费用。
    sources: tuple = ()  # 决定物理输入的源码文件列表。

    cost_unit = "相对投资单位"  # 规划费用的单位由具体网架设置。
    line_options = ()  # 普通固定网架没有规划选项；Case33 在母网中提供逐线路型号。
    budgets = (0., 1., 2., np.inf)  # 网架的默认测试预算，可在 Notebook 中覆盖。

    @property  # 固定拓扑以根向支路作为规划输入。
    def planning_branches(self):  # 返回内部节点索引，−1 表示根节点。
        return tuple((int(parent),i) for i,parent in enumerate(self.parent))  # 可重构算例覆写为包含新建走廊的端点表。

    def __post_init__(self):  # 根据无向支路建立径向顺序并统一所有数组的索引。
        adjacency = {i: [] for i in (self.root, *self.nodes)}  # 先构造无向邻接表，再确定根向拓扑。
        for k, (a, b) in enumerate(self.branches):  # 支路编号用于稍后重排对应的设备参数。
            adjacency[a].append((b, k))  # 向端点 a 记录邻接节点 b 及支路编号。
            adjacency[b].append((a, k))  # 同时记录反方向，使输入走廊无须预先定向。
        parent, edge, order = {self.root: None}, {}, [self.root]  # 记录父节点、入边和从根到叶的遍历次序。
        for node in order:  # 沿从根开始的队列遍历网络。
            for child, k in adjacency[node]:  # 查看当前节点相邻的每条支路。
                if child not in parent:  # 只访问尚未确定父节点的新节点。
                    parent[child], edge[child] = node, k  # 记录该节点的父节点和入边编号。
                    order.append(child)  # 将新节点加入从根到叶的处理顺序。
        assert len(order) == self.n+1 and len(self.branches) == self.n  # 支路递推只适用于连通径向网。
        self.parent = np.array([-1 if parent[i] == self.root else self.nodes.index(parent[i])  # 根端映射为 −1，其余父节点映射为内部索引。
                                for i in self.nodes])  # 内部索引中 -1 表示接电源根节点。
        self.order = [self.nodes.index(i) for i in order[1:]]  # 将真实节点号转换为数组索引。
        indices = [edge[i] for i in self.nodes]  # 第 i 条模型支路必须以第 i 个非根节点为受端。
        self.r, self.reactance = (np.asarray(a)[indices] for a in (self.r, self.reactance))  # 按受端节点重排支路阻抗。
        self.capacity = np.broadcast_to(self.capacity, (self.n,))[indices]  # 支路有功上限使用相同重排。
        self.vmin, self.vmax = (np.broadcast_to(a, (self.n,)) for a in (self.vmin, self.vmax))  # 电压限值按节点给定，不按边重排。

    @property  # 非根节点数按当前节点表读取。
    def n(self):  # 返回径向网络的支路数及非根节点数。
        return len(self.nodes)  # 树中每个非根节点恰有一条入边。

    @cached_property  # 缓存独立负荷在完整节点数组中的位置。
    def selected(self):  # 返回三维负荷坐标对应的内部节点索引。
        return np.array([self.nodes.index(i) for i in self.load_nodes])  # 三个独立负荷节点在全网数组中的位置。

    @cached_property  # 接根支路列表与查询负荷无关，只计算一次。
    def roots(self):  # 找出计算电源送出功率所需的支路。
        return np.flatnonzero(self.parent < 0)  # 电源直接送出的各条支路。

    @cached_property  # 固定拓扑的直接子支路关系只计算一次。
    def children(self):  # 返回各支路受端的下游入边索引。
        return [np.flatnonzero(self.parent == i) for i in range(self.n)]  # 各支路受端节点直接连接的下游支路。

    @cached_property  # 缓存固定拓扑的下游关联矩阵。
    def D(self):  # 构造支路功率汇总与路径压降使用的矩阵 D。
        """D[e,j]=1 表示 j 在支路 e 下游；D.T 表示各节点到根的路径。"""
        matrix = np.eye(self.n)  # 每条支路首先包含自己的受端节点。
        for i in reversed(self.order):  # 先处理叶节点，保证子树信息已完整累加。
            if self.parent[i] >= 0:  # 接根支路没有需要更新的非根父支路。
                matrix[self.parent[i]] += matrix[i]  # 从叶到根汇总子树节点，形成下游关联矩阵 D。
        return matrix  # 返回元素为零或一的下游关联矩阵。

    @cached_property  # 固定的负荷嵌入矩阵只计算一次。
    def E(self):  # 把独立负荷坐标映射到全网节点。
        return np.eye(self.n)[:, self.selected]  # E(n×3) 把三个独立负荷坐标嵌入全网节点向量。

    @cached_property  # 固定背景有功只从原始工况提取一次。
    def fixed_p(self):  # 返回除独立负荷节点以外的有功负荷。
        values = self.original_p.copy()  # 复制原始数组，避免清空可变节点时修改输入工况。
        values[self.selected] = 0.  # 清空可变节点的基准值，防止随后叠加独立负荷时重复计入。
        return values  # 返回单位仍为 kW 的固定有功背景。

    @cached_property  # 固定背景无功只从原始工况提取一次。
    def fixed_q(self):  # 返回除独立负荷节点以外的无功负荷。
        values = self.original_q.copy()  # 复制原始数组，保持原工况可供交叉核验。
        values[self.selected] = 0.  # 这些节点的无功由独立有功和固定 Q/P 比例重新给定。
        return values  # 返回单位仍为 kvar 的固定无功背景。

    def loads(self, power):  # 在不改变背景工况的前提下组装查询负荷。
        """只替换独立坐标，其余负荷不变；返回全网 P/Q（p.u.）。"""
        power = np.asarray(power).reshape(-1, len(self.load_nodes))  # 支持一次传入多个三维负荷样本。
        return ((self.fixed_p+power@self.E.T)/self.base,  # P_node=P_fixed+E*p，kW 转为标幺。
                (self.fixed_q+(power*self.q_ratio)@self.E.T)/self.base)  # Q_node=Q_fixed+E*diag(Q/P)*p。

    def ppc(self, power=None):  # 构造供独立节点潮流程序读取的 MATPOWER 格式数据。
        """标准 MATPOWER 数据，供独立节点导纳矩阵潮流交叉核验。"""
        power = self.original_p[self.selected] if power is None else power  # 未指定查询时采用原始独立节点负荷。
        p, q = self.loads(power)  # 得到含固定背景的完整标幺 P/Q。
        bus = np.zeros((self.n+1, 13))  # MATPOWER 节点表含根节点和全部负荷节点，共 13 列。
        bus[:, 0], bus[:, 1] = (self.root, *self.nodes), 1  # 写入真实节点编号，默认节点类型为 PQ。
        bus[0, 1] = 3  # 根节点设置为平衡节点。
        bus[1:, 2], bus[1:, 3] = p[0]*self.base/1000, q[0]*self.base/1000  # 节点功率从标幺恢复为 MATPOWER 要求的 MW/Mvar。
        bus[:, 6:8], bus[:, 9:11] = 1., [self.voltage_kv, 1.]  # 设置区域、初始电压、电压基准和分区编号。
        bus[:, 11:13] = 1.  # 根节点电压上下限固定为 1 p.u.。
        bus[1:, 11], bus[1:, 12] = np.sqrt(self.vmax), np.sqrt(self.vmin)  # 非根节点电压上下限从平方值恢复为幅值。
        gen = np.zeros((1, 21))  # 建立一台根节点电源的标准发电机表。
        gen[0, [0, 3, 4, 5, 6, 7, 8]] = [self.root, self.source_qmax*self.base/1000,  # 填写源端节点及无功上限，单位为 Mvar。
            -self.source_qmax*self.base/1000, 1., self.base/1000, 1., self.source_pmax*self.base/1000]  # 补入无功下限、电压设定、功率基准、状态和有功上限。
        branch = np.zeros((self.n, 13))  # 为每条在运支路建立标准参数行。
        branch[:, 0] = [self.root if i < 0 else self.nodes[i] for i in self.parent]  # 用定向后的父节点恢复支路真实送端编号。
        branch[:, 1], branch[:, 2], branch[:, 3] = self.nodes, self.r, self.reactance  # 填写受端节点及标幺电阻、电抗。
        branch[:, 10:13] = [1., -360., 360.]  # 所有固定方案支路均闭合，相角差不另加限制。
        return dict(version="2", baseMVA=self.base/1000, bus=bus, gen=gen, branch=branch)  # 返回可直接由独立潮流软件读取的标准数据结构。
