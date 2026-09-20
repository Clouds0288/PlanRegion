"""唯一网架数据源；原始负荷用 kW/kvar，径向模型阻抗及限值用标幺值。

四节点费用为元；33 节点合成改造费用为相对投资单位，由 cost_unit 标识。
"""
from dataclasses import dataclass  # 用数据类声明网架与设备参数，避免重复编写初始化赋值。
from functools import cached_property  # 静态派生矩阵和方案表仅计算一次。
from itertools import combinations, product  # 仅为四节点小算例生成树结构和型号组合。

import numpy as np  # 按统一顺序存储物理参数和关联矩阵。


@dataclass(frozen=True)  # 型号参数创建后不可修改。
class LineType:  # 四节点算例中的线路设备型号。
    name: str  # 用于审阅和展示的型号名称。
    r_ohm_km: float  # 单位长度电阻，Ω/km。
    x_ohm_km: float  # 单位长度电抗，Ω/km。
    capacity_kw: float  # 给定的送端有功上限，kW。
    cable_cny_m: float  # 电缆单位长度造价，元/m。


@dataclass(frozen=True)  # 走廊基础数据创建后不可修改。
class Corridor:  # 两节点之间可选的线路走廊。
    name: str  # 走廊编号。
    start: int  # 走廊一个端点的真实节点号。
    end: int  # 走廊另一个端点的真实节点号。
    length_m: float  # 走廊实际长度，m。
    existing: bool  # 是否为既有走廊，用于区分新增通道费用。


@dataclass(frozen=True)  # 投资项目的配置保持不可变。
class Project:  # 固定拓扑下的单条支路改造项目。
    """同走廊增设一回相同线路；建成后的等值阻抗为原值的一半。"""
    name: str  # 改造项目名称。
    branch: tuple[int, int]  # 原支路的两个真实端点。
    cost: float  # 项目的增量投资，单位由母网定义。


@dataclass(frozen=True)  # 每条候选支路的型号表不可变。
class LineOptions:  # 紧凑模型直接读取的逐线路型号配置。
    """一条在运支路的可选型号；branch 是受端节点在 nodes 中的索引。

    r/reactance 为标幺阻抗，cost 为各型号的增量投资；第 0 型保持原状。
    """
    branch: int  # 支路受端在非根节点数组中的索引。
    r: tuple  # 各型号的标幺电阻。
    reactance: tuple  # 各型号的标幺电抗。
    cost: tuple  # 各型号相对基础网架的新增费用。


@dataclass(frozen=True)  # 四节点母网的输入配置不可变。
class Network:  # 含候选走廊和设备型号的小规模母网配置。
    name: str  # 算例名称，也用于选择网架文件及输出目录。
    nodes: tuple[int, ...]  # 包含电源根节点的真实节点号。
    root: int  # 固定电压电源的节点号。
    corridors: tuple[Corridor, ...]  # 已有和候选走廊列表。
    lines: tuple[LineType, ...]  # 各走廊可选的线路型号列表。
    power_factor: float  # 独立负荷的统一功率因数。
    voltage_kv: float  # 线电压基准，kV。
    voltage_min_pu: float  # 节点电压幅值下限，p.u.。
    transformer_kva: float  # 电源变压器视在容量，kVA；本小算例也用作标幺基准。
    new_corridor_cny_m: float  # 新建走廊的额外单位长度费用，元/m。

    @cached_property  # 缓存不会随运行查询变化的负荷节点列表。
    def load_nodes(self):  # 返回四节点算例的全部非根节点。
        return tuple(i for i in self.nodes if i != self.root)  # 四节点算例中，三个非根节点均为独立负荷。

    @cached_property  # 无损总负荷界只需由固定设备参数计算一次。
    def power_limit(self):  # 返回独立有功负荷的共同外界。
        return self.transformer_kva * self.power_factor  # 无损有功总负荷外界，kVA×功率因数得到 kW。

    @cached_property  # 走廊型号费用只由网架数据生成一次。
    def cost(self):  # 按建设向量相同次序计算各选项的增量投资。
        # x 按“走廊优先、线型次之”展开；保留现有走廊原型号的费用为零。
        return np.array([  # 返回与走廊和型号展开顺序一致的费用向量。
            0. if edge.existing and k == 0 else edge.length_m * (  # 保留既有走廊原型号免费，其他选择按长度计价。
                line.cable_cny_m + (0. if edge.existing else self.new_corridor_cny_m))  # 新建走廊额外计入通道费用。
            for edge in self.corridors for k, line in enumerate(self.lines)  # 先遍历走廊，再遍历该走廊的所有型号。
        ])  # 完成一维建设费用数组。

    cost_unit = "元"  # 明确四节点实验投资的物理单位。

    @property  # 源码列表由算例名称即时组成，无需独立配置副本。
    def sources(self):  # 列出决定本算例物理输入的文件。
        return ("Network/__init__.py", f"Network/{self.name}.py")  # 实验结果记录这些文件的指纹以便复现。

    @cached_property  # 完整方案表仅为小算例对照生成一次。
    def designs(self):  # 枚举可连接全部节点的树及其线路型号。
        """枚举连通树及设备选型；每个方案转换为统一的径向网架数据。"""
        designs = []  # 收集合规的固定径向方案。
        for edges in combinations(range(len(self.corridors)), len(self.nodes)-1):  # 树必须恰有 |V|-1 条边。
            reached = {self.root}  # 从电源根节点检查选中走廊能否连通全网。
            for _ in self.load_nodes:  # 至多传播非根节点数轮即可完成连通性检查。
                for e in edges:  # 每轮沿所有选中走廊扩展已到达节点。
                    edge = self.corridors[e]  # 读取当前走廊的端点。
                    if edge.start in reached or edge.end in reached:  # 至少一个端点可达时，另一个端点也可达。
                        reached.update((edge.start, edge.end))  # 同时将两端加入已到达集合。
            if len(reached) != len(self.nodes):  # 在 |V|-1 条边条件下，全连通等价于生成树。
                continue  # 不连通组合不属于合法径向方案。
            for types in product(range(len(self.lines)), repeat=len(edges)):  # 每条选中走廊独立选择一种设备型号。
                choice = tuple(zip(edges, types))  # (走廊编号, 线型编号) 的离散建设方案。
                x = np.zeros(len(self.cost))  # 建设向量与费用向量采用相同展开次序。
                branches, r, reactance, capacity = [], [], [], []  # 收集该方案的支路、阻抗和给定有功容量。
                zbase = self.voltage_kv**2/(self.transformer_kva/1000)  # Zbase=Ubase²/Sbase，kV²/MVA 得到 Ω。
                for e, k in choice:  # 按选中的走廊与型号实例化物理参数。
                    edge, line = self.corridors[e], self.lines[k]  # 同时读取走廊长度和设备单位参数。
                    x[e*len(self.lines)+k] = 1  # 本走廊只启用当前选中线型。
                    branches.append((edge.start, edge.end))  # 无向走廊随后统一定向为从根到叶。
                    r.append(line.r_ohm_km*edge.length_m/1000/zbase)  # 长度换成 km，再把 Ω 换成标幺电阻。
                    reactance.append(line.x_ohm_km*edge.length_m/1000/zbase)  # 同样转换电抗。
                    capacity.append(line.capacity_kw/self.transformer_kva)  # 本算例给定的是送端有功容量，不是电流热限。
                designs.append(RadialNetwork(  # 将一个完整建设组合转换为统一径向数据。
                    self.name, self.root, self.load_nodes, branches, r, reactance,  # 传入根节点、非根节点、选中支路及标幺阻抗。
                    self.transformer_kva, self.voltage_kv, np.zeros(len(self.load_nodes)),  # 指定功率基准、电压基准和初始有功负荷。
                    np.zeros(len(self.load_nodes)), self.load_nodes,  # 该小算例没有固定背景无功，三个节点均独立变化。
                    np.full(len(self.load_nodes), np.tan(np.arccos(self.power_factor))),  # 由功率因数计算统一的 Q/P 比例。
                    self.voltage_min_pu**2, 1., self.power_limit, capacity=capacity,  # 电压幅值限值先平方，支路有功容量沿用给定值。
                    source_smax=1., choice=choice, x=x, cost=float(self.cost@x), sources=self.sources))  # 源端视在上限为 1 p.u.，同时记录本方案投资与输入来源。
        return tuple(sorted(designs, key=lambda d: (d.cost, tuple(d.x))))  # 按投资递增排列，为 AC 最小投资扫描提供顺序。


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
    choice: tuple = ()  # 四节点小算例的走廊与型号组合。
    x: object = ()  # 本固定方案的离散建设记录。
    cost: float = 0.  # 本方案增量投资，基础固定方案默认零费用。
    sources: tuple = ()  # 决定物理输入的源码文件列表。

    cost_unit = "相对投资单位"  # 固定网架默认使用相对投资单位；四节点母网另行标为元。
    line_options = ()  # 普通固定网架没有规划选项；Case33 在母网中提供逐线路型号。

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

    @property  # 固定网架只有当前一个合法方案。
    def designs(self):  # 提供小算例基准所需的统一方案接口。
        return (self,)  # 普通固定网架不产生额外建设组合。

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
