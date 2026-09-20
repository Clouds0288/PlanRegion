"""唯一网架数据源；原始负荷用 kW/kvar，径向模型阻抗及限值用标幺值。

四节点费用为元；33 节点合成改造费用为相对投资单位，由 cost_unit 标识。
"""
from dataclasses import dataclass
from functools import cached_property
from itertools import combinations, product

import numpy as np


@dataclass(frozen=True)
class LineType:
    name: str
    r_ohm_km: float
    x_ohm_km: float
    capacity_kw: float
    cable_cny_m: float


@dataclass(frozen=True)
class Corridor:
    name: str
    start: int
    end: int
    length_m: float
    existing: bool


@dataclass(frozen=True)
class Project:
    """同走廊增设一回相同线路；建成后的等值阻抗为原值的一半。"""
    name: str
    branch: tuple[int, int]
    cost: float


@dataclass(frozen=True)
class LineOptions:
    """一条在运支路的可选型号；branch 是受端节点在 nodes 中的索引。

    r/reactance 为标幺阻抗，cost 为各型号的增量投资；第 0 型保持原状。
    """
    branch: int
    r: tuple
    reactance: tuple
    cost: tuple


@dataclass(frozen=True)
class Network:
    name: str
    nodes: tuple[int, ...]
    root: int
    corridors: tuple[Corridor, ...]
    lines: tuple[LineType, ...]
    power_factor: float
    voltage_kv: float
    voltage_min_pu: float
    transformer_kva: float
    new_corridor_cny_m: float

    @cached_property
    def load_nodes(self):
        return tuple(i for i in self.nodes if i != self.root)  # 四节点算例中，三个非根节点均为独立负荷。

    @cached_property
    def power_limit(self):
        return self.transformer_kva * self.power_factor  # 无损有功总负荷外界，kVA×功率因数得到 kW。

    @cached_property
    def cost(self):
        # x 按“走廊优先、线型次之”展开；保留现有走廊原型号的费用为零。
        return np.array([
            0. if edge.existing and k == 0 else edge.length_m * (
                line.cable_cny_m + (0. if edge.existing else self.new_corridor_cny_m))
            for edge in self.corridors for k, line in enumerate(self.lines)
        ])

    cost_unit = "元"

    @property
    def sources(self):
        return ("Network/__init__.py", f"Network/{self.name}.py")

    @cached_property
    def designs(self):
        """枚举连通树及设备选型；每个方案转换为统一的径向网架数据。"""
        designs = []
        for edges in combinations(range(len(self.corridors)), len(self.nodes)-1):  # 树必须恰有 |V|-1 条边。
            reached = {self.root}  # 从电源根节点检查选中走廊能否连通全网。
            for _ in self.load_nodes:
                for e in edges:
                    edge = self.corridors[e]
                    if edge.start in reached or edge.end in reached:
                        reached.update((edge.start, edge.end))
            if len(reached) != len(self.nodes):  # 在 |V|-1 条边条件下，全连通等价于生成树。
                continue
            for types in product(range(len(self.lines)), repeat=len(edges)):  # 每条选中走廊独立选择一种设备型号。
                choice = tuple(zip(edges, types))  # (走廊编号, 线型编号) 的离散建设方案。
                x = np.zeros(len(self.cost))  # 建设向量与费用向量采用相同展开次序。
                branches, r, reactance, capacity = [], [], [], []
                zbase = self.voltage_kv**2/(self.transformer_kva/1000)  # Zbase=Ubase²/Sbase，kV²/MVA 得到 Ω。
                for e, k in choice:
                    edge, line = self.corridors[e], self.lines[k]
                    x[e*len(self.lines)+k] = 1  # 本走廊只启用当前选中线型。
                    branches.append((edge.start, edge.end))  # 无向走廊随后统一定向为从根到叶。
                    r.append(line.r_ohm_km*edge.length_m/1000/zbase)  # 长度换成 km，再把 Ω 换成标幺电阻。
                    reactance.append(line.x_ohm_km*edge.length_m/1000/zbase)  # 同样转换电抗。
                    capacity.append(line.capacity_kw/self.transformer_kva)  # 本算例给定的是送端有功容量，不是电流热限。
                designs.append(RadialNetwork(
                    self.name, self.root, self.load_nodes, branches, r, reactance,
                    self.transformer_kva, self.voltage_kv, np.zeros(len(self.load_nodes)),
                    np.zeros(len(self.load_nodes)), self.load_nodes,
                    np.full(len(self.load_nodes), np.tan(np.arccos(self.power_factor))),
                    self.voltage_min_pu**2, 1., self.power_limit, capacity=capacity,
                    source_smax=1., choice=choice, x=x, cost=float(self.cost@x), sources=self.sources))
        return tuple(sorted(designs, key=lambda d: (d.cost, tuple(d.x))))  # 按投资递增排列，为 AC 最小投资扫描提供顺序。


@dataclass
class RadialNetwork:
    """固定方案；阻抗和限值用 p.u.，负荷用 kW/kvar。

    nodes 只列非根节点，load_nodes 是独立变化的负荷节点。
    capacity 是支路送端有功上限；未给定的限值用 inf，不补造额定值。
    """

    name: str
    root: int
    nodes: tuple
    branches: object
    r: object
    reactance: object
    base: float  # 数值为 kVA；P/Q 用相同基准转换成标幺值。
    voltage_kv: float
    original_p: np.ndarray
    original_q: np.ndarray
    load_nodes: tuple
    q_ratio: np.ndarray
    vmin: object  # 节点电压幅值下限的平方（p.u.²）。
    vmax: object  # 节点电压幅值上限的平方（p.u.²）。
    power_limit: float  # 三个独立节点的总负荷外界（kW），不直接等同于源端容量。
    capacity: object = np.inf
    source_pmax: float = np.inf
    source_qmax: float = np.inf
    source_smax: float = np.inf
    choice: tuple = ()
    x: object = ()
    cost: float = 0.
    sources: tuple = ()

    cost_unit = "相对投资单位"

    def __post_init__(self):
        adjacency = {i: [] for i in (self.root, *self.nodes)}  # 先构造无向邻接表，再确定根向拓扑。
        for k, (a, b) in enumerate(self.branches):
            adjacency[a].append((b, k))
            adjacency[b].append((a, k))
        parent, edge, order = {self.root: None}, {}, [self.root]  # 记录父节点、入边和从根到叶的遍历次序。
        for node in order:
            for child, k in adjacency[node]:
                if child not in parent:
                    parent[child], edge[child] = node, k
                    order.append(child)
        assert len(order) == self.n+1 and len(self.branches) == self.n  # 支路递推只适用于连通径向网。
        self.parent = np.array([-1 if parent[i] == self.root else self.nodes.index(parent[i])
                                for i in self.nodes])  # 内部索引中 -1 表示接电源根节点。
        self.order = [self.nodes.index(i) for i in order[1:]]  # 将真实节点号转换为数组索引。
        indices = [edge[i] for i in self.nodes]  # 第 i 条模型支路必须以第 i 个非根节点为受端。
        self.r, self.reactance = (np.asarray(a)[indices] for a in (self.r, self.reactance))  # 按受端节点重排支路阻抗。
        self.capacity = np.broadcast_to(self.capacity, (self.n,))[indices]  # 支路有功上限使用相同重排。
        self.vmin, self.vmax = (np.broadcast_to(a, (self.n,)) for a in (self.vmin, self.vmax))  # 电压限值按节点给定，不按边重排。

    @property
    def n(self):
        return len(self.nodes)

    @property
    def designs(self):
        return (self,)

    @cached_property
    def selected(self):
        return np.array([self.nodes.index(i) for i in self.load_nodes])  # 三个独立负荷节点在全网数组中的位置。

    @cached_property
    def roots(self):
        return np.flatnonzero(self.parent < 0)  # 电源直接送出的各条支路。

    @cached_property
    def children(self):
        return [np.flatnonzero(self.parent == i) for i in range(self.n)]  # 各支路受端节点直接连接的下游支路。

    @cached_property
    def D(self):
        """D[e,j]=1 表示 j 在支路 e 下游；D.T 表示各节点到根的路径。"""
        matrix = np.eye(self.n)  # 每条支路首先包含自己的受端节点。
        for i in reversed(self.order):
            if self.parent[i] >= 0:
                matrix[self.parent[i]] += matrix[i]  # 从叶到根汇总子树节点，形成下游关联矩阵 D。
        return matrix

    @cached_property
    def E(self):
        return np.eye(self.n)[:, self.selected]  # E(n×3) 把三个独立负荷坐标嵌入全网节点向量。

    @cached_property
    def fixed_p(self):
        values = self.original_p.copy()
        values[self.selected] = 0.  # 清空可变节点的基准值，防止随后叠加独立负荷时重复计入。
        return values

    @cached_property
    def fixed_q(self):
        values = self.original_q.copy()
        values[self.selected] = 0.  # 这些节点的无功由独立有功和固定 Q/P 比例重新给定。
        return values

    def loads(self, power):
        """只替换独立坐标，其余负荷不变；返回全网 P/Q（p.u.）。"""
        power = np.asarray(power).reshape(-1, len(self.load_nodes))  # 支持一次传入多个三维负荷样本。
        return ((self.fixed_p+power@self.E.T)/self.base,  # P_node=P_fixed+E*p，kW 转为标幺。
                (self.fixed_q+(power*self.q_ratio)@self.E.T)/self.base)  # Q_node=Q_fixed+E*diag(Q/P)*p。

    def ppc(self, power=None):
        """标准 MATPOWER 数据，供独立节点导纳矩阵潮流交叉核验。"""
        power = self.original_p[self.selected] if power is None else power
        p, q = self.loads(power)
        bus = np.zeros((self.n+1, 13))
        bus[:, 0], bus[:, 1] = (self.root, *self.nodes), 1
        bus[0, 1] = 3
        bus[1:, 2], bus[1:, 3] = p[0]*self.base/1000, q[0]*self.base/1000
        bus[:, 6:8], bus[:, 9:11] = 1., [self.voltage_kv, 1.]
        bus[:, 11:13] = 1.
        bus[1:, 11], bus[1:, 12] = np.sqrt(self.vmax), np.sqrt(self.vmin)
        gen = np.zeros((1, 21))
        gen[0, [0, 3, 4, 5, 6, 7, 8]] = [self.root, self.source_qmax*self.base/1000,
            -self.source_qmax*self.base/1000, 1., self.base/1000, 1., self.source_pmax*self.base/1000]
        branch = np.zeros((self.n, 13))
        branch[:, 0] = [self.root if i < 0 else self.nodes[i] for i in self.parent]
        branch[:, 1], branch[:, 2], branch[:, 3] = self.nodes, self.r, self.reactance
        branch[:, 10:13] = [1., -360., 360.]
        return dict(version="2", baseMVA=self.base/1000, bus=bus, gen=gen, branch=branch)
