"""候选网架与选定运行树；节点、走廊、型号分别使用唯一索引。"""
from dataclasses import asdict, dataclass
from functools import cached_property
import hashlib
import json

import numpy as np
from scipy import sparse


@dataclass(frozen=True)
class Project:
    """同走廊并联一回相同线路；费用为增量投资。"""
    name: str
    branch: tuple[int, int]
    cost: float


@dataclass(frozen=True)
class TypeParameters:
    """阻抗、送端有功容量用 p.u.；费用单位由网架定义。"""
    id: str
    r: float
    reactance: float
    capacity: float
    investment_cost: float


@dataclass(frozen=True)
class Corridor:
    id: str
    endpoints: tuple
    existing_type: str | None
    initial_active: bool
    types: tuple[TypeParameters, ...]


@dataclass
class Network:
    """完整候选图。初始开合状态只用于生成初始方案，不限制规划。

    nodes 仅含非根节点；required 指定必须接入的节点。
    型号按走廊顺序连续展开，x/P/Q/ell 与 r/reactance/capacity/cost 共用该顺序。
    负荷用 kW/kvar，base 用 kVA，vmin/vmax 为电压标幺值的平方。
    """
    name: str
    root: object
    nodes: tuple
    corridors: tuple[Corridor, ...]
    base: float
    voltage_kv: float
    original_p: np.ndarray
    original_q: np.ndarray
    load_nodes: tuple
    q_ratio: np.ndarray
    vmin: object
    vmax: object
    power_limit: float
    required: object = True
    source_pmax: float = np.inf
    source_qmax: float = np.inf
    source_smax: float = np.inf
    sources: tuple = ()
    road_allowed: object = True

    cost_unit = '相对投资单位'
    budgets = (0., 1., 2., np.inf)

    def __post_init__(self):
        self.nodes, self.corridors, self.load_nodes = tuple(self.nodes), tuple(self.corridors), tuple(self.load_nodes)
        if self.root in self.nodes or len(set(self.nodes)) != self.n:
            raise ValueError('nodes must be unique and exclude the root')
        self.node_index = {node: i for i, node in enumerate(self.nodes)}
        self.node_index[self.root] = -1
        self.selected = np.array([self.node_index[node] for node in self.load_nodes], dtype=int)
        if np.any(self.selected < 0) or len(set(self.load_nodes)) != len(self.load_nodes):
            raise ValueError('load_nodes must be distinct nonroot nodes')
        for name in ('original_p', 'original_q', 'vmin', 'vmax'):
            setattr(self, name, np.broadcast_to(getattr(self, name), (self.n,)).astype(float))
        self.q_ratio = np.broadcast_to(self.q_ratio, (len(self.load_nodes),)).astype(float)
        self.required = np.broadcast_to(self.required, (self.n,)).astype(bool)
        self.required |= (self.original_p != 0.) | (self.original_q != 0.)
        self.required[self.selected] = True
        if len({c.id for c in self.corridors}) != self.n_corridors:
            raise ValueError('corridor IDs must be unique')
        self.road_allowed = np.broadcast_to(self.road_allowed, (self.n_corridors,)).astype(bool).copy()
        senders, receivers, keys, blocks, parameters = [], [], [], [], []
        for corridor in self.corridors:
            a, b = (self.node_index[node] for node in corridor.endpoints)
            if a == b or not corridor.types or len({t.id for t in corridor.types}) != len(corridor.types):
                raise ValueError(f'invalid corridor: {corridor.id}')
            if corridor.initial_active and corridor.existing_type not in {t.id for t in corridor.types}:
                raise ValueError(f'missing initial type: {corridor.id}')
            a, b = (b, a) if b < 0 else (a, b)
            senders.append(a)
            receivers.append(b)
            blocks.append(slice(len(keys), len(keys)+len(corridor.types)))
            keys.extend((corridor.id, t.id) for t in corridor.types)
            parameters.extend(corridor.types)
        self.senders, self.receivers = np.array(senders), np.array(receivers)
        self.type_keys, self.type_slices = tuple(keys), tuple(blocks)
        self.type_corridor = np.repeat(np.arange(self.n_corridors), [len(c.types) for c in self.corridors])
        self.r, self.reactance, self.capacity, self.cost = (
            np.array([getattr(t, name) for t in parameters], dtype=float)
            for name in ('r', 'reactance', 'capacity', 'investment_cost'))

    @property
    def n(self):
        return len(self.nodes)

    @property
    def n_corridors(self):
        return len(self.corridors)

    @property
    def n_types(self):
        return len(self.type_keys)

    @cached_property
    def corridor_types(self):
        """走廊投入状态 z = corridor_types @ x。"""
        return sparse.csr_matrix((np.ones(self.n_types), (self.type_corridor, np.arange(self.n_types))),
                                 shape=(self.n_corridors, self.n_types))

    @cached_property
    def receiving(self):
        return sparse.csr_matrix((np.ones(self.n_corridors), (self.receivers, np.arange(self.n_corridors))),
                                 shape=(self.n, self.n_corridors))

    @cached_property
    def sending(self):
        edges = np.flatnonzero(self.senders >= 0)
        return sparse.csr_matrix((np.ones(len(edges)), (self.senders[edges], edges)),
                                 shape=(self.n, self.n_corridors))

    @cached_property
    def incidence(self):
        return self.receiving-self.sending

    @cached_property
    def E(self):
        return sparse.csr_matrix((np.ones(len(self.selected)), (self.selected, np.arange(len(self.selected)))),
                                 shape=(self.n, len(self.selected)))

    @cached_property
    def fixed_p(self):
        values = self.original_p.copy()
        values[self.selected] = 0.
        return values

    @cached_property
    def fixed_q(self):
        values = self.original_q.copy()
        values[self.selected] = 0.
        return values

    def loads(self, power):
        power = np.asarray(power).reshape(-1, len(self.load_nodes))
        return ((self.fixed_p+power@self.E.T)/self.base,
                (self.fixed_q+(power*self.q_ratio)@self.E.T)/self.base)

    @property
    def initial_plan(self):
        return {c.id: c.existing_type if c.initial_active else None for c in self.corridors}

    def encode_plan(self, plan):
        """I/O 边界：完整走廊方案转为型号向量。"""
        if set(plan) != {c.id for c in self.corridors}:
            raise ValueError('plan must specify every corridor')
        for corridor in self.corridors:
            if plan[corridor.id] not in {None, *(t.id for t in corridor.types)}:
                raise ValueError(f'unknown type for corridor {corridor.id}')
        return np.array([plan[c] == t for c, t in self.type_keys], dtype=int)

    def decode_plan(self, x):
        x = self._binary_selection(x)
        return {c.id: next((t.id for t, chosen in zip(c.types, x[s]) if chosen), None)
                for c, s in zip(self.corridors, self.type_slices)}

    def _binary_selection(self, x):
        x = np.asarray(x)
        if x.shape != (self.n_types,) or not np.isin(x, (0, 1)).all():
            raise ValueError('x must be a binary vector in type order')
        if np.any(self.corridor_types@x > 1):
            raise ValueError('at most one type may be selected per corridor')
        return x

    def tree(self, x):
        return OperatingTree(self, self._binary_selection(x))

    @property
    def fingerprint(self):
        payload = asdict(self)
        payload['cost_unit'] = self.cost_unit
        encoded = json.dumps(payload, sort_keys=True, default=lambda v: v.tolist(), separators=(',', ':'))
        return hashlib.sha256(encoded.encode()).hexdigest()


class OperatingTree:
    """选定方案的根向视图；节点和型号索引指回 Network，不复制候选图。"""

    def __init__(self, network, x):
        self.network = net = network
        chosen = np.flatnonzero(x)
        adjacency = [[] for _ in range(net.n+1)]  # root 的内部索引 -1 同时指向末项。
        for k in chosen:
            e = net.type_corridor[k]
            a, b = net.senders[e], net.receivers[e]
            adjacency[a].append((b, k))
            adjacency[b].append((a, k))
        parents, edges, order = {-1: -1}, {}, [-1]
        for node in order:
            for child, k in adjacency[node]:
                if child not in parents:
                    parents[child], edges[child] = node, k
                    order.append(child)
        if len(chosen) != len(order)-1 or any(i not in parents for i in np.flatnonzero(net.required)):
            raise ValueError('selection must be a rooted tree covering every required node')
        self.node_indices = np.array(sorted(order[1:]), dtype=int)
        local = {node: i for i, node in enumerate(self.node_indices)}
        local[-1] = -1
        self.parent = np.array([local[parents[i]] for i in self.node_indices], dtype=int)
        self.order = np.array([local[i] for i in order[1:]], dtype=int)
        self.type_indices = np.array([edges[i] for i in self.node_indices], dtype=int)
        self.direction = np.where(net.receivers[net.type_corridor[self.type_indices]] == self.node_indices, 1, -1)
        self.nodes = tuple(net.nodes[i] for i in self.node_indices)

    @property
    def n(self):
        return len(self.nodes)

    @property
    def r(self):
        return self.network.r[self.type_indices]

    @property
    def reactance(self):
        return self.network.reactance[self.type_indices]

    @property
    def capacity(self):
        return self.network.capacity[self.type_indices]

    @property
    def vmin(self):
        return self.network.vmin[self.node_indices]

    @property
    def vmax(self):
        return self.network.vmax[self.node_indices]

    @property
    def base(self):
        return self.network.base

    @property
    def load_nodes(self):
        return self.network.load_nodes

    @property
    def q_ratio(self):
        return self.network.q_ratio

    @property
    def source_pmax(self):
        return self.network.source_pmax

    @property
    def source_qmax(self):
        return self.network.source_qmax

    @property
    def source_smax(self):
        return self.network.source_smax

    @property
    def power_limit(self):
        return self.network.power_limit

    @property
    def cost(self):
        return float(self.network.cost[self.type_indices].sum())

    @property
    def fixed_p(self):
        return self.network.fixed_p[self.node_indices]

    @property
    def fixed_q(self):
        return self.network.fixed_q[self.node_indices]

    @cached_property
    def E(self):
        return self.network.E[self.node_indices].toarray()

    @cached_property
    def roots(self):
        return np.flatnonzero(self.parent < 0)

    @cached_property
    def children(self):
        return [np.flatnonzero(self.parent == i) for i in range(self.n)]

    @cached_property
    def D(self):
        """D[i,j]=1 表示节点 j 在入边 i 的下游。"""
        matrix = np.eye(self.n)
        for i in reversed(self.order):
            if self.parent[i] >= 0:
                matrix[self.parent[i]] += matrix[i]
        return matrix

    def loads(self, power):
        return tuple(values[:, self.node_indices] for values in self.network.loads(power))

    def ppc(self, power=None):
        """MATPOWER 数据供独立节点导纳潮流验证；内部节点统一编号 0..n。"""
        net = self.network
        power = net.original_p[net.selected] if power is None else power
        p, q = self.loads(power)
        bus = np.zeros((self.n+1, 13))
        bus[:, 0], bus[:, 1] = np.arange(self.n+1), 1
        bus[0, 1] = 3
        bus[1:, 2], bus[1:, 3] = p[0]*net.base/1000, q[0]*net.base/1000
        bus[:, 6:8], bus[:, 9:11], bus[:, 11:13] = 1., [net.voltage_kv, 1.], 1.
        bus[1:, 11], bus[1:, 12] = np.sqrt(self.vmax), np.sqrt(self.vmin)
        gen = np.zeros((1, 21))
        gen[0, [0, 3, 4, 5, 6, 7, 8]] = [0, net.source_qmax*net.base/1000,
            -net.source_qmax*net.base/1000, 1., net.base/1000, 1., net.source_pmax*net.base/1000]
        branch = np.zeros((self.n, 13))
        branch[:, 0], branch[:, 1] = self.parent+1, np.arange(1, self.n+1)
        branch[:, 2], branch[:, 3], branch[:, 10:13] = self.r, self.reactance, [1., -360., 360.]
        return dict(version='2', baseMVA=net.base/1000, bus=bus, gen=gen, branch=branch)
