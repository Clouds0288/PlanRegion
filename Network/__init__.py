"""网架数据；节点、走廊、设备和费用统一使用 kW、kV、m、元。"""
from dataclasses import dataclass
from functools import cached_property

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
        return tuple(i for i in self.nodes if i != self.root)

    @cached_property
    def power_limit(self):
        return self.transformer_kva * self.power_factor

    @cached_property
    def cost(self):
        return np.array([
            0. if edge.existing and k == 0 else edge.length_m * (
                line.cable_cny_m + (0. if edge.existing else self.new_corridor_cny_m))
            for edge in self.corridors for k, line in enumerate(self.lines)
        ])
