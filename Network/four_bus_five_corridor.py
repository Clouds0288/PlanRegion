"""四节点五走廊算例；线路参数与造价沿用原始数据。"""
from dataclasses import dataclass

import numpy as np

from . import Corridor, Network, TypeParameters


@dataclass(frozen=True)
class LineType:
    name: str
    r_ohm_km: float
    x_ohm_km: float
    capacity_kw: float
    cable_cny_m: float


class FourBus(Network):
    corridor_data = (
        ('01', (0, 1), 220., 'L'),
        ('12', (1, 2), 180., 'L'),
        ('13', (1, 3), 240., 'L'),
        ('02', (0, 2), 300., None),
        ('23', (2, 3), 160., None),
    )
    lines = (
        LineType('L', 1.15, .08, 35., 28.),
        LineType('M', .62, .08, 65., 43.),
        LineType('H', .32, .08, 100., 65.),
    )
    power_factor, voltage_kv, voltage_min_pu = .95, .4, .93
    transformer_kva, new_corridor_cny_m = 150., 45.
    cost_unit, budgets = '元', (20000., 40000., 60000., np.inf)

    def __init__(self):
        base = self.transformer_kva
        zbase = self.voltage_kv**2/(base/1000.)
        corridors = []
        for name, endpoints, length, existing in self.corridor_data:
            types = tuple(TypeParameters(t.name, t.r_ohm_km*length/1000./zbase,
                t.x_ohm_km*length/1000./zbase, t.capacity_kw/base,
                0. if t.name == existing else length*(t.cable_cny_m+
                    (self.new_corridor_cny_m if existing is None else 0.))) for t in self.lines)
            corridors.append(Corridor(name, endpoints, existing, existing is not None, types))
        super().__init__('four_bus_five_corridor', 0, (1, 2, 3), tuple(corridors),
            base, self.voltage_kv, np.zeros(3), np.zeros(3), (1, 2, 3),
            np.full(3, np.tan(np.arccos(self.power_factor))), self.voltage_min_pu**2, 1.,
            base*self.power_factor, source_smax=1.,
            sources=('Network/__init__.py', 'Network/four_bus_five_corridor.py'))


network = FourBus()
