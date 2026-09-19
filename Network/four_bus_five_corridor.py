"""四节点五走廊合成算例：根节点 0，负荷节点 1/2/3，设备 L/M/H。"""
from . import Corridor, LineType, Network


network = Network(
    name="four_bus_five_corridor",
    nodes=(0, 1, 2, 3),
    root=0,
    corridors=(
        Corridor("01", 0, 1, 220., True),
        Corridor("12", 1, 2, 180., True),
        Corridor("13", 1, 3, 240., True),
        Corridor("02", 0, 2, 300., False),
        Corridor("23", 2, 3, 160., False),
    ),
    lines=(
        LineType("L", 1.15, .08, 35., 28.),
        LineType("M", .62, .08, 65., 43.),
        LineType("H", .32, .08, 100., 65.),
    ),
    power_factor=.95,
    voltage_kv=.4,
    voltage_min_pu=.93,
    transformer_kva=150.,
    new_corridor_cny_m=45.,
)
