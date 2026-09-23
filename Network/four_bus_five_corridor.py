"""原四节点五走廊基础算例；保留线路、型号与造价，不预先生成建设组合。"""
from dataclasses import dataclass  # 走廊与型号参数使用具名字段，便于审阅。
from functools import cached_property  # 逐线路型号表只从原始数据计算一次。
import numpy as np  # 标幺转换与固定方案数组。
from . import Corridor, TypeParameters, RadialNetwork


@dataclass(frozen=True)  # 原始设备参数不可变。
class LineType:  # 全部走廊共用的 L/M/H 型号。
    name: str  # 型号名称。
    r_ohm_km: float  # 单位长度电阻，Ω/km。
    x_ohm_km: float  # 单位长度电抗，Ω/km。
    capacity_kw: float  # 送端有功上限，kW。
    cable_cny_m: float  # 电缆造价，元/m。


class FourBus(RadialNetwork):  # 基础测试与 case33 都直接交给同一套求解流程。
    corridor_data = (  # ID、端点、长度(m)、既有型号；仅保存原始物理数据。
        ('01', (0, 1), 220., 'L'),
        ('12', (1, 2), 180., 'L'),
        ('13', (1, 3), 240., 'L'),
        ('02', (0, 2), 300., None),
        ('23', (2, 3), 160., None),
    )
    lines = (  # 型号、阻抗、送端有功上限和造价沿用原始算例。
        LineType('L',1.15,.08,35.,28.),  # 基础 L 型。
        LineType('M',.62,.08,65.,43.),  # 中等 M 型。
        LineType('H',.32,.08,100.,65.),  # 加强 H 型。
    )
    power_factor, voltage_kv, voltage_min_pu = .95, .4, .93  # 原功率因数和电压设置。
    transformer_kva, new_corridor_cny_m = 150., 45.  # 原配变容量和新通道单位造价。
    cost_unit, budgets = '元', (20000.,40000.,60000.,np.inf)  # 恢复原四档投资预算。

    def __init__(self):  # 基础运行网是三条既有走廊，规划输入另含全部五条走廊。
        existing = [row for row in self.corridor_data if row[3] is not None]
        zbase = self.voltage_kv**2/(self.transformer_kva/1000.)  # kV²/MVA 得到基准阻抗 Ω。
        line = self.lines[0]  # 既有线路均为 L 型。
        super().__init__(  # 与 33 节点使用相同的基础负荷、限值及运行接口。
            'four_bus_five_corridor',0,(1,2,3),[row[1] for row in existing],
            [line.r_ohm_km*row[2]/1000./zbase for row in existing],
            [line.x_ohm_km*row[2]/1000./zbase for row in existing],
            self.transformer_kva,self.voltage_kv,np.zeros(3),np.zeros(3),(1,2,3),  # 三个独立负荷，无固定背景。
            np.full(3,np.tan(np.arccos(self.power_factor))),self.voltage_min_pu**2,1.,  # Q/P 和电压平方限值。
            self.transformer_kva*self.power_factor,capacity=line.capacity_kw/self.transformer_kva,  # 无损负荷外界及基础线路有功上限。
            source_smax=1.,sources=('Network/__init__.py','Network/four_bus_five_corridor.py'))  # 配变视在容量为 1 p.u.。

    @cached_property
    def corridors(self):
        zbase = self.voltage_kv**2/(self.base/1000.)
        result = []
        for name, endpoints, length, existing_type in self.corridor_data:
            types = tuple(TypeParameters(t.name, t.r_ohm_km*length/1000./zbase,
                t.x_ohm_km*length/1000./zbase, t.capacity_kw/self.base,
                0. if t.name == existing_type else length*(t.cable_cny_m+
                    (self.new_corridor_cny_m if existing_type is None else 0.))) for t in self.lines)
            result.append(Corridor(name, endpoints, existing_type, existing_type is not None,
                                   False, types))
        return tuple(result)


network = FourBus()  # 保留原网架模块可直接导入的入口。
