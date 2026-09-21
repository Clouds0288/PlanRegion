"""原四节点五走廊基础算例；保留线路、型号与造价，不预先生成建设组合。"""
from dataclasses import dataclass  # 走廊与型号参数使用具名字段，便于审阅。
from functools import cached_property  # 逐线路型号表只从原始数据计算一次。
import numpy as np  # 标幺转换与固定方案数组。
from . import LineOptions, RadialNetwork  # 与 case33 共用规划输入和径向运行接口。


@dataclass(frozen=True)  # 原始设备参数不可变。
class LineType:  # 全部走廊共用的 L/M/H 型号。
    name: str  # 型号名称。
    r_ohm_km: float  # 单位长度电阻，Ω/km。
    x_ohm_km: float  # 单位长度电抗，Ω/km。
    capacity_kw: float  # 送端有功上限，kW。
    cable_cny_m: float  # 电缆造价，元/m。


@dataclass(frozen=True)  # 原始走廊参数不可变。
class Corridor:  # 既有与新建走廊使用同一数据结构。
    name: str  # 走廊编号。
    start: int  # 一个端点的真实节点号。
    end: int  # 另一个端点的真实节点号。
    length_m: float  # 长度，m。
    existing: bool  # 既有走廊不收新增通道费。


class FourBus(RadialNetwork):  # 基础测试与 case33 都直接交给同一套求解流程。
    corridors = (  # 原五条走廊完整保留；既有线路也可退出所选运行树。
        Corridor('01',0,1,220.,True),  # 既有 0–1。
        Corridor('12',1,2,180.,True),  # 既有 1–2。
        Corridor('13',1,3,240.,True),  # 既有 1–3。
        Corridor('02',0,2,300.,False),  # 候选新建 0–2。
        Corridor('23',2,3,160.,False),  # 候选新建 2–3。
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
        existing = [e for e in self.corridors if e.existing]  # 只提取既有网架，不展开建设组合。
        zbase = self.voltage_kv**2/(self.transformer_kva/1000.)  # kV²/MVA 得到基准阻抗 Ω。
        line = self.lines[0]  # 既有线路均为 L 型。
        super().__init__(  # 与 33 节点使用相同的基础负荷、限值及运行接口。
            'four_bus_five_corridor',0,(1,2,3),[(e.start,e.end) for e in existing],
            [line.r_ohm_km*e.length_m/1000./zbase for e in existing],  # 电阻转标幺。
            [line.x_ohm_km*e.length_m/1000./zbase for e in existing],  # 电抗转标幺。
            self.transformer_kva,self.voltage_kv,np.zeros(3),np.zeros(3),(1,2,3),  # 三个独立负荷，无固定背景。
            np.full(3,np.tan(np.arccos(self.power_factor))),self.voltage_min_pu**2,1.,  # Q/P 和电压平方限值。
            self.transformer_kva*self.power_factor,capacity=line.capacity_kw/self.transformer_kva,  # 无损负荷外界及基础线路有功上限。
            source_smax=1.,sources=('Network/__init__.py','Network/four_bus_five_corridor.py'))  # 配变视在容量为 1 p.u.。

    @property  # 模型使用内部索引，原始数据保留真实节点号。
    def planning_branches(self):  # 所有五条走廊均可进入规划，而不是固定原来的三条树边。
        index = {self.root:-1,**{node:i for i,node in enumerate(self.nodes)}}  # 根索引 −1，其余按 nodes 排列。
        return tuple((index[e.start],index[e.end]) for e in self.corridors)  # 方向由主问题决定。

    @cached_property  # 每条走廊的三个型号仅生成一次。
    def line_options(self):  # 与 case33 共用逐线路输入；optional 允许开断或不建设。
        zbase = self.voltage_kv**2/(self.base/1000.)  # 所有型号使用同一标幺基准。
        return tuple(LineOptions(i,  # branch 对应 planning_branches 中的走廊位置。
            tuple(t.r_ohm_km*e.length_m/1000./zbase for t in self.lines),  # L/M/H 电阻。
            tuple(t.x_ohm_km*e.length_m/1000./zbase for t in self.lines),  # L/M/H 电抗。
            tuple(0. if e.existing and k==0 else e.length_m*(t.cable_cny_m+(0. if e.existing else self.new_corridor_cny_m)) for k,t in enumerate(self.lines)),  # 保留既有 L 型免费，其他选择计原造价。
            tuple(t.capacity_kw/self.base for t in self.lines),True)  # 型号有功上限；每条走廊可不投入。
            for i,e in enumerate(self.corridors))  # 只展开逐线路选项，不生成组合表。

    def design(self, choice):  # 将一次优化结果转换为独立 AC 可读取的固定径向网。
        selected = [(e,int(k)) for e,k in enumerate(choice) if k>=0]  # −1 表示走廊未投入，0/1/2 表示 L/M/H。
        options = self.line_options  # 费用、阻抗和容量均从唯一型号表读取。
        return RadialNetwork(self.name,self.root,self.nodes,  # 径向接口会检查连通性并从根重新定向。
            [(self.corridors[e].start,self.corridors[e].end) for e,k in selected],
            [options[e].r[k] for e,k in selected],[options[e].reactance[k] for e,k in selected],  # 当前选型的阻抗。
            self.base,self.voltage_kv,self.original_p,self.original_q,self.load_nodes,self.q_ratio,  # 同一基础负荷和坐标。
            self.vmin,self.vmax,self.power_limit,capacity=[options[e].capacity[k] for e,k in selected],  # 同一节点限值和逐型号容量。
            source_smax=self.source_smax,x=np.asarray(choice,dtype=int),  # 同一配变容量与本次逐线路选择。
            cost=sum(options[e].cost[k] for e,k in selected),sources=self.sources)  # 只记录这一方案的费用与数据来源。


network = FourBus()  # 保留原网架模块可直接导入的入口。
