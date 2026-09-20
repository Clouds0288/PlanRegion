"""MATPOWER case33bw 原始径向网架；三个实际节点的负荷独立变化。"""
from pathlib import Path
from copy import copy
from functools import cached_property
from itertools import product
import re

import numpy as np

from . import LineOptions, Project, RadialNetwork


class Case33(RadialNetwork):
    # 合成投资试验：费用为相对单位，不是 MATPOWER 提供的工程造价。
    projects = (Project("A", (2, 3), 2.), Project("B", (12, 13), 1.),
                Project("C", (23, 24), 1.), Project("D", (27, 28), 1.))

    def __init__(self, load_nodes=(18, 25, 33)):
        path = Path(__file__).parent/"data"/"case33bw.m"
        source = path.read_text(encoding="utf-8")
        tables = {}
        for name in ("bus", "gen", "branch"):
            block = re.search(rf"mpc\.{name}\s*=\s*\[(.*?)\];", source, re.S)[1]
            block = re.sub(r"%[^\n]*", "", block)
            tables[name] = np.array([[float(x) for x in row.split()]
                                     for row in block.split(";") if row.strip()])
        bus, gen, branch = (tables[k] for k in ("bus", "gen", "branch"))
        base_mva = float(re.search(r"mpc.baseMVA\s*=\s*([\d.]+)", source)[1])
        root = int(bus[bus[:, 1] == 3, 0][0])
        rows = bus[bus[:, 0] != root]
        nodes = tuple(map(int, rows[:, 0]))
        selected = [nodes.index(i) for i in load_nodes]  # 独立变化的真实负荷节点对应的数组位置。
        active = branch[branch[:, 10] == 1]  # 五条常开联络线保持断开。
        zbase = bus[0, 9]**2/base_mva  # 原始线路阻抗单位是 Ω，用 Ubase²/Sbase 换成标幺值。
        power_limit = gen[0, 8]*1000-rows[:, 2].sum()+rows[selected, 2].sum()  # 电源有功上限减固定背景负荷，得到独立节点无损总量外界。
        # rateA=0：原始数据未启用线路热限；baseMVA 不是配变容量。
        super().__init__(
            "case33bw", root, nodes, active[:, :2].astype(int), active[:, 2]/zbase,
            active[:, 3]/zbase, base_mva*1000, bus[0, 9], rows[:, 2], rows[:, 3],
            tuple(load_nodes), rows[selected, 3]/rows[selected, 2], rows[:, 12]**2,
            rows[:, 11]**2, power_limit, source_pmax=gen[0, 8]/base_mva,  # 电压幅值限值先平方；Pmax 从 MW 转标幺。
            source_qmax=gen[0, 3]/base_mva,  # Qmax 从 Mvar 转标幺，保留原始电源运行限制。
            sources=("Network/__init__.py", "Network/case33bw.py", "Network/data/case33bw.m"))

    @cached_property
    def line_options(self):
        """紧凑模型只读取逐线路选项，不生成完整建设组合。"""
        options = []
        for project in self.projects:
            e = self.nodes.index(project.branch[1])
            options.append(LineOptions(e, (self.r[e], self.r[e]/2),
                                       (self.reactance[e], self.reactance[e]/2),
                                       (0., project.cost)))  # 同走廊并联一回，R/X 减半；费用仍由 projects 唯一定义。
        return tuple(options)

    def design(self, choice):
        """按一次求解给出的线路型号实例化网架；不遍历其他组合。"""
        design = copy(self)
        design.r, design.reactance = self.r.copy(), self.reactance.copy()
        for k, options in zip(choice, self.line_options):
            design.r[options.branch] = options.r[k]
            design.reactance[options.branch] = options.reactance[k]
        design.x = np.asarray(choice, dtype=int)  # 四条线路分别取 0/1 型，与历史建设向量一致。
        design.cost = sum(o.cost[k] for o, k in zip(self.line_options, choice))
        return design

    @cached_property
    def designs(self):
        """仅供小算例枚举对照；紧凑 MILP/MISOCP 与联合割不读取此属性。"""
        designs = [self.design(choice) for choice in product(
            *(range(len(o.cost)) for o in self.line_options))]
        return tuple(sorted(designs, key=lambda d: (d.cost, tuple(d.x))))


network = Case33()
