"""MATPOWER case33bw：32 条常闭线路和 5 条联络线均参与径向重构。"""
from pathlib import Path
import re

import numpy as np

from . import Corridor, Network, Project, TypeParameters


class Case33(Network):
    # 并联升级的合成投资，单位并非 MATPOWER 提供的工程造价。
    projects = (  # 前 4、8、16、32 项形成嵌套候选集；费用与线路只在此定义。
        Project("A", (2, 3), 2.),     # 原候选：公共上游线路。
        Project("B", (12, 13), 1.),   # 原候选：节点 18 所在支路。
        Project("C", (23, 24), 1.),   # 原候选：节点 25 所在支路。
        Project("D", (27, 28), 1.),   # 原候选：节点 33 所在支路。
        Project("E", (5, 6), 1.),     # 新候选：下游分支的上游线路。
        Project("F", (13, 14), 1.),   # 新候选：进一步改善节点 18 方向。
        Project("G", (24, 25), 1.),   # 新候选：节点 25 的入线。
        Project("H", (28, 29), 1.),   # 新候选：进一步改善节点 33 方向。
        Project("I", (3, 4), 1.),     # 十六线路新增：节点 18、33 方向的公共上游。
        Project("J", (6, 7), 1.),     # 十六线路新增：节点 18 所在分支的入口。
        Project("K", (9, 10), 1.),    # 十六线路新增：节点 18 方向的中段。
        Project("L", (16, 17), 1.),   # 十六线路新增：节点 18 方向的末段。
        Project("M", (3, 23), 1.),    # 十六线路新增：节点 25 所在分支的入口。
        Project("N", (6, 26), 1.),    # 十六线路新增：节点 33 所在分支的入口。
        Project("O", (29, 30), 1.),   # 十六线路新增：节点 33 方向的中下游。
        Project("P", (31, 32), 1.),   # 十六线路新增：节点 33 方向的末段。
        Project("Q", (1, 2), 1.),
        Project("R", (4, 5), 1.),
        Project("S", (7, 8), 1.),
        Project("T", (8, 9), 1.),
        Project("U", (10, 11), 1.),
        Project("V", (11, 12), 1.),
        Project("W", (14, 15), 1.),
        Project("X", (15, 16), 1.),
        Project("Y", (17, 18), 1.),
        Project("Z", (2, 19), 1.),
        Project("AA", (19, 20), 1.),
        Project("AB", (20, 21), 1.),
        Project("AC", (21, 22), 1.),
        Project("AD", (26, 27), 1.),
        Project("AE", (30, 31), 1.),
        Project("AF", (32, 33), 1.),
    )

    def __init__(self, load_nodes=(18, 25, 33), upgrade_count=4):
        if not 0 <= upgrade_count <= len(type(self).projects):
            raise ValueError(f'upgrade_count must be between 0 and {len(type(self).projects)}')
        self.upgrade_count = upgrade_count
        self.projects = type(self).projects[:upgrade_count]
        source = (Path(__file__).parent/'data'/'case33bw.m').read_text(encoding='utf-8')
        tables = {}
        for name in ('bus', 'gen', 'branch'):
            block = re.search(rf"mpc\.{name}\s*=\s*\[(.*?)\];", source, re.S)[1]
            block = re.sub(r"%[^\n]*", "", block)
            tables[name] = np.array([[float(x) for x in row.split()]
                                     for row in block.split(';') if row.strip()])
        bus, gen, branch = (tables[k] for k in ('bus', 'gen', 'branch'))
        base_mva = float(re.search(r"mpc.baseMVA\s*=\s*([\d.]+)", source)[1])
        root = int(bus[bus[:, 1] == 3, 0][0])
        rows = bus[bus[:, 0] != root]
        nodes = tuple(map(int, rows[:, 0]))
        selected = [nodes.index(i) for i in load_nodes]
        zbase = bus[0, 9]**2/base_mva
        projects = {frozenset(p.branch): p for p in self.projects}
        corridors = []
        for row in branch:
            a, b = map(int, row[:2])
            # rateA=0 表示未提供线路限额；baseMVA 也不是变压器容量。
            original = TypeParameters('existing', row[2]/zbase, row[3]/zbase,
                                      row[5]/base_mva if row[5] else np.inf, 0.)
            project = projects.get(frozenset((a, b)))
            types = (original,) if project is None else (original,
                TypeParameters('parallel', original.r/2, original.reactance/2, original.capacity, project.cost))
            corridors.append(Corridor(f'{a}-{b}', (a, b), 'existing', bool(row[10]), types))
        power_limit = gen[0, 8]*1000-rows[:, 2].sum()+rows[selected, 2].sum()
        super().__init__('case33bw', root, nodes, tuple(corridors),
            base_mva*1000, bus[0, 9], rows[:, 2], rows[:, 3], tuple(load_nodes),
            rows[selected, 3]/rows[selected, 2], rows[:, 12]**2, rows[:, 11]**2, power_limit,
            source_pmax=gen[0, 8]/base_mva, source_qmax=gen[0, 3]/base_mva,
            sources=('Network/__init__.py', 'Network/case33bw.py', 'Network/data/case33bw.m'))


network = Case33()
