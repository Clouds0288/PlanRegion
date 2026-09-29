"""MATPOWER case33bw：限定七条开关，原始线路参数不升级，保持全网径向连通。"""
from pathlib import Path
import re

import numpy as np

from . import Corridor, Network, TypeParameters


LOAD_NODES = (18, 25)             # 二维：(18, 25)；三维：(18, 25, 30)


class Case33(Network):
    # 1. 唯一预算入口：相对原始状态，断开或闭合一条线路各计一次。
    switch_budget = 7
    cost_unit = '次开合变动'
    budgets = (switch_budget,)
    switchable_branches = ((21, 8), (7, 8), (22, 12), (11, 12),
                           (9, 15), (33, 18), (25, 29))

    def __init__(self, load_nodes=LOAD_NODES):
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
        switchable = {frozenset(edge) for edge in self.switchable_branches}
        corridors = []
        for row in branch:
            a, b = map(int, row[:2])
            allowed = frozenset((a, b)) in switchable
            active = bool(row[10])
            # 2. 预算表达式为 c0+c@x：常闭系数 -1，常开系数 +1。
            cost = (1-2*int(active)) if allowed else 0.
            # rateA=0 表示未提供线路限额；baseMVA 也不是变压器容量。
            original = TypeParameters('existing', row[2]/zbase, row[3]/zbase,
                                      row[5]/base_mva if row[5] else np.inf, cost)
            corridors.append(Corridor(f'{a}-{b}', (a, b), 'existing', active,
                                      (original,), switchable=allowed))
        power_limit = gen[0, 8]*1000-rows[:, 2].sum()+rows[selected, 2].sum()
        super().__init__('case33bw', root, nodes, tuple(corridors),
            base_mva*1000, bus[0, 9], rows[:, 2], rows[:, 3], tuple(load_nodes),
            rows[selected, 3]/rows[selected, 2], rows[:, 12]**2, rows[:, 11]**2, power_limit,
            cost_offset=sum(c.initial_active and c.switchable for c in corridors),
            source_pmax=gen[0, 8]/base_mva, source_qmax=gen[0, 3]/base_mva,
            sources=('Network/__init__.py', 'Network/case33bw.py', 'Network/data/case33bw.m'))


network = Case33()
