"""MATPOWER case33bw 的算例，原始线路参数不升级，保持全网径向连通。

Case33：限定七条开关，预算为相对原始状态的开合次数。
Case33Plan：扩展规划，S1–S5 为可开断的既有线路（基态闭合、开断不计费），C1–C5 为基态不建的候选走廊
（即原五条联络线，阻抗不变），预算为所建候选的相对建设费之和。
Case33S：勘察算例，网架与费用同 Case33Plan，全部线路 250 A，按信息状态屏蔽不可用的候选路。
"""
from pathlib import Path
import re

import numpy as np

from . import Corridor, Network, TypeParameters


LOAD_NODES = (18, 25)             # 二维：(18, 25)；三维：(18, 25, 30)
CURRENT_LIMIT = 200.             # A；主线、SOCP 和 AC 共用，可改为逐线路数组


class _Case33bw(Network):
    """按 case33bw.m 建网：price(a, b, active) 给可变线路的费用系数，None 为固定线路。"""

    def _build(self, name, load_nodes, current_limit, price):
        # 1. MATPOWER 的 bus / gen / branch 表与 baseMVA
        source = (Path(__file__).parent/'data'/'case33bw.m').read_text(encoding='utf-8')
        tables = {}
        for table in ('bus', 'gen', 'branch'):
            block = re.search(rf"mpc\.{table}\s*=\s*\[(.*?)\];", source, re.S)[1]
            block = re.sub(r"%[^\n]*", "", block)
            tables[table] = np.array([[float(x) for x in row.split()] for row in block.split(';') if row.strip()])
        bus, gen, branch = (tables[k] for k in ('bus', 'gen', 'branch'))
        base_mva = float(re.search(r"mpc.baseMVA\s*=\s*([\d.]+)", source)[1])
        root = int(bus[bus[:, 1] == 3, 0][0])
        rows = bus[bus[:, 0] != root]
        nodes = tuple(map(int, rows[:, 0]))
        selected = [nodes.index(i) for i in load_nodes]
        zbase = bus[0, 9]**2/base_mva
        # 2. 输入为线路电流有效值 A；额定上限独立于原 rateA 和派生电流界。
        current_limit = np.broadcast_to(current_limit, (len(branch),)).astype(float)
        if np.any(np.isnan(current_limit)) or np.any(current_limit <= 0.):
            raise ValueError('current_limit must be positive amperes or infinity')
        ell_limit = (current_limit/(base_mva*1000/(np.sqrt(3)*bus[0, 9])))**2
        # 3. 走廊：rateA=0 表示未提供线路限额；baseMVA 也不是变压器容量。
        corridors = []
        for row, limit in zip(branch, ell_limit):
            a, b = map(int, row[:2])
            cost = price(a, b, bool(row[10]))
            original = TypeParameters('existing', row[2]/zbase, row[3]/zbase,
                                      row[5]/base_mva if row[5] else np.inf, 0. if cost is None else cost, limit)
            corridors.append(Corridor(f'{a}-{b}', (a, b), 'existing', bool(row[10]), (original,),
                                      switchable=cost is not None))
        # 4. 预算表达式为 c0+c@x，c0=-c@x0 使原始方案的费用为 0。
        power_limit = gen[0, 8]*1000-rows[:, 2].sum()+rows[selected, 2].sum()
        super().__init__(name, root, nodes, tuple(corridors),
            base_mva*1000, bus[0, 9], rows[:, 2], rows[:, 3], tuple(load_nodes),
            rows[selected, 3]/rows[selected, 2], rows[:, 12]**2, rows[:, 11]**2, power_limit,
            cost_offset=-sum(c.types[0].investment_cost for c in corridors if c.initial_active and c.switchable),
            source_pmax=gen[0, 8]/base_mva, source_qmax=gen[0, 3]/base_mva,
            sources=('Network/__init__.py', 'Network/case33bw.py', 'Network/data/case33bw.m'))


class Case33(_Case33bw):
    """唯一预算入口：相对原始状态，断开或闭合一条线路各计一次（常闭系数 -1，常开系数 +1）。"""
    switch_budget = 7
    cost_unit = '次开合变动'
    budgets = (switch_budget,)
    switchable_branches = ((21, 8), (7, 8), (22, 12), (11, 12),
                           (9, 15), (33, 18), (25, 29))

    def __init__(self, load_nodes=LOAD_NODES, *, current_limit=np.inf):
        switchable = {frozenset(edge) for edge in self.switchable_branches}
        self._build('case33bw', load_nodes, current_limit,
                    lambda a, b, active: (1-2*int(active)) if frozenset((a, b)) in switchable else None)


class Case33Plan(_Case33bw):
    """扩展规划：基态与原 Case33 一致（S1–S5 闭合，C1–C5 不建）；开断 S 不计费，建设 C 计相对建设费。"""
    switches = ((7, 8), (11, 12), (14, 15), (28, 29), (32, 33))                          # S1–S5
    candidates = {(8, 21): 4., (9, 15): 4., (12, 22): 4., (18, 33): 1., (25, 29): 1.}   # C1–C5：端点 → 相对建设费
    plan_budget = 14.   # 全部候选都可建：共 87 个径向方案
    cost_unit = '相对建设费'
    budgets = (plan_budget,)

    def __init__(self, load_nodes=LOAD_NODES, *, current_limit=np.inf):
        prices = {**{frozenset(edge): 0 for edge in self.switches},
                  **{frozenset(edge): cost for edge, cost in self.candidates.items()}}
        self._build('case33bw_plan', load_nodes, current_limit, lambda a, b, active: prices.get(frozenset((a, b))))


class Case33S(Case33Plan):
    """勘察算例 Case33-S：网架与费用同 Case33Plan（C1–C5 为候选路，1 单位为 C4 造价），全部线路 250 A，电压 0.9–1.1 p.u.
    （原数据）。available 为可用候选路的序号（信息状态下的道路掩码），其余候选路的 road_allowed 为 False。"""
    current_limit = 250.        # A，全部线路
    status_quo = (90., 420.)    # 现状点 z^0（kW）：端口 18、25 的原负荷

    def __init__(self, load_nodes=LOAD_NODES, *, available=range(5)):
        super().__init__(load_nodes, current_limit=self.current_limit)
        self.name = 'case33bw_s'
        ids = {frozenset(c.endpoints): c.id for c in self.corridors}
        self.roads = tuple(ids[frozenset(edge)] for edge in self.candidates)   # C1–C5 的走廊 ID
        self.road_allowed = np.array([c.id not in self.roads or self.roads.index(c.id) in available
                                      for c in self.corridors])


network = Case33()
