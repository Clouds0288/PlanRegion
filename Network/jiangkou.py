"""江口数据：Jiangkou() 接入既有径向网；load_jiangkou('full') 读取完整候选图。

full 中非负荷中间节点允许不投入，供后续拓扑模型适配；本模块不运行求解器。
可用 Jiangkou(load_nodes=('B000025', 'B000078', 'B000083')) 指定独立负荷建筑。
"""
from functools import cached_property
import json
from pathlib import Path

import numpy as np

from . import Corridor, TypeParameters, RadialNetwork


def load_jiangkou(mode="existing"):
    """返回 existing/full 数据；负荷为 kW，options 的阻抗/容量为 p.u.，费用为元。"""
    kinds = {"existing": {"ele", "service_ele"},
             "full": {"ele", "road", "service_ele", "service_road"}}[mode]
    data = json.loads((Path(__file__).parent/"data"/"jiangkou.json").read_text(encoding="utf-8"))
    for name, columns in (("nodes", "node_columns"), ("corridors", "corridor_columns")):
        fields = data.pop(columns)
        data[name] = [dict(zip(fields, row, strict=True)) for row in data[name]]
    data["mode"] = mode
    data["corridors"] = [e for e in data["corridors"] if e["kind"] in kinds
                         and (mode == "full" or not e["normally_open"])]
    endpoints = {e[k] for e in data["corridors"] for k in ("start", "end")}
    data["nodes"] = [n for n in data["nodes"] if n["id"] in endpoints]
    for node in data["nodes"]:
        node["required"] = mode == "existing" or node["kind"] == "building" or node["id"] == data["root"]
    base, voltage, lines = data["transformer_kva"], data["voltage_kv"], data["line_types"]
    zbase = voltage**2/(base/1000.)
    capacity_factor = np.sqrt(3)*voltage*1000*data["power_factor"]*data["line_loading_limit"]/base
    for edge in data["corridors"]:
        original = edge["line_type"]
        edge["initial_active"] = original is not None and not edge["normally_open"]
        edge["must_use"] = mode == "existing"
        types = ([original] if original else []) + [t for t in data["planning_line_types"]
                if original is None or lines[t]["max_i_ka"] > lines[original]["max_i_ka"]]
        factor = 1. if original else data["new_corridor_cost_factor"]
        edge["options"] = [dict(line_type=t,
            r=lines[t]["r_ohm_per_km"]*edge["length_m"]/1000./zbase,
            reactance=lines[t]["x_ohm_per_km"]*edge["length_m"]/1000./zbase,
            capacity=lines[t]["max_i_ka"]*capacity_factor,
            cost=0. if t == original else edge["length_m"]*lines[t]["cost_yuan_per_m"]*factor)
            for t in types]
    return data


class Jiangkou(RadialNetwork):
    """只升级已有线路的网架；第 0 选项保持原状，其余选项不降低原线路容量。"""
    cost_unit, budgets = "元", (0., np.inf)

    def __init__(self, load_nodes=None):
        data = self.data = load_jiangkou()
        self.node_ids = {n["source_id"]: n["id"] for n in data["nodes"]}
        load_nodes = tuple(self.node_ids[n] for n in (data["default_load_nodes"] if load_nodes is None else load_nodes))
        rows = [n for n in data["nodes"] if n["id"] != data["root"]]
        nodes = tuple(n["id"] for n in rows)
        p = np.array([n["p_kw"] for n in rows])
        ratio = np.tan(np.arccos(data["power_factor"]))
        corridors = data["corridors"]
        baseline = [e["options"][0] for e in corridors]
        base = data["transformer_kva"]
        power_limit = base*data["power_factor"]-p.sum()+sum(p[nodes.index(n)] for n in load_nodes)
        super().__init__("jiangkou_existing", data["root"], nodes,
            [(e["start"], e["end"]) for e in corridors],
            [o["r"] for o in baseline], [o["reactance"] for o in baseline],
            base, data["voltage_kv"], p, p*ratio, load_nodes, np.full(len(load_nodes), ratio),
            data["voltage_min_pu"]**2, data["voltage_max_pu"]**2, power_limit,
            capacity=[o["capacity"] for o in baseline], source_smax=1.,
            sources=("Network/__init__.py", "Network/jiangkou.py", "Network/data/jiangkou.json"))
        by_endpoints = {frozenset((e["start"], e["end"])): e for e in corridors}
        self.corridor_data = tuple(by_endpoints[frozenset((self.root if a < 0 else self.nodes[a], n))]
                                  for n, a in zip(self.nodes, self.parent))

    @cached_property
    def corridors(self):
        return tuple(Corridor(e['id'],
            (self.root if a < 0 else self.nodes[a], n), e['line_type'],
            e['initial_active'], e['must_use'],
            tuple(TypeParameters(o['line_type'], o['r'], o['reactance'], o['capacity'], o['cost'])
                  for o in e['options']))
            for e, n, a in zip(self.corridor_data, self.nodes, self.parent))


network = Jiangkou()
