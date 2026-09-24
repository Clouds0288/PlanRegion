"""江口完整候选图：建筑必须接入，非负荷中间节点可随选定路径投入。"""
import json
from pathlib import Path

import numpy as np

from . import Corridor, Network, TypeParameters


def load_jiangkou(mode="full"):
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


class Jiangkou(Network):
    """既有线路可重构、升级，新走廊按给定工程费用计价。"""
    cost_unit, budgets = "元", (0., np.inf)

    def __init__(self, load_nodes=None, mode='full'):
        data = self.data = load_jiangkou(mode)
        self.node_ids = {n["source_id"]: n["id"] for n in data["nodes"]}
        load_nodes = tuple(self.node_ids[n] for n in (data["default_load_nodes"] if load_nodes is None else load_nodes))
        rows = [n for n in data["nodes"] if n["id"] != data["root"]]
        nodes = tuple(n["id"] for n in rows)
        p = np.array([n["p_kw"] for n in rows])
        ratio = np.tan(np.arccos(data["power_factor"]))
        corridors = data["corridors"]
        base = data["transformer_kva"]
        power_limit = base*data["power_factor"]-p.sum()+sum(p[nodes.index(n)] for n in load_nodes)
        candidates = tuple(Corridor(e['id'], (e['start'], e['end']), e['line_type'], e['initial_active'],
            tuple(TypeParameters(o['line_type'], o['r'], o['reactance'], o['capacity'], o['cost'])
                  for o in e['options'])) for e in corridors)
        super().__init__('jiangkou_'+mode, data["root"], nodes, candidates,
            base, data["voltage_kv"], p, p*ratio, load_nodes, np.full(len(load_nodes), ratio),
            data["voltage_min_pu"]**2, data["voltage_max_pu"]**2, power_limit,
            required=[n['required'] for n in rows], source_smax=1.,
            sources=("Network/__init__.py", "Network/jiangkou.py", "Network/data/jiangkou.json"))


network = Jiangkou()
