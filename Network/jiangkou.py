"""江口低压台区（0.4 kV、1000 kVA）现状网架的算例：只换导线，端口的根→端口路径上的线段可升级。

Jiangkou：全网 342 段径向树。端口路径上的走廊可换为载流更大的规划型号（原型号免费，造价按长度×单价，万元），
按段整段换线：路径并集上原型号相同、中间不分叉的连续走廊为一段（Network.sections），共用一个决策；其余线路固定；
预算为所换导线的造价之和。负荷为 Network/data/jiangkou.json 的基准负荷（原始楼宇峰值 ×0.6）。
reduced=True 为等值网：只保留端口路径（主干），路径外的子树并入挂接点——负荷取基态（端口为零、现状导线）下
子树首段的送端功率（含子树网损），挂接点的电压下限抬高子树内相对挂接点的最大压降；挂在根上的子树保留首段并入其首节点，
零负荷子树删去。
"""
import json
from pathlib import Path

import numpy as np

from . import Corridor, Network, TypeParameters


LOAD_NODES = ('B000078', 'B000042')   # 二维：左下、右上；三维再加左上 B000011
BUDGET = 13.                          # 万元；二维测试预算：够换右上整条路径（11.1）或左下大部分（15.4），换不了两条
YUAN = 1e4                            # 元 / 万元


def _flow(parent, r, reactance, p, q):
    """径向树的精确潮流（DistFlow 不动点）：parent 为父节点下标（根为 -1，父先于子），返回各入边送端 P/Q 与节点电压平方，标幺。"""
    n = len(parent)
    below = np.eye(n)   # below[i,j]=1 表示节点 j 在入边 i 的下游
    for i in reversed(range(n)):
        if parent[i] >= 0:
            below[parent[i]] += below[i]
    ell, v = np.zeros(n), np.ones(n)
    for _ in range(100):
        P, Q = below@(p+r*ell), below@(q+reactance*ell)
        sending = np.where(parent >= 0, v[parent], 1.)
        ell = (P**2+Q**2)/sending
        v = 1.-below.T@(2*(r*P+reactance*Q)-(r**2+reactance**2)*ell)
    return P, Q, v, below


class Jiangkou(Network):
    """load_nodes 为楼宇编号；reduced 为按端口路径构造的等值网。"""
    cost_unit = '万元'
    budgets = (BUDGET, np.inf)

    def __init__(self, load_nodes=LOAD_NODES, *, reduced=False):
        # 1. 数据：根向树（走廊为 (父, 子)，父先于子）、型号表与基值
        data = json.loads((Path(__file__).parent/'data'/'jiangkou.json').read_text(encoding='utf-8'))
        rows = [dict(zip(data['node_columns'], row)) for row in data['nodes']]
        edges = [dict(zip(data['corridor_columns'], row)) for row in data['corridors']]
        lines, base, voltage_kv, root = data['line_types'], data['transformer_kva'], data['voltage_kv'], data['root']
        zbase, ibase = voltage_kv**2/(base/1000.), base/(np.sqrt(3)*voltage_kv)
        ratio = np.tan(np.arccos(data['power_factor']))
        nodes = tuple(row['id'] for row in rows[1:])
        index = {node: i for i, node in enumerate(nodes)}
        index[root] = -1
        incoming = {e['child']: e for e in edges}
        parent = np.array([index[incoming[node]['parent']] for node in nodes])
        p = np.array([row['p_kw'] for row in rows[1:]])
        q = p*ratio
        vmin = np.full(len(nodes), data['voltage_min_pu']**2)

        def option(line, length, cost):
            t = lines[line]
            return TypeParameters(line, t['r_ohm_per_km']*length/1000./zbase, t['x_ohm_per_km']*length/1000./zbase,
                                  np.inf, cost, (t['max_i_ka']*1000.*data['line_loading_limit']/ibase)**2)

        # 2. 端口路径：根→各端口的走廊可升级，其余走廊固定为原型号
        path = set()
        for node in load_nodes:
            while node != root:
                path.add(node)
                node = incoming[node]['parent']
        corridors = {}
        for node in nodes:
            e = incoming[node]
            upgrades = [t for t in data['planning_line_types'] if node in path
                        and lines[t]['max_i_ka'] > lines[e['line_type']]['max_i_ka']]
            corridors[node] = Corridor(e['id'], (e['parent'], node), e['line_type'], True,
                (option(e['line_type'], e['length_m'], 0.),
                 *(option(t, e['length_m'], e['length_m']*lines[t]['cost_yuan_per_m']/YUAN) for t in upgrades)),
                switchable=bool(upgrades))
        # 段：路径并集上原型号相同、中间不分叉的连续走廊共用一个升级决策（整段换线）
        section_of, sections = {}, []
        for node in nodes:
            above = incoming[node]['parent']
            if node not in path:
                continue
            if (above != root and incoming[above]['line_type'] == incoming[node]['line_type']
                    and sum(incoming[other]['parent'] == above for other in path) == 1):
                section_of[node] = section_of[above]
                sections[section_of[node]].append(corridors[node].id)
            else:
                section_of[node] = len(sections)
                sections.append([corridors[node].id])

        # 3. 等值网：基态潮流给出路径外子树的首段送端功率与相对挂接点的最大压降
        kept = list(nodes)
        if reduced:
            existing = [corridors[node].types[0] for node in nodes]
            fixed = np.where(np.isin(nodes, load_nodes), 0., p)/base
            P, Q, v, below = _flow(parent, np.array([t.r for t in existing]), np.array([t.reactance for t in existing]),
                                   fixed, fixed*ratio)
            drop = {i: (v[parent[i]] if parent[i] >= 0 else 1.)-v[below[i] > 0].min() for i in range(len(nodes))}
            kept = [node for node in nodes if node in path]
            for i, node in enumerate(nodes):
                a = parent[i]
                if node in path or (a >= 0 and nodes[a] not in path):
                    continue
                if a >= 0:
                    # 3a. 挂在主干节点上的子树整体并入挂接点
                    p[a] += P[i]*base
                    q[a] += Q[i]*base
                    vmin[a] = max(vmin[a], data['voltage_min_pu']**2+drop[i])
                elif fixed[below[i] > 0].sum() > 0.:
                    # 3b. 挂在根上的子树保留首段，其余并入首节点
                    children = [j for j in range(len(nodes)) if parent[j] == i]
                    p[i] += sum(P[j] for j in children)*base
                    q[i] += sum(Q[j] for j in children)*base
                    vmin[i] = max(vmin[i], data['voltage_min_pu']**2+v[i]-v[below[i] > 0].min())
                    kept.append(node)
        keep = [index[node] for node in kept]

        # 4. 网架：源端 1000 kVA；载流上限折为 ell_limit，送端有功不另设上限
        power_limit = base*data['power_factor']-p[keep].sum()+sum(p[index[node]] for node in load_nodes)
        super().__init__('jiangkou_reduced' if reduced else 'jiangkou', root, tuple(kept),
            tuple(corridors[node] for node in kept), base, voltage_kv, p[keep], q[keep], tuple(load_nodes),
            np.full(len(load_nodes), ratio), vmin[keep], data['voltage_max_pu']**2, power_limit, source_smax=1.,
            sources=('Network/__init__.py', 'Network/jiangkou.py', 'Network/data/jiangkou.json'))
        self.sections = tuple(map(tuple, sections))


network = Jiangkou()
