"""五节点六候选道路概念案例；参数与 2026-09-25 已冻结实验一致。"""
import numpy as np

from . import Corridor, Network, TypeParameters


EXISTING = [('03', (0, 3), 85.), ('31', (3, 1), 55.),
            ('04', (0, 4), 65.), ('42', (4, 2), 55.)]
ROUTES = {
    'A': {'corridors': ['A'], 'endpoints': [0, 1], 'capacity_kw': 90., 'investment_cost': 2.5},
    'B': {'corridors': ['B'], 'endpoints': [0, 2], 'capacity_kw': 90., 'investment_cost': 2.5},
    'C': {'corridors': ['C'], 'endpoints': [1, 2], 'capacity_kw': 90., 'investment_cost': 1.},
    'D': {'corridors': ['D'], 'endpoints': [1, 4], 'capacity_kw': 70., 'investment_cost': 1.},
    'E': {'corridors': ['E'], 'endpoints': [2, 3], 'capacity_kw': 90., 'investment_cost': 2.},
    'F': {'corridors': ['F'], 'endpoints': [3, 4], 'capacity_kw': 40., 'investment_cost': 1.},
}
for data in ROUTES.values():
    data.update(r_ohm=.6, reactance_ohm=.3, survey_cost=1., road_usable_probability=.8)


class Concept5(Network):
    """所有节点必接入；既有线路可开断，候选道路可用不等于必须建设。"""
    budgets = (4.,)
    cost_unit = '相对投资单位'

    def __init__(self, allowed_roads=None):
        allowed_roads = set(ROUTES) if allowed_roads is None else set(allowed_roads)
        if not allowed_roads <= set(ROUTES):
            raise ValueError('unknown candidate road')
        base, voltage_kv = 100., 11.
        zbase = voltage_kv**2/(base/1000.)
        corridors = []
        for name, endpoints, capacity_kw in EXISTING:
            parameters = TypeParameters('existing', .6/zbase, .3/zbase, capacity_kw/base, 0.)
            corridors.append(Corridor(name, endpoints, 'existing', True, (parameters,)))
        for route, data in ROUTES.items():
            parameters = TypeParameters('new', data['r_ohm']/zbase, data['reactance_ohm']/zbase,
                                        data['capacity_kw']/base, data['investment_cost'])
            corridors.append(Corridor(route, tuple(data['endpoints']), None, False, (parameters,)))
        super().__init__('concept_five_bus', 0, (1, 2, 3, 4), tuple(corridors), base, voltage_kv,
                         np.array([40., 40., 0., 0.]), np.array([12., 12., 0., 0.]), (1, 2), [.3, .3],
                         .95**2, 1.05**2, 135., source_pmax=1.35, source_qmax=1.,
                         sources=('Network/__init__.py', 'Network/concept5.py'),
                         road_allowed=[c.existing_type is not None or c.id in allowed_roads for c in corridors])


network = Concept5()
