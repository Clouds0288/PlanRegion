"""统一型号索引、径向重构、固定方案和可选节点。"""
from dataclasses import fields, replace

import numpy as np
import pytest
from gurobipy import GRB

from Network import Corridor, Network, TypeParameters
from Network.case33bw import Case33
from Network.four_bus_five_corridor import FourBus
from model import GridPhysics, MasterProblem, SubProblem


def modified(network, **changes):
    if 'corridors' in changes and 'road_allowed' not in changes:
        allowed = {c.id: network.road_allowed[i] for i, c in enumerate(network.corridors)}
        changes['road_allowed'] = [allowed[c.id] for c in changes['corridors']]
    return Network(**({f.name: getattr(network, f.name) for f in fields(Network)} | changes))


@pytest.mark.parametrize('network', [FourBus(), Case33()])
def test_one_index_for_types_flows_and_costs(network):
    e = GridPhysics(network, 'socp')
    x = network.encode_plan(network.initial_plan)
    assert x.shape == network.r.shape == network.cost.shape == (network.n_types,)
    assert network.decode_plan(x) == network.initial_plan
    assert e.slack_slice.stop-e.slack_slice.start == 2*network.n_corridors
    problem = MasterProblem(e, fixed_plan=network.initial_plan, threads=1)
    with problem.model:
        problem.model.update()
        np.testing.assert_array_equal(problem.x.LB, x)
        np.testing.assert_array_equal(problem.x.UB, x)
        assert problem.x.shape == (network.n_types,)
        for corridor in network.corridors:
            row = problem.model.getRow(problem.model.getConstrByName(f'one_type[{corridor.id}]'))
            assert {row.getVar(i).VarName for i in range(row.size())} == {
                f'x[{corridor.id},{k.id}]' for k in corridor.types}
        balance_rows = [c for c in problem.model.getConstrs()
                        if c.ConstrName.startswith(('active_balance[', 'reactive_balance[', 'voltage_drop['))]
        assert len(balance_rows) == 2*network.n+network.n_corridors


def test_fixed_single_types_still_count_investment():
    net = FourBus()
    corridors = tuple(replace(c, types=(replace(c.types[0], investment_cost=7.),))
                      for c in net.corridors if c.initial_active)
    net = modified(net, corridors=corridors)
    for budget, feasible in ((20., False), (21., True)):
        problem = MasterProblem(GridPhysics(net, 'linear'), budget=budget, power=[3., 4., 5.], threads=1)
        with problem.model:
            answer = problem.solve()
        assert (answer is not None) == feasible
        if feasible:
            assert answer['objective'] == 21.
            np.testing.assert_array_equal(answer['x'], np.ones(3))


def optional_network():
    line = TypeParameters('line', .01, .005, 1., 0.)
    edges = ((0, 1), (1, 2), (0, 3), (3, 2), (3, 4))
    return Network('optional', 0, (1, 2, 3, 4),
        tuple(Corridor(str(i), ends, 'line', i < 2, (line,)) for i, ends in enumerate(edges)),
        100., .4, np.zeros(4), np.zeros(4), (1, 2), [.3, .3], .9**2, 1., 100.,
        required=[True, True, False, False], source_smax=1.)


@pytest.mark.parametrize('method', ['linear', 'socp'])
def test_optional_nodes_follow_selected_paths(method):
    net = optional_network()
    for active, node_count in (((0, 1), 2), ((0, 2, 3), 3)):
        plan = {c.id: 'line' if int(c.id) in active else None for c in net.corridors}
        x = net.encode_plan(plan)
        tree = net.tree(x)
        assert tree.n == node_count
        assert 4 not in tree.nodes
        e = GridPhysics(net, method)
        problem = MasterProblem(e, fixed_plan=plan, power=[10., 10.], threads=1)
        with problem.model:
            answer = problem.solve()
            np.testing.assert_array_equal(problem.active_nodes.X > .5, np.isin(net.nodes, tree.nodes))
        assert answer['feasible']
        assert SubProblem(e, threads=1).solve(x, np.array([10., 10.]))['feasible']
    missing = net.initial_plan | {'1': None}
    with pytest.raises(ValueError, match='required'):
        net.tree(net.encode_plan(missing))
    problem = MasterProblem(GridPhysics(net, method), fixed_plan=missing,
                            power=np.zeros(2), threads=1)
    with problem.model:
        assert problem.solve() is None


@pytest.mark.parametrize('method', ['linear', 'socp'])
def test_joint_cut_is_valid_with_optional_nodes(method):
    net = optional_network()
    e = GridPhysics(net, method)
    x, power = net.encode_plan(net.initial_plan), np.array([80., 80.])
    cut = SubProblem(e, threads=1).solve(x, power)['cut']
    assert cut is not None and cut[0]+cut[1:3]@power+cut[3:]@x < -1e-9
    problem = MasterProblem(e, threads=1)
    with problem.model:
        problem.model.setObjective(cut[0]+cut[1:3]@problem.power+cut[3:]@problem.x, GRB.MINIMIZE)
        problem.model.optimize()
        assert problem.model.Status == GRB.OPTIMAL
        assert problem.model.ObjBound >= -1e-7


def test_fingerprint_covers_physics_and_initial_status():
    net = FourBus()
    assert net.fingerprint == FourBus().fingerprint
    assert net.fingerprint != modified(net, vmin=.8).fingerprint
    changed = replace(net.corridors[0], initial_active=False)
    assert net.fingerprint != modified(net, corridors=(changed, *net.corridors[1:])).fingerprint


def test_invalid_plan_or_tree_fails_at_boundary():
    net = FourBus()
    with pytest.raises(ValueError, match='every corridor'):
        net.encode_plan({'01': 'L'})
    with pytest.raises(ValueError, match='unknown type'):
        net.encode_plan(net.initial_plan | {'01': 'missing'})
    cycle = net.initial_plan | {'02': 'L'}
    with pytest.raises(ValueError, match='rooted tree'):
        net.tree(net.encode_plan(cycle))
    x = net.encode_plan(net.initial_plan)
    x[1] = 1
    with pytest.raises(ValueError, match='one type'):
        net.tree(x)
