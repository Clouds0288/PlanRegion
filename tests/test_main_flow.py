"""RB 主线：锥几何与夹逼测度、细分规则、面判据几何、网架支撑认证与一轮 B，以及一个分区的端到端构域。"""
from collections import Counter
from types import SimpleNamespace

import numpy as np
import pytest
from threadpoolctl import threadpool_limits

import main
from monitor import RunMonitor
from Network.case33bw import Case33, CURRENT_LIMIT
from region import (CERTIFIED_FACES, Cone, Network, Radial, Support, allowed_offsets, boundary_faces, box_vertices,
                    build_partition, chebyshev_center, choose_criterion, clip_box, cone_outer, contains, face_measures,
                    halfspaces, piece, polytope_volume, sandwich)


def settings(d):
    return SimpleNamespace(threads=1, obbt_workers=1, tau=main.REGION_TAU, discovery_eps=main.DISCOVERY_EPS,
                           discovery_share=main.DISCOVERY_SHARE, mip_seconds=main.MIP_SECONDS, mip_gap=main.MIP_GAP,
                           min_width=main.MIN_WIDTH[d], max_cones=main.MAX_CONES[d], network_eps=d*main.REGION_TAU/2.)


def radial(code, seconds=120.):
    network = Case33(load_nodes=(18, 25, 30)[:len(code)], current_limit=CURRENT_LIMIT)
    monitor = RunMonitor()
    monitor.time_limit = seconds
    part = Radial(network, tuple(1 if c == 'p' else -1 for c in code), Case33.switch_budget, monitor,
                  settings(len(code)))
    part.deadline = monitor.clock()+seconds
    return part


def support(normal, ub, point):
    normal = np.asarray(normal, float)
    return dict(lb=-np.inf if point is None else float(normal@point), ub=ub, point=point, status=2, normal=normal)


def test_cone_outer_and_sandwich_clip_network_sets_to_the_outer_bound():
    for d, cube_volume, outer_volume in ((2, .64, .18), (3, .512, .036)):
        U, V = np.eye(d), .5*np.eye(d)   # 内三角顶点 0.5·e_j；c=2·1，外界 c@xi<=1.2 即腿长 0.6 的单纯形
        halfspace = (np.full(d, 2.), 1.2)
        np.testing.assert_allclose(sorted(np.linalg.norm(cone_outer(U, halfspace), axis=1)), [0.]+[.6]*d)
        cube = box_vertices(d)*.8
        assert polytope_volume(cube) == pytest.approx(cube_volume)
        assert polytope_volume(piece(cube, U, halfspace)) == pytest.approx(outer_volume)
        inner, outer = sandwich([(U, V, halfspace)], [('cube', cube)], d, {})
        assert inner == pytest.approx(outer_volume) and outer == pytest.approx(outer_volume)
        small = box_vertices(d)*.1
        inner, _ = sandwich([(U, V, halfspace)], [('small', small)], d, {})
        assert inner == pytest.approx(polytope_volume(np.vstack([np.zeros(d), V])))
        assert len(piece(box_vertices(d)+2., U, halfspace)) == 0
        np.testing.assert_allclose(clip_box(cone_outer(U, None)).max(axis=0), np.ones(d))


def test_split_options_star_edge_and_longest_edge_midpoint():
    for code in ('nn', 'nnn'):
        part = radial(code)
        d = part.d
        cone = Cone(0, tuple(range(d)), (), .5*np.eye(d))
        cone.final = dict(point=np.full(d, 1./d)*part.bounds, x=None)
        star, fallback = part.options(cone)
        assert len(star) == d and all(d in keys for keys in star)                 # 内部解点：星形剖分
        assert all(len(set(keys)) == d for keys in [*star, *fallback])
        assert len(fallback) == 2 and part.midpoints == {(0, 1): d+1}            # 最长棱（并列取首条）中点二分
        np.testing.assert_allclose(part.directions[d+1][:2], [2**-.5]*2)
        cone.final = dict(point=None, x=None)
        assert len(part.options(cone)) == 1
    near_edge = np.array([.48, .48, .04])*part.bounds   # 三维：一个锥坐标偏小，在对边上二分
    cone.final = dict(point=near_edge, x=None)
    edge, _ = part.options(cone)
    assert len(edge) == 2 and {keys[2] for keys in edge} == {2}


def test_face_criteria_boundary_faces_measures_and_center():
    corner = np.array([[0., 0.], [.5, 0.], [.5, .5], [0., .5]])
    faces = halfspaces(corner)
    boundary = boundary_faces(faces)
    assert boundary.sum() == 2                                      # xi_1>=0、xi_2>=0 由分区盒认证
    assert choose_criterion(faces, boundary) == 'origin'            # 原点在 P_x 上且离两个非边界面 0.5
    np.testing.assert_allclose(sorted(allowed_offsets(faces, 'origin', None, .005)[~boundary]), [.5/.995]*2)
    np.testing.assert_allclose(face_measures(corner, faces), .5)
    inside = halfspaces(corner+.2)
    center = chebyshev_center(inside)
    np.testing.assert_allclose(center, [.45, .45])
    assert choose_criterion(inside, boundary_faces(inside)) == 'center'
    upper = inside[np.argmax(inside[:, 0]), :]                      # 面 xi_1<=0.7：β'=β+tau·(β-n@c)
    assert allowed_offsets(upper[None], 'center', center, .005)[0] == pytest.approx(.7+.005*.25)


def test_network_certifies_faces_clips_its_outer_and_grows_with_violating_points():
    corner = np.array([[0., 0.], [.5, 0.], [.5, .5], [0., .5]])
    network = Network((0,), 'x', corner, 1., 2, main.REGION_TAU, .01)
    assert list(network.face_status).count('BOUNDARY') == 2 and list(network.face_status).count('PENDING') == 2
    for _ in range(2):   # 两个非边界面的支撑上界恰为 β：面认证，O_x 被裁到 P_x
        normal = network.next_normal()
        network.add(normal, support(normal, .5, .5*np.abs(normal)))
    assert network.status == 'certified' and np.isin(network.face_status, CERTIFIED_FACES).all()
    assert polytope_volume(network.outer) == pytest.approx(.25) and network.ratio() == pytest.approx(0.)
    grown = Network((1,), 'y', corner, 1., 2, main.REGION_TAU, .01)
    grown.add(np.array([0., 1.]), support([0., 1.], .9, np.array([.1, .8])))   # 审计点越过允许值：补进 V_x
    assert grown.version == 2 and any(np.allclose(v, [.1, .8]) for v in grown.vertices)
    assert grown.outer[:, 1].max() <= .9 and grown.status == 'active'
    empty = Network((2,), 'z', np.empty((0, 2)), 1., 2, main.REGION_TAU, .01)
    np.testing.assert_allclose(empty.next_normal(), [1., 0.])     # 尚无认证点：先查种子方向
    empty.add(np.array([1., 0.]), support([1., 0.], .3, np.array([.3, .1])))
    assert len(empty.vertices) == 1 and not empty.full


def test_support_pass_keeps_certified_points_inside_the_support_outer():
    part = radial('nn')
    x = tuple(int(v) for v in part.network.encode_plan(part.network.initial_plan))   # 原点可行：射线顶点作种子
    builder = Support(part)
    with threadpool_limits(limits=1):
        assert part.origin_feasible(x) and all(part.vertex(x, key) is not None for key in range(2))
        builder.run([x])
    network = builder.networks[x]
    assert network.status in ('certified', 'eps_B') and builder.count > 0 and x in part.equations.boxes
    assert network.full and contains(network.vertices, halfspaces(network.outer), 1e-7).all()   # P_x ⊆ O_x
    assert [key for key, _ in builder.sets()] == [(x, network.version)]


def test_one_partition_certifies_with_support_networks_and_records_frames():
    network = Case33(load_nodes=(18, 25), current_limit=CURRENT_LIMIT)
    monitor = RunMonitor()
    monitor.time_limit = 120.
    with threadpool_limits(limits=1):
        result = build_partition(network, (1, 1), Case33.switch_budget, monitor, settings(2))
    assert result['status'] == 'certified' and result['certified']
    assert result['cones'] >= 1 and len(result['outer']) == result['cones']
    assert result['gap'] <= 2*main.REGION_TAU
    events = Counter(item['patch']['event'] for item in monitor.history)
    assert events['phase_start'] == events['partition_end'] == 1 and events['cone'] >= 2
    assert all(np.all(np.asarray(v) >= -1e-9) for v in result['inner'])   # 分区内为非负幅值
