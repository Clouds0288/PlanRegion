"""RCUT 主线：锥几何与夹逼测度、细分与停滞规则、带 OBBT 紧化行的割平面循环，以及一个分区的端到端构域。"""
from collections import Counter
from types import SimpleNamespace

import numpy as np
import pytest
from threadpoolctl import threadpool_limits

import main
from model import new_model
from monitor import RunMonitor
from Network.case33bw import Case33, CURRENT_LIMIT
from region import (CUT_ACCEPTED, Cone, Cutting, Radial, box_vertices, build_partition, clip_box, cone_outer,
                    contains, halfspaces, piece, polytope_volume, sandwich, stagnated)


def settings(d):
    return SimpleNamespace(threads=1, obbt_workers=1, tau=main.REGION_TAU, discovery_eps=main.DISCOVERY_EPS,
                           discovery_share=main.DISCOVERY_SHARE, mip_seconds=main.MIP_SECONDS, mip_gap=main.MIP_GAP,
                           min_width=main.MIN_WIDTH[d], max_cones=main.MAX_CONES[d], threshold=main.CUT_THRESHOLD,
                           patience=main.CUT_PATIENCE, point_tol=main.POINT_TOL)


def radial(code, seconds=120.):
    network = Case33(load_nodes=(18, 25, 30)[:len(code)], current_limit=CURRENT_LIMIT)
    monitor = RunMonitor()
    monitor.time_limit = seconds
    part = Radial(network, tuple(1 if c == 'p' else -1 for c in code), Case33.switch_budget, monitor,
                  settings(len(code)))
    part.deadline = monitor.clock()+seconds
    return part


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


def test_stagnation_needs_patience_consecutive_small_cuts():
    assert stagnated([.5, .3, .001, .005, .009], .01, 3)
    assert not stagnated([.5, .001, .02, .005, .009], .01, 3)   # 中间一次 >= threshold 重新计数
    assert not stagnated([.001, .001], .01, 3)
    assert not stagnated([.5, .3, .01, .005, .009], .01, 3)     # 等于 threshold 不算小割（同主线）


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


def test_cut_loop_carries_obbt_rows_and_keeps_certified_ray_points():
    part = radial('nn')
    x = tuple(int(v) for v in part.network.encode_plan(part.network.initial_plan))   # 原点可行：两条射线顶点都存在
    cutting = Cutting(part)

    def rows():
        with new_model('rows', 1) as model:
            part.oracle._build(model, np.asarray(x), np.zeros(2))
            model.update()
            linear = model.NumConstrs
            model.remove(model.getQConstrs())   # 与 SubProblem._cut 相同：割 LP 只去掉锥约束
            model.update()
            return linear, model.NumConstrs
    with threadpool_limits(limits=1):
        plain = rows()
        cutting.run([x])
        tight = rows()
        added = 2*len(part.equations.boxes[x])+2*sum(x)
        assert tight[0]-plain[0] == tight[1]-plain[1] == added   # SP 与割 LP 都带全部盒约束行与反向锥包络行
        state = cutting.networks[x]
        assert state['status'] in CUT_ACCEPTED and cutting.count > 0
        points = np.array([part.vertex(x, key) for key in range(2)])   # 紧化射线认证点：有效割不会切掉
        assert contains(points, halfspaces(state['vertices']), 1e-7).all()
        assert [key for key, _ in cutting.sets()] == [(x, 0)]


def test_one_partition_certifies_with_cut_networks_and_records_frames():
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
