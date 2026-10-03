"""维数贯通：三维并集、独立扫描、误差分母与原生过程回放。"""
from itertools import product
from unittest.mock import patch

import numpy as np
import pytest

import main
from Network.case33bw import Case33, LOAD_NODES
from main import SCAN_DIVISIONS
from geometry import polytope_volume, union_measure
from monitor import COMPARISONS, NativeWindow, RunMonitor, grid_comparison
from plot import cut_polygon, voxel_faces
from tests.planning_checks import mainline_stubs, recorded_monitor, synthetic_reference


CUBE = np.array(list(product((0., 1.), repeat=3)))


def finished(monitor, d):
    """合成回放的终态结果：各锥的内外域行（带符号 kW）。"""
    rows = [row for row in monitor.state['cones'].values() if row]
    result = dict(status='time_limit', certified=False, inner=[dict(vertices=row['inner']) for row in rows],
                  outer=[dict(vertices=row['outer']) for row in rows])
    monitor._emit('region_end', phase='构域停止', result=result, partition=None)
    return result


@pytest.mark.parametrize('dimension', [2, 3])
def test_three_comparison_modes_switch_colors_and_show_six_metrics(dimension):
    monitor, _, _ = recorded_monitor(dimension)
    result = finished(monitor, dimension)
    states = np.ones((2,)*dimension, np.int8)
    states.flat[0] = states.flat[-1] = -1
    socp = states.copy()
    socp.flat[0] = 1
    bounds = np.full(dimension, 100.)
    monitor.validation(synthetic_reference(states, socp, bounds=np.zeros(dimension), axis_lower=-bounds), result)
    window = NativeWindow(monitor)
    try:
        window.root.withdraw()
        meshes = []
        for label in COMPARISONS.values():
            window.comparison_mode.set(label)
            window._refresh_validation()
            window.canvases['C'].draw()
            axis = window.axes['C']
            assert label in axis.get_title()
            assert [item.get_text() for item in axis.get_legend().get_texts()][1:] == ['遗漏', '多余']
            lines = window.comparison_text.get().splitlines()
            assert len(lines) == 5 and all('遗漏' in line and '多余' in line for line in lines[:3])
            assert lines[3].startswith('本组遗漏按分区') and lines[4].startswith('本组多余按分区')
            if dimension == 2:
                meshes.append(np.asarray(axis.collections[0].get_array()).copy())
            else:
                assert {item.get_gid() for item in axis.collections} == {'reference', 'missed-cells', 'extra-cells'}
        if dimension == 2:
            for i in range(3):
                for j in range(i):
                    assert not np.array_equal(meshes[i], meshes[j])
    finally:
        window.close()


def test_union_volume_nearly_collinear_overlap_preserves_small_growth():
    import json
    from pathlib import Path
    samples = json.loads((Path(__file__).parent/'data/union_volume_regression.json').read_text())
    previous = None
    for sample in samples:
        polytopes = [np.asarray(p) for p in sample['polytopes']]
        volume = union_measure(polytopes, 3)
        assert volume == pytest.approx(sample['expected'], abs=1e-11)
        assert union_measure(polytopes[::-1], 3) == pytest.approx(volume, abs=1e-11)
        if previous is not None:
            # 只有 A 增加一个边界点，并集增量不能超过 A 本身的增量。
            delta = polytope_volume(polytopes[0])-polytope_volume(previous[0])
            assert -1e-11 <= volume-previous[1] <= delta+1e-11
        previous = polytopes[0], volume


def test_case33_node_configuration_reaches_all_dimensions(tmp_path):
    assert Case33().load_nodes == LOAD_NODES
    with mainline_stubs() as calls:
        main.main('case33', 3, show=False, output=tmp_path)
    network, monitor = calls.build.call_args.args[0], calls.build.call_args.kwargs['monitor']
    assert network.load_nodes == (18, 25, 30)
    assert monitor.output.name == 'case33_18_25_30.json.gz'
    selected = np.isin(network.nodes, network.load_nodes)
    np.testing.assert_array_equal(network.fixed_p[~selected], network.original_p[~selected])
    assert np.all(network.fixed_p[selected] == 0.)


@pytest.mark.parametrize('case,nodes', [('fourbus', (1, 2, 3)), ('case33', (18, 25, 30))])
@pytest.mark.parametrize('dimension', [2, 3])
def test_dimension_setting_controls_nodes_and_recording(case, nodes, dimension, tmp_path):
    with mainline_stubs() as calls:
        main.main(case, dimension, show=False, output=tmp_path)
    network, monitor = calls.build.call_args.args[0], calls.build.call_args.kwargs['monitor']
    assert network.load_nodes == nodes[:dimension]
    assert calls.build.call_args.kwargs['settings'].min_width == main.MIN_WIDTH[dimension]
    assert monitor.output == tmp_path/'mode_1'/f"{case}_{'_'.join(map(str, nodes[:dimension]))}.json.gz"
    calls.validate.assert_called_once_with(monitor)   # 构域结束后校验同一份记录，再由记录画收敛过程图
    calls.curves.assert_called_once_with(monitor)
    calls.draw.assert_called_once_with(calls.curves.return_value, dimension*main.REGION_TAU, monitor.output)


@pytest.mark.parametrize('dimension', [2, 3])
def test_validation_reads_network_and_scan_divisions_from_the_record(dimension, tmp_path):
    nodes = (18, 25, 30)[:dimension]
    result = dict(inner=[], outer=[], axis_lower=-np.ones(dimension), axis_bounds=np.ones(dimension))
    monitor = RunMonitor(output=tmp_path/'record.json.gz')
    monitor.state.update(network='case33bw', load_nodes=list(nodes), budget=7, result=result)
    with patch('main.reference_grid') as grid, patch('main.export_comparison') as export, \
         patch.object(monitor, 'validation') as validation:
        assert main.validate(monitor) is export.return_value
    network, budget, scanned, output = grid.call_args.args
    assert network.load_nodes == nodes and budget == 7 and scanned is result and output == main.SCAN_OUTPUT
    assert np.isfinite(network.ell_limit).all()   # 与构域相同的主线电流限额
    assert grid.call_args.kwargs['divisions'] == SCAN_DIVISIONS[dimension]
    validation.assert_called_once_with(grid.return_value, result)
    assert export.call_args.args[-1] == monitor.output   # 对比导出到记录旁


def test_three_dimensional_measure_keeps_overlap_gaps_and_face_ownership():
    assert union_measure([CUBE, CUBE+[.5, .5, 0.]], 3) == pytest.approx(1.75)
    assert union_measure([CUBE, CUBE+[2., 0., 0.]], 3) == pytest.approx(2.)
    assert union_measure([CUBE, CUBE+[1., 0., 0.]], 3) == pytest.approx(2.)
    assert union_measure([CUBE, CUBE.copy(), .5*CUBE+.25], 3) == pytest.approx(1.)
    # 三个盒子围出一个空隙，不能用全体顶点的凸包替代并集。
    pieces = [CUBE*[.2, 1., 1.], CUBE*[.2, 1., 1.]+[.8, 0., 0.], CUBE*[1., .2, 1.]]
    assert union_measure(pieces, 3) == pytest.approx(.52)
    basis, _ = np.linalg.qr([[1., 2., 3.], [3., 1., 4.], [2., 5., 1.]])
    rotated = [p@basis.T+[.3, .4, .5] for p in pieces]
    assert union_measure(rotated, 3) == pytest.approx(.52, abs=1e-9)


def test_three_dimensional_metrics_count_cells_and_declared_denominators():
    truth = np.ones((2, 2, 2), dtype=int)
    truth[0, 0, 0] = -1
    result = dict(inner=[dict(vertices=CUBE*[1., 2., 2.])], outer=[dict(vertices=CUBE*2.)])
    metrics = grid_comparison(dict(axis_lower=np.zeros(3), bounds=np.full(3, 2.), states=truth, socp_states=truth),
                              result)['metrics']
    assert metrics['inner'] == dict(mr_percent=400/7, fr_percent=25., missed_cells=4,
                                  extra_cells=1, reference_cells=7, computed_cells=4)
    assert metrics['outer']['mr_percent'] == 0.
    assert metrics['outer']['fr_percent'] == 12.5
    assert metrics['outer']['computed_cells'] == 8


def test_voxel_surface_keeps_disconnected_regions_and_physical_scale():
    states = -np.ones((3, 3, 3), dtype=int)
    states[0, 0, 0] = states[2, 2, 2] = 1
    faces = voxel_faces(states, [3., 6., 9.])
    assert faces.shape == (12, 4, 3)
    np.testing.assert_array_equal(faces.min(axis=(0, 1)), [0., 0., 0.])
    np.testing.assert_array_equal(faces.max(axis=(0, 1)), [3., 6., 9.])
    assert not any(np.all((face.mean(axis=0) > [1., 2., 3.]) & (face.mean(axis=0) < [2., 4., 6.])) for face in faces)
    assert voxel_faces(np.ones((2, 2, 2)), [2., 2., 2.]).shape == (24, 4, 3)


def test_three_dimensional_replay_tracks_cones_cuts_and_validation(tmp_path):
    monitor, x, cut = recorded_monitor(3, tmp_path/'three.json.gz')
    polygon = cut_polygon(cut, x, np.full(3, 100.))
    np.testing.assert_allclose(polygon.sum(axis=1), 150., atol=1e-9)
    assert len(polygon) == 6
    cut_index = next(i for i, item in enumerate(monitor.history) if item['patch']['event'] == 'cut')
    result = finished(monitor, 3)
    monitor.validation(synthetic_reference(np.ones((2,)*3), bounds=np.zeros(3), axis_lower=-np.full(3, 100.)), result)
    monitor.save()
    restored = RunMonitor()
    restored.load_recording(monitor.output)
    window = NativeWindow(restored)
    try:
        window.root.withdraw()
        window.seek(cut_index)
        scheme_axis = window.scheme_views['---:1'][1]
        assert scheme_axis.name == '3d' and 'p_{3}' in scheme_axis.get_zlabel()
        removed = window._removed_3d(restored.frame(cut_index), '---:1')
        assert sum(polytope_volume(p) for p in removed) == pytest.approx(500000.)
        window.seek(cut_index-1)   # 割之前的评分帧：网架已出现，尚无割
        assert not any(item.get_gid() == 'cut----:1----:1' for item in window.scheme_views['---:1'][1].collections)
        window.seek(0)
        assert not window.scheme_views['---:1'][0].winfo_manager()
        assert any(item.get_gid() == 'reference' for item in window.axes['C'].collections)
        for canvas in window.canvases.values():
            canvas.draw()
        for _, _, canvas in window.scheme_views.values():
            canvas.draw()
    finally:
        window.close()


@pytest.mark.parametrize('shape', [(1, 3), (2, 3, 4)])
def test_ac_panel_renders_rectangular_grid_and_single_row(shape):
    monitor, _, _ = recorded_monitor(len(shape))
    result = finished(monitor, len(shape))
    states = np.ones(shape, dtype=np.int8)
    states.flat[0] = -1
    bounds = 10.*np.arange(1, len(shape)+1)
    monitor.validation(synthetic_reference(states, bounds=bounds, axis_lower=-bounds), result)
    window = NativeWindow(monitor)
    try:
        window.root.withdraw()
        window.show()
        window.canvases['C'].draw()
        title = window.axes['C'].get_title()
        assert '×'.join(map(str, shape)) in title and 'AC' in title
        if len(shape) == 2:
            coordinates = window.axes['C'].collections[0].get_coordinates()
            np.testing.assert_allclose(coordinates.max(axis=(0, 1)), bounds)
    finally:
        window.close()
