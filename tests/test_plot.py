"""算法顶点与展示网格的边界、缓存和离线回放回归。"""
from io import StringIO
from itertools import product
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
import json
import unittest

import numpy as np

from monitor import RunMonitor
from plot import region_geometry, region_view, surface, json_value, pack_replay
from region import RegionState
from plot import sample_region


CUBE = np.array(list(product((0., 1.), repeat=3)))


class PresentationTests(unittest.TestCase):
    def setUp(self):
        self.region = RegionState([2., 3., 5.], 10., .002)
        self.region.add_scheme([0], {'line': 'type0'}, 0.)
        self.region.add_point([0], CUBE*.5)

    def test_view_preserves_raw_vertices_certificates_and_classification(self):
        result = self.region.finish(certified=True)
        original = json.dumps(json_value(result))
        view = region_view(result, self.region.bounds)
        self.assertNotIn('geometry', result)
        self.assertNotIn('faces', result['inner'][0])
        self.assertTrue(view['inner'][0]['faces'])
        np.testing.assert_array_equal(view['inner'][0]['vertices'], result['inner'][0]['vertices'])
        self.assertEqual(view['inner'][0]['choice'], result['inner'][0]['choice'])
        self.assertNotIn('certificates', result)
        self.assertNotIn('cuts', result)
        self.assertEqual(region_view(view, self.region.bounds), view)
        points = np.array([[.2, .3, .4], [.9, .9, .9], [.5, .5, .5005]])*self.region.bounds
        np.testing.assert_array_equal(sample_region(result, points, self.region.bounds),
                                      sample_region(view, points, self.region.bounds))
        self.assertEqual(json.dumps(json_value(result)), original)

    def test_surface_preserves_degenerate_and_thin_domains(self):
        for poly in ([], [[.2, .3, .4]], CUBE*[1., 0., 0.], CUBE*[1., 1., 0.]):
            mesh = surface(poly, self.region.bounds)
            self.assertEqual(mesh['faces'], [])
            np.testing.assert_array_equal(np.asarray(mesh['vertices']).reshape(-1, 3),
                                          np.asarray(poly).reshape(-1, 3)*self.region.bounds)
        thin = CUBE*[1., .5, 2e-9]
        self.assertTrue(surface(thin, self.region.bounds)['faces'])

    def test_surface_cache_refreshes_after_cut_and_bound_change(self):
        cache = {}
        old = region_geometry(self.region.records.values(), self.region.bounds, cache)
        self.region.apply_cut(np.array([.75, -.5, 0., 0., 0.]))
        updated = region_geometry(self.region.records.values(), self.region.bounds, cache)
        self.assertEqual(old[0]['inner'], updated[0]['inner'])
        self.assertNotEqual(old[0]['outer'], updated[0]['outer'])
        scaled = region_geometry(self.region.records.values(), 2*self.region.bounds, cache)
        np.testing.assert_array_equal(scaled[0]['outer']['vertices'],
                                      2*np.asarray(updated[0]['outer']['vertices']))
        self.assertEqual(len(cache), 1)

    def test_recorded_geometry_is_independent_of_later_growth(self):
        monitor = RunMonitor()
        monitor._emit('geometry', **monitor._geometry(self.region))
        first = monitor.frame(0)
        self.region.add_point([0], CUBE*.8)
        monitor._emit('geometry', **monitor._geometry(self.region))
        self.assertEqual(monitor.frame(0), first)
        self.assertNotEqual(monitor.state['schemes'], first['schemes'])


    def test_terminal_monitor_does_not_build_surfaces(self):
        with patch('plot.surface', side_effect=AssertionError('unexpected plotting')):
            monitor = RunMonitor()
            monitor._emit('geometry', **monitor._geometry(self.region))
            self.assertNotIn('faces', monitor.state['schemes']['A'])


    def test_export_and_reload_without_event_recording(self):
        # 新协议统一保存帧，不再提供另一套静态 HTML 导出。
        with TemporaryDirectory() as folder:
            monitor = RunMonitor(output=Path(folder)/'monitor.json.gz')
            monitor._emit('completed', result=self.region.finish(True))
            monitor.save()
            restored = RunMonitor()
            restored.load_recording(monitor.output)
            self.assertEqual(restored.state, monitor.state)
            self.assertEqual([p.name for p in Path(folder).iterdir()], ['monitor.json.gz'])


    def test_large_html_recording_restores_packed_history(self):
        # 2.1：大记录同样只写原生窗口压缩帧，保留小数精度。
        with TemporaryDirectory() as folder:
            monitor = RunMonitor(output=Path(folder)/'monitor.json.gz')
            for i in range(1001):
                monitor._emit('point', point=[i*1e-12, 0.])
            monitor.save()
            restored = RunMonitor()
            restored.load_recording(monitor.output)
            self.assertEqual(restored.history, monitor.history)



if __name__ == '__main__':
    unittest.main()
