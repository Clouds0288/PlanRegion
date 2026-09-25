"""算法顶点与展示网格的边界、缓存和离线回放回归。"""
from io import StringIO
from itertools import product
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
import json
import unittest

import numpy as np

from plot import RunMonitor, region_geometry, region_view, surface, json_value, pack_replay
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
        with RunMonitor(record=True, stream=StringIO()) as monitor:
            monitor('geometry', records=self.region.records.values(), bounds=self.region.bounds)
            first = json.loads(monitor.snapshot())['geometry']
            self.region.add_point([0], CUBE*.8)
            monitor('geometry', records=self.region.records.values(), bounds=self.region.bounds)
            self.assertEqual(monitor.history[0]['patch']['geometry'], first)
            self.assertNotEqual(monitor.state['geometry'], first)
            self.assertTrue(monitor.state['geometry'][0]['inner']['faces'])
            self.assertNotIn('records', monitor.state)

    def test_terminal_monitor_does_not_build_surfaces(self):
        with patch('plot.surface', side_effect=AssertionError('unexpected plotting')):
            result = self.region.finish(certified=False)
            with RunMonitor(record=False, stream=StringIO()) as monitor:
                monitor('geometry', records=self.region.records.values(), bounds=self.region.bounds)
                monitor('region_end', region=result, bounds=self.region.bounds)
                self.assertEqual(monitor.geometry_cache, {})
                self.assertEqual(monitor.history, [])

    def test_export_and_reload_without_event_recording(self):
        result = self.region.finish(certified=True)
        original = json.dumps(json_value(result))
        with TemporaryDirectory() as folder:
            with RunMonitor(record=False, output=folder, stream=StringIO()) as monitor:
                monitor('completed', results=[result], bounds=self.region.bounds)
                monitor.save_snapshot()
            self.assertEqual([path.name for path in Path(folder).iterdir()], ['live_view.html'])
            self.assertIn('window.SAVED_REPLAY_GZIP=',
                          (Path(folder)/'live_view.html').read_text(encoding='utf-8'))
            with RunMonitor(record=False, stream=StringIO()) as restored:
                restored.load_recording(Path(folder)/'live_view.html')
                saved = json.loads(restored.snapshot())
                self.assertEqual(saved['history'], [])
                self.assertTrue(saved['results'][0]['inner'][0]['faces'])
                self.assertEqual(saved['results'], monitor.state['results'])
        self.assertEqual(json.dumps(json_value(result)), original)

    def test_large_html_recording_restores_packed_history(self):
        with TemporaryDirectory() as folder:
            with RunMonitor(record=True, output=folder, stream=StringIO()) as monitor:
                monitor('phase_start', geometry=[{'vertices': [[0., -0., 1e-12]]}])
                monitor('point', point=[1., 2., 3.], query=1)
                monitor('completed', message='done')
                monitor.state['padding'] = 'x'*10_000_001
                with patch('plot.pack_replay', wraps=pack_replay) as packed:
                    monitor.save_snapshot()
                    packed.assert_called_once()
            with RunMonitor(record=False, stream=StringIO()) as restored:
                restored.load_recording(Path(folder)/'live_view.html')
                self.assertEqual(restored.history, monitor.history)
                self.assertEqual(restored.events, monitor.events)
                self.assertEqual(restored.state, monitor.state | {
                    'history_total': len(monitor.history), 'recording_version': 1})


if __name__ == '__main__':
    unittest.main()
