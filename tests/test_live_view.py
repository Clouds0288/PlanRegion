"""展示字段及浏览器端纯几何/回放回归，不修改求解状态。"""
from pathlib import Path
import subprocess
import unittest

from plot import region_geometry
from region import RegionState


class LiveViewTests(unittest.TestCase):
    def test_selection_is_exported_as_an_independent_list(self):
        region = RegionState([2., 3., 5.], 10., 0.)
        region.add_scheme([0, 1], {'e': 'k'}, 0.)
        shown = region_geometry(region.records.values(), region.bounds)
        self.assertEqual(shown[0]['x'], [0, 1])
        shown[0]['x'][0] = 1
        self.assertEqual(region.records[(0, 1)]['x'].tolist(), [0, 1])

    def test_viewer_geometry_and_replay(self):
        result = subprocess.run(['node', str(Path(__file__).with_name('live_view.test.cjs'))],
                                capture_output=True, text=True, encoding='utf-8', timeout=30)
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)


if __name__ == '__main__':
    unittest.main()
