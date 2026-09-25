"""勘察图从结果生成，轮数、观测标签和输出位置不依赖全局状态。"""
from copy import deepcopy
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

from Network.concept5 import ROUTES
from survey import BUDGET, STOPPING_AREA_TOLERANCE, simulate_surveys
from plot import render_survey
from shapely.geometry import shape


class SurveyPlotTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = Path(__file__).parent/'data/concept5_reference.json'
        cls.reference = json.loads(path.read_text(encoding='utf-8'))
        cls.protocol = dict(budget=BUDGET, routes=ROUTES,
                            stopping_area_tolerance=STOPPING_AREA_TOLERANCE)

    def test_labels_and_matrix_follow_another_observation_sequence(self):
        domains = {frozenset(row['allowed']): dict(inner=shape(row['inner']), outer=shape(row['outer']),
                                                  plan_ids=row['plan_ids'])
                   for row in self.reference['domains']}
        result = simulate_surveys(domains, dict(A=False, B=True, C=False, D=True, E=True, F=False))
        result['protocol'] = self.protocol
        result['baseline_share'] = result['states'][0]['confirmed_area']/result['states'][0]['optimistic_area']
        # 此处只检查观测标签和评分；非凸面板使用固定绘图夹具。
        result['nonconvexity_witness'] = self.reference['result']['nonconvexity_witness']
        self.assertNotEqual(result['trace'], self.reference['result']['trace'])
        figures = {}
        def capture(fig, path, **kwargs):
            figures[Path(path).stem] = fig
        with TemporaryDirectory() as folder, patch.object(Figure, 'savefig', autospec=True, side_effect=capture):
            render_survey(result, folder)
        fig = figures['conditional_road_values']
        ax = next(ax for ax in fig.axes if ax.get_label() == 'conditional_value')
        labels = [label.get_text() for label in ax.get_xticklabels()]
        self.assertEqual(len(labels), len(result['states']))
        self.assertEqual(labels[1], '1 (A−)')
        for route in ROUTES:
            line = next(line for line in ax.lines if line.get_label() == route)
            scores = [row for row in result['all_candidates'] if row['route'] == route]
            self.assertEqual(list(line.get_xdata()), [row['step'] for row in scores])
            self.assertEqual(list(line.get_ydata()), [row['information_efficiency'] for row in scores])
        heat = next(ax for ax in fig.axes if ax.get_label() == 'score_matrix')
        self.assertEqual(heat.images[0].get_array().shape, (len(ROUTES), len(result['states'])))

    def test_render_preserves_input_and_global_style(self):
        result = deepcopy(self.reference['result'])
        result['protocol'] = self.protocol
        original, style = deepcopy(result), mpl.rcParams.copy()
        with TemporaryDirectory() as folder:
            output = Path(folder)/'figures'
            render_survey(result, output)
            self.assertEqual({path.name for path in output.iterdir()},
                             {f'{name}.{suffix}' for name in ('concept5_overview', 'conditional_road_values')
                              for suffix in ('pdf', 'svg', 'png')})
            self.assertTrue(all(path.stat().st_size > 1000 for path in output.iterdir()))
        self.assertEqual(result, original)
        self.assertEqual(dict(mpl.rcParams), dict(style))
        self.assertEqual(plt.get_fignums(), [])


if __name__ == '__main__':
    unittest.main()
