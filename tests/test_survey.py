"""Formal survey entry and conditional information value, without experiment imports."""
import ast
from contextlib import redirect_stdout
from io import StringIO
from itertools import product
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from shapely.geometry import shape

import main
import survey
from Network.concept5 import ROUTES
from region import GEOMETRY_TOL


ROOT = Path(__file__).resolve().parents[1]


class SurveyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reference = json.loads((ROOT/'tests/data/concept5_reference.json').read_text(encoding='utf-8'))
        cls.domains = {frozenset(row['allowed']): dict(inner=shape(row['inner']), outer=shape(row['outer']),
                                                      plan_ids=row['plan_ids'])
                       for row in cls.reference['domains']}

    def test_expected_gap_reduction_in_every_information_state(self):
        checked = 0
        for statuses in product((0, 1, -1), repeat=len(ROUTES)):
            yes = frozenset(r for r, s in zip(ROUTES, statuses) if s == 1)
            no = frozenset(r for r, s in zip(ROUTES, statuses) if s == -1)
            before = survey.information_state(self.domains, yes, no)
            for row in survey.candidate_values(self.domains, yes, no):
                route, probability = row['route'], row['road_usable_probability']
                plus = survey.information_state(self.domains, yes | {route}, no)
                minus = survey.information_state(self.domains, yes, no | {route})
                reduction = before['information_gap']-probability*plus['information_gap']-(1-probability)*minus['information_gap']
                self.assertAlmostEqual(reduction, row['marginal_information_value'], delta=1e-7)
                self.assertLessEqual(row['marginal_information_value_lower'], reduction+1e-5)
                self.assertGreaterEqual(row['marginal_information_value_upper'], reduction-1e-5)
                checked += 1
        self.assertEqual(checked, 1458)

    def test_first_choice_does_not_use_hidden_observations(self):
        for flip in (False, True):
            truth = {'A': flip, 'B': not flip, 'C': flip, 'D': not flip, 'E': True, 'F': False}
            result = survey.simulate_surveys(self.domains, truth)
            self.assertEqual(result['trace'][0]['route'], 'A')
            self.assertAlmostEqual(result['trace'][0]['marginal_information_value'],
                                   self.reference['result']['trace'][0]['marginal_information_value'])

    def test_default_main_dispatches_to_formal_survey(self):
        # 新默认入口迁为二维 FourBus；勘察从 survey.run_survey 单独启动。
        with patch.object(main, 'run', return_value='finished') as run:
            self.assertEqual(main.main(), 'finished')
        self.assertIsInstance(run.call_args.args[0], main.FourBus)
        self.assertEqual(run.call_args.args[0].load_nodes, (1, 2))


    def test_production_modules_do_not_import_experiments_or_tests(self):
        paths = [ROOT/name for name in ('main.py', 'model.py', 'region.py', 'plot.py', 'vertify.py',
                                       'survey.py', 'monitor.py')]
        paths.extend((ROOT/'Network').glob('*.py'))
        for path in paths:
            tree = ast.parse(path.read_text(encoding='utf-8-sig'))
            for node in ast.walk(tree):
                if isinstance(node, ast.ImportFrom):
                    names = [node.module or '']
                elif isinstance(node, ast.Import):
                    names = [alias.name for alias in node.names]
                else:
                    continue
                self.assertFalse(any(name.split('.')[0] in {'tests', 'experiments'} for name in names), path.name)
        self.assertNotIn('.codex', (ROOT/'plot.py').read_text(encoding='utf-8'))

    def test_standalone_run_matches_frozen_reference(self):
        with TemporaryDirectory() as directory:
            output = Path(directory)/'fresh'
            original, caches = survey.DomainCache, []
            def create_cache(**settings):
                caches.append(original(**settings))
                return caches[-1]
            with redirect_stdout(StringIO()), patch.object(survey, 'DomainCache', side_effect=create_cache):
                result = survey.run_survey(threads=1, time_limit=60.)
            cache = caches[0]
            saved = result  # 独立入口直接返回同一份结果，不再经 main 写网页记录。
            expected = self.reference['result']
            self.assertEqual([(r['route'], r['survey_observation']) for r in result['trace']],
                             [(r['route'], r['survey_observation']) for r in expected['trace']])
            self.assertEqual(result['uninspected'], ['F'])
            self.assertEqual(len(result['all_candidates']), 21)
            self.assertLessEqual(result['states'][-1]['information_gap_upper'], survey.STOPPING_AREA_TOLERANCE)
            self.assertLessEqual(result['remaining_max_marginal_upper'], survey.STOPPING_AREA_TOLERANCE)
            self.assertTrue(all(row['ranking_certified'] for row in result['trace']))
            for actual, old in zip(result['states'], expected['states'], strict=True):
                for name in ('confirmed_area', 'optimistic_area'):
                    self.assertAlmostEqual(actual[name], old[name], delta=.5)
            for allowed, row in cache.items():
                old = self.domains[allowed]
                inner, outer = row['inner'], row['outer']
                # 归一化几何容差换回 kW；方向上界会改变微小的边界舍入。
                tolerance = 135.*GEOMETRY_TOL
                self.assertLess(inner.difference(old['outer'].buffer(tolerance)).area, 1e-6)
                self.assertLess(old['inner'].difference(outer.buffer(tolerance)).area, 1e-6)
                self.assertLess(inner.symmetric_difference(old['inner']).area, .5)
                for subset, smaller in cache.items():
                    if subset < allowed:
                        for kind in ('inner', 'outer'):
                            self.assertLess(smaller[kind].difference(row[kind]).area, 1e-6)
            audit = saved['audit']
            self.assertTrue(audit['passed'])
            self.assertEqual(audit['ac_inner_vertices'], audit['ac_feasible_vertices'])
            self.assertTrue(audit['midpoint_infeasible'])
            self.assertEqual(saved['protocol']['schema'], 'survey-v2')
            scheme_ids = {row['id'] for row in saved['schemes']}
            self.assertIn(saved['nonconvexity_witness']['plan_a'], scheme_ids)
            self.assertIn(saved['nonconvexity_witness']['plan_b'], scheme_ids)
            self.assertNotIn('reference_directory', result)


if __name__ == '__main__':
    unittest.main()
