"""完整 AC、区域差集和唯一结果文件的回归验证。"""
import tempfile
import unittest
import numpy as np
import gurobipy as gp
from Network.four_bus_five_corridor import network
from plot import voxel_surface
from model import ACPowerFlow
from vertify import BenchmarkResult, disagreement


class ACReferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.models = [ACPowerFlow(d) for d in network.designs]
        cls.environment = gp.Env(empty=True)
        cls.environment.setParam('OutputFlag', 0)
        cls.environment.start()

    @classmethod
    def tearDownClass(cls):
        for model in cls.models:
            model.close()
        cls.environment.dispose()

    def test_all_designs_against_explicit_global_ac(self):
        rng = np.random.default_rng(20260919)
        counts = {-1: 0, 1: 0}
        for model in self.models:
            power = rng.dirichlet(np.ones(4))[:3]*network.power_limit
            actual = int(model.classify(power)[0])
            expected = model.global_status(power, self.environment)
            self.assertNotEqual(expected, 0, 'QCQP 未完成判定')
            self.assertEqual(actual, expected, (model.network.choice, power))
            counts[actual] += 1
        self.assertEqual(len(self.models), 216)
        self.assertGreater(counts[-1], 0)
        self.assertGreater(counts[1], 0)

    def test_ac_witness_satisfies_complex_nodal_power_flow(self):
        power = np.array([5., 7., 4.])
        for model in self.models[::9]:
            c = model.network
            status, ell = model.classify(power, return_currents=True)
            self.assertEqual(status[0], 1)
            P, Q, v, _ = model.state(power, ell)
            voltage, current = np.ones(c.n, dtype=complex), np.zeros(c.n, dtype=complex)
            for i in c.order:
                parent = 1+0j if c.parent[i] < 0 else voltage[c.parent[i]]
                current[i] = np.conj((P[0, i]+1j*Q[0, i])/parent)
                voltage[i] = parent-(c.r[i]+1j*c.reactance[i])*current[i]
            demand = voltage*np.conj(current-np.array([current[j].sum() for j in c.children]))
            p, q = c.loads(power)
            np.testing.assert_allclose(demand, (p+1j*q)[0], atol=1e-10, rtol=0)
            np.testing.assert_allclose(abs(voltage)**2, v[0], atol=1e-10, rtol=0)

    def test_transformer_checks_sending_end_power(self):
        power = np.full(3, network.power_limit/3)
        for model in self.models:
            self.assertEqual(model.classify(power)[0], -1)


class RegionComparisonTests(unittest.TestCase):
    def test_symmetric_difference_counts_both_sides(self):
        labels = np.array([0, 0, 1, 2, 2, 3, 3, 3])
        row = disagreement((labels & 1)>0, (labels & 2)>0)
        self.assertAlmostEqual(row['region_error_percent'], 50.)
        self.assertAlmostEqual(row['fr_percent'], 25.)
        self.assertAlmostEqual(row['mr_percent'], 40.)

    def test_empty_region_has_undefined_conditional_rate(self):
        for label, fr, mr in [(0, None, None), (1, 100., None), (2, None, 100.)]:
            labels = np.full(8, label)
            row = disagreement((labels & 1)>0, (labels & 2)>0)
            self.assertEqual(row['fr_percent'], fr)
            self.assertEqual(row['mr_percent'], mr)

    def test_surface_preserves_cavity_and_volume(self):
        mask = np.ones((3, 3, 3), dtype=bool)
        mask[1, 1, 1] = False
        h = np.array([.75, 1.5, 2.])
        mesh = voxel_surface(mask, h)
        points = np.asarray(mesh['vertices'])
        triangles = points[np.asarray(mesh['triangles'])]
        volume = np.einsum('ij,ij->i', triangles[:, 0],
                           np.cross(triangles[:, 1], triangles[:, 2])).sum()/6
        self.assertAlmostEqual(volume, mask.sum()*h.prod())

    def test_result_roundtrip_preserves_masks_and_deduplicates_geometry(self):
        masks = np.random.default_rng(4).random((4, 2, 3, 3, 3)) > .5
        ac_cost = np.random.default_rng(5).choice([0., 20000., 40000., np.inf], (3,3,3))
        masks[2,0] = ac_cost<=20000.
        masks[2,1] = np.isfinite(ac_cost)
        record = dict(inner=[[0, 0, 0]], outer=[[1, 2, 3]])
        result = BenchmarkResult(masks, dict(divisions=3, budgets=[20000., None]),
                                 {'socp': [[record], [record]]}, ac_cost)
        with tempfile.TemporaryDirectory(prefix='planregion-test-') as folder:
            result.save(folder)
            restored = BenchmarkResult.load(folder)
            np.testing.assert_array_equal(restored.masks, masks)
            self.assertEqual(restored.regions, result.regions)
            import json
            with np.load(f'{folder}/result.npz') as data:
                self.assertEqual(set(data.files), {'masks', 'ac_cost', 'record'})
                np.testing.assert_array_equal(data['ac_cost'], ac_cost)
                self.assertEqual(len(json.loads(str(data['record']))['geometry']), 1)


if __name__ == '__main__':
    unittest.main()
