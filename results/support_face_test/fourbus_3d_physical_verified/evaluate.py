"""Audit the frozen physical run on the unchanged independent reference grid."""
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from threadpoolctl import threadpool_limits

folder = Path(__file__).resolve().parent
workspace = folder.parents[2]
sys.path.insert(0, str(folder/'source'))
from monitor import RunMonitor
from model import PortPhysics, PLANNING_TOL
from Network.four_bus_five_corridor import FourBus
sys.path.append(str(workspace))
from tests.planning_checks import margin

comparison = json.loads((folder/'comparison.json').read_text(encoding='utf-8'))
scan = workspace/'results/scans/fourbus_signed/mode_1/1_2_3/pf_0.95_1_+1_budget_20000.npz'
with np.load(scan) as saved:
    reference = {k: saved[k].copy() for k in ('axis_lower', 'bounds', 'states')}
comparison['reference'] = dict(path=str(scan), sha256=hashlib.sha256(scan.read_bytes()).hexdigest(),
    shape=list(reference['states'].shape), axis_lower=reference['axis_lower'].tolist(),
    bounds=reference['bounds'].tolist(), scope='Existing independent reference, metrics restricted to this grid')
with threadpool_limits(limits=1):
    run = comparison['runs'][0]
    monitor = RunMonitor(output=folder/'support_1/monitor.json.gz')
    monitor.load_recording(monitor.output)
    monitor.validation(reference, run['result'])
    monitor.save()
    run['metrics'] = monitor.validation_state['validation']['metrics']
    run['frames'] = len(monitor.history)
    run['partitions'] = []
    residuals = []
    for path in sorted((folder/'support_1').glob('partition_*.json')):
        saved = json.loads(path.read_text(encoding='utf-8'))
        sign = np.asarray(saved['sign'])
        equations = PortPhysics(FourBus(load_nodes=(1, 2, 3)), sign)
        run['partitions'].append(dict(sign=saved['sign'], **saved['coverage']))
        for row in saved['schemes']:
            for point in row['certificates']:
                residuals.append(-margin(equations, np.asarray(row['x']), np.asarray(point['p'])*sign,
                                         np.asarray(point['state'])))
        for step in saved['discovery_history']:
            if step.get('candidate_feasible'):
                residuals.append(-margin(equations, np.asarray(step['x']), np.asarray(step['p'])*sign,
                                         np.asarray(step['state'])))
    comparison['witness_audit'] = dict(witnesses=len(residuals), maximum_violation=max(residuals),
        threshold=PLANNING_TOL, passed=max(residuals) <= PLANNING_TOL)
assert comparison['witness_audit']['passed']
(folder/'comparison.json').write_text(json.dumps(comparison, indent=2, ensure_ascii=False), encoding='utf-8')
print(run['result']['timing'], run['result']['counts'], flush=True)
print(run['metrics'], flush=True)
print(comparison['witness_audit'], flush=True)
