"""Reproduce the isolated Jiangkou existing / infinite-budget 2D experiment."""
from pathlib import Path
from time import perf_counter
import hashlib
import json
import sys
import traceback

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'source'))
sys.stdout.reconfigure(encoding='utf-8', line_buffering=True)
import numpy as np
from threadpoolctl import threadpool_limits
from Network.jiangkou import Jiangkou
from main import ContinuousRegion, RegionTimeout
from model import PLANNING_TOL
from region import initial_polytope, polytope_volume, halfspaces, contains, clip_polytope
from vertify import ACPowerFlow
from experiment_runtime import save, serial

LIMIT = 3600.
BOUNDS = None
solver = None
started = perf_counter()
elapsed_before = 0.


def observe(event, **data):
    row = dict(seconds=elapsed_before+perf_counter()-started, event=event,
               **{k:v for k,v in data.items() if k not in ('state', 'region', 'geometry')})
    with (ROOT/'results'/'events.jsonl').open('a', encoding='utf-8') as stream:
        stream.write(json.dumps(serial(row), ensure_ascii=False)+'\n')
    if event in ('geometry', 'mode_switch', 'maximum', 'region_end') and solver is not None:
        checkpoint = solver.region.finish(False)
        checkpoint.update(status='running', counts=solver.counts, seconds=elapsed_before+perf_counter()-started,
                          residual_mode=solver.residual_mode, bounds=BOUNDS)
        save('domain_checkpoint.json', checkpoint)
    if event in ('geometry', 'query_end', 'coverage', 'mode_switch', 'sp_end', 'mp_start'):
        if event!='geometry' or solver is None or solver.counts['cuts']%20==0:
            print(f'{elapsed_before+perf_counter()-started:.1f}s {event}: '+str(data.get('message', ''))[:180])


def ac_seed(equations, choice, power):
    design = equations.network.design(choice)
    ac = ACPowerFlow(design)
    status, ell = ac.classify(power, return_currents=True)
    if status[0] != 1:
        return None
    x = equations.selection(choice)
    state = np.zeros(len(equations.upper))
    active = equations.selected[np.flatnonzero(x)]
    state[equations.ell.start+active] = ell[0, equations.receivers[active]]
    state = equations.restore(x, power, state)
    margin = equations.margin(x, power, state)
    if margin < -PLANNING_TOL:
        return None
    np.savez_compressed(ROOT/'results'/'bootstrap_certificate.npz',
                        x=x, power=power, state=state, ell=ell, choice=choice, margin=margin)
    return x, state


def main():
    global solver, BOUNDS, elapsed_before
    # Exercise genuine 2D corner cases before any long computation.
    triangle = initial_polytope(np.ones(2), 1.)
    assert np.isclose(polytope_volume(triangle), .5)
    assert contains([[.2,.3],[1.,1.]], halfspaces(triangle)).tolist() == [True,False]
    assert np.isclose(polytope_volume(clip_polytope(triangle,.4,[-1,0])), .32)
    assert contains([[0,0]],halfspaces(np.array([[0.,0.]]))).all()
    network = Jiangkou(load_nodes=('B000025','B000078'))
    BOUNDS = np.repeat(network.power_limit, 2)
    choice = np.array([len(o.cost)-1 for o in network.line_options])
    nominal = network.original_p[network.selected]
    preflight = dict(nodes=network.n+1, branches=network.n,
                     building_ids=['B000025','B000078'], load_node_ids=network.load_nodes,
                     nominal_kw=nominal, fixed_load_kw=network.fixed_p.sum(),
                     transformer_kva=network.base, bounds_kw=BOUNDS,
                     binary_options=sum(len(o.cost) for o in network.line_options),
                     budget=None, time_limit_seconds=LIMIT, tau=.002,
                     baseline_ac=ACPowerFlow(network).classify([nominal,[0,0]]),
                     all_max_ac=ACPowerFlow(network.design(choice)).classify([nominal,[0,0]]))
    save('preflight.json', preflight)
    print(json.dumps(serial(preflight),ensure_ascii=False))
    final_path=ROOT/'results'/'domain.json'
    prior=json.loads(final_path.read_text(encoding='utf-8')) if final_path.exists() else {}
    if 'envs\\methods\\' in prior.get('python',''):
        elapsed_before=prior['timing']['total_seconds']
    solver = ContinuousRegion(network,'socp',np.inf,BOUNDS,tau=.002,
                              observer=observe,residual_mode='physical',time_limit=max(0.,LIMIT-elapsed_before))
    solver.started-=elapsed_before
    solver.time_limit=LIMIT
    if elapsed_before:
        solver.times.update({k:prior['timing'][k] for k in solver.times})
        solver.phase_seconds.update(prior['phases'])
        solver.counts.update(prior['counts'])
    print(f'{perf_counter()-started:.1f}s equations_ready {solver.equations.G.shape}')
    try:
        previous_path = ROOT/'results'/'domain_checkpoint.json'
        previous = json.loads(previous_path.read_text(encoding='utf-8')) if previous_path.exists() else None
        if previous and previous['certificates']:
            solver.residual_mode=previous.get('residual_mode','physical')
            solver.region.cuts = [np.asarray(cut) for cut in previous['cuts']]
            for record in previous['certificates']:
                x = np.asarray(record['x'],int)
                solver.region.add_scheme(x,record['choice'],record['cost'])
                solver.region.add_point(x,record['inner'])
            solver.counts['reused_cuts'] = len(solver.region.cuts)
            solver.counts['reused_schemes'] = len(previous['certificates'])
            solver.geometry(message='继承已保存的切割域，直接继续全局剩余域搜索')
        else:
            seed = ac_seed(solver.equations, choice, nominal)
            if seed is not None:
                x, _ = seed
                solver.refine(x, nominal)
        # Domain coverage needs a valid total upper bound, not a separate proof
        # of maximum total load or minimum investment at an infinite budget.
        while True:
            witness = solver.residual()
            observe('coverage',bound=witness['bound'],complete=witness['complete'])
            if witness['complete']:
                result = solver.finish('certified',witness['bound'],None)
                break
            if witness['x'] is None:
                if solver.residual_mode=='physical':
                    solver.residual_mode='light'
                    observe('mode_switch',message='物理剩余域搜索较慢，继承全部结果转轻量模式')
                    continue
                result = solver.finish('unknown',witness['bound'],None)
                break
            before = solver.region.progress
            ok = solver.refine(witness['x'],
                seed=witness['p'] if witness.get('feasible') else None,witness=witness['p'])
            if not ok or before==solver.region.progress:
                result = solver.finish('unknown',witness['bound'],None)
                break
    except RegionTimeout:
        result = solver.finish('time_limit',None,None)
    except Exception:
        save('error.json',dict(traceback=traceback.format_exc()))
        result = solver.finish('partial_error',None,None)
    result.update(bounds=BOUNDS, building_ids=preflight['building_ids'],
                  nominal_kw=nominal, dimension=2, measure_unit='kW^2',
                  python=sys.executable,
                  experiment_seconds=elapsed_before+perf_counter()-started)
    save('domain.json',result)
    original = json.loads((ROOT/'source_manifest.json').read_text(encoding='utf-8'))
    changed = [p for p,h in original.items()
               if hashlib.sha256((ROOT.parents[1]/p).read_bytes()).hexdigest()!=h]
    save('source_integrity.json',dict(original_sources_unchanged=not changed,
        externally_changed_files=changed,test_writes='tests/jiangkou_2d only; never restore production files'))
    print(json.dumps(serial({k:result[k] for k in ('status','schemes','counts','inner_volume','outer_volume','volume_gap','timing')}),ensure_ascii=False))


if __name__ == '__main__':
    with threadpool_limits(limits=1):
        main()
