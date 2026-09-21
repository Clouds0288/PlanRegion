"""优化前/独立/增量三组结果对比，以及按原始时间顺序核对去重证书。"""
from argparse import ArgumentParser
from pathlib import Path
import json

import numpy as np

from main import BenchmarkResult
from plot import RunMonitor
from region import contains, halfspaces, GEOMETRY_TOL
from vertify import METHODS


def audit_events(folder):
    state, domains = {}, []
    counts = dict(sp=0, ordinary_sp=0, ordinary_covered=0, skipped=0,
                  support_sp=0, covered_support_sp=0, unsupported_skips=0)
    for line in (Path(folder)/'events.jsonl').read_text(encoding='utf-8').splitlines():
        frame = json.loads(line)
        state.update(frame['patch'])
        if 'geometry' in frame['patch']:
            domains = [(row, halfspaces(np.asarray(row['inner']['vertices'])/state['bounds']))
                       for row in state['geometry'] if row['inner']['vertices']]
        event = state['event']
        if event not in ('sp_start', 'sp_skip'):
            continue
        p = np.asarray(state['point'])/state['bounds']
        covered = [(row, eq) for row, eq in domains if contains([p], eq)[0]]
        if event == 'sp_skip':
            counts['skipped'] += 1
            owners = [row for row, _ in covered if row['choice'] == state['covered_by']]
            assert owners, ('skip without a previous certificate', frame['id'])
            budget = np.inf if state['budget'] is None else state['budget']
            assert owners[0]['cost'] <= budget+1e-8
        else:
            counts['sp'] += 1
            mode = state.get('mode')
            if mode == 'SP-vertex':
                counts['ordinary_sp'] += 1
                counts['ordinary_covered'] += bool(covered)
                assert not covered, ('unnecessary ordinary SP', frame['id'])
            if mode == 'SP-gap-support':
                counts['support_sp'] += 1
                counts['covered_support_sp'] += bool(covered)
                witness = np.asarray(state['uncovered_witness'])/state['bounds']
                assert not any(contains([witness], eq)[0] for _, eq in domains), ('covered witness', frame['id'])
                support = np.asarray(state['support_points'])/state['bounds']
                assert contains([witness], halfspaces(support))[0], ('invalid simplex', frame['id'])
                assert any(np.allclose(p, vertex, rtol=0., atol=1e-12) for vertex in support)
                assert not any(row['choice'] == state['choice'] for row, _ in covered), ('own certificate repeated', frame['id'])
    result = BenchmarkResult.load(folder)
    assert counts['sp'] == sum(r['counts']['sp'] for r in result.metadata['continuous'])
    assert counts['skipped'] == sum(r['counts']['sp_skipped'] for r in result.metadata['continuous'])
    return dict(passed=True, geometry_tolerance=GEOMETRY_TOL, **counts)


def compare(baseline, independent, incremental):
    folders = dict(baseline=Path(baseline), independent=Path(independent), incremental=Path(incremental))
    results = {key: BenchmarkResult.load(folder) for key, folder in folders.items()}
    event_audits = {key: audit_events(folders[key]) for key in ('independent', 'incremental')}
    old = results['baseline']
    rows_by_key = {name: {(r['method'], r['budget']): r for r in result.metadata['continuous']}
                   for name, result in results.items()}
    totals, rows = {}, []
    for name, result in results.items():
        assert result.metadata['audit']['passed'], (name, 'missing independent audit')
        np.testing.assert_array_equal(result.states[METHODS.index('ac')], old.states[METHODS.index('ac')])
        conflict = (result.states*old.states == -1)
        assert not conflict.any(), (name, 'contradictory definite classifications')
        totals[name] = dict(seconds=sum(r['timing']['total_seconds'] for r in result.metadata['continuous']),
                            sp=sum(r['counts']['sp'] for r in result.metadata['continuous']),
                            mp_seconds=sum(r['timing']['mp_seconds'] for r in result.metadata['continuous']),
                            sp_seconds=sum(r['timing']['sp_seconds'] for r in result.metadata['continuous']),
                            geometry_seconds=sum(r['timing']['geometry_seconds'] for r in result.metadata['continuous']))
    for key, base in rows_by_key['baseline'].items():
        row = dict(method=key[0], budget=key[1], passed=True)
        for name in results:
            current = rows_by_key[name][key]
            assert current['status'] == 'certified', (name, key)
            assert current['coverage_bound'] is None or current['coverage_bound'] <= GEOMETRY_TOL
            assert abs(current['max_total']-base['max_total']) <= .003, (name, key, 'maximum')
            assert current['inner_volume'] <= base['outer_volume']*(1+1e-8)
            assert base['inner_volume'] <= current['outer_volume']*(1+1e-8)
            row[name] = dict(seconds=current['timing']['total_seconds'], sp=current['counts']['sp'],
                             max_total=current['max_total'], gap=current['volume_gap'])
        rows.append(row)
    report = dict(passed=True, folders={k:str(v.resolve()) for k,v in folders.items()},
                  totals=totals, rows=rows, event_audits=event_audits,
                  note='相同网架、预算、精度、单线程和完整记录；独立与增量计时分开；AC 与审核另计。')
    for name in ('independent', 'incremental'):
        target, result = folders[name], results[name]
        (target/'optimization.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
        result.metadata['optimization'] = report
        result.save(target)
        with RunMonitor(record=False, output=target, open_browser=False) as monitor:
            monitor.load_recording(target/'replay.json')
            monitor.state['metadata'] = result.metadata
            monitor.save_snapshot()
    return report


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('baseline')
    parser.add_argument('independent')
    parser.add_argument('incremental')
    args = parser.parse_args()
    print(json.dumps(compare(args.baseline, args.independent, args.incremental), ensure_ascii=False, indent=2))
