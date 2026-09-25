"""五节点道路勘察：共享对偶割、条件信息价值与停止证书。"""
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
from datetime import datetime
import hashlib
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from shapely.geometry import GeometryCollection, MultiPoint, Point, mapping
from shapely.ops import unary_union
from threadpoolctl import threadpool_limits

from Network.concept5 import Concept5, ROUTES
from main import build_continuous_region
from model import PlanningEquations, PlanningModel
from plot import json_value
from vertify import ACPowerFlow

ROOT = Path(__file__).resolve().parent
BUDGET = 4.
STOPPING_AREA_TOLERANCE = 5.
DOMAIN_AREA_TOLERANCE = .5
# 合成观测场景；F 没有预设可用或不可用真值。
SURVEY_REALIZATION = {'A': True, 'B': False, 'C': True, 'D': False, 'E': True}



def write_json(path, value):
    path.write_text(json.dumps(json_value(value), ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')


def polygon_union(rows):
    return unary_union([MultiPoint(row['vertices']).convex_hull for row in rows]) if rows else GeometryCollection()


class DomainCache(dict):
    """Lazy road-mask queries with a global physical cut pool; no construction enumeration."""

    def __init__(self, output, *, threads=1, time_limit=60., share_cuts=True):
        super().__init__()
        self.output = Path(output)
        self.output.mkdir(parents=True, exist_ok=True)
        self.threads, self.time_limit, self.share_cuts = threads, time_limit, share_cuts
        self.cuts, self.queries, self.regions = [], [], {}

    def __missing__(self, allowed):
        allowed = frozenset(allowed)
        network = Concept5(allowed)
        name = ''.join(sorted(allowed)) or 'baseline'
        started = perf_counter()
        initial_cuts = len(self.cuts) if self.share_cuts else 0
        local_cuts = list(self.cuts) if self.share_cuts else []
        attempts = []

        def progress(event, **data):
            if event == 'cut':
                local_cuts.append(np.asarray(data['cut']).copy())

        # Identical geometry and tolerance; physical residual search is a recorded fallback.
        for mode in ('light', 'physical'):
            result = build_continuous_region(
                network, 'socp', BUDGET, np.array([135., 135.]), tau=1e-5,
                residual_mode=mode, time_limit=self.time_limit, threads=self.threads,
                progress=progress, cuts=local_cuts,
            )
            inner, outer = (polygon_union(result[key]) for key in ('inner', 'outer'))
            attempts.append({'mode': mode, 'status': result['status'], 'counts': result['counts'],
                             'seconds': result['timing']['total_seconds'], 'area_gap': outer.area-inner.area})
            write_json(self.output/(name+'_'+mode+'.json'), result)
            if result['status'] == 'certified' and outer.area-inner.area <= DOMAIN_AREA_TOLERANCE:
                break
        else:
            raise RuntimeError(f'Uncertified road mask {name}: {attempts}')
        # Independent approximations need not be nested. Preserve inclusion using
        # certified subsets/supersets, without changing any previously scored domain.
        smaller = [value for key, value in self.items() if key < allowed]
        larger = [value for key, value in self.items() if allowed < key]
        inner = unary_union([inner, *(value['inner'] for value in smaller)])
        for value in larger:
            inner = inner.intersection(value['inner'])
            outer = outer.intersection(value['outer'])
        outer = unary_union([outer, *(value['outer'] for value in smaller)])
        if outer.area-inner.area > DOMAIN_AREA_TOLERANCE:
            raise RuntimeError(f'Monotone area interval too wide for {name}: {outer.area-inner.area}')
        if self.share_cuts:
            self.cuts = local_cuts
        self.regions[allowed] = result
        domain = {'inner': inner, 'outer': outer, 'plan_ids': None}
        self[allowed] = domain
        query = {'allowed': sorted(allowed), 'seconds': perf_counter()-started,
                 'initial_cuts': initial_cuts, 'new_cuts': len(local_cuts)-initial_cuts,
                 'scheme_count': len(result['inner']), 'inner_area': inner.area,
                 'outer_area': outer.area, 'area_gap': outer.area-inner.area, 'attempts': attempts}
        self.queries.append(query)
        write_json(self.output/'queries.json', self.queries)
        print(json.dumps(query, ensure_ascii=False), flush=True)
        return domain


def increase_bounds(before, after):
    if before['plan_ids'] is not None and before['plan_ids'] == after['plan_ids']:
        return (0., 0., 0.)
    if before['inner'].difference(after['inner']).area > 1e-5:
        raise RuntimeError('Expected monotonic domain inclusion failed')
    return (after['inner'].difference(before['inner']).area,
            after['inner'].difference(before['outer']).area,
            after['outer'].difference(before['inner']).area)


def candidate_values(domains, known_usable, known_unusable):
    allowed = frozenset(ROUTES)-known_unusable
    scores = []
    for route, data in ROUTES.items():
        if route in known_usable or route in known_unusable:
            continue
        gain = increase_bounds(domains[known_usable], domains[known_usable | {route}])
        removal = increase_bounds(domains[allowed-{route}], domains[allowed])
        probability, cost = data['road_usable_probability'], data['survey_cost']
        row = {'route': route, 'road_usable_probability': probability, 'survey_cost': cost}
        for index, suffix in enumerate(('', '_lower', '_upper')):
            row['gain_if_usable'+suffix], row['removal_if_unusable'+suffix] = gain[index], removal[index]
            row['expected_expansion'+suffix] = probability*gain[index]
            row['expected_exclusion'+suffix] = (1-probability)*removal[index]
            row['marginal_information_value'+suffix] = probability*gain[index]+(1-probability)*removal[index]
            row['information_efficiency'+suffix] = row['marginal_information_value'+suffix]/cost
        scores.append(row)
    return sorted(scores, key=lambda row: (-row['information_efficiency'], row['route']))


def information_state(domains, known_usable, known_unusable):
    confirmed, optimistic = domains[known_usable], domains[frozenset(ROUTES)-known_unusable]
    return {'known_usable': sorted(known_usable), 'known_unusable': sorted(known_unusable),
            'confirmed': mapping(confirmed['inner']), 'optimistic': mapping(optimistic['inner']),
            'confirmed_outer': mapping(confirmed['outer']), 'optimistic_outer': mapping(optimistic['outer']),
            'confirmed_area': confirmed['inner'].area, 'confirmed_area_upper': confirmed['outer'].area,
            'optimistic_area': optimistic['inner'].area, 'optimistic_area_upper': optimistic['outer'].area,
            'information_gap': optimistic['inner'].difference(confirmed['inner']).area,
            'information_gap_lower': optimistic['inner'].difference(confirmed['outer']).area,
            'information_gap_upper': optimistic['outer'].difference(confirmed['inner']).area}


def simulate_surveys(domains, survey_realization):
    known_usable, known_unusable = frozenset(), frozenset()
    states, trace, all_candidates, total_cost = [], [], [], 0.
    for step in range(len(ROUTES)+1):
        before = information_state(domains, known_usable, known_unusable)
        before.update(step=step, survey_total_cost=total_cost)
        states.append(before)
        scores = candidate_values(domains, known_usable, known_unusable)
        all_candidates.extend({'step': step, **row} for row in scores)
        maximum_upper = max((row['marginal_information_value_upper'] for row in scores), default=0.)
        if maximum_upper <= STOPPING_AREA_TOLERANCE and before['information_gap_upper'] <= STOPPING_AREA_TOLERANCE:
            return {'states': states, 'trace': trace, 'all_candidates': all_candidates,
                    'stop_reason': 'residual_domain_gap_and_all_marginals_below_tolerance',
                    'uninspected': [row['route'] for row in scores],
                    'remaining_max_marginal_upper': maximum_upper}
        if not scores or scores[0]['marginal_information_value'] <= 1e-9:
            raise RuntimeError('Residual gap needs combination look-ahead, not a zero-value stop')
        choice = scores[0]
        certified = all(choice['information_efficiency_lower'] >= row['information_efficiency_upper']-1e-10
                        for row in scores[1:])
        observation = bool(survey_realization[choice['route']])
        if observation:
            known_usable |= {choice['route']}
        else:
            known_unusable |= {choice['route']}
        total_cost += choice['survey_cost']
        after = information_state(domains, known_usable, known_unusable)
        gain = after['confirmed_area']-before['confirmed_area']
        removal = before['optimistic_area']-after['optimistic_area']
        if abs(before['information_gap']-after['information_gap']-gain-removal) > 1e-3:
            raise RuntimeError('Gap decomposition failed')
        trace.append({'step': step+1, **choice, 'survey_observation': observation,
                      'ranking_certified': certified, 'realized_gain': gain, 'realized_removal': removal,
                      'confirmed_area_after': after['confirmed_area'], 'optimistic_area_after': after['optimistic_area'],
                      'information_gap_after': after['information_gap'], 'survey_total_cost': total_cost})
    raise RuntimeError('No stopping certificate')


def nonconvexity_witness(equations, solved, allowed, domain):
    candidates = []
    for entry in solved:
        if not set(entry['routes']) <= allowed:
            continue
        for answer in entry['answers']:
            candidates.append((np.asarray(answer['p'])*(1.-1e-5), entry['id']))
    best = None
    for i, (a, plan_a) in enumerate(candidates):
        for b, plan_b in candidates[i+1:]:
            if plan_a == plan_b:
                continue
            middle = (a+b)/2.
            distance = domain['outer'].distance(Point(middle))
            if best is None or distance > best['midpoint_distance_lower']:
                best = {'p_a':a.tolist(), 'p_b':b.tolist(), 'p_mid':middle.tolist(),
                        'plan_a':plan_a, 'plan_b':plan_b, 'midpoint_distance_lower':distance}
    if best is None or best['midpoint_distance_lower'] < 1.:
        raise RuntimeError('No clearly certified nonconvexity witness')
    problem = PlanningModel(equations, power=np.asarray(best['p_mid']), budget=BUDGET, threads=1)
    for route in ROUTES:
        if route not in allowed:
            problem.choices[route, 'new'].UB = 0.
    try:
        problem.solve(time_limit=30.)
        best['midpoint_full_model_status'] = problem.model.Status
        if problem.model.Status != 3:
            raise RuntimeError('Integer feasibility check did not certify midpoint infeasible')
    finally:
        problem.model.dispose()
    best['allowed'] = sorted(allowed)
    best['convex_hull_excess_area'] = domain['inner'].convex_hull.area-domain['inner'].area
    return best


def audit_results(cache):
    """独立 AC 点检查；不使用旧实验结果作为生产依赖。"""
    plans = {}
    network = Concept5()
    for result in cache.regions.values():
        for row in result['inner']:
            selection = tuple(network.encode_plan(row['choice']))
            plans.setdefault(selection, set()).update(tuple(p) for p in row['vertices'])
    ac_total = ac_feasible = 0
    for selection, vertices in plans.items():
        solver = ACPowerFlow(network.tree(np.array(selection)))
        try:
            labels = solver.classify(np.array(sorted(vertices)))
            ac_total += len(labels)
            ac_feasible += int(np.count_nonzero(labels == 1))
        finally:
            solver.close()
    return {'ac_inner_vertices': ac_total, 'ac_feasible_vertices': ac_feasible,
            'max_domain_area_uncertainty': max(row['area_gap'] for row in cache.queries)}


def run_survey(*, output=None, threads=1, time_limit=60.):
    """从网架和观测场景独立计算，目录必须是新目录。"""
    output = Path(output) if output is not None else ROOT/'results'/'concept5'/datetime.now().strftime('%Y-%m-%d_%H%M%S')
    output.mkdir(parents=True, exist_ok=False)
    network = Concept5()
    protocol = dict(schema='survey-v1', algorithm='dual_SP_with_shared_cuts', model='SOCP',
                    budget=BUDGET, cost_unit=network.cost_unit, survey_cost_unit='relative survey unit',
                    reference_max=network.power_limit, load_nodes=list(network.load_nodes),
                    network=asdict(network), network_fingerprint=network.fingerprint,
                    routes=ROUTES, synthetic_realization=SURVEY_REALIZATION,
                    survey_cost_deducted_from_budget=False,
                    stopping_area_tolerance=STOPPING_AREA_TOLERANCE,
                    domain_area_tolerance=DOMAIN_AREA_TOLERANCE, tau=1e-5,
                    threads=threads, time_limit_per_attempt=time_limit,
                    sources={name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                             for name in ('model.py', 'main.py', 'region.py', 'Network/__init__.py',
                                          'Network/concept5.py', 'survey.py')})
    write_json(output/'protocol.json', protocol)
    with threadpool_limits(limits=1):
        cache = DomainCache(output/'domains', threads=threads, time_limit=time_limit)
        started = perf_counter()
        result = simulate_surveys(cache, SURVEY_REALIZATION)
        elapsed = perf_counter()-started
        result.update(baseline_area=result['states'][0]['confirmed_area'],
                      baseline_share=result['states'][0]['confirmed_area']/result['states'][0]['optimistic_area'],
                      seconds=elapsed, query_count=len(cache), joint_cut_count=len(cache.cuts))
        allowed = frozenset(result['states'][-1]['known_usable'])
        solved = [dict(id=index, plan=row['choice'], cost=row['cost'],
                       routes=[route for route in ROUTES if row['choice'][route] is not None],
                       answers=[{'p': point} for point in row['vertices']])
                  for index, row in enumerate(cache.regions[allowed]['inner'])]
        audit_started = perf_counter()
        result['nonconvexity_witness'] = nonconvexity_witness(
            PlanningEquations(network, 'socp'), solved, allowed, cache[allowed])
        audit = audit_results(cache)
        audit.update(seconds=perf_counter()-audit_started,
                     midpoint_infeasible=result['nonconvexity_witness']['midpoint_full_model_status'] == 3,
                     all_rankings_certified=all(row['ranking_certified'] for row in result['trace']))
        audit['passed'] = (audit['ac_inner_vertices'] == audit['ac_feasible_vertices']
                           and audit['midpoint_infeasible'] and audit['all_rankings_certified']
                           and audit['max_domain_area_uncertainty'] <= DOMAIN_AREA_TOLERANCE)
    write_json(output/'results.json', result)
    write_json(output/'domains.json', [dict(allowed=sorted(allowed), inner=mapping(value['inner']),
                                           outer=mapping(value['outer']), plan_ids=None)
                                      for allowed, value in cache.items()])
    write_json(output/'cuts.json', cache.cuts)
    write_json(output/'schemes.json', solved)
    write_json(output/'audit.json', audit)
    for name in ('trace', 'all_candidates'):
        rows = result[name]
        with (output/(name+'.csv')).open('w', encoding='utf-8-sig', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    lines = ['# 五节点道路勘察结果', '',
             f'计算耗时 {elapsed:.3f} 秒；按需查询 {len(cache)} 个道路集合，共享 {len(cache.cuts)} 条物理联合割。',
             f'观测顺序：{", ".join(row["route"]+("+" if row["survey_observation"] else "−") for row in result["trace"])}；未勘察：{", ".join(result["uninspected"])}。', '',
             '|步数|确认域 / kW²|乐观域 / kW²|信息间隙上界 / kW²|', '|---|---:|---:|---:|']
    lines.extend(f'|{row["step"]}|{row["confirmed_area"]:.3f}|{row["optimistic_area"]:.3f}|{row["information_gap_upper"]:.3f}|'
                 for row in result['states'])
    lines.extend(['', f'独立 AC 顶点检查：{audit["ac_feasible_vertices"]}/{audit["ac_inner_vertices"]}；这不是连续全域的 AC 认证。',
                  f'核查通过：{audit["passed"]}。参数与源文件指纹见 protocol.json，逐轮评分见 all_candidates.csv。'])
    (output/'report.md').write_text('\n'.join(lines)+'\n', encoding='utf-8')
    print(json.dumps({'output': str(output), 'seconds': elapsed, 'passed': audit['passed']}, ensure_ascii=False), flush=True)
    if not audit['passed']:
        raise RuntimeError('Survey validation failed; see audit.json')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--time-limit', type=float, default=60.)
    args = parser.parse_args()
    run_survey(output=args.output, threads=args.threads, time_limit=args.time_limit)


if __name__ == '__main__':
    main()
