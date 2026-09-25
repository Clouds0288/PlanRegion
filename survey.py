"""五节点道路勘察：共享对偶割、条件信息价值与停止证书。"""
from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime
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


def polygon_union(rows):
    return unary_union([MultiPoint(row['vertices']).convex_hull for row in rows]) if rows else GeometryCollection()


class DomainCache(dict):
    """按道路集合缓存认证域，共享物理联合割；计算期间不写文件。"""

    def __init__(self, *, threads=1, time_limit=60., share_cuts=True):
        super().__init__()
        self.threads, self.time_limit, self.share_cuts = threads, time_limit, share_cuts
        self.cuts, self.regions = [], {}

    def __missing__(self, allowed):
        # 1. 设置本次允许的道路，继承同一物理模型的有效割
        allowed = frozenset(allowed)
        network = Concept5(allowed)
        local_cuts = list(self.cuts) if self.share_cuts else []

        def progress(event, **data):
            if event == 'cut':
                local_cuts.append(np.asarray(data['cut']).copy())

        # 2. 对偶构域；light 未达到认证精度时用 physical 继续，保留有效割
        for mode in ('light', 'physical'):
            result = build_continuous_region(
                network, 'socp', BUDGET, np.array([135., 135.]), tau=1e-5,
                residual_mode=mode, time_limit=self.time_limit, threads=self.threads,
                progress=progress, cuts=local_cuts,
            )
            inner, outer = (polygon_union(result[key]) for key in ('inner', 'outer'))
            if result['status'] == 'certified' and outer.area-inner.area <= DOMAIN_AREA_TOLERANCE:
                break
        else:
            raise RuntimeError(f'Roads {sorted(allowed)}: {result["status"]}, area gap {outer.area-inner.area:g}')

        # 3. 用已认证的子集/超集维持道路域包含关系，不修改此前评分的域
        smaller = [value for key, value in self.items() if key < allowed]
        larger = [value for key, value in self.items() if allowed < key]
        inner = unary_union([inner, *(value['inner'] for value in smaller)])
        for value in larger:
            inner = inner.intersection(value['inner'])
            outer = outer.intersection(value['outer'])
        outer = unary_union([outer, *(value['outer'] for value in smaller)])
        if outer.area-inner.area > DOMAIN_AREA_TOLERANCE:
            raise RuntimeError(f'Roads {sorted(allowed)}: area gap {outer.area-inner.area:g}')

        # 4. 缓存本次域、构域结果及共享割
        if self.share_cuts:
            self.cuts = local_cuts
        self.regions[allowed] = result
        domain = {'inner': inner, 'outer': outer, 'plan_ids': None}
        self[allowed] = domain
        return domain


def increase_bounds(before, after):
    """返回新增面积的估计、下界和上界，单位 kW²。"""
    if before['plan_ids'] is not None and before['plan_ids'] == after['plan_ids']:
        return (0., 0., 0.)
    return (after['inner'].difference(before['inner']).area,
            after['inner'].difference(before['outer']).area,
            after['outer'].difference(before['inner']).area)


def candidate_values(domains, known_usable, known_unusable):
    """计算两种观测下的条件价值，按单位勘察费的信息价值排序。"""
    # 1. 对每条未知道路分别计算确认域扩展和乐观域排除
    allowed = frozenset(ROUTES)-known_unusable
    scores = []
    for route, data in ROUTES.items():
        if route in known_usable or route in known_unusable:
            continue
        gain = increase_bounds(domains[known_usable], domains[known_usable | {route}])
        removal = increase_bounds(domains[allowed-{route}], domains[allowed])
        probability, cost = data['road_usable_probability'], data['survey_cost']
        row = {'route': route, 'road_usable_probability': probability, 'survey_cost': cost}
        # 2. 将面积区间传播到期望价值及单位费用价值
        for index, suffix in enumerate(('', '_lower', '_upper')):
            row['gain_if_usable'+suffix], row['removal_if_unusable'+suffix] = gain[index], removal[index]
            row['expected_expansion'+suffix] = probability*gain[index]
            row['expected_exclusion'+suffix] = (1-probability)*removal[index]
            row['marginal_information_value'+suffix] = probability*gain[index]+(1-probability)*removal[index]
            row['information_efficiency'+suffix] = row['marginal_information_value'+suffix]/cost
        scores.append(row)
    # 3. 同值时按道路名称确定顺序
    return sorted(scores, key=lambda row: (-row['information_efficiency'], row['route']))


def information_state(domains, known_usable, known_unusable):
    """由道路信息生成确认/乐观域及面积区间，坐标 kW、面积 kW²。"""
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
    """逐轮评分、观测、更新；仅记录决策所需的状态和候选评分。"""
    known_usable, known_unusable = frozenset(), frozenset()
    states, trace, all_candidates, total_cost = [], [], [], 0.
    for step in range(len(ROUTES)+1):
        # 1. 计算当前信息上下域及全部未知路线的条件价值
        before = information_state(domains, known_usable, known_unusable)
        before.update(step=step, survey_total_cost=total_cost)
        states.append(before)
        scores = candidate_values(domains, known_usable, known_unusable)
        all_candidates.extend({'step': step, **row} for row in scores)
        # 2. 同时检查单路边际值和整个剩余信息间隙的停止证书
        maximum_upper = max((row['marginal_information_value_upper'] for row in scores), default=0.)
        if maximum_upper <= STOPPING_AREA_TOLERANCE and before['information_gap_upper'] <= STOPPING_AREA_TOLERANCE:
            return {'states': states, 'trace': trace, 'all_candidates': all_candidates,
                    'stop_reason': 'residual_domain_gap_and_all_marginals_below_tolerance',
                    'uninspected': [row['route'] for row in scores],
                    'remaining_max_marginal_upper': maximum_upper}
        if not scores or scores[0]['marginal_information_value'] <= 1e-9:
            raise RuntimeError('Residual gap needs combination look-ahead, not a zero-value stop')
        # 3. 选择信息效率最高的道路，再读取该道路的观测
        choice = scores[0]
        certified = all(choice['information_efficiency_lower'] >= row['information_efficiency_upper']-1e-10
                        for row in scores[1:])
        observation = bool(survey_realization[choice['route']])
        if observation:
            known_usable |= {choice['route']}
        else:
            known_unusable |= {choice['route']}
        # 4. 更新勘察费及实际新增/排除面积，保存本轮决策
        total_cost += choice['survey_cost']
        after = information_state(domains, known_usable, known_unusable)
        gain = after['confirmed_area']-before['confirmed_area']
        removal = before['optimistic_area']-after['optimistic_area']
        trace.append({'step': step+1, **choice, 'survey_observation': observation,
                      'ranking_certified': certified, 'realized_gain': gain, 'realized_removal': removal,
                      'confirmed_area_after': after['confirmed_area'], 'optimistic_area_after': after['optimistic_area'],
                      'information_gap_after': after['information_gap'], 'survey_total_cost': total_cost})
    raise RuntimeError('No stopping certificate')


def nonconvexity_witness(equations, solved, allowed, domain):
    """寻找两个可行端点及其不可行中点；所有负荷坐标为 kW。"""
    # 1. 收集不同建设方案的已认证端点
    candidates = []
    for entry in solved:
        if not set(entry['routes']) <= allowed:
            continue
        for answer in entry['answers']:
            candidates.append((np.asarray(answer['p'])*(1.-1e-5), entry['id']))
    # 2. 选择距数值外域最远的中点
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
    # 3. 用完整整数模型独立证明中点不可行
    problem = PlanningModel(equations, power=np.asarray(best['p_mid']), budget=BUDGET, threads=1)
    for route in ROUTES:
        if route not in allowed:
            problem.choices[route, 'new'].UB = 0.
    with problem.model:
        problem.solve(time_limit=30.)
        best['midpoint_full_model_status'] = problem.model.Status
        if problem.model.Status != 3:
            raise RuntimeError('Integer feasibility check did not certify midpoint infeasible')
    best['allowed'] = sorted(allowed)
    best['convex_hull_excess_area'] = domain['inner'].convex_hull.area-domain['inner'].area
    return best


def audit_results(cache):
    """独立 AC 点检查；不使用旧实验结果作为生产依赖。"""
    # 1. 按完整选型合并重复顶点
    plans = {}
    network = Concept5()
    for result in cache.regions.values():
        for row in result['inner']:
            selection = tuple(network.encode_plan(row['choice']))
            plans.setdefault(selection, set()).update(tuple(p) for p in row['vertices'])
    # 2. 对每个方案的认证顶点执行独立 AC 校核
    ac_total = ac_feasible = 0
    for selection, vertices in plans.items():
        solver = ACPowerFlow(network.tree(np.array(selection)))
        try:
            labels = solver.classify(np.array(sorted(vertices)))
            ac_total += len(labels)
            ac_feasible += int(np.count_nonzero(labels == 1))
        finally:
            solver.close()
    # 3. 汇总已检验点及最大数值面积区间
    return {'ac_inner_vertices': ac_total, 'ac_feasible_vertices': ac_feasible,
            'max_domain_area_uncertainty': max(row['outer'].area-row['inner'].area for row in cache.values())}


def run_survey(*, output=None, threads=1, time_limit=60.):
    """计算、核查、保存；全部科研结果集中在一个 results.json。"""
    # 1. 设置案例与复现参数
    output = Path(output) if output is not None else ROOT/'results'/'concept5'/datetime.now().strftime('%Y-%m-%d_%H%M%S')
    output.mkdir(parents=True, exist_ok=True)
    network = Concept5()
    protocol = dict(schema='survey-v2', algorithm='dual_SP_with_shared_cuts', model='SOCP',
                    budget=BUDGET, cost_unit=network.cost_unit, survey_cost_unit='relative survey unit',
                    reference_max=network.power_limit, load_nodes=list(network.load_nodes),
                    network=asdict(network), network_fingerprint=network.fingerprint,
                    routes=ROUTES, synthetic_realization=SURVEY_REALIZATION,
                    survey_cost_deducted_from_budget=False,
                    stopping_area_tolerance=STOPPING_AREA_TOLERANCE,
                    domain_area_tolerance=DOMAIN_AREA_TOLERANCE, tau=1e-5,
                    threads=threads, time_limit_per_attempt=time_limit)
    with threadpool_limits(limits=1):
        # 2. 按需求域并完成逐轮勘察；缓存及割只驻留内存
        cache = DomainCache(threads=threads, time_limit=time_limit)
        started = perf_counter()
        result = simulate_surveys(cache, SURVEY_REALIZATION)
        elapsed = perf_counter()-started
        result.update(baseline_area=result['states'][0]['confirmed_area'],
                      baseline_share=result['states'][0]['confirmed_area']/result['states'][0]['optimistic_area'],
                      seconds=elapsed, query_count=len(cache), joint_cut_count=len(cache.cuts))
        # 3. 整理最终方案，生成非凸见证并独立核查
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
    # 4. 统一保存一次；绘图只需读取这个文件
    result.update(protocol=protocol, schemes=solved, audit=audit)
    (output/'results.json').write_text(
        json.dumps(json_value(result), ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
    print(f'勘察完成：{elapsed:.3f}s，{len(cache)} 个道路集合，{len(cache.cuts)} 条共享割；结果：{output}', flush=True)
    if not audit['passed']:
        raise RuntimeError('Survey validation failed; see results.json audit')
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
