"""独立 FourBus SOCP 外域实验；不调用或修改 main.py。

运行：python experiments/fourbus_outer.py --budget 20000

G = max_{(x,p) in outer} min_{y,eta>=0} eta。
用现有物理模型的完整锥对偶消去内层 min；p 与对偶乘子的乘积
使 G 成为非凸混合整数二次约束问题，交给 Gurobi 的全局算法。
SP 保留现有实现：固定 G 返回的 x,p，认证或返回全局联合割。
只有 G 的全局上界 <= epsilon 才停止认证；超时目标值不是上界。
"""
import argparse
import json
from pathlib import Path
import sys
from time import perf_counter
from types import SimpleNamespace

import gurobipy as gp
from gurobipy import GRB
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from Network.four_bus_five_corridor import FourBus
from model import GridPhysics, MasterProblem, PLANNING_TOL, SubProblem, new_model


def export_sp(equations):
    """只导出现有 SP 的线性行、状态界和锥；p 为 kW，y 为原始状态。"""
    with new_model('SP_coefficients', 1) as model:
        x = model.addVars(equations.keys, ub=1., name='x')
        p = model.addVars(equations.network.load_nodes, lb=-GRB.INFINITY, name='p_kw')
        eta = model.addVar(name='violation')
        operation = equations.add_operation(model, x, p, eta)
        model.update()
        rows = model.getConstrs()
        coefficients = model.getA().toarray()
        sign = np.array([1. if row.Sense == '<' else -1. for row in rows])
        assert all(row.Sense in ('<', '>') for row in rows)
        coefficients *= sign[:, None]
        state_ids = [var.index for var in operation.state.tolist()]
        size = len(state_ids)
        sp_y = np.vstack([coefficients[:, state_ids], np.eye(size), -np.eye(size)])
        sp_eta = np.r_[coefficients[:, eta.index], np.zeros(2*size)]
        sp_x = np.vstack([coefficients[:, [var.index for var in x.values()]],
                          np.zeros((2*size, len(x)))])
        sp_p = np.vstack([coefficients[:, [var.index for var in p.values()]],
                          np.zeros((2*size, len(p)))])
        sp_rhs = np.r_[sign*np.array(model.getAttr('RHS', rows)),
                       equations.y_ub_global, -equations.y_lb_global]
        cone_y, cone_constant, cone_slices = [], [], []
        position = {index: j for j, index in enumerate(state_ids)}
        for head, tail in operation.cones:
            start = len(cone_constant)
            for expression in [head, *tail]:
                row = np.zeros(size)
                for j in range(expression.size()):
                    row[position[expression.getVar(j).index]] += expression.getCoeff(j)
                cone_y.append(row)
                cone_constant.append(expression.getConstant())
            cone_slices.append(slice(start, len(cone_constant)))
    return SimpleNamespace(sp_y=sp_y, sp_eta=sp_eta, sp_x=sp_x, sp_p=sp_p,
                           sp_rhs=sp_rhs, cone_y=np.array(cone_y),
                           cone_constant=np.array(cone_constant), cone_slices=cone_slices)


class GlobalViolation:
    """全局搜索所有合法 x,p；完整 SP 对偶的可行域不依赖已知方案列表。"""

    def __init__(self, equations, budget, axis_bounds, total_bound, *, threads=4):
        self.equations = equations
        self.problem = MasterProblem(equations, budget=budget, cuts_only=True, threads=threads)
        self.model = model = self.problem.model
        model.ModelName = 'G_max_min_eta'
        model.Params.NonConvex = 2
        model.Params.Aggregate = 0
        model.Params.ScaleFlag = 0
        model.Params.MIPGapAbs = 1e-7
        self.problem.power.UB = axis_bounds
        model.addConstr(self.problem.power.sum() <= total_bound, name='MP2_total_bound')
        template = export_sp(equations)
        self.template = template

        # 1. <= 行乘子非负；eta >= 0 对应 sum(放松等式两侧乘子) <= 1。
        upper = np.where(template.sp_eta < 0., 1., GRB.INFINITY)
        self.row_dual = model.addMVar(len(template.sp_rhs), ub=upper, name='row_dual')
        model.addConstr(-template.sp_eta @ self.row_dual <= 1., name='eta_stationarity')

        # 2. Lorentz 锥自对偶；头 >= 0、尾自由，保留每个完整锥。
        lower = np.full(len(template.cone_constant), -GRB.INFINITY)
        for block in template.cone_slices:
            lower[block.start] = 0.
        self.cone_dual = model.addMVar(len(lower), lb=lower, name='cone_dual')
        for j, block in enumerate(template.cone_slices):
            head = self.cone_dual[block.start].item()
            tail = self.cone_dual[block.start+1:block.stop].tolist()
            model.addQConstr(gp.quicksum(item*item for item in tail) <= head*head,
                             name=f'dual_cone[{j}]')
        model.addConstr(template.sp_y.T @ self.row_dual == template.cone_y.T @ self.cone_dual,
                        name='y_stationarity')

        # 3. x*乘子用指示约束精确表达；不设置任意乘子上界或枚举建设组合。
        self.choice_term = model.addMVar(equations.network.n_types,
                                        lb=-GRB.INFINITY, name='choice_term')
        delta = template.sp_x.T @ self.row_dual
        for j, choice in enumerate(self.problem.choices.values()):
            term = self.choice_term[j].item()
            model.addGenConstrIndicator(choice, 0, term == 0.)
            model.addGenConstrIndicator(choice, 1, term == delta[j].item())

        # 4. 目标是内层 min 的对偶值；唯一连续双线性项为 p*平衡乘子。
        dual_value = (self.choice_term.sum()-template.sp_rhs @ self.row_dual
                      -template.cone_constant @ self.cone_dual).item()
        loads = list(self.problem.loads.values())
        for i, j in zip(*np.nonzero(template.sp_p)):
            dual_value += template.sp_p[i, j]*self.row_dual[i].item()*loads[j]
        # FourBus 的零潮流、v=1 提供 eta* <= max(p)/base 的解析上界。
        self.violation = model.addVar(ub=float(max(axis_bounds)/equations.network.base), name='G_violation')
        model.addQConstr(self.violation <= dual_value, name='dual_objective')
        model.setObjective(self.violation, GRB.MAXIMIZE)

    def candidate_lower(self, x, power):
        """用锥可行乘子和有限 y 盒计算保守下界，不把平方残差当作范数残差。"""
        template = self.template
        row_dual = np.maximum(self.row_dual.X, 0.)
        cone_dual = self.cone_dual.X.copy()
        for block in template.cone_slices:
            cone_dual[block.start] = max(cone_dual[block.start],
                                        np.linalg.norm(cone_dual[block.start+1:block.stop]))
        scale = max(1., -template.sp_eta@row_dual)
        row_dual, cone_dual = row_dual/scale, cone_dual/scale
        h = template.sp_y.T@row_dual-template.cone_y.T@cone_dual
        value = (row_dual@(template.sp_x@x+template.sp_p@power-template.sp_rhs)
                 -cone_dual@template.cone_constant
                 +np.maximum(h, 0.)@self.equations.y_lb_global
                 +np.minimum(h, 0.)@self.equations.y_ub_global)
        return max(0., float(value))

    def solve(self, time_limit, epsilon):
        """返回 G 的候选及上下界；只有 bound 可以给出全局停止证书。"""
        model = self.model
        model.Params.TimeLimit = time_limit
        model.Params.BestBdStop = epsilon
        model.optimize()
        if model.Status not in (GRB.OPTIMAL, GRB.TIME_LIMIT, GRB.USER_OBJ_LIMIT):
            raise RuntimeError(f'G status={model.Status}')
        if model.SolCount and model.MaxVio > PLANNING_TOL:
            raise RuntimeError(f'G MaxVio={model.MaxVio:g} > {PLANNING_TOL:g}')
        bound = min(float(model.ObjBound), self.violation.UB)
        x = np.rint(self.problem.x.X).astype(int) if model.SolCount else None
        power = self.problem.power.X if model.SolCount else None
        answer = dict(objective=self.candidate_lower(x, power) if model.SolCount else None,
                      bound=bound, status=int(model.Status), seconds=float(model.Runtime),
                      x=x, p=power)
        # 外域只会收缩，上轮有效上界继续成立；不让超时产生的较松界覆盖已有证书。
        self.violation.UB = bound
        return answer


def run_experiment(budget=20000., epsilon=1e-4, iterations=100, time_limit=30., threads=4,
                   output='results/fourbus_outer', cuts=()):
    """五步主循环在此处，方便逐行调试；结果与主线输出分开。"""
    started = perf_counter()
    equations = GridPhysics(FourBus(), 'socp')
    net = equations.network
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)

    # 1. 完整 MP2 初始化：三个负荷方向 + 总负荷，x/p/y 都允许变化。
    initial = []
    for direction in [*np.eye(len(net.load_nodes)), np.ones(len(net.load_nodes))]:
        problem = MasterProblem(equations, budget=budget, direction=direction, threads=threads)
        with problem.model:
            answer = problem.solve(time_limit=60.)
        initial.append(dict(direction=direction.tolist(), x=answer['x'].tolist(),
                            p=answer['p'].tolist(), bound=float(answer['bound'])))
        print(f"MP2 direction={direction.tolist()} bound={answer['bound']:.6f} kW", flush=True)
    axis_bounds = np.array([item['bound'] for item in initial[:-1]])
    total_bound = initial[-1]['bound']
    result = dict(budget=budget, epsilon=epsilon, axis_bounds=axis_bounds.tolist(),
                  total_bound=total_bound, initial=initial, trace=[],
                  cuts=[np.asarray(cut).tolist() for cut in cuts], initial_cut_count=len(cuts),
                  status='running', certified=False)
    sp = SubProblem(equations, threads=threads)
    search = GlobalViolation(equations, budget, axis_bounds, total_bound, threads=threads)
    with search.model:
        for cut in cuts:
            search.problem.add_cut(np.asarray(cut))
        for step in range(1, iterations+1):
            # 2. G 同时选 x,p：对所有预算内网架最大化“各自最优运行后的违反量”。
            answer = search.solve(time_limit, epsilon)
            record = dict(step=step, objective=answer['objective'], bound=answer['bound'],
                          status=answer['status'], seconds=answer['seconds'])
            result['trace'].append(record)
            print(f"G {step}: objective={answer['objective']} upper={answer['bound']:.8g}", flush=True)

            # 3. 只有全局上界达标才认证整个联合外域；无反例但有间隙则尚未证明。
            if answer['bound'] <= epsilon:
                result.update(status='residual_certified', certified=True)
                break
            if answer['objective'] is None or answer['objective'] <= PLANNING_TOL:
                result['status'] = 'global_search_unresolved'
                break
            x, power = answer['x'], answer['p']
            record.update(x=x.tolist(), p=power.tolist(), plan=net.decode_plan(x))

            # 4. 固定 G 返回的同一对 x,p，SP 只调整 y；失败仅排除这对组合。
            checked = sp.solve(x, power)
            record['sp_feasible'] = checked['feasible']
            if checked['feasible']:
                raise RuntimeError('Positive G violation disagrees with the original SP certificate')
            cut = checked['cut']
            record['cut_at_candidate'] = float(cut[0]+cut[1:4]@power+cut[4:]@x)
            print(f"  SP cut={record['cut_at_candidate']:.8g} p={np.round(power, 4).tolist()} plan={record['plan']}", flush=True)

            # 5. 加入含 x,p 的全局有效割，再回到步骤 2，允许重新选网架和负荷。
            search.problem.add_cut(cut)
            result['cuts'].append(cut.tolist())
            result['seconds'] = perf_counter()-started
            (output/'results.json').write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding='utf-8')
        else:
            result['status'] = 'iteration_limit'
        result['seconds'] = perf_counter()-started
        (output/'results.json').write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding='utf-8')
    # 保存可继续查询的外域本身；不包含 G 的对偶变量，也不混入完整电气模型。
    outer = MasterProblem(equations, budget=budget, cuts_only=True,
                          cuts=map(np.array, result['cuts']), threads=threads)
    with outer.model:
        outer.power.UB = axis_bounds
        outer.model.addConstr(outer.power.sum() <= total_bound, name='MP2_total_bound')
        outer.model.setObjective(0.)
        outer.model.write(str(output/'outer.lp'))
    print(f"{result['status']}: {len(result['cuts'])} cuts, {result['seconds']:.2f} s -> {output/'results.json'}", flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--budget', type=float, default=20000.)
    parser.add_argument('--epsilon', type=float, default=1e-4)
    parser.add_argument('--iterations', type=int, default=100)
    parser.add_argument('--time-limit', type=float, default=30., help='每次 G 的全局搜索时限，秒')
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--output', default='results/fourbus_outer')
    parser.add_argument('--cuts-from', help='同一 FourBus 模型已有的 results.json，仅复用有效割')
    arguments = vars(parser.parse_args())
    cuts_from = arguments.pop('cuts_from')
    if cuts_from:
        arguments['cuts'] = json.loads(Path(cuts_from).read_text(encoding='utf-8'))['cuts']
    run_experiment(**arguments)
