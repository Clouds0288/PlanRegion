"""独立正负功率实验的连续 SP；原物理装配与 Gurobi 取割 LP 保持一致。"""
from time import perf_counter

import clarabel
import gurobipy as gp
from gurobipy import GRB
import numpy as np
from scipy import sparse

from model import SubProblem, new_model, PLANNING_TOL, SP_TIME_LIMIT


def solve_conic(model, equations, operation, x, fixed, threads, deadline, tolerance, margin):
    """原 SOCP 的固定变量消元、连续求解及原始残差核验。"""
    model.update()
    variables, rows = model.getVars(), model.getConstrs()
    n = len(variables)

    def affine(expression):
        expression = gp.LinExpr(expression)
        coefficients = np.zeros(n)
        for j in range(expression.size()):
            coefficients[expression.getVar(j).index] += expression.getCoeff(j)
        return coefficients, expression.getConstant()

    # 1. 原始等式、不等式和变量界，统一写成 Az+s=b。
    matrix = model.getA().tocsc()
    rhs = np.asarray(model.getAttr('RHS', rows))
    senses = np.asarray(model.getAttr('Sense', rows))
    lb, ub = np.array(model.getAttr('LB', variables)), np.array(model.getAttr('UB', variables))
    for key, selected in zip(equations.keys, x):
        if not selected:
            fixed.update({flow[key].index: 0. for flow in (operation.P, operation.Q, operation.ell)})
    for e in equations.types:
        if sum(x[j] for j, key in enumerate(equations.keys) if key[0] == e):
            fixed.update({operation.plus[e].index: 0., operation.minus[e].index: 0.})
    free = np.array([j for j in range(n) if j not in fixed])
    solution = np.zeros(n)
    solution[list(fixed)] = list(fixed.values())
    reduced = matrix[:, free]
    keep = np.asarray(abs(reduced).sum(axis=1)).ravel() != 0.
    a, b = reduced[keep], (rhs-matrix@solution)[keep]
    equal = senses[keep] == '='
    direction = np.where(senses[keep] == '>', -1., 1.)
    eye = sparse.eye(len(free), format='csc')
    lower, upper = lb[free] > -GRB.INFINITY, ub[free] < GRB.INFINITY
    blocks = [a[equal], sparse.diags(direction[~equal])@a[~equal], -eye[lower], eye[upper]]
    values = [b[equal], direction[~equal]*b[~equal], -lb[free][lower], ub[free][upper]]
    cones = [clarabel.ZeroConeT(int(equal.sum())),
             clarabel.NonnegativeConeT(int((~equal).sum()+lower.sum()+upper.sum()))]

    # 2. 每个原二阶锥保存为仿射向量，直接交给连续锥求解器。
    cone_rows = []
    for j, (head, tail) in enumerate(operation.cones):
        expressions = [affine(item) for item in (head, *tail)]
        a, b = np.array([item[0] for item in expressions]), np.array([item[1] for item in expressions])
        cone_rows.append((a, b))
        if j < len(x) and not x[j]:
            continue  # 开断型号 P=Q=ell=0，原电流锥恒成立。
        if j < len(x) and equations.network.name == 'case33bw':
            # Case33 小电流锥等价缩放，避免 v 与 ell 数量级悬殊。
            key = equations.keys[j]
            v = operation.v[equations.ends[key[0]][0]]
            ell, P, Q = operation.ell[key], operation.P[key], operation.Q[key]
            expressions = [affine(item) for item in (v+100*ell, 20*P, 20*Q, v-100*ell)]
            a, b = np.array([item[0] for item in expressions]), np.array([item[1] for item in expressions])
        blocks.append(sparse.csc_matrix(-a[:, free]))
        constant = b+a@solution
        constant[0] -= margin*(100. if j < len(x) and equations.network.name == 'case33bw' else 1.)
        values.append(constant)
        cones.append(clarabel.SecondOrderConeT(len(b)))
    objective, _ = affine(model.getObjective())
    objective *= model.ModelSense
    # 射线的线性等式先精确消元，避免小阻抗压降行在锥求解中损失精度。
    if equal.any():
        left, singular, right = np.linalg.svd(blocks.pop(0).toarray(), full_matrices=True)
        rank = int((singular > singular[0]*max(len(free), int(equal.sum()))*np.finfo(float).eps).sum())
        offset = right[:rank].T@((left[:, :rank].T@values.pop(0))/singular[:rank])
        basis = right[rank:].T
        cones.pop(0)
        matrix_free = sparse.vstack(blocks, format='csc')
        conic_matrix = sparse.csc_matrix(matrix_free@basis)
        conic_rhs = np.concatenate(values)-matrix_free@offset
        conic_objective = basis.T@objective[free]
    else:
        conic_matrix = sparse.vstack(blocks, format='csc')
        conic_rhs, conic_objective = np.concatenate(values), objective[free]
    settings = clarabel.DefaultSettings()
    settings.verbose = False
    settings.max_threads = threads
    settings.static_regularization_constant = 1e-10
    settings.tol_gap_abs = settings.tol_gap_rel = tolerance
    settings.tol_feas = tolerance
    settings.time_limit = max(0., deadline-perf_counter())
    count = len(conic_objective)
    answer = clarabel.DefaultSolver(sparse.csc_matrix((count, count)), conic_objective,
        conic_matrix, conic_rhs, cones, settings).solve()
    if answer.status == clarabel.SolverStatus.MaxTime:
        raise TimeoutError('连续子问题达到总时限')
    if answer.status not in (clarabel.SolverStatus.Solved, clarabel.SolverStatus.AlmostSolved):
        raise RuntimeError(f'Continuous SOCP: {answer.status}, x={list(x)}')
    solution[free] = offset+basis@answer.x if equal.any() else answer.x

    # 3. 以原始约束重新核验；锥求解器的缩放残差不直接作为物理证书。
    residual = matrix@solution-rhs
    violation = max(0., np.max(np.where(senses == '=', np.abs(residual),
                    np.where(senses == '>', -residual, residual))),
                    np.max(lb-solution), np.max(solution-ub))
    for constraint in model.getQConstrs():
        expression = model.getQCRow(constraint)
        a, b = affine(expression.getLinExpr())
        value = a@solution+b+sum(expression.getCoeff(j)*solution[expression.getVar1(j).index]
            *solution[expression.getVar2(j).index] for j in range(expression.size()))
        violation = max(violation, value-constraint.QCRHS)
    normals = []
    for a, b in cone_rows:
        tail = (a@solution+b)[1:]
        length = np.linalg.norm(tail)
        normals.append(tail/length if length else tail)
    return solution, violation, normals


class PortSubProblem(SubProblem):
    def solve(self, x, power, time_limit=None, *, score_only=False):
        self.calls += 1
        deadline = perf_counter()+(SP_TIME_LIMIT[self.equations.method] if time_limit is None else time_limit)
        with new_model('port_SP', self.threads) as model:
            choice, p, eta, operation, _ = self._build(model, x, power)
            model.update()
            fixed = {v.index: float(a) for v, a in zip([*choice.values(), *p.values()], [*x, *power])}
            solution, violation, normals = solve_conic(model, self.equations, operation, x,
                fixed, self.threads, deadline, 1e-10, 0.)
            value = float(solution[eta.index])
            if max(0., value)+violation <= PLANNING_TOL:
                state = solution[[v.index for v in operation.state.tolist()]]
                return dict(cut=None, state=state, feasible=True, eta=value)
            if value <= PLANNING_TOL:
                raise RuntimeError(f'SP: eta={value:g}, MaxVio={violation:g}; certificate exceeds tolerance')
        if score_only:
            return dict(cut=None, state=None, feasible=False, eta=value, cone_normals=normals)
        cut = self.generate_cut(x, power, normals, time_limit=deadline-perf_counter())
        return dict(cut=cut, state=None, feasible=False, eta=value)
