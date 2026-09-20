"""原线性 Benders 规划流程；从 NodePower Notebook 提取，供独立重复计时。"""
from itertools import combinations, product
from math import gcd
from pathlib import Path
from time import sleep

import numpy as np
from gurobipy import GRB
from IPython.display import clear_output

from Network.four_bus_five_corridor import network
from model import MP1, MP2, SP, TOL, add_cut
from region import simplex, update_outer, add_certificate, ResidualSearch
from plot import show_benders
from replay import Replay


def advance(solver, steps, pause, draw):
    for step in range(steps):
        if solver.finished:
            break
        solver.step()
        if draw is not None:
            clear_output(wait=True)
            draw(solver)
        if pause and not solver.finished and step + 1 < steps:
            sleep(pause)
    return solver.result


class MP1SP:
    """给定节点负荷求最低费用；total=True 只限定总量，允许自由分配。"""

    def __init__(self, power, net=network, *, total=False, fixed_x=None, log=None, sp=None):
        mode = "vertex" if fixed_x is not None else "MP1-total" if total else "MP1"
        self._initialize(net, mode, power, log, sp)
        self.master, self.x, self.p = MP1(net, power, self.cuts, total=total, fixed_x=fixed_x)

    def _initialize(self, net, mode, target, log, sp):
        self.network, self.target = net, target
        self.log = {"queries": [], "history": []} if log is None else log
        self.sp = SP(net) if sp is None else sp
        self.query_id = len(self.queries)
        query = {"mode": mode}
        if mode != "vertex":
            query["target"] = None if np.ndim(target) == 0 and np.isinf(target) else target
        self.queries.append(query)
        self.finished, self.feasible = False, None

    @property
    def queries(self):
        return self.log["queries"]

    @property
    def history(self):
        return self.log["history"]

    @property
    def cuts(self):
        return [r["cut"] for r in self.history if r["cut"] is not None]

    def step(self):
        self.master.optimize()
        if self.master.Status == GRB.INFEASIBLE:
            row = dict(query=self.query_id, x=None, p=None, eta=None, cut=None)
            self.finished, self.feasible = True, False
        else:
            x = np.rint([v.X for v in self.x.values()])
            p = np.array([v.X for v in self.p.values()])
            eta, cut = self.sp.solve(x, p)
            row = dict(query=self.query_id, x=x, p=p, eta=eta, cut=cut)
            self.finished = self.feasible = eta <= TOL
            if cut is not None:
                add_cut(self.master, self.x, self.p, cut)
        self.history.append(row)
        if self.finished:
            self.master.dispose()
        return row

    @property
    def result(self):
        if not self.finished:
            return {"status": "迭代中"}
        if not self.feasible:
            return {"status": "不可行"}
        row = self.history[-1]
        return dict(status="最优", cost=float(self.network.cost@row["x"]),
                    power=row["p"], total=float(sum(row["p"])), x=row["x"])

    def run(self, steps=10, pause=1, show=True):
        return advance(self, steps, pause, show_benders if show else None)


class MP2SP(MP1SP):
    """给定预算求最大总负荷，节点负荷自由分配。"""

    def __init__(self, budget=np.inf, net=network, *, log=None, sp=None):
        self._initialize(net, "MP2", budget, log, sp)
        self.master, self.x, self.p = MP2(net, budget, self.cuts)


class Planning:
    """恢复费用前沿；优先认证尚未被全局并集覆盖的局部凸区域。"""

    def __init__(self, budget=np.inf, net=network):
        self.network, self.budget = net, budget
        self.log = {"queries": [], "history": []}
        self.sp = SP(net)
        self.quantum = gcd(*np.rint(net.cost[net.cost > 0]).astype(int))
        designs = []
        ne, nk = len(net.corridors), len(net.lines)
        for edges in combinations(range(ne), len(net.nodes)-1):
            reached = {net.root}
            for _ in net.load_nodes:
                for e in edges:
                    edge = net.corridors[e]
                    if edge.start in reached or edge.end in reached:
                        reached.update((edge.start, edge.end))
            if len(reached) != len(net.nodes):
                continue
            for types in product(range(nk), repeat=len(edges)):
                x = np.zeros(ne*nk)
                for e, k in zip(edges, types):
                    x[e*nk+k] = 1
                if net.cost@x <= budget:
                    designs.append(x)
        self.designs = np.array(sorted(designs, key=lambda x: (net.cost@x, tuple(x)))).reshape(-1, ne*nk)
        self.costs = self.designs@net.cost
        self.polytopes = [simplex(net) for _ in self.designs]
        self.certified, self.search = {}, None
        self.session = MP2SP(budget, net, log=self.log, sp=self.sp)
        self.phase, self.finished = "MP2", False
        self.replay = None

    @property
    def history(self):
        return self.log["history"]

    @property
    def queries(self):
        return self.log["queries"]

    @property
    def frontier(self):
        return [(self.queries[r["query"]]["target"], float(self.network.cost@r["x"]))
                for r in self.history if self.queries[r["query"]]["mode"] == "MP1-total"
                and r["x"] is not None and r["cut"] is None]

    def _next_vertex(self):
        self.phase = "vertex"
        if self.search is None:
            self.search = ResidualSearch(self.designs, self.polytopes)
        self.candidate = self.search.next(self.polytopes, self.certified)
        self.finished = self.candidate is None

    def step(self):
        if self.phase == "vertex":
            # x、p 均已确定，直接调用同一个 SP；无需重建固定变量的 MP1。
            x, p = self.candidate
            eta, cut = self.sp.solve(x, p)
            row = dict(query=len(self.queries), x=x, p=p, eta=eta, cut=cut)
            self.queries.append(dict(mode="vertex"))
            self.history.append(row)
        else:
            row = self.session.step()
        self.polytopes = update_outer(self.costs, self.designs, self.polytopes,
                                     row, self.queries[row["query"]])
        if row["eta"] is not None and row["cut"] is None:
            add_certificate(self.certified, row)
        if self.phase == "vertex":
            self._next_vertex()
        elif self.session.finished:
            result = self.session.result
            if not self.session.feasible:
                self._next_vertex()
            elif self.phase == "MP2":
                self.session = MP1SP(result["total"], self.network, total=True, log=self.log, sp=self.sp)
                self.phase = "MP1-total"
            else:
                budget = round(result["cost"])-self.quantum
                if budget >= 0:
                    self.session = MP2SP(budget, self.network, log=self.log, sp=self.sp)
                    self.phase = "MP2"
                else:
                    self._next_vertex()
        return row

    @property
    def result(self):
        capacity = next((float(sum(r["p"])) for r in self.history
                         if r["query"] == 0 and r["eta"] is not None and r["cut"] is None), None)
        return dict(status="完整规划—可调度域已认证" if self.finished else "迭代中：蓝色为已认证内域",
                    iterations=len(self.history), cuts=sum(r["cut"] is not None for r in self.history),
                    max_total=capacity, frontier=sorted(self.frontier))

    def run(self, steps=10, pause=1, show=True):
        result = advance(self, steps, pause, (lambda s: print(s.result)) if show else None)
        if show:
            if self.replay is None:
                self.replay = Replay([self])
            self.replay.show(Path("results") / f"{self.network.name}.html")
        return result
