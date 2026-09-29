"""所有正负分区的初始化MP均须提供合格原始物理证书。"""
from itertools import product
import unittest

import numpy as np
from threadpoolctl import threadpool_limits

from Network.case33bw import Case33
from experiments.fourbus_signed import PortPhysics
from model import MasterProblem


class Case33SignedSupportTests(unittest.TestCase):
    def test_all_orthant_initial_supports(self):
        with threadpool_limits(limits=1):
            for nodes in ((18, 25), (18, 25, 30)):
                d = len(nodes)
                for sign in product((1, -1), repeat=d):
                    equations = PortPhysics(Case33(load_nodes=nodes), sign)
                    for direction in (*np.eye(d), np.ones(d)):
                        with self.subTest(nodes=nodes, sign=sign, direction=direction):
                            problem = MasterProblem(equations, budget=7, direction=direction, threads=4)
                            with problem.model:
                                answer = problem.solve()
                                self.assertTrue(answer['feasible'])
                                self.assertLessEqual(problem.model.MaxVio, 1e-8)


if __name__ == '__main__':
    unittest.main()
