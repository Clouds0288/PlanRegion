"""Same-query comparison of Gurobi proof strategies; numerical model unchanged."""
import json
from pathlib import Path
import numpy as np
from threadpoolctl import threadpool_limits
from tests.case33_3d import PlanningOracle3D
from tests.case33_two_stage import save


def probe():
    root = Path('results/case33_3d/20260924')
    result = json.loads((root/'budget_4/result.json').read_text(encoding='utf-8'))
    event = next(e for e in reversed(result['events']) if e['mode'] == 'distance' and e['seconds'] > 40
                 and e['status'] == 2 and e['bound'] > 20.)
    answers = []
    for focus in (0,2,3):
        oracle = PlanningOracle3D(4.,threads=2)
        for cut in result['supports']:
            oracle.add_support(cut['weights'],cut['bound'])
        oracle.certificates = [{**c,**{k:np.asarray(c[k]) for k in ('x','p','state')}} for c in result['certificates']]
        oracle.problem.model.Params.MIPFocus = focus
        answer = oracle.solve(target=event['target'],time_limit=90.,gap_kw=4.,warm_count=8,
                              label=f'focus_{focus}')
        answers.append(dict(focus=focus,**{k:answer[k] for k in ('status','bound','objective','seconds','nodes','feasible')}))
        save(root/'parameter_probe.json',dict(target=event['target'],threads=2,results=answers,
            source_hashes=result['source_hashes'],official_parameter_reference=
            'https://docs.gurobi.com/projects/optimizer/en/current/reference/parameters.html#parameter-MIPFocus'))
        oracle.close()


if __name__ == '__main__':
    with threadpool_limits(limits=1):
        probe()
