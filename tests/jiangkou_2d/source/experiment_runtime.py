"""Experiment-only optimizer heartbeat and stalled physical-search fallback."""
from pathlib import Path
from time import perf_counter, sleep
import json
import numpy as np
from gurobipy import GRB

OUTPUT = Path(__file__).resolve().parents[1] / 'results'
PHYSICAL_STALL_SECONDS = 180.


def serial(value):
    if isinstance(value, dict):
        return {str(k): serial(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [serial(v) for v in value]
    if isinstance(value, np.ndarray):
        return serial(value.tolist())
    if isinstance(value, np.generic):
        return serial(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def save(name, value):
    path = OUTPUT / name
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(serial(value), ensure_ascii=False, indent=2), encoding='utf-8')
    for attempt in range(20):
        try:
            temp.replace(path)
            break
        except PermissionError:
            if attempt==19:
                raise
            sleep(.025)


def optimize(model, mode):
    started, last = perf_counter(), [0.]
    progress = dict(mode=mode, physical_stall_seconds=PHYSICAL_STALL_SECONDS,
                    time_limit=model.Params.TimeLimit)

    def callback(m, where):
        now = perf_counter()
        if where == GRB.Callback.MIP:
            progress.update(bound=m.cbGet(GRB.Callback.MIP_OBJBND),
                            incumbent=m.cbGet(GRB.Callback.MIP_OBJBST),
                            nodes=m.cbGet(GRB.Callback.MIP_NODCNT),
                            solutions=m.cbGet(GRB.Callback.MIP_SOLCNT))
        if now-last[0] >= 15.:
            progress['seconds'] = now-started
            save('solver_progress.json', progress)
            last[0] = now
        # A global bound remains valid after interruption; no feasibility inference.
        if mode == 'physical' and now-started >= PHYSICAL_STALL_SECONDS:
            progress['stop_reason'] = 'physical_search_stalled_switch_to_light'
            m.terminate()

    model.optimize(callback)
    progress.update(seconds=perf_counter()-started, status=int(model.Status),
                    solutions=model.SolCount,
                    bound=float(model.ObjBound),
                    incumbent=float(model.ObjVal) if model.SolCount else None)
    save('solver_progress.json', progress)
    with (OUTPUT/'solver_calls.jsonl').open('a', encoding='utf-8') as stream:
        stream.write(json.dumps(serial(progress), ensure_ascii=False)+'\n')
