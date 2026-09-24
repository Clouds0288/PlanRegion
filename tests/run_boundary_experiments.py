"""Serial fresh-process experiment schedule; stops immediately on a failed audit."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--stage', choices=('screen', 'formal'), required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    output = Path(args.output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    schedule = []
    if args.stage == 'screen':
        schedule = [(variant, 2., 0, 120.) for variant in ('A', 'A-L', 'B', 'C')]
    else:
        for seed in range(3):
            for budget in (0., 2., float('inf')):
                variants = ['A', 'A-L', 'B', 'C'] if budget == 2. else ['A', 'B', 'C']
                variants = variants[seed:]+variants[:seed]
                schedule.extend((variant, budget, seed, 300.) for variant in variants)
    inputs = ['tests/benchmark_boundary_search.py', 'tests/test_boundary_search.py',
              'tests/run_boundary_experiments.py', 'docs/notation.md', 'docs/case33_boundary_test_plan.md',
              'main.py', 'model.py', 'region.py', 'vertify.py', 'Network/__init__.py',
              'Network/case33bw.py', 'Network/data/case33bw.m', 'tests/planning_checks.py',
              'tests/benchmark_math_acceleration.py']
    for name in inputs:
        target = output/'source'/name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT/name, target)
    (output/'schedule.json').write_text(json.dumps([
        dict(variant=v, budget=None if b == float('inf') else b, seed=s, seconds=t)
        for v,b,s,t in schedule], indent=2), encoding='utf-8')
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    for index, (variant, budget, seed, seconds) in enumerate(schedule):
        name = f'{index+1:02d}_{variant}_b{budget:g}_s{seed}'
        directory = output/name
        if directory.exists():
            raise FileExistsError(f'Refusing to overwrite an experiment: {directory}')
        directory.mkdir()
        command = [sys.executable, '-X', 'utf8', '-m', 'tests.benchmark_boundary_search',
                   '--variant', variant, '--budget', str(budget), '--seconds', str(seconds),
                   '--seed', str(seed), '--threads', '20', '--output', str(directory)]
        print(json.dumps(dict(event='start', case=name, utc=datetime.now(timezone.utc).isoformat())), flush=True)
        with (directory/'console.log').open('w', encoding='utf-8') as log:
            finished = subprocess.run(command, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
        if finished.returncode:
            print(json.dumps(dict(event='failed', case=name, returncode=finished.returncode)), flush=True)
            sys.exit(finished.returncode)
        data = json.loads((directory/'result.json').read_text(encoding='utf-8'))
        print(json.dumps(dict(event='done', case=name, status=data['status'], seconds=data['seconds'],
                              counts=data['counts'], grid=data['grid_progress'])), flush=True)


if __name__ == '__main__':
    main()
