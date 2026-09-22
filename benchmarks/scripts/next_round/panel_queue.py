"""Run the frozen initial depth/library panel, waiting only for qualified inputs."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time
from launch import launch
from access import digest
from prepare import ROOT,dump


def execute(out,attempt=1,workers=8):
    libraries=[f'{c}_rep{r}' for c in ['HepG2','K562'] for r in [1,2,3]]
    jobs=[(library,depth) for library in libraries for depth in ['10M','full']]
    dump(out/f'provenance/initial_panel_queue.attempt{attempt}.freeze.json',dict(jobs=jobs,
        score_code_sha256=digest(Path(__file__).with_name('run.py')),
        posterior_code_sha256=digest(Path(__file__).with_name('shrinkage.py')),
        maximum_jobs=2,workers_per_job=workers,seed=2026,robustness_matrix=False,
        scheduling='start each fixed job when its own library and depth calibration are ready; all12 jobs remain required',
        screen='full only, all three biological libraries; no available-subset replacement'))
    active=[];completed=[]
    while jobs or active:
        for library,depth in list(jobs):
            if len(active)>=2:break
            required=[out/f'logs/{library}.fragments.attempt1/exit.json',out/f'logs/calibration_{depth}.attempt2/exit.json']
            for receipt in required:
                if receipt.exists() and json.loads(receipt.read_text())['returncode']!=0:
                    raise RuntimeError('prerequisite failed; no automatic scientific repair: '+receipt.parent.name)
            if not all(p.exists() for p in required):continue
            inventory=json.loads((out/f'provenance/{library}.fragment_inventory.json').read_text())
            if inventory['complete_fragments']<10000000:raise RuntimeError('10M unavailable; no silent workload substitution')
            jobs.remove((library,depth));name=f'panel_{library}_{depth}.attempt1'
            command=[sys.executable,'-B',str(Path(__file__).with_name('run.py')),'--outdir',str(out),
                     '--stage','score','--library',library,'--depth',depth,'--calibration-depth',depth,'--workers',str(workers)]
            launch(out,name,command);active.append(name)
        for name in list(active):
            receipt=out/'logs'/name/'exit.json'
            if receipt.exists():
                result=json.loads(receipt.read_text())
                if result['returncode']!=0:raise RuntimeError('panel job failed; preserved '+name)
                completed.append(name);active.remove(name)
        if jobs or active:time.sleep(20)
    dump(out/'provenance/initial_panel_execution.json',dict(completed_jobs=completed,exit_codes_zero=True,
        independently_verified=False,scientific_screen_passed=False))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,required=True);p.add_argument('--attempt',type=int,default=1);p.add_argument('--workers',type=int,default=8)
    a=p.parse_args();execute(a.outdir.resolve(),a.attempt,a.workers)
