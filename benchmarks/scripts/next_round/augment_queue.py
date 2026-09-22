"""Complete the additional fragment-dependence ledger for every library."""
import argparse
import json
from pathlib import Path
import sys
import time
from launch import launch
from access import digest
from prepare import dump


def execute(out):
    waiting=[f'{c}_rep{r}' for c in ['HepG2','K562'] for r in [1,2,3]]
    dump(out/'provenance/fragment_link_queue.freeze.json',dict(libraries=waiting,maximum_jobs=2,
        code_sha256=digest(Path(__file__).with_name('fragment_links.py')),
        scientific_scores_or_thresholds_changed=False))
    active=[];completed=[]
    while waiting or active:
        for library in list(waiting):
            ready=out/f'logs/{library}.fragments.attempt1/exit.json'
            if len(active)>=2:break
            if ready.exists():
                if json.loads(ready.read_text())['returncode']!=0:raise RuntimeError('fragment prerequisite failed')
                name=f'links_{library}.attempt1'
                launch(out,name,[sys.executable,'-B',str(Path(__file__).with_name('fragment_links.py')),'--outdir',str(out),'--library',library])
                waiting.remove(library);active.append(name)
        for name in list(active):
            receipt=out/'logs'/name/'exit.json'
            if receipt.exists():
                if json.loads(receipt.read_text())['returncode']!=0:raise RuntimeError('fragment link job failed: '+name)
                completed.append(name);active.remove(name)
        if waiting or active:time.sleep(20)
    dump(out/'provenance/fragment_link_queue.completed.json',dict(completed=completed,scientific_scores_changed=False))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,required=True);a=p.parse_args();execute(a.outdir.resolve())
