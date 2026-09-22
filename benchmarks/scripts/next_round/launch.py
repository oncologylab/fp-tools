"""Launch an immutable supervised attempt; never infer completion from launch."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
from access import digest
from prepare import ROOT,dump


def launch(out,name,command,trace=False):
    out=Path(out).resolve()
    spec=out/'provenance'/f'{name}.job.json'
    if trace:command=['/usr/bin/strace','-f','-e','trace=openat','-o',str(out/'logs'/f'{name}.openat.local.log')]+command
    value=dict(schema='fp-tools-logged-experiment-job-v1',name=name,cwd=str(ROOT),
        outdir=str(out/'logs'/name),environment={k:'1' for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','PYTHONUNBUFFERED']},
        command=command)
    dump(spec,value)
    with (out/'logs'/f'{name}.supervisor.log').open('xb') as log:
        child=subprocess.Popen([sys.executable,'-B',str(ROOT/'benchmarks/scripts/run_logged_experiment.py'),
                                '--spec',str(spec),'--spec-sha256',digest(spec)],
                                stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True,
                                env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'})
    print('Launched supervisor',child.pid,'for',name,'; outputs not yet verified',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,required=True);p.add_argument('--name',required=True)
    p.add_argument('--trace',action='store_true');p.add_argument('command',nargs=argparse.REMAINDER)
    a=p.parse_args();command=a.command[1:] if a.command[0]=='--' else a.command
    launch(a.outdir,a.name,command,a.trace)
