"""Wait for every frozen initial endpoint, then summarize without refitting."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time


def main(out):
    receipt=out/'logs/initial_panel_queue.attempt2/exit.json'
    while not receipt.exists():time.sleep(20)
    if json.loads(receipt.read_text())['returncode']!=0:
        raise RuntimeError('Initial panel failed; no partial-panel advancement')
    here=Path(__file__).parent
    for stage in ['metrics','calibration','shape']:
        subprocess.run([sys.executable,'-B',str(here/'analyse.py'),'--outdir',str(out),'--stage',stage],check=True)
        print('Summary completed',stage,flush=True)
    subprocess.run([sys.executable,'-B',str(here/'confounds.py'),'--outdir',str(out)],check=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,required=True)
    main(p.parse_args().outdir.resolve())
