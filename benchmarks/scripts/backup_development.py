"""Explicit migration backup; content-blind inventory and protected-path denylist.

Not a scientific data reader. Copy/check are delegated to rclone. Completion
requires an unchanged source inventory and a successful remote checksum check.
Never deletes scientific inputs or remote objects. No credentials are exported.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
import time

PROJECTS=['fp-tools','fp-tools-research','fp-tools-next-round','fp-tools-calibration','fp-tools-nuclear-bias']
SKIP={'.git','.venv','.agents','.claude','.codex','.pytest_cache','.ruff_cache','__pycache__',
      '.ipynb_checkpoints','.mypy_cache','.cache','node_modules','development_backup_20261019',
      'build','site','poetry.toml','.env','.netrc','.npmrc','.pypirc','secrets.toml','rclone.conf'}
DENIED=['encsr587trp','encsr485tlp','sk-n-sh','gm23338','nakeddna_rep2','nakeddna_rep3',
        'naked_rep2','naked_rep3','srr13450941','srr13450942','frozen_parametric_holdout']


def write(path,value):
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value,indent=2)+'\n');temporary.replace(path)


def inventory(base,out):
    out.mkdir(parents=True,exist_ok=False)
    denied=list(DENIED)
    # Authorized existing identity metadata only; do not open protected labels.
    manifest=base/'fp-tools-research/benchmarks/manifests/frozen_parametric_factorization_v1.tsv'
    with manifest.open() as f:
        for row in csv.DictReader(f,delimiter='\t'):
            if row['cell_type'] in ['GM23338','SK-N-SH']:
                denied.extend([row['file_accession'].lower(),row['experiment_accession'].lower()])
    records=[];excluded=[];links=[]
    for project in PROJECTS:
        root=base/project
        if not root.is_dir():raise FileNotFoundError(root)
        for directory,dirs,files in os.walk(root,followlinks=False):
            parent=Path(directory)
            for name in list(dirs):
                p=parent/name;relative=p.relative_to(base).as_posix()
                if p.is_symlink():
                    links.append(dict(path=relative,target=os.readlink(p)));dirs.remove(name)
                elif name in SKIP or any(t in relative.lower() for t in denied):
                    excluded.append(dict(path=relative,reason='protected identity' if any(t in relative.lower() for t in denied) else 'rebuildable environment, credential/config or backup state'))
                    dirs.remove(name)
            for name in files:
                p=parent/name;relative=p.relative_to(base).as_posix()
                if p.is_symlink():links.append(dict(path=relative,target=os.readlink(p)));continue
                if name in SKIP or name.startswith('.env.') or name.endswith(('.pem','.key','.pyc')) or any(t in relative.lower() for t in denied):
                    excluded.append(dict(path=relative,reason='protected identity or credential/cache policy'));continue
                st=p.stat()
                if not p.is_file():excluded.append(dict(path=relative,reason='not regular file'));continue
                if '\n' in relative or '\r' in relative:raise ValueError('unsafe file-list name')
                records.append(dict(path=relative,bytes=st.st_size,mtime_ns=st.st_mtime_ns))
    # Smaller reports, manifests and models are secured before large raw inputs.
    records.sort(key=lambda r:(r['bytes'],r['path']))
    with (out/'inventory.jsonl').open('w') as f:
        for row in records:f.write(json.dumps(row)+'\n')
    (out/'files.txt').write_text(''.join(r['path']+'\n' for r in records))
    write(out/'exclusions.json',excluded);write(out/'symlinks.json',links)
    summary=dict(files=len(records),bytes=sum(r['bytes'] for r in records),
        largest_files=sorted(records,key=lambda r:r['bytes'],reverse=True)[:10],
        excluded_entries=len(excluded),symlinks=len(links),protected_contents_copied=False,
        scope='five named project trees; credentials, environments, caches and protected identities excluded; symlinks inventoried but not followed')
    write(out/'inventory_summary.json',summary)
    print(json.dumps(summary,indent=2))


def run(base,out,remote):
    import fcntl
    lock=(out/'job.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    records=[json.loads(line) for line in (out/'inventory.jsonl').read_text().splitlines()]
    def unchanged():
        changed=[]
        for r in records:
            p=base/r['path']
            if p.is_symlink() or not p.is_file():changed.append(r['path']);continue
            st=p.stat()
            if (st.st_size,st.st_mtime_ns)!=(r['bytes'],r['mtime_ns']):changed.append(r['path'])
        if changed:
            write(out/'changed_sources.json',changed)
            raise RuntimeError(f'{len(changed)} source files changed; completion withheld')
    status=dict(state='RUNNING',started_unix=time.time(),pid=os.getpid(),
        boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),
        process_start_ticks=Path('/proc/self/stat').read_text().split()[21],remote=remote,
        files=len(records),bytes=sum(r['bytes'] for r in records),protected_contents_copied=False)
    write(out/'status.json',status)
    def call(command):
        with (out/'commands.jsonl').open('a') as log:
            log.write(json.dumps(dict(time=time.time(),argv=command))+'\n')
        result=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
        with (out/'console.log').open('a') as log:log.write(result.stdout)
        with (out/'commands.jsonl').open('a') as log:
            log.write(json.dumps(dict(time=time.time(),exit_code=result.returncode))+'\n')
        result.check_returncode()
    try:
        unchanged()
        call(['rclone','copy',str(out),remote+'/recovery','--include','*.json','--include','*.jsonl','--include','files.txt','--checksum'])
        call(['rclone','copy',str(base),remote+'/workspace','--files-from-raw',str(out/'files.txt'),
              '--checksum','--immutable','--transfers','4','--checkers','8','--order-by','size,ascending',
              '--retries','8','--low-level-retries','20','--stats','60s','--log-level','INFO','--log-file',str(out/'transfer.log')])
        unchanged();status['state']='VERIFYING';write(out/'status.json',status)
        call(['rclone','check',str(base),remote+'/workspace','--files-from-raw',str(out/'files.txt'),
              '--one-way','--checkers','8','--combined',str(out/'checksum_check.txt'),
              '--log-file',str(out/'check.log'),'--log-level','INFO'])
        unchanged();status['state']='COMPLETE_VERIFIED';status['finished_unix']=time.time()
    except Exception as error:
        status['state']='FAILED_INCOMPLETE';status['error']=str(error);status['finished_unix']=time.time()
        raise
    finally:
        write(out/'status.json',status)
        subprocess.run(['rclone','copy',str(out),remote+'/recovery','--include','*.json','--include','*.jsonl',
            '--include','*.txt','--include','*.log','--checksum'],check=False)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['inventory','run'])
    p.add_argument('--base',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--remote',default='yilab:Yaoxiang/fp-tools/development_backup_2026-10-19')
    a=p.parse_args()
    if a.action=='inventory':inventory(a.base.resolve(),a.out.resolve())
    else:run(a.base.resolve(),a.out.resolve(),a.remote)
