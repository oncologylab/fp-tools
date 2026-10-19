"""Bounded-space, checksum-verified tar-part transport for a frozen inventory.

Uses the content-blind inventory made by backup_development.py. Two-GiB parts
avoid per-file Box API overhead and large-file upload limits. Resume regenerates
the deterministic tar stream and verifies existing parts; no source is deleted.
"""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time


def save(path,data):
    temp=path.with_suffix('.tmp');temp.write_text(json.dumps(data,indent=2)+'\n');temp.replace(path)


def run(base,out,remote,part_bytes=2*1024**3):
    lock=(out/'job.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    records=[json.loads(x) for x in (out/'inventory.jsonl').read_text().splitlines()]
    def unchanged():
        bad=[]
        for r in records:
            p=base/r['path']
            if p.is_symlink() or not p.is_file():bad.append(r['path']);continue
            st=p.stat()
            if (st.st_size,st.st_mtime_ns)!=(r['bytes'],r['mtime_ns']):bad.append(r['path'])
        if bad:
            save(out/'changed_sources.json',bad)
            raise RuntimeError(f'{len(bad)} changed sources; snapshot cannot be certified')
    def command(argv):
        with (out/'commands.jsonl').open('a') as f:f.write(json.dumps(dict(time=time.time(),argv=argv))+'\n')
        result=subprocess.run(argv,capture_output=True,text=True)
        with (out/'commands.jsonl').open('a') as f:f.write(json.dumps(dict(time=time.time(),exit_code=result.returncode))+'\n')
        if result.returncode:raise RuntimeError(result.stderr[-3000:])
        return result.stdout
    def sync_metadata():
        command(['rclone','copy',str(out),remote+'/recovery','--include','*.json','--include','*.jsonl',
                 '--include','*.txt','--include','*.md','--include','*.py','--checksum','--transfers','4'])
    status=dict(state='RUNNING',format='ordered GNU tar stream split into parts',started_unix=time.time(),
        pid=os.getpid(),boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),
        process_start_ticks=Path('/proc/self/stat').read_text().split()[21],remote=remote,
        files=len(records),source_bytes=sum(r['bytes'] for r in records),verified_bytes=0,verified_parts=0,
        protected_contents_copied=False)
    parts=[];tar=None
    try:
        unchanged();save(out/'status.json',status);sync_metadata()
        names=out/'files.null';names.write_bytes(b''.join(r['path'].encode()+b'\0' for r in records))
        stderr=(out/'tar.stderr.log').open('a')
        tar=subprocess.Popen(['tar','--format=gnu','--numeric-owner','--owner=0','--group=0',
            '--no-recursion','-C',str(base),'--null','-T',str(names),'-cf','-'],stdout=subprocess.PIPE,stderr=stderr)
        index=0
        while True:
            scratch=out/'part.incomplete'
            if scratch.exists():scratch.rename(out/f'interrupted_part_{time.time_ns()}.bin')
            sha1=hashlib.sha1();sha256=hashlib.sha256();size=0
            with scratch.open('xb') as f:
                while size<part_bytes:
                    data=tar.stdout.read(min(8*1024**2,part_bytes-size))
                    if not data:break
                    f.write(data);sha1.update(data);sha256.update(data);size+=len(data)
            if not size:
                scratch.unlink();break
            name=f'workspace.tar.part-{index:06d}'
            remote_path=remote+'/parts/'+name
            command(['rclone','copyto',str(scratch),remote_path,'--immutable','--checksum',
                '--retries','8','--low-level-retries','20','--stats','60s','--log-file',str(out/'transfer.log')])
            info=json.loads(command(['rclone','lsjson',remote_path,'--stat','--hash']))
            if info['Size']!=size or info['Hashes'].get('sha1')!=sha1.hexdigest():raise RuntimeError('remote part checksum mismatch')
            parts.append(dict(name=name,bytes=size,sha1=sha1.hexdigest(),sha256=sha256.hexdigest()))
            save(out/'parts.json',parts)
            status.update(verified_parts=len(parts),verified_bytes=sum(p['bytes'] for p in parts),updated_unix=time.time())
            save(out/'status.json',status);sync_metadata()
            # Delete only the successfully verified transport spool, never an input.
            scratch.unlink();index+=1
        if tar.wait()!=0:raise RuntimeError('tar reported an incomplete or changing source')
        unchanged();status['state']='VERIFYING';save(out/'status.json',status);sync_metadata()
        listing=json.loads(command(['rclone','lsjson',remote+'/parts','--files-only','--hash']))
        actual={r['Name']:r for r in listing}
        if set(actual)!={p['name'] for p in parts}:raise RuntimeError('remote part-set mismatch')
        for p in parts:
            r=actual[p['name']]
            if r['Size']!=p['bytes'] or r['Hashes'].get('sha1')!=p['sha1']:raise RuntimeError('final part checksum mismatch')
        status['state']='COMPLETE_VERIFIED';status['finished_unix']=time.time()
    except Exception as error:
        status.update(state='FAILED_INCOMPLETE',error=str(error),finished_unix=time.time())
        raise
    finally:
        if tar is not None and tar.poll() is None:tar.terminate();tar.wait()
        save(out/'status.json',status)
        try:sync_metadata()
        except Exception as error:
            (out/'metadata_upload_error.log').write_text(str(error))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--base',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--remote',required=True);p.add_argument('--part-bytes',type=int,default=2*1024**3)
    a=p.parse_args();run(a.base.resolve(),a.out.resolve(),a.remote,a.part_bytes)
