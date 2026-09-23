"""Package an existing evidence folder; never run or repair science."""
import argparse
from hashlib import sha256
from pathlib import Path,PurePosixPath
import stat
import zipfile


def digest(path):
    h=sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()


def execute(root,destination):
    if root.name!='fp_tools_next_round_evidence':raise ValueError('one required top-level folder')
    if destination.exists():raise ValueError('refuse to overwrite an evidence archive')
    listed=[]
    for line in (root/'SHA256SUMS').read_text().splitlines():
        expected,name=line.split('  ',1);p=PurePosixPath(name)
        if p.is_absolute() or '..' in p.parts or (root/name).is_symlink():raise ValueError('unsafe payload')
        if digest(root/name)!=expected:raise ValueError('checksum failure '+name)
        listed.append(name)
    actual={p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file() and p.name!='SHA256SUMS'}
    if actual!=set(listed) or len(listed)!=len(set(listed)):raise ValueError('unmanifested or duplicate payload')
    with zipfile.ZipFile(destination,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as archive:
        for path in sorted(root.rglob('*')):
            if path.is_symlink():raise ValueError('symlinks are not payloads')
            if path.is_file():archive.write(path,root.name+'/'+path.relative_to(root).as_posix())
    with zipfile.ZipFile(destination) as archive:
        if archive.testzip() is not None:raise AssertionError('ZIP CRC failure')
        for item in archive.infolist():
            p=PurePosixPath(item.filename)
            if p.is_absolute() or '..' in p.parts or p.parts[0]!=root.name or stat.S_ISLNK(item.external_attr>>16):raise AssertionError('unsafe archive entry')
    print('ZIP',destination,'bytes',destination.stat().st_size,'SHA256',digest(destination),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('--destination',type=Path,required=True)
    a=p.parse_args();execute(a.root.resolve(),a.destination.resolve())
