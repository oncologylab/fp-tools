"""Explicit input identities; no recursive scientific discovery."""
import hashlib
import json
from pathlib import Path

DENIED = ('encsr587trp','encsr485tlp','sk-n-sh','gm23338','nakeddna_rep2',
          'nakeddna_rep3','naked_rep2','naked_rep3','srr13450941','srr13450942',
          'nutrient')


def deny(path, identity):
    value = f'{path} {Path(path).resolve()} {identity}'.lower()
    if any(token in value for token in DENIED):
        raise PermissionError('protected identity or resolved path')


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for data in iter(lambda: handle.read(8*1024*1024), b''):
            h.update(data)
    return h.hexdigest()


def dump(path, value):
    with Path(path).open('x') as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write('\n')


class Access:
    def __init__(self, out):
        self.out = Path(out)
        self.inputs = json.loads((self.out/'inputs.local.json').read_text())

    def path(self, key):
        record = self.inputs[key]
        path = Path(record['path'])
        deny(path, record['identity'])
        if str(path.resolve()) != record['resolved_path']:
            raise PermissionError('resolved path changed')
        if digest(path) != record['sha256']:
            raise ValueError('input checksum changed: '+key)
        with (self.out/'logs/access.jsonl').open('a') as log:
            log.write(json.dumps(dict(input_id=key,sha256=record['sha256'],identity=record['identity']))+'\n')
        return path
