"""Explicit resolved-path and identity allowlist, with an append-only receipt."""
from hashlib import sha256
import json
from pathlib import Path

DENY = ['ENCSR587TRP','ENCSR485TLP','SK-N-SH','GM23338','NakedDNA_rep3','naked_rep3',
        'naked_dna_10m_rep3','SRR13450942']


def digest(path):
    h=sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(8*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def deny(path, identity):
    text=str(path)+' '+str(Path(path).resolve())+' '+identity
    if any(value.lower() in text.lower() for value in DENY):
        raise PermissionError('protected dataset identity or path')


class Access:
    def __init__(self, registry, log):
        self.entries=json.loads(Path(registry).read_text())['inputs']
        self.log=Path(log)

    def path(self, key, verify=True):
        entry=self.entries[key]
        path=Path(entry['path'])
        deny(path,entry['identity'])
        if str(path.resolve()) != entry['resolved_path']:
            raise PermissionError('resolved input changed')
        if verify and digest(path) != entry['sha256']:
            raise ValueError('input checksum changed: '+key)
        with self.log.open('a') as handle:
            handle.write(json.dumps({'input_id':key,'identity':entry['identity'],'sha256':entry['sha256'],
                                     'checksum_verified':verify})+'\n')
        return path
