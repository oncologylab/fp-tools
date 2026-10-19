"""Permitted/protected input lists, identity checks and the immutable plan freeze."""
import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import platform
import sys

from protocol import DEPTHS, LIBRARIES, POLICY, SOURCE

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
NEXT = ROOT.parent / 'fp-tools-next-round'
SOURCE_DIR = NEXT / 'benchmarks/results' / SOURCE
ARCHIVE = NEXT / 'benchmarks/results/fp_tools_next_round_evidence_2026-09-23_97ec5a4.zip'
ARCHIVE_SHA256 = '236aa5e744611a81d263cd90d52d30ef4ed80ae4d9315c7ef01076dad94cf0b6'
PROTECTED = POLICY['protected']


def digest(path):
    h = sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def deny(path, identity=''):
    text = f'{path} {Path(path).resolve()} {identity}'.lower()
    if any(value.lower() in text for value in PROTECTED):
        raise PermissionError(f'protected dataset identity or path: {path}')


def permitted_relative_paths():
    paths = ['tables/control_counts.full.tsv.gz', 'arrays/NakedDNA_rep1.propensities.npz',
             'tables/nested_sites.tsv.gz', 'tables/task_ledger.tsv', 'tables/per_site_scores.tsv.gz',
             'tables/covariate_baseline_all_libraries.tsv.gz', 'tables/empirical_power.tsv',
             'models/frozen_LOG21.npz']
    paths += [f'tables/counts.{library}.{depth}.tsv.gz' for library in LIBRARIES for depth in DEPTHS]
    paths += [f'arrays/empirical_simulations/block_{block:03d}.npz' for block in range(138)]
    return paths


def build_registry():
    entries = {}
    for relative in permitted_relative_paths():
        path = SOURCE_DIR / relative
        identity = f'next-round derived development data: {relative}'
        deny(path, identity)
        entries[relative] = dict(path=str(path), resolved_path=str(path.resolve()), identity=identity,
                                 bytes=path.stat().st_size, sha256=digest(path))
    return dict(schema='fp-tools-calibration-round-inputs-v1', source=SOURCE,
                source_archive=dict(path=str(ARCHIVE), sha256=ARCHIVE_SHA256),
                permitted=entries,
                protected_identities=PROTECTED,
                protected_description=['SK-N-SH ENCSR587TRP labels', 'GM23338 ENCSR485TLP labels',
                                       'final naked-DNA replicate 3 (SRR13450942)', 'nutrient-stress application data'],
                examined_not_used=['naked-DNA replicate 2', 'chr16-22/X transfer labels', 'raw BAM/FASTA/bigWig'],
                policy='only permitted entries are opened; each is re-hashed before use')


class Access:
    """Re-verify identity and checksum each time a permitted input is opened."""

    def __init__(self, registry_path, log_path):
        self.registry = json.loads(Path(registry_path).read_text())['permitted']
        self.log = Path(log_path)

    def path(self, relative):
        entry = self.registry[relative]
        path = Path(entry['path'])
        deny(path, entry['identity'])
        if str(path.resolve()) != entry['resolved_path']:
            raise PermissionError('resolved input changed: ' + relative)
        if digest(path) != entry['sha256']:
            raise ValueError('input checksum changed: ' + relative)
        with self.log.open('a') as handle:
            handle.write(json.dumps(dict(input=relative, sha256=entry['sha256'],
                                         time_utc=datetime.now(timezone.utc).isoformat())) + '\n')
        return path


def dump(path, value):
    with Path(path).open('x') as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write('\n')


def code_identity():
    # Scientific modules only; reporting, verification and packaging code is
    # written later and records its own identity in its outputs.
    files = [HERE / name for name in ['protocol.py', 'core.py', 'inputs.py', 'run.py', 'test_core.py']]
    files.append(HERE.parent / 'conditional_protection_score.py')
    return {str(f.relative_to(ROOT)): digest(f) for f in files}


def freeze(out):
    out = Path(out)
    if out.exists():
        raise FileExistsError('fresh attempt directory required')
    if digest(ARCHIVE) != ARCHIVE_SHA256:
        raise ValueError('next-round archive identity changed')
    out.mkdir(parents=True)
    registry = build_registry()
    dump(out / 'inputs.json', registry)
    import numpy, pandas, scipy, sklearn, matplotlib
    plan = dict(policy=POLICY, code=code_identity(), inputs_sha256=digest(out / 'inputs.json'),
                created_utc=datetime.now(timezone.utc).isoformat(),
                runtime=dict(python=sys.version, platform=platform.platform(), numpy=numpy.__version__,
                             pandas=pandas.__version__, scipy=scipy.__version__, sklearn=sklearn.__version__,
                             matplotlib=matplotlib.__version__),
                status='frozen before any new scoring')
    dump(out / 'PLAN.freeze.json', plan)
    print('frozen', out, digest(out / 'PLAN.freeze.json'))


def check_frozen_code(out):
    plan = json.loads((Path(out) / 'PLAN.freeze.json').read_text())
    current = code_identity()
    if current != plan['code']:
        changed = sorted(k for k in set(current) | set(plan['code']) if current.get(k) != plan['code'].get(k))
        raise RuntimeError('code changed after freeze: ' + ', '.join(changed))
    if plan['policy'] != json.loads(json.dumps(POLICY)):
        raise RuntimeError('policy changed after freeze')
    return plan


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--outdir', type=Path, required=True)
    freeze(parser.parse_args().outdir)
