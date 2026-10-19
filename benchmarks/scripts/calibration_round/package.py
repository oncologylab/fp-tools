"""Assemble the self-contained calibration-round evidence folder and ZIP (no science)."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
import shutil
import zipfile

from inputs import HERE, ROOT, SOURCE_DIR

PERSONAL = re.compile(r'/home/[A-Za-z0-9_.-]+|' r'/tmp/' r'claude[^\s"\']*', re.I)
CODE = ['protocol.py', 'core.py', 'inputs.py', 'run.py', 'test_core.py', 'figures.py', 'ledgers.py',
        'reanalyse.py', 'verify.py', 'package.py', 'slides.py', 'README.md']
DOCS = ['CALIBRATION_ROUND_BASELINE_20261019.md', 'CALIBRATION_ROUND_PLAN_20261019.md',
        'CALIBRATION_ROUND_20261019.md', 'CALIBRATION_ROUND_MODEL_20261019.md']


def sanitize(text):
    return PERSONAL.sub('<LOCAL>', text)


def assemble(result, staging, extras):
    if staging.exists():
        raise FileExistsError(staging)
    staging.mkdir(parents=True)
    registry = json.loads((result / 'inputs.json').read_text())
    for rel in registry['permitted']:
        target = staging / 'data' / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(SOURCE_DIR / rel, target)
    (staging / 'inputs.json').write_text(sanitize((result / 'inputs.json').read_text()))
    shutil.copy2(result / 'PLAN.freeze.json', staging / 'PLAN.freeze.json')
    for sub in ['tables', 'models', 'figures']:
        shutil.copytree(result / sub, staging / 'results' / sub)
    shutil.copy2(result / 'result.json', staging / 'results/result.json')
    (staging / 'results/logs').mkdir()
    for log in list((result / 'logs').iterdir()) + [result.parent / (result.name + '.attempt1.log')]:
        (staging / 'results/logs' / log.name).write_text(sanitize(log.read_text()))
    (staging / 'code').mkdir()
    for name in CODE:
        shutil.copy2(HERE / name, staging / 'code' / name)
    shutil.copy2(HERE.parent / 'conditional_protection_score.py', staging / 'code/conditional_protection_score.py')
    (staging / 'code/requirements.txt').write_text('numpy==2.4.6\npandas==2.3.3\nscipy==1.17.1\nscikit-learn==1.9.0\n'
                                                   'matplotlib==3.10.9\npytest\n')
    (staging / 'docs').mkdir()
    for name in DOCS:
        shutil.copy2(ROOT / 'benchmarks' / name, staging / 'docs' / name)
    for src, rel in extras:
        target = staging / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, target)
    origin = dict(original_inputs_json_sha256=sha256((result / 'inputs.json').read_bytes()).hexdigest(),
                  packaged_inputs_json_note='local paths replaced by <LOCAL>; checksums of inputs unchanged',
                  source_archive=registry['source_archive']['sha256'])
    (staging / 'ORIGIN_INDEX.json').write_text(json.dumps(origin, indent=2) + '\n')
    lines = []
    for path in sorted(p for p in staging.rglob('*') if p.is_file()):
        lines.append(sha256(path.read_bytes()).hexdigest() + '  ' + path.relative_to(staging).as_posix())
    (staging / 'SHA256SUMS').write_text('\n'.join(lines) + '\n')


def privacy(staging):
    hits = []
    for path in staging.rglob('*'):
        if path.is_file() and path.suffix in ('.md', '.json', '.tsv', '.txt', '.py', '.log', '.jsonl'):
            if PERSONAL.search(path.read_text(errors='ignore')):
                hits.append(path.relative_to(staging).as_posix())
    if hits:
        raise ValueError('personal path patterns remain: ' + ', '.join(hits))


def archive(staging, zip_path):
    if zip_path.exists():
        raise FileExistsError(zip_path)
    root = staging.name
    with zipfile.ZipFile(zip_path, 'x', zipfile.ZIP_DEFLATED, compresslevel=9) as z:
        for path in sorted(p for p in staging.rglob('*') if p.is_file()):
            z.write(path, f'{root}/{path.relative_to(staging).as_posix()}')
    print(zip_path, zip_path.stat().st_size, sha256(zip_path.read_bytes()).hexdigest())


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--result', type=Path, required=True)
    parser.add_argument('--staging', type=Path, required=True)
    parser.add_argument('--zip', type=Path)
    parser.add_argument('--extra', nargs=2, action='append', default=[], metavar=('SOURCE', 'RELATIVE'))
    a = parser.parse_args()
    assemble(a.result.resolve(), a.staging.resolve(), [(Path(s), r) for s, r in a.extra])
    privacy(a.staging.resolve())
    if a.zip:
        archive(a.staging.resolve(), a.zip.resolve())
