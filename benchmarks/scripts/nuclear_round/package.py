"""Explicit permitted payloads only; no raw-data or recursive data discovery."""
import argparse
from datetime import date
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import zipfile

ROOT=Path(__file__).resolve().parents[3]


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()


def build(out,attempt=None):
    out.mkdir(parents=True,exist_ok=False)
    old=ROOT.parent/'fp-tools-calibration/benchmarks/results/fp_tools_calibration_evidence_2026-10-19'
    # This previously verified, privacy-filtered review package contains only
    # authorized derived training data. Its original checksums stay intact.
    shutil.copytree(old,out/'historical_calibration')
    reporting=ROOT/'benchmarks/results/calibration_reporting_20261019_v2'
    shutil.copytree(reporting,out/'corrected_reporting')
    for name in ['figures.py','slides.py','curve_reporting.py','test_curve_reporting.py','protocol.py']:
        destination=out/'reporting_code'/name;destination.parent.mkdir(exist_ok=True)
        shutil.copyfile(ROOT/'benchmarks/scripts/calibration_round'/name,destination)
    shutil.copyfile(ROOT/'benchmarks/CALIBRATION_REPORTING_CORRECTION_20261019.md',out/'REPORTING_CORRECTIONS.md')
    shutil.copyfile(ROOT/'benchmarks/NUCLEAR_BIAS_ROUND_20261019.md',out/'NUCLEAR_EXPERIMENT_PLAN.md')
    status='REPORTING_CORRECTION_ONLY; nuclear scientific experiment not included'
    if attempt:
        if not (attempt/'control.decision.json').exists():raise ValueError('control stage incomplete')
        target=out/'nuclear';target.mkdir()
        for name in ['arrays','models','tables']:
            shutil.copytree(attempt/name,target/name)
        for name in ['EXPERIMENT_PLAN.freeze.json','control.decision.json','extraction.complete.json','fit.complete.json']:
            shutil.copyfile(attempt/name,target/name)
        for name in ['diagnostics.freeze.json','diagnostics.complete.json','EXECUTIVE_FINDINGS.md','NEXT_DEEP_RESEARCH_PROMPT.md']:
            if (attempt/name).exists():shutil.copyfile(attempt/name,target/name)
        if (attempt/'figures').exists():shutil.copytree(attempt/'figures',target/'figures')
        for name in ['EXECUTIVE_FINDINGS.md','NEXT_DEEP_RESEARCH_PROMPT.md']:
            if (attempt/name).exists():shutil.copyfile(attempt/name,out/name)
        (target/'logs').mkdir()
        for name in ['fits.jsonl','access.jsonl','extract.console.log','fit.console.log','control.console.log']:
            if (attempt/'logs'/name).exists():
                text=(attempt/'logs'/name).read_text().replace(str(ROOT.parent),'<LOCAL_ROOT>')
                (target/'logs'/name).write_text(text)
        for name in ['extract.launch.json','extract.exit.json','fit.launch.json','fit.exit.json','control.launch.json','control.exit.json']:
            if (attempt/'logs'/name).exists():
                text=(attempt/'logs'/name).read_text().replace(str(ROOT.parent),'<LOCAL_ROOT>')
                (target/'logs'/name).write_text(text)
        # No local input registry, raw collated BAM, or absolute-path launch logs.
        records=json.loads((attempt/'inputs.local.json').read_text())
        inputs=out/'inputs';inputs.mkdir()
        for key,name in [('controls','control_counts.full.tsv.gz'),('old_scores','control_reference_scores.tsv.gz'),('ledger','task_ledger.tsv')]:
            p=Path(records[key]['path'])
            if digest(p)!=records[key]['sha256']:raise ValueError('input changed')
            shutil.copyfile(p,inputs/name)
        for name in ['verify.py','run.py','access.py','protocol.py','test_access.py','test_run.py','diagnostics.py','finalise.py','package.py']:
            dest=out/'code'/name;dest.parent.mkdir(exist_ok=True)
            shutil.copyfile(Path(__file__).parent/name,dest)
        source=out/'source';source.mkdir()
        for name in ['src/fp_tools/tools/parametric_bias.py','benchmarks/scripts/evaluate_parametric_bias.py',
                     'benchmarks/scripts/conditional_protection_score.py']:
            destination=source/name;destination.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(ROOT/name,destination)
        status=json.loads((attempt/'control.decision.json').read_text())['decision']
    (out/'README_FIRST.md').write_text(
        '# Nuclear-bias follow-up and corrected calibration evidence\n\n'
        f'Status: **{status}**. No production change or promotion.\n\n'
        'Read REPORTING_CORRECTIONS.md before historical material. The historical\n'
        'archive is preserved for provenance; use corrected_reporting/figures for\n'
        'the revised plots and editable slides. Historical scientific tables are unchanged.\n\n'
        'Historical independent verification (read-only):\n\n'
        '```bash\npython -I -B historical_calibration/code/verify.py historical_calibration\n```\n\n'
        'If nuclear/ is present, verify all new regional counts, scores, calls and\n'
        'control gates independently with `python -B code/verify.py .`. This does\n'
        'not train, extract or repair results. Dependencies: NumPy, pandas, SciPy;\n'
        'historical verification additionally uses scikit-learn. Raw-data training\n'
        'requires the original allowlisted environment; source/ preserves its code.\n\n'
        'Both rounds use already examined development data. Protected labels and\n'
        'final naked-DNA controls remain closed. DWM means dinucleotide weight matrix;\n'
        'its propensity factorial arm is not a native TOBIAS benchmark.\n')
    (out/'ORIGIN_INDEX.json').write_text(json.dumps(dict(
        historical_commit='690978f',historical_archive_sha256='ef1e2367c37f610441cf8e719f311538b9b30c954bd3ec079f66b1fd250b07c9',
        code_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        scope='authorized derived training data only; no raw BAM/FASTA or protected content',
        original_vs_packaged='historical_calibration retains original packaged hashes; nuclear freeze records original scientific inputs; root SHA256SUMS records this package'),indent=2)+'\n')
    files=sorted(p for p in out.rglob('*') if p.is_file())
    for p in files:
        if p.is_symlink() or p.suffix in {'.bam','.fastq','.fa','.fasta'}:raise ValueError('unsafe payload')
    (out/'SHA256SUMS').write_text(''.join(f'{digest(p)}  {p.relative_to(out)}\n' for p in files))
    archive=out.with_suffix('.zip')
    if archive.exists():raise FileExistsError(archive)
    with zipfile.ZipFile(archive,'x',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(out.rglob('*')):
            if p.is_file():z.write(p,Path(out.name)/p.relative_to(out))
    with zipfile.ZipFile(archive) as z:
        if z.testzip() is not None:raise ValueError('CRC failure')
        assert len(z.namelist())==len(set(z.namelist()))
        assert all(not Path(n).is_absolute() and '..' not in Path(n).parts for n in z.namelist())
    print(json.dumps(dict(path=str(archive),bytes=archive.stat().st_size,sha256=digest(archive)),indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,required=True);p.add_argument('--attempt',type=Path)
    a=p.parse_args();build(a.outdir.resolve(),a.attempt.resolve() if a.attempt else None)
