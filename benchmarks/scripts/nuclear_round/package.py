"""Explicit permitted payloads only; no raw-data or recursive data discovery."""
import argparse
from datetime import date
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import zipfile
from access import deny

ROOT=Path(__file__).resolve().parents[3]


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()


def build(out,attempt=None):
    out.mkdir(parents=True,exist_ok=False)
    sanitized=[];input_identities={}
    old=ROOT.parent/'fp-tools-calibration/benchmarks/results/fp_tools_calibration_evidence_2026-10-19'
    deny(old,'verified historical derived-data package')
    if old.is_symlink() or any(p.is_symlink() for p in old.rglob('*')):raise ValueError('symlink in historical payload')
    # This previously verified, privacy-filtered review package contains only
    # authorized derived training data. Its original checksums stay intact.
    shutil.copytree(old,out/'historical_calibration')
    reporting=ROOT/'benchmarks/results/calibration_reporting_20261019_v2'
    if reporting.is_symlink() or any(p.is_symlink() for p in reporting.rglob('*')):raise ValueError('symlink in reporting payload')
    shutil.copytree(reporting,out/'corrected_reporting')
    for name in ['figures.py','slides.py','curve_reporting.py','test_curve_reporting.py','protocol.py']:
        destination=out/'reporting_code'/name;destination.parent.mkdir(exist_ok=True)
        shutil.copyfile(ROOT/'benchmarks/scripts/calibration_round'/name,destination)
    shutil.copyfile(ROOT/'benchmarks/CALIBRATION_REPORTING_CORRECTION_20261019.md',out/'REPORTING_CORRECTIONS.md')
    shutil.copyfile(ROOT/'benchmarks/NUCLEAR_BIAS_ROUND_20261019.md',out/'NUCLEAR_EXPERIMENT_PLAN.md')
    plan=(ROOT/'benchmarks/NUCLEAR_BIAS_ROUND_20261019.md').read_text()
    (out/'ORIGINAL_PROMOTION_REQUIREMENTS.md').write_text('# Original promotion requirements\n\n'+
        plan.split('## Original promotion requirements (unchanged)\n\n',1)[1].split('## Resources and delivery',1)[0])
    (out/'BASELINE_STATE.md').write_text('''# Baseline state

Base commit: 690978f, rejected calibration round. Its source and archive remain
unchanged. This branch first corrected reporting, then froze nuclear sequence
training at source commit 28823cc before any new fit or scoring. See
EXPERIMENT_PLAN.freeze.json when included; it records exact source/input hashes.

The prior nuclear LOG21 grid used the older first-10-million preparation without
a blacklist and was not leave-one-training-chromosome-out. It cannot be reused
as held-chromosome evidence for this round. New windows are genome-grid selected,
without labels or count-based selection, with matched hg38 blacklist exclusion.

Existing corrected extraction, original control arrays, the 6,020-site nested
cohort, complete 21-task ledger and historical calibration scores were already
available. No historical scientific stage was rerun to search for a new result.
All 19 old numerical tables remain byte-identical. The independent old verifier
checked 32 fits, 33,200 exact-tail rows, 16,086 calls and 6,300 fold metrics with
zero failures. Corrected plotting summaries match all 80 available task/arm rows.

Protected labels/controls remain closed. The earlier examined test/transfer
results and ChromBPNet comparisons are historical, not new confirmation.
''')
    status='REPORTING_CORRECTION_ONLY; nuclear scientific experiment not included'
    if attempt:
        if not (attempt/'control.decision.json').exists():raise ValueError('control stage incomplete')
        target=out
        for name in ['arrays','models','tables']:
            shutil.copytree(attempt/name,target/name)
        for name in ['EXPERIMENT_PLAN.freeze.json','control.decision.json','extraction.complete.json','fit.complete.json']:
            shutil.copyfile(attempt/name,target/name)
        for name in ['diagnostics.freeze.json','diagnostics.complete.json','IMPLEMENTATION_CORRECTION.freeze.json','DIAGNOSTICS_IMPLEMENTATION_CORRECTION.freeze.json','context_export.json','environment.json','EXECUTIVE_FINDINGS.md','NEXT_DEEP_RESEARCH_PROMPT.md']:
            if (attempt/name).exists():shutil.copyfile(attempt/name,target/name)
        if (attempt/'provenance').exists():shutil.copytree(attempt/'provenance',target/'provenance')
        if (attempt/'figures').exists():shutil.copytree(attempt/'figures',target/'figures')
        for name in ['EXECUTIVE_FINDINGS.md','NEXT_DEEP_RESEARCH_PROMPT.md']:
            if (attempt/name).exists():shutil.copyfile(attempt/name,out/name)
        (target/'logs').mkdir()
        for name in ['fits.jsonl','access.jsonl','extract.console.log','fit.console.log','control.console.log']:
            if (attempt/'logs'/name).exists():
                text=(attempt/'logs'/name).read_text().replace(str(ROOT.parent),'<LOCAL_ROOT>')
                (target/'logs'/name).write_text(text)
                sanitized.append(dict(file='logs/'+name,original_sha256=digest(attempt/'logs'/name),
                    packaged_sha256=digest(target/'logs'/name),transformation='replace local workspace root by <LOCAL_ROOT>'))
        for name in ['extract.launch.json','extract.exit.json','fit.launch.json','fit.exit.json','control.launch.json','control.exit.json']:
            if (attempt/'logs'/name).exists():
                text=(attempt/'logs'/name).read_text().replace(str(ROOT.parent),'<LOCAL_ROOT>')
                (target/'logs'/name).write_text(text)
                sanitized.append(dict(file='logs/'+name,original_sha256=digest(attempt/'logs'/name),
                    packaged_sha256=digest(target/'logs'/name),transformation='replace local workspace root by <LOCAL_ROOT>'))
        # No local input registry, raw collated BAM, or absolute-path launch logs.
        records=json.loads((attempt/'inputs.local.json').read_text())
        input_identities={k:dict(identity=v['identity'],original_sha256=v['sha256']) for k,v in records.items()}
        inputs=out/'inputs';inputs.mkdir()
        for key,name in [('controls','control_counts.full.tsv.gz'),('control_arrays','control_observed.full.npz'),
                         ('blacklist','hg38-blacklist.v2.bed'),('old_scores','control_reference_scores.tsv.gz'),('ledger','task_ledger.tsv')]:
            p=Path(records[key]['path'])
            deny(p,records[key]['identity'])
            if str(p.resolve())!=records[key]['resolved_path']:raise ValueError('input resolved path changed')
            if digest(p)!=records[key]['sha256']:raise ValueError('input changed')
            shutil.copyfile(p,inputs/name)
        for name in ['verify.py','reanalyse.py','replot.py','run.py','access.py','protocol.py','test_access.py','test_run.py','diagnostics.py','finalise.py','package.py','export_contexts.py','chromatin.py','chromatin_metrics.py','verify_chromatin.py','report_chromatin.py','replot_chromatin.py']:
            dest=out/'code'/name;dest.parent.mkdir(exist_ok=True)
            shutil.copyfile(Path(__file__).parent/name,dest)
        versions=json.loads((attempt/'environment.json').read_text())['packages']
        (out/'code/requirements.txt').write_text(''.join(f'{name}=={versions[name]}\n' for name in
            ['numpy','pandas','scipy','scikit-learn','matplotlib']))
        source=out/'source';source.mkdir()
        for name in ['src/fp_tools/tools/parametric_bias.py','benchmarks/scripts/evaluate_parametric_bias.py',
                     'benchmarks/scripts/conditional_protection_score.py','benchmarks/scripts/fixed_shape_normalization.py','src/fp_tools/tools/functional_footprints.py']:
            destination=source/name;destination.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(ROOT/name,destination)
        status=json.loads((attempt/'control.decision.json').read_text())['decision']
        (out/'LIMITATIONS_AND_BLOCKERS.md').write_text('''# Limits and unrun work

This package does not establish a footprint-detection improvement or promotion.
The primary control screen and native-DWM support prerequisite govern the stop
decision; see control.decision.json. No new chromatin discrimination,
canonical difficult-family rescue, depth/seed replication or end-to-end resource
qualification is inferred from a control-only experiment.

Controls are repeatedly examined naked-DNA replicate 1, shared across tasks.
Repeated physical sites and fragment support induce dependence. Wilson bounds
are descriptive; chromosome bootstrap intervals condition on fixed predictions
and do not establish full refit uncertainty. Historical ChIP-negative calls are
not a true false-positive estimate. Protein-free-DNA calibration alone does not
establish chromatin calibration. Context 21 bases is distinct from one-base
prediction resolution. Training source, depth, preparation and window selection
differ from the frozen mitochondrial reference, limiting causal attribution.

Native production DWM coverage/model-fitting scope remains unqualified. The
same-detector DWM propensity arm is a diagnostic control, not direct TOBIAS.
All locked-control, functional uncertainty and production workflow gates remain
unmet. No protected content was opened and no main/release/manuscript changes
were made. Any further experiment requires a new scientific question and freeze.
''')
    (out/'README_FIRST.md').write_text(
        '# Nuclear-bias follow-up and corrected calibration evidence\n\n'
        f'Status: **{status}**. No production change or promotion.\n\n'
        'Read REPORTING_CORRECTIONS.md before historical material. The historical\n'
        'archive is preserved for provenance; use corrected_reporting/figures for\n'
        'the revised plots and editable slides. Historical scientific tables are unchanged.\n\n'
        'Editable slide text specifies Arial. The Linux PDF preview uses compatible\n'
        'font substitution; both pages were rendered and visually checked.\n\n'
        'Historical independent verification (read-only):\n\n'
        '```bash\npython -I -B historical_calibration/code/verify.py historical_calibration\n```\n\n'
        'If EXPERIMENT_PLAN.freeze.json is present, verify all new counts, scores, calls and\n'
        'control gates independently with `python -B code/verify.py .`. This does\n'
        'not train, extract or repair results. Dependencies: NumPy, pandas, SciPy;\n'
        'write independently recomputed tables elsewhere with\n'
        '`python -B code/reanalyse.py . --outdir ../independent_reanalysis`.\n'
        'Recreate a representative figure with `python -B code/replot.py .\n'
        '--output ../representative_control.png` (also needs Matplotlib).\n'
        'historical verification additionally uses scikit-learn. Raw-data training\n'
        'requires the original allowlisted environment; source/ preserves its code.\n\n'
        'Both rounds use already examined development data. Protected labels and\n'
        'final naked-DNA controls remain closed. DWM means dinucleotide weight matrix;\n'
        'its propensity factorial arm is not a native TOBIAS benchmark.\n')
    (out/'ORIGIN_INDEX.json').write_text(json.dumps(dict(
        historical_commit='690978f',historical_archive_sha256='ef1e2367c37f610441cf8e719f311538b9b30c954bd3ec079f66b1fd250b07c9',
        code_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        scope='authorized derived training data only; no raw BAM/FASTA or protected content',
        inputs=input_identities,sanitized_payloads=sanitized,
        source_mapping={'benchmarks/scripts/nuclear_round/':'code/',
            'original_runner':'provenance/run.pre_control_correction.py',
            'original_diagnostics':'provenance/diagnostics.pre_missingness_correction.py',
            'original_parametric_and_score_modules':'source/ preserves their repository-relative paths',
            'raw_inputs':'intentionally omitted; original hashes retained; derived arrays and short contexts enable reanalysis'},
        original_vs_packaged='historical_calibration retains original packaged hashes; nuclear freeze records original scientific inputs; root SHA256SUMS records this package'),indent=2)+'\n')
    (out/'DATA_DICTIONARY.md').write_text('''# Data dictionary

Historical fields and denominators retain their archived definitions. All nuclear
tables are tab-separated; compressed tables use gzip. Missing values are not zero.

- `site_id`: original stable physical/task-site identifier; `physical_id` and
  `deduplicated_fit_representative` distinguish repeated task windows from unique
  control observations. Repeated task rows are not independent controls.
- `TFBS_chr`: held training chromosome (chr1–15); `cell` denotes the motif-task
  panel, not a separate naked-DNA biological library. All controls use replicate 1.
- `plus/minus_center`, `left`, `right`, `total`: motif-oriented integer cuts;
  center is −7..7, shoulders −40..−12 and 12..40, total their sum per strand.
- `q_plus/q_minus`: sequence-expected conditional central share; not binding
  probabilities. Expected central counts are strand total multiplied by q.
- `score`: capped −log10 conditional two-strand count-tail probability, rounded
  to ten decimal places. Strict call threshold is >1.3010299957. `ceiling` is the
  same score with zero central counts; it is not a biological power estimate.
- `status`: CALL, NO_CALL, ABSTAIN_UNREACHABLE, NO_INFORMATION, UNAVAILABLE or EXCLUDED.
  Abstention is not a negative occupancy label. Eligibility/support counts and
  informative denominators are retained separately in control calibration tables.
- `control.l2_*.npz`: ordered site IDs; plus/minus observed arrays and log bias
  arrays, each 4200 × 201, offsets −100..100. No pickle or raw sequencing payload.
- `training_windows.npz`: chromosome, start, sequence (200+82 bases including
  context margins), and integer counts with axes window × strand × 200 bases.
  `models/folds.json` lists selected row indices; every held chromosome is absent.
- Model NPZ: main base coefficients and pair coefficients at distances 1/5/9.
  JSON sidecars retain hashes, training histories, seeds and fitting provenance.
- `power.tsv`: whole-fragment resampling; r is an endpoint-count weighting
  parameter. `realised_depletion` is the corresponding marginal central-cut
  depletion, undefined for a zero-source-center window. Draws from one window
  are not independent biological libraries.
- Curve plot sources: thresholds are summarized within held chromosomes and
  averaged equally over chromosomes, libraries and available tasks. Mean-curve
  area is not average precision; scalar AP is calculated from original scores.

`FAIL`, `INCOMPLETE`, `NOT_TESTED_LOCKED` and `DEVELOPMENT_CRITERION_MET` are scoped
statuses. No status in this package establishes complete promotion.
''')
    if attempt:
        shutil.copyfile(ROOT/'benchmarks/manifests/end_to_end_resource_protocol_v1.json',out/'end_to_end_resource_protocol_v1.json')
        if (attempt/'chromatin/complete.json').exists():
            sub=out/'chromatin';sub.mkdir()
            for name in ['arrays','tables','figures','provenance']:
                shutil.copytree(attempt/'chromatin'/name,sub/name)
            for name in ['EXPERIMENT_PLAN.freeze.json','complete.json']:
                shutil.copyfile(attempt/'chromatin'/name,sub/name)
            (sub/'logs').mkdir()
            for p in sorted((attempt/'chromatin/logs').iterdir()):
                if not p.is_file() or p.is_symlink():raise ValueError('unexpected log')
                text=p.read_text().replace(str(ROOT.parent),'<LOCAL_ROOT>')
                dest=sub/'logs'/p.name;dest.write_text(text)
                sanitized.append(dict(file=str(dest.relative_to(out)),original_sha256=digest(p),
                    packaged_sha256=digest(dest),transformation='replace local workspace root'))
            shutil.copyfile(attempt/'final.decision.json',out/'final.decision.json')
            (out/'LIMITATIONS_AND_BLOCKERS.md').write_text('''# Limits and blockers

The control component passed; the conditional fixed-cohort chromatin experiment
ran and failed the accuracy screen. See EXECUTIVE_FINDINGS.md and final.decision.json.
All 6 libraries, full/10M seed2026 were evaluated; further depths/seeds and full
resource timing were intentionally stopped. Native production DWM strand and
complete-fragment support remains incomplete; the comparator is DWM propensity,
not a direct TOBIAS run. Locked safety and functional uncertainty remain unqualified.
Previously constructed nested cohorts are reused development evidence, not fresh
confirmation. Bootstrap conditions on fitted predictions. Naked-DNA task windows
share one biological source and dependent fragments. ChIP-negative calls are not
known false positives. Source, depth, preparation and window selection differ
from mitochondrial training. No main, release or manuscript changes were made.

All-task unstratified shape sources and arrays are complete. Additional
count/expected-count-stratified shape panels and full functional-band calibration
were not completed; no canonical difficult-family rescue is claimed. The
chromatin independent verifier checks every new count-tail/call and rank/AP row;
its block aggregation replay reuses the frozen contract. This is explicitly not
an independent reimplementation of cohort matching or raw-data extraction.
''')
            with (out/'README_FIRST.md').open('a') as f:
                f.write('\n## Final conditional-stage decision\n\nThe control receipt above predates the conditional chromatin diagnostic. The final decision is **REJECTED_FOR_FOOTPRINT_IMPROVEMENT**. Read final.decision.json and EXECUTIVE_FINDINGS.md. Verify new chromatin counts, calls and rankings using `python -B code/verify_chromatin.py .`; independently replay headline aggregation with `python -B code/verify_chromatin.py . --output ../chromatin_replay`. The complete derived arrays and all-task plot sources are in chromatin/. Historical slides remain explicitly historical.\n')
            with (out/'BASELINE_STATE.md').open('a') as f:
                f.write('\nThe successful control component triggered the already authorized fixed-cohort chromatin diagnostic, frozen at ff23cdf before scoring. It completed 216,720 rows, including 72,240 new model rows, and failed the accuracy screen. No model selection or new matching followed these results.\n')
            origin=json.loads((out/'ORIGIN_INDEX.json').read_text());origin['sanitized_payloads']=sanitized
            origin['chromatin_input_identities']=json.loads((sub/'EXPERIMENT_PLAN.freeze.json').read_text())['inputs']
            (out/'ORIGIN_INDEX.json').write_text(json.dumps(origin,indent=2)+'\n')
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
