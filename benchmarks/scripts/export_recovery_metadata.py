"""Export explicit provenance and small aggregate tables, never raw inputs."""
import hashlib
import json
from pathlib import Path
import shutil

ROOT=Path(__file__).resolve().parents[2]


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()


def export():
    target=ROOT/'benchmarks/recovery/20261019';target.mkdir(parents=True,exist_ok=True)
    attempt=ROOT/'benchmarks/results/nuclear_bias_20261019_v1';records=[]
    metadata={
        'nuclear_input_provenance.json':attempt/'inputs.local.json',
        'chromatin_input_provenance.json':attempt/'chromatin/inputs.local.json',
        'calibration_input_provenance.json':ROOT.parent/'fp-tools-calibration/benchmarks/results/calibration_round_20261019_v1/inputs.json',
        'next_round_experiment.freeze.json':ROOT.parent/'fp-tools-next-round/benchmarks/results/next_round_20260922_v1/EXPERIMENT_PLAN.freeze.json',
        'scientific_environment.json':attempt/'environment.json',
        'control_decision.json':attempt/'control.decision.json',
        'final_decision.json':attempt/'final.decision.json',
    }
    for name,source in metadata.items():
        original=source.read_text();json.loads(original)
        (target/name).write_text(original.replace(str(ROOT.parent),'<PROJECTS>'))
        records.append(dict(file=name,source=str(source.relative_to(ROOT.parent)),
            original_sha256=sha(source),packaged_sha256=sha(target/name),
            transformation='replace local project parent with <PROJECTS>; preserve artifact hashes'))
    tables={
        'chromatin_screen_summary.tsv':attempt/'chromatin/tables/screen_summary.tsv',
        'chromatin_task_metrics.tsv':attempt/'chromatin/tables/task_metrics.tsv',
        'chromatin_shape_metrics.tsv':attempt/'chromatin/tables/shape_metrics.tsv',
        'task_ledger.tsv':attempt/'tables/task_ledger.tsv',
        'promotion_gates.tsv':attempt/'tables/promotion_gates.tsv',
        'control_task_means.tsv':attempt/'tables/control_task_means.tsv',
        'control_calls.tsv':attempt/'tables/control_calls.tsv',
        'control_likelihood.tsv':attempt/'tables/paired_likelihood.tsv',
        'power_curve_summary.tsv':attempt/'tables/power_curve_summary.tsv',
    }
    for name,source in tables.items():
        if source.stat().st_size>500000:raise ValueError('not a compact summary')
        shutil.copyfile(source,target/name)
        records.append(dict(file=name,source=str(source.relative_to(ROOT.parent)),
            original_sha256=sha(source),packaged_sha256=sha(target/name),transformation='unchanged aggregate table'))
    archives=[
        ROOT/'benchmarks/results/fp_tools_next_round_evidence_2026-10-19_1473f2d.zip',
        ROOT.parent/'fp-tools-calibration/benchmarks/results/fp_tools_calibration_evidence_2026-10-19_eaa31b1.zip',
        ROOT.parent/'fp-tools-next-round/benchmarks/results/fp_tools_next_round_evidence_2026-09-23_97ec5a4.zip',
        ROOT.parent/'fp-tools-research/benchmarks/results/deep_research_review_20260922_v1/fp_tools_research_review_2026-09-22.zip']
    (target/'archive_index.json').write_text(json.dumps([dict(name=p.name,relative_original_path=str(p.relative_to(ROOT.parent)),
        bytes=p.stat().st_size,sha256=sha(p),storage='original host only; not uploaded') for p in archives],indent=2)+'\n')
    (target/'ORIGIN_INDEX.json').write_text(json.dumps(records,indent=2)+'\n')
    payload=sorted(p for p in target.iterdir() if p.is_file() and p.name!='SHA256SUMS')
    (target/'SHA256SUMS').write_text(''.join(f'{sha(p)}  {p.name}\n' for p in payload))
    print(json.dumps(dict(files=len(payload)+1,bytes=sum(p.stat().st_size for p in target.iterdir() if p.is_file()),
        raw_inputs_read=False,protected_contents_read=False),indent=2))


if __name__=='__main__':export()
