#!/usr/bin/env python3
"""Apply the unchanged descriptive gate ledger to verified corrected counts."""

import argparse
import csv
import gzip
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_direct_bias_gate_ledger as ledger
from audit_chrombpnet_comparison import check_identity, checked_path, canonical_id
from audit_direct_bias_safety import record, write_immutable

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "benchmarks/results/corrected_count_reassessment_20260922_v1/reassessment_manifest.json"
SOURCE_SHA = "71235831ce3559b5e5c81a5b30e4564879cf7afb11593cc69c221e9eb3d8b979"
VERIFICATION = ROOT / "benchmarks/results/corrected_count_jobs_20260922/verify_attempt1/exit.json"


def build(out):
    manifest = json.loads(checked_path(dict(path=str(SOURCE), sha256=SOURCE_SHA), ROOT).read_text())
    check_identity(manifest, "stage_id")
    if (manifest.get("schema") != "fp-tools-corrected-count-reassessment-v1" or len(manifest.get("stages", [])) != 96
            or any(manifest.get(k) is not False for k in
                   ("candidate_retuned", "new_holdout_labels_read", "final_untouched_control_read", "promotion_eligible"))):
        raise ValueError("incorrect corrected-count evidence scope")
    verified = json.loads(VERIFICATION.read_text())
    check_identity(verified, "exit_id")
    if verified["returncode"] != 0:
        raise ValueError("independent verification did not succeed")
    launch = json.loads(checked_path(verified["launch"], ROOT).read_text())
    if "--verify-only" not in " ".join(launch["command"]) or str(SOURCE.parent) not in " ".join(launch["command"]):
        raise ValueError("verification receipt refers to another command")
    tables = {}
    for name in ("metrics", "shapes", "coverage", "bootstrap"):
        with gzip.open(checked_path(manifest["outputs"][name], ROOT), "rt") as stream:
            tables[name] = list(csv.DictReader(stream, delimiter="\t"))
    tasks, exceptions = ledger.summarize(tables)
    result = dict(schema="fp-tools-corrected-count-gate-ledger-v1", source=record(SOURCE),
        independent_verification=record(VERIFICATION), code=[record(Path(__file__)), record(Path(ledger.__file__))],
        evidence="post_hoc_fixed_candidate_corrected_count_reassessment", task_count=len(tasks),
        tasks=tasks, shape_exceptions=exceptions, promotion_eligible=False,
        scope="Unchanged descriptive ledger on corrected counts. Not the locked mean-gain or untouched promotion test. Safety, functional uncertainty, output compatibility, and end-to-end resources remain separate.")
    result["stage_id"] = canonical_id(result)
    out.mkdir(parents=True, exist_ok=True)
    write_immutable(out / "gate_ledger.json", result)
    print(f"Verified corrected ledger: {len(tasks)} tasks, {len(exceptions)} shape exceptions, no promotion")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    out = args.outdir.resolve()
    if not out.is_relative_to(ROOT / "benchmarks/results"):
        parser.error("use research results")
    build(out)
