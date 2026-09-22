#!/usr/bin/env python3
"""Frozen, label-free chr1 extraction audit; does not rescore any candidate."""

import argparse
import ast
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys

import numpy as np
import pandas as pd
import pysam

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_strand_functional_profiles as current

LEGACY_REF = "4df283a"
EXTRACTOR = "benchmarks/scripts/build_strand_functional_profiles.py"


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(4 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def record(path):
    return {"path": str(Path(path).resolve()), "sha256": sha(path)}


def check(item):
    p = Path(item["path"])
    if sha(p) != item["sha256"]:
        raise ValueError(f"changed frozen input: {p}")
    return p


def oracle(sites, bam):
    """Independent per-site counting from sorted cut positions, not window queries."""
    plus = np.zeros((len(sites), 201), dtype=np.uint32)
    minus = np.zeros_like(plus)
    positions = [[], []]
    with pysam.AlignmentFile(str(bam), "rb") as source:
        length = source.get_reference_length("chr1")
        for read in source.fetch("chr1"):
            if read.flag & (4 | 256 | 512 | 1024 | 2048) or read.mapping_quality < 30:
                continue
            qlen = read.query_length or read.infer_query_length() or 0
            cut = (read.reference_end + qlen - (read.query_alignment_end or qlen) - 5
                   if read.is_reverse else read.reference_start - (read.query_alignment_start or 0) + 4)
            positions[int(read.is_reverse)].append(cut)
    positions = [np.sort(np.asarray(p, dtype=np.int64)) for p in positions]
    valid = np.zeros(len(sites), dtype=bool)
    for i, row in enumerate(sites.itertuples(index=False)):
        if row.TFBS_chr != "chr1":
            raise ValueError("oracle may only inspect training chromosome chr1")
        center = (int(row.TFBS_start) + int(row.TFBS_end)) // 2
        start, end = center - 100, center + 101
        if start < 0 or end > length:
            continue
        valid[i] = True
        for cuts, target in zip(positions, (plus, minus)):
            left, right = np.searchsorted(cuts, [start, end], side="left")
            target[i] = np.bincount(cuts[left:right] - start, minlength=201)
    return plus, minus, valid


def legacy_extractor(blob):
    tree = ast.parse(blob)
    names = {"_usable_read", "_merged_ranges", "extract_strand_cut_profiles"}
    nodes = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name in names]
    if len(nodes) != len(names):
        raise ValueError("historical extractor definitions missing")
    module = ast.Module(body=nodes, type_ignores=[])
    namespace = {"np": np, "pd": pd, "Path": Path,
                 "REQUIRED_COLUMNS": current.REQUIRED_COLUMNS,
                 "cut_position_from_alignment": current.cut_position_from_alignment}
    exec(compile(module, "frozen_legacy_extractor", "exec"), namespace)
    return namespace["extract_strand_cut_profiles"]


def runtime():
    return {"python": platform.python_version(), "numpy": np.__version__,
            "pandas": pd.__version__, "pysam": pysam.__version__}


def freeze(root, output):
    output.mkdir(parents=True, exist_ok=False)
    registry = ROOT / "benchmarks/results/direct_bias_stability_96_20260908_v1/configuration.freeze.json"
    if sha(registry) != "9e870f9a8b66bf3e37f7c330e097ccbecef8c83eb36f4b0985b97e9cb7f1d599":
        raise ValueError("unexpected source matrix")
    matrix = json.loads(registry.read_text())["matrix"]
    inputs = {}
    for cell in ("HepG2", "K562"):
        row = next(r for r in matrix if r["cell"] == cell and r["sample"] == f"{cell}_rep1"
                   and r["depth"] == "10m" and r["seed"] == 2026)
        source = root / f"benchmarks/results/frozen_parametric_factorization_v1_run2/development_inputs/mito_log21_l2_1e-4_shift_4_-5/{cell}.sites.tsv.gz"
        cols = ["cell", "tf", "TFBS_chr", "TFBS_start", "TFBS_end", "TFBS_strand", "chromosome_split"]
        sites = pd.read_csv(source, sep="\t", usecols=cols)
        sites = sites.loc[sites.TFBS_chr.eq("chr1")].reset_index(drop=True)
        if sites.empty or not sites.chromosome_split.eq("train").all() or not sites.cell.eq(cell).all():
            raise ValueError("invalid training-only subset")
        subset = output / f"{cell}.chr1.sites.tsv"
        sites.to_csv(subset, sep="\t", index=False)
        entries = {"source_sites": record(source), "sites": record(subset),
                   "bam": {"path": str(root / row["bam"]), "sha256": row["bam_sha256"]},
                   "bai": {"path": str(root / row["bai"]), "sha256": row["bai_sha256"]}}
        for entry in entries.values():
            check(entry)
        inputs[cell] = entries
    blob = subprocess.check_output(["git", "show", f"{LEGACY_REF}:{EXTRACTOR}"], cwd=ROOT)
    document = {"schema": "fp-tools-strand-extraction-audit-v1", "chromosomes": ["chr1"],
                "scope": "training_only_implementation_diagnostic_no_labels_or_scoring",
                "flank": 100, "read_shift": [4, -5], "minimum_mapq": 30,
                "keep_duplicates": False, "inputs": inputs, "runtime": runtime(),
                "legacy_ref": LEGACY_REF, "legacy_sha256": hashlib.sha256(blob).hexdigest(),
                "code": [record(ROOT / p) for p in ("benchmarks/scripts/audit_strand_extraction.py",
                         EXTRACTOR, "src/fp_tools/tools/parametric_bias.py")],
                "source_registry": record(registry), "new_holdout_labels_read": False}
    target = output / "configuration.freeze.json"
    target.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")
    (output / "configuration.sha256").write_text(sha(target) + "\n")
    print("Frozen chr1 coordinates, BAM identities, legacy/current code and runtime before extraction", flush=True)


def run(output):
    target = output / "configuration.freeze.json"
    if sha(target) != (output / "configuration.sha256").read_text().strip():
        raise ValueError("configuration changed")
    config = json.loads(target.read_text())
    if config["runtime"] != runtime():
        raise ValueError("runtime changed")
    for entry in config["code"] + [config["source_registry"]]:
        check(entry)
    blob = subprocess.check_output(["git", "show", f"{config['legacy_ref']}:{EXTRACTOR}"], cwd=ROOT)
    if hashlib.sha256(blob).hexdigest() != config["legacy_sha256"]:
        raise ValueError("legacy code changed")
    if (output / "audit.json").exists():
        raise ValueError("completed audit exists; do not rerun")
    legacy = legacy_extractor(blob)
    rows = []
    for cell, inputs in config["inputs"].items():
        for entry in inputs.values():
            check(entry)
        sites = pd.read_csv(inputs["sites"]["path"], sep="\t")
        if not sites.TFBS_chr.eq("chr1").all() or not sites.chromosome_split.eq("train").all():
            raise ValueError("training-only firewall failed")
        bam = inputs["bam"]["path"]
        kwargs = dict(flank=100, read_shift=(4, -5), minimum_mapq=30, keep_duplicates=False)
        old, new, independent = legacy(sites, bam, **kwargs), current.extract_strand_cut_profiles(sites, bam, **kwargs), oracle(sites, bam)
        if any(not np.array_equal(a, b) for a, b in zip(new, independent)):
            raise ValueError("corrected extractor differs from independent oracle")
        if not np.array_equal(old[2], new[2]):
            raise ValueError("validity support changed")
        delta = np.stack([new[s].astype(np.int64) - old[s] for s in (0, 1)], axis=1)
        for tf, group in sites.groupby("tf", sort=True):
            ix = group.index.to_numpy()
            part = delta[ix]
            rows.append({"cell": cell, "tf": tf, "sites": len(ix),
                         "sites_changed": int(np.any(part != 0, axis=(1, 2)).sum()),
                         "legacy_total_cuts": int(old[0][ix].sum() + old[1][ix].sum()),
                         "corrected_total_cuts": int(new[0][ix].sum() + new[1][ix].sum()),
                         "removed_duplicate_counts": int(-part[part < 0].sum()),
                         "recovered_counts": int(part[part > 0].sum())})
        np.savez_compressed(output / f"{cell}.comparison.npz", site_hash=current.site_hashes(sites),
                            old_plus=old[0], old_minus=old[1], new_plus=new[0], new_minus=new[1], valid=new[2])
        print(f"{cell}: independent oracle matched all {len(sites)} training sites", flush=True)
    table = output / "per_tf_counts.tsv"
    pd.DataFrame(rows).to_csv(table, sep="\t", index=False)
    result = {"configuration_sha256": sha(target), "oracle_parity": True,
              "new_holdout_labels_read": False, "candidate_scoring_performed": False,
              "outputs": [record(table)] + [record(output / f"{c}.comparison.npz") for c in config["inputs"]]}
    (output / "audit.json").write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--freeze-only", action="store_true")
    args = parser.parse_args()
    output = args.outdir.resolve()
    if not output.is_relative_to(ROOT / "benchmarks/results"):
        parser.error("output must stay in the research worktree's benchmark results")
    if args.freeze_only:
        freeze(args.data_root.resolve(), output)
    else:
        run(output)
