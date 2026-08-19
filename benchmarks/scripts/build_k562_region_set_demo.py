#!/usr/bin/env python3
"""Build the compact K562 CTCF-bound versus matched-control region demo."""

from __future__ import annotations

import gzip
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pyBigWig


ROOT = Path(__file__).resolve().parents[2]
PROJECT = ROOT / "data/public/processed/encode_atac_8cell_20260813"
ATAC_PEAKS = PROJECT / "input/encode_conservative_idr_peaks/ENCFF695IGF.bed.gz"
CHIP_PEAKS = ROOT / "data/public/raw/encode/K562.CTCF.TF_ChIP-seq.ENCFF362OPG.bed.gz"
FOOTPRINTS = [
    PROJECT / "samples/K562_rep1/footprints/K562_rep1_footprints.bw",
    PROJECT / "samples/K562_rep2/footprints/K562_rep2_footprints.bw",
]
CORRECTED = [
    PROJECT / "samples/K562_rep1/normalize/K562_rep1_corrected_q95_scaled.bw",
    PROJECT / "samples/K562_rep2/normalize/K562_rep2_corrected_q95_scaled.bw",
]
GENOME = ROOT / "data/public/raw/genome/hg38.fa"
WORK = ROOT / "benchmarks/results/k562_ctcf_region_set_demo"
DOC_REPORT = ROOT / "docs/demos/reports/region_set_K562_CTCF.html"
DOC_DATA = ROOT / "docs/demos/data/region_set_K562_CTCF_results.tsv"
DOC_QC = ROOT / "docs/demos/data/region_set_K562_CTCF_matching_qc.tsv"


def read_bed(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    rows = []
    with opener(path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith(('#', 'track', 'browser')):
                continue
            fields = line.rstrip("\n").split("\t")
            rows.append((fields[0], int(fields[1]), int(fields[2])))
    return rows


def overlaps_by_chrom(query, reference):
    refs = {}
    for chrom, start, end in reference:
        refs.setdefault(chrom, []).append((start, end))
    for chrom in refs:
        refs[chrom].sort()
    result = []
    cursors = {chrom: 0 for chrom in refs}
    for chrom, start, end in sorted(query):
        intervals = refs.get(chrom, [])
        cursor = cursors.get(chrom, 0)
        while cursor < len(intervals) and intervals[cursor][1] <= start:
            cursor += 1
        cursors[chrom] = cursor
        result.append(bool(cursor < len(intervals) and intervals[cursor][0] < end))
    return result


def remove_overlapping_candidates(frame):
    keep = []
    for _chrom, subset in frame.groupby("chrom", sort=False):
        previous_end = -1
        for index, row in subset.sort_values(["start", "end"]).iterrows():
            if int(row.start) >= previous_end:
                keep.append(index)
                previous_end = int(row.end)
    return frame.loc[keep].copy()


def build_matched_beds(seed=1, strata=50, per_stratum=20):
    WORK.mkdir(parents=True, exist_ok=True)
    peaks = sorted(read_bed(ATAC_PEAKS))
    chip = read_bed(CHIP_PEAKS)
    frame = pd.DataFrame(peaks, columns=["chrom", "start", "end"])
    frame["ctcf_bound"] = overlaps_by_chrom(peaks, chip)
    handles = [pyBigWig.open(str(path)) for path in CORRECTED]
    try:
        baseline = []
        for row in frame.itertuples(index=False):
            values = [handle.stats(row.chrom, row.start, row.end, type="mean")[0] for handle in handles]
            finite = [float(value) for value in values if value is not None and np.isfinite(value)]
            baseline.append(float(np.mean(finite)) if finite else np.nan)
    finally:
        for handle in handles:
            handle.close()
    frame["baseline_signal"] = baseline
    frame = frame[np.isfinite(frame["baseline_signal"])].copy()
    bound = frame[frame["ctcf_bound"]]
    control = frame[~frame["ctcf_bound"]]
    lower = max(bound["baseline_signal"].min(), control["baseline_signal"].min())
    upper = min(bound["baseline_signal"].max(), control["baseline_signal"].max())
    frame = frame[frame["baseline_signal"].between(lower, upper, inclusive="both")].copy()
    frame = remove_overlapping_candidates(frame)
    frame["stratum"] = pd.qcut(
        frame["baseline_signal"], q=strata, labels=False, duplicates="drop"
    ).astype(int)
    rng = np.random.default_rng(seed)
    selected = []
    for stratum, subset in frame.groupby("stratum"):
        group_1 = subset[subset["ctcf_bound"]]
        group_2 = subset[~subset["ctcf_bound"]]
        n = min(len(group_1), len(group_2), per_stratum)
        if n < 2:
            continue
        selected.append(group_1.loc[rng.choice(group_1.index, size=n, replace=False)])
        selected.append(group_2.loc[rng.choice(group_2.index, size=n, replace=False)])
    matched = pd.concat(selected).sort_values(["chrom", "start", "end"]).copy()
    paths = {}
    for is_bound, label in ((True, "CTCF_bound"), (False, "matched_control")):
        output = WORK / f"{label}.bed"
        subset = matched[matched["ctcf_bound"] == is_bound]
        subset[["chrom", "start", "end", "stratum", "baseline_signal"]].to_csv(
            output, sep="\t", header=False, index=False, float_format="%.8g"
        )
        paths[label] = output
    qc = matched.groupby(["stratum", "ctcf_bound"])["baseline_signal"].agg(["count", "mean", "median"])
    qc.to_csv(WORK / "matching_qc.tsv", sep="\t")
    return paths


def main():
    inputs = [ATAC_PEAKS, CHIP_PEAKS, GENOME, *FOOTPRINTS, *CORRECTED]
    missing = [str(path) for path in inputs if not path.exists()]
    if missing:
        raise SystemExit("Missing source file(s): " + ", ".join(missing))
    beds = build_matched_beds()
    output = WORK / "output"
    command = [
        str(ROOT / ".venv/bin/diff-footprints"),
        "--comparison-axis", "regions",
        "--signals", *map(str, FOOTPRINTS),
        "--sample-names", "K562_rep1", "K562_rep2",
        "--cond-names", "K562", "K562",
        "--regions", str(beds["CTCF_bound"]), str(beds["matched_control"]),
        "--region-labels", "CTCF_bound", "matched_control",
        "--region-strata-column", "4",
        "--genome", str(GENOME),
        "--motif-db", "jaspar2026_vertebrates",
        "--aggregate-signals", *map(str, CORRECTED),
        "--plot-aggregate", "top",
        "--plot-aggregate-top-n", "12",
        "--aggregate-flank", "100",
        "--min-regions-per-set", "10",
        "--outdir", str(output),
        "--prefix", "K562_CTCF_regions",
        "--cores", "4",
    ]
    subprocess.run(command, cwd=ROOT, check=True)
    report = output / "K562_CTCF_regions_CTCF_bound_vs_matched_control.html"
    results = output / "K562_CTCF_regions_results.txt"
    DOC_REPORT.parent.mkdir(parents=True, exist_ok=True)
    DOC_DATA.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(report, DOC_REPORT)
    shutil.copyfile(results, DOC_DATA)
    shutil.copyfile(WORK / "matching_qc.tsv", DOC_QC)
    for label, bed in beds.items():
        shutil.copyfile(bed, ROOT / f"docs/demos/data/region_set_K562_{label}.bed")
    print(f"Wrote {DOC_REPORT}")
    print(f"Wrote {DOC_DATA}")


if __name__ == "__main__":
    main()
