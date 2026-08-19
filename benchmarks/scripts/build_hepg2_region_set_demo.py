#!/usr/bin/env python3
"""Build the compact HepG2 HNF4A/FOXA2 region-set demonstration."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import itertools
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pyBigWig


ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data/public/raw/encode"
ALL_SITE = ROOT / "data/public/processed/encode_k562_hepg2_all_site_aggregate_20260817/work"
ATAC_PEAKS = RAW / "HepG2.ATAC.ENCFF609BSU.bed.gz"
HNF4A_PEAKS = RAW / "HepG2.HNF4A.ENCFF704BPD.bed.gz"
FOXA2_PEAKS = RAW / "HepG2.FOXA2.ENCFF656PGC.bed.gz"
FOOTPRINTS = [ALL_SITE / f"footprints/HepG2_rep{rep}_footprints.bw" for rep in (1, 2, 3)]
CORRECTED = [
    ALL_SITE / f"normalized_corrected_q95/HepG2_rep{rep}_corrected.background_scale_q95.bw"
    for rep in (1, 2, 3)
]
GENOME = ROOT / "data/public/raw/genome/hg38.fa"
WORK = ROOT / "benchmarks/results/hepg2_hnf4a_foxa2_region_demo"
DOC_BROWSER = ROOT / "docs/demos/reports/region_set_HepG2_HNF4A_FOXA2"
DOC_DATA = ROOT / "docs/demos/data"

GROUPS = (
    ("HNF4A + FOXA2", "HNF4A_FOXA2"),
    ("HNF4A only", "HNF4A_only"),
    ("FOXA2 only", "FOXA2_only"),
    ("No HNF4A/FOXA2", "No_HNF4A_FOXA2"),
)
DEFAULT_MOTIFS = (
    "MA1494.2",  # HNF4A
    "MA0484.3",  # HNF4G
    "MA0047.4",  # FOXA2
    "MA0148.5",  # FOXA1
    "MA0046.3",  # HNF1A
    "MA0153.2",  # HNF1B
    "MA0102.5",  # CEBPA
    "MA0466.4",  # CEBPB
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--strata", type=int, default=50)
    parser.add_argument("--per-stratum", type=int, default=100)
    parser.add_argument("--cores", type=int, default=4)
    return parser.parse_args()


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_bed(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    rows = []
    with opener(path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            if chrom.startswith("chr") and "_" not in chrom and chrom not in {"chrM", "chrMT"}:
                rows.append((chrom, start, end))
    return sorted(set(rows))


def overlap_flags(query, reference):
    refs = {}
    for chrom, start, end in reference:
        refs.setdefault(chrom, []).append((start, end))
    for intervals in refs.values():
        intervals.sort()
    flags = []
    cursors = {chrom: 0 for chrom in refs}
    for chrom, start, end in query:
        intervals = refs.get(chrom, [])
        cursor = cursors.get(chrom, 0)
        while cursor < len(intervals) and intervals[cursor][1] <= start:
            cursor += 1
        cursors[chrom] = cursor
        flags.append(cursor < len(intervals) and intervals[cursor][0] < end)
    return flags


def remove_overlapping_regions(frame):
    keep = []
    for _chrom, subset in frame.groupby("chrom", sort=False):
        previous_end = -1
        for index, row in subset.sort_values(["start", "end"]).iterrows():
            if int(row.start) >= previous_end:
                keep.append(index)
                previous_end = int(row.end)
    return frame.loc[keep].copy()


def baseline_signal(frame):
    handles = [pyBigWig.open(str(path)) for path in CORRECTED]
    try:
        output = []
        for row in frame.itertuples(index=False):
            values = [handle.stats(row.chrom, row.start, row.end, type="mean")[0] for handle in handles]
            finite = [float(value) for value in values if value is not None and np.isfinite(value)]
            output.append(float(np.mean(finite)) if finite else np.nan)
        return output
    finally:
        for handle in handles:
            handle.close()


def build_matched_regions(seed, strata, per_stratum):
    WORK.mkdir(parents=True, exist_ok=True)
    peaks = read_bed(ATAC_PEAKS)
    frame = pd.DataFrame(peaks, columns=["chrom", "start", "end"])
    frame["hnf4a"] = overlap_flags(peaks, read_bed(HNF4A_PEAKS))
    frame["foxa2"] = overlap_flags(peaks, read_bed(FOXA2_PEAKS))
    frame["group"] = np.select(
        [
            frame["hnf4a"] & frame["foxa2"],
            frame["hnf4a"] & ~frame["foxa2"],
            ~frame["hnf4a"] & frame["foxa2"],
        ],
        [GROUPS[0][0], GROUPS[1][0], GROUPS[2][0]],
        default=GROUPS[3][0],
    )
    frame["baseline_signal"] = baseline_signal(frame)
    frame = frame[np.isfinite(frame["baseline_signal"])].copy()
    frame = remove_overlapping_regions(frame)

    ranges = frame.groupby("group")["baseline_signal"].agg(["min", "max"])
    lower, upper = float(ranges["min"].max()), float(ranges["max"].min())
    frame = frame[frame["baseline_signal"].between(lower, upper, inclusive="both")].copy()
    frame["stratum"] = pd.qcut(
        frame["baseline_signal"], q=strata, labels=False, duplicates="drop"
    ).astype(int)

    rng = np.random.default_rng(seed)
    selected = []
    group_names = [label for label, _slug in GROUPS]
    for stratum, subset in frame.groupby("stratum", sort=True):
        candidates = {label: subset[subset["group"] == label] for label in group_names}
        n = min(per_stratum, *(len(candidates[label]) for label in group_names))
        if n < 2:
            continue
        for label in group_names:
            chosen = rng.choice(candidates[label].index.to_numpy(), size=n, replace=False)
            selected.append(candidates[label].loc[chosen])
    matched = pd.concat(selected, ignore_index=True)

    counts = matched.groupby(["stratum", "group"]).size().unstack(fill_value=0)
    if counts.nunique(axis=1).ne(1).any():
        raise RuntimeError("Region counts are not equal within every retained stratum")
    paths = []
    for label, slug in GROUPS:
        output = WORK / f"{slug}.bed"
        subset = matched[matched["group"] == label].sort_values(["chrom", "start", "end"])
        subset[["chrom", "start", "end", "stratum", "baseline_signal"]].to_csv(
            output, sep="\t", header=False, index=False, float_format="%.8g"
        )
        paths.append(output)

    qc = (
        matched.groupby(["stratum", "group"], sort=True)["baseline_signal"]
        .agg(n_regions="size", mean="mean", median="median", sd="std")
        .reset_index()
    )
    qc.to_csv(WORK / "matching_qc.tsv", sep="\t", index=False, float_format="%.8g")
    summary = matched.groupby("group")["baseline_signal"].agg(n_regions="size", mean="mean", sd="std")
    pooled_sd = float(matched["baseline_signal"].std(ddof=1))
    max_smd = max(
        abs(float(summary.loc[first, "mean"] - summary.loc[second, "mean"])) / pooled_sd
        for first, second in itertools.combinations(summary.index, 2)
    )
    summary = summary.reset_index()
    summary["maximum_pairwise_standardized_mean_difference"] = max_smd
    summary.to_csv(
        WORK / "matching_summary.tsv", sep="\t", index=False, float_format="%.8g"
    )
    if max_smd >= 0.1:
        raise RuntimeError(f"Baseline matching failed: maximum pairwise SMD = {max_smd:.3f}")
    return paths, len(matched) // len(GROUPS), max_smd


def write_source_manifest():
    rows = [
        ("ATAC-seq", "HepG2", "ENCSR291GJU", "ENCFF624SON", "footprint replicate 1", FOOTPRINTS[0]),
        ("ATAC-seq", "HepG2", "ENCSR291GJU", "ENCFF926KFU", "footprint replicate 2", FOOTPRINTS[1]),
        ("ATAC-seq", "HepG2", "ENCSR291GJU", "ENCFF990VCP", "footprint replicate 3", FOOTPRINTS[2]),
        ("ATAC-seq", "HepG2", "ENCSR291GJU", "ENCFF624SON", "aggregate signal replicate 1", CORRECTED[0]),
        ("ATAC-seq", "HepG2", "ENCSR291GJU", "ENCFF926KFU", "aggregate signal replicate 2", CORRECTED[1]),
        ("ATAC-seq", "HepG2", "ENCSR291GJU", "ENCFF990VCP", "aggregate signal replicate 3", CORRECTED[2]),
        ("ATAC-seq", "HepG2", "ENCSR291GJU", "ENCFF609BSU", "accessible regions", ATAC_PEAKS),
        ("TF ChIP-seq", "HNF4A", "ENCSR469FBY", "ENCFF704BPD", "region definition", HNF4A_PEAKS),
        ("TF ChIP-seq", "FOXA2", "ENCSR490AMH", "ENCFF656PGC", "region definition", FOXA2_PEAKS),
    ]
    table = pd.DataFrame(
        [
            {
                "assay": assay,
                "target": target,
                "experiment_accession": experiment,
                "file_accession": accession,
                "role": role,
                "artifact": Path(path).name,
                "artifact_sha256": sha256(path),
            }
            for assay, target, experiment, accession, role, path in rows
        ]
    )
    table.to_csv(WORK / "source_manifest.tsv", sep="\t", index=False)


def gzip_copy(source, destination):
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(gzip.compress(Path(source).read_bytes(), compresslevel=9, mtime=0))


def publish_compact_files(region_paths, output):
    DOC_DATA.mkdir(parents=True, exist_ok=True)
    for path, (_label, slug) in zip(region_paths, GROUPS):
        gzip_copy(path, DOC_DATA / f"region_set_HepG2_{slug}.bed.gz")
    shutil.copyfile(WORK / "matching_qc.tsv", DOC_DATA / "region_set_HepG2_matching_qc.tsv")
    shutil.copyfile(WORK / "matching_summary.tsv", DOC_DATA / "region_set_HepG2_matching_summary.tsv")
    shutil.copyfile(WORK / "source_manifest.tsv", DOC_DATA / "region_set_HepG2_source_manifest.tsv")
    gzip_copy(
        output / "HepG2_HNF4A_FOXA2_regions_results.txt",
        DOC_DATA / "region_set_HepG2_HNF4A_FOXA2_results.tsv.gz",
    )


def main():
    args = parse_args()
    required = [ATAC_PEAKS, HNF4A_PEAKS, FOXA2_PEAKS, GENOME, *FOOTPRINTS, *CORRECTED]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise SystemExit("Missing source file(s): " + ", ".join(missing))

    region_paths, regions_per_group, max_smd = build_matched_regions(
        args.seed, args.strata, args.per_stratum
    )
    write_source_manifest()
    output = WORK / "output"
    command = [
        str(ROOT / ".venv/bin/diff-footprints"),
        "--comparison-axis", "regions",
        "--signals", *map(str, FOOTPRINTS),
        "--sample-names", "HepG2 rep 1", "HepG2 rep 2", "HepG2 rep 3",
        "--cond-names", "HepG2", "HepG2", "HepG2",
        "--regions", *map(str, region_paths),
        "--region-labels", *(label for label, _slug in GROUPS),
        "--region-strata-column", "4",
        "--genome", str(GENOME),
        "--motif-db", "jaspar2026_vertebrates",
        "--aggregate-signals", *map(str, CORRECTED),
        "--plot-aggregate", "top",
        "--plot-aggregate-motifs", *DEFAULT_MOTIFS,
        "--default-aggregate-plots", "8",
        "--aggregate-site-set", "all",
        "--aggregate-flank", "100",
        "--min-regions-per-set", "10",
        "--skip-excel",
        "--outdir", str(output),
        "--prefix", "HepG2_HNF4A_FOXA2_regions",
        "--cores", str(args.cores),
    ]
    subprocess.run(command, cwd=ROOT, check=True)

    reports = sorted(output.glob("HepG2_HNF4A_FOXA2_regions_*_vs_*.html"))
    if len(reports) != 6:
        raise RuntimeError(f"Expected six pairwise reports, found {len(reports)}")
    browser_command = [
        str(ROOT / ".venv/bin/review-multi-comparisons"),
        "--inputs", *map(str, reports),
        "--output-dir", str(DOC_BROWSER),
        "--title", "HepG2 HNF4A/FOXA2 region footprints",
        "--default-comparison", GROUPS[0][0], GROUPS[3][0],
        "--default-aggregate-motifs", *DEFAULT_MOTIFS,
        "--default-aggregate-plots", "8",
        "--documentation-url", "../../../",
    ]
    subprocess.run(browser_command, cwd=ROOT, check=True)
    publish_compact_files(region_paths, output)
    print(f"Regions per group: {regions_per_group:,}")
    print(f"Maximum pairwise baseline SMD: {max_smd:.4f}")
    print(f"Wrote {DOC_BROWSER / 'index.html'}")


if __name__ == "__main__":
    main()
