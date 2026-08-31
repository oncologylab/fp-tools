#!/usr/bin/env python3
"""Scan preregistered holdout motifs inside label-free ATAC peak regions."""

from __future__ import annotations

import argparse
from collections import defaultdict
from hashlib import sha256
import gzip
import json
from pathlib import Path
import sys
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))

from fp_tools.utils.fasta import open_fasta  # noqa: E402
from fp_tools.utils.motifs import MotifList  # noqa: E402
from fp_tools.utils.regions import OneRegion, RegionList  # noqa: E402


OUTPUT_COLUMNS = [
    "cell",
    "tf",
    "motif_id",
    "motif_family",
    "TFBS_chr",
    "TFBS_start",
    "TFBS_end",
    "TFBS_strand",
    "motif_score",
    "peak_start",
    "peak_end",
    "chromosome_split",
]


def locked_holdout_tasks(study: dict) -> pd.DataFrame:
    """Return tasks from either the main matrix or a single-task freeze."""

    if "tasks" in study:
        tasks = pd.DataFrame(study["tasks"])
    elif "task" in study and "cell" in study:
        task = dict(study["task"])
        task.setdefault("cell", study["cell"])
        task.setdefault("split", "locked_holdout")
        tasks = pd.DataFrame([task])
    else:
        raise ValueError("study must contain tasks or a cell plus one task")
    required = {"cell", "tf", "motif_id", "motif_family", "split"}
    missing = required.difference(tasks.columns)
    if missing:
        raise ValueError("study tasks lack: " + ", ".join(sorted(missing)))
    return tasks


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def read_and_merge_peaks(path: str | Path, chromosomes: set[str]) -> list[tuple[str, int, int]]:
    path = Path(path)
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with _open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_number} has fewer than three BED columns")
            chromosome = fields[0]
            if chromosome not in chromosomes:
                continue
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number} has non-integer coordinates") from exc
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{line_number} has an invalid interval")
            intervals[chromosome].append((start, end))
    merged: list[tuple[str, int, int]] = []
    for chromosome in sorted(intervals):
        ordered = sorted(intervals[chromosome])
        current_start, current_end = ordered[0]
        for start, end in ordered[1:]:
            if start <= current_end:
                current_end = max(current_end, end)
            else:
                merged.append((chromosome, current_start, current_end))
                current_start, current_end = start, end
        merged.append((chromosome, current_start, current_end))
    if not merged:
        raise ValueError(f"no supported peak intervals were found in {path}")
    return merged


def count_input_peaks(path: str | Path) -> int:
    with _open_text(Path(path)) as handle:
        return sum(
            1
            for line in handle
            if line.strip() and not line.startswith(("#", "track", "browser"))
        )


def chromosome_split(chromosome: str, study: dict) -> str:
    for split, chromosomes in study["chromosome_split"].items():
        if chromosome in chromosomes:
            return split
    return "excluded"


def prepare_motifs(
    motif_path: str | Path,
    motif_ids: Iterable[str],
    background: np.ndarray,
    pvalue: float,
) -> MotifList:
    wanted = set(map(str, motif_ids))
    all_motifs = MotifList().from_file(str(motif_path))
    motifs = MotifList([motif for motif in all_motifs if motif.id in wanted])
    observed = {motif.id for motif in motifs}
    missing = sorted(wanted.difference(observed))
    if missing:
        raise ValueError("motif database lacks: " + ", ".join(missing))
    for motif in motifs:
        motif.bg = np.asarray(background, dtype=float)
        motif.set_prefix("id")
        motif.get_pssm()
        motif.get_threshold(pvalue)
    motifs.set_background()
    motifs.setup_moods_scanner()
    return motifs


def fetch_peak_sequences(
    genome: str | Path,
    peaks: Sequence[tuple[str, int, int]],
) -> tuple[list[tuple[str, int, int, str]], float]:
    records: list[tuple[str, int, int, str]] = []
    gc_bases = 0
    total_bases = 0
    with open_fasta(genome) as fasta:
        bounds = dict(zip(fasta.references, fasta.lengths))
        for chromosome, start, end in peaks:
            if chromosome not in bounds or end > int(bounds[chromosome]):
                raise ValueError(f"peak is outside FASTA bounds: {chromosome}:{start}-{end}")
            sequence = fasta.fetch(chromosome, start, end).upper()
            records.append((chromosome, start, end, sequence))
            gc_bases += sequence.count("G") + sequence.count("C")
            total_bases += sum(sequence.count(base) for base in "ACGT")
    gc_fraction = gc_bases / total_bases if total_bases else 0.5
    return records, float(gc_fraction)


def scan_records(
    records: Sequence[tuple[str, int, int, str]],
    motifs: MotifList,
    tasks: pd.DataFrame,
    study: dict,
    cell: str,
) -> pd.DataFrame:
    task_by_motif = {
        str(row.motif_id): row for row in tasks.itertuples(index=False)
    }
    rows = []
    for chromosome, peak_start, peak_end, sequence in records:
        region = OneRegion([chromosome, peak_start, peak_end])
        by_motif: dict[str, RegionList] = defaultdict(RegionList)
        for site in motifs.scan_sequence(sequence, region):
            by_motif[str(site.name)].append(site)
        for motif_id, motif_sites in by_motif.items():
            task = task_by_motif[motif_id]
            for site in motif_sites.resolve_overlaps(priority="higher"):
                rows.append(
                    {
                        "cell": cell,
                        "tf": str(task.tf),
                        "motif_id": motif_id,
                        "motif_family": str(task.motif_family),
                        "TFBS_chr": str(site.chrom),
                        "TFBS_start": int(site.start),
                        "TFBS_end": int(site.end),
                        "TFBS_strand": str(site.strand),
                        "motif_score": float(site.score),
                        "peak_start": int(peak_start),
                        "peak_end": int(peak_end),
                        "chromosome_split": chromosome_split(str(site.chrom), study),
                    }
                )
    output = pd.DataFrame(rows, columns=OUTPUT_COLUMNS)
    if output.empty:
        raise ValueError(f"motif scan produced no sites for {cell}")
    return output.sort_values(
        ["tf", "TFBS_chr", "TFBS_start", "TFBS_end", "TFBS_strand"],
        kind="mergesort",
    ).reset_index(drop=True)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cell", required=True)
    parser.add_argument("--peaks", type=Path, required=True)
    parser.add_argument("--genome", type=Path, required=True)
    parser.add_argument(
        "--motifs",
        type=Path,
        default=Path(
            "src/fp_tools/resources/motifs/"
            "JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
        ),
    )
    parser.add_argument(
        "--study",
        type=Path,
        default=Path("benchmarks/manifests/footprint_functional_v1.spec.json"),
    )
    parser.add_argument("--motif-pvalue", type=float, default=1e-4)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    if not 0 < args.motif_pvalue < 1:
        raise SystemExit("--motif-pvalue must be between zero and one")
    study = json.loads(args.study.read_text(encoding="utf-8"))
    tasks = locked_holdout_tasks(study)
    tasks = tasks[
        tasks["split"].eq("locked_holdout") & tasks["cell"].astype(str).eq(args.cell)
    ].copy()
    if tasks.empty:
        raise SystemExit(f"study has no locked holdout tasks for {args.cell}")
    chromosomes = set(
        chromosome
        for values in study["chromosome_split"].values()
        for chromosome in values
    )
    peaks = read_and_merge_peaks(args.peaks, chromosomes)
    records, gc_fraction = fetch_peak_sequences(args.genome, peaks)
    background = np.array(
        [(1.0 - gc_fraction) / 2.0, gc_fraction / 2.0, gc_fraction / 2.0, (1.0 - gc_fraction) / 2.0]
    )
    motifs = prepare_motifs(
        args.motifs,
        tasks["motif_id"].astype(str),
        background,
        args.motif_pvalue,
    )
    sites = scan_records(records, motifs, tasks, study, args.cell)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    sites.to_csv(args.out, sep="\t", index=False)
    counts = (
        sites.groupby(["cell", "tf", "motif_id", "motif_family", "chromosome_split"], sort=True)
        .size()
        .rename("sites")
        .reset_index()
    )
    counts_path = args.out.with_name(args.out.name.replace(".tsv.gz", ".counts.tsv"))
    counts.to_csv(counts_path, sep="\t", index=False)
    manifest = {
        "schema": "fp-tools-locked-holdout-motif-scan-v1",
        "locked_holdout_labels_read": False,
        "cell": args.cell,
        "study": str(args.study),
        "study_sha256": file_sha256(args.study),
        "peaks": str(args.peaks),
        "peaks_sha256": file_sha256(args.peaks),
        "genome": str(args.genome),
        "genome_sha256": file_sha256(args.genome),
        "motifs": str(args.motifs),
        "motifs_sha256": file_sha256(args.motifs),
        "motif_pvalue": float(args.motif_pvalue),
        "gc_fraction": gc_fraction,
        "background": background.tolist(),
        "input_peaks": count_input_peaks(args.peaks),
        "merged_peaks": int(len(peaks)),
        "motifs_scanned": int(len(motifs)),
        "sites": int(len(sites)),
        "outputs": {
            args.out.name: {"path": str(args.out), "sha256": file_sha256(args.out)},
            counts_path.name: {"path": str(counts_path), "sha256": file_sha256(counts_path)},
        },
    }
    manifest_path = args.out.with_name(args.out.name.replace(".tsv.gz", ".manifest.json"))
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(counts.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
