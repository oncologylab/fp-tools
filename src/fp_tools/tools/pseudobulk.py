#!/usr/bin/env python
"""Pseudobulk fragment grouping utility."""

from __future__ import annotations

import argparse
import gzip
import multiprocessing
import re
import shlex
from concurrent.futures import ProcessPoolExecutor
from collections import defaultdict
from pathlib import Path

import pandas as pd
import yaml

try:
    import pyBigWig
except ImportError:  # pragma: no cover - dependency is declared, keep import-time failure friendly
    pyBigWig = None

try:
    import pysam
except ImportError:  # pragma: no cover
    pysam = None


def _open_text(path: str | Path):
    path = Path(path)
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def _annotation_separator(path: str | Path) -> str:
    if Path(path).suffix.lower() == ".csv":
        return ","
    with _open_text(path) as handle:
        sample = handle.readline()
    return "," if sample.count(",") > sample.count("\t") else "\t"


def _normalize_barcode(barcode: str, strip_suffix: bool = True) -> str:
    barcode = str(barcode).strip()
    return re.sub(r"-\d+$", "", barcode) if strip_suffix else barcode


def _safe_group_label(label: str) -> str:
    label = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(label).strip())
    return label.strip("_") or "group"


def _parse_chrom_list(values: str | None) -> set[str] | None:
    if not values:
        return None
    return {value.strip() for value in values.split(",") if value.strip()}


def load_genome_sizes(path: str | Path) -> list[tuple[str, int]]:
    rows = []
    with _open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, size, *_ = line.rstrip("\n").split("\t")
            rows.append((chrom, int(size)))
    if not rows:
        raise ValueError(f"No genome sizes found in {path}")
    return rows


def load_annotations(
    path: str | Path,
    barcode_column: str,
    group_by: list[str],
    strip_barcode_suffix: bool = True,
) -> dict[str, str]:
    """Map barcodes to pseudobulk group labels."""

    frame = pd.read_csv(path, sep=_annotation_separator(path))
    missing = [column for column in [barcode_column, *group_by] if column not in frame.columns]
    if missing:
        raise ValueError(f"Missing annotation columns: {missing}")
    mapping = {}
    for _, row in frame.iterrows():
        group = "__".join(_safe_group_label(row[column]) for column in group_by)
        barcode = _normalize_barcode(row[barcode_column], strip_barcode_suffix)
        mapping[barcode] = group
    return mapping


def _compress_no_index(path: Path) -> Path:
    gz_path = path.with_suffix(path.suffix + ".gz")
    with path.open("rb") as source, gzip.open(gz_path, "wb") as target:
        target.writelines(source)
    path.unlink()
    return gz_path


def _compress_and_index(path: Path) -> Path:
    if pysam is None:
        raise RuntimeError("pysam is required for --index-output")
    gz_path = path.with_suffix(path.suffix + ".gz")
    pysam.tabix_compress(str(path), str(gz_path), force=True)
    pysam.tabix_index(str(gz_path), preset="bed", force=True)
    path.unlink()
    return gz_path


def write_cutsite_bigwig(
    fragment_file: str | Path,
    output: str | Path,
    genome_sizes: str | Path,
    cpm: bool = True,
) -> Path:
    """Write a sparse cut-site bigWig from 10x-style fragment intervals."""

    if pyBigWig is None:
        raise RuntimeError("pyBigWig is required for --write-cutsite-bigwigs")
    genome = load_genome_sizes(genome_sizes)
    chrom_sizes = dict(genome)
    counts: dict[str, dict[int, int]] = defaultdict(lambda: defaultdict(int))
    total_cuts = 0

    with _open_text(fragment_file) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            chrom = fields[0]
            if chrom not in chrom_sizes:
                continue
            start = max(0, int(fields[1]))
            end = min(chrom_sizes[chrom], int(fields[2]))
            if end <= start:
                continue
            multiplicity = int(fields[4]) if len(fields) > 4 and fields[4].isdigit() else 1
            counts[chrom][start] += multiplicity
            counts[chrom][end - 1] += multiplicity
            total_cuts += 2 * multiplicity

    scale = (1_000_000.0 / total_cuts) if cpm and total_cuts else 1.0
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    bw = pyBigWig.open(str(output), "w")
    try:
        bw.addHeader(genome)
        for chrom, _size in genome:
            positions = sorted(counts.get(chrom, {}))
            if not positions:
                continue
            bw.addEntries(
                [chrom] * len(positions),
                positions,
                ends=[position + 1 for position in positions],
                values=[counts[chrom][position] * scale for position in positions],
            )
    finally:
        bw.close()
    return output


def _compress_group(args: tuple[str, str, bool]) -> tuple[str, str]:
    group, path_text, index_output = args
    path = Path(path_text)
    compressed = _compress_and_index(path) if index_output else _compress_no_index(path)
    return group, str(compressed)


def _write_group_bigwig(args: tuple[str, str, str, str, bool]) -> tuple[str, str]:
    group, fragment_file, output, genome_sizes, cpm_normalize = args
    write_cutsite_bigwig(fragment_file, output, genome_sizes, cpm=cpm_normalize)
    return group, output


def group_fragments(
    fragments: str | Path,
    annotations: str | Path,
    outdir: str | Path,
    group_by: list[str],
    barcode_column: str = "barcode",
    min_cells: int = 1,
    min_fragments: int = 1,
    strip_barcode_suffix: bool = True,
    include_chroms: set[str] | None = None,
    exclude_chroms: set[str] | None = None,
    compress_output: bool = False,
    index_output: bool = False,
    genome_sizes: str | Path | None = None,
    write_cutsite_bigwigs: bool = False,
    cpm_normalize: bool = True,
    cores: int | None = None,
) -> pd.DataFrame:
    """Split fragments into pseudobulk group files and write a manifest/QC table."""

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    if index_output and not compress_output:
        compress_output = True
    if write_cutsite_bigwigs and genome_sizes is None:
        raise ValueError("--genome-sizes is required with --write-cutsite-bigwigs")
    cores = multiprocessing.cpu_count() if cores is None else max(1, int(cores))

    barcode_to_group = load_annotations(annotations, barcode_column, group_by, strip_barcode_suffix=strip_barcode_suffix)
    handles = {}
    cells_by_group: dict[str, set[str]] = defaultdict(set)
    fragments_by_group: dict[str, int] = defaultdict(int)
    paths_by_group: dict[str, Path] = {}
    matched_fragments = 0
    skipped_unmatched = 0
    skipped_chrom = 0

    try:
        with _open_text(fragments) as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 4:
                    continue
                chrom = fields[0]
                if include_chroms is not None and chrom not in include_chroms:
                    skipped_chrom += 1
                    continue
                if exclude_chroms is not None and chrom in exclude_chroms:
                    skipped_chrom += 1
                    continue
                barcode = _normalize_barcode(fields[3], strip_barcode_suffix)
                group = barcode_to_group.get(barcode)
                if group is None:
                    skipped_unmatched += 1
                    continue
                if group not in handles:
                    path = outdir / f"{group}.fragments.tsv"
                    paths_by_group[group] = path
                    handles[group] = path.open("w", encoding="utf-8")
                handles[group].write(line)
                cells_by_group[group].add(barcode)
                count = int(fields[4]) if len(fields) > 4 and fields[4].isdigit() else 1
                fragments_by_group[group] += count
                matched_fragments += count
    finally:
        for handle in handles.values():
            handle.close()

    if compress_output:
        tasks = [(group, str(path), index_output) for group, path in paths_by_group.items()]
        if cores > 1 and len(tasks) > 1:
            with ProcessPoolExecutor(max_workers=min(cores, len(tasks))) as executor:
                for group, compressed in executor.map(_compress_group, tasks):
                    paths_by_group[group] = Path(compressed)
        else:
            for task in tasks:
                group, compressed = _compress_group(task)
                paths_by_group[group] = Path(compressed)

    keep_by_group: dict[str, bool] = {}
    bigwigs_by_group: dict[str, str] = defaultdict(str)
    for group in sorted(paths_by_group):
        cells = len(cells_by_group[group])
        fragments_count = fragments_by_group[group]
        keep_by_group[group] = cells >= min_cells and fragments_count >= min_fragments

    if write_cutsite_bigwigs:
        tasks = [
            (
                group,
                str(paths_by_group[group]),
                str(outdir / f"{group}.cutsites.cpm.bw"),
                str(genome_sizes),
                cpm_normalize,
            )
            for group in sorted(paths_by_group)
            if keep_by_group[group]
        ]
        if cores > 1 and len(tasks) > 1:
            with ProcessPoolExecutor(max_workers=min(cores, len(tasks))) as executor:
                for group, bigwig in executor.map(_write_group_bigwig, tasks):
                    bigwigs_by_group[group] = bigwig
        else:
            for task in tasks:
                group, bigwig = _write_group_bigwig(task)
                bigwigs_by_group[group] = bigwig

    rows = []
    for group in sorted(paths_by_group):
        rows.append(
            {
                "group": group,
                "fragment_file": str(paths_by_group[group]),
                "n_cells": len(cells_by_group[group]),
                "n_fragments": fragments_by_group[group],
                "cutsite_bigwig": bigwigs_by_group[group],
                "passes_filters": keep_by_group[group],
            }
        )
    manifest = pd.DataFrame(rows)
    manifest_path = outdir / "pseudobulk_manifest.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False)

    config = {
        "version": 1,
        "pseudobulk_manifest": str(manifest_path),
        "group_by": group_by,
        "matched_fragments": matched_fragments,
        "skipped_unmatched_fragments": skipped_unmatched,
        "skipped_chromosome_fragments": skipped_chrom,
        "samples": [row for row in manifest.to_dict(orient="records") if row["passes_filters"]],
    }
    with (outdir / "fp_tools_manifest.yml").open("w", encoding="utf-8") as handle:
        yaml.safe_dump(config, handle, sort_keys=False)
    return manifest


def write_downstream_commands(
    manifest: pd.DataFrame,
    output: str | Path,
    genome_sizes: str | Path | None = None,
    cores: int | None = None,
) -> Path:
    """Write a reproducible shell plan for pseudobulk BAM/bigWig generation."""

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    cores = max(1, int(cores))
    genome_arg = shlex.quote(str(genome_sizes)) if genome_sizes else "GENOME_SIZES.txt"
    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        "# Generated by pseudobulk-fragments. Review paths and install bedtools, samtools, and bedGraphToBigWig before running.",
    ]
    if genome_sizes is None:
        lines.append("# Replace GENOME_SIZES.txt with a two-column chromosome sizes file before running bigWig/BAM steps.")
    lines.append("")

    runnable = manifest[manifest["passes_filters"].astype(bool)] if "passes_filters" in manifest.columns else manifest
    for _, row in runnable.iterrows():
        group = str(row["group"])
        fragment_file = Path(str(row["fragment_file"]))
        prefix = fragment_file.with_suffix("")
        bed = prefix.with_suffix(".bed")
        bedgraph = prefix.with_suffix(".bedGraph")
        unsorted_bam = prefix.with_suffix(".bam")
        sorted_bam = prefix.with_suffix(".sorted.bam")
        bigwig = prefix.with_suffix(".bw")
        lines.extend(
            [
                f"# {group}",
                f"awk 'BEGIN{{OFS=\"\\t\"}} !/^#/ {{print $1,$2,$3}}' {shlex.quote(str(fragment_file))} | sort -k1,1 -k2,2n > {shlex.quote(str(bed))}",
                f"bedtools genomecov -bg -i {shlex.quote(str(bed))} -g {genome_arg} > {shlex.quote(str(bedgraph))}",
                f"bedGraphToBigWig {shlex.quote(str(bedgraph))} {genome_arg} {shlex.quote(str(bigwig))}",
                f"bedToBam -i {shlex.quote(str(bed))} -g {genome_arg} > {shlex.quote(str(unsorted_bam))}",
                f"samtools sort -@ {cores} -o {shlex.quote(str(sorted_bam))} {shlex.quote(str(unsorted_bam))}",
                f"samtools index -@ {cores} {shlex.quote(str(sorted_bam))}",
                "",
            ]
        )
    output.write_text("\n".join(lines), encoding="utf-8")
    output.chmod(0o755)
    return output


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Group single-cell ATAC fragments into pseudobulk fragment files.")
    parser.add_argument("--fragments", required=True, help="10x-style fragments TSV/TSV.GZ with barcode in column 4.")
    parser.add_argument("--annotations", required=True, help="Cell annotation TSV or CSV.")
    parser.add_argument("--group-by", required=True, help="Comma-separated annotation columns to group by, e.g. donor,cell_type.")
    parser.add_argument("--barcode-column", default="barcode", help="Annotation barcode column (default: barcode).")
    parser.add_argument("--no-strip-barcode-suffix", action="store_true", help="Require exact barcode matches instead of matching AAAC-1 to AAAC.")
    parser.add_argument("--include-chroms", help="Comma-separated chromosomes to keep, e.g. chr1,chr2,chrX.")
    parser.add_argument("--exclude-chroms", help="Comma-separated chromosomes to skip, e.g. chrM,chrY.")
    parser.add_argument("--min-cells", type=int, default=1, help="Minimum cells for passes_filters (default: 1).")
    parser.add_argument("--min-fragments", type=int, default=1, help="Minimum fragments for passes_filters (default: 1).")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument("--compress-output", action="store_true", help="Write grouped fragments as .tsv.gz files.")
    parser.add_argument("--index-output", action="store_true", help="BGZF-compress and tabix-index grouped fragments for random access.")
    parser.add_argument("--write-cutsite-bigwigs", action="store_true", help="Write one sparse cut-site bigWig per kept pseudobulk group.")
    parser.add_argument("--no-cpm-normalize", action="store_true", help="Write raw cut counts instead of CPM-normalized bigWig values.")
    parser.add_argument("--write-downstream-commands", action="store_true", help="Write a shell script for BED/BAM/bigWig generation from kept pseudobulk groups.")
    parser.add_argument("--genome-sizes", help="Two-column chromosome sizes file used by generated bedtools/UCSC commands and cut-site bigWigs.")
    parser.add_argument("--cores", type=int, default=None, help="Cores for compression, bigWig writing, and generated samtools commands (default: all available cores).")
    args = parser.parse_args(argv)

    manifest = group_fragments(
        args.fragments,
        args.annotations,
        args.outdir,
        group_by=[column.strip() for column in args.group_by.split(",") if column.strip()],
        barcode_column=args.barcode_column,
        min_cells=args.min_cells,
        min_fragments=args.min_fragments,
        strip_barcode_suffix=not args.no_strip_barcode_suffix,
        include_chroms=_parse_chrom_list(args.include_chroms),
        exclude_chroms=_parse_chrom_list(args.exclude_chroms),
        compress_output=args.compress_output,
        index_output=args.index_output,
        genome_sizes=args.genome_sizes,
        write_cutsite_bigwigs=args.write_cutsite_bigwigs,
        cpm_normalize=not args.no_cpm_normalize,
        cores=args.cores,
    )
    print(f"Wrote {len(manifest)} pseudobulk groups to {args.outdir}")
    if args.write_downstream_commands:
        command_path = write_downstream_commands(
            manifest,
            Path(args.outdir) / "pseudobulk_downstream_commands.sh",
            genome_sizes=args.genome_sizes,
            cores=args.cores,
        )
        print(f"Wrote downstream command plan to {command_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
