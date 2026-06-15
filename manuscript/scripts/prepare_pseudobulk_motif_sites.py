#!/usr/bin/env python
"""Scan exact motif-centered sites for the 10x PBMC pseudobulk example."""

from __future__ import annotations

import argparse
from collections import defaultdict
import heapq
from pathlib import Path

import pysam

from fp_tools.utils.motifs import MotifList
from fp_tools.utils.regions import OneRegion


DEFAULT_CANDIDATES = {
    "B_cell": ["PAX5", "EBF1", "POU2F2", "SPIB", "MEF2C", "BHLHE41", "BACH2"],
    "T_NK": ["TCF7", "TCF7L2", "LEF1", "Gata3", "RUNX3", "RORA", "TBX21", "EOMES", "IRF4"],
    "Myeloid": ["Spi1", "CEBPB", "CEBPA", "CEBPD", "JUNB", "FOS", "BATF", "MAF", "REL"],
    "Control": ["CTCF"],
}


def parse_candidates(text: str | None) -> dict[str, list[str]]:
    if not text:
        return DEFAULT_CANDIDATES
    out: dict[str, list[str]] = {}
    for block in text.split(";"):
        if not block.strip():
            continue
        label, values = block.split(":", 1)
        out[label.strip()] = [value.strip() for value in values.split(",") if value.strip()]
    return out


def load_peaks(path: Path, chroms: set[str] | None, max_peaks: int | None) -> list[OneRegion]:
    peaks = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, *_ = line.rstrip("\n").split("\t")
            if chroms is not None and chrom not in chroms:
                continue
            peaks.append(OneRegion([chrom, int(start), int(end)]))
            if max_peaks is not None and len(peaks) >= max_peaks:
                break
    return peaks


def select_motifs(motif_file: Path, candidates: dict[str, list[str]], pvalue: float, naming: str) -> tuple[MotifList, dict[str, str], dict[str, str]]:
    wanted = {name.upper(): lineage for lineage, names in candidates.items() for name in names}
    selected = MotifList()
    prefix_to_tf: dict[str, str] = {}
    prefix_to_lineage: dict[str, str] = {}
    for motif in MotifList().from_file(str(motif_file)):
        if motif.name.upper() not in wanted:
            continue
        motif.set_prefix(naming)
        motif.get_threshold(pvalue)
        selected.append(motif)
        prefix_to_tf[motif.prefix] = motif.name
        prefix_to_lineage[motif.prefix] = wanted[motif.name.upper()]
    selected.set_background()
    return selected, prefix_to_tf, prefix_to_lineage


def scan_sites(
    peaks: list[OneRegion],
    genome: Path,
    motifs: MotifList,
    prefix_to_tf: dict[str, str],
    prefix_to_lineage: dict[str, str],
    outdir: Path,
    plot_sites_per_tf: int,
    site_selection: str,
) -> dict[str, dict[str, int | float | str]]:
    outdir.mkdir(parents=True, exist_ok=True)
    motifs.setup_moods_scanner()
    fasta = pysam.FastaFile(str(genome))
    total_counts: dict[str, int] = defaultdict(int)
    motif_counts: dict[str, int] = defaultdict(int)
    seen: dict[str, set[tuple[str, int, int, str, str]]] = defaultdict(set)
    selected: dict[str, list[tuple[float, str, int, int, str, str]]] = defaultdict(list)
    metadata: dict[str, dict[str, int | float | str]] = {}
    try:
        chrom_sizes = dict(zip(fasta.references, fasta.lengths))
        for region in peaks:
            if region.chrom not in chrom_sizes or region.end > chrom_sizes[region.chrom]:
                continue
            seq = fasta.fetch(region.chrom, region.start, region.end)
            for site in motifs.scan_sequence(seq, region):
                tf = prefix_to_tf.get(site.name)
                if tf is None:
                    continue
                key = (site.chrom, int(site.start), int(site.end), str(site.strand), site.name)
                if key in seen[tf]:
                    continue
                seen[tf].add(key)
                total_counts[tf] += 1
                motif_counts[site.name] += 1
                score = float(site.score)
                item = (score, site.chrom, int(site.start), int(site.end), site.name, str(site.strand))
                if site_selection == "all" or plot_sites_per_tf <= 0:
                    selected[tf].append(item)
                else:
                    heap = selected[tf]
                    if len(heap) < plot_sites_per_tf:
                        heapq.heappush(heap, item)
                    elif score > heap[0][0]:
                        heapq.heapreplace(heap, item)
                metadata[tf] = {
                    "tf": tf,
                    "lineage": prefix_to_lineage[site.name],
                    "motif_prefixes": ",".join(sorted(prefix for prefix, value in prefix_to_tf.items() if value == tf)),
                }
    finally:
        fasta.close()

    for tf, items in selected.items():
        items_sorted = sorted(items, key=lambda row: (-row[0], row[1], row[2], row[3], row[4], row[5]))
        with (outdir / f"{tf}.motif_hits.bed").open("w", encoding="utf-8") as handle:
            for score, chrom, start, end, motif_name, strand in items_sorted:
                handle.write(f"{chrom}\t{start}\t{end}\t{motif_name}\t{score}\t{strand}\n")

    for tf, row in metadata.items():
        plotted_sites = len(selected.get(tf, []))
        scores = [item[0] for item in selected.get(tf, [])]
        row['total_motif_hits'] = total_counts[tf]
        row['plotted_sites'] = plotted_sites
        row['n_sites'] = plotted_sites
        row['selection_method'] = "all" if site_selection == "all" or plot_sites_per_tf <= 0 else f"top_score_{plot_sites_per_tf}"
        row['score_min_selected'] = min(scores) if scores else float("nan")
        row['motif_site_counts'] = ",".join(f"{prefix}:{motif_counts[prefix]}" for prefix in str(row['motif_prefixes']).split(","))
    return metadata


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--peaks", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--motifs", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--summary", default=None)
    parser.add_argument("--candidates", default=None, help="Semicolon-delimited lineage:TF,TF specification.")
    parser.add_argument("--chroms", default="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX")
    parser.add_argument("--max-peaks", type=int, default=None)
    parser.add_argument("--plot-sites-per-tf", type=int, default=2000, help="Number of top-scoring motif hits to write per TF for plotting (default: 2000; <=0 writes all hits).")
    parser.add_argument("--site-selection", choices=["top-score", "all"], default="top-score", help="How to choose plotted motif hits after counting all hits (default: top-score).")
    parser.add_argument("--max-sites-per-tf", type=int, default=None, help="Deprecated alias for --plot-sites-per-tf.")
    parser.add_argument("--motif-pvalue", type=float, default=1e-4)
    parser.add_argument("--naming", choices=["id", "name", "name_id", "id_name"], default="name_id")
    args = parser.parse_args(argv)

    candidates = parse_candidates(args.candidates)
    chroms = {chrom.strip() for chrom in args.chroms.split(",") if chrom.strip()} if args.chroms else None
    peaks = load_peaks(Path(args.peaks), chroms, args.max_peaks)
    motifs, prefix_to_tf, prefix_to_lineage = select_motifs(Path(args.motifs), candidates, args.motif_pvalue, args.naming)
    if not motifs:
        raise SystemExit("No requested motifs were found in the motif file.")
    plot_sites = args.max_sites_per_tf if args.max_sites_per_tf is not None else args.plot_sites_per_tf
    selection = "top-score" if args.site_selection == "top-score" else "all"
    metadata = scan_sites(peaks, Path(args.genome), motifs, prefix_to_tf, prefix_to_lineage, Path(args.outdir), plot_sites, selection)

    summary_path = Path(args.summary) if args.summary else Path(args.outdir) / "motif_centered_site_summary.tsv"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with summary_path.open("w", encoding="utf-8") as handle:
        handle.write(
            "tf\tlineage\ttotal_motif_hits\tplotted_sites\tn_sites\tselection_method\t"
            "score_min_selected\tmotif_prefixes\tmotif_site_counts\n"
        )
        for tf in sorted(metadata):
            row = metadata[tf]
            handle.write(
                f"{tf}\t{row['lineage']}\t{row['total_motif_hits']}\t{row['plotted_sites']}\t"
                f"{row['n_sites']}\t{row['selection_method']}\t{row['score_min_selected']}\t"
                f"{row['motif_prefixes']}\t{row['motif_site_counts']}\n"
            )
    print(f"Wrote motif-centered sites for {len(metadata)} TFs to {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
