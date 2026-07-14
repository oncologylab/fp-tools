#!/usr/bin/env python3
"""Assemble LCMV RNA quantifications into gene-level count matrices."""

from __future__ import annotations

import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path


def attributes(value: str) -> dict[str, str]:
    result = {}
    for item in value.rstrip(";").split(";"):
        fields = item.strip().split(" ", 1)
        if len(fields) == 2:
            result[fields[0]] = fields[1].strip().strip('"')
    return result


def transcript_map(gtf: Path) -> tuple[dict[str, tuple[str, str]], dict[str, str]]:
    opener = gzip.open if gtf.suffix == ".gz" else open
    tx, genes = {}, {}
    with opener(gtf, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            attr = attributes(fields[8])
            gene = attr.get("gene_id", "").split(".")[0]
            symbol = attr.get("gene_name", gene)
            if gene:
                genes[gene] = symbol
            transcript = attr.get("transcript_id", "").split(".")[0]
            if transcript and gene:
                tx[transcript] = (gene, symbol)
    return tx, genes


def write_tx2gene(gtf: Path, path: Path) -> None:
    """Write the versioned transcript mapping expected by tximport."""
    opener = gzip.open if gtf.suffix == ".gz" else open
    seen = set()
    rows = []
    with opener(gtf, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "transcript":
                continue
            attr = attributes(fields[8])
            transcript = attr.get("transcript_id", "")
            gene = attr.get("gene_id", "").split(".")[0]
            if transcript and gene and transcript not in seen:
                seen.add(transcript)
                rows.append((transcript, gene, attr.get("gene_name", gene)))
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["transcript_id", "gene_id", "gene_symbol"])
        writer.writerows(rows)


def read_samples(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_matrix(path: Path, values: dict[str, dict[str, float]], samples: list[str], genes: dict[str, str], rounded: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id", "gene_symbol", *samples])
        for gene in sorted(values):
            row = [values[gene].get(sample, 0.0) for sample in samples]
            writer.writerow([gene, genes.get(gene, gene), *([int(round(x)) for x in row] if rounded else [f"{x:.6f}" for x in row])])


def kallisto(project: Path, rows: list[dict[str, str]], txmap: dict[str, tuple[str, str]], genes: dict[str, str]) -> None:
    sample_ids = [row["sample"] for row in rows]
    counts, tpms = defaultdict(dict), defaultdict(dict)
    for sample in sample_ids:
        abundance = project / "rna/uniform_kallisto" / sample / "abundance.tsv"
        if not abundance.is_file():
            raise FileNotFoundError(abundance)
        with abundance.open(encoding="utf-8", newline="") as handle:
            for record in csv.DictReader(handle, delimiter="\t"):
                key = record["target_id"].split(".")[0]
                if key not in txmap:
                    continue
                gene, symbol = txmap[key]
                genes.setdefault(gene, symbol)
                counts[gene][sample] = counts[gene].get(sample, 0.0) + float(record["est_counts"])
                tpms[gene][sample] = tpms[gene].get(sample, 0.0) + float(record["tpm"])
    out = project / "rna/counts/uniform_kallisto"
    write_matrix(out / "gene_estimated_counts.tsv", counts, sample_ids, genes)
    write_matrix(out / "gene_counts_rounded.tsv", counts, sample_ids, genes, rounded=True)
    write_matrix(out / "gene_tpm.tsv", tpms, sample_ids, genes)


def paper(project: Path, rows: list[dict[str, str]], txmap: dict[str, tuple[str, str]], genes: dict[str, str]) -> None:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["author"]].append(row)
    for author, subset in grouped.items():
        sample_ids = [row["sample"] for row in subset]
        values = defaultdict(dict)
        for row in subset:
            sample = row["sample"]
            if row["paper_method"] == "tophat2_htseq":
                source = project / "rna/paper_specific" / sample / "htseq_counts.tsv"
                with source.open(encoding="utf-8") as handle:
                    for line in handle:
                        gene, count = line.rstrip("\n").split("\t")[:2]
                        gene = gene.split(".")[0]
                        if not gene.startswith("__"):
                            values[gene][sample] = float(count)
            else:
                source = project / "rna/paper_specific" / sample / "abundance.tsv"
                with source.open(encoding="utf-8", newline="") as handle:
                    for record in csv.DictReader(handle, delimiter="\t"):
                        mapped = txmap.get(record["target_id"].split(".")[0])
                        if mapped:
                            gene, symbol = mapped
                            genes.setdefault(gene, symbol)
                            values[gene][sample] = values[gene].get(sample, 0.0) + float(record["est_counts"])
        write_matrix(project / "rna/counts/paper_specific" / f"{author.lower().replace('-', '_')}_gene_counts.tsv", values, sample_ids, genes, rounded=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    args = parser.parse_args()
    txmap, genes = transcript_map(args.gtf)
    rows = read_samples(args.project / "rna/metadata/samples.tsv")
    write_tx2gene(args.gtf, args.project / "rna/metadata/tx2gene.tsv")
    kallisto(args.project, rows, txmap, genes)
    paper(args.project, rows, txmap, genes)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
