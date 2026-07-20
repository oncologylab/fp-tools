#!/usr/bin/env python3
"""Build versioned multimodal metadata for the LCMV CD8 ATAC/RNA analysis."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in fields} for row in rows)


def joined(values) -> str:
    return ",".join(sorted({str(value) for value in values if str(value)}))


def slug(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_")


def load_atac_outputs(paths: list[Path]) -> dict[str, dict[str, str]]:
    outputs: dict[str, dict[str, str]] = {}
    for path in paths:
        for row in read_tsv(path):
            sample = row["sample"]
            if sample in outputs:
                raise ValueError(f"Duplicate ATAC output sample {sample} in {path}")
            outputs[sample] = row
    return outputs


def analysis_units(atac: list[dict[str, str]]) -> list[list[dict[str, str]]]:
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in atac:
        grouped[(row["condition_pair_id"], row["replicate"])].append(row)
    units: list[list[dict[str, str]]] = []
    for group in grouped.values():
        partitions = [row["technical_partition"] for row in group]
        if len(group) > 1 and all(partitions):
            units.append(sorted(group, key=lambda row: row["technical_partition"]))
        else:
            units.extend([[row] for row in group])
    return sorted(
        units,
        key=lambda group: (
            group[0]["collection"],
            group[0]["author"],
            group[0]["condition_pair_id"],
            int(group[0]["replicate"]),
        ),
    )


def expected_atac_output(root: Path, run: str) -> dict[str, str]:
    sample = root / "samples" / run
    return {
        "sample": run,
        "bam": str(sample / "alignment" / f"{run}.filtered.bam"),
        "peaks": str(sample / "peaks" / f"{run}.narrowPeak"),
        "bigwig": str(sample / "tracks" / f"{run}.rp10m.bw"),
    }


def validate_comparisons(
    comparisons: list[dict[str, str]], resolved: list[dict[str, str]]
) -> None:
    conditions: dict[str, set[str]] = defaultdict(set)
    authors: dict[str, set[str]] = defaultdict(set)
    for row in resolved:
        conditions[row["condition_pair_id"]].add(row["assay"])
        authors[row["condition_pair_id"]].add(row["author"])
    seen: set[str] = set()
    for row in comparisons:
        name = row["comparison"]
        if name in seen:
            raise ValueError(f"Duplicate v2 comparison: {name}")
        seen.add(name)
        required = set(row["assay_scope"].split("+"))
        for condition in (row["cond1"], row["cond2"]):
            if condition not in conditions:
                raise ValueError(f"Unknown condition {condition} in {name}")
            if not required.issubset(conditions[condition]):
                raise ValueError(
                    f"{name} requires {sorted(required)} but {condition} has "
                    f"{sorted(conditions[condition])}"
                )
            if authors[condition] != {row["author"]}:
                raise ValueError(f"Cross-study or mislabeled comparison: {name}")


def build(
    resolved_path: Path,
    project: Path,
    comparison_path: Path,
    atac_output_paths: list[Path],
    new_atac_root: Path,
) -> dict[str, int]:
    resolved = read_tsv(resolved_path)
    comparisons = read_tsv(comparison_path)
    validate_comparisons(comparisons, resolved)
    outputs = load_atac_outputs(atac_output_paths)
    atac = [row for row in resolved if row["assay"] == "ATAC"]
    rna = [row for row in resolved if row["assay"] == "RNA"]
    rna_by_condition: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rna:
        rna_by_condition[row["condition_pair_id"]].append(row)

    units = analysis_units(atac)
    inventory_rows: list[dict[str, object]] = []
    primary_fp: list[dict[str, object]] = []
    supporting_fp: list[dict[str, object]] = []
    descriptive_fp: list[dict[str, object]] = []
    for group in units:
        first = group[0]
        merged = len(group) > 1
        if merged:
            sample = (
                f"{slug(first['author'])}_{first['harmonized_condition']}_"
                f"rep{first['replicate']}"
            )
            merge = project / "inputs" / "technical_merges" / sample
            output = {
                "bam": str(merge / f"{sample}.filtered.bam"),
                "peaks": str(merge / f"{sample}.narrowPeak"),
                "bigwig": "",
            }
        else:
            sample = first["gsm_accession"]
            output = outputs.get(first["run_accession"])
            if output is None:
                output = expected_atac_output(new_atac_root, first["run_accession"])
        rna_matches = rna_by_condition.get(first["condition_pair_id"], [])
        inventory_rows.append(
            {
                "ID": sample,
                "Author": first["author"],
                "Year": first["year"],
                "Series": first["series"],
                "Subseries": first["subseries"],
                "Infection": first["infection"],
                "Day": first["day"],
                "Tissue": first["tissue"],
                "Original_condition": joined(item["original_condition"] for item in group),
                "Harmonized_condition": first["harmonized_condition"],
                "Broader_state": first["broader_state"],
                "Condition_pair_ID": first["condition_pair_id"],
                "Biological_replicate": first["replicate"],
                "ATAC_GSM": joined(item["gsm_accession"] for item in group),
                "ATAC_SRR": joined(item["run_accession"] for item in group),
                "RNA_GSM": joined(item["gsm_accession"] for item in rna_matches),
                "RNA_SRR": joined(item["run_accession"] for item in rna_matches),
                "RNA_biological_replicates": joined(item["replicate"] for item in rna_matches),
                "RNA_match_scope": first["rna_match_status"],
                "Pairing_level": "condition" if rna_matches else "none",
                "Same_biological_specimen": "false",
                "Collection": first["collection"],
                "Include_primary_paired_analysis": first["include_in_primary_paired_analysis"],
                "ATAC_BAM": output["bam"],
                "ATAC_peaks": output["peaks"],
                "ATAC_bigwig": output["bigwig"],
            }
        )
        fp_row = {
            "sample": sample,
            "condition": first["condition_pair_id"],
            "bam": output["bam"],
            "peaks": output["peaks"],
        }
        descriptive_fp.append(
            {**fp_row, "broad_condition": first["broader_state"], "collection": first["collection"]}
        )
        (primary_fp if first["collection"] == "primary" else supporting_fp).append(fp_row)

    metadata = project / "metadata"
    write_tsv(metadata / "multimodal_inventory.tsv", inventory_rows, list(inventory_rows[0]))
    sample_fields = ["sample", "condition", "bam", "peaks"]
    write_tsv(metadata / "samples.tsv", primary_fp, sample_fields)
    write_tsv(metadata / "supporting_samples.tsv", supporting_fp, sample_fields)
    write_tsv(
        metadata / "supporting_analysis_samples.tsv",
        [*primary_fp, *supporting_fp],
        sample_fields,
    )
    write_tsv(
        metadata / "samples_descriptive.tsv",
        descriptive_fp,
        [*sample_fields, "broad_condition", "collection"],
    )

    primary = [row for row in comparisons if row["analysis_tier"].startswith("primary_")]
    supporting_atac = [row for row in comparisons if row["analysis_tier"] == "supporting_atac_only"]
    supporting_rna = [row for row in comparisons if row["analysis_tier"] == "supporting_rna_only"]
    comparison_fields = list(comparisons[0])
    write_tsv(metadata / "comparison_design.tsv", comparisons, comparison_fields)
    write_tsv(metadata / "comparisons.tsv", primary, comparison_fields)
    write_tsv(metadata / "rna_comparisons.tsv", primary, comparison_fields)
    write_tsv(metadata / "supporting_atac_comparisons.tsv", supporting_atac, comparison_fields)
    write_tsv(metadata / "supporting_rna_comparisons.tsv", supporting_rna, comparison_fields)

    rna_by_gsm: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rna:
        rna_by_gsm[row["gsm_accession"]].append(row)
    rna_rows: list[dict[str, object]] = []
    for gsm, group in sorted(rna_by_gsm.items()):
        first = group[0]
        rna_rows.append(
            {
                "sample": gsm,
                "condition": first["condition_pair_id"],
                "broad_condition": first["broader_state"],
                "collection": first["collection"],
                "author": first["author"],
                "series": first["series"],
                "gsm_accession": gsm,
                "srr_accessions": joined(item["run_accession"] for item in group),
                "replicate": first["replicate"],
                "library_layout": first["library_layout"],
                "fastq_1": ";".join(item["fastq_1"] for item in group),
                "fastq_2": ";".join(item["fastq_2"] for item in group if item["fastq_2"]),
                "paper_method": "tophat2_htseq" if first["author"] in {"Milner", "Scott-Browne"} else "kallisto",
            }
        )
    rna_metadata = project / "rna" / "metadata"
    write_tsv(rna_metadata / "samples.tsv", rna_rows, list(rna_rows[0]))
    write_tsv(
        rna_metadata / "samples_fine.tsv",
        [
            {"sample": row["sample"], "condition": row["condition"]}
            for row in rna_rows
            if row["collection"] == "primary"
        ],
        ["sample", "condition"],
    )
    write_tsv(
        rna_metadata / "samples_descriptive.tsv",
        [
            {
                "sample": row["sample"],
                "condition": row["broad_condition"],
                "collection": row["collection"],
            }
            for row in rna_rows
        ],
        ["sample", "condition", "collection"],
    )

    counts = {
        "atac_runs": len(atac),
        "atac_units": len(units),
        "primary_atac_units": len(primary_fp),
        "supporting_atac_units": len(supporting_fp),
        "rna_samples": len(rna_rows),
        "primary_rna_samples": sum(row["collection"] == "primary" for row in rna_rows),
        "primary_comparisons": len(primary),
        "supporting_atac_comparisons": len(supporting_atac),
        "supporting_rna_comparisons": len(supporting_rna),
    }
    print("Built LCMV v2 metadata: " + ", ".join(f"{key}={value}" for key, value in counts.items()))
    return counts


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolved-runs", type=Path, required=True)
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--comparisons", type=Path, required=True)
    parser.add_argument("--atac-outputs", type=Path, action="append", default=[])
    parser.add_argument("--new-atac-root", type=Path, required=True)
    args = parser.parse_args()
    build(
        args.resolved_runs.resolve(),
        args.project.resolve(),
        args.comparisons.resolve(),
        [path.resolve() for path in args.atac_outputs],
        args.new_atac_root.resolve(),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
