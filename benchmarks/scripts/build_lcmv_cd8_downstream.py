#!/usr/bin/env python3
"""Build multimodal metadata for the LCMV CD8 ATAC/RNA downstream analysis."""

from __future__ import annotations

import argparse
import csv
import itertools
from collections import defaultdict
from pathlib import Path


FINE_COMPARISONS = [
    ("guan_MPEC_vs_NAIVE", "guan_MPEC_D8", "guan_NAIVE_D0"),
    ("guan_SLEC_vs_NAIVE", "guan_SLEC_D8", "guan_NAIVE_D0"),
    ("guan_SLEC_vs_MPEC", "guan_SLEC_D8", "guan_MPEC_D8"),
    ("scott_browne_MPEC_vs_NAIVE", "scott_browne_MPEC_D8_ATAC_ONLY", "scott_browne_NAIVE_D0"),
    ("scott_browne_SLEC_vs_NAIVE", "scott_browne_SLEC_D8_ATAC_ONLY", "scott_browne_NAIVE_D0"),
    ("scott_browne_SLEC_vs_MPEC", "scott_browne_SLEC_D8_ATAC_ONLY", "scott_browne_MPEC_D8_ATAC_ONLY"),
    ("beltra_TexProg2_vs_TexProg1", "beltra_TexProg2", "beltra_TexProg1"),
    ("beltra_TexInt_vs_TexProg1", "beltra_TexInt", "beltra_TexProg1"),
    ("beltra_TexTerm_vs_TexProg1", "beltra_TexTerm", "beltra_TexProg1"),
    ("beltra_TexInt_vs_TexProg2", "beltra_TexInt", "beltra_TexProg2"),
    ("beltra_TexTerm_vs_TexProg2", "beltra_TexTerm", "beltra_TexProg2"),
    ("beltra_TexTerm_vs_TexInt", "beltra_TexTerm", "beltra_TexInt"),
]
RNA_FINE_COMPARISONS = [row for row in FINE_COMPARISONS if not row[0].startswith("scott_browne_")]
BROAD_STATES = ["NAIVE", "MPEC", "SLEC", "TRM", "TEXprog", "TEXeff", "TEXterm"]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def joined(values) -> str:
    return ",".join(sorted({str(value) for value in values if str(value)}))


def build(root: Path, project: Path) -> None:
    metadata = root / "metadata"
    resolved = read_tsv(metadata / "resolved_runs.tsv")
    outputs = {row["sample"]: row for row in read_tsv(metadata / "all_atac_samples.tsv")}
    atac = [row for row in resolved if row["assay"] == "ATAC"]
    rna = [row for row in resolved if row["assay"] == "RNA"]
    rna_by_condition: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rna:
        rna_by_condition[row["condition_pair_id"]].append(row)

    merge_srrs = {"SRR4435505", "SRR4435506"}
    units: list[list[dict[str, str]]] = [[row for row in atac if row["run_accession"] in merge_srrs]]
    units.extend([[row] for row in atac if row["run_accession"] not in merge_srrs])

    paired_rows: list[dict[str, object]] = []
    fp_rows: list[dict[str, object]] = []
    broad_rows: list[dict[str, object]] = []
    for group in units:
        first = group[0]
        is_merged = len(group) > 1
        sample = "ScottBrowne_NAIVE_D0_rep1" if is_merged else first["gsm_accession"]
        if is_merged:
            bam = project / "inputs" / "technical_merges" / sample / f"{sample}.filtered.bam"
            peaks = project / "inputs" / "technical_merges" / sample / f"{sample}.narrowPeak"
            bigwig = ""
        else:
            output = outputs[first["run_accession"]]
            bam, peaks, bigwig = output["bam"], output["peaks"], output["bigwig"]
        rna_matches = rna_by_condition.get(first["condition_pair_id"], [])
        row = {
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
            "ATAC_BAM": str(bam),
            "ATAC_peaks": str(peaks),
            "ATAC_bigwig": str(bigwig),
        }
        paired_rows.append(row)
        fp_rows.append({"sample": sample, "condition": first["condition_pair_id"], "bam": str(bam), "peaks": str(peaks)})
        broad_rows.append({"sample": sample, "condition": first["broader_state"], "bam": str(bam), "peaks": str(peaks)})

    write_tsv(root / "LCMV_CD8_mm10_ATAC_RNA.txt", paired_rows, list(paired_rows[0]))
    condition_rows = []
    for condition in sorted({row["condition_pair_id"] for row in resolved}):
        atac_group = [row for row in atac if row["condition_pair_id"] == condition]
        rna_group = rna_by_condition.get(condition, [])
        first = (atac_group or rna_group)[0]
        condition_rows.append({
            "Condition_pair_ID": condition,
            "Author": first["author"],
            "Harmonized_condition": first["harmonized_condition"],
            "Broader_state": first["broader_state"],
            "Collection": first["collection"],
            "ATAC_GSM": joined(row["gsm_accession"] for row in atac_group),
            "ATAC_SRR": joined(row["run_accession"] for row in atac_group),
            "RNA_GSM": joined(row["gsm_accession"] for row in rna_group),
            "RNA_SRR": joined(row["run_accession"] for row in rna_group),
            "Pairing_level": "condition" if atac_group and rna_group else "none",
            "Same_biological_specimen": "false",
            "Include_primary_paired_analysis": first["include_in_primary_paired_analysis"],
        })
    write_tsv(root / "metadata" / "condition_assay_pairs.tsv", condition_rows, list(condition_rows[0]))
    write_tsv(project / "metadata" / "samples.tsv", fp_rows, ["sample", "condition", "bam", "peaks"])
    write_tsv(project / "metadata" / "samples_broad.tsv", broad_rows, ["sample", "condition", "bam", "peaks"])
    write_tsv(project / "metadata" / "comparisons.tsv", [{"comparison": a, "cond1": b, "cond2": c} for a, b, c in FINE_COMPARISONS], ["comparison", "cond1", "cond2"])
    write_tsv(project / "metadata" / "rna_comparisons.tsv", [{"comparison": a, "cond1": b, "cond2": c} for a, b, c in RNA_FINE_COMPARISONS], ["comparison", "cond1", "cond2"])
    broad_comparisons = [{"comparison": f"{b}_vs_{a}", "cond1": b, "cond2": a} for a, b in itertools.combinations(BROAD_STATES, 2)]
    write_tsv(project / "metadata" / "comparisons_broad.tsv", broad_comparisons, ["comparison", "cond1", "cond2"])

    rna_by_gsm: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rna:
        rna_by_gsm[row["gsm_accession"]].append(row)
    rna_rows = []
    for gsm, group in sorted(rna_by_gsm.items()):
        first = group[0]
        rna_rows.append({
            "sample": gsm,
            "condition": first["condition_pair_id"],
            "broad_condition": first["broader_state"],
            "author": first["author"],
            "series": first["series"],
            "gsm_accession": gsm,
            "srr_accessions": joined(item["run_accession"] for item in group),
            "replicate": first["replicate"],
            "library_layout": first["library_layout"],
            "fastq_1": ";".join(item["fastq_1"] for item in group),
            "fastq_2": ";".join(item["fastq_2"] for item in group if item["fastq_2"]),
            "paper_method": "tophat2_htseq" if first["author"] in {"Milner", "Scott-Browne"} else "kallisto",
        })
    write_tsv(project / "rna" / "metadata" / "samples.tsv", rna_rows, list(rna_rows[0]))
    write_tsv(project / "rna" / "metadata" / "samples_fine.tsv", [{"sample": row["sample"], "condition": row["condition"]} for row in rna_rows], ["sample", "condition"])
    write_tsv(project / "rna" / "metadata" / "samples_broad.tsv", [{"sample": row["sample"], "condition": row["broad_condition"]} for row in rna_rows], ["sample", "condition"])

    if (len(atac), len(units), len(rna_by_gsm), len(FINE_COMPARISONS), len(broad_comparisons)) != (33, 32, 25, 12, 21):
        raise RuntimeError("unexpected collection or comparison count")
    print(f"Wrote {len(paired_rows)} paired ATAC rows, {len(rna_rows)} RNA rows, 12 fine and 21 broad comparisons")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("data/public/raw/lcmv_cd8_bulk"))
    parser.add_argument("--project", type=Path, default=Path("data/public/processed/lcmv_cd8_bulk_fp_rna"))
    args = parser.parse_args()
    build(args.root.resolve(), args.project.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
