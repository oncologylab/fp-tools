#!/usr/bin/env python3
"""Validate a versioned LCMV ATAC/RNA collection and completed outputs."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd
import pyBigWig
import pysam


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def split_ids(series: pd.Series) -> set[str]:
    return {
        value
        for cell in series.fillna("")
        for value in str(cell).split(",")
        if value
    }


def assert_bigwig(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise AssertionError(f"Missing bigWig: {path}")
    handle = pyBigWig.open(str(path))
    try:
        if "chr1" not in handle.chroms():
            raise AssertionError(f"Unexpected bigWig chromosomes: {path}")
    finally:
        handle.close()


def assert_bed(path: Path, require_sorted: bool = False) -> None:
    previous = None
    rows = 0
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise AssertionError(f"BED row has fewer than three fields: {path}")
            start, end = int(fields[1]), int(fields[2])
            if start < 0 or end <= start:
                raise AssertionError(f"Invalid BED interval in {path}: {line.rstrip()}")
            key = (fields[0], start, end)
            if require_sorted and previous is not None and key < previous:
                raise AssertionError(f"BED is not sorted: {path}")
            previous, rows = key, rows + 1
    if rows == 0:
        raise AssertionError(f"Empty BED: {path}")


def validate(
    raw_metadata: Path,
    project: Path,
    selection_path: Path,
    comparison_path: Path,
    atac_qc_paths: list[Path],
    verify_checksums: bool,
    expected_motifs: int,
) -> dict[str, object]:
    selected = pd.read_csv(selection_path, sep="\t", dtype=str).fillna("")
    resolved = pd.read_csv(raw_metadata / "resolved_runs.tsv", sep="\t", dtype=str).fillna("")
    downloads = pd.read_csv(raw_metadata / "download_manifest.tsv", sep="\t", dtype=str).fillna("")
    inventory = pd.read_csv(project / "metadata/multimodal_inventory.tsv", sep="\t", dtype=str).fillna("")
    rna_meta = pd.read_csv(project / "rna/metadata/samples.tsv", sep="\t", dtype=str).fillna("")
    design = pd.read_csv(comparison_path, sep="\t", dtype=str).fillna("")

    if len(selected) != selected.gsm_accession.nunique():
        raise AssertionError("Duplicate curated GSM accession")
    if set(selected.gsm_accession) != set(resolved.gsm_accession):
        raise AssertionError("Resolved runs do not cover the exact curated GSM set")
    if resolved.run_accession.duplicated().any():
        raise AssertionError("Duplicate resolved run accession")
    expected_fastqs = sum(1 + bool(row.fastq_2) for row in resolved.itertuples())
    if len(downloads) != expected_fastqs:
        raise AssertionError("Download manifest does not cover every resolved FASTQ")

    primary_selection = selected.query("collection == 'primary'")
    primary_conditions = set(primary_selection.condition_pair_id)
    for condition in primary_conditions:
        assays = set(primary_selection.query("condition_pair_id == @condition").assay)
        if assays != {"ATAC", "RNA"}:
            raise AssertionError(f"Primary condition is not paired: {condition}")
    if (primary_selection.include_in_primary_paired_analysis != "true").any():
        raise AssertionError("Primary selection contains an excluded library")
    supporting = selected.query("collection == 'supplemental'")
    if (supporting.include_in_primary_paired_analysis != "false").any():
        raise AssertionError("Supporting library entered primary paired analysis")

    if split_ids(inventory.ATAC_GSM) != set(resolved.query("assay == 'ATAC'").gsm_accession):
        raise AssertionError("ATAC inventory does not cover every curated ATAC GSM")
    if split_ids(inventory.RNA_GSM) != set(resolved.query("assay == 'RNA' and collection == 'primary'").gsm_accession):
        raise AssertionError("Multimodal inventory does not cover every primary RNA GSM")
    if (inventory.Same_biological_specimen != "false").any():
        raise AssertionError("Condition-level pairing was mislabeled as specimen pairing")

    for row in downloads.itertuples(index=False):
        path = Path(row.local_path)
        if not path.is_file() or path.stat().st_size != int(row.expected_bytes):
            raise AssertionError(f"FASTQ size mismatch: {path}")
        if verify_checksums and md5(path) != row.checksum:
            raise AssertionError(f"FASTQ checksum mismatch: {path}")

    for row in inventory.itertuples(index=False):
        bam = Path(row.ATAC_BAM)
        if not bam.is_file() or not Path(str(bam) + ".bai").is_file():
            raise AssertionError(f"Missing BAM/BAI: {bam}")
        pysam.quickcheck(str(bam))
        assert_bed(Path(row.ATAC_peaks))

    primary_samples = pd.read_csv(project / "metadata/samples.tsv", sep="\t", dtype=str)
    supporting_samples = pd.read_csv(project / "metadata/supporting_samples.tsv", sep="\t", dtype=str)
    assert_bed(project / "peaks/merged_peaks_filtered.bed", require_sorted=True)
    for collection, samples in (("primary", primary_samples), ("supporting", supporting_samples)):
        base = project if collection == "primary" else project / "supporting_assay_only"
        for sample in samples["sample"]:
            root = base / "samples" / sample
            assert_bigwig(root / "atac_correct" / f"{sample}_corrected.bw")
            assert_bigwig(root / "normalize" / f"{sample}_corrected_q95_scaled.bw")
            assert_bigwig(root / "footprints" / f"{sample}_footprints.bw")
            motif_root = root / "match_motifs"
            completed = list(motif_root.glob("*/beds/.done"))
            all_beds = list(motif_root.glob("*/beds/*_all.bed"))
            if len(completed) != expected_motifs or len(all_beds) != expected_motifs:
                raise AssertionError(f"Incomplete per-motif outputs for {sample}")

    primary_design = design[design.analysis_tier.str.startswith("primary_")]
    supporting_atac_design = design.query("analysis_tier == 'supporting_atac_only'")
    supporting_rna_design = design.query("analysis_tier == 'supporting_rna_only'")
    for root, expected in (
        (project / "comparisons", set(primary_design.comparison)),
        (project / "supporting_assay_only/comparisons", set(supporting_atac_design.comparison)),
    ):
        observed = {path.parent.name for path in root.glob("*/diff_footprints_results.txt")}
        if observed != expected:
            raise AssertionError(f"Comparison result set mismatch under {root}")
        for name in expected:
            for filename in ("diff_footprints_results.txt", "diff_footprints_replicate_report.tsv"):
                path = root / name / filename
                if not path.is_file() or path.stat().st_size == 0:
                    raise AssertionError(f"Missing comparison output: {path}")

    for sample in rna_meta["sample"]:
        info = json.loads((project / "rna/uniform_kallisto" / sample / "run_info.json").read_text())
        if info["k-mer length"] != 21 or info["n_processed"] <= 0 or info["n_pseudoaligned"] <= 0:
            raise AssertionError(f"Invalid uniform Kallisto result: {sample}")
    for sample in rna_meta.query("paper_method == 'kallisto'")["sample"]:
        info = json.loads((project / "rna/paper_specific" / sample / "run_info.json").read_text())
        if info["k-mer length"] != 31 or info["n_pseudoaligned"] <= 0:
            raise AssertionError(f"Invalid paper-specific Kallisto result: {sample}")
    for sample in rna_meta.query("paper_method == 'tophat2_htseq'")["sample"]:
        root = project / "rna/paper_specific" / sample
        if not (root / "htseq_counts.done").is_file() or not (root / "htseq_counts.tsv").is_file():
            raise AssertionError(f"Incomplete HTSeq result: {sample}")

    matrix = pd.read_csv(
        project / "rna/counts/uniform_kallisto/gene_counts_tximport_length_scaled.tsv",
        sep="\t",
    )
    if matrix.shape[1] != len(rna_meta) + 2 or matrix.isna().any().any():
        raise AssertionError("Invalid uniform tximport count matrix")
    for root, expected in (
        (project / "rna/deseq2/within_study", set(primary_design.comparison)),
        (project / "rna/deseq2/supporting_rna_only", set(supporting_rna_design.comparison)),
    ):
        observed = {path.stem for path in root.glob("*.tsv")}
        if observed != expected:
            raise AssertionError(f"RNA comparison result set mismatch under {root}")
        for path in root.glob("*.tsv"):
            table = pd.read_csv(path, sep="\t")
            required = {"gene_id", "gene_symbol", "log2FoldChange", "pvalue", "padj"}
            if len(table) < 1000 or not required.issubset(table.columns):
                raise AssertionError(f"Invalid RNA result: {path}")
    for path in (
        project / "rna/deseq2/descriptive_sample_correlations.tsv",
        project / "rna/deseq2/descriptive_pca.tsv",
    ):
        if not path.is_file() or path.stat().st_size == 0:
            raise AssertionError(f"Missing descriptive RNA output: {path}")
    if (project / "rna/deseq2/exploratory_broad_state").exists():
        raise AssertionError("v2 must not contain pooled cross-study differential tests")

    qc_tables = [pd.read_csv(path, sep="\t") for path in atac_qc_paths if path.is_file()]
    flags: list[dict[str, str]] = []
    if qc_tables:
        qc = pd.concat(qc_tables, ignore_index=True)
        for row in qc.itertuples(index=False):
            reasons = []
            if row.frip < 0.05:
                reasons.append("FRiP<0.05")
            if row.tss_enrichment < 5:
                reasons.append("TSS_enrichment<5")
            if row.mitochondrial_fraction > 0.2:
                reasons.append("mitochondrial_fraction>0.2")
            if reasons:
                flags.append({"assay": "ATAC", "sample": row.sample, "flags": ";".join(reasons)})
    for sample in rna_meta["sample"]:
        info = json.loads((project / "rna/uniform_kallisto" / sample / "run_info.json").read_text())
        reasons = []
        if info["p_pseudoaligned"] < 50:
            reasons.append("uniform_k21_pseudoalignment<50%")
        paper_info = project / "rna/paper_specific" / sample / "run_info.json"
        if paper_info.is_file() and json.loads(paper_info.read_text())["p_pseudoaligned"] < 50:
            reasons.append("paper_k31_pseudoalignment<50%")
        if reasons:
            flags.append({"assay": "RNA", "sample": sample, "flags": ";".join(reasons)})

    validation = project / "validation"
    validation.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(flags, columns=["assay", "sample", "flags"]).to_csv(
        validation / "qc_flags.tsv", sep="\t", index=False
    )
    summary = {
        "libraries": len(selected),
        "runs": len(resolved),
        "fastqs": len(downloads),
        "primary_conditions": len(primary_conditions),
        "atac_runs": int((resolved.assay == "ATAC").sum()),
        "atac_units": len(inventory),
        "primary_atac_units": len(primary_samples),
        "supporting_atac_units": len(supporting_samples),
        "rna_samples": len(rna_meta),
        "primary_rna_samples": int((rna_meta.collection == "primary").sum()),
        "motifs": expected_motifs,
        "primary_atac_comparisons": len(primary_design),
        "primary_rna_comparisons": len(primary_design),
        "supporting_atac_comparisons": len(supporting_atac_design),
        "supporting_rna_comparisons": len(supporting_rna_design),
        "checksums_verified": verify_checksums,
        "qc_flags": len(flags),
        "status": "passed",
    }
    (validation / "audit_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-metadata", type=Path, required=True)
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--selection", type=Path, required=True)
    parser.add_argument("--comparisons", type=Path, required=True)
    parser.add_argument("--atac-qc", type=Path, action="append", default=[])
    parser.add_argument("--verify-checksums", action="store_true")
    parser.add_argument("--expected-motifs", type=int, default=1019)
    args = parser.parse_args()
    summary = validate(
        args.raw_metadata.resolve(),
        args.project.resolve(),
        args.selection.resolve(),
        args.comparisons.resolve(),
        [path.resolve() for path in args.atac_qc],
        args.verify_checksums,
        args.expected_motifs,
    )
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
