#!/usr/bin/env python3
"""Validate the curated LCMV ATAC/RNA collection and completed outputs."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd
import pyBigWig
import pysam


EXPECTED = {
    "libraries": 58,
    "runs": 60,
    "fastqs": 105,
    "atac_runs": 33,
    "atac_units": 32,
    "rna_samples": 25,
    "motifs": 1019,
    "fine_atac": 12,
    "broad_atac": 21,
    "fine_rna": 9,
    "broad_rna": 21,
}


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def split_ids(series: pd.Series) -> set[str]:
    return {value for cell in series.fillna("") for value in str(cell).split(",") if value}


def assert_bigwig(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise AssertionError(f"Missing bigWig: {path}")
    handle = pyBigWig.open(str(path))
    try:
        if "chr1" not in handle.chroms():
            raise AssertionError(f"bigWig does not use the expected mm10 chromosomes: {path}")
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


def validate(raw: Path, project: Path, selection: Path, verify_checksums: bool) -> dict[str, object]:
    selected = pd.read_csv(selection, sep="\t", dtype=str).fillna("")
    resolved = pd.read_csv(raw / "metadata/resolved_runs.tsv", sep="\t", dtype=str).fillna("")
    downloads = pd.read_csv(raw / "metadata/download_manifest.tsv", sep="\t", dtype=str).fillna("")
    paired = pd.read_csv(raw / "LCMV_CD8_mm10_ATAC_RNA.txt", sep="\t", dtype=str).fillna("")
    condition_pairs = pd.read_csv(raw / "metadata/condition_assay_pairs.tsv", sep="\t", dtype=str).fillna("")
    rna_meta = pd.read_csv(project / "rna/metadata/samples.tsv", sep="\t", dtype=str).fillna("")

    if (len(selected), len(resolved), len(downloads)) != (EXPECTED["libraries"], EXPECTED["runs"], EXPECTED["fastqs"]):
        raise AssertionError("Unexpected curated library/run/FASTQ count")
    if resolved.run_accession.duplicated().any() or selected.gsm_accession.duplicated().any():
        raise AssertionError("Duplicate curated GSM or resolved run")
    if resolved.query("assay == 'ATAC'").run_accession.nunique() != EXPECTED["atac_runs"]:
        raise AssertionError("Unexpected ATAC run count")
    if len(paired) != EXPECTED["atac_units"] or len(rna_meta) != EXPECTED["rna_samples"]:
        raise AssertionError("Unexpected downstream ATAC or RNA sample count")
    if split_ids(paired.ATAC_GSM) != set(resolved.query("assay == 'ATAC'").gsm_accession):
        raise AssertionError("Paired table does not cover every ATAC GSM")
    if split_ids(paired.RNA_GSM) != set(resolved.query("assay == 'RNA'").gsm_accession):
        raise AssertionError("Paired table does not cover every matched RNA GSM")
    primary = condition_pairs.query("Include_primary_paired_analysis == 'true'")
    if len(primary) != 9 or (primary.Pairing_level != "condition").any() or (primary.Same_biological_specimen != "false").any():
        raise AssertionError("Primary condition pairing semantics are incorrect")

    for row in downloads.itertuples(index=False):
        path = Path(row.local_path)
        if not path.is_file() or path.stat().st_size != int(row.expected_bytes):
            raise AssertionError(f"FASTQ size mismatch: {path}")
        if verify_checksums and md5(path) != row.checksum:
            raise AssertionError(f"FASTQ checksum mismatch: {path}")

    atac_outputs = pd.read_csv(raw / "metadata/all_atac_samples.tsv", sep="\t", dtype=str).fillna("")
    if len(atac_outputs) != EXPECTED["atac_runs"]:
        raise AssertionError("Unexpected ATAC output manifest count")
    for row in atac_outputs.itertuples(index=False):
        bam = Path(row.bam)
        if not bam.is_file() or not Path(str(bam) + ".bai").is_file():
            raise AssertionError(f"Missing BAM/BAI: {bam}")
        pysam.quickcheck(str(bam))
        assert_bed(Path(row.peaks))
        assert_bigwig(Path(row.bigwig))

    samples = pd.read_csv(project / "metadata/samples.tsv", sep="\t", dtype=str)
    assert_bed(project / "peaks/merged_peaks_filtered.bed", require_sorted=True)
    for sample in samples["sample"]:
        root = project / "samples" / sample
        assert_bigwig(root / "atac_correct" / f"{sample}_corrected.bw")
        assert_bigwig(root / "normalize" / f"{sample}_corrected_q95_scaled.bw")
        assert_bigwig(root / "footprints" / f"{sample}_footprints.bw")
        motif_root = root / "match_motifs"
        completed = list(motif_root.glob("*/beds/.done"))
        all_beds = list(motif_root.glob("*/beds/*_all.bed"))
        if len(completed) != EXPECTED["motifs"] or len(all_beds) != EXPECTED["motifs"]:
            raise AssertionError(f"Incomplete per-motif outputs for {sample}")

    for root, expected in ((project / "comparisons", EXPECTED["fine_atac"]), (project / "cross_study_atlas/comparisons", EXPECTED["broad_atac"])):
        results = list(root.glob("*/diff_footprints_results.txt"))
        reports = list(root.glob("*/diff_footprints_replicate_report.tsv"))
        if len(results) != expected or len(reports) != expected or any(path.stat().st_size == 0 for path in results + reports):
            raise AssertionError(f"Incomplete ATAC comparisons under {root}")

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
        assigned = 0
        with (root / "htseq_counts.tsv").open(encoding="utf-8") as handle:
            for line in handle:
                gene, count = line.rstrip("\n").split("\t")[:2]
                if not gene.startswith("__"):
                    assigned += int(count)
        if assigned == 0:
            raise AssertionError(f"HTSeq assigned no reads: {sample}")

    matrix = pd.read_csv(project / "rna/counts/uniform_kallisto/gene_counts_tximport_length_scaled.tsv", sep="\t")
    if matrix.shape[1] != EXPECTED["rna_samples"] + 2 or matrix.isna().any().any():
        raise AssertionError("Invalid tximport count matrix")
    for root, expected in ((project / "rna/deseq2/within_study", EXPECTED["fine_rna"]), (project / "rna/deseq2/exploratory_broad_state", EXPECTED["broad_rna"])):
        files = list(root.glob("*.tsv"))
        if len(files) != expected:
            raise AssertionError(f"Unexpected RNA comparison count under {root}")
        for path in files:
            table = pd.read_csv(path, sep="\t")
            required = {"gene_id", "gene_symbol", "log2FoldChange", "pvalue", "padj"}
            if len(table) < 1000 or not required.issubset(table.columns):
                raise AssertionError(f"Invalid DESeq2 result: {path}")
    if len(list((project / "rna/deseq2/beltra_tmm_voom_limma").glob("*.tsv"))) != 6:
        raise AssertionError("Incomplete Beltra TMM/voom/limma sensitivity results")
    if len(list((project / "rna/deseq2/sensitivity_exclude_GSM3045265").glob("*.tsv"))) != 3:
        raise AssertionError("Incomplete Guan RNA sensitivity results")
    if len(list((project / "sensitivity/exclude_GSM3045301/comparisons").glob("*/diff_footprints_results.txt"))) != 3:
        raise AssertionError("Incomplete Guan ATAC sensitivity results")

    qc = pd.concat([
        pd.read_csv(raw / "atac_primary/reports/qc_summary.tsv", sep="\t"),
        pd.read_csv(raw / "atac_supplemental/reports/qc_summary.tsv", sep="\t"),
    ], ignore_index=True)
    flags = []
    for row in qc.itertuples(index=False):
        reasons = []
        if row.frip < 0.05: reasons.append("FRiP<0.05")
        if row.tss_enrichment < 5: reasons.append("TSS_enrichment<5")
        if row.mitochondrial_fraction > 0.2: reasons.append("mitochondrial_fraction>0.2")
        if reasons: flags.append({"assay": "ATAC", "sample": row.sample, "flags": ";".join(reasons)})
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
    pd.DataFrame(flags, columns=["assay", "sample", "flags"]).to_csv(validation / "qc_flags.tsv", sep="\t", index=False)
    prior_checksum_status = False
    prior_summary = validation / "audit_summary.json"
    if prior_summary.is_file():
        prior_checksum_status = bool(json.loads(prior_summary.read_text()).get("checksums_verified"))
    summary = {
        **EXPECTED,
        "checksums_verified": verify_checksums or prior_checksum_status,
        "qc_flags": len(flags),
        "status": "passed",
    }
    (validation / "audit_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw", type=Path, default=Path("data/public/raw/lcmv_cd8_bulk"))
    parser.add_argument("--project", type=Path, default=Path("data/public/processed/lcmv_cd8_bulk_fp_rna"))
    parser.add_argument("--selection", type=Path, default=Path("benchmarks/manifests/compact/lcmv_cd8_libraries.tsv"))
    parser.add_argument("--verify-checksums", action="store_true")
    args = parser.parse_args()
    summary = validate(args.raw.resolve(), args.project.resolve(), args.selection.resolve(), args.verify_checksums)
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
