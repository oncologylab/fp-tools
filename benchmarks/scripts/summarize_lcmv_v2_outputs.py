#!/usr/bin/env python3
"""Create descriptive integration, inventory, and transfer assets for LCMV v2."""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pyBigWig
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def directory_size(path: Path) -> int:
    return sum(item.stat().st_size for item in path.rglob("*") if item.is_file())


def file_count(path: Path) -> int:
    return sum(item.is_file() for item in path.rglob("*"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_peaks(path: Path) -> list[tuple[str, int, int]]:
    peaks = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            chrom, start, end = line.rstrip("\n").split("\t")[:3]
            peaks.append((chrom, int(start), int(end)))
    if len(peaks) > 10_000:
        indices = np.linspace(0, len(peaks) - 1, 10_000, dtype=int)
        peaks = [peaks[index] for index in indices]
    return peaks


def atac_descriptive(project: Path, validation: Path) -> None:
    primary = pd.read_csv(project / "metadata/samples.tsv", sep="\t", dtype=str)
    supporting = pd.read_csv(project / "metadata/supporting_samples.tsv", sep="\t", dtype=str)
    primary["collection"] = "primary"
    supporting["collection"] = "supporting"
    samples = pd.concat([primary, supporting], ignore_index=True)
    peaks = read_peaks(project / "peaks/merged_peaks_filtered.bed")
    matrix = []
    for row in samples.itertuples(index=False):
        base = project if row.collection == "primary" else project / "supporting_assay_only"
        path = base / "samples" / row.sample / "normalize" / f"{row.sample}_corrected_q95_scaled.bw"
        bw = pyBigWig.open(str(path))
        try:
            values = []
            chroms = bw.chroms()
            for chrom, start, end in peaks:
                if chrom not in chroms:
                    values.append(0.0)
                    continue
                value = bw.stats(chrom, start, min(end, chroms[chrom]), type="mean")[0]
                values.append(0.0 if value is None or not np.isfinite(value) else value)
        finally:
            bw.close()
        matrix.append(np.log1p(values))
    values = np.asarray(matrix, dtype=float)
    correlations = np.corrcoef(values)
    pd.DataFrame(correlations, index=samples["sample"], columns=samples["sample"]).to_csv(
        validation / "atac_descriptive_sample_correlations.tsv", sep="\t"
    )
    centered = values - values.mean(axis=0, keepdims=True)
    u, singular, _ = np.linalg.svd(centered, full_matrices=False)
    scores = u[:, :5] * singular[:5]
    pca = samples[["sample", "condition", "collection"]].copy()
    for index in range(scores.shape[1]):
        pca[f"PC{index + 1}"] = scores[:, index]
    pca.to_csv(validation / "atac_descriptive_pca.tsv", sep="\t", index=False)


def correlation_heatmap(table: Path, output: Path, title: str) -> None:
    frame = pd.read_csv(table, sep="\t", index_col=0)
    values = frame.apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    size = max(8.0, min(18.0, 4.0 + len(frame) * 0.22))
    fig, axis = plt.subplots(figsize=(size, size))
    image = axis.imshow(values, vmin=-1, vmax=1, cmap="coolwarm", aspect="auto")
    axis.set_title(title)
    axis.set_xticks(range(len(frame.columns)), frame.columns, rotation=90, fontsize=6)
    axis.set_yticks(range(len(frame.index)), frame.index, fontsize=6)
    fig.colorbar(image, ax=axis, label="Pearson correlation", shrink=0.75)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def linked_results(project: Path, validation: Path) -> None:
    design = pd.read_csv(project / "metadata/comparisons.tsv", sep="\t", dtype=str)
    rows = []
    for item in design.itertuples(index=False):
        fp_path = project / "comparisons" / item.comparison / "diff_footprints_results.txt"
        rna_path = project / "rna/deseq2/within_study" / f"{item.comparison}.tsv"
        fp = pd.read_csv(fp_path, sep="\t")
        rna = pd.read_csv(rna_path, sep="\t")
        change_columns = [column for column in fp if column.endswith("_change")]
        pvalue_columns = [column for column in fp if column.endswith("_pvalue")]
        change = change_columns[0] if change_columns else "change"
        pvalue = pvalue_columns[0] if pvalue_columns else "pvalue"
        motif_column = next((column for column in ("TF", "name", "motif_name") if column in fp), fp.columns[0])
        fp_sig = fp[pd.to_numeric(fp[pvalue], errors="coerce") < 0.05].copy()
        rna_sig = rna[pd.to_numeric(rna["padj"], errors="coerce") < 0.05].copy()
        left = fp[[motif_column, change]].rename(columns={motif_column: "gene_symbol", change: "fp_change"})
        joined = left.merge(rna[["gene_symbol", "log2FoldChange"]], on="gene_symbol", how="inner").dropna()
        correlation = joined["fp_change"].corr(joined["log2FoldChange"], method="spearman") if len(joined) >= 3 else np.nan
        rows.append(
            {
                "comparison": item.comparison,
                "analysis_tier": item.analysis_tier,
                "significant_motifs_p_lt_0_05": len(fp_sig),
                "significant_genes_fdr_lt_0_05": len(rna_sig),
                "matched_tf_genes": len(joined),
                "tf_gene_effect_spearman": correlation,
            }
        )
    pd.DataFrame(rows).to_csv(
        validation / "linked_atac_rna_comparisons.tsv", sep="\t", index=False
    )


def inventory(raw: Path, project: Path, validation: Path) -> None:
    locations = [
        ("raw_fastq_atac", "raw", raw / "fastq/atac", "ATAC FASTQs"),
        ("raw_fastq_rna", "raw", raw / "fastq/rna", "RNA FASTQs"),
        ("raw_atac_v1", "raw", raw / "atac_primary", "Reusable v1 ATAC BAM/peak/bigWig outputs"),
        ("raw_atac_v2_additions", "raw", raw / "atac_v2_additions", "New v2 ATAC preprocessing outputs"),
        ("v2_primary_samples", "downstream", project / "samples", "Primary paired footprint outputs"),
        ("v2_primary_comparisons", "downstream", project / "comparisons", "Primary paired differential footprints"),
        ("v2_supporting", "downstream", project / "supporting_assay_only", "Supporting ATAC-only outputs"),
        ("v2_rna_counts", "downstream", project / "rna/counts", "RNA gene count matrices"),
        ("v2_rna_results", "downstream", project / "rna/deseq2", "Within-study RNA results and descriptive integration"),
        ("v2_metadata", "downstream", project / "metadata", "Sample and comparison metadata"),
        ("v2_reports", "downstream", project / "reports", "Human-readable reports"),
        ("v2_validation", "downstream", project / "validation", "QC, audit, and transfer manifests"),
    ]
    rows = []
    for key, tier, path, description in locations:
        rows.append(
            {
                "data_family": key,
                "transfer_tier": tier,
                "absolute_path": str(path.resolve()),
                "files": file_count(path) if path.exists() else 0,
                "bytes": directory_size(path) if path.exists() else 0,
                "description": description,
            }
        )
    pd.DataFrame(rows).to_csv(validation / "data_inventory.tsv", sep="\t", index=False)


def transfer_manifests(raw: Path, project: Path, validation: Path, checksums: bool) -> None:
    repo = project.parents[3]
    downstream = sorted(
        path
        for path in project.rglob("*")
        if path.is_file()
        and "/state/" not in str(path)
        and "/logs/" not in str(path)
        and "superseded" not in str(path)
        and "invalid_k31" not in str(path)
    )
    resolved = pd.read_csv(raw / "metadata/v2/resolved_runs.tsv", sep="\t", dtype=str).fillna("")
    reproducibility = set(downstream)
    for row in resolved.itertuples(index=False):
        for value in (row.fastq_1, row.fastq_2):
            if value:
                reproducibility.add(Path(value))
    inventory_table = pd.read_csv(project / "metadata/multimodal_inventory.tsv", sep="\t", dtype=str)
    for row in inventory_table.itertuples(index=False):
        for value in (row.ATAC_BAM, f"{row.ATAC_BAM}.bai", row.ATAC_peaks):
            reproducibility.add(Path(value))
    transfer = validation / "transfer"
    transfer.mkdir(parents=True, exist_ok=True)
    for name, paths in (("downstream_files.txt", downstream), ("full_reproducibility_files.txt", sorted(reproducibility))):
        (transfer / name).write_text(
            "\n".join(str(path.resolve().relative_to(repo.resolve())) for path in paths) + "\n",
            encoding="utf-8",
        )
    if checksums:
        with (transfer / "downstream_sha256.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["sha256", "bytes", "path"])
            for path in downstream:
                writer.writerow([sha256(path), path.stat().st_size, path.resolve().relative_to(repo.resolve())])
    (transfer / "README.md").write_text(
        "# LCMV CD8 v2 transfer\n\n"
        "Run from the repository root on the destination server. Replace `<remote-host>` with the source server.\n\n"
        "```bash\n"
        "rsync -a --partial --append-verify --files-from=data/public/processed/lcmv_cd8_bulk_fp_rna_v2/validation/transfer/downstream_files.txt <remote-host>:/home/exouser/projects/fp-tools/ ./\n"
        "```\n\n"
        "Use `full_reproducibility_files.txt` instead to retrieve FASTQs and raw ATAC inputs as well.\n",
        encoding="utf-8",
    )


def report(project: Path, validation: Path) -> None:
    inventory_table = pd.read_csv(validation / "data_inventory.tsv", sep="\t")
    audit = json.loads((validation / "audit_summary.json").read_text())
    design = pd.read_csv(project / "metadata/comparison_design.tsv", sep="\t")
    condition_rows = pd.read_csv(project / "metadata/multimodal_inventory.tsv", sep="\t")
    condition_summary = condition_rows.groupby(
        ["Author", "Condition_pair_ID", "Infection", "Day", "Tissue", "Collection"], as_index=False
    ).agg(ATAC_units=("ID", "nunique"), RNA_GSM=("RNA_GSM", "first"))
    condition_summary["RNA_samples"] = condition_summary["RNA_GSM"].fillna("").map(lambda value: len([item for item in str(value).split(",") if item]))
    body = f"""<!doctype html><html><head><meta charset="utf-8"><title>LCMV CD8 v2 audit</title>
<style>body{{font:15px system-ui;margin:0;color:#172033;background:#f4f7fb}}main{{max-width:1200px;margin:auto;padding:32px}}h1{{color:#113b68}}.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:12px}}.card,section{{background:white;border:1px solid #dbe4ef;border-radius:12px;padding:16px;margin:16px 0}}.n{{font-size:28px;font-weight:800;color:#0b6e69}}table{{border-collapse:collapse;width:100%;font-size:13px}}th,td{{padding:7px;border-bottom:1px solid #e6ebf2;text-align:left}}th{{background:#edf4fa;position:sticky;top:0}}.scroll{{overflow:auto;max-height:520px}}</style></head><body><main>
<h1>LCMV CD8 ATAC + RNA v2</h1><p>Audited paired within-study analysis with assay-only supporting data and descriptive cross-study integration.</p>
<div class="cards">{''.join(f'<div class="card"><div class="n">{html.escape(str(value))}</div>{html.escape(key.replace("_", " "))}</div>' for key, value in audit.items() if isinstance(value, (int, float)) and key not in {"motifs"})}</div>
<section><h2>Condition coverage</h2><div class="scroll">{condition_summary.drop(columns="RNA_GSM").to_html(index=False, escape=True)}</div></section>
<section><h2>Comparison design</h2><div class="scroll">{design.to_html(index=False, escape=True)}</div></section>
<section><h2>Data locations</h2><div class="scroll">{inventory_table.to_html(index=False, escape=True)}</div></section>
</main></body></html>"""
    reports = project / "reports"
    reports.mkdir(parents=True, exist_ok=True)
    (reports / "lcmv_cd8_v2_overview.html").write_text(body, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw", type=Path, default=Path("data/public/raw/lcmv_cd8_bulk"))
    parser.add_argument("--project", type=Path, default=Path("data/public/processed/lcmv_cd8_bulk_fp_rna_v2"))
    parser.add_argument("--checksums", action="store_true")
    args = parser.parse_args()
    raw, project = args.raw.resolve(), args.project.resolve()
    validation = project / "validation"
    validation.mkdir(parents=True, exist_ok=True)
    atac_descriptive(project, validation)
    correlation_heatmap(
        validation / "atac_descriptive_sample_correlations.tsv",
        validation / "atac_descriptive_sample_correlations.png",
        "LCMV v2 ATAC sample correlations",
    )
    correlation_heatmap(
        project / "rna/deseq2/descriptive_sample_correlations.tsv",
        validation / "rna_descriptive_sample_correlations.png",
        "LCMV v2 RNA sample correlations",
    )
    linked_results(project, validation)
    inventory(raw, project, validation)
    transfer_manifests(raw, project, validation, args.checksums)
    report(project, validation)
    print(f"Wrote LCMV v2 reports and transfer assets under {validation}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
