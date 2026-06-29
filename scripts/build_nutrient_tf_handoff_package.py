#!/usr/bin/env python3
"""Build a compact nutrient-stress TF handoff package.

The package is designed for upload to ChatGPT/Deep Research-like tools.  It
combines motif-level fp-tools differential footprint summaries with the RNA
annotations already overlaid on the all-motif aggregate-grid source tables.
"""

from __future__ import annotations

import argparse
import gzip
import math
import re
import shutil
import zipfile
from dataclasses import dataclass
from datetime import date
from pathlib import Path

import pandas as pd


LINEAGE = {
    "hpafii": "more epithelial-like pancreatic cancer cell line",
    "aspc1": "intermediate pancreatic cancer cell line",
    "panc1": "more mesenchymal-like pancreatic cancer cell line",
}

DISPLAY_CELL = {
    "hpafii": "HPAFII",
    "aspc1": "AsPC1",
    "panc1": "Panc1",
}

STRESS_ORDER = {
    "FBS": 0,
    "Glc": 1,
    "Met.Cys": 2,
    "Gln.Arg": 3,
    "Gln": 3,
    "Arg": 3,
    "BCAA": 4,
    "Trp": 5,
    "Lys": 6,
}


@dataclass(frozen=True)
class SourceSet:
    cell_key: str
    source_tsv: Path
    source_pdf: Path


def parse_args() -> argparse.Namespace:
    today = date.today().strftime("%Y%m%d")
    return argparse.ArgumentParser(
        description="Build nutrient-stress TF handoff tables, prompts, and zips."
    ).parse_args()


def find_sources(root: Path) -> list[SourceSet]:
    found = []
    for cell_key in ("hpafii", "aspc1", "panc1"):
        project = root / f"nutrient_{cell_key}_ctrl_vs_10fbs" / "reports"
        source_tsv = project / f"{cell_key}_all_motif_aggregate_grid_v6_rna_source.tsv"
        source_pdf = project / f"{cell_key}_all_motif_aggregate_grid_v6_rna.pdf"
        if not source_tsv.exists():
            raise FileNotFoundError(source_tsv)
        if not source_pdf.exists():
            raise FileNotFoundError(source_pdf)
        found.append(SourceSet(cell_key=cell_key, source_tsv=source_tsv, source_pdf=source_pdf))
    return found


def parse_rna_pairs(value: object) -> list[tuple[str, float]]:
    if not isinstance(value, str) or not value.strip() or value.strip().lower() == "nan":
        return []
    pairs = []
    for item in value.split(";"):
        if ":" not in item:
            continue
        gene, raw_val = item.split(":", 1)
        gene = gene.strip()
        if not gene:
            continue
        try:
            val = float(raw_val)
        except ValueError:
            continue
        if math.isfinite(val):
            pairs.append((gene, val))
    return pairs


def parse_condition(condition: str) -> dict[str, object]:
    stress = "Other"
    for key in ("Met.Cys", "Gln.Arg", "BCAA", "FBS", "Glc", "Gln", "Arg", "Trp", "Lys"):
        if key in condition:
            stress = key
            break
    dose_match = re.match(r"^([0-9]+(?:\.[0-9]+)?)", condition)
    dose = float(dose_match.group(1)) if dose_match else float("nan")
    return {
        "stress_type": stress,
        "dose": dose,
        "stress_order": STRESS_ORDER.get(stress, 99),
        "condition_sort": STRESS_ORDER.get(stress, 99) * 100000 - (dose if math.isfinite(dose) else -1),
    }


def evidence_class(delta_fp: float, rna_log2fc: float) -> str:
    fp_state = "footprint_gain" if delta_fp > 0 else "footprint_loss" if delta_fp < 0 else "footprint_neutral"
    if rna_log2fc >= 0.5:
        rna_state = "RNA_up"
    elif rna_log2fc <= -0.5:
        rna_state = "RNA_down"
    else:
        rna_state = "RNA_neutral"
    return f"{fp_state}_{rna_state}"


def load_motif_tables(sources: list[SourceSet]) -> pd.DataFrame:
    frames = []
    for source in sources:
        df = pd.read_csv(source.source_tsv, sep="\t")
        df.insert(0, "cell_line", DISPLAY_CELL[source.cell_key])
        df.insert(1, "cell_line_key", source.cell_key)
        df.insert(2, "cell_line_context", LINEAGE[source.cell_key])
        meta = pd.DataFrame([parse_condition(c) for c in df["condition"].astype(str)])
        df = pd.concat([df, meta], axis=1)
        frames.append(df)
    out = pd.concat(frames, ignore_index=True)
    out["abs_delta_fp"] = out["delta_fp"].abs()
    out["neg_log10_fdr"] = -out["fdr"].clip(lower=1e-300).map(math.log10)
    return out


def build_tf_long(motif_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for row in motif_df.itertuples(index=False):
        for gene, rna_fc in parse_rna_pairs(row.rna_log2fc):
            rows.append(
                {
                    "cell_line": row.cell_line,
                    "cell_line_key": row.cell_line_key,
                    "cell_line_context": row.cell_line_context,
                    "condition": row.condition,
                    "comparison": row.comparison,
                    "stress_type": row.stress_type,
                    "dose": row.dose,
                    "tf_gene": gene,
                    "motif_prefix": row.motif_prefix,
                    "motif_name": row.motif_name,
                    "motif_id": row.motif_id,
                    "delta_fp": row.delta_fp,
                    "abs_delta_fp": abs(row.delta_fp),
                    "pvalue": row.pvalue,
                    "fdr": row.fdr,
                    "neg_log10_fdr": row.neg_log10_fdr,
                    "n_sites": row.n_sites,
                    "rna_log2fc": rna_fc,
                    "abs_rna_log2fc": abs(rna_fc),
                    "evidence_class": evidence_class(row.delta_fp, rna_fc),
                    "aggregate_profile": row.aggregate_profile,
                    "profile_source": row.profile_source,
                }
            )
    return pd.DataFrame(rows)


def summarize_tf_conditions(tf_long: pd.DataFrame) -> pd.DataFrame:
    grouped = []
    keys = ["cell_line", "cell_line_key", "cell_line_context", "condition", "stress_type", "dose", "tf_gene"]
    for key, g in tf_long.groupby(keys, dropna=False):
        best_fdr_idx = g["fdr"].idxmin()
        best_abs_idx = g["abs_delta_fp"].idxmax()
        best_fdr = g.loc[best_fdr_idx]
        best_abs = g.loc[best_abs_idx]
        grouped.append(
            {
                **dict(zip(keys, key)),
                "n_motif_rows": len(g),
                "n_motif_models": g["motif_prefix"].nunique(),
                "best_motif_by_fdr": best_fdr["motif_prefix"],
                "best_motif_name": best_fdr["motif_name"],
                "best_motif_id": best_fdr["motif_id"],
                "best_fdr": best_fdr["fdr"],
                "best_pvalue": best_fdr["pvalue"],
                "max_abs_delta_fp": best_abs["abs_delta_fp"],
                "delta_fp_at_max_abs": best_abs["delta_fp"],
                "mean_delta_fp": g["delta_fp"].mean(),
                "rna_log2fc_mean": g["rna_log2fc"].mean(),
                "rna_log2fc_max_abs": g.loc[g["abs_rna_log2fc"].idxmax(), "rna_log2fc"],
                "dominant_evidence_class": g["evidence_class"].value_counts().index[0],
                "profile_sources": ";".join(sorted(map(str, g["profile_source"].dropna().unique()))),
            }
        )
    summary = pd.DataFrame(grouped)
    summary["neg_log10_best_fdr"] = -summary["best_fdr"].clip(lower=1e-300).map(math.log10)
    summary["priority_score"] = (
        summary["neg_log10_best_fdr"].clip(upper=300)
        + summary["max_abs_delta_fp"].rank(pct=True)
        + summary["rna_log2fc_max_abs"].abs().rank(pct=True)
    )
    return summary.sort_values(
        ["cell_line", "condition", "priority_score", "max_abs_delta_fp"],
        ascending=[True, True, False, False],
    )


def top_candidates(summary: pd.DataFrame, n: int = 50) -> pd.DataFrame:
    return (
        summary.sort_values(["cell_line", "condition", "priority_score"], ascending=[True, True, False])
        .groupby(["cell_line", "condition"], group_keys=False)
        .head(n)
        .reset_index(drop=True)
    )


def cross_cellline_recurrence(summary: pd.DataFrame) -> pd.DataFrame:
    significant = summary[summary["best_fdr"] <= 0.05].copy()
    rows = []
    for key, g in significant.groupby(["stress_type", "tf_gene"], dropna=False):
        cell_lines = sorted(g["cell_line"].unique())
        if len(cell_lines) < 2:
            continue
        rows.append(
            {
                "stress_type": key[0],
                "tf_gene": key[1],
                "n_cell_lines": len(cell_lines),
                "cell_lines": ";".join(cell_lines),
                "conditions": ";".join(sorted(g["condition"].unique())),
                "best_fdr": g["best_fdr"].min(),
                "max_abs_delta_fp": g["max_abs_delta_fp"].max(),
                "max_abs_rna_log2fc": g["rna_log2fc_max_abs"].abs().max(),
                "evidence_classes": ";".join(sorted(g["dominant_evidence_class"].unique())),
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["n_cell_lines", "best_fdr", "max_abs_delta_fp"], ascending=[False, True, False]
    )


def write_gzip_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        df.to_csv(handle, sep="\t", index=False)


def write_plain_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def write_docs(outdir: Path) -> None:
    meta = outdir / "metadata"
    prompts = outdir / "prompts"
    meta.mkdir(parents=True, exist_ok=True)
    prompts.mkdir(parents=True, exist_ok=True)
    (meta / "data_dictionary.md").write_text(
        """# Data Dictionary

This handoff summarizes nutrient-stress footprint and RNA evidence from fp-tools.

## Main tables

- `motif_rna_long.tsv.gz`: motif-level differential footprint rows from the RNA-annotated aggregate-grid source tables.
- `tf_condition_long.tsv.gz`: one row per TF gene, motif model, cell line, and nutrient-stress condition.
- `tf_condition_summary.tsv`: one row per TF gene, cell line, and condition after aggregation across motif models.
- `top_candidates_by_cellline_condition.tsv`: top ranked TFs per cell line and condition.
- `cross_cellline_recurrence.tsv`: TFs recurring across at least two cell lines within a stress type.

## Key columns

- `delta_fp`: fp-tools differential footprint value for nutrient stress versus `10_FBS_Ctrl`.
- `fdr` / `best_fdr`: Benjamini-Hochberg adjusted motif-level FDR from fp-tools.
- `rna_log2fc`: RNA expression log2 fold-change from the RUVr k=20 corrected DESeq2 log2-normalized RNA matrix.
- `evidence_class`: footprint gain/loss combined with RNA up/down/neutral.
- `cell_line_context`: lineage framing used for biological interpretation.

## Interpretation notes

Motif-associated footprint evidence supports regulatory motif activity but does not prove direct TF-specific occupancy. Some motif models map to TF families or dimers. Motif scanning used JASPAR 2026; motif-to-gene labels were assigned using the available JASPAR 2024 mapping and motif names.
""",
        encoding="utf-8",
    )
    (meta / "methods_summary.md").write_text(
        """# Methods Summary For Interpretation

The comparisons are nutrient stress versus `10_FBS_Ctrl` within each pancreatic cancer cell line. HPAFII is treated as more epithelial-like, AsPC1 as intermediate, and Panc1 as more mesenchymal-like.

Footprint evidence comes from fp-tools differential footprint analysis and aggregate-grid source tables. RNA evidence comes from the RUVr k=20 corrected DESeq2 log2-normalized RNA matrix. Raw counts were used upstream only to avoid reporting TFs without expression evidence.

Candidate ranking is intentionally broad: TFs are prioritized by footprint FDR, absolute footprint change, and absolute RNA log2FC. Concordant and discordant RNA/footprint patterns are both retained because TF expression and motif-associated footprint changes can be biologically decoupled.
""",
        encoding="utf-8",
    )
    (prompts / "master_prompt.md").write_text(
        """# Prompt: Integrated Nutrient-Stress TF Discovery

You are analyzing fp-tools footprint and RNA-seq summary tables from three pancreatic cancer cell lines: HPAFII, AsPC1, and Panc1. HPAFII is more epithelial-like, AsPC1 is intermediate, and Panc1 is more mesenchymal-like. All comparisons are nutrient stress versus `10_FBS_Ctrl`.

Use the attached tables to identify transcription factors with both epigenetic evidence from motif-associated footprint changes and transcriptional evidence from RNA log2FC. Analyze results separately by cell line, nutrient stress type, and dose, then identify recurrent and lineage-specific patterns across the three cell lines.

Important interpretation rules:

1. Treat `delta_fp` as motif-associated differential footprint evidence, not direct TF occupancy.
2. Treat `rna_log2fc` as TF gene expression evidence from RUVr k=20 corrected DESeq2 log2-normalized RNA.
3. Do not discard discordant footprint/RNA patterns; explain plausible regulatory interpretations.
4. Be cautious with motif-family and dimer motifs because motif scanning used JASPAR 2026 while motif-to-gene labels use the available JASPAR 2024 mapping.
5. Prioritize TFs that are biologically meaningful in pancreatic cancer, nutrient stress, metabolism, epithelial/mesenchymal lineage, stress response, apoptosis, ferroptosis, autophagy, or therapy resistance.

Please produce:

- A ranked table of candidate TFs by cell line and nutrient stress condition.
- A cross-cell-line recurrence summary.
- A lineage-specific interpretation comparing HPAFII, AsPC1, and Panc1.
- A short list of the most promising TF targets for follow-up experiments.
- Literature-supported rationale with citations.
- Caveats and suggested validation experiments.
""",
        encoding="utf-8",
    )
    (prompts / "target_prioritization_prompt.md").write_text(
        """# Prompt: TF Target Prioritization

Using the attached `top_candidates_by_cellline_condition.tsv` and `cross_cellline_recurrence.tsv`, prioritize TFs that could be candidate regulatory targets after nutrient stress in pancreatic cancer cells.

Rank candidates by:

1. strength and significance of footprint evidence;
2. RNA expression change;
3. recurrence across nutrient stresses or cell lines;
4. pancreatic cancer relevance;
5. plausibility as an intervention target;
6. lineage specificity in epithelial-like HPAFII, intermediate AsPC1, or mesenchymal-like Panc1.

For each prioritized TF, report the cell line, stress condition, footprint direction, RNA direction, biological interpretation, relevant literature, and recommended validation experiment.
""",
        encoding="utf-8",
    )
    (prompts / "literature_validation_prompt.md").write_text(
        """# Prompt: Literature Validation

Validate the highest-priority TF candidates from the attached tables using current literature. Focus on pancreatic cancer, nutrient stress, amino-acid/glucose/serum limitation, EMT or epithelial lineage state, AP-1, ATF/ISR, MYC/MAX, TEAD/YAP, HNF, NF-kB, interferon/IRF, and metabolic stress pathways when relevant.

For each TF, determine whether published evidence supports:

- altered TF expression under nutrient stress;
- chromatin accessibility or footprint changes;
- pancreatic cancer growth, survival, lineage state, invasion, therapy resistance, or metabolic adaptation;
- druggability or practical perturbation strategies.

Return a concise evidence table with citations and a final shortlist of TFs for experimental follow-up.
""",
        encoding="utf-8",
    )


def make_zip(zip_path: Path, source_dir: Path, include_suffixes: tuple[str, ...] | None = None) -> None:
    if zip_path.exists():
        zip_path.unlink()
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in sorted(source_dir.rglob("*")):
            if path.is_dir():
                continue
            if include_suffixes and not path.name.endswith(include_suffixes):
                continue
            zf.write(path, path.relative_to(source_dir))


def main() -> None:
    _ = parse_args()
    stamp = date.today().strftime("%Y%m%d")
    processed = Path("data/public/processed")
    outdir = processed / f"nutrient_tf_handoff_{stamp}"
    tables = outdir / "tables"
    figures = outdir / "figures"
    metadata = outdir / "metadata"
    outdir.mkdir(parents=True, exist_ok=True)
    tables.mkdir(exist_ok=True)
    figures.mkdir(exist_ok=True)
    metadata.mkdir(exist_ok=True)

    sources = find_sources(processed)
    motif_df = load_motif_tables(sources)
    tf_long = build_tf_long(motif_df)
    tf_summary = summarize_tf_conditions(tf_long)
    top = top_candidates(tf_summary)
    recurring = cross_cellline_recurrence(tf_summary)

    condition_meta = (
        motif_df[
            [
                "cell_line",
                "cell_line_key",
                "cell_line_context",
                "condition",
                "comparison",
                "stress_type",
                "dose",
                "stress_order",
                "condition_sort",
            ]
        ]
        .drop_duplicates()
        .sort_values(["cell_line_key", "condition_sort", "condition"])
    )

    write_gzip_tsv(motif_df, tables / "motif_rna_long.tsv.gz")
    write_gzip_tsv(tf_long, tables / "tf_condition_long.tsv.gz")
    write_plain_tsv(tf_summary, tables / "tf_condition_summary.tsv")
    write_plain_tsv(top, tables / "top_candidates_by_cellline_condition.tsv")
    write_plain_tsv(recurring, tables / "cross_cellline_recurrence.tsv")
    write_plain_tsv(condition_meta, metadata / "condition_metadata.tsv")

    manifest_rows = []
    for source in sources:
        dest = figures / source.source_pdf.name
        shutil.copy2(source.source_pdf, dest)
        manifest_rows.append(
            {
                "cell_line": DISPLAY_CELL[source.cell_key],
                "figure_pdf": f"figures/{dest.name}",
                "source_table": str(source.source_tsv),
                "local_pdf_source": str(source.source_pdf),
            }
        )
    write_plain_tsv(pd.DataFrame(manifest_rows), metadata / "figure_manifest.tsv")

    write_docs(outdir)

    (outdir / "README.md").write_text(
        f"""# fp-tools Nutrient TF Handoff

Created: {stamp}

This directory contains compact tables and prompts for identifying transcription factors with both motif-associated footprint evidence and RNA expression evidence after nutrient stress in HPAFII, AsPC1, and Panc1 pancreatic cancer cell lines.

Recommended upload to ChatGPT:

1. Upload `fp_tools_nutrient_tf_handoff_tables_prompts_{stamp}.zip`.
2. Start with `prompts/master_prompt.md`.
3. Upload `fp_tools_nutrient_tf_handoff_figures_{stamp}.zip` only if visual aggregate profiles are needed.

The compact tables are the primary source for analysis. The PDFs are optional visual references.
""",
        encoding="utf-8",
    )

    table_zip = outdir / f"fp_tools_nutrient_tf_handoff_tables_prompts_{stamp}.zip"
    figure_zip = outdir / f"fp_tools_nutrient_tf_handoff_figures_{stamp}.zip"
    with zipfile.ZipFile(table_zip, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for sub in ("README.md", "tables", "metadata", "prompts"):
            path = outdir / sub
            if path.is_file():
                zf.write(path, path.relative_to(outdir))
            else:
                for item in sorted(path.rglob("*")):
                    if item.is_file() and item.parent.name != "figures":
                        zf.write(item, item.relative_to(outdir))
    with zipfile.ZipFile(figure_zip, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for item in sorted((outdir / "figures").rglob("*")):
            if item.is_file():
                zf.write(item, item.relative_to(outdir))
        zf.write(metadata / "figure_manifest.tsv", (metadata / "figure_manifest.tsv").relative_to(outdir))

    print(f"Output directory: {outdir}")
    print(f"Table/prompt zip: {table_zip}")
    print(f"Figure zip: {figure_zip}")
    print(f"Rows: motif={len(motif_df):,}, tf_long={len(tf_long):,}, tf_summary={len(tf_summary):,}")


if __name__ == "__main__":
    main()
