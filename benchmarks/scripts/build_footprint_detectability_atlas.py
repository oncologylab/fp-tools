#!/usr/bin/env python3
"""Build an expression-aware atlas of consistently weak motif footprints.

The atlas is intentionally descriptive.  A low motif-level footprint score is
reported as a weak-shape hypothesis, not as evidence that a transcription
factor is absent or that footprint calling failed.  Orthogonal binding labels
are required for those conclusions.

Repeated pairwise reports are first ranked within their own analysis and then
collapsed to independent biological contexts.  This prevents one ENCODE cell
line, or a repeated nutrient control, from being counted once for every
comparison in which it appears.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


RESULT_FILENAME = "diff_footprints_results.txt"
MOTIF_COLUMNS = ("motif_id", "name", "output_prefix")
OPTIONAL_MOTIF_COLUMNS = ("cluster",)


@dataclass(frozen=True)
class AtlasThresholds:
    low_percentile: float = 0.10
    minimum_low_context_fraction: float = 0.70
    minimum_expression: float = 4.0
    minimum_encode_contexts: int = 5
    minimum_nutrient_contexts: int = 10


@dataclass(frozen=True)
class AtlasArtifacts:
    context_scores: Path
    motif_summary: Path
    candidates: Path
    input_manifest: Path
    report: Path
    metadata: Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def portable_path(path: Path, path_root: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(path_root.resolve()).as_posix()
    except ValueError:
        return f"external/{resolved.name}"


def parse_project_spec(spec: str) -> tuple[str, Path]:
    if "=" not in spec:
        raise ValueError(f"Nutrient project must be CELL=DIR, received: {spec}")
    cell, raw_path = spec.split("=", 1)
    cell = cell.strip()
    path = Path(raw_path).expanduser()
    if not cell or not raw_path.strip():
        raise ValueError(f"Nutrient project must be CELL=DIR, received: {spec}")
    return cell, path


def _score_columns(frame: pd.DataFrame) -> list[str]:
    columns = [
        column
        for column in frame.columns
        if column.endswith("_mean_score") and "_rep" not in column
    ]
    if not columns:
        raise ValueError("Result table has no condition-level *_mean_score columns")
    return columns


def _validate_result_table(frame: pd.DataFrame, path: Path) -> None:
    missing = [column for column in MOTIF_COLUMNS if column not in frame.columns]
    if missing:
        raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")
    if frame["motif_id"].duplicated().any():
        duplicate = frame.loc[frame["motif_id"].duplicated(), "motif_id"].iloc[0]
        raise ValueError(f"{path} contains duplicate motif_id {duplicate}")


def read_result_contexts(
    path: Path,
    cohort: str,
    context_prefix: str,
    source_analysis: str,
) -> pd.DataFrame:
    """Read one differential result and return ranked motif-context rows."""

    frame = pd.read_csv(path, sep="\t")
    _validate_result_table(frame, path)
    metadata_columns = list(MOTIF_COLUMNS) + [
        column for column in OPTIONAL_MOTIF_COLUMNS if column in frame.columns
    ]
    rows: list[pd.DataFrame] = []
    for score_column in _score_columns(frame):
        condition = score_column[: -len("_mean_score")]
        scores = pd.to_numeric(frame[score_column], errors="coerce")
        ranked = frame[metadata_columns].copy()
        if "cluster" not in ranked:
            ranked["cluster"] = ""
        ranked["cohort"] = cohort
        ranked["biological_context"] = (
            condition if not context_prefix else f"{context_prefix}|{condition}"
        )
        ranked["source_analysis"] = source_analysis
        ranked["score"] = scores
        ranked["percentile"] = scores.rank(method="average", pct=True)
        rows.append(ranked)
    return pd.concat(rows, ignore_index=True)


def collect_encode_results(root: Path) -> tuple[pd.DataFrame, list[Path]]:
    paths = sorted(root.glob(f"pairs/*/results/{RESULT_FILENAME}"))
    if not paths:
        raise FileNotFoundError(f"No ENCODE pairwise result tables found under {root}")
    frames = [
        read_result_contexts(
            path,
            cohort="ENCODE",
            context_prefix="",
            source_analysis=path.parents[1].name,
        )
        for path in paths
    ]
    return pd.concat(frames, ignore_index=True), paths


def collect_nutrient_results(
    projects: Iterable[tuple[str, Path]],
) -> tuple[pd.DataFrame, list[Path]]:
    frames: list[pd.DataFrame] = []
    paths: list[Path] = []
    for cell, root in projects:
        project_paths = sorted(root.glob(f"comparisons/*/{RESULT_FILENAME}"))
        if not project_paths:
            raise FileNotFoundError(f"No nutrient result tables found under {root}")
        paths.extend(project_paths)
        for path in project_paths:
            frames.append(
                read_result_contexts(
                    path,
                    cohort="NUTRIENT",
                    context_prefix=cell,
                    source_analysis=f"{cell}|{path.parent.name}",
                )
            )
    if not frames:
        raise ValueError("At least one nutrient project is required")
    return pd.concat(frames, ignore_index=True), paths


def validate_motif_metadata(frame: pd.DataFrame) -> None:
    for motif_id, group in frame.groupby("motif_id", sort=False):
        for column in ("name", "output_prefix"):
            values = group[column].fillna("").astype(str).unique()
            if len(values) > 1:
                rendered = ", ".join(sorted(values)[:5])
                raise ValueError(
                    f"Inconsistent {column} metadata for {motif_id}: {rendered}"
                )


def join_unique(values: pd.Series) -> str:
    return ";".join(
        sorted(
            {
                item
                for value in values
                if pd.notna(value)
                for item in str(value).split(";")
                if item
            }
        )
    )


def collapse_biological_contexts(
    records: pd.DataFrame,
    low_percentile: float,
) -> pd.DataFrame:
    """Collapse repeated pair/control rows to one row per biological context."""

    validate_motif_metadata(records)
    group_columns = [
        "cohort",
        "biological_context",
        "motif_id",
        "name",
        "output_prefix",
    ]
    collapsed = (
        records.groupby(group_columns, as_index=False, dropna=False)
        .agg(
            n_source_analyses=("source_analysis", "nunique"),
            cluster=("cluster", join_unique),
            median_score=("score", "median"),
            median_percentile=("percentile", "median"),
        )
        .sort_values(["cohort", "biological_context", "median_percentile", "motif_id"])
        .reset_index(drop=True)
    )
    collapsed["bottom_percentile"] = collapsed["median_percentile"] <= low_percentile
    return collapsed


def load_expression(path: Path | None) -> pd.DataFrame:
    if path is None:
        return pd.DataFrame(columns=["expression_gene", "expression_median", "expression_samples"])
    frame = pd.read_csv(path, sep="\t")
    if "gene_key" not in frame:
        raise ValueError(f"{path} is missing required gene_key column")
    numeric = frame.drop(columns=["gene_key"]).apply(pd.to_numeric, errors="coerce")
    result = pd.DataFrame(
        {
            "expression_gene": frame["gene_key"].astype(str).str.upper(),
            "expression_median": numeric.median(axis=1, skipna=True),
            "expression_samples": numeric.notna().sum(axis=1),
        }
    )
    return (
        result.groupby("expression_gene", as_index=False)
        .agg(
            expression_median=("expression_median", "median"),
            expression_samples=("expression_samples", "max"),
        )
    )


def _cohort_summary(context_scores: pd.DataFrame, cohort: str) -> pd.DataFrame:
    subset = context_scores[context_scores["cohort"] == cohort]
    prefix = cohort.lower()
    return (
        subset.groupby(["motif_id", "name", "output_prefix"], as_index=False)
        .agg(
            **{
                f"{prefix}_clusters": ("cluster", join_unique),
                f"{prefix}_contexts": ("biological_context", "nunique"),
                f"{prefix}_median_score": ("median_score", "median"),
                f"{prefix}_median_percentile": ("median_percentile", "median"),
                f"{prefix}_bottom_fraction": ("bottom_percentile", "mean"),
            }
        )
    )


def summarize_motifs(
    context_scores: pd.DataFrame,
    expression: pd.DataFrame,
    thresholds: AtlasThresholds,
) -> pd.DataFrame:
    encode = _cohort_summary(context_scores, "ENCODE")
    nutrient = _cohort_summary(context_scores, "NUTRIENT")
    motif_columns = ["motif_id", "name", "output_prefix"]
    summary = encode.merge(nutrient, on=motif_columns, how="outer", validate="one_to_one")
    summary["cluster"] = summary[["encode_clusters", "nutrient_clusters"]].apply(
        lambda row: ";".join(
            sorted(
                {
                    item
                    for value in row
                    if pd.notna(value)
                    for item in str(value).split(";")
                    if item
                }
            )
        ),
        axis=1,
    )

    summary["expression_gene"] = summary["name"].where(
        summary["name"].astype(str).str.fullmatch(r"[A-Za-z0-9-]+"), ""
    ).str.upper()
    summary = summary.merge(expression, on="expression_gene", how="left", validate="many_to_one")

    weak_both = (
        (summary["encode_contexts"] >= thresholds.minimum_encode_contexts)
        & (summary["nutrient_contexts"] >= thresholds.minimum_nutrient_contexts)
        & (summary["encode_median_percentile"] <= thresholds.low_percentile)
        & (summary["nutrient_median_percentile"] <= thresholds.low_percentile)
        & (summary["encode_bottom_fraction"] >= thresholds.minimum_low_context_fraction)
        & (summary["nutrient_bottom_fraction"] >= thresholds.minimum_low_context_fraction)
    )
    expression_known = summary["expression_median"].notna()
    expression_pass = summary["expression_median"] >= thresholds.minimum_expression
    summary["candidate_status"] = np.select(
        [weak_both & expression_pass, weak_both & ~expression_known, weak_both],
        [
            "weak_shape_expressed",
            "weak_shape_expression_unknown",
            "weak_shape_low_expression",
        ],
        default="not_consistently_weak",
    )
    summary["weak_shape_candidate"] = summary["candidate_status"] == "weak_shape_expressed"
    summary["combined_percentile"] = np.sqrt(
        summary["encode_median_percentile"] * summary["nutrient_median_percentile"]
    )
    summary["interpretation"] = np.where(
        summary["weak_shape_candidate"],
        "Expression-supported weak aggregate shape; orthogonal occupancy is required to call a detection failure.",
        "Not selected as an expression-supported weak-shape candidate under the configured thresholds.",
    )
    return summary.sort_values(
        ["weak_shape_candidate", "combined_percentile", "motif_id"],
        ascending=[False, True, True],
        na_position="last",
    ).reset_index(drop=True)


def build_input_manifest(paths: Iterable[Path], path_root: Path) -> pd.DataFrame:
    rows = []
    for path in sorted(set(Path(item) for item in paths)):
        rows.append(
            {
                "path": portable_path(path, path_root),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    return pd.DataFrame(rows)


def _format_report_table(frame: pd.DataFrame) -> str:
    display = frame.copy()
    numeric_columns = [
        column
        for column in display.columns
        if column.endswith(("_percentile", "_fraction", "_median"))
        or column == "combined_percentile"
    ]
    for column in numeric_columns:
        if column in display:
            display[column] = display[column].map(
                lambda value: "" if pd.isna(value) else f"{float(value):.4f}"
            )
    return display.to_html(index=False, escape=True, table_id="atlas-table", border=0)


def write_html_report(
    summary: pd.DataFrame,
    context_scores: pd.DataFrame,
    thresholds: AtlasThresholds,
    output: Path,
) -> None:
    candidate_count = int(summary["weak_shape_candidate"].sum())
    context_counts = context_scores.groupby("cohort")["biological_context"].nunique().to_dict()
    columns = [
        "motif_id",
        "name",
        "cluster",
        "candidate_status",
        "expression_median",
        "encode_median_percentile",
        "encode_bottom_fraction",
        "nutrient_median_percentile",
        "nutrient_bottom_fraction",
        "combined_percentile",
    ]
    table = _format_report_table(summary[columns])
    document = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>fp-tools footprint detectability atlas</title>
<style>
body {{ color:#172033; background:#f5f7fb; font-family:Inter,Arial,sans-serif; margin:0; }}
main {{ max-width:1700px; margin:0 auto; padding:32px; }}
h1 {{ margin:0 0 8px; }}
.notice {{ background:#fff7df; border:1px solid #e5c468; border-radius:8px; padding:14px 16px; }}
.summary {{ display:flex; gap:12px; flex-wrap:wrap; margin:18px 0; }}
.card {{ background:white; border:1px solid #d8deea; border-radius:8px; padding:12px 16px; }}
input {{ width:min(560px,100%); box-sizing:border-box; padding:10px; margin:8px 0 16px; }}
.table-wrap {{ overflow:auto; background:white; border:1px solid #d8deea; border-radius:8px; }}
table {{ border-collapse:collapse; width:100%; font-size:13px; white-space:nowrap; }}
th,td {{ border-bottom:1px solid #e5e9f1; padding:8px 10px; text-align:left; }}
th {{ background:#eef2f8; position:sticky; top:0; }}
tr:hover {{ background:#f7faff; }}
code {{ background:#eef2f8; padding:2px 4px; }}
</style>
</head>
<body><main>
<h1>Footprint detectability atlas</h1>
<p class="notice"><strong>Interpretation:</strong> this report identifies weak aggregate-shape hypotheses. A low score does not establish TF absence, failed Tn5 correction, or failed footprint detection. Those conclusions require matched occupancy controls and correction ablations.</p>
<div class="summary">
<div class="card"><strong>{len(summary):,}</strong><br>motifs</div>
<div class="card"><strong>{context_counts.get('ENCODE', 0):,}</strong><br>ENCODE cell contexts</div>
<div class="card"><strong>{context_counts.get('NUTRIENT', 0):,}</strong><br>nutrient cell-condition contexts</div>
<div class="card"><strong>{candidate_count:,}</strong><br>expression-supported candidates</div>
</div>
<p>Candidates require median percentile ≤ {thresholds.low_percentile:.2f} and bottom-percentile membership in at least {thresholds.minimum_low_context_fraction:.0%} of both cohorts, with median expression ≥ {thresholds.minimum_expression:g}.</p>
<label for="filter"><strong>Filter motifs</strong></label><br>
<input id="filter" type="search" placeholder="Search motif, TF, cluster, or status">
<div class="table-wrap">{table}</div>
</main>
<script>
const input=document.getElementById('filter');
const rows=[...document.querySelectorAll('#atlas-table tbody tr')];
input.addEventListener('input',()=>{{const q=input.value.toLowerCase();rows.forEach(row=>{{row.hidden=!row.textContent.toLowerCase().includes(q);}});}});
</script></body></html>"""
    output.write_text(document, encoding="utf-8")


def build_atlas(
    encode_root: Path,
    nutrient_projects: Iterable[tuple[str, Path]],
    expression_path: Path | None,
    outdir: Path,
    thresholds: AtlasThresholds = AtlasThresholds(),
    path_root: Path = Path("."),
) -> AtlasArtifacts:
    encode, encode_paths = collect_encode_results(encode_root)
    nutrient, nutrient_paths = collect_nutrient_results(nutrient_projects)
    records = pd.concat([encode, nutrient], ignore_index=True)
    context_scores = collapse_biological_contexts(records, thresholds.low_percentile)
    expression = load_expression(expression_path)
    summary = summarize_motifs(context_scores, expression, thresholds)
    candidates = summary[summary["weak_shape_candidate"]].copy()

    outdir.mkdir(parents=True, exist_ok=True)
    artifacts = AtlasArtifacts(
        context_scores=outdir / "detectability_context_scores.tsv.gz",
        motif_summary=outdir / "detectability_motif_summary.tsv",
        candidates=outdir / "detectability_candidates.tsv",
        input_manifest=outdir / "detectability_input_manifest.tsv",
        report=outdir / "detectability_atlas.html",
        metadata=outdir / "detectability_metadata.json",
    )
    context_scores.to_csv(artifacts.context_scores, sep="\t", index=False, compression="gzip")
    summary.to_csv(artifacts.motif_summary, sep="\t", index=False)
    candidates.to_csv(artifacts.candidates, sep="\t", index=False)

    input_paths = encode_paths + nutrient_paths
    if expression_path is not None:
        input_paths.append(expression_path)
    manifest = build_input_manifest(input_paths, path_root)
    manifest.to_csv(artifacts.input_manifest, sep="\t", index=False)
    write_html_report(summary, context_scores, thresholds, artifacts.report)

    metadata = {
        "schema_version": 1,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "interpretation": "Weak aggregate-shape hypotheses; not occupancy or failure calls.",
        "motifs": int(len(summary)),
        "weak_shape_candidates": int(len(candidates)),
        "contexts": {
            key: int(value)
            for key, value in context_scores.groupby("cohort")["biological_context"].nunique().items()
        },
        "thresholds": thresholds.__dict__,
        "artifacts": {
            key: value.name for key, value in artifacts.__dict__.items() if key != "metadata"
        },
    }
    artifacts.metadata.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return artifacts


def default_nutrient_projects(root: Path) -> list[tuple[str, Path]]:
    return [
        ("ASPC1", root / "data/public/processed/nutrient_aspc1_ctrl_vs_10fbs"),
        ("HPAFII", root / "data/public/processed/nutrient_hpafii_ctrl_vs_10fbs"),
        ("PANC1", root / "data/public/processed/nutrient_panc1_ctrl_vs_10fbs"),
    ]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--encode-root",
        type=Path,
        default=Path("data/public/processed/encode_cancer_pairwise_q95_20260814"),
    )
    parser.add_argument(
        "--nutrient-project",
        action="append",
        default=[],
        metavar="CELL=DIR",
        help="Repeat for each nutrient project. Defaults to the three local projects.",
    )
    parser.add_argument(
        "--expression",
        type=Path,
        default=Path("data/public/raw/nutrient_rna/nutrient_rna_deseq2_log2norm_ruvr_k20.tsv.gz"),
    )
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--path-root", type=Path, default=Path("."))
    parser.add_argument("--low-percentile", type=float, default=0.10)
    parser.add_argument("--minimum-low-context-fraction", type=float, default=0.70)
    parser.add_argument("--minimum-expression", type=float, default=4.0)
    parser.add_argument("--minimum-encode-contexts", type=int, default=5)
    parser.add_argument("--minimum-nutrient-contexts", type=int, default=10)
    args = parser.parse_args(argv)

    if not 0 < args.low_percentile <= 1:
        parser.error("--low-percentile must be in (0, 1]")
    if not 0 <= args.minimum_low_context_fraction <= 1:
        parser.error("--minimum-low-context-fraction must be in [0, 1]")

    projects = (
        [parse_project_spec(spec) for spec in args.nutrient_project]
        if args.nutrient_project
        else default_nutrient_projects(args.path_root)
    )
    thresholds = AtlasThresholds(
        low_percentile=args.low_percentile,
        minimum_low_context_fraction=args.minimum_low_context_fraction,
        minimum_expression=args.minimum_expression,
        minimum_encode_contexts=args.minimum_encode_contexts,
        minimum_nutrient_contexts=args.minimum_nutrient_contexts,
    )
    artifacts = build_atlas(
        args.encode_root,
        projects,
        args.expression,
        args.outdir,
        thresholds=thresholds,
        path_root=args.path_root,
    )
    summary = pd.read_csv(artifacts.motif_summary, sep="\t")
    candidates = summary[summary["weak_shape_candidate"]]
    print(
        f"wrote {len(summary):,} motifs and {len(candidates):,} expression-supported "
        f"weak-shape candidates to {args.outdir}"
    )
    if not candidates.empty:
        columns = [
            "name",
            "motif_id",
            "expression_median",
            "encode_median_percentile",
            "nutrient_median_percentile",
        ]
        print(candidates[columns].head(20).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
