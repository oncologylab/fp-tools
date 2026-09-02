#!/usr/bin/env python3
"""Render a concise frozen-parametric versus ChromBPNet comparison report."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


COMPARISON_SCHEMA = "fp-tools-chrombpnet-frozen-comparison-v1"
REPORT_SCHEMA = "fp-tools-parametric-chrombpnet-report-v1"
METHODS = {
    "DWM_conventional_geometry": ("Conventional DWM", "#6F6F6F"),
    "frozen_parametric_bias_conventional_geometry": (
        "Frozen CPU parametric",
        "#D55E00",
    ),
    "ChromBPNet_bias_conventional_geometry": ("ChromBPNet bias", "#0072B2"),
}
PARAMETRIC = "frozen_parametric_bias_conventional_geometry"
DEEP_BIAS = "ChromBPNet_bias_conventional_geometry"
DWM = "DWM_conventional_geometry"


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def checked_output(document: dict, name: str, manifest_path: Path) -> Path:
    record = document.get("outputs", {}).get(name, {})
    path = Path(record.get("path", ""))
    if not path.is_file():
        raise ValueError(f"{manifest_path} lacks declared {name} output")
    if file_sha256(path) != record.get("sha256"):
        raise ValueError(f"{manifest_path} {name} checksum mismatch")
    return path


def select_report_rows(
    metrics: pd.DataFrame,
    bootstrap: pd.DataFrame,
    safety: pd.DataFrame,
    *,
    cell: str,
    tf: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series]:
    selected = metrics[
        metrics["cell"].astype(str).eq(cell)
        & metrics["tf"].astype(str).eq(tf)
        & metrics["method"].astype(str).isin(METHODS)
    ].copy()
    if set(selected["method"].astype(str)) != set(METHODS):
        raise ValueError(f"comparison does not contain all report methods for {cell} {tf}")
    if selected["status"].astype(str).ne("eligible").any():
        raise ValueError(f"{cell} {tf} is not eligible on common finite support")
    if selected.duplicated("method").any():
        raise ValueError(f"comparison contains duplicate methods for {cell} {tf}")

    paired = bootstrap[
        bootstrap["cell"].astype(str).eq(cell)
        & bootstrap["tf"].astype(str).eq(tf)
        & bootstrap["method"].astype(str).eq(PARAMETRIC)
        & bootstrap["baseline"].astype(str).isin([DWM, DEEP_BIAS])
    ].copy()
    if set(paired["baseline"].astype(str)) != {DWM, DEEP_BIAS}:
        raise ValueError("paired bootstrap table lacks DWM or ChromBPNet comparison")
    if paired.duplicated("baseline").any():
        raise ValueError("paired bootstrap table contains duplicate comparisons")

    safe = safety[
        safety["cell"].astype(str).eq(cell)
        & safety["tf"].astype(str).eq(tf)
        & safety["residual"].astype(str).eq("deviance")
    ]
    if len(safe) != 1:
        raise ValueError(f"expected one naked-DNA deviance safety row for {cell} {tf}")
    return selected, paired, safe.iloc[0]


def build_metrics_table(
    metrics: pd.DataFrame,
    bootstrap: pd.DataFrame,
    safety: pd.Series,
) -> pd.DataFrame:
    table = metrics.copy()
    table["display_method"] = table["method"].map(
        {method: label for method, (label, _color) in METHODS.items()}
    )
    table["auroc_gain_lower_95_vs_dwm"] = np.nan
    table["auroc_gain_upper_95_vs_dwm"] = np.nan
    table["relative_auprc_gain_lower_95_vs_dwm"] = np.nan
    table["relative_auprc_gain_upper_95_vs_dwm"] = np.nan
    table["auroc_gain_lower_95_vs_chrombpnet_bias"] = np.nan
    table["auroc_gain_upper_95_vs_chrombpnet_bias"] = np.nan
    table["relative_auprc_gain_lower_95_vs_chrombpnet_bias"] = np.nan
    table["relative_auprc_gain_upper_95_vs_chrombpnet_bias"] = np.nan
    candidate_index = table.index[table["method"].astype(str).eq(PARAMETRIC)][0]
    for row in bootstrap.to_dict("records"):
        suffix = "dwm" if row["baseline"] == DWM else "chrombpnet_bias"
        for metric in ("auroc_gain", "relative_auprc_gain"):
            for bound in ("lower_95", "upper_95"):
                table.loc[candidate_index, f"{metric}_{bound}_vs_{suffix}"] = row[
                    f"{metric}_{bound}"
                ]
    for key in (
        "naked_sites",
        "finite_support",
        "false_positive_calls",
        "false_positive_rate",
        "false_positive_rate_upper_95",
    ):
        table[key] = safety[key]
    order = {method: index for index, method in enumerate(METHODS)}
    table["_order"] = table["method"].map(order)
    return table.sort_values("_order").drop(columns="_order").reset_index(drop=True)


def render_report(
    metrics: pd.DataFrame,
    bootstrap: pd.DataFrame,
    profiles: pd.DataFrame,
    safety: pd.Series,
    *,
    cell: str,
    tf: str,
    output_pdf: Path,
    output_png: Path,
) -> None:
    task_profiles = profiles[
        profiles["cell"].astype(str).eq(cell)
        & profiles["tf"].astype(str).eq(tf)
        & profiles["method"].astype(str).isin(METHODS)
    ].copy()
    if set(task_profiles["method"].astype(str)) != set(METHODS):
        raise ValueError("aggregate table lacks one or more report methods")

    by_method = metrics.set_index("method")
    dwm_pair = bootstrap[bootstrap["baseline"].astype(str).eq(DWM)].iloc[0]
    deep_pair = bootstrap[bootstrap["baseline"].astype(str).eq(DEEP_BIAS)].iloc[0]
    candidate = by_method.loc[PARAMETRIC]
    n_positive = int(candidate["n_positive"])
    n_negative = int(candidate["n_negative"])

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    figure = plt.figure(figsize=(11.69, 8.27), facecolor="white")
    figure.suptitle(
        f"{cell} {tf}: frozen CPU sequence-bias model approaches ChromBPNet bias",
        x=0.04,
        y=0.965,
        ha="left",
        fontsize=17,
        weight="bold",
    )
    figure.text(
        0.04,
        0.922,
        f"Locked chr19–22/X comparison on identical finite motif sites "
        f"({n_positive:,} ChIP-positive; {n_negative:,} matched negative)",
        ha="left",
        fontsize=10,
        color="#444444",
    )

    axis = figure.add_axes([0.07, 0.42, 0.59, 0.42])
    for method, (label, color) in METHODS.items():
        data = task_profiles[task_profiles["method"].astype(str).eq(method)].sort_values(
            "position"
        )
        x = data["position"].to_numpy(dtype=float)
        y = data["positive_minus_negative"].to_numpy(dtype=float)
        lower = data["lower_95"].to_numpy(dtype=float)
        upper = data["upper_95"].to_numpy(dtype=float)
        axis.plot(x, y, color=color, linewidth=2.0, label=label)
        axis.fill_between(x, lower, upper, color=color, alpha=0.09, linewidth=0)
    axis.axhline(0, color="#888888", linewidth=0.8)
    axis.axvspan(-10, 10, color="#CCCCCC", alpha=0.18, linewidth=0)
    axis.set(
        xlabel="Position relative to motif center (bp)",
        ylabel="ChIP-positive − matched-negative aggregate residual",
        title="Aggregate footprint separation (95% chromosome-bootstrap bands)",
        xlim=(-100, 100),
    )
    axis.legend(frameon=False, ncol=3, loc="upper right", fontsize=8)
    axis.grid(axis="y", color="#E6E6E6", linewidth=0.6)

    table_axis = figure.add_axes([0.70, 0.44, 0.27, 0.37])
    table_axis.axis("off")
    table_rows = []
    for method, (label, _color) in METHODS.items():
        row = by_method.loc[method]
        table_rows.append(
            [
                label,
                f"{row['auroc']:.3f}",
                f"{row['auprc']:.3f}",
                f"{row['functional_separation']:.3f}",
            ]
        )
    table = table_axis.table(
        cellText=table_rows,
        colLabels=["Method", "AUROC", "AUPRC", "Shape sep."],
        loc="center",
        cellLoc="center",
        colWidths=[0.43, 0.18, 0.18, 0.21],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(7.5)
    table.scale(1.0, 1.65)
    for (row_index, column), cell_box in table.get_celld().items():
        cell_box.set_edgecolor("#C7C7C7")
        cell_box.set_linewidth(0.5)
        if row_index == 0:
            cell_box.set_facecolor("#EAEAEA")
            cell_box.set_text_props(weight="bold")
        elif row_index == 2:
            cell_box.set_facecolor("#FDE9DE")
            if column == 0:
                cell_box.set_text_props(weight="bold", color="#A33F00")
    table_axis.set_title("Same geometry, same sites", fontsize=10, weight="bold", pad=8)

    gain_text = (
        "CPU parametric vs DWM\n"
        f"AUROC {candidate['auroc_gain_over_dwm']:+.3f} "
        f"[{dwm_pair['auroc_gain_lower_95']:+.3f}, {dwm_pair['auroc_gain_upper_95']:+.3f}]\n"
        f"Relative AUPRC {100 * candidate['relative_auprc_gain_over_dwm']:+.1f}% "
        f"[{100 * dwm_pair['relative_auprc_gain_lower_95']:+.1f}%, "
        f"{100 * dwm_pair['relative_auprc_gain_upper_95']:+.1f}%]"
    )
    figure.text(
        0.04,
        0.315,
        gain_text,
        ha="left",
        va="top",
        fontsize=10,
        weight="bold",
        bbox={"boxstyle": "round,pad=0.55", "facecolor": "#FDE9DE", "edgecolor": "#D55E00"},
    )
    figure.text(
        0.365,
        0.315,
        "Head-to-head with ChromBPNet bias\n"
        f"AUROC difference CI [{deep_pair['auroc_gain_lower_95']:+.3f}, "
        f"{deep_pair['auroc_gain_upper_95']:+.3f}]\n"
        f"Relative AUPRC difference CI "
        f"[{100 * deep_pair['relative_auprc_gain_lower_95']:+.1f}%, "
        f"{100 * deep_pair['relative_auprc_gain_upper_95']:+.1f}%]",
        ha="left",
        va="top",
        fontsize=9.5,
        bbox={"boxstyle": "round,pad=0.55", "facecolor": "#E5F1F8", "edgecolor": "#0072B2"},
    )
    figure.text(
        0.69,
        0.315,
        "Independent enzyme safety\n"
        f"Naked-DNA replicate 2: {int(safety['false_positive_calls'])}/"
        f"{int(safety['finite_support'])} false calls\n"
        f"FPR {100 * safety['false_positive_rate']:.1f}%; Wilson upper "
        f"{100 * safety['false_positive_rate_upper_95']:.2f}%",
        ha="left",
        va="top",
        fontsize=9.5,
        bbox={"boxstyle": "round,pad=0.55", "facecolor": "#E8F3E8", "edgecolor": "#3A7D44"},
    )
    figure.text(
        0.04,
        0.125,
        "Conclusion",
        fontsize=10,
        weight="bold",
    )
    figure.text(
        0.04,
        0.086,
        "A 21-bp reverse-complement-tied log-linear sequence model recovers most of the deep bias model's CTCF gain "
        "while remaining CPU-only and kilobyte-sized.",
        fontsize=9.5,
    )
    figure.text(
        0.04,
        0.045,
        "Limitation: internal K562 CTCF result at 10M observed fragments. MYC regresses, HepG2 deep-reference validation is pending, "
        "new promotion holdouts remain unopened, and fp-tools defaults are unchanged.",
        fontsize=8.5,
        color="#555555",
    )
    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_pdf, format="pdf", bbox_inches="tight")
    figure.savefig(output_png, format="png", dpi=170, bbox_inches="tight")
    plt.close(figure)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison-manifest", type=Path, required=True)
    parser.add_argument("--safety", type=Path, required=True)
    parser.add_argument("--cell", default="K562")
    parser.add_argument("--tf", default="CTCF")
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)

    document = json.loads(args.comparison_manifest.read_text(encoding="utf-8"))
    if document.get("schema") != COMPARISON_SCHEMA:
        raise ValueError("unsupported frozen ChromBPNet comparison manifest")
    metrics_path = checked_output(document, "metrics", args.comparison_manifest)
    bootstrap_path = checked_output(document, "bootstrap", args.comparison_manifest)
    profiles_path = checked_output(document, "profiles", args.comparison_manifest)
    metrics = pd.read_csv(metrics_path, sep="\t")
    bootstrap = pd.read_csv(bootstrap_path, sep="\t")
    profiles = pd.read_csv(profiles_path, sep="\t")
    safety = pd.read_csv(args.safety, sep="\t")
    selected, paired, safe = select_report_rows(
        metrics, bootstrap, safety, cell=args.cell, tf=args.tf
    )
    summary = build_metrics_table(selected, paired, safe)

    args.outdir.mkdir(parents=True, exist_ok=True)
    stem = f"{args.cell}_{args.tf}_parametric_vs_DWM_and_ChromBPNet"
    pdf_path = args.outdir / f"{stem}.pdf"
    png_path = args.outdir / f"{stem}.png"
    metrics_output = args.outdir / "metrics.tsv"
    readme_path = args.outdir / "README.md"
    render_report(
        selected,
        paired,
        profiles,
        safe,
        cell=args.cell,
        tf=args.tf,
        output_pdf=pdf_path,
        output_png=png_path,
    )
    summary.to_csv(metrics_output, sep="\t", index=False)
    readme_path.write_text(
        f"# {args.cell} {args.tf}: frozen parametric bias comparison\n\n"
        "The one-page report compares conventional DWM, a frozen CPU-only "
        "21-bp log-linear sequence-bias model, and the official pinned "
        "ChromBPNet bias network on identical chr19–22/X motif sites and the "
        "same conventional footprint geometry. The parametric gains over DWM "
        "are supported by chromosome-block bootstrap intervals and independent "
        "naked-DNA replicate-2 safety.\n\n"
        "Limitations: this is an internal 10M-fragment CTCF result, not a general "
        "TF correction. MYC regresses, HepG2 deep-reference validation is pending, "
        "and the new promotion holdouts remain unopened. It does not justify a "
        "package-default change.\n",
        encoding="utf-8",
    )
    outputs = {
        "pdf": pdf_path,
        "preview": png_path,
        "metrics": metrics_output,
        "readme": readme_path,
    }
    manifest = {
        "schema": REPORT_SCHEMA,
        "comparison_manifest": {
            "path": str(args.comparison_manifest),
            "sha256": file_sha256(args.comparison_manifest),
        },
        "safety": {"path": str(args.safety), "sha256": file_sha256(args.safety)},
        "cell": args.cell,
        "tf": args.tf,
        "outputs": {
            name: {"path": str(path), "sha256": file_sha256(path)}
            for name, path in outputs.items()
        },
    }
    canonical = json.dumps(manifest, sort_keys=True, separators=(",", ":"))
    manifest["report_id"] = sha256(canonical.encode()).hexdigest()
    manifest_path = args.outdir / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(summary[["display_method", "auroc", "auprc"]].to_string(index=False))
    print(pdf_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
