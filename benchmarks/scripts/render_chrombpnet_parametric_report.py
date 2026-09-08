#!/usr/bin/env python3
"""Render a concise frozen-parametric versus ChromBPNet comparison report."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from audit_chrombpnet_comparison import check_identity, checked_path, canonical_id  # noqa: E402
from audit_direct_bias_safety import OPERATOR, SCHEMA as SAFETY_SCHEMA, read_thresholds  # noqa: E402
from audit_direct_bias_safety_coverage import SCHEMA_COVERAGE  # noqa: E402


COMPARISON_SCHEMA = "fp-tools-chrombpnet-frozen-comparison-v1"
REPORT_SCHEMA = "fp-tools-parametric-chrombpnet-report-v2"
METHODS = {
    "DWM_conventional_geometry": ("DWM expectations", "#6F6F6F"),
    "frozen_parametric_bias_conventional_geometry": (
        "Frozen CPU parametric",
        "#D55E00",
    ),
    "ChromBPNet_bias_conventional_geometry": ("ChromBPNet bias", "#0072B2"),
}
PARAMETRIC = "frozen_parametric_bias_conventional_geometry"
DEEP_BIAS = "ChromBPNet_bias_conventional_geometry"
DWM = "DWM_conventional_geometry"
RAW = "raw_geometry"
SAFETY_COLUMNS = (
    "support", "informative_support", "false_positive_calls", "false_positive_rate",
    "false_positive_rate_upper_95", "informative_false_positive_rate",
    "informative_false_positive_rate_upper_95", "paired_support", "paired_false_positive_calls",
    "paired_dwm_calls", "paired_false_positive_rate", "paired_false_positive_rate_upper_95",
    "paired_false_positive_rate_increase_over_dwm", "dwm_unavailable_informative_sites",
)


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def checked_output(document: dict, name: str, manifest_path: Path, data_root: Path = Path(".")) -> Path:
    record = document.get("outputs", {}).get(name, {})
    if not record.get("path"):
        raise ValueError(f"{manifest_path} lacks declared {name} output")
    return checked_path(record, data_root)


def validate_safety_link(coverage, source, configuration, comparison, *, cell):
    if (coverage.get("schema") != SCHEMA_COVERAGE or source.get("schema") != SAFETY_SCHEMA
            or coverage.get("operator") != OPERATOR or source.get("operator") != OPERATOR
            or coverage.get("replicate") != "NakedDNA_rep2" or source.get("replicate") != "NakedDNA_rep2"
            or coverage.get("final_untouched_safety") is not False
            or coverage.get("new_holdout_labels_read") is not False):
        raise ValueError("report requires coverage-qualified direct-operator safety")
    for document, key in ((coverage, "coverage_id"), (source, "safety_id"), (configuration, "configuration_id")):
        check_identity(document, key)
    if (cell not in {"K562", "HepG2"} or cell not in configuration.get("cells", {})
            or configuration["cells"][cell]["comparison_result_id"] != comparison["result_id"]
            or coverage["configuration_id"] != configuration["configuration_id"]
            or source["configuration_id"] != configuration["configuration_id"]
            or coverage["threshold_id"] != source["threshold_id"]):
        raise ValueError("safety operator configuration does not match the reported comparison")


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
        & bootstrap["baseline"].astype(str).isin([DWM, RAW, DEEP_BIAS])
    ].copy()
    if set(paired["baseline"].astype(str)) != {DWM, RAW, DEEP_BIAS}:
        raise ValueError("paired bootstrap table lacks DWM, raw or ChromBPNet comparison")
    if paired.duplicated("baseline").any():
        raise ValueError("paired bootstrap table contains duplicate comparisons")

    if not set(SAFETY_COLUMNS + ("passes_candidate_screen", "dispersion")).issubset(safety.columns):
        raise ValueError("report requires coverage-qualified direct-operator safety rows")
    safe = safety[
        safety["cell"].astype(str).eq(cell)
        & safety["tf"].astype(str).eq(tf)
    ]
    if len(safe) != 1:
        raise ValueError(f"expected one direct-operator safety row for {cell} {tf}")
    if safe.iloc[0].passes_candidate_screen != True:
        raise ValueError("candidate lacks a passing exact-operator safety screen")
    gains = paired[paired.baseline.isin([DWM, RAW])]
    if not (gains.auroc_gain_lower_95.gt(0).all() and gains.relative_auprc_gain_lower_95.gt(0).all()):
        raise ValueError("report requires significant AUROC and AUPRC gains over DWM and raw")
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
        suffix = {DWM: "dwm", RAW: "raw", DEEP_BIAS: "chrombpnet_bias"}[row["baseline"]]
        for metric in ("auroc_gain", "relative_auprc_gain"):
            for bound in ("lower_95", "upper_95"):
                table.loc[candidate_index, f"{metric}_{bound}_vs_{suffix}"] = row[
                    f"{metric}_{bound}"
                ]
    for key in SAFETY_COLUMNS:
        table[key] = np.nan
        table.loc[candidate_index, key] = safety[key]
    order = {method: index for index, method in enumerate(METHODS)}
    table["_order"] = table["method"].map(order)
    return table.sort_values("_order").drop(columns="_order").reset_index(drop=True)


def render_report(
    metrics: pd.DataFrame,
    bootstrap: pd.DataFrame,
    profiles: pd.DataFrame,
    safety: pd.Series,
    balance: pd.Series,
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
        f"{cell} {tf}: frozen direct-bias footprint comparison",
        x=0.04,
        y=0.965,
        ha="left",
        fontsize=17,
        weight="bold",
    )
    figure.text(
        0.04,
        0.922,
        f"Development chr19–22/X; identical finite motif sites "
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
        title="Protection shape (pointwise 95% chromosome-bootstrap bands)",
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
        colLabels=["Method", "AUROC", "AUPRC", "RMS sep."],
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
        "Replicate-2 candidate safety\n"
        f"Absolute: {int(safety.false_positive_calls)}/{int(safety.informative_support)} informative sites\n"
        f"DWM-comparable: {int(safety.paired_false_positive_calls)}/{int(safety.paired_support)}; "
        f"upper {100 * safety.paired_false_positive_rate_upper_95:.2f}%\n"
        f"DWM calls: {int(safety.paired_dwm_calls)}; unavailable: {int(safety.dwm_unavailable_informative_sites)}",
        ha="left",
        va="top",
        fontsize=9.5,
        bbox={"boxstyle": "round,pad=0.55", "facecolor": "#E8F3E8", "edgecolor": "#3A7D44"},
    )
    figure.text(
        0.04,
        0.125,
        "Scope and limits",
        fontsize=10,
        weight="bold",
    )
    figure.text(
        0.04,
        0.086,
        f"Observed cuts: replicate 1, 10M, seed 2026. Signed deviance, dispersion {safety.dispersion:.6g}; "
        "same center/shoulder geometry for all three arms.\n"
        "DWM is an expected-bias reference. Scores and displayed residual profiles use different normalization; "
        "this is not a production-command benchmark.",
        fontsize=8.5,
    )
    figure.text(
        0.04,
        0.045,
        f"Motif-score imbalance remains (SMD {balance.motif_score_smd_after:.3f}); accessibility SMD {balance.accessibility_smd_after:.3f}. "
        "Previously examined chromosomes; sparse control counts.\n"
        "No general TF claim, no calibrated functional uncertainty, no final untouched safety. "
        "New holdouts remain unopened; production defaults are unchanged.",
        fontsize=8,
        color="#555555",
    )
    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_pdf, format="pdf", metadata={"CreationDate": None, "ModDate": None})
    figure.savefig(output_png, format="png", dpi=170)
    plt.close(figure)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison-manifest", type=Path, required=True)
    parser.add_argument("--safety-manifest", "--safety", dest="safety_manifest", type=Path, required=True,
                        help="Coverage-qualified direct-operator JSON manifest; legacy safety TSVs are rejected")
    parser.add_argument("--audit-manifest", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, default=Path("."))
    parser.add_argument("--cell", default="K562")
    parser.add_argument("--tf", default="CTCF")
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)

    if args.cell not in {"K562", "HepG2"}:
        raise ValueError("report accepts previously examined development cells only")
    document = json.loads(args.comparison_manifest.read_text(encoding="utf-8"))
    if document.get("schema") != COMPARISON_SCHEMA:
        raise ValueError("unsupported frozen ChromBPNet comparison manifest")
    check_identity(document, "result_id")
    coverage = json.loads(args.safety_manifest.read_text())
    if coverage.get("schema") != SCHEMA_COVERAGE:
        raise ValueError("report requires coverage-qualified direct-operator safety")
    source = json.loads(checked_path(coverage["source_safety"], args.data_root).read_text())
    configuration = json.loads(checked_path(source["configuration"], args.data_root).read_text())
    validate_safety_link(coverage, source, configuration, document, cell=args.cell)
    threshold_freeze, _ = read_thresholds(checked_path(source["thresholds"], args.data_root),
                                          configuration["configuration_id"], args.data_root)
    if threshold_freeze["threshold_id"] != source["threshold_id"]:
        raise ValueError("report threshold identity mismatch")
    audit = json.loads(args.audit_manifest.read_text())
    check_identity(audit, "audit_id")
    if (audit.get("schema") != "fp-tools-chrombpnet-development-audit-v1"
            or audit.get("new_holdout_labels_read") is not False
            or not any(item["result_id"] == document["result_id"]
                       and item["sha256"] == file_sha256(args.comparison_manifest) for item in audit["comparisons"])):
        raise ValueError("shape/balance audit does not match reported comparison")
    tables = {name: pd.read_csv(checked_output(audit, name, args.audit_manifest, args.data_root), sep="\t")
              for name in ("shape_audit", "site_balance")}
    shape = tables["shape_audit"].query("cell == @args.cell and tf == @args.tf and method == @PARAMETRIC")
    if len(shape) != 1 or shape.canonical_protection.iloc[0] != True:
        raise ValueError("report requires canonical central protection")
    balance_rows = tables["site_balance"].query("cell == @args.cell and tf == @args.tf")
    if len(balance_rows) != 1:
        raise ValueError("report requires exact-task covariate balance evidence")
    balance = balance_rows.iloc[0]
    metrics_path = checked_output(document, "metrics", args.comparison_manifest, args.data_root)
    bootstrap_path = checked_output(document, "bootstrap", args.comparison_manifest, args.data_root)
    profiles_path = checked_output(document, "profiles", args.comparison_manifest, args.data_root)
    metrics = pd.read_csv(metrics_path, sep="\t")
    bootstrap = pd.read_csv(bootstrap_path, sep="\t")
    profiles = pd.read_csv(profiles_path, sep="\t")
    safety = pd.read_csv(checked_output(coverage, "qualified_safety", args.safety_manifest, args.data_root), sep="\t")
    selected, paired, safe = select_report_rows(
        metrics, bootstrap, safety, cell=args.cell, tf=args.tf
    )
    if not np.isclose(safe.dispersion, configuration["cells"][args.cell]["dispersion"], rtol=0, atol=1e-14):
        raise ValueError("safety row dispersion differs from the frozen operator")
    summary = build_metrics_table(selected, paired, safe)

    if args.outdir.exists() and any(args.outdir.iterdir()):
        raise ValueError("use a new report directory; historical reports are immutable")
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
        balance,
        cell=args.cell,
        tf=args.tf,
        output_pdf=pdf_path,
        output_png=png_path,
    )
    summary.to_csv(metrics_output, sep="\t", index=False)
    readme_path.write_text(
        f"# {args.cell} {args.tf}: frozen parametric bias comparison\n\n"
        "The one-page report compares DWM expectations, a frozen CPU-compatible "
        "21-bp log-linear sequence-bias model, and the official pinned "
        "ChromBPNet bias network on identical chr19–22/X motif sites and the "
        "same conventional footprint geometry. AUROC and relative AUPRC gains "
        "over both DWM and raw have positive chromosome-block bootstrap intervals, "
        "and the aggregate passes the signed central-protection diagnostic.\n\n"
        f"All expected profiles are matched to each site's observed total before signed deviance "
        f"(dispersion {safe.dispersion:.16g}). Displayed residuals additionally use outer-flank "
        "normalization. DWM here is an expected-bias reference, not the production corrected-track scorer. "
        "ChromBPNet regulatory/no-bias predictions have a different role and are covered in the complete TF audit.\n\n"
        f"Exact-method candidate screening uses validation-only frozen thresholds and existing naked-DNA replicate 2. "
        f"Absolute calls: {int(safe.false_positive_calls)}/{int(safe.informative_support)} informative sites "
        f"(Wilson upper {100 * safe.informative_false_positive_rate_upper_95:.3f}%). "
        f"Paired calls: {int(safe.paired_false_positive_calls)}/{int(safe.paired_support)} DWM-comparable sites "
        f"(upper {100 * safe.paired_false_positive_rate_upper_95:.3f}%, DWM calls {int(safe.paired_dwm_calls)}, "
        f"increase {100 * safe.paired_false_positive_rate_increase_over_dwm:.3f} percentage points). "
        f"{int(safe.dwm_unavailable_informative_sites)} informative site(s) lack a DWM expectation and are excluded "
        "only from the paired denominator. Coverage qualification is post hoc and preserves all original scores and thresholds.\n\n"
        f"Limits: {args.cell} {args.tf}, observed replicate 1 at 10M fragments, seed 2026; "
        "previously examined development chromosomes. Motif-score/accessibility imbalance remains "
        f"(SMD {balance.motif_score_smd_after:.3f}/{balance.accessibility_smd_after:.3f}). "
        "Sparse naked-DNA counts and no large residual flag do not establish absence of occupancy leakage. "
        "Pointwise bootstrap bands do not establish calibrated functional uncertainty. "
        "The two recorded cell comparisons have different dispersions; do not pool them as one common-operator experiment. "
        "No general TF replacement passes promotion gates. Depth, seed, biological-replicate and resource gates "
        "for this exact direct operator remain incomplete. Final naked-DNA replicate 3 and new cell holdout labels "
        "remain unopened. Production defaults are unchanged.\n",
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
        "comparison": {"result_id": document["result_id"], "sha256": file_sha256(args.comparison_manifest)},
        "safety": {"coverage_id": coverage["coverage_id"], "sha256": file_sha256(args.safety_manifest),
                   "configuration_id": configuration["configuration_id"], "threshold_id": source["threshold_id"]},
        "audit": {"audit_id": audit["audit_id"], "sha256": file_sha256(args.audit_manifest)},
        "script_sha256": file_sha256(Path(__file__)),
        "evidence": "development_diagnostic_with_coverage_qualified_candidate_safety",
        "promotion_eligible": False,
        "dispersion": float(safe.dispersion),
        "cell": args.cell,
        "tf": args.tf,
        "outputs": {
            name: {"path": path.name, "sha256": file_sha256(path)}
            for name, path in outputs.items()
        },
    }
    manifest["report_id"] = canonical_id(manifest)
    manifest_path = args.outdir / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (args.outdir / "SHA256SUMS").write_text("".join(f"{file_sha256(path)}  {path.name}\n"
                                                   for path in [*outputs.values(), manifest_path]))
    print(summary[["display_method", "auroc", "auprc"]].to_string(index=False))
    print(pdf_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
