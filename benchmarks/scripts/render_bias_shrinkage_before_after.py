#!/usr/bin/env python3
"""Render concise reports for safety-qualified frozen bias-shrinkage gains."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd


TEST_SCHEMA = "fp-tools-frozen-bias-shrinkage-test-v1"
NAKED_SCHEMA = "fp-tools-frozen-bias-shrinkage-naked-dna-v1"
POLICY_SCHEMA = "fp-tools-frozen-bias-shrinkage-policy-v1"
REPORT_SCHEMA = "fp-tools-bias-shrinkage-before-after-v1"
BASELINE_METHOD = "DWM_conventional_deviance"
CANDIDATE_METHOD = "frozen_tf_specific_shrinkage"


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def checked_output(document: dict, name: str, source: Path) -> Path:
    record = document.get("outputs", {}).get(name, {})
    path = Path(record.get("path", ""))
    if not path.is_file():
        raise ValueError(f"{source} lacks the declared {name} output")
    if file_sha256(path) != record.get("sha256"):
        raise ValueError(f"{source} {name} checksum mismatch")
    return path


def load_inputs(
    test_manifest: Path,
    naked_manifest: Path,
    policy_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict, dict, dict]:
    test = json.loads(test_manifest.read_text(encoding="utf-8"))
    naked = json.loads(naked_manifest.read_text(encoding="utf-8"))
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    if test.get("schema") != TEST_SCHEMA:
        raise ValueError("unsupported frozen bias-shrinkage test manifest")
    if naked.get("schema") != NAKED_SCHEMA:
        raise ValueError("unsupported frozen bias-shrinkage naked-DNA manifest")
    if policy.get("schema") != POLICY_SCHEMA:
        raise ValueError("unsupported frozen bias-shrinkage policy")
    policy_ids = {test.get("policy_id"), naked.get("policy_id"), policy.get("policy_id")}
    if len(policy_ids) != 1:
        raise ValueError("test, naked-DNA, and policy artifacts do not share one policy")
    if test.get("models_refitted_on_test") is not False:
        raise ValueError("test manifest does not certify no-refit evaluation")
    if test.get("raw_guardrail") is not True:
        raise ValueError("test manifest does not certify the raw-signal guardrail")
    if naked.get("models_refitted_on_naked_dna") is not False:
        raise ValueError("naked-DNA manifest does not certify frozen models")
    if naked.get("thresholds_changed_on_naked_dna") is not False:
        raise ValueError("naked-DNA thresholds were not held frozen")
    metrics = pd.read_csv(checked_output(test, "metrics", test_manifest), sep="\t")
    bootstrap = pd.read_csv(
        checked_output(test, "bootstrap", test_manifest), sep="\t"
    )
    profiles = pd.read_csv(checked_output(test, "profiles", test_manifest), sep="\t")
    safety = pd.read_csv(checked_output(naked, "rates", naked_manifest), sep="\t")
    return metrics, bootstrap, profiles, safety, test, naked, policy


def qualified_candidates(
    metrics: pd.DataFrame,
    bootstrap: pd.DataFrame,
    safety: pd.DataFrame,
    *,
    minimum_contexts: int,
) -> pd.DataFrame:
    if minimum_contexts < 1:
        raise ValueError("minimum contexts must be positive")
    required_metrics = {
        "cell",
        "tf",
        "motif_family",
        "method",
        "status",
        "positive_sites",
        "negative_sites",
        "auroc",
        "auprc",
        "raw_auroc",
        "raw_auprc",
        "dwm_auroc",
        "dwm_auprc",
        "auroc_gain_over_dwm",
        "relative_auprc_gain_over_dwm",
    }
    missing = required_metrics.difference(metrics.columns)
    if missing:
        raise ValueError("test metrics lack columns: " + ", ".join(sorted(missing)))
    candidate = metrics[metrics["method"].astype(str).eq(CANDIDATE_METHOD)].copy()
    if candidate.duplicated(["cell", "tf"]).any():
        raise ValueError("test metrics contain duplicate TF-specific candidates")

    required_bootstrap = {
        "cell",
        "tf",
        "method",
        "baseline",
        "auroc_gain_lower_95",
        "auroc_gain_upper_95",
        "relative_auprc_gain_lower_95",
        "relative_auprc_gain_upper_95",
    }
    missing = required_bootstrap.difference(bootstrap.columns)
    if missing:
        raise ValueError("bootstrap table lacks columns: " + ", ".join(sorted(missing)))
    paired_bootstrap = bootstrap[
        bootstrap["method"].astype(str).eq(CANDIDATE_METHOD)
        & bootstrap["baseline"].astype(str).eq(BASELINE_METHOD)
    ].copy()
    if paired_bootstrap.duplicated(["cell", "tf"]).any():
        raise ValueError("bootstrap table contains duplicate DWM comparisons")
    candidate = candidate.merge(
        paired_bootstrap[
            [
                "cell",
                "tf",
                "auroc_gain_lower_95",
                "auroc_gain_upper_95",
                "relative_auprc_gain_lower_95",
                "relative_auprc_gain_upper_95",
            ]
        ],
        on=["cell", "tf"],
        how="left",
        validate="one_to_one",
    )

    required_safety = {
        "cell",
        "tf",
        "method",
        "finite_sites",
        "calls",
        "false_positive_rate",
        "false_positive_rate_upper_95",
        "false_positive_rate_increase_over_dwm",
        "passes_safety",
    }
    missing = required_safety.difference(safety.columns)
    if missing:
        raise ValueError("naked-DNA table lacks columns: " + ", ".join(sorted(missing)))
    candidate_safety = safety[
        safety["method"].astype(str).eq(CANDIDATE_METHOD)
    ].copy()
    if candidate_safety.duplicated(["cell", "tf"]).any():
        raise ValueError("naked-DNA table contains duplicate TF-specific safety rows")
    candidate = candidate.merge(
        candidate_safety[
            [
                "cell",
                "tf",
                "finite_sites",
                "calls",
                "false_positive_rate",
                "false_positive_rate_upper_95",
                "false_positive_rate_increase_over_dwm",
                "passes_safety",
            ]
        ],
        on=["cell", "tf"],
        how="left",
        validate="one_to_one",
    )
    candidate["significant_dwm_gain"] = (
        (candidate["auroc_gain_lower_95"] > 0.0)
        & (candidate["relative_auprc_gain_lower_95"] > 0.0)
    )
    candidate["qualified_context"] = (
        candidate["status"].astype(str).eq("eligible")
        & candidate["passes_safety"].eq(True)
        & candidate["significant_dwm_gain"]
    )
    counts = (
        candidate.groupby("tf")["qualified_context"]
        .sum()
        .astype(int)
        .rename("qualified_contexts")
    )
    candidate = candidate.merge(counts, on="tf", how="left", validate="many_to_one")
    candidate["report_qualified"] = (
        candidate["qualified_context"]
        & (candidate["qualified_contexts"] >= minimum_contexts)
    )
    return candidate.sort_values(["tf", "cell"]).reset_index(drop=True)


def safe_token(value: object) -> str:
    token = re.sub(r"[^A-Za-z0-9_-]+", "_", str(value)).strip("_-")
    if not token:
        raise ValueError("empty report filename token")
    return token


def profile_curve(
    profiles: pd.DataFrame,
    *,
    cell: str,
    tf: str,
    method: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    selected = profiles[
        profiles["cell"].astype(str).eq(cell)
        & profiles["tf"].astype(str).eq(tf)
        & profiles["method"].astype(str).eq(method)
    ].sort_values("position")
    if selected.empty or selected["position"].duplicated().any():
        raise ValueError(f"profile table is incomplete for {cell}/{tf}/{method}")
    position = selected["position"].to_numpy(dtype=float)
    curve = selected["positive_minus_negative"].to_numpy(dtype=float)
    lower = selected["lower_95"].to_numpy(dtype=float)
    upper = selected["upper_95"].to_numpy(dtype=float)
    flank = (np.abs(position) >= 50) & (np.abs(position) <= 100)
    center = float(np.nanmean(curve[flank]))
    scale = float(np.sqrt(np.nanmean(np.square(curve[flank] - center))))
    if not np.isfinite(scale) or scale <= 1e-12:
        scale = 1.0
    return position, (curve - center) / scale, (lower - center) / scale, (
        upper - center
    ) / scale


def render_report(
    rows: pd.DataFrame,
    profiles: pd.DataFrame,
    *,
    alpha: float,
    source: str,
    output: Path,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = rows.sort_values("cell").reset_index(drop=True)
    tf = str(rows.iloc[0]["tf"])
    figure, axes = plt.subplots(
        1,
        len(rows),
        figsize=(11.0, 8.5),
        squeeze=False,
    )
    axes = axes[0]
    for axis, row in zip(axes, rows.to_dict("records"), strict=True):
        for method, label, color in (
            (BASELINE_METHOD, "Before: conventional DWM", "#666666"),
            (CANDIDATE_METHOD, "After: frozen parametric", "#6A3D9A"),
        ):
            position, curve, lower, upper = profile_curve(
                profiles,
                cell=str(row["cell"]),
                tf=tf,
                method=method,
            )
            axis.fill_between(position, lower, upper, color=color, alpha=0.12)
            axis.plot(position, curve, color=color, linewidth=1.8, label=label)
        axis.axvspan(-15, 15, color="#E6AB02", alpha=0.08, linewidth=0)
        axis.axvline(0, color="#444444", linewidth=0.7, linestyle="--")
        axis.axhline(0, color="#999999", linewidth=0.6)
        axis.set_title(str(row["cell"]), fontsize=12)
        axis.set_xlabel("Position from motif center (bp)")
        axis.grid(axis="y", alpha=0.18, linewidth=0.6)
    axes[0].set_ylabel("Flank-standardized ChIP-positive − matched-negative")
    axes[0].legend(frameon=False, fontsize=8, loc="lower left")

    figure.suptitle(
        f"{tf} footprint detection: conventional DWM vs frozen parametric correction",
        fontsize=15,
        y=0.965,
    )
    source_label = {
        "parametric_lambda": "sample-calibrated",
        "parametric_direct": "direct local-total-scaled",
    }.get(source, source.replace("parametric_", "").replace("_", " "))
    figure.text(
        0.5,
        0.925,
        f"TF-specific validation choice: subtract {100.0 * alpha:.0f}% of the "
        f"{source_label} expected-bias signal; coefficients remain frozen",
        ha="center",
        fontsize=9,
    )

    table_rows = []
    for row in rows.to_dict("records"):
        table_rows.append(
            [
                str(row["cell"]),
                f"{row['dwm_auroc']:.3f} → {row['auroc']:.3f}",
                f"{row['auroc_gain_over_dwm']:+.3f} "
                f"[{row['auroc_gain_lower_95']:+.3f}, {row['auroc_gain_upper_95']:+.3f}]",
                f"{row['dwm_auprc']:.3f} → {row['auprc']:.3f}",
                f"{100.0 * row['relative_auprc_gain_over_dwm']:+.1f}% "
                f"[{100.0 * row['relative_auprc_gain_lower_95']:+.1f}, "
                f"{100.0 * row['relative_auprc_gain_upper_95']:+.1f}]",
                f"{row['raw_auroc']:.3f} / {row['raw_auprc']:.3f}",
                f"{100.0 * row['false_positive_rate']:.1f}% "
                f"(upper {100.0 * row['false_positive_rate_upper_95']:.1f}%)",
            ]
        )
    table = figure.add_axes([0.025, 0.14, 0.95, 0.16])
    table.axis("off")
    rendered = table.table(
        cellText=table_rows,
        colLabels=[
            "Cell",
            "AUROC DWM → new",
            "AUROC Δ [95% CI]",
            "AUPRC DWM → new",
            "AUPRC relative Δ [95% CI]",
            "Raw AUROC / AUPRC",
            "Naked-DNA FPR",
        ],
        loc="center",
        cellLoc="center",
    )
    rendered.auto_set_font_size(False)
    rendered.set_fontsize(7.0)
    rendered.scale(1.0, 1.65)
    for (row_index, _column), cell in rendered.get_celld().items():
        if row_index == 0:
            cell.set_facecolor("#EAEAEA")
            cell.set_text_props(weight="bold")
        cell.set_edgecolor("#BBBBBB")
        cell.set_linewidth(0.5)

    figure.text(
        0.5,
        0.075,
        "Result: significant improvement over conventional DWM in both internal cell contexts; "
        "independent naked-DNA replicate-2 safety passed.",
        ha="center",
        fontsize=9,
        weight="bold",
    )
    figure.text(
        0.5,
        0.04,
        "Limitation: CTCF-specific internal chromosome-holdout result—not a general method, "
        "not external promotion evidence, and not a package-default change.",
        ha="center",
        fontsize=8,
        color="#555555",
    )
    figure.subplots_adjust(top=0.88, bottom=0.34, left=0.08, right=0.98, wspace=0.20)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, format="pdf")
    plt.close(figure)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--test-manifest", type=Path, required=True)
    parser.add_argument("--naked-manifest", type=Path, required=True)
    parser.add_argument("--policy", type=Path, required=True)
    parser.add_argument("--minimum-contexts", type=int, default=2)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)

    metrics, bootstrap, profiles, safety, test, naked, policy = load_inputs(
        args.test_manifest,
        args.naked_manifest,
        args.policy,
    )
    summary = qualified_candidates(
        metrics,
        bootstrap,
        safety,
        minimum_contexts=args.minimum_contexts,
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    summary_path = args.outdir / "metrics.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    reports = []
    qualifying_tfs = sorted(summary.loc[summary["report_qualified"], "tf"].unique())
    for tf in qualifying_tfs:
        rows = summary[
            summary["tf"].astype(str).eq(str(tf)) & summary["report_qualified"]
        ]
        choice = policy.get("per_tf_choices", {}).get(str(tf))
        if not choice:
            raise ValueError(f"policy lacks a TF-specific choice for {tf}")
        output = args.outdir / f"{safe_token(tf)}_before_after.pdf"
        render_report(
            rows,
            profiles,
            alpha=float(choice["alpha"]),
            source=str(choice["source"]),
            output=output,
        )
        reports.append(
            {
                "tf": str(tf),
                "contexts": sorted(rows["cell"].astype(str).tolist()),
                "path": str(output),
                "sha256": file_sha256(output),
            }
        )

    readme = args.outdir / "README.md"
    readme.write_text(
        "# CTCF frozen parametric bias-shrinkage result\n\n"
        "The report compares the conventional DWM residual with a TF-specific, "
        "validation-frozen partial subtraction of a control-trained parametric "
        "expected-bias signal. It is emitted only when chromosome-block bootstrap "
        "intervals for both AUROC and AUPRC exclude zero in at least two cell "
        "contexts and the exact frozen policy passes independent naked-DNA "
        "replicate-2 safety.\n\n"
        "This is an internal K562/HepG2 chromosome-holdout result for CTCF. It "
        "does not establish a general correction method, it has not been tested "
        "on the unopened SK-N-SH or GM23338 promotion holdouts, and it does not "
        "justify changing the fp-tools default.\n",
        encoding="utf-8",
    )
    manifest = {
        "schema": REPORT_SCHEMA,
        "policy_id": test["policy_id"],
        "test_manifest": {
            "path": str(args.test_manifest),
            "sha256": file_sha256(args.test_manifest),
        },
        "naked_manifest": {
            "path": str(args.naked_manifest),
            "sha256": file_sha256(args.naked_manifest),
        },
        "policy": {"path": str(args.policy), "sha256": file_sha256(args.policy)},
        "qualification": {
            "minimum_contexts": args.minimum_contexts,
            "auroc_gain_lower_95_gt": 0.0,
            "relative_auprc_gain_lower_95_gt": 0.0,
            "independent_naked_dna_safety": True,
        },
        "reports": reports,
        "outputs": {
            "metrics": {"path": str(summary_path), "sha256": file_sha256(summary_path)},
            "readme": {"path": str(readme), "sha256": file_sha256(readme)},
        },
    }
    canonical = json.dumps(manifest, sort_keys=True, separators=(",", ":"))
    manifest["report_id"] = sha256(canonical.encode()).hexdigest()
    manifest_path = args.outdir / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(summary[["cell", "tf", "report_qualified"]].to_string(index=False))
    print(f"reports emitted: {len(reports)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
