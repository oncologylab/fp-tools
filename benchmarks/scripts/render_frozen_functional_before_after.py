#!/usr/bin/env python3
"""Render concise before/after reports only for safety-qualified TF gains."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import re

import pandas as pd


TEST_SCHEMA = "fp-tools-frozen-functional-test-results-v1"
NAKED_SCHEMA = "fp-tools-frozen-functional-naked-dna-v1"
REPORT_SCHEMA = "fp-tools-frozen-functional-before-after-v1"
DETECTOR_EVIDENCE_SCHEMA = "fp-tools-frozen-detector-evidence-v1"
BASELINE_METHOD = "DWM_conventional_geometry"
RAW_GUARDED_REPORT_CLASSES = {
    "robust_tf_specific_gain",
    "depth_dependent_tf_specific_gain",
}


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
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict, dict]:
    test = json.loads(test_manifest.read_text(encoding="utf-8"))
    naked = json.loads(naked_manifest.read_text(encoding="utf-8"))
    if test.get("schema") != TEST_SCHEMA:
        raise ValueError("unsupported frozen functional test manifest")
    if naked.get("schema") != NAKED_SCHEMA:
        raise ValueError("unsupported frozen functional naked-DNA manifest")
    if test.get("policy_id") != naked.get("policy_id"):
        raise ValueError("test and naked-DNA manifests use different policies")
    if test.get("models_refitted_on_test") is not False:
        raise ValueError("test manifest does not certify no-refit evaluation")
    if naked.get("models_refitted_on_naked_dna") is not False:
        raise ValueError("naked-DNA manifest does not certify frozen models")
    metrics = pd.read_csv(checked_output(test, "metrics", test_manifest), sep="\t")
    bootstrap = pd.read_csv(checked_output(test, "bootstrap", test_manifest), sep="\t")
    profiles = pd.read_csv(checked_output(test, "profiles", test_manifest), sep="\t")
    safety = pd.read_csv(checked_output(naked, "rates", naked_manifest), sep="\t")
    return metrics, bootstrap, profiles, safety, test, naked


def load_detector_evidence(manifest: Path) -> tuple[pd.DataFrame, dict, Path]:
    document = json.loads(manifest.read_text(encoding="utf-8"))
    if document.get("schema") != DETECTOR_EVIDENCE_SCHEMA:
        raise ValueError("unsupported frozen detector evidence manifest")
    evidence_path = checked_output(document, "per_tf_evidence", manifest)
    return pd.read_csv(evidence_path, sep="\t"), document, evidence_path


def apply_raw_guardrail(
    summary: pd.DataFrame,
    evidence: pd.DataFrame,
) -> pd.DataFrame:
    if evidence.duplicated(["cell", "tf"]).any():
        raise ValueError("detector evidence contains duplicate TF tasks")
    keep = [
        "cell",
        "tf",
        "detector_classification",
        "test_raw_auroc",
        "test_raw_auprc",
        "test_auroc",
        "test_auprc",
        "test_auroc_gain_over_raw",
        "test_relative_auprc_gain_over_raw",
        "test_raw_bootstrap_auroc_gain_lower_95",
        "test_raw_bootstrap_auroc_gain_upper_95",
        "test_raw_bootstrap_relative_auprc_gain_lower_95",
        "test_raw_bootstrap_relative_auprc_gain_upper_95",
        "replicate_samples",
        "replicate_auroc_gain_over_raw_positive_fraction",
        "replicate_auprc_gain_over_raw_positive_fraction",
        "depth_both_gain_over_raw_fraction",
        "depth_high_both_gain_over_raw_fraction",
    ]
    missing = {"cell", "tf", "detector_classification"}.difference(evidence.columns)
    if missing:
        raise ValueError("detector evidence lacks columns: " + ", ".join(sorted(missing)))
    output = summary.merge(
        evidence[[column for column in keep if column in evidence]],
        on=["cell", "tf"],
        how="left",
        validate="one_to_one",
    )
    output["passes_raw_guardrail"] = output["detector_classification"].isin(
        RAW_GUARDED_REPORT_CLASSES
    )
    output["report_qualified"] = (
        output["report_qualified"] & output["passes_raw_guardrail"]
    )
    return output


def qualified_candidates(
    metrics: pd.DataFrame,
    bootstrap: pd.DataFrame,
    safety: pd.DataFrame,
) -> pd.DataFrame:
    required_metrics = {
        "cell",
        "tf",
        "motif_family",
        "candidate_id",
        "method",
        "status",
        "n_positive",
        "n_negative",
        "dwm_auroc",
        "dwm_auprc",
        "auroc",
        "auprc",
        "auroc_gain_over_dwm",
        "relative_auprc_gain_over_dwm",
        "functional_separation_relative_change_over_dwm",
    }
    missing = required_metrics.difference(metrics.columns)
    if missing:
        raise ValueError("test metrics lack columns: " + ", ".join(sorted(missing)))
    candidate = metrics[
        metrics["method"].astype(str).str.startswith("frozen_")
    ].copy()
    if candidate.duplicated(["cell", "tf"]).any():
        raise ValueError("test metrics contain duplicate frozen TF candidates")

    required_bootstrap = {
        "cell",
        "tf",
        "method",
        "auroc_gain_lower_95",
        "auroc_gain_upper_95",
        "relative_auprc_gain_lower_95",
        "relative_auprc_gain_upper_95",
    }
    missing = required_bootstrap.difference(bootstrap.columns)
    if missing:
        raise ValueError("bootstrap table lacks columns: " + ", ".join(sorted(missing)))
    frozen_bootstrap = bootstrap[
        bootstrap["method"].astype(str).str.startswith("frozen_")
    ].copy()
    if frozen_bootstrap.duplicated(["cell", "tf"]).any():
        raise ValueError("bootstrap table contains duplicate frozen TF candidates")
    candidate = candidate.merge(
        frozen_bootstrap[
            [
                "cell",
                "tf",
                "method",
                "auroc_gain_lower_95",
                "auroc_gain_upper_95",
                "relative_auprc_gain_lower_95",
                "relative_auprc_gain_upper_95",
            ]
        ],
        on=["cell", "tf", "method"],
        how="left",
        validate="one_to_one",
    )

    required_safety = {
        "cell",
        "tf",
        "method",
        "candidate_false_positive_rate",
        "candidate_wilson_upper_95",
        "candidate_minus_dwm",
        "passes_safety",
    }
    missing = required_safety.difference(safety.columns)
    if missing:
        raise ValueError("naked-DNA table lacks columns: " + ", ".join(sorted(missing)))
    paired = safety[safety["method"].astype(str).eq("paired_safety")].copy()
    if paired.duplicated(["cell", "tf"]).any():
        raise ValueError("naked-DNA table contains duplicate paired TF safety rows")
    candidate = candidate.merge(
        paired[
            [
                "cell",
                "tf",
                "candidate_false_positive_rate",
                "candidate_wilson_upper_95",
                "candidate_minus_dwm",
                "passes_safety",
            ]
        ],
        on=["cell", "tf"],
        how="left",
        validate="one_to_one",
    )
    candidate["passes_safety"] = candidate["passes_safety"].eq(True)
    candidate["significant_auroc_gain"] = candidate["auroc_gain_lower_95"] > 0.0
    candidate["significant_auprc_gain"] = (
        candidate["relative_auprc_gain_lower_95"] > 0.0
    )
    candidate["eligible"] = candidate["status"].astype(str).eq("eligible")
    candidate["report_qualified"] = (
        candidate["eligible"]
        & candidate["passes_safety"]
        & candidate["significant_auroc_gain"]
        & candidate["significant_auprc_gain"]
    )
    return candidate.sort_values(["cell", "tf"]).reset_index(drop=True)


def safe_token(value: object) -> str:
    token = re.sub(r"[^A-Za-z0-9_-]+", "_", str(value)).strip("_-")
    if not token:
        raise ValueError("empty report filename token")
    return token


def _profile_for(
    profiles: pd.DataFrame,
    *,
    cell: str,
    tf: str,
    method: str,
) -> pd.DataFrame:
    selected = profiles[
        profiles["cell"].astype(str).eq(cell)
        & profiles["tf"].astype(str).eq(tf)
        & profiles["method"].astype(str).eq(method)
    ].sort_values("position")
    if selected.empty:
        raise ValueError(f"profile table lacks {cell}/{tf}/{method}")
    if selected["position"].duplicated().any():
        raise ValueError(f"profile table duplicates positions for {cell}/{tf}/{method}")
    return selected


def render_report(row: pd.Series, profiles: pd.DataFrame, output: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cell = str(row["cell"])
    tf = str(row["tf"])
    methods = [BASELINE_METHOD, str(row["method"])]
    titles = ["Conventional DWM", "Frozen parametric + functional"]
    figure, axes = plt.subplots(
        2,
        2,
        figsize=(10.8, 6.7),
        sharex="col",
        gridspec_kw={"height_ratios": [2.0, 1.1]},
    )
    for column, (method, title) in enumerate(zip(methods, titles, strict=True)):
        values = _profile_for(profiles, cell=cell, tf=tf, method=method)
        position = values["position"].to_numpy(dtype=float)
        positive = values["positive_mean"].to_numpy(dtype=float)
        negative = values["negative_mean"].to_numpy(dtype=float)
        difference = values["positive_minus_negative"].to_numpy(dtype=float)
        lower = values["lower_95"].to_numpy(dtype=float)
        upper = values["upper_95"].to_numpy(dtype=float)
        axes[0, column].plot(position, positive, color="#C23B33", label="ChIP positive")
        axes[0, column].plot(position, negative, color="#2A6FBB", label="Matched negative")
        axes[1, column].fill_between(
            position,
            lower,
            upper,
            color="#6A3D9A",
            alpha=0.18,
            linewidth=0,
        )
        axes[1, column].plot(position, difference, color="#6A3D9A", linewidth=1.6)
        for axis in axes[:, column]:
            axis.axvline(0, color="#666666", linewidth=0.7, linestyle="--")
            axis.axhline(0, color="#999999", linewidth=0.6)
        axes[0, column].set_title(title)
        axes[1, column].set_xlabel("Position from motif center (bp)")
    axes[0, 0].set_ylabel("Normalized residual")
    axes[1, 0].set_ylabel("Positive − negative")
    axes[0, 0].legend(frameon=False, fontsize=8)
    metrics = (
        f"AUROC {row['dwm_auroc']:.3f} → {row['auroc']:.3f} "
        f"(Δ {row['auroc_gain_over_dwm']:+.3f}; 95% CI "
        f"{row['auroc_gain_lower_95']:+.3f} to {row['auroc_gain_upper_95']:+.3f})\n"
        f"AUPRC {row['dwm_auprc']:.3f} → {row['auprc']:.3f} "
        f"(relative Δ {100.0 * row['relative_auprc_gain_over_dwm']:+.1f}%; 95% CI "
        f"{100.0 * row['relative_auprc_gain_lower_95']:+.1f}% to "
        f"{100.0 * row['relative_auprc_gain_upper_95']:+.1f}%)\n"
        f"Independent naked-DNA FPR {100.0 * row['candidate_false_positive_rate']:.1f}% "
        f"(Wilson upper {100.0 * row['candidate_wilson_upper_95']:.1f}%)"
    )
    if pd.notna(row.get("test_raw_auroc")):
        metrics += (
            "\nRaw-signal guardrail on exact common support: "
            f"AUROC {row['test_raw_auroc']:.3f} → {row['test_auroc']:.3f} "
            f"(Δ {row['test_auroc_gain_over_raw']:+.3f}; 95% CI "
            f"{row['test_raw_bootstrap_auroc_gain_lower_95']:+.3f} to "
            f"{row['test_raw_bootstrap_auroc_gain_upper_95']:+.3f}); "
            f"relative AUPRC Δ "
            f"{100.0 * row['test_relative_auprc_gain_over_raw']:+.1f}% "
            f"(95% CI "
            f"{100.0 * row['test_raw_bootstrap_relative_auprc_gain_lower_95']:+.1f}% "
            f"to {100.0 * row['test_raw_bootstrap_relative_auprc_gain_upper_95']:+.1f}%)"
        )
    if pd.notna(row.get("replicate_samples")):
        replicate_fraction = min(
            float(row["replicate_auroc_gain_over_raw_positive_fraction"]),
            float(row["replicate_auprc_gain_over_raw_positive_fraction"]),
        )
        metrics += (
            f"\nRaw-guard stability: both metrics positive in "
            f"{100.0 * replicate_fraction:.0f}% of {int(row['replicate_samples'])} "
            f"biological replicates and "
            f"{100.0 * row['depth_high_both_gain_over_raw_fraction']:.0f}% of "
            "25M/50M depth-seed runs."
        )
    figure.suptitle(f"{cell} — {tf}: internal frozen before/after result", y=0.985)
    figure.text(0.5, 0.885, metrics, ha="center", va="top", fontsize=8.2)
    figure.text(
        0.5,
        0.018,
        "Research result on K562/HepG2 chromosome holdout; not external promotion evidence.",
        ha="center",
        fontsize=8,
        color="#555555",
    )
    top = 0.69 if pd.notna(row.get("test_raw_auroc")) else 0.76
    figure.subplots_adjust(top=top, bottom=0.12, left=0.08, right=0.98, hspace=0.10)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, format="pdf")
    plt.close(figure)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--test-manifest", type=Path, required=True)
    parser.add_argument("--naked-manifest", type=Path, required=True)
    parser.add_argument(
        "--detector-evidence-manifest",
        type=Path,
        help="Require a raw-guarded robust/depth-dependent detector classification.",
    )
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)

    metrics, bootstrap, profiles, safety, test, naked = load_inputs(
        args.test_manifest,
        args.naked_manifest,
    )
    summary = qualified_candidates(metrics, bootstrap, safety)
    evidence_document = None
    evidence_path = None
    if args.detector_evidence_manifest is not None:
        evidence, evidence_document, evidence_path = load_detector_evidence(
            args.detector_evidence_manifest
        )
        if evidence_document.get("policy_id") != test.get("policy_id"):
            raise ValueError("report and detector evidence use different policies")
        summary = apply_raw_guardrail(summary, evidence)
    args.outdir.mkdir(parents=True, exist_ok=True)
    summary_path = args.outdir / "before_after_metrics.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    reports = []
    for _index, row in summary[summary["report_qualified"]].iterrows():
        output = args.outdir / (
            f"{safe_token(row['cell'])}_{safe_token(row['tf'])}_before_after.pdf"
        )
        render_report(row, profiles, output)
        reports.append(
            {
                "cell": str(row["cell"]),
                "tf": str(row["tf"]),
                "path": str(output),
                "sha256": file_sha256(output),
            }
        )
    readme = args.outdir / "README.md"
    raw_guard_text = (
        " A raw-signal guardrail was also required: only robust or "
        "depth-dependent gains with biological-replicate and depth/seed "
        "evidence were emitted."
        if evidence_document is not None
        else ""
    )
    readme.write_text(
        "# Frozen parametric before/after reports\n\n"
        "PDFs are emitted only when both chromosome-block bootstrap intervals "
        "exclude zero and the frozen candidate passes independent naked-DNA "
        f"replicate-2 safety.{raw_guard_text} These are internal K562/HepG2 chromosome-holdout "
        "results, not evidence of transfer to the unopened promotion holdouts. "
        "The current DWM package default is unchanged.\n",
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
        "detector_evidence_manifest": (
            {
                "path": str(args.detector_evidence_manifest),
                "sha256": file_sha256(args.detector_evidence_manifest),
                "evidence_path": str(evidence_path),
                "evidence_sha256": file_sha256(evidence_path),
            }
            if evidence_document is not None
            else None
        ),
        "qualification": {
            "eligible": True,
            "auroc_gain_lower_95_gt": 0.0,
            "relative_auprc_gain_lower_95_gt": 0.0,
            "independent_naked_dna_safety": True,
            "raw_signal_guardrail": evidence_document is not None,
        },
        "tasks": int(len(summary)),
        "reports_emitted": int(len(reports)),
        "reports": reports,
        "outputs": {
            "metrics": {"path": str(summary_path), "sha256": file_sha256(summary_path)},
            "readme": {"path": str(readme), "sha256": file_sha256(readme)},
        },
    }
    canonical = json.dumps(manifest, sort_keys=True, separators=(",", ":"))
    manifest["report_id"] = sha256(canonical.encode()).hexdigest()
    manifest_path = args.outdir / "before_after_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(summary[["cell", "tf", "report_qualified"]].to_string(index=False))
    print(f"reports emitted: {len(reports)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
