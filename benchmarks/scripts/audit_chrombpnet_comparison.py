#!/usr/bin/env python3
"""Audit completed development comparisons without refitting or unlocking labels.

All six existing methods retain their original scores, sites, dispersion and
bootstrap results. This is a diagnostic report, never a promotion evaluator.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))
sys.path.insert(0, str(SCRIPT_DIR))

PINNED_COMMIT = "09938fdb4397ec0006510e5251e48920a505d4de"
KEYS = ["cell", "tf", "method"]
SITE_KEYS = ["cell", "tf", "TFBS_chr", "TFBS_start", "TFBS_end", "TFBS_strand"]
PARAMETRIC = "frozen_parametric_bias_conventional_geometry"
DWM = "DWM_conventional_geometry"
METHODS = {
    "raw_geometry": ("Raw observed cuts", "cuts", "#333333"),
    DWM: ("DWM expected + deviance", "normalized deviance", "#777777"),
    PARAMETRIC: ("Frozen parametric bias + deviance", "normalized deviance", "#D55E00"),
    "ChromBPNet_bias_conventional_geometry": (
        "ChromBPNet bias + deviance", "normalized deviance", "#0072B2"
    ),
    "ChromBPNet_regulatory_geometry": (
        "ChromBPNet regulatory prediction", "total-matched predicted cuts", "#009E73"
    ),
    "ChromBPNet_nobias_geometry": (
        "ChromBPNet no-bias prediction", "total-matched predicted cuts", "#CC79A7"
    ),
}


def canonical_id(document: dict) -> str:
    return sha256(json.dumps(document, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def data_path(value: str | Path, data_root: Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    candidate = data_root / path
    # Versioned study specifications are in the research worktree, while the
    # large immutable artifacts can remain in the original checkout.
    return candidate if candidate.is_file() else REPOSITORY / path


def checked_path(record: dict, data_root: Path) -> Path:
    path = data_path(record["path"], data_root)
    if not path.is_file() or file_sha256(path) != record.get("sha256"):
        raise ValueError(f"artifact checksum mismatch: {path}")
    if "bytes" in record and path.stat().st_size != record["bytes"]:
        raise ValueError(f"artifact size mismatch: {path}")
    return path


def check_identity(document: dict, field: str) -> None:
    value = {key: item for key, item in document.items() if key != field}
    if canonical_id(value) != document.get(field):
        raise ValueError(f"invalid canonical {field}")


def load_comparison(path: Path, data_root: Path) -> tuple[dict, dict, dict]:
    document = json.loads(path.read_text())
    freeze = json.loads(path.with_name("chrombpnet_comparison_inputs.freeze.json").read_text())
    if (not freeze.get("cells") or not set(freeze["cells"]).issubset({"K562", "HepG2"})
            or freeze.get("test_labels_previously_opened_for_other_models") is not True):
        raise ValueError("audit accepts only previously opened development comparisons")
    if document.get("schema") != "fp-tools-chrombpnet-frozen-comparison-v1":
        raise ValueError("unsupported comparison schema")
    if freeze.get("schema") != "fp-tools-chrombpnet-comparison-inputs-v1":
        raise ValueError("unsupported comparison input schema")
    check_identity(document, "result_id")
    check_identity(freeze, "input_id")
    if document["input_id"] != freeze["input_id"]:
        raise ValueError("comparison input identity mismatch")
    if any(item.get("pinned_chrombpnet_commit") != PINNED_COMMIT for item in (document, freeze)):
        raise ValueError("comparison does not use pinned ChromBPNet source")
    for name in ("common_finite_support", "prediction_totals_matched_to_observed",
                 "frozen_parametric_expected_from_test_artifact"):
        if freeze.get(name) is not True:
            raise ValueError(f"comparison does not certify {name}")
    for record in freeze["inputs"]:
        checked_path(record, data_root)
    paths = {key: checked_path(record, data_root) for key, record in document["outputs"].items()}
    required = {"metrics", "profiles", "bootstrap", "site_scores", "coverage"}
    if not required.issubset(paths):
        raise ValueError("comparison lacks required outputs")
    tables = {key: pd.read_csv(paths[key], sep="\t") for key in required}
    scores = tables["site_scores"]
    metrics = tables["metrics"]
    if scores.duplicated(SITE_KEYS).any() or metrics.duplicated(KEYS).any():
        raise ValueError("duplicate comparison task or site")
    if set(scores["cell"]) != set(freeze["cells"]):
        raise ValueError("comparison output cell differs from frozen input")
    if not np.isfinite(scores[list(METHODS)].to_numpy(dtype=float)).all():
        raise ValueError("comparison scores are not finite on common support")
    if not scores["chip_label"].isin([0, 1]).all():
        raise ValueError("comparison labels are not binary")
    for (cell, tf), task in scores.groupby(["cell", "tf"]):
        rows = metrics[metrics.cell.eq(cell) & metrics.tf.eq(tf)]
        if set(rows.method) != set(METHODS):
            raise ValueError(f"missing comparison methods for {cell} {tf}")
        if (not rows.n_sites.eq(len(task)).all()
                or not rows.n_positive.eq(task.chip_label.eq(1).sum()).all()
                or not rows.n_negative.eq(task.chip_label.eq(0).sum()).all()):
            raise ValueError("metric support differs from common site support")
    return document, freeze, tables


def standardized_difference(first, second) -> float:
    first = np.asarray(first, dtype=float)
    second = np.asarray(second, dtype=float)
    first, second = first[np.isfinite(first)], second[np.isfinite(second)]
    if len(first) < 2 or len(second) < 2:
        return np.nan
    difference = float(first.mean() - second.mean())
    scale = np.sqrt((first.var(ddof=1) + second.var(ddof=1)) / 2)
    if scale == 0:
        return 0.0 if difference == 0 else float(np.copysign(np.inf, difference))
    return float(difference / scale)


def balance_audit(source: pd.DataFrame, retained: pd.DataFrame, totals: np.ndarray) -> pd.DataFrame:
    if source.duplicated(SITE_KEYS).any() or retained.duplicated(SITE_KEYS).any():
        raise ValueError("duplicate source or retained site")
    if len(totals) != len(source):
        raise ValueError("raw totals do not align to source sites")
    source = source.copy()
    source["raw_total"] = totals
    merged = retained[SITE_KEYS + ["chip_label"]].merge(
        source, on=SITE_KEYS, how="left", validate="one_to_one",
        suffixes=("_retained", ""), indicator=True,
    )
    if not merged["_merge"].eq("both").all():
        raise ValueError("retained sites missing from frozen input")
    if not merged.chip_label.eq(merged.chip_label_retained).all():
        raise ValueError("retained label differs from frozen input label")
    rows = []
    for (cell, tf), before in source.groupby(["cell", "tf"], sort=True):
        after = merged[merged.cell.eq(cell) & merged.tf.eq(tf)]
        row = {"cell": cell, "tf": tf}
        for label, name in ((1, "positive"), (0, "negative")):
            original = before.chip_label.eq(label).sum()
            remaining = after.chip_label.eq(label).sum()
            row[f"{name}_original"] = int(original)
            row[f"{name}_retained"] = int(remaining)
            row[f"{name}_retained_fraction"] = float(remaining / original) if original else np.nan
        row["retention_difference"] = row["positive_retained_fraction"] - row["negative_retained_fraction"]
        for covariate in ("motif_score", "accessibility", "central_accessibility", "raw_total"):
            if covariate not in source:
                continue
            for name, frame in (("before", before), ("after", after)):
                pos = frame.loc[frame.chip_label.eq(1), covariate].to_numpy(dtype=float)
                neg = frame.loc[frame.chip_label.eq(0), covariate].to_numpy(dtype=float)
                row[f"{covariate}_smd_{name}"] = standardized_difference(pos, neg)
                if name == "after":
                    row[f"{covariate}_positive_mean"] = float(np.mean(pos)) if len(pos) else np.nan
                    row[f"{covariate}_negative_mean"] = float(np.mean(neg)) if len(neg) else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def shape_audit(profiles: pd.DataFrame) -> pd.DataFrame:
    if not np.isfinite(profiles[["position", "positive_mean", "negative_mean",
                                 "positive_minus_negative"]].to_numpy(dtype=float)).all():
        raise ValueError("aggregate profiles must contain finite values")
    rows = []
    for key, group in profiles.groupby(KEYS, sort=True):
        group = group.sort_values("position")
        if group.position.duplicated().any():
            raise ValueError("duplicate aggregate position")
        center = group[group.position.abs() <= 5]
        left = group[group.position.between(-40, -15)]
        right = group[group.position.between(15, 40)]
        if len(center) != 11 or len(left) != 26 or len(right) != 26:
            raise ValueError("aggregate does not cover the fixed shape windows")
        c = float(center.positive_minus_negative.mean())
        l = float(left.positive_minus_negative.mean())
        r = float(right.positive_minus_negative.mean())
        depletion = (l + r) / 2 - c
        significant = int((center.upper_95 < 0).sum())
        denominator = abs(l) + abs(r)
        rows.append({
            **dict(zip(KEYS, key)), "center_contrast": c,
            "left_shoulder": l, "right_shoulder": r, "depletion": depletion,
            "shoulder_asymmetry": (r - l) / denominator if denominator > 0 else 0.0,
            "significant_center_points": significant,
            # This reproduces the existing diagnostic strong shape screen;
            # pointwise bands do not establish simultaneous functional coverage.
            "canonical_protection": bool(depletion >= 0.02 and c < 0 and significant >= 2),
            "both_shoulders_positive": bool(l > 0 and r > 0),
            "band_interpretation": "pointwise_chromosome_bootstrap_diagnostic",
        })
    return pd.DataFrame(rows)


def task_decisions(metrics: pd.DataFrame, shapes: pd.DataFrame, bootstrap: pd.DataFrame) -> pd.DataFrame:
    selected = metrics[metrics.method.eq(PARAMETRIC)].merge(
        shapes, on=KEYS, how="left", validate="one_to_one", indicator=True
    )
    if not selected["_merge"].eq("both").all():
        raise ValueError("missing candidate shape evidence")
    selected = selected.drop(columns="_merge")
    rows = []
    for row in selected.to_dict("records"):
        paired = bootstrap[bootstrap.cell.eq(row["cell"]) & bootstrap.tf.eq(row["tf"])
                           & bootstrap.method.eq(PARAMETRIC)]
        for baseline, label in ((DWM, "dwm"), ("raw_geometry", "raw")):
            pair = paired[paired.baseline.eq(baseline)]
            if len(pair) > 1:
                raise ValueError("duplicate paired bootstrap comparison")
            row[f"significant_vs_{label}"] = bool(
                len(pair) == 1 and pair.iloc[0].auroc_gain_lower_95 > 0
                and pair.iloc[0].relative_auprc_gain_lower_95 > 0
            )
        if row["status"] != "eligible":
            interpretation = "underpowered"
        elif not row["canonical_protection"]:
            interpretation = "no_canonical_protection"
        elif not row["significant_vs_dwm"]:
            interpretation = "no_significant_gain_over_dwm"
        elif not row["significant_vs_raw"]:
            interpretation = "no_significant_gain_over_raw"
        else:
            interpretation = "internal_shape_and_discrimination_gain_safety_unverified"
        row.update(interpretation=interpretation, evidence="development_diagnostic",
                   exact_operator_safety="not_established", promotion_eligible=False)
        rows.append(row)
    return pd.DataFrame(rows)


def load_source_support(freeze: dict, data_root: Path) -> tuple[pd.DataFrame, np.ndarray]:
    from build_strand_functional_profiles import site_hashes

    sources, totals = [], []
    for cell in freeze["cells"]:
        records = [r for r in freeze["inputs"] if r["purpose"] == f"{cell}-test-manifest"]
        if len(records) != 1:
            raise ValueError("expected one frozen test artifact per cell")
        document = json.loads(checked_path(records[0], data_root).read_text())
        if (document.get("schema") != "fp-tools-strand-functional-profiles-v1"
                or document.get("metadata", {}).get("labels_used") is not False
                or document.get("metadata", {}).get("cell") != cell):
            raise ValueError("invalid label-free development profile artifact")
        source_path = checked_path({"path": document["sites"], "sha256": document["sites_sha256"]}, data_root)
        array_path = checked_path({"path": document["profiles_npz"], "sha256": document["profiles_sha256"]}, data_root)
        source = pd.read_csv(source_path, sep="\t")
        if (set(source.cell) != {cell} or not source.chip_label.isin([0, 1]).all()
                or not set(source.TFBS_chr).issubset({"chr19", "chr20", "chr21", "chr22", "chrX"})):
            raise ValueError("source sites are not the recorded development chromosome test")
        with np.load(array_path, allow_pickle=False) as arrays:
            if not np.array_equal(site_hashes(source), arrays["site_hash"]):
                raise ValueError("source site identities do not align to profile arrays")
            observed = arrays["plus_observed"] + arrays["minus_observed"]
            totals.append(observed.sum(axis=1))
        sources.append(source)
    return pd.concat(sources, ignore_index=True), np.concatenate(totals)


def render_panels(profiles: pd.DataFrame, metrics: pd.DataFrame, dispersion: float, path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    with PdfPages(path) as pdf:
        for (cell, tf), task in profiles.groupby(["cell", "tf"], sort=True):
            figure, axes = plt.subplots(2, 3, figsize=(11.69, 8.27))
            figure.suptitle(f"{cell} {tf}: completed development comparison (dispersion {dispersion:.6g})", fontsize=14)
            for axis, (method, (label, unit, color)) in zip(axes.flat, METHODS.items()):
                curve = task[task.method.eq(method)].sort_values("position")
                metric = metrics[metrics.cell.eq(cell) & metrics.tf.eq(tf) & metrics.method.eq(method)].iloc[0]
                x = curve.position.to_numpy()
                axis.plot(x, curve.positive_minus_negative, color=color, linewidth=1.4)
                axis.fill_between(x, curve.lower_95, curve.upper_95, color=color, alpha=0.18)
                axis.axhline(0, color="#999999", linewidth=0.6)
                axis.axvspan(-5, 5, color="#CCCCCC", alpha=0.2)
                axis.set_title(f"{label}\nAUROC {metric.auroc:.3f} · AUPRC {metric.auprc:.3f}", fontsize=9)
                axis.set(xlim=(-100, 100), xlabel="Position relative to motif center (bp)",
                         ylabel=f"Positive − negative ({unit})")
                axis.tick_params(labelsize=8)
                axis.spines[["top", "right"]].set_visible(False)
            figure.text(0.05, 0.04,
                        f"Same sites: {int(metric.n_positive):,} positive / {int(metric.n_negative):,} negative; {metric.status}. "
                        "Shading: pointwise 95% chromosome-bootstrap bands.\n"
                        "Panels use their stated units and independent y scales. Regulatory/no-bias panels are predictions. "
                        "No exact-method safety or promotion claim.", fontsize=8)
            figure.subplots_adjust(left=0.08, right=0.98, top=0.87, bottom=0.15, wspace=0.38, hspace=0.44)
            pdf.savefig(figure)
            plt.close(figure)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison-manifest", type=Path, action="append", required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.outdir.exists():
        raise ValueError("use a new output directory; completed audits are immutable")
    loaded = [load_comparison(path, args.data_root) for path in args.comparison_manifest]
    cells = [cell for _, freeze, _ in loaded for cell in freeze["cells"]]
    if len(cells) != len(set(cells)):
        raise ValueError("multiple comparisons for the same cell must be audited separately")
    args.outdir.mkdir(parents=True)
    all_metrics, all_shapes, all_balance, all_decisions = [], [], [], []
    outputs = {}
    for document, freeze, tables in loaded:
        source, totals = load_source_support(freeze, args.data_root)
        balance = balance_audit(source, tables["site_scores"], totals)
        shapes = shape_audit(tables["profiles"])
        metrics = tables["metrics"].copy()
        metrics["dispersion"] = float(freeze["dispersion"])
        metrics["input_id"] = freeze["input_id"]
        metrics["result_id"] = document["result_id"]
        metrics["profile_units"] = metrics.method.map({key: values[1] for key, values in METHODS.items()})
        decisions = task_decisions(metrics, shapes, tables["bootstrap"])
        all_metrics.append(metrics)
        all_shapes.append(shapes)
        all_balance.append(balance)
        all_decisions.append(decisions)
        name = "_".join(freeze["cells"]) + "_all_TF_aggregate_audit"
        outputs[name] = args.outdir / f"{name}.pdf"
        render_panels(tables["profiles"], metrics, freeze["dispersion"], outputs[name])
    for name, frames in (("per_tf_metrics", all_metrics), ("shape_audit", all_shapes),
                         ("site_balance", all_balance), ("parametric_decisions", all_decisions)):
        outputs[name] = args.outdir / f"{name}.tsv"
        pd.concat(frames, ignore_index=True).to_csv(outputs[name], sep="\t", index=False)
    outputs["limitations"] = args.outdir / "README.md"
    dispersions = {cell: freeze["dispersion"] for _, freeze, _ in loaded for cell in freeze["cells"]}
    outputs["limitations"].write_text(
        "# Completed ChromBPNet comparison audit\n\n"
        "This is a diagnostic audit of previously opened K562/HepG2 development "
        "chromosomes chr19–22/X. All six methods use identical finite motif sites "
        "within each cell. Original inputs, outputs and canonical manifest hashes "
        "were verified. Scores and bootstrap results were not recomputed.\n\n"
        f"Recorded dispersions: {json.dumps(dispersions, sort_keys=True)}. "
        "Different dispersions must not be pooled as a common correction experiment.\n\n"
        "DWM here means DWM expected cuts, total-matched to the observed profile, "
        "then signed-deviance residuals and conventional center/shoulder geometry. "
        "It is a controlled bias-model comparison, not a direct run of the "
        "production corrected-track footprint command. Parametric and ChromBPNet "
        "bias arms use the same operator within each comparison. Regulatory and "
        "no-bias arms score total-matched model predictions, not corrected observed "
        "ATAC. They are scientific references and change more than the bias model.\n\n"
        "The raw/prediction aggregate units are cuts; correction aggregates are "
        "outer-flank normalized residuals. Each panel labels its units. Functional "
        "separation is an unsigned RMS contrast and can rise with central enrichment. "
        "The separate shape audit checks signed center, shoulders and the existing "
        "diagnostic protection criterion. Pointwise bands do not certify calibrated "
        "simultaneous uncertainty. Balance diagnostics report positive-minus-negative "
        "standardized differences; no covariate matching or outcome selection is "
        "performed by this audit. Raw total cuts also reflect biological occupancy.\n\n"
        "The earlier one-page report cited factorization-residual naked-DNA safety. "
        "That score uses a fitted sample bias strength and per-site flank background; "
        "it is not the direct parametric residual used here. Exact-operator safety "
        "is therefore not established by that citation. No method is promoted by "
        "this report. Untouched naked-DNA replicate 3 and new SK-N-SH/GM23338 "
        "labels were not read. Depth/replicate stability, CPU/memory limits, model "
        "size, occupancy leakage, uncertainty and every other locked gate remain "
        "separate required evidence.\n"
    )
    manifest = {
        "schema": "fp-tools-chrombpnet-development-audit-v1",
        "evidence": "development_diagnostic", "promotion_eligible": False,
        "new_holdout_labels_read": False,
        "comparisons": [{"path": str(path.resolve()), "sha256": file_sha256(path),
                         "result_id": document["result_id"], "input_id": freeze["input_id"]}
                        for path, (document, freeze, _) in zip(args.comparison_manifest, loaded)],
        "script": {"path": str(Path(__file__).resolve()), "sha256": file_sha256(Path(__file__))},
        "outputs": {key: {"path": str(path.resolve()), "sha256": file_sha256(path)}
                    for key, path in outputs.items()},
    }
    manifest["audit_id"] = canonical_id(manifest)
    (args.outdir / "audit_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(pd.concat(all_decisions)[["cell", "tf", "interpretation"]].to_string(index=False))
    print(args.outdir / "audit_manifest.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
