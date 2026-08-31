#!/usr/bin/env python3
"""Decompose strand-footprint gains into bias-model and cut-shift effects.

The four prespecified arms form a balanced 2 x 2 experiment.  Detector
settings are selected by their mean score across all arms, then held common
while model, shift, and interaction effects are calculated for every TF.
"""

from __future__ import annotations

import argparse
from hashlib import blake2b, sha256
import json
from pathlib import Path

import numpy as np
import pandas as pd


ARMS = (
    "MT_SELMA10_4m4",
    "MT_SELMA10_4m5",
    "MT_LOG81_4m4",
    "MT_LOG81_4m5",
)
GROUP_KEYS = ("cell", "tf", "motif_family", "training_scope")
SETTING_KEYS = ("channel_set", "smoother", "window_limit")
METRICS = ("auroc", "auprc", "brier")


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_seed(*values: object, seed: int = 2026) -> int:
    digest = blake2b(digest_size=8)
    for value in (seed, *values):
        digest.update(str(value).encode())
        digest.update(b"\0")
    return int.from_bytes(digest.digest(), "little") % (2**32 - 1)


def select_common_settings(metrics: pd.DataFrame) -> pd.DataFrame:
    """Select one detector setting without favoring any factorial arm."""

    required = set(GROUP_KEYS + SETTING_KEYS + ("bias_configuration", "status", "selection_score"))
    missing = sorted(required.difference(metrics.columns))
    if missing:
        raise ValueError(f"factorial metrics lack columns: {missing}")
    passing = metrics[
        metrics["status"].eq("ok") & metrics["bias_configuration"].isin(ARMS)
    ].copy()
    duplicate_keys = list(GROUP_KEYS + SETTING_KEYS) + ["bias_configuration"]
    if passing.duplicated(duplicate_keys).any():
        raise ValueError("factorial metrics contain duplicate arm/setting rows")
    wide = passing.pivot(
        index=list(GROUP_KEYS + SETTING_KEYS),
        columns="bias_configuration",
        values="selection_score",
    )
    complete = wide.dropna(subset=list(ARMS)).reset_index()
    if complete.empty:
        return complete
    complete["mean_arm_selection_score"] = complete[list(ARMS)].mean(axis=1)
    complete["minimum_arm_selection_score"] = complete[list(ARMS)].min(axis=1)
    return (
        complete.sort_values(
            list(GROUP_KEYS)
            + ["mean_arm_selection_score", "minimum_arm_selection_score"]
            + list(SETTING_KEYS),
            ascending=[True] * len(GROUP_KEYS) + [False, False, True, True, True],
            kind="mergesort",
        )
        .groupby(list(GROUP_KEYS), as_index=False, sort=True)
        .head(1)
        .reset_index(drop=True)
    )


def factorial_effects(metrics: pd.DataFrame, settings: pd.DataFrame) -> pd.DataFrame:
    """Calculate Shapley-balanced model and shift contributions."""

    if settings.empty:
        return pd.DataFrame()
    selected = metrics.merge(
        settings[list(GROUP_KEYS + SETTING_KEYS)],
        on=list(GROUP_KEYS + SETTING_KEYS),
        how="inner",
        validate="many_to_one",
    )
    selected = selected[
        selected["status"].eq("ok") & selected["bias_configuration"].isin(ARMS)
    ].copy()
    rows = []
    for keys, group in selected.groupby(list(GROUP_KEYS + SETTING_KEYS), sort=True):
        identity = dict(zip(GROUP_KEYS + SETTING_KEYS, keys))
        indexed = group.set_index("bias_configuration")
        if not set(ARMS).issubset(indexed.index):
            continue
        for metric in METRICS:
            values = {arm: float(indexed.loc[arm, metric]) for arm in ARMS}
            if not all(np.isfinite(list(values.values()))):
                continue
            direction = -1.0 if metric == "brier" else 1.0
            oriented = {arm: direction * value for arm, value in values.items()}
            shift_contribution = 0.5 * (
                oriented["MT_SELMA10_4m4"]
                - oriented["MT_SELMA10_4m5"]
                + oriented["MT_LOG81_4m4"]
                - oriented["MT_LOG81_4m5"]
            )
            model_contribution = 0.5 * (
                oriented["MT_SELMA10_4m4"]
                - oriented["MT_LOG81_4m4"]
                + oriented["MT_SELMA10_4m5"]
                - oriented["MT_LOG81_4m5"]
            )
            interaction = (
                oriented["MT_SELMA10_4m4"]
                - oriented["MT_SELMA10_4m5"]
                - oriented["MT_LOG81_4m4"]
                + oriented["MT_LOG81_4m5"]
            )
            total = (
                oriented["MT_SELMA10_4m4"] - oriented["MT_LOG81_4m5"]
            )
            if not np.isclose(total, shift_contribution + model_contribution, atol=1e-12):
                raise AssertionError("factorial contributions do not sum to diagonal contrast")
            threshold = 0.01
            if total < threshold:
                attribution = "no_material_gain"
            elif abs(interaction) >= max(0.02, 0.5 * abs(total)):
                attribution = "interaction_sensitive"
            elif shift_contribution > 0 and model_contribution > 0:
                ratio = shift_contribution / max(model_contribution, np.finfo(float).eps)
                if ratio >= 2:
                    attribution = "shift_dominant"
                elif ratio <= 0.5:
                    attribution = "bias_model_dominant"
                else:
                    attribution = "joint"
            elif shift_contribution > model_contribution:
                attribution = "shift_dominant"
            else:
                attribution = "bias_model_dominant"
            best_arm = max(oriented, key=lambda arm: (oriented[arm], arm))
            rows.append(
                {
                    **identity,
                    "metric": metric,
                    "reference_value_log81_4m5": values["MT_LOG81_4m5"],
                    "candidate_value_selma10_4m4": values["MT_SELMA10_4m4"],
                    "candidate_minus_reference_raw": (
                        values["MT_SELMA10_4m4"] - values["MT_LOG81_4m5"]
                    ),
                    "total_improvement": total,
                    "shift_contribution": shift_contribution,
                    "bias_model_contribution": model_contribution,
                    "model_by_shift_interaction": interaction,
                    "best_arm": best_arm,
                    "attribution": attribution,
                    **{f"value_{arm.lower()}": values[arm] for arm in ARMS},
                }
            )
    return pd.DataFrame(rows)


def bootstrap_mean_interval(
    values: np.ndarray,
    *,
    iterations: int,
    seed: int,
) -> tuple[float, float, float]:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if len(finite) == 0:
        return np.nan, np.nan, np.nan
    mean = float(np.mean(finite))
    if len(finite) < 2 or iterations < 2:
        return mean, mean, mean
    rng = np.random.default_rng(seed)
    samples = finite[rng.integers(0, len(finite), size=(iterations, len(finite)))]
    lower, upper = np.quantile(np.mean(samples, axis=1), [0.025, 0.975])
    return mean, float(lower), float(upper)


def summarize_effects(effects: pd.DataFrame, *, bootstraps: int, seed: int) -> pd.DataFrame:
    rows = []
    components = (
        "total_improvement",
        "shift_contribution",
        "bias_model_contribution",
        "model_by_shift_interaction",
    )
    for (scope, metric), group in effects.groupby(["training_scope", "metric"], sort=True):
        for component in components:
            mean, lower, upper = bootstrap_mean_interval(
                group[component].to_numpy(dtype=float),
                iterations=bootstraps,
                seed=stable_seed(scope, metric, component, seed=seed),
            )
            rows.append(
                {
                    "training_scope": scope,
                    "metric": metric,
                    "component": component,
                    "cell_tf_pairs": int(len(group)),
                    "mean": mean,
                    "bootstrap_lower": lower,
                    "bootstrap_upper": upper,
                    "median": float(group[component].median()),
                    "positive_fraction": float(np.mean(group[component] > 0)),
                }
            )
    return pd.DataFrame(rows)


def render_heatmaps(effects: pd.DataFrame, output: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    output.parent.mkdir(parents=True, exist_ok=True)
    components = ["shift_contribution", "bias_model_contribution", "total_improvement"]
    labels = ["Cut shift", "Bias model", "Total: SELMA +4/−4 vs log81 +4/−5"]
    with PdfPages(output) as pdf:
        for metric in ("auroc", "auprc"):
            for scope in sorted(effects["training_scope"].unique()):
                subset = effects[
                    effects["metric"].eq(metric) & effects["training_scope"].eq(scope)
                ].copy()
                if subset.empty:
                    continue
                subset["cell_tf"] = subset["cell"].astype(str) + " — " + subset["tf"].astype(str)
                subset = subset.sort_values("total_improvement", ascending=False)
                matrix = subset[components].to_numpy(dtype=float)
                limit = max(float(np.nanquantile(np.abs(matrix), 0.98)), 0.01)
                height = max(4.5, 0.32 * len(subset) + 1.8)
                figure, axis = plt.subplots(figsize=(9.0, height))
                image = axis.imshow(matrix, aspect="auto", cmap="RdBu_r", vmin=-limit, vmax=limit)
                axis.set_xticks(np.arange(len(labels)), labels=labels)
                axis.set_yticks(np.arange(len(subset)), labels=subset["cell_tf"].tolist())
                axis.set_title(f"{metric.upper()} factorial decomposition — {scope}")
                for row in range(len(subset)):
                    for column in range(len(components)):
                        value = matrix[row, column]
                        axis.text(column, row, f"{value:+.3f}", ha="center", va="center", fontsize=7)
                figure.colorbar(image, ax=axis, label="Improvement (positive is better)")
                figure.tight_layout()
                pdf.savefig(figure)
                plt.close(figure)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metrics", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--bootstraps", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args(argv)
    if args.bootstraps < 2:
        raise SystemExit("--bootstraps must be at least two")
    metrics = pd.read_csv(args.metrics, sep="\t")
    observed_arms = set(metrics["bias_configuration"].dropna().astype(str))
    missing_arms = sorted(set(ARMS).difference(observed_arms))
    if missing_arms:
        raise SystemExit(f"factorial metrics are missing arms: {missing_arms}")
    settings = select_common_settings(metrics)
    effects = factorial_effects(metrics, settings)
    summary = summarize_effects(effects, bootstraps=args.bootstraps, seed=args.seed)
    args.outdir.mkdir(parents=True, exist_ok=True)
    settings_path = args.outdir / "strand_factorial_common_settings.tsv"
    effects_path = args.outdir / "strand_factorial_effects.tsv"
    summary_path = args.outdir / "strand_factorial_summary.tsv"
    figure_path = args.outdir / "strand_factorial_heatmaps.pdf"
    settings.to_csv(settings_path, sep="\t", index=False)
    effects.to_csv(effects_path, sep="\t", index=False)
    summary.to_csv(summary_path, sep="\t", index=False)
    render_heatmaps(effects, figure_path)
    manifest = {
        "schema": "fp-tools-strand-bias-factorial-v1",
        "locked_test_labels_read": False,
        "training_labels_used_by_input_detector": True,
        "selection_rule": "maximize mean selection score across all four arms at a common detector setting",
        "arms": list(ARMS),
        "metrics": str(args.metrics),
        "metrics_sha256": file_sha256(args.metrics),
        "bootstraps": args.bootstraps,
        "seed": args.seed,
        "outputs": {
            path.name: {"path": str(path), "sha256": file_sha256(path)}
            for path in (settings_path, effects_path, summary_path, figure_path)
        },
    }
    (args.outdir / "strand_factorial_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
