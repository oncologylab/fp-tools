#!/usr/bin/env python3
"""Score fixed ChIP-labeled motif sites from one or more signal bigWigs."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

from fp_tools.utils import bigwig as pyBigWig


SITE_COLUMNS = {
    "cell", "tf", "TFBS_chr", "TFBS_start", "TFBS_end", "chip_label",
}
SIGNAL_COLUMNS = {"cell", "method", "signal"}


def read_sites(paths: list[Path]) -> pd.DataFrame:
    frame = pd.concat([pd.read_csv(path, sep="\t") for path in paths], ignore_index=True)
    missing = sorted(SITE_COLUMNS.difference(frame.columns))
    if missing:
        raise ValueError("site table is missing columns: " + ", ".join(missing))
    key = ["cell", "tf", "TFBS_chr", "TFBS_start", "TFBS_end"]
    frame = frame.drop_duplicates(key, keep="first").copy()
    frame["chip_label"] = pd.to_numeric(frame["chip_label"], errors="raise").astype(int)
    return frame


def read_signals(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t")
    missing = sorted(SIGNAL_COLUMNS.difference(frame.columns))
    if missing:
        raise ValueError("signal manifest is missing columns: " + ", ".join(missing))
    if frame.duplicated(["cell", "method"]).any():
        raise ValueError("signal manifest has duplicate cell/method rows")
    absent = [str(path) for path in frame["signal"] if not Path(path).is_file()]
    if absent:
        raise FileNotFoundError("signal files do not exist: " + ", ".join(absent[:3]))
    return frame


def score_centers(sites: pd.DataFrame, signal: Path) -> np.ndarray:
    """Read the base-resolution value at each motif center."""

    scores = np.full(len(sites), np.nan, dtype=float)
    handle = pyBigWig.open(str(signal))
    try:
        chromosomes = handle.chroms()
        for output_index, row in enumerate(sites.itertuples(index=False)):
            chrom = str(row.TFBS_chr)
            center = (int(row.TFBS_start) + int(row.TFBS_end)) // 2
            if chrom not in chromosomes or center < 0 or center >= int(chromosomes[chrom]):
                continue
            value = handle.values(chrom, center, center + 1)[0]
            if value is not None and np.isfinite(value):
                scores[output_index] = float(value)
    finally:
        handle.close()
    return scores


def evaluate(sites: pd.DataFrame, signals: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    prediction_frames = []
    metric_rows = []
    for cell, cell_signals in signals.groupby("cell", sort=True):
        subset = sites[sites["cell"] == cell].copy().reset_index(drop=True)
        if subset.empty:
            continue
        method_scores = {
            str(signal_row.method): score_centers(subset, Path(signal_row.signal))
            for signal_row in cell_signals.itertuples(index=False)
        }
        common = np.logical_and.reduce(
            [np.isfinite(scores) for scores in method_scores.values()]
        )
        matched = subset.loc[common].copy().reset_index(drop=True)
        for method, scores in method_scores.items():
            predictions = matched.copy()
            predictions["method"] = method
            predictions["score"] = scores[common]
            prediction_frames.append(predictions)
            for tf, group in predictions.groupby("tf", sort=True):
                labels = group["chip_label"].to_numpy(dtype=int)
                tf_scores = group["score"].to_numpy(dtype=float)
                if len(np.unique(labels)) != 2:
                    continue
                metric_rows.append(
                    {
                        "cell": str(cell),
                        "tf": str(tf),
                        "method": method,
                        "n_sites": int(len(group)),
                        "positive_sites": int(labels.sum()),
                        "auroc": float(roc_auc_score(labels, tf_scores)),
                        "auprc": float(average_precision_score(labels, tf_scores)),
                    }
                )
    predictions = pd.concat(prediction_frames, ignore_index=True) if prediction_frames else pd.DataFrame()
    return predictions, pd.DataFrame(metric_rows)


def paired_chromosome_bootstrap(
    predictions: pd.DataFrame,
    baseline_method: str,
    n_bootstrap: int,
    seed: int,
) -> pd.DataFrame:
    """Bootstrap paired method-minus-baseline metrics by chromosome."""

    if n_bootstrap <= 0 or predictions.empty:
        return pd.DataFrame()
    index = [
        "cell", "tf", "TFBS_chr", "TFBS_start", "TFBS_end", "chip_label",
    ]
    wide = predictions.pivot(index=index, columns="method", values="score").reset_index()
    if baseline_method not in wide:
        raise ValueError(f"baseline method is absent from predictions: {baseline_method}")
    rng = np.random.default_rng(seed)
    rows = []
    for (cell, tf), group in wide.groupby(["cell", "tf"], sort=True):
        if group[baseline_method].isna().any():
            raise ValueError(f"baseline method {baseline_method} is missing for {cell}/{tf}")
        methods = [
            method
            for method in sorted(predictions["method"].unique())
            if method != baseline_method and method in group and group[method].notna().all()
        ]
        blocks = [
            np.flatnonzero(group["TFBS_chr"].to_numpy() == chromosome)
            for chromosome in sorted(group["TFBS_chr"].unique())
        ]
        labels = group["chip_label"].to_numpy(dtype=int)
        baseline = group[baseline_method].to_numpy(dtype=float)
        for method in methods:
            candidate = group[method].to_numpy(dtype=float)
            distributions = {"auroc": [], "auprc": []}
            for _ in range(n_bootstrap):
                selected = rng.integers(0, len(blocks), size=len(blocks))
                positions = np.concatenate([blocks[block_index] for block_index in selected])
                sample_labels = labels[positions]
                if len(np.unique(sample_labels)) != 2:
                    continue
                sample_baseline = baseline[positions]
                sample_candidate = candidate[positions]
                distributions["auroc"].append(
                    roc_auc_score(sample_labels, sample_candidate)
                    - roc_auc_score(sample_labels, sample_baseline)
                )
                distributions["auprc"].append(
                    average_precision_score(sample_labels, sample_candidate)
                    - average_precision_score(sample_labels, sample_baseline)
                )
            for metric, values in distributions.items():
                array = np.asarray(values, dtype=float)
                rows.append(
                    {
                        "cell": str(cell),
                        "tf": str(tf),
                        "baseline_method": baseline_method,
                        "method": method,
                        "metric": metric,
                        "ci_low": float(np.quantile(array, 0.025)) if len(array) else np.nan,
                        "ci_high": float(np.quantile(array, 0.975)) if len(array) else np.nan,
                        "probability_delta_gt_zero": float(np.mean(array > 0)) if len(array) else np.nan,
                        "successful_bootstraps": int(len(array)),
                        "requested_bootstraps": int(n_bootstrap),
                        "resampling_unit": "chromosome",
                        "n_blocks": int(len(blocks)),
                    }
                )
    return pd.DataFrame(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sites", nargs="+", type=Path, required=True)
    parser.add_argument("--signals", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--chromosomes", nargs="*", default=[])
    parser.add_argument("--write-predictions", action="store_true")
    parser.add_argument("--baseline-method", default="raw")
    parser.add_argument("--bootstrap", type=int, default=0)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args(argv)

    sites = read_sites(args.sites)
    if args.chromosomes:
        sites = sites[sites["TFBS_chr"].isin(args.chromosomes)].copy()
    predictions, metrics = evaluate(sites, read_signals(args.signals))
    args.outdir.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(args.outdir / "bigwig_site_metrics.tsv", sep="\t", index=False)
    bootstrap = paired_chromosome_bootstrap(
        predictions,
        baseline_method=args.baseline_method,
        n_bootstrap=args.bootstrap,
        seed=args.seed,
    )
    if not bootstrap.empty:
        bootstrap.to_csv(args.outdir / "bigwig_site_bootstrap.tsv", sep="\t", index=False)
    if args.write_predictions:
        predictions.to_csv(
            args.outdir / "bigwig_site_predictions.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )
    print(metrics.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
