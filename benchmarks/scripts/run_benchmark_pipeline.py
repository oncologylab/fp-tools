#!/usr/bin/env python3
"""Run the fp-tools benchmark summary pipeline from labeled prediction TSVs."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PAPER_SCRIPTS = SCRIPT_DIR.parent.parent / "paper" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))
sys.path.insert(0, str(PAPER_SCRIPTS))

from compute_binary_metrics import bootstrap_confidence_intervals, compute_metrics  # noqa: E402
from compute_calibration import compute_calibration  # noqa: E402
from plot_benchmark_panels import plot_metrics  # noqa: E402
from plot_calibration_panels import plot_calibration  # noqa: E402


def _read_prediction_tables(paths: list[str | Path]) -> pd.DataFrame:
    frames = []
    for path in paths:
        frame = pd.read_csv(path, sep="\t")
        frame["source_table"] = str(path)
        frames.append(frame)
    if not frames:
        raise ValueError("At least one prediction table is required.")
    return pd.concat(frames, ignore_index=True)


def write_run_summary(outdir: Path, outputs: dict[str, Path | list[Path]]) -> Path:
    summary = outdir / "benchmark_run_summary.md"
    with summary.open("w", encoding="utf-8") as handle:
        handle.write("# fp-tools benchmark run summary\n\n")
        for label, value in outputs.items():
            if isinstance(value, list):
                handle.write(f"- {label}:\n")
                for item in value:
                    handle.write(f"  - `{item}`\n")
            else:
                handle.write(f"- {label}: `{value}`\n")
    return summary


def run_benchmark_pipeline(
    predictions: list[str | Path],
    outdir: str | Path,
    label_col: str = "label",
    score_col: str = "score",
    group_cols: list[str] | None = None,
    bins: int = 10,
    bootstrap: int = 0,
    seed: int = 2026,
    title: str = "fp-tools public benchmark",
) -> dict[str, Path | list[Path]]:
    """Run metrics, calibration, and figure generation from labeled predictions."""

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    figures = outdir / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    group_cols = group_cols or ["tf", "cell", "method"]

    combined = _read_prediction_tables(predictions)
    group_cols = [column for column in group_cols if column in combined.columns]

    combined_path = outdir / "combined_predictions.tsv"
    metrics_path = outdir / "binary_metrics.tsv"
    calibration_bins_path = outdir / "calibration_bins.tsv"
    calibration_summary_path = outdir / "calibration_summary.tsv"
    combined.to_csv(combined_path, sep="\t", index=False)

    metrics = compute_metrics(combined, label_col, score_col, group_cols)
    metrics.to_csv(metrics_path, sep="\t", index=False)
    outputs: dict[str, Path | list[Path]] = {
        "combined_predictions": combined_path,
        "binary_metrics": metrics_path,
    }

    if bootstrap > 0:
        bootstrap_path = outdir / "binary_metrics_bootstrap.tsv"
        ci = bootstrap_confidence_intervals(
            combined,
            label_col,
            score_col,
            group_cols,
            n_bootstrap=bootstrap,
            seed=seed,
        )
        ci.to_csv(bootstrap_path, sep="\t", index=False)
        outputs["binary_metrics_bootstrap"] = bootstrap_path

    bins_df, calibration_summary = compute_calibration(combined, label_col, score_col, group_cols, bins=bins)
    bins_df.to_csv(calibration_bins_path, sep="\t", index=False)
    calibration_summary.to_csv(calibration_summary_path, sep="\t", index=False)
    outputs["calibration_bins"] = calibration_bins_path
    outputs["calibration_summary"] = calibration_summary_path

    outputs["benchmark_figures"] = plot_metrics(metrics, figures / "benchmark_summary", title=title)
    outputs["calibration_figures"] = plot_calibration(
        bins_df,
        calibration_summary,
        figures / "calibration_summary",
        title=f"{title} calibration",
    )
    outputs["run_summary"] = write_run_summary(outdir, outputs)
    return outputs


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", nargs="+", required=True, help="One or more labeled prediction TSVs.")
    parser.add_argument("--outdir", required=True, help="Output benchmark result directory.")
    parser.add_argument("--label-col", default="label")
    parser.add_argument("--score-col", default="score")
    parser.add_argument("--group-cols", nargs="*", default=["tf", "cell", "method"])
    parser.add_argument("--bins", type=int, default=10)
    parser.add_argument("--bootstrap", type=int, default=0, help="Optional number of bootstrap resamples.")
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--title", default="fp-tools public benchmark")
    args = parser.parse_args()

    outputs = run_benchmark_pipeline(
        args.predictions,
        args.outdir,
        label_col=args.label_col,
        score_col=args.score_col,
        group_cols=args.group_cols,
        bins=args.bins,
        bootstrap=args.bootstrap,
        seed=args.seed,
        title=args.title,
    )
    for label, value in outputs.items():
        if isinstance(value, list):
            for item in value:
                print(f"{label}\t{item}")
        else:
            print(f"{label}\t{value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
