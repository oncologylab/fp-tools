#!/usr/bin/env python3
"""Compare fixed TF candidates across independent read-sampling seeds."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


KEY = ["cell", "tf", "correction"]


def compare_seeds(first: pd.DataFrame, second: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    columns = KEY + ["n_sites", "positive_sites", "auroc", "auprc"]
    merged = first[columns].merge(
        second[columns], on=KEY, suffixes=("_seed_a", "_seed_b"), validate="one_to_one"
    )
    merged["delta_auroc_seed_b_minus_a"] = merged["auroc_seed_b"] - merged["auroc_seed_a"]
    merged["delta_auprc_seed_b_minus_a"] = merged["auprc_seed_b"] - merged["auprc_seed_a"]
    summaries = []
    for (cell, tf), group in merged.groupby(["cell", "tf"], sort=True):
        winner_a_auroc = group.loc[group["auroc_seed_a"].idxmax(), "correction"]
        winner_b_auroc = group.loc[group["auroc_seed_b"].idxmax(), "correction"]
        winner_a_auprc = group.loc[group["auprc_seed_a"].idxmax(), "correction"]
        winner_b_auprc = group.loc[group["auprc_seed_b"].idxmax(), "correction"]
        summaries.append(
            {
                "cell": cell,
                "tf": tf,
                "corrections_compared": len(group),
                "winner_auroc_seed_a": winner_a_auroc,
                "winner_auroc_seed_b": winner_b_auroc,
                "stable_auroc_winner": winner_a_auroc == winner_b_auroc,
                "winner_auprc_seed_a": winner_a_auprc,
                "winner_auprc_seed_b": winner_b_auprc,
                "stable_auprc_winner": winner_a_auprc == winner_b_auprc,
                "max_abs_seed_delta_auroc": group["delta_auroc_seed_b_minus_a"].abs().max(),
                "max_abs_seed_delta_auprc": group["delta_auprc_seed_b_minus_a"].abs().max(),
            }
        )
    return merged, pd.DataFrame(summaries)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed-a", type=Path, required=True)
    parser.add_argument("--seed-b", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args(argv)
    rows, summary = compare_seeds(
        pd.read_csv(args.seed_a, sep="\t"), pd.read_csv(args.seed_b, sep="\t")
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    rows.to_csv(args.output, sep="\t", index=False)
    summary.to_csv(args.summary, sep="\t", index=False)
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
