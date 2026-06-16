#!/usr/bin/env python3
"""Compare significant motifs between two diff-footprints result tables."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


KEY_COLUMNS = ["output_prefix", "name", "motif_id"]


def _read_results(path: Path, comparison: str, method: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    sig_col = f"{comparison}_significant_fdr05"
    required = KEY_COLUMNS + [
        f"{comparison}_change",
        f"{comparison}_pvalue",
        f"{comparison}_qvalue_bh",
        f"{comparison}_mean_delta_fp",
        f"{comparison}_mean_log2fc",
        sig_col,
    ]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")

    keep = list(required)
    for col in df.columns:
        if col.endswith("_bound") or col.endswith("_mean_score") or col.endswith("_score_sd"):
            keep.append(col)
    keep = list(dict.fromkeys(keep))
    out = df[keep].copy()
    out[sig_col] = out[sig_col].map(_as_bool)
    out[f"{comparison}_direction"] = out[f"{comparison}_change"].map(
        lambda value: "method_group1_up" if pd.to_numeric(value, errors="coerce") > 0 else "method_group2_up"
    )

    rename = {}
    for col in out.columns:
        if col not in KEY_COLUMNS:
            rename[col] = f"{method}_{col}"
    return out.rename(columns=rename)


def _as_bool(value) -> bool:
    if isinstance(value, bool):
        return value
    if pd.isna(value):
        return False
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def compare_results(method_a_results: Path, method_b_results: Path, method_a_name: str, method_b_name: str, comparison: str) -> pd.DataFrame:
    a = _read_results(method_a_results, comparison, method_a_name)
    b = _read_results(method_b_results, comparison, method_b_name)
    merged = a.merge(b, on=KEY_COLUMNS, how="outer", validate="one_to_one")

    a_sig = f"{method_a_name}_{comparison}_significant_fdr05"
    b_sig = f"{method_b_name}_{comparison}_significant_fdr05"
    a_change = f"{method_a_name}_{comparison}_change"
    b_change = f"{method_b_name}_{comparison}_change"
    a_direction = f"{method_a_name}_{comparison}_direction"
    b_direction = f"{method_b_name}_{comparison}_direction"

    merged[a_sig] = merged[a_sig].fillna(False).map(_as_bool)
    merged[b_sig] = merged[b_sig].fillna(False).map(_as_bool)
    merged["significance_class"] = "neither"
    merged.loc[merged[a_sig] & merged[b_sig], "significance_class"] = "shared"
    merged.loc[merged[a_sig] & ~merged[b_sig], "significance_class"] = f"{method_a_name}_only"
    merged.loc[~merged[a_sig] & merged[b_sig], "significance_class"] = f"{method_b_name}_only"

    a_numeric = pd.to_numeric(merged[a_change], errors="coerce")
    b_numeric = pd.to_numeric(merged[b_change], errors="coerce")
    merged["direction_flip"] = (merged[a_sig] & merged[b_sig] & (a_numeric * b_numeric < 0)).fillna(False)
    merged["manual_review_priority"] = merged["significance_class"].map(
        {
            f"{method_a_name}_only": 1,
            f"{method_b_name}_only": 1,
            "shared": 2,
            "neither": 3,
        }
    )
    merged = merged.sort_values(
        ["manual_review_priority", "direction_flip", a_direction, b_direction, "name", "motif_id"],
        ascending=[True, False, True, True, True, True],
        na_position="last",
    )
    return merged


def write_outputs(df: pd.DataFrame, outdir: Path, method_a_name: str, method_b_name: str, comparison: str) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    a_sig = f"{method_a_name}_{comparison}_significant_fdr05"
    b_sig = f"{method_b_name}_{comparison}_significant_fdr05"
    a_change = f"{method_a_name}_{comparison}_change"
    b_change = f"{method_b_name}_{comparison}_change"

    df.to_csv(outdir / f"{method_a_name}_vs_{method_b_name}_all_motifs.csv", index=False)
    df.loc[df[a_sig] & ~df[b_sig]].to_csv(outdir / f"{method_a_name}_only_significant.csv", index=False)
    df.loc[~df[a_sig] & df[b_sig]].to_csv(outdir / f"{method_b_name}_only_significant.csv", index=False)
    df.loc[df[a_sig] & df[b_sig]].to_csv(outdir / "shared_significant.csv", index=False)

    summary = pd.DataFrame(
        [
            {"metric": "total_motifs", "count": len(df)},
            {"metric": f"{method_a_name}_significant", "count": int(df[a_sig].sum())},
            {"metric": f"{method_b_name}_significant", "count": int(df[b_sig].sum())},
            {"metric": f"{method_a_name}_only", "count": int((df[a_sig] & ~df[b_sig]).sum())},
            {"metric": f"{method_b_name}_only", "count": int((~df[a_sig] & df[b_sig]).sum())},
            {"metric": "shared", "count": int((df[a_sig] & df[b_sig]).sum())},
            {"metric": "neither", "count": int((~df[a_sig] & ~df[b_sig]).sum())},
            {"metric": "direction_flip_shared", "count": int(df["direction_flip"].sum())},
            {"metric": f"{method_a_name}_positive_change", "count": int((pd.to_numeric(df[a_change], errors="coerce") > 0).sum())},
            {"metric": f"{method_a_name}_negative_change", "count": int((pd.to_numeric(df[a_change], errors="coerce") < 0).sum())},
            {"metric": f"{method_b_name}_positive_change", "count": int((pd.to_numeric(df[b_change], errors="coerce") > 0).sum())},
            {"metric": f"{method_b_name}_negative_change", "count": int((pd.to_numeric(df[b_change], errors="coerce") < 0).sum())},
        ]
    )
    summary.to_csv(outdir / "significance_summary.csv", index=False)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--method-a-name", required=True, help="Name for the first method, for example corrected_q95.")
    parser.add_argument("--method-a-results", type=Path, required=True, help="diff-footprints results table for the first method.")
    parser.add_argument("--method-b-name", required=True, help="Name for the second method, for example none.")
    parser.add_argument("--method-b-results", type=Path, required=True, help="diff-footprints results table for the second method.")
    parser.add_argument("--comparison", required=True, help="Comparison prefix in result columns, for example K562_HepG2.")
    parser.add_argument("--outdir", type=Path, required=True, help="Output directory for comparison CSV files.")
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    df = compare_results(args.method_a_results, args.method_b_results, args.method_a_name, args.method_b_name, args.comparison)
    write_outputs(df, args.outdir, args.method_a_name, args.method_b_name, args.comparison)


if __name__ == "__main__":
    main()
