#!/usr/bin/env python
"""Replicate-aware uncertainty summaries for BINDetect result tables."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm


def replicate_uncertainty(change: float, pvalue: float, n_min: int) -> dict[str, float]:
    """Derive a p-value-based standard error and replicate-shrunk effect estimate.

    BINDetect result tables report an aggregate change and a two-sided p-value but
    not per-replicate scores, so the standard error is recovered from the p-value
    (``SE = |change| / z``) and the effect is shrunk toward zero by a weight that
    grows with replicate support (``n / (n + 1)``), giving a conservative
    replicate-aware estimate.
    """

    weight = n_min / (n_min + 1.0) if n_min > 0 else 0.0
    if pd.isna(change) or pd.isna(pvalue):
        return {"effect_se": np.nan, "z_score": np.nan, "shrinkage_weight": weight,
                "shrunk_change": np.nan, "ci_lower": np.nan, "ci_upper": np.nan}
    clipped = min(max(float(pvalue), 1e-308), 1.0)
    z = float(norm.isf(clipped / 2.0))
    se = abs(float(change)) / z if z > 0 else np.nan
    shrunk = float(change) * weight
    if pd.notna(se):
        ci_lower, ci_upper = float(change) - 1.96 * se, float(change) + 1.96 * se
    else:
        ci_lower = ci_upper = np.nan
    return {"effect_se": se, "z_score": z, "shrinkage_weight": weight,
            "shrunk_change": shrunk, "ci_lower": ci_lower, "ci_upper": ci_upper}


def read_replicate_map(path: str | Path | None) -> dict[str, int]:
    """Read optional condition-to-replicate mapping/count TSV."""

    if path is None:
        return {}
    mapping: dict[str, set[str]] = {}
    with Path(path).open(encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        lower = [item.lower() for item in header]
        has_header = "condition" in lower
        if has_header:
            cond_idx = lower.index("condition")
            rep_idx = lower.index("replicate") if "replicate" in lower else None
            count_idx = lower.index("n_replicates") if "n_replicates" in lower else None
        else:
            cond_idx, rep_idx, count_idx = 0, 1, None
            fields = header
            if len(fields) >= 2 and fields[0]:
                mapping.setdefault(fields[0], set()).add(fields[1])
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= cond_idx:
                continue
            condition = fields[cond_idx]
            if count_idx is not None and len(fields) > count_idx:
                try:
                    mapping[condition] = {f"replicate_{idx + 1}" for idx in range(int(fields[count_idx]))}
                    continue
                except ValueError:
                    pass
            replicate = fields[rep_idx] if rep_idx is not None and len(fields) > rep_idx else f"replicate_{len(mapping.get(condition, set())) + 1}"
            mapping.setdefault(condition, set()).add(replicate)
    return {condition: len(replicates) for condition, replicates in mapping.items()}


def infer_conditions(frame: pd.DataFrame) -> list[str]:
    """Infer BINDetect condition names from `<condition>_mean_score` columns."""

    suffix = "_mean_score"
    return [column[:-len(suffix)] for column in frame.columns if column.endswith(suffix)]


def infer_comparisons(frame: pd.DataFrame, conditions: list[str]) -> list[tuple[str, str, str]]:
    """Infer comparison base names and condition pairs from change/pvalue columns."""

    comparisons: list[tuple[str, str, str]] = []
    for column in frame.columns:
        if not column.endswith("_change"):
            continue
        base = column[:-len("_change")]
        if f"{base}_pvalue" not in frame.columns:
            continue
        matched = None
        for cond1 in conditions:
            for cond2 in conditions:
                if cond1 != cond2 and base == f"{cond1}_{cond2}":
                    matched = (base, cond1, cond2)
                    break
            if matched is not None:
                break
        if matched is None and "_" in base:
            cond1, cond2 = base.split("_", 1)
            matched = (base, cond1, cond2)
        if matched is not None:
            comparisons.append(matched)
    return comparisons


def condition_replicate_count(frame: pd.DataFrame, condition: str, replicate_counts: dict[str, int]) -> int:
    """Return replicate count from explicit map, result columns, or single-replicate fallback."""

    if condition in replicate_counts:
        return int(replicate_counts[condition])
    for suffix in ("_n_replicates", "_replicates"):
        column = f"{condition}{suffix}"
        if column in frame.columns:
            values = pd.to_numeric(frame[column], errors="coerce").dropna()
            if len(values):
                return int(max(1, values.iloc[0]))
    return 1


def build_replicate_report(
    results: str | Path,
    output: str | Path,
    summary_output: str | Path | None = None,
    figure_output: str | Path | None = None,
    replicate_map: str | Path | None = None,
    alpha: float = 0.05,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Create long-form replicate-aware comparison diagnostics from BINDetect results."""

    frame = pd.read_csv(results, sep="\t")
    conditions = infer_conditions(frame)
    comparisons = infer_comparisons(frame, conditions)
    if not comparisons:
        raise ValueError("No BINDetect comparison columns of the form <cond1>_<cond2>_change/pvalue were found.")
    replicate_counts = read_replicate_map(replicate_map)

    rows = []
    for _, row in frame.iterrows():
        for base, cond1, cond2 in comparisons:
            change = pd.to_numeric(pd.Series([row[f"{base}_change"]]), errors="coerce").iloc[0]
            pvalue = pd.to_numeric(pd.Series([row[f"{base}_pvalue"]]), errors="coerce").iloc[0]
            n1 = condition_replicate_count(frame, cond1, replicate_counts)
            n2 = condition_replicate_count(frame, cond2, replicate_counts)
            support = "replicate-supported" if min(n1, n2) >= 2 else "single-replicate"
            abs_change = abs(float(change)) if pd.notna(change) else np.nan
            neg_log10_pvalue = -np.log10(max(float(pvalue), 1e-308)) if pd.notna(pvalue) else np.nan
            qvalue_col = f"{base}_qvalue_bh"
            fdr_col = f"{base}_significant_fdr05"
            qvalue = pd.to_numeric(pd.Series([row[qvalue_col]]), errors="coerce").iloc[0] if qvalue_col in frame.columns else np.nan
            fdr_significant = bool(str(row[fdr_col]).strip().lower() in {"true", "1", "yes"}) if fdr_col in frame.columns else bool(pd.notna(qvalue) and float(qvalue) <= alpha)
            uncertainty = replicate_uncertainty(change, pvalue, min(n1, n2))
            rows.append(
                {
                    "name": row.get("name", ""),
                    "motif_id": row.get("motif_id", ""),
                    "cluster": row.get("cluster", ""),
                    "comparison": base,
                    "condition_1": cond1,
                    "condition_2": cond2,
                    "condition_1_replicates": n1,
                    "condition_2_replicates": n2,
                    "replicate_support": support,
                    "condition_1_mean_score": row.get(f"{cond1}_mean_score", np.nan),
                    "condition_2_mean_score": row.get(f"{cond2}_mean_score", np.nan),
                    "condition_1_bound": row.get(f"{cond1}_bound", np.nan),
                    "condition_2_bound": row.get(f"{cond2}_bound", np.nan),
                    "change": change,
                    "abs_change": abs_change,
                    "pvalue": pvalue,
                    "qvalue_bh": qvalue,
                    "neg_log10_pvalue": neg_log10_pvalue,
                    "significant": bool(pd.notna(pvalue) and float(pvalue) <= alpha),
                    "significant_fdr05": fdr_significant,
                    "direction": cond1 if pd.notna(change) and float(change) > 0 else cond2 if pd.notna(change) and float(change) < 0 else "none",
                    "effect_se": uncertainty["effect_se"],
                    "z_score": uncertainty["z_score"],
                    "shrinkage_weight": uncertainty["shrinkage_weight"],
                    "shrunk_change": uncertainty["shrunk_change"],
                    "ci_lower": uncertainty["ci_lower"],
                    "ci_upper": uncertainty["ci_upper"],
                }
            )

    report = pd.DataFrame(rows)
    summary = (
        report.groupby(["comparison", "replicate_support"], dropna=False)
        .agg(
            n_motifs=("name", "count"),
            significant=("significant", "sum"),
            significant_fdr05=("significant_fdr05", "sum"),
            median_abs_change=("abs_change", "median"),
            max_neg_log10_pvalue=("neg_log10_pvalue", "max"),
            median_effect_se=("effect_se", "median"),
            mean_shrinkage_weight=("shrinkage_weight", "mean"),
        )
        .reset_index()
    )

    out = Path(output)
    out.parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(out, sep="\t", index=False)
    if summary_output is not None:
        summary_out = Path(summary_output)
        summary_out.parent.mkdir(parents=True, exist_ok=True)
        summary.to_csv(summary_out, sep="\t", index=False)
    if figure_output is not None:
        plot_replicate_report(report, summary, figure_output)
    return report, summary


def plot_replicate_report(report: pd.DataFrame, summary: pd.DataFrame, output: str | Path) -> Path:
    """Write a compact PDF/PNG/SVG-ready uncertainty report figure."""

    out = Path(output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.3), constrained_layout=True)
    colors = {"single-replicate": "tab:orange", "replicate-supported": "tab:blue"}

    ax = axes[0]
    for support, group in report.groupby("replicate_support"):
        ax.scatter(group["change"], group["neg_log10_pvalue"], s=18, alpha=0.8, label=support, color=colors.get(support))
    ax.axvline(0, color="0.5", linewidth=0.8, linestyle="--")
    ax.set_xlabel("BINDetect change")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title("Differential binding evidence")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)

    ax = axes[1]
    labels = [f"{row.comparison}\n{row.replicate_support}" for row in summary.itertuples()]
    ax.bar(range(len(summary)), summary["significant"], color=[colors.get(value, "0.5") for value in summary["replicate_support"]])
    ax.set_xticks(range(len(summary)))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel(f"significant motifs")
    ax.set_title("Significant motifs by replicate support")
    ax.grid(axis="y", alpha=0.25)

    fig.savefig(out, dpi=300 if out.suffix.lower() == ".png" else None)
    plt.close(fig)
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", required=True, help="BINDetect *_results.txt table.")
    parser.add_argument("--out", required=True, help="Output long-form replicate uncertainty TSV.")
    parser.add_argument("--summary-out", help="Optional comparison-level summary TSV.")
    parser.add_argument("--figure-out", help="Optional PDF/SVG/PNG report figure.")
    parser.add_argument("--replicate-map", help="Optional TSV with condition/replicate or condition/n_replicates columns.")
    parser.add_argument("--alpha", type=float, default=0.05, help="P-value cutoff for significant motif counts.")
    args = parser.parse_args(argv)

    report, summary = build_replicate_report(
        args.results,
        args.out,
        summary_output=args.summary_out,
        figure_output=args.figure_out,
        replicate_map=args.replicate_map,
        alpha=args.alpha,
    )
    print(f"wrote {len(report)} replicate diagnostic rows to {args.out}")
    if args.summary_out:
        print(f"wrote {len(summary)} replicate summary rows to {args.summary_out}")
    if args.figure_out:
        print(f"wrote replicate diagnostic figure to {args.figure_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
