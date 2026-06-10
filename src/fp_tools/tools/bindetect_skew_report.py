#!/usr/bin/env python3
"""
Skewness and shift report generator for BINDetect summary tables.

This is an fp-tools-specific reporting module. It scans BINDetect result tables
for comparison columns, computes direction/shift/skew statistics for each
comparison, and writes a compact multi-page PDF plus optional JSON summary.
"""

from __future__ import annotations
import argparse
import json
import math
import re
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import FancyBboxPatch
from fp_tools.utils.plotting_style import PDF_FONT_SIZE, apply_pdf_style, ascii_tick_formatter

# ---------- Utilities ----------

def infer_comp_from_filename(path: str | Path) -> str:
    """Infer comparison label from filename (e.g., D_vs_HY_results.txt → D_vs_HY)."""
    stem = Path(path).stem  # e.g., "D_vs_HY_results"
    m = re.search(r"(.*?)(_results?)$", stem, flags=re.IGNORECASE)
    return m.group(1) if m else stem

def hodges_lehmann_onesample(x: np.ndarray) -> float:
    """One-sample Hodges–Lehmann estimator: median of Walsh averages (xi+xj)/2, i<=j."""
    x = np.asarray(x, dtype=float)
    n = x.size
    walsh = []
    for i in range(n):
        walsh.append((x[i] + x[i:]) / 2.0)
    return float(np.median(np.concatenate(walsh, axis=0)))

def fisher_skewness_g1(x: np.ndarray) -> float:
    """Fisher–Pearson standardized third central moment (bias-corrected)."""
    return float(stats.skew(x, bias=False))

def fisher_skewness_se(n: int) -> float:
    """SE of Fisher skewness under normality (valid for n>8)."""
    return math.sqrt(6 * n * (n - 1) / ((n - 2) * (n + 1) * (n + 3)))

def bowley_skewness(x: np.ndarray) -> float:
    """
    Bowley (quartile) skewness: gB = (Q3 + Q1 - 2*Q2) / (Q3 - Q1).
    Robust to outliers; undefined if Q3 == Q1 (return 0).
    """
    q1, q2, q3 = np.quantile(x, [0.25, 0.5, 0.75], method="linear")
    denom = (q3 - q1)
    if denom <= 0:
        return 0.0
    return float((q3 + q1 - 2.0*q2) / denom)

def binomial_sign_test(pos: int, n: int) -> float:
    """Two-sided binomial test for proportion of positives (H0: p=0.5)."""
    try:
        return stats.binomtest(pos, n, p=0.5, alternative="two-sided").pvalue
    except AttributeError:  # SciPy < 1.7
        return stats.binom_test(pos, n, p=0.5, alternative="two-sided")

def sign_randomization_perm_p(x: np.ndarray, stat_fn, n_perm: int, rng: np.random.Generator) -> tuple[float, float]:
    """
    Permutation (sign-randomization) test for symmetry about 0 using statistic stat_fn(x).
    Returns (observed_stat, two-sided p) with small-sample correction (Phipson & Smyth, 2010).
    """
    x = np.asarray(x, dtype=float)
    obs = stat_fn(x)
    n = x.size
    extreme = 0
    for _ in range(int(n_perm)):
        signs = rng.choice([-1.0, 1.0], size=n, replace=True)
        val = stat_fn(signs * x)
        if abs(val) >= abs(obs) - 1e-15:
            extreme += 1
    pval = (extreme + 1.0) / (n_perm + 1.0)
    return float(obs), float(pval)

def autodetect_pairs(df: pd.DataFrame) -> List[Tuple[str, str]]:
    """Find all <prefix>_change + <prefix>_pvalue pairs."""
    change_cols = [c for c in df.columns if c.endswith("_change")]
    pairs = []
    for ch in change_cols:
        prefix = ch[:-len("_change")]
        pv = prefix + "_pvalue"
        if pv in df.columns:
            pairs.append((ch, pv))
    if not pairs:
        raise ValueError("No <prefix>_change + <prefix>_pvalue column pairs found.")
    return pairs

def read_numeric_tsv(path: str | Path) -> pd.DataFrame:
    """Read TSV and coerce numeric-looking columns to numeric (NaNs for non-numeric cells)."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    for col in df.columns:
        ser = pd.to_numeric(df[col], errors="coerce")
        if ser.notna().sum() >= max(5, 0.5 * len(ser)):
            df[col] = ser
    return df

# ---------- Plotting helpers ----------

FIGSIZE = (11, 4.5)  # uniform across all pages
TITLE_FONTSIZE = PDF_FONT_SIZE
AX_FONTSIZE = PDF_FONT_SIZE
TABLE_FONT = PDF_FONT_SIZE
apply_pdf_style()


def _apply_ascii_minus(ax, x_decimals: int | None = None, y_decimals: int | None = None) -> None:
    """Force ASCII hyphen-minus in tick labels for PDF output."""
    ax.xaxis.set_major_formatter(ascii_tick_formatter(x_decimals))
    ax.yaxis.set_major_formatter(ascii_tick_formatter(y_decimals))

def page_hist_density(ax, x: np.ndarray, title: str):
    ax.hist(x, bins="auto", density=True, alpha=0.35, label="Changes")
    try:
        from scipy.stats import gaussian_kde
        kde = gaussian_kde(x[~np.isnan(x)])
        xs = np.linspace(np.nanmin(x), np.nanmax(x), 400)
        ax.plot(xs, kde(xs), lw=2, label="KDE")
    except Exception:
        pass
    ax.axvline(0, color="k", ls="--", lw=1)
    ax.set_xlabel("Change", fontsize=AX_FONTSIZE)
    ax.set_ylabel("Density", fontsize=AX_FONTSIZE)
    ax.set_title(title, fontsize=TITLE_FONTSIZE, pad=8)
    ax.legend(loc="best", fontsize=PDF_FONT_SIZE)
    _apply_ascii_minus(ax)

def page_volcano(ax, change: np.ndarray, pval: np.ndarray, highlight: np.ndarray | None, title: str):
    y = -np.log10(np.clip(pval, 1e-308, None))
    if highlight is None:
        ax.scatter(change, y, s=8, alpha=0.6, edgecolor="none", c="tab:gray")
    else:
        colors = np.full(change.shape, "tab:gray", dtype=object)
        neg_hi = (change < 0) & highlight
        pos_hi = (change > 0) & highlight
        colors[neg_hi] = "tab:green"
        colors[pos_hi] = "tab:red"
        ax.scatter(change, y, s=8, alpha=0.65, edgecolor="none", c=colors)
    ax.axvline(0, color="k", ls="--", lw=1)
    ax.set_xlabel("Change", fontsize=AX_FONTSIZE)
    ax.set_ylabel("-log10(pvalue)", fontsize=AX_FONTSIZE)
    ax.set_title(title, fontsize=TITLE_FONTSIZE, pad=8)
    _apply_ascii_minus(ax)

def page_mirrored_ecdf(ax, x: np.ndarray, title: str):
    x = x[np.isfinite(x)]
    pos = np.sort(x[x > 0]); neg = np.sort(np.abs(x[x < 0]))
    if pos.size:
        ax.plot(pos, np.arange(1, pos.size + 1) / pos.size, label="ECDF(+change)")
    if neg.size:
        ax.plot(neg, np.arange(1, neg.size + 1) / neg.size, label="ECDF(|-change|)")
    ax.set_xlabel("Magnitude", fontsize=AX_FONTSIZE)
    ax.set_ylabel("ECDF", fontsize=AX_FONTSIZE)
    ax.set_title(title + " (Mirrored tails)", fontsize=TITLE_FONTSIZE, pad=8)
    ax.legend(loc="best", fontsize=PDF_FONT_SIZE)
    _apply_ascii_minus(ax)

def page_sign_bar(ax, pos: int, neg: int, n: int, prop_pos: float):
    ax.bar(["Positive", "Negative"], [pos, neg])
    ax.set_ylabel("Count (non-zero)", fontsize=AX_FONTSIZE)
    ax.set_title(f"Sign counts (n={n}, prop+={prop_pos:.3f})", fontsize=TITLE_FONTSIZE, pad=8)
    _apply_ascii_minus(ax)

# ---------- Summary page (original language, clean layout) ----------

def page_summary(pdf: PdfPages, key: str, res: dict, params: dict):
    """Append a concise text summary page with interpretation + parameters."""
    import textwrap

    fig = plt.figure(figsize=FIGSIZE)
    fig.suptitle(f"{key} — Summary & Methods", fontsize=TITLE_FONTSIZE, weight="bold", y=0.96)

    # Two columns: left explanation, right results/parameters card
    gs = fig.add_gridspec(
        nrows=1, ncols=2, left=0.05, right=0.98, top=0.90, bottom=0.10, wspace=0.12,
        width_ratios=[0.57, 0.43]
    )
    ax_left  = fig.add_subplot(gs[0, 0]); ax_left.axis("off")
    ax_right = fig.add_subplot(gs[0, 1]); ax_right.axis("off")

    def wrap_bullet(prefix: str, text: str, width: int) -> str:
        return textwrap.fill(text, width=width, initial_indent=prefix, subsequent_indent=" " * len(prefix))

    left_width = 54

    expl_lines = [
        "How to read this page:",
        wrap_bullet("• Direction: ", "Sign test asks whether positive and negative TF changes are balanced.", left_width),
        "",
        wrap_bullet("• Shift: ", "Wilcoxon tests whether the typical change differs from 0. "
                    "HL is the typical effect size on the same scale as the change score.", left_width),
        "",
        wrap_bullet("• Asymmetry: ", f"{params['skew_stat'].capitalize()} skewness asks whether one tail is heavier "
                    "than the other after accounting for sign flips.", left_width),
        "",
        "Quick guide:",
        wrap_bullet("• ", "Positive HL means stronger binding in the first condition; negative HL means stronger "
                    "binding in the second.", left_width),
        wrap_bullet("• ", "A small skewness p-value supports tail imbalance, not just a center shift.", left_width),
        wrap_bullet("• ", "Bowley skewness is robust; Fisher skewness is more sensitive to outliers.", left_width),
    ]
    ax_left.text(0.0, 1.0, "\n".join(expl_lines), ha="left", va="top", fontsize=AX_FONTSIZE, color="#1E1E1E")

    # Right column: results + parameters monospace card
    def kv_block(pairs, keyw=16, width=40):
        lines = []
        for k, v in pairs:
            kfix = (k + ":").ljust(keyw)
            lines.append(textwrap.fill(str(v), width=width, initial_indent=kfix, subsequent_indent=" " * keyw))
        return "\n".join(lines)

    results_pairs = [
        ("n TFs (finite)",      f"{res['n_TFs']}"),
        ("+/- (nonzero)",       f"{res['positives']} / {res['negatives']}  (prop+ = {res['prop_positive']:.3f})"),
        ("Sign test p",         f"{res['sign_test_p']:.3e}"),
        ("Wilcoxon p",          f"{res['wilcoxon_p']:.3e}"),
        ("HL estimate",         f"{res['HL_estimate']:.4f}"),
        ("Skewness value",      f"{res['skewness_value']:.4f}"),
        ("Skewness p",          f"{res['skewness_p']:.3e}"),
    ]
    params_pairs = [
        ("skew_method", params["skew_method"]),
        ("skew_stat",   params["skew_stat"]),
        ("n_perm",      params["n_perm"]),
        ("seed",        params["seed"]),
        ("file",        Path(params["file"]).name),
    ]

    results_txt = kv_block(results_pairs)
    params_txt  = kv_block(params_pairs)

    card = FancyBboxPatch((0.00, 0.00), 1.00, 1.00, transform=ax_right.transAxes,
                          boxstyle="round,pad=0.014,rounding_size=8",
                          linewidth=0.6, edgecolor="#C7CED8", facecolor="#F9FAFB", zorder=0)
    ax_right.add_patch(card)

    ax_right.text(0.04, 0.94, "Results", ha="left", va="top",
                  fontsize=AX_FONTSIZE, fontweight="bold", color="#2A2F36", transform=ax_right.transAxes)
    ax_right.text(0.04, 0.90, results_txt, ha="left", va="top",
                  fontsize=AX_FONTSIZE, fontweight="bold", color="#2A2F36", transform=ax_right.transAxes)
    ax_right.text(0.04, 0.46, "Parameters", ha="left", va="top",
                  fontsize=AX_FONTSIZE, fontweight="bold", color="#2A2F36", transform=ax_right.transAxes)
    ax_right.text(0.04, 0.42, params_txt, ha="left", va="top",
                  fontsize=AX_FONTSIZE, fontweight="bold", color="#2A2F36", transform=ax_right.transAxes)

    pdf.savefig(fig)
    plt.close(fig)

# ---------- Core analysis ----------

def analyze_one(x: np.ndarray,
                skew_method: str,
                skew_stat: str,
                n_perm: int,
                rng) -> Dict[str, float]:
    x = x[np.isfinite(x)]
    n = x.size
    if n < 8:
        raise ValueError("Need at least 8 finite changes for stable tests.")

    pos = int(np.sum(x > 0)); neg = int(np.sum(x < 0)); nz = pos + neg
    prop_pos = float(pos / nz) if nz > 0 else float("nan")
    p_sign = binomial_sign_test(pos, nz) if nz > 0 else float("nan")

    wres = stats.wilcoxon(x, zero_method="wilcox", correction=False,
                          alternative="two-sided", mode="auto")
    wilcoxon_p = float(wres.pvalue)
    HL = hodges_lehmann_onesample(x)

    # Select skewness statistic
    stat_fn = bowley_skewness if skew_stat == "bowley" else fisher_skewness_g1

    if skew_method == "z":
        if stat_fn is not fisher_skewness_g1:
            raise ValueError("--skew-method z requires --skew-stat fisher")
        g1 = fisher_skewness_g1(x)
        se = fisher_skewness_se(n)
        z = float(g1 / se)
        p_skew = float(2.0 * stats.norm.sf(abs(z)))
        skew_val = g1
    else:
        obs, p_skew = sign_randomization_perm_p(x, stat_fn=stat_fn, n_perm=n_perm, rng=rng)
        skew_val = obs

    return dict(
        n_TFs=int(n),
        positives=pos, negatives=neg, nonzero=nz, prop_positive=prop_pos,
        sign_test_p=p_sign,
        wilcoxon_p=wilcoxon_p, HL_estimate=HL,
        skewness_value=skew_val, skewness_p=p_skew
    )
# ---------- Core runner (re-usable by CLI and BINDetect) ----------

def _default_output_paths(results_tsv: str | Path,
                          out_json: str | None,
                          out_pdf: str | None) -> tuple[str, str]:
    """
    If out_json / out_pdf are None, create defaults in the same directory as
    results_tsv:
        <stem>_skew.json
        <stem>_skew.pdf
    where <stem> is the filename without extension.
    """
    results_path = Path(results_tsv)
    stem = results_path.stem  # e.g. "All_motifs_results"
    if out_json is None:
        out_json = str(results_path.with_name(f"{stem}_skew.json"))
    if out_pdf is None:
        out_pdf = str(results_path.with_name(f"{stem}_skew.pdf"))
    return out_json, out_pdf


def run_skew_report(results_tsv: str | Path,
                    out_json: str | None = None,
                    out_pdf: str | None = None,
                    skew_method: str = "perm",
                    skew_stat: str = "bowley",
                    n_perm: int = 5000,
                    seed: int = 1) -> dict:
    """
    Core entry point: do the skew/shift analysis and optionally write JSON/PDF.

    Returns a dict {comparison_key -> summary_dict}.
    Does NOT print to stdout; printing is left to the CLI wrapper.
    """
    results_tsv = str(results_tsv)
    out_json, out_pdf = _default_output_paths(results_tsv, out_json, out_pdf)

    rng = np.random.default_rng(seed)
    df = read_numeric_tsv(results_tsv)
    pairs = autodetect_pairs(df)
    inferred_label = infer_comp_from_filename(results_tsv)

    summaries: Dict[str, Dict[str, float]] = {}
    pdf = PdfPages(out_pdf) if out_pdf else None

    for change_col, pval_col in pairs:
        key_from_cols = change_col[:-len("_change")]  # e.g., "D_HY"
        display_label = f"{key_from_cols}  ({inferred_label})"

        x = pd.to_numeric(df[change_col], errors="coerce").to_numpy()
        res = analyze_one(
            x,
            skew_method=skew_method,
            skew_stat=skew_stat,
            n_perm=n_perm,
            rng=rng,
        )

        res_store = {"change_col": change_col, "pval_col": pval_col, **res}
        summaries[key_from_cols] = res_store

        if pdf is not None:
            # Page 1
            fig, axs = plt.subplots(1, 2, figsize=FIGSIZE, constrained_layout=True)
            page_hist_density(axs[0], x, title=f"{display_label}: Change distribution")

            p = pd.to_numeric(df[pval_col], errors="coerce").to_numpy()
            highlight = None
            hi_col = f"{key_from_cols}_highlighted"
            if hi_col in df.columns:
                h = df[hi_col].astype(str).str.upper() == "TRUE"
                highlight = h.to_numpy()
            if np.isfinite(p).any():
                page_volcano(
                    axs[1], x, p, highlight, title=f"{display_label}: Volcano"
                )
            else:
                axs[1].axis("off")
                axs[1].set_title(f"{display_label}: Volcano (no p-values)")
            pdf.savefig(fig)
            plt.close(fig)

            # Page 2
            fig, axs = plt.subplots(1, 2, figsize=FIGSIZE, constrained_layout=True)
            page_mirrored_ecdf(axs[0], x, title=f"{display_label}: Tails")
            page_sign_bar(
                axs[1],
                res["positives"],
                res["negatives"],
                res["nonzero"],
                res["prop_positive"],
            )
            capsule = (
                f"n={res['n_TFs']}  |  sign p={res['sign_test_p']:.3e}\n"
                f"wilcoxon p={res['wilcoxon_p']:.3e}  |  HL={res['HL_estimate']:.4f}\n"
                f"skew={res['skewness_value']:.4f}  |  p={res['skewness_p']:.3e}\n"
                f"({skew_method}, {skew_stat})"
            )
            axs[0].text(
                0.02,
                0.02,
                capsule,
                transform=axs[0].transAxes,
                fontsize=PDF_FONT_SIZE,
                ha="left",
                va="bottom",
                bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"),
            )
            pdf.savefig(fig)
            plt.close(fig)

            # Page 3 — summary (original language)
            params = dict(
                file=str(results_tsv),
                skew_method=skew_method,
                skew_stat=skew_stat,
                n_perm=int(n_perm),
                seed=int(seed),
            )
            page_summary(
                pdf,
                key=f"{key_from_cols}  ({inferred_label})",
                res={**res_store},
                params=params,
            )

    if out_json:
        with open(out_json, "w", encoding="utf-8") as fh:
            json.dump(summaries, fh, indent=2)

    if pdf is not None:
        pdf.close()

    return summaries

# ---------- CLI ----------

def main():
    ap = argparse.ArgumentParser(description="Skew/shift report for BINDetect results")
    ap.add_argument("results_tsv", help="Path to BINDetect *_results.txt (tab-delimited)")
    ap.add_argument("--out-json", default=None, help="Write JSON summary (all comparisons)")
    ap.add_argument("--out-pdf", default=None, help="Write multi-page PDF per comparison")
    ap.add_argument("--skew-method", choices=["perm", "z"], default="perm",
                    help="Skewness test: 'perm' (sign-randomization; robust) or 'z' (Fisher).")
    ap.add_argument("--skew-stat", choices=["bowley", "fisher"], default="bowley",
                    help="Skewness statistic: robust Bowley (default) or Fisher (moment-based).")
    ap.add_argument("--n-perm", type=int, default=5000, help="Permutations for 'perm' method.")
    ap.add_argument("--seed", type=int, default=1, help="RNG seed.")
    args = ap.parse_args()

    # Always get defaults in same directory if user did not specify paths
    out_json, out_pdf = _default_output_paths(args.results_tsv, args.out_json, args.out_pdf)

    summaries = run_skew_report(
        results_tsv=args.results_tsv,
        out_json=out_json,
        out_pdf=out_pdf,
        skew_method=args.skew_method,
        skew_stat=args.skew_stat,
        n_perm=args.n_perm,
        seed=args.seed,
    )

    # Console summary (same as before)
    print("\n=== BINDetect Skew/Shift Report ===")
    print(f"Input file : {args.results_tsv}")
    for key, res in summaries.items():
        print(f"\n[{key}]")
        print(f"  change_col     : {res['change_col']}")
        print(f"  pval_col       : {res['pval_col']}")
        print(f"  n (finite)     : {res['n_TFs']}")
        print(f"  + / - (nonzero): {res['positives']} / {res['negatives']}  (prop+={res['prop_positive']:.3f})")
        print(f"  Sign test p    : {res['sign_test_p']:.3e}")
        print(f"  Wilcoxon p     : {res['wilcoxon_p']:.3e}")
        print(f"  HL estimate    : {res['HL_estimate']:.4f}")
        print(f"  Skew value     : {res['skewness_value']:.4f}")
        print(f"  Skew p         : {res['skewness_p']:.3e}  "
              f"(method={args.skew_method}, stat={args.skew_stat})")

    print(f"\nWrote JSON → {out_json}")
    print(f"Wrote PDF  → {out_pdf}\n")

# ---------- Programmatic entry points for BINDetect ----------

def generate_skew_report(results_tsv: str | Path,
                         out_json: str | None = None,
                         out_pdf: str | None = None,
                         skew_method: str = "perm",
                         skew_stat: str = "bowley",
                         n_perm: int = 5000,
                         seed: int = 1,
                         **kwargs) -> dict:
    """
    Entry point used by BINDetect’s driver.

    - results_tsv: path to *_results.txt
    - out_json / out_pdf: optional; if None, default names are created in the
      same directory as results_tsv.

    Returns the summaries dict from run_skew_report().
    """
    # Any extra kwargs are ignored to keep this robust to caller changes.
    return run_skew_report(
        results_tsv=results_tsv,
        out_json=out_json,
        out_pdf=out_pdf,
        skew_method=skew_method,
        skew_stat=skew_stat,
        n_perm=n_perm,
        seed=seed,
    )


def main_from_kwargs(**kwargs) -> dict:
    """
    Alternate programmatic entry: accepts keyword args similar to the CLI.
    This is here because some BINDetect versions look for main_from_kwargs().
    """
    results_tsv = kwargs.pop("results_tsv")
    out_json = kwargs.pop("out_json", None)
    out_pdf = kwargs.pop("out_pdf", None)
    skew_method = kwargs.pop("skew_method", "perm")
    skew_stat = kwargs.pop("skew_stat", "bowley")
    n_perm = int(kwargs.pop("n_perm", 5000))
    seed = int(kwargs.pop("seed", 1))

    return generate_skew_report(
        results_tsv=results_tsv,
        out_json=out_json,
        out_pdf=out_pdf,
        skew_method=skew_method,
        skew_stat=skew_stat,
        n_perm=n_perm,
        seed=seed,
    )

if __name__ == "__main__":
    main()
