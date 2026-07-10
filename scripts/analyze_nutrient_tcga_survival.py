#!/usr/bin/env python3
"""TCGA-PAAD survival analysis for nutrient-stress motif-linked TFs."""

from __future__ import annotations

import argparse
import gzip
import math
import shutil
import urllib.request
import zipfile
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
from scipy import optimize, stats


PAAD_EXPR_URL = "https://tcga.xenahubs.net/download/TCGA.PAAD.sampleMap/HiSeqV2.gz"
PANCAN_SURVIVAL_URL = (
    "https://pancanatlas.xenahubs.net/download/"
    "Survival_SupplementalTable_S1_20171025_xena_sp"
)


@dataclass
class Endpoint:
    name: str
    time_col: str
    event_col: str


ENDPOINTS = [
    Endpoint("OS", "OS.time", "OS"),
    Endpoint("DSS", "DSS.time", "DSS"),
    Endpoint("PFI", "PFI.time", "PFI"),
    Endpoint("DFI", "DFI.time", "DFI"),
]


KM_FONT = {
    "family": "sans-serif",
    "sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "size": 9,
    "weight": "bold",
}


def apply_km_style() -> None:
    plt.rcParams.update(
        {
            "font.family": KM_FONT["family"],
            "font.sans-serif": KM_FONT["sans-serif"],
            "font.size": KM_FONT["size"],
            "font.weight": KM_FONT["weight"],
            "axes.labelweight": "bold",
            "axes.titleweight": "bold",
            "axes.linewidth": 0.9,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def bold_axis_text(ax) -> None:
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(9)
        label.set_fontweight("bold")
    ax.xaxis.label.set_fontsize(9)
    ax.xaxis.label.set_fontweight("bold")
    ax.yaxis.label.set_fontsize(9)
    ax.yaxis.label.set_fontweight("bold")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Analyze TCGA-PAAD survival for nutrient-stress motif-linked TFs."
    )
    p.add_argument(
        "--handoff-dir",
        default="data/public/processed/nutrient_tf_handoff_20260629",
        help="Nutrient TF handoff directory with tables/ subdirectory.",
    )
    p.add_argument(
        "--rna-raw",
        default="data/public/raw/nutrient_rna/nutrient_rna_raw_counts_ruvr_k20_gene_universe.tsv.gz",
        help="Local nutrient raw RNA count matrix.",
    )
    p.add_argument(
        "--outdir",
        default=None,
        help="Output directory. Defaults to data/public/processed/nutrient_tcga_survival_YYYYMMDD_HHMMSS.",
    )
    p.add_argument("--min-local-raw-mean", type=float, default=1.0)
    p.add_argument("--min-group-size", type=int, default=10)
    p.add_argument("--min-events", type=int, default=8)
    p.add_argument("--plot-top", type=int, default=40)
    p.add_argument("--km-dpi", type=int, default=150)
    p.add_argument("--skip-km-plots", action="store_true")
    p.add_argument("--force-download", action="store_true")
    return p.parse_args()


def safe_gene(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if not text:
        return ""
    if "|" in text:
        text = text.split("|", 1)[0]
    return text.upper()


def safe_filename(text: str) -> str:
    keep = []
    for ch in str(text):
        if ch.isalnum() or ch in "._-":
            keep.append(ch)
        else:
            keep.append("_")
    return "".join(keep).strip("_") or "item"


def download(url: str, path: Path, force: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.stat().st_size > 0 and not force:
        return
    tmp = path.with_suffix(path.suffix + ".tmp")
    with urllib.request.urlopen(url, timeout=120) as response, tmp.open("wb") as out:
        shutil.copyfileobj(response, out)
    tmp.replace(path)


def bh_fdr(pvalues: pd.Series) -> pd.Series:
    p = pd.to_numeric(pvalues, errors="coerce")
    out = pd.Series(np.nan, index=p.index, dtype=float)
    ok = p.notna() & np.isfinite(p) & (p >= 0) & (p <= 1)
    if not ok.any():
        return out
    vals = p[ok].to_numpy(float)
    order = np.argsort(vals)
    ranked = vals[order]
    n = len(vals)
    adj = ranked * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    restored = np.empty_like(adj)
    restored[order] = adj
    out.loc[ok] = restored
    return out


def load_nutrient_tf_universe(handoff_dir: Path, rna_raw: Path, min_raw_mean: float):
    summary = pd.read_csv(handoff_dir / "tables" / "tf_condition_summary.tsv", sep="\t")
    recurrence = pd.read_csv(handoff_dir / "tables" / "cross_cellline_recurrence.tsv", sep="\t")
    tf_genes = sorted({safe_gene(x) for x in summary["tf_gene"].dropna() if safe_gene(x)})

    raw = pd.read_csv(rna_raw, sep="\t", compression="gzip")
    sample_cols = [
        c for c in raw.columns if c not in {"gene_key", "ensembl_gene_id", "HGNC"}
    ]
    raw["gene_symbol"] = raw["HGNC"].map(safe_gene)
    raw = raw[raw["gene_symbol"] != ""].copy()
    raw["local_raw_mean"] = raw[sample_cols].mean(axis=1)
    raw["local_raw_max"] = raw[sample_cols].max(axis=1)
    expr = (
        raw.groupby("gene_symbol", as_index=False)[["local_raw_mean", "local_raw_max"]]
        .max()
        .sort_values("gene_symbol")
    )
    universe = pd.DataFrame({"gene": tf_genes})
    universe = universe.merge(expr, left_on="gene", right_on="gene_symbol", how="left")
    universe.drop(columns=["gene_symbol"], inplace=True)
    universe["locally_expressed"] = universe["local_raw_mean"].fillna(0) >= min_raw_mean

    nutrient_gene = (
        summary.groupby("tf_gene")
        .agg(
            nutrient_best_fdr=("best_fdr", "min"),
            nutrient_max_abs_delta_fp=("max_abs_delta_fp", "max"),
            nutrient_max_abs_rna_log2fc=("rna_log2fc_max_abs", "max"),
            nutrient_n_cell_lines=("cell_line", "nunique"),
            nutrient_n_conditions=("condition", "nunique"),
        )
        .reset_index()
    )
    nutrient_gene["gene"] = nutrient_gene["tf_gene"].map(safe_gene)
    nutrient_gene = nutrient_gene.drop(columns=["tf_gene"]).groupby("gene", as_index=False).max()

    recurrence_gene = (
        recurrence.groupby("tf_gene")
        .agg(
            recurrence_n_stress_types=("stress_type", "nunique"),
            recurrence_max_cell_lines=("n_cell_lines", "max"),
            recurrence_best_fdr=("best_fdr", "min"),
        )
        .reset_index()
    )
    recurrence_gene["gene"] = recurrence_gene["tf_gene"].map(safe_gene)
    recurrence_gene = (
        recurrence_gene.drop(columns=["tf_gene"]).groupby("gene", as_index=False).max()
    )
    return universe, nutrient_gene, recurrence_gene


def load_tcga_expression(expr_gz: Path) -> pd.DataFrame:
    expr = pd.read_csv(expr_gz, sep="\t", compression="gzip")
    first_col = expr.columns[0]
    expr = expr.rename(columns={first_col: "gene"})
    expr["gene"] = expr["gene"].map(safe_gene)
    expr = expr[expr["gene"] != ""].copy()
    sample_cols = [c for c in expr.columns if c != "gene"]
    for c in sample_cols:
        expr[c] = pd.to_numeric(expr[c], errors="coerce")
    expr = expr.groupby("gene", as_index=True)[sample_cols].mean()
    tumor_cols = [
        c for c in expr.columns if len(c.split("-")) >= 4 and c.split("-")[3][:2] == "01"
    ]
    return expr[tumor_cols]


def load_survival(survival_path: Path) -> pd.DataFrame:
    surv = pd.read_csv(survival_path, sep="\t")
    surv = surv[surv["cancer type abbreviation"] == "PAAD"].copy()
    surv["sample"] = surv["sample"].astype(str)
    for ep in ENDPOINTS:
        surv[ep.time_col] = pd.to_numeric(surv[ep.time_col], errors="coerce")
        surv[ep.event_col] = pd.to_numeric(surv[ep.event_col], errors="coerce")
    return surv


def kaplan_meier(time: np.ndarray, event: np.ndarray):
    order = np.argsort(time)
    time = time[order]
    event = event[order]
    xs = [0.0]
    ys = [1.0]
    s = 1.0
    for t in np.unique(time[event == 1]):
        at_risk = np.sum(time >= t)
        events = np.sum((time == t) & (event == 1))
        if at_risk <= 0:
            continue
        xs.extend([t, t])
        ys.extend([s, s * (1.0 - events / at_risk)])
        s = ys[-1]
    if len(time):
        xs.append(float(np.max(time)))
        ys.append(s)
    return np.asarray(xs), np.asarray(ys)


def median_survival(time: np.ndarray, event: np.ndarray) -> float:
    xs, ys = kaplan_meier(time, event)
    below = np.where(ys <= 0.5)[0]
    if len(below) == 0:
        return np.nan
    return xs[below[0]]


def logrank_test(time: np.ndarray, event: np.ndarray, group: np.ndarray) -> tuple[float, float]:
    group = group.astype(int)
    observed_minus_expected = 0.0
    variance = 0.0
    for t in np.unique(time[event == 1]):
        risk = time >= t
        evt = (time == t) & (event == 1)
        n = risk.sum()
        n1 = (risk & (group == 1)).sum()
        n0 = n - n1
        d = evt.sum()
        d1 = (evt & (group == 1)).sum()
        if n <= 1 or n1 == 0 or n0 == 0 or d == 0:
            continue
        expected = d * n1 / n
        var = d * (n - d) * n1 * n0 / (n * n * (n - 1))
        observed_minus_expected += d1 - expected
        variance += var
    if variance <= 0:
        return np.nan, np.nan
    chi2 = observed_minus_expected * observed_minus_expected / variance
    return float(chi2), float(stats.chi2.sf(chi2, 1))


def cox_univariate(time: np.ndarray, event: np.ndarray, x: np.ndarray) -> dict[str, float]:
    ok = np.isfinite(time) & np.isfinite(event) & np.isfinite(x)
    time = time[ok].astype(float)
    event = event[ok].astype(int)
    x = x[ok].astype(float)
    if len(time) < 20 or event.sum() < 8 or np.nanstd(x) <= 0:
        return {}
    x = (x - np.mean(x)) / np.std(x)
    event_times = np.unique(time[event == 1])

    def neg_loglik(beta: float) -> float:
        bx = beta * x
        total = float(np.sum(bx[event == 1]))
        for t in event_times:
            d = np.sum((time == t) & (event == 1))
            risk_sum = np.sum(np.exp(bx[time >= t]))
            if risk_sum <= 0:
                return np.inf
            total -= d * math.log(risk_sum)
        return -total

    result = optimize.minimize_scalar(neg_loglik, bounds=(-8, 8), method="bounded")
    if not result.success or not np.isfinite(result.x):
        return {}
    beta = float(result.x)
    bx = beta * x
    hess = 0.0
    for t in event_times:
        d = np.sum((time == t) & (event == 1))
        risk = time >= t
        w = np.exp(bx[risk])
        if w.sum() <= 0:
            continue
        xm = np.sum(w * x[risk]) / np.sum(w)
        x2m = np.sum(w * x[risk] * x[risk]) / np.sum(w)
        hess -= d * (x2m - xm * xm)
    if hess >= 0 or not np.isfinite(hess):
        return {}
    se = math.sqrt(-1.0 / hess)
    z = beta / se
    p = 2 * stats.norm.sf(abs(z))
    return {
        "cox_beta_zexpr": beta,
        "cox_hr_zexpr": math.exp(beta),
        "cox_ci95_low_zexpr": math.exp(beta - 1.96 * se),
        "cox_ci95_high_zexpr": math.exp(beta + 1.96 * se),
        "cox_p_zexpr": p,
    }


def analyze_endpoint(
    expr: pd.DataFrame,
    survival: pd.DataFrame,
    genes: list[str],
    endpoint: Endpoint,
    min_group_size: int,
    min_events: int,
) -> tuple[pd.DataFrame, dict[tuple[str, str], pd.DataFrame]]:
    rows = []
    plot_data = {}
    available_samples = [s for s in expr.columns if s in set(survival["sample"])]
    surv = survival.set_index("sample").loc[available_samples].copy()
    time = pd.to_numeric(surv[endpoint.time_col], errors="coerce")
    event = pd.to_numeric(surv[endpoint.event_col], errors="coerce")
    valid_surv = time.notna() & event.notna() & (time > 0)
    for gene in genes:
        if gene not in expr.index:
            continue
        vals = expr.loc[gene, available_samples]
        df = pd.DataFrame(
            {
                "sample": available_samples,
                "expr": pd.to_numeric(vals, errors="coerce").to_numpy(float),
                "time_days": time.to_numpy(float),
                "event": event.to_numpy(float),
            }
        )
        df = df[valid_surv.to_numpy() & np.isfinite(df["expr"])].copy()
        if df.empty or df["expr"].nunique() < 2:
            rows.append({"gene": gene, "endpoint": endpoint.name, "status": "no_expression_variance"})
            continue
        median_expr = float(df["expr"].median())
        df["group"] = np.where(df["expr"] > median_expr, "High", "Low")
        n_high = int((df["group"] == "High").sum())
        n_low = int((df["group"] == "Low").sum())
        events = int(df["event"].sum())
        if n_high < min_group_size or n_low < min_group_size or events < min_events:
            rows.append(
                {
                    "gene": gene,
                    "endpoint": endpoint.name,
                    "status": "insufficient_group_or_events",
                    "n": len(df),
                    "n_high": n_high,
                    "n_low": n_low,
                    "events": events,
                    "median_expr_cutoff": median_expr,
                }
            )
            continue
        group_binary = (df["group"] == "High").to_numpy(int)
        chi2, lr_p = logrank_test(
            df["time_days"].to_numpy(float), df["event"].to_numpy(int), group_binary
        )
        cox = cox_univariate(
            df["time_days"].to_numpy(float), df["event"].to_numpy(int), df["expr"].to_numpy(float)
        )
        low = df[df["group"] == "Low"]
        high = df[df["group"] == "High"]
        row = {
            "gene": gene,
            "endpoint": endpoint.name,
            "status": "tested",
            "n": len(df),
            "n_high": n_high,
            "n_low": n_low,
            "events": events,
            "events_high": int(high["event"].sum()),
            "events_low": int(low["event"].sum()),
            "median_expr_cutoff": median_expr,
            "expr_low_mean": float(low["expr"].mean()),
            "expr_high_mean": float(high["expr"].mean()),
            "median_survival_high_days": median_survival(
                high["time_days"].to_numpy(float), high["event"].to_numpy(int)
            ),
            "median_survival_low_days": median_survival(
                low["time_days"].to_numpy(float), low["event"].to_numpy(int)
            ),
            "logrank_chi2": chi2,
            "logrank_p": lr_p,
        }
        row.update(cox)
        rows.append(row)
        plot_data[(gene, endpoint.name)] = df
    return pd.DataFrame(rows), plot_data


def plot_km(gene: str, endpoint: str, df: pd.DataFrame, row: pd.Series, out: Path, dpi: int):
    apply_km_style()
    fig, ax = plt.subplots(figsize=(3.6, 3.6))
    colors = {"Low": "#2f6fb0", "High": "#c43c39"}
    for group in ["Low", "High"]:
        sub = df[df["group"] == group]
        xs, ys = kaplan_meier(sub["time_days"].to_numpy(float), sub["event"].to_numpy(int))
        ax.step(xs / 30.4375, ys, where="post", color=colors[group], lw=2.2, label=f"{group} ({len(sub)})")
    p = row.get("logrank_p", np.nan)
    hr = row.get("cox_hr_zexpr", np.nan)
    fdr = row.get("logrank_fdr", np.nan)
    text = f"log-rank p={p:.2g}\nFDR={fdr:.2g}\nHR/z={hr:.2f}" if np.isfinite(p) else ""
    ax.text(0.98, 0.05, text, ha="right", va="bottom", transform=ax.transAxes, fontsize=9, fontweight="bold")
    ax.set_title(f"{gene} TCGA-PAAD {endpoint}", fontsize=9, fontweight="bold")
    ax.set_xlabel("Time (months)", fontsize=9, fontweight="bold")
    ax.set_ylabel("Survival probability", fontsize=9, fontweight="bold")
    ax.set_ylim(-0.03, 1.03)
    ax.set_box_aspect(1)
    ax.tick_params(axis="both", labelsize=9, width=0.9)
    bold_axis_text(ax)
    ax.legend(frameon=False, loc="upper right", prop={"size": 9, "weight": "bold"})
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=dpi)
    plt.close(fig)


def plot_summary(results: pd.DataFrame, merged: pd.DataFrame, outdir: Path, plot_top: int):
    outdir.mkdir(parents=True, exist_ok=True)
    tested = results[results["status"] == "tested"].copy()
    if tested.empty:
        return
    os_res = tested[tested["endpoint"] == "OS"].sort_values(["logrank_fdr", "logrank_p"]).head(plot_top)
    if not os_res.empty:
        fig, ax = plt.subplots(figsize=(8.5, max(5, 0.18 * len(os_res))))
        y = np.arange(len(os_res))
        vals = -np.log10(os_res["logrank_fdr"].clip(lower=1e-300))
        colors = np.where(os_res["cox_hr_zexpr"] > 1, "#c43c39", "#2f6fb0")
        ax.barh(y, vals, color=colors)
        ax.set_yticks(y)
        ax.set_yticklabels(os_res["gene"], fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel("-log10 OS log-rank FDR")
        ax.set_title("Top TCGA-PAAD survival-associated nutrient TFs", fontweight="bold")
        ax.spines[["top", "right"]].set_visible(False)
        fig.tight_layout()
        fig.savefig(outdir / "top_os_survival_tfs.png", dpi=180)
        fig.savefig(outdir / "top_os_survival_tfs.pdf")
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.5, 5.2))
    os_all = tested[tested["endpoint"] == "OS"].copy()
    if not os_all.empty:
        x = np.log2(os_all["cox_hr_zexpr"].clip(lower=1e-6))
        y = -np.log10(os_all["logrank_p"].clip(lower=1e-300))
        sig = os_all["logrank_fdr"] < 0.1
        ax.scatter(x[~sig], y[~sig], s=18, color="#999999", alpha=0.65, label="FDR >= 0.1")
        ax.scatter(x[sig], y[sig], s=30, color="#c43c39", alpha=0.9, label="FDR < 0.1")
        for _, r in os_all.sort_values(["logrank_fdr", "logrank_p"]).head(15).iterrows():
            ax.text(np.log2(r["cox_hr_zexpr"]), -np.log10(max(r["logrank_p"], 1e-300)), r["gene"], fontsize=7)
        ax.axvline(0, color="#333333", lw=0.8)
        ax.set_xlabel("log2 Cox HR per 1 SD expression")
        ax.set_ylabel("-log10 log-rank p")
        ax.set_title("TCGA-PAAD OS survival scan", fontweight="bold")
        ax.legend(frameon=False, fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)
        fig.tight_layout()
        fig.savefig(outdir / "os_survival_volcano.png", dpi=180)
        fig.savefig(outdir / "os_survival_volcano.pdf")
        plt.close(fig)

    heat = merged[merged["endpoint"] == "OS"].copy()
    heat = heat[heat["status"] == "tested"].copy()
    if not heat.empty:
        heat["combined_rank_score"] = (
            -np.log10(heat["logrank_p"].clip(lower=1e-300))
            + np.log10(heat["nutrient_max_abs_delta_fp"].fillna(0) + 1.0)
            + heat["nutrient_n_cell_lines"].fillna(0)
        )
        top = heat.sort_values("combined_rank_score", ascending=False).head(min(plot_top, 50))
        cols = [
            "nutrient_max_abs_delta_fp",
            "nutrient_max_abs_rna_log2fc",
            "nutrient_n_cell_lines",
            "recurrence_n_stress_types",
            "cox_hr_zexpr",
            "logrank_fdr",
        ]
        mat = top.set_index("gene")[cols].copy()
        mat["logrank_fdr"] = -np.log10(mat["logrank_fdr"].clip(lower=1e-300))
        mat["cox_hr_zexpr"] = np.log2(mat["cox_hr_zexpr"].clip(lower=1e-6))
        labels = [
            "max |delta FP|",
            "max |RNA log2FC|",
            "cell lines",
            "stress types",
            "log2 HR",
            "-log10 FDR",
        ]
        scaled = mat.copy()
        for c in scaled.columns:
            v = scaled[c].astype(float)
            if v.max() != v.min():
                scaled[c] = (v - v.min()) / (v.max() - v.min())
            else:
                scaled[c] = 0.5
        fig, ax = plt.subplots(figsize=(8, max(5, 0.18 * len(scaled))))
        im = ax.imshow(scaled.to_numpy(), aspect="auto", cmap="viridis")
        ax.set_xticks(np.arange(len(labels)))
        ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=8)
        ax.set_yticks(np.arange(len(scaled)))
        ax.set_yticklabels(scaled.index, fontsize=8)
        ax.set_title("Integrated nutrient footprint/RNA and TCGA OS evidence", fontweight="bold")
        fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02, label="column-scaled value")
        fig.tight_layout()
        fig.savefig(outdir / "integrated_nutrient_tcga_heatmap.png", dpi=180)
        fig.savefig(outdir / "integrated_nutrient_tcga_heatmap.pdf")
        plt.close(fig)


def plot_km_axis(ax, gene: str, endpoint: str, df: pd.DataFrame, row: pd.Series) -> None:
    apply_km_style()
    colors = {"Low": "#2f6fb0", "High": "#c43c39"}
    for group in ["Low", "High"]:
        sub = df[df["group"] == group]
        xs, ys = kaplan_meier(sub["time_days"].to_numpy(float), sub["event"].to_numpy(int))
        ax.step(
            xs / 30.4375,
            ys,
            where="post",
            color=colors[group],
            lw=1.7,
            label=f"{group} ({len(sub)})",
        )
    p = row.get("logrank_p", np.nan)
    fdr = row.get("logrank_fdr", np.nan)
    hr = row.get("cox_hr_zexpr", np.nan)
    stats_text = f"p={p:.1g}  FDR={fdr:.1g}\nHR/z={hr:.2f}" if np.isfinite(p) else ""
    ax.text(0.98, 0.05, stats_text, ha="right", va="bottom", transform=ax.transAxes, fontsize=9, fontweight="bold")
    ax.set_title(f"{gene} {endpoint}", fontsize=9, fontweight="bold")
    ax.set_xlabel("Months", fontsize=9, fontweight="bold")
    ax.set_ylabel("Survival", fontsize=9, fontweight="bold")
    ax.set_ylim(-0.03, 1.03)
    ax.set_box_aspect(1)
    ax.tick_params(axis="both", labelsize=9, width=0.9)
    bold_axis_text(ax)
    ax.legend(frameon=False, loc="upper right", prop={"size": 9, "weight": "bold"})
    ax.spines[["top", "right"]].set_visible(False)


def plot_km_pdf_panels(
    results: pd.DataFrame,
    plot_data: dict[tuple[str, str], pd.DataFrame],
    outdir: Path,
    panels_per_page: int = 12,
) -> None:
    apply_km_style()
    panel_dir = outdir / "figures" / "km_panel_pdfs"
    panel_dir.mkdir(parents=True, exist_ok=True)
    lookup = results.set_index(["gene", "endpoint"])
    endpoint_paths = []
    for endpoint in sorted(results["endpoint"].dropna().unique()):
        ranked = (
            results[(results["endpoint"] == endpoint) & (results["status"] == "tested")]
            .sort_values(["logrank_fdr", "logrank_p", "gene"])
            .copy()
        )
        pdf_path = panel_dir / f"{safe_filename(endpoint)}_KM_panels_ranked.pdf"
        endpoint_paths.append(pdf_path)
        with PdfPages(pdf_path) as pdf:
            for start in range(0, len(ranked), panels_per_page):
                chunk = ranked.iloc[start : start + panels_per_page]
                fig, axes = plt.subplots(4, 3, figsize=(9.5, 12.5))
                axes = axes.ravel()
                for ax in axes:
                    ax.axis("off")
                for ax, (_, row) in zip(axes, chunk.iterrows()):
                    ax.axis("on")
                    key = (row["gene"], row["endpoint"])
                    plot_km_axis(ax, row["gene"], row["endpoint"], plot_data[key], lookup.loc[key])
                page = start // panels_per_page + 1
                fig.suptitle(
                    f"TCGA-PAAD {endpoint} Kaplan-Meier curves ranked by log-rank FDR - page {page}",
                    fontsize=12,
                    fontweight="bold",
                    y=0.995,
                )
                fig.tight_layout(rect=[0, 0, 1, 0.975], h_pad=0.7, w_pad=0.5)
                pdf.savefig(fig)
                plt.close(fig)

    combined_path = panel_dir / "all_endpoints_KM_panels_ranked.pdf"
    with PdfPages(combined_path) as pdf:
        for pdf_path in endpoint_paths:
            # Re-open endpoint PDFs is not portable without extra dependencies, so
            # regenerate into the combined file in endpoint order.
            endpoint = pdf_path.name.split("_", 1)[0]
            ranked = (
                results[(results["endpoint"] == endpoint) & (results["status"] == "tested")]
                .sort_values(["logrank_fdr", "logrank_p", "gene"])
                .copy()
            )
            for start in range(0, len(ranked), panels_per_page):
                chunk = ranked.iloc[start : start + panels_per_page]
                fig, axes = plt.subplots(4, 3, figsize=(9.5, 12.5))
                axes = axes.ravel()
                for ax in axes:
                    ax.axis("off")
                for ax, (_, row) in zip(axes, chunk.iterrows()):
                    ax.axis("on")
                    key = (row["gene"], row["endpoint"])
                    plot_km_axis(ax, row["gene"], row["endpoint"], plot_data[key], lookup.loc[key])
                page = start // panels_per_page + 1
                fig.suptitle(
                    f"TCGA-PAAD {endpoint} Kaplan-Meier curves ranked by log-rank FDR - page {page}",
                    fontsize=12,
                    fontweight="bold",
                    y=0.995,
                )
                fig.tight_layout(rect=[0, 0, 1, 0.975], h_pad=0.7, w_pad=0.5)
                pdf.savefig(fig)
                plt.close(fig)


def write_readme(outdir: Path, args: argparse.Namespace, counts: dict[str, int], upload_hint: str):
    text = f"""# Nutrient-Stress TF TCGA Survival Package

Generated: {datetime.now().isoformat(timespec='seconds')}

This package tests motif-linked, locally expressed nutrient-stress transcription factors
against TCGA-PAAD survival endpoints from UCSC Xena.

## Inputs

- Nutrient TF handoff: `{args.handoff_dir}`
- Local raw RNA expression: `{args.rna_raw}`
- TCGA-PAAD expression: `{PAAD_EXPR_URL}`
- PanCanAtlas survival endpoints: `{PANCAN_SURVIVAL_URL}`

## Filters

- Motif-linked TFs found in the nutrient footprint/RNA handoff tables.
- Local expression filter: raw RNA mean >= {args.min_local_raw_mean:g}.
- TCGA filter: gene present in the PAAD expression matrix and variable across samples.
- Survival test filter: each high/low expression group has at least {args.min_group_size} samples and the endpoint has at least {args.min_events} events.

## Statistics

For each TF and endpoint, tumor samples were split into high and low expression groups
using the median TCGA expression value. Survival curves use the Kaplan-Meier estimator.
The high-vs-low comparison uses a two-sided log-rank test. A univariate Cox proportional
hazards model was also fit using z-scored continuous TF expression. P-values are adjusted
within each endpoint using Benjamini-Hochberg FDR correction.

## Output Counts

- Motif-linked TF genes: {counts.get('motif_linked_tfs', 0)}
- Locally expressed motif-linked TF genes: {counts.get('locally_expressed_tfs', 0)}
- TF genes available in TCGA-PAAD expression: {counts.get('tcga_available_tfs', 0)}
- Tested endpoint rows: {counts.get('tested_endpoint_rows', 0)}

## Main Files

- `tables/all_tf_survival_results.tsv`: all endpoint-level survival test rows.
- `tables/tf_survival_ranked_by_endpoint.tsv`: tested rows ranked by endpoint FDR.
- `tables/nutrient_tf_tcga_merged_evidence.tsv`: nutrient footprint/RNA evidence merged with TCGA survival statistics.
- `tables/tf_universe_filter_audit.tsv`: inclusion and exclusion reasons for every motif-linked TF.
- `figures/km_curves/`: Kaplan-Meier curves for tested TF-endpoint pairs.
- `figures/km_panel_pdfs/`: multi-page PDF panels containing the same curves ranked by endpoint FDR.
- `figures/summary/`: ranked and integrated summary figures.

## Upload

Suggested Box upload target:

`{upload_hint}`
"""
    (outdir / "README.md").write_text(text)
    (outdir / "methods_survival_analysis.md").write_text(
        "\n".join(text.splitlines()[10:42]) + "\n"
    )


def main() -> None:
    args = parse_args()
    handoff_dir = Path(args.handoff_dir)
    outdir = Path(
        args.outdir
        or f"data/public/processed/nutrient_tcga_survival_{datetime.now():%Y%m%d_%H%M%S}"
    )
    tables_dir = outdir / "tables"
    tcga_dir = outdir / "tcga_downloads"
    figures_dir = outdir / "figures"
    km_dir = figures_dir / "km_curves"
    summary_fig_dir = figures_dir / "summary"
    for d in [tables_dir, tcga_dir, km_dir, summary_fig_dir]:
        d.mkdir(parents=True, exist_ok=True)

    expr_path = tcga_dir / "TCGA_PAAD_HiSeqV2.gz"
    survival_path = tcga_dir / "Survival_SupplementalTable_S1_20171025_xena_sp.tsv"
    download(PAAD_EXPR_URL, expr_path, force=args.force_download)
    download(PANCAN_SURVIVAL_URL, survival_path, force=args.force_download)

    universe, nutrient_gene, recurrence_gene = load_nutrient_tf_universe(
        handoff_dir, Path(args.rna_raw), args.min_local_raw_mean
    )
    tcga_expr = load_tcga_expression(expr_path)
    survival = load_survival(survival_path)
    available = set(tcga_expr.index)
    universe["tcga_available"] = universe["gene"].isin(available)
    universe["testable"] = universe["locally_expressed"] & universe["tcga_available"]
    universe["exclusion_reason"] = ""
    universe.loc[~universe["locally_expressed"], "exclusion_reason"] = "not_locally_expressed"
    universe.loc[
        universe["locally_expressed"] & ~universe["tcga_available"], "exclusion_reason"
    ] = "not_in_tcga_paad_expression"
    universe.loc[universe["testable"], "exclusion_reason"] = "included"
    universe.to_csv(tables_dir / "tf_universe_filter_audit.tsv", sep="\t", index=False)

    genes = universe.loc[universe["testable"], "gene"].sort_values().tolist()
    all_results = []
    all_plot_data = {}
    for ep in ENDPOINTS:
        res, pdata = analyze_endpoint(
            tcga_expr, survival, genes, ep, args.min_group_size, args.min_events
        )
        all_results.append(res)
        all_plot_data.update(pdata)
    results = pd.concat(all_results, ignore_index=True)
    for ep_name, idx in results.groupby("endpoint").groups.items():
        mask = results.index.isin(idx) & (results["status"] == "tested")
        results.loc[mask, "logrank_fdr"] = bh_fdr(results.loc[mask, "logrank_p"])
        results.loc[mask, "cox_fdr_zexpr"] = bh_fdr(results.loc[mask, "cox_p_zexpr"])
    results.to_csv(tables_dir / "all_tf_survival_results.tsv", sep="\t", index=False)

    ranked = (
        results[results["status"] == "tested"]
        .sort_values(["endpoint", "logrank_fdr", "logrank_p", "gene"])
        .copy()
    )
    ranked.to_csv(tables_dir / "tf_survival_ranked_by_endpoint.tsv", sep="\t", index=False)

    merged = results.merge(nutrient_gene, on="gene", how="left").merge(
        recurrence_gene, on="gene", how="left"
    )
    merged.to_csv(tables_dir / "nutrient_tf_tcga_merged_evidence.tsv", sep="\t", index=False)

    if not args.skip_km_plots:
        result_lookup = results.set_index(["gene", "endpoint"])
        for key, df in all_plot_data.items():
            gene, ep = key
            row = result_lookup.loc[key]
            plot_km(
                gene,
                ep,
                df,
                row,
                km_dir / f"{safe_filename(gene)}_{safe_filename(ep)}_KM.png",
                args.km_dpi,
            )
        plot_km_pdf_panels(results, all_plot_data, outdir)

    plot_summary(results, merged, summary_fig_dir, args.plot_top)

    counts = {
        "motif_linked_tfs": int(len(universe)),
        "locally_expressed_tfs": int(universe["locally_expressed"].sum()),
        "tcga_available_tfs": int(universe["testable"].sum()),
        "tested_endpoint_rows": int((results["status"] == "tested").sum()),
    }
    upload_hint = (
        "yilab:Yaoxiang/thesis_project/manuscript/review_nutrient_stress_tfs/"
        f"{outdir.name}"
    )
    write_readme(outdir, args, counts, upload_hint)

    zip_path = outdir.with_suffix(".zip")
    if zip_path.exists():
        zip_path.unlink()
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in outdir.rglob("*"):
            if path.is_file():
                zf.write(path, path.relative_to(outdir.parent))

    print(f"Output: {outdir}")
    print(f"Archive: {zip_path}")
    print(f"Tested endpoint rows: {counts['tested_endpoint_rows']}")


if __name__ == "__main__":
    main()
