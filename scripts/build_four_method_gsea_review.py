#!/usr/bin/env python3
"""Build a four-method ENCODE normalization-comparison review folder.

This supersedes the earlier GSEA-like review. Methods 2--4 now compare motif
level score matrices directly, and all methods use the same strict significance
gate for volcano coloring and aggregate-profile selection.
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import scipy.stats

from fp_tools.tools.bindetect import _existing_result_motifs
from fp_tools.tools.bindetect_functions import (
    build_bindetect_aggregate_payload,
    plot_interactive_bindetect,
)


SAMPLES = ["K562_rep1", "K562_rep2", "K562_rep3", "HepG2_rep1", "HepG2_rep2", "HepG2_rep3"]
GROUP1 = SAMPLES[:3]
GROUP2 = SAMPLES[3:]
COMPARISON = "K562_HepG2"
COMPARISON_TUPLE = ("K562", "HepG2")
SIG_P_CUTOFF = 0.001
SIG_DELTA_CUTOFF = 0.1
SIG_BOUND_CUTOFF = 500
AGGREGATE_TOP_N = 500
PSEUDOCOUNT_SCALE = 1000.0


def bh(pvalues: np.ndarray) -> np.ndarray:
    pvals = np.asarray(pvalues, dtype=float)
    out = np.full(pvals.shape, np.nan)
    ok = np.isfinite(pvals)
    if not ok.any():
        return out
    clipped = np.clip(pvals[ok], 0.0, 1.0)
    order = np.argsort(clipped)
    ranked = clipped[order]
    adj = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    back = np.empty_like(adj)
    back[order] = np.clip(adj, 0.0, 1.0)
    out[ok] = back
    return out


def bonferroni(pvalues: pd.Series) -> pd.Series:
    vals = pd.to_numeric(pvalues, errors="coerce").fillna(1.0)
    return (vals * len(vals)).clip(upper=1.0)


def read_results(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    for col, default in [
        (f"{COMPARISON}_pvalue", 1.0),
        (f"{COMPARISON}_qvalue_bh", 1.0),
        (f"{COMPARISON}_change", 0.0),
        (f"{COMPARISON}_mean_delta_fp", 0.0),
        ("K562_bound", 0),
        ("HepG2_bound", 0),
        ("total_tfbs", 0),
    ]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").fillna(default)
    return df


def bed_score_columns(n_fields: int, n_samples: int = 6, n_conditions: int = 2) -> tuple[int, int]:
    sample_start = n_fields - n_samples - (2 * n_conditions)
    if sample_start < 6:
        raise ValueError(f"Cannot infer sample score columns from {n_fields} fields")
    return sample_start, sample_start + n_samples


def iter_score_chunks(run_dir: Path, motif: str, chunksize: int = 200_000):
    bed = run_dir / motif / "beds" / f"{motif}_all.bed"
    if not bed.exists():
        return
    with bed.open("r", encoding="utf-8") as handle:
        first = handle.readline()
    if not first.strip():
        return
    n_fields = len(first.rstrip("\n").split("\t"))
    start, end = bed_score_columns(n_fields)
    for chunk in pd.read_csv(
        bed,
        sep="\t",
        header=None,
        usecols=list(range(start, end)),
        names=SAMPLES,
        chunksize=chunksize,
        comment="#",
    ):
        yield chunk.apply(pd.to_numeric, errors="coerce")


def build_matrix_from_beds(run_dir: Path, motifs: list[str], out_prefix: Path) -> pd.DataFrame:
    rows = []
    site_summary = []
    for motif in motifs:
        sums = np.zeros(len(SAMPLES), dtype=float)
        count = 0
        for chunk in iter_score_chunks(run_dir, motif):
            vals = chunk.to_numpy(dtype=float)
            finite = np.isfinite(vals).all(axis=1)
            vals = vals[finite]
            if vals.size:
                sums += vals.sum(axis=0)
                count += vals.shape[0]
        rows.append(pd.Series(sums / count if count else np.nan, index=SAMPLES, name=motif))
        site_summary.append({"output_prefix": motif, "n_sites": count})
    matrix = pd.DataFrame(rows)
    matrix.index.name = "output_prefix"
    matrix.to_csv(out_prefix.with_suffix(".matrix.tsv"), sep="\t")
    pd.DataFrame(site_summary).to_csv(out_prefix.with_suffix(".site_counts.tsv"), sep="\t", index=False)
    return matrix


def build_sum_matrix_from_beds(run_dir: Path, motifs: list[str], out_prefix: Path) -> pd.DataFrame:
    rows = []
    site_summary = []
    for motif in motifs:
        sums = np.zeros(len(SAMPLES), dtype=float)
        count = 0
        for chunk in iter_score_chunks(run_dir, motif):
            vals = chunk.to_numpy(dtype=float)
            finite = np.isfinite(vals).all(axis=1)
            vals = vals[finite]
            if vals.size:
                sums += np.clip(vals, 0.0, None).sum(axis=0)
                count += vals.shape[0]
        rows.append(pd.Series(sums if count else np.nan, index=SAMPLES, name=motif))
        site_summary.append({"output_prefix": motif, "n_sites": count})
    matrix = pd.DataFrame(rows)
    matrix.index.name = "output_prefix"
    matrix.to_csv(out_prefix.with_suffix(".summed_scores.tsv"), sep="\t")
    pd.DataFrame(site_summary).to_csv(out_prefix.with_suffix(".site_counts.tsv"), sep="\t", index=False)
    return matrix


def score_to_pseudocount_matrix(matrix: pd.DataFrame, scale: float = PSEUDOCOUNT_SCALE) -> pd.DataFrame:
    values = np.nan_to_num(matrix.to_numpy(dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    counts = np.rint(np.clip(values, 0.0, None) * float(scale)).astype(np.int64)
    return pd.DataFrame(counts, index=matrix.index, columns=matrix.columns)


def summed_score_to_count_matrix(sum_matrix: pd.DataFrame) -> pd.DataFrame:
    values = np.nan_to_num(sum_matrix.to_numpy(dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    counts = np.rint(np.clip(values, 0.0, None)).astype(np.int64)
    return pd.DataFrame(counts, index=sum_matrix.index, columns=sum_matrix.columns)


def run_pydeseq2(counts_by_motif: pd.DataFrame, out_native: Path, cores: int) -> pd.DataFrame:
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except Exception as exc:
        raise ImportError('PyDESeq2 is required to rebuild Method 3. Install with: pip install "fp-tools-bio[deseq2]"') from exc

    metadata = pd.DataFrame({"condition": ["K562", "K562", "K562", "HepG2", "HepG2", "HepG2"]}, index=SAMPLES)
    dds = DeseqDataSet(
        counts=counts_by_motif.T.loc[SAMPLES].astype(int),
        metadata=metadata,
        design="~condition",
        ref_level=["condition", "HepG2"],
        min_replicates=2,
        n_cpus=max(1, int(cores)),
        quiet=True,
    )
    dds.deseq2()
    stats = DeseqStats(dds, contrast=["condition", "K562", "HepG2"], n_cpus=max(1, int(cores)), quiet=True)
    stats.summary()
    native = stats.results_df.copy()
    native.index.name = "output_prefix"
    native = native.reset_index()
    native.to_csv(out_native, sep="\t", index=False)
    return native


def metadata_for(metadata: pd.DataFrame) -> pd.DataFrame:
    keep = ["output_prefix", "name", "motif_id", "cluster", "total_tfbs", "K562_bound", "HepG2_bound"]
    return metadata[[c for c in keep if c in metadata.columns]].drop_duplicates("output_prefix")


def bound_pass(df: pd.DataFrame) -> pd.Series:
    bounds = df[[c for c in ["K562_bound", "HepG2_bound"] if c in df.columns]].apply(pd.to_numeric, errors="coerce").fillna(0)
    if bounds.empty:
        return pd.Series(False, index=df.index)
    return bounds.max(axis=1) > SIG_BOUND_CUTOFF


def direction(change: pd.Series) -> np.ndarray:
    return np.where(pd.to_numeric(change, errors="coerce").fillna(0.0) >= 0, "K562_up", "HepG2_up")


def apply_method1_significance(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out[f"{COMPARISON}_bonferroni_pvalue"] = bonferroni(out[f"{COMPARISON}_pvalue"])
    change = pd.to_numeric(out[f"{COMPARISON}_change"], errors="coerce").fillna(0.0)
    sig = (
        (out[f"{COMPARISON}_bonferroni_pvalue"] < SIG_P_CUTOFF)
        & (change.abs() > SIG_DELTA_CUTOFF)
        & bound_pass(out)
    )
    out[f"{COMPARISON}_significant_fdr05"] = sig
    out[f"{COMPARISON}_highlighted"] = sig
    return out


def apply_fdr_significance(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    qval = pd.to_numeric(out[f"{COMPARISON}_qvalue_bh"], errors="coerce").fillna(1.0)
    change = pd.to_numeric(out[f"{COMPARISON}_change"], errors="coerce").fillna(0.0)
    sig = (qval < SIG_P_CUTOFF) & (change.abs() > SIG_DELTA_CUTOFF) & bound_pass(out)
    out[f"{COMPARISON}_significant_fdr05"] = sig
    out[f"{COMPARISON}_highlighted"] = sig
    return out


def moderated_score_matrix(matrix: pd.DataFrame, metadata: pd.DataFrame, out_native: Path) -> pd.DataFrame:
    rows = []
    variances = []
    raw = []
    for motif, row in matrix.iterrows():
        x = pd.to_numeric(row[GROUP1], errors="coerce").to_numpy(dtype=float)
        y = pd.to_numeric(row[GROUP2], errors="coerce").to_numpy(dtype=float)
        x = x[np.isfinite(x)]
        y = y[np.isfinite(y)]
        vx = np.var(x, ddof=1) if len(x) > 1 else 0.0
        vy = np.var(y, ddof=1) if len(y) > 1 else 0.0
        pooled = ((len(x) - 1) * vx + (len(y) - 1) * vy) / max(1, len(x) + len(y) - 2) if len(x) and len(y) else np.nan
        raw.append((motif, x, y, pooled))
        if np.isfinite(pooled) and pooled > 0:
            variances.append(float(pooled))
    prior_var = float(np.median(variances)) if variances else 1.0
    if not np.isfinite(prior_var) or prior_var <= 0:
        prior_var = 1.0
    resid_df = 4
    prior_df = 4.0
    total_df = resid_df + prior_df
    for motif, x, y, pooled in raw:
        if len(x) == 0 or len(y) == 0:
            effect, se, tstat, pval, welch_p = np.nan, np.nan, np.nan, 1.0, 1.0
        else:
            effect = float(np.mean(x) - np.mean(y))
            pooled = prior_var if not np.isfinite(pooled) else max(float(pooled), 0.0)
            mod_var = (prior_df * prior_var + resid_df * pooled) / total_df
            se = float(np.sqrt(mod_var * (1.0 / len(x) + 1.0 / len(y))))
            tstat = effect / se if se > 0 else 0.0
            pval = float(2.0 * scipy.stats.t.sf(abs(tstat), total_df))
            welch = scipy.stats.ttest_ind(x, y, equal_var=False, nan_policy="omit")
            welch_p = float(welch.pvalue) if np.isfinite(welch.pvalue) else 1.0
        rows.append(
            {
                "output_prefix": motif,
                "footprint_score_delta": effect,
                "moderated_t": tstat,
                "moderated_se": se,
                "moderated_df": total_df,
                "pvalue": pval,
                "welch_pvalue": welch_p,
                "K562_mean_score": float(np.mean(x)) if len(x) else np.nan,
                "HepG2_mean_score": float(np.mean(y)) if len(y) else np.nan,
            }
        )
    native = pd.DataFrame(rows)
    native["padj"] = bh(native["pvalue"].to_numpy(dtype=float))
    native = metadata_for(metadata).merge(native, on="output_prefix", how="right")
    native.to_csv(out_native, sep="\t", index=False)
    return native


def empirical_bayes_log_matrix(matrix: pd.DataFrame, metadata: pd.DataFrame, out_native: Path) -> pd.DataFrame:
    """Match the original Method 2 effect scale: moderated log2(score * 1000 + 1)."""
    values = np.nan_to_num(matrix.to_numpy(dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    raw_matrix = pd.DataFrame(values, index=matrix.index, columns=matrix.columns)
    log_matrix = pd.DataFrame(
        np.log2(np.clip(values, 0.0, None) * PSEUDOCOUNT_SCALE + 1.0),
        index=matrix.index,
        columns=matrix.columns,
    )
    native = moderated_score_matrix(log_matrix, metadata, out_native)
    native["effect"] = native["footprint_score_delta"]
    native["qvalue"] = native["padj"]
    native["K562_mean_log_score"] = native["K562_mean_score"]
    native["HepG2_mean_log_score"] = native["HepG2_mean_score"]
    native["K562_mean_raw_score"] = native["output_prefix"].map(raw_matrix[GROUP1].mean(axis=1))
    native["HepG2_mean_raw_score"] = native["output_prefix"].map(raw_matrix[GROUP2].mean(axis=1))
    native["raw_score_delta"] = native["K562_mean_raw_score"] - native["HepG2_mean_raw_score"]
    native.to_csv(out_native, sep="\t", index=False)
    return native


def eb_to_diff(native: pd.DataFrame, out_path: Path, method_label: str) -> pd.DataFrame:
    out = metadata_for(native).copy()
    data = native.set_index("output_prefix")
    out[f"{COMPARISON}_change"] = out["output_prefix"].map(data["footprint_score_delta"]).astype(float)
    out[f"{COMPARISON}_pvalue"] = out["output_prefix"].map(data["pvalue"]).astype(float).clip(lower=1e-308)
    out[f"{COMPARISON}_qvalue_bh"] = out["output_prefix"].map(data["padj"]).astype(float).fillna(1.0)
    out[f"{COMPARISON}_mean_delta_fp"] = out[f"{COMPARISON}_change"]
    out[f"{COMPARISON}_mean_log2fc"] = out[f"{COMPARISON}_change"]
    if "raw_score_delta" in data.columns:
        out[f"{COMPARISON}_raw_score_delta"] = out["output_prefix"].map(data["raw_score_delta"]).astype(float)
    out["method_label"] = method_label
    out = apply_fdr_significance(out)
    out.to_csv(out_path, sep="\t", index=False)
    return out


def deseq_to_diff(native: pd.DataFrame, raw_matrix: pd.DataFrame, metadata: pd.DataFrame, out_path: Path, method_label: str) -> pd.DataFrame:
    means = pd.DataFrame(
        {
            "output_prefix": raw_matrix.index,
            "K562_mean_score": raw_matrix[GROUP1].mean(axis=1).to_numpy(dtype=float),
            "HepG2_mean_score": raw_matrix[GROUP2].mean(axis=1).to_numpy(dtype=float),
        }
    )
    means["score_delta"] = means["K562_mean_score"] - means["HepG2_mean_score"]
    out = metadata_for(metadata).merge(native, on="output_prefix", how="right").merge(means, on="output_prefix", how="left")
    out = out.rename(
        columns={
            "log2FoldChange": "deseq2_log2FoldChange",
            "pvalue": "deseq2_pvalue",
            "padj": "deseq2_padj",
        }
    )
    out[f"{COMPARISON}_change"] = pd.to_numeric(out["deseq2_log2FoldChange"], errors="coerce").fillna(0.0)
    out[f"{COMPARISON}_pvalue"] = pd.to_numeric(out["deseq2_pvalue"], errors="coerce").fillna(1.0).clip(lower=1e-308)
    out[f"{COMPARISON}_qvalue_bh"] = pd.to_numeric(out["deseq2_padj"], errors="coerce").fillna(1.0)
    out[f"{COMPARISON}_raw_score_delta"] = pd.to_numeric(out["score_delta"], errors="coerce").fillna(0.0)
    out[f"{COMPARISON}_mean_delta_fp"] = out[f"{COMPARISON}_raw_score_delta"]
    out[f"{COMPARISON}_mean_log2fc"] = out[f"{COMPARISON}_change"]
    out["method_label"] = method_label
    out = apply_fdr_significance(out)
    out.to_csv(out_path, sep="\t", index=False)
    return out


def native_bindetect_to_diff(native: pd.DataFrame, audit_native: pd.DataFrame, out_path: Path, method_label: str) -> pd.DataFrame:
    out = metadata_for(native).copy()
    data = native.set_index("output_prefix")
    out[f"{COMPARISON}_change"] = out["output_prefix"].map(data[f"{COMPARISON}_change"]).astype(float)
    out[f"{COMPARISON}_pvalue"] = out["output_prefix"].map(data[f"{COMPARISON}_pvalue"]).astype(float).clip(lower=1e-308)
    out[f"{COMPARISON}_qvalue_bh"] = out["output_prefix"].map(data[f"{COMPARISON}_qvalue_bh"]).astype(float).fillna(1.0)
    for col in [f"{COMPARISON}_mean_delta_fp", f"{COMPARISON}_mean_log2fc", f"{COMPARISON}_delta_fp_se", f"{COMPARISON}_log2fc_se"]:
        if col in data.columns:
            out[col] = out["output_prefix"].map(data[col]).astype(float)
    if audit_native is not None and not audit_native.empty:
        audit = audit_native.set_index("output_prefix")
        audit_cols = {
            "footprint_score_delta": f"{COMPARISON}_matrix_log_score_delta",
            "raw_score_delta": f"{COMPARISON}_raw_score_delta",
            "pvalue": f"{COMPARISON}_matrix_log_score_pvalue",
            "padj": f"{COMPARISON}_matrix_log_score_qvalue_bh",
            "moderated_se": f"{COMPARISON}_matrix_log_score_se",
        }
        for source, dest in audit_cols.items():
            if source in audit.columns:
                out[dest] = out["output_prefix"].map(audit[source]).astype(float)
    out["method_label"] = method_label
    out = apply_fdr_significance(out)
    out.to_csv(out_path, sep="\t", index=False)
    return out


def read_or_build_method2(norm_sample: Path, motifs: list[str], metadata: pd.DataFrame, outdir: Path) -> pd.DataFrame:
    native_path = outdir / "02_method2_sample_quantile_scores_limma_native.tsv"
    result_path = outdir / "02_method2_sample_quantile_scores_limma_results.tsv"
    matrix = build_matrix_from_beds(norm_sample, motifs, outdir / "02_method2_sample_quantile_scores")
    native = empirical_bayes_log_matrix(matrix, metadata, native_path)
    return eb_to_diff(native, result_path, "Sample-quantile motif footprint-score matrix, empirical-Bayes log-score test")


def read_or_build_method3(norm_none: Path, motifs: list[str], metadata: pd.DataFrame, fp: Path, outdir: Path, cores: int) -> pd.DataFrame:
    native_path = outdir / "03_method3_raw_score_pseudocount_deseq2_native.tsv"
    result_path = outdir / "03_method3_raw_score_pseudocount_deseq2_results.tsv"
    raw_matrix = build_matrix_from_beds(norm_none, motifs, outdir / "03_method3_raw_scores")
    summed_scores = build_sum_matrix_from_beds(norm_none, motifs, outdir / "03_method3_raw_score_summed")
    pseudo = summed_score_to_count_matrix(summed_scores)
    pseudo.to_csv(outdir / "03_method3_raw_score_summed_count_matrix.tsv", sep="\t")
    native = run_pydeseq2(pseudo, native_path, cores)
    return deseq_to_diff(native, raw_matrix, metadata, result_path, "Raw summed motif footprint-score DESeq2")


def read_or_build_method4(norm_q95: Path, motifs: list[str], metadata: pd.DataFrame, outdir: Path) -> pd.DataFrame:
    native_path = outdir / "04_method4_q95_scores_limma_native.tsv"
    result_path = outdir / "04_method4_q95_scores_limma_results.tsv"
    matrix = build_matrix_from_beds(norm_q95, motifs, outdir / "04_method4_q95_scores")
    audit_native = empirical_bayes_log_matrix(matrix, metadata, native_path)
    return native_bindetect_to_diff(metadata, audit_native, result_path, "Q95 corrected BINDetect results, no internal normalization")


def unified_summary(methods: dict[str, pd.DataFrame], out_path: Path) -> pd.DataFrame:
    keys = ["output_prefix", "name", "motif_id", "cluster", "total_tfbs"]
    merged = None
    for method, df in methods.items():
        extra = [f"{COMPARISON}_bonferroni_pvalue"] if f"{COMPARISON}_bonferroni_pvalue" in df.columns else []
        keep = df[keys + [c for c in ["K562_bound", "HepG2_bound"] if c in df.columns] + [
            f"{COMPARISON}_change",
            f"{COMPARISON}_pvalue",
            f"{COMPARISON}_qvalue_bh",
            f"{COMPARISON}_significant_fdr05",
            *extra,
        ]].copy()
        keep[f"{method}_direction"] = direction(keep[f"{COMPARISON}_change"])
        rename = {
            f"{COMPARISON}_change": f"{method}_delta_score",
            f"{COMPARISON}_pvalue": f"{method}_pvalue",
            f"{COMPARISON}_qvalue_bh": f"{method}_qvalue_bh",
            f"{COMPARISON}_significant_fdr05": f"{method}_significant",
            "K562_bound": f"{method}_K562_bound",
            "HepG2_bound": f"{method}_HepG2_bound",
            f"{COMPARISON}_bonferroni_pvalue": f"{method}_bonferroni_pvalue",
        }
        keep = keep.rename(columns=rename)
        if merged is None:
            merged = keep
        else:
            merged = merged.merge(keep.drop(columns=[c for c in keys[1:] if c in keep.columns]), on="output_prefix", how="outer")
    assert merged is not None
    sig_cols = [f"{method}_significant" for method in methods]
    for col in sig_cols:
        merged[col] = merged[col].fillna(False).astype(bool)
    merged["n_methods_significant"] = merged[sig_cols].sum(axis=1)
    merged["methods_significant"] = merged.apply(lambda r: ";".join([m for m in methods if bool(r.get(f"{m}_significant", False))]), axis=1)
    dir_cols = [f"{method}_direction" for method in methods]
    merged["direction_agreement"] = merged.apply(
        lambda r: len({r[c] for c in dir_cols if pd.notna(r.get(c)) and bool(r.get(c.replace("_direction", "_significant"), False))}) <= 1,
        axis=1,
    )
    merged.to_csv(out_path, index=False)
    return merged


def write_pairwise_summaries(summary: pd.DataFrame, methods: list[str], outdir: Path) -> None:
    for i, left in enumerate(methods):
        for right in methods[i + 1:]:
            left_sig = summary[f"{left}_significant"].fillna(False).astype(bool)
            right_sig = summary[f"{right}_significant"].fillna(False).astype(bool)
            prefix = f"{left}_vs_{right}"
            subset = summary.copy()
            subset.to_csv(outdir / f"{prefix}_all_motifs.csv", index=False)
            subset.loc[left_sig & ~right_sig].to_csv(outdir / f"{prefix}_{left}_only_significant.csv", index=False)
            subset.loc[~left_sig & right_sig].to_csv(outdir / f"{prefix}_{right}_only_significant.csv", index=False)
            subset.loc[left_sig & right_sig].to_csv(outdir / f"{prefix}_shared_significant.csv", index=False)


def find_signal_files(fp_dir: Path, subdir: str, pattern_suffix: str) -> list[str]:
    root = fp_dir / subdir
    paths = []
    for sample in SAMPLES:
        matches = sorted(root.rglob(f"{sample}*{pattern_suffix}"))
        if not matches:
            raise FileNotFoundError(f"Could not find {sample}*{pattern_suffix} under {root}")
        paths.append(str(matches[0]))
    return paths


def make_aggregate_args(source_dir: Path, aggregate_signals: list[str], cores: int, aggregate_max_sites: int | None = None) -> SimpleNamespace:
    return SimpleNamespace(
        outdir=str(source_dir),
        signals=list(aggregate_signals),
        aggregate_signals=list(aggregate_signals),
        cond_groups={"K562": [0, 1, 2], "HepG2": [3, 4, 5]},
        sample_names=list(SAMPLES),
        normalization="none",
        aggregate_normalization="none",
        aggregate_site_set="all",
        plot_aggregate="sig",
        plot_aggregate_top_n=AGGREGATE_TOP_N,
        aggregate_pvalue_threshold=SIG_P_CUTOFF,
        aggregate_sig_only=True,
        aggregate_sig_no_fallback=True,
        aggregate_flank=100,
        aggregate_max_sites=aggregate_max_sites,
        cores=max(1, int(cores)),
    )


def make_logo_args(source_dir: Path) -> SimpleNamespace:
    return SimpleNamespace(outdir=str(source_dir))


def write_full_html_report(
    result_table: pd.DataFrame,
    source_dir: Path,
    aggregate_signals: list[str],
    html_out: Path,
    report_label: str,
    cores: int,
    change_label: str = "Delta differential footprint score",
    aggregate_max_sites: int | None = None,
) -> None:
    report_table = result_table.copy()
    sig = report_table[f"{COMPARISON}_significant_fdr05"].astype(bool)
    report_table[f"{COMPARISON}_highlighted"] = sig
    logo_args = make_logo_args(source_dir)
    motifs = _existing_result_motifs(report_table, COMPARISON_TUPLE, logo_args)
    aggregate_args = make_aggregate_args(source_dir, aggregate_signals, cores, aggregate_max_sites=aggregate_max_sites)
    aggregate_data = build_bindetect_aggregate_payload(motifs, report_table, COMPARISON_TUPLE, aggregate_args)
    plot_interactive_bindetect(
        motifs,
        list(COMPARISON_TUPLE),
        str(html_out),
        aggregate_data=aggregate_data,
        title="Differential footprint report",
        report_label=report_label,
        change_label=change_label,
    )


def organize_review_folder(outdir: Path) -> None:
    csv_dir = outdir / "csv"
    plots_dir = outdir / "plots"
    csv_dir.mkdir(exist_ok=True)
    plots_dir.mkdir(exist_ok=True)
    final_root = {
        "00_all_methods_motif_comparison.csv",
        "01_method1_tobias_sample_quantile.html",
        "02_method2_qnorm_limma.html",
        "03_method3_raw_score_deseq2.html",
        "04_method4_q95_limma.html",
    }
    for path in sorted(outdir.iterdir()):
        if path.is_dir():
            continue
        if path.name in final_root:
            continue
        if path.suffix == ".html":
            path.unlink()
        elif path.suffix == ".pdf":
            shutil.move(str(path), str(plots_dir / path.name))
        elif path.suffix in {".csv", ".tsv", ".gz"} or path.name.endswith(".tsv.gz"):
            shutil.move(str(path), str(csv_dir / path.name))
        else:
            path.unlink()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fp-dir", type=Path, default=Path("data/public/processed/encode_k562_hepg2_atac_replicates/fp_tools"))
    ap.add_argument("--outdir", type=Path, default=None)
    ap.add_argument("--html-cores", type=int, default=4, help="Cores to use while rebuilding aggregate payloads for HTML reports.")
    args = ap.parse_args()

    fp = args.fp_dir
    outdir = args.outdir or fp / "review_normalization_comparison_20260618"
    outdir.mkdir(parents=True, exist_ok=True)

    norm_sample = fp / "diff_footprints_jaspar2026_vertebrates_norm_sample_quantile"
    norm_none = fp / "diff_footprints_jaspar2026_vertebrates_norm_none"
    norm_q95 = fp / "diff_footprints_jaspar2026_vertebrates_norm_corrected_q95"
    metadata_sample = read_results(norm_sample / "diff_footprints_results.txt")
    metadata_none = read_results(norm_none / "diff_footprints_results.txt")
    metadata_q95 = read_results(norm_q95 / "diff_footprints_results.txt")
    motifs = metadata_sample["output_prefix"].astype(str).tolist()
    raw_corrected_signals = find_signal_files(fp, "atac_correct", "_corrected.bw")
    q95_corrected_signals = find_signal_files(fp, "normalized_corrected_bigwigs/peak_q95", "_corrected.background_scale_q95.bw")

    method1 = apply_method1_significance(metadata_sample.copy())
    method1.to_csv(outdir / "01_method1_tobias_sample_quantile_results.tsv", sep="\t", index=False)
    write_full_html_report(
        method1,
        norm_sample,
        q95_corrected_signals,
        outdir / "01_method1_tobias_sample_quantile.html",
        "Method 1: TOBIAS/BINDetect sample-quantile results; significant = Bonferroni p < 0.001, |delta score| > 0.1, and max bound sites > 500.",
        args.html_cores,
        aggregate_max_sites=AGGREGATE_TOP_N,
    )

    method2 = read_or_build_method2(norm_sample, motifs, metadata_sample, outdir)
    write_full_html_report(
        method2,
        norm_sample,
        q95_corrected_signals,
        outdir / "02_method2_qnorm_limma.html",
        "Method 2: sample-quantile motif footprint-score matrix with empirical-Bayes test on log2(score * 1000 + 1); significant = FDR < 0.001, |log-score effect| > 0.1, and max bound sites > 500.",
        args.html_cores,
        change_label="Mean log2 footprint-score difference (K562 - HepG2)",
        aggregate_max_sites=AGGREGATE_TOP_N,
    )

    method3 = read_or_build_method3(norm_none, motifs, metadata_none, fp, outdir, args.html_cores)
    write_full_html_report(
        method3,
        norm_none,
        raw_corrected_signals,
        outdir / "03_method3_raw_score_deseq2.html",
        "Method 3: raw summed motif footprint-score counts with default DESeq2 size-factor normalization; significant = FDR < 0.001, |DESeq2 log2 fold change| > 0.1, and max bound sites > 500.",
        args.html_cores,
        change_label="DESeq2 log2 fold change (K562 / HepG2)",
        aggregate_max_sites=AGGREGATE_TOP_N,
    )

    method4 = read_or_build_method4(norm_q95, motifs, metadata_q95, outdir)
    write_full_html_report(
        method4,
        norm_q95,
        q95_corrected_signals,
        outdir / "04_method4_q95_limma.html",
        "Method 4: Q95-corrected BINDetect results with no internal normalization; significant = FDR < 0.001, |native BINDetect change| > 0.1, and max bound sites > 500.",
        args.html_cores,
        change_label="Q95 BINDetect differential footprint score",
        aggregate_max_sites=AGGREGATE_TOP_N,
    )

    methods = {
        "method1_tobias_qnorm": method1,
        "method2_qnorm_limma": method2,
        "method3_raw_score_deseq2": method3,
        "method4_q95_limma": method4,
    }
    summary = unified_summary(methods, outdir / "00_all_methods_motif_comparison.csv")
    write_pairwise_summaries(summary, list(methods), outdir)
    pd.DataFrame(
        [{"method": name, "significant": int(df[f"{COMPARISON}_significant_fdr05"].astype(bool).sum())} for name, df in methods.items()]
    ).to_csv(outdir / "00_significance_summary.csv", index=False)
    organize_review_folder(outdir)
    print(f"Wrote {outdir}")
    print(summary[["output_prefix", "n_methods_significant", "methods_significant"]].head().to_string(index=False))


if __name__ == "__main__":
    main()
