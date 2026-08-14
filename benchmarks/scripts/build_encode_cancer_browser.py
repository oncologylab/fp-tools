#!/usr/bin/env python3
"""Run, validate, and export the seven-line ENCODE cancer footprint resource."""

from __future__ import annotations

import argparse
import concurrent.futures
import gzip
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import subprocess

import numpy as np
import pandas as pd
import pyBigWig

from fp_tools.utils.empirical_bayes import fit_moderated_contrast
from fp_tools.utils.project_layout import normalized_bigwig_path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_PROJECT = ROOT / "data/public/processed/encode_cancer_7line_20260814"
DEFAULT_MANIFEST = ROOT / "benchmarks/manifests/encode_cancer_7line_20260814.tsv"
DEFAULT_SPEC = ROOT / "benchmarks/manifests/encode_cancer_7line_20260814.spec.json"
DEFAULT_COMPARISONS = ROOT / "benchmarks/manifests/encode_cancer_7line_20260814_comparisons.tsv"
DEFAULT_GENOME = ROOT / "data/public/raw/genome/hg38.fa"
DEFAULT_SITE = ROOT / "docs/ENCODE-Cancer-Cell-lines-Footprinting"
EXPECTED_MOTIFS = 1019
EXPECTED_COMPARISONS = 21
RELEASE_DATE = "2026-08-14"


def load_spec(path: Path) -> dict:
    payload = json.loads(path.read_text(encoding="utf-8"))
    conditions = payload.get("conditions") or {}
    if len(conditions) != 7:
        raise ValueError(f"Expected seven conditions, found {len(conditions)}")
    samples = [sample for values in conditions.values() for sample in values]
    if len(samples) != 15 or len(samples) != len(set(samples)):
        raise ValueError("Expected exactly 15 unique biological-replicate samples")
    if int(payload.get("motifs_per_sample", 0)) != EXPECTED_MOTIFS:
        raise ValueError(f"Expected {EXPECTED_MOTIFS} motifs per sample")
    peak = payload.get("peak_universe")
    if not peak or set(peak) != {"regions", "covered_bp", "md5"}:
        raise ValueError("The project specification must contain a locked peak universe")
    return payload


def read_design(manifest_path: Path, spec_path: Path, comparisons_path: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    spec = load_spec(spec_path)
    manifest = pd.read_csv(manifest_path, sep="\t", dtype=str, keep_default_na=False)
    comparisons = pd.read_csv(comparisons_path, sep="\t", dtype=str, keep_default_na=False)
    expected_samples = {sample for values in spec["conditions"].values() for sample in values}
    if len(manifest) != 15 or set(manifest["sample"]) != expected_samples:
        raise ValueError("Manifest does not contain the exact 15 samples in the project specification")
    if set(manifest["condition"]) != set(spec["conditions"]):
        raise ValueError("Manifest conditions do not match the project specification")
    if set(manifest["condition"]).intersection({"GM12878", "IMR-90", "DND-41"}):
        raise ValueError("Reference or excluded cell lines are present in the cancer-only design")
    if set(manifest.loc[manifest["condition"].eq("A549"), "peak_accession"]) != {"ENCFF876UEM"}:
        raise ValueError("A549 must use the conservative-IDR peak file ENCFF876UEM")
    expected_pairs = {frozenset(pair) for pair in itertools.combinations(sorted(spec["conditions"]), 2)}
    observed_pairs = {frozenset((row.cond1, row.cond2)) for row in comparisons.itertuples(index=False)}
    if len(comparisons) != EXPECTED_COMPARISONS or observed_pairs != expected_pairs:
        raise ValueError("Comparison table does not contain every unordered seven-line pair exactly once")
    if comparisons["comparison"].duplicated().any():
        raise ValueError("Comparison IDs must be unique")
    return manifest, comparisons, spec


def _result_path(project: Path, comparison: str | None = None) -> Path:
    return project / "comparisons" / "all_21_pairwise" / "diff_footprints_results.txt"


def _comparison_columns(cond1: str, cond2: str) -> dict[str, str]:
    base = f"{cond1}_{cond2}_ebayes"
    return {
        "mean1": f"{cond1}_mean_score",
        "sd1": f"{cond1}_score_sd",
        "mean2": f"{cond2}_mean_score",
        "sd2": f"{cond2}_score_sd",
        "effect": f"{base}_effect",
        "ci_lower": f"{base}_ci_lower",
        "ci_upper": f"{base}_ci_upper",
        "moderated_t": f"{base}_moderated_t",
        "moderated_df": f"{base}_moderated_df",
        "pvalue": f"{base}_pvalue",
        "qvalue": f"{base}_qvalue_bh",
        "significant": f"{base}_significant_fdr05",
    }


def validate_result(path: Path, cond1: str, cond2: str) -> pd.DataFrame:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Missing comparison result: {path}")
    table = pd.read_csv(path, sep="\t")
    columns = _comparison_columns(cond1, cond2)
    required = {"output_prefix", "name", "motif_id", "cluster", "total_tfbs", *columns.values()}
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(f"{path} is missing empirical-Bayes columns: {sorted(missing)}")
    if len(table) != EXPECTED_MOTIFS or table["output_prefix"].nunique() != EXPECTED_MOTIFS:
        raise ValueError(f"{path} does not contain exactly {EXPECTED_MOTIFS} motifs")
    for label in ("effect", "ci_lower", "ci_upper", "moderated_t", "pvalue", "qvalue"):
        values = pd.to_numeric(table[columns[label]], errors="coerce")
        if values.isna().any() or not np.isfinite(values).all():
            raise ValueError(f"{path} contains non-finite {label} values")
    moderated_df = pd.to_numeric(table[columns["moderated_df"]], errors="coerce").to_numpy(dtype=float)
    if np.isnan(moderated_df).any() or np.isneginf(moderated_df).any() or (moderated_df <= 0).any():
        raise ValueError(f"{path} contains invalid moderated_df values")
    pvalues = pd.to_numeric(table[columns["pvalue"]], errors="raise")
    qvalues = pd.to_numeric(table[columns["qvalue"]], errors="raise")
    if ((pvalues < 0) | (pvalues > 1) | (qvalues < 0) | (qvalues > 1)).any():
        raise ValueError(f"{path} contains p- or q-values outside [0,1]")
    return table


def recover_result_from_replicate_matrix(
    project: Path,
    manifest: pd.DataFrame,
    comparisons: pd.DataFrame,
) -> Path:
    """Rebuild the compact comparison table from completed motif summaries.

    ``diff-footprints`` writes the per-sample motif-score matrix before its
    presentation outputs. That matrix is the sufficient input for every
    replicate empirical-Bayes contrast in this resource, so a presentation-
    table validation failure must not require repeating the site-level run.
    """

    output = _result_path(project)
    matrix_path = output.with_name("diff_footprints_replicate_motif_score_matrix.tsv")
    if not matrix_path.is_file() or matrix_path.stat().st_size == 0:
        raise ValueError("The completed replicate motif-score matrix is unavailable")
    ordered = manifest.sort_values(["condition", "biological_replicate"]).reset_index(drop=True)
    samples = ordered["sample"].tolist()
    score_table = pd.read_csv(matrix_path, sep="\t")
    if len(score_table) != EXPECTED_MOTIFS or set(samples).difference(score_table.columns):
        raise ValueError("The replicate motif-score matrix is incomplete")
    if score_table["motif"].duplicated().any():
        raise ValueError("The replicate motif-score matrix contains duplicate motifs")
    scores = score_table.set_index("motif")[samples].apply(pd.to_numeric, errors="coerce")
    if scores.isna().any().any() or not np.isfinite(scores.to_numpy(dtype=float)).all():
        raise ValueError("The replicate motif-score matrix contains non-finite scores")

    first_sample = samples[0]
    summary_path = project / "samples" / first_sample / "match_motifs" / "motif_matches_results.txt"
    summary = pd.read_csv(
        summary_path,
        sep="\t",
        usecols=["output_prefix", "name", "motif_id", "cluster", "total_tfbs"],
    )
    if len(summary) != EXPECTED_MOTIFS or summary["output_prefix"].duplicated().any():
        raise ValueError("The motif annotation summary is incomplete")
    summary = summary.set_index("output_prefix").loc[scores.index]
    summary.index.name = "output_prefix"
    summary = summary.reset_index()
    n_sites = pd.to_numeric(score_table.set_index("motif").loc[scores.index, "n_sites"], errors="raise")
    total_tfbs = pd.to_numeric(summary["total_tfbs"], errors="raise")
    if not np.array_equal(n_sites.to_numpy(dtype=int), total_tfbs.to_numpy(dtype=int)):
        raise ValueError("Motif-site counts differ between completed summaries")

    derived_columns: dict[str, np.ndarray] = {}
    sample_conditions = dict(zip(ordered["sample"], ordered["condition"], strict=True))
    for condition, condition_samples in ordered.groupby("condition", sort=False)["sample"]:
        values = scores[condition_samples.tolist()]
        derived_columns[f"{condition}_mean_score"] = values.mean(axis=1).to_numpy()
        derived_columns[f"{condition}_score_sd"] = values.std(axis=1, ddof=1).to_numpy()
    for comparison in comparisons.itertuples(index=False):
        fitted = fit_moderated_contrast(
            scores,
            sample_conditions,
            comparison.cond1,
            comparison.cond2,
        )
        base = f"{comparison.cond1}_{comparison.cond2}_ebayes"
        for source, destination in (
            ("effect", f"{base}_effect"),
            ("ci_lower", f"{base}_ci_lower"),
            ("ci_upper", f"{base}_ci_upper"),
            ("moderated_t", f"{base}_moderated_t"),
            ("moderated_df", f"{base}_moderated_df"),
            ("pvalue", f"{base}_pvalue"),
            ("qvalue_bh", f"{base}_qvalue_bh"),
            ("significant_fdr05", f"{base}_significant_fdr05"),
        ):
            derived_columns[destination] = fitted.loc[scores.index, source].to_numpy()

    result = pd.concat(
        [summary.reset_index(drop=True), pd.DataFrame(derived_columns)],
        axis=1,
    )

    temporary = output.with_name(output.name + ".part")
    result.to_csv(temporary, sep="\t", index=False, float_format="%.10g")
    temporary.replace(output)
    for comparison in comparisons.itertuples(index=False):
        validate_result(output, comparison.cond1, comparison.cond2)
    return output


def run_comparisons(project: Path, manifest: pd.DataFrame, comparisons: pd.DataFrame, genome: Path, cores: int) -> None:
    peaks = project / "peaks" / "merged_peaks_filtered.bed"
    if not peaks.is_file() or not genome.is_file():
        raise ValueError("Verified project peaks and GRCh38 FASTA are required")
    output = _result_path(project)
    if output.is_file():
        try:
            for comparison in comparisons.itertuples(index=False):
                validate_result(output, comparison.cond1, comparison.cond2)
            return
        except ValueError:
            pass
    matrix = output.with_name("diff_footprints_replicate_motif_score_matrix.tsv")
    if matrix.is_file() and matrix.stat().st_size > 0:
        recover_result_from_replicate_matrix(project, manifest, comparisons)
        return
    ordered = manifest.sort_values(["condition", "biological_replicate"]).reset_index(drop=True)
    selected = ordered["sample"].tolist()
    conditions = ordered["condition"].tolist()
    outdir = output.parent
    outdir.mkdir(parents=True, exist_ok=True)
    command = [
        str(ROOT / ".venv/bin/diff-footprints"),
        "--layout", "custom",
        "--sample-dirs", *[str(project / "samples" / sample) for sample in selected],
        "--sample-names", *selected,
        "--cond-names", *conditions,
        "--peaks", str(peaks),
        "--genome", str(genome),
        "--motif-db", "jaspar2026_vertebrates",
        "--outdir", str(outdir),
        "--prefix", "diff_footprints",
        "--normalization", "sample-quantile",
        "--replicate-report", "off",
        "--plot-aggregate", "off",
        "--motif-outputs", "summary",
        "--skip-excel",
        "--cores", str(cores),
        "--verbosity", "3",
    ]
    log = project / "logs" / "all_21_pairwise.diff_footprints.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("a", encoding="utf-8") as handle:
        handle.write("\n$ " + " ".join(command) + "\n")
        handle.flush()
        result = subprocess.run(command, cwd=ROOT, stdout=handle, stderr=subprocess.STDOUT, check=False)
    if result.returncode:
        raise RuntimeError(f"Pairwise analysis failed; see {log}")
    for comparison in comparisons.itertuples(index=False):
        validate_result(output, comparison.cond1, comparison.cond2)


def standardized_comparisons(project: Path, comparisons: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    frames = []
    summaries = []
    for comparison in comparisons.itertuples(index=False):
        table = validate_result(_result_path(project, comparison.comparison), comparison.cond1, comparison.cond2)
        columns = _comparison_columns(comparison.cond1, comparison.cond2)
        out = pd.DataFrame({
            "comparison": comparison.comparison,
            "condition1": comparison.cond1,
            "condition2": comparison.cond2,
            "prefix": table["output_prefix"].astype(str),
            "name": table["name"].astype(str),
            "motif_id": table["motif_id"].astype(str),
            "cluster": table["cluster"].astype(str),
            "n_sites": pd.to_numeric(table["total_tfbs"], errors="raise").astype(int),
            "mean1": pd.to_numeric(table[columns["mean1"]], errors="raise"),
            "sd1": pd.to_numeric(table[columns["sd1"]], errors="raise"),
            "mean2": pd.to_numeric(table[columns["mean2"]], errors="raise"),
            "sd2": pd.to_numeric(table[columns["sd2"]], errors="raise"),
            "effect": pd.to_numeric(table[columns["effect"]], errors="raise"),
            "ci_lower": pd.to_numeric(table[columns["ci_lower"]], errors="raise"),
            "ci_upper": pd.to_numeric(table[columns["ci_upper"]], errors="raise"),
            "moderated_t": pd.to_numeric(table[columns["moderated_t"]], errors="raise"),
            "moderated_df": pd.to_numeric(table[columns["moderated_df"]], errors="raise"),
            "pvalue": pd.to_numeric(table[columns["pvalue"]], errors="raise"),
            "qvalue": pd.to_numeric(table[columns["qvalue"]], errors="raise"),
            "significant": table[columns["significant"]].astype(str).str.lower().eq("true"),
        })
        frames.append(out)
        positive = int((out["significant"] & out["effect"].gt(0)).sum())
        negative = int((out["significant"] & out["effect"].lt(0)).sum())
        summaries.append({
            "comparison": comparison.comparison,
            "condition1": comparison.cond1,
            "condition2": comparison.cond2,
            "motifs_tested": len(out),
            "condition1_enriched_q05": positive,
            "condition2_enriched_q05": negative,
            "significant_total_q05": positive + negative,
        })
    combined = pd.concat(frames, ignore_index=True)
    if len(combined) != EXPECTED_MOTIFS * EXPECTED_COMPARISONS:
        raise ValueError("Combined result does not contain 21,399 motif-comparison rows")
    return combined, pd.DataFrame(summaries)


def coordinate_digest(path: Path) -> str:
    value = hashlib.sha256()
    columns = ["motif", "TFBS_chr", "TFBS_start", "TFBS_end", "TFBS_strand"]
    for chunk in pd.read_csv(path, sep="\t", usecols=columns, chunksize=100_000):
        text = chunk.to_csv(sep="\t", index=False, header=False, lineterminator="\n")
        value.update(text.encode("utf-8"))
    return value.hexdigest()


def _select_profile_sites(sites: pd.DataFrame, max_sites_per_motif: int) -> tuple[pd.DataFrame, dict[str, int]]:
    counts = sites.groupby("motif", sort=False, observed=True).size().astype(int).to_dict()
    selected = []
    for _motif, group in sites.groupby("motif", sort=False, observed=True):
        if max_sites_per_motif > 0 and len(group) > max_sites_per_motif:
            indices = np.linspace(0, len(group) - 1, max_sites_per_motif, dtype=int)
            group = group.iloc[np.unique(indices)]
        selected.append(group)
    return pd.concat(selected, ignore_index=True), {str(key): int(value) for key, value in counts.items()}


def prepare_shared_profile_sites(
    project: Path,
    manifest: pd.DataFrame,
    max_sites_per_motif: int,
) -> tuple[Path, dict[str, int]]:
    """Stream one shared-scan cache into a bounded coordinate selection.

    The motif coordinates are identical across samples because match-motifs ran
    in shared-scan mode. Reading one coordinate cache once avoids loading or
    rescanning the 50-million-row cache independently for every profile worker.
    """

    first_sample = str(manifest.iloc[0]["sample"])
    source = project / "samples" / first_sample / "match_motifs" / "cache" / "motif_sites.tsv.gz"
    summary_path = project / "samples" / first_sample / "match_motifs" / "motif_matches_results.txt"
    upstream_marker = project / "state" / "match_motifs.json"
    if not source.is_file() or not summary_path.is_file() or not upstream_marker.is_file():
        raise ValueError("Verified shared motif cache, summary, and marker are required")

    summary = pd.read_csv(summary_path, sep="\t", usecols=["output_prefix", "total_tfbs"])
    total_counts = {
        str(row.output_prefix): int(row.total_tfbs)
        for row in summary.itertuples(index=False)
    }
    if len(total_counts) != EXPECTED_MOTIFS or any(value <= 0 for value in total_counts.values()):
        raise ValueError("Motif summary does not contain 1,019 positive site counts")
    expected_selected = sum(min(value, max_sites_per_motif) for value in total_counts.values())
    selection_dir = project / "reports" / "browser"
    selection_dir.mkdir(parents=True, exist_ok=True)
    destination = selection_dir / f"profile_sites_max{max_sites_per_motif}.tsv.gz"
    marker = destination.with_suffix(destination.suffix + ".json")
    signature = {
        "schema": "fp-tools.encode-cancer-profile-sites.v1",
        "source_sample": first_sample,
        "source_size": source.stat().st_size,
        "source_mtime_ns": source.stat().st_mtime_ns,
        "upstream_marker_sha256": hashlib.sha256(upstream_marker.read_bytes()).hexdigest(),
        "max_sites_per_motif": max_sites_per_motif,
        "motifs": len(total_counts),
        "source_rows": sum(total_counts.values()),
        "selected_rows": expected_selected,
    }
    if destination.is_file() and marker.is_file():
        cached = json.loads(marker.read_text(encoding="utf-8"))
        if (
            all(cached.get(key) == value for key, value in signature.items())
            and cached.get("selection_sha256") == hashlib.sha256(destination.read_bytes()).hexdigest()
        ):
            return destination, total_counts

    targets = {
        motif: np.unique(np.linspace(0, count - 1, min(count, max_sites_per_motif), dtype=np.int64))
        for motif, count in total_counts.items()
    }
    seen = {motif: 0 for motif in total_counts}
    written = 0
    columns = ["motif", "TFBS_chr", "TFBS_start", "TFBS_end", "TFBS_strand"]
    temporary = destination.with_name(destination.name + ".part")
    with gzip.open(temporary, "wt", encoding="utf-8", compresslevel=6) as handle:
        write_header = True
        for chunk in pd.read_csv(source, sep="\t", usecols=columns, chunksize=1_000_000):
            selected_groups = []
            for motif, group in chunk.groupby("motif", sort=False, observed=True):
                motif = str(motif)
                if motif not in targets:
                    raise ValueError(f"Unexpected motif in shared coordinate cache: {motif}")
                start = seen[motif]
                stop = start + len(group)
                target = targets[motif]
                local = target[(target >= start) & (target < stop)] - start
                if len(local):
                    selected_groups.append(group.iloc[local])
                seen[motif] = stop
            if selected_groups:
                selected = pd.concat(selected_groups, ignore_index=True)
                selected.to_csv(handle, sep="\t", index=False, header=write_header, lineterminator="\n")
                write_header = False
                written += len(selected)
    if seen != total_counts:
        mismatched = {motif: (seen[motif], total_counts[motif]) for motif in total_counts if seen[motif] != total_counts[motif]}
        temporary.unlink(missing_ok=True)
        raise ValueError(f"Shared coordinate cache counts differ from motif summary: {list(mismatched.items())[:3]}")
    if written != expected_selected:
        temporary.unlink(missing_ok=True)
        raise ValueError(f"Selected {written} profile sites instead of {expected_selected}")
    temporary.replace(destination)
    signature["selection_sha256"] = hashlib.sha256(destination.read_bytes()).hexdigest()
    marker.write_text(json.dumps(signature, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return destination, total_counts


def _profile_sums_by_motif(
    bw: pyBigWig.pyBigWig,
    sites: pd.DataFrame,
    flank: int,
    tile_bp: int = 2_000_000,
) -> tuple[dict[str, np.ndarray], dict[str, int]]:
    """Aggregate strand-oriented windows with tiled bigWig reads."""

    width = flank * 2
    sums = {str(motif): np.zeros(width, dtype=np.float64) for motif in sites["motif"].unique()}
    used = {motif: 0 for motif in sums}
    chrom_sizes = bw.chroms()
    working = sites.copy()
    working["center"] = (
        pd.to_numeric(working["TFBS_start"], errors="raise").to_numpy(dtype=np.int64)
        + pd.to_numeric(working["TFBS_end"], errors="raise").to_numpy(dtype=np.int64)
    ) // 2
    for chrom, chrom_sites in working.groupby("TFBS_chr", sort=False, observed=True):
        chrom = str(chrom)
        chrom_size = chrom_sizes.get(chrom)
        if chrom_size is None:
            continue
        chrom_sites = chrom_sites.loc[
            chrom_sites["center"].ge(flank) & chrom_sites["center"].le(chrom_size - flank)
        ].copy()
        if chrom_sites.empty:
            continue
        chrom_sites["tile"] = (chrom_sites["center"] // tile_bp).astype(int)
        for tile, tile_sites in chrom_sites.groupby("tile", sort=True, observed=True):
            tile_start = int(tile) * tile_bp
            tile_end = min(chrom_size, tile_start + tile_bp)
            read_start = max(0, tile_start - flank)
            read_end = min(chrom_size, tile_end + flank)
            signal = np.asarray(bw.values(chrom, read_start, read_end, numpy=True), dtype=float)
            signal = np.nan_to_num(signal, nan=0.0, posinf=0.0, neginf=0.0)
            for motif, motif_sites in tile_sites.groupby("motif", sort=False, observed=True):
                motif = str(motif)
                centers = motif_sites["center"].to_numpy(dtype=np.int64)
                strands = motif_sites["TFBS_strand"].astype(str).to_numpy()
                for start in range(0, len(centers), 4096):
                    batch_centers = centers[start:start + 4096]
                    offsets = batch_centers[:, None] - read_start + np.arange(-flank, flank)
                    values = signal[offsets]
                    reverse = strands[start:start + 4096] == "-"
                    if reverse.any():
                        values[reverse] = values[reverse, ::-1]
                    sums[motif] += values.sum(axis=0, dtype=np.float64)
                    used[motif] += len(values)
    return sums, used


def _compute_sample_profiles(task: tuple[str, str, str, str, int, int, dict[str, int]]) -> tuple[str, str]:
    sample, bigwig_path, sites_path, output_path, flank, max_sites_per_motif, total_counts = task
    sites = pd.read_csv(
        sites_path,
        sep="\t",
        usecols=["motif", "TFBS_chr", "TFBS_start", "TFBS_end", "TFBS_strand"],
    )
    selected_counts = sites.groupby("motif", sort=False, observed=True).size().astype(int).to_dict()
    expected_selected_counts = {
        motif: min(count, max_sites_per_motif)
        for motif, count in total_counts.items()
    }
    if selected_counts != expected_selected_counts:
        raise ValueError(f"{sample} profile-site selection counts are incomplete")
    records = {}
    with pyBigWig.open(bigwig_path) as bw:
        sums, used = _profile_sums_by_motif(bw, sites, flank)
    for motif in sorted(total_counts):
        count = used.get(motif, 0)
        profile = sums[motif] / count if count else np.zeros(flank * 2, dtype=float)
        records[motif] = {
            "n_sites": total_counts[motif],
            "n_profile_sites": count,
            "profile": [round(float(value), 6) for value in profile],
        }
    if len(records) != EXPECTED_MOTIFS:
        raise ValueError(f"{sample} produced {len(records)} profiles instead of {EXPECTED_MOTIFS}")
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(destination, "wt", encoding="utf-8", compresslevel=6) as handle:
        json.dump(
            {
                "sample": sample,
                "flank": flank,
                "max_sites_per_motif": max_sites_per_motif,
                "motifs": records,
            },
            handle,
            separators=(",", ":"),
            allow_nan=False,
        )
    return sample, str(destination)


def compute_profiles(
    project: Path,
    manifest: pd.DataFrame,
    flank: int,
    workers: int,
    max_sites_per_motif: int,
) -> None:
    local_dir = project / "reports" / "browser" / "sample_profiles"
    selected_sites, total_counts = prepare_shared_profile_sites(project, manifest, max_sites_per_motif)
    tasks = []
    for row in manifest.itertuples(index=False):
        bigwig = normalized_bigwig_path(project, row.sample)
        if not bigwig.is_file():
            raise ValueError(f"Missing normalized track for {row.sample}")
        output = local_dir / f"{row.sample}.json.gz"
        if output.is_file() and output.stat().st_size > 0:
            with gzip.open(output, "rt", encoding="utf-8") as handle:
                cached = json.load(handle)
            if (
                cached.get("sample") == row.sample
                and cached.get("flank") == flank
                and cached.get("max_sites_per_motif") == max_sites_per_motif
                and len(cached.get("motifs") or {}) == EXPECTED_MOTIFS
            ):
                continue
        tasks.append(
            (
                row.sample,
                str(bigwig),
                str(selected_sites),
                str(output),
                flank,
                max_sites_per_motif,
                total_counts,
            )
        )
    (project / "reports" / "browser").mkdir(parents=True, exist_ok=True)
    selection_sha256 = hashlib.sha256(selected_sites.read_bytes()).hexdigest()
    pd.DataFrame([
        {
            "sample": sample,
            "shared_coordinate_selection": str(selected_sites),
            "selection_sha256": selection_sha256,
        }
        for sample in manifest["sample"]
    ]).to_csv(project / "reports" / "browser" / "motif_coordinate_audit.tsv", sep="\t", index=False)
    if tasks:
        with concurrent.futures.ProcessPoolExecutor(max_workers=max(1, min(workers, len(tasks)))) as pool:
            for sample, output in pool.map(_compute_sample_profiles, tasks):
                print(f"profiles complete: {sample} -> {output}", flush=True)


def _json_number(value: object) -> float | int | None:
    if value is None or pd.isna(value):
        return None
    number = float(value)
    if not math.isfinite(number):
        return None
    return int(number) if number.is_integer() else number


def _motif_record(row: object) -> dict:
    moderated_df = float(row.moderated_df)
    return {
        "prefix": row.prefix,
        "name": row.name,
        "motif_id": row.motif_id,
        "cluster": row.cluster,
        "n_sites": int(row.n_sites),
        "mean1": _json_number(row.mean1),
        "sd1": _json_number(row.sd1),
        "mean2": _json_number(row.mean2),
        "sd2": _json_number(row.sd2),
        "effect": _json_number(row.effect),
        "ci_lower": _json_number(row.ci_lower),
        "ci_upper": _json_number(row.ci_upper),
        "moderated_t": _json_number(row.moderated_t),
        "moderated_df": None if np.isposinf(moderated_df) else _json_number(moderated_df),
        "normal_limit": bool(np.isposinf(moderated_df)),
        "pvalue": _json_number(row.pvalue),
        "qvalue": _json_number(row.qvalue),
        "significant": bool(row.significant),
    }


def _write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, separators=(",", ":"), ensure_ascii=False, allow_nan=False), encoding="utf-8")


def build_site(
    project: Path,
    site: Path,
    manifest: pd.DataFrame,
    comparisons: pd.DataFrame,
    spec: dict,
    flank: int,
    max_sites_per_motif: int,
) -> None:
    combined, summary = standardized_comparisons(project, comparisons)
    report_dir = project / "reports" / "browser"
    report_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(report_dir / "comparison_summary.tsv", sep="\t", index=False)
    data_dir = site / "data"
    comparison_dir = data_dir / "comparisons"
    profile_dir = data_dir / "profiles"
    for directory in (data_dir, comparison_dir, profile_dir):
        directory.mkdir(parents=True, exist_ok=True)
    for comparison in comparisons.itertuples(index=False):
        subset = combined.loc[combined["comparison"].eq(comparison.comparison)]
        payload = {
            "schema": "fp-tools.encode-cancer-comparison.v1",
            "id": comparison.comparison,
            "condition1": comparison.cond1,
            "condition2": comparison.cond2,
            "motifs": [_motif_record(row) for row in subset.itertuples(index=False)],
        }
        _write_json(comparison_dir / f"{comparison.comparison}.json", payload)
    combined.to_csv(data_dir / "all_pairwise_results.tsv.gz", sep="\t", index=False, compression="gzip", float_format="%.10g")
    summary.to_csv(data_dir / "comparison_summary.tsv", sep="\t", index=False)
    manifest.to_csv(data_dir / "samples.tsv", sep="\t", index=False)

    profiles_by_sample = {}
    for sample in manifest["sample"]:
        path = project / "reports" / "browser" / "sample_profiles" / f"{sample}.json.gz"
        if not path.is_file():
            raise ValueError(f"Missing precomputed profile cache: {path}")
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            profiles_by_sample[sample] = json.load(handle)["motifs"]
    motif_prefixes = sorted(set(combined["prefix"]))
    shards: dict[str, dict] = {f"{idx:02x}": {} for idx in range(32)}
    motif_index = {}
    for prefix in motif_prefixes:
        shard = f"{int(hashlib.sha256(prefix.encode('utf-8')).hexdigest()[:2], 16) % 32:02x}"
        site_counts = {profiles_by_sample[sample][prefix]["n_sites"] for sample in profiles_by_sample}
        if len(site_counts) != 1:
            raise ValueError(f"Motif-site counts differ across samples for {prefix}")
        profile_site_counts = {
            profiles_by_sample[sample][prefix]["n_profile_sites"]
            for sample in profiles_by_sample
        }
        if len(profile_site_counts) != 1:
            raise ValueError(f"Aggregate-profile site counts differ across samples for {prefix}")
        shards[shard][prefix] = {
            "n_sites": site_counts.pop(),
            "n_profile_sites": profile_site_counts.pop(),
            "samples": {
                sample: profiles_by_sample[sample][prefix]["profile"]
                for sample in manifest["sample"]
            },
        }
        motif_index[prefix] = shard
    for shard, motifs in shards.items():
        _write_json(profile_dir / f"{shard}.json", {"schema": "fp-tools.encode-cancer-profiles.v1", "motifs": motifs})
    _write_json(data_dir / "motif_index.json", motif_index)

    peak = json.loads((project / "peaks" / "peak_universe_qc.json").read_text(encoding="utf-8"))
    if {key: peak[key] for key in ("regions", "covered_bp", "md5")} != spec["peak_universe"]:
        raise ValueError("Project peak universe no longer matches the pinned specification")
    metadata = {
        "schema": "fp-tools.encode-cancer-browser.v1",
        "release_date": RELEASE_DATE,
        "genome": "GRCh38",
        "motif_database": "JASPAR 2026 CORE vertebrates non-redundant",
        "differential_normalization": "Sample-quantile normalization against shared background-score distributions",
        "motif_count": EXPECTED_MOTIFS,
        "comparison_count": EXPECTED_COMPARISONS,
        "sample_count": len(manifest),
        "profile_axis": list(range(-flank, flank)),
        "max_profile_sites_per_motif": max_sites_per_motif,
        "profile_site_selection": "Deterministic evenly spaced sites in genomic scan order; negative-strand windows are reverse-oriented.",
        "peak_universe": peak,
        "conditions": [
            {
                "name": condition,
                "samples": list(samples),
                "experiment": manifest.loc[manifest["condition"].eq(condition), "experiment"].iloc[0],
            }
            for condition, samples in spec["conditions"].items()
        ],
        "comparisons": [
            {
                "id": row.comparison,
                "condition1": row.cond1,
                "condition2": row.cond2,
                "file": f"data/comparisons/{row.comparison}.json",
            }
            for row in comparisons.itertuples(index=False)
        ],
        "downloads": {
            "all_results": "data/all_pairwise_results.tsv.gz",
            "comparison_summary": "data/comparison_summary.tsv",
            "samples": "data/samples.tsv",
        },
    }
    _write_json(data_dir / "metadata.json", metadata)
    size = sum(path.stat().st_size for path in site.rglob("*") if path.is_file())
    if size > 50 * 1024 * 1024:
        raise ValueError(f"Static browser exceeds the 50 MiB budget: {size / 1024 / 1024:.1f} MiB")
    print(json.dumps({"site": str(site), "bytes": size, "motifs": len(motif_prefixes), "comparisons": len(comparisons)}, indent=2))


def clean_intermediate_html(project: Path, comparisons: pd.DataFrame) -> None:
    directory = _result_path(project).parent
    for path in directory.glob("diff_footprints_*.html"):
        path.unlink()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("comparisons", "profiles", "site", "verify", "clean-html", "all"))
    parser.add_argument("--project", type=Path, default=DEFAULT_PROJECT)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--spec", type=Path, default=DEFAULT_SPEC)
    parser.add_argument("--comparisons", type=Path, default=DEFAULT_COMPARISONS)
    parser.add_argument("--genome", type=Path, default=DEFAULT_GENOME)
    parser.add_argument("--site", type=Path, default=DEFAULT_SITE)
    parser.add_argument("--cores", type=int, default=min(28, os.cpu_count() or 1))
    parser.add_argument("--profile-workers", type=int, default=4)
    parser.add_argument("--flank", type=int, default=100)
    parser.add_argument("--max-profile-sites-per-motif", type=int, default=1500)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    project = args.project.expanduser().resolve()
    manifest, comparisons, spec = read_design(args.manifest, args.spec, args.comparisons)
    if args.mode in {"comparisons", "all"}:
        run_comparisons(project, manifest, comparisons, args.genome.expanduser().resolve(), args.cores)
    if args.mode in {"profiles", "all"}:
        compute_profiles(
            project,
            manifest,
            args.flank,
            args.profile_workers,
            args.max_profile_sites_per_motif,
        )
    if args.mode in {"site", "all"}:
        build_site(
            project,
            args.site.expanduser().resolve(),
            manifest,
            comparisons,
            spec,
            args.flank,
            args.max_profile_sites_per_motif,
        )
    if args.mode in {"verify", "all"}:
        combined, summary = standardized_comparisons(project, comparisons)
        if len(combined) != 21_399 or len(summary) != 21:
            raise ValueError("Final result dimensions are incomplete")
        print(f"verified {len(combined)} motif-comparison rows across {len(summary)} comparisons")
    if args.mode == "clean-html":
        clean_intermediate_html(project, comparisons)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
