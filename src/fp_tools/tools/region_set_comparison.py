"""Differential footprint analysis between genomic region sets.

The region sets are the experimental groups. Biological samples are repeated
measurements of the same groups, and genomic regions (not motif hits) are the
unit of aggregation and resampling.
"""

from __future__ import annotations

import itertools
import os
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import pyBigWig
import pysam

from fp_tools.tools.diff_footprint_helpers import plot_interactive_diff_footprints
from fp_tools.utils.empirical_bayes import (
    benjamini_hochberg,
    fit_moderated_paired_contrast,
)
from fp_tools.utils.logger import FpToolsLogger
from fp_tools.utils.motifs import MotifList
from fp_tools.utils.regions import OneRegion, RegionList
from fp_tools.utils.utilities import expand_dirs, make_directory


@dataclass(frozen=True)
class RegionRecord:
    set_index: int
    set_label: str
    region_id: str
    chrom: str
    start: int
    end: int
    stratum: str
    fields: tuple[str, ...]


def _read_region_sets(paths, labels, strata_column=None):
    records = []
    for set_index, (path, label) in enumerate(zip(paths, labels)):
        seen = set()
        with open(path, encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                if not line.strip() or line.startswith(('#', 'track', 'browser')):
                    continue
                fields = tuple(line.rstrip("\n").split("\t"))
                if len(fields) < 3:
                    raise ValueError(f"{path}:{line_number} has fewer than three BED columns")
                try:
                    start, end = int(fields[1]), int(fields[2])
                except ValueError as exc:
                    raise ValueError(f"{path}:{line_number} has non-integer BED coordinates") from exc
                if start < 0 or end <= start:
                    raise ValueError(f"{path}:{line_number} has an invalid BED interval")
                coordinate = (fields[0], start, end)
                if coordinate in seen:
                    raise ValueError(f"Duplicate region in {path}: {fields[0]}:{start}-{end}")
                seen.add(coordinate)
                if strata_column is not None:
                    if strata_column < 1 or strata_column > len(fields):
                        raise ValueError(
                            f"{path}:{line_number} does not contain requested BED column {strata_column}"
                        )
                    stratum = fields[strata_column - 1].strip()
                    if not stratum:
                        raise ValueError(f"{path}:{line_number} has an empty matching stratum")
                else:
                    stratum = "all"
                records.append(RegionRecord(
                    set_index=set_index,
                    set_label=label,
                    region_id=f"{label}:{fields[0]}:{start}-{end}",
                    chrom=fields[0],
                    start=start,
                    end=end,
                    stratum=stratum,
                    fields=fields,
                ))
    if any(not any(record.set_index == idx for record in records) for idx in range(len(paths))):
        raise ValueError("Every --regions BED file must contain at least one interval")
    _reject_overlapping_regions(records)
    return records


def _reject_overlapping_regions(records):
    by_chrom = defaultdict(list)
    for record in records:
        by_chrom[record.chrom].append(record)
    for chrom, chrom_records in by_chrom.items():
        ordered = sorted(chrom_records, key=lambda item: (item.start, item.end, item.set_index))
        previous = ordered[0]
        for current in ordered[1:]:
            if current.start < previous.end:
                raise ValueError(
                    "Region sets must be mutually non-overlapping; found "
                    f"{previous.region_id} and {current.region_id}"
                )
            if current.end > previous.end:
                previous = current


def _validate_strata(records, labels):
    table = pd.DataFrame({
        "region_set": [record.set_label for record in records],
        "stratum": [record.stratum for record in records],
    })
    counts = table.groupby(["stratum", "region_set"], sort=True).size().unstack(fill_value=0)
    missing = [stratum for stratum, row in counts.iterrows() if any(row.get(label, 0) == 0 for label in labels)]
    if missing:
        preview = ", ".join(map(str, missing[:8]))
        raise ValueError(f"Every matching stratum must occur in every region set; missing set(s) in: {preview}")
    return counts


def _prepare_motifs(paths, sequences, naming, pvalue):
    gc_numerator = sum(sequence.upper().count("G") + sequence.upper().count("C") for sequence in sequences)
    gc_denominator = sum(len(sequence) for sequence in sequences)
    gc = gc_numerator / gc_denominator if gc_denominator else 0.5
    background = np.array([(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2], dtype=float)
    motifs = MotifList()
    for path in expand_dirs(paths):
        motifs += MotifList().from_file(path)
    if not motifs:
        raise ValueError("No motifs were found in the selected motif database")
    for motif in motifs:
        motif.bg = background
        motif.set_prefix(naming)
        motif.get_pssm()
        motif.get_threshold(pvalue)
    duplicated = {name for name, count in Counter(motif.prefix.upper() for motif in motifs).items() if count > 1}
    counters = {name: 1 for name in duplicated}
    for motif in motifs:
        key = motif.prefix.upper()
        if key in duplicated:
            motif.prefix = f"{motif.prefix}_{counters[key]}"
            counters[key] += 1
    motifs.setup_moods_scanner()
    return motifs


def _stratified_effect(frame, label_1, label_2, value_column):
    effects, weights = [], []
    for _stratum, subset in frame.groupby("stratum", sort=False):
        values_1 = pd.to_numeric(
            subset.loc[subset["region_set"] == label_1, value_column], errors="coerce"
        ).dropna().to_numpy(dtype=float)
        values_2 = pd.to_numeric(
            subset.loc[subset["region_set"] == label_2, value_column], errors="coerce"
        ).dropna().to_numpy(dtype=float)
        if not len(values_1) or not len(values_2):
            continue
        effects.append(float(np.mean(values_1) - np.mean(values_2)))
        weights.append(float(len(values_1) * len(values_2) / (len(values_1) + len(values_2))))
    if not weights:
        return np.nan
    return float(np.average(effects, weights=weights))


def _permutation_pvalue(frame, label_1, label_2, value_column, observed, n_permutations, rng):
    if not np.isfinite(observed) or n_permutations <= 0:
        return np.nan
    strata = []
    for _stratum, subset in frame.groupby("stratum", sort=False):
        values_1 = pd.to_numeric(
            subset.loc[subset["region_set"] == label_1, value_column], errors="coerce"
        ).dropna().to_numpy(dtype=float)
        values_2 = pd.to_numeric(
            subset.loc[subset["region_set"] == label_2, value_column], errors="coerce"
        ).dropna().to_numpy(dtype=float)
        if len(values_1) and len(values_2):
            strata.append((np.concatenate([values_1, values_2]), len(values_1), len(values_2)))
    if not strata:
        return np.nan
    weights = np.array([n1 * n2 / (n1 + n2) for _values, n1, n2 in strata], dtype=float)
    extreme = 0
    for _ in range(int(n_permutations)):
        effects = []
        for values, n1, _n2 in strata:
            shuffled = rng.permutation(values)
            effects.append(float(np.mean(shuffled[:n1]) - np.mean(shuffled[n1:])))
        permuted = float(np.average(effects, weights=weights))
        extreme += int(abs(permuted) >= abs(observed))
    return float((extreme + 1) / (int(n_permutations) + 1))


def _bootstrap_ci(frame, label_1, label_2, value_column, n_bootstrap, rng):
    if n_bootstrap <= 0:
        return (np.nan, np.nan)
    prepared = []
    for stratum, subset in frame.groupby("stratum", sort=False):
        values = {}
        for label in (label_1, label_2):
            values[label] = pd.to_numeric(
                subset.loc[subset["region_set"] == label, value_column], errors="coerce"
            ).dropna().to_numpy(dtype=float)
        if len(values[label_1]) and len(values[label_2]):
            prepared.append((stratum, values))
    estimates = []
    for _ in range(int(n_bootstrap)):
        rows = []
        for stratum, values in prepared:
            for label in (label_1, label_2):
                sampled = rng.choice(values[label], size=len(values[label]), replace=True)
                rows.extend({"stratum": stratum, "region_set": label, value_column: value} for value in sampled)
        estimates.append(_stratified_effect(pd.DataFrame(rows), label_1, label_2, value_column))
    finite = np.asarray(estimates, dtype=float)
    finite = finite[np.isfinite(finite)]
    return tuple(np.quantile(finite, [0.025, 0.975])) if len(finite) else (np.nan, np.nan)


def _scan_region_scores(records, motifs, genome, signals, sample_names):
    fasta = pysam.FastaFile(genome)
    bigwigs = [pyBigWig.open(path, "rb") for path in signals]
    motif_rows = defaultdict(list)
    try:
        fasta_bounds = dict(zip(fasta.references, fasta.lengths))
        signal_bounds = [handle.chroms() for handle in bigwigs]
        for record in records:
            if record.chrom not in fasta_bounds or record.end > fasta_bounds[record.chrom]:
                raise ValueError(f"Region is outside the genome FASTA: {record.region_id}")
            if any(record.chrom not in bounds or record.end > bounds[record.chrom] for bounds in signal_bounds):
                raise ValueError(f"Region is outside a signal bigWig: {record.region_id}")
            region = OneRegion([record.chrom, record.start, record.end])
            sequence = fasta.fetch(record.chrom, record.start, record.end)
            sites = motifs.scan_sequence(sequence, region)
            by_motif = defaultdict(RegionList)
            for site in sites:
                by_motif[site.name].append(site)
            for prefix, motif_sites in by_motif.items():
                motif_sites = motif_sites.resolve_overlaps()
                sample_site_scores = [[] for _ in signals]
                centers = []
                for site in motif_sites:
                    center = int(site.start + (site.end - site.start) / 2)
                    centers.append((site.chrom, center, site.strand))
                    for sample_index, handle in enumerate(bigwigs):
                        value = handle.values(site.chrom, center, center + 1, numpy=True)
                        value = float(np.nan_to_num(value[0], nan=0.0))
                        sample_site_scores[sample_index].append(value)
                row = {
                    "region_set": record.set_label,
                    "set_index": record.set_index,
                    "region_id": record.region_id,
                    "chrom": record.chrom,
                    "start": record.start,
                    "end": record.end,
                    "stratum": record.stratum,
                    "n_sites": len(centers),
                    "centers": centers,
                }
                for sample_name, values in zip(sample_names, sample_site_scores):
                    row[sample_name] = float(np.mean(values))
                motif_rows[prefix].append(row)
    finally:
        fasta.close()
        for handle in bigwigs:
            handle.close()
    return motif_rows


def _comparison_rows(motifs, motif_rows, labels, sample_names, args):
    results, replicate_rows = [], []
    rng = np.random.default_rng(int(args.random_seed))
    total_by_set = Counter(record.set_label for record in args.region_records)
    pairs = list(itertools.combinations(labels, 2))
    for label_1, label_2 in pairs:
        comparison = f"{label_1}_vs_{label_2}"
        pair_results = []
        effects_by_motif = {}
        row_lookup = {}
        for motif in motifs:
            frame = pd.DataFrame(motif_rows.get(motif.prefix, []))
            if frame.empty:
                frame = pd.DataFrame(columns=["region_set", "stratum", "n_sites", *sample_names])
            frame = frame[frame["region_set"].isin([label_1, label_2])]
            counts = frame.groupby("region_set").size().to_dict()
            site_counts = frame.groupby("region_set")["n_sites"].sum().to_dict() if not frame.empty else {}
            if frame.empty:
                analysis_frame = frame
            else:
                stratum_counts = frame.groupby(["stratum", "region_set"]).size().unstack(fill_value=0)
                shared_strata = stratum_counts.index[
                    (stratum_counts.get(label_1, 0) > 0) & (stratum_counts.get(label_2, 0) > 0)
                ]
                analysis_frame = frame[frame["stratum"].isin(shared_strata)]
            analyzed_counts = analysis_frame.groupby("region_set").size().to_dict()
            sufficient = min(
                analyzed_counts.get(label_1, 0), analyzed_counts.get(label_2, 0)
            ) >= int(args.min_regions_per_set)
            replicate_effects = {}
            if sufficient:
                for sample in sample_names:
                    replicate_effects[sample] = _stratified_effect(
                        analysis_frame, label_1, label_2, sample
                    )
                sufficient = all(np.isfinite(value) for value in replicate_effects.values())
            if sufficient:
                for sample in sample_names:
                    replicate_rows.append({
                        "comparison": comparison,
                        "motif": motif.prefix,
                        "sample": sample,
                        "effect": replicate_effects[sample],
                    })
            effects_by_motif[motif.prefix] = replicate_effects
            row = {
                "comparison": comparison,
                "region_set_1": label_1,
                "region_set_2": label_2,
                "output_prefix": motif.prefix,
                "name": motif.name,
                "motif_id": motif.id,
                "effect": np.nan,
                "ci_lower": np.nan,
                "ci_upper": np.nan,
                "pvalue": np.nan,
                "qvalue_bh": np.nan,
                "significant_fdr05": False,
                "status": "tested" if sufficient else "insufficient motif-containing regions",
                "statistical_method": (
                    "paired empirical-Bayes moderated t" if len(sample_names) >= 2
                    else "within-stratum label permutation"
                ),
                "n_replicates": len(sample_names),
                "n_regions_set_1": total_by_set[label_1],
                "n_regions_set_2": total_by_set[label_2],
                "n_motif_regions_set_1": counts.get(label_1, 0),
                "n_motif_regions_set_2": counts.get(label_2, 0),
                "n_analyzed_regions_set_1": analyzed_counts.get(label_1, 0),
                "n_analyzed_regions_set_2": analyzed_counts.get(label_2, 0),
                "motif_prevalence_set_1": counts.get(label_1, 0) / total_by_set[label_1],
                "motif_prevalence_set_2": counts.get(label_2, 0) / total_by_set[label_2],
                "n_motif_sites_set_1": int(site_counts.get(label_1, 0)),
                "n_motif_sites_set_2": int(site_counts.get(label_2, 0)),
            }
            if sufficient and len(sample_names) == 1:
                sample = sample_names[0]
                observed = replicate_effects[sample]
                row["effect"] = observed
                row["pvalue"] = _permutation_pvalue(
                    analysis_frame, label_1, label_2, sample, observed,
                    int(args.region_permutations), rng,
                )
                row["ci_lower"], row["ci_upper"] = _bootstrap_ci(
                    analysis_frame, label_1, label_2, sample,
                    int(args.region_bootstrap), rng,
                )
            pair_results.append(row)
            row_lookup[motif.prefix] = row

        if len(sample_names) >= 2:
            matrix = pd.DataFrame.from_dict(effects_by_motif, orient="index", columns=sample_names)
            complete = matrix.dropna()
            if len(complete) >= 2:
                model = fit_moderated_paired_contrast(matrix)
                for prefix, model_row in model.iterrows():
                    if prefix not in row_lookup or not np.isfinite(model_row["pvalue"]):
                        continue
                    row = row_lookup[prefix]
                    for column in ("effect", "ci_lower", "ci_upper", "pvalue", "qvalue_bh"):
                        row[column] = float(model_row[column])
                    row["significant_fdr05"] = bool(model_row["significant_fdr05"])
        else:
            tested = [row for row in pair_results if np.isfinite(row["pvalue"])]
            qvalues = benjamini_hochberg([row["pvalue"] for row in tested])
            for row, qvalue in zip(tested, qvalues):
                row["qvalue_bh"] = float(qvalue)
                row["significant_fdr05"] = bool(qvalue <= 0.05)
        results.extend(pair_results)
    return pd.DataFrame(results), pd.DataFrame(replicate_rows)


def _aggregate_payload(motifs, motif_rows, result_rows, labels, signals, sample_names, args):
    mode = getattr(args, "plot_aggregate", "sig")
    if mode == "off":
        return {"x": [], "motifs": []}
    ranked = result_rows.copy()
    ranked["p_sort"] = pd.to_numeric(ranked["pvalue"], errors="coerce").fillna(1.0)
    ranked["effect_sort"] = pd.to_numeric(ranked["effect"], errors="coerce").abs().fillna(0.0)
    ranked = ranked.sort_values(["p_sort", "effect_sort"], ascending=[True, False])
    top_n = max(1, int(getattr(args, "plot_aggregate_top_n", 20)))
    if mode == "all":
        selected = ranked
    elif mode == "top":
        selected = ranked.head(top_n)
    else:
        selected = ranked[pd.to_numeric(ranked["qvalue_bh"], errors="coerce") <= 0.05].head(top_n)
        if selected.empty:
            selected = ranked.head(top_n)
    flank = max(1, int(getattr(args, "aggregate_flank", 100)))
    x = list(range(-flank, flank))
    handles = [pyBigWig.open(path, "rb") for path in signals]
    motif_lookup = {motif.prefix: motif for motif in motifs}
    payloads = []
    try:
        for prefix in selected["output_prefix"]:
            motif = motif_lookup[prefix]
            rows = motif_rows.get(prefix, [])
            conditions = []
            for label in labels:
                label_rows = [row for row in rows if row["region_set"] == label]
                sample_payloads = []
                for sample, handle in zip(sample_names, handles):
                    enhancer_profiles = []
                    for row in label_rows:
                        site_profiles = []
                        for chrom, center, strand in row["centers"]:
                            start, end = center - flank, center + flank
                            if start < 0 or end > handle.chroms(chrom):
                                continue
                            values = np.nan_to_num(handle.values(chrom, start, end, numpy=True), nan=0.0)
                            if strand == "-":
                                values = values[::-1]
                            site_profiles.append(values)
                        if site_profiles:
                            enhancer_profiles.append(np.mean(np.asarray(site_profiles), axis=0))
                    profile = (
                        np.mean(np.asarray(enhancer_profiles), axis=0)
                        if enhancer_profiles else np.zeros(len(x), dtype=float)
                    )
                    sample_payloads.append({
                        "name": f"{label}_{sample}",
                        "profile": [round(float(value), 6) for value in profile],
                        "fp_score": round(float(np.mean(profile)), 6),
                    })
                mean_profile = np.mean(
                    np.asarray([sample["profile"] for sample in sample_payloads], dtype=float), axis=0
                )
                conditions.append({
                    "name": label,
                    "profile": [round(float(value), 6) for value in mean_profile],
                    "samples": sample_payloads,
                    "n_sites": len(label_rows),
                    "fp_score": round(float(np.mean(mean_profile)), 6),
                })
            result = selected.loc[selected["output_prefix"] == prefix].iloc[0]
            payloads.append({
                "prefix": prefix,
                "name": motif.name,
                "motif_id": motif.id,
                "change": float(result["effect"]) if np.isfinite(result["effect"]) else 0.0,
                "pvalue": float(result["pvalue"]) if np.isfinite(result["pvalue"]) else 1.0,
                "n_sites": sum(condition["n_sites"] for condition in conditions),
                "site_set": "all motif-containing regions",
                "conditions": conditions,
            })
    finally:
        for handle in handles:
            handle.close()
    return {
        "x": x,
        "motifs": payloads,
        "comparison": " / ".join(labels),
        "normalization": "none",
        "site_set": "equal enhancer weight",
        "max_sites_per_motif": None,
        "x_label": "Distance from motif center (bp)",
        "y_label": "Corrected cut-site signal" if getattr(args, "aggregate_signals", None) else "Footprint score signal",
    }


def _write_reports(results, motifs, motif_rows, labels, args):
    motif_lookup = {motif.prefix: motif for motif in motifs}
    aggregate_signals = list(getattr(args, "aggregate_signals", None) or args.signals)
    for comparison, pair_rows in results.groupby("comparison", sort=False):
        label_1 = pair_rows.iloc[0]["region_set_1"]
        label_2 = pair_rows.iloc[0]["region_set_2"]
        report_motifs = []
        for _, row in pair_rows.iterrows():
            motif = motif_lookup[row["output_prefix"]]
            motif.change = float(row["effect"]) if np.isfinite(row["effect"]) else 0.0
            motif.pvalue = float(row["pvalue"]) if np.isfinite(row["pvalue"]) else 1.0
            motif.qvalue = float(row["qvalue_bh"]) if np.isfinite(row["qvalue_bh"]) else 1.0
            motif.group = (
                f"{label_1}_up" if motif.change > 0 and motif.qvalue <= 0.05
                else f"{label_2}_up" if motif.change < 0 and motif.qvalue <= 0.05
                else "n.s."
            )
            motif.region_stats = row.to_dict()
            report_motifs.append(motif)
        aggregate = _aggregate_payload(
            motifs, motif_rows, pair_rows, [label_1, label_2], aggregate_signals,
            args.sample_names, args,
        )
        html_path = os.path.join(args.outdir, f"{args.prefix}_{comparison}.html")
        plot_interactive_diff_footprints(
            report_motifs,
            [label_1, label_2],
            html_path,
            aggregate_data=aggregate,
            title="Region-set footprint report",
            report_label="Equal enhancer weighting; matching strata preserved",
            change_label="Stratum-adjusted footprint-score difference",
            results_table=pair_rows,
        )


def run_region_set_comparison(args):
    """Run ``diff-footprints --comparison-axis regions``."""

    logger = FpToolsLogger("diff-footprints", getattr(args, "verbosity", 1))
    logger.begin()
    if not args.genome:
        raise ValueError("--genome is required with --comparison-axis regions")
    if getattr(args, "normalization", "none") != "none" and not getattr(args, "norm_off", False):
        raise ValueError("Region-set comparisons currently require --normalization none")
    if int(args.region_permutations) < 1 or int(args.region_bootstrap) < 1:
        raise ValueError("--region-permutations and --region-bootstrap must be positive")
    if int(args.min_regions_per_set) < 2:
        raise ValueError("--min-regions-per-set must be at least 2")

    args.outdir = os.path.abspath(args.outdir)
    make_directory(args.outdir)
    paths = [os.path.abspath(path) for path in args.regions]
    labels = list(args.region_labels or [Path(path).stem for path in paths])
    if len(paths) < 2:
        raise ValueError("At least two --regions BED files are required")
    if len(labels) != len(paths):
        raise ValueError("--region-labels must contain one label per --regions BED file")
    if len(set(labels)) != len(labels):
        raise ValueError("--region-labels must be unique")
    sample_names = list(args.sample_names or [Path(path).stem for path in args.signals])
    if len(sample_names) != len(args.signals) or len(set(sample_names)) != len(sample_names):
        raise ValueError("--sample-names must provide one unique name per --signals bigWig")
    args.sample_names = sample_names
    if getattr(args, "aggregate_signals", None) and len(args.aggregate_signals) != len(args.signals):
        raise ValueError("--aggregate-signals must contain one bigWig per --signals bigWig")

    logger.info("Reading and validating region sets")
    records = _read_region_sets(paths, labels, getattr(args, "region_strata_column", None))
    args.region_records = records
    balance = _validate_strata(records, labels)
    balance_long = balance.reset_index().melt(id_vars="stratum", var_name="region_set", value_name="n_regions")
    balance_path = os.path.join(args.outdir, args.prefix + "_region_balance.tsv")
    balance_long.to_csv(balance_path, sep="\t", index=False)

    logger.info("Reading motifs and scanning region sequences")
    fasta = pysam.FastaFile(args.genome)
    try:
        sequences = [fasta.fetch(record.chrom, record.start, record.end) for record in records]
    finally:
        fasta.close()
    motifs = _prepare_motifs(args.motifs, sequences, args.naming, args.motif_pvalue)
    motif_rows = _scan_region_scores(records, motifs, args.genome, args.signals, sample_names)

    logger.info("Calculating region-set contrasts")
    results, replicate_effects = _comparison_rows(motifs, motif_rows, labels, sample_names, args)
    results_path = os.path.join(args.outdir, args.prefix + "_results.txt")
    effects_path = os.path.join(args.outdir, args.prefix + "_region_replicate_effects.tsv")
    results.to_csv(results_path, sep="\t", index=False, na_rep="NA")
    replicate_effects.to_csv(effects_path, sep="\t", index=False, na_rep="NA")
    if not getattr(args, "skip_excel", False):
        excel_path = os.path.join(args.outdir, args.prefix + "_results.xlsx")
        with pd.ExcelWriter(excel_path, engine="xlsxwriter") as writer:
            results.to_excel(writer, sheet_name="Results", index=False)
            replicate_effects.to_excel(writer, sheet_name="Replicate effects", index=False)
            balance_long.to_excel(writer, sheet_name="Region balance", index=False)

    logger.info("Writing interactive region-set report(s)")
    _write_reports(results, motifs, motif_rows, labels, args)
    logger.end()
    return results
