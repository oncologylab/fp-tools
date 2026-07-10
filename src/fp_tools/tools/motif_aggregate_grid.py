"""Static multi-page motif aggregate grids from multi-comparison reports."""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import math
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

from fp_tools.tools.diff_footprint_helpers import _aggregate_fp_score, _mean_profile, _read_bed_centers
from fp_tools.tools.plot_aggregate_batch import read_embedded_payload
from fp_tools.utils.project_layout import corrected_bigwig_path, is_project_layout, match_motifs_dir, project_root, read_sample_table, reports_dir


NUTRIENT_GROUP_ORDER = {
    "FBS": 0,
    "Glc": 1,
    "Met.Cys": 2,
    "Gln.Arg": 3,
    "Gln": 4,
    "Arg": 5,
    "BCAA": 6,
    "Trp": 7,
    "Lys": 8,
}


@dataclass(frozen=True)
class ComparisonView:
    label: str
    payload: dict
    condition: str
    index: int


@dataclass(frozen=True)
class MotifView:
    prefix: str
    name: str
    motif_id: str
    sort_score: float


@dataclass(frozen=True)
class SampleView:
    sample: str
    condition: str
    corrected_bigwig: Path
    match_dir: Path


@dataclass(frozen=True)
class RnaFcView:
    label: str
    values: tuple[tuple[str, float], ...]


def _as_float(value, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    if not math.isfinite(out):
        return default
    return out


def _comparison_condition(label: str, payload: dict | None = None) -> str:
    text = str(label or "").strip()
    if " vs " in text:
        return text.split(" vs ", 1)[0].strip()
    if payload:
        conditions = payload.get("conditions") or []
        if conditions:
            return str(conditions[0])
    return text


def _split_condition(condition: str) -> tuple[float, str]:
    match = re.match(r"^\s*([0-9]+(?:\.[0-9]+)?)_(.+?)\s*$", condition)
    if match:
        return float(match.group(1)), match.group(2)
    return float("-inf"), condition


def nutrient_sort_key(condition: str) -> tuple[int, float, str]:
    concentration, nutrient = _split_condition(condition)
    group_rank = NUTRIENT_GROUP_ORDER.get(nutrient, len(NUTRIENT_GROUP_ORDER))
    return (group_rank, -concentration, condition)


def ordered_comparisons(review_payload: dict) -> list[ComparisonView]:
    views = []
    for idx, item in enumerate(review_payload.get("comparisons") or []):
        payload = item.get("payload") or {}
        label = str(item.get("label") or f"Comparison {idx + 1}")
        condition = _comparison_condition(label, payload)
        views.append(ComparisonView(label=label, payload=payload, condition=condition, index=idx))
    return sorted(views, key=lambda item: (*nutrient_sort_key(item.condition), item.index))


def _motif_label(item: dict | MotifView) -> str:
    name = str(item.name if isinstance(item, MotifView) else item.get("name") or "")
    motif_id = str(item.motif_id if isinstance(item, MotifView) else item.get("motif_id") or "")
    return f"{name} ({motif_id})" if motif_id else name


def _point_map(payload: dict) -> dict[str, dict]:
    return {str(point.get("prefix")): point for point in payload.get("points") or [] if point.get("prefix")}


def _aggregate_map(payload: dict) -> dict[str, dict]:
    aggregate = payload.get("aggregate") or {}
    return {str(motif.get("prefix")): motif for motif in aggregate.get("motifs") or [] if motif.get("prefix")}


def _split_tf_symbols(value: str) -> list[str]:
    text = str(value or "").strip()
    if not text:
        return []
    parts = re.split(r"::|/|;|,", text)
    out = []
    seen = set()
    for part in parts:
        gene = re.sub(r"[^A-Za-z0-9-]", "", part).upper()
        if gene and gene not in seen:
            seen.add(gene)
            out.append(gene)
    return out


def _motif_id_base(value: str) -> str:
    match = re.search(r"(MA\d+)", str(value or ""))
    return match.group(1) if match else ""


def _read_motif_gene_map(path: str | Path | None) -> tuple[dict[str, list[str]], dict[str, list[str]]]:
    if not path:
        return {}, {}
    table = pd.read_csv(path, sep="\t")
    if "motif" not in table.columns or "gene_symbol" not in table.columns:
        raise ValueError("motif gene map must contain motif and gene_symbol columns")
    by_motif: dict[str, list[str]] = {}
    by_base: dict[str, list[str]] = {}
    for _, row in table.iterrows():
        motif = str(row.get("motif") or "").strip()
        genes = _split_tf_symbols(str(row.get("gene_symbol") or ""))
        if not motif or not genes:
            continue
        by_motif[motif] = genes
        base = _motif_id_base(motif)
        if base:
            by_base.setdefault(base, genes)
    return by_motif, by_base


def _genes_for_motif(motif: MotifView, motif_map: dict[str, list[str]], base_map: dict[str, list[str]]) -> list[str]:
    if motif.prefix in motif_map:
        return motif_map[motif.prefix]
    if motif.motif_id in motif_map:
        return motif_map[motif.motif_id]
    base = _motif_id_base(motif.motif_id) or _motif_id_base(motif.prefix)
    if base and base in base_map:
        return base_map[base]
    return _split_tf_symbols(motif.name)


def _format_rna_fc(items: list[tuple[str, float]]) -> str:
    if not items:
        return ""
    first = f"RNA {items[0][0]}={items[0][1]:+.2f}"
    if len(items) == 1:
        return first
    return first + f"\n{items[1][0]}={items[1][1]:+.2f}"


def collect_motifs(comparisons: Iterable[ComparisonView]) -> list[MotifView]:
    records: dict[str, dict] = {}
    max_changes: dict[str, float] = {}
    for comparison in comparisons:
        points = _point_map(comparison.payload)
        for prefix, motif in _aggregate_map(comparison.payload).items():
            point = points.get(prefix) or motif
            entry = records.setdefault(
                prefix,
                {
                    "prefix": prefix,
                    "name": str(point.get("name") or motif.get("name") or prefix),
                    "motif_id": str(point.get("motif_id") or motif.get("motif_id") or ""),
                },
            )
            if not entry.get("motif_id") and (point.get("motif_id") or motif.get("motif_id")):
                entry["motif_id"] = str(point.get("motif_id") or motif.get("motif_id"))
            max_changes[prefix] = max(max_changes.get(prefix, 0.0), abs(_as_float(point.get("change"))))
    motifs = [
        MotifView(prefix=prefix, name=record["name"], motif_id=record["motif_id"], sort_score=max_changes.get(prefix, 0.0))
        for prefix, record in records.items()
    ]
    return sorted(motifs, key=lambda item: (-item.sort_score, _motif_label(item).lower(), item.prefix))


def collect_motifs_from_payloads(payloads: Iterable[dict]) -> list[MotifView]:
    records: dict[str, dict] = {}
    max_changes: dict[str, float] = {}
    for payload in payloads:
        for comparison in ordered_comparisons(payload):
            points = _point_map(comparison.payload)
            aggregate = _aggregate_map(comparison.payload)
            for prefix, point in points.items():
                motif = aggregate.get(prefix, {})
                records.setdefault(
                    prefix,
                    {
                        "prefix": prefix,
                        "name": str(point.get("name") or motif.get("name") or prefix),
                        "motif_id": str(point.get("motif_id") or motif.get("motif_id") or ""),
                    },
                )
                max_changes[prefix] = max(max_changes.get(prefix, 0.0), abs(_as_float(point.get("change"))))
    motifs = [
        MotifView(prefix=prefix, name=record["name"], motif_id=record["motif_id"], sort_score=max_changes.get(prefix, 0.0))
        for prefix, record in records.items()
    ]
    return sorted(motifs, key=lambda item: (-item.sort_score, _motif_label(item).lower(), item.prefix))


def _ordered_motifs_for_payload(review_payload: dict, motif_order: list[MotifView] | None = None) -> list[MotifView]:
    if motif_order is None:
        return collect_motifs_from_payloads([review_payload])
    point_prefixes = {prefix for comparison in ordered_comparisons(review_payload) for prefix in _point_map(comparison.payload)}
    return [motif for motif in motif_order if motif.prefix in point_prefixes]


def _trim_profile(x_values: list, profile: list, flank: int) -> tuple[list[float], list[float]]:
    pairs = []
    for x, y in zip(x_values, profile):
        xf = _as_float(x, default=float("nan"))
        yf = _as_float(y, default=float("nan"))
        if math.isfinite(xf) and math.isfinite(yf) and -flank <= xf <= flank:
            pairs.append((xf, yf))
    if not pairs:
        for x, y in zip(x_values, profile):
            xf = _as_float(x, default=float("nan"))
            yf = _as_float(y, default=float("nan"))
            if math.isfinite(xf) and math.isfinite(yf):
                pairs.append((xf, yf))
    if not pairs:
        return [], []
    return [p[0] for p in pairs], [p[1] for p in pairs]


def _condition_color(payload: dict, condition: str, index: int) -> str:
    colors = payload.get("colors") or {}
    fallback = ["#dc2626", "#2563eb", "#16a34a", "#9333ea", "#f97316", "#0891b2"]
    return str(colors.get(f"{condition}_up") or colors.get(condition) or fallback[index % len(fallback)])


def _samples_for_motif(payload: dict, motif: dict, flank: int) -> list[dict]:
    x_values = motif.get("x") or (payload.get("aggregate") or {}).get("x") or []
    rows = []
    for cond_idx, condition in enumerate(motif.get("conditions") or []):
        condition_name = str(condition.get("name") or f"condition{cond_idx + 1}")
        color = _condition_color(payload, condition_name, cond_idx)
        for sample in condition.get("samples") or []:
            x, y = _trim_profile(x_values, sample.get("profile") or [], flank)
            if not x:
                continue
            rows.append(
                {
                    "condition": condition_name,
                    "sample": str(sample.get("name") or condition_name),
                    "x": x,
                    "y": y,
                    "color": color,
                    "score": _as_float(sample.get("fp_score")),
                }
            )
    return sorted(rows, key=lambda row: (row["condition"], row["sample"]))


def _read_project_samples(project: str | Path) -> list[SampleView]:
    project = Path(project)
    sample_table = _project_sample_table(project)
    if sample_table is None:
        return []
    rows = []
    for record in read_sample_table(sample_table):
        rows.append(SampleView(sample=record.sample, condition=record.condition, corrected_bigwig=corrected_bigwig_path(project, record.sample), match_dir=match_motifs_dir(project, record.sample)))
    return rows


def _project_sample_table(project: str | Path) -> Path | None:
    project = Path(project)
    candidates = [
        project / "metadata" / "samples.tsv",
        project / "metadata" / "sample_table.tsv",
        project / "samples.tsv",
        project / "sample_table.tsv",
    ]
    for path in candidates:
        if path.exists() and path.stat().st_size > 0:
            return path
    return None


def _read_sample_conditions(project: str | Path) -> dict[str, str]:
    table = _project_sample_table(project)
    if table is None:
        return {}
    return {record.sample: record.condition for record in read_sample_table(table)}


def _mean_by_condition(matrix: pd.DataFrame, sample_conditions: dict[str, str], samples: list[str]) -> dict[str, pd.Series]:
    available = [sample for sample in samples if sample in matrix.columns and sample in sample_conditions]
    out = {}
    for condition in sorted({sample_conditions[sample] for sample in available}):
        cols = [sample for sample in available if sample_conditions[sample] == condition]
        if cols:
            out[condition] = matrix[cols].mean(axis=1)
    return out


def compute_rna_fc_map(
    review_payload: dict,
    project: str | Path | None,
    motif_order: list[MotifView] | None,
    rna_log2norm: str | Path | None,
    rna_raw_counts: str | Path | None,
    motif_gene_map: str | Path | None,
    min_raw_mean: float = 1.0,
) -> dict[tuple[int, str], RnaFcView]:
    if not project or not rna_log2norm or not motif_gene_map:
        return {}
    sample_conditions = _read_sample_conditions(project)
    if not sample_conditions:
        return {}
    expr = pd.read_csv(rna_log2norm, sep="\t")
    if "gene_key" not in expr.columns:
        raise ValueError("RNA log2-normalized matrix must contain a gene_key column")
    expr["gene_key_upper"] = expr["gene_key"].astype(str).str.upper()
    expr = expr.drop_duplicates("gene_key_upper").set_index("gene_key_upper")
    sample_cols = [col for col in expr.columns if col in sample_conditions]
    expr = expr[sample_cols].apply(pd.to_numeric, errors="coerce")
    raw = None
    raw_means = None
    if rna_raw_counts:
        raw = pd.read_csv(rna_raw_counts, sep="\t")
        if "gene_key" not in raw.columns:
            raise ValueError("RNA raw count matrix must contain a gene_key column")
        raw["gene_key_upper"] = raw["gene_key"].astype(str).str.upper()
        raw = raw.drop_duplicates("gene_key_upper").set_index("gene_key_upper")
        raw_cols = [col for col in raw.columns if col in sample_conditions]
        raw = raw[raw_cols].apply(pd.to_numeric, errors="coerce")
        raw_means = _mean_by_condition(raw, sample_conditions, raw_cols)
    expr_means = _mean_by_condition(expr, sample_conditions, sample_cols)
    motif_map, base_map = _read_motif_gene_map(motif_gene_map)
    comparisons = ordered_comparisons(review_payload)
    motifs = _ordered_motifs_for_payload(review_payload, motif_order)
    motif_by_prefix = {motif.prefix: motif for motif in motifs}
    rna_map: dict[tuple[int, str], RnaFcView] = {}
    for comparison in comparisons:
        conditions = _wanted_conditions(comparison)
        if len(conditions) < 2 or conditions[0] not in expr_means or conditions[1] not in expr_means:
            continue
        cond1, cond2 = conditions[:2]
        points = _point_map(comparison.payload)
        for prefix in points:
            motif = motif_by_prefix.get(prefix)
            if not motif:
                point = points[prefix]
                motif = MotifView(prefix=prefix, name=str(point.get("name") or prefix), motif_id=str(point.get("motif_id") or ""), sort_score=0.0)
            values = []
            for gene in _genes_for_motif(motif, motif_map, base_map):
                gene_key = gene.upper()
                if gene_key not in expr.index:
                    continue
                if raw_means is not None and cond1 in raw_means and cond2 in raw_means:
                    max_raw = max(_as_float(raw_means[cond1].get(gene_key), float("nan")), _as_float(raw_means[cond2].get(gene_key), float("nan")))
                    if not math.isfinite(max_raw) or max_raw < min_raw_mean:
                        continue
                fc = _as_float(expr_means[cond1].get(gene_key), float("nan")) - _as_float(expr_means[cond2].get(gene_key), float("nan"))
                if math.isfinite(fc):
                    values.append((gene, fc))
            if values:
                rna_map[(comparison.index, prefix)] = RnaFcView(label=_format_rna_fc(values), values=tuple(values))
    return rna_map


def _motif_all_bed(match_dir: Path, prefix: str, cache: dict[tuple[str, str], Path | None] | None = None) -> Path | None:
    cache_key = (str(match_dir), prefix)
    if cache is not None and cache_key in cache:
        return cache[cache_key]
    candidates = [
        match_dir / prefix / "beds" / f"{prefix}_all.bed",
        match_dir / prefix / f"{prefix}_all.bed",
    ]
    for path in candidates:
        if path.exists() and path.stat().st_size > 0:
            if cache is not None:
                cache[cache_key] = path
            return path
    matches = sorted(match_dir.glob(f"*/beds/{prefix}_all.bed")) + sorted(match_dir.glob(f"*/{prefix}_all.bed"))
    for path in matches:
        if path.exists() and path.stat().st_size > 0:
            if cache is not None:
                cache[cache_key] = path
            return path
    if cache is not None:
        cache[cache_key] = None
    return None


def _profile_cache_from_payload(review_payload: dict) -> dict[tuple[str, str], dict]:
    cache: dict[tuple[str, str], dict] = {}
    for comparison in ordered_comparisons(review_payload):
        x_values = (comparison.payload.get("aggregate") or {}).get("x") or []
        for prefix, aggregate in _aggregate_map(comparison.payload).items():
            for condition in aggregate.get("conditions") or []:
                condition_name = str(condition.get("name") or "")
                if condition_name:
                    condition_record = dict(condition)
                    condition_record["_x"] = x_values
                    cache.setdefault((prefix, condition_name), condition_record)
    return cache


def _assemble_cached_aggregate(prefix: str, point: dict, comparison: ComparisonView, condition_cache: dict[tuple[str, str], dict], flank: int) -> dict | None:
    wanted_conditions = list((comparison.payload.get("conditions") or [])[:2])
    if len(wanted_conditions) < 2:
        wanted_conditions = [comparison.condition, "10_FBS_Ctrl"]
    conditions = []
    x_values = None
    for condition in wanted_conditions:
        cached = condition_cache.get((prefix, condition))
        if not cached:
            return None
        if x_values is None and cached.get("_x"):
            x_values = cached.get("_x")
        conditions.append(cached)
    n_sites = max([_as_float(condition.get("n_sites"), 0.0) for condition in conditions] or [0.0])
    return {
        "prefix": prefix,
        "name": str(point.get("name") or prefix),
        "motif_id": str(point.get("motif_id") or ""),
        "change": _as_float(point.get("change")),
        "pvalue": _as_float(point.get("pvalue"), 1.0),
        "fdr": point.get("fdr", ""),
        "n_sites": int(n_sites),
        "site_set": "all",
        "x": x_values or list(range(-flank, flank)),
        "conditions": conditions,
        "profile_source": "assembled",
    }


def _recompute_aggregate_motif(
    prefix: str,
    point: dict,
    comparison: ComparisonView,
    project_samples: list[SampleView],
    flank: int,
    bed_cache: dict[tuple[str, str], Path | None] | None = None,
    centers_cache: dict[str, list] | None = None,
    profile_cache: dict[tuple[str, str], list[float]] | None = None,
) -> dict | None:
    wanted_conditions = list((comparison.payload.get("conditions") or [])[:2])
    if len(wanted_conditions) < 2:
        wanted_conditions = [comparison.condition, "10_FBS_Ctrl"]
    sample_rows_by_condition = {cond: [sample for sample in project_samples if sample.condition == cond] for cond in wanted_conditions}
    conditions = []
    all_centers = []
    seen_centers = set()
    for condition in wanted_conditions:
        samples = []
        condition_profiles = []
        condition_centers = []
        for sample in sample_rows_by_condition.get(condition, []):
            bed = _motif_all_bed(sample.match_dir, prefix, cache=bed_cache)
            if bed is None or not sample.corrected_bigwig.exists():
                continue
            bed_key = str(bed)
            if centers_cache is not None and bed_key in centers_cache:
                centers = centers_cache[bed_key]
            else:
                centers = _read_bed_centers(str(bed))
                if centers_cache is not None:
                    centers_cache[bed_key] = centers
            if not centers:
                continue
            for center in centers:
                if center not in seen_centers:
                    seen_centers.add(center)
                    all_centers.append(center)
            condition_centers.extend(centers)
            profile_key = (str(sample.corrected_bigwig), prefix)
            if profile_cache is not None and profile_key in profile_cache:
                profile = profile_cache[profile_key]
            else:
                profile = _mean_profile(str(sample.corrected_bigwig), centers, flank)
                if profile_cache is not None:
                    profile_cache[profile_key] = profile
            condition_profiles.append(profile)
            samples.append({"name": sample.sample, "profile": profile, "fp_score": round(float(_aggregate_fp_score(profile)), 6), "source": "recomputed"})
        if condition_profiles:
            mean_profile = [round(float(v), 6) for v in np.nanmean(np.asarray(condition_profiles, dtype=float), axis=0)]
            fp_score = round(float(_aggregate_fp_score(mean_profile)), 6)
        else:
            mean_profile = [0.0] * (flank * 2)
            fp_score = 0.0
        conditions.append({"name": condition, "profile": mean_profile, "samples": samples, "n_sites": len(condition_centers), "fp_score": fp_score})
    if not all_centers or not any(condition.get("samples") for condition in conditions):
        return None
    return {
        "prefix": prefix,
        "name": str(point.get("name") or prefix),
        "motif_id": str(point.get("motif_id") or ""),
        "change": _as_float(point.get("change")),
        "pvalue": _as_float(point.get("pvalue"), 1.0),
        "fdr": point.get("fdr", ""),
        "n_sites": len(all_centers),
        "site_set": "all",
        "x": list(range(-flank, flank)),
        "conditions": conditions,
        "profile_source": "recomputed",
    }


def _wanted_conditions(comparison: ComparisonView) -> list[str]:
    wanted_conditions = list((comparison.payload.get("conditions") or [])[:2])
    if len(wanted_conditions) < 2:
        wanted_conditions = [comparison.condition, "10_FBS_Ctrl"]
    return wanted_conditions


def _compute_condition_profile_worker(task: tuple[str, str, list[tuple[str, str, str]], int]) -> tuple[str, str, dict | None]:
    prefix, condition, sample_rows, flank = task
    samples = []
    condition_profiles = []
    condition_centers = []
    seen_centers = set()
    for sample_name, corrected_bigwig, match_dir in sample_rows:
        bed = _motif_all_bed(Path(match_dir), prefix)
        if bed is None or not Path(corrected_bigwig).exists():
            continue
        centers = _read_bed_centers(str(bed))
        if not centers:
            continue
        unique_centers = []
        for center in centers:
            if center not in seen_centers:
                seen_centers.add(center)
                unique_centers.append(center)
        condition_centers.extend(unique_centers)
        profile = _mean_profile(str(corrected_bigwig), unique_centers, flank)
        condition_profiles.append(profile)
        samples.append({"name": sample_name, "profile": profile, "fp_score": round(float(_aggregate_fp_score(profile)), 6), "source": "recomputed"})
    if not condition_profiles:
        return prefix, condition, None
    mean_profile = [round(float(v), 6) for v in np.nanmean(np.asarray(condition_profiles, dtype=float), axis=0)]
    record = {
        "name": condition,
        "profile": mean_profile,
        "samples": samples,
        "n_sites": len(condition_centers),
        "fp_score": round(float(_aggregate_fp_score(mean_profile)), 6),
    }
    return prefix, condition, record


def _compute_missing_condition_profiles(
    missing_pairs: set[tuple[str, str]],
    project_samples: list[SampleView],
    flank: int,
    cores: int | None = None,
) -> dict[tuple[str, str], dict]:
    if not missing_pairs or not project_samples:
        return {}
    by_condition: dict[str, list[tuple[str, str, str]]] = {}
    for sample in project_samples:
        by_condition.setdefault(sample.condition, []).append((sample.sample, str(sample.corrected_bigwig), str(sample.match_dir)))
    tasks = []
    for prefix, condition in sorted(missing_pairs):
        sample_rows = by_condition.get(condition) or []
        if sample_rows:
            tasks.append((prefix, condition, sample_rows, flank))
    if not tasks:
        return {}
    workers = cores or (os.cpu_count() or 1)
    workers = max(1, min(workers, len(tasks)))
    out: dict[tuple[str, str], dict] = {}
    if workers == 1:
        iterator = map(_compute_condition_profile_worker, tasks)
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as pool:
            iterator = pool.map(_compute_condition_profile_worker, tasks)
            for prefix, condition, record in iterator:
                if record:
                    out[(prefix, condition)] = record
            return out
    for prefix, condition, record in iterator:
        if record:
            out[(prefix, condition)] = record
    return out


def _assemble_recomputed_aggregate(prefix: str, point: dict, conditions: list[str], condition_profiles: dict[tuple[str, str], dict], flank: int) -> dict | None:
    records = []
    for condition in conditions:
        record = condition_profiles.get((prefix, condition))
        if not record:
            return None
        records.append(record)
    return {
        "prefix": prefix,
        "name": str(point.get("name") or prefix),
        "motif_id": str(point.get("motif_id") or ""),
        "change": _as_float(point.get("change")),
        "pvalue": _as_float(point.get("pvalue"), 1.0),
        "fdr": point.get("fdr", ""),
        "n_sites": sum(int(record.get("n_sites") or 0) for record in records),
        "site_set": "all",
        "x": list(range(-flank, flank)),
        "conditions": records,
        "profile_source": "recomputed",
    }


def prepare_aggregate_maps(
    review_payload: dict,
    project: str | Path | None = None,
    fill_missing: bool = False,
    recompute_missing: bool = False,
    flank: int = 60,
    cores: int | None = None,
) -> dict[tuple[int, str], tuple[dict | None, str]]:
    project_samples = _read_project_samples(project) if project and recompute_missing else []
    condition_cache = _profile_cache_from_payload(review_payload) if fill_missing else {}
    prepared = {}
    pending: list[tuple[ComparisonView, str, dict, list[str]]] = []
    missing_condition_pairs: set[tuple[str, str]] = set()
    for comparison in ordered_comparisons(review_payload):
        points = _point_map(comparison.payload)
        aggregates = _aggregate_map(comparison.payload)
        for prefix, point in points.items():
            aggregate = aggregates.get(prefix)
            if aggregate:
                prepared[(comparison.index, prefix)] = (aggregate, "html")
                continue
            assembled = _assemble_cached_aggregate(prefix, point, comparison, condition_cache, flank) if condition_cache else None
            if assembled:
                prepared[(comparison.index, prefix)] = (assembled, "assembled")
                continue
            conditions = _wanted_conditions(comparison)
            if project_samples:
                pending.append((comparison, prefix, point, conditions))
                missing_condition_pairs.update((prefix, condition) for condition in conditions)
            else:
                prepared[(comparison.index, prefix)] = (None, "missing")
    recomputed_profiles = _compute_missing_condition_profiles(missing_condition_pairs, project_samples, flank, cores=cores) if pending else {}
    for comparison, prefix, point, conditions in pending:
        recomputed = _assemble_recomputed_aggregate(prefix, point, conditions, recomputed_profiles, flank)
        prepared[(comparison.index, prefix)] = (recomputed, "recomputed" if recomputed else "missing")
    return prepared


def _nice_ticks(vmin: float, vmax: float, n: int = 4) -> list[float]:
    if not math.isfinite(vmin) or not math.isfinite(vmax) or vmin == vmax:
        return [0.0]
    span = vmax - vmin
    raw = span / max(1, n - 1)
    power = 10 ** math.floor(math.log10(abs(raw)))
    frac = raw / power
    step = (1 if frac <= 1 else 2 if frac <= 2 else 2.5 if frac <= 2.5 else 5 if frac <= 5 else 10) * power
    start = math.ceil(vmin / step) * step
    end = math.floor(vmax / step) * step
    ticks = []
    value = start
    while value <= end + step / 2:
        ticks.append(float(f"{value:.12g}"))
        value += step
    return ticks or [vmin, vmax]


def _format_number(value: float) -> str:
    if not math.isfinite(value):
        return "NA"
    av = abs(value)
    if av >= 1:
        return f"{value:.2f}"
    if av >= 0.01:
        return f"{value:.3f}"
    if av >= 0.001:
        return f"{value:.4f}"
    return f"{value:.1e}"


def _add_fig_line(fig, xs: list[float], ys: list[float], **kwargs) -> None:
    fig.add_artist(Line2D(xs, ys, transform=fig.transFigure, **kwargs))


def _add_fig_text(fig, x: float, y: float, text: str, **kwargs) -> None:
    fig.text(x, y, text, transform=fig.transFigure, **kwargs)


def _fig_rect(fig, x0: float, y0: float, width: float, height: float, **kwargs) -> None:
    fig.add_artist(Rectangle((x0, y0), width, height, transform=fig.transFigure, **kwargs))


def _add_trace_legend(fig, comparisons: list[ComparisonView], page_width: float, page_height: float) -> None:
    if not comparisons:
        return
    first = comparisons[0]
    conditions = list((first.payload.get("conditions") or [])[:2])
    stress_label = "Column condition"
    ctrl_label = conditions[1] if len(conditions) > 1 else "Control"
    stress_color = _condition_color(first.payload, conditions[0] if conditions else first.condition, 0)
    ctrl_color = _condition_color(first.payload, ctrl_label, 1)
    x0 = 0.70
    y = 1 - 0.18 / page_height
    line_w = 0.026
    text_dx = 0.006
    _add_fig_line(fig, [x0, x0 + line_w], [y, y], color=stress_color, lw=1.4)
    _add_fig_text(fig, x0 + line_w + text_dx, y, f"{stress_label} (nutrient stress)", ha="left", va="center", fontsize=5.8, fontweight="bold")
    x1 = x0 + 0.17
    _add_fig_line(fig, [x1, x1 + line_w], [y, y], color=ctrl_color, lw=1.4)
    _add_fig_text(fig, x1 + line_w + text_dx, y, ctrl_label, ha="left", va="center", fontsize=5.8, fontweight="bold")


def _domain_for_row(comparisons: list[ComparisonView], motif: MotifView, aggregate_maps: dict[tuple[int, str], tuple[dict | None, str]], flank: int) -> tuple[float, float]:
    values = []
    for comparison in comparisons:
        aggregate = aggregate_maps.get((comparison.index, motif.prefix), (None, "missing"))[0]
        if not aggregate:
            continue
        for row in _samples_for_motif(comparison.payload, aggregate, flank):
            values.extend(row["y"])
    if not values:
        return 0.0, 1.0
    vmin = min(values)
    vmax = max(values)
    pad = max((vmax - vmin) * 0.06, 1e-6)
    return vmin - pad, vmax + pad


def write_source_table(
    review_payload: dict,
    output: str | Path,
    flank: int = 60,
    motif_order: list[MotifView] | None = None,
    aggregate_maps: dict[tuple[int, str], tuple[dict | None, str]] | None = None,
    rna_fc_map: dict[tuple[int, str], RnaFcView] | None = None,
) -> int:
    comparisons = ordered_comparisons(review_payload)
    motifs = _ordered_motifs_for_payload(review_payload, motif_order)
    aggregate_maps = aggregate_maps or prepare_aggregate_maps(review_payload, flank=flank)
    rna_fc_map = rna_fc_map or {}
    rows = []
    for motif in motifs:
        for comparison in comparisons:
            points = _point_map(comparison.payload)
            point = points.get(motif.prefix) or {}
            aggregate, profile_source = aggregate_maps.get((comparison.index, motif.prefix), ({}, "missing"))
            aggregate = aggregate or {}
            rna = rna_fc_map.get((comparison.index, motif.prefix))
            rows.append(
                {
                    "motif_prefix": motif.prefix,
                    "motif_name": motif.name,
                    "motif_id": motif.motif_id,
                    "comparison": comparison.label,
                    "condition": comparison.condition,
                    "delta_fp": point.get("change", aggregate.get("change", "")),
                    "pvalue": point.get("pvalue", aggregate.get("pvalue", "")),
                    "fdr": point.get("fdr", ""),
                    "n_sites": aggregate.get("n_sites", aggregate.get("sites", "")),
                    "aggregate_profile": bool(aggregate),
                    "profile_source": profile_source,
                    "rna_tf_genes": ";".join(gene for gene, _ in rna.values) if rna else "",
                    "rna_log2fc": ";".join(f"{gene}:{value:.6g}" for gene, value in rna.values) if rna else "",
                }
            )
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()) if rows else ["motif_prefix"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def plot_grid_pdf(
    review_payload: dict,
    output: str | Path,
    rows_per_page: int = 16,
    flank: int = 60,
    title: str | None = None,
    motif_order: list[MotifView] | None = None,
    aggregate_maps: dict[tuple[int, str], tuple[dict | None, str]] | None = None,
    rna_fc_map: dict[tuple[int, str], RnaFcView] | None = None,
    repeat_column_labels: str = "none",
) -> tuple[int, int]:
    comparisons = ordered_comparisons(review_payload)
    motifs = _ordered_motifs_for_payload(review_payload, motif_order)
    aggregate_maps = aggregate_maps or prepare_aggregate_maps(review_payload, flank=flank)
    rna_fc_map = rna_fc_map or {}
    if repeat_column_labels not in {"none", "row"}:
        raise ValueError("repeat_column_labels must be 'none' or 'row'")
    if not comparisons:
        raise ValueError("No comparisons were found in the review payload")
    if not motifs:
        raise ValueError("No aggregate motif profiles were found in the review payload")
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"], "font.size": 6, "axes.linewidth": 0.6, "pdf.fonttype": 42, "ps.fonttype": 42})
    n_cols = len(comparisons)
    panel = 1.0
    gap = 0.045
    label_width = 1.65
    left_margin = 0.08
    right_margin = 0.08
    top_margin = 0.42
    header_height = 0.34
    bottom_margin = 0.18
    plot_left = 0.22
    plot_right = 0.045
    row_column_labels = repeat_column_labels == "row"
    plot_top = 0.16 if row_column_labels else 0.055
    plot_bottom = 0.16
    page_width = left_margin + label_width + n_cols * panel + (n_cols - 1) * gap + right_margin
    page_height = top_margin + header_height + rows_per_page * panel + (rows_per_page - 1) * gap + bottom_margin
    with PdfPages(output) as pdf:
        for page_idx, start in enumerate(range(0, len(motifs), rows_per_page), start=1):
            page_motifs = motifs[start : start + rows_per_page]
            fig = plt.figure(figsize=(page_width, page_height), constrained_layout=False)
            _add_fig_text(fig, left_margin / page_width, 1 - 0.18 / page_height, title or str(review_payload.get("title") or "Motif aggregate grid"), ha="left", va="center", fontsize=8, fontweight="bold")
            _add_fig_text(fig, 0.5, 1 - 0.18 / page_height, "Corrected cut-site signal; all motif sites; global order by max |dFP|" if motif_order else "Corrected cut-site signal; all motif sites; order by max |dFP|", ha="center", va="center", fontsize=6, fontweight="bold")
            _add_trace_legend(fig, comparisons, page_width, page_height)
            _add_fig_text(fig, left_margin / page_width, 1 - (top_margin + header_height * 0.45) / page_height, f"Page {page_idx}", ha="left", va="center", fontsize=6, fontweight="bold")
            header_y = page_height - top_margin - header_height * 0.55
            for col_idx, comparison in enumerate(comparisons):
                panel_x = left_margin + label_width + col_idx * (panel + gap)
                _add_fig_text(fig, (panel_x + panel / 2) / page_width, header_y / page_height, comparison.condition, ha="center", va="center", fontsize=5.7, fontweight="bold")
            for row_idx, motif in enumerate(page_motifs, start=1):
                vmin, vmax = _domain_for_row(comparisons, motif, aggregate_maps, flank)
                ticks = _nice_ticks(vmin, vmax, 4)
                row_top = page_height - top_margin - header_height - (row_idx - 1) * (panel + gap)
                panel_y = row_top - panel
                _add_fig_text(fig, (left_margin + label_width - 0.05) / page_width, (panel_y + panel / 2) / page_height, _motif_label(motif), ha="right", va="center", fontsize=5.6, fontweight="bold")
                for col_idx, comparison in enumerate(comparisons):
                    panel_x = left_margin + label_width + col_idx * (panel + gap)
                    if row_column_labels:
                        label_y = panel_y + panel - 0.055
                        _fig_rect(
                            fig,
                            (panel_x + 0.055) / page_width,
                            (label_y - 0.035) / page_height,
                            (panel - 0.11) / page_width,
                            0.07 / page_height,
                            facecolor="white",
                            edgecolor="#d7dee8",
                            linewidth=0.28,
                        )
                        _add_fig_text(
                            fig,
                            (panel_x + panel / 2) / page_width,
                            label_y / page_height,
                            comparison.condition,
                            ha="center",
                            va="center",
                            fontsize=5.2,
                            fontweight="bold",
                            color="#111827",
                        )
                    plot_x0 = panel_x + plot_left
                    plot_y0 = panel_y + plot_bottom
                    plot_w = panel - plot_left - plot_right
                    plot_h = panel - plot_top - plot_bottom
                    fx0 = plot_x0 / page_width
                    fy0 = plot_y0 / page_height
                    fw = plot_w / page_width
                    fh = plot_h / page_height
                    aggregate, profile_source = aggregate_maps.get((comparison.index, motif.prefix), (None, "missing"))
                    point = _point_map(comparison.payload).get(motif.prefix) or aggregate or {}
                    change = _as_float(point.get("change"), default=float("nan"))
                    rna = rna_fc_map.get((comparison.index, motif.prefix))
                    _fig_rect(fig, fx0, fy0, fw, fh, facecolor="white", edgecolor="black", linewidth=0.45)
                    sx = lambda value: (plot_x0 + ((value + flank) / (2 * flank)) * plot_w) / page_width
                    sy = lambda value: (plot_y0 + ((value - vmin) / (vmax - vmin or 1.0)) * plot_h) / page_height
                    if aggregate:
                        for sample in _samples_for_motif(comparison.payload, aggregate, flank):
                            _add_fig_line(fig, [sx(x) for x in sample["x"]], [sy(y) for y in sample["y"]], color=sample["color"], lw=0.45, alpha=0.9)
                    else:
                        _add_fig_text(fig, (plot_x0 + plot_w / 2) / page_width, (plot_y0 + plot_h / 2) / page_height, "No\nprofile", ha="center", va="center", fontsize=5.5, fontweight="bold")
                    _add_fig_text(fig, (plot_x0 + 0.02) / page_width, (plot_y0 + plot_h - 0.035) / page_height, f"dFP={_format_number(change)}", ha="left", va="top", fontsize=5.2, fontweight="bold", bbox={"boxstyle": "round,pad=0.08", "facecolor": "white", "edgecolor": "none", "alpha": 0.82})
                    if rna and rna.label:
                        _add_fig_text(fig, (plot_x0 + 0.02) / page_width, (plot_y0 + plot_h - 0.13) / page_height, rna.label, ha="left", va="top", fontsize=4.5, fontweight="bold", color="#111827", bbox={"boxstyle": "round,pad=0.06", "facecolor": "white", "edgecolor": "none", "alpha": 0.75})
                    for x_tick in (-flank, 0, flank):
                        x_fig = sx(x_tick)
                        _add_fig_line(fig, [x_fig, x_fig], [fy0, fy0 + fh], color="#e6edf5" if x_tick else "#999999", lw=0.3, ls=(0, (2, 2)) if x_tick == 0 else "solid")
                        _add_fig_text(fig, x_fig, (plot_y0 - 0.035) / page_height, str(x_tick), ha="center", va="top", fontsize=4.8)
                    for y_tick in ticks:
                        y_fig = sy(y_tick)
                        _add_fig_line(fig, [fx0, fx0 + fw], [y_fig, y_fig], color="#edf2f7", lw=0.25)
                        if col_idx == 0:
                            _add_fig_text(fig, (plot_x0 - 0.025) / page_width, y_fig, _format_number(y_tick), ha="right", va="center", fontsize=4.8)
            pdf.savefig(fig)
            plt.close(fig)
    return len(motifs), math.ceil(len(motifs) / rows_per_page)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="plot-motif-aggregate-grid", description="Create a multi-page motif-by-comparison aggregate PDF from review-multi-comparisons HTML.")
    parser.add_argument("--input-html", help="review_multi_comparisons.html with embedded comparison payloads.")
    parser.add_argument("--order-htmls", nargs="*", default=[], help="Optional review HTML files used only to build one shared motif row order by max |delta FP|.")
    parser.add_argument("--outdir", help="Project directory used with --layout project.")
    parser.add_argument("--layout", choices=["project", "custom"], default="project", help="Use fp-tools project layout under --outdir (default: project).")
    parser.add_argument("--output", help="Output multi-page PDF. In project mode, defaults to reports/motif_aggregate_grid.pdf.")
    parser.add_argument("--source-tsv", help="Output source TSV. Defaults to <output stem>_source.tsv.")
    parser.add_argument("--rows-per-page", type=int, default=16, help="Motif rows per PDF page (default: 16).")
    parser.add_argument("--flank", type=int, default=60, help="Distance from motif center shown in each subplot (default: 60 bp).")
    parser.add_argument("--fill-missing-profiles", action="store_true", help="Fill missing motif aggregate panels from condition profiles already embedded elsewhere in the review report.")
    parser.add_argument("--recompute-missing-profiles", action="store_true", help="Slower fallback: recompute still-missing motif aggregate profiles from project sample bigWigs and motif BEDs.")
    parser.add_argument("--cores", type=int, default=None, help="Worker processes for --recompute-missing-profiles (default: all available cores).")
    parser.add_argument("--rna-log2norm", help="Optional DESeq2/RUVr log2-normalized RNA matrix with gene_key and sample columns.")
    parser.add_argument("--rna-raw-counts", help="Optional raw RNA count matrix used to filter unexpressed TFs.")
    parser.add_argument("--motif-gene-map", help="Optional motif-to-gene map with motif and gene_symbol columns.")
    parser.add_argument("--rna-min-raw-mean", type=float, default=1.0, help="Minimum mean raw count in either compared condition for a TF to be shown (default: 1.0).")
    parser.add_argument("--repeat-column-labels", choices=["none", "row"], default="none", help="Repeat comparison labels inside each motif row panel instead of showing them only in the page header (default: none).")
    parser.add_argument("--title", help="PDF title.")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    project = None
    if is_project_layout(args.layout):
        if not args.outdir and not args.input_html:
            parser.error("--layout project requires --outdir or --input-html")
        if args.outdir:
            project = project_root(args.outdir)
            if not args.input_html:
                args.input_html = str(reports_dir(project) / "review_multi_comparisons.html")
            if not args.output:
                args.output = str(reports_dir(project) / "motif_aggregate_grid.pdf")
    if not args.input_html:
        parser.error("provide --input-html or use --layout project with --outdir")
    if not args.output:
        parser.error("provide --output or use --layout project with --outdir")
    if args.recompute_missing_profiles and project is None:
        parser.error("--recompute-missing-profiles requires --outdir in project layout")
    if args.rows_per_page < 1:
        parser.error("--rows-per-page must be at least 1")
    if args.rna_log2norm and not args.motif_gene_map:
        parser.error("--rna-log2norm requires --motif-gene-map")
    if args.rna_raw_counts and not args.rna_log2norm:
        parser.error("--rna-raw-counts requires --rna-log2norm")
    payload = read_embedded_payload(args.input_html)
    if payload.get("schema") != "fp-tools.review-multi-comparisons.v1":
        raise SystemExit(f"{args.input_html} is not a review-multi-comparisons HTML payload")
    order_payloads = []
    for html in args.order_htmls:
        order_payload = read_embedded_payload(html)
        if order_payload.get("schema") != "fp-tools.review-multi-comparisons.v1":
            raise SystemExit(f"{html} is not a review-multi-comparisons HTML payload")
        order_payloads.append(order_payload)
    motif_order = collect_motifs_from_payloads(order_payloads) if order_payloads else None
    aggregate_maps = prepare_aggregate_maps(payload, project=project, fill_missing=args.fill_missing_profiles or args.recompute_missing_profiles, recompute_missing=args.recompute_missing_profiles, flank=args.flank, cores=args.cores)
    rna_fc_map = compute_rna_fc_map(payload, project, motif_order, args.rna_log2norm, args.rna_raw_counts, args.motif_gene_map, min_raw_mean=args.rna_min_raw_mean)
    output = Path(args.output)
    source_tsv = Path(args.source_tsv) if args.source_tsv else output.with_name(f"{output.stem}_source.tsv")
    motifs, pages = plot_grid_pdf(payload, output, rows_per_page=args.rows_per_page, flank=args.flank, title=args.title, motif_order=motif_order, aggregate_maps=aggregate_maps, rna_fc_map=rna_fc_map, repeat_column_labels=args.repeat_column_labels)
    row_count = write_source_table(payload, source_tsv, flank=args.flank, motif_order=motif_order, aggregate_maps=aggregate_maps, rna_fc_map=rna_fc_map)
    print(f"Wrote {output} ({motifs} motifs, {pages} pages)")
    print(f"Wrote {source_tsv} ({row_count} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
