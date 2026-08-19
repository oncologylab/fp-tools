"""Interactive aggregate HTML reports for plot-aggregate."""

from __future__ import annotations

import argparse
import base64
import csv
import gzip
import html
import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
from fp_tools.utils import bigwig as pyBigWig


DEFAULT_COLORS = ["#2563eb", "#dc2626", "#16a34a", "#9333ea", "#f97316", "#0891b2", "#7c3aed", "#64748b"]


def _read_manifest(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [dict(row) for row in reader]
    fields = set(reader.fieldnames or [])
    required = {"sample", "signal"}
    missing = required - fields
    if "match_dir" not in fields and "sample_dir" not in fields:
        missing.add("match_dir/sample_dir")
    if missing:
        raise SystemExit(f"Manifest is missing required column(s): {', '.join(sorted(missing))}")
    for row in rows:
        if not row.get("match_dir") and row.get("sample_dir"):
            row["match_dir"] = row["sample_dir"]
        row["_groups_defined"] = "condition" in fields and bool(row.get("condition"))
    return rows


def _read_bed_centers(path: Path) -> list[tuple[str, int]]:
    centers = []
    if not path.exists():
        return centers
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            try:
                centers.append((fields[0], (int(fields[1]) + int(fields[2])) // 2))
            except ValueError:
                continue
    return centers


def _read_bed_sites(path: Path) -> set[tuple[str, int, int]]:
    sites = set()
    if not path.exists():
        return sites
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            try:
                sites.add((fields[0], int(fields[1]), int(fields[2])))
            except ValueError:
                continue
    return sites


def _motif_bed_path(match_dir: str | Path, prefix: str, condition: str | None = None, sample: str | None = None, site_set: str = "bound") -> Path:
    """Return the requested motif-site BED, with sensible fallbacks."""
    bed_dir = Path(match_dir) / prefix / "beds"
    candidates = []
    site_set = (site_set or "bound").lower()
    if site_set in {"bound", "unbound"}:
        for label in (condition, sample):
            if label:
                candidates.append(bed_dir / f"{prefix}_{label}_{site_set}.bed")
        candidates.extend(sorted(bed_dir.glob(f"{prefix}_*_{site_set}.bed")))
    if site_set == "all":
        candidates.append(bed_dir / f"{prefix}_all.bed")
    candidates.append(bed_dir / f"{prefix}_all.bed")
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return bed_dir / f"{prefix}_all.bed"


def _bed_source_label(path: str | Path | None) -> str:
    if path is None:
        return "motif-site set"
    name = Path(path).name
    if name.endswith("_bound.bed"):
        return "bound.bed"
    if name.endswith("_all.bed"):
        return "all.bed"
    if name.endswith("_unbound.bed"):
        return "unbound.bed"
    return name or "motif-site set"


def _summarize_site_sets(values: list[str]) -> str:
    clean = sorted({v for v in values if v})
    if not clean:
        return "motif-site set"
    if len(clean) == 1:
        return clean[0]
    return "mixed beds"


def _discover_motifs(match_dir: str | Path) -> list[dict[str, str | int | float]]:
    root = Path(match_dir)
    result_files = sorted(root.glob("*_results.txt"))
    motifs = []
    if result_files:
        with result_files[0].open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                prefix = row.get("output_prefix") or row.get("name") or row.get("motif_id")
                if not prefix:
                    continue
                total = float(row.get("total_tfbs") or 0)
                motifs.append({
                    "prefix": prefix,
                    "name": row.get("name") or prefix,
                    "motif_id": row.get("motif_id") or "",
                    "score": total,
                    "sites": int(total),
                })
    if not motifs:
        for bed in sorted(root.glob("*/beds/*_all.bed")):
            prefix = bed.parent.parent.name
            centers = _read_bed_centers(bed)
            motifs.append({"prefix": prefix, "name": prefix, "motif_id": "", "score": len(centers), "sites": len(centers)})
    motifs.sort(key=lambda row: (-float(row["score"]), str(row["name"])))
    return motifs


def _mean_profile(signal: str | Path, centers: list[tuple[str, int]], flank: int) -> list[float]:
    profiles = []
    with pyBigWig.open(str(signal)) as bw:
        chroms = bw.chroms()
        for chrom, center in centers:
            if chrom not in chroms:
                continue
            start = center - flank
            end = center + flank
            if start < 0 or end > chroms[chrom] or end <= start:
                continue
            values = np.asarray(bw.values(chrom, start, end, numpy=True), dtype=float)
            if values.size != flank * 2:
                continue
            profiles.append(np.nan_to_num(values, nan=0.0, posinf=0.0, neginf=0.0))
    if not profiles:
        return [0.0] * (flank * 2)
    return [round(float(v), 6) for v in np.nanmean(np.vstack(profiles), axis=0)]


def _x_values(flank: int) -> list[int]:
    return list(range(-flank, 0)) + list(range(1, flank + 1))


def _condition_colors(conditions: list[str]) -> dict[str, str]:
    return {cond: DEFAULT_COLORS[idx % len(DEFAULT_COLORS)] for idx, cond in enumerate(conditions)}


def _normalize_profile(profile: list[float], mode: str, baseline: float | None = None, scale: float | None = None) -> list[float]:
    arr = np.asarray(profile, dtype=float)
    if mode == "none" or arr.size == 0:
        return [round(float(v), 6) for v in arr]
    center = float(np.nanmedian(arr)) if baseline is None else baseline
    spread = float(np.nanpercentile(arr, 90) - np.nanpercentile(arr, 10)) if scale is None else scale
    if not np.isfinite(spread) or abs(spread) < 1e-12:
        spread = 1.0
    return [round(float(v), 6) for v in ((arr - center) / spread)]


def _profile_mean(profiles: list[list[float]]) -> list[float]:
    if not profiles:
        return []
    return [round(float(v), 6) for v in np.nanmean(np.asarray(profiles, dtype=float), axis=0)]


def _select_motif_prefixes(motifs: list[dict[str, str | int | float]], requested: list[str] | None, top_n: int) -> list[str]:
    if not requested:
        return [str(motif["prefix"]) for motif in motifs[:top_n]]
    lookup = {}
    for motif in motifs:
        for key in ("prefix", "name", "motif_id"):
            value = str(motif.get(key) or "").lower()
            if value:
                lookup[value] = str(motif["prefix"])
    selected = []
    missing = []
    for value in requested:
        prefix = lookup.get(str(value).lower())
        if prefix is None:
            missing.append(str(value))
            continue
        if prefix not in selected:
            selected.append(prefix)
    if missing:
        raise SystemExit(f"Requested motif(s) not found in match directory: {', '.join(missing)}")
    return selected


def _logo_data_uri(path: Path) -> str:
    if not path.exists():
        return ""
    return "data:image/png;base64," + base64.b64encode(path.read_bytes()).decode("ascii")


def _logo_map_from_match_dirs(match_dirs: list[str | Path], prefixes: list[str]) -> dict[str, dict[str, str]]:
    logos = {}
    for prefix in prefixes:
        for match_dir in match_dirs:
            png = Path(match_dir) / prefix / f"{prefix}.png"
            if png.exists():
                logos[prefix] = {"png": _logo_data_uri(png)}
                break
    return logos


def build_payload(rows: list[dict[str, str]], flank: int, top_n: int, normalization: str = "none", motif_names: list[str] | None = None, site_set: str = "bound") -> dict:
    """Build the compact aggregate payload from manifest/sample-dir rows."""

    x = _x_values(flank)
    normalization = (normalization or "none").replace("_", "-")
    sample_rows = []
    motif_meta: dict[str, dict[str, str | int | float]] = {}
    motif_scores: dict[str, float] = defaultdict(float)

    for idx, row in enumerate(rows):
        sample = row["sample"]
        condition = row.get("condition") or sample
        label = row.get("label") or sample
        sample_rows.append({"sample": sample, "condition": condition, "label": label, "row": row, "idx": idx})
        for motif in _discover_motifs(row["match_dir"]):
            prefix = str(motif["prefix"])
            motif_meta.setdefault(prefix, motif)
            motif_scores[prefix] = max(motif_scores[prefix], float(motif.get("score") or 0.0))

    ranked_prefixes = sorted(motif_scores, key=lambda pfx: (-motif_scores[pfx], str(motif_meta[pfx].get("name", pfx))))
    ranked_motifs = [motif_meta[prefix] for prefix in ranked_prefixes]
    selected_prefixes = _select_motif_prefixes(ranked_motifs, motif_names, top_n)
    raw_profiles: dict[tuple[str, str], list[float]] = {}
    profile_values_for_norm = []
    for sample_info in sample_rows:
        row = sample_info["row"]
        sample = sample_info["sample"]
        for prefix in selected_prefixes:
            bed = _motif_bed_path(row["match_dir"], prefix, condition=sample_info["condition"], sample=sample, site_set=site_set)
            centers = _read_bed_centers(bed)
            profile = _mean_profile(row["signal"], centers, flank)
            raw_profiles[(sample, prefix)] = profile
            sample_info.setdefault("site_counts", {})[prefix] = len(centers)
            sample_info.setdefault("site_sets", {})[prefix] = _read_bed_sites(bed)
            sample_info.setdefault("bed_sources", {})[prefix] = _bed_source_label(bed)
            profile_values_for_norm.extend(profile)

    global_baseline = global_scale = None
    if normalization == "sample-quantile" and profile_values_for_norm:
        arr = np.asarray(profile_values_for_norm, dtype=float)
        global_baseline = float(np.nanmedian(arr))
        global_scale = float(np.nanpercentile(arr, 90) - np.nanpercentile(arr, 10))

    conditions = []
    for sample_info in sample_rows:
        if sample_info["condition"] not in conditions:
            conditions.append(sample_info["condition"])
    colors = _condition_colors(conditions)
    groups_defined = any(bool(sample_info["row"].get("_groups_defined")) for sample_info in sample_rows)

    motifs = []
    for prefix in selected_prefixes:
        meta = motif_meta[prefix]
        series = []
        condition_profiles: dict[str, list[list[float]]] = defaultdict(list)
        for sample_info in sample_rows:
            sample = sample_info["sample"]
            condition = sample_info["condition"]
            profile = raw_profiles.get((sample, prefix), [0.0] * len(x))
            if normalization == "condition-quantile":
                cond_values = []
                for other in sample_rows:
                    if other["condition"] == condition:
                        cond_values.extend(raw_profiles.get((other["sample"], prefix), []))
                arr = np.asarray(cond_values, dtype=float) if cond_values else np.asarray(profile, dtype=float)
                baseline = float(np.nanmedian(arr))
                scale = float(np.nanpercentile(arr, 90) - np.nanpercentile(arr, 10))
                profile = _normalize_profile(profile, normalization, baseline, scale)
            else:
                profile = _normalize_profile(profile, normalization, global_baseline, global_scale)
            condition_profiles[condition].append(profile)
            avg_score = float(np.nanmean(np.asarray(profile, dtype=float))) if profile else 0.0
            series.append({
                "id": f"sample::{sample}",
                "label": sample_info["label"],
                "kind": "sample",
                "condition": condition,
                "profile": profile,
                "avg_score": round(avg_score, 6),
                "sites": int(sample_info.get("site_counts", {}).get(prefix, 0)),
                "bed_source": str(sample_info.get("bed_sources", {}).get(prefix, "motif-site set")),
            })
        for condition, profiles in condition_profiles.items():
            mean_profile = _profile_mean(profiles)
            avg_score = float(np.nanmean(np.asarray(mean_profile, dtype=float))) if mean_profile else 0.0
            condition_site_coords = set()
            for sample_info in sample_rows:
                if sample_info["condition"] == condition:
                    condition_site_coords.update(sample_info.get("site_sets", {}).get(prefix, set()))
            condition_sites = len(condition_site_coords)
            condition_site_set = _summarize_site_sets([str(s.get("bed_source") or "") for s in series if s.get("kind") == "sample" and s.get("condition") == condition])
            series.append({
                "id": f"condition::{condition}",
                "label": f"{condition} mean",
                "kind": "condition",
                "condition": condition,
                "profile": mean_profile,
                "avg_score": round(avg_score, 6),
                "sites": condition_sites,
                "bed_source": condition_site_set,
            })
        motif_site_set = _summarize_site_sets([str(s.get("bed_source") or "") for s in series])
        union_site_coords = set()
        for sample_info in sample_rows:
            union_site_coords.update(sample_info.get("site_sets", {}).get(prefix, set()))
        union_sites = len(union_site_coords) if union_site_coords else int(meta.get("sites") or 0)
        motifs.append({"prefix": prefix, "name": str(meta.get("name") or prefix), "motif_id": str(meta.get("motif_id") or ""), "score": round(float(motif_scores[prefix]), 6), "sites": union_sites, "site_set": motif_site_set, "series": series})

    logos = _logo_map_from_match_dirs([row["match_dir"] for row in rows], selected_prefixes)
    return {"schema": "fp-tools.aggregate.batch.v2", "x": x, "motifs": motifs, "conditions": conditions, "colors": colors, "logos": logos, "groups_defined": groups_defined, "normalization": normalization, "site_set": _summarize_site_sets([str(m.get("site_set") or "") for m in motifs]), "x_label": "Distance from motif center (bp)", "y_label": "Corrected cut-site signal (a.u.)" if normalization == "none" else "Normalized corrected cut-site signal (a.u.)"}


def build_payload_from_tfbs(tfbs: list[str | Path], signals: list[str | Path], labels: list[str], conditions: list[str], flank: int, normalization: str = "none", motif_labels: list[str] | None = None, groups_defined: bool = False) -> dict:
    """Build an aggregate HTML payload directly from BED files and bigWigs."""

    x = _x_values(flank)
    normalization = (normalization or "none").replace("_", "-")
    motifs = []
    colors = _condition_colors(list(dict.fromkeys(conditions)))
    for idx, bed in enumerate(tfbs):
        path = Path(bed)
        prefix = path.stem
        name = motif_labels[idx] if motif_labels and idx < len(motif_labels) else prefix
        centers = _read_bed_centers(path)
        raw_profiles = [_mean_profile(signal, centers, flank) for signal in signals]
        flat = [value for profile in raw_profiles for value in profile]
        baseline = scale = None
        if normalization == "sample-quantile" and flat:
            arr = np.asarray(flat, dtype=float)
            baseline = float(np.nanmedian(arr))
            scale = float(np.nanpercentile(arr, 90) - np.nanpercentile(arr, 10))
        series = []
        condition_profiles: dict[str, list[list[float]]] = defaultdict(list)
        for sample_label, condition, profile in zip(labels, conditions, raw_profiles):
            if normalization == "condition-quantile":
                cond_values = [v for other_cond, other_profile in zip(conditions, raw_profiles) if other_cond == condition for v in other_profile]
                arr = np.asarray(cond_values, dtype=float) if cond_values else np.asarray(profile, dtype=float)
                profile = _normalize_profile(profile, normalization, float(np.nanmedian(arr)), float(np.nanpercentile(arr, 90) - np.nanpercentile(arr, 10)))
            else:
                profile = _normalize_profile(profile, normalization, baseline, scale)
            condition_profiles[condition].append(profile)
            avg_score = float(np.nanmean(np.asarray(profile, dtype=float))) if profile else 0.0
            series.append({"id": f"sample::{sample_label}", "label": sample_label, "kind": "sample", "condition": condition, "profile": profile, "avg_score": round(avg_score, 6), "sites": len(centers), "bed_source": _bed_source_label(path)})
        for condition, profiles in condition_profiles.items():
            mean_profile = _profile_mean(profiles)
            avg_score = float(np.nanmean(np.asarray(mean_profile, dtype=float))) if mean_profile else 0.0
            series.append({"id": f"condition::{condition}", "label": f"{condition} mean", "kind": "condition", "condition": condition, "profile": mean_profile, "avg_score": round(avg_score, 6), "sites": len(centers), "bed_source": _bed_source_label(path)})
        motifs.append({"prefix": prefix, "name": name, "motif_id": "", "score": len(centers), "sites": len(centers), "site_set": _bed_source_label(path), "series": series})
    return {"schema": "fp-tools.aggregate.batch.v2", "x": x, "motifs": motifs, "conditions": list(dict.fromkeys(conditions)), "colors": colors, "logos": {}, "groups_defined": bool(groups_defined), "normalization": normalization, "site_set": _summarize_site_sets([str(m.get("site_set") or "") for m in motifs]), "x_label": "Distance from motif center (bp)", "y_label": "Corrected cut-site signal (a.u.)" if normalization == "none" else "Normalized corrected cut-site signal (a.u.)"}


def _compressed_json_b64(payload: dict) -> str:
    text = json.dumps(payload, separators=(",", ":"), ensure_ascii=False)
    return base64.b64encode(gzip.compress(text.encode("utf-8"), compresslevel=9)).decode("ascii")


def _decode_payload_b64(payload_b64: str) -> dict:
    return json.loads(gzip.decompress(base64.b64decode(payload_b64)).decode("utf-8"))


def read_embedded_payload(path: str | Path) -> dict:
    text = Path(path).read_text(encoding="utf-8")
    match = re.search(r'reportPayloadB64="([^"]+)"', text)
    if not match:
        match = re.search(r"const\s+reportPayloadB64\s*=\s*['\"]([^'\"]+)['\"]", text)
    if not match:
        raise ValueError(f"Could not find reportPayloadB64 in {path}")
    return _decode_payload_b64(match.group(1))


def _series_from_diff_payload(payload: dict, source_label: str) -> dict:
    aggregate = payload.get("aggregate") or {}
    motifs = []
    conditions_seen = []
    colors = dict(payload.get("colors") or {})
    source_logos = dict(payload.get("logos") or {})
    logos = {}
    site_set = str(aggregate.get("site_set") or payload.get("site_set") or "motif-site set")
    bed_source = "bound.bed" if site_set == "bound" else ("all.bed" if site_set == "all" else site_set)
    for motif in aggregate.get("motifs") or []:
        series = []
        prefix = str(motif.get("prefix") or motif.get("name") or "motif")
        if prefix in source_logos:
            logos[prefix] = source_logos[prefix]
        for cond in motif.get("conditions") or []:
            cond_name = str(cond.get("name") or "condition")
            if cond_name not in conditions_seen:
                conditions_seen.append(cond_name)
            if cond_name not in colors:
                colors[cond_name] = colors.get(f"{cond_name}_up") or DEFAULT_COLORS[len(colors) % len(DEFAULT_COLORS)]
            for sample in cond.get("samples") or []:
                label = str(sample.get("name") or f"{source_label} {cond_name}")
                profile = sample.get("profile") or []
                avg_score = float(np.nanmean(np.asarray(profile, dtype=float))) if profile else 0.0
                series.append({"id": f"{source_label}::sample::{label}", "label": label, "kind": "sample", "condition": cond_name, "profile": profile, "avg_score": round(avg_score, 6), "fp_score": float(sample.get("fp_score") or avg_score or 0.0), "sites": int(sample.get("sites") or sample.get("n_sites") or 0), "bed_source": bed_source})
            if cond.get("profile") is not None:
                profile = cond.get("profile") or []
                avg_score = float(np.nanmean(np.asarray(profile, dtype=float))) if profile else 0.0
                series.append({"id": f"{source_label}::condition::{cond_name}", "label": f"{source_label} {cond_name} mean" if source_label else f"{cond_name} mean", "kind": "condition", "condition": cond_name, "profile": profile, "avg_score": round(avg_score, 6), "fp_score": float(cond.get("fp_score") or avg_score or 0.0), "sites": int(cond.get("sites") or cond.get("n_sites") or 0), "bed_source": bed_source})
        motifs.append({"prefix": prefix, "name": str(motif.get("name") or motif.get("prefix") or "motif"), "motif_id": str(motif.get("motif_id") or ""), "score": abs(float(motif.get("change") or 0.0)), "sites": int(motif.get("n_sites") or motif.get("sites") or 0), "site_set": bed_source, "series": series})
    return {"schema": "fp-tools.aggregate.batch.v2", "x": aggregate.get("x") or [], "motifs": motifs, "conditions": conditions_seen, "colors": {cond: colors.get(cond) or colors.get(f"{cond}_up") or DEFAULT_COLORS[idx % len(DEFAULT_COLORS)] for idx, cond in enumerate(conditions_seen)}, "logos": logos, "groups_defined": True, "normalization": aggregate.get("normalization") or payload.get("normalization") or "none", "site_set": bed_source, "x_label": aggregate.get("x_label") or "Distance from motif center (bp)", "y_label": aggregate.get("y_label") or "Corrected cut-site signal (a.u.)"}


def _ensure_batch_payload(payload: dict, source_label: str = "") -> dict:
    if payload.get("schema") == "fp-tools.aggregate.batch.v2":
        return payload
    if "aggregate" in payload:
        return _series_from_diff_payload(payload, source_label)
    if "samples" in payload:
        motifs_by_prefix: dict[str, dict] = {}
        conditions = []
        for sample in payload.get("samples") or []:
            condition = str(sample.get("condition") or sample.get("sample") or "sample")
            if condition not in conditions:
                conditions.append(condition)
            for motif in sample.get("motifs") or []:
                prefix = str(motif.get("prefix") or motif.get("name") or "motif")
                entry = motifs_by_prefix.setdefault(prefix, {"prefix": prefix, "name": str(motif.get("name") or prefix), "motif_id": str(motif.get("motif_id") or ""), "score": float(motif.get("score") or 0.0), "sites": int(motif.get("sites") or 0), "series": []})
                profile = motif.get("profile") or []
                avg_score = float(np.nanmean(np.asarray(profile, dtype=float))) if profile else 0.0
                bed_source = str(motif.get("bed_source") or motif.get("site_set") or "motif-site set")
                entry["series"].append({"id": f"sample::{sample.get('sample')}", "label": str(sample.get("label") or sample.get("sample") or condition), "kind": "sample", "condition": condition, "profile": profile, "avg_score": round(avg_score, 6), "sites": int(motif.get("sites") or 0), "bed_source": bed_source})
                entry["site_set"] = _summarize_site_sets([str(entry.get("site_set") or ""), bed_source])
        return {"schema": "fp-tools.aggregate.batch.v2", "x": payload.get("x") or [], "motifs": list(motifs_by_prefix.values()), "conditions": conditions, "colors": _condition_colors(conditions), "groups_defined": bool(payload.get("groups_defined")), "normalization": payload.get("normalization") or "none", "site_set": payload.get("site_set") or "motif-site set", "x_label": payload.get("x_label") or "Distance from motif center (bp)", "y_label": payload.get("y_label") or "Corrected cut-site signal (a.u.)"}
    raise ValueError("Unsupported aggregate HTML payload schema")


def merge_payloads(payloads: list[dict]) -> dict:
    if not payloads:
        raise ValueError("No aggregate payloads were provided")
    normalized = [_ensure_batch_payload(payload, f"report{idx + 1}") for idx, payload in enumerate(payloads)]
    merged = {"schema": "fp-tools.aggregate.batch.v2", "x": normalized[0].get("x") or [], "motifs": [], "conditions": [], "colors": {}, "logos": {}, "groups_defined": any(bool(p.get("groups_defined")) for p in normalized), "normalization": ", ".join(sorted({str(p.get("normalization") or "none") for p in normalized})), "site_set": _summarize_site_sets([str(p.get("site_set") or "") for p in normalized]), "x_label": normalized[0].get("x_label") or "Distance from motif center (bp)", "y_label": normalized[0].get("y_label") or "Corrected cut-site signal (a.u.)"}
    motifs_by_prefix: dict[str, dict] = {}
    for payload in normalized:
        for cond in payload.get("conditions") or []:
            if cond not in merged["conditions"]:
                merged["conditions"].append(cond)
        merged["colors"].update(payload.get("colors") or {})
        merged["logos"].update(payload.get("logos") or {})
        if not merged["x"] and payload.get("x"):
            merged["x"] = payload["x"]
        for motif in payload.get("motifs") or []:
            prefix = str(motif.get("prefix") or motif.get("name") or "motif")
            entry = motifs_by_prefix.setdefault(prefix, {"prefix": prefix, "name": str(motif.get("name") or prefix), "motif_id": str(motif.get("motif_id") or ""), "score": float(motif.get("score") or 0.0), "sites": int(motif.get("sites") or 0), "site_set": str(motif.get("site_set") or payload.get("site_set") or ""), "series": []})
            entry["score"] = max(float(entry.get("score") or 0.0), float(motif.get("score") or 0.0))
            entry["sites"] = max(int(entry.get("sites") or 0), int(motif.get("sites") or 0))
            entry["site_set"] = _summarize_site_sets([str(entry.get("site_set") or ""), str(motif.get("site_set") or payload.get("site_set") or "")])
            entry["series"].extend(motif.get("series") or [])
    if not merged["colors"]:
        merged["colors"] = _condition_colors(merged["conditions"])
    else:
        for idx, cond in enumerate(merged["conditions"]):
            merged["colors"].setdefault(cond, DEFAULT_COLORS[idx % len(DEFAULT_COLORS)])
    merged["motifs"] = sorted(motifs_by_prefix.values(), key=lambda m: (str(m.get("name") or ""), str(m.get("prefix") or "")))
    return merged


def write_html(payload: dict, output: str | Path, title: str, default_layout: str = "2x2", show_summary: bool = True) -> None:
    payload = _ensure_batch_payload(payload)
    payload["default_layout"] = default_layout
    escaped_title = html.escape(title)
    payload_b64 = _compressed_json_b64(payload)
    summary_css = "" if show_summary else ".main-layout{grid-template-columns:minmax(0,1fr)}.waterfall-card{display:none}"
    document = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>{escaped_title}</title><style>
:root{{--ink:#152133;--muted:#596579;--line:#d9e2ec;--grid:#e8eef5;--bg:#eef3f8;--accent:#173b73}}*{{box-sizing:border-box}}body{{margin:0;font-family:Arial,Helvetica,sans-serif;background:var(--bg);color:var(--ink);font-weight:700}}.wrap{{max-width:min(1880px,calc(100vw - 16px));margin:5px auto;padding:0 5px}}.panel,.plot-card,.waterfall-card,.card{{background:#fff;border:1px solid var(--line);border-radius:7px}}.panel{{box-shadow:0 14px 34px rgba(21,33,51,.10);overflow:hidden}}.head{{padding:7px 14px 5px;border-bottom:1px solid var(--line);background:linear-gradient(180deg,#fff 0%,#f7fafc 100%)}}h1{{margin:0;font-size:19px;line-height:1.1;font-weight:900}}.sub{{margin:2px 0 0;color:var(--muted);font-size:11px}}.top-row{{display:grid;grid-template-columns:280px minmax(540px,1fr) 210px;gap:7px;padding:6px 8px;border-bottom:1px solid var(--line);background:#fbfdff}}.card{{padding:5px;min-height:48px}}.section-title{{font-size:10px;line-height:1.05;text-transform:uppercase;letter-spacing:.08em;color:#728197;margin:0 0 4px;font-weight:900}}.controls,.summary-controls{{display:flex;flex-wrap:wrap;align-items:center;gap:5px}}label{{font-size:10px;color:#52606d;text-transform:uppercase;letter-spacing:.06em;font-weight:900}}select,input{{border:1px solid #cbd5e1;border-radius:6px;background:white;color:var(--ink);font-family:Arial,Helvetica,sans-serif;font-size:11px;font-weight:800;padding:4px 6px}}input[type=checkbox]{{width:13px;height:13px;padding:0}}button,summary{{border:1px solid #b8c5d6;background:#fff;color:var(--accent);border-radius:6px;padding:4px 7px;font-family:Arial,Helvetica,sans-serif;font-size:10px;font-weight:900;cursor:pointer}}button:hover,summary:hover{{background:#f2f6fb}}.color-row,.pill{{display:flex;align-items:center;justify-content:space-between;gap:6px;font-size:11px;color:#334e68;font-weight:800;border:1px solid #e6edf5;border-radius:999px;padding:3px 6px;background:#fbfdff}}.color-row input{{width:24px;height:17px;padding:0}}.advanced{{margin:6px 8px}}.options-grid{{display:grid;grid-template-columns:280px minmax(520px,1fr) 210px;gap:8px;margin-top:8px}}.advanced summary{{display:inline-flex;list-style:none}}.advanced summary::-webkit-details-marker{{display:none}}.style-panel{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:8px;margin-top:8px}}.sample-style-row{{display:grid;grid-template-columns:minmax(0,1fr) 28px 58px 58px 64px;gap:5px;align-items:center;border:1px solid #dbe5f0;border-radius:6px;padding:5px;background:#fff}}.sample-style-name{{font-size:10px;color:#334e68;font-weight:900;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}}.sample-style-row input[type=color]{{width:28px;height:22px;padding:1px}}.sample-style-row input[type=number]{{width:58px;min-width:58px}}.sample-style-row select{{width:64px}}.main-layout{{display:grid;grid-template-columns:330px minmax(0,1fr);align-items:start;gap:7px;padding:7px;background:#f8fbff}}.waterfall-card{{padding:6px;min-width:0;max-width:330px}}.waterfall-card svg{{width:318px;max-width:100%;margin-top:5px}}.grid{{display:grid;gap:7px;align-items:start}}.grid.g1x1{{grid-template-columns:340px}}.grid.g1x2{{grid-template-columns:repeat(2,340px)}}.grid.g2x2{{grid-template-columns:repeat(2,340px)}}.grid.g2x3{{grid-template-columns:repeat(3,340px)}}.plot-card{{padding:5px;min-width:0;width:340px;cursor:pointer}}.plot-card.active{{border-color:#173b73;box-shadow:0 0 0 2px rgba(23,59,115,.10)}}.panel-tools{{display:grid;grid-template-columns:54px minmax(0,1fr);align-items:center;gap:4px;margin-bottom:4px}}.panel-label{{font-size:10px;color:#52606d;text-transform:uppercase;letter-spacing:.06em;font-weight:900;white-space:nowrap}}.panel-actions{{grid-column:1/3;display:flex;align-items:center;justify-content:space-between;gap:5px}}.panel-tf{{width:100%;min-width:0}}.sample-picker{{position:relative;display:inline-block}}.sample-picker summary{{min-width:118px;max-width:150px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}}.sample-menu{{display:grid;gap:2px;position:absolute;right:0;z-index:30;width:210px;max-height:220px;overflow:auto;border:1px solid var(--line);border-radius:6px;background:#fff;box-shadow:0 10px 24px rgba(21,33,51,.16);padding:5px}}.sample-menu label{{display:flex;align-items:center;gap:5px;text-transform:none;letter-spacing:0;color:#334e68;font-size:11px}}svg{{width:100%;height:auto;display:block;background:#fff}}.axis{{stroke:#3b4552;stroke-width:1.2}}.grid-line{{stroke:var(--grid);stroke-width:1}}.zero{{stroke:#677386;stroke-width:1.25;stroke-dasharray:4 4}}.tick{{font-family:Arial,Helvetica,sans-serif;font-size:10px;fill:var(--muted);font-weight:800}}.summary-label{{fill:#334e68;stroke:#fff;stroke-width:2.5;paint-order:stroke;stroke-linejoin:round}}.axis-label{{font-family:Arial,Helvetica,sans-serif;font-size:10px;fill:var(--ink);font-weight:900}}.plot-title{{font-family:Arial,Helvetica,sans-serif;font-size:12px;fill:var(--ink);font-weight:900}}@media(max-width:1260px){{.options-grid{{grid-template-columns:1fr 1fr}}.style-panel{{grid-template-columns:repeat(2,minmax(0,1fr))}}.grid.g1x2,.grid.g2x2,.grid.g2x3{{grid-template-columns:repeat(2,340px)}}}}@media(max-width:980px){{.main-layout{{grid-template-columns:1fr}}.waterfall-card{{max-width:330px}}.grid.g1x2,.grid.g2x2,.grid.g2x3{{grid-template-columns:1fr}}}}@media(max-width:760px){{.options-grid,.style-panel{{grid-template-columns:1fr}}.grid.g1x1,.grid.g1x2,.grid.g2x2,.grid.g2x3{{grid-template-columns:minmax(0,1fr)}}.plot-card{{width:auto}}}}
{summary_css}</style></head><body><div class="wrap"><div class="panel"><div class="head"><h1>{escaped_title}</h1><p class="sub" id="report-detail">Standalone multi-sample, multi-TF aggregate report</p></div><details class="advanced"><summary>Plot options</summary><div class="options-grid"><div><p class="section-title">Groups</p><div class="card controls" id="color-controls"></div></div><div><p class="section-title">Layout</p><div class="card controls"><select id="layout" aria-label="Grid"><option value="1x1">1x1</option><option value="1x2">1x2</option><option value="2x2">2x2</option><option value="2x3">2x3</option></select><label class="pill"><input id="show-mean" type="checkbox">Mean</label><label>Mean width <input id="mean-width" type="number" min="0.2" max="6" step="0.1" value="1.05"></label><label>Mean type <select id="mean-type"><option value="solid">Solid</option><option value="dash">Dash</option><option value="dot">Dot</option></select></label><button id="download-grid">Download grid SVG</button></div></div><div><p class="section-title">Status</p><div class="card"><p id="status-detail" class="sub">Loading report</p></div></div></div><div><p class="section-title">Sample line styles</p><div class="style-panel" id="sample-style-controls"></div></div></div></details><div class="main-layout"><div class="waterfall-card"><p class="section-title">TF site summary</p><div class="summary-controls"><label>Sort by <select id="summary-sort"></select></label><label>Rows <select id="summary-rows"><option value="20">20</option><option value="50">50</option><option value="100">100</option><option value="all">All</option></select></label></div><svg id="waterfall-chart" viewBox="0 0 320 520"></svg></div><div id="plot-grid" class="grid g2x2"></div></div></div></div><script>
const DEFAULT_COLORS={json.dumps(DEFAULT_COLORS)};const reportPayloadB64="{payload_b64}";let payload=null,slotPrefixes=[],slotSamples=[],activeSlot=0,sampleLineStyles={{}};const reportDetail=document.getElementById('report-detail'),statusDetail=document.getElementById('status-detail'),colorControls=document.getElementById('color-controls'),layoutSel=document.getElementById('layout'),plotGrid=document.getElementById('plot-grid'),waterfallChart=document.getElementById('waterfall-chart'),sampleStyleControls=document.getElementById('sample-style-controls'),showMean=document.getElementById('show-mean'),meanWidth=document.getElementById('mean-width'),meanType=document.getElementById('mean-type'),summarySort=document.getElementById('summary-sort'),summaryRows=document.getElementById('summary-rows');
function escText(v){{return String(v??'').replace(/[&<>\"']/g,c=>({{'&':'&amp;','<':'&lt;','>':'&gt;','\"':'&quot;',"'":'&#39;'}}[c]))}}
function b64ToBytes(b64){{return Uint8Array.from(atob(b64),c=>c.charCodeAt(0))}}
async function decodePayload(){{if(!('DecompressionStream'in window))throw new Error('This report needs a modern browser with gzip DecompressionStream support.');const ds=new DecompressionStream('gzip');const stream=new Blob([b64ToBytes(reportPayloadB64)]).stream().pipeThrough(ds);return JSON.parse(await new Response(stream).text())}}
function motifLabel(m){{return m?(m.motif_id?`${{m.name}} (${{m.motif_id}})`:m.name):''}}
function motifByPrefix(prefix){{return payload.motifs.find(m=>m.prefix===prefix)||payload.motifs[0]}}
function slotCount(layout){{return layout==='1x1'?1:layout==='1x2'?2:layout==='2x2'?4:6}}
function currentColors(){{const out={{...payload.colors}};(payload.conditions||[]).forEach((c,i)=>{{out[c]=out[c]||out[`${{c}}_up`]||DEFAULT_COLORS[i%DEFAULT_COLORS.length]}});document.querySelectorAll('[data-color-condition]').forEach(inp=>out[inp.dataset.colorCondition]=inp.value);return out}}
function renderColors(){{colorControls.innerHTML=(payload.conditions||[]).map(c=>`<label class="color-row"><span>${{escText(c)}}</span><input type="color" data-color-condition="${{escText(c)}}" value="${{currentColors()[c]||'#64748b'}}"></label>`).join('');colorControls.querySelectorAll('input').forEach(inp=>inp.addEventListener('input',renderAll))}}
function allSamples(){{const seen=new Map();(payload.motifs||[]).forEach(m=>(m.series||[]).filter(s=>s.kind==='sample').forEach(s=>{{if(!seen.has(s.label))seen.set(s.label,{{label:s.label,condition:s.condition}})}}));return [...seen.values()]}}
function allSampleLabels(){{return allSamples().map(s=>s.label)}}
function ensureSlots(refill=false){{const n=slotCount(layoutSel.value),labels=allSampleLabels();if(refill||!slotPrefixes.length)slotPrefixes=payload.motifs.slice(0,n).map(m=>m.prefix);while(slotPrefixes.length<n)slotPrefixes.push(payload.motifs[slotPrefixes.length%payload.motifs.length].prefix);slotPrefixes=slotPrefixes.slice(0,n);while(slotSamples.length<n)slotSamples.push(new Set(labels));slotSamples=slotSamples.slice(0,n).map(set=>{{const clean=new Set([...set].filter(label=>labels.includes(label)));return clean.size?clean:new Set(labels)}});activeSlot=Math.min(activeSlot,n-1)}}
function lineDash(type){{return type==='dash'?'6 4':(type==='dot'?'1.2 3':'')}}
function dashAttr(type){{const dash=lineDash(type);return dash?` stroke-dasharray="${{dash}}"`:''}}
function lineWidthValue(v,fallback){{const n=Number(v);return Number.isFinite(n)&&n>0?Math.min(8,Math.max(.1,n)):fallback}}
function alphaValue(v,fallback){{const n=Number(v);return Number.isFinite(n)?Math.min(1,Math.max(.05,n)):fallback}}
function sampleStyle(label,defaults={{}}){{const stored=sampleLineStyles[label]||{{}};return{{color:stored.color||defaults.color||'#2563eb',alpha:stored.alpha??0.9,width:stored.width||2,type:stored.type||'solid'}}}}
function renderSampleStyles(){{const colors=currentColors();sampleStyleControls.innerHTML=allSamples().map(s=>{{const style=sampleStyle(s.label,{{color:colors[s.condition]}});return `<label class="sample-style-row" title="Adjust ${{escText(s.label)}}"><span class="sample-style-name">${{escText(s.label)}}</span><input data-sample-color="${{escText(s.label)}}" type="color" aria-label="Color for ${{escText(s.label)}}" value="${{style.color}}"><input data-sample-alpha="${{escText(s.label)}}" type="number" min="0.05" max="1" step="0.05" value="${{style.alpha}}"><input data-sample-width="${{escText(s.label)}}" type="number" min="0.2" max="5" step="0.1" value="${{style.width}}"><select data-sample-type="${{escText(s.label)}}"><option value="solid"${{style.type==='solid'?' selected':''}}>Solid</option><option value="dash"${{style.type==='dash'?' selected':''}}>Dash</option><option value="dot"${{style.type==='dot'?' selected':''}}>Dot</option></select></label>`}}).join('');sampleStyleControls.querySelectorAll('[data-sample-color]').forEach(el=>el.addEventListener('input',()=>{{sampleLineStyles[el.dataset.sampleColor]={{...(sampleLineStyles[el.dataset.sampleColor]||{{}}),color:el.value}};renderAll(false)}}));sampleStyleControls.querySelectorAll('[data-sample-alpha]').forEach(el=>el.addEventListener('input',()=>{{sampleLineStyles[el.dataset.sampleAlpha]={{...(sampleLineStyles[el.dataset.sampleAlpha]||{{}}),alpha:alphaValue(el.value,.9)}};renderAll(false)}}));sampleStyleControls.querySelectorAll('[data-sample-width]').forEach(el=>el.addEventListener('input',()=>{{sampleLineStyles[el.dataset.sampleWidth]={{...(sampleLineStyles[el.dataset.sampleWidth]||{{}}),width:lineWidthValue(el.value,2)}};renderAll(false)}}));sampleStyleControls.querySelectorAll('[data-sample-type]').forEach(el=>el.addEventListener('change',()=>{{sampleLineStyles[el.dataset.sampleType]={{...(sampleLineStyles[el.dataset.sampleType]||{{}}),type:el.value}};renderAll(false)}}))}}
function samplesForMotif(motif,idx=null){{const allowed=idx===null?new Set(allSampleLabels()):slotSamples[idx];return (motif.series||[]).filter(s=>s.kind==='sample'&&allowed.has(s.label))}}
function conditionMeans(motif,samples){{const conditionSeries=(motif.series||[]).filter(s=>s.kind==='condition'),by={{}};samples.forEach(s=>{{(by[s.condition]||(by[s.condition]=[])).push(s)}});return Object.entries(by).map(([condition,rows])=>{{const len=Math.max(...rows.map(r=>r.profile.length),0),profile=[];for(let i=0;i<len;i++)profile.push(rows.reduce((acc,r)=>acc+(Number(r.profile[i])||0),0)/rows.length);const avg=profile.length?profile.reduce((a,b)=>a+b,0)/profile.length:0,condMeta=conditionSeries.find(s=>s.condition===condition)||{{}},sites=Number(condMeta.sites)||rows.reduce((acc,r)=>acc+(Number(r.sites)||0),0),sampleFp=rows.map(r=>Number(r.fp_score)).filter(Number.isFinite),fpScore=Number.isFinite(Number(condMeta.fp_score))?Number(condMeta.fp_score):(sampleFp.length?sampleFp.reduce((a,b)=>a+b,0)/sampleFp.length:avg);return{{id:`condition::${{condition}}`,label:`${{condition}} mean`,kind:'condition',condition,profile,avg_score:avg,fp_score:fpScore,sites,bed_source:condMeta.bed_source||bedSummaryFromRows(rows)}}}})}}
function bedSummaryFromRows(rows){{const vals=[...new Set(rows.map(r=>r.bed_source||r.site_set).filter(Boolean))];return vals.length===1?vals[0]:(vals.length>1?'mixed beds':(payload.site_set||'motif-site set'))}}
function subplotSubtitle(motif,samples,means){{const bed=bedSummaryFromRows([...means,...samples])||motif.site_set||payload.site_set||'motif-site set';return `${{samples.length}} sample${{samples.length===1?'':'s'}} - ${{motif.sites||0}} union sites - ${{bed}}`}}
function niceTicks(min,max,n){{const out=[];for(let i=0;i<n;i++)out.push(min+(max-min)*(i/Math.max(1,n-1)));return out}}
function fmt(v){{const a=Math.abs(v);if(!Number.isFinite(v))return'';if(a===0)return'0';if(a>=1)return v.toFixed(1).replace('-0.0','0.0');if(a>=.01)return v.toFixed(2).replace('-0.00','0.00');if(a>=.001)return v.toFixed(3).replace('-0.000','0.000');return v.toExponential(1).replace('-0.0e+0','0')}}
function pathD(profile,x,sx,sy){{return profile.map((y,i)=>`${{i?'L':'M'}}${{sx(x[i]).toFixed(2)}},${{sy(y).toFixed(2)}}`).join(' ')}}
function updateStatus(){{statusDetail.textContent=`${{slotPrefixes.length}} panels - ${{payload.motifs.length}} TFs available - ${{allSamples().length}} samples`}}
function samplePickerHtml(idx){{const selected=slotSamples[idx]||new Set(allSampleLabels()),summary=`Samples: ${{selected.size}}`;return `<details class="sample-picker" data-picker="${{idx}}"><summary>${{summary}}</summary><div class="sample-menu">${{allSamples().map(s=>`<label><input type="checkbox" data-slot-sample="${{idx}}" data-sample="${{escText(s.label)}}" ${{selected.has(s.label)?'checked':''}}> ${{escText(s.label)}} <small>${{escText(s.condition)}}</small></label>`).join('')}}</div></details>`}}
function panelTfHtml(idx,prefix){{return `<span class="panel-label">Panel ${{idx+1}}</span><select class="panel-tf" data-panel-tf="${{idx}}">${{payload.motifs.map(m=>`<option value="${{escText(m.prefix)}}" ${{m.prefix===prefix?'selected':''}}>${{escText(motifLabel(m))}}</option>`).join('')}}</select>`}}
function renderSummaryControls(){{summarySort.innerHTML=`<option value="__union__">Union sites</option>`+(payload.conditions||[]).map(c=>`<option value="${{escText(c)}}">${{escText(c)}} sites</option>`).join('');if(!summarySort.value)summarySort.value='__union__'}}
function conditionSitesForMotif(motif){{const out={{}};(payload.conditions||[]).forEach(c=>out[c]=0);(motif.series||[]).filter(s=>s.kind==='condition').forEach(s=>{{out[s.condition]=Number(s.sites)||0}});(motif.series||[]).filter(s=>s.kind==='sample').forEach(s=>{{if(!out[s.condition])out[s.condition]=0;if(!(motif.series||[]).some(r=>r.kind==='condition'&&r.condition===s.condition&&Number(r.sites)))out[s.condition]+=Number(s.sites)||0}});return out}}
function drawSummary(){{
const conditions=payload.conditions||[],sortCond=summarySort.value||'__union__',limitValue=summaryRows.value||'20';
let rows=(payload.motifs||[]).map(m=>{{const sites=conditionSitesForMotif(m),unionSites=Number(m.sites)||conditions.reduce((acc,c)=>Math.max(acc,Number(sites[c])||0),0),sortSites=sortCond==='__union__'?unionSites:(Number(sites[sortCond])||0);return{{motif:m,sites,unionSites,sortSites,total:conditions.reduce((acc,c)=>acc+(Number(sites[c])||0),0)}}}}).sort((a,b)=>(b.sortSites-a.sortSites)||(b.unionSites-a.unionSites)||(b.total-a.total)||motifLabel(a.motif).localeCompare(motifLabel(b.motif)));
const totalRows=rows.length;if(limitValue!=='all')rows=rows.slice(0,Number(limitValue)||20);
const width=330,rowH=13,gap=3,margin={{top:62,bottom:26}},height=Math.max(520,margin.top+rows.length*(rowH+gap)+margin.bottom),tfX=6,colors=currentColors(),plotX=82,plotW=76,colGap=4,shownConds=conditions.slice(0,Math.max(1,Math.min(2,conditions.length))),siteColumns=[{{key:'__union__',label:'Union sites',color:'#475569'}}].concat(shownConds.map((cond,i)=>({{key:cond,label:`${{cond}} sites`,color:colors[cond]||DEFAULT_COLORS[i%DEFAULT_COLORS.length]}}))),maxSites=Math.max(...rows.flatMap(r=>siteColumns.map(col=>col.key==='__union__'?r.unionSites:(Number(r.sites[col.key])||0))),1),sortLabel=sortCond==='__union__'?'Union sites':`${{sortCond}} sites`;
waterfallChart.setAttribute('viewBox',`0 0 ${{width}} ${{height}}`);
let parts=[`<rect width="${{width}}" height="${{height}}" fill="#fff"/><text x="${{width/2}}" y="17" class="plot-title" text-anchor="middle">TF site summary</text><text x="${{width/2}}" y="32" class="tick" text-anchor="middle">Top ${{rows.length}} of ${{totalRows}} motifs by ${{escText(sortLabel)}}</text><text x="${{tfX}}" y="52" class="tick">TF</text>`];
siteColumns.forEach((col,i)=>{{const x=plotX+i*(plotW+colGap);parts.push(`<text x="${{x}}" y="52" class="tick">${{escText(col.label)}}</text><line x1="${{x}}" y1="56" x2="${{x+plotW}}" y2="56" class="grid-line"/>`)}});
rows.forEach((r,i)=>{{const y=margin.top+i*(rowH+gap),tfLabel=motifLabel(r.motif);parts.push(`<text x="${{tfX}}" y="${{y+rowH-2}}" class="tick summary-label">${{escText(tfLabel).slice(0,12)}}</text>`);siteColumns.forEach((col,j)=>{{const x=plotX+j*(plotW+colGap),sites=col.key==='__union__'?r.unionSites:(Number(r.sites[col.key])||0),barW=sites/maxSites*plotW,siteLabel=sites>=1000?`${{Math.round(sites/100)/10}}k`:String(sites);parts.push(`<rect x="${{x}}" y="${{y}}" width="${{barW}}" height="${{rowH}}" fill="${{col.color}}" fill-opacity="${{col.key==='__union__'?'0.50':'0.72'}}"><title>${{escText(tfLabel)}} ${{escText(col.label)}}=${{sites}}</title></rect><text x="${{x+plotW-2}}" y="${{y+rowH-2}}" class="tick summary-label" text-anchor="end">${{siteLabel}}</text>`)}})}});
parts.push(`<text x="${{plotX}}" y="${{height-8}}" class="tick">site counts use the plotted site set</text>`);
waterfallChart.innerHTML=parts.join('')
}}
function drawAggregate(motif,idx){{const samples=samplesForMotif(motif,idx),means=conditionMeans(motif,samples),x=payload.x,width=340,height=340,margin={{top:36,right:16,bottom:40,left:52}},innerW=width-margin.left-margin.right,innerH=height-margin.top-margin.bottom,colors=currentColors(),sampleSeries=[...samples].sort((a,b)=>(Number(b.fp_score??b.avg_score??0)-Number(a.fp_score??a.avg_score??0))),series=showMean.checked?sampleSeries.concat(means):sampleSeries;const allY=[...series.flatMap(s=>s.profile),...means.flatMap(s=>s.profile)].filter(Number.isFinite);let rawMin=Math.min(...allY,0),rawMax=Math.max(...allY,1e-9),pad=Math.max((rawMax-rawMin||1)*.06,1e-6),step=niceStep((rawMax-rawMin+2*pad)/4),ymin=Math.floor((rawMin-pad)/step)*step,ymax=Math.ceil((rawMax+pad)/step)*step;const sx=v=>margin.left+((v-x[0])/(x[x.length-1]-x[0]||1))*innerW,sy=v=>margin.top+innerH-((v-ymin)/(ymax-ymin||1))*innerH,xTicks=[x[0],0,x[x.length-1]],yTicks=niceTicks(ymin,ymax,4);let parts=[`<svg class="aggregate-panel" data-panel="${{idx}}" viewBox="0 0 ${{width}} ${{height}}"><rect width="${{width}}" height="${{height}}" fill="#fff"/><text x="${{width/2}}" y="16" text-anchor="middle" class="plot-title">${{escText(motifLabel(motif))}}</text><text x="${{width/2}}" y="30" text-anchor="middle" class="tick">${{escText(subplotSubtitle(motif,samples,means))}}</text>`];yTicks.forEach(v=>parts.push(`<text x="${{margin.left-7}}" y="${{sy(v)+3}}" class="tick" text-anchor="end">${{fmt(v)}}</text>`));xTicks.forEach(v=>parts.push(`<line x1="${{sx(v)}}" y1="${{margin.top}}" x2="${{sx(v)}}" y2="${{margin.top+innerH}}" class="grid-line"/><text x="${{sx(v)}}" y="${{margin.top+innerH+17}}" class="tick" text-anchor="middle">${{v}}</text>`));parts.push(`<line x1="${{sx(0)}}" y1="${{margin.top}}" x2="${{sx(0)}}" y2="${{margin.top+innerH}}" class="zero"/><line x1="${{margin.left}}" y1="${{margin.top+innerH}}" x2="${{margin.left+innerW}}" y2="${{margin.top+innerH}}" class="axis"/><line x1="${{margin.left}}" y1="${{margin.top}}" x2="${{margin.left}}" y2="${{margin.top+innerH}}" class="axis"/>`);sampleSeries.forEach((s,i)=>{{const style=sampleStyle(s.label,{{color:colors[s.condition]||DEFAULT_COLORS[i%DEFAULT_COLORS.length]}}),dash=dashAttr(style.type);parts.push(`<path d="${{pathD(s.profile,x,sx,sy)}}" fill="none" stroke="${{style.color}}" stroke-width="${{lineWidthValue(style.width,2)}}"${{dash}} stroke-opacity="${{alphaValue(style.alpha,.9)}}"><title>${{escText(s.label)}} - ${{escText(s.condition)}}</title></path>`)}});if(showMean.checked)means.forEach((s,i)=>{{const dash=dashAttr(meanType.value),color=colors[s.condition]||DEFAULT_COLORS[i%DEFAULT_COLORS.length];parts.push(`<path d="${{pathD(s.profile,x,sx,sy)}}" fill="none" stroke="${{color}}" stroke-width="${{lineWidthValue(meanWidth.value,1.05)}}"${{dash}} stroke-opacity="0.95" stroke-linecap="round"/><text x="${{margin.left+8}}" y="${{margin.top+11+i*11}}" font-family="Arial,Helvetica,sans-serif" font-size="8.5" font-weight="900" fill="${{color}}">${{escText(s.condition)}} mean</text>`)}});parts.push(`<text x="${{margin.left+innerW/2}}" y="${{height-10}}" class="axis-label" text-anchor="middle">${{escText(payload.x_label)}}</text><text x="14" y="${{margin.top+innerH/2}}" class="axis-label" text-anchor="middle" transform="rotate(-90 14 ${{margin.top+innerH/2}})">${{escText(payload.y_label)}}</text></svg>`);return parts.join('')}}
function renderPlots(){{ensureSlots();updateStatus();plotGrid.className=`grid g${{layoutSel.value}}`;plotGrid.innerHTML=slotPrefixes.map((prefix,idx)=>{{const motif=motifByPrefix(prefix);return `<div class="plot-card${{idx===activeSlot?' active':''}}" data-card="${{idx}}"><div class="panel-tools">${{panelTfHtml(idx,prefix)}}<div class="panel-actions">${{samplePickerHtml(idx)}}<button data-download-panel="${{idx}}">Download SVG</button></div></div>${{drawAggregate(motif,idx)}}</div>`}}).join('');plotGrid.querySelectorAll('.plot-card').forEach(card=>card.addEventListener('click',ev=>{{if(ev.target.closest('button,details,input,label,summary,select'))return;activeSlot=Number(card.dataset.card);renderAll(false)}}));plotGrid.querySelectorAll('[data-panel-tf]').forEach(sel=>sel.addEventListener('change',()=>{{const idx=Number(sel.dataset.panelTf);slotPrefixes[idx]=sel.value;activeSlot=idx;renderAll(false)}}));plotGrid.querySelectorAll('[data-slot-sample]').forEach(inp=>inp.addEventListener('change',()=>{{const idx=Number(inp.dataset.slotSample),set=slotSamples[idx];if(inp.checked)set.add(inp.dataset.sample);else set.delete(inp.dataset.sample);if(!set.size)allSampleLabels().forEach(label=>set.add(label));renderAll(false)}}));plotGrid.querySelectorAll('[data-download-panel]').forEach(btn=>btn.addEventListener('click',ev=>{{ev.stopPropagation();const idx=Number(btn.dataset.downloadPanel),svg=document.querySelector(`.aggregate-panel[data-panel="${{idx}}"]`);if(svg)downloadBlob(svgBlob(svg),`plot_aggregate_batch_panel_${{idx+1}}.svg`)}}));drawSummary()}}
function renderAll(refreshStyles=true){{renderColors();if(refreshStyles)renderSampleStyles();renderPlots()}}
function svgBlob(svgNode){{const clone=svgNode.cloneNode(true);clone.setAttribute('xmlns','http://www.w3.org/2000/svg');return new Blob([new XMLSerializer().serializeToString(clone)],{{type:'image/svg+xml;charset=utf-8'}})}}
function downloadBlob(blob,filename){{const url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download=filename;document.body.appendChild(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(url),1000)}}
function downloadGrid(){{const svgs=[...document.querySelectorAll('.aggregate-panel')];if(!svgs.length)return;const w=340,h=340,cols=layoutSel.value==='2x3'?3:(layoutSel.value==='1x1'?1:2),rows=Math.ceil(svgs.length/cols);let parts=[`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${{cols*w}} ${{rows*h}}">`];svgs.forEach((svg,i)=>parts.push(`<g transform="translate(${{(i%cols)*w}},${{Math.floor(i/cols)*h}})">${{svg.innerHTML}}</g>`));parts.push('</svg>');downloadBlob(new Blob(parts,{{type:'image/svg+xml;charset=utf-8'}}),'plot_aggregate_batch_grid.svg')}}
document.getElementById('download-grid').addEventListener('click',downloadGrid);[layoutSel,showMean,meanType].forEach(el=>el.addEventListener('change',()=>{{ensureSlots(layoutSel===el);renderAll(false)}}));[summarySort,summaryRows].forEach(el=>el.addEventListener('change',()=>renderAll(false)));meanWidth.addEventListener('input',()=>renderAll(false));decodePayload().then(data=>{{payload=data;layoutSel.value=payload.default_layout||'2x2';ensureSlots(true);renderSummaryControls();reportDetail.textContent=`${{payload.conditions.length}} groups - ${{payload.motifs.length}} TFs - ${{allSamples().length}} samples`;renderAll()}}).catch(err=>{{statusDetail.textContent=`Could not open report payload: ${{err.message}}`}});
</script></body></html>"""
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    Path(output).write_text(document, encoding="utf-8")


def write_html(payload: dict, output: str | Path, title: str, default_layout: str = "2x2", show_summary: bool = True) -> None:
    payload = _ensure_batch_payload(payload)
    payload["default_layout"] = default_layout
    escaped_title = html.escape(title)
    payload_b64 = _compressed_json_b64(payload)
    summary_css = "" if show_summary else ".main-layout{grid-template-columns:minmax(0,1fr)}.waterfall-card{display:none}"
    document = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>{escaped_title}</title><style>
:root{{--ink:#152133;--muted:#596579;--line:#d9e2ec;--grid:#e8eef5;--bg:#eef3f8;--accent:#173b73}}*{{box-sizing:border-box}}body{{margin:0;font-family:Arial,Helvetica,sans-serif;background:var(--bg);color:var(--ink);font-weight:700}}.wrap{{max-width:min(1900px,calc(100vw - 12px));margin:4px auto;padding:0 4px}}.panel,.plot-card,.waterfall-card,.card,.motif-card{{background:#fff;border:1px solid var(--line);border-radius:7px}}.panel{{box-shadow:0 12px 30px rgba(21,33,51,.10);overflow:hidden}}.head{{padding:8px 14px 6px;border-bottom:1px solid var(--line);background:#fff}}h1{{margin:0;font-size:22px;line-height:1.08;font-weight:900}}.sub{{margin:2px 0 0;color:var(--muted);font-size:11px}}.options{{padding:7px 8px;border-bottom:1px solid var(--line);background:#fbfdff}}.option-grid{{display:grid;grid-template-columns:260px 420px minmax(0,1fr);gap:8px;align-items:stretch}}.card{{padding:5px 6px;min-width:0}}.top-controls .card{{display:flex;align-items:center;gap:8px;min-height:34px}}.top-controls .card:first-child{{flex:2 1 620px}}.top-controls .card:nth-child(2){{flex:0 0 390px}}.top-controls .card:nth-child(2) .control-row{{flex-wrap:nowrap}}.top-controls .card:nth-child(3){{flex:1 1 260px}}.section-title{{font-size:10px;line-height:1.05;text-transform:uppercase;letter-spacing:.08em;color:#728197;margin:0;font-weight:900;white-space:nowrap}}.controls{{display:flex;flex-wrap:wrap;align-items:center;gap:6px}}label{{font-size:10px;color:#52606d;text-transform:uppercase;letter-spacing:.06em;font-weight:900}}select,input{{border:1px solid #cbd5e1;border-radius:6px;background:white;color:var(--ink);font-family:Arial,Helvetica,sans-serif;font-size:11px;font-weight:800;padding:4px 6px}}input[type=checkbox]{{width:13px;height:13px;padding:0}}button{{border:1px solid #b8c5d6;background:#fff;color:var(--accent);border-radius:6px;padding:4px 7px;font-family:Arial,Helvetica,sans-serif;font-size:10px;font-weight:900;cursor:pointer}}button:hover{{background:#f2f6fb}}.export-stack{{display:grid;gap:5px}}.export-stack button{{text-align:left}}.color-grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(82px,1fr));gap:5px;margin-top:6px}}.color-row{{display:flex;align-items:center;justify-content:space-between;gap:5px;font-size:11px;color:#334e68;font-weight:800;border:1px solid #e6edf5;border-radius:999px;padding:3px 6px;background:#fbfdff;min-width:0}}.color-row span{{white-space:nowrap;overflow:hidden;text-overflow:ellipsis}}.color-row input{{width:24px;height:17px;padding:0}}.sample-style-panel{{display:grid;grid-template-columns:1fr;gap:5px;max-height:238px;overflow:auto}}.sample-style-row{{display:grid;grid-template-columns:17px minmax(90px,1fr) 30px 54px 54px 64px;gap:5px;align-items:center;border:1px solid #dbe5f0;border-radius:6px;padding:4px;background:#fff}}.sample-style-name{{font-size:10px;color:#334e68;font-weight:900;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}}.sample-style-row input[type=color]{{width:30px;height:21px;padding:1px}}.sample-style-row input[type=number]{{width:54px;min-width:54px}}.sample-style-row select{{width:64px}}.motif-toolbar{{display:flex;align-items:center;justify-content:space-between;gap:8px;margin-bottom:5px}}.motif-grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(165px,1fr));gap:6px}}.motif-card{{padding:5px;min-width:0}}.motif-card.active{{outline:2px solid #93c5fd}}.motif-card select{{width:100%;min-width:0;margin-bottom:4px}}.motif-logo{{height:90px;display:flex;align-items:center;justify-content:center;overflow:hidden}}.motif-logo img{{max-width:100%;max-height:88px}}.logo-empty{{font-size:11px;color:#94a3b8}}.main-layout{{display:grid;grid-template-columns:330px minmax(0,1fr);align-items:start;gap:7px;padding:7px;background:#f8fbff}}.waterfall-card{{padding:6px;min-width:0;max-width:330px}}.waterfall-card svg{{width:318px;max-width:100%;margin-top:5px}}.grid{{display:grid;gap:7px;align-items:start;position:relative}}.grid.cols1{{grid-template-columns:340px}}.grid.cols2{{grid-template-columns:repeat(2,340px)}}.grid.cols3{{grid-template-columns:repeat(3,340px)}}.grid.cols4{{grid-template-columns:repeat(4,340px)}}.plot-card{{padding:5px;min-width:0;width:340px;cursor:pointer}}.plot-card.active{{border-color:#173b73;box-shadow:0 0 0 2px rgba(23,59,115,.10)}}.plot-card.grouped{{border-color:#0f766e;box-shadow:0 0 0 2px rgba(15,118,110,.12)}}.panel-tools{{display:grid;grid-template-columns:54px minmax(0,1fr);align-items:center;gap:4px;margin-bottom:4px}}.panel-label{{font-size:10px;color:#52606d;text-transform:uppercase;letter-spacing:.06em;font-weight:900;white-space:nowrap}}.panel-actions{{grid-column:1/3;display:flex;align-items:center;justify-content:space-between;gap:5px}}.panel-tf{{width:100%;min-width:0}}.sample-picker{{position:relative;display:inline-block}}.sample-picker summary{{list-style:none;border:1px solid #b8c5d6;background:#fff;color:var(--accent);border-radius:6px;padding:4px 7px;font-size:10px;font-weight:900;min-width:118px;max-width:150px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;cursor:pointer}}.sample-picker summary::-webkit-details-marker{{display:none}}.sample-menu{{display:grid;gap:2px;position:absolute;right:0;z-index:30;width:210px;max-height:220px;overflow:auto;border:1px solid var(--line);border-radius:6px;background:#fff;box-shadow:0 10px 24px rgba(21,33,51,.16);padding:5px}}.sample-menu label{{display:flex;align-items:center;gap:5px;text-transform:none;letter-spacing:0;color:#334e68;font-size:11px}}svg{{width:100%;height:auto;display:block;background:#fff}}.axis{{stroke:#3b4552;stroke-width:1.2}}.grid-line{{stroke:var(--grid);stroke-width:1}}.zero{{stroke:#677386;stroke-width:1.25;stroke-dasharray:4 4}}.tick{{font-family:Arial,Helvetica,sans-serif;font-size:10px;fill:var(--muted);font-weight:800}}.summary-label{{fill:#334e68;stroke:#fff;stroke-width:2.5;paint-order:stroke;stroke-linejoin:round}}.axis-label{{font-family:Arial,Helvetica,sans-serif;font-size:10px;fill:var(--ink);font-weight:900}}.plot-title{{font-family:Arial,Helvetica,sans-serif;font-size:12px;fill:var(--ink);font-weight:900}}.legend-box{{position:absolute;top:6px;right:6px;z-index:5;background:rgba(255,255,255,.92);border:1px solid #d9e2ec;border-radius:6px;padding:5px;display:grid;gap:3px;max-width:220px}}.legend-row{{display:grid;grid-template-columns:32px minmax(0,1fr);gap:5px;align-items:center;font-size:10px;color:#334e68;font-weight:900}}.legend-row span{{overflow:hidden;text-overflow:ellipsis;white-space:nowrap}}@media(max-width:1500px){{.option-grid{{grid-template-columns:250px minmax(360px,1fr)}}.selected-card{{grid-column:1/3}}.grid.cols3,.grid.cols4{{grid-template-columns:repeat(2,340px)}}}}@media(max-width:980px){{.option-grid,.main-layout{{grid-template-columns:1fr}}.selected-card{{grid-column:auto}}.waterfall-card{{max-width:none}}.grid.cols3,.grid.cols4{{grid-template-columns:repeat(2,340px)}}}}@media(max-width:760px){{.motif-grid{{grid-template-columns:1fr}}.plot-card{{width:100%}}.grid.cols1,.grid.cols2,.grid.cols3,.grid.cols4{{grid-template-columns:minmax(0,1fr)}}}}{summary_css}
</style></head><body><div class="wrap"><section class="panel"><header class="head"><h1>{escaped_title}</h1><p class="sub" id="report-detail"></p></header><div class="options"><div class="option-grid"><div class="card"><p class="section-title">Export editable SVG</p><div class="export-stack"><button id="download-logo">Download motif logo panel</button><button id="download-summary">Download bar plot panel</button><button id="download-grid">Download motif aggregate panel</button></div><p class="section-title" style="margin-top:8px">Groups</p><div class="color-grid" id="color-controls"></div></div><div class="card"><p class="section-title">Sample line styles</p><div id="sample-style-controls" class="sample-style-panel"></div></div><div class="card selected-card"><div class="motif-toolbar"><p class="section-title">Selected motifs</p><div class="controls"><label>Plots <select id="plot-count"></select></label><label>Columns <select id="plot-cols"><option value="1">1</option><option value="2" selected>2</option><option value="3">3</option><option value="4">4</option></select></label><button id="group-autoscale">Group autoscale</button><button id="reset-autoscale">Reset autoscale</button></div></div><div id="motif-grid" class="motif-grid"></div></div></div></div><div class="main-layout"><aside class="waterfall-card"><p class="section-title">TF site summary</p><div class="controls"><label>Sort <select id="summary-sort"></select></label><label>Rows <select id="summary-rows"><option value="20">20</option><option value="40">40</option><option value="80">80</option><option value="all">All</option></select></label></div><svg id="waterfall-chart" viewBox="0 0 330 520"></svg></aside><main id="plot-grid" class="grid cols2"></main></div></section></div><script>
const reportPayloadB64="{payload_b64}";
const DEFAULT_COLORS={json.dumps(DEFAULT_COLORS)};
let payload=null,slotPrefixes=[],slotSamples=[],sampleLineStyles={{}},activeSlot=0,groupAutoscalePanels=new Set(),groupAutoscaleDomain=null;
const plotGrid=document.getElementById('plot-grid'),motifGrid=document.getElementById('motif-grid'),plotCountSel=document.getElementById('plot-count'),plotColsSel=document.getElementById('plot-cols'),sampleStyleControls=document.getElementById('sample-style-controls'),colorControls=document.getElementById('color-controls'),summarySort=document.getElementById('summary-sort'),summaryRows=document.getElementById('summary-rows'),waterfallChart=document.getElementById('waterfall-chart'),reportDetail=document.getElementById('report-detail');
function escText(v){{return String(v??'').replace(/[&<>"']/g,c=>({{'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}}[c]))}}
async function decodePayload(){{const bin=Uint8Array.from(atob(reportPayloadB64),c=>c.charCodeAt(0));if('DecompressionStream' in window){{const stream=new Blob([bin]).stream().pipeThrough(new DecompressionStream('gzip'));return JSON.parse(await new Response(stream).text())}}throw new Error('This browser cannot decompress the embedded report payload.')}}
function currentColors(){{const out={{...(payload.colors||{{}})}};colorControls.querySelectorAll('input').forEach(inp=>out[inp.dataset.colorCondition]=inp.value);return out}}
function motifLabel(m){{return m?(m.motif_id?`${{m.name}} (${{m.motif_id}})`:m.name||m.prefix):''}}
function motifByPrefix(prefix){{return (payload.motifs||[]).find(m=>m.prefix===prefix)||payload.motifs[0]}}
function allSamples(){{const seen=new Map();(payload.motifs||[]).forEach(m=>(m.series||[]).filter(s=>s.kind==='sample').forEach(s=>{{if(!seen.has(s.label))seen.set(s.label,{{label:s.label,condition:s.condition}})}}));return [...seen.values()]}}
function allSampleLabels(){{return allSamples().map(s=>s.label)}}
function initialPlotCount(){{return payload.default_layout==='1x1'?1:(payload.default_layout==='1x2'?2:(payload.default_layout==='2x3'?6:4))}}
function ensureSlots(refill=false){{const n=Number(plotCountSel.value)||initialPlotCount(),labels=allSampleLabels();if(refill||!slotPrefixes.length)slotPrefixes=(payload.motifs||[]).slice(0,n).map(m=>m.prefix);while(slotPrefixes.length<n)slotPrefixes.push((payload.motifs||[])[slotPrefixes.length%Math.max(1,payload.motifs.length)]?.prefix);slotPrefixes=slotPrefixes.slice(0,n);while(slotSamples.length<n)slotSamples.push(new Set(labels));slotSamples=slotSamples.slice(0,n).map(set=>{{const clean=new Set([...set].filter(label=>labels.includes(label)));return clean.size?clean:new Set(labels)}});activeSlot=Math.max(0,Math.min(activeSlot,n-1));groupAutoscalePanels=new Set([...groupAutoscalePanels].filter(i=>i<n))}}
function lineDash(type){{return type==='dash'?'6 4':(type==='dot'?'1.2 3':'')}}
function dashAttr(type){{const dash=lineDash(type);return dash?` stroke-dasharray="${{dash}}"`:''}}
function lineWidthValue(v,fallback){{const n=Number(v);return Number.isFinite(n)&&n>0?Math.min(8,Math.max(.1,n)):fallback}}
function alphaValue(v,fallback){{const n=Number(v);return Number.isFinite(n)?Math.min(1,Math.max(.05,n)):fallback}}
function sampleStyle(label,defaults={{}}){{const stored=sampleLineStyles[label]||{{}};return{{visible:stored.visible!==false,color:stored.color||defaults.color||'#2563eb',alpha:stored.alpha??0.9,width:stored.width||2,type:stored.type||'solid'}}}}
function renderColors(){{colorControls.innerHTML=(payload.conditions||[]).map(c=>`<label class="color-row"><span>${{escText(c)}}</span><input type="color" data-color-condition="${{escText(c)}}" value="${{currentColors()[c]||'#64748b'}}"></label>`).join('');colorControls.querySelectorAll('input').forEach(inp=>inp.addEventListener('input',()=>{{renderAll(false)}}))}}
function renderSampleStyles(){{const colors=currentColors();sampleStyleControls.innerHTML=allSamples().map(s=>{{const style=sampleStyle(s.label,{{color:colors[s.condition]}});return `<label class="sample-style-row" title="Adjust ${{escText(s.label)}}"><input data-sample-visible="${{escText(s.label)}}" type="checkbox" ${{style.visible?'checked':''}}><span class="sample-style-name">${{escText(s.label)}}</span><input data-sample-color="${{escText(s.label)}}" type="color" value="${{style.color}}"><input data-sample-alpha="${{escText(s.label)}}" type="number" min="0.05" max="1" step="0.05" value="${{style.alpha}}"><input data-sample-width="${{escText(s.label)}}" type="number" min="0.2" max="5" step="0.1" value="${{style.width}}"><select data-sample-type="${{escText(s.label)}}"><option value="solid"${{style.type==='solid'?' selected':''}}>Solid</option><option value="dash"${{style.type==='dash'?' selected':''}}>Dash</option><option value="dot"${{style.type==='dot'?' selected':''}}>Dot</option></select></label>`}}).join('');sampleStyleControls.querySelectorAll('[data-sample-visible]').forEach(el=>el.addEventListener('change',()=>{{sampleLineStyles[el.dataset.sampleVisible]={{...(sampleLineStyles[el.dataset.sampleVisible]||{{}}),visible:el.checked}};renderAll(false)}}));sampleStyleControls.querySelectorAll('[data-sample-color]').forEach(el=>el.addEventListener('input',()=>{{sampleLineStyles[el.dataset.sampleColor]={{...(sampleLineStyles[el.dataset.sampleColor]||{{}}),color:el.value}};renderAll(false)}}));sampleStyleControls.querySelectorAll('[data-sample-alpha]').forEach(el=>el.addEventListener('input',()=>{{sampleLineStyles[el.dataset.sampleAlpha]={{...(sampleLineStyles[el.dataset.sampleAlpha]||{{}}),alpha:alphaValue(el.value,.9)}};renderAll(false)}}));sampleStyleControls.querySelectorAll('[data-sample-width]').forEach(el=>el.addEventListener('input',()=>{{sampleLineStyles[el.dataset.sampleWidth]={{...(sampleLineStyles[el.dataset.sampleWidth]||{{}}),width:lineWidthValue(el.value,2)}};renderAll(false)}}));sampleStyleControls.querySelectorAll('[data-sample-type]').forEach(el=>el.addEventListener('change',()=>{{sampleLineStyles[el.dataset.sampleType]={{...(sampleLineStyles[el.dataset.sampleType]||{{}}),type:el.value}};renderAll(false)}}))}}
function samplesForMotif(motif,idx=null){{const allowed=idx===null?new Set(allSampleLabels()):slotSamples[idx];return (motif.series||[]).filter(s=>s.kind==='sample'&&allowed.has(s.label)&&sampleStyle(s.label).visible)}}
function bedSummaryFromRows(rows){{const vals=[...new Set(rows.map(r=>r.bed_source||r.site_set).filter(Boolean))];return vals.length===1?vals[0]:(vals.length>1?'mixed beds':(payload.site_set||'motif-site set'))}}
function niceTicks(min,max,n){{const out=[];for(let i=0;i<n;i++)out.push(min+(max-min)*(i/Math.max(1,n-1)));return out}}
function fmt(v){{const a=Math.abs(v);if(!Number.isFinite(v))return'';if(a===0)return'0';if(a>=1)return v.toFixed(1).replace('-0.0','0.0');if(a>=.01)return v.toFixed(2).replace('-0.00','0.00');if(a>=.001)return v.toFixed(3).replace('-0.000','0.000');return v.toExponential(1).replace('-0.0e+0','0')}}
function pathD(profile,x,sx,sy){{return profile.map((y,i)=>`${{i?'L':'M'}}${{sx(x[i]).toFixed(2)}},${{sy(y).toFixed(2)}}`).join(' ')}}
function yDomainForPanel(idx){{if(groupAutoscaleDomain&&groupAutoscalePanels.has(idx))return groupAutoscaleDomain;const motif=motifByPrefix(slotPrefixes[idx]),samples=samplesForMotif(motif,idx),allY=samples.flatMap(s=>s.profile).filter(Number.isFinite);let rawMin=Math.min(...allY,0),rawMax=Math.max(...allY,1e-9),pad=Math.max((rawMax-rawMin||1)*.06,1e-6),step=niceStep((rawMax-rawMin+2*pad)/4),ymin=Math.floor((rawMin-pad)/step)*step,ymax=Math.ceil((rawMax+pad)/step)*step;return [ymin,ymax]}}
function samplePickerHtml(idx){{const selected=slotSamples[idx]||new Set(allSampleLabels()),summary=`Samples: ${{selected.size}}`;return `<details class="sample-picker"><summary>${{summary}}</summary><div class="sample-menu">${{allSamples().map(s=>`<label><input type="checkbox" data-slot-sample="${{idx}}" data-sample="${{escText(s.label)}}" ${{selected.has(s.label)?'checked':''}}> ${{escText(s.label)}} <small>${{escText(s.condition)}}</small></label>`).join('')}}</div></details>`}}
function panelTfHtml(idx,prefix){{return `<span class="panel-label">Plot ${{idx+1}}</span><select class="panel-tf" data-panel-tf="${{idx}}">${{(payload.motifs||[]).map(m=>`<option value="${{escText(m.prefix)}}" ${{m.prefix===prefix?'selected':''}}>${{escText(motifLabel(m))}}</option>`).join('')}}</select>`}}
function logoHtml(prefix){{const logo=(payload.logos||{{}})[prefix]||{{}};return logo.png?`<img alt="Motif logo" src="${{logo.png}}">`:'<span class="logo-empty">Motif logo unavailable</span>'}}
function renderSelectedMotifs(){{motifGrid.innerHTML=slotPrefixes.map((prefix,idx)=>{{const motif=motifByPrefix(prefix);return `<div class="motif-card${{idx===activeSlot?' active':''}}" data-motif-card="${{idx}}"><select data-card-tf="${{idx}}">${{(payload.motifs||[]).map(m=>`<option value="${{escText(m.prefix)}}" ${{m.prefix===prefix?'selected':''}}>${{escText(motifLabel(m))}}</option>`).join('')}}</select><div class="motif-logo">${{logoHtml(prefix)}}</div></div>`}}).join('');motifGrid.querySelectorAll('[data-motif-card]').forEach(card=>card.addEventListener('click',ev=>{{if(ev.target.closest('select'))return;activeSlot=Number(card.dataset.motifCard);renderAll(false)}}));motifGrid.querySelectorAll('[data-card-tf]').forEach(sel=>sel.addEventListener('change',()=>{{slotPrefixes[Number(sel.dataset.cardTf)]=sel.value;activeSlot=Number(sel.dataset.cardTf);groupAutoscaleDomain=null;renderAll(false)}}))}}
function renderSummaryControls(){{summarySort.innerHTML=`<option value="__union__">Union sites</option>`+(payload.conditions||[]).map(c=>`<option value="${{escText(c)}}">${{escText(c)}} sites</option>`).join('');if(!summarySort.value)summarySort.value='__union__'}}
function conditionSitesForMotif(motif){{const out={{}};(payload.conditions||[]).forEach(c=>out[c]=0);(motif.series||[]).filter(s=>s.kind==='condition').forEach(s=>{{out[s.condition]=Number(s.sites)||0}});(motif.series||[]).filter(s=>s.kind==='sample').forEach(s=>{{if(!out[s.condition])out[s.condition]=0;if(!(motif.series||[]).some(r=>r.kind==='condition'&&r.condition===s.condition&&Number(r.sites)))out[s.condition]+=Number(s.sites)||0}});return out}}
function drawSummary(){{const conditions=payload.conditions||[],sortCond=summarySort.value||'__union__',limitValue=summaryRows.value||'20';let rows=(payload.motifs||[]).map(m=>{{const sites=conditionSitesForMotif(m),unionSites=Number(m.sites)||conditions.reduce((acc,c)=>Math.max(acc,Number(sites[c])||0),0),sortSites=sortCond==='__union__'?unionSites:(Number(sites[sortCond])||0);return{{motif:m,sites,unionSites,sortSites,total:conditions.reduce((acc,c)=>acc+(Number(sites[c])||0),0)}}}}).sort((a,b)=>(b.sortSites-a.sortSites)||(b.unionSites-a.unionSites)||(b.total-a.total)||motifLabel(a.motif).localeCompare(motifLabel(b.motif)));const totalRows=rows.length;if(limitValue!=='all')rows=rows.slice(0,Number(limitValue)||20);const width=330,rowH=13,gap=3,margin={{top:62,bottom:26}},height=Math.max(520,margin.top+rows.length*(rowH+gap)+margin.bottom),tfX=6,colors=currentColors(),plotX=86,plotW=72,colGap=4,shownConds=conditions.slice(0,Math.max(1,Math.min(2,conditions.length))),siteColumns=[{{key:'__union__',label:'Union',color:'#475569'}}].concat(shownConds.map((cond,i)=>({{key:cond,label:cond,color:colors[cond]||DEFAULT_COLORS[i%DEFAULT_COLORS.length]}}))),maxSites=Math.max(...rows.flatMap(r=>siteColumns.map(col=>col.key==='__union__'?r.unionSites:(Number(r.sites[col.key])||0))),1),sortLabel=sortCond==='__union__'?'Union sites':`${{sortCond}} sites`;waterfallChart.setAttribute('viewBox',`0 0 ${{width}} ${{height}}`);let parts=[`<rect width="${{width}}" height="${{height}}" fill="#fff"/><text x="${{width/2}}" y="17" class="plot-title" text-anchor="middle">TF site summary</text><text x="${{width/2}}" y="32" class="tick" text-anchor="middle">Top ${{rows.length}} of ${{totalRows}} motifs by ${{escText(sortLabel)}}</text><text x="${{tfX}}" y="52" class="tick">TF</text>`];siteColumns.forEach((col,i)=>{{const x=plotX+i*(plotW+colGap);parts.push(`<text x="${{x}}" y="52" class="tick">${{escText(col.label)}}</text><line x1="${{x}}" y1="56" x2="${{x+plotW}}" y2="56" class="grid-line"/>`)}});rows.forEach((r,i)=>{{const y=margin.top+i*(rowH+gap),tfLabel=motifLabel(r.motif);parts.push(`<text x="${{tfX}}" y="${{y+rowH-2}}" class="tick summary-label">${{escText(tfLabel).slice(0,13)}}</text>`);siteColumns.forEach((col,j)=>{{const x=plotX+j*(plotW+colGap),sites=col.key==='__union__'?r.unionSites:(Number(r.sites[col.key])||0),barW=sites/maxSites*plotW,siteLabel=sites>=1000?`${{Math.round(sites/100)/10}}k`:String(sites);parts.push(`<rect x="${{x}}" y="${{y}}" width="${{barW}}" height="${{rowH}}" fill="${{col.color}}" fill-opacity="${{col.key==='__union__'?'0.50':'0.72'}}"><title>${{escText(tfLabel)}} ${{escText(col.label)}}=${{sites}}</title></rect><text x="${{x+plotW-2}}" y="${{y+rowH-2}}" class="tick summary-label" text-anchor="end">${{siteLabel}}</text>`)}})}});parts.push(`<text x="${{plotX}}" y="${{height-8}}" class="tick">site counts use plotted BED</text>`);waterfallChart.innerHTML=parts.join('')}}
function drawAggregate(motif,idx){{const samples=samplesForMotif(motif,idx),x=payload.x,width=340,height=340,margin={{top:28,right:14,bottom:38,left:50}},innerW=width-margin.left-margin.right,innerH=height-margin.top-margin.bottom,colors=currentColors(),sampleSeries=[...samples].sort((a,b)=>(Number(b.fp_score??b.avg_score??0)-Number(a.fp_score??a.avg_score??0))),domain=yDomainForPanel(idx),ymin=domain[0],ymax=domain[1],sx=v=>margin.left+((v-x[0])/(x[x.length-1]-x[0]||1))*innerW,sy=v=>margin.top+innerH-((v-ymin)/(ymax-ymin||1))*innerH,xTicks=[x[0],0,x[x.length-1]],yTicks=niceTicks(ymin,ymax,4),bed=bedSummaryFromRows(sampleSeries)||motif.site_set||payload.site_set||'bed';let parts=[`<svg class="aggregate-panel" data-panel="${{idx}}" viewBox="0 0 ${{width}} ${{height}}"><style>.axis{{stroke:#3b4552;stroke-width:1.2}}.grid-line{{stroke:#e8eef5;stroke-width:1}}.zero{{stroke:#677386;stroke-width:1.25;stroke-dasharray:4 4}}.tick{{font-family:Arial,Helvetica,sans-serif;font-size:10px;fill:#596579;font-weight:800}}.axis-label{{font-family:Arial,Helvetica,sans-serif;font-size:10px;fill:#152133;font-weight:900}}.plot-title{{font-family:Arial,Helvetica,sans-serif;font-size:12px;fill:#152133;font-weight:900}}</style><rect width="${{width}}" height="${{height}}" fill="#fff"/><text x="${{width/2}}" y="16" text-anchor="middle" class="plot-title">${{escText(motifLabel(motif))}}</text>`];yTicks.forEach(v=>parts.push(`<text x="${{margin.left-7}}" y="${{sy(v)+3}}" class="tick" text-anchor="end">${{fmt(v)}}</text>`));xTicks.forEach(v=>parts.push(`<line x1="${{sx(v)}}" y1="${{margin.top}}" x2="${{sx(v)}}" y2="${{margin.top+innerH}}" class="grid-line"/><text x="${{sx(v)}}" y="${{margin.top+innerH+17}}" class="tick" text-anchor="middle">${{v}}</text>`));parts.push(`<line x1="${{sx(0)}}" y1="${{margin.top}}" x2="${{sx(0)}}" y2="${{margin.top+innerH}}" class="zero"/><line x1="${{margin.left}}" y1="${{margin.top+innerH}}" x2="${{margin.left+innerW}}" y2="${{margin.top+innerH}}" class="axis"/><line x1="${{margin.left}}" y1="${{margin.top}}" x2="${{margin.left}}" y2="${{margin.top+innerH}}" class="axis"/>`);sampleSeries.forEach((s,i)=>{{const style=sampleStyle(s.label,{{color:colors[s.condition]||DEFAULT_COLORS[i%DEFAULT_COLORS.length]}}),dash=dashAttr(style.type);parts.push(`<path d="${{pathD(s.profile,x,sx,sy)}}" fill="none" stroke="${{style.color}}" stroke-width="${{lineWidthValue(style.width,2)}}"${{dash}} stroke-opacity="${{alphaValue(style.alpha,.9)}}"><title>${{escText(s.label)}} - ${{escText(s.condition)}}</title></path>`)}});parts.push(`<text x="${{margin.left+6}}" y="${{margin.top+innerH-8}}" class="tick" fill="#94a3b8">${{motif.sites||0}}</text><text x="${{margin.left+innerW/2}}" y="${{height-9}}" class="axis-label" text-anchor="middle">${{escText(payload.x_label)}}</text><text x="14" y="${{margin.top+innerH/2}}" class="axis-label" text-anchor="middle" transform="rotate(-90 14 ${{margin.top+innerH/2}})">${{escText(payload.y_label)}}</text></svg>`);return parts.join('')}}
function legendHtml(){{return `<div class="legend-box">${{allSamples().filter(s=>sampleStyle(s.label).visible).map(s=>{{const style=sampleStyle(s.label,{{color:currentColors()[s.condition]}}),dash=lineDash(style.type);return `<div class="legend-row"><svg viewBox="0 0 32 8"><line x1="1" y1="4" x2="31" y2="4" stroke="${{style.color}}" stroke-width="${{lineWidthValue(style.width,2)}}" stroke-opacity="${{alphaValue(style.alpha,.9)}}" ${{dash?`stroke-dasharray="${{dash}}"`:''}}/></svg><span>${{escText(s.label)}}</span></div>`}}).join('')}}</div>`}}
function renderPlots(){{ensureSlots();plotGrid.className=`grid cols${{plotColsSel.value||2}}`;plotGrid.innerHTML=legendHtml()+slotPrefixes.map((prefix,idx)=>{{const motif=motifByPrefix(prefix);return `<div class="plot-card${{idx===activeSlot?' active':''}}${{groupAutoscalePanels.has(idx)?' grouped':''}}" data-card="${{idx}}"><div class="panel-tools">${{panelTfHtml(idx,prefix)}}<div class="panel-actions"><label><input type="checkbox" data-group-panel="${{idx}}" ${{groupAutoscalePanels.has(idx)?'checked':''}}> Autoscale group</label>${{samplePickerHtml(idx)}}<button data-download-panel="${{idx}}">Download SVG</button></div></div>${{drawAggregate(motif,idx)}}</div>`}}).join('');plotGrid.querySelectorAll('.plot-card').forEach(card=>card.addEventListener('click',ev=>{{if(ev.target.closest('button,details,input,label,summary,select'))return;activeSlot=Number(card.dataset.card);renderAll(false)}}));plotGrid.querySelectorAll('[data-panel-tf]').forEach(sel=>sel.addEventListener('change',()=>{{slotPrefixes[Number(sel.dataset.panelTf)]=sel.value;activeSlot=Number(sel.dataset.panelTf);groupAutoscaleDomain=null;renderAll(false)}}));plotGrid.querySelectorAll('[data-slot-sample]').forEach(inp=>inp.addEventListener('change',()=>{{const idx=Number(inp.dataset.slotSample),set=slotSamples[idx];if(inp.checked)set.add(inp.dataset.sample);else set.delete(inp.dataset.sample);if(!set.size)allSampleLabels().forEach(label=>set.add(label));groupAutoscaleDomain=null;renderAll(false)}}));plotGrid.querySelectorAll('[data-group-panel]').forEach(inp=>inp.addEventListener('change',()=>{{const idx=Number(inp.dataset.groupPanel);if(inp.checked)groupAutoscalePanels.add(idx);else groupAutoscalePanels.delete(idx);groupAutoscaleDomain=null;renderAll(false)}}));plotGrid.querySelectorAll('[data-download-panel]').forEach(btn=>btn.addEventListener('click',ev=>{{ev.stopPropagation();const svg=document.querySelector(`.aggregate-panel[data-panel="${{btn.dataset.downloadPanel}}"]`);if(svg)downloadBlob(svgBlob(svg),`plot_aggregate_panel_${{Number(btn.dataset.downloadPanel)+1}}.svg`)}}));drawSummary();renderSelectedMotifs()}}
function renderAll(refreshStyles=true){{renderColors();if(refreshStyles)renderSampleStyles();renderPlots();reportDetail.textContent=`${{payload.conditions.length}} groups - ${{payload.motifs.length}} TFs - ${{allSamples().length}} samples`}}
function styledSvgClone(svgNode){{const clone=svgNode.cloneNode(true);clone.setAttribute('xmlns','http://www.w3.org/2000/svg');clone.setAttribute('font-family','Arial,Helvetica,sans-serif');return clone}}
function svgBlob(svgNode){{return new Blob([new XMLSerializer().serializeToString(styledSvgClone(svgNode))],{{type:'image/svg+xml;charset=utf-8'}})}}
function downloadBlob(blob,filename){{const url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download=filename;document.body.appendChild(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(url),1000)}}
function legendSvg(x,y){{const rows=allSamples().filter(s=>sampleStyle(s.label).visible),h=Math.max(24,rows.length*15+10),w=150;let parts=[`<g transform="translate(${{x}},${{y}})"><rect width="${{w}}" height="${{h}}" rx="6" fill="#fff" stroke="#d9e2ec"/>`];rows.forEach((s,i)=>{{const style=sampleStyle(s.label,{{color:currentColors()[s.condition]}}),dash=lineDash(style.type),yy=13+i*15;parts.push(`<line x1="8" y1="${{yy}}" x2="42" y2="${{yy}}" stroke="${{style.color}}" stroke-width="${{lineWidthValue(style.width,2)}}" stroke-opacity="${{alphaValue(style.alpha,.9)}}" ${{dash?`stroke-dasharray="${{dash}}"`:''}}/><text x="48" y="${{yy+3}}" font-family="Arial,Helvetica,sans-serif" font-size="10" font-weight="900" fill="#334e68">${{escText(s.label)}}</text>`)}});parts.push('</g>');return parts.join('')}}
function downloadGrid(){{const svgs=[...document.querySelectorAll('.aggregate-panel')];if(!svgs.length)return;const w=340,h=340,cols=Number(plotColsSel.value)||2,rows=Math.ceil(svgs.length/cols),legendW=160,totalW=cols*w+legendW,totalH=Math.max(rows*h,120);let parts=[`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${{totalW}} ${{totalH}}" font-family="Arial,Helvetica,sans-serif"><rect width="100%" height="100%" fill="#fff"/>`];svgs.forEach((svg,i)=>{{const clone=styledSvgClone(svg);parts.push(`<g transform="translate(${{(i%cols)*w}},${{Math.floor(i/cols)*h}})">${{clone.innerHTML}}</g>`)}});parts.push(legendSvg(cols*w+6,6));parts.push('</svg>');downloadBlob(new Blob(parts,{{type:'image/svg+xml;charset=utf-8'}}),'plot_aggregate_grid.svg')}}
function downloadLogoPanel(){{const cards=[...motifGrid.querySelectorAll('.motif-card')],w=230,h=132,gap=8,cols=Math.min(4,Math.max(1,cards.length)),rows=Math.ceil(cards.length/cols);let parts=[`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${{cols*w+(cols-1)*gap}} ${{rows*h+(rows-1)*gap}}" font-family="Arial,Helvetica,sans-serif"><rect width="100%" height="100%" fill="#fff"/>`];slotPrefixes.forEach((prefix,i)=>{{const motif=motifByPrefix(prefix),x=(i%cols)*(w+gap),y=Math.floor(i/cols)*(h+gap),logo=(payload.logos||{{}})[prefix]||{{}};parts.push(`<g transform="translate(${{x}},${{y}})"><rect width="${{w}}" height="${{h}}" rx="7" fill="#fff" stroke="#d9e2ec"/><text x="8" y="18" font-size="12" font-weight="900" fill="#152133">${{escText(motifLabel(motif)).slice(0,30)}}</text>`);if(logo.png)parts.push(`<image x="12" y="28" width="${{w-24}}" height="90" preserveAspectRatio="xMidYMid meet" href="${{logo.png}}"/>`);else parts.push(`<text x="${{w/2}}" y="76" text-anchor="middle" font-size="11" font-weight="800" fill="#94a3b8">Motif logo unavailable</text>`);parts.push('</g>')}});parts.push('</svg>');downloadBlob(new Blob(parts,{{type:'image/svg+xml;charset=utf-8'}}),'plot_aggregate_motif_logo_panel.svg')}}
function updateGroupAutoscale(){{const panels=[...groupAutoscalePanels];if(!panels.length){{groupAutoscaleDomain=null;renderAll(false);return}}let vals=[];panels.forEach(idx=>{{const motif=motifByPrefix(slotPrefixes[idx]);samplesForMotif(motif,idx).forEach(s=>vals=vals.concat(s.profile.filter(Number.isFinite)))}});let rawMin=Math.min(...vals,0),rawMax=Math.max(...vals,1e-9),pad=Math.max((rawMax-rawMin||1)*.06,1e-6),step=niceStep((rawMax-rawMin+2*pad)/4),ymin=Math.floor((rawMin-pad)/step)*step,ymax=Math.ceil((rawMax+pad)/step)*step;groupAutoscaleDomain=[ymin,ymax];renderAll(false)}}
document.getElementById('download-grid').addEventListener('click',downloadGrid);document.getElementById('download-summary').addEventListener('click',()=>downloadBlob(svgBlob(waterfallChart),'plot_aggregate_barplot.svg'));document.getElementById('download-logo').addEventListener('click',downloadLogoPanel);document.getElementById('group-autoscale').addEventListener('click',updateGroupAutoscale);document.getElementById('reset-autoscale').addEventListener('click',()=>{{groupAutoscaleDomain=null;groupAutoscalePanels.clear();renderAll(false)}});[summarySort,summaryRows,plotColsSel].forEach(el=>el.addEventListener('change',()=>renderAll(false)));plotCountSel.addEventListener('change',()=>{{ensureSlots(false);groupAutoscaleDomain=null;renderAll(false)}});decodePayload().then(data=>{{payload=data;for(let i=1;i<=12;i++)plotCountSel.insertAdjacentHTML('beforeend',`<option value="${{i}}">${{i}}</option>`);plotCountSel.value=String(Math.min(12,initialPlotCount()));plotColsSel.value=plotCountSel.value==='1'?'1':(plotCountSel.value==='2'?'2':(plotCountSel.value==='6'?'3':'2'));ensureSlots(true);renderSummaryControls();renderAll()}}).catch(err=>{{reportDetail.textContent=`Could not open report payload: ${{err.message}}`}});
</script></body></html>"""
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    Path(output).write_text(document, encoding="utf-8")


def write_html(payload: dict, output: str | Path, title: str, default_layout: str = "2x3", show_summary: bool = True) -> None:
    """Write the column-oriented plot-aggregate HTML report."""

    payload = merge_payloads([payload])
    payload["default_layout"] = default_layout
    payload_b64 = _compressed_json_b64(payload)
    escaped_title = html.escape(title or "Aggregate motif footprint browser")
    document = f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>{escaped_title}</title><style>
:root{{--ink:#152133;--muted:#596579;--line:#d9e2ec;--grid:#e8eef5;--bg:#eef3f8;--accent:#173b73}}*{{box-sizing:border-box}}body{{margin:0;font-family:Arial,Helvetica,sans-serif;background:var(--bg);color:var(--ink);font-weight:700}}.wrap{{max-width:min(2400px,calc(100vw - 12px));margin:4px auto;padding:0 4px}}.panel,.card,.motif-col{{background:#fff;border:1px solid var(--line);border-radius:7px}}.panel{{box-shadow:0 12px 30px rgba(21,33,51,.10);overflow:hidden}}.head{{padding:8px 14px 6px;border-bottom:1px solid var(--line);background:#fff}}h1{{margin:0;font-size:22px;line-height:1.08;font-weight:900}}.sub{{margin:2px 0 0;color:var(--muted);font-size:11px}}.top-controls{{display:flex;flex-wrap:wrap;gap:6px;padding:5px 8px;border-bottom:1px solid var(--line);background:#fbfdff;align-items:center}}.card{{padding:5px 6px;min-width:0}}.top-controls .card{{display:flex;align-items:center;gap:8px;min-height:34px}}.top-controls .card:first-child{{flex:2 1 620px}}.top-controls .card:nth-child(2){{flex:0 0 390px}}.top-controls .card:nth-child(2) .control-row{{flex-wrap:nowrap}}.top-controls .card:nth-child(3){{flex:1 1 260px}}.section-title{{font-size:10px;line-height:1.05;text-transform:uppercase;letter-spacing:.08em;color:#728197;margin:0;font-weight:900;white-space:nowrap}}.button-stack{{display:flex;flex-wrap:wrap;align-items:center;gap:5px}}button{{border:1px solid #b8c5d6;background:#fff;color:var(--accent);border-radius:6px;padding:4px 7px;font-family:Arial,Helvetica,sans-serif;font-size:10px;font-weight:900;cursor:pointer}}button:hover{{background:#f2f6fb}}label{{font-size:10px;color:#52606d;text-transform:uppercase;letter-spacing:.06em;font-weight:900}}select,input{{border:1px solid #cbd5e1;border-radius:6px;background:white;color:var(--ink);font-family:Arial,Helvetica,sans-serif;font-size:11px;font-weight:800;padding:4px 6px}}input[type=checkbox]{{width:13px;height:13px;padding:0}}.control-row,.color-grid{{display:flex;flex-wrap:wrap;align-items:center;gap:6px}}.color-row{{display:flex;align-items:center;gap:5px;font-size:11px;color:#334e68;font-weight:800;border:1px solid #e6edf5;border-radius:999px;padding:3px 6px;background:#fbfdff}}.color-row input{{width:26px;height:18px;padding:0}}.motif-stage{{padding:7px;background:#f8fbff;overflow-x:auto}}.motif-columns{{display:grid;grid-auto-flow:column;grid-auto-columns:minmax(280px,1fr);gap:7px;align-items:start;width:100%;min-width:0}}.motif-col{{display:grid;grid-template-rows:auto auto auto auto auto;gap:6px;padding:6px;min-width:280px}}.motif-col.active{{border-color:#173b73;box-shadow:0 0 0 2px rgba(23,59,115,.10)}}.motif-col.grouped{{border-color:#0f766e;box-shadow:0 0 0 2px rgba(15,118,110,.12)}}.motif-select{{width:100%;min-width:0}}.motif-logo{{height:96px;display:flex;align-items:center;justify-content:center;overflow:hidden;border:1px solid #edf2f7;border-radius:6px;background:#fff}}.motif-logo img{{max-width:100%;max-height:92px}}.logo-empty{{font-size:11px;color:#94a3b8}}.plot-controls{{display:grid;grid-template-columns:1fr 1fr;gap:5px;align-items:center}}.plot-controls button{{padding:4px 6px}}.sample-panel{{border:1px solid #e6edf5;border-radius:6px;padding:5px;background:#fff;display:grid;gap:5px}}.sample-group{{display:grid;gap:3px}}.sample-group-title{{display:flex;align-items:center;gap:5px;font-size:10px;color:#334e68;font-weight:900}}.sample-dot{{width:8px;height:8px;border-radius:99px;display:inline-block}}.sample-style-row{{display:grid;grid-template-columns:14px minmax(56px,1fr) 24px 42px 42px 48px;gap:3px;align-items:center}}.sample-style-name{{font-size:10px;color:#334e68;font-weight:900;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}}.sample-style-row input[type=color]{{width:24px;height:19px;padding:1px}}.sample-style-row input[type=number]{{width:42px;min-width:42px;padding:3px}}.sample-style-row select{{width:48px;padding:3px}}.legend-panel{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:2px;border:1px solid #e6edf5;border-radius:6px;padding:5px;background:#fff}}.legend-row{{display:grid;grid-template-columns:34px minmax(0,1fr);gap:5px;align-items:center;font-size:10px;color:#334e68;font-weight:900}}.legend-row span{{overflow:hidden;text-overflow:ellipsis;white-space:nowrap}}svg{{width:100%;height:auto;display:block;background:#fff}}.aggregate-panel{{border:1px solid #edf2f7;border-radius:6px}}.axis{{stroke:#3b4552;stroke-width:1.2}}.grid-line{{stroke:var(--grid);stroke-width:1}}.zero{{stroke:#677386;stroke-width:1.25;stroke-dasharray:4 4}}.tick{{font-family:Arial,Helvetica,sans-serif;font-size:10px;fill:#596579;font-weight:800}}.axis-label{{font-family:Arial,Helvetica,sans-serif;font-size:10px;fill:#152133;font-weight:900}}.plot-title{{font-family:Arial,Helvetica,sans-serif;font-size:12px;fill:#152133;font-weight:900}}(max-width:760px){{.top-controls .card{{flex:1 1 100%}}.motif-columns{{grid-auto-columns:minmax(280px,92vw)}}}}
</style></head><body><div class="wrap"><section class="panel"><header class="head"><h1>{escaped_title}</h1><p class="sub" id="report-detail"></p></header><div class="top-controls"><div class="card"><p class="section-title">Export editable SVG</p><div class="button-stack"><button id="download-logo">Download motif logo panel</button><button id="download-grid">Download motif aggregate panel</button><button id="download-combined">Download combined panel</button></div></div><div class="card"><p class="section-title">Layout</p><div class="control-row"><label>Plots <select id="plot-count"></select></label><button id="group-autoscale">Group autoscale</button><button id="reset-autoscale">Reset autoscale</button></div></div><div class="card"><p class="section-title">Groups</p><div id="color-controls" class="color-grid"></div></div></div><main class="motif-stage"><div id="motif-columns" class="motif-columns"></div></main></section></div><script>
const reportPayloadB64="{payload_b64}",DEFAULT_COLORS={json.dumps(DEFAULT_COLORS)};let payload=null,slotPrefixes=[],slotSampleStyles=[],activeSlot=0,groupAutoscalePanels=new Set(),groupAutoscaleDomain=null;const reportDetail=document.getElementById('report-detail'),plotCountSel=document.getElementById('plot-count'),motifColumns=document.getElementById('motif-columns'),colorControls=document.getElementById('color-controls');
function decodePayload(){{return new Promise((resolve,reject)=>{{try{{const binary=atob(reportPayloadB64),bytes=new Uint8Array(binary.length);for(let i=0;i<binary.length;i++)bytes[i]=binary.charCodeAt(i);new Response(new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip'))).text().then(t=>resolve(JSON.parse(t))).catch(reject)}}catch(e){{reject(e)}}}})}}
function escText(v){{return String(v??'').replace(/[&<>"']/g,c=>({{'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}}[c]))}}
function motifLabel(m){{return m.motif_id?`${{m.name}} (${{m.motif_id}})`:m.name}}
function motifByPrefix(prefix){{return (payload.motifs||[]).find(m=>m.prefix===prefix)||payload.motifs[0]}}
function currentColors(){{const colors={{...(payload.colors||{{}})}};colorControls.querySelectorAll('input').forEach(inp=>colors[inp.dataset.condition]=inp.value);return colors}}
function allSamples(){{const seen=new Map();(payload.motifs||[]).forEach(m=>(m.series||[]).filter(s=>s.kind==='sample').forEach(s=>{{if(!seen.has(s.label))seen.set(s.label,{{label:s.label,condition:s.condition}})}}));return [...seen.values()]}}
function allSampleLabels(){{return allSamples().map(s=>s.label)}}
function initialPlotCount(){{return Math.min(6,Math.max(1,(payload.motifs||[]).length))}}
function defaultStyle(sample){{return {{visible:true,color:currentColors()[sample.condition]||'#2563eb',alpha:.9,width:.7,type:'solid'}}}}
function ensureSlots(refill=false){{const n=Number(plotCountSel.value)||initialPlotCount(),labels=allSampleLabels();if(refill||!slotPrefixes.length)slotPrefixes=(payload.motifs||[]).slice(0,n).map(m=>m.prefix);while(slotPrefixes.length<n)slotPrefixes.push((payload.motifs||[])[slotPrefixes.length%Math.max(1,payload.motifs.length)]?.prefix);slotPrefixes=slotPrefixes.slice(0,n);while(slotSampleStyles.length<n)slotSampleStyles.push(Object.fromEntries(allSamples().map(s=>[s.label,defaultStyle(s)])));slotSampleStyles=slotSampleStyles.slice(0,n);slotSampleStyles.forEach(map=>allSamples().forEach(s=>{{if(!map[s.label])map[s.label]=defaultStyle(s)}}));groupAutoscalePanels=new Set([...groupAutoscalePanels].filter(i=>i<n));activeSlot=Math.max(0,Math.min(activeSlot,n-1))}}
function lineDash(type){{return type==='dash'?'6 4':(type==='dot'?'1.2 3':'')}}
function dashAttr(type){{const dash=lineDash(type);return dash?` stroke-dasharray="${{dash}}"`:''}}
function lineWidthValue(v,fallback){{const n=Number(v);return Number.isFinite(n)&&n>0?Math.min(8,Math.max(.1,n)):fallback}}
function alphaValue(v,fallback){{const n=Number(v);return Number.isFinite(n)?Math.min(1,Math.max(.05,n)):fallback}}
function sampleStyle(idx,label){{const s=slotSampleStyles[idx]?.[label]||{{}};return {{visible:s.visible!==false,color:s.color||'#2563eb',alpha:s.alpha??.9,width:s.width||.7,type:s.type||'solid'}}}}
function samplesForMotif(motif,idx){{return (motif.series||[]).filter(s=>s.kind==='sample'&&sampleStyle(idx,s.label).visible)}}
function renderColors(){{const card=colorControls.closest('.card');if(!payload.groups_defined){{colorControls.innerHTML='';if(card)card.style.display='none';return}}if(card)card.style.display='flex';colorControls.innerHTML=(payload.conditions||[]).map(c=>`<label class="color-row"><span>${{escText(c)}}</span><input type="color" data-condition="${{escText(c)}}" value="${{(payload.colors||{{}})[c]||'#64748b'}}"></label>`).join('');colorControls.querySelectorAll('input').forEach(inp=>inp.addEventListener('input',()=>{{slotSampleStyles.forEach(map=>allSamples().forEach(s=>{{if(s.condition===inp.dataset.condition)map[s.label].color=inp.value}}));renderAll(false)}}))}}
function logoHtml(prefix){{const logo=(payload.logos||{{}})[prefix]||{{}};return logo.png?`<img alt="Motif logo" src="${{logo.png}}">`:'<span class="logo-empty">Motif logo unavailable</span>'}}
function sampleControls(idx){{const colors=currentColors(),samples=allSamples(),groups=payload.groups_defined?(payload.conditions||[]):['__all__'];return groups.map(cond=>{{const groupSamples=cond==='__all__'?samples:samples.filter(s=>s.condition===cond),title=cond==='__all__'?'Samples':`${{escText(cond)}} samples`,dot=cond==='__all__'?'#64748b':(colors[cond]||'#64748b');return `<div class="sample-group"><div class="sample-group-title"><span class="sample-dot" style="background:${{dot}}"></span>${{title}}</div>${{groupSamples.map(s=>{{const st=sampleStyle(idx,s.label);return `<label class="sample-style-row"><input data-visible="${{idx}}:${{escText(s.label)}}" type="checkbox" ${{st.visible?'checked':''}}><span class="sample-style-name">${{escText(s.label)}}</span><input data-color="${{idx}}:${{escText(s.label)}}" type="color" value="${{st.color}}"><input data-alpha="${{idx}}:${{escText(s.label)}}" type="number" min="0.05" max="1" step="0.05" value="${{st.alpha}}"><input data-width="${{idx}}:${{escText(s.label)}}" type="number" min="0.2" max="5" step="0.1" value="${{st.width}}"><select data-type="${{idx}}:${{escText(s.label)}}"><option value="solid"${{st.type==='solid'?' selected':''}}>Solid</option><option value="dash"${{st.type==='dash'?' selected':''}}>Dash</option><option value="dot"${{st.type==='dot'?' selected':''}}>Dot</option></select></label>`}}).join('')}}</div>`}}).join('')}}
function legendHtml(idx){{return allSamples().filter(s=>sampleStyle(idx,s.label).visible).map(s=>{{const st=sampleStyle(idx,s.label),dash=lineDash(st.type);return `<div class="legend-row"><svg viewBox="0 0 34 8"><line x1="1" y1="4" x2="33" y2="4" stroke="${{st.color}}" stroke-width="${{lineWidthValue(st.width,2)}}" stroke-opacity="${{alphaValue(st.alpha,.9)}}" ${{dash?`stroke-dasharray="${{dash}}"`:''}}/></svg><span>${{escText(s.label)}}</span></div>`}}).join('')}}
function niceTicks(min,max,n){{const out=[];for(let i=0;i<n;i++)out.push(min+(max-min)*(i/Math.max(1,n-1)));return out}}
function fmt(v){{const a=Math.abs(v);if(!Number.isFinite(v))return'';if(a===0)return'0';if(a>=1)return v.toFixed(1).replace('-0.0','0.0');if(a>=.01)return v.toFixed(2).replace('-0.00','0.00');if(a>=.001)return v.toFixed(3).replace('-0.000','0.000');return v.toExponential(1).replace('-0.0e+0','0')}}
function pathD(profile,x,sx,sy){{return profile.map((y,i)=>`${{i?'L':'M'}}${{sx(x[i]).toFixed(2)}},${{sy(y).toFixed(2)}}`).join(' ')}}
function yDomainForPanel(idx){{if(groupAutoscaleDomain&&groupAutoscalePanels.has(idx))return groupAutoscaleDomain;const motif=motifByPrefix(slotPrefixes[idx]),vals=samplesForMotif(motif,idx).flatMap(s=>s.profile).filter(Number.isFinite);let rawMin=Math.min(...vals,0),rawMax=Math.max(...vals,1e-9),pad=Math.max((rawMax-rawMin||1)*.06,1e-6),step=niceStep((rawMax-rawMin+2*pad)/4),ymin=Math.floor((rawMin-pad)/step)*step,ymax=Math.ceil((rawMax+pad)/step)*step;return [ymin,ymax]}}
function drawAggregate(motif,idx){{const samples=samplesForMotif(motif,idx),x=payload.x,width=340,height=340,margin={{top:28,right:14,bottom:38,left:50}},innerW=width-margin.left-margin.right,innerH=height-margin.top-margin.bottom,sampleSeries=[...samples].sort((a,b)=>(Number(b.fp_score??b.avg_score??0)-Number(a.fp_score??a.avg_score??0))),domain=yDomainForPanel(idx),ymin=domain[0],ymax=domain[1],sx=v=>margin.left+((v-x[0])/(x[x.length-1]-x[0]||1))*innerW,sy=v=>margin.top+innerH-((v-ymin)/(ymax-ymin||1))*innerH,xTicks=[x[0],0,x[x.length-1]],yTicks=niceTicks(ymin,ymax,4);let parts=[`<svg class="aggregate-panel" data-panel="${{idx}}" viewBox="0 0 ${{width}} ${{height}}"><style>.axis{{stroke:#3b4552;stroke-width:1.2}}.grid-line{{stroke:#e8eef5;stroke-width:1}}.zero{{stroke:#677386;stroke-width:1.25;stroke-dasharray:4 4}}.tick{{font-family:Arial,Helvetica,sans-serif;font-size:10px;fill:#596579;font-weight:800}}.axis-label{{font-family:Arial,Helvetica,sans-serif;font-size:10px;fill:#152133;font-weight:900}}.plot-title{{font-family:Arial,Helvetica,sans-serif;font-size:12px;fill:#152133;font-weight:900}}</style><rect width="${{width}}" height="${{height}}" fill="#fff"/><text x="${{width/2}}" y="16" text-anchor="middle" class="plot-title">${{escText(motifLabel(motif))}}</text>`];yTicks.forEach(v=>parts.push(`<text x="${{margin.left-7}}" y="${{sy(v)+3}}" class="tick" text-anchor="end">${{fmt(v)}}</text>`));xTicks.forEach(v=>parts.push(`<line x1="${{sx(v)}}" y1="${{margin.top}}" x2="${{sx(v)}}" y2="${{margin.top+innerH}}" class="grid-line"/><text x="${{sx(v)}}" y="${{margin.top+innerH+17}}" class="tick" text-anchor="middle">${{v}}</text>`));parts.push(`<line x1="${{sx(0)}}" y1="${{margin.top}}" x2="${{sx(0)}}" y2="${{margin.top+innerH}}" class="zero"/><line x1="${{margin.left}}" y1="${{margin.top+innerH}}" x2="${{margin.left+innerW}}" y2="${{margin.top+innerH}}" class="axis"/><line x1="${{margin.left}}" y1="${{margin.top}}" x2="${{margin.left}}" y2="${{margin.top+innerH}}" class="axis"/>`);sampleSeries.forEach(s=>{{const st=sampleStyle(idx,s.label);parts.push(`<path d="${{pathD(s.profile,x,sx,sy)}}" fill="none" stroke="${{st.color}}" stroke-width="${{lineWidthValue(st.width,2)}}"${{dashAttr(st.type)}} stroke-opacity="${{alphaValue(st.alpha,.9)}}"><title>${{escText(s.label)}} - ${{escText(s.condition)}}</title></path>`)}});parts.push(`<text x="${{margin.left+6}}" y="${{margin.top+innerH-8}}" class="tick" fill="#94a3b8">${{motif.sites||0}}</text><text x="${{margin.left+innerW/2}}" y="${{height-9}}" class="axis-label" text-anchor="middle">${{escText(payload.x_label)}}</text><text x="14" y="${{margin.top+innerH/2}}" class="axis-label" text-anchor="middle" transform="rotate(-90 14 ${{margin.top+innerH/2}})">${{escText(payload.y_label)}}</text></svg>`);return parts.join('')}}
function renderColumns(){{ensureSlots();motifColumns.style.gridTemplateColumns=`repeat(${{slotPrefixes.length}}, minmax(280px, 1fr))`;motifColumns.innerHTML=slotPrefixes.map((prefix,idx)=>{{const motif=motifByPrefix(prefix);return `<section class="motif-col${{idx===activeSlot?' active':''}}${{groupAutoscalePanels.has(idx)?' grouped':''}}" data-col="${{idx}}"><select class="motif-select" data-motif="${{idx}}">${{(payload.motifs||[]).map(m=>`<option value="${{escText(m.prefix)}}" ${{m.prefix===prefix?'selected':''}}>${{escText(motifLabel(m))}}</option>`).join('')}}</select><div class="motif-logo">${{logoHtml(prefix)}}</div><div class="plot-controls"><label><input type="checkbox" data-group="${{idx}}" ${{groupAutoscalePanels.has(idx)?'checked':''}}> Autoscale group</label><button data-download-panel="${{idx}}">Download SVG</button></div><div class="sample-panel">${{sampleControls(idx)}}</div><div class="legend-panel">${{legendHtml(idx)}}</div>${{drawAggregate(motif,idx)}}</section>`}}).join('');bindColumnEvents()}}
function bindColumnEvents(){{motifColumns.querySelectorAll('[data-col]').forEach(col=>col.addEventListener('click',ev=>{{if(ev.target.closest('select,input,button'))return;activeSlot=Number(col.dataset.col);renderAll(false)}}));motifColumns.querySelectorAll('[data-motif]').forEach(sel=>sel.addEventListener('change',()=>{{slotPrefixes[Number(sel.dataset.motif)]=sel.value;activeSlot=Number(sel.dataset.motif);groupAutoscaleDomain=null;renderAll(false)}}));motifColumns.querySelectorAll('[data-group]').forEach(inp=>inp.addEventListener('change',()=>{{const idx=Number(inp.dataset.group);if(inp.checked)groupAutoscalePanels.add(idx);else groupAutoscalePanels.delete(idx);groupAutoscaleDomain=null;renderAll(false)}}));motifColumns.querySelectorAll('[data-visible]').forEach(inp=>inp.addEventListener('change',()=>{{const [idx,label]=inp.dataset.visible.split(':');slotSampleStyles[Number(idx)][label]={{...(slotSampleStyles[Number(idx)][label]||{{}}),visible:inp.checked}};renderAll(false)}}));motifColumns.querySelectorAll('[data-color]').forEach(inp=>inp.addEventListener('input',()=>{{const [idx,label]=inp.dataset.color.split(':');slotSampleStyles[Number(idx)][label]={{...(slotSampleStyles[Number(idx)][label]||{{}}),color:inp.value}};renderAll(false)}}));motifColumns.querySelectorAll('[data-alpha]').forEach(inp=>inp.addEventListener('input',()=>{{const [idx,label]=inp.dataset.alpha.split(':');slotSampleStyles[Number(idx)][label]={{...(slotSampleStyles[Number(idx)][label]||{{}}),alpha:alphaValue(inp.value,.9)}};renderAll(false)}}));motifColumns.querySelectorAll('[data-width]').forEach(inp=>inp.addEventListener('input',()=>{{const [idx,label]=inp.dataset.width.split(':');slotSampleStyles[Number(idx)][label]={{...(slotSampleStyles[Number(idx)][label]||{{}}),width:lineWidthValue(inp.value,2)}};renderAll(false)}}));motifColumns.querySelectorAll('[data-type]').forEach(sel=>sel.addEventListener('change',()=>{{const [idx,label]=sel.dataset.type.split(':');slotSampleStyles[Number(idx)][label]={{...(slotSampleStyles[Number(idx)][label]||{{}}),type:sel.value}};renderAll(false)}}));motifColumns.querySelectorAll('[data-download-panel]').forEach(btn=>btn.addEventListener('click',ev=>{{ev.stopPropagation();const svg=document.querySelector(`.aggregate-panel[data-panel="${{btn.dataset.downloadPanel}}"]`);if(svg)downloadBlob(svgBlob(svg),`plot_aggregate_panel_${{Number(btn.dataset.downloadPanel)+1}}.svg`)}}))}}
function renderAll(refreshColors=true){{if(refreshColors)renderColors();renderColumns();const groupText=payload.groups_defined?`${{payload.conditions.length}} groups`:'ungrouped samples';reportDetail.textContent=`${{groupText}} - ${{payload.motifs.length}} TFs - ${{allSamples().length}} samples`}}
function styledSvgClone(svgNode){{const clone=svgNode.cloneNode(true);clone.setAttribute('xmlns','http://www.w3.org/2000/svg');clone.setAttribute('font-family','Arial,Helvetica,sans-serif');return clone}}
function svgBlob(svgNode){{return new Blob([new XMLSerializer().serializeToString(styledSvgClone(svgNode))],{{type:'image/svg+xml;charset=utf-8'}})}}
function downloadBlob(blob,filename){{const url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download=filename;document.body.appendChild(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(url),1000)}}
function downloadGrid(){{const svgs=[...document.querySelectorAll('.aggregate-panel')],w=340,h=340,gap=8,totalW=svgs.length*w+(svgs.length-1)*gap,totalH=h;let parts=[`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${{totalW}} ${{totalH}}" font-family="Arial,Helvetica,sans-serif"><rect width="100%" height="100%" fill="#fff"/>`];svgs.forEach((svg,i)=>{{const clone=styledSvgClone(svg);parts.push(`<g transform="translate(${{i*(w+gap)}},0)">${{clone.innerHTML}}</g>`)}});parts.push('</svg>');downloadBlob(new Blob(parts,{{type:'image/svg+xml;charset=utf-8'}}),'plot_aggregate_grid.svg')}}
function downloadLogoPanel(){{const w=230,h=132,gap=8,totalW=slotPrefixes.length*w+(slotPrefixes.length-1)*gap;let parts=[`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${{totalW}} ${{h}}" font-family="Arial,Helvetica,sans-serif"><rect width="100%" height="100%" fill="#fff"/>`];slotPrefixes.forEach((prefix,i)=>{{const motif=motifByPrefix(prefix),x=i*(w+gap),logo=(payload.logos||{{}})[prefix]||{{}};parts.push(`<g transform="translate(${{x}},0)"><rect width="${{w}}" height="${{h}}" rx="7" fill="#fff" stroke="#d9e2ec"/><text x="8" y="18" font-size="12" font-weight="900" fill="#152133">${{escText(motifLabel(motif)).slice(0,30)}}</text>`);if(logo.png)parts.push(`<image x="12" y="28" width="${{w-24}}" height="90" preserveAspectRatio="xMidYMid meet" href="${{logo.png}}"/>`);else parts.push(`<text x="${{w/2}}" y="76" text-anchor="middle" font-size="11" font-weight="800" fill="#94a3b8">Motif logo unavailable</text>`);parts.push('</g>')}});parts.push('</svg>');downloadBlob(new Blob(parts,{{type:'image/svg+xml;charset=utf-8'}}),'plot_aggregate_motif_logo_panel.svg')}}
function updateGroupAutoscale(){{const panels=[...groupAutoscalePanels];if(!panels.length){{groupAutoscaleDomain=null;renderAll(false);return}}let vals=[];panels.forEach(idx=>{{const motif=motifByPrefix(slotPrefixes[idx]);samplesForMotif(motif,idx).forEach(s=>vals=vals.concat(s.profile.filter(Number.isFinite)))}});let rawMin=Math.min(...vals,0),rawMax=Math.max(...vals,1e-9),pad=Math.max((rawMax-rawMin||1)*.06,1e-6),step=niceStep((rawMax-rawMin+2*pad)/4),ymin=Math.floor((rawMin-pad)/step)*step,ymax=Math.ceil((rawMax+pad)/step)*step;groupAutoscaleDomain=[ymin,ymax];renderAll(false)}}
document.getElementById('download-grid').addEventListener('click',downloadGrid);document.getElementById('download-logo').addEventListener('click',downloadLogoPanel);document.getElementById('download-combined').addEventListener('click',downloadGrid);document.getElementById('group-autoscale').addEventListener('click',updateGroupAutoscale);document.getElementById('reset-autoscale').addEventListener('click',()=>{{groupAutoscaleDomain=null;groupAutoscalePanels.clear();renderAll(false)}});plotCountSel.addEventListener('change',()=>{{ensureSlots(false);groupAutoscaleDomain=null;renderAll(false)}});decodePayload().then(data=>{{payload=data;for(let i=1;i<=12;i++)plotCountSel.insertAdjacentHTML('beforeend',`<option value="${{i}}">${{i}}</option>`);plotCountSel.value=String(initialPlotCount());ensureSlots(true);renderAll()}}).catch(err=>{{reportDetail.textContent=`Could not open report payload: ${{err.message}}`}});
</script></body></html>"""
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    Path(output).write_text(document, encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Create an interactive aggregate HTML report from match-motifs or embedded diff-footprints outputs.")
    parser.add_argument("--manifest", help="TSV with sample, signal, and match_dir columns.")
    parser.add_argument("--input-html", nargs="*", default=[], help="Existing aggregate/diff-footprints HTML report(s) with embedded reportPayloadB64 payloads.")
    parser.add_argument("--output", required=True, help="Output self-contained HTML file.")
    parser.add_argument("--flank", type=int, default=100, help="Flank around motif centers for aggregate profiles (default: 100).")
    parser.add_argument("--top-n", type=int, default=30, help="Number of motifs to preload from manifest mode (default: 30).")
    parser.add_argument("--motifs", nargs="*", help="Motif prefixes, names, or IDs to preload from manifest mode.")
    parser.add_argument("--site-set", choices=["bound", "all", "unbound"], default="bound", help="Motif-site BED set to use from match directories in manifest mode (default: bound).")
    parser.add_argument("--normalization", choices=["none", "sample-quantile", "condition-quantile"], default="none", help="Profile scaling for manifest mode (default: none).")
    parser.add_argument("--default-layout", choices=["1x1", "1x2", "2x2", "2x3"], default="2x2", help="Initial panel grid layout (default: 2x2).")
    parser.add_argument("--title", default="Aggregate motif footprint browser")
    parser.add_argument("--hide-summary", action="store_true", help="Hide the TF site summary sidebar in the HTML report.")
    args = parser.parse_args(argv)
    payloads = []
    if args.manifest:
        payloads.append(
            build_payload(
                _read_manifest(args.manifest),
                flank=max(1, args.flank),
                top_n=max(1, args.top_n),
                normalization=args.normalization,
                motif_names=args.motifs,
                site_set=args.site_set,
            )
        )
    for path in args.input_html:
        payloads.append(read_embedded_payload(path))
    if not payloads:
        parser.error("provide --manifest and/or --input-html")
    payload = merge_payloads(payloads)
    write_html(payload, args.output, args.title, default_layout=args.default_layout, show_summary=not args.hide_summary)
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
