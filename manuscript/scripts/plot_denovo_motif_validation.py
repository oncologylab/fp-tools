#!/usr/bin/env python
"""Plot de novo motif validation from a two-condition replicate experiment."""

from __future__ import annotations

import argparse
import base64
import gzip
import json
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.font_manager import FontProperties
from matplotlib.patches import PathPatch
from matplotlib.textpath import TextPath
from matplotlib.transforms import Affine2D
import numpy as np
import pandas as pd
import pyBigWig

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

REPO_ROOT = SCRIPT_DIR.parents[1]
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from figure_style import apply_style, bold_all_text  # noqa: E402
from fp_tools.utils.motifs import MotifList  # noqa: E402


DEFAULT_VALIDATION_DIR = (
    REPO_ROOT
    / "data/public/processed/encode_k562_hepg2_atac_replicates/fp_tools/denovo_motif_validation_maxcover_n250"
)
DEFAULT_JASPAR = (
    REPO_ROOT
    / "data/public/raw/jaspar/2026/JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
)
CONDITION_COLORS = {
    "K562": "#1f77b4",
    "HepG2": "#d62728",
    "Bcell": "#1f77b4",
    "Tcell": "#d62728",
}
DEFAULT_CONDITIONS = ("K562", "HepG2")
ACTIVE_CONDITIONS = DEFAULT_CONDITIONS
SITE_WINDOW_BP = 10
N_AGGREGATE_EXAMPLES = 8
TN5_ADAPTER_CORE = "CTGTCTCTTATACACATCT"
DISPLAY_SAMPLE = "K562_rep1"
DISPLAY_SAMPLE_BASENAME = "K562_rep1.ENCFF077FBI_corrected.bw"
MAX_RANKING_SITES_PER_MOTIF = 1000
N_AGGREGATE_PER_CONDITION = 4
LOGO_COLORS = {"A": "#2f9e44", "C": "#1971c2", "G": "#f08c00", "T": "#d6336c"}
LOGO_FONT = FontProperties(family="DejaVu Sans", weight="bold")
rcParams["svg.fonttype"] = "none"


def load_report_payload(report_html: Path) -> dict:
    text = report_html.read_text(encoding="utf-8")
    match = re.search(r'reportPayloadB64="([^"]+)"', text)
    if match is None:
        raise ValueError(f"Could not find reportPayloadB64 in {report_html}")
    return json.loads(gzip.decompress(base64.b64decode(match.group(1))).decode("utf-8"))


def jaspar_name_map(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    motifs = MotifList().from_file(str(path))
    return {motif.id: motif.name for motif in motifs}


def load_streme_motifs(validation_dir: Path, id_to_name: dict[str, str]) -> pd.DataFrame:
    rows = []
    directions = [
        (
            f"{ACTIVE_CONDITIONS[0]}_vs_{ACTIVE_CONDITIONS[1]}_streme",
            ACTIVE_CONDITIONS[0],
            f"{ACTIVE_CONDITIONS[0]} candidates",
        ),
        (
            f"{ACTIVE_CONDITIONS[1]}_vs_{ACTIVE_CONDITIONS[0]}_streme",
            ACTIVE_CONDITIONS[1],
            f"{ACTIVE_CONDITIONS[1]} candidates",
        ),
    ]
    for direction, source_prefix, label in directions:
        path = validation_dir / "motifs" / direction / "motif_summary.tsv"
        if not path.exists():
            continue
        df = pd.read_csv(path, sep="\t")
        discovered = df[df["source"].eq("MEME")].copy()
        tomtom = df[df["source"].eq("Tomtom")].copy()
        tomtom["q_value"] = pd.to_numeric(tomtom.get("q_value"), errors="coerce")
        best_by_motif = tomtom.sort_values("q_value").drop_duplicates("motif_id").set_index("motif_id")
        for _, row in discovered.iterrows():
            motif_id = str(row["motif_id"])
            motif_number = motif_id.split("-", 1)[0]
            name = f"{source_prefix}_denovo_{motif_number}"
            best = best_by_motif.loc[motif_id] if motif_id in best_by_motif.index else None
            if best is not None and pd.notna(best["q_value"]):
                target_id = str(best["target_id"])
                q_value = float(best["q_value"])
                target_name = id_to_name.get(target_id, target_id)
                match = f"{target_name} ({target_id})" if q_value <= 0.05 else "no confident match"
            else:
                target_id = ""
                target_name = ""
                q_value = float("nan")
                match = "no confident match"
            rows.append(
                {
                    "name": name,
                    "direction": label,
                    "source_condition": source_prefix,
                    "de_novo_motif": motif_id,
                    "consensus": str(row["consensus"]),
                    "sites": int(float(row["sites"])),
                    "e_value": float(row["e_value"]),
                    "target_id": target_id,
                    "target_name": target_name,
                    "tomtom_q_value": q_value,
                    "tomtom_label": match,
                    "confident_tomtom": bool(pd.notna(q_value) and q_value <= 0.05),
                }
            )
    return pd.DataFrame(rows)


def parse_streme_probability_matrices(path: Path) -> dict[str, np.ndarray]:
    if not path.exists():
        return {}
    lines = path.read_text(encoding="utf-8").splitlines()
    matrices: dict[str, np.ndarray] = {}
    motif_id: str | None = None
    idx = 0
    while idx < len(lines):
        line = lines[idx].strip()
        if line.startswith("MOTIF "):
            parts = line.split()
            motif_id = parts[1] if len(parts) > 1 else None
            idx += 1
            continue
        if motif_id and line.startswith("letter-probability matrix:"):
            match = re.search(r"\bw=\s*(\d+)", line)
            width = int(match.group(1)) if match else 0
            matrix = []
            for row_line in lines[idx + 1 : idx + 1 + width]:
                values = [float(value) for value in row_line.split()[:4]]
                if len(values) == 4:
                    matrix.append(values)
            if matrix:
                matrices[motif_id] = np.asarray(matrix, dtype=float)
            idx += width + 1
            continue
        idx += 1
    return matrices


def load_streme_probability_matrices(validation_dir: Path) -> dict[str, np.ndarray]:
    matrices: dict[str, np.ndarray] = {}
    directions = [
        (f"{ACTIVE_CONDITIONS[0]}_vs_{ACTIVE_CONDITIONS[1]}_streme", ACTIVE_CONDITIONS[0]),
        (f"{ACTIVE_CONDITIONS[1]}_vs_{ACTIVE_CONDITIONS[0]}_streme", ACTIVE_CONDITIONS[1]),
    ]
    for direction, source_prefix in directions:
        path = validation_dir / "motifs" / direction / "streme" / "streme.txt"
        for motif_id, matrix in parse_streme_probability_matrices(path).items():
            motif_number = motif_id.split("-", 1)[0]
            matrices[f"{source_prefix}_denovo_{motif_number}"] = matrix
    return matrices


def center_minus_flank_score(profile: np.ndarray, xvals: np.ndarray, center_bp: int = 10, edge_bp: int = 20) -> float:
    center = profile[np.abs(xvals) <= center_bp]
    left = profile[xvals <= xvals.min() + edge_bp]
    right = profile[xvals >= xvals.max() - edge_bp]
    flanks = np.concatenate([left, right])
    if center.size == 0 or flanks.size == 0:
        return float("nan")
    return float(np.nanmean(center) - np.nanmean(flanks))


def center_flank_score(profile: np.ndarray, xvals: np.ndarray, center_bp: int = 10, edge_bp: int = 20) -> float:
    return center_minus_flank_score(profile, xvals, center_bp=center_bp, edge_bp=edge_bp)


def footprint_shape_score(profile: np.ndarray, xvals: np.ndarray, center_bp: int = 10, flank_min_bp: int = 20, flank_max_bp: int = 50) -> float:
    center = profile[np.abs(xvals) <= center_bp]
    flanks = profile[(np.abs(xvals) >= flank_min_bp) & (np.abs(xvals) <= flank_max_bp)]
    if center.size == 0 or flanks.size == 0:
        return float("nan")
    return float(np.nanmean(flanks) - np.nanmean(center))


def motif_lookup(payload: dict) -> dict[str, dict]:
    return {motif["name"]: motif for motif in payload.get("aggregate", {}).get("motifs", [])}


def aggregate_scores(payload: dict, xvals: np.ndarray) -> pd.DataFrame:
    rows = []
    for motif in payload.get("aggregate", {}).get("motifs", []):
        if "_denovo_" not in str(motif.get("name", "")):
            continue
        scores = {}
        for condition in motif.get("conditions", []):
            profile = np.asarray(condition.get("profile", []), dtype=float)
            if profile.size:
                scores[str(condition.get("name", ""))] = center_minus_flank_score(profile, xvals)
        if not scores:
            continue
        stronger = min(scores, key=scores.get)
        rows.append(
            {
                "name": motif.get("name", ""),
                "motif_id": motif.get("motif_id", ""),
                "n_sites_payload": int(motif.get("n_sites", 0)),
                "min_center_minus_flank": min(scores.values()),
                "max_center_minus_flank": max(scores.values()),
                "spread_center_minus_flank": max(scores.values()) - min(scores.values()),
                "deeper_condition": stronger,
            }
        )
    return pd.DataFrame(rows)


def low_complexity(consensus: str) -> bool:
    seq = re.sub(r"[^ACGT]", "", str(consensus).upper())
    if not seq:
        return True
    return max(seq.count(base) for base in "ACGT") / len(seq) >= 0.55 or len(set(seq)) <= 2


def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def max_window_identity(query: str, target: str) -> float:
    if not query or not target:
        return 0.0
    shorter, longer = (query, target) if len(query) <= len(target) else (target, query)
    if len(shorter) < 8:
        return 0.0
    best = 0.0
    for offset in range(len(longer) - len(shorter) + 1):
        window = longer[offset : offset + len(shorter)]
        matches = sum(a == b for a, b in zip(shorter, window))
        best = max(best, matches / len(shorter))
    return best


def tn5_bias_like(consensus: str) -> bool:
    seq = re.sub(r"[^ACGT]", "", str(consensus).upper())
    if len(seq) < 8:
        return False
    adapter = TN5_ADAPTER_CORE
    return max(max_window_identity(seq, adapter), max_window_identity(reverse_complement(seq), adapter)) >= 0.75


def quality_flags(row: pd.Series) -> str:
    flags = []
    if bool(row.get("low_complexity", False)):
        flags.append("low_complexity")
    if bool(row.get("very_broad", False)):
        flags.append("broad_site_set")
    if bool(row.get("tn5_bias_like", False)):
        flags.append("tn5_bias_like")
    if not bool(row.get("center_depleted", False)):
        flags.append("no_center_depletion")
    return ";".join(flags) if flags else "pass"


def choose_aggregate_examples(table: pd.DataFrame, n: int = N_AGGREGATE_EXAMPLES) -> list[str]:
    def rank_rows(rows: pd.DataFrame) -> pd.DataFrame:
        rows = rows.copy()
        empty_numeric = pd.Series(np.nan, index=rows.index)
        shape_values = rows["k562_rep1_footprint_shape_score"] if "k562_rep1_footprint_shape_score" in rows else empty_numeric
        site_values = rows["sites"] if "sites" in rows else empty_numeric
        rows["_shape"] = pd.to_numeric(shape_values, errors="coerce").fillna(-np.inf)
        rows["_sites"] = pd.to_numeric(site_values, errors="coerce").fillna(0)
        return rows.sort_values(
            ["_shape", "_sites"],
            ascending=[False, False],
        )

    base = table[
        table["source_condition"].astype(str).eq(ACTIVE_CONDITIONS[0])
        & ~table["confident_tomtom"]
        & ~table["low_complexity"]
        & ~table["very_broad"]
        & ~table["tn5_bias_like"]
        & table["k562_rep1_footprint_shape_score"].notna()
    ].copy()
    with_min_sites = base[pd.to_numeric(base["sites"], errors="coerce").fillna(0).ge(20)]
    if len(with_min_sites) >= n:
        base = with_min_sites
    if base.empty:
        base = table[
            table["source_condition"].astype(str).eq(ACTIVE_CONDITIONS[0])
            & ~table["low_complexity"]
            & ~table["very_broad"]
            & ~table["tn5_bias_like"]
            & table["k562_rep1_footprint_shape_score"].notna()
        ].copy()
    if base.empty:
        base = table[table["k562_rep1_footprint_shape_score"].notna()].copy()
    return rank_rows(base)["name"].astype(str).head(n).tolist()


def choose_condition_aggregate_examples(
    table: pd.DataFrame,
    condition: str,
    score_column: str,
    n: int = N_AGGREGATE_PER_CONDITION,
) -> list[str]:
    rows = table[
        table["source_condition"].astype(str).eq(condition)
        & ~table["confident_tomtom"]
        & ~table["low_complexity"]
        & ~table["very_broad"]
        & ~table["tn5_bias_like"]
        & table[score_column].notna()
    ].copy()
    with_min_sites = rows[pd.to_numeric(rows["sites"], errors="coerce").fillna(0).ge(20)]
    if len(with_min_sites) >= n:
        rows = with_min_sites
    if rows.empty:
        rows = table[
            table["source_condition"].astype(str).eq(condition)
            & ~table["low_complexity"]
            & ~table["very_broad"]
            & ~table["tn5_bias_like"]
            & table[score_column].notna()
        ].copy()
    if rows.empty:
        rows = table[table["source_condition"].astype(str).eq(condition)].copy()
    rows["_shape"] = pd.to_numeric(rows.get(score_column), errors="coerce").fillna(-np.inf)
    rows["_sites"] = pd.to_numeric(rows.get("sites"), errors="coerce").fillna(0)
    return rows.sort_values(["_shape", "_sites"], ascending=[False, False])["name"].astype(str).head(n).tolist()


def selected_motif_table(streme: pd.DataFrame, denovo_results: pd.DataFrame, scores: pd.DataFrame) -> pd.DataFrame:
    table = streme.merge(denovo_results, on="name", how="left", suffixes=("", "_result"))
    if not scores.empty:
        table = table.merge(scores, on="name", how="left")
    empty_numeric = pd.Series(np.nan, index=table.index)
    for column in [
        "k562_rep1_footprint_shape_score",
        "k562_rep1_center_minus_flank",
        "k562_rep1_k562_bound_sites",
    ]:
        if column not in table:
            table[column] = np.nan
    table["low_complexity"] = table["consensus"].map(low_complexity)
    table["very_broad"] = pd.to_numeric(table["total_tfbs"] if "total_tfbs" in table else empty_numeric, errors="coerce").fillna(0).gt(30000)
    table["center_depleted"] = pd.to_numeric(
        table["k562_rep1_center_minus_flank"],
        errors="coerce",
    ).lt(0)
    table["tn5_bias_like"] = table["consensus"].map(tn5_bias_like)
    highlighted_col = f"{ACTIVE_CONDITIONS[0]}_{ACTIVE_CONDITIONS[1]}_highlighted"
    if highlighted_col in table:
        table["highlighted_de_novo_only"] = table[highlighted_col].astype(str).eq("True")
    else:
        table["highlighted_de_novo_only"] = False
    table["quality_flags"] = table.apply(quality_flags, axis=1)
    selected_motifs = choose_aggregate_examples(table)
    panel_by_name = {name: chr(ord("C") + idx) for idx, name in enumerate(selected_motifs)}
    table["display_panel"] = table["name"].map(panel_by_name).fillna("")
    table["selected_for_aggregate"] = table["display_panel"].ne("")
    sort_key = {name: idx for idx, name in enumerate(selected_motifs)}
    table["_sort"] = table["name"].map(sort_key).fillna(1000)
    table = table.sort_values(
        ["selected_for_aggregate", "highlighted_de_novo_only", "confident_tomtom", "_sort", "tomtom_q_value"],
        ascending=[False, False, False, True, True],
    )
    return table.drop(columns=["_sort"])


def candidate_motif_table(streme: pd.DataFrame, denovo_results: pd.DataFrame) -> pd.DataFrame:
    table = streme.merge(denovo_results, on="name", how="left", suffixes=("", "_result"))
    empty_numeric = pd.Series(np.nan, index=table.index)
    table["low_complexity"] = table["consensus"].map(low_complexity)
    table["very_broad"] = pd.to_numeric(table["total_tfbs"] if "total_tfbs" in table else empty_numeric, errors="coerce").fillna(0).gt(30000)
    table["tn5_bias_like"] = table["consensus"].map(tn5_bias_like)
    return table


def read_site_windows(paths: list[Path], half_width: int = SITE_WINDOW_BP) -> dict[str, list[tuple[int, int]]]:
    windows: dict[str, list[tuple[int, int]]] = {}
    for path in paths:
        if not path.exists():
            continue
        with path.open() as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                chrom = fields[0]
                try:
                    start = int(fields[1])
                    end = int(fields[2])
                except ValueError:
                    continue
                center = (start + end) // 2
                windows.setdefault(chrom, []).append((max(0, center - half_width), center + half_width + 1))
    return merge_interval_dict(windows)


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [intervals[0]]
    for start, end in intervals[1:]:
        prev_start, prev_end = merged[-1]
        if start <= prev_end:
            merged[-1] = (prev_start, max(prev_end, end))
        else:
            merged.append((start, end))
    return merged


def merge_interval_dict(data: dict[str, list[tuple[int, int]]]) -> dict[str, list[tuple[int, int]]]:
    return {chrom: merge_intervals(intervals) for chrom, intervals in data.items() if intervals}


def interval_count(data: dict[str, list[tuple[int, int]]]) -> int:
    return sum(len(intervals) for intervals in data.values())


def interval_bp(data: dict[str, list[tuple[int, int]]]) -> int:
    return sum(end - start for intervals in data.values() for start, end in intervals)




def subtract_intervals(
    data_a: dict[str, list[tuple[int, int]]],
    data_b: dict[str, list[tuple[int, int]]],
) -> dict[str, list[tuple[int, int]]]:
    out: dict[str, list[tuple[int, int]]] = {}
    for chrom, intervals in data_a.items():
        blockers = data_b.get(chrom, [])
        if not blockers:
            out[chrom] = list(intervals)
            continue
        j = 0
        for start, end in intervals:
            pieces = [(start, end)]
            while j < len(blockers) and blockers[j][1] <= start:
                j += 1
            k = j
            while k < len(blockers) and blockers[k][0] < end:
                block_start, block_end = blockers[k]
                next_pieces = []
                for piece_start, piece_end in pieces:
                    if block_end <= piece_start or block_start >= piece_end:
                        next_pieces.append((piece_start, piece_end))
                        continue
                    if piece_start < block_start:
                        next_pieces.append((piece_start, block_start))
                    if block_end < piece_end:
                        next_pieces.append((block_end, piece_end))
                pieces = next_pieces
                if not pieces:
                    break
                k += 1
            out.setdefault(chrom, []).extend(pieces)
    return merge_interval_dict(out)


def intersection(data_a: dict[str, list[tuple[int, int]]], data_b: dict[str, list[tuple[int, int]]]) -> dict[str, list[tuple[int, int]]]:
    out: dict[str, list[tuple[int, int]]] = {}
    for chrom in set(data_a).intersection(data_b):
        a = data_a[chrom]
        b = data_b[chrom]
        i = j = 0
        while i < len(a) and j < len(b):
            start = max(a[i][0], b[j][0])
            end = min(a[i][1], b[j][1])
            if start < end:
                out.setdefault(chrom, []).append((start, end))
            if a[i][1] < b[j][1]:
                i += 1
            else:
                j += 1
    return merge_interval_dict(out)


def coverage_rows_for_label(label: str, de_paths: list[Path], jaspar_paths: list[Path], half_width: int) -> list[dict[str, object]]:
    de_windows = read_site_windows(de_paths, half_width=half_width)
    jaspar_windows = read_site_windows(jaspar_paths, half_width=half_width)
    shared = intersection(de_windows, jaspar_windows)
    jaspar_only = subtract_intervals(jaspar_windows, de_windows)
    denovo_only = subtract_intervals(de_windows, jaspar_windows)
    de_bp = interval_bp(de_windows)
    jaspar_bp = interval_bp(jaspar_windows)
    shared_bp = interval_bp(shared)
    jaspar_only_bp = interval_bp(jaspar_only)
    denovo_only_bp = interval_bp(denovo_only)
    rows = [
        {
            "site_set": label,
            "class": "shared",
            "merged_site_windows": interval_count(shared),
            "bp": shared_bp,
            "fraction_of_jaspar_bp": shared_bp / jaspar_bp if jaspar_bp else 0.0,
            "fraction_of_denovo_bp": shared_bp / de_bp if de_bp else 0.0,
            "window_half_width_bp": half_width,
        },
        {
            "site_set": label,
            "class": "JASPAR-only",
            "merged_site_windows": interval_count(jaspar_only),
            "bp": jaspar_only_bp,
            "fraction_of_jaspar_bp": jaspar_only_bp / jaspar_bp if jaspar_bp else 0.0,
            "fraction_of_denovo_bp": 0.0,
            "window_half_width_bp": half_width,
        },
        {
            "site_set": label,
            "class": "de novo-only",
            "merged_site_windows": interval_count(denovo_only),
            "bp": denovo_only_bp,
            "fraction_of_jaspar_bp": 0.0,
            "fraction_of_denovo_bp": denovo_only_bp / de_bp if de_bp else 0.0,
            "window_half_width_bp": half_width,
        },
    ]
    return rows


def collect_bound_beds(root: Path, set_name: str, condition: str | None = None, include_denovo: bool = True) -> list[Path]:
    base = root / "diff_footprints" / set_name
    pattern = f"*_{condition}_bound.bed" if condition else "*_bound.bed"
    paths = sorted(base.glob(f"*/beds/{pattern}"))
    if include_denovo:
        return paths
    return [p for p in paths if "_denovo_" not in p.name and "_denovo_" not in str(p.parent.parent.name)]


def corrected_signal_paths(validation_dir: Path) -> dict[str, list[tuple[str, Path]]]:
    paths_by_condition: dict[str, list[tuple[str, Path]]] = {condition: [] for condition in ACTIVE_CONDITIONS}
    atac_dir = validation_dir.parent / "atac_correct"
    for path in sorted(atac_dir.glob("*/*_corrected.bw")):
        sample = path.parent.name
        for condition in ACTIVE_CONDITIONS:
            if sample.startswith(f"{condition}_"):
                paths_by_condition.setdefault(condition, []).append((sample, path))
                break
    return paths_by_condition


def display_signal_path(validation_dir: Path) -> tuple[str, Path] | None:
    for sample, path in corrected_signal_paths(validation_dir).get(ACTIVE_CONDITIONS[0], []):
        if sample == DISPLAY_SAMPLE and path.name == DISPLAY_SAMPLE_BASENAME:
            return sample, path
    for sample, path in corrected_signal_paths(validation_dir).get(ACTIVE_CONDITIONS[0], []):
        if sample == DISPLAY_SAMPLE:
            return sample, path
    return None


def rep1_signal_path(validation_dir: Path, condition: str) -> tuple[str, Path] | None:
    preferred_sample = f"{condition}_rep1"
    for sample, path in corrected_signal_paths(validation_dir).get(condition, []):
        if sample == preferred_sample and path.name.endswith("_corrected.bw"):
            return sample, path
    for sample, path in corrected_signal_paths(validation_dir).get(condition, []):
        if sample.endswith("_rep1") and path.name.endswith("_corrected.bw"):
            return sample, path
    paths = corrected_signal_paths(validation_dir).get(condition, [])
    return paths[0] if paths else None


def motif_bed_dir(validation_dir: Path, motif_name: str) -> Path | None:
    base = validation_dir / "diff_footprints" / "denovo_only"
    if not base.exists():
        return None
    matches = sorted(path for path in base.iterdir() if path.is_dir() and path.name.startswith(f"{motif_name}_"))
    return matches[0] / "beds" if matches else None


def read_bed_centers(
    path: Path,
    flank: int,
    chrom_sizes: dict[str, int],
    max_sites: int | None = 12000,
) -> list[tuple[str, int, int]]:
    sites: list[tuple[str, int, int]] = []
    if not path.exists():
        return sites
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            chrom = fields[0]
            if chrom not in chrom_sizes:
                continue
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                continue
            center = (start + end) // 2
            window_start = center - flank
            window_end = center + flank + 1
            if window_start < 0 or window_end > chrom_sizes[chrom]:
                continue
            sites.append((chrom, window_start, window_end))
            if max_sites is not None and len(sites) >= max_sites:
                break
    return sites


def mean_bigwig_profile(path: Path, sites: list[tuple[str, int, int]]) -> list[float]:
    if not sites:
        return []
    profiles = []
    with pyBigWig.open(str(path)) as bw:
        for chrom, start, end in sites:
            values = np.asarray(bw.values(chrom, start, end, numpy=True), dtype=float)
            if values.size == end - start:
                profiles.append(values)
    if not profiles:
        return []
    return np.nanmean(np.vstack(profiles), axis=0).tolist()


def compute_selected_aggregates(
    validation_dir: Path,
    selected: pd.DataFrame,
    xvals: np.ndarray,
    max_motifs: int = 16,
) -> dict[str, dict]:
    if xvals.size == 0:
        xvals = np.arange(-100, 101, dtype=float)
    flank = int(max(abs(float(xvals.min())), abs(float(xvals.max()))))
    signal_paths = corrected_signal_paths(validation_dir)
    all_signal_paths = [path for paths in signal_paths.values() for _, path in paths]
    if not all_signal_paths:
        return {}
    with pyBigWig.open(str(all_signal_paths[0])) as bw:
        chrom_sizes = bw.chroms()

    motifs: dict[str, dict] = {}
    rows = selected[selected["selected_for_aggregate"]].copy()
    if len(rows) < max_motifs:
        extra = selected[
            ~selected["selected_for_aggregate"]
            & ~selected["low_complexity"]
            & ~selected["very_broad"]
            & ~selected["tn5_bias_like"]
        ].head(max_motifs - len(rows))
        rows = pd.concat([rows, extra], ignore_index=True)
    for _, row in rows.iterrows():
        motif_name = str(row["name"])
        beds_dir = motif_bed_dir(validation_dir, motif_name)
        if beds_dir is None:
            continue
        motif = {
            "name": motif_name,
            "motif_id": str(row.get("de_novo_motif", "")),
            "n_sites": 0,
            "conditions": [],
        }
        for condition in ACTIVE_CONDITIONS:
            bed = beds_dir / f"{beds_dir.parent.name}_{condition}_bound.bed"
            sites = read_bed_centers(bed, flank=flank, chrom_sizes=chrom_sizes)
            motif["n_sites"] += len(sites)
            samples = []
            for sample, signal_path in signal_paths.get(condition, []):
                profile = mean_bigwig_profile(signal_path, sites)
                if profile:
                    samples.append({"name": sample, "profile": profile})
            if samples:
                condition_profiles = np.vstack([np.asarray(sample["profile"], dtype=float) for sample in samples])
                motif["conditions"].append(
                    {
                        "name": condition,
                        "profile": np.nanmean(condition_profiles, axis=0).tolist(),
                        "samples": samples,
                    }
                )
        if motif["conditions"]:
            motifs[motif_name] = motif
    return motifs


def compute_condition_rep1_aggregates(
    validation_dir: Path,
    candidates: pd.DataFrame,
    xvals: np.ndarray,
    condition: str,
    max_sites_per_motif: int | None = MAX_RANKING_SITES_PER_MOTIF,
) -> tuple[dict[str, dict], pd.DataFrame]:
    if xvals.size == 0:
        xvals = np.arange(-100, 101, dtype=float)
    flank = int(max(abs(float(xvals.min())), abs(float(xvals.max()))))
    display_signal = rep1_signal_path(validation_dir, condition)
    if display_signal is None:
        return {}, pd.DataFrame()
    sample_name, signal_path = display_signal
    with pyBigWig.open(str(signal_path)) as bw:
        chrom_sizes = bw.chroms()

    candidate_rows = candidates[
        candidates["source_condition"].astype(str).eq(condition)
        & ~candidates["confident_tomtom"]
        & ~candidates["low_complexity"]
        & ~candidates["very_broad"]
        & ~candidates["tn5_bias_like"]
    ].copy()
    with_min_sites = candidate_rows[pd.to_numeric(candidate_rows["sites"], errors="coerce").fillna(0).ge(20)]
    if len(with_min_sites) >= N_AGGREGATE_PER_CONDITION:
        candidate_rows = with_min_sites
    if candidate_rows.empty:
        candidate_rows = candidates[candidates["source_condition"].astype(str).eq(condition)].copy()

    motifs: dict[str, dict] = {}
    score_rows: list[dict[str, object]] = []
    for _, row in candidate_rows.iterrows():
        motif_name = str(row["name"])
        beds_dir = motif_bed_dir(validation_dir, motif_name)
        if beds_dir is None:
            continue
        bed = beds_dir / f"{beds_dir.parent.name}_{condition}_bound.bed"
        sites = read_bed_centers(
            bed,
            flank=flank,
            chrom_sizes=chrom_sizes,
            max_sites=max_sites_per_motif,
        )
        profile = mean_bigwig_profile(signal_path, sites)
        if not profile:
            continue
        profile_array = np.asarray(profile, dtype=float)
        shape_score = footprint_shape_score(profile_array, xvals)
        center_minus_flank = center_minus_flank_score(profile_array, xvals)
        motifs[motif_name] = {
            "name": motif_name,
            "motif_id": str(row.get("de_novo_motif", "")),
            "n_sites": len(sites),
            "conditions": [
                {
                    "name": condition,
                    "profile": profile,
                    "samples": [{"name": sample_name, "profile": profile}],
                }
            ],
        }
        prefix = condition.lower()
        score_rows.append(
            {
                "name": motif_name,
                f"{prefix}_rep1_footprint_shape_score": shape_score,
                f"{prefix}_rep1_center_minus_flank": center_minus_flank,
                f"{prefix}_rep1_bound_sites": len(sites),
            }
        )
    return motifs, pd.DataFrame(score_rows)


def compute_k562_rep1_aggregates(
    validation_dir: Path,
    candidates: pd.DataFrame,
    xvals: np.ndarray,
    max_sites_per_motif: int | None = MAX_RANKING_SITES_PER_MOTIF,
) -> tuple[dict[str, dict], pd.DataFrame]:
    motifs, scores = compute_condition_rep1_aggregates(
        validation_dir,
        candidates,
        xvals,
        ACTIVE_CONDITIONS[0],
        max_sites_per_motif=max_sites_per_motif,
    )
    if not scores.empty:
        prefix = ACTIVE_CONDITIONS[0].lower()
        scores = scores.rename(
            columns={
                f"{prefix}_rep1_footprint_shape_score": "k562_rep1_footprint_shape_score",
                f"{prefix}_rep1_center_minus_flank": "k562_rep1_center_minus_flank",
                f"{prefix}_rep1_bound_sites": "k562_rep1_k562_bound_sites",
            }
        )
    return motifs, scores


def compute_global_coverage(validation_dir: Path, half_width: int = SITE_WINDOW_BP) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    denovo_set = "denovo_only"
    jaspar_set = "jaspar2026_plus_denovo"
    rows.extend(
        coverage_rows_for_label(
            "all_bound",
            collect_bound_beds(validation_dir, denovo_set, include_denovo=True),
            collect_bound_beds(validation_dir, jaspar_set, include_denovo=False),
            half_width,
        )
    )
    for condition in ACTIVE_CONDITIONS:
        rows.extend(
            coverage_rows_for_label(
                f"{condition}_bound",
                collect_bound_beds(validation_dir, denovo_set, condition=condition, include_denovo=True),
                collect_bound_beds(validation_dir, jaspar_set, condition=condition, include_denovo=False),
                half_width,
            )
        )
    return pd.DataFrame(rows)


def plot_coverage_bar(ax, coverage: pd.DataFrame, motif_sets: pd.DataFrame) -> None:
    ax.set_title(
        "A. De novo-only site coverage of JASPAR-supported footprints",
        loc="left",
        pad=9,
        fontsize=9.6,
        fontweight="bold",
    )
    plot_sets = [f"{ACTIVE_CONDITIONS[0]}_bound", f"{ACTIVE_CONDITIONS[1]}_bound"]
    labels = list(ACTIVE_CONDITIONS)
    colors = {"shared": "#4b9b62", "de novo-only": "#4f8dcf"}
    classes = ["shared", "de novo-only"]
    y = np.arange(len(plot_sets))
    left = np.zeros(len(plot_sets), dtype=float)
    for cls in classes:
        values = []
        for site_set in plot_sets:
            match = coverage[(coverage["site_set"].eq(site_set)) & (coverage["class"].eq(cls))]
            values.append(float(match["bp"].iloc[0]) if not match.empty else 0.0)
        ax.barh(y, values, left=left, color=colors[cls], edgecolor="white", linewidth=0.5, label=cls, height=0.58)
        left += np.asarray(values, dtype=float)
    xmax = max(float(np.nanmax(left)) * 1.06, 1.0)
    tick_step = 10_000_000
    xticks = np.arange(0, np.ceil(xmax / tick_step) * tick_step + 1, tick_step)
    ax.set_xlim(0, xticks[-1] if xticks.size else xmax)
    ax.set_xticks(xticks)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel(f"Merged site-window bp (+/-{SITE_WINDOW_BP} bp around motif centers)")
    ax.ticklabel_format(axis="x", style="plain", useOffset=False)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda value, _pos: f"{int(value):,}"))
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color="0.9", linewidth=0.55)
    ax.legend(loc="upper right", bbox_to_anchor=(1.0, 1.13), frameon=False, fontsize=7.1, ncol=2)
    bold_all_text(ax)


def format_tomtom(row: pd.Series) -> str:
    q = row.get("tomtom_q_value")
    label = str(row.get("tomtom_label", ""))
    if pd.notna(q) and bool(row.get("confident_tomtom", False)):
        return f"{label.split('(', 1)[0].strip()}-like (q={float(q):.1e})"
    if label and label != "nan":
        return label
    return "no confident match"


def motif_short_label(name: object) -> str:
    match = re.search(r"_denovo_(\d+)", str(name))
    if match:
        return f"dn{match.group(1)}"
    return str(name).replace("_", " ")


def draw_sequence_logo(
    ax,
    matrix: np.ndarray | None,
    consensus: str,
    bounds: tuple[float, float, float, float],
    *,
    border: bool = False,
) -> None:
    logo_ax = ax.inset_axes(bounds)
    if border:
        logo_ax.set_facecolor("white")
        logo_ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
        for spine in logo_ax.spines.values():
            spine.set_visible(True)
            spine.set_color("black")
            spine.set_linewidth(1.15)
    else:
        logo_ax.set_axis_off()
    if matrix is None or matrix.size == 0:
        for idx, base in enumerate(str(consensus).upper()):
            logo_ax.text(
                idx + 0.5,
                0.5,
                base,
                ha="center",
                va="center",
                fontsize=7.5,
                fontweight="bold",
                color=LOGO_COLORS.get(base, "0.35"),
            )
        width = max(len(str(consensus)), 1)
        if border:
            logo_ax.set_xlim(-0.45, width + 0.45)
            logo_ax.set_ylim(-0.12, 1.12)
        else:
            logo_ax.set_xlim(0, width)
            logo_ax.set_ylim(0, 1)
        return

    bases = ["A", "C", "G", "T"]
    width = matrix.shape[0]
    if border:
        logo_ax.set_xlim(-0.45, width + 0.45)
        logo_ax.set_ylim(-0.15, 2.15)
    else:
        logo_ax.set_xlim(0, width)
        logo_ax.set_ylim(0, 2.0)
    for pos, probs in enumerate(matrix):
        probs = np.asarray(probs, dtype=float)
        probs = probs / probs.sum() if probs.sum() > 0 else np.full(4, 0.25)
        entropy = -float(np.sum([p * np.log2(p) for p in probs if p > 0]))
        information = max(0.0, 2.0 - entropy)
        y_offset = 0.0
        for base_idx in np.argsort(probs):
            base = bases[int(base_idx)]
            height = float(probs[base_idx]) * information
            if height <= 0.012:
                continue
            text_path = TextPath((0, 0), base, size=1.0, prop=LOGO_FONT)
            bbox = text_path.get_extents()
            scale_x = 0.82 / max(bbox.width, 1e-6)
            scale_y = height / max(bbox.height, 1e-6)
            transform = (
                Affine2D()
                .scale(scale_x, scale_y)
                .translate(pos + 0.09 - bbox.x0 * scale_x, y_offset - bbox.y0 * scale_y)
                + logo_ax.transData
            )
            logo_ax.add_patch(PathPatch(text_path, transform=transform, color=LOGO_COLORS[base], linewidth=0))
            y_offset += height


def plot_discovery_table(ax, selected: pd.DataFrame, logo_matrices: dict[str, np.ndarray]) -> pd.DataFrame:
    ax.axis("off")
    ax.set_title(
        "B. K562-derived de novo motifs used for K562 aggregate examples",
        loc="left",
        pad=5,
        fontsize=9.4,
        fontweight="bold",
    )
    display = selected[selected["selected_for_aggregate"]].copy()
    if display.empty:
        display = selected.head(N_AGGREGATE_EXAMPLES).copy()
    display = display.head(N_AGGREGATE_EXAMPLES)
    blocks = [(0.00, 0.485), (0.515, 1.00)]
    col_fracs = [0.00, 0.18, 0.39, 0.87, 1.00]
    headers = ["motif", "K562\nfootprints", "sequence logo", "score"]
    top = 0.84
    header_h = 0.16
    row_h = 0.165
    for block_left, block_right in blocks:
        block_width = block_right - block_left
        xs = [block_left + frac * block_width for frac in col_fracs]
        ax.add_patch(
            plt.Rectangle(
                (block_left, top - header_h),
                block_width,
                header_h,
                transform=ax.transAxes,
                facecolor="#f4f7fb",
                edgecolor="0.72",
                linewidth=0.65,
            )
        )
        for idx, h in enumerate(headers):
            ax.text(
                xs[idx] + 0.008 * block_width,
                top - 0.035,
                h,
                transform=ax.transAxes,
                fontsize=6.9,
                fontweight="bold",
                va="top",
            )
        for x in xs:
            ax.plot([x, x], [top - header_h - 4 * row_h, top], color="0.78", linewidth=0.55, transform=ax.transAxes)
        for row_idx in range(5):
            y = top - header_h - row_idx * row_h
            ax.plot([block_left, block_right], [y, y], color="0.78", linewidth=0.55, transform=ax.transAxes)
        ax.plot([block_left, block_right], [top, top], color="0.62", linewidth=0.8, transform=ax.transAxes)
    summary_rows = []
    for idx, (_, row) in enumerate(display.iterrows()):
        block_left, block_right = blocks[idx // 4]
        block_width = block_right - block_left
        xs = [block_left + frac * block_width for frac in col_fracs]
        row_top = top - header_h - (idx % 4) * row_h
        y = row_top - 0.035
        footprint_n = pd.to_numeric(pd.Series([row.get("K562_bound", np.nan)]), errors="coerce").iloc[0]
        score = pd.to_numeric(
            pd.Series([row.get("k562_rep1_footprint_shape_score", np.nan)]),
            errors="coerce",
        ).iloc[0]
        ax.text(xs[0] + 0.008 * block_width, y, motif_short_label(row["name"]), transform=ax.transAxes, fontsize=7.0, va="top", fontweight="bold")
        ax.text(xs[1] + 0.008 * block_width, y, f"{int(footprint_n)}" if pd.notna(footprint_n) else "NA", transform=ax.transAxes, fontsize=6.8, va="top")
        ax.text(xs[3] + 0.008 * block_width, y, f"{score:.2f}" if pd.notna(score) else "NA", transform=ax.transAxes, fontsize=6.8, va="top")
        draw_sequence_logo(
            ax,
            logo_matrices.get(str(row["name"])),
            str(row["consensus"]),
            (xs[2] + 0.010 * block_width, row_top - row_h + 0.023, (xs[3] - xs[2]) - 0.020 * block_width, row_h - 0.046),
        )
        summary_rows.append(row.to_dict())
    bold_all_text(ax)
    return pd.DataFrame(summary_rows)


def aggregate_title(motif: dict) -> tuple[str, str]:
    prefix = str(motif.get("name", ""))
    if prefix.startswith(ACTIVE_CONDITIONS[0]):
        source = f"{ACTIVE_CONDITIONS[0]} candidates"
    elif prefix.startswith(ACTIVE_CONDITIONS[1]):
        source = f"{ACTIVE_CONDITIONS[1]} candidates"
    else:
        source = "de novo candidates"
    motif_name = str(motif.get("name", "de novo motif")).replace("_", " ")
    motif_name = re.sub(r" K562 denovo \d+ \d+-.*$", "", motif_name)
    return motif_name, source


def plot_aggregate_panel(
    ax,
    motif: dict,
    selected_row: pd.Series | None,
    xvals: np.ndarray,
    panel_label: str,
    summary_rows: list[dict[str, object]],
    *,
    show_xlabel: bool,
    show_ylabel: bool,
):
    for condition in motif.get("conditions", []):
        condition_name = str(condition["name"])
        if condition_name != ACTIVE_CONDITIONS[0]:
            continue
        for sample in condition.get("samples", []):
            sample_name = str(sample.get("name", "sample"))
            if sample_name != DISPLAY_SAMPLE:
                continue
            profile = np.asarray(sample["profile"], dtype=float)
            sample_score = center_minus_flank_score(profile, xvals)
            shape_score = footprint_shape_score(profile, xvals)
            display_mask = (xvals >= -60) & (xvals <= 60)
            ax.plot(
                xvals[display_mask],
                profile[display_mask],
                color=CONDITION_COLORS.get(ACTIVE_CONDITIONS[0], "#1f77b4"),
                linestyle="solid",
                linewidth=0.85,
                alpha=0.98,
                zorder=2,
            )
            summary_rows.append(
                {
                    "panel": panel_label,
                    "motif": motif.get("name", ""),
                    "motif_id": motif.get("motif_id", ""),
                    "streme_sites": int(float(selected_row.get("sites", 0))) if selected_row is not None else np.nan,
                    "k562_bound_footprints": int(float(selected_row.get("K562_bound", 0))) if selected_row is not None else np.nan,
                    "profile_sites": int(motif.get("n_sites", 0)),
                    "sample": sample_name,
                    "condition": condition_name,
                    "k562_rep1_footprint_shape_score": shape_score,
                    "k562_rep1_center_minus_flank": sample_score,
                    "center_minus_flank": sample_score,
                    "selection_role": "K562_rep1 aggregate example",
                }
            )
    ax.axvline(0, color="0.35", linewidth=0.75)
    ax.axhline(0, color="0.78", linewidth=0.55, zorder=0)
    ax.set_xlim(-60, 60)
    footprint_n = int(float(selected_row.get("K562_bound", motif.get("n_sites", 0)))) if selected_row is not None else int(motif.get("n_sites", 0))
    shape_score = pd.to_numeric(
        pd.Series([selected_row.get("k562_rep1_footprint_shape_score", np.nan) if selected_row is not None else np.nan]),
        errors="coerce",
    ).iloc[0]
    ax.set_title(
        f"{panel_label}. K562 {motif_short_label(motif.get('name', ''))}, N={footprint_n:,}\nscore={shape_score:.2f}",
        fontsize=7.1,
        fontweight="bold",
    )
    ax.set_xlabel("Distance from motif center (bp)" if show_xlabel else "")
    ax.set_ylabel("Normalized cut-site signal" if show_ylabel else "")
    ax.set_box_aspect(1)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="0.9", linewidth=0.5)
    bold_all_text(ax)


def plot_condition_aggregate_mini(
    ax,
    motif: dict,
    selected_row: pd.Series,
    xvals: np.ndarray,
    condition: str,
    logo_matrices: dict[str, np.ndarray],
    summary_rows: list[dict[str, object]],
) -> None:
    motif_name = str(motif.get("name", ""))
    matrix = logo_matrices.get(motif_name)
    width = int(matrix.shape[0]) if matrix is not None and matrix.size else len(str(selected_row.get("consensus", "")))
    profile = None
    sample_name = ""
    for condition_row in motif.get("conditions", []):
        if str(condition_row.get("name")) != condition:
            continue
        for sample in condition_row.get("samples", []):
            sample_name = str(sample.get("name", ""))
            profile = np.asarray(sample.get("profile", []), dtype=float)
            break
    if profile is None or profile.size == 0:
        ax.axis("off")
        return

    display_mask = (xvals >= -60) & (xvals <= 60)
    ax.plot(
        xvals[display_mask],
        profile[display_mask],
        color=CONDITION_COLORS.get(condition, "#1f77b4"),
        linestyle="solid",
        linewidth=1.0,
        alpha=0.98,
        zorder=2,
    )
    for boundary in (-width / 2.0, width / 2.0):
        ax.axvline(boundary, color="0.55", linestyle=(0, (1.2, 1.8)), linewidth=0.85, zorder=1)
    ax.axhline(0, color="0.78", linewidth=0.55, zorder=0)
    ax.set_xlim(-60, 60)
    ax.set_box_aspect(1)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="0.9", linewidth=0.5)

    bound_n = pd.to_numeric(pd.Series([selected_row.get(f"{condition}_bound", motif.get("n_sites", 0))]), errors="coerce").iloc[0]
    draw_sequence_logo(
        ax,
        matrix,
        str(selected_row.get("consensus", "")),
        (0.0, 1.06, 1.0, 0.23),
        border=True,
    )
    ax.text(
        0.04,
        0.06,
        f"{int(bound_n)}" if pd.notna(bound_n) else "NA",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=6.5,
        fontweight="normal",
        bbox={
            "boxstyle": "round,pad=0.16",
            "facecolor": "white",
            "edgecolor": "none",
            "alpha": 0.78,
        },
    )

    shape_score = footprint_shape_score(profile, xvals)
    center_minus_flank = center_minus_flank_score(profile, xvals)
    summary_rows.append(
        {
            "panel_group": condition,
            "motif": motif_name,
            "motif_id": motif.get("motif_id", ""),
            "streme_sites": int(float(selected_row.get("sites", 0))),
            "condition_bound_footprints": int(bound_n) if pd.notna(bound_n) else np.nan,
            "profile_sites": int(motif.get("n_sites", 0)),
            "motif_width": width,
            "sample": sample_name,
            "condition": condition,
            "rep1_footprint_shape_score": shape_score,
            "rep1_center_minus_flank": center_minus_flank,
            "selection_role": f"{condition}_rep1 aggregate example",
        }
    )
    bold_all_text(ax)


def plot_condition_aggregate_group(
    fig,
    subspec,
    panel_label: str,
    condition: str,
    motif_names: list[str],
    motifs_by_name: dict[str, dict],
    selected_lookup: pd.DataFrame,
    xvals: np.ndarray,
    logo_matrices: dict[str, np.ndarray],
    summary_rows: list[dict[str, object]],
) -> None:
    nested = subspec.subgridspec(3, 2, height_ratios=[0.13, 1.0, 1.0], hspace=0.95, wspace=0.42)
    title_ax = fig.add_subplot(nested[0, :])
    title_ax.axis("off")
    title_ax.text(0.0, 0.5, f"{panel_label}. {condition}", ha="left", va="center", fontsize=9.2, fontweight="bold")
    for idx, motif_name in enumerate(motif_names[:N_AGGREGATE_PER_CONDITION]):
        ax = fig.add_subplot(nested[1 + idx // 2, idx % 2])
        row = selected_lookup.loc[motif_name]
        plot_condition_aggregate_mini(
            ax,
            motifs_by_name[motif_name],
            row,
            xvals,
            condition,
            logo_matrices,
            summary_rows,
        )


def plot_validation(validation_dir: Path, jaspar: Path, out_prefix: Path, conditions: tuple[str, str] = DEFAULT_CONDITIONS) -> None:
    global ACTIVE_CONDITIONS
    ACTIVE_CONDITIONS = conditions
    motif_sets = pd.read_csv(validation_dir / "motifs" / "motif_set_summary.tsv", sep="\t")
    id_to_name = jaspar_name_map(jaspar)
    streme = load_streme_motifs(validation_dir, id_to_name)
    logo_matrices = load_streme_probability_matrices(validation_dir)
    denovo_results = pd.read_csv(validation_dir / "diff_footprints" / "denovo_only" / "diff_footprints_results.txt", sep="\t")
    payload = load_report_payload(
        validation_dir
        / "diff_footprints"
        / "denovo_only"
        / f"diff_footprints_{ACTIVE_CONDITIONS[0]}_{ACTIVE_CONDITIONS[1]}.html"
    )
    xvals = np.asarray(payload.get("aggregate", {}).get("x", []), dtype=float)
    if xvals.size == 0:
        xvals = np.arange(-100, 101, dtype=float)
    motifs_by_name = motif_lookup(payload)
    candidate_table = candidate_motif_table(streme, denovo_results)
    score_frames = [aggregate_scores(payload, xvals)]
    computed_motifs, k562_scores = compute_k562_rep1_aggregates(validation_dir, candidate_table, xvals)
    hepg2_motifs, hepg2_scores = compute_condition_rep1_aggregates(validation_dir, candidate_table, xvals, ACTIVE_CONDITIONS[1])
    motifs_by_name.update(computed_motifs)
    motifs_by_name.update(hepg2_motifs)
    for score_frame in (k562_scores, hepg2_scores):
        if not score_frame.empty:
            score_frames.append(score_frame)
    scores = pd.DataFrame()
    for score_frame in score_frames:
        if score_frame.empty:
            continue
        scores = score_frame if scores.empty else scores.merge(score_frame, on="name", how="outer")
    selected = selected_motif_table(streme, denovo_results, scores)

    k562_score_col = f"{ACTIVE_CONDITIONS[0].lower()}_rep1_footprint_shape_score"
    if k562_score_col not in selected and "k562_rep1_footprint_shape_score" in selected:
        k562_score_col = "k562_rep1_footprint_shape_score"
    hepg2_score_col = f"{ACTIVE_CONDITIONS[1].lower()}_rep1_footprint_shape_score"
    k562_names = choose_condition_aggregate_examples(selected, ACTIVE_CONDITIONS[0], k562_score_col)
    hepg2_names = choose_condition_aggregate_examples(selected, ACTIVE_CONDITIONS[1], hepg2_score_col)

    full_k562_motifs, full_k562_scores = compute_k562_rep1_aggregates(
        validation_dir,
        selected[selected["name"].astype(str).isin(k562_names)].copy(),
        xvals,
        max_sites_per_motif=None,
    )
    full_hepg2_motifs, full_hepg2_scores = compute_condition_rep1_aggregates(
        validation_dir,
        selected[selected["name"].astype(str).isin(hepg2_names)].copy(),
        xvals,
        ACTIVE_CONDITIONS[1],
        max_sites_per_motif=None,
    )
    motifs_by_name.update(full_k562_motifs)
    motifs_by_name.update(full_hepg2_motifs)

    full_score_frames = [frame for frame in (full_k562_scores, full_hepg2_scores) if not frame.empty]
    if full_score_frames:
        full_scores = full_score_frames[0]
        for frame in full_score_frames[1:]:
            full_scores = full_scores.merge(frame, on="name", how="outer")
        selected = selected.set_index("name")
        full_scores = full_scores.set_index("name")
        for column in full_scores.columns:
            selected.loc[full_scores.index, column] = full_scores[column]
        selected = selected.reset_index()

    selected_lookup = selected.set_index("name") if not selected.empty else pd.DataFrame()
    aggregate_group_rows: list[dict[str, object]] = []
    for panel_label, condition, names in [
        ("C", ACTIVE_CONDITIONS[0], k562_names),
        ("D", ACTIVE_CONDITIONS[1], hepg2_names),
    ]:
        for order, name in enumerate(names, start=1):
            if name not in selected_lookup.index:
                continue
            row = selected_lookup.loc[name].to_dict()
            matrix = logo_matrices.get(str(name))
            row.update(
                {
                    "name": name,
                    "panel_group": panel_label,
                    "aggregate_condition": condition,
                    "aggregate_order": order,
                    "motif_width": int(matrix.shape[0]) if matrix is not None and matrix.size else len(str(row.get("consensus", ""))),
                    "condition_bound_footprints": row.get(f"{condition}_bound", np.nan),
                }
            )
            aggregate_group_rows.append(row)
    coverage = compute_global_coverage(validation_dir, half_width=SITE_WINDOW_BP)

    apply_style(base_size=8.8)
    fig = plt.figure(figsize=(6.45, 7.05))
    gs = fig.add_gridspec(4, 4, height_ratios=[0.78, 1.05, 1.12, 1.12], hspace=0.66, wspace=0.46)

    ax_coverage = fig.add_subplot(gs[0, :])
    plot_coverage_bar(ax_coverage, coverage, motif_sets)

    ax_table = fig.add_subplot(gs[1, :])
    discovery_rows = plot_discovery_table(ax_table, selected, logo_matrices)

    summary_rows: list[dict[str, object]] = []
    plot_condition_aggregate_group(
        fig,
        gs[2:, :2],
        "C",
        ACTIVE_CONDITIONS[0],
        [name for name in k562_names if name in motifs_by_name and name in selected_lookup.index],
        motifs_by_name,
        selected_lookup,
        xvals,
        logo_matrices,
        summary_rows,
    )
    plot_condition_aggregate_group(
        fig,
        gs[2:, 2:],
        "D",
        ACTIVE_CONDITIONS[1],
        [name for name in hepg2_names if name in motifs_by_name and name in selected_lookup.index],
        motifs_by_name,
        selected_lookup,
        xvals,
        logo_matrices,
        summary_rows,
    )

    fig.text(0.54, 0.035, "Distance from motif center (bp)", ha="center", va="center", fontsize=8.8, fontweight="bold")
    fig.text(0.035, 0.31, "Normalized cut-site signal", ha="center", va="center", rotation="vertical", fontsize=8.8, fontweight="bold")
    fig.subplots_adjust(left=0.09, right=0.99, top=0.965, bottom=0.08, hspace=0.66, wspace=0.46)

    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    svg_path = out_prefix.with_suffix(".svg")
    fig.savefig(svg_path, bbox_inches="tight")
    svg_path.write_text("\n".join(line.rstrip() for line in svg_path.read_text().splitlines()) + "\n")
    print(f"Wrote {svg_path}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validation-dir", type=Path, default=DEFAULT_VALIDATION_DIR)
    parser.add_argument("--jaspar", type=Path, default=DEFAULT_JASPAR)
    parser.add_argument("--out-prefix", type=Path, default=REPO_ROOT / "manuscript/figures/denovo_motif_validation")
    parser.add_argument("--conditions", nargs=2, default=DEFAULT_CONDITIONS, metavar=("COND1", "COND2"))
    args = parser.parse_args(argv)
    plot_validation(args.validation_dir, args.jaspar, args.out_prefix, tuple(args.conditions))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
