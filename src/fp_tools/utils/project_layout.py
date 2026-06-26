"""Shared project-layout helpers for command-oriented fp-tools workflows."""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path


MITO_CHROMS = ("chrM", "chrMT", "M", "MT", "Mito")


@dataclass(frozen=True)
class SampleRecord:
    sample: str
    condition: str
    bam: str = ""
    peaks: str = ""


@dataclass(frozen=True)
class ComparisonRecord:
    comparison: str
    cond1: str
    cond2: str
    cond1_samples: tuple[str, ...] = ()
    cond2_samples: tuple[str, ...] = ()


def is_project_layout(value: str | None) -> bool:
    return str(value or "").strip().lower() == "project"


def project_root(outdir: str | Path | None) -> Path:
    if not outdir:
        raise ValueError("--layout project requires --outdir")
    return Path(outdir).expanduser().resolve()


def samples_root(project: str | Path) -> Path:
    return Path(project) / "samples"


def peaks_dir(project: str | Path) -> Path:
    return Path(project) / "peaks"


def comparisons_dir(project: str | Path) -> Path:
    return Path(project) / "comparisons"


def reports_dir(project: str | Path) -> Path:
    return Path(project) / "reports"


def logs_dir(project: str | Path) -> Path:
    return Path(project) / "logs"


def merged_peaks_path(project: str | Path) -> Path:
    return peaks_dir(project) / "merged_peaks.bed"


def analysis_peaks_path(project: str | Path) -> Path:
    return peaks_dir(project) / "merged_peaks.analysis.bed"


def corrected_bigwig_path(project: str | Path, sample: str) -> Path:
    return samples_root(project) / sample / "atac_correct" / f"{sample}_corrected.bw"


def normalized_bigwig_path(project: str | Path, sample: str) -> Path:
    return samples_root(project) / sample / "normalize" / f"{sample}_corrected_q95_scaled.bw"


def footprint_bigwig_path(project: str | Path, sample: str) -> Path:
    return samples_root(project) / sample / "footprints" / f"{sample}_footprints.bw"


def candidate_footprints_path(project: str | Path, sample: str) -> Path:
    return samples_root(project) / sample / "footprints" / f"{sample}_candidate_footprints.bed"


def match_motifs_dir(project: str | Path, sample: str) -> Path:
    return samples_root(project) / sample / "match_motifs"


def sample_dir(project: str | Path, sample: str) -> Path:
    return samples_root(project) / sample


def comparison_dir(project: str | Path, comparison: str) -> Path:
    return comparisons_dir(project) / safe_name(comparison)


def normalize_qc_dir(project: str | Path) -> Path:
    return logs_dir(project) / "normalize_q95"


def review_output_path(project: str | Path, name: str = "review_multi_comparisons.html") -> Path:
    return reports_dir(project) / name


def safe_name(value: str) -> str:
    text = str(value).strip()
    text = re.sub(r"\s+", "_", text)
    text = re.sub(r"[^A-Za-z0-9_.:+-]+", "_", text)
    text = text.strip("_")
    return text or "comparison"


def _first(row: dict[str, str], names: tuple[str, ...], default: str = "") -> str:
    lower_map = {key.lower(): key for key in row}
    for name in names:
        key = lower_map.get(name.lower())
        if key is not None:
            return str(row.get(key, "")).strip()
    return default


def _read_tsv(path: str | Path) -> list[dict[str, str]]:
    with Path(path).expanduser().open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} is missing a header row")
        return [dict(row) for row in reader]


def read_sample_table(path: str | Path) -> list[SampleRecord]:
    rows = []
    for line_no, row in enumerate(_read_tsv(path), start=2):
        sample = _first(row, ("sample", "sample_id", "id"))
        condition = _first(row, ("condition", "condition_label", "cond", "group"))
        bam = _first(row, ("bam", "bam_path", "input_bam"))
        peaks = _first(row, ("peaks", "peak_bed", "peak", "bed"))
        if not sample or not condition:
            raise ValueError(f"{path}:{line_no} requires sample and condition columns")
        rows.append(SampleRecord(sample=sample, condition=condition, bam=bam, peaks=peaks))
    if not rows:
        raise ValueError(f"{path} does not contain any sample rows")
    seen = set()
    duplicates = []
    for row in rows:
        if row.sample in seen:
            duplicates.append(row.sample)
        seen.add(row.sample)
    if duplicates:
        raise ValueError(f"{path} contains duplicate sample IDs: {', '.join(sorted(set(duplicates)))}")
    return rows


def _split_samples(value: str) -> tuple[str, ...]:
    if not value:
        return ()
    return tuple(part.strip() for part in re.split(r"[,;]", value) if part.strip())


def read_comparison_table(path: str | Path) -> list[ComparisonRecord]:
    records = []
    for line_no, row in enumerate(_read_tsv(path), start=2):
        comparison = _first(row, ("comparison", "comparison_id", "name"))
        cond1 = _first(row, ("cond1", "condition1", "target", "target_label"))
        cond2 = _first(row, ("cond2", "condition2", "baseline", "baseline_label"))
        cond1_samples = _split_samples(_first(row, ("cond1_samples", "condition1_samples", "target_samples", "target_id")))
        cond2_samples = _split_samples(_first(row, ("cond2_samples", "condition2_samples", "baseline_samples", "baseline_ids")))
        if not cond2_samples:
            legacy_1 = _first(row, ("baseline_id_1",))
            legacy_2 = _first(row, ("baseline_id_2",))
            cond2_samples = tuple(value for value in (legacy_1, legacy_2) if value)
        if not cond1 or not cond2:
            raise ValueError(f"{path}:{line_no} requires cond1 and cond2 columns")
        if not comparison:
            comparison = f"{cond1}_vs_{cond2}"
        records.append(
            ComparisonRecord(
                comparison=comparison,
                cond1=cond1,
                cond2=cond2,
                cond1_samples=cond1_samples,
                cond2_samples=cond2_samples,
            )
        )
    if not records:
        raise ValueError(f"{path} does not contain any comparison rows")
    seen = set()
    duplicates = []
    for record in records:
        if record.comparison in seen:
            duplicates.append(record.comparison)
        seen.add(record.comparison)
    if duplicates:
        raise ValueError(f"{path} contains duplicate comparisons: {', '.join(sorted(set(duplicates)))}")
    return records


def samples_for_condition(samples: list[SampleRecord], condition: str, explicit_samples: tuple[str, ...] = ()) -> list[SampleRecord]:
    sample_by_id = {row.sample: row for row in samples}
    if explicit_samples:
        missing = [sample for sample in explicit_samples if sample not in sample_by_id]
        if missing:
            raise ValueError(f"Unknown sample(s) in comparison table: {', '.join(missing)}")
        selected = [sample_by_id[sample] for sample in explicit_samples]
        wrong = [row.sample for row in selected if row.condition != condition]
        if wrong:
            raise ValueError(
                "Comparison table sample(s) do not match condition {0}: {1}".format(
                    condition,
                    ", ".join(wrong),
                )
            )
        return selected
    selected = [row for row in samples if row.condition == condition]
    if not selected:
        raise ValueError(f"No samples found for condition: {condition}")
    return selected


def write_analysis_peaks(merged_peaks: str | Path, output: str | Path, drop_chroms: tuple[str, ...] = MITO_CHROMS) -> Path:
    merged = Path(merged_peaks).expanduser()
    output_path = Path(output).expanduser()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    drop = set(drop_chroms)
    with merged.open("r", encoding="utf-8") as source, output_path.open("w", encoding="utf-8") as target:
        for line in source:
            if not line.strip() or line.startswith("#"):
                continue
            chrom = line.split("\t", 1)[0]
            if chrom in drop:
                continue
            target.write(line)
    return output_path
