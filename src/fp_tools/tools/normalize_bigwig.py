#!/usr/bin/env python
"""Normalize bigWig signal tracks using shared background regions."""

from __future__ import annotations

import argparse
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pyBigWig


@dataclass
class BackgroundStats:
    sample: str
    input_bigwig: str
    output_bigwig: str
    background_median: float
    background_q90: float
    background_q95: float
    background_q97_5: float
    background_q99: float
    background_mad: float
    background_iqr: float
    scaling_stat: str
    scaling_value: float
    target_scaling_value: float
    scale_factor: float


def _safe_stem(path: str | Path) -> str:
    name = Path(path).name
    for suffix in (".bigWig", ".bigwig", ".bw"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    label = re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("_")
    return label or "track"


def _read_background(path: str | Path) -> list[tuple[str, int, int]]:
    regions: list[tuple[str, int, int]] = []
    with Path(path).expanduser().open("r", encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"Background BED line {line_no} has fewer than 3 columns")
            chrom, start_text, end_text = fields[:3]
            start = int(start_text)
            end = int(end_text)
            if start < 0 or end <= start:
                raise ValueError(f"Background BED line {line_no} has invalid coordinates")
            regions.append((chrom, start, end))
    if not regions:
        raise ValueError(f"No usable background regions found in {path}")
    return regions


def _read_chrom_sizes(path: str | Path) -> dict[str, int]:
    sizes: dict[str, int] = {}
    with Path(path).expanduser().open("r", encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"Chrom sizes line {line_no} has fewer than 2 columns")
            sizes[fields[0]] = int(fields[1])
    if not sizes:
        raise ValueError(f"No chromosome sizes found in {path}")
    return sizes


def _parse_quantile_stat(stat: str) -> float | None:
    if not stat.startswith("q"):
        return None
    try:
        quantile = float(stat[1:])
    except ValueError as exc:
        raise ValueError(f"Unsupported statistic: {stat}") from exc
    if quantile <= 0 or quantile >= 100:
        raise ValueError(f"Quantile statistic must be between q0 and q100, got {stat}")
    return quantile / 100.0


def _format_stat_for_filename(stat: str) -> str:
    return stat.replace(".", "_")


def _stat(values: np.ndarray, stat: str) -> float:
    if stat == "median":
        return float(np.median(values))
    if stat == "iqr":
        q25, q75 = np.quantile(values, [0.25, 0.75])
        return float(q75 - q25)
    quantile = _parse_quantile_stat(stat)
    if quantile is not None:
        return float(np.quantile(values, quantile))
    raise ValueError(f"Unsupported statistic: {stat}")


def _background_values(bigwig: str | Path, regions: list[tuple[str, int, int]]) -> np.ndarray:
    values: list[float] = []
    with pyBigWig.open(str(bigwig)) as bw:
        chroms = bw.chroms()
        for chrom, start, end in regions:
            if chrom not in chroms:
                continue
            clipped_end = min(end, chroms[chrom])
            if clipped_end <= start:
                continue
            summary = bw.stats(chrom, start, clipped_end, type="mean")
            value = summary[0] if summary else None
            values.append(0.0 if value is None or not math.isfinite(value) else float(value))
    if not values:
        raise ValueError(f"No background signal values could be read from {bigwig}")
    return np.asarray(values, dtype=float)


def _fit_background_stats(
    bigwigs: list[str],
    background_regions: list[tuple[str, int, int]],
    method: str,
    stat: str,
    target: str,
) -> list[dict[str, float]]:
    rows = []
    stat_values = []
    for bigwig in bigwigs:
        values = _background_values(bigwig, background_regions)
        median = float(np.median(values))
        q25, q75 = np.quantile(values, [0.25, 0.75])
        mad = float(np.median(np.abs(values - median)))
        iqr = float(q75 - q25)
        chosen = _stat(values, stat)
        rows.append(
            {
                "background_median": median,
                "background_q90": float(np.quantile(values, 0.9)),
                "background_q95": float(np.quantile(values, 0.95)),
                "background_q97_5": float(np.quantile(values, 0.975)),
                "background_q99": float(np.quantile(values, 0.99)),
                "background_mad": mad,
                "background_iqr": iqr,
                "scaling_stat": stat,
                "chosen_stat": chosen,
            }
        )
        stat_values.append(chosen)

    if method == "none":
        for row in rows:
            row["scale_factor"] = 1.0
        return rows

    if method == "background-zscore":
        for row in rows:
            row["scale_factor"] = 1.0
        return rows

    finite_stats = np.asarray([value for value in stat_values if math.isfinite(value)], dtype=float)
    if finite_stats.size == 0:
        raise ValueError("No finite background statistics were available for scaling")
    target_value = float(np.median(finite_stats) if target == "median" else np.mean(finite_stats))
    if target_value <= 0:
        raise ValueError("Target background statistic is <= 0; cannot calculate scale factors")
    for row in rows:
        sample_value = row["chosen_stat"]
        if not math.isfinite(sample_value) or sample_value <= 0:
            raise ValueError("Sample background statistic is <= 0; cannot calculate scale factor")
        row["target_scaling_value"] = target_value
        row["scale_factor"] = target_value / sample_value
    return rows


def _output_path(input_bigwig: str | Path, outdir: Path, method: str, stat: str) -> Path:
    suffix = method.replace("-", "_")
    if method == "background-scale":
        suffix = f"{suffix}_{_format_stat_for_filename(stat)}"
    return outdir / f"{_safe_stem(input_bigwig)}.{suffix}.bw"


def _transform_values(values: np.ndarray, method: str, stats: dict[str, float]) -> np.ndarray:
    values = np.nan_to_num(values.astype(float, copy=False), nan=0.0, posinf=0.0, neginf=0.0)
    if method == "none":
        return values
    if method == "background-scale":
        return values * stats["scale_factor"]
    if method == "background-zscore":
        mad = stats["background_mad"]
        if not math.isfinite(mad) or mad <= 0:
            mad = 1.0
        return (values - stats["background_median"]) / mad
    raise ValueError(f"Unsupported normalization method: {method}")


def _header_from_bigwig(input_bigwig: str | Path, chrom_sizes: dict[str, int] | None = None) -> list[tuple[str, int]]:
    with pyBigWig.open(str(input_bigwig)) as bw:
        input_chroms = bw.chroms()
    if chrom_sizes is None:
        return list(input_chroms.items())
    missing = sorted(set(input_chroms) - set(chrom_sizes))
    if missing:
        raise ValueError(f"--chrom-sizes is missing contigs present in {input_bigwig}: {', '.join(missing[:5])}")
    return [(chrom, chrom_sizes[chrom]) for chrom in input_chroms]


def _write_normalized_bigwig(
    input_bigwig: str | Path,
    output_bigwig: str | Path,
    method: str,
    stats: dict[str, float],
    chrom_sizes: dict[str, int] | None = None,
) -> None:
    header = _header_from_bigwig(input_bigwig, chrom_sizes)
    with pyBigWig.open(str(input_bigwig)) as source, pyBigWig.open(str(output_bigwig), "w") as target:
        target.addHeader(header)
        for chrom, _size in header:
            intervals = source.intervals(chrom)
            if not intervals:
                continue
            starts = [int(start) for start, _end, _value in intervals]
            ends = [int(end) for _start, end, _value in intervals]
            values = _transform_values(np.asarray([value for _start, _end, value in intervals], dtype=float), method, stats)
            target.addEntries([chrom] * len(starts), starts, ends=ends, values=values.tolist())


def normalize_bigwigs(
    bigwigs: list[str],
    background: str | Path,
    outdir: str | Path,
    method: str = "background-scale",
    stat: str = "q90",
    target: str = "median",
    chrom_sizes: str | Path | None = None,
    warn_scale_low: float = 0.5,
    warn_scale_high: float = 2.0,
) -> list[BackgroundStats]:
    if not bigwigs:
        raise ValueError("--bigwigs requires at least one input bigWig")
    outdir_path = Path(outdir).expanduser()
    outdir_path.mkdir(parents=True, exist_ok=True)
    background_regions = _read_background(background)
    chrom_size_dict = _read_chrom_sizes(chrom_sizes) if chrom_sizes else None
    fitted = _fit_background_stats(bigwigs, background_regions, method, stat, target)

    rows: list[BackgroundStats] = []
    for bigwig, stats in zip(bigwigs, fitted):
        output = _output_path(bigwig, outdir_path, method, stat)
        _write_normalized_bigwig(bigwig, output, method, stats, chrom_size_dict)
        scale = stats["scale_factor"]
        if method == "background-scale" and (scale < warn_scale_low or scale > warn_scale_high):
            print(
                f"WARNING: scale factor for {bigwig} is {scale:.4g}; check library quality or background regions.",
                file=sys.stderr,
            )
        rows.append(
            BackgroundStats(
                sample=_safe_stem(bigwig),
                input_bigwig=str(bigwig),
                output_bigwig=str(output),
                background_median=stats["background_median"],
                background_q90=stats["background_q90"],
                background_q95=stats["background_q95"],
                background_q97_5=stats["background_q97_5"],
                background_q99=stats["background_q99"],
                background_mad=stats["background_mad"],
                background_iqr=stats["background_iqr"],
                scaling_stat=str(stats["scaling_stat"]),
                scaling_value=stats["chosen_stat"],
                target_scaling_value=stats.get("target_scaling_value", stats["chosen_stat"]),
                scale_factor=scale,
            )
        )

    _write_qc_tables(rows, outdir_path)
    return rows


def _write_qc_tables(rows: list[BackgroundStats], outdir: Path) -> None:
    header = [
        "sample",
        "input_bigwig",
        "output_bigwig",
        "background_median",
        "background_q90",
        "background_q95",
        "background_q97_5",
        "background_q99",
        "background_mad",
        "background_iqr",
        "scaling_stat",
        "scaling_value",
        "target_scaling_value",
        "scale_factor",
    ]
    qc_path = outdir / "normalize_bigwig_qc.tsv"
    manifest_path = outdir / "normalize_bigwig_manifest.tsv"
    with qc_path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write(
                "\t".join(
                    [
                        row.sample,
                        row.input_bigwig,
                        row.output_bigwig,
                        f"{row.background_median:.10g}",
                        f"{row.background_q90:.10g}",
                        f"{row.background_q95:.10g}",
                        f"{row.background_q97_5:.10g}",
                        f"{row.background_q99:.10g}",
                        f"{row.background_mad:.10g}",
                        f"{row.background_iqr:.10g}",
                        row.scaling_stat,
                        f"{row.scaling_value:.10g}",
                        f"{row.target_scaling_value:.10g}",
                        f"{row.scale_factor:.10g}",
                    ]
                )
                + "\n"
            )
    with manifest_path.open("w", encoding="utf-8") as handle:
        handle.write("sample\tinput_bigwig\toutput_bigwig\n")
        for row in rows:
            handle.write(f"{row.sample}\t{row.input_bigwig}\t{row.output_bigwig}\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Normalize bigWig tracks using robust statistics from shared background BED regions. "
            "For corrected cut-site bigWigs, the recommended method is background-scale."
        )
    )
    parser.add_argument("--bigwigs", nargs="+", required=True, help="Input bigWig files to normalize together.")
    parser.add_argument("--background", required=True, help="Shared background BED used to estimate sample statistics.")
    parser.add_argument("--outdir", required=True, help="Output directory for normalized bigWigs and QC tables.")
    parser.add_argument(
        "--method",
        choices=["background-scale", "background-zscore", "none"],
        default="background-scale",
        help="Normalization method (default: background-scale).",
    )
    parser.add_argument(
        "--stat",
        default="q90",
        help=(
            "Background statistic used by background-scale (default: q90). "
            "Use median, iqr, or quantiles such as q90, q95, q97.5, or q99."
        ),
    )
    parser.add_argument(
        "--target",
        choices=["median", "mean"],
        default="median",
        help="Across-sample target statistic for background-scale (default: median).",
    )
    parser.add_argument("--chrom-sizes", help="Optional chromosome sizes file for output validation/header.")
    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    normalize_bigwigs(
        bigwigs=args.bigwigs,
        background=args.background,
        outdir=args.outdir,
        method=args.method,
        stat=args.stat,
        target=args.target,
        chrom_sizes=args.chrom_sizes,
    )


if __name__ == "__main__":
    main()
