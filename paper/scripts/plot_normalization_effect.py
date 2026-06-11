#!/usr/bin/env python
"""Build a raw-vs-normalized aggregate figure for the manuscript."""

from __future__ import annotations

import argparse
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pyBigWig

from fp_tools.tools.plot_aggregate import (
    apply_quantile_normalization_to_signal_dict,
    build_condition_groups,
    calculate_group_aggregates,
    plot_normalization_comparison,
)
from fp_tools.utils.regions import OneRegion, RegionList


def read_signal_dict(tfbs: Path, signal_paths: dict[str, Path], flank: int, limit: int) -> tuple[dict, dict, list, dict]:
    regions = RegionList().from_bed(str(tfbs))[:limit]
    for region in regions:
        OneRegion.set_width(region, flank * 2)
    regions_dict = {"IRF1 sites": regions}
    region_names = ["IRF1 sites"]
    signal_dict: dict[str, dict[tuple, np.ndarray]] = {}
    for name, path in signal_paths.items():
        bw = pyBigWig.open(str(path))
        try:
            signal_dict[name] = {}
            for region in regions:
                signal_dict[name][region.tup()] = np.asarray(region.get_signal(bw), dtype=float)
        finally:
            bw.close()
    motif_widths = {"IRF1 sites": 10}
    return signal_dict, regions_dict, region_names, motif_widths


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tfbs", default="test_data/IRF1_all.bed")
    parser.add_argument("--bcell", default="test_data/Bcell_footprints.bw")
    parser.add_argument("--tcell", default="test_data/Tcell_footprints.bw")
    parser.add_argument("--out-prefix", default="paper/manuscript/figures/normalization_effect")
    parser.add_argument("--flank", type=int, default=80)
    parser.add_argument("--limit", type=int, default=1200)
    args = parser.parse_args(argv)

    signal_paths = {
        "Bcell_rep1": Path(args.bcell),
        "Bcell_rep2": Path(args.bcell),
        "Tcell_rep1": Path(args.tcell),
        "Tcell_rep2": Path(args.tcell),
    }
    signal_dict, regions_dict, region_names, motif_widths = read_signal_dict(
        Path(args.tfbs), signal_paths, args.flank, args.limit
    )
    scale_factors = {
        "Bcell_rep1": 0.60,
        "Bcell_rep2": 1.40,
        "Tcell_rep1": 0.70,
        "Tcell_rep2": 1.60,
    }
    for sample, factor in scale_factors.items():
        for key in signal_dict[sample]:
            signal_dict[sample][key] = signal_dict[sample][key] * factor

    sample_names = list(signal_dict)
    condition_names, condition_groups = build_condition_groups(
        sample_names,
        ["Bcell", "Bcell", "Tcell", "Tcell"],
    )
    run_args = SimpleNamespace(
        width=args.flank * 2,
        remove_outliers=0.99,
        log_transform=False,
        normalize=False,
        smooth=3,
        flank=args.flank,
    )
    raw, _, _, raw_stats = calculate_group_aggregates(
        signal_dict, regions_dict, region_names, condition_names, condition_groups, motif_widths, run_args
    )
    normalized_signal = apply_quantile_normalization_to_signal_dict(
        signal_dict,
        region_names,
        regions_dict,
        sample_names,
        condition_names,
        condition_groups,
        "sample-quantile",
        logger=None,
    )
    normalized, _, _, norm_stats = calculate_group_aggregates(
        normalized_signal, regions_dict, region_names, condition_names, condition_groups, motif_widths, run_args
    )

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    for suffix in (".png", ".pdf"):
        plot_normalization_comparison(
            raw,
            normalized,
            condition_names,
            region_names,
            out_prefix.with_suffix(suffix),
            title="IRF1 aggregate profiles before and after sample-quantile normalization",
        )
    stats_path = out_prefix.with_suffix(".tsv")
    with stats_path.open("w", encoding="utf-8") as handle:
        handle.write("normalization\tcondition\tregions\tn_replicates\tmean_profile\tmean_profile_sd\tmean_flank\tmean_center\taggregate_fp_score\n")
        for label, rows in (("raw", raw_stats), ("sample_quantile", norm_stats)):
            for row in rows:
                handle.write(label + "\t" + "\t".join(str(row[col]) for col in [
                    "condition", "regions", "n_replicates", "mean_profile",
                    "mean_profile_sd", "mean_flank", "mean_center", "aggregate_fp_score",
                ]) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
