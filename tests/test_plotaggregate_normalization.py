import argparse
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from fp_tools.parsers import add_aggregate_arguments
from fp_tools.tools.plot_aggregate import (
    apply_quantile_normalization_to_signal_dict,
    build_condition_groups,
    calculate_group_aggregates,
    plot_normalization_comparison,
)


class DummyRegion:
    def __init__(self, name):
        self.name = name

    def tup(self):
        return self.name


class PlotAggregateNormalizationTest(unittest.TestCase):
    def test_parser_accepts_replicate_normalization_options(self):
        parser = add_aggregate_arguments(argparse.ArgumentParser())
        args = parser.parse_args([
            "--cond-names", "B", "B", "T", "T",
            "--normalization", "sample-quantile",
            "--normalization-comparison-output", "norm.png",
            "--output_aggregated_stats", "stats.csv",
            "--show-replicate-sd",
        ])
        self.assertEqual(args.cond_names, ["B", "B", "T", "T"])
        self.assertEqual(args.normalization, "sample-quantile")
        self.assertTrue(args.show_replicate_sd)

    def test_repeated_condition_names_group_replicates(self):
        conditions, groups = build_condition_groups(
            ["B_rep1", "B_rep2", "T_rep1"],
            ["B", "B", "T"],
        )
        self.assertEqual(conditions, ["B", "T"])
        self.assertEqual(groups["B"], ["B_rep1", "B_rep2"])

    def test_sample_quantile_normalization_changes_scaled_aggregates(self):
        regions = {"sites": [DummyRegion("r1"), DummyRegion("r2")]}
        region_names = ["sites"]
        signal_dict = {
            "B_rep1": {"r1": np.array([1.0, 2.0, 1.0]), "r2": np.array([2.0, 3.0, 2.0])},
            "B_rep2": {"r1": np.array([2.0, 4.0, 2.0]), "r2": np.array([4.0, 6.0, 4.0])},
            "T_rep1": {"r1": np.array([1.5, 2.5, 1.5]), "r2": np.array([2.5, 3.5, 2.5])},
            "T_rep2": {"r1": np.array([6.0, 10.0, 6.0]), "r2": np.array([10.0, 14.0, 10.0])},
        }
        condition_names, groups = build_condition_groups(
            list(signal_dict),
            ["B", "B", "T", "T"],
        )
        args = SimpleNamespace(width=3, remove_outliers=1, log_transform=False, normalize=False, smooth=1, flank=1)
        raw, _, _, _ = calculate_group_aggregates(signal_dict, regions, region_names, condition_names, groups, {"sites": 1}, args)
        normalized_signal = apply_quantile_normalization_to_signal_dict(
            signal_dict, region_names, regions, list(signal_dict), condition_names, groups, "sample-quantile", logger=None
        )
        norm, sd, _, stats = calculate_group_aggregates(normalized_signal, regions, region_names, condition_names, groups, {"sites": 1}, args)
        before = abs(float(np.mean(raw["B"]["sites"])) - float(np.mean(raw["T"]["sites"])))
        after = abs(float(np.mean(norm["B"]["sites"])) - float(np.mean(norm["T"]["sites"])))
        self.assertLess(after, before)
        self.assertTrue(np.isfinite(sd["B"]["sites"]).all())
        self.assertEqual({row["condition"] for row in stats}, {"B", "T"})

    def test_comparison_figure_writes_file(self):
        profile = np.array([1.0, 2.0, 1.0])
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "comparison.png"
            plot_normalization_comparison(
                {"B": {"sites": profile}},
                {"B": {"sites": profile * 0.5}},
                ["B"],
                ["sites"],
                output,
            )
            self.assertTrue(output.exists())
            self.assertGreater(output.stat().st_size, 1000)


if __name__ == "__main__":
    unittest.main()
