import argparse
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from fp_tools.parsers import add_aggregate_arguments
from fp_tools.utils.normalization import normalize_arrays
from fp_tools.tools.plot_aggregate import (
    apply_quantile_normalization_to_signal_dict,
    build_condition_groups,
    calculate_group_aggregates,
    normalize_site_profiles,
    profile_shape_diagnostics,
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

    def test_parser_accepts_site_shape_options(self):
        parser = add_aggregate_arguments(argparse.ArgumentParser())
        args = parser.parse_args([
            "--site-normalization", "flank-rms",
            "--site-normalization-clip", "4",
            "--show-site-ci",
            "--shape-diagnostics",
        ])
        self.assertEqual(args.site_normalization, "flank-rms")
        self.assertEqual(args.site_normalization_clip, 4.0)
        self.assertTrue(args.show_site_ci)
        self.assertTrue(args.shape_diagnostics)

    def test_flank_rms_normalization_reveals_shared_shape(self):
        profile = np.r_[np.ones(20) * 2, np.zeros(20), np.ones(20) * 2]
        matrix = np.vstack([profile, profile, profile, profile * 10])
        normalized = normalize_site_profiles(matrix, mode="flank-rms", clip=5)
        np.testing.assert_allclose(normalized[0], normalized[-1])
        self.assertLess(float(np.mean(normalized[:, 20:40])), float(np.mean(normalized[:, :20])))

    def test_shape_diagnostics_separate_strong_and_underpowered(self):
        rng = np.random.default_rng(7)
        strong = rng.normal(0, 0.1, size=(100, 120))
        strong[:, 52:68] -= 1.0
        diagnostics = profile_shape_diagnostics(strong, motif_width=16, flank=60)
        self.assertEqual(diagnostics["detectability"], "strong")
        self.assertGreater(diagnostics["depletion_z"], 5)
        underpowered = profile_shape_diagnostics(strong[:10], motif_width=16, flank=60)
        self.assertEqual(underpowered["detectability"], "underpowered")

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
        raw, _, _, _, _ = calculate_group_aggregates(signal_dict, regions, region_names, condition_names, groups, {"sites": 1}, args)
        normalized_signal = apply_quantile_normalization_to_signal_dict(
            signal_dict, region_names, regions, list(signal_dict), condition_names, groups, "sample-quantile", logger=None
        )
        norm, sd, _, _, stats = calculate_group_aggregates(normalized_signal, regions, region_names, condition_names, groups, {"sites": 1}, args)
        before = abs(float(np.mean(raw["B"]["sites"])) - float(np.mean(raw["T"]["sites"])))
        after = abs(float(np.mean(norm["B"]["sites"])) - float(np.mean(norm["T"]["sites"])))
        self.assertLess(after, before)
        self.assertTrue(np.isfinite(sd["B"]["sites"]).all())
        self.assertEqual({row["condition"] for row in stats}, {"B", "T"})

    def test_sample_quantile_matches_shared_normalization_helper(self):
        regions = {"sites": [DummyRegion("r1"), DummyRegion("r2")]}
        region_names = ["sites"]
        signal_dict = {
            "B_rep1": {"r1": np.array([1.0, 2.0, 1.0]), "r2": np.array([2.0, 3.0, 2.0])},
            "B_rep2": {"r1": np.array([1.5, 2.5, 1.5]), "r2": np.array([2.5, 3.5, 2.5])},
            "T_rep1": {"r1": np.array([6.0, 8.0, 6.0]), "r2": np.array([8.0, 10.0, 8.0])},
            "T_rep2": {"r1": np.array([7.0, 9.0, 7.0]), "r2": np.array([9.0, 11.0, 9.0])},
        }
        sample_names = list(signal_dict)
        condition_names, groups = build_condition_groups(sample_names, ["B", "B", "T", "T"])
        normalized_signal = apply_quantile_normalization_to_signal_dict(
            signal_dict, region_names, regions, sample_names, condition_names, groups, "sample-quantile", logger=None
        )
        arrays = [np.concatenate([signal_dict[name][reg.tup()] for rid in region_names for reg in regions[rid]]) for name in sample_names]
        expected_arrays, _, _ = normalize_arrays(arrays, sample_names, mode="sample-quantile", logger=None)
        for name, expected in zip(sample_names, expected_arrays):
            observed = np.concatenate([normalized_signal[name][reg.tup()] for rid in region_names for reg in regions[rid]])
            np.testing.assert_allclose(observed, expected, rtol=1e-7, atol=1e-7)

    def test_condition_quantile_matches_condition_level_normalizers(self):
        from fp_tools.utils.normalization import fit_quantile_normalizers

        regions = {"sites": [DummyRegion("r1"), DummyRegion("r2")]}
        region_names = ["sites"]
        signal_dict = {
            "B_rep1": {"r1": np.array([1.0, 2.0, 1.0]), "r2": np.array([2.0, 3.0, 2.0])},
            "B_rep2": {"r1": np.array([1.5, 2.5, 1.5]), "r2": np.array([2.5, 3.5, 2.5])},
            "T_rep1": {"r1": np.array([6.0, 8.0, 6.0]), "r2": np.array([8.0, 10.0, 8.0])},
            "T_rep2": {"r1": np.array([7.0, 9.0, 7.0]), "r2": np.array([9.0, 11.0, 9.0])},
        }
        sample_names = list(signal_dict)
        condition_names, groups = build_condition_groups(sample_names, ["B", "B", "T", "T"])
        normalized_signal = apply_quantile_normalization_to_signal_dict(
            signal_dict, region_names, regions, sample_names, condition_names, groups, "condition-quantile", logger=None
        )
        condition_arrays = []
        for condition in condition_names:
            sample_arrays = [np.concatenate([signal_dict[name][reg.tup()] for rid in region_names for reg in regions[rid]]) for name in groups[condition]]
            condition_arrays.append(np.mean(np.vstack(sample_arrays), axis=0))
        norm_objects, _ = fit_quantile_normalizers(condition_arrays, condition_names, logger=None)
        for condition in condition_names:
            for name in groups[condition]:
                original = np.concatenate([signal_dict[name][reg.tup()] for rid in region_names for reg in regions[rid]])
                expected = np.maximum(0.0, norm_objects[condition].normalize(original))
                observed = np.concatenate([normalized_signal[name][reg.tup()] for rid in region_names for reg in regions[rid]])
                np.testing.assert_allclose(observed, expected, rtol=1e-7, atol=1e-7)

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
