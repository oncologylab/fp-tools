import argparse
import tempfile
import unittest
from pathlib import Path

import numpy as np

from fp_tools.parsers import add_aggregate_arguments, add_scorebigwig_arguments
from fp_tools.tools.plot_aggregate import default_multiscale_output, plot_multiscale_aggregate_npz
from fp_tools.utils.multiscale import (
    aggregate_multiscale_tensor, load_multiscale_npz, multiscale_depletion,
    hybrid_footprint_score, parse_scales, summarize_multiscale, symmetric_depletion,
    trim_multiscale_features, write_multiscale_npz,
)


class MultiscaleScoringTest(unittest.TestCase):
    def test_parse_scales_sorts_deduplicates_and_validates(self):
        self.assertEqual(parse_scales([32, 8, 8, 16]), (8, 16, 32))
        with self.assertRaises(ValueError):
            parse_scales([2])

    def test_synthetic_depletion_scores_peak_at_center(self):
        signal = np.ones(101, dtype=float) * 5.0
        signal[48:53] = 0.0
        features = multiscale_depletion(signal, [8, 16])
        summary = summarize_multiscale(features)
        self.assertGreater(summary[50], 2.5)
        self.assertAlmostEqual(float(summary[50]), float(np.max(summary)), places=6)


    def test_multiscale_npz_roundtrip_and_aggregate(self):
        records = [
            (("chr1", 10, 13), {8: np.array([1.0, 2.0, 3.0]), 16: np.array([4.0, 5.0, 6.0])}),
            (("chr1", 20, 22), {8: np.array([3.0, 5.0]), 16: np.array([7.0, 9.0])}),
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir) / "multiscale.npz"
            write_multiscale_npz(str(out), records, [16, 8], "max")
            data = load_multiscale_npz(str(out))

        self.assertEqual(data["format_version"].item(), "fp-tools-multiscale-npz-v1")
        self.assertEqual(data["scales"].tolist(), [8, 16])
        self.assertEqual(data["offsets"].tolist(), [0, 3, 5])
        self.assertEqual(data["tensor"].shape, (2, 5))
        aggregate = aggregate_multiscale_tensor(data)
        self.assertEqual(aggregate.shape, (2, 3))
        self.assertAlmostEqual(float(aggregate[0, 0]), 2.0)
        self.assertAlmostEqual(float(aggregate[0, 2]), 3.0)

    def test_trim_multiscale_features(self):
        features = {8: np.arange(8), 16: np.arange(8) + 10}
        trimmed = trim_multiscale_features(features, 2)
        self.assertEqual(trimmed[8].tolist(), [2, 3, 4, 5])

    def test_multiscale_parser_options_are_available(self):
        parser = add_scorebigwig_arguments(argparse.ArgumentParser())
        args = parser.parse_args(["--score", "multiscale", "--scales", "8", "16", "--multiscale-summary", "mean"])
        self.assertEqual(args.score, "multiscale")
        self.assertEqual(args.scales, [8, 16])
        self.assertEqual(args.multiscale_summary, "mean")

    def test_symmetric_depletion_and_hybrid_recover_wide_shape(self):
        signal = np.ones(201, dtype=float) * 4.0
        signal[84:117] = -2.0
        symmetric = symmetric_depletion(signal, center_width=33, flank_width=32)
        self.assertEqual(int(np.argmax(symmetric)), 100)
        self.assertGreater(symmetric[100], 2.0)
        legacy = np.zeros_like(signal)
        hybrid = hybrid_footprint_score(signal, legacy, weight=0.2)
        self.assertAlmostEqual(float(hybrid[100]), float(0.2 * symmetric[100]))

    def test_hybrid_parser_options_are_available(self):
        parser = add_scorebigwig_arguments(argparse.ArgumentParser())
        args = parser.parse_args([
            "--score", "hybrid", "--hybrid-center-width", "35",
            "--hybrid-flank-width", "30", "--hybrid-weight", "0.1",
        ])
        self.assertEqual(args.score, "hybrid")
        self.assertEqual(args.hybrid_center_width, 35)
        self.assertEqual(args.hybrid_flank_width, 30)
        self.assertEqual(args.hybrid_weight, 0.1)

    def test_plotaggregate_multiscale_parser_options_are_available(self):
        parser = add_aggregate_arguments(argparse.ArgumentParser())
        args = parser.parse_args([
            "--multiscale-npz",
            "tensor.npz",
            "--output-multiscale-aggregate",
            "aggregate_multiscale.pdf",
        ])
        self.assertEqual(args.multiscale_npz, "tensor.npz")
        self.assertEqual(args.output_multiscale_aggregate, "aggregate_multiscale.pdf")
        self.assertEqual(default_multiscale_output("aggregate.pdf"), "aggregate_multiscale.pdf")

    def test_plot_multiscale_aggregate_npz_writes_figure(self):
        records = [
            (("chr1", 10, 13), {8: np.array([1.0, 2.0, 3.0]), 16: np.array([4.0, 5.0, 6.0])}),
            (("chr1", 20, 23), {8: np.array([2.0, 4.0, 6.0]), 16: np.array([8.0, 10.0, 12.0])}),
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            npz_path = tmp / "multiscale.npz"
            output = tmp / "plotaggregate_multiscale.pdf"
            write_multiscale_npz(str(npz_path), records, [8, 16], "max")
            returned = plot_multiscale_aggregate_npz(npz_path, output)
            self.assertEqual(returned, output)
            self.assertTrue(output.exists())
            self.assertGreater(output.stat().st_size, 1000)


if __name__ == "__main__":
    unittest.main()
