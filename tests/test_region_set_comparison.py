import argparse
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
import pyBigWig
import pysam

from fp_tools.parsers import add_diff_footprints_arguments
from fp_tools.tools.region_set_comparison import (
    _read_region_sets,
    _stratified_effect,
    run_region_set_comparison,
)
from fp_tools.tools.diff_footprint_helpers import select_aggregate_rows


class RegionSetComparisonTest(unittest.TestCase):
    def test_parser_accepts_region_axis_inputs(self):
        parser = add_diff_footprints_arguments(argparse.ArgumentParser(prog="diff-footprints"))
        args = parser.parse_args([
            "--comparison-axis", "regions",
            "--signals", "sample.bw",
            "--regions", "bound.bed", "control.bed",
            "--region-labels", "bound", "control",
            "--region-strata-column", "4",
            "--genome", "genome.fa",
        ])
        self.assertEqual(args.comparison_axis, "regions")
        self.assertEqual(args.region_labels, ["bound", "control"])
        self.assertEqual(args.region_strata_column, 4)
        self.assertIsNone(args.plot_aggregate_motifs)
        self.assertEqual(args.default_aggregate_plots, 4)

    def test_parser_accepts_ordered_aggregate_motifs(self):
        parser = add_diff_footprints_arguments(argparse.ArgumentParser(prog="diff-footprints"))
        args = parser.parse_args([
            "--plot-aggregate-motifs", "MA1494.2", "MA0047.4",
            "--default-aggregate-plots", "8",
        ])
        self.assertEqual(args.plot_aggregate_motifs, ["MA1494.2", "MA0047.4"])
        self.assertEqual(args.default_aggregate_plots, 8)

    def test_explicit_aggregate_motifs_preserve_order_and_reject_ambiguity(self):
        rows = pd.DataFrame([
            {"output_prefix": "HNF4A_MA1494.2", "motif_id": "MA1494.2", "name": "HNF4A"},
            {"output_prefix": "HNF4A_MA0114.5", "motif_id": "MA0114.5", "name": "HNF4A"},
            {"output_prefix": "FOXA2_MA0047.4", "motif_id": "MA0047.4", "name": "FOXA2"},
        ])
        selected = select_aggregate_rows(rows, ["MA0047.4", "MA1494.2"])
        self.assertEqual(selected["motif_id"].tolist(), ["MA0047.4", "MA1494.2"])
        with self.assertRaisesRegex(ValueError, "Ambiguous aggregate motif"):
            select_aggregate_rows(rows, ["HNF4A"])

    def test_stratified_effect_respects_strata(self):
        frame = pd.DataFrame({
            "stratum": ["low", "low", "high", "high"],
            "region_set": ["A", "B", "A", "B"],
            "rep1": [2.0, 1.0, 20.0, 10.0],
        })
        self.assertAlmostEqual(_stratified_effect(frame, "A", "B", "rep1"), 5.5)

    def test_overlapping_region_sets_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "a.bed").write_text("chr1\t10\t20\n")
            (root / "b.bed").write_text("chr1\t19\t30\n")
            with self.assertRaisesRegex(ValueError, "mutually non-overlapping"):
                _read_region_sets([root / "a.bed", root / "b.bed"], ["A", "B"])

    def test_end_to_end_single_sample_writes_quantification_and_html(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta = root / "genome.fa"
            fasta.write_text(">chr1\n" + "ACGT" * 600 + "\n")
            pysam.faidx(str(fasta))
            motif = root / "motif.pfm"
            motif.write_text(
                ">M1 TEST\n"
                "A [10 0 0 0]\n"
                "C [0 10 0 0]\n"
                "G [0 0 10 0]\n"
                "T [0 0 0 10]\n"
                ">M2 TEST2\n"
                "A [10 0 0 0]\n"
                "C [0 10 0 0]\n"
                "G [0 0 10 0]\n"
                "T [0 0 0 10]\n"
            )
            regions_a = root / "a.bed"
            regions_b = root / "b.bed"
            regions_a.write_text("".join(f"chr1\t{40+i*20}\t{56+i*20}\t{i%2}\n" for i in range(10)))
            regions_b.write_text("".join(f"chr1\t{500+i*20}\t{516+i*20}\t{i%2}\n" for i in range(10)))
            signal = root / "sample.bw"
            values = np.ones(2400, dtype=float)
            for i in range(10):
                values[40+i*20:56+i*20] = 5.0
            with pyBigWig.open(str(signal), "w") as bw:
                bw.addHeader([("chr1", 2400)])
                bw.addEntries("chr1", 0, values=values.tolist(), span=1, step=1)
            outdir = root / "out"
            args = argparse.Namespace(
                genome=str(fasta), signals=[str(signal)], aggregate_signals=None,
                sample_names=["rep1"], cond_names=None, regions=[str(regions_a), str(regions_b)],
                region_labels=["bound", "control"], region_strata_column=4,
                region_permutations=49, region_bootstrap=40, min_regions_per_set=5,
                random_seed=3, motifs=[str(motif)], naming="name_id", motif_pvalue=0.01,
                normalization="none", norm_off=False, outdir=str(outdir), prefix="test",
                skip_excel=True, plot_aggregate="top", plot_aggregate_top_n=1,
                plot_aggregate_motifs=None, default_aggregate_plots=1,
                aggregate_flank=5, verbosity=0,
            )
            result = run_region_set_comparison(args)
            tested = result.loc[result["output_prefix"].str.startswith("TEST_M1")].iloc[0]
            self.assertGreater(tested["effect"], 3.0)
            self.assertEqual(tested["status"], "tested")
            self.assertTrue((outdir / "test_results.txt").exists())
            report = outdir / "test_bound_vs_control.html"
            self.assertTrue(report.exists())
            html = report.read_text()
            self.assertIn("Region-set footprint report", html)
            self.assertIn("Download results TSV", html)

            signal_2 = root / "sample_2.bw"
            values_2 = np.ones(2400, dtype=float)
            for i in range(10):
                values_2[40+i*20:56+i*20] = 4.5
            with pyBigWig.open(str(signal_2), "w") as bw:
                bw.addHeader([("chr1", 2400)])
                bw.addEntries("chr1", 0, values=values_2.tolist(), span=1, step=1)
            replicated_outdir = root / "replicated_out"
            replicated_args = argparse.Namespace(**vars(args))
            replicated_args.signals = [str(signal), str(signal_2)]
            replicated_args.sample_names = ["rep1", "rep2"]
            replicated_args.outdir = str(replicated_outdir)
            replicated = run_region_set_comparison(replicated_args)
            replicated_tested = replicated.loc[
                replicated["output_prefix"].str.startswith("TEST_M1")
            ].iloc[0]
            self.assertEqual(
                replicated_tested["statistical_method"],
                "paired empirical-Bayes moderated t",
            )
            self.assertTrue(np.isfinite(replicated_tested["pvalue"]))
            replicate_effects = pd.read_csv(
                replicated_outdir / "test_region_replicate_effects.tsv", sep="\t"
            )
            self.assertEqual(set(replicate_effects["sample"]), {"rep1", "rep2"})


if __name__ == "__main__":
    unittest.main()
