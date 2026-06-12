import csv
import importlib.util
import pathlib
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

from fp_tools.utils.multiscale import write_multiscale_npz


ROOT = pathlib.Path(__file__).resolve().parents[1]


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


download_manifest = load_module("download_manifest", ROOT / "benchmarks" / "scripts" / "download_manifest.py")
compute_binary_metrics = load_module("compute_binary_metrics", ROOT / "benchmarks" / "scripts" / "compute_binary_metrics.py")
compute_calibration = load_module("compute_calibration", ROOT / "benchmarks" / "scripts" / "compute_calibration.py")
plot_benchmark_panels = load_module("plot_benchmark_panels", ROOT / "manuscript" / "scripts" / "plot_benchmark_panels.py")
plot_calibration_panels = load_module("plot_calibration_panels", ROOT / "manuscript" / "scripts" / "plot_calibration_panels.py")
plot_multiscale_npz = load_module("plot_multiscale_npz", ROOT / "manuscript" / "scripts" / "plot_multiscale_npz.py")
build_encode_manifest = load_module("build_encode_manifest", ROOT / "benchmarks" / "scripts" / "build_encode_manifest.py")
build_motif_removal_benchmark = load_module("build_motif_removal_benchmark", ROOT / "benchmarks" / "scripts" / "build_motif_removal_benchmark.py")
build_label_overlap_benchmark = load_module("build_label_overlap_benchmark", ROOT / "benchmarks" / "scripts" / "build_label_overlap_benchmark.py")
run_benchmark_pipeline = load_module("run_benchmark_pipeline", ROOT / "benchmarks" / "scripts" / "run_benchmark_pipeline.py")
score_peaks_with_pwm = load_module("score_peaks_with_pwm", ROOT / "benchmarks" / "scripts" / "score_peaks_with_pwm.py")
footprint_from_bam = load_module("footprint_from_bam", ROOT / "benchmarks" / "scripts" / "footprint_from_bam.py")
footprint_occupancy_score = load_module("footprint_occupancy_score", ROOT / "benchmarks" / "scripts" / "footprint_occupancy_score.py")
build_tf_feature_table = load_module("build_tf_feature_table", ROOT / "benchmarks" / "scripts" / "build_tf_feature_table.py")
evaluate_methods = load_module("evaluate_methods", ROOT / "benchmarks" / "scripts" / "evaluate_methods.py")
plot_method_comparison = load_module("plot_method_comparison", ROOT / "manuscript" / "scripts" / "plot_method_comparison.py")


class FeatureTableAndMethodsTest(unittest.TestCase):

    def test_gc_content(self):
        self.assertAlmostEqual(build_tf_feature_table.gc_content("GGCC"), 1.0)
        self.assertAlmostEqual(build_tf_feature_table.gc_content("ATAT"), 0.0)
        self.assertAlmostEqual(build_tf_feature_table.gc_content("ACGT"), 0.5)

    def test_read_accessibility(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            bed = pathlib.Path(tmpdir) / "a.bed"
            bed.write_text("#chrom\tstart\tend\tname\tscore\nchr1\t0\t10\tp1\t3.5\nchr1\t20\t30\tp2\t7.0\n", encoding="utf-8")
            acc = build_tf_feature_table.read_accessibility(bed)
            self.assertEqual(acc, {"p1": 3.5, "p2": 7.0})

    def test_integrated_oof_scores_separable(self):
        import numpy as np

        rng = np.random.default_rng(0)
        n = 200
        motif = np.concatenate([rng.normal(2, 1, n // 2), rng.normal(-2, 1, n // 2)])
        frame = pd.DataFrame({
            "accessibility": rng.normal(0, 1, n),
            "motif": motif,
            "gc": rng.normal(0.5, 0.1, n),
            "label": np.array([1] * (n // 2) + [0] * (n // 2)),
        })
        scores = evaluate_methods.integrated_oof_scores(frame, ["accessibility", "motif", "gc"], seed=1)
        self.assertEqual(len(scores), n)
        # Positive class should get higher mean score than negative class.
        self.assertGreater(scores[frame["label"] == 1].mean(), scores[frame["label"] == 0].mean())


class FootprintFromBamTest(unittest.TestCase):

    def test_footprint_score_positive_when_center_depleted(self):
        counts = np.array([5, 5, 5, 0, 0, 0, 5, 5, 5], dtype=float)
        # center = indices 3:6 (all zero), flanks of width 3 on each side (all 5s)
        score = footprint_from_bam.footprint_score(counts, 3, 6, 3)
        self.assertGreater(score, 0)

    def test_footprint_score_negative_when_center_enriched(self):
        counts = np.array([0, 0, 0, 9, 9, 9, 0, 0, 0], dtype=float)
        self.assertLess(footprint_from_bam.footprint_score(counts, 3, 6, 3), 0)

    def test_build_cutsites_counts_tn5_insertions(self):
        import pysam

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = pathlib.Path(tmpdir)
            bam_path = tmp / "mini.bam"
            header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chrT", "LN": 200}]}
            with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
                a = pysam.AlignedSegment()
                a.query_name = "r1"
                a.query_sequence = "A" * 20
                a.flag = 0  # forward
                a.reference_id = 0
                a.reference_start = 50
                a.mapping_quality = 60
                a.cigarstring = "20M"
                a.query_qualities = pysam.qualitystring_to_array("I" * 20)
                bam.write(a)
            pysam.index(str(bam_path))
            counts = footprint_from_bam.build_cutsites(str(bam_path), "chrT", 200)
            # forward read start 50 + shift 4 = insertion at 54
            self.assertEqual(int(counts[54]), 1)
            self.assertEqual(int(counts.sum()), 1)


class FootprintOccupancyScoreTest(unittest.TestCase):

    def test_best_match_finds_motif_offset(self):
        from fp_tools.tools.variants import read_pwm_motifs

        with tempfile.TemporaryDirectory() as tmpdir:
            meme = pathlib.Path(tmpdir) / "m.meme"
            meme.write_text(
                "MEME version 4\n\nALPHABET= ACGT\n\nMOTIF m M\n"
                "letter-probability matrix: alength= 4 w= 4 nsites= 1 E= 0\n"
                "0.91 0.03 0.03 0.03\n0.91 0.03 0.03 0.03\n0.91 0.03 0.03 0.03\n0.91 0.03 0.03 0.03\n",
                encoding="utf-8",
            )
            motif = read_pwm_motifs(meme)[0]
            score, off = footprint_occupancy_score.best_match("CCCCAAAACCCC", motif)
            self.assertEqual(off, 4)
            self.assertGreater(score, 0)

    def test_footprint_score_sign(self):
        import numpy as np

        depleted = np.array([4.0, 4.0, 0.0, 0.0, 4.0, 4.0])
        self.assertGreater(footprint_occupancy_score.footprint_score(depleted, 2, 4, 2), 0)


class PlotMethodComparisonTest(unittest.TestCase):

    def test_plot_writes_three_formats(self):
        metrics = pd.DataFrame(
            {
                "group": ["CTCF/accessibility", "CTCF/fp-tools-motif", "global"],
                "auroc": [0.53, 0.87, 0.70],
                "auprc": [0.36, 0.81, 0.58],
            }
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = pathlib.Path(tmpdir) / "cmp"
            outputs = plot_method_comparison.plot_method_comparison(metrics, prefix)
            self.assertEqual(len(outputs), 3)
            for path in outputs:
                self.assertTrue(path.exists())


class ScorePeaksWithPwmTest(unittest.TestCase):

    def test_best_pwm_match_ranks_motif_bearing_peak_higher(self):
        import pysam

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = pathlib.Path(tmpdir)
            # Synthetic contig: a strong A-run motif planted in peak 1, absent in peak 2.
            seq = ("A" * 40) + ("AAAAAAAA") + ("C" * 40) + ("CGCGCGCG") + ("C" * 40)
            fa = tmp / "mini.fa"
            fa.write_text(f">chrT\n{seq}\n", encoding="utf-8")
            pysam.faidx(str(fa))

            # MEME motif favouring poly-A (length 8).
            meme = tmp / "polyA.meme"
            rows = "\n".join("0.97 0.01 0.01 0.01" for _ in range(8))
            meme.write_text(
                "MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\n"
                "Background letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n"
                "MOTIF polyA TEST\nletter-probability matrix: alength= 4 w= 8 nsites= 1 E= 0\n"
                + rows + "\n",
                encoding="utf-8",
            )

            peaks = tmp / "peaks.bed"
            peaks.write_text("chrT\t0\t48\tpeak_A\nchrT\t88\t136\tpeak_C\n", encoding="utf-8")
            out = tmp / "scored.bed"
            n = score_peaks_with_pwm.score_peaks_with_pwm(peaks, fa, meme, out)
            self.assertEqual(n, 2)

            scored = pd.read_csv(out, sep="\t")
            by_name = dict(zip(scored["name"], scored["score"]))
            self.assertGreater(by_name["peak_A"], by_name["peak_C"])

    def test_chroms_not_in_genome_are_skipped(self):
        import pysam

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = pathlib.Path(tmpdir)
            fa = tmp / "mini.fa"
            fa.write_text(">chrT\n" + "ACGT" * 20 + "\n", encoding="utf-8")
            pysam.faidx(str(fa))
            meme = tmp / "m.meme"
            meme.write_text(
                "MEME version 4\n\nALPHABET= ACGT\n\nMOTIF m M\n"
                "letter-probability matrix: alength= 4 w= 4 nsites= 1 E= 0\n"
                "0.7 0.1 0.1 0.1\n0.1 0.7 0.1 0.1\n0.1 0.1 0.7 0.1\n0.1 0.1 0.1 0.7\n",
                encoding="utf-8",
            )
            peaks = tmp / "peaks.bed"
            peaks.write_text("chrT\t0\t40\tkeep\nchrZ\t0\t40\tdrop\n", encoding="utf-8")
            out = tmp / "scored.bed"
            n = score_peaks_with_pwm.score_peaks_with_pwm(peaks, fa, meme, out)
            self.assertEqual(n, 1)
            scored = pd.read_csv(out, sep="\t")
            self.assertEqual(list(scored["name"]), ["keep"])


class BenchmarkScriptsTest(unittest.TestCase):

    def test_encode_manifest_choose_file_prefers_requested_assembly(self):
        files = [
            {"status": "released", "output_type": "alignments", "file_format": "bam", "assembly": "hg19", "accession": "OLD"},
            {"status": "released", "output_type": "alignments", "file_format": "bam", "assembly": "GRCh38", "accession": "NEW"},
        ]
        selected = build_encode_manifest.choose_file(files, "alignments", "bam", assembly="GRCh38")
        self.assertEqual(selected["accession"], "NEW")


    def test_encode_smoke_manifest_is_versioned_and_downloadable(self):
        manifest = ROOT / "benchmarks" / "manifests" / "encode_k562_ctcf_smoke.tsv"
        rows = download_manifest.read_manifest(manifest)
        self.assertEqual(len(rows), 3)
        self.assertEqual({row["source"] for row in rows}, {"ENCODE"})
        self.assertEqual({row["assembly"] for row in rows}, {"GRCh38"})
        self.assertTrue(any(row["file_format"] == "bam" for row in rows))
        self.assertTrue(all(row["url"].startswith("https://www.encodeproject.org/") for row in rows))

    def test_download_manifest_dry_run_writes_report(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = pathlib.Path(tmpdir)
            manifest = tmp / "manifest.tsv"
            report = tmp / "download_report.tsv"
            with manifest.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=["file_accession", "url", "local_path", "checksum"], delimiter="	")
                writer.writeheader()
                writer.writerow({"file_accession": "FILE1", "url": "https://example.org/file.bam", "local_path": str(tmp / "file.bam"), "checksum": ""})
            rows = download_manifest.read_manifest(manifest)
            results = [download_manifest.download_one(row, dry_run=True) for row in rows]
            download_manifest.write_report(results, report)
            text = report.read_text(encoding="utf-8")
            self.assertIn("dry_run", text)
            self.assertIn("FILE1", text)

    def test_compute_binary_metrics_global_and_grouped(self):
        df = pd.DataFrame(
            {
                "label": [1, 0, 1, 0, 1, 0],
                "score": [0.95, 0.05, 0.80, 0.10, 0.70, 0.20],
                "tf": ["A", "A", "A", "A", "B", "B"],
                "cell": ["C", "C", "C", "C", "C", "C"],
                "method": ["fp", "fp", "fp", "fp", "fp", "fp"],
            }
        )
        metrics = compute_binary_metrics.compute_metrics(df, "label", "score", ["tf", "cell", "method"])
        self.assertIn("global", set(metrics["group"]))
        self.assertGreaterEqual(metrics.loc[metrics["group"] == "global", "auprc"].iloc[0], 0.99)

    def test_build_motif_removal_benchmark_zeroes_baseline_and_summarizes(self):
        df = pd.DataFrame(
            {
                "label": [1, 0, 1, 0],
                "tf": ["CTCF", "CTCF", "GATA1", "GATA1"],
                "cell": ["K562", "K562", "K562", "K562"],
                "motif_family": ["CTCF", "CTCF", "GATA", "GATA"],
                "motif_score": [12.0, 8.0, 10.0, 7.0],
                "rank_score": [0.91, 0.15, 0.77, 0.20],
                "candidate_score": [5.0, 1.0, 4.0, 1.5],
            }
        )
        table = build_motif_removal_benchmark.build_motif_removal_table(
            df,
            remove_col="motif_family",
            remove_values=["CTCF"],
            baseline_score_col="motif_score",
            recovery_score_cols=["rank_score", "candidate_score"],
        )
        self.assertEqual(set(table["method"]), {"motif_removed_baseline", "rank_score", "candidate_score"})
        self.assertEqual(len(table), 6)
        baseline = table[table["method"] == "motif_removed_baseline"]
        self.assertTrue((baseline["score"] == 0.0).all())
        self.assertTrue(table["removed"].all())
        self.assertEqual(set(table["removal_target"]), {"CTCF"})

        summary = build_motif_removal_benchmark.summarize_removal_table(table)
        self.assertEqual(set(summary["method"]), {"motif_removed_baseline", "rank_score", "candidate_score"})
        self.assertTrue((summary["n"] == 2).all())

    def test_motif_removal_table_can_include_controls(self):
        df = pd.DataFrame(
            {
                "label": [1, 0, 1],
                "motif_family": ["ETS", "AP1", "ETS"],
                "rank_score": [0.8, 0.2, 0.7],
            }
        )
        table = build_motif_removal_benchmark.build_motif_removal_table(
            df,
            remove_col="motif_family",
            remove_values=["ETS"],
            recovery_score_cols=["rank_score"],
            include_controls=True,
        )
        self.assertEqual(len(table), 3)
        self.assertEqual(table["removed"].tolist(), [True, False, True])

    def test_build_label_overlap_benchmark_outputs_metrics_ready_table(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = pathlib.Path(tmpdir)
            predictions = tmp / "predictions.bed"
            labels = tmp / "labels.bed"
            out = tmp / "benchmark.tsv"
            predictions.write_text(
                "#chrom\tstart\tend\tname\tscore\tmotif_family\n"
                "chr1\t10\t20\tsite1\t0.9\tCTCF\n"
                "chr1\t30\t40\tsite2\t0.2\tCTCF\n"
                "chr2\t5\t15\tsite3\t0.7\tGATA\n",
                encoding="utf-8",
            )
            labels.write_text("chr1\t15\t25\nchr2\t0\t6\n", encoding="utf-8")

            table = build_label_overlap_benchmark.build_label_overlap_table(
                predictions,
                labels,
                out,
                min_overlap_bp=2,
                method="reranked",
                tf="CTCF",
                cell="K562",
                metadata_cols=["name", "motif_family"],
            )
            self.assertTrue(out.exists())

        self.assertEqual(table["label"].tolist(), [1, 0, 0])
        self.assertEqual(table["method"].unique().tolist(), ["reranked"])
        self.assertEqual(table["tf"].unique().tolist(), ["CTCF"])
        self.assertIn("motif_family", table.columns)

    def test_compute_binary_metrics_bootstrap_confidence_intervals(self):
        df = pd.DataFrame(
            {
                "label": [1, 0, 1, 0, 1, 0],
                "score": [0.95, 0.05, 0.80, 0.10, 0.70, 0.20],
                "method": ["fp", "fp", "fp", "fp", "fp", "fp"],
            }
        )
        boot = compute_binary_metrics.bootstrap_confidence_intervals(
            df,
            "label",
            "score",
            ["method"],
            n_bootstrap=25,
            seed=7,
        )
        auprc = boot[(boot["group"] == "global") & (boot["metric"] == "auprc")].iloc[0]
        self.assertEqual(int(auprc["n_bootstrap"]), 25)
        self.assertGreaterEqual(float(auprc["ci_high"]), float(auprc["ci_low"]))
        self.assertGreater(int(auprc["successful_bootstraps"]), 0)

    def test_compute_binary_metrics_allows_raw_non_probability_scores(self):
        df = pd.DataFrame(
            {
                "label": [1, 0, 1, 0],
                "score": [5.0, 1.0, 4.0, 2.0],
                "method": ["raw", "raw", "raw", "raw"],
            }
        )
        metrics = compute_binary_metrics.compute_metrics(df, "label", "score", ["method"])
        global_row = metrics[metrics["group"] == "global"].iloc[0]
        self.assertGreaterEqual(global_row["auprc"], 0.99)
        self.assertTrue(np.isnan(global_row["brier"]))

    def test_run_benchmark_pipeline_writes_tables_figures_and_summary(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = pathlib.Path(tmpdir)
            predictions = tmp / "predictions.tsv"
            outdir = tmp / "benchmark_run"
            pd.DataFrame(
                {
                    "label": [1, 0, 1, 0],
                    "score": [0.9, 0.1, 0.8, 0.2],
                    "tf": ["CTCF", "CTCF", "IRF1", "IRF1"],
                    "cell": ["K562", "K562", "K562", "K562"],
                    "method": ["fp-tools", "fp-tools", "fp-tools", "fp-tools"],
                }
            ).to_csv(predictions, sep="\t", index=False)

            outputs = run_benchmark_pipeline.run_benchmark_pipeline(
                [predictions],
                outdir,
                bins=5,
                bootstrap=5,
                title="synthetic benchmark",
            )

            self.assertTrue(outputs["combined_predictions"].exists())
            self.assertTrue(outputs["binary_metrics"].exists())
            self.assertTrue(outputs["binary_metrics_bootstrap"].exists())
            self.assertTrue(outputs["run_summary"].exists())
            for key in ("benchmark_figures", "calibration_figures"):
                self.assertEqual(len(outputs[key]), 3)
                self.assertTrue(all(path.exists() for path in outputs[key]))

    def test_compute_and_plot_calibration_reports(self):
        df = pd.DataFrame(
            {
                "label": [1, 0, 1, 0, 1, 0, 1, 0],
                "score": [0.95, 0.10, 0.80, 0.15, 0.70, 0.30, 0.60, 0.40],
                "tf": ["A", "A", "A", "A", "B", "B", "B", "B"],
                "cell": ["K562"] * 8,
                "method": ["model"] * 8,
            }
        )
        bins, summary = compute_calibration.compute_calibration(df, "label", "score", ["tf", "cell", "method"], bins=5)
        self.assertIn("global", set(summary["group"]))
        self.assertTrue((summary["ece"] >= 0).all())
        self.assertEqual(len(bins[bins["group"] == "global"]), 5)

        with tempfile.TemporaryDirectory() as tmpdir:
            out_prefix = pathlib.Path(tmpdir) / "figure_calibration"
            outputs = plot_calibration_panels.plot_calibration(bins, summary, out_prefix)
            self.assertEqual({path.suffix for path in outputs}, {".pdf", ".svg", ".png"})
            for output in outputs:
                self.assertTrue(output.exists())

    def test_plot_multiscale_npz_writes_all_formats(self):
        records = [
            (("chr1", 0, 5), {8: np.array([1, 2, 3, 2, 1]), 16: np.array([0, 1, 2, 1, 0])}),
            (("chr1", 10, 15), {8: np.array([2, 3, 4, 3, 2]), 16: np.array([1, 2, 3, 2, 1])}),
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            npz_path = pathlib.Path(tmpdir) / "multiscale.npz"
            out_prefix = pathlib.Path(tmpdir) / "figure_multiscale"
            write_multiscale_npz(str(npz_path), records, [8, 16], "max")
            outputs = plot_multiscale_npz.plot_multiscale_npz(npz_path, out_prefix)
            self.assertEqual([path.suffix for path in outputs], [".pdf", ".svg", ".png"])
            for output in outputs:
                self.assertTrue(output.exists())
                self.assertGreater(output.stat().st_size, 0)

    def test_plot_benchmark_panels_writes_all_formats(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out_prefix = pathlib.Path(tmpdir) / "figure_benchmark_summary"
            metrics = pd.DataFrame(
                {
                    "group": ["global", "A/C/fp", "B/C/fp"],
                    "n": [6, 4, 2],
                    "positives": [3, 2, 1],
                    "auroc": [1.0, 1.0, 1.0],
                    "auprc": [1.0, 1.0, 1.0],
                    "recall_at_1pct_fdr": [1.0, 1.0, 1.0],
                    "recall_at_5pct_fdr": [1.0, 1.0, 1.0],
                    "recall_at_10pct_fdr": [1.0, 1.0, 1.0],
                    "brier": [0.05, 0.04, 0.06],
                }
            )
            outputs = plot_benchmark_panels.plot_metrics(metrics, out_prefix)
            self.assertEqual({path.suffix for path in outputs}, {".pdf", ".svg", ".png"})
            for path in outputs:
                self.assertTrue(path.exists())


if __name__ == "__main__":
    unittest.main()
