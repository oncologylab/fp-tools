import os
import pathlib
import subprocess
import tempfile
import unittest

import numpy as np
import pandas as pd
from fp_tools.utils import bigwig as pyBigWig
from fp_tools.utils.signals import (
    add_bias_prediction_window,
    atac_correct_arrays,
    fast_rolling_math,
    finalize_bias_prediction,
    footprint_score_array,
    footprint_score_array_fast,
    local_maxima_indices,
)


ROOT = pathlib.Path(__file__).resolve().parents[1]
BIN = ROOT / ".venv" / "bin"


def max_cores() -> str:
    return str(max(1, os.cpu_count() or 1))


def bigwig_window_summary(path, chrom="chr4", start=74000, end=75000):
    bw = pyBigWig.open(str(path))
    try:
        chroms = bw.chroms()
        intervals = bw.intervals(chrom, start, end) or []
        total = round(sum(float(item[2]) for item in intervals), 6)
        mean = round(total / len(intervals), 6) if intervals else 0.0
        return chroms, len(intervals), total, mean
    finally:
        bw.close()

def run_command(command, timeout=90, env=None):
    merged_env = os.environ.copy()
    if env:
        merged_env.update(env)
    result = subprocess.run(
        [str(item) for item in command],
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=merged_env,
        timeout=timeout,
        check=False,
    )
    if result.returncode != 0:
        raise AssertionError(f"Command failed: {' '.join(map(str, command))}\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}")
    return result


class CliGoldenRegressionTest(unittest.TestCase):
    def test_fast_footprint_kernel_matches_legacy_kernel(self):
        rng = np.random.default_rng(7)
        signal = rng.normal(loc=0.0, scale=2.0, size=500).astype("float64")
        signal[::11] = 0.0
        legacy = footprint_score_array(signal, 10, 30, 20, 50)
        fast = footprint_score_array_fast(signal, 10, 30, 20, 50)
        np.testing.assert_allclose(fast, legacy, rtol=1e-6, atol=1e-6)

    def test_atac_correct_cython_kernel_matches_numpy_path(self):
        rng = np.random.default_rng(13)
        length = 421
        window = 100
        k_flank = 3
        reg_end = length - k_flank
        overlaps = int(window / 10)
        prediction_matrix = np.zeros((overlaps, length), dtype="float64")
        prediction_matrix[:, k_flank:reg_end] = np.nan
        prediction_sum = np.zeros(length, dtype="float64")
        prediction_count = np.zeros(length, dtype=np.int64)
        row = 0
        for start in range(k_flank, reg_end - window, 10):
            end = min(start + window, reg_end)
            prediction = rng.random(end - start).astype("float64")
            prediction_matrix[row, start:end] = prediction
            add_bias_prediction_window(prediction_sum, prediction_count, prediction, start, end)
            row = row + 1 if row < overlaps - 1 else 0

        with np.errstate(invalid="ignore"):
            expected_bias = np.nanmean(prediction_matrix, axis=0)
        observed_bias = finalize_bias_prediction(prediction_sum, prediction_count, k_flank, reg_end)
        np.testing.assert_allclose(observed_bias, expected_bias, rtol=1e-12, atol=1e-12, equal_nan=True)

        uncorrected = rng.normal(loc=0.2, scale=2.0, size=length).astype("float64")
        uncorrected[::17] = 0.0
        correction_factor = 1.37
        observed_uncorrected, observed_expected, observed_corrected = atac_correct_arrays(
            uncorrected.copy(),
            expected_bias.astype("float64", copy=True),
            window,
            correction_factor,
        )

        signal_sum = fast_rolling_math(uncorrected.copy(), window, "sum")
        signal_sum[np.isnan(signal_sum)] = 0
        bias_sum = fast_rolling_math(expected_bias.copy(), window, "sum")
        nulls = np.logical_or(np.isclose(bias_sum, 0), np.isnan(bias_sum))
        bias_sum[nulls] = 1
        bias_probas = expected_bias / bias_sum
        bias_probas[nulls] = 0
        expected = signal_sum * bias_probas
        expected_uncorrected = uncorrected.copy() * correction_factor
        expected *= correction_factor
        expected_corrected = expected_uncorrected - expected
        uncorrected_sum = fast_rolling_math(expected_uncorrected, window, "sum")
        uncorrected_sum[np.isnan(uncorrected_sum)] = 0
        corrected_sum = fast_rolling_math(np.abs(expected_corrected), window, "sum")
        corrected_sum[np.isnan(corrected_sum)] = 0
        corrected_pos = np.copy(expected_corrected)
        corrected_pos[corrected_pos < 0] = 0
        corrected_pos_sum = fast_rolling_math(corrected_pos, window, "sum")
        corrected_pos_sum[np.isnan(corrected_pos_sum)] = 0
        corrected_neg_sum = corrected_sum - corrected_pos_sum
        zero_sum = corrected_pos_sum == 0
        corrected_pos_sum[zero_sum] = np.nan
        scale_factor = (uncorrected_sum - corrected_neg_sum) / corrected_pos_sum
        scale_factor[zero_sum] = 1
        scale_factor[scale_factor < 1] = 1
        pos_bool = expected_corrected > 0
        expected_corrected[pos_bool] *= scale_factor[pos_bool]

        np.testing.assert_allclose(observed_uncorrected, expected_uncorrected, rtol=1e-12, atol=1e-12, equal_nan=True)
        np.testing.assert_allclose(observed_expected, expected, rtol=1e-12, atol=1e-12, equal_nan=True)
        np.testing.assert_allclose(observed_corrected, expected_corrected, rtol=1e-12, atol=1e-12, equal_nan=True)

    def test_local_maxima_cython_kernel_matches_numpy_path(self):
        values = np.array([0.0, 1.0, 1.0, 0.0, -np.inf, 2.0, 1.0, -np.inf, 0.0], dtype="float64")
        left = np.empty_like(values)
        right = np.empty_like(values)
        left[0] = -np.inf
        left[1:] = values[:-1]
        right[-1] = -np.inf
        right[:-1] = values[1:]
        mask = np.isfinite(values) & (values >= left) & (values >= right) & ((values > left) | (values > right))
        self.assertEqual(local_maxima_indices(values), np.flatnonzero(mask).tolist())

    def test_footprint_scores_sum_bigwig_summary_is_stable(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output = pathlib.Path(tmpdir) / "footprints_sum.bw"
            run_command(
                [
                    BIN / "call-footprints",
                    "--signal",
                    "test_data/Bcell_corrected.bw",
                    "--regions",
                    "test_data/merged_peaks.bed",
                    "--output",
                    output,
                    "--score",
                    "sum",
                    "--window",
                    "20",
                    "--cores",
                    max_cores(),
                    "--verbosity",
                    "1",
                ]
            )

            self.assertTrue(output.exists())
            self.assertGreater(output.stat().st_size, 0)
            chroms, count, total, mean = bigwig_window_summary(output)

            self.assertEqual(chroms, {"chr4": 190214555})
            self.assertEqual(count, 186)
            self.assertAlmostEqual(total, 106.694972, delta=1e-4)
            self.assertAlmostEqual(mean, 0.573629, places=6)

    def test_footprint_scores_sum_is_stable_across_core_counts(self):
        summaries = []
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = pathlib.Path(tmpdir)
            for cores in ("1", "2"):
                output = tmp / f"footprints_sum_cores_{cores}.bw"
                run_command(
                    [
                        BIN / "call-footprints",
                        "--signal",
                        "test_data/Bcell_corrected.bw",
                        "--regions",
                        "test_data/merged_peaks.bed",
                        "--output",
                        output,
                        "--score",
                        "sum",
                        "--window",
                        "20",
                        "--cores",
                        cores,
                        "--verbosity",
                        "1",
                    ]
                )
                self.assertTrue(output.exists())
                self.assertGreater(output.stat().st_size, 0)
                summaries.append(bigwig_window_summary(output))

        self.assertEqual(summaries[0], summaries[1])
        self.assertEqual(summaries[0][0], {"chr4": 190214555})
        self.assertGreater(summaries[0][1], 0)
        self.assertGreater(summaries[0][2], 0.0)

    def test_batch_footprint_scoring_with_nested_workers_is_stable(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = pathlib.Path(tmpdir)
            out_a = tmp / "Bcell_batch_sum.bw"
            out_b = tmp / "Tcell_batch_sum.bw"
            run_command(
                [
                    BIN / "call-footprints",
                    "--signals",
                    "test_data/Bcell_corrected.bw",
                    "test_data/Tcell_corrected.bw",
                    "--outputs",
                    out_a,
                    out_b,
                    "--regions",
                    "test_data/merged_peaks.bed",
                    "--score",
                    "sum",
                    "--window",
                    "20",
                    "--cores",
                    "16",
                    "--verbosity",
                    "1",
                ]
            )
            for output in (out_a, out_b):
                self.assertTrue(output.exists())
                self.assertGreater(output.stat().st_size, 0)
                chroms, count, total, _ = bigwig_window_summary(output)
                self.assertEqual(chroms, {"chr4": 190214555})
                self.assertGreater(count, 0)
                self.assertGreater(total, 0.0)

    def test_plot_aggregate_text_summary_is_stable(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = pathlib.Path(tmpdir)
            output = tmp / "aggregate.pdf"
            output_txt = tmp / "aggregate.txt"
            run_command(
                [
                    BIN / "plot-aggregate",
                    "--TFBS",
                    "test_data/IRF1_all.bed",
                    "--signals",
                    "test_data/Bcell_footprints.bw",
                    "--output",
                    output,
                    "--output-txt",
                    output_txt,
                    "--flank",
                    "20",
                    "--verbosity",
                    "1",
                ]
            )
            self.assertTrue(output_txt.exists())
            rows = output_txt.read_text(encoding="utf-8").splitlines()
            values = [float(value) for value in rows[2].split("\t")[2].split(",")]

        self.assertEqual(rows[0], "### AGGREGATE")
        self.assertEqual(len(values), 40)
        self.assertEqual(round(sum(values[:5]), 4), 53.5609)
        self.assertEqual(round(values[-1], 4), 10.3948)

    def test_diff_footprints_one_motif_summary_is_stable(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = pathlib.Path(tmpdir) / "diff_footprints"
            run_command(
                [
                    BIN / "diff-footprints",
                    "--signals",
                    "test_data/Bcell_footprints.bw",
                    "test_data/Tcell_footprints.bw",
                    "--motifs",
                    "test_data/individual_motifs/MA0050.2.jaspar",
                    "--genome",
                    "test_data/genome.fa.gz",
                    "--peaks",
                    "test_data/merged_peaks.bed",
                    "--cond-names",
                    "Bcell",
                    "Tcell",
                    "--outdir",
                    outdir,
                    "--prefix",
                    "diff_footprints_probe",
                    "--cores",
                    max_cores(),
                    "--skip-excel",
                    "--verbosity",
                    "1",
                ],
                timeout=120,
            )
            results = pd.read_csv(outdir / "diff_footprints_probe_results.txt", sep="\t")

        self.assertEqual(len(results), 1)
        row = results.iloc[0]
        self.assertEqual(row["name"], "IRF1")
        self.assertEqual(row["motif_id"], "MA0050.2")
        self.assertEqual(int(row["total_tfbs"]), 3269)
        self.assertEqual(int(row["Bcell_bound"]), 1367)
        self.assertEqual(int(row["Tcell_bound"]), 672)
        self.assertAlmostEqual(float(row["Bcell_mean_score"]), 10.58374, places=5)
        self.assertAlmostEqual(float(row["Tcell_mean_score"]), 7.52570, places=5)
        self.assertAlmostEqual(float(row["Bcell_Tcell_change"]), 0.35168, places=5)
        self.assertIn("Bcell_Tcell_qvalue_bh", results.columns)
        self.assertIn("Bcell_Tcell_significant_fdr05", results.columns)
        for column in (
            "Bcell_score_sd",
            "Tcell_score_sd",
            "Bcell_Tcell_delta_fp_se",
            "Bcell_Tcell_log2fc_se",
        ):
            self.assertNotIn(column, results.columns)
        self.assertGreaterEqual(float(row["Bcell_Tcell_qvalue_bh"]), 0.0)
        self.assertLessEqual(float(row["Bcell_Tcell_qvalue_bh"]), 1.0)

    def test_diff_footprints_summary_mode_skips_per_motif_exports(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = pathlib.Path(tmpdir) / "diff_footprints_summary"
            completed = run_command(
                [
                    BIN / "diff-footprints",
                    "--signals",
                    "test_data/Bcell_footprints.bw",
                    "test_data/Tcell_footprints.bw",
                    "--motifs",
                    "test_data/individual_motifs/MA0050.2.jaspar",
                    "--genome",
                    "test_data/genome.fa.gz",
                    "--peaks",
                    "test_data/merged_peaks.bed",
                    "--cond-names",
                    "Bcell",
                    "Tcell",
                    "--outdir",
                    outdir,
                    "--prefix",
                    "diff_footprints_probe",
                    "--cores",
                    max_cores(),
                    "--skip-excel",
                    "--plot-aggregate",
                    "off",
                    "--motif-outputs",
                    "summary",
                    "--verbosity",
                    "2",
                ],
                timeout=120,
            )
            results = pd.read_csv(outdir / "diff_footprints_probe_results.txt", sep="\t")

        self.assertIn("diff-footprints (run started", completed.stdout)
        self.assertIn("Creating diff-footprints plot(s)", completed.stdout)
        self.assertIn("Finished diff-footprints run", completed.stdout)
        self.assertEqual(len(results), 1)
        self.assertEqual(int(results.iloc[0]["total_tfbs"]), 3269)
        self.assertAlmostEqual(float(results.iloc[0]["Bcell_Tcell_change"]), 0.35168, places=5)
        self.assertFalse((outdir / "IRF1_MA0050.2" / "beds" / "IRF1_MA0050.2_all.bed").exists())
        self.assertFalse((outdir / "IRF1_MA0050.2" / "IRF1_MA0050.2_overview.txt").exists())

    def test_match_motifs_logs_use_public_command_name(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = pathlib.Path(tmpdir) / "match_motifs_logging"
            completed = run_command(
                [
                    BIN / "match-motifs",
                    "--signals",
                    "test_data/Bcell_footprints.bw",
                    "test_data/Tcell_footprints.bw",
                    "--motifs",
                    "test_data/individual_motifs/MA0050.2.jaspar",
                    "--genome",
                    "test_data/genome.fa.gz",
                    "--peaks",
                    "test_data/merged_peaks.bed",
                    "--sample-names",
                    "Bcell",
                    "Tcell",
                    "--cond-names",
                    "Bcell",
                    "Tcell",
                    "--outdir",
                    outdir,
                    "--cores",
                    max_cores(),
                    "--skip-excel",
                    "--plot-aggregate",
                    "off",
                    "--motif-outputs",
                    "summary",
                    "--verbosity",
                    "2",
                ],
                timeout=120,
            )

        self.assertIn("match-motifs (run started", completed.stdout)
        self.assertIn("Creating match-motifs summary output(s)", completed.stdout)
        self.assertIn("Finished match-motifs run", completed.stdout)
        self.assertNotIn("diff-footprints", completed.stdout)

    def test_match_motifs_auto_writes_compact_cache_and_per_motif_beds(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = pathlib.Path(tmpdir) / "match_motifs_auto"
            run_command(
                [
                    BIN / "match-motifs",
                    "--signals",
                    "test_data/Bcell_footprints.bw",
                    "--motifs",
                    "test_data/individual_motifs/MA0050.2.jaspar",
                    "--genome",
                    "test_data/genome.fa.gz",
                    "--peaks",
                    "test_data/merged_peaks.bed",
                    "--sample-names",
                    "Bcell",
                    "--outdir",
                    outdir,
                    "--prefix",
                    "motif_matches",
                    "--cores",
                    max_cores(),
                    "--skip-excel",
                    "--verbosity",
                    "1",
                ],
                timeout=120,
            )
            results = pd.read_csv(outdir / "motif_matches_results.txt", sep="\t")
            self.assertEqual(len(results), 1)
            self.assertEqual(int(results.iloc[0]["total_tfbs"]), 3269)
            self.assertEqual(int(results.iloc[0]["Bcell_bound"]), 891)
            self.assertAlmostEqual(float(results.iloc[0]["Bcell_mean_score"]), 10.58374, places=5)
            self.assertTrue((outdir / "cache" / "motif_sites.tsv.gz").exists())
            self.assertFalse((outdir / "cache" / "motif_sites.zip").exists())
            self.assertTrue((outdir / "cache" / "background_scores.tsv.gz").exists())
            bed_dir = outdir / "IRF1_MA0050.2" / "beds"
            self.assertTrue((bed_dir / "IRF1_MA0050.2_all.bed").exists())
            self.assertTrue((bed_dir / "IRF1_MA0050.2_Bcell_bound.bed").exists())
            self.assertTrue((bed_dir / "IRF1_MA0050.2_Bcell_unbound.bed").exists())
            self.assertTrue((bed_dir / ".done").exists())
            self.assertFalse((outdir / "IRF1_MA0050.2" / "IRF1_MA0050.2_overview.txt").exists())

    def test_match_motifs_shared_project_scan_matches_per_sample_scan(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = pathlib.Path(tmpdir)
            projects = {"shared": root / "shared", "per": root / "per"}
            for project in projects.values():
                (project / "peaks").mkdir(parents=True)
                (project / "peaks" / "merged_peaks_filtered.bed").write_text(
                    (ROOT / "test_data" / "merged_peaks.bed").read_text(encoding="utf-8"),
                    encoding="utf-8",
                )
                (project / "samples.tsv").write_text("sample\tcondition\nA\tA\nB\tB\n", encoding="utf-8")
                for sample, fixture in [("A", "Bcell_footprints.bw"), ("B", "Tcell_footprints.bw")]:
                    footprint_dir = project / "samples" / sample / "footprints"
                    footprint_dir.mkdir(parents=True)
                    os.symlink(ROOT / "test_data" / fixture, footprint_dir / f"{sample}_footprints.bw")

            base = [
                BIN / "match-motifs",
                "--genome",
                "test_data/genome.fa.gz",
                "--motifs",
                "test_data/individual_motifs/MA0050.2.jaspar",
                "--cores",
                "4",
                "--skip-excel",
                "--verbosity",
                "1",
            ]
            run_command(
                base
                + [
                    "--sample-table",
                    projects["shared"] / "samples.tsv",
                    "--peaks",
                    projects["shared"] / "peaks" / "merged_peaks_filtered.bed",
                    "--outdir",
                    projects["shared"],
                ],
                timeout=120,
                env={"FP_TOOLS_SYNC_MATCH_BEDS": "1"},
            )
            run_command(
                base
                + [
                    "--sample-table",
                    projects["per"] / "samples.tsv",
                    "--peaks",
                    projects["per"] / "peaks" / "merged_peaks_filtered.bed",
                    "--outdir",
                    projects["per"],
                    "--match-scan-mode",
                    "per-sample",
                ],
                timeout=120,
            )

            for sample in ["A", "B"]:
                shared_dir = projects["shared"] / "samples" / sample / "match_motifs"
                per_dir = projects["per"] / "samples" / sample / "match_motifs"
                self.assertEqual(
                    (shared_dir / "motif_matches_results.txt").read_text(encoding="utf-8"),
                    (per_dir / "motif_matches_results.txt").read_text(encoding="utf-8"),
                )
                for name in [
                    "IRF1_MA0050.2_all.bed",
                    f"IRF1_MA0050.2_{sample}_bound.bed",
                    f"IRF1_MA0050.2_{sample}_unbound.bed",
                ]:
                    self.assertEqual(
                        (shared_dir / "IRF1_MA0050.2" / "beds" / name).read_text(encoding="utf-8"),
                        (per_dir / "IRF1_MA0050.2" / "beds" / name).read_text(encoding="utf-8"),
                    )

    def test_match_motifs_shared_project_summary_skips_per_motif_beds(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            project = pathlib.Path(tmpdir) / "summary"
            (project / "peaks").mkdir(parents=True)
            (project / "peaks" / "merged_peaks_filtered.bed").write_text(
                (ROOT / "test_data" / "merged_peaks.bed").read_text(encoding="utf-8"),
                encoding="utf-8",
            )
            (project / "samples.tsv").write_text("sample\tcondition\nA\tA\nB\tB\n", encoding="utf-8")
            for sample, fixture in [("A", "Bcell_footprints.bw"), ("B", "Tcell_footprints.bw")]:
                footprint_dir = project / "samples" / sample / "footprints"
                footprint_dir.mkdir(parents=True)
                os.symlink(ROOT / "test_data" / fixture, footprint_dir / f"{sample}_footprints.bw")

            run_command(
                [
                    BIN / "match-motifs",
                    "--sample-table",
                    project / "samples.tsv",
                    "--genome",
                    "test_data/genome.fa.gz",
                    "--peaks",
                    project / "peaks" / "merged_peaks_filtered.bed",
                    "--motifs",
                    "test_data/individual_motifs/MA0050.2.jaspar",
                    "--outdir",
                    project,
                    "--motif-outputs",
                    "summary",
                    "--cores",
                    "4",
                    "--skip-excel",
                    "--verbosity",
                    "1",
                ],
                timeout=120,
                env={"FP_TOOLS_SYNC_MATCH_BEDS": "1"},
            )

            for sample in ["A", "B"]:
                match_dir = project / "samples" / sample / "match_motifs"
                self.assertTrue((match_dir / "motif_matches_results.txt").is_file())
                self.assertTrue((match_dir / "cache" / "motif_sites.tsv.gz").is_file())
                self.assertTrue((match_dir / "cache" / "background_scores.tsv.gz").is_file())
                self.assertFalse((match_dir / "IRF1_MA0050.2").exists())

    def test_diff_footprints_replicate_grouping_writes_diagnostics(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = pathlib.Path(tmpdir) / "diff_footprints_reps"
            run_command(
                [
                    BIN / "diff-footprints",
                    "--signals",
                    "test_data/demo_Bcell_rep1_footprints.bw",
                    "test_data/demo_Bcell_rep2_footprints.bw",
                    "test_data/demo_Tcell_rep1_footprints.bw",
                    "test_data/demo_Tcell_rep2_footprints.bw",
                    "--motifs",
                    "test_data/individual_motifs/MA0050.2.jaspar",
                    "--genome",
                    "test_data/genome.fa.gz",
                    "--peaks",
                    "test_data/merged_peaks.bed",
                    "--cond-names",
                    "Bcell",
                    "Bcell",
                    "Tcell",
                    "Tcell",
                    "--outdir",
                    outdir,
                    "--prefix",
                    "diff_footprints_probe",
                    "--cores",
                    max_cores(),
                    "--skip-excel",
                    "--verbosity",
                    "1",
                    "--normalization",
                    "sample-quantile",
                    "--replicate-report",
                    "on",
                ],
                timeout=120,
            )
            results = pd.read_csv(outdir / "diff_footprints_probe_results.txt", sep="	")
            report = pd.read_csv(outdir / "diff_footprints_probe_replicate_report.tsv", sep="	")
            replicate_matrix_exists = (
                outdir / "diff_footprints_probe_replicate_motif_score_matrix.tsv"
            ).exists()

        row = results.iloc[0]
        for column in (
            "Bcell_n_replicates",
            "Bcell_score_sd",
            "Tcell_n_replicates",
            "Tcell_score_sd",
            "Bcell_Tcell_mean_delta_fp",
            "Bcell_Tcell_mean_log2fc",
            "Bcell_Tcell_delta_fp_se",
            "Bcell_Tcell_log2fc_se",
            "Bcell_Tcell_qvalue_bh",
            "Bcell_Tcell_significant_fdr05",
        ):
            self.assertIn(column, results.columns)
        self.assertEqual(int(row["Bcell_n_replicates"]), 2)
        self.assertEqual(int(row["Tcell_n_replicates"]), 2)
        self.assertGreater(float(row["Bcell_score_sd"]), 0.0)
        self.assertGreater(float(row["Bcell_mean_score"]), float(row["Tcell_mean_score"]))
        self.assertGreater(float(row["Bcell_Tcell_mean_delta_fp"]), 0.0)
        self.assertGreater(float(row["Bcell_Tcell_mean_log2fc"]), 0.0)
        self.assertTrue((report["replicate_support"] == "replicate-counts-only").all())
        self.assertTrue(replicate_matrix_exists)

    @unittest.skipUnless(os.environ.get("FP_TOOLS_RUN_SLOW_REGRESSIONS") == "1", "slow atac-correct regression is opt-in")
    def test_atacorrect_fixture_smoke_outputs_corrected_bigwig(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = pathlib.Path(tmpdir)
            run_command(
                [
                    BIN / "atac-correct",
                    "--bams",
                    "test_data/Bcell.bam",
                    "--genome",
                    "test_data/genome.fa.gz",
                    "--peaks",
                    "test_data/merged_peaks.bed",
                    "--blacklist",
                    "test_data/blacklist.bed",
                    "--outdir",
                    outdir,
                    "--prefix",
                    "Bcell_ci",
                    "--cores",
                    max_cores(),
                    "--track-off",
                    "bias",
                    "expected",
                    "uncorrected",
                    "--verbosity",
                    "1",
                ],
                timeout=300,
            )
            self.assertTrue((outdir / "Bcell_ci_corrected.bw").exists())


if __name__ == "__main__":
    unittest.main()
