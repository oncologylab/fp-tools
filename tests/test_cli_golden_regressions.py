import os
import pathlib
import subprocess
import tempfile
import unittest

import pandas as pd
import pyBigWig


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

def run_command(command, timeout=90):
    result = subprocess.run(
        [str(item) for item in command],
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=timeout,
        check=False,
    )
    if result.returncode != 0:
        raise AssertionError(f"Command failed: {' '.join(map(str, command))}\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}")
    return result


class CliGoldenRegressionTest(unittest.TestCase):
    def test_footprint_scores_sum_bigwig_summary_is_stable(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output = pathlib.Path(tmpdir) / "footprints_sum.bw"
            run_command(
                [
                    BIN / "FootprintScores",
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
        self.assertEqual((count, total, mean), (186, 106.694972, 0.573629))

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

    def test_plot_aggregate_text_summary_is_stable(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = pathlib.Path(tmpdir)
            output = tmp / "aggregate.pdf"
            output_txt = tmp / "aggregate.txt"
            run_command(
                [
                    BIN / "PlotAggregate",
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

    def test_bindetect_one_motif_summary_is_stable(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = pathlib.Path(tmpdir) / "bindetect"
            run_command(
                [
                    BIN / "BINDetect",
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
                    "bindetect_probe",
                    "--cores",
                    max_cores(),
                    "--skip-excel",
                    "--verbosity",
                    "1",
                ],
                timeout=120,
            )
            results = pd.read_csv(outdir / "bindetect_probe_results.txt", sep="\t")

        self.assertEqual(len(results), 1)
        row = results.iloc[0]
        self.assertEqual(row["name"], "IRF1")
        self.assertEqual(row["motif_id"], "MA0050.2")
        self.assertEqual(int(row["total_tfbs"]), 3269)
        self.assertEqual(int(row["Bcell_bound"]), 1099)
        self.assertEqual(int(row["Tcell_bound"]), 620)
        self.assertAlmostEqual(float(row["Bcell_mean_score"]), 10.68326, places=5)
        self.assertAlmostEqual(float(row["Tcell_mean_score"]), 7.58304, places=5)
        self.assertAlmostEqual(float(row["Bcell_Tcell_change"]), 0.34019, places=5)
        self.assertIn("Bcell_Tcell_qvalue_bh", results.columns)
        self.assertIn("Bcell_Tcell_significant_fdr05", results.columns)
        self.assertGreaterEqual(float(row["Bcell_Tcell_qvalue_bh"]), 0.0)
        self.assertLessEqual(float(row["Bcell_Tcell_qvalue_bh"]), 1.0)

    def test_bindetect_replicate_grouping_writes_diagnostics(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = pathlib.Path(tmpdir) / "bindetect_reps"
            run_command(
                [
                    BIN / "BINDetect",
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
                    "bindetect_probe",
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
            results = pd.read_csv(outdir / "bindetect_probe_results.txt", sep="	")
            report = pd.read_csv(outdir / "bindetect_probe_replicate_report.tsv", sep="	")

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
        self.assertTrue((report["replicate_support"] == "replicate-supported").all())

    @unittest.skipUnless(os.environ.get("FP_TOOLS_RUN_SLOW_REGRESSIONS") == "1", "slow ATACorrect regression is opt-in")
    def test_atacorrect_fixture_smoke_outputs_corrected_bigwig(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = pathlib.Path(tmpdir)
            run_command(
                [
                    BIN / "ATACorrect",
                    "--bam",
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
