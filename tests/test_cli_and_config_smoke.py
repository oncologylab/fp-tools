import contextlib
import io
import argparse
import pathlib
import subprocess
import tempfile
import unittest
from types import SimpleNamespace

from fp_tools.cli_batch import run_config_file
from fp_tools.gui_config import expand_jobs, load_yaml_config, normalize_config
from fp_tools.parsers import add_aggregate_arguments, add_atacorrect_arguments, add_bindetect_arguments, add_scorebigwig_arguments
from fp_tools.tools import bindetect
from fp_tools.tools.bindetect import _prepare_condition_metadata


ROOT = pathlib.Path(__file__).resolve().parents[1]
CONFIG_DIR = ROOT / "examples" / "gui_configs"


class CliAndConfigSmokeTest(unittest.TestCase):
    def test_all_example_yaml_configs_expand_to_jobs(self):
        config_paths = sorted(CONFIG_DIR.glob("*.yml"))
        self.assertGreaterEqual(len(config_paths), 1)

        for path in config_paths:
            with self.subTest(config=path.name):
                config = normalize_config(load_yaml_config(path))
                jobs = expand_jobs(config)
                self.assertGreaterEqual(len(jobs), 1)
                for job in jobs:
                    self.assertTrue(job.job_id)
                    self.assertIn(
                        job.tool,
                        {
                            "atac-correct",
                            "call-footprints",
                            "diff-footprints",
                            "normalize-bigwig",
                            "plot-aggregate",
                            "plot-aggregate-batch",
                            "motif-discovery",
                            "motif-summary",
                            "fp-tools-score-variants",
                            "pseudobulk-fragments",
                            "find-signature-fp",
                            "pseudobulk-footprints",
                        },
                    )
                    self.assertEqual(job.command[0], job.tool)

    def test_all_example_yaml_configs_support_dry_run(self):
        for path in sorted(CONFIG_DIR.glob("*.yml")):
            with self.subTest(config=path.name):
                stdout = io.StringIO()
                with contextlib.redirect_stdout(stdout):
                    code = run_config_file(path, dry_run=True)
                self.assertEqual(code, 0)
                self.assertIn("[", stdout.getvalue())

    def test_packaged_entry_points_print_help(self):
        commands = [
            "atac-correct",
            "call-footprints",
            "match-motifs",
            "diff-footprints",
            "normalize-bigwig",
            "plot-aggregate",
            "plot-aggregate-batch",
            "run-workflow",
            "motif-discovery",
            "motif-summary",
            "fp-tools-score-variants",
            "pseudobulk-fragments",
            "find-signature-fp",
            "ATACorrect",
            "FootprintScores",
            "ScoreBigwig",
            "BINDetect",
            "PlotAggregate",
        ]
        for command in commands:
            exe = ROOT / ".venv" / "bin" / command
            if not exe.exists():
                self.skipTest(f"{exe} is not available in this checkout")
            with self.subTest(command=command):
                result = subprocess.run(
                    [str(exe), "--help"],
                    cwd=ROOT,
                    text=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    timeout=30,
                    check=False,
                )
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn("usage:", result.stdout.lower())

    def test_call_footprints_accepts_batch_signal_arguments(self):
        parser = add_scorebigwig_arguments(argparse.ArgumentParser())
        args = parser.parse_args(
            [
                "--signals",
                "A_corrected.bw",
                "B_corrected.bw",
                "--outputs",
                "A_footprints.bw",
                "B_footprints.bw",
                "--regions",
                "peaks.bed",
                "--output-beds",
                "A_candidates.bed",
                "B_candidates.bed",
            ]
        )
        self.assertEqual(args.signals, ["A_corrected.bw", "B_corrected.bw"])
        self.assertEqual(args.outputs, ["A_footprints.bw", "B_footprints.bw"])
        self.assertEqual(args.output_beds, ["A_candidates.bed", "B_candidates.bed"])

    def test_atac_correct_accepts_batch_bam_and_peak_arguments(self):
        parser = add_atacorrect_arguments(argparse.ArgumentParser())
        args = parser.parse_args(
            [
                "--bams",
                "A.bam",
                "B.bam",
                "--sample-names",
                "A",
                "B",
                "--genome",
                "genome.fa",
                "--peaks",
                "A_peaks.bed",
                "B_peaks.bed",
            ]
        )
        self.assertEqual(args.bams, ["A.bam", "B.bam"])
        self.assertEqual(args.sample_names, ["A", "B"])
        self.assertEqual(args.peaks, ["A_peaks.bed", "B_peaks.bed"])

    def test_plot_aggregate_accepts_match_dir_html_arguments(self):
        parser = add_aggregate_arguments(argparse.ArgumentParser())
        args = parser.parse_args(
            [
                "--match-dir",
                "results/motif_matches/sample",
                "--signals",
                "sample_corrected.bw",
                "--motifs",
                "SPIB",
                "CEBPB",
                "--format",
                "html",
                "--output",
                "aggregate.html",
            ]
        )
        self.assertEqual(args.match_dir, ["results/motif_matches/sample"])
        self.assertEqual(args.motifs, ["SPIB", "CEBPB"])
        self.assertEqual(args.format, "html")

    def test_diff_footprints_accepts_sample_names_separate_from_conditions(self):
        parser = add_bindetect_arguments(argparse.ArgumentParser(prog="diff-footprints"))
        args = parser.parse_args(
            [
                "--signals",
                "K562_rep1_footprints.bw",
                "K562_rep2_footprints.bw",
                "HepG2_rep1_footprints.bw",
                "HepG2_rep2_footprints.bw",
                "--sample-names",
                "K562_R1",
                "K562_R2",
                "HepG2_R1",
                "HepG2_R2",
                "--cond-names",
                "K562",
                "K562",
                "HepG2",
                "HepG2",
                "--genome",
                "hg38.fa.gz",
                "--peaks",
                "peaks.bed",
            ]
        )
        prepared = _prepare_condition_metadata(args)
        self.assertEqual(prepared.cond_names, ["K562", "HepG2"])
        self.assertEqual(prepared.sample_names, ["K562_R1", "K562_R2", "HepG2_R1", "HepG2_R2"])
        self.assertEqual(prepared.condition_samples["K562"], ["K562_R1", "K562_R2"])
        self.assertEqual(prepared.condition_samples["HepG2"], ["HepG2_R1", "HepG2_R2"])
        self.assertEqual(prepared.sample_to_condition["K562_R2"], "K562")

    def test_diff_footprints_accepts_folder_inputs(self):
        parser = add_bindetect_arguments(argparse.ArgumentParser(prog="diff-footprints"))
        args = parser.parse_args(
            [
                "--sample-dirs",
                "results/samples/A",
                "results/samples/B",
                "--cond-names",
                "A",
                "B",
                "--genome",
                "hg38.fa.gz",
                "--peaks",
                "peaks.bed",
            ]
        )
        self.assertEqual(args.sample_dirs, ["results/samples/A", "results/samples/B"])
        self.assertIsNone(args.project_dir)

    def test_diff_footprints_disambiguates_sample_condition_name_collisions(self):
        parser = add_bindetect_arguments(argparse.ArgumentParser(prog="diff-footprints"))
        args = parser.parse_args(
            [
                "--signals",
                "Bcell_footprints.bw",
                "Tcell_footprints.bw",
                "--sample-names",
                "Bcell",
                "Tcell",
                "--cond-names",
                "Bcell",
                "Tcell",
                "--genome",
                "hg38.fa.gz",
                "--peaks",
                "peaks.bed",
            ]
        )
        prepared = _prepare_condition_metadata(args)
        self.assertEqual(prepared.sample_names, ["Bcell_sample", "Tcell_sample"])
        self.assertEqual(prepared.condition_samples["Bcell"], ["Bcell_sample"])
        self.assertEqual(prepared.condition_samples["Tcell"], ["Tcell_sample"])

    def test_cached_match_dirs_merge_to_tfbs_tmp_files(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = pathlib.Path(tmpdir)
            for sample, score in [("A", "1.25000"), ("B", "2.50000")]:
                bed_dir = root / sample / "match_motifs" / "TF1_MA0001.1" / "beds"
                bed_dir.mkdir(parents=True)
                (bed_dir / "TF1_MA0001.1_all.bed").write_text(
                    "chr1\t10\t20\tTF1\t8.0\t+\tchr1\t1\t100\t" + score + "\n",
                    encoding="utf-8",
                )
            outdir = root / "out"
            (outdir / "TF1_MA0001.1" / "beds").mkdir(parents=True)
            args = SimpleNamespace(
                cached_match_dirs=[
                    str(root / "A" / "match_motifs"),
                    str(root / "B" / "match_motifs"),
                ],
                sample_names=["A_sample", "B_sample"],
                outdir=str(outdir),
                peak_header_list=["peak_chr", "peak_start", "peak_end"],
            )
            logger = bindetect.FpToolsLogger("test", 0)
            bindetect._write_cached_tfbs_tmp_files(args, ["TF1_MA0001.1"], logger)
            tmp_file = outdir / "TF1_MA0001.1" / "beds" / "TF1_MA0001.1.tmp"
            self.assertEqual(
                tmp_file.read_text(encoding="utf-8").strip(),
                "chr1\t10\t20\tTF1\t8.0\t+\tchr1\t1\t100\t1.25000\t2.50000",
            )

    def test_match_motifs_sample_names_label_match_only_outputs(self):
        parser = add_bindetect_arguments(argparse.ArgumentParser(prog="match-motifs"), command_name="match-motifs")
        args = parser.parse_args(
            [
                "--signals",
                "file_stem_A.bw",
                "file_stem_B.bw",
                "--sample-names",
                "SampleA",
                "SampleB",
                "--genome",
                "hg38.fa.gz",
                "--peaks",
                "peaks.bed",
            ]
        )
        args.match_only = True
        prepared = _prepare_condition_metadata(args)
        self.assertEqual(prepared.cond_names, ["SampleA", "SampleB"])
        self.assertEqual(prepared.sample_names, ["SampleA", "SampleB"])
        self.assertEqual(prepared.condition_samples["SampleA"], ["SampleA"])


if __name__ == "__main__":
    unittest.main()
