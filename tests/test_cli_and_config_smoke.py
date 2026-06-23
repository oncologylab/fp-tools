import contextlib
import io
import argparse
import pathlib
import subprocess
import unittest

from fp_tools.cli_batch import run_config_file
from fp_tools.gui_config import expand_jobs, load_yaml_config, normalize_config
from fp_tools.parsers import add_aggregate_arguments, add_scorebigwig_arguments


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


if __name__ == "__main__":
    unittest.main()
