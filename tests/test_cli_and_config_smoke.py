import contextlib
import io
import pathlib
import subprocess
import unittest

from fp_tools.cli_batch import run_config_file
from fp_tools.gui_config import expand_jobs, load_yaml_config, normalize_config


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
                    self.assertIn(job.tool, {"ATACorrect", "FootprintScores", "BINDetect", "PlotAggregate"})
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
            "ATACorrect",
            "FootprintScores",
            "BINDetect",
            "PlotAggregate",
            "fp-tools-run",
            "fp-tools-gui",
            "fp-tools-build-tfbs-features",
            "fp-tools-train-tfbs-model",
            "fp-tools-predict-tfbs",
            "fp-tools-generate-candidates",
            "fp-tools-rerank-candidates",
            "fp-tools-export-candidate-fasta",
            "fp-tools-meme-command",
            "fp-tools-motif-discovery",
            "fp-tools-summarize-motifs",
            "fp-tools-score-variants",
            "fp-tools-pseudobulk",
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


if __name__ == "__main__":
    unittest.main()
