import argparse
import pathlib
import subprocess
import unittest

from fp_tools.gui_config import expand_jobs, normalize_config
from fp_tools.tools.diff_footprints import _resolve_motif_arguments
from fp_tools.utils.motif_databases import (
    DEFAULT_MOTIF_DB,
    MOTIF_DATABASES,
    motif_db_path,
    normalize_motif_db_key,
    resolve_motif_inputs,
)
from fp_tools.utils.motifs import MotifList


ROOT = pathlib.Path(__file__).resolve().parents[1]


class MotifDatabaseTest(unittest.TestCase):
    def test_default_resolves_to_jaspar2026_vertebrates(self):
        paths = resolve_motif_inputs(None, None)
        self.assertEqual(len(paths), 1)
        self.assertTrue(paths[0].endswith("JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt"))
        self.assertEqual(normalize_motif_db_key("jaspar2026"), DEFAULT_MOTIF_DB)

    def test_user_motifs_are_not_augmented_unless_motif_db_is_requested(self):
        self.assertEqual(resolve_motif_inputs(["custom.jaspar"], None), ["custom.jaspar"])
        paths = resolve_motif_inputs(["custom.jaspar"], "hocomoco14")
        self.assertEqual(paths[-1], "custom.jaspar")
        self.assertTrue(paths[0].endswith("H14CORE_jaspar_format.txt"))

    def test_optional_resolver_does_not_add_default_database(self):
        self.assertEqual(resolve_motif_inputs(None, None, use_default=False), [])
        paths = resolve_motif_inputs(None, "jaspar2026", use_default=False)
        self.assertEqual(len(paths), 1)
        self.assertTrue(paths[0].endswith("JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt"))

    def test_all_bundled_motif_files_exist_and_parse(self):
        expected_counts = {
            "jaspar2026_vertebrates": 1019,
            "jaspar2026_core": 2633,
            "hocomoco14_core": 1595,
            "hocomoco14_invivo": 1595,
            "hocomoco14_invitro": 1579,
            "hocomoco14_rsnp": 1595,
        }
        for key, db in MOTIF_DATABASES.items():
            with self.subTest(key=key):
                path = motif_db_path(key)
                self.assertTrue(path.exists(), path)
                motifs = MotifList().from_file(str(path))
                self.assertGreater(len(motifs), 0)
                self.assertEqual(len(motifs[0].counts), 4)
                if key in expected_counts:
                    self.assertEqual(len(motifs), expected_counts[key])
                self.assertIn(db.license.split(";")[0], db.license)

    def test_diff_footprints_resolver_sets_default_motifs_on_args(self):
        args = argparse.Namespace(motifs=None, motif_db=None)
        paths = _resolve_motif_arguments(args)
        self.assertEqual(paths, args.motifs)
        self.assertEqual(len(paths), 1)
        self.assertTrue(paths[0].endswith("JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt"))

    def test_cli_lists_motif_databases_without_required_inputs(self):
        exe = ROOT / ".venv" / "bin" / "diff-footprints"
        if not exe.exists():
            self.skipTest(f"{exe} is not available in this checkout")
        result = subprocess.run(
            [str(exe), "--list-motif-dbs"],
            cwd=ROOT,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=30,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("jaspar2026_vertebrates", result.stdout)
        self.assertIn("hocomoco14_core", result.stdout)

    def test_optional_motif_clis_list_motif_databases_without_required_inputs(self):
        commands = ["sc-footprinting", "discover-motifs"]
        for command in commands:
            exe = ROOT / ".venv" / "bin" / command
            if not exe.exists():
                self.skipTest(f"{exe} is not available in this checkout")
            with self.subTest(command=command):
                result = subprocess.run(
                    [str(exe), "--list-motif-dbs"],
                    cwd=ROOT,
                    text=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    timeout=30,
                    check=False,
                )
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn("jaspar2026_vertebrates", result.stdout)
                self.assertIn("hocomoco14_core", result.stdout)

    def test_config_no_longer_requires_motifs_for_diff_footprints(self):
        config = normalize_config(
            {
                "version": 1,
                "run_mode": "batch",
                "comparisons": [
                    {
                        "comparison_id": "demo",
                        "tool": "diff-footprints",
                        "signals": ["a.bw", "b.bw"],
                        "genome": "genome.fa",
                        "peaks": "peaks.bed",
                        "cond_names": ["A", "B"],
                        "outdir": "out",
                    }
                ],
            }
        )
        jobs = expand_jobs(config)
        self.assertEqual(len(jobs), 1)
        self.assertNotIn("--motifs", jobs[0].command)

    def test_config_can_pass_builtin_motif_db(self):
        config = normalize_config(
            {
                "version": 1,
                "run_mode": "batch",
                "comparisons": [
                    {
                        "comparison_id": "demo",
                        "tool": "diff-footprints",
                        "signals": ["a.bw", "b.bw"],
                        "genome": "genome.fa",
                        "peaks": "peaks.bed",
                        "motif_db": "hocomoco14_core",
                        "cond_names": ["A", "B"],
                        "outdir": "out",
                    }
                ],
            }
        )
        command = expand_jobs(config)[0].command
        self.assertIn("--motif-db", command)
        self.assertIn("hocomoco14_core", command)

    def test_config_can_pass_known_motif_db(self):
        config = normalize_config(
            {
                "version": 1,
                "run_mode": "batch",
                "samples": [
                    {
                        "sample_id": "demo",
                        "tool": "discover-motifs",
                        "fasta": "sites.fa",
                        "outdir": "motifs",
                        "known_motif_db": "jaspar2026_vertebrates",
                    }
                ],
            }
        )
        command = expand_jobs(config)[0].command
        self.assertIn("--known-motif-db", command)
        self.assertIn("jaspar2026_vertebrates", command)


if __name__ == "__main__":
    unittest.main()
