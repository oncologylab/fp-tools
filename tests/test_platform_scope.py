from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest import mock

from fp_tools.tools import bulk_footprinting, prepare_atac


class PlatformScopeTest(unittest.TestCase):
    def test_prepare_atac_rejects_non_linux_before_writes_or_runtime_setup(self):
        for host in ("Darwin", "Windows"):
            with self.subTest(host=host), tempfile.TemporaryDirectory() as tmpdir:
                output = Path(tmpdir) / "project"
                config = Path(tmpdir) / "defaults.yml"
                with mock.patch(
                    "fp_tools.platform_support.platform.system", return_value=host
                ), mock.patch.object(prepare_atac, "prepare_command_runtime") as runtime:
                    with self.assertRaises(SystemExit) as raised:
                        prepare_atac.main(["--write-default-config", str(config)])
                self.assertEqual(raised.exception.code, 2)
                self.assertFalse(output.exists())
                self.assertFalse(config.exists())
                runtime.assert_not_called()

    def test_prepare_atac_help_remains_available_on_non_linux(self):
        with mock.patch(
            "fp_tools.platform_support.platform.system", return_value="Darwin"
        ), self.assertRaises(SystemExit) as raised:
            prepare_atac.main(["--help"])
        self.assertEqual(raised.exception.code, 0)

    def test_bulk_raw_mode_rejects_non_linux_for_every_runtime(self):
        for host in ("Darwin", "Windows"):
            for runtime_mode in ("auto", "managed", "system", "container"):
                with self.subTest(host=host, runtime=runtime_mode), tempfile.TemporaryDirectory() as tmpdir:
                    output = Path(tmpdir) / "project"
                    with mock.patch(
                        "fp_tools.platform_support.platform.system", return_value=host
                    ), mock.patch.object(
                        bulk_footprinting, "prepare_command_runtime"
                    ) as runtime:
                        with self.assertRaises(SystemExit) as raised:
                            bulk_footprinting.main(
                                [
                                    "--reads-table",
                                    "reads.tsv",
                                    "--comparison-table",
                                    "comparisons.tsv",
                                    "--genome",
                                    "hg38",
                                    "--outdir",
                                    str(output),
                                    "--runtime",
                                    runtime_mode,
                                ]
                            )
                    self.assertEqual(raised.exception.code, 2)
                    self.assertFalse(output.exists())
                    runtime.assert_not_called()

    def test_linux_bulk_raw_mode_reaches_the_existing_workflow(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            reads = root / "reads.tsv"
            comparisons = root / "comparisons.tsv"
            reads.write_text("sample\tcondition\tfastq_1\n", encoding="utf-8")
            comparisons.write_text("comparison\tcond1\tcond2\n", encoding="utf-8")
            with mock.patch(
                "fp_tools.platform_support.platform.system", return_value="Linux"
            ), mock.patch.object(
                bulk_footprinting, "run_bulk_footprinting", return_value=0
            ) as run:
                result = bulk_footprinting.main(
                    [
                        "--reads-table",
                        str(reads),
                        "--comparison-table",
                        str(comparisons),
                        "--genome",
                        "hg38",
                        "--outdir",
                        str(root / "project"),
                        "--dry-run",
                    ]
                )
        self.assertEqual(result, 0)
        run.assert_called_once()

    def test_bam_mode_rejects_raw_only_options(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "project"
            with self.assertRaises(SystemExit) as raised:
                bulk_footprinting.main(
                    [
                        "--sample-table",
                        "samples.tsv",
                        "--comparison-table",
                        "comparisons.tsv",
                        "--genome",
                        "genome.fa.gz",
                        "--outdir",
                        str(output),
                        "--profile",
                        "modern",
                    ]
                )
            self.assertEqual(raised.exception.code, 2)
            self.assertFalse(output.exists())

    def test_missing_raw_input_is_rejected_before_runtime_setup(self):
        with tempfile.TemporaryDirectory() as tmpdir, mock.patch.object(
            bulk_footprinting, "prepare_command_runtime"
        ) as runtime:
            root = Path(tmpdir)
            comparisons = root / "comparisons.tsv"
            comparisons.write_text("comparison\tcond1\tcond2\n", encoding="utf-8")
            with self.assertRaises(SystemExit) as raised:
                bulk_footprinting.main(
                    [
                        "--reads-table",
                        str(root / "missing.tsv"),
                        "--comparison-table",
                        str(comparisons),
                        "--genome",
                        "hg38",
                        "--outdir",
                        str(root / "project"),
                    ]
                )
            self.assertEqual(raised.exception.code, 2)
            runtime.assert_not_called()
            self.assertFalse((root / "project").exists())


if __name__ == "__main__":
    unittest.main()
