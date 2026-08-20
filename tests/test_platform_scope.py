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

    def test_bulk_wrapper_rejects_raw_input_and_preparation_options(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "project"
            for extra in (("--reads-table", "reads.tsv"), ("--profile", "modern")):
                with self.subTest(flag=extra[0]), self.assertRaises(SystemExit) as raised:
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
                            *extra,
                        ]
                    )
                self.assertEqual(raised.exception.code, 2)
                self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
