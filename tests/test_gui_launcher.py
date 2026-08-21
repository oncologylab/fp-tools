from __future__ import annotations

import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

from fp_tools import gui_app
from fp_tools.command_registry import COMMAND_TARGETS, dispatch_command
from fp_tools.desktop import INTERNAL_LIST_EXAMPLES_FLAG, main as desktop_main
from fp_tools.gui_app import (
    GENERIC_TOOL_DEFAULTS,
    _all_example_files,
    _format_extra_args,
    _parse_extra_args,
    _validation_errors_markup,
)
from fp_tools.cli_gui import _access_urls, _startup_messages
from fp_tools.utils.subprocess_commands import (
    fp_tools_subprocess_command,
    python_script_subprocess_command,
    resolve_fp_tools_subprocess,
)


class GuiLauncherTests(unittest.TestCase):
    def test_validation_error_markup_wraps_paths_and_escapes_html(self):
        markup = _validation_errors_markup(
            [
                r"Missing fragments: C:\Users\Researcher\a very long project\inputs\cells.fragments.tsv.gz",
                "Invalid <sample> & comparison",
            ]
        )

        self.assertIn('class="fp-validation-errors"', markup)
        self.assertEqual(markup.count("<li>"), 2)
        self.assertIn("&lt;sample&gt; &amp; comparison", markup)
        self.assertNotIn("<sample>", markup)

    def test_run_control_state_tracks_validation_and_keeps_launch_guard(self):
        normalized = {
            "version": 1,
            "run_mode": "single",
            "defaults": {},
            "samples": [
                {
                    "tool": "bulk-footprinting",
                    "sample_table": "samples.tsv",
                    "comparison_table": "comparisons.tsv",
                    "genome": "genome.fa",
                    "outdir": "results",
                }
            ],
            "comparisons": [],
        }
        fake_streamlit = MagicMock()
        fake_streamlit.session_state.current_config = normalized
        fake_streamlit.button.return_value = True

        with (
            patch.object(gui_app, "st", fake_streamlit),
            patch.object(gui_app, "normalize_config", return_value=normalized),
            patch.object(gui_app, "validate_gui_config", return_value=["missing input"]),
            patch.object(gui_app, "_current_config_tool", return_value="bulk-footprinting"),
            patch.object(gui_app, "materialize_run_config") as materialize,
            patch.object(gui_app, "launch_config_async") as launch,
        ):
            gui_app._render_run_controls(Path("runs"), label="bulk_footprinting")

        self.assertTrue(fake_streamlit.button.call_args.kwargs["disabled"])
        materialize.assert_not_called()
        launch.assert_not_called()

        fake_streamlit.reset_mock()
        fake_streamlit.session_state.current_config = normalized
        fake_streamlit.button.return_value = False
        with (
            patch.object(gui_app, "st", fake_streamlit),
            patch.object(gui_app, "normalize_config", return_value=normalized),
            patch.object(gui_app, "validate_gui_config", return_value=[]),
            patch.object(gui_app, "_current_config_tool", return_value="bulk-footprinting"),
        ):
            gui_app._render_run_controls(Path("runs"), label="bulk_footprinting")

        self.assertFalse(fake_streamlit.button.call_args.kwargs["disabled"])

    def test_local_startup_explains_browser_and_ssh_access(self):
        text = "\n".join(_startup_messages("127.0.0.1", 8891, Path("runs")))
        self.assertIn("http://127.0.0.1:8891", text)
        self.assertIn("ssh -N -L 8891:127.0.0.1:8891 USER@SERVER", text)
        self.assertIn("Ctrl+C", text)

    def test_network_mode_warns_about_authentication(self):
        text = "\n".join(_startup_messages("0.0.0.0", 8891, Path("runs")))
        self.assertIn("http://SERVER_IP:8891", text)
        self.assertIn("does not add authentication", text)
        self.assertEqual(
            _access_urls("0.0.0.0", 8891),
            ["http://127.0.0.1:8891", "http://SERVER_IP:8891"],
        )

    def test_desktop_registry_covers_bam_first_commands(self):
        expected = {
            "bulk-footprinting",
            "atac-correct",
            "call-footprints",
            "match-motifs",
            "diff-footprints",
            "normalize-bigwig",
            "plot-aggregate",
            "review-multi-comparisons",
            "run-yaml-workflow",
            "fp-tools-gui",
            "fp-tools-runtime",
            "discover-motifs",
            "summarize-motifs",
            "pseudobulk-fragments",
            "find-signature-fp",
            "sc-footprinting",
        }
        self.assertTrue(expected.issubset(COMMAND_TARGETS))
        self.assertNotIn("prepare-atac", COMMAND_TARGETS)

    def test_frozen_child_command_routes_through_desktop_executable(self):
        with patch("fp_tools.utils.subprocess_commands.is_frozen", return_value=True):
            command = fp_tools_subprocess_command("plot-aggregate", ["--help"])
        self.assertEqual(command[0], __import__("sys").executable)
        self.assertEqual(command[1:3], ["--fp-tools-internal-command", "plot-aggregate"])
        self.assertEqual(command[-1], "--help")

    def test_non_fp_tools_child_command_is_not_rewritten(self):
        self.assertEqual(resolve_fp_tools_subprocess(["samtools", "--version"]), ["samtools", "--version"])

    def test_frozen_python_helper_routes_through_desktop_executable(self):
        with patch("fp_tools.utils.subprocess_commands.is_frozen", return_value=True):
            command = python_script_subprocess_command("helper.py", ["--value", "1"])
        self.assertEqual(command[0], __import__("sys").executable)
        self.assertEqual(command[1:3], ["--fp-tools-internal-python-script", "helper.py"])
        self.assertEqual(command[-2:], ["--value", "1"])

    def test_dispatch_restores_sys_argv(self):
        import sys

        previous = sys.argv
        with self.assertRaises(SystemExit) as raised:
            dispatch_command("summarize-motifs", ["--help"])
        self.assertEqual(raised.exception.code, 0)
        self.assertIs(sys.argv, previous)

    def test_desktop_lists_packaged_bam_first_examples(self):
        with patch("builtins.print") as printer:
            self.assertEqual(desktop_main([INTERNAL_LIST_EXAMPLES_FLAG]), 0)
        names = {str(call.args[0]) for call in printer.call_args_list}
        self.assertIn("bulk_footprinting_bam.yml", names)

    def test_packaged_examples_are_available_outside_the_source_directory(self):
        packaged = {path.name: path.read_text(encoding="utf-8") for path in _all_example_files()}
        names = set(packaged)
        self.assertIn("bulk_footprinting_bam.yml", names)
        self.assertIn("call_footprints_single.yml", names)
        source_dir = Path(__file__).resolve().parents[1] / "examples" / "gui_configs"
        for source in source_dir.glob("*.yml"):
            self.assertEqual(packaged[source.name], source.read_text(encoding="utf-8"))

    def test_generic_forms_keep_supported_blank_fields(self):
        self.assertIn("sample_table", GENERIC_TOOL_DEFAULTS["bulk-footprinting"])
        self.assertNotIn("reads_table", GENERIC_TOOL_DEFAULTS["bulk-footprinting"])
        self.assertTrue(
            {"meme_txt", "tomtom_tsv"}.issubset(GENERIC_TOOL_DEFAULTS["summarize-motifs"])
        )
        self.assertTrue(
            {"labels", "output_dir", "output_html"}.issubset(
                GENERIC_TOOL_DEFAULTS["review-multi-comparisons"]
            )
        )
        self.assertTrue(
            {"candidates", "fasta", "script", "execute"}.issubset(
                GENERIC_TOOL_DEFAULTS["discover-motifs"]
            )
        )

    def test_extra_args_round_trip_windows_path_with_spaces(self):
        arguments = ["--meme-txt", r"C:\Research data\meme.txt", "--title", "B cell motifs"]
        self.assertEqual(_parse_extra_args(_format_extra_args(arguments)), arguments)


if __name__ == "__main__":
    unittest.main()
