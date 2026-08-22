from __future__ import annotations

import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

from fp_tools import gui_app
from fp_tools.command_registry import COMMAND_TARGETS, dispatch_command
from fp_tools.desktop import (
    INTERNAL_GUI_SERVER_FLAG,
    INTERNAL_LIST_EXAMPLES_FLAG,
    INTERNAL_NATIVE_SMOKE_FLAG,
    main as desktop_main,
)
from fp_tools.desktop_window import DesktopLaunchError, _parse_desktop_args, _server_command
from fp_tools.gui_app import (
    GENERIC_TOOL_DEFAULTS,
    _all_example_files,
    _config_form_mode,
    _config_widget_key,
    _format_extra_args,
    _parse_extra_args,
    _set_config,
    _updated_single_config,
    _validation_errors_markup,
)
from fp_tools.cli_gui import _access_urls, _startup_messages
from fp_tools.utils.subprocess_commands import (
    fp_tools_subprocess_command,
    python_script_subprocess_command,
    resolve_fp_tools_subprocess,
)


class GuiLauncherTests(unittest.TestCase):
    class _SessionState(dict):
        def __getattr__(self, name):
            try:
                return self[name]
            except KeyError as exc:
                raise AttributeError(name) from exc

        def __setattr__(self, name, value):
            self[name] = value

    def test_desktop_default_opens_native_window(self):
        with patch("fp_tools.desktop_window.launch_native_gui", return_value=0) as launch:
            self.assertEqual(desktop_main(["--run-dir", "runs"]), 0)
        launch.assert_called_once_with(["--run-dir", "runs"], auto_close=False)

    def test_desktop_native_smoke_uses_auto_close(self):
        with patch("fp_tools.desktop_window.launch_native_gui", return_value=0) as launch:
            self.assertEqual(
                desktop_main([INTERNAL_NATIVE_SMOKE_FLAG, "--run-dir", "runs"]),
                0,
            )
        launch.assert_called_once_with(["--run-dir", "runs"], auto_close=True)

    def test_desktop_internal_server_never_opens_browser(self):
        with patch("fp_tools.cli_gui.main") as launch:
            self.assertEqual(
                desktop_main([INTERNAL_GUI_SERVER_FLAG, "--port", "8898"]),
                0,
            )
        launch.assert_called_once_with(["--port", "8898", "--no-browser"])

    def test_desktop_window_is_local_and_spawns_private_server(self):
        args = _parse_desktop_args(["--port", "8899", "--run-dir", "runs"])
        self.assertEqual(args.port, 8899)
        command = _server_command(args.port, Path(args.run_dir))
        self.assertEqual(command[1], INTERNAL_GUI_SERVER_FLAG)
        self.assertIn("--no-browser", command)
        with self.assertRaises(DesktopLaunchError):
            _parse_desktop_args(["--host", "0.0.0.0"])

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

        fake_streamlit = MagicMock()
        with patch.object(gui_app, "st", fake_streamlit):
            gui_app._apply_page_style()
        style = fake_streamlit.markdown.call_args.args[0]
        self.assertIn("overflow-wrap: anywhere !important", style)
        self.assertIn("word-break: normal !important", style)
        self.assertIn("button:disabled", style)
        self.assertIn("padding-top: 4.25rem !important", style)
        self.assertIn('[data-testid="stTextInputRootElement"]', style)
        self.assertIn('[data-testid="stTextAreaRootElement"]', style)
        self.assertIn('[data-testid="stNumberInputContainer"]', style)
        self.assertIn(".stFormSubmitButton > button", style)
        self.assertIn("@media (max-width: 1500px)", style)

    def test_page_heading_reuses_the_home_hero_card(self):
        fake_streamlit = MagicMock()
        with patch.object(gui_app, "st", fake_streamlit):
            gui_app._render_page_heading("call-footprints", "Score corrected signal.")
        markup = fake_streamlit.markdown.call_args.args[0]
        self.assertIn('class="fp-hero fp-page-heading"', markup)
        self.assertIn("call-footprints", markup)
        self.assertIn("Score corrected signal.", markup)

    def test_fresh_gui_defaults_do_not_assume_repository_paths(self):
        self.assertEqual(GENERIC_TOOL_DEFAULTS["match-motifs"]["signals"], "")
        self.assertEqual(GENERIC_TOOL_DEFAULTS["bulk-footprinting"]["sample_table"], "")
        self.assertEqual(GENERIC_TOOL_DEFAULTS["sc-footprinting"]["fragments"], "")
        self.assertFalse(GENERIC_TOOL_DEFAULTS["sc-footprinting"]["dry_run"])

    def test_validation_messages_are_humanized_without_losing_context(self):
        self.assertEqual(
            gui_app._friendly_validation_message("bulk_run: 'sample_table' file does not exist: /work/samples.tsv"),
            "Samples TSV File not found: /work/samples.tsv",
        )
        self.assertEqual(
            gui_app._friendly_validation_message("Missing fragments: C:\\work\\cells.tsv.gz"),
            "Missing fragments: C:\\work\\cells.tsv.gz",
        )

    def test_mobile_navigation_and_compact_sidebar_styles_are_present(self):
        fake_streamlit = MagicMock()
        with patch.object(gui_app, "st", fake_streamlit):
            gui_app._apply_page_style()
        style = fake_streamlit.markdown.call_args.args[0]
        self.assertIn('[data-testid="stExpandSidebarButton"]', style)
        self.assertIn('content: "Open navigation"', style)
        self.assertIn('content: "Close navigation"', style)
        self.assertIn(".fp-workspace-path", style)
        self.assertIn('[data-testid="stMain"] [data-testid="stWidgetLabel"]', style)
        self.assertIn("-webkit-text-fill-color: var(--fp-text) !important", style)
        self.assertIn("flex-direction: column !important", style)
        self.assertIn(
            '[data-testid="stHorizontalBlock"] > [data-testid="stColumn"]',
            style,
        )
        self.assertNotIn(".fp-sidebar-brand-subtitle", style)
        self.assertNotIn(".fp-run-dir-pill", style)

    def test_set_config_advances_widget_revision_and_reruns(self):
        state = self._SessionState(config_revision=3)
        fake_streamlit = MagicMock(session_state=state)
        config = {
            "version": 1,
            "run_mode": "single",
            "defaults": {},
            "samples": [{"sample_id": "loaded", "tool": "bulk-footprinting"}],
            "comparisons": [],
        }
        with patch.object(gui_app, "st", fake_streamlit):
            _set_config(config)
            self.assertEqual(_config_widget_key("sample_table"), "cfg_4_sample_table")
        self.assertEqual(state["current_config"]["samples"], config["samples"])
        self.assertIsNone(state["current_config"]["run_root"])
        self.assertEqual(state["config_revision"], 4)
        self.assertEqual(state["config_update_notice"], "Current config updated.")
        fake_streamlit.rerun.assert_called_once_with()

    def test_config_update_notice_uses_streamlit_default_icon(self):
        state = self._SessionState(config_update_notice="Current config updated.")
        fake_streamlit = MagicMock(session_state=state)
        with patch.object(gui_app, "st", fake_streamlit):
            gui_app._render_config_update_notice()
        fake_streamlit.toast.assert_called_once_with("Current config updated.")
        self.assertNotIn("config_update_notice", state)

    def test_query_page_overrides_stale_session_navigation(self):
        state = self._SessionState(gui_page="bulk-footprinting")
        fake_streamlit = MagicMock(
            session_state=state,
            query_params={"page": "normalize-bigwig"},
        )
        with patch.object(gui_app, "st", fake_streamlit):
            page = gui_app._current_page_from_query()
        self.assertEqual(page, "normalize-bigwig")
        self.assertEqual(state["gui_page"], "normalize-bigwig")

    def test_loaded_single_config_preserves_metadata_defaults_and_hidden_fields(self):
        state = self._SessionState(
            config_revision=2,
            current_config={
                "version": 1,
                "run_mode": "single",
                "run_root": "runs",
                "defaults": {"cores": 2, "normalization": "sample-quantile"},
                "samples": [
                    {
                        "sample_id": "loaded_bulk",
                        "tool": "bulk-footprinting",
                        "sample_table": r"C:\\project path\\samples.tsv",
                        "hidden_future_option": "keep-me",
                    }
                ],
                "comparisons": [],
            },
        )
        fake_streamlit = MagicMock(session_state=state)
        with patch.object(gui_app, "st", fake_streamlit):
            self.assertEqual(_config_form_mode("bulk-footprinting"), "Single run")
            self.assertEqual(gui_app._current_single_params("bulk-footprinting")["cores"], 2)
            updated = _updated_single_config(
                "bulk-footprinting",
                {"outdir": r"C:\\project path\\changed", "cores": 6},
                job_id="fallback",
            )
        item = updated["samples"][0]
        self.assertEqual(item["sample_id"], "loaded_bulk")
        self.assertEqual(item["hidden_future_option"], "keep-me")
        self.assertEqual(item["outdir"], r"C:\\project path\\changed")
        self.assertEqual(updated["run_root"], "runs")
        self.assertEqual(updated["defaults"]["normalization"], "sample-quantile")

    def test_loaded_batch_config_selects_matching_editor_mode(self):
        fake_streamlit = MagicMock(
            session_state=self._SessionState(
                current_config={
                    "version": 1,
                    "run_mode": "batch",
                    "defaults": {},
                    "samples": [],
                    "comparisons": [
                        {"comparison_id": "a_vs_b", "tool": "diff-footprints"}
                    ],
                }
            )
        )
        with patch.object(gui_app, "st", fake_streamlit):
            self.assertEqual(
                _config_form_mode("diff-footprints"),
                "Batch comparison list",
            )

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
