import pathlib
import tomllib
import unittest

from fp_tools.gui_config import canonical_tool_name


class PublicApiAliasesTest(unittest.TestCase):
    def test_console_scripts_are_registered(self):
        pyproject = pathlib.Path(__file__).resolve().parents[1] / "pyproject.toml"
        data = tomllib.loads(pyproject.read_text(encoding="utf-8"))

        scripts = data["project"]["scripts"]
        poetry_scripts = data["tool"]["poetry"]["scripts"]
        expected = {
            "atac-correct": "fp_tools.cli:main",
            "call-footprints": "fp_tools.cli_scorebigwig:main",
            "match-motifs": "fp_tools.tools.match_motifs:main",
            "diff-footprints": "fp_tools.tools.diff_footprints:diff_footprints_cli",
            "normalize-bigwig": "fp_tools.tools.normalize_bigwig:main",
            "plot-aggregate": "fp_tools.cli_plotaggregate:main",
            "review-multi-comparisons": "fp_tools.tools.review_multi_comparisons:main",
            "plot-motif-aggregate-grid": "fp_tools.cli_compat:plot_motif_aggregate_grid_main",
            "prepare-atac": "fp_tools.tools.prepare_atac:main",
            "bulk-footprinting": "fp_tools.tools.bulk_footprinting:main",
            "run-yaml-workflow": "fp_tools.cli_batch:main",
            "run-workflow": "fp_tools.cli_compat:run_workflow_main",
            "fp-tools-gui": "fp_tools.cli_gui:main",
            "discover-motifs": "fp_tools.tools.motif_discovery:motif_discovery_plan_main",
            "summarize-motifs": "fp_tools.tools.motif_discovery:motif_report_main",
            "motif-discovery": "fp_tools.cli_compat:motif_discovery_main",
            "motif-summary": "fp_tools.cli_compat:motif_summary_main",
            "pseudobulk-fragments": "fp_tools.tools.pseudobulk:main",
            "find-signature-fp": "fp_tools.tools.find_signature_fp:main",
            "sc-footprinting": "fp_tools.tools.pseudobulk_footprints:main",
            "pseudobulk-footprints": "fp_tools.cli_compat:pseudobulk_footprints_main",
        }

        self.assertEqual(scripts, expected)
        self.assertEqual(poetry_scripts, expected)

    def test_config_accepts_public_names(self):
        aliases = {
            "atac-correct": "atac-correct",
            "call-footprints": "call-footprints",
            "match-motifs": "match-motifs",
            "diff-footprints": "diff-footprints",
            "normalize-bigwig": "normalize-bigwig",
            "plot-aggregate": "plot-aggregate",
            "review-multi-comparisons": "review-multi-comparisons",
            "plot-motif-aggregate-grid": "plot-aggregate",
            "prepare-atac": "prepare-atac",
            "bulk-footprinting": "bulk-footprinting",
            "run-workflow": "run-yaml-workflow",
            "discover-motifs": "discover-motifs",
            "motif-discovery": "discover-motifs",
            "summarize-motifs": "summarize-motifs",
            "motif-summary": "summarize-motifs",
            "pseudobulk-fragments": "pseudobulk-fragments",
            "find-signature-fp": "find-signature-fp",
            "sc-footprinting": "sc-footprinting",
            "pseudobulk-footprints": "sc-footprinting",
        }

        for alias, canonical in aliases.items():
            self.assertEqual(canonical_tool_name(alias), canonical)


if __name__ == "__main__":
    unittest.main()
