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
            "match-motifs": "fp_tools.tools.bindetect:match_motifs_cli",
            "diff-footprints": "fp_tools.tools.bindetect:diff_footprints_cli",
            "plot-aggregate": "fp_tools.cli_plotaggregate:main",
            "plot-aggregate-batch": "fp_tools.tools.plot_aggregate_batch:main",
            "run-workflow": "fp_tools.cli_batch:main",
            "motif-discovery": "fp_tools.tools.motif_discovery:motif_discovery_plan_main",
            "motif-summary": "fp_tools.tools.motif_discovery:motif_report_main",
            "pseudobulk-fragments": "fp_tools.tools.pseudobulk:main",
            "ATACorrect": "fp_tools.cli:main",
            "FootprintScores": "fp_tools.cli_scorebigwig:main",
            "ScoreBigwig": "fp_tools.cli_scorebigwig:main",
            "BINDetect": "fp_tools.tools.bindetect:run_cli",
            "PlotAggregate": "fp_tools.cli_plotaggregate:main",
        }

        self.assertEqual(scripts, expected)
        self.assertEqual(poetry_scripts, expected)

    def test_config_accepts_public_names_and_tobias_aliases(self):
        aliases = {
            "atac-correct": "atac-correct",
            "ATACorrect": "atac-correct",
            "call-footprints": "call-footprints",
            "FootprintScores": "call-footprints",
            "ScoreBigwig": "call-footprints",
            "match-motifs": "match-motifs",
            "diff-footprints": "diff-footprints",
            "BINDetect": "diff-footprints",
            "plot-aggregate": "plot-aggregate",
            "PlotAggregate": "plot-aggregate",
        }

        for alias, canonical in aliases.items():
            self.assertEqual(canonical_tool_name(alias), canonical)


if __name__ == "__main__":
    unittest.main()
