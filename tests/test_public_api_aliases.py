import pathlib
import tomllib
import unittest

from fp_tools.gui_config import canonical_tool_name


class PublicApiAliasesTest(unittest.TestCase):
    def test_new_console_scripts_are_registered(self):
        pyproject = pathlib.Path(__file__).resolve().parents[1] / "pyproject.toml"
        data = tomllib.loads(pyproject.read_text(encoding="utf-8"))

        scripts = data["project"]["scripts"]
        poetry_scripts = data["tool"]["poetry"]["scripts"]
        expected = {
            "atac-correct": "fp_tools.cli:main",
            "score-footprints": "fp_tools.cli_scorebigwig:main",
            "detect-tf-binding": "fp_tools.tools.bindetect:run_cli",
            "plot-aggregate": "fp_tools.cli_plotaggregate:main",
            "fp-tools-build-tfbs-features": "fp_tools.tools.tfbs_features:main",
            "fp-tools-train-tfbs-model": "fp_tools.tools.tfbs_model:train_main",
            "fp-tools-predict-tfbs": "fp_tools.tools.tfbs_model:predict_main",
            "fp-tools-generate-candidates": "fp_tools.tools.candidates:main",
            "fp-tools-rerank-candidates": "fp_tools.tools.rerank:main",
            "fp-tools-export-candidate-fasta": "fp_tools.tools.motif_discovery:export_fasta_main",
            "fp-tools-meme-command": "fp_tools.tools.motif_discovery:meme_command_main",
            "fp-tools-motif-discovery": "fp_tools.tools.motif_discovery:motif_discovery_plan_main",
            "fp-tools-summarize-motifs": "fp_tools.tools.motif_discovery:motif_report_main",
            "fp-tools-score-variants": "fp_tools.tools.variants:main",
            "fp-tools-pseudobulk": "fp_tools.tools.pseudobulk:main",
            "fp-tools-bindetect-replicate-report": "fp_tools.tools.bindetect_replicate_report:main",
            "fp-tools-decompose-competition": "fp_tools.tools.footprint_competition:main",
        }

        for name, target in expected.items():
            self.assertEqual(scripts[name], target)
            self.assertEqual(poetry_scripts[name], target)

        self.assertNotIn("fp-tools-motif-discovery-plan", scripts)
        self.assertNotIn("fp-tools-motif-discovery-plan", poetry_scripts)

    def test_config_accepts_new_aliases_and_legacy_names(self):
        aliases = {
            "atac-correct": "ATACorrect",
            "ATACorrect": "ATACorrect",
            "score-footprints": "FootprintScores",
            "FootprintScores": "FootprintScores",
            "detect-tf-binding": "BINDetect",
            "BINDetect": "BINDetect",
            "plot-aggregate": "PlotAggregate",
            "PlotAggregate": "PlotAggregate",
        }

        for alias, canonical in aliases.items():
            self.assertEqual(canonical_tool_name(alias), canonical)


if __name__ == "__main__":
    unittest.main()
