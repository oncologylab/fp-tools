"""Contract tests that keep documentation in sync with packaged entry points.

These guard against README/MkDocs/PyPI drift where the documented command
surface diverges from the console scripts declared in ``pyproject.toml``.
"""

import pathlib
import re
import tomllib
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[1]

REMOVED_ALIASES = {
    "ATACorrect",
    "FootprintScores",
    "ScoreBigwig",
    "BINDetect",
    "PlotAggregate",
}


def _load_pyproject():
    return tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))


def _verify_help_commands(readme: str):
    """Return the set of commands that have a ``<cmd> --help`` line in the README."""
    return set(re.findall(r"^([\w.-]+) --help$", readme, flags=re.MULTILINE))


class DocsEntryPointContractTest(unittest.TestCase):
    def setUp(self):
        self.data = _load_pyproject()
        self.project_scripts = self.data["project"]["scripts"]
        self.readme = (ROOT / "README.md").read_text(encoding="utf-8")
        self.site_docs = "\n".join(
            (ROOT / path).read_text(encoding="utf-8")
            for path in [
                "docs/index.md",
                "docs/api.md",
            ]
        )
        self.api_reference = (ROOT / "docs" / "api.md").read_text(encoding="utf-8")

    def test_setuptools_and_poetry_scripts_match(self):
        poetry_scripts = self.data["tool"]["poetry"]["scripts"]
        self.assertEqual(
            self.project_scripts,
            poetry_scripts,
            "[project.scripts] and [tool.poetry.scripts] have drifted apart.",
        )

    def test_primary_entry_points_are_documented_in_readme_and_manual(self):
        for command in self.project_scripts:
            self.assertIn(command, self.readme, f"{command} is missing from README.md")
            self.assertIn(
                command, self.site_docs, f"{command} is missing from MkDocs pages"
            )

    def test_help_block_exactly_covers_non_alias_commands(self):
        documented = _verify_help_commands(self.readme)
        expected = set(self.project_scripts)
        # Every primary command must appear in the README --help verification block ...
        self.assertEqual(
            expected - documented,
            set(),
            "Commands declared in pyproject but missing a `--help` check in README.md.",
        )
        # ... and the verification block must not invent commands that do not exist.
        self.assertEqual(
            documented - set(self.project_scripts),
            set(),
            "README `--help` block references commands that are not entry points.",
        )

    def test_tobias_compatible_aliases_are_removed(self):
        for alias in REMOVED_ALIASES:
            self.assertNotIn(alias, self.project_scripts)
        documented = _verify_help_commands(self.readme)
        self.assertFalse(documented & REMOVED_ALIASES)

    def test_gui_extra_is_declared_and_documented(self):
        extras = self.data["project"].get("optional-dependencies", {})
        self.assertIn(
            "gui", extras, "Expected a [project.optional-dependencies] gui extra."
        )
        self.assertTrue(
            any("streamlit" in dep for dep in extras["gui"]),
            "The gui extra should provide streamlit.",
        )
        self.assertNotIn(
            "streamlit",
            "\n".join(self.data["project"]["dependencies"]),
            "streamlit should be an optional extra, not a core dependency.",
        )
        self.assertIn("fp-tools-bio[gui]", self.readme)

    def test_api_reference_is_command_manual(self):
        for command in self.project_scripts:
            self.assertIn(f"### `{command}`", self.api_reference)
            self.assertIn(f"usage: {command}", self.api_reference)
        self.assertNotIn("::: fp_tools", self.api_reference)

    def test_public_site_docs_use_current_wording(self):
        public_docs = "\n".join(
            (ROOT / path).read_text(encoding="utf-8")
            for path in [
                "README.md",
                "docs/index.md",
                "docs/reports.md",
                "docs/api.md",
            ]
        )
        for stale in [
            "BINDetect",
            "BINDetect-style",
            "Compatibility Aliases",
            "Compatibility Aliases",
            "Open static report",
            "GUI Preview",
            "The PyPI package is",
            "The Python import name is",
            "ATACCorrect",
        ]:
            self.assertNotIn(stale, public_docs)

    def test_prepare_atac_profiles_are_explained_in_plain_language(self):
        public_atac_docs = "\n".join(
            (ROOT / path).read_text(encoding="utf-8")
            for path in ["README.md", "docs/index.md", "docs/api.md"]
        )
        for required in [
            "fastp",
            "mapping quality 30",
            "MACS3",
            "Trim Galore",
            "Picard",
            "HOMER",
            "reads Bowtie2 reports at more than one genomic location",
        ]:
            self.assertIn(required, public_atac_docs)
        for rejected in [
            "CUT&Tag branches",
            "`XS:i:` filtering",
            "legacy-compatible",
            "Bowtie2/Picard",
        ]:
            self.assertNotIn(rejected, public_atac_docs)

    def test_local_markdown_links_exist(self):
        markdown_files = [
            ROOT / "README.md",
            ROOT / "AGENTS.md",
            ROOT / "DEV_PLAN.md",
            ROOT / "RELEASE_CHECKLIST.md",
            *sorted((ROOT / "docs").glob("*.md")),
            *sorted((ROOT / "benchmarks").glob("*.md")),
        ]
        pattern = re.compile(r"\[[^\]]+\]\((?!https?://|mailto:|#)([^)#]+)")
        missing = []
        for source in markdown_files:
            for target in pattern.findall(source.read_text(encoding="utf-8")):
                path = (source.parent / target).resolve()
                if not path.exists():
                    missing.append(f"{source.relative_to(ROOT)} -> {target}")
        self.assertEqual(missing, [])


if __name__ == "__main__":
    unittest.main()
