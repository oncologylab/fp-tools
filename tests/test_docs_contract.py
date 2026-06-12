"""Contract tests that keep documentation in sync with the packaged entry points.

These guard against the README/MANUAL/PyPI drift where the documented command
surface diverges from the console scripts declared in ``pyproject.toml``.
"""

import pathlib
import re
import tomllib
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[1]

# Backward-compatibility aliases for the classical TOBIAS-style command names.
# They are intentionally not part of the primary ``--help`` verification block.
LEGACY_ALIASES = {
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
        self.manual = (ROOT / "MANUAL.md").read_text(encoding="utf-8")

    def test_setuptools_and_poetry_scripts_match(self):
        poetry_scripts = self.data["tool"]["poetry"]["scripts"]
        self.assertEqual(
            self.project_scripts,
            poetry_scripts,
            "[project.scripts] and [tool.poetry.scripts] have drifted apart.",
        )

    def test_primary_entry_points_are_documented_in_readme_and_manual(self):
        primary = set(self.project_scripts) - LEGACY_ALIASES
        for command in primary:
            self.assertIn(command, self.readme, f"{command} is missing from README.md")
            self.assertIn(command, self.manual, f"{command} is missing from MANUAL.md")

    def test_help_block_exactly_covers_non_alias_commands(self):
        documented = _verify_help_commands(self.readme)
        expected = set(self.project_scripts) - LEGACY_ALIASES
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

    def test_tobias_compatible_aliases_are_registered_but_not_primary_help_commands(self):
        for alias in LEGACY_ALIASES:
            self.assertIn(alias, self.project_scripts)
        documented = _verify_help_commands(self.readme)
        self.assertFalse(documented & LEGACY_ALIASES)

    def test_gui_extra_is_declared_and_documented(self):
        extras = self.data["project"].get("optional-dependencies", {})
        self.assertIn("gui", extras, "Expected a [project.optional-dependencies] gui extra.")
        self.assertTrue(
            any("streamlit" in dep for dep in extras["gui"]),
            "The gui extra should provide streamlit.",
        )
        self.assertNotIn(
            "streamlit",
            "\n".join(self.data["project"]["dependencies"]),
            "streamlit should be an optional extra, not a core dependency.",
        )
        self.assertIn('fp-tools-bio[gui]', self.readme)


if __name__ == "__main__":
    unittest.main()
