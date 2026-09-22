"""Regressions for executable documentation and installed-package presentation."""

import argparse
from pathlib import Path
import re
import tomllib

from fp_tools.parsers import add_aggregate_arguments

ROOT = Path(__file__).resolve().parents[1]


def test_aggregate_guide_options_are_accepted_including_output_prose():
    parser = add_aggregate_arguments(argparse.ArgumentParser())
    guide = (ROOT / "docs/get-started/commands/plot-aggregate.md").read_text()
    options = set(re.findall(r"--[A-Za-z][A-Za-z0-9_-]*", guide))
    assert options <= set(parser._option_string_actions)
    assert "--share_y" not in parser.format_help()


def test_package_description_has_public_logo_and_reader_summary():
    readme = (ROOT / "README.md").read_text()
    logo = re.search(r'<img src="([^"]+)"', readme).group(1)
    assert logo == "https://oncologylab.github.io/fp-tools/assets/fp_tools_logo_horizontal.svg"
    meta = tomllib.loads((ROOT / "pyproject.toml").read_text())
    assert meta["project"]["description"] == meta["tool"]["poetry"]["description"]
    assert "Cython" not in meta["project"]["description"]


def test_command_reference_preserves_original_title_anchor():
    reference = (ROOT / "docs/api.md").read_text(encoding="utf-8")
    generator = (ROOT / "scripts/generate_api_reference.py").read_text(encoding="utf-8")
    assert 'id="api-reference"' in reference
    assert 'id="api-reference"' in generator


def test_local_container_examples_are_pinned_and_loopback_only():
    for file in ("README.md", "docs/get-started/installation.md"):
        text = (ROOT / file).read_text()
        assert "fp-tools.git#main" not in text
        assert "Complete reproducible environment" not in text
    text = (ROOT / "docs/get-started/installation.md").read_text()
    assert "-p 127.0.0.1:8891:8891" in text


def test_single_cell_help_describes_count_matrix(capsys):
    import pytest
    from fp_tools.tools.find_signature_fp import main as signature_main
    from fp_tools.tools.pseudobulk_footprints import build_parser as sc_parser
    with pytest.raises(SystemExit) as result:
        signature_main(["--help"])
    assert result.value.code == 0
    signature_help = " ".join(capsys.readouterr().out.split())
    for help_text in (signature_help, sc_parser()._option_string_actions["--h5ad"].help):
        assert "genomic-bin counts" in help_text
        assert "selected" in help_text
