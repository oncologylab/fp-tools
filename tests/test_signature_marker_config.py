"""Marker lists from GUI/YAML must retain the single-value CLI contract."""

import argparse
import subprocess
import sys
from copy import deepcopy
from unittest.mock import patch

import pytest

from fp_tools.gui_config import (
    build_cli_command, config_to_yaml_text, expand_jobs, parse_yaml_text, validate_config,
)


@pytest.fixture
def signature_parser():
    from fp_tools.tools.find_signature_fp import main

    captured = []

    class Captured(Exception):
        pass

    def capture(parser, *args, **kwargs):
        captured.append(parser)
        raise Captured

    with patch.object(argparse.ArgumentParser, "parse_args", capture):
        with pytest.raises(Captured):
            main([])
    return captured[0]


@pytest.mark.parametrize("markers", [
    ["STAT6", "FOSB", "CEBPA", "IRF8", "RELA", "ZNF683", "NR4A1", "SMAD3"],
    ["STAT6", "CEBPA"], ["STAT6"], "STAT6,CEBPA", [], None,
])
@pytest.mark.parametrize("section", ["samples", "comparisons"])
def test_marker_yaml_roundtrip_and_parser(signature_parser, markers, section):
    from fp_tools.tools.find_signature_fp import MARKERS

    config = {
        "defaults": {"tool": "find-signature-fp", "annotations": "annotations.tsv",
                     "fragments": "fragments.tsv", "h5ad": "cells.h5ad",
                     "outdir": "output", "tf_site_dir": "motif_sites", "markers": markers},
        section: [{"job_id": "signature"}],
    }
    original = deepcopy(config)
    assert validate_config(config) == []
    job = expand_jobs(parse_yaml_text(config_to_yaml_text(config)))[0]
    parsed = signature_parser.parse_args(job.command[1:])
    expected = ",".join(markers) if isinstance(markers, list) else markers
    assert parsed.markers == (expected or ",".join(MARKERS))
    assert config == original


@pytest.mark.parametrize("override", [["ZNF683", "STAT6"], [], "CEBPA,STAT6"])
def test_explicit_marker_override(override):
    config = {"defaults": {"tool": "find-signature-fp", "markers": ["FOSB", "IRF8"]},
              "samples": [{"markers": override}]}
    command = expand_jobs(config)[0].command
    expected = ",".join(override) if isinstance(override, list) else override
    assert command == (["find-signature-fp", "--markers", expected] if expected
                       else ["find-signature-fp"])


def test_gui_default_editor_produces_runnable_marker_argument():
    from fp_tools.gui_app import GENERIC_TOOL_DEFAULTS, _prepare_generic_params

    params = _prepare_generic_params("find-signature-fp", {
        "markers": GENERIC_TOOL_DEFAULTS["find-signature-fp"]["markers"],
    })
    assert len(params["markers"]) == 8
    command = build_cli_command("find-signature-fp", params)
    assert command == ["find-signature-fp", "--markers", ",".join(params["markers"])]


def test_other_list_arguments_keep_separate_values():
    assert build_cli_command("match-motifs", {"motifs": ["one.meme", "two.meme"]}) == [
        "match-motifs", "--motifs", "one.meme", "two.meme",
    ]
    assert build_cli_command("another-tool", {"markers": ["A", "B"]}) == [
        "another-tool", "--markers", "A", "B",
    ]


def test_list_yaml_matches_direct_signature_analysis(tmp_path):
    from fp_tools.cli_batch import run_config_file
    from scripts.smoke_desktop_bundle import (
        assert_signature_outputs, signature_config, write_signature_fixture,
    )

    inputs = write_signature_fixture(tmp_path / "inputs")
    direct_output = tmp_path / "direct"
    direct = signature_config(inputs, direct_output)
    direct["markers"] = ",".join(direct["markers"])
    command = expand_jobs(direct)[0].command
    result = subprocess.run(
        [sys.executable, "-m", "fp_tools.tools.find_signature_fp", *command[1:]],
        capture_output=True, text=True, timeout=120,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    yaml_output = tmp_path / "yaml"
    path = tmp_path / "signature.yml"
    path.write_text(config_to_yaml_text(signature_config(inputs, yaml_output)), encoding="utf-8")
    assert run_config_file(path, run_root=tmp_path / "runs") == 0
    assert_signature_outputs(yaml_output, direct_output)
