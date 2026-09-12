"""Keep the static GUI demo's copyable YAML compatible with real commands."""

import argparse
import contextlib
import html
import io
import json
from pathlib import Path
import re
from unittest.mock import patch

import pytest

from fp_tools.cli_batch import run_config_file
from fp_tools.gui_config import expand_jobs, parse_yaml_text, validate_config
from fp_tools.parsers import (
    add_aggregate_arguments,
    add_atacorrect_arguments,
    add_diff_footprints_arguments,
    add_scorebigwig_arguments,
)


DEMO = Path(__file__).resolve().parents[1] / "docs/demos/gui/fp-tools-gui-static-demo.html"


def _examples():
    source = DEMO.read_text(encoding="utf-8")
    literals = re.findall(r'^\s*yaml: ("(?:[^"\\]|\\.)*"),?$', source, re.MULTILINE)
    examples = [json.loads(value) for value in literals]
    examples.extend(html.unescape(value) for value in re.findall(r"<pre>(tool: .*?)</pre>", source, re.DOTALL))
    assert len(examples) == 11, "Check every command example and the initial home preview"
    return examples


def _main_parser(main):
    """Capture parsers built inside main without reading inputs or running jobs."""
    captured = []

    class ParserCaptured(Exception):
        pass

    def capture(parser, *args, **kwargs):
        captured.append(parser)
        raise ParserCaptured

    with patch.object(argparse.ArgumentParser, "parse_args", capture):
        with pytest.raises(ParserCaptured):
            main([])
    return captured[0]


def _parser(tool):
    if tool == "atac-correct":
        return add_atacorrect_arguments(argparse.ArgumentParser())
    if tool == "call-footprints":
        return add_scorebigwig_arguments(argparse.ArgumentParser())
    if tool in {"match-motifs", "diff-footprints"}:
        return add_diff_footprints_arguments(argparse.ArgumentParser(), command_name=tool)
    if tool == "plot-aggregate":
        return add_aggregate_arguments(argparse.ArgumentParser())
    if tool == "normalize-bigwig":
        from fp_tools.tools.normalize_bigwig import build_parser

        return build_parser()
    if tool == "sc-footprinting":
        from fp_tools.tools.pseudobulk_footprints import build_parser

        return build_parser()
    if tool == "pseudobulk-fragments":
        from fp_tools.tools.pseudobulk import main

        return _main_parser(main)
    if tool == "find-signature-fp":
        from fp_tools.tools.find_signature_fp import main

        return _main_parser(main)
    raise AssertionError(f"No parser check for demo tool: {tool}")


@pytest.mark.parametrize("yaml_text", _examples())
def test_demo_yaml_validates_and_expands_to_supported_arguments(yaml_text):
    config = parse_yaml_text(yaml_text)
    assert validate_config(config) == []
    jobs = expand_jobs(config)
    assert len(jobs) == 1
    job = jobs[0]
    args = _parser(job.tool).parse_args(job.command[1:])
    if job.tool == "sc-footprinting" and args.fragments:
        assert args.genome_sizes, "Fragment-based workflows require chromosome sizes"


@pytest.mark.parametrize("yaml_text", _examples())
def test_demo_yaml_dry_run_needs_no_example_input_files(yaml_text, tmp_path):
    path = tmp_path / "example.yml"
    path.write_text(yaml_text, encoding="utf-8")
    output = io.StringIO()
    with contextlib.redirect_stdout(output):
        assert run_config_file(path, run_root=tmp_path / "runs", dry_run=True) == 0
    assert expand_jobs(parse_yaml_text(yaml_text))[0].tool in output.getvalue()
    assert not (tmp_path / "runs").exists()
