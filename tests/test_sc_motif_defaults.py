from unittest.mock import patch
import shlex

import pandas as pd
import pytest

from fp_tools.gui_app import GENERIC_TOOL_DEFAULTS
from fp_tools.gui_config import build_cli_command, config_to_yaml_text, expand_jobs, parse_yaml_text
from fp_tools.tools.pseudobulk_footprints import build_parser, run_pseudobulk_footprints


@pytest.mark.parametrize("custom,database,expected", [
    (False, None, "jaspar2026_vertebrates"),
    (True, None, None),
    (False, "hocomoco14_core", "hocomoco14_core"),
    (True, "hocomoco14_core", "hocomoco14_core"),
])
def test_sc_motif_selection_matches_generated_commands(tmp_path, custom, database, expected):
    annotations = tmp_path / "annotations.tsv"
    annotations.write_text("barcode\tcell_type\tsnap_cell_type\tumap_1\tumap_2\n")
    command = ["--fragments", "fragments.tsv", "--annotations", str(annotations),
               "--h5ad", "cells.h5ad", "--genome", "genome.fa", "--genome-sizes", "genome.sizes",
               "--peaks", "peaks.bed", "--group-by", "cell_type", "--outdir", str(tmp_path / "output"), "--dry-run"]
    if custom:
        command.extend(["--motifs", "custom.jaspar"])
    if database:
        command.extend(["--motif-db", database])
    args = build_parser().parse_args(command)
    with patch("fp_tools.tools.pseudobulk_footprints._filter_bed_by_chroms", return_value="peaks.bed"), patch(
        "fp_tools.tools.pseudobulk_footprints._group_inputs", return_value=pd.DataFrame([
            {"group": "B", "passes_filters": True, "source_type": "fragments", "fragment_file": "B.tsv"}
        ])
    ), patch("fp_tools.tools.pseudobulk_footprints.resolve_motif_inputs", return_value=["resolved.meme"]) as resolve:
        assert run_pseudobulk_footprints(args) == 0
        resolve.assert_called_once_with(["custom.jaspar"] if custom else None, expected, use_default=False)
    lines = (tmp_path / "output" / "pseudobulk_footprint_commands.sh").read_text().splitlines()
    diff = next(shlex.split(line) for line in lines if line.startswith("diff-footprints "))
    assert ("--motifs" in diff) == custom
    assert ("--motif-db" in diff) == bool(expected)
    if expected:
        assert diff[diff.index("--motif-db") + 1] == expected
    assert GENERIC_TOOL_DEFAULTS["sc-footprinting"]["motif_db"] == ""
    assert "motifs" in GENERIC_TOOL_DEFAULTS["sc-footprinting"]


@pytest.mark.parametrize("motifs,database", [(None, None), (["custom.jaspar"], None),
                                           (None, "hocomoco14_core"), (["custom.jaspar"], "hocomoco14_core")])
def test_sc_yaml_roundtrip_preserves_explicit_motif_choices(motifs, database):
    item = {"tool": "sc-footprinting", "fragments": "fragments.tsv", "annotations": "annotations.tsv",
            "h5ad": "cells.h5ad", "group_by": "cell_type", "genome": "genome.fa",
            "genome_sizes": "genome.sizes", "peaks": "peaks.bed", "outdir": "output"}
    defaults = {"motifs": motifs, "motif_db": database}
    config = {"version": 1, "defaults": defaults, "samples": [item]}
    restored = parse_yaml_text(config_to_yaml_text(config))
    job = expand_jobs(restored)[0]
    command = build_cli_command(job.tool, job.params)
    assert ("--motifs" in command) == bool(motifs)
    assert ("--motif-db" in command) == bool(database)
