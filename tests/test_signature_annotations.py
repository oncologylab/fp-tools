from argparse import Namespace
from unittest.mock import patch

import pytest

from fp_tools.gui_config import _validate_gui_input_paths
from fp_tools.tools.find_signature_fp import read_annotations
from fp_tools.tools.pseudobulk_footprints import run_pseudobulk_footprints


@pytest.mark.parametrize("separator,suffix", [("\t", ".tsv"), (",", ".csv")])
def test_missing_signature_annotation_rejected_by_reader_and_gui(tmp_path, separator, suffix):
    path = tmp_path / ("annotations" + suffix)
    path.write_text(separator.join(["barcode", "cell_type", "umap_1", "umap_2"]) + "\n")
    with pytest.raises(SystemExit, match="snap_cell_type"):
        read_annotations(path)
    for tool in ("find-signature-fp", "sc-footprinting"):
        assert any("snap_cell_type" in error for error in _validate_gui_input_paths(tool, {"annotations": str(path)}, "test"))
    path.write_text(separator.join(["barcode", "cell_type", "snap_cell_type", "umap_1", "umap_2"]) + "\n")
    assert "snap_cell_type" in read_annotations(path).columns
    assert not _validate_gui_input_paths("find-signature-fp", {"annotations": str(path)}, "test")


def test_sc_annotation_preflight_precedes_output_and_grouping(tmp_path):
    path = tmp_path / "annotations.tsv"
    path.write_text("barcode\tcell_type\tumap_1\tumap_2\n")
    output = tmp_path / "output"
    args = Namespace(annotations=str(path), outdir=str(output), single_cell_signature_script=None)
    with patch("fp_tools.tools.pseudobulk_footprints._group_inputs") as group:
        with pytest.raises(SystemExit, match="snap_cell_type"):
            run_pseudobulk_footprints(args)
        group.assert_not_called()
    assert not output.exists()


def test_signature_cli_rejects_before_output_creation(tmp_path):
    from fp_tools.tools.find_signature_fp import main

    path = tmp_path / "annotations.tsv"
    path.write_text("barcode\tcell_type\tumap_1\tumap_2\n")
    output = tmp_path / "output"
    with pytest.raises(SystemExit, match="snap_cell_type"):
        main(["--annotations", str(path), "--fragments", "absent.tsv", "--h5ad", "absent.h5ad", "--outdir", str(output)])
    assert not output.exists()
