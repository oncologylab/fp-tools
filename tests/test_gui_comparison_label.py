from fp_tools.gui_app import GUI_FIELD_LABELS
from fp_tools.gui_config import validate_config


def test_comparisons_label_agrees_with_required_validation():
    assert GUI_FIELD_LABELS["comparison_table"] == "Comparisons TSV"
    config = {"samples": [{"tool": "bulk-footprinting", "sample_table": "samples.tsv",
                           "genome": "hg38", "outdir": "results"}]}
    assert any("comparison_table" in error for error in validate_config(config))
    config["samples"][0]["comparison_table"] = "comparisons.tsv"
    assert validate_config(config) == []
