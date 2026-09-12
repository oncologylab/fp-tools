from fp_tools.gui_app import GUI_FIELD_LABELS
from fp_tools.gui_config import validate_config


def test_desktop_audit_uses_current_comparisons_label():
    import ast
    from pathlib import Path

    path = Path(__file__).resolve().parents[1] / "scripts/audit_desktop_gui.py"
    tree = ast.parse(path.read_text(encoding="utf-8"))
    labels = [
        node.value for node in ast.walk(tree)
        if isinstance(node, ast.Constant) and isinstance(node.value, str)
        and node.value.startswith("Comparisons TSV")
    ]
    assert labels, "The desktop audit must exercise the comparisons input"
    assert set(labels) == {GUI_FIELD_LABELS["comparison_table"]}


def test_comparisons_label_agrees_with_required_validation():
    assert GUI_FIELD_LABELS["comparison_table"] == "Comparisons TSV"
    config = {"samples": [{"tool": "bulk-footprinting", "sample_table": "samples.tsv",
                           "genome": "hg38", "outdir": "results"}]}
    assert any("comparison_table" in error for error in validate_config(config))
    config["samples"][0]["comparison_table"] = "comparisons.tsv"
    assert validate_config(config) == []
