from copy import deepcopy

import pytest

from fp_tools.gui_config import expand_jobs, validate_config, validate_gui_config


@pytest.mark.parametrize("section", ["samples", "comparisons"])
def test_defaults_validate_like_expanded_jobs(tmp_path, section):
    signal = tmp_path / "signal.bw"
    regions = tmp_path / "regions.bed"
    signal.touch()
    regions.touch()
    config = {
        "defaults": {"tool": "call-footprints", "signal": str(signal),
                     "regions": str(regions), "output": str(tmp_path / "out.bw")},
        section: [{"job_id": "inherited"}],
    }
    original = deepcopy(config)
    assert validate_config(config) == []
    assert validate_gui_config(config) == []
    assert expand_jobs(config)[0].params["signal"] == str(signal)
    assert config == original
    config[section][0]["signal"] = ""
    assert any("signal" in error for error in validate_config(config))
    assert expand_jobs(config)[0].params["signal"] == ""


@pytest.mark.parametrize("field,value,expected", [
    ("signals", ["/nonexistent/defaults-signal.bw"], "does not exist"),
    ("score", "invalid-score", "unsupported 'score'"),
])
def test_invalid_inherited_gui_fields(field, value, expected):
    config = {"defaults": {field: value}, "samples": [{"tool": "call-footprints"}]}
    assert any(expected in error for error in validate_gui_config(config))


def test_inherited_gui_boundary_and_numeric_validation():
    config = {"defaults": {"tool": "bulk-footprinting", "reads_table": "reads.tsv"},
              "samples": [{"sample_id": "one"}]}
    assert any("FASTQ preparation" in error for error in validate_gui_config(config))
    config = {"defaults": {"tool": "normalize-bigwig", "stat": "q101"},
              "comparisons": [{"comparison_id": "one"}]}
    assert any("unsupported 'stat'" in error for error in validate_gui_config(config))
