from copy import deepcopy
import pytest
from streamlit.testing.v1 import AppTest


def config_page():
    import streamlit as st
    from pathlib import Path
    from fp_tools import gui_app
    from io import BytesIO
    from unittest.mock import patch

    if "test_upload" in st.session_state:
        with patch.object(st, "file_uploader", return_value=BytesIO(st.session_state.test_upload)):
            gui_app._render_config_page(Path(st.session_state.test_run_dir))
    else:
        gui_app._render_config_page(Path(st.session_state.test_run_dir))


def page(tmp_path):
    signal = tmp_path / "signal.bw"
    regions = tmp_path / "regions.bed"
    signal.touch()
    regions.touch()
    app = AppTest.from_function(config_page)
    app.session_state.test_run_dir = str(tmp_path)
    app.session_state.current_config = {
        "samples": [{"tool": "call-footprints", "signal": str(signal),
                     "regions": str(regions), "output": str(tmp_path / "out.bw")}],
    }
    return app.run()


@pytest.mark.parametrize("text", ["not_a_mapping", "[one, two]", "samples: [", "samples: [1]",
                                  "samples: 5", "version: []\nsamples: [{}]", "false"])
def test_invalid_yaml_preserves_config_and_controls(tmp_path, text):
    app = page(tmp_path)
    assert not app.exception
    original = deepcopy(app.session_state.current_config)
    assert not app.button(key="run_config").disabled
    app.text_area[0].set_value(text)
    app.button(key="config_apply_text").click().run()
    assert not app.exception
    assert app.error
    assert app.session_state.current_config == original
    assert app.button(key="config_save_btn")
    assert app.button(key="run_config").disabled
    from fp_tools.gui_config import config_to_yaml_text
    app.text_area[0].set_value(config_to_yaml_text(original))
    app.button(key="config_apply_text").click().run()
    assert not app.exception
    assert not app.button(key="run_config").disabled


def test_missing_path_preserves_config_and_controls(tmp_path):
    app = page(tmp_path)
    original = deepcopy(app.session_state.current_config)
    app.text_input(key="config_load_path").set_value(str(tmp_path / "missing.yml"))
    app.button(key="config_load_path_btn").click().run()
    assert not app.exception
    assert app.error
    assert app.session_state.current_config == original
    assert app.button(key="config_save_btn")
    assert app.button(key="run_config").disabled


@pytest.mark.parametrize("contents", [b"\xff", b"not_a_mapping", b"samples: ["])
def test_bad_upload_preserves_config_and_controls(tmp_path, contents):
    app = page(tmp_path)
    original = deepcopy(app.session_state.current_config)
    app.session_state.test_upload = contents
    app.run()
    app.button(key="config_apply_upload").click().run()
    assert not app.exception
    assert app.error
    assert app.session_state.current_config == original
    assert app.button(key="config_save_btn")
    assert app.button(key="run_config").disabled


def test_unexpected_loader_error_is_not_hidden(monkeypatch):
    from fp_tools import gui_app

    def broken_loader():
        raise RuntimeError("unexpected implementation error")

    with pytest.raises(RuntimeError, match="unexpected implementation"):
        gui_app._apply_config_input(broken_loader)
