import io
import sys
import subprocess
from unittest import mock

import pytest

from fp_tools.tools import bulk_footprinting, pseudobulk_footprints


def test_wrapper_defaults_are_automatic():
    for module in (bulk_footprinting, pseudobulk_footprints):
        parser = module.build_parser()
        assert parser.get_default("cores") is None


def test_bulk_child_output_reaches_console_and_logs(tmp_path, capsys):
    stdout = tmp_path / "stdout.log"
    stderr = tmp_path / "stderr.log"
    code = bulk_footprinting._run(
        [sys.executable, "-c", "import sys; print('stage output'); print('stage error', file=sys.stderr); sys.exit(7)"],
        stdout, stderr,
    )
    captured = capsys.readouterr()
    assert code == 7
    assert "stage output" in captured.out
    assert "stage error" in captured.err
    assert str(stderr) in captured.err
    assert "stage output" in stdout.read_text()
    assert "stage error" in stderr.read_text()


def test_bulk_commands_resolve_auto_budget():
    parser = bulk_footprinting.build_parser()
    args = parser.parse_args([
        "--sample-table", "samples.tsv", "--comparison-table", "comparisons.tsv",
        "--genome", "hg38", "--outdir", "project",
    ])
    from fp_tools.utils import resources
    with mock.patch.object(resources, "available_cores", return_value=36):
        commands = bulk_footprinting.build_commands(args)
    for _, command in commands:
        if "--cores" in command:
            assert command[command.index("--cores") + 1] == "36"
    assert args.cores is None


def test_core_budget_respects_affinity_and_overrides():
    from fp_tools.utils import resources
    with mock.patch.object(resources.os, "process_cpu_count", return_value=None, create=True), \
         mock.patch.object(resources.os, "sched_getaffinity", return_value={0, 1, 2, 3}, create=True):
        assert resources.available_cores() == 4
        assert resources.resolve_cores(None) == 4
        assert resources.resolve_cores(2) == 2
        warnings = []
        assert resources.resolve_cores(36, warn=warnings.append) == 4
        assert resources.resolve_cores(0, warn=warnings.append) == 4
        assert len(warnings) == 2


def test_core_budget_falls_back_to_one():
    from fp_tools.utils import resources
    with mock.patch.object(resources.os, "process_cpu_count", return_value=None, create=True), \
         mock.patch.object(resources.os, "sched_getaffinity", side_effect=OSError, create=True), \
         mock.patch.object(resources.os, "cpu_count", return_value=None):
        assert resources.available_cores() == 1


def test_output_is_visible_before_child_exits(tmp_path, monkeypatch):
    from fp_tools.utils.workflow_execution import run_logged
    release = tmp_path / "release"

    class Console(io.StringIO):
        def write(self, text):
            result = super().write(text)
            if "ready" in self.getvalue():
                release.touch()
            return result

    console = Console()
    monkeypatch.setattr(sys, "stdout", console)
    script = (
        "import os,time,sys; from pathlib import Path; "
        "os.write(1,b'ready'); deadline=time.monotonic()+10; "
        "p=Path(sys.argv[1])\n"
        "while not p.exists() and time.monotonic()<deadline: time.sleep(.01)\n"
        "sys.exit(0 if p.exists() else 9)"
    )
    result = run_logged([sys.executable, "-c", script, str(release)])
    assert result.returncode == 0
    assert console.getvalue() == "ready"


def test_simultaneous_large_streams_are_preserved(tmp_path, capsys):
    from fp_tools.utils.workflow_execution import run_logged
    script = (
        "import os,threading\n"
        "def emit(fd, byte):\n"
        " for _ in range(128): os.write(fd, byte*8192)\n"
        "t=threading.Thread(target=emit,args=(2,b'E')); t.start(); emit(1,b'O'); t.join()"
    )
    out, err = tmp_path / "out", tmp_path / "err"
    result = run_logged([sys.executable, "-c", script], stdout_log=out, stderr_log=err)
    captured = capsys.readouterr()
    assert result.returncode == 0
    assert out.read_bytes() == b"O" * 1048576
    assert err.read_bytes() == b"E" * 1048576
    assert captured.out == "O" * 1048576
    assert "E" * 1048576 in captured.err


def test_unicode_partial_and_invalid_bytes(tmp_path, capsys):
    from fp_tools.utils.workflow_execution import run_logged
    out = tmp_path / "out"
    script = "import os,time; os.write(1,b'\\xc3'); time.sleep(.05); os.write(1,b'\\xa9\\xff')"
    assert run_logged([sys.executable, "-c", script], stdout_log=out).returncode == 0
    assert out.read_bytes() == b"\xc3\xa9\xff"
    assert capsys.readouterr().out == "é\ufffd"


def test_pipeline_data_is_not_echoed(tmp_path, capsys):
    from fp_tools.utils.workflow_execution import run_logged
    first = [sys.executable, "-c", "import os; os.write(1,bytes(range(256))); os.write(2,b'producer note\\n')"]
    second = [sys.executable, "-c", "import sys; b=sys.stdin.buffer.read(); print(len(b))"]
    result = run_logged(second, input_command=first, capture_stdout=True, stderr_log=tmp_path / "err")
    captured = capsys.readouterr()
    assert result.returncode == 0
    assert result.stdout.strip() == b"256"
    assert captured.out == ""
    assert "producer note" in captured.err


def test_binary_file_output_is_preserved(tmp_path, capsys):
    from fp_tools.utils.workflow_execution import run_logged
    output = tmp_path / "data"
    with output.open("wb") as handle:
        result = run_logged([sys.executable, "-c", "import os; os.write(1,bytes(range(256)))"], stdout_target=handle)
    assert result.returncode == 0
    assert output.read_bytes() == bytes(range(256))
    assert capsys.readouterr().out == ""


def test_interruption_terminates_child(tmp_path, monkeypatch):
    from fp_tools.utils import workflow_execution
    original = subprocess.Popen.wait
    interrupted = []

    def wait(process, *args, **kwargs):
        if not interrupted:
            interrupted.append(process)
            raise KeyboardInterrupt
        return original(process, *args, **kwargs)

    monkeypatch.setattr(subprocess.Popen, "wait", wait)
    with pytest.raises(KeyboardInterrupt):
        workflow_execution.run_logged([sys.executable, "-c", "import time; time.sleep(60)"], stdout_log=tmp_path / "out")
    assert interrupted[0].poll() is not None


def test_startup_failure_records_failed_job(tmp_path, monkeypatch):
    import json
    from fp_tools import cli_batch
    from fp_tools.gui_config import expand_jobs
    job = expand_jobs({"samples": [{"tool": "call-footprints", "job_id": "broken"}]})[0]
    monkeypatch.setattr(cli_batch, "resolve_fp_tools_subprocess", lambda _: [str(tmp_path / "missing")])
    assert cli_batch.run_job(job, tmp_path) == 1
    status = json.loads((tmp_path / "broken" / "status.json").read_text())
    assert status["status"] == "failed"
    assert status["exit_code"] == 1


def test_sample_concurrency_never_exceeds_budget():
    from fp_tools.tools import atacorrect, score_bigwig, diff_footprints
    from fp_tools.utils import resources
    with mock.patch.object(resources, "available_cores", return_value=4):
        for module in (atacorrect, score_bigwig, diff_footprints):
            assert module._sample_worker_plan(12, None, 8) == (4, 1)
            assert module._sample_worker_plan(12, 2, 8) == (2, 1)


def test_gui_automatic_override_survives_roundtrip():
    from fp_tools import gui_app
    from fp_tools.gui_config import config_to_yaml_text, parse_yaml_text, expand_jobs
    config = {"defaults": {"cores": 2}, "samples": [{"tool": "bulk-footprinting", "cores": None}]}
    restored = parse_yaml_text(config_to_yaml_text(config))
    assert "--cores" not in expand_jobs(restored)[0].command
    restored["samples"][0]["cores"] = 3
    command = expand_jobs(restored)[0].command
    assert command[command.index("--cores") + 1] == "3"
    for tool in ("bulk-footprinting", "sc-footprinting"):
        assert gui_app.GENERIC_TOOL_DEFAULTS[tool]["cores"] is None


@pytest.mark.parametrize("assembly", ["hg38", "mm10"])
def test_managed_assembly_is_propagated_to_bulk_stages(tmp_path, monkeypatch, capsys, assembly):
    from fp_tools.utils import references
    samples, comparisons = tmp_path / "samples.tsv", tmp_path / "comparisons.tsv"
    samples.write_text("sample\tcondition\tbam\tpeaks\nA1\tA\tA.bam\tA.bed\nB1\tB\tB.bam\tB.bed\n")
    comparisons.write_text("comparison\tcond1\tcond2\nA_B\tA\tB\n")
    cache = tmp_path / "cache"
    args = bulk_footprinting.build_parser().parse_args([
        "--sample-table", str(samples), "--comparison-table", str(comparisons),
        "--genome", assembly, "--reference-dir", str(cache),
        "--outdir", str(tmp_path / "project"), "--dry-run",
    ])
    with mock.patch.object(references, "verified_urlopen") as download:
        assert bulk_footprinting.run_bulk_footprinting(args) == 0
        download.assert_not_called()
    output = capsys.readouterr().out
    assert str(cache / assembly / f"{assembly}.fa") in output
    assert str(cache / assembly / f"{assembly}.blacklist.bed") in output
    assert not cache.exists()
    assert not (tmp_path / "project").exists()
    manifest = references.REFERENCE_MANIFEST[assembly]
    assert f"/{assembly}/" in manifest["fasta_url"]
    assert manifest["blacklist_url"].endswith(f"/{assembly}-blacklist.v2.bed.gz")


def test_preprocessing_config_keeps_auto_portable(tmp_path):
    from fp_tools.tools.prepare_atac import write_default_config, load_settings
    from fp_tools.utils import resources
    import yaml
    path = write_default_config(tmp_path / "prepare.yml")
    assert yaml.safe_load(path.read_text())["resources"]["cores"] is None
    with mock.patch.object(resources, "available_cores", return_value=36):
        assert load_settings(path)["resources"]["cores"] == 36


def test_pseudobulk_generated_commands_resolve_automatic_cores(tmp_path):
    import pandas as pd
    from fp_tools.tools.pseudobulk import write_downstream_commands
    from fp_tools.utils import resources
    with mock.patch.object(resources, "available_cores", return_value=36):
        path = write_downstream_commands(pd.DataFrame([
            {"group": "A", "fragment_file": "A.fragments.tsv"},
        ]), tmp_path / "commands.sh")
    assert "samtools sort -@ 36" in path.read_text()


def test_core_form_can_clear_a_saved_limit():
    from streamlit.testing.v1 import AppTest

    script = '''
import streamlit as st
from fp_tools.gui_app import _render_core_limit
with st.form("settings"):
    cores = _render_core_limit(2, "cores_test")
    if st.form_submit_button("Apply"):
        st.session_state.saved_cores = cores
'''

    original = sys.modules["__main__"]
    try:
        app = AppTest.from_string(script, default_timeout=30).run()
        assert not app.exception
        assert app.number_input[0].value == 2
        app.number_input[0].set_value(None)
        app.button[0].click().run()
        assert not app.exception
        assert app.session_state.saved_cores is None
        app.number_input[0].set_value(3)
        app.button[0].click().run()
        assert not app.exception
        assert app.session_state.saved_cores == 3
    finally:
        sys.modules["__main__"] = original


def test_nested_wrapper_relays_child_once(tmp_path, capsys):
    from fp_tools.utils.workflow_execution import run_logged
    inner_out, inner_err = tmp_path / "inner.out", tmp_path / "inner.err"
    script = (
        "import sys; from fp_tools.utils.workflow_execution import run_logged; "
        f"run_logged([sys.executable,'-c',\"print('nested result')\"], stdout_log={str(inner_out)!r}, stderr_log={str(inner_err)!r})"
    )
    result = run_logged([sys.executable, "-c", script], stdout_log=tmp_path / "outer.out")
    assert result.returncode == 0
    assert capsys.readouterr().out == "nested result\n"
    assert inner_out.read_text() == "nested result\n"
    assert (tmp_path / "outer.out").read_text() == "nested result\n"


def test_batch_core_columns_keep_integer_limits_and_null_defaults(monkeypatch):
    import pandas as pd
    from fp_tools import gui_app
    frame = pd.DataFrame([{"sample_id": "A", "cores": 2.0}, {"sample_id": "B", "cores": None}, {"sample_id": None, "cores": None}])
    monkeypatch.setattr(gui_app.st, "subheader", lambda _: None)
    monkeypatch.setattr(gui_app.st, "data_editor", lambda *a, **k: frame)
    monkeypatch.setattr(gui_app, "_config_widget_key", lambda key: key)
    rows = gui_app._data_editor("samples", frame.to_dict("records"), "test")
    assert type(rows[0]["cores"]) is int
    assert rows[1]["cores"] is None
    assert len(rows) == 2
