from hashlib import sha256
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

import pytest

if sys.platform != 'linux':
    pytest.skip('research process supervision requires Linux /proc', allow_module_level=True)

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'benchmarks/scripts'))
import run_logged_experiment as m


def make_spec(tmp_path, code="print('completed fixture')", **updates):
    spec = {'schema': m.SCHEMA, 'name': 'fixture', 'command': [sys.executable, '-c', code],
            'cwd': str(tmp_path), 'outdir': str(tmp_path / 'execution'),
            'environment': {'OPENBLAS_NUM_THREADS': '1'}}
    spec.update(updates)
    path = tmp_path / 'job.json'
    path.write_text(json.dumps(spec))
    return path, sha256(path.read_bytes()).hexdigest()


def command(path, digest):
    return [sys.executable, str(Path(m.__file__)), '--spec', str(path), '--spec-sha256', digest]


def verify_receipt(path, field):
    document = json.loads(path.read_text())
    value = {key: item for key, item in document.items() if key != field}
    assert sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest() == document[field]
    return document


@pytest.mark.parametrize('returncode', [0, 7])
def test_real_child_exit_log_and_identities_are_preserved(tmp_path, returncode):
    path, digest = make_spec(tmp_path, f"print('literal log'); raise SystemExit({returncode})")
    done = subprocess.run(command(path, digest), capture_output=True, text=True)
    assert done.returncode == returncode, done.stderr
    out = tmp_path / 'execution'
    launch = verify_receipt(out / 'launch.json', 'launch_id')
    result = verify_receipt(out / 'exit.json', 'exit_id')
    assert launch['child']['pid'] != launch['supervisor']['pid']
    assert launch['child']['start_ticks'] > 0 and launch['child']['boot_id'] == launch['supervisor']['boot_id']
    assert launch['specification'] == m.file_record(out / 'specification.json')
    assert (out / 'specification.json').read_bytes() == path.read_bytes()
    assert launch['specification']['sha256'] == digest
    assert result['returncode'] == returncode and result['scientific_outputs_verified'] is False
    assert result['signals_received_by_supervisor'] == []
    assert result['stdout'] == m.file_record(out / 'stdout.log')
    assert (out / 'stdout.log').read_text() == 'literal log\n'


def test_specification_corruption_rejects_before_output_creation(tmp_path):
    path, digest = make_spec(tmp_path)
    path.write_text(path.read_text() + ' ')
    with pytest.raises(ValueError, match='checksum mismatch'):
        m.supervise(path, digest)
    assert not (tmp_path / 'execution').exists()


@pytest.mark.parametrize('updates', [
    {'environment': {'PRIVATE_TOKEN': 'not-a-secret-fixture'}},
    {'command': 'python -c malicious-shell-string'},
    {'command': ['python', '-c', 'pass']},
    {'outdir': 'relative-output'},
])
def test_invalid_specification_rejects_before_output_creation(tmp_path, updates):
    path, digest = make_spec(tmp_path, **updates)
    with pytest.raises(ValueError):
        m.supervise(path, digest)
    assert not (tmp_path / 'execution').exists()


def test_existing_execution_is_not_repeated_or_overwritten(tmp_path):
    path, digest = make_spec(tmp_path, "from pathlib import Path; p=Path('counter'); p.write_text(p.read_text()+'x' if p.exists() else 'x')")
    assert subprocess.run(command(path, digest), capture_output=True).returncode == 0
    out = tmp_path / 'execution'
    before = {p.name: p.read_bytes() for p in out.iterdir()}
    assert subprocess.run(command(path, digest), capture_output=True).returncode != 0
    assert {p.name: p.read_bytes() for p in out.iterdir()} == before
    assert (tmp_path / 'counter').read_text() == 'x'


def test_missing_executable_records_spawn_failure(tmp_path):
    path, digest = make_spec(tmp_path, command=[str(tmp_path / 'missing-program')])
    done = subprocess.run(command(path, digest), capture_output=True)
    assert done.returncode == 127
    launch = verify_receipt(tmp_path / 'execution/launch.json', 'launch_id')
    result = verify_receipt(tmp_path / 'execution/exit.json', 'exit_id')
    assert launch['child'] is None and launch['spawn_error']
    assert result['returncode'] == 127 and result['stdout']['size_bytes'] == 0


def test_sigterm_reaches_owned_child_and_has_a_durable_exit_record(tmp_path):
    path, digest = make_spec(tmp_path, 'import time; time.sleep(60)')
    supervisor = subprocess.Popen(command(path, digest), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                  start_new_session=True)
    launch_path = tmp_path / 'execution/launch.json'
    try:
        deadline = time.monotonic() + 10
        while not launch_path.exists() and time.monotonic() < deadline:
            assert supervisor.poll() is None
            time.sleep(.02)
        launch = verify_receipt(launch_path, 'launch_id')
        supervisor.send_signal(signal.SIGTERM)
        assert supervisor.wait(timeout=10) == 143
        result = verify_receipt(tmp_path / 'execution/exit.json', 'exit_id')
        assert result['returncode'] == -signal.SIGTERM
        assert [item['signal'] for item in result['signals_received_by_supervisor']] == ['SIGTERM']
        assert not Path(f"/proc/{launch['child']['pid']}").exists()
    finally:
        if supervisor.poll() is None:
            supervisor.kill()
            supervisor.wait()
        if launch_path.exists():
            identity = json.loads(launch_path.read_text())['child']
            try:
                if m.process_identity(identity['pid']) == identity:
                    os.killpg(identity['pid'], signal.SIGKILL)
            except (ProcessLookupError, FileNotFoundError):
                pass


def test_receipt_write_failure_stops_the_unrecorded_child(tmp_path, monkeypatch):
    path, digest = make_spec(tmp_path, 'import time; time.sleep(60)')
    children = []
    original = m.subprocess.Popen

    def capture(*args, **kwargs):
        child = original(*args, **kwargs)
        children.append(child)
        return child

    def failed_write(*args, **kwargs):
        raise OSError('simulated receipt storage failure')

    monkeypatch.setattr(m.subprocess, 'Popen', capture)
    monkeypatch.setattr(m, 'write_receipt', failed_write)
    with pytest.raises(OSError, match='simulated receipt'):
        m.supervise(path, digest)
    assert len(children) == 1 and children[0].poll() == -signal.SIGTERM


def test_detached_supervisor_survives_launcher_exit_with_a_frozen_specification(tmp_path):
    path, digest = make_spec(tmp_path, "from pathlib import Path; import time\nwhile not Path('release').exists(): time.sleep(.02)\nprint('released')")
    launcher = "import json,subprocess,sys; p=subprocess.Popen(json.loads(sys.argv[1]),start_new_session=True,stdin=subprocess.DEVNULL,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL); print(p.pid)"
    parent = subprocess.run([sys.executable, '-c', launcher, json.dumps(command(path, digest))],
                            capture_output=True, text=True, check=True)
    supervisor_pid = int(parent.stdout.strip())
    out = tmp_path / 'execution'
    launch = None
    try:
        deadline = time.monotonic() + 10
        while not (out / 'launch.json').exists() and time.monotonic() < deadline:
            time.sleep(.02)
        launch = verify_receipt(out / 'launch.json', 'launch_id')
        assert launch['supervisor']['pid'] == supervisor_pid
        assert os.getsid(supervisor_pid) == supervisor_pid
        assert os.getsid(launch['child']['pid']) == launch['child']['pid']
        assert m.process_identity(supervisor_pid) == launch['supervisor']
        path.write_text(path.read_text() + ' ')
        (tmp_path / 'release').touch()
        deadline = time.monotonic() + 10
        while not (out / 'exit.json').exists() and time.monotonic() < deadline:
            time.sleep(.02)
        result = verify_receipt(out / 'exit.json', 'exit_id')
        assert result['returncode'] == 0
        assert launch['specification']['sha256'] == digest
        assert m.file_record(out / 'specification.json')['sha256'] == digest
        assert m.file_record(path)['sha256'] != digest
    finally:
        (tmp_path / 'release').touch()
        if launch is not None:
            try:
                if m.process_identity(supervisor_pid) == launch['supervisor']:
                    os.kill(supervisor_pid, signal.SIGTERM)
            except (ProcessLookupError, FileNotFoundError):
                pass
