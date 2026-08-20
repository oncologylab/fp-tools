# [`fp-tools-gui`](../../api.md#fp-tools-gui)

Launch the browser interface for configuring and running fp-tools commands.

Bulk GUI workflows start from coordinate-sorted BAM/BAI files and matching
peak BED files. The GUI does not perform FASTQ-to-BAM preprocessing.
Missing inputs and unsupported options are reported before a run starts.

The GUI is available through the Python package, the complete container, and
the self-contained desktop downloads on the
[release page](https://github.com/oncologylab/fp-tools/releases).

## Example command

```bash
fp-tools-gui --host 127.0.0.1 --port 8891 --run-dir project/gui_runs --no-browser
```

## Primary inputs

- `--host` — interface on which the GUI listens (default: `127.0.0.1`).
- `--port` — fixed browser port.
- `--run-dir` — directory for GUI-managed configurations and runs.
- `--no-browser` — start the server without opening a local browser.

## Main outputs

- `{run_dir}/{timestamp}_{label}/config.yml` — reusable command-compatible YAML saved for a configured run.
- `{run_dir}/{timestamp}_{label}/status.json`, `launcher_stdout.log`, and `launcher_stderr.log` — launcher state and captured batch-runner output.
- `{run_dir}/{timestamp}_{label}/{job_id}/status.json`, `command.txt`, `stdout.log`, and `stderr.log` — per-job state, exact command, and analysis logs.
- The exact analysis files documented by the selected command; the GUI does not introduce GUI-only scientific outputs.

Files under `{run_dir}` are local run state. A saved YAML remains runnable with
`run-yaml-workflow --config {run_dir}/{timestamp}_{label}/config.yml`.

Open the [GUI Demo](../../gui.md), or see the
[complete `fp-tools-gui` reference](../../api.md#fp-tools-gui).

## Local computer

Open the desktop executable or run `fp-tools-gui`. A browser opens after the
server is ready. If it does not, open the local URL printed in the terminal.

## Remote Linux server

Start fp-tools on the server without exposing a network port:

```bash
fp-tools-gui --no-browser --port 8891
```

On your computer, create an SSH tunnel and keep that terminal open:

```bash
ssh -N -L 8891:127.0.0.1:8891 USER@SERVER
```

Open `http://127.0.0.1:8891`. Binding with `--host 0.0.0.0` is also supported,
but fp-tools does not add authentication; protect direct network access with a
firewall, VPN, or reverse proxy.
