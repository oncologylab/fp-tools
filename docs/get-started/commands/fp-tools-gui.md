# [`fp-tools-gui`](../../api.md#fp-tools-gui)

Launch the optional browser interface for configuring and running fp-tools
commands.

## Example command

```bash
fp-tools-gui \
  --host 0.0.0.0 \
  --port 8891 \
  --run-dir project/gui_runs
```

## Primary inputs

- `--host` — interface on which the GUI listens.
- `--port` — fixed browser port.
- `--run-dir` — directory for GUI-managed configurations and runs.

## Main outputs

- `{run_dir}/{timestamp}_{label}/config.yml` — reusable command-compatible YAML saved for a configured run.
- `{run_dir}/{timestamp}_{label}/status.json`, `launcher_stdout.log`, and `launcher_stderr.log` — launcher state and captured batch-runner output.
- `{run_dir}/{timestamp}_{label}/{job_id}/status.json`, `command.txt`, `stdout.log`, and `stderr.log` — per-job state, exact command, and analysis logs.
- The exact analysis files documented by the selected command; the GUI does not introduce GUI-only scientific outputs.

Files under `{run_dir}` are local run state. A saved YAML remains runnable with
`run-yaml-workflow --config {run_dir}/{timestamp}_{label}/config.yml`.

Open the [GUI Demo](../../gui.md), or see the
[complete `fp-tools-gui` reference](../../api.md#fp-tools-gui).
