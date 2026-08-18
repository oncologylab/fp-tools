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

- Reusable YAML configurations.
- Runs and outputs from the corresponding command implementation.
- Browser views for monitoring configured jobs.

Open the [GUI Demo](../../gui.md), or see the
[complete `fp-tools-gui` reference](../../api.md#fp-tools-gui).
