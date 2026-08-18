# `fp-tools-gui`

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

- Optional host and port settings.
- A directory for GUI-managed runs.
- User-selected command inputs and parameters in the browser.

## Main outputs

- Reusable YAML configurations.
- Runs and outputs from the corresponding command implementation.
- Browser views for monitoring configured jobs.

Open the [GUI Demo](../../gui.md), or see the
[complete `fp-tools-gui` reference](../../api.md#fp-tools-gui).
