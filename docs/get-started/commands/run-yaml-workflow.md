# [`run-yaml-workflow`](../../api.md#run-yaml-workflow)

Run one or more fp-tools jobs from a reusable YAML configuration.

## Example command

```bash
run-yaml-workflow --config examples/gui_configs/diff_footprints_single.yml
```

## Primary inputs

- `--config` — command-compatible YAML configuration exported by the GUI or written directly.

## Main outputs

- The exact files documented for each command named in the YAML; YAML does not create a separate analysis format.
- Standard output containing the expanded command lines when `--dry-run` is used.
- `{run_root}/{job_id}/config.yml` and `command.txt` — normalized per-job configuration and exact command.
- `{run_root}/{job_id}/status.json`, `stdout.log`, and `stderr.log` — completion state and captured command output.
- `{run_root}/batch_index.tsv` — one-row-per-job batch status index.

Paths are resolved according to the YAML runner and remain independent of GUI
state. Inspect the dry-run expansion before starting a long workflow.

The YAML remains command-compatible and does not require GUI state. See the
[complete `run-yaml-workflow` reference](../../api.md#run-yaml-workflow).
