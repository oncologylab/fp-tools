# [`run-yaml-workflow`](../../api.md#run-yaml-workflow)

Run one or more fp-tools commands from a saved YAML configuration. Use this to
repeat a GUI run or apply the same settings to several samples or comparisons.

## Example command

```bash
run-yaml-workflow --config workflow.yml --dry-run
run-yaml-workflow --config workflow.yml --run-root project/yaml_runs
```

## Primary inputs

- `--config` — YAML configuration exported by the GUI or written directly.
- `--run-root` — optional directory for job status and logs; it does not replace each command's analysis output directory.
- `--dry-run` — print the commands without starting the analysis.

Replace `workflow.yml` with your saved configuration. The first command prints
the jobs without running them; check their inputs and output paths before
running the second command. Jobs run sequentially.

## Main outputs

- The analysis files documented for each command named in the YAML.
- Standard output containing the expanded command lines when `--dry-run` is used.
- `{run_root}/{job_id}/config.yml` and `command.txt` — saved settings and exact command for each job.
- `{run_root}/{job_id}/status.json`, `stdout.log`, and `stderr.log` — completion state and captured command output.
- `{run_root}/batch_index.tsv` — one-row-per-job batch status index.

Relative paths in the YAML are interpreted from the directory where you launch
the command. Run from the same directory each time, or use absolute paths. If
neither `--run-root` nor a YAML `run_root` is supplied, logs go to a new
`fp-tools-batch-{timestamp}` folder in that directory.

Open `batch_index.tsv` to check which jobs succeeded. For a failed job, inspect
its `stderr.log`. Add `--fail-fast` to stop after the first failed job. See the
[complete `run-yaml-workflow` reference](../../api.md#run-yaml-workflow).
