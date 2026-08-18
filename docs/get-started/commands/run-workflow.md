# `run-workflow`

Run one or more fp-tools jobs from a reusable YAML configuration.

## Example command

```bash
run-workflow \
  --config examples/gui_configs/diff_footprints_single.yml
```

## Primary inputs

- A YAML configuration exported by the GUI or written directly.
- Optional tool filters and run-metadata directory.

## Main outputs

- The same outputs produced by the configured fp-tools commands.
- Expanded-command preview in dry-run mode.
- Run metadata and logs when a run root is supplied.

The YAML remains command-compatible and does not require GUI state. See the
[complete `run-workflow` reference](../../api.md#run-workflow).
