# GUI Status

## Goal
Provide an optional browser-based GUI for `fp-tools` that runs as a per-user process on a Linux server without changing the core command-line workflows.

## Current Status
The first working GUI layer is implemented.

Available commands:
- `fp-tools-gui`
- `fp-tools-run --config ...`

Core packaged commands remain primary and unchanged:
- `ATACorrect`
- `FootprintScores`
- `BINDetect`
- `PlotAggregate`

## Implemented

### Architecture
- GUI remains isolated from `src/fp_tools/tools/`.
- Core logic was not moved out of the packaged commands.
- YAML is shared between GUI and optional config-driven CLI.
- Plain CLI does not require YAML.

### Files Added
- `src/fp_tools/gui_app.py`
- `src/fp_tools/gui_jobs.py`
- `src/fp_tools/gui_forms.py`
- `src/fp_tools/gui_config.py`
- `src/fp_tools/cli_batch.py`
- `src/fp_tools/cli_gui.py`

### Config Model
- GUI can run without a preexisting YAML file.
- GUI can load YAML.
- GUI can save YAML.
- Every GUI run materializes a normalized `config.yml`.
- GUI-saved YAML can be rerun with:
  - `fp-tools-run --config <file>.yml`

### Batch Support
- sample-list batch support for:
  - `ATACorrect`
  - `FootprintScores`
  - `PlotAggregate`
  - single-condition `BINDetect`
- comparison-list batch support for:
  - multi-condition `BINDetect`
  - replicate grouping by repeated `cond_names`

### GUI Pages
- `Home`
- `Run History`
- `ATACorrect`
- `FootprintScores`
- `BINDetect`
- `PlotAggregate`
- `Config`

### GUI Behavior
- supports direct form-driven runs
- supports loading example YAML configs
- supports loading uploaded YAML
- supports loading YAML from a path
- supports lightweight pre-launch config validation
- supports background job launch
- supports run-history inspection
- detects primary output paths for completed child jobs
- supports changing the GUI run directory from the sidebar

### Example Assets
Ready-to-load GUI YAML:
- `examples/gui_configs/`

GUI run metadata and logs:
- `examples/gui_runs/`

GUI example output files:
- `examples/gui_demo_outputs/`

### Validation Completed
- local syntax checks for new GUI/config modules
- local `fp-tools-run --help`
- local `fp-tools-gui --help`
- local real YAML-driven `PlotAggregate` run through `fp-tools-run`
- local `fp-tools-gui` launch validation
- remote install on `cy232`
- remote `fp-tools-gui --help`
- remote `fp-tools-run --help`
- remote GUI HTTP response check
- remote example YAML run through `fp-tools-run`

## Current Limitations

### GUI UX
- page forms currently cover the most important options, not every advanced CLI flag
- no inline preview yet for generated PDFs, HTML, or tables
- run history is functional but still basic beyond log/status/output-path inspection
- Home page does not yet provide rich quick links into recent runs

### Execution Model
- GUI jobs are launched in the background, but there is no explicit cancel button yet
- output discovery is still path-based rather than a richer results browser
- top-level run status handling was improved, but the status layer should still be hardened with more edge-case testing

### Validation
- the GUI has been smoke-tested locally and remotely
- not every page and every option combination has been exhaustively exercised through the GUI yet
- most deep validation still exists at the CLI level rather than the GUI level

## Confirmed Design Rules
- one GUI process per Linux user
- no shared multi-user service
- no login in the first version
- direct CLI remains primary
- YAML on CLI is optional, not required
- GUI and config runner must not break existing direct CLI workflows

## Recommended Usage

Launch with an auto-selected port:

```bash
fp-tools-gui --host 0.0.0.0 --run-dir examples/gui_runs
```

Or use a fixed port:

```bash
fp-tools-gui --host 0.0.0.0 --port 8891 --run-dir examples/gui_runs
```

Run a saved YAML config directly:

```bash
fp-tools-run --config examples/gui_configs/plotaggregate_single.yml
```

## Next Tasks

### Highest Priority
1. Expand per-page validation and user feedback beyond the current required-field checks.
2. Improve output browsing from simple detected paths into a richer results view.
3. Harden run-status refresh and batch-status aggregation.

### Medium Priority
1. Add broader advanced-option coverage on the tool pages.
2. Improve Run History with better recent-run navigation and filtering.
3. Add cleaner success/failure summaries after background launch.

### Optional Enhancements
1. Inline PDF/HTML previews.
2. Cancel running jobs.
3. Re-run from previous run folder.
4. Preset templates for common workflows.

## Keep This File
This file is still needed.

It is no longer just a future design draft. It now serves as:
- implementation status
- remaining-work tracker
- boundary document for keeping the GUI isolated from core `fp-tools` logic
