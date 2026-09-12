# [`fp-tools-gui`](../../api.md#fp-tools-gui)

Launch the browser interface for configuring and running fp-tools commands.
The Windows and Apple Silicon desktop downloads present the same interface in
a native fp-tools application window.

Bulk GUI workflows start from coordinate-sorted BAM/BAI files and matching
peak BED files. The GUI does not perform FASTQ-to-BAM preprocessing.
Missing inputs and unsupported options are reported before a run starts.

The GUI is available through the Python package, the complete container, and
the self-contained desktop downloads on the
[release page](https://github.com/oncologylab/fp-tools/releases).

## Example command

```bash
fp-tools-gui --host 127.0.0.1 --port 8891 --run-dir project/gui_runs
```

## Primary inputs

- `--host` — interface on which the GUI listens (default: `127.0.0.1`).
- `--port` — fixed browser port.
- `--run-dir` — directory for GUI-managed configurations and runs.

## Start an analysis

1. Select the workflow or command from the sidebar.
2. Enter your input paths and output directory, then select **Update page config**.
3. Review the displayed command and resolve any validation errors before starting the run.
4. Open **Run History** to check progress, read logs, and find the output files.

Use the **Config** page to save the settings as YAML or load a previous run's
configuration.

## Main outputs

- `{run_dir}/{timestamp}_{label}/config.yml` — saved YAML settings for the run.
- `{run_dir}/{timestamp}_{label}/status.json`, `launcher_stdout.log`, and `launcher_stderr.log` — overall run status and logs.
- `{run_dir}/{timestamp}_{label}/{job_id}/status.json`, `command.txt`, `stdout.log`, and `stderr.log` — each job's status, exact command, and analysis logs.
- The analysis files documented by the selected command, written to the output directory you chose.

A saved YAML can also be run from the command line with
`run-yaml-workflow --config {run_dir}/{timestamp}_{label}/config.yml`.

Open the [GUI Demo](../../gui.md), or see the
[complete `fp-tools-gui` reference](../../api.md#fp-tools-gui).

## Local computer

Open the desktop executable to use the native application window. When using
the Python package, run `fp-tools-gui`; a browser opens after the server is
ready. If it does not, open the local URL printed in the terminal.

## Remote Linux server

Start fp-tools on the server. `--no-browser` prevents it from opening a browser
on the server, and the default host setting limits access to that server:

```bash
fp-tools-gui --no-browser --port 8891
```

On your computer, create an SSH tunnel and keep that terminal open:

```bash
ssh -N -L 8891:127.0.0.1:8891 USER@SERVER
```

Open `http://127.0.0.1:8891` on your computer to use the server's GUI.
