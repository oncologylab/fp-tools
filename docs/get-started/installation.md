---
hide:
  - toc
---

# Installation

Choose the recommended option for your computer. Bulk workflows in the graphical
user interface (GUI)
and desktop app start from coordinate-sorted BAM/BAI files and matching peak
BED files. Single-cell workflows start from fragments and cell annotations.

| Computer | Recommended installation |
| --- | --- |
| Windows 10/11 x64 | Desktop app |
| Mac with Apple silicon | Desktop app |
| Intel Mac or Linux | Python package |

## Desktop app

Download the app for your computer, then open it. fp-tools opens in its own
application window; no browser or Python installation is required.

[Download for Windows](https://github.com/oncologylab/fp-tools/releases/download/v0.2.7/fp-tools-gui-windows-x64.exe){ .md-button }
[Download for Apple silicon](https://github.com/oncologylab/fp-tools/releases/download/v0.2.7/fp-tools-gui-macos-apple-silicon.dmg){ .md-button }

Windows may ask you to confirm the unsigned app download. The macOS app
is unsigned and has not been notarized by Apple, so Gatekeeper may report that
Apple cannot verify the developer. Download it only from the official
OncologyLab GitHub release page and verify the published SHA-256 checksum.

### Verify your download

Download [SHA256SUMS.txt](https://github.com/oncologylab/fp-tools/releases/download/v0.2.7/SHA256SUMS.txt)
from the same release as your app. In the folder containing your download, run:

=== "macOS Terminal"

    ```bash
    shasum -a 256 fp-tools-gui-macos-apple-silicon.dmg
    ```

=== "Windows PowerShell"

    ```powershell
    Get-FileHash .\fp-tools-gui-windows-x64.exe -Algorithm SHA256
    ```

Compare the complete hash with the line for that filename in `SHA256SUMS.txt`.
The values must match; uppercase and lowercase hexadecimal letters are equivalent.

On macOS, drag `fp-tools.app` to Applications and try to open it once. If macOS
blocks it, open **System Settings > Privacy & Security**, find the fp-tools
message, and select **Open Anyway**. On a managed Mac, an administrator may need
to approve the app.

As an advanced fallback, remove the quarantine attribute in Terminal and open
the app:

```bash
xattr -dr com.apple.quarantine /Applications/fp-tools.app && open /Applications/fp-tools.app
```

Use this command only after downloading fp-tools from the official OncologyLab
GitHub release page and verifying its checksum.

Optional de novo motif discovery prepares its external tools on first use.

## Python package

Use Python 3.11–3.13. Create an isolated environment in your analysis folder:

=== "macOS / Linux"

    ```bash
    python3 -m venv .venv
    source .venv/bin/activate
    python -m pip install --upgrade fp-tools-bio
    bulk-footprinting --help
    ```

=== "Windows PowerShell"

    ```powershell
    py -m venv .venv
    .\.venv\Scripts\python.exe -m pip install --upgrade fp-tools-bio
    .\.venv\Scripts\bulk-footprinting.exe --help
    ```

This installs all fp-tools commands, including the GUI. Use `<command> --help`
to view a command's options, or follow a workflow linked below.

If a command is not found, reactivate the environment in that terminal, or run
its executable directly from `.venv/bin/` (macOS/Linux) or `.venv\Scripts\`
(Windows). The desktop app does not require this Python setup.

Multi-line commands ending in `\` in this manual use bash or zsh, the usual
macOS/Linux shells. In Windows PowerShell, put the command on one line, removing
each trailing `\`, or use a backtick for continuation. File paths containing
spaces must be quoted. For example, this template works on one line:

```text
bulk-footprinting --sample-table samples.tsv --comparison-table comparisons.tsv --genome hg38 --outdir project
```

Replace these paths with your files and run from their containing folder. In an
unactivated Windows environment, use `.\.venv\Scripts\bulk-footprinting.exe`.

To use the browser interface, run:

```bash
fp-tools-gui
```

The GUI normally opens in your browser. If it does not, open the URL printed in
the terminal. The default port is 8891; the launcher can choose another available
port. In an unactivated Windows environment, run `.\.venv\Scripts\fp-tools-gui.exe`.

### Running on a remote Linux server

On the server, start the GUI:

```bash
fp-tools-gui --port 8891 --no-browser
```

On your computer, open an SSH tunnel, replacing `USER` and `SERVER` with your
login name and server address:

```bash
ssh -N -L 8891:127.0.0.1:8891 USER@SERVER
```

Keep both commands running and open `http://127.0.0.1:8891` on your computer.

??? note "Optional Docker installation"

    Docker provides a versioned command-line and graphical environment. With a
    current Docker installation, build the released source and start the local
    interface from the folder containing your data:

    ```bash
    docker build -t fp-tools:0.2.8 https://github.com/oncologylab/fp-tools.git#v0.2.8
    docker run --rm -p 127.0.0.1:8891:8891 -v "${PWD}:/work" fp-tools:0.2.8
    ```

    Open `http://127.0.0.1:8891`. Your current folder is mounted as `/work`, so
    select input files under `/work` in the app. This starts a local interface,
    not an authenticated multi-user service. For remote use, keep the loopback
    binding and use the SSH tunnel above.

    The release tag fixes the fp-tools source version. Dependency resolution can
    still change a rebuild; to reuse the packaged environment, download the
    matching architecture's container archive and checksum from the release,
    then load it with `docker load -i <archive.tar.gz>`. Published images use
    the tag `fp-tools:v0.2.8`.

    The Linux container also supports FASTQ-to-BAM preparation with
    `prepare-atac`; native Windows and macOS installations do not.

??? note "Optional FASTQ-to-BAM preparation"

    The footprinting workflows start from BAM/BAI and peak BED files. Linux
    users who need read preprocessing can run
    [`prepare-atac`](commands/prepare-atac.md) separately before starting the
    bulk workflow.

## Start an analysis

1. Choose an installation above.
2. Obtain matched inputs: your BAM/BAI and peaks, or the complete small example
   in the [single-cell workflow](workflows/single-cell.md).
3. Choose a new output folder and run the workflow or load its YAML settings file
   in the app.
4. Follow progress in the terminal or the app's run view; inspect logs if a stage
   fails.
5. Open the result named in the workflow guide and read its interpretation notes.

- [Bulk ATAC-seq workflow](workflows/bulk-atac-seq.md)
- [Single-cell ATAC-seq workflow](workflows/single-cell.md)
- [De novo motif discovery](workflows/de-novo-motif-discovery.md)
