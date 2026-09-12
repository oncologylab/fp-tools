---
hide:
  - toc
---

# Installation

Choose the recommended option for your computer. Bulk workflows in the GUI
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

[Download for Windows](https://github.com/oncologylab/fp-tools/releases/download/v0.2.3/fp-tools-gui-windows-x64.exe){ .md-button }
[Download for Apple silicon](https://github.com/oncologylab/fp-tools/releases/download/v0.2.3/fp-tools-gui-macos-apple-silicon.dmg){ .md-button }

Windows may ask you to confirm the unsigned app download. The macOS app
is unsigned and has not been notarized by Apple, so Gatekeeper may report that
Apple cannot verify the developer. Download it only from the official
OncologyLab GitHub release page and verify the published SHA-256 checksum.

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

Use Python 3.11–3.13:

```bash
python -m pip install --upgrade fp-tools-bio
bulk-footprinting --help
```

This installs all fp-tools commands, including the GUI. Use `<command> --help`
to view a command's options, or follow a workflow linked below.

On Windows, use `py` instead of `python` if needed.

To use the browser interface, run:

```bash
fp-tools-gui
```

The GUI normally opens in your browser. If it does not, open
`http://127.0.0.1:8891`.

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

    Docker provides the same complete command-line and GUI environment:

    ```bash
    docker build -t fp-tools:latest https://github.com/oncologylab/fp-tools.git#main
    docker run --rm -p 8891:8891 -v "${PWD}:/work" fp-tools:latest
    ```

    Open `http://127.0.0.1:8891`. Your current folder is available as `/work`.

    The Linux container also supports FASTQ-to-BAM preparation with
    `prepare-atac`; native Windows and macOS installations do not.

??? note "Optional FASTQ-to-BAM preparation"

    The footprinting workflows start from BAM/BAI and peak BED files. Linux
    users who need read preprocessing can run
    [`prepare-atac`](commands/prepare-atac.md) separately before starting the
    bulk workflow.

## Start an analysis

- [Bulk ATAC-seq workflow](workflows/bulk-atac-seq.md)
- [Single-cell ATAC-seq workflow](workflows/single-cell.md)
- [De novo motif discovery](workflows/de-novo-motif-discovery.md)
