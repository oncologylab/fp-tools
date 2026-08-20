---
hide:
  - toc
---

# Installation

Choose the recommended option for your computer. Every GUI and desktop route
starts from coordinate-sorted BAM/BAI files and matching peak BED files.

| Computer | Recommended installation |
| --- | --- |
| Windows 10/11 x64 | Desktop app |
| Mac with Apple silicon | Desktop app |
| Intel Mac or Linux | Python package |

## Desktop app

Download the app for your computer, then open it. The GUI starts in your
browser.

[Download for Windows](https://github.com/oncologylab/fp-tools/releases/download/v0.2.0rc10/fp-tools-gui-windows-x64.exe){ .md-button }
[Download for Apple silicon](https://github.com/oncologylab/fp-tools/releases/download/v0.2.0rc10/fp-tools-gui-macos-apple-silicon.dmg){ .md-button }

The preview apps are unsigned. Windows may ask you to confirm the download. On
macOS, Control-click the app and choose **Open** the first time.

Optional de novo motif discovery prepares its external tools on first use.

## Python package

Use Python 3.11–3.13:

```bash
python -m pip install --upgrade --pre fp-tools-bio
fp-tools-gui
```

The GUI normally opens automatically. If it does not, open
`http://127.0.0.1:8891`.

On Windows, use `py` instead of `python` if needed.

### Running on a remote Linux server

```bash
fp-tools-gui --host 0.0.0.0 --port 8891 --no-browser
```

Open `http://SERVER_IP:8891` from your computer, replacing `SERVER_IP` with the
server address. The port must be permitted by the server firewall.

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
