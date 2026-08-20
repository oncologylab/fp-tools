# Installation

Choose the desktop executable, complete container, or Python package.

## Desktop executable

The desktop download includes fp-tools, Python, the scientific libraries, and
the GUI. Python does not need to be installed.

=== "Windows"

    Download `fp-tools-gui-windows-x64.exe` from the
    [release page](https://github.com/oncologylab/fp-tools/releases),
    then open it. This preview is unsigned, so Windows may ask you to confirm.

=== "macOS Apple Silicon"

    Download `fp-tools-gui-macos-apple-silicon.dmg`, open it, and launch
    `fp-tools`. This preview is unsigned; on first launch, Control-click the app
    and choose **Open**.

The app downloads its pinned analysis runtime on the first FASTQ or motif-
discovery run. Windows uses a private WSL2 distribution; fp-tools enables it
with permission if needed, and Windows may request one restart.

## Complete container

The container includes fp-tools and its external genomics programs. Docker
automatically downloads the image when this one command starts the GUI.

=== "Windows PowerShell"

    ```powershell
    docker run --rm -p 8891:8891 -v "${PWD}:/work" ghcr.io/oncologylab/fp-tools-bio:latest
    ```

=== "macOS"

    ```bash
    docker run --rm -p 8891:8891 -v "$PWD:/work" ghcr.io/oncologylab/fp-tools-bio:latest
    ```

=== "Linux"

    ```bash
    docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp -p 8891:8891 -v "$PWD:/work" ghcr.io/oncologylab/fp-tools-bio:latest
    ```

Open `http://127.0.0.1:8891`. Files in the current directory appear under
`/work` in the GUI. Replace the default GUI command after the image name to run
a command directly, for example:

```bash
docker run --rm -v "$PWD:/work" ghcr.io/oncologylab/fp-tools-bio:latest atac-correct --help
```

## Command-line package

=== "Windows PowerShell"

    ```powershell
    py -m pip install --upgrade --pre fp-tools-bio
    ```

=== "macOS / Linux"

    ```bash
    python3 -m pip install --upgrade --pre fp-tools-bio
    ```

Confirm that the commands are available:

```bash
atac-correct --help
diff-footprints --help
```

## Managed analysis runtime

Raw-read and motif-discovery commands automatically install their pinned tools
in the fp-tools cache. Inspect or prepare that cache with:

```bash
fp-tools-runtime status
fp-tools-runtime install core
```

Use `--runtime system` to use programs already on `PATH`, or `--runtime
container` to run the complete public image. Workflows starting from BAM,
bigWig, BED, or fragment inputs do not download the raw-read runtime.

## GUI

```bash
fp-tools-gui
```

The standard installation includes the GUI. It opens a local browser and
prints the URL. See [`fp-tools-gui`](commands/fp-tools-gui.md) for remote-server
access.

The GUI saves YAML configurations that can be run directly with
[`run-yaml-workflow`](commands/run-yaml-workflow.md).

## Source installation

```bash
git clone https://github.com/oncologylab/fp-tools.git
cd fp-tools
python -m pip install -e .
```

See the [API Reference](../api.md) for complete command options.

## Platform support

| Platform | Desktop executable | Python package | Complete container |
| --- | --- | --- | --- |
| Windows 10/11 x64 | Available | Available | Docker Desktop |
| macOS Apple Silicon | Available | Available | Docker Desktop |
| macOS Intel | — | Available | Docker Desktop |
| Linux x86_64/ARM64 | — | Available | Docker Engine |

`prepare-atac` orchestrates tools such as Bowtie 2, SAMtools, BEDTools, fastp,
and MACS3. De novo discovery can call STREME and Tomtom. The complete container
includes these programs.

Next, choose the [bulk ATAC-seq workflow](workflows/bulk-atac-seq.md), the
[single-cell workflow](workflows/single-cell.md), or an individual command in
the [tool overview](tool-overview.md).
