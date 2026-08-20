# Installation

Choose the desktop executable, complete container, or Python package.

## Desktop executable

The desktop download includes fp-tools, Python, the scientific libraries, and
the GUI. Python does not need to be installed.

=== "Windows"

    Download `fp-tools-gui-windows-x64.exe` from the
    [release page](https://github.com/oncologylab/fp-tools/releases),
    then open it.

=== "macOS Apple Silicon"

    Download `fp-tools-gui-macos-apple-silicon.tar.gz`, extract it, and open
    `fp-tools-gui`.

=== "macOS Intel"

    Download `fp-tools-gui-macos-intel.tar.gz`, extract it, and open
    `fp-tools-gui`.

=== "Linux x64 / ARM64"

    Download the matching `fp-tools-gui-linux-*.tar.gz`, extract it, and run:

    ```bash
    ./fp-tools-gui
    ```

The desktop executable covers the GUI and native fp-tools analyses from BAM,
fragment, bigWig, BED, and motif inputs. Use the complete container when a run
also needs raw-read programs or MEME Suite.

## Complete container

The container includes fp-tools and its external genomics programs. Docker
automatically downloads the image when this one command starts the GUI.

=== "Windows PowerShell"

    ```powershell
    docker run --rm -p 8891:8891 -v "${PWD}:/work" ghcr.io/oncologylab/fp-tools:latest
    ```

=== "macOS"

    ```bash
    docker run --rm -p 8891:8891 -v "$PWD:/work" ghcr.io/oncologylab/fp-tools:latest
    ```

=== "Linux"

    ```bash
    docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp -p 8891:8891 -v "$PWD:/work" ghcr.io/oncologylab/fp-tools:latest
    ```

Open `http://127.0.0.1:8891`. Files in the current directory appear under
`/work` in the GUI. Replace the default GUI command after the image name to run
a command directly, for example:

```bash
docker run --rm -v "$PWD:/work" ghcr.io/oncologylab/fp-tools:latest atac-correct --help
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

## External analysis programs

Raw-read preparation uses external genomics programs. The default `modern`
profile requires Fastp, Bowtie2 (including `bowtie2-build`), Samtools, Bedtools,
and MACS3. Check the active environment before downloading data:

```bash
prepare-atac --doctor --profile modern
```

The alternative `homer-atac` profile uses Trim Galore, Bowtie2, Picard,
Samtools, Bedtools, and HOMER. The doctor report identifies every missing
program. Commands that start from existing bigWigs do not require the complete
raw-read toolchain.

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

| Platform | Desktop / Python | Complete container |
| --- | --- | --- |
| Windows 10/11 x64 | Native | Docker Desktop |
| macOS Intel/Apple Silicon | Native | Docker Desktop |
| Linux x86_64/ARM64 | Native | Docker Engine |

`prepare-atac` orchestrates tools such as Bowtie 2, SAMtools, BEDTools, fastp,
and MACS3. De novo discovery can call STREME and Tomtom. The complete container
includes these programs.

Next, choose the [bulk ATAC-seq workflow](workflows/bulk-atac-seq.md), the
[single-cell workflow](workflows/single-cell.md), or an individual command in
the [tool overview](tool-overview.md).
