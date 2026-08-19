# Installation

fp-tools supports Python 3.11–3.13 on current Windows, macOS, and Linux systems.

## Command-line package

=== "Windows PowerShell"

    ```powershell
    py -m pip install --upgrade fp-tools-bio
    ```

=== "macOS / Linux"

    ```bash
    python3 -m pip install --upgrade fp-tools-bio
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

| Platform | Core analysis and GUI | Raw FASTQ preparation | External MEME tools |
| --- | --- | --- | --- |
| Windows 10/11 | Native | WSL | WSL |
| macOS Intel/Apple Silicon | Native | Install command-line tools | Install MEME Suite |
| Linux x86_64/ARM64 | Native | Install command-line tools | Install MEME Suite |

`prepare-atac` orchestrates tools such as Bowtie 2, SAMtools, BEDTools, fastp,
and MACS3. De novo discovery can call STREME and Tomtom. These external programs
are not required for BAM-, fragment-, bigWig-, or BED-based fp-tools analyses.

Next, choose the [bulk ATAC-seq workflow](workflows/bulk-atac-seq.md), the
[single-cell workflow](workflows/single-cell.md), or an individual command in
the [tool overview](tool-overview.md).
