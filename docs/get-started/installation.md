# Installation

fp-tools requires Python 3.12 or later.

## Command-line package

```bash
pip install fp-tools-bio
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

## Optional GUI

```bash
pip install "fp-tools-bio[gui]"
fp-tools-gui
```

The GUI saves YAML configurations that can be run directly with
[`run-yaml-workflow`](commands/run-yaml-workflow.md).

## Source installation

```bash
git clone https://github.com/oncologylab/fp-tools.git
cd fp-tools
python -m pip install -e .
```

See the [API Reference](../api.md) for complete command options.

Next, choose the [bulk ATAC-seq workflow](workflows/bulk-atac-seq.md), the
[single-cell workflow](workflows/single-cell.md), or an individual command in
the [tool overview](tool-overview.md).
