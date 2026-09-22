# Single-cell ATAC-seq workflow

`sc-footprinting` groups fragments into pseudobulk samples, runs the footprint
analysis, and calculates per-cell footprint signatures. A pseudobulk sample
combines fragments from cells that share an annotation, such as a cell type.

## Main commands

<div class="fp-command-chain" markdown="1">

[`sc-footprinting`](../commands/sc-footprinting.md)

</div>

- [`sc-footprinting`](../commands/sc-footprinting.md) runs the complete workflow.
- [`pseudobulk-fragments`](../commands/pseudobulk-fragments.md) and [`find-signature-fp`](../commands/find-signature-fp.md) remain available as focused utilities.

## Prepare your inputs

For your first run, use the [complete PBMC example](#small-pbmc-example) below.
The general command and annotation rows in this section are templates for your
own data. Run them from a folder containing the named input files.

| Input | Required content |
| --- | --- |
| `fragments.tsv.gz` | Fragment records with chromosome, start, end, and barcode in the first four columns. |
| `cell_annotations.tsv` | Cell barcodes, group labels, and UMAP coordinates, as shown below. |
| `genomic_bin_counts.h5ad` | AnnData with genomic-bin counts and matching cell names; an embedding alone is insufficient. |
| `hg38.chrom.sizes` | Two columns: chromosome name and length. |
| `hg38.fa.gz` | Reference FASTA matching the fragments and peaks. |
| `merged_peaks.bed` | Accessible regions in the same assembly and chromosome naming scheme. |

The annotation table needs these columns:

```tsv
barcode	cell_type	snap_cell_type	umap_1	umap_2
AAACGAAAGAAACGCC-1	B_cell	B_cell	1.2	-0.8
AAACGAAAGAAAGCAG-1	T_cell	T_cell	-2.1	0.4
```

Use your actual barcodes and UMAP coordinates. `cell_type` contains broad
labels and `snap_cell_type` contains detailed labels; they can be identical
if you use one annotation level. `--group-by cell_type` combines cells with
the same `cell_type` value.

The AnnData cell names must match the annotation barcodes. Its features must
be genomic bins such as `chr1:0-500`, with a boolean `var['selected']` column
and a count matrix. The companion motif-activity calculation uses 500-base
bins. The workflow uses `obsm['X_spectral']` for nearest-neighbor smoothing
when available, otherwise it uses the annotation UMAP coordinates.

## Run the workflow

```bash
sc-footprinting \
  --fragments fragments.tsv.gz \
  --annotations cell_annotations.tsv \
  --h5ad genomic_bin_counts.h5ad \
  --group-by cell_type \
  --genome-sizes hg38.chrom.sizes \
  --genome hg38.fa.gz \
  --peaks merged_peaks.bed \
  --motif-db jaspar2026_vertebrates \
  --outdir project/pseudobulk
```

## Review the results

Check `project/pseudobulk/pseudobulk_footprint_manifest.tsv` for completed
groups. Per-cell score tables, heatmaps, and UMAPs are written under
`project/pseudobulk/plots/single_cell_footprinting/`. The
[command guide](../commands/sc-footprinting.md) lists the intermediate tracks
and motif comparison files.

Per-cell signatures use neighboring cells to reduce sparse signal. Compare
patterns across annotated groups without treating each score as an independent
binding measurement. See the
[Single-cell output example](../output-examples/single-cell-atac-seq.md) for a
visual guide.

## Small PBMC example

This real peripheral blood mononuclear cell (PBMC) example contains 300 cells:
100 B cells, 100 monocytes, and 100 T or natural killer cells. It retains human
hg38 chromosome 22, 101,637 genomic bins (5,669 selected), and 110,888 fragment
records. Counts and cell coordinates come from the annotated 10x PBMC5k
dataset distributed by SnapATAC2; coordinates were not recomputed for this
subset. Accessible regions are the source's selected bins, not newly called
peaks.

Download `fp-tools-pbmc-chr22-demo-v1.zip` and its `.sha256` file from the
[v0.2.8 release](https://github.com/oncologylab/fp-tools/releases/tag/v0.2.8).
The ZIP is approximately 14.3 MB. In macOS Terminal, run
`shasum -a 256 fp-tools-pbmc-chr22-demo-v1.zip`; on Linux use `sha256sum` instead
of `shasum -a 256`. In Windows PowerShell, run
`Get-FileHash .\fp-tools-pbmc-chr22-demo-v1.zip -Algorithm SHA256`.
Compare the complete hash with the `.sha256` file, then extract the entire
folder. It includes fragments and their index, annotations,
count AnnData, the chromosome FASTA and index, chromosome sizes, a matched
blacklist, selected regions, and `workflow.yml`. `manifest.json` records the
source hashes and subset definition; `SHA256SUMS.txt` records every input hash.

Open a terminal **inside the extracted `fp-tools-pbmc-chr22-demo-v1` folder**.
With the Python package installed, run this supplied example:

```bash
run-yaml-workflow --config workflow.yml
```

For the Windows or macOS desktop app, open the GUI's **Config** page and load
`workflow.yml`. Replace the relative input filenames in the editor with their
full paths inside the extracted folder, and set `outdir` to a new results
folder using its full path. Apply the configuration, resolve any validation
errors, then start the run. Saved YAML remains runnable from the command line.
The command uses all
available cores. The first use also prepares the motif database.

Check `results/pseudobulk_footprint_manifest.tsv` for the three completed groups.
Open the figures under `results/plots/single_cell_footprinting/` and inspect the
group-level signal tracks and motif comparisons alongside them. This small
subset teaches the workflow; it does not reproduce the full PBMC example
figures or validate binding, biological differences, or full-scale performance.
There are no biological replicates in this example.

## Prepare the full PBMC5k inputs

For a source installation, the repository includes
[`prepare_10x_pbmc5k_scatac.py`](https://github.com/oncologylab/fp-tools/blob/v0.2.8/benchmarks/scripts/prepare_10x_pbmc5k_scatac.py).
Run the following from the root of a v0.2.8 source checkout in a separate
Python environment with SnapATAC2, AnnData, pandas, and Matplotlib installed:

```bash
python benchmarks/scripts/prepare_10x_pbmc5k_scatac.py --chroms chr22
```

The script retrieves the matching fragment file and annotated AnnData through
`snapatac2.datasets.pbmc5k`, preserving the source count matrix and boolean
`var['selected']`. It writes barcode-matched annotations for all 4,437 cells
(422 B cells, 1,702 monocytes, and 2,313 T or natural killer cells), selected-bin
BED files, and a preparation summary. The `--chroms` setting limits the demo
region list and chromosome-size table; it does not subset the source cells or
count matrix. BED coordinates use zero-based starts and exclude the end.

The annotated AnnData comes from
[the SnapATAC2 PBMC5k archive](https://osf.io/download/e9vc3/); its SHA-256 is
`592f1551c27d0cfe4d81e7febad624d6b7d3ebf977b0c3ea64e06b3f3d76f078`.
The matching [10x fragment file](https://cf.10xgenomics.com/samples/cell-atac/2.0.0/atac_pbmc_5k_nextgem/atac_pbmc_5k_nextgem_fragments.tsv.gz)
has SHA-256 `5fe44c0f8f76ce1534c1ae418cf0707ca5ef712004eee77c3d98d2d4b35ceaec`.
Use an hg38 reference and matching blacklist when adapting the command above.
For the exact 300-cell subset, the versioned
[`build_pbmc_chr22_demo.py`](https://github.com/oncologylab/fp-tools/blob/v0.2.8/benchmarks/scripts/build_pbmc_chr22_demo.py)
accepts these prepared files plus the reference and blacklist; its `--help`
lists the required paths. It selects the first 100 sorted barcodes per broad
group and writes the complete ZIP and checksums.
