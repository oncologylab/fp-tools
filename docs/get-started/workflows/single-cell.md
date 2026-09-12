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

| Input | Required content |
| --- | --- |
| `fragments.tsv.gz` | Fragment records with chromosome, start, end, and barcode in the first four columns. |
| `cell_annotations.tsv` | Cell barcodes, group labels, and UMAP coordinates, as shown below. |
| `cell_embedding.h5ad` | AnnData with genomic-bin counts and matching cell names; an embedding alone is insufficient. |
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
  --h5ad cell_embedding.h5ad \
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
