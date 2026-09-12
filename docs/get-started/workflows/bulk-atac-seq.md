# Bulk ATAC-seq workflow

`bulk-footprinting` runs a complete analysis from coordinate-sorted BAM files,
their adjacent BAI indexes, and peak BED files. The workflow corrects Tn5 bias,
scores footprints, matches motifs, compares conditions, and creates interactive
reports.

## 1. Prepare the sample table

Create `samples.tsv` with one row per biological sample and four required
columns. Here is a minimal two-condition, two-replicate design:

```tsv
sample	condition	bam	peaks
control_rep1	control	data/bams/control_rep1.sorted.bam	data/peaks/control_rep1.peaks.bed
control_rep2	control	data/bams/control_rep2.sorted.bam	data/peaks/control_rep2.peaks.bed
treated_rep1	treated	data/bams/treated_rep1.sorted.bam	data/peaks/treated_rep1.peaks.bed
treated_rep2	treated	data/bams/treated_rep2.sorted.bam	data/peaks/treated_rep2.peaks.bed
```

[Download this sample table](../../demos/data/bulk/samples.tsv). Each `sample`
must be unique. Samples with the same `condition` are treated as biological
replicates. Every BAM must be coordinate-sorted and have an adjacent index such
as `control_rep1.sorted.bam.bai`. BAMs and peak BEDs must use the same genome
assembly and chromosome names.

| Column | What to enter |
| --- | --- |
| `sample` | A unique name for one biological sample, such as `control_rep1`. |
| `condition` | The group to compare, such as `control` or `treated`. |
| `bam` | Path to that sample's coordinate-sorted BAM file. |
| `peaks` | Path to the matching peak BED file. |

## 2. Choose the comparison

Create `comparisons.tsv` to define which conditions to compare:

```tsv
comparison	cond1	cond2
treated_vs_control	treated	control
```

[Download this comparison table](../../demos/data/bulk/comparisons.tsv). The
reported change is `cond1` relative to `cond2`. Add another row for each
additional comparison. Use a unique `comparison` name for the output folder;
`cond1` and `cond2` must match values in the sample table's `condition` column.

## 3. Run the workflow

For human hg38 data, run:

```bash
bulk-footprinting \
  --sample-table samples.tsv \
  --comparison-table comparisons.tsv \
  --genome hg38 \
  --outdir project \
  --cores 8
```

On first use, fp-tools downloads the hg38 FASTA and blacklist, verifies their
checksums, builds the FASTA index, and saves them in the managed reference
cache. Later runs reuse the verified files. The default motif database is
`jaspar2026_vertebrates` when no motif option is supplied.

## Managed and custom references

fp-tools has managed support for these assemblies:

| `--genome` | Species | Automatically managed files |
| --- | --- | --- |
| `hg38` | Human | hg38 FASTA, FASTA index, and hg38 blacklist |
| `mm10` | Mouse | mm10 FASTA, FASTA index, and mm10 blacklist |

Managed files are stored under `~/.cache/fp-tools/references` by default. Use
`--reference-dir /shared/fp-tools-references` to place the cache elsewhere.
Use `--blacklist custom.blacklist.bed` to replace the managed blacklist, or
`--no-blacklist` to disable blacklist filtering.

For another assembly, provide the FASTA path explicitly. fp-tools never
guesses a blacklist for a custom reference:

```bash
bulk-footprinting \
  --sample-table samples.tsv \
  --comparison-table comparisons.tsv \
  --genome /references/custom.fa \
  --blacklist /references/custom.blacklist.bed \
  --outdir project \
  --cores 8
```

Omit `--blacklist` when the custom assembly has no blacklist.

## Choose motifs

Motif selection follows these rules:

| Options | Motifs used |
| --- | --- |
| No motif options | `jaspar2026_vertebrates` |
| `--motifs custom.jaspar` | Only motifs in `custom.jaspar` |
| `--motif-db DB --motifs custom.jaspar` | Motifs from `DB` plus the custom file |

Run `bulk-footprinting --list-motif-dbs` to list the packaged databases without
supplying workflow inputs.

## Review the results

The main outputs are the per-sample corrected cut-site and footprint-score
bigWigs, motif matches, differential statistics, aggregate profiles, and an
interactive comparison report under `{project}`. See the
[`bulk-footprinting` command guide](../commands/bulk-footprinting.md) for exact
file patterns and the [bulk output example](../output-examples/bulk-atac-seq.md)
for a visual tour.

Here, `{project}` means the directory supplied to `--outdir` (`project` in the
example). `{sample}` in an output filename means a name from your sample table.
Review the differential statistics together with the corrected cut-site
profiles and agreement between biological replicates.

## ENCODE example

To practice with public data, use the compact
[ENCODE BAM/peak sample sheet](../../demos/data/encode/encode_hepg2_k562_bams.tsv),
[comparison table](../../demos/data/encode/encode_hepg2_k562_comparisons.tsv),
and [download helper](../../demos/data/encode/download_encode_hepg2_k562.sh).
You can also explore the finished
[ENCODE cancer-cell-line reports](../../reports.md).

## Starting from FASTQ files

On the Linux CLI or in the Linux container, run
[`prepare-atac`](../commands/prepare-atac.md) first. It writes filtered BAM/BAI
files, peak BED files, and `metadata/samples.tsv`, which can then be passed to
`bulk-footprinting`. FASTQ preparation is not part of the BAM-first bulk
workflow and is not available in the native Windows or macOS applications.
