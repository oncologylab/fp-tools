# Bulk ATAC-seq workflow

This guide starts with either paired FASTQ files or aligned BAM and peak files,
then compares motif-associated footprint scores between biological conditions.
The worked example uses three HepG2 and three K562 biological replicates from
ENCODE experiments [ENCSR291GJU](https://www.encodeproject.org/experiments/ENCSR291GJU/)
and [ENCSR868FGK](https://www.encodeproject.org/experiments/ENCSR868FGK/).

!!! warning "Plan storage and runtime before downloading"
    The six released BAM files total about 52 GB; their indexes and analysis
    outputs require additional space. The FASTQ example is also a substantial
    download. Start with `--dry-run` when using aligned data and use a compute
    system appropriate for whole-genome ATAC-seq.

## Before you start

Install fp-tools and verify the external programs used by raw-read processing:

```bash
prepare-atac --doctor --profile modern
```

For a custom genome, provide a matching FASTA, Bowtie2 index, chromosome names,
blacklist, and peak files. All BAMs, peaks, and the FASTA must use the same
assembly. The examples below use hg38.

## Starting from FASTQ files

[`prepare-atac`](../commands/prepare-atac.md) accepts a tab-separated or
comma-separated sample sheet. This is the six-sample ENCODE subset used here:

```text
run_accession	sample	condition	replicate	fastq_1	fastq_2	fastq_1_md5	fastq_2_md5
HepG2_bio1_run1	HepG2_rep3	HepG2	1	https://www.encodeproject.org/files/ENCFF529IMC/@@download/ENCFF529IMC.fastq.gz	https://www.encodeproject.org/files/ENCFF835SCS/@@download/ENCFF835SCS.fastq.gz	2d6a39e5386b5ed406373675f8f5b002	30e02a6230fcf312e1b346a8f904283a
HepG2_bio2_run1	HepG2_rep1	HepG2	2	https://www.encodeproject.org/files/ENCFF134LRP/@@download/ENCFF134LRP.fastq.gz	https://www.encodeproject.org/files/ENCFF734TQO/@@download/ENCFF734TQO.fastq.gz	4510b6f8383934206e1818f1d3dfa48a	a404d7362e67a008de8ac9c84ce327d6
HepG2_bio3_run1	HepG2_rep2	HepG2	3	https://www.encodeproject.org/files/ENCFF167KHV/@@download/ENCFF167KHV.fastq.gz	https://www.encodeproject.org/files/ENCFF317QCM/@@download/ENCFF317QCM.fastq.gz	156555346b1e676e8f2f1ba9f3e5f767	68fe7de2489647c61a6dbd59e0fee64a
K562_bio1_run1	K562_rep3	K562	1	https://www.encodeproject.org/files/ENCFF260KVA/@@download/ENCFF260KVA.fastq.gz	https://www.encodeproject.org/files/ENCFF761EDP/@@download/ENCFF761EDP.fastq.gz	0abafff85fc6023f0e7f1f56df17f645	d35a86e928b28d3884edcb38ef84cc8f
K562_bio2_run1	K562_rep2	K562	2	https://www.encodeproject.org/files/ENCFF098UCE/@@download/ENCFF098UCE.fastq.gz	https://www.encodeproject.org/files/ENCFF703BGR/@@download/ENCFF703BGR.fastq.gz	c0835be0503b8044de379f3018af30e8	ea88127b6dca5943fc711d6bfefbdea3
K562_bio3_run1	K562_rep1	K562	3	https://www.encodeproject.org/files/ENCFF354EXH/@@download/ENCFF354EXH.fastq.gz	https://www.encodeproject.org/files/ENCFF575ZTZ/@@download/ENCFF575ZTZ.fastq.gz	211c3a7e4244f06de4f22734a97f6371	a77fcb1cc2310ef81d2663de1e605af6
```

[Download the ENCODE FASTQ sheet](../../demos/data/encode/encode_hepg2_k562_fastq_urls.tsv)
or use the [local FASTQ template](../../demos/data/encode/local_fastq_template.tsv).
The downloadable ENCODE sheet selects one released paired sequencing run per
biological replicate, making the subset explicit and reproducible.

Each row is one paired sequencing run. Rows sharing `sample`, `condition`, and
`replicate` are treated as technical runs and combined into one biological
sample. `fastq_1_md5` and `fastq_2_md5` are optional for local files but strongly
recommended for downloaded files. Relative paths are interpreted from the
directory where the command is run.

```bash
prepare-atac --samples encode_hepg2_k562_fastq_urls.tsv --genome hg38 --outdir project --cores 8
```

The downstream aligned-data table is written to
`project/metadata/samples.tsv`; use that generated table for the remaining
commands.

## Starting from aligned ENCODE data

The sample table requires `sample`, `condition`, `bam`, and `peaks`. Additional
provenance columns are allowed and ignored by the analysis reader.

```text
sample	condition	bam	peaks
HepG2_rep1	HepG2	encode_data/bams/ENCFF624SON.bam	encode_data/peaks/ENCFF536RJV.bed
HepG2_rep2	HepG2	encode_data/bams/ENCFF926KFU.bam	encode_data/peaks/ENCFF536RJV.bed
HepG2_rep3	HepG2	encode_data/bams/ENCFF990VCP.bam	encode_data/peaks/ENCFF536RJV.bed
K562_rep1	K562	encode_data/bams/ENCFF077FBI.bam	encode_data/peaks/ENCFF855PCP.bed
K562_rep2	K562	encode_data/bams/ENCFF128WZG.bam	encode_data/peaks/ENCFF855PCP.bed
K562_rep3	K562	encode_data/bams/ENCFF534DCE.bam	encode_data/peaks/ENCFF855PCP.bed
```

Download the complete [HepG2/K562 sample table](../../demos/data/encode/encode_hepg2_k562_bams.tsv)
and [download/checksum helper](../../demos/data/encode/download_encode_hepg2_k562.sh),
then run the helper from the directory containing the sample table. It creates
the relative paths shown above and indexes each BAM with Samtools.
For your own already-aligned data, copy the
[local BAM/peak template](../../demos/data/encode/local_bam_peak_template.tsv)
and replace its portable relative paths.

Repeated `condition` values define biological replicates. Every BAM needs a
`.bai` index next to it. Peak files may be shared by samples from the same
condition; fp-tools merges all supplied peak sets for the project analysis.

## Define the comparison

The comparison direction is `cond1 - cond2`. For HepG2 versus K562:

```text
comparison	cond1	cond2
HepG2_vs_K562	HepG2	K562
```

[Download the HepG2/K562 comparison table](../../demos/data/encode/encode_hepg2_k562_comparisons.tsv).
Condition names must exactly match the sample table; the `comparison` value
becomes the output-directory name.

## Run the complete aligned-data workflow

```bash
bulk-footprinting --sample-table encode_hepg2_k562_bams.tsv --comparison-table encode_hepg2_k562_comparisons.tsv \
  --genome hg38.fa.gz --blacklist hg38.blacklist.bed --motif-db jaspar2026_vertebrates --outdir project --cores 8
```

Validate table contents and see the expanded commands before beginning the
long run:

```bash
bulk-footprinting --sample-table encode_hepg2_k562_bams.tsv --comparison-table encode_hepg2_k562_comparisons.tsv \
  --genome hg38.fa.gz --blacklist hg38.blacklist.bed --outdir project --dry-run
```

The wrapper runs `atac-correct`, `call-footprints`, `match-motifs`,
`diff-footprints`, and `review-multi-comparisons`. See the
[`bulk-footprinting` guide](../commands/bulk-footprinting.md) for its exact
output tree.

## Full seven-cell-line design

The complete locked study contains 17 biological samples across A549, HCT116,
HepG2, K562, MCF-7, PC-3, and Panc1, with 21 unordered pairwise comparisons:

- [17-sample ENCODE BAM/peak table](../../demos/data/encode/encode_cancer_7line_bams.tsv)
- [21-comparison table](../../demos/data/encode/encode_cancer_7line_comparisons.tsv)
- [Interactive ENCODE cancer-cell-line output](../../reports.md)

The downloadable full sample table records accessions, URLs, and MD5 checksums,
but its `bam` and `peaks` columns intentionally use portable local paths. It is
a study design and provenance file; download and index the listed files before
running it. ENCODE peak URLs end in `.bed.gz`: verify `peak_md5` against the
downloaded compressed file, then decompress it to the plain-text `.bed` path
shown in the `peaks` column.

!!! note "Why the wrapper output is not byte-for-byte identical to the demo"
    The interactive seven-cell-line demo uses the same locked ENCODE source
    accessions, but each comparison uses its own union of released IDR peaks and
    Q95 background scaling. `bulk-footprinting` uses one project-wide merged
    peak set and does not insert that pair-specific scaling step. Use the wrapper
    for the documented standard workflow; use the published demo as the exact
    record of the specialized pairwise analysis.

## Run the core commands separately

The same analysis can be run step by step:

<div class="fp-command-chain" markdown="1">

[`atac-correct`](../commands/atac-correct.md)
<span>→</span>
[`call-footprints`](../commands/call-footprints.md)
<span>→</span>
[`match-motifs`](../commands/match-motifs.md)
<span>→</span>
[`diff-footprints`](../commands/diff-footprints.md)

</div>

Use the exact project-layout paths documented on each command page. The
[Bulk output example](../output-examples/bulk-atac-seq.md) shows how to inspect
the final interactive report.

## Representative ENCODE QC files

These compact artifacts were produced from the locked A549 cancer-cell-line
inputs and illustrate the file formats emitted by the core commands. They are
representative QC examples, not a reproduction of the pair-specific demo:

- [ATACorrect diagnostic PDF](../../demos/qc/encode/A549_rep1_atacorrect.pdf)
- [Normalization QC table](../../demos/qc/encode/A549_normalize_bigwig_qc.tsv)
- [Replicate motif-score matrix](../../demos/qc/encode/A549_motif_score_matrix.tsv)
