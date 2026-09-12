# [`prepare-atac`](../../api.md#prepare-atac)

**Linux CLI or Linux container only.** Download public ATAC-seq reads or use
local FASTQ files, then trim, align, filter, call peaks, calculate alignment
coverage, and write QC files. The GUI and native macOS/Windows installations
start from coordinate-sorted BAM/BAI and matching peak BED files.

## Example command

```bash
prepare-atac --samples metadata.tsv --genome hg38 --outdir project
```

## Primary inputs

- `--samples` — TSV or CSV sample sheet. Provide `sample`, `condition`, and `fastq_1` paths or URLs; add `fastq_2` for paired-end reads. A public sequencing run can instead be supplied in `run_accession`.
- `--genome` — managed `hg38` or `mm10` reference label, or a custom label used with explicit reference options.
- `--outdir` — project directory represented by `{project}` below.

For paired-end local files, `metadata.tsv` can contain:

```tsv
sample	condition	fastq_1	fastq_2
control_1	control	reads/control_1_R1.fastq.gz	reads/control_1_R2.fastq.gz
treated_1	treated	reads/treated_1_R1.fastq.gz	reads/treated_1_R2.fastq.gz
```

Repeated rows with the same `sample`, `condition`, and `replicate` combine
technical sequencing runs. Different `sample` values sharing a `condition` are
biological replicates.

## Main outputs

For each `{sample}`, the default modern profile writes:

| Path | Meaning |
| --- | --- |
| `{project}/samples/{sample}/alignment/{sample}.filtered.bam` | Coordinate-sorted, filtered ATAC-seq alignment used downstream. |
| `{project}/samples/{sample}/alignment/{sample}.filtered.bam.bai` | Samtools index for the filtered BAM. |
| `{project}/samples/{sample}/peaks/{sample}.narrowPeak` | MACS3 narrow-peak calls before project-level merging. |
| `{project}/samples/{sample}/tracks/{sample}.rp10m.bw` | Sequencing-depth-normalized alignment coverage bigWig for viewing read coverage. |
| `{project}/samples/{sample}/qc/{sample}.fastp.html` | Interactive Fastp read-trimming QC report. |
| `{project}/samples/{sample}/qc/{sample}.fastp.json` | Machine-readable Fastp metrics. |
| `{project}/samples/{sample}/qc/flagstat.tsv` | Samtools alignment and filtering counts. |
| `{project}/samples/{sample}/qc/fragment_lengths.tsv` | Fragment-length distribution used to inspect ATAC-seq periodicity. |
| `{project}/samples/{sample}/qc/metrics.json` | Consolidated per-sample QC metrics. |
| `{project}/samples/{sample}/qc/commands.log` | External commands used for that sample. |

Project-level files include:

| Path | Meaning |
| --- | --- |
| `{project}/peaks/merged_peaks.bed` | Union of sample peak intervals. |
| `{project}/peaks/merged_peaks_filtered.bed` | Analysis peak set after excluded chromosomes are removed. |
| `{project}/metadata/resolved_runs.tsv` | Resolved local/downloaded FASTQ files and run grouping. |
| `{project}/metadata/samples.tsv` | `sample`, `condition`, `bam`, and `peaks` table ready for `bulk-footprinting`. |
| `{project}/reports/qc_summary.tsv` | Cross-sample QC summary. |

## Reference options

With `--genome hg38` or `--genome mm10`, the command downloads and verifies the
matching reference and blacklist as needed. Use `--reference-dir` to choose
where references are cached.

For a custom genome label, provide `--fasta` and `--macs-genome-size`. Supply
`--bowtie2-index` if you already have an index; otherwise the command builds one
from the FASTA. `--blacklist` replaces the managed hg38/mm10 blacklist, while
`--no-blacklist` disables it. Custom FASTA inputs use no blacklist unless you
provide one.

Review `reports/qc_summary.tsv` and the per-sample QC reports, then pass
`metadata/samples.tsv` to [`bulk-footprinting`](bulk-footprinting.md) with your
comparison table. The [bulk workflow guide](../workflows/bulk-atac-seq.md) shows
the required comparison columns. See the
[complete `prepare-atac` reference](../../api.md#prepare-atac).
