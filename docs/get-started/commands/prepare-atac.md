# [`prepare-atac`](../../api.md#prepare-atac)

**Linux CLI or Linux container only.** Download public ATAC-seq reads or use
local FASTQ files, then trim, align, filter, call peaks, calculate alignment
coverage, and write QC files. The GUI and native macOS/Windows installations
start from filtered BAM/BAI and peak BED files.

## Example command

```bash
prepare-atac --samples metadata.tsv --genome hg38 --outdir project
```

## Primary inputs

- `--samples` — TSV or CSV sample sheet containing `sample`, `condition`, and either paired `fastq_1`/`fastq_2` paths or URLs. See the [ENCODE and local examples](../workflows/bulk-atac-seq.md#optional-fastq-preparation-on-linux).
- `--genome` — packaged `hg38` or `mm10` reference label, or a custom label used with explicit reference options.
- `--outdir` — project directory represented by `{project}` below.

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
| `{project}/samples/{sample}/tracks/{sample}.rp10m.bw` | Sequencing-depth-normalized alignment coverage bigWig. `rp10m` is retained only as the historical filename suffix. |
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
| `{project}/metadata/samples.tsv` | Downstream `sample`, `condition`, `bam`, and `peaks` table accepted by core commands. |
| `{project}/reports/qc_summary.tsv` | Cross-sample QC summary. |

Continue with [`atac-correct`](atac-correct.md), or see the
[complete `prepare-atac` reference](../../api.md#prepare-atac).
