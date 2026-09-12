---
core_nav:
  previous:
    title: Tool overview
    url: get-started/tool-overview/
  next:
    title: call-footprints
    url: get-started/commands/call-footprints/
---

# [`atac-correct`](../../api.md#atac-correct)

Correct ATAC-seq cut-site signal for Tn5 sequence bias. Start with
coordinate-sorted BAM files, adjacent BAI indexes, and matching peak BED files.
Use the corrected signal as the input to `call-footprints`.

## Example command

```bash
atac-correct \
  --sample-table project/metadata/samples.tsv \
  --genome hg38.fa.gz \
  --blacklist hg38.blacklist.bed \
  --outdir project
```

## Primary inputs

- `--sample-table` — tab-separated file with one row per sample and the columns `sample`, `condition`, `bam`, and `peaks`. Give each sample a unique name and use the same condition label for biological replicates.
- `--genome` — reference FASTA whose chromosome names and assembly match every BAM and peak BED.
- `--blacklist` — optional assembly-matched BED of regions to exclude from bias estimation and corrected output.
- `--outdir` — project directory represented by `{project}` below.

Use the [sample-table example](../workflows/bulk-atac-seq.md#1-prepare-the-sample-table)
to prepare your inputs. Replace the FASTA and blacklist paths with your own
files. A custom FASTA does not automatically select a blacklist.

## Main outputs

For each `{sample}` in the table, this command writes:

| Path | Meaning |
| --- | --- |
| `{project}/samples/{sample}/atac_correct/{sample}_corrected.bw` | Bias-corrected cut-site signal. Positive positions have more observed cuts than expected; negative positions have fewer. |
| `{project}/samples/{sample}/atac_correct/{sample}_atacorrect.pdf` | Diagnostic plots comparing learned Tn5 sequence bias before and after correction. Omitted with `--skip-qc`. |
| `{project}/samples/{sample}/atac_correct/{sample}_AtacBias.pickle` | Saved bias model for reuse. Downstream commands use the corrected bigWig and do not need this file. |

With `--write-tracks all`, the same directory also contains:

| Path | Meaning |
| --- | --- |
| `{sample}_uncorrected.bw` | Observed base-resolution cut-site signal after the configured forward/reverse read shifts and sequencing-depth normalization. |
| `{sample}_bias.bw` | Tn5 sequence-bias score predicted from the reference sequence. |
| `{sample}_expected.bw` | Expected cut-site signal after the sequence-bias score is scaled to local observed cuts. |

Project-level peak outputs are `{project}/peaks/merged_peaks.bed` and
`{project}/peaks/merged_peaks_filtered.bed`. A direct single-BAM run writes the
same `{prefix}_*.bw`, `{prefix}_atacorrect.pdf`, and
`{prefix}_AtacBias.pickle` patterns directly under `{outdir}`.

Open the QC PDF to inspect the sequence-bias diagnostics, then use the
corrected bigWig for footprint scoring. A negative corrected value means fewer
cuts than the bias model expected; it does not by itself identify a bound TF.

[Open a representative ENCODE ATACCorrect PDF](../../demos/qc/encode/A549_rep1_atacorrect.pdf).
Continue with [`call-footprints`](call-footprints.md), or see the
[complete `atac-correct` reference](../../api.md#atac-correct).
