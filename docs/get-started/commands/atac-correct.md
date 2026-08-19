---
core_nav:
  previous:
    title: prepare-atac
    url: get-started/commands/prepare-atac/
  next:
    title: call-footprints
    url: get-started/commands/call-footprints/
---

# [`atac-correct`](../../api.md#atac-correct)

Estimate Tn5 sequence bias from aligned ATAC-seq fragments and subtract the
expected bias contribution from the observed cut-site signal. Run this before
footprint scoring.

## Example command

```bash
atac-correct \
  --sample-table project/metadata/samples.tsv \
  --genome hg38.fa.gz \
  --blacklist hg38.blacklist.bed \
  --outdir project
```

## Primary inputs

- `--sample-table` — TSV with `sample`, `condition`, `bam`, and `peaks`; BAM indexes must be adjacent to the BAMs.
- `--genome` — reference FASTA whose chromosome names and assembly match every BAM and peak BED.
- `--blacklist` — BED intervals excluded from bias estimation and corrected output.
- `--outdir` — project directory represented by `{project}` below.

## Main outputs

For each `{sample}`, project layout writes:

| Path | Meaning |
| --- | --- |
| `{project}/samples/{sample}/atac_correct/{sample}_corrected.bw` | Bias-corrected cut-site signal. Positive positions have more observed cuts than expected; negative positions have fewer. |
| `{project}/samples/{sample}/atac_correct/{sample}_atacorrect.pdf` | Diagnostic plots comparing learned Tn5 sequence bias before and after correction. Omitted with `--skip-qc`. |
| `{project}/samples/{sample}/atac_correct/{sample}_AtacBias.pickle` | Serialized learned bias model for reuse or advanced debugging. It is not required by downstream commands. |

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

[Open a representative ENCODE ATACCorrect PDF](../../demos/qc/encode/A549_rep1_atacorrect.pdf).
Continue with [`call-footprints`](call-footprints.md), or see the
[complete `atac-correct` reference](../../api.md#atac-correct).
