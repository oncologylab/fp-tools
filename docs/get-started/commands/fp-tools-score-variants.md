# `fp-tools-score-variants`

Annotate genomic variants with footprint overlap, local sequence changes, and
optional motif or model score changes.

## Example command

```bash
fp-tools-score-variants \
  --variants variants.vcf \
  --genome hg38.fa.gz \
  --out project/variants/variant_scores.tsv
```

## Primary inputs

- Variants in VCF or BED-like form.
- Reference genome FASTA.
- Optional footprint tracks, candidate intervals, motifs, or a trained model.

## Main outputs

- A variant-level TSV containing the requested overlap and score-change columns.

See the [complete `fp-tools-score-variants` reference](../../api.md#fp-tools-score-variants).
