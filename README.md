# fp-tools

`fp-tools` is a Python package for ATAC-seq footprint analysis. It helps users turn ATAC-seq data into bias-corrected cut-site tracks, footprint scores, motif-centered comparisons, aggregate plots, and static HTML reports.

Documentation: <https://oncologylab.github.io/fp-tools/>

## Install

```bash
pip install fp-tools-bio
```

To use the browser interface:

```bash
pip install "fp-tools-bio[gui]"
fp-tools-gui
```

The GUI opens in a browser and writes the same YAML configs that can be run from the command line.

## What You Can Do

- Correct ATAC-seq cut-site signal for Tn5 sequence bias.
- Call footprint scores from corrected bigWig tracks.
- Scan known motif databases, including bundled JASPAR 2026 and HOCOMOCO files.
- Compare footprint scores across conditions or replicates.
- Generate volcano-style differential footprint HTML reports.
- Plot motif-centered aggregate footprints.
- Group single-cell ATAC fragments into pseudobulk cell-type profiles.

## Typical Workflow

```text
atac-correct
  -> call-footprints
  -> match-motifs
  -> diff-footprints
  -> plot-aggregate or plot-aggregate-batch
```

Use `match-motifs` to inspect motif sites and bound/unbound calls for a single footprint track. For two or more conditions, use `diff-footprints`; it can scan the same motif database, compare conditions, and write an interactive HTML report.

## Minimal Example

```bash
atac-correct \
  --bam sample.bam \
  --genome hg38.fa.gz \
  --peaks peaks.bed \
  --blacklist hg38.blacklist.bed \
  --outdir results/atac_correct/sample

call-footprints \
  --signal results/atac_correct/sample/sample_corrected.bw \
  --regions peaks.bed \
  --output results/footprints/sample_footprints.bw

match-motifs \
  --signals results/footprints/sample_footprints.bw \
  --genome hg38.fa.gz \
  --peaks peaks.bed \
  --motif-db jaspar2026_vertebrates \
  --outdir results/motif_matches/sample

diff-footprints \
  --signals \
    conditionA_rep1_footprints.bw conditionA_rep2_footprints.bw \
    conditionB_rep1_footprints.bw conditionB_rep2_footprints.bw \
  --aggregate-signals \
    conditionA_rep1_corrected.bw conditionA_rep2_corrected.bw \
    conditionB_rep1_corrected.bw conditionB_rep2_corrected.bw \
  --genome hg38.fa.gz \
  --peaks peaks.bed \
  --cond-names conditionA conditionA conditionB conditionB \
  --motif-db jaspar2026_vertebrates \
  --normalization none \
  --plot-aggregate sig \
  --outdir results/diff_footprints/conditionA_vs_conditionB
```

Repeated names in `--cond-names` define biological replicates. The main report will be written inside the output folder as a standalone HTML file.

## Pseudobulk Workflow

Use `pseudobulk-footprints` when starting from single-cell ATAC fragments and cell annotations. With a motif database, marker motif-site BEDs, and an h5ad file containing the cell embedding, the standard pseudobulk route also writes Fig4-style single-cell footprinting plots.

```bash
pseudobulk-footprints \
  --fragments pbmc_fragments.tsv.gz \
  --annotations cell_annotations.tsv \
  --group-by cell_type \
  --genome-sizes hg38.chrom.sizes \
  --genome hg38.fa.gz \
  --peaks peaks.bed \
  --motif-db jaspar2026_vertebrates \
  --tf-site-dir marker_motif_sites \
  --single-cell-signature-h5ad pbmc_embedding.h5ad \
  --outdir results/pseudobulk
```

Key outputs include pseudobulk fragments, pseudo-BAMs, corrected bigWigs, footprint-score bigWigs, differential footprint reports, aggregate plots, and `plots/single_cell_footprinting/single_cell_footprinting.svg`.

## De Novo Motif Workflow

Use candidate footprints from `call-footprints --output-bed` to prepare de novo motif discovery. This route writes a reproducible MEME/STREME/DREME script and can compare discovered motifs to a built-in motif database.

```bash
motif-discovery \
  --candidates results/footprints/sample_candidates.bed \
  --genome hg38.fa.gz \
  --flank 75 \
  --method streme \
  --known-motif-db jaspar2026_vertebrates \
  --outdir results/de_novo/sample
```

Use discovered motifs alone for a de novo-only run, or add them to a standard database run:

```bash
diff-footprints \
  --signals conditionA_footprints.bw conditionB_footprints.bw \
  --genome hg38.fa.gz \
  --peaks peaks.bed \
  --cond-names conditionA conditionB \
  --motif-db jaspar2026_vertebrates \
  --motifs results/de_novo/sample/streme/streme.txt \
  --outdir results/diff_footprints/database_plus_de_novo
```

## Main Commands

| Command | Use |
| --- | --- |
| `atac-correct` | Bias-correct ATAC-seq cut-site signal. |
| `call-footprints` | Create footprint score tracks. |
| `match-motifs` | Scan motifs for one sample. |
| `diff-footprints` | Compare motif footprints across conditions. |
| `normalize-bigwig` | Scale bigWigs before aggregate plotting. |
| `plot-aggregate` | Make static aggregate footprint plots. |
| `plot-aggregate-batch` | Make interactive aggregate HTML reports. |
| `motif-discovery` | Prepare candidate-centered de novo motif discovery. |
| `motif-summary` | Summarize motif discovery outputs. |
| `fp-tools-score-variants` | Score variants against candidate footprints and motifs. |
| `pseudobulk-fragments` | Group single-cell fragments by annotation. |
| `pseudobulk-footprints` | Run the full pseudobulk footprint workflow, including optional Fig4-style single-cell plots. |
| `run-workflow` | Run a saved YAML config. |
| `fp-tools-gui` | Open the optional browser GUI. |

Check any command with `--help`:

```bash
atac-correct --help
call-footprints --help
match-motifs --help
diff-footprints --help
normalize-bigwig --help
plot-aggregate --help
plot-aggregate-batch --help
run-workflow --help
fp-tools-gui --help
motif-discovery --help
motif-summary --help
fp-tools-score-variants --help
pseudobulk-fragments --help
pseudobulk-footprints --help
```

## GUI

```bash
fp-tools-gui --host 0.0.0.0 --port 8891
```

Open the printed URL in a browser. If running on a server or cloud VM, make sure the selected port is allowed by the firewall or security group.

## More Examples

Example YAML configs are in `examples/gui_configs/`. They can be loaded in the GUI or run directly:

```bash
run-workflow --config examples/gui_configs/footprintscores_single.yml
```

Static report demos and GUI screenshots are available in the documentation site.
