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
- Plot per-cell footprint-signature heatmaps and UMAP reports from single-cell fragments.

## Typical Workflow

[`atac-correct`](https://oncologylab.github.io/fp-tools/api/#atac-correct)
-> [`call-footprints`](https://oncologylab.github.io/fp-tools/api/#call-footprints)
-> [`match-motifs`](https://oncologylab.github.io/fp-tools/api/#match-motifs) or [`motif-discovery`](https://oncologylab.github.io/fp-tools/api/#motif-discovery)
-> [`diff-footprints`](https://oncologylab.github.io/fp-tools/api/#diff-footprints)
-> [`plot-aggregate`](https://oncologylab.github.io/fp-tools/api/#plot-aggregate)

Use [`match-motifs`](https://oncologylab.github.io/fp-tools/api/#match-motifs) to inspect motif sites and bound/unbound calls for one or more footprint tracks. Use [`motif-discovery`](https://oncologylab.github.io/fp-tools/api/#motif-discovery) when candidate footprint intervals should be searched for de novo motifs. For two or more conditions, use [`diff-footprints`](https://oncologylab.github.io/fp-tools/api/#diff-footprints); it can scan the same motif database, compare conditions, and write an interactive HTML report.

## Minimal Example

Create a simple sample table:

```text
sample	condition	bam	peaks
A	conditionA	A.bam	A_peaks.bed
B	conditionB	B.bam	B_peaks.bed
```

For condition comparisons, create a comparison table:

```text
comparison	cond1	cond2
conditionA_vs_conditionB	conditionA	conditionB
```

```bash
atac-correct \
  --sample-table project/metadata/samples.tsv \
  --genome hg38.fa.gz \
  --blacklist hg38.blacklist.bed \
  --outdir project

normalize-bigwig \
  --sample-table project/metadata/samples.tsv \
  --background project/peaks/merged_peaks_filtered.bed \
  --outdir project \
  --method background-scale \
  --stat q95 \
  --target median

call-footprints \
  --sample-table project/metadata/samples.tsv \
  --regions project/peaks/merged_peaks_filtered.bed \
  --outdir project

match-motifs \
  --sample-table project/metadata/samples.tsv \
  --genome hg38.fa.gz \
  --peaks project/peaks/merged_peaks_filtered.bed \
  --motif-db jaspar2026_vertebrates \
  --outdir project

diff-footprints \
  --sample-table project/metadata/samples.tsv \
  --comparison-table project/metadata/comparisons.tsv \
  --genome hg38.fa.gz \
  --peaks project/peaks/merged_peaks_filtered.bed \
  --motif-db jaspar2026_vertebrates \
  --outdir project
```

With a sample table and `--outdir project`, fp-tools uses the recommended project layout by default: merged peaks in `project/peaks`, per-sample outputs in `project/samples/<sample>/`, differential reports in `project/comparisons`, and review pages in `project/reports`. `atac-correct` writes both the raw merged peak set and a filtered `project/peaks/merged_peaks_filtered.bed` with mitochondrial chromosomes excluded; downstream project commands use this filtered BED when peak/background regions are omitted or when the project raw merged BED is passed. Custom paths remain available with `--layout custom`.

After `match-motifs`, project-mode `diff-footprints` reuses each sample folder's motif-site and background-score caches instead of rescanning motifs or rereading footprint bigWigs. The first cached comparison may create internal per-motif shard caches for faster reuse; later comparisons with the same sample folders reuse those shards. In project mode, `match-motifs` uses one shared motif scan across samples by default and then writes standard per-sample folders. Repeated samples with the same `condition` are treated as biological replicates. Per-motif BED folders are written by default in the background after report-ready outputs; use `match-motifs --motif-outputs summary` only when you want cache-only output. HTML aggregate profiles in `sig` and `top` modes are capped by `--plot-aggregate-top-n`; increase this value to show more motif profiles.

For a portable multi-condition project script, see
`examples/nutrient_stress_project/run_ctrl_vs_10fbs.sh`. It documents the
expected raw-data layout, explains how to prepare clean `samples.tsv` and
`comparisons.tsv` files from an `ATAC_Nutrients_hg38_*.txt` table, supports
`CHECK_ONLY=1` input validation, and then runs the full project workflow.

## Methodological Improvements

fp-tools keeps the interpretable TOBIAS-style center-versus-flank footprint score, but improves the workflow around it. Multi-sample projects can q95-scale corrected cut-site tracks over shared background regions before footprint scoring, reducing sample-level signal shifts without forcing full distributions to match. Footprint scoring and candidate detection use optimized Cython-backed kernels by default, with a legacy kernel available for exact historical comparisons. Known-motif analysis uses one shared motif scan across project samples, writes compact motif-site/background caches, and lets `diff-footprints` reuse those caches for replicate-aware comparisons and interactive reports instead of rescanning motifs for every contrast.

## Pseudobulk Workflow

For single-cell ATAC data, start with [`pseudobulk-fragments`](https://oncologylab.github.io/fp-tools/api/#pseudobulk-fragments) to group fragments by cell annotation. Then run the standard footprint workflow on the grouped pseudobulk samples. After motif-aware comparisons are available, use [`find-signature-fp`](https://oncologylab.github.io/fp-tools/api/#find-signature-fp) to write marker footprint-signature heatmaps and UMAP reports. [`pseudobulk-footprints`](https://oncologylab.github.io/fp-tools/api/#pseudobulk-footprints) is the convenience wrapper that can run grouping, correction, scoring, motif reports, aggregate plots, and optional signature reporting in one command.

```bash
pseudobulk-footprints \
  --fragments pbmc_fragments.tsv.gz \
  --annotations cell_annotations.tsv \
  --group-by cell_type \
  --genome-sizes hg38.chrom.sizes \
  --genome hg38.fa.gz \
  --peaks merged_peaks.bed \
  --motif-db jaspar2026_vertebrates \
  --tf-site-dir marker_motif_sites \
  --single-cell-signature-h5ad pbmc_embedding.h5ad \
  --outdir project/pseudobulk
```

To run the signature report as a standalone step:

```bash
find-signature-fp \
  --annotations cell_annotations.tsv \
  --fragments pbmc_fragments.tsv.gz \
  --h5ad pbmc_embedding.h5ad \
  --tf-site-dir marker_motif_sites \
  --all-motif-results project/pseudobulk/pseudobulk_diff_footprints_results.txt \
  --all-motif-diff-dir project/pseudobulk/diff_footprints \
  --outdir project/pseudobulk/signature_fp
```

Key outputs include pseudobulk fragments, pseudo-BAMs, corrected bigWigs, footprint-score bigWigs, differential footprint reports, aggregate plots, and single-cell footprint-signature heatmaps/UMAPs.

## De Novo Motif Workflow

Use candidate footprints from [`call-footprints`](https://oncologylab.github.io/fp-tools/api/#call-footprints) `--call-candidates` in project mode, or `--output-bed`/`--output-beds` in custom mode, to prepare de novo motif discovery. This route writes a reproducible MEME/STREME/DREME script and can compare discovered motifs to a built-in motif database.

```bash
motif-discovery \
  --candidates project/samples/sample/footprints/sample_candidate_footprints.bed \
  --genome hg38.fa.gz \
  --flank 75 \
  --method streme \
  --known-motif-db jaspar2026_vertebrates \
  --outdir project/de_novo/sample
```

Use discovered motifs alone for a de novo-only run, or add them to a standard database run:

```bash
diff-footprints \
  --signals conditionA_footprints.bw conditionB_footprints.bw \
  --sample-names conditionA conditionB \
  --genome hg38.fa.gz \
  --peaks merged_peaks.bed \
  --cond-names conditionA conditionB \
  --motif-db jaspar2026_vertebrates \
  --motifs project/de_novo/sample/streme/streme.txt \
  --outdir project/comparisons/database_plus_de_novo
```

## Main Commands

| Command | Use |
| --- | --- |
| [`atac-correct`](https://oncologylab.github.io/fp-tools/api/#atac-correct) | Bias-correct ATAC-seq cut-site signal. |
| [`call-footprints`](https://oncologylab.github.io/fp-tools/api/#call-footprints) | Create footprint score tracks from one or more bigWigs. |
| [`match-motifs`](https://oncologylab.github.io/fp-tools/api/#match-motifs) | Scan motifs for one or more footprint tracks. |
| [`diff-footprints`](https://oncologylab.github.io/fp-tools/api/#diff-footprints) | Compare motif footprints across conditions. |
| [`normalize-bigwig`](https://oncologylab.github.io/fp-tools/api/#normalize-bigwig) | Scale bigWigs before aggregate plotting. |
| [`plot-aggregate`](https://oncologylab.github.io/fp-tools/api/#plot-aggregate) | Make aggregate footprint plots as PDF/SVG-style output or HTML. |
| [`review-multi-comparisons`](https://oncologylab.github.io/fp-tools/api/#review-multi-comparisons) | Review multiple differential-footprint HTML reports in one page. |
| [`plot-motif-aggregate-grid`](https://oncologylab.github.io/fp-tools/api/#plot-motif-aggregate-grid) | Export multi-page motif-by-comparison aggregate PDFs from review reports. |
| [`motif-discovery`](https://oncologylab.github.io/fp-tools/api/#motif-discovery) | Prepare candidate-centered de novo motif discovery. |
| [`motif-summary`](https://oncologylab.github.io/fp-tools/api/#motif-summary) | Summarize motif discovery outputs. |
| [`pseudobulk-fragments`](https://oncologylab.github.io/fp-tools/api/#pseudobulk-fragments) | Group single-cell fragments by annotation. |
| [`find-signature-fp`](https://oncologylab.github.io/fp-tools/api/#find-signature-fp) | Plot per-cell footprint-signature heatmaps and UMAP reports. |
| [`pseudobulk-footprints`](https://oncologylab.github.io/fp-tools/api/#pseudobulk-footprints) | Run the full pseudobulk footprint workflow, including optional signature reporting. |
| [`run-workflow`](https://oncologylab.github.io/fp-tools/api/#run-workflow) | Run a saved YAML config. |
| [`fp-tools-gui`](https://oncologylab.github.io/fp-tools/api/#fp-tools-gui) | Open the optional browser GUI. |
| [`fp-tools-score-variants`](https://oncologylab.github.io/fp-tools/api/#fp-tools-score-variants) | Optional variant-annotation utility; not part of the standard footprint workflow. |

Check any command with `--help`:

```bash
atac-correct --help
call-footprints --help
match-motifs --help
diff-footprints --help
normalize-bigwig --help
plot-aggregate --help
review-multi-comparisons --help
plot-motif-aggregate-grid --help
run-workflow --help
fp-tools-gui --help
motif-discovery --help
motif-summary --help
fp-tools-score-variants --help
pseudobulk-fragments --help
find-signature-fp --help
pseudobulk-footprints --help
```

## GUI

```bash
fp-tools-gui --host 0.0.0.0 --port 8891
```

Open the printed URL in a browser. If running on a server or cloud VM, make sure the selected port is allowed by the firewall or security group.

## More Examples

Example YAML configs are in `examples/gui_configs/`. They can be loaded in the GUI or run directly with [`run-workflow`](https://oncologylab.github.io/fp-tools/api/#run-workflow):

```bash
run-workflow --config examples/gui_configs/call_footprints_single.yml
```

Static report demos and GUI screenshots are available in the documentation site.
