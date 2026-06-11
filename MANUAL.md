# fp-tools Manual

`fp-tools` is a standalone footprinting package for ATAC-seq style workflows. It provides command-first tools for bias correction, footprint scoring, differential binding detection, and aggregate signal plotting.

The PyPI distribution is named `fp-tools-bio`; the installed Python package is `fp_tools`.

## Install

```bash
pip install fp-tools-bio
```

The core install is intentionally light for cluster and container use. The optional
Streamlit GUI is installed separately through the `gui` extra:

```bash
pip install "fp-tools-bio[gui]"
```

## Release status legend

The same three-way status model used in `README.md` applies here:

- **Released** — shipped in the current PyPI distribution (`fp-tools-bio`, version 0.1.7).
- **Development branch** — implemented and callable from the current source tree / `main`, with unit-test coverage, but biological validation is still expanding and the public release is pending.
- **Planned** — described in the validation roadmap (`DEV_PLAN.md`); not yet implemented or not yet benchmarked.

The PyPI page is authoritative for what `pip install fp-tools-bio` provides today.

## Commands

### Released (PyPI 0.1.7)

- `atac-correct`: correct ATAC-seq cutsite signal for Tn5 sequence bias.
- `score-footprints`: calculate footprint, multiscale, sum, mean, or pass-through scores from bigWig signal.
- `detect-tf-binding`: scan motifs, infer bound sites, and compare TF binding across conditions.
- `plot-aggregate`: plot aggregate signal around TFBS or region sets.
- `fp-tools-run`: run optional YAML batch configs.
- `fp-tools-gui`: launch the optional Streamlit GUI wrapper.

Legacy aliases remain available for compatibility: `ATACorrect`, `FootprintScores`, `ScoreBigwig`, `BINDetect`, and `PlotAggregate`.

### Development branch (current source tree)

Additional opt-in commands support tabular TFBS feature/model workflows, candidate generation and reranking, de novo motif-discovery orchestration, variant scoring, 10x-style pseudobulk fragment grouping with optional indexed fragments and cut-site bigWigs, replicate-aware BINDetect reports, and multiscale competition decomposition. These commands are present and callable from the current source tree and carry unit-test coverage; biological validation benchmarks for several of them are still being expanded, and they are not all guaranteed to be exposed in the released wheel until the next tagged release.

Direct CLI usage is the primary interface. YAML configs and the GUI are optional wrapper paths and do not replace the plain command-line tools.

## Verify

```bash
atac-correct --help
score-footprints --help
detect-tf-binding --help
plot-aggregate --help
fp-tools-run --help
fp-tools-gui --help
fp-tools-build-tfbs-features --help
fp-tools-train-tfbs-model --help
fp-tools-predict-tfbs --help
fp-tools-generate-candidates --help
fp-tools-rerank-candidates --help
fp-tools-export-candidate-fasta --help
fp-tools-meme-command --help
fp-tools-motif-discovery-plan --help
fp-tools-summarize-motifs --help
fp-tools-score-variants --help
fp-tools-pseudobulk --help
fp-tools-bindetect-replicate-report --help
fp-tools-decompose-competition --help
```

## Minimal Workflow

### 1. atac-correct

```bash
atac-correct \
  --bam test_data/Bcell.bam \
  --genome test_data/genome.fa.gz \
  --peaks test_data/merged_peaks.bed \
  --blacklist test_data/blacklist.bed \
  --outdir examples/atacorrect/atac-correct_test2 \
  --cores 1
```

### 2. score-footprints

```bash
score-footprints \
  --signal examples/atacorrect/atac-correct_test2/Bcell_corrected.bw \
  --regions test_data/merged_peaks.bed \
  --output examples/scorebigwig/score-footprints_test2/Bcell_footprints.bw \
  --cores 1
```

### 3. detect-tf-binding

```bash
detect-tf-binding \
  --motifs test_data/motifs.jaspar \
  --signals test_data/Bcell_footprints.bw test_data/Tcell_footprints.bw \
  --genome test_data/genome.fa.gz \
  --peaks test_data/merged_peaks_annotated.bed \
  --peak-header test_data/merged_peaks_annotated_header.txt \
  --outdir examples/bindetect/detect-tf-binding_output_htmlfix_014 \
  --cond-names Bcell Tcell \
  --cores 1
```

### 4. plot-aggregate

```bash
plot-aggregate \
  --TFBS test_data/IRF1_all.bed \
  --signals test_data/Bcell_corrected.bw \
  --output examples/reports/plotaggregate_control_mode_test.pdf \
  --output_aggregated_scores examples/reports/plotaggregate_control_mode_test_scores.csv
```

## GUI

Start the GUI on a Linux server:

```bash
fp-tools-gui --host 0.0.0.0 --run-dir examples/gui_runs
```

If `--port` is omitted, the launcher picks a free port and prints the exact URL. A fixed port can also be supplied:

```bash
fp-tools-gui --host 0.0.0.0 --port 8891 --run-dir examples/gui_runs
```

The GUI can run directly from forms, load YAML, save YAML, and inspect run history. GUI run metadata and logs are written under `examples/gui_runs/`; ready-to-load YAML examples are in `examples/gui_configs/`.

## YAML Runner

Run a saved config directly from the command line:

```bash
fp-tools-run --config examples/gui_configs/plotaggregate_single.yml
```

GUI-saved YAML files can be rerun the same way. YAML is optional for normal CLI use.

## Extra Features

### detect-tf-binding Replicate Grouping

```bash
detect-tf-binding --motifs test_data/motifs.jaspar \
  --signals test_data/demo_Bcell_rep1_footprints.bw test_data/demo_Bcell_rep2_footprints.bw test_data/demo_Tcell_rep1_footprints.bw test_data/demo_Tcell_rep2_footprints.bw \
  --genome test_data/genome.fa.gz \
  --peaks test_data/merged_peaks_annotated.bed \
  --peak-header test_data/merged_peaks_annotated_header.txt \
  --outdir examples/bindetect/detect-tf-binding_output_synthetic_replicates_demo \
  --cond-names Bcell Bcell Tcell Tcell \
  --cores 40
```

Grouped results are written to `bindetect_results.txt` under the output directory.

### detect-tf-binding Skewness Report

For multi-condition runs, the skewness report is written automatically to:

```text
<outdir>/bindetect_results_skewness_report.pdf
```

Single-condition runs do not produce this report.

### plot-aggregate Control Overlay

```bash
plot-aggregate --TFBS test_data/IRF1_all.bed \
  --signals test_data/Bcell_corrected.bw test_data/Tcell_corrected.bw \
  --signal-labels Bcell Tcell \
  --control-label Bcell \
  --output examples/reports/plotaggregate_control_mode_test.pdf \
  --output_aggregated_scores examples/reports/plotaggregate_control_mode_test_scores.csv
```

### plot-aggregate Directory Input

```bash
plot-aggregate --TFBS examples/plotaggregate_tfbs_dir \
  --signals test_data/Bcell_corrected.bw \
  --output examples/reports/plotaggregate_dirinput_test.pdf \
  --output_aggregated_scores examples/reports/plotaggregate_dirinput_test_scores.csv
```

### plot-aggregate Fixed Grid Layout

```bash
plot-aggregate --TFBS examples/plotaggregate_tfbs_grid/CTCF_Bcell_bound.bed examples/plotaggregate_tfbs_grid/IRF1_Bcell_bound.bed examples/plotaggregate_tfbs_grid/ETS1_Bcell_bound.bed examples/plotaggregate_tfbs_grid/GATA3_Bcell_bound.bed examples/plotaggregate_tfbs_grid/RUNX1_Bcell_bound.bed \
  --signals test_data/Bcell_corrected.bw test_data/Tcell_corrected.bw \
  --signal-labels Bcell_corrected Tcell_corrected \
  --grid 2x5 \
  --output examples/reports/plotaggregate_grid_2x5_test.pdf \
  --output_aggregated_scores examples/reports/plotaggregate_grid_2x5_test_scores.csv
```

### plot-aggregate CSV Exports

```bash
plot-aggregate --TFBS test_data/IRF1_all.bed \
  --signals test_data/Bcell_corrected.bw \
  --output examples/reports/plotaggregate_signals_export_test.pdf \
  --output_aggregated_signals examples/reports/plotaggregate_signals_export_test_signals.csv \
  --output_aggregated_scores examples/reports/plotaggregate_signals_export_test_scores.csv
```

The legacy `--output-csv` alias writes the aggregated signal CSV.

### Larger Motif Databases

Large JASPAR and HOCOMOCO-scale motif databases are supported through the same detect-tf-binding command path. A JASPAR2026-scale run completed with 1019 motifs.

*Yi Lab, 2026.*

## Real-data pseudobulk example

The paper workflow validates `fp-tools-pseudobulk` on public 10x PBMC Multiome ATAC fragments. Rebuild the prepared annotations and TF site sets from official 10x analysis tables, then run pseudobulk grouping with high core counts for per-group compression and bigWig generation:

```bash
python benchmarks/scripts/prepare_10x_pbmc_pseudobulk.py --write-example-archive

fp-tools-pseudobulk \
  --fragments data/public/raw/10x_pbmc/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz \
  --annotations data/public/processed/pseudobulk_pbmc/pbmc_10x_cell_annotations.tsv \
  --group-by cell_type \
  --min-cells 300 \
  --min-fragments 50000 \
  --index-output \
  --write-cutsite-bigwigs \
  --genome-sizes data/public/processed/pseudobulk_pbmc/hg38.chrom.sizes \
  --cores 32 \
  --outdir data/public/processed/pseudobulk_pbmc/run

python paper/scripts/plot_pseudobulk_tf_aggregates.py \
  --manifest data/public/processed/pseudobulk_pbmc/run/pseudobulk_manifest.tsv \
  --tf-site-dir data/public/processed/pseudobulk_pbmc/tf_sites \
  --out-prefix paper/manuscript/figures/supp_pseudobulk_tf_aggregates
```

Keep raw 10x files and generated pseudobulk fragments/bigWigs out of the main repository. Only scripts, manifests, compact source tables, and manuscript figures are tracked here; reusable large example-data archives belong in `oncologylab/fp-tools-data` release assets (https://github.com/oncologylab/fp-tools-data/releases/tag/pbmc-pseudobulk-v1).
