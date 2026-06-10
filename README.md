# fp-tools

`fp-tools` is a standalone footprinting package for ATAC-seq style workflows. It provides command-first tools for bias correction, footprint scoring, differential binding detection, and aggregate signal plotting.

The PyPI distribution is named `fp-tools-bio`; the installed Python package is `fp_tools`.

## Install

```bash
pip install fp-tools-bio
```

## Commands

- `atac-correct`: correct ATAC-seq cutsite signal for Tn5 sequence bias.
- `score-footprints`: calculate footprint, multiscale, sum, mean, or pass-through scores from bigWig signal.
- `detect-tf-binding`: scan motifs, infer bound sites, and compare TF binding across conditions.
- `plot-aggregate`: plot aggregate signal around TFBS or region sets.
- `fp-tools-run`: run optional YAML batch configs.
- `fp-tools-gui`: launch the optional Streamlit GUI wrapper.
- `fp-tools-build-tfbs-features`: build tabular supervised TFBS features from candidates, genome, signals, and labels.
- `fp-tools-train-tfbs-model`: train an optional motif-centric supervised TFBS model from tabular features.
- `fp-tools-predict-tfbs`: apply a saved tabular TFBS model to feature tables.
- `fp-tools-generate-candidates`: nominate motif-free footprint candidates from score bigWigs and peak BEDs.
- `fp-tools-rerank-candidates`: combine candidate, motif, family, and model scores into a ranked site table.
- `fp-tools-export-candidate-fasta`: export candidate-centered FASTA for de novo motif discovery.
- `fp-tools-meme-command`: print a MEME/DREME command for exported candidate FASTA.
- `fp-tools-motif-discovery-plan`: write or run a MEME/DREME, Tomtom, and fp-tools summary command plan.
- `fp-tools-summarize-motifs`: summarize MEME/Tomtom results as TSV and HTML reports.
- `fp-tools-score-variants`: annotate variants with allele checks and footprint/candidate overlaps.
- `fp-tools-pseudobulk`: group single-cell ATAC fragments into pseudobulk fragment files and a manifest.
- `fp-tools-bindetect-replicate-report`: summarize BINDetect comparison effects, p-values, replicate support, and uncertainty.
- `fp-tools-decompose-competition`: decompose multiscale footprint signal into competing TF-scale and nucleosome-scale components.

Direct CLI usage is the primary interface. YAML configs and the GUI are optional wrapper paths and do not replace the plain command-line tools.

Legacy aliases remain available for compatibility: `ATACorrect`, `FootprintScores`, `ScoreBigwig`, `BINDetect`, and `PlotAggregate`.

## Feature Comparison Across the Field

`fp-tools current` describes released functionality. `fp-tools roadmap` describes planned opt-in modules under development, not current released functionality.

Symbols: ✅ native first-class support, ⚠️ partial or indirect support, ❌ absent.

| Feature | fp-tools<br>current | fp-tools<br>roadmap | TOBIAS | HINT-<br>ATAC | PRINT /<br>scPrinter | TAMC | maxATAC | ChromBPNet<br>Suite | REUNION | DENIS |
|---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| Bulk ATAC<br>footprinting | ✅ | ✅ | ✅ | ✅ | ✅ | ⚠️ | ❌ | ❌ | ⚠️ | ⚠️ |
| Tn5 bias<br>correction | ✅ | ✅ | ✅ | ✅ | ✅ | ❌ | ❌ | ✅ | ❌ | ❌ |
| Classical<br>footprint scoring | ✅ | ✅ | ✅ | ✅ | ⚠️ | ⚠️ | ❌ | ❌ | ❌ | ⚠️ |
| Multiscale /<br>nucleosome-aware | ❌ | ✅ | ❌ | ⚠️ | ✅ | ❌ | ❌ | ⚠️ | ❌ | ❌ |
| Supervised<br>TFBS prediction | ❌ | ✅ | ❌ | ❌ | ✅ | ✅ | ✅ | ✅ | ✅ | ❌ |
| Motif-relaxed /<br>motif-free recovery | ❌ | ✅ | ❌ | ❌ | ❌ | ❌ | ❌ | ❌ | ✅ | ✅ |
| De novo<br>motif discovery | ❌ | ✅ | ❌ | ❌ | ✅ | ❌ | ❌ | ✅ | ❌ | ✅ |
| Variant<br>scoring | ❌ | ✅ | ❌ | ❌ | ❌ | ❌ | ✅ | ✅ | ❌ | ❌ |
| scATAC /<br>pseudobulk support | ⚠️ | ✅ | ⚠️ | ⚠️ | ✅ | ⚠️ | ✅ | ⚠️ | ✅ | ✅ |
| Visualization /<br>reporting | ✅ | ✅ | ✅ | ⚠️ | ✅ | ⚠️ | ⚠️ | ✅ | ⚠️ | ✅ |
| GUI / YAML /<br>batch execution | ✅ | ✅ | ❌ | ❌ | ❌ | ❌ | ❌ | ❌ | ❌ | ❌ |

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

Development builds also include an opt-in multiscale depletion mode for comparing narrow TF-scale and broader nucleosome-scale signal structure:

```bash
score-footprints \
  --signal examples/atacorrect/atac-correct_test2/Bcell_corrected.bw \
  --regions test_data/merged_peaks.bed \
  --score multiscale \
  --scales 8 16 24 32 64 100 147 \
  --multiscale-summary max \
  --output examples/scorebigwig/score-footprints_multiscale/Bcell_multiscale_max.bw \
  --output-multiscale-npz examples/scorebigwig/score-footprints_multiscale/Bcell_multiscale_tensor.npz \
  --cores 1
```

The optional `.npz` sidecar stores scale-by-position scores for downstream heatmaps and paper figures:

```bash
python paper/scripts/plot_multiscale_npz.py \
  --multiscale-npz examples/scorebigwig/score-footprints_multiscale/Bcell_multiscale_tensor.npz \
  --out-prefix paper/figures/figure_multiscale_summary
```

`plot-aggregate` can also render the sidecar as a companion scale-by-position aggregate figure while keeping its standard bigWig aggregate plot:

```bash
plot-aggregate --TFBS test_data/IRF1_all.bed \
  --signals examples/scorebigwig/score-footprints_multiscale/Bcell_multiscale_max.bw \
  --multiscale-npz examples/scorebigwig/score-footprints_multiscale/Bcell_multiscale_tensor.npz \
  --output examples/reports/plotaggregate_multiscale.pdf \
  --output-multiscale-aggregate examples/reports/plotaggregate_multiscale_tensor.pdf
```

### Optional supervised TFBS model

Build a feature table from candidate sites, then train a motif-centric tabular classifier with a binary `label` column and numeric features such as motif score, footprint score, and multiscale summaries:

```bash
fp-tools-build-tfbs-features \
  --candidates examples/bindetect/motif_free_candidates.bed \
  --genome test_data/genome.fa.gz \
  --signals examples/scorebigwig/score-footprints_multiscale/Bcell_multiscale_max.bw \
  --signal-labels footprint \
  --labels-bed benchmarks/labels/ctcf_chip_positive.bed \
  --out benchmarks/features/ctcf_train.tsv

fp-tools-train-tfbs-model \
  --train-table benchmarks/features/ctcf_train.tsv \
  --model-out models/ctcf_tabular.pkl \
  --metrics-out models/ctcf_tabular_metrics.tsv

fp-tools-predict-tfbs \
  --model models/ctcf_tabular.pkl \
  --features benchmarks/features/ctcf_test.tsv \
  --out benchmarks/predictions/ctcf_test_predictions.tsv
```

### Optional motif-free candidates

Nominate high-scoring candidate binding intervals from footprint or multiscale summary signal without requiring motif hits:

```bash
fp-tools-generate-candidates \
  --signal examples/scorebigwig/score-footprints_multiscale/Bcell_multiscale_max.bw \
  --peaks test_data/merged_peaks.bed \
  --out examples/bindetect/motif_free_candidates.bed \
  --min-score 1.0 \
  --top-n-per-region 5
```

To include motif-relaxed candidates, scan motifs externally with a weaker PWM threshold, then merge the resulting BED with the signal-local maxima. The motif intervals are filtered to peaks and scored from the same signal track:

```bash
fp-tools-generate-candidates \
  --signal examples/scorebigwig/score-footprints_multiscale/Bcell_multiscale_max.bw \
  --peaks test_data/merged_peaks.bed \
  --motif-sites benchmarks/motifs/ctcf_relaxed_p1e-3.bed \
  --motif-signal-math max \
  --out examples/bindetect/motif_relaxed_candidates.bed
```

### Optional motif-relaxed reranking

Combine motif-free candidates, motif/family annotations, and supervised probabilities into a ranked table:

```bash
fp-tools-rerank-candidates \
  --sites benchmarks/predictions/ctcf_test_predictions.tsv \
  --score-columns binding_probability candidate_score footprint_mean \
  --weights 2.0 1.0 1.0 \
  --family-map benchmarks/motif_families.tsv \
  --motif-column motif_id \
  --out benchmarks/predictions/ctcf_test_reranked.tsv
```

### Optional de novo motif prep

Export candidate-centered sequences for external motif discovery tools:

```bash
fp-tools-export-candidate-fasta \
  --candidates examples/bindetect/motif_free_candidates.bed \
  --genome test_data/genome.fa.gz \
  --out examples/bindetect/motif_free_candidates.fa \
  --flank 50

fp-tools-meme-command \
  --fasta examples/bindetect/motif_free_candidates.fa \
  --outdir examples/bindetect/denovo_motifs \
  --method meme \
  --extra-args -dna -mod zoops

fp-tools-motif-discovery-plan \
  --fasta examples/bindetect/motif_free_candidates.fa \
  --outdir examples/bindetect/denovo_motifs \
  --method meme \
  --known-motifs test_data/individual_motifs/MA0050.2.jaspar \
  --extra-args -dna -mod zoops -nmotifs 10
```

The generated plan ends by running `fp-tools-summarize-motifs`, which writes TSV plus self-contained HTML reports with inline consensus-logo panels.

### Optional variant scoring scaffold

Annotate variants against the genome and motif-free candidate intervals as a first footprint-aware variant table:

```bash
fp-tools-score-variants \
  --variants variants.bed \
  --genome test_data/genome.fa.gz \
  --candidate-scores examples/bindetect/motif_free_candidates.bed \
  --sequence-flank 20 \
  --kmer-size 3 \
  --motifs test_data/individual_motifs/MA0050.2.jaspar \
  --motif-flank 30 \
  --tfbs-model models/ctcf_tabular.pkl \
  --out examples/bindetect/scored_variants.tsv
```

The output includes deterministic ref/alt sequence-context features such as GC shift, allele length delta, exact k-mer disruption, optional best-motif PWM score deltas, and optional tabular model probability deltas.

### Optional pseudobulk fragments

Group single-cell ATAC fragments by annotation columns and write a manifest for downstream bulk processing:

```bash
fp-tools-pseudobulk \
  --fragments pbmc_fragments.tsv.gz \
  --annotations cell_annotations.tsv \
  --group-by donor,cell_type \
  --min-cells 200 \
  --min-fragments 50000 \
  --write-downstream-commands \
  --genome-sizes refs/hg38.chrom.sizes \
  --cores 16 \
  --outdir pseudobulk_pbmc
```

The optional downstream command script contains reproducible `bedtools`, `bedGraphToBigWig`, and `samtools` steps for kept pseudobulk groups.

### Optional replicate-aware BINDetect report

Summarize a BINDetect `*_results.txt` table into per-comparison effect sizes, p-values, replicate support, and replicate-aware uncertainty:

```bash
fp-tools-bindetect-replicate-report \
  --results BINDetect_output/bindetect_results.txt \
  --replicate-map replicate_counts.tsv \
  --out replicate_report.tsv \
  --summary-out replicate_summary.tsv \
  --figure-out replicate_report.png
```

Replicate counts can come from a `condition`/`n_replicates` (or `condition`/`replicate`) map, from `<condition>_n_replicates` columns in the results table, or fall back to single-replicate. Each comparison is labelled `single-replicate` or `replicate-supported`, and reported with a p-value-derived standard error, shrinkage weight, and confidence interval.

### Optional competition-aware footprint decomposition

Decompose a multiscale NPZ sidecar (`score-footprints --output-multiscale-npz`) into competing TF-scale and nucleosome-scale footprint components per region:

```bash
fp-tools-decompose-competition \
  --npz sample_multiscale.npz \
  --tf-band 3,30 \
  --nucleosome-band 120,200 \
  --out competition.tsv \
  --summary-out competition_summary.tsv \
  --figure-out competition.png
```

Each region is partitioned into TF-only, nucleosome-only, and shared (competing) signal, with a `competition_index` and a `dominant_component` label (`tf`, `nucleosome`, `competing`, or `none`).

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
