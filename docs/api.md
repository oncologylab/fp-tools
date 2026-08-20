---
hide:
  - navigation
---

# API Reference

Direct CLI commands are the primary interface. Each reference includes a method summary, practical example, primary inputs, outputs, and the complete command options.

| Command | Purpose |
| --- | --- |
| [`prepare-atac`](#prepare-atac) | **Linux CLI or Linux container only.** Download public ATAC-seq reads or use local FASTQ files, then trim, align, filter, call peaks, calculate alignment coverage, and write QC files. The GUI and native macOS/Windows installations start from filtered BAM/BAI and peak BED files. |
| [`bulk-footprinting`](#bulk-footprinting) | Run bulk ATAC-seq from BAM/BAI and peak BED inputs through interactive reports. |
| [`atac-correct`](#atac-correct) | Estimate Tn5 sequence bias from aligned ATAC-seq fragments and subtract the expected bias contribution from the observed cut-site signal. Run this before footprint scoring. |
| [`call-footprints`](#call-footprints) | Calculate a continuous footprint score from bias-corrected cut-site signal within accessible regions. Higher local depletion relative to flanking signal produces stronger footprint evidence. |
| [`match-motifs`](#match-motifs) | Scan accessible regions for motif instances, measure the footprint score at each instance, and classify sample-specific bound and unbound sites. Run this when motif locations and per-sample motif summaries are needed. |
| [`diff-footprints`](#diff-footprints) | Compare motif-associated footprint scores across conditions or between user-defined region sets measured in the same sample(s). |
| [`normalize-bigwig`](#normalize-bigwig) | Scale corrected cut-site signals using statistics measured over the same background regions. Use this optional step when samples require an explicitly shared signal scale before downstream scoring or plotting. |
| [`plot-aggregate`](#plot-aggregate) | Plot average signal around motif sites or other genomic regions as a static figure or interactive HTML report. |
| [`review-multi-comparisons`](#review-multi-comparisons) | Combine differential-footprint reports as a scalable browser bundle or one self-contained HTML report. |
| [`run-yaml-workflow`](#run-yaml-workflow) | Run one or more fp-tools jobs from a reusable YAML configuration. |
| [`fp-tools-gui`](#fp-tools-gui) | Launch the browser interface for configuring and running fp-tools commands. |
| [`fp-tools-runtime`](#fp-tools-runtime) | Inspect, install, or repair the private external-tool runtime. Linux provides raw-read and de novo motif components; macOS and Windows provide the optional de novo motif component. |
| [`discover-motifs`](#discover-motifs) | Prepare or run de novo motif discovery from candidate footprint intervals or an existing FASTA file. |
| [`summarize-motifs`](#summarize-motifs) | Summarize MEME, STREME, DREME, and Tomtom results in a compact report. |
| [`pseudobulk-fragments`](#pseudobulk-fragments) | Group single-cell ATAC fragments by a cell-annotation column to create pseudobulk inputs. |
| [`find-signature-fp`](#find-signature-fp) | Calculate and plot per-cell footprint signatures from completed pseudobulk or motif analyses. |
| [`sc-footprinting`](#sc-footprinting) | Run grouping, bias correction, footprint scoring, motif analysis, and per-cell signature reporting for single-cell ATAC-seq data. |

## `prepare-atac`

<div class="fp-api-card" markdown="1">

**Linux CLI or Linux container only.** Download public ATAC-seq reads or use
local FASTQ files, then trim, align, filter, call peaks, calculate alignment
coverage, and write QC files. The GUI and native macOS/Windows installations
start from filtered BAM/BAI and peak BED files.

**Example command**

```bash
prepare-atac --samples metadata.tsv --genome hg38 --outdir project
```

**Primary inputs**

- `--samples` — TSV or CSV sample sheet containing `sample`, `condition`, and either paired `fastq_1`/`fastq_2` paths or URLs. See the [bulk workflow guide](get-started/workflows/bulk-atac-seq.md).
- `--genome` — packaged `hg38` or `mm10` reference label, or a custom label used with explicit reference options.
- `--outdir` — project directory represented by `{project}` below.

Repeated rows with the same `sample`, `condition`, and `replicate` combine
technical sequencing runs. Different `sample` values sharing a `condition` are
biological replicates.

**Main outputs**

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

**Complete options**

```text
usage: prepare-atac [-h] [--samples SAMPLES] [--genome GENOME]
                    [--outdir OUTDIR] [--config CONFIG]
                    [--profile {modern,homer-atac}] [--id-column ID_COLUMN]
                    [--sample-column SAMPLE_COLUMN]
                    [--condition-column CONDITION_COLUMN]
                    [--include INCLUDE [INCLUDE ...]]
                    [--reference-dir REFERENCE_DIR] [--fasta FASTA]
                    [--bowtie2-index BOWTIE2_INDEX] [--blacklist BLACKLIST]
                    [--tss TSS] [--macs-genome-size MACS_GENOME_SIZE]
                    [--cores CORES]
                    [--max-parallel-samples MAX_PARALLEL_SAMPLES]
                    [--memory-gb MEMORY_GB] [--keep-intermediates]
                    [--no-resume] [--fail-fast] [--dry-run] [--doctor]
                    [--write-default-config PATH]
                    [--runtime {auto,managed,system,container}]

Download single- or paired-end ATAC-seq reads and prepare filtered BAM, peak
BED, sequencing-depth-normalized alignment coverage bigWig, and QC files.

options:
  -h, --help            show this help message and exit
  --samples SAMPLES     TSV/CSV sample metadata with an ID/run_accession or
                        fastq_1 column.
  --genome GENOME       Named genome (hg38/mm10) or custom genome label.
  --outdir OUTDIR       Output project directory.
  --config CONFIG       Optional preprocessing YAML; CLI values override
                        packaged defaults.
  --profile {modern,homer-atac}
                        Processing method: modern uses fastp, samtools, and
                        MACS3; homer-atac uses Trim Galore, Picard, and HOMER
                        (default: modern).
  --id-column ID_COLUMN
                        Explicit accession column name.
  --sample-column SAMPLE_COLUMN
                        Explicit sample-name column.
  --condition-column CONDITION_COLUMN
                        Explicit condition column.
  --include INCLUDE [INCLUDE ...]
                        Only process these sample names or accessions.
  --reference-dir REFERENCE_DIR
                        Reference cache root (default: ~/.cache/fp-
                        tools/references).
  --fasta FASTA         Custom reference FASTA.
  --bowtie2-index BOWTIE2_INDEX
                        Existing Bowtie2 index prefix.
  --blacklist BLACKLIST
                        Custom blacklist BED.
  --tss TSS             Optional TSS BED for enrichment QC.
  --macs-genome-size MACS_GENOME_SIZE
                        MACS3 genome size or hs/mm shorthand for custom
                        genomes.
  --cores CORES         Total core budget.
  --max-parallel-samples MAX_PARALLEL_SAMPLES
                        Maximum samples processed concurrently (default:
                        config value).
  --memory-gb MEMORY_GB
                        Enforced total memory budget in GiB; reserves 8 GiB
                        for the host.
  --keep-intermediates  Keep trimmed FASTQs and intermediate alignment files
                        under each sample .work directory.
  --no-resume           Recompute completed samples even when fingerprints
                        match.
  --fail-fast           Stop after the first failed sample.
  --dry-run             Validate and list the planned samples without
                        downloading or processing.
  --doctor              Report external preprocessing dependencies and exit.
  --write-default-config PATH
                        Write the fully documented default YAML and exit.
  --runtime {auto,managed,system,container}
                        External-tool runtime: auto/managed provisions the
                        pinned fp-tools runtime, system uses PATH, and
                        container uses the complete image (default: auto).
```

</div>

## `bulk-footprinting`

<div class="fp-api-card" markdown="1">

Run bulk ATAC-seq from BAM/BAI and peak BED inputs through interactive reports.

The [bulk workflow guide](get-started/workflows/bulk-atac-seq.md) provides a runnable
HepG2-versus-K562 ENCODE example.

**Example command**

```bash
bulk-footprinting --sample-table samples.tsv --comparison-table comparisons.tsv --genome hg38.fa.gz \
  --outdir project --cores 8
```

**Primary inputs**

- `--sample-table` — sample, condition, coordinate-sorted BAM, and peak BED columns.
- `--comparison-table` — comparison, condition 1, and condition 2 columns.
- `--genome` — reference FASTA matching the BAM and peak coordinates.
- `--outdir` — project output directory.
- `--cores` — total worker cores.

**Main outputs**

`{project}` is the `--outdir`, `{sample}` comes from the sample table, and
`{comparison}` comes from the comparison table:

| Path | Meaning |
| --- | --- |
| `{project}/samples/{sample}/atac_correct/{sample}_corrected.bw` | Bias-corrected cut-site signal. |
| `{project}/samples/{sample}/footprints/{sample}_footprints.bw` | Footprint score signal. |
| `{project}/samples/{sample}/match_motifs/motif_matches_results.txt` | Per-sample motif summary and binding calls. |
| `{project}/comparisons/{comparison}/diff_footprints_results.txt` | Motif-level differential statistics. |
| `{project}/comparisons/{comparison}/diff_footprints_{cond1}_{cond2}.html` | Portable interactive comparison report. |
| `{project}/reports/review_multi_comparisons/index.html` | Static browser combining every requested comparison. |
| `{project}/reports/review_multi_comparisons.html` | Aggregate-free portable review written when standalone HTML review mode is selected. |
| `{project}/logs/bulk_footprinting/bulk_footprinting_commands.sh` | Exact commands generated for the workflow stages. |
| `{project}/logs/bulk_footprinting/{stage}.stdout.log` and `{stage}.stderr.log` | Stage-specific logs for troubleshooting. |

FASTQ preprocessing is intentionally separate. Linux users can run
[`prepare-atac`](#prepare-atac) first, then provide its generated
`metadata/samples.tsv` to this command.

**Complete options**

```text
usage: bulk-footprinting [-h] --sample-table SAMPLE_TABLE --comparison-table
                         COMPARISON_TABLE --genome GENOME
                         [--blacklist BLACKLIST] [--motifs [MOTIFS ...]]
                         [--motif-db MOTIF_DB]
                         [--normalization {none,condition-quantile,sample-quantile}]
                         [--plot-aggregate {sig,all,top,off}]
                         [--review-format {auto,bundle,standalone,none}]
                         --outdir OUTDIR [--cores CORES] [--resume] [--force]
                         [--dry-run] [--fail-fast]
                         [--runtime {auto,managed,system,container}]

Run the complete bulk ATAC-seq footprinting workflow for explicit comparisons.

options:
  -h, --help            show this help message and exit
  --sample-table SAMPLE_TABLE
                        TSV with sample, condition, BAM, and peak BED columns.
  --comparison-table COMPARISON_TABLE
                        TSV with comparison, cond1, and cond2 columns.
  --genome GENOME       Reference FASTA matching the BAM and peak coordinates.
  --blacklist BLACKLIST
                        Optional blacklist BED used during bias correction.
  --motifs [MOTIFS ...]
                        Optional motif files.
  --motif-db MOTIF_DB   Built-in motif database (default:
                        jaspar2026_vertebrates).
  --normalization {none,condition-quantile,sample-quantile}
                        Differential-stage normalization (default: none).
  --plot-aggregate {sig,all,top,off}
                        Aggregate profiles generated by diff-footprints
                        (default: all).
  --review-format {auto,bundle,standalone,none}
                        Combined-review output; auto uses standalone when
                        aggregation is off and bundle otherwise (default:
                        auto).
  --outdir OUTDIR       Project output directory.
  --cores CORES         Total worker cores passed to each stage (default: 1).
  --resume              Skip stages whose expected outputs are complete.
  --force               Rerun stages even when outputs already exist.
  --dry-run             Validate inputs and print the commands without running
                        them.
  --fail-fast           Stop at the first failed stage.
  --runtime {auto,managed,system,container}
                        External-tool runtime: auto/managed provisions the
                        pinned fp-tools runtime, system uses PATH, and
                        container uses the complete image (default: auto).
```

</div>

## `atac-correct`

<div class="fp-api-card" markdown="1">

Estimate Tn5 sequence bias from aligned ATAC-seq fragments and subtract the
expected bias contribution from the observed cut-site signal. Run this before
footprint scoring.

**Example command**

```bash
atac-correct --sample-table project/metadata/samples.tsv --genome hg38.fa.gz --blacklist hg38.blacklist.bed --outdir project
```

**Primary inputs**

- `--sample-table` — TSV with `sample`, `condition`, `bam`, and `peaks`; BAM indexes must be adjacent to the BAMs.
- `--genome` — reference FASTA whose chromosome names and assembly match every BAM and peak BED.
- `--blacklist` — BED intervals excluded from bias estimation and corrected output.
- `--outdir` — project directory represented by `{project}` below.

**Main outputs**

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

**Complete options**

```text
usage: atac-correct [-h] [--bams [<bam> ...]] [--fragments [<fragments.tsv.gz> ...]]
                    [-g <fasta>] [-p [<bed> ...]] [--regions-in <bed>]
                    [--regions-out <bed>] [--blacklist <bed>] [--extend <int>]
                    [--split-strands] [--norm-off] [--write-tracks [<track> ...]]
                    [--track-off [<track> ...]] [--skip-qc]
                    [--scale-corrected {auto,none,q95}] [--scale-background <bed>]
                    [--scale-corrected-bigwigs [<bigwig> ...]]
                    [--scale-target {median,mean}] [--scale-chrom-sizes <chrom.sizes>]
                    [--merged-peaks-out <bed>] [--drop-chroms [<chrom> ...]]
                    [--k_flank <int>] [--read_shift <int> <int>] [--bg_shift <int>]
                    [--window <int>] [--score_mat <mat>] [--bias-pkl <obj>]
                    [--prefix <prefix>] [--sample-names [<name> ...]]
                    [--sample-table <tsv>] [--layout {custom,project}]
                    [--sample-output-root <directory>] [--outdir <directory>]
                    [--cores <int>] [--sample-workers <int>] [--split <int>]
                    [--verbosity <int>]

__________________________________________________________________________________________

                                  fp-tools atac-correct
__________________________________________________________________________________________

atac-correct corrects ATAC-seq cutsite signal for Tn5 sequence bias.

Usage:
atac-correct --bams <reads.bam> [<more_reads.bam> ...] --genome <genome.fa> --peaks
<merged_peaks.bed> [<sample_peaks.bed> ...]

Output files:
- <outdir>/<sample>/<sample>_corrected.bw for multi-BAM runs
- optional auxiliary tracks with --write-tracks
- <outdir>/<prefix>_atacorrect.pdf

------------------------------------------------------------------------------------------

Required arguments:
  --bams [<bam> ...]               One or more .bam files containing reads to be corrected
  --fragments [<fragments.tsv.gz> ...]
                                   One or more 10x-style fragment files; uses cut sites
                                   directly without creating pseudo-BAMs
  -g <fasta>, --genome <fasta>     A .fasta-file containing whole genomic sequence
  -p [<bed> ...], --peaks [<bed> ...]
                                   One shared merged peak BED, or multiple per-sample peak
                                   BEDs to merge internally

Optional arguments:
  --regions-in <bed>               Input regions for estimating bias (default: regions not
                                   in peaks.bed)
  --regions-out <bed>              Output regions (default: peaks.bed)
  --blacklist <bed>                Blacklisted regions in .bed-format (default: None)
  --extend <int>                   Extend output regions with basepairs
                                   upstream/downstream (default: 100)
  --split-strands                  Write out tracks per strand
  --norm-off                       Switches off normalization based on number of reads
  --write-tracks [<track> ...]     Cut-site signal bigWigs to write (default: corrected;
                                   use all for corrected, observed/uncorrected, sequence-
                                   bias, and expected signals)
  --track-off [<track> ...]        Compatibility option to switch off individual bigWig
                                   tracks after --write-tracks is resolved
  --skip-qc                        Skip atac-correct diagnostic PDF and pre/post bias
                                   verification counts. Corrected bigWig output is
                                   unchanged.
  --scale-corrected {auto,none,q95}
                                   Optionally q95-scale corrected bigWigs after
                                   correction. In auto mode this runs only when --scale-
                                   corrected-bigwigs has more than one track (default:
                                   auto)
  --scale-background <bed>         Shared BED regions used to estimate q95 scaling for
                                   --scale-corrected
  --scale-corrected-bigwigs [<bigwig> ...]
                                   Corrected bigWigs to q95-scale together. Include the
                                   current sample's corrected bigWig or omit to scale only
                                   the current output
  --scale-target {median,mean}     Across-sample q95 target for --scale-corrected
                                   (default: median)
  --scale-chrom-sizes <chrom.sizes>
                                   Optional chromosome sizes file for scaled bigWig output
                                   validation
  --merged-peaks-out <bed>         Path for internally merged peak BED when multiple
                                   --peaks files are supplied (default:
                                   <outdir>/merged_peaks.bed)
  --drop-chroms [<chrom> ...]      Drop any chromosomes in the list from the correction.
                                   The default is to drop the mitochrondrial chromosome.
                                   Default: ['chrM', 'chrMT', 'M', 'MT', 'Mito']

Advanced atac-correct arguments (no need to touch):
  --k_flank <int>                  Flank +/- of cutsite to estimate bias from (default:
                                   12)
  --read_shift <int> <int>         Read shift for forward and reverse reads (default: 4
                                   -5)
  --bg_shift <int>                 Read shift for estimation of background frequencies
                                   (default: 100)
  --window <int>                   Window size for calculating expected signal (default:
                                   100)
  --score_mat <mat>                Type of matrix to use for bias estimation (PWM/DWM)
                                   (default: DWM)
  --bias-pkl <obj>                 Path to a pre-calculated AtacBias.pkl-object, as output
                                   from a previous atac-correct run (default: None). Can
                                   be used to bypass the internal bias estimation.

Run arguments:
  --prefix <prefix>                Prefix for output files in single-BAM runs (default:
                                   BAM filename stem)
  --sample-names [<name> ...]      Sample labels for --bams (default: BAM filename stems)
  --sample-table <tsv>             Project sample table with sample, condition, bam, and
                                   peaks columns
  --layout {custom,project}        Use fp-tools standard project output layout under
                                   --outdir (default: project when --sample-table is
                                   provided)
  --sample-output-root <directory>
                                   Sample output root; writes each sample under
                                   <root>/<sample>/atac_correct, typically
                                   <project>/samples
  --outdir <directory>             Output directory for files (default: current working
                                   directory)
  --cores <int>                    Number of cores to use for computation (default: all
                                   available cores)
  --sample-workers <int>           Number of samples to process concurrently for multi-BAM
                                   runs (default: auto when --cores is set)
  --split <int>                    Split of multiprocessing jobs (default: 100)
  --verbosity <int>                Level of output logging (0: silent, 1: errors/warnings,
                                   2: info, 3: stats, 4: debug, 5: spam) (default: 3)
```

</div>

## `call-footprints`

<div class="fp-api-card" markdown="1">

Calculate a continuous footprint score from bias-corrected cut-site signal
within accessible regions. Higher local depletion relative to flanking signal
produces stronger footprint evidence.

**Example command**

```bash
call-footprints --signals A_corrected.bw B_corrected.bw --sample-names A B \
  --regions merged_peaks.bed --sample-output-root project/samples
```

**Primary inputs**

- `--signals` — one bias-corrected cut-site signal bigWig per sample.
- `--sample-names` — labels in the same order as `--signals`.
- `--regions` — BED intervals in which scores are calculated; normally the project merged, filtered peaks.
- `--sample-output-root` — root represented by `{sample_root}` below.

**Main outputs**

| Path | Meaning |
| --- | --- |
| `{sample_root}/{sample}/footprints/{sample}_footprints.bw` | Base-resolution footprint score bigWig used by `match-motifs` and `diff-footprints`. |
| `{sample_root}/{sample}/footprints/{sample}_candidate_footprints.bed` | Optional local score maxima for de novo motif discovery; written with `--call-candidates`. |
| user-selected `*.npz` | Optional compressed scale-by-position score arrays when `--score multiscale` is used with an NPZ output option. |

In direct mode, `--output result.bw` writes exactly `result.bw`; multiple
signals written through `--outdir {outdir}` use
`{outdir}/{signal_stem}_footprints.bw`.

**Complete options**

```text
usage: call-footprints [-h] [-s <bigwig>] [--signals [<bigwig> ...]] [-o <bigwig>]
                       [--outputs [<bigwig> ...]] [-r <bed>] [--score <score>]
                       [--absolute] [--extend <int>] [--smooth <int>]
                       [--min-limit <float>] [--max-limit <float>] [--scales [<int> ...]]
                       [--multiscale-summary <method>] [--output-multiscale-npz <npz>]
                       [--output-multiscale-npzs [<npz> ...]] [--output-bed <bed>]
                       [--output-beds [<bed> ...]] [--output-bed-dir <directory>]
                       [--call-candidates] [--top-n <int>] [--min-score <float>]
                       [--call-width <bp>] [--min-distance <bp>] [--fp-min <int>]
                       [--fp-max <int>] [--flank-min <int>] [--flank-max <int>]
                       [--footprint-kernel {fast,reference}] [--window <int>]
                       [--sample-names [<name> ...]] [--sample-table <tsv>]
                       [--layout {custom,project}] [--sample-output-root <directory>]
                       [--outdir <directory>] [--cores <int>] [--sample-workers <int>]
                       [--split <int>] [--verbosity <int>]

__________________________________________________________________________________________

                                 fp-tools call-footprints
__________________________________________________________________________________________

call-footprints calculates footprint, sum, mean, or pass-through scores from one or more
bigWig signals and can optionally call ranked footprint candidate intervals.

Usage: call-footprints --signals <cutsites.bw> [<more_cutsites.bw> ...] --regions
<regions.bed> --outdir <output_dir>
   or: call-footprints --signal <cutsites.bw> --regions <regions.bed> --output <output.bw>

Output:
- footprint score bigWig(s)
- optional candidate BED from --output-bed/--output-beds for de novo motif discovery

------------------------------------------------------------------------------------------

Required arguments:
  -s <bigwig>, --signal <bigwig>        A .bw file of ATAC-seq cutsite signal
  --signals [<bigwig> ...]              One or more .bw files of ATAC-seq cutsite signal
  -o <bigwig>, --output <bigwig>        Full path to output bigwig
  --outputs [<bigwig> ...]              Output bigWig path per --signals input
  -r <bed>, --regions <bed>             Genomic regions to run footprinting within

Optional arguments:
  --score <score>                       Type of scoring to perform on cutsites
                                        (footprint/sum/mean/none/multiscale) (default:
                                        footprint)
  --absolute                            Convert bigwig signal to absolute values before
                                        calculating score
  --extend <int>                        Extend input regions with bp (default: 100)
  --smooth <int>                        Smooth output signal by mean in <bp> windows
                                        (default: no smoothing)
  --min-limit <float>                   Limit input bigwig value range (default: no lower
                                        limit)
  --max-limit <float>                   Limit input bigwig value range (default: no upper
                                        limit)

Parameters for score == multiscale:
  --scales [<int> ...]                  Window sizes for multiscale depletion scoring
                                        (default: 8 16 24 32 64 100 147)
  --multiscale-summary <method>         How to collapse scale-specific scores into the
                                        output bigWig (default: max)
  --output-multiscale-npz <npz>         Optional compressed NumPy sidecar with per-region
                                        scale-by-position multiscale scores (only for
                                        --score multiscale)
  --output-multiscale-npzs [<npz> ...]  Output multiscale NPZ sidecar per --signals input

Optional footprint candidate BED calling:
  --output-bed <bed>                    Optional BED-like file of genomic coordinates for
                                        footprint peaks used by de novo motif discovery
  --output-beds [<bed> ...]             Candidate BED path per --signals input
  --output-bed-dir <directory>          Directory for candidate BED files derived from
                                        --signals names
  --call-candidates                     In project/sample-output-root mode, also write
                                        candidate footprint BEDs for de novo motif
                                        discovery
  --top-n <int>                         Keep only the top N footprint calls by score
                                        (default: keep all)
  --min-score <float>                   Minimum footprint score for candidate BED calls
                                        (default: no threshold)
  --call-width <bp>                     Width of candidate BED intervals centered on local
                                        maxima (default: 50)
  --min-distance <bp>                   Minimum distance between retained local footprint
                                        centers within a region (default: 20)

Parameters for score == footprint:
  --fp-min <int>                        Minimum footprint width (default: 20)
  --fp-max <int>                        Maximum footprint width (default: 50)
  --flank-min <int>                     Minimum range of flanking regions (default: 10)
  --flank-max <int>                     Maximum range of flanking regions (default: 30)
  --footprint-kernel {fast,reference}   Footprint scoring kernel (default: fast; use
                                        reference for the original scalar implementation)

Parameters for score == sum:
  --window <int>                        The window for calculation of sum (default: 100)

Run arguments:
  --sample-names [<name> ...]           Sample labels for --signals when using project
                                        layout
  --sample-table <tsv>                  Project sample table with sample, condition, bam,
                                        and peaks columns
  --layout {custom,project}             Use fp-tools standard project output layout under
                                        --outdir (default: project when --sample-table is
                                        provided)
  --sample-output-root <directory>      Sample output root; writes each sample under
                                        <root>/<sample>/footprints, typically
                                        <project>/samples
  --outdir <directory>                  Output directory used with --signals when
                                        --outputs is not supplied
  --cores <int>                         Number of cores to use for computation (default:
                                        all available cores)
  --sample-workers <int>                Number of input signals to process concurrently
                                        for multi-signal runs (default: auto when --cores
                                        is set)
  --split <int>                         Split of multiprocessing jobs (default: 100)
  --verbosity <int>                     Level of output logging (0: silent, 1:
                                        errors/warnings, 2: info, 3: stats, 4: debug, 5:
                                        spam) (default: 3)
```

</div>

## `match-motifs`

<div class="fp-api-card" markdown="1">

Scan accessible regions for motif instances, measure the footprint score at
each instance, and classify sample-specific bound and unbound sites. Run this
when motif locations and per-sample motif summaries are needed.

**Example command**

```bash
match-motifs --signals A_footprints.bw B_footprints.bw --sample-names A B --genome hg38.fa.gz \
  --peaks merged_peaks.bed --motif-db jaspar2026_vertebrates --sample-output-root project/samples
```

**Primary inputs**

- `--signals` — one footprint score bigWig per sample.
- `--sample-names` — sample labels in the same order as `--signals`.
- `--genome` — assembly-matched reference FASTA used to scan motif sequences.
- `--peaks` — accessible-region BED searched for motif instances.
- `--motif-db` — packaged motif collection; the example uses JASPAR 2026 vertebrates.
- `--sample-output-root` — root represented by `{sample_root}` below.

**Main outputs**

For each `{sample}`, the default output directory is
`{sample_root}/{sample}/match_motifs/`:

| Path | Meaning |
| --- | --- |
| `motif_matches_results.txt` | Tab-separated motif summary with site counts and per-sample mean scores. |
| `motif_matches_results.xlsx` | Excel copy of the motif summary unless `--skip-excel` is used. |
| `motif_matches_distances.txt` | Motif-similarity distances used for motif clustering. |
| `motif_matches_replicate_motif_score_matrix.tsv` | Motif-by-sample footprint score matrix when multiple samples are analyzed together. |
| `cache/motif_sites.tsv.gz` | Compact scanned motif-site cache reusable by differential analysis. |
| `cache/background_scores.tsv.gz` | Compact background-score cache. |
| `{motif}/beds/{motif}_{sample}_all.bed` | All scanned instances for one motif. |
| `{motif}/beds/{motif}_{sample}_bound.bed` | Instances classified as bound in the sample. |
| `{motif}/beds/{motif}_{sample}_unbound.bed` | Instances classified as unbound in the sample. |

`{motif}` follows the selected `--naming` convention, such as
`CTCF_MA0139.2`. `--motif-outputs summary` omits the per-motif BED files but
keeps the summary and reusable caches.

**Complete options**

```text
usage: match-motifs [-h] [--signals [<bigwig> ...]] [--peaks <bed>] [--genome <fasta>]
                    [--motifs [<motifs> ...]] [--motif-db <name>] [--list-motif-dbs]
                    [--sample-names [<name> ...]] [--cond-names [<name> ...]]
                    [--sample-table <tsv>] [--layout {custom,project}]
                    [--match-scan-mode {auto,shared,per-sample}]
                    [--sample-output-root <directory>] [--peak-header <file>]
                    [--naming <string>] [--motif-pvalue <float>] [--bound-pvalue <float>]
                    [--cluster-threshold <float>] [--pseudo <float>] [--skip-excel]
                    [--output-peaks <bed>] [--norm-off]
                    [--normalization {condition-quantile,sample-quantile,none}]
                    [--aggregate-signals [<bigwig> ...]]
                    [--plot-aggregate {sig,all,top,off}] [--plot-aggregate-top-n <int>]
                    [--plot-aggregate-motifs <motif> [<motif> ...]]
                    [--default-aggregate-plots <int>]
                    [--aggregate-pvalue-threshold <float>] [--aggregate-flank <bp>]
                    [--aggregate-normalization {match,none,sample-quantile,size-factor}]
                    [--aggregate-site-set {all,bound}]
                    [--motif-outputs {auto,summary,full}] [--report-label <text>]
                    [--outdir <directory>] [--prefix <prefix>] [--cores <int>]
                    [--sample-workers <int>] [--split <int>] [--debug] [--verbosity <int>]

__________________________________________________________________________________________

                                  fp-tools match-motifs
__________________________________________________________________________________________

match-motifs scans motifs in open chromatin regions for one or more footprint score tracks
and infers sample-specific bound and unbound motif sites.

Usage:
match-motifs --signals <footprints.bw> [<more_footprints.bw> ...] --genome <genome.fasta>
--peaks <peaks.bed> [--motif-db jaspar2026_vertebrates | --motifs <motifs.txt>]

Output files:
- <outdir>/<prefix>_results.{txt,xlsx}
- <outdir>/<prefix>_distances.txt
- <outdir>/cache/* compact reuse caches
- <outdir>/<TF>/beds/*_all.bed, *_bound.bed, and *_unbound.bed by default
- optional <outdir>/<TF>/<TF>_overview.{txt,xlsx} with --motif-outputs full

------------------------------------------------------------------------------------------

Required arguments:
  --signals [<bigwig> ...]         One or more footprint score bigWigs (.bigwig format)
  --peaks <bed>                    Peaks.bed containing open chromatin regions
  --genome <fasta>                 Genome .fasta file

Optional arguments:
  --motifs [<motifs> ...]          Motif file(s) in pfm/jaspar/meme/transfac format; if
                                   omitted, the built-in JASPAR 2026 vertebrates set is
                                   used
  --motif-db <name>                Built-in motif database to use or add to --motifs
                                   (default when --motifs is omitted:
                                   jaspar2026_vertebrates)
  --list-motif-dbs                 List available built-in motif databases and exit
  --sample-names [<name> ...]      Sample labels for --signals (default: prefix of each
                                   --signals file)
  --cond-names [<name> ...]        Optional condition labels for --signals (default:
                                   prefix of each --signals file)
  --sample-table <tsv>             Project sample table with sample and condition columns
                                   plus bam/peaks for upstream steps
  --layout {custom,project}        Use fp-tools standard project output layout under
                                   --outdir (default: project when --sample-table is
                                   provided)
  --match-scan-mode {auto,shared,per-sample}
                                   Project match-motifs scan mode. auto uses one shared
                                   motif scan for multi-sample project runs; per-sample
                                   preserves independent sample scans.
  --sample-output-root <directory>
                                   Sample output root; writes one match_motifs folder
                                   under <root>/<sample> for each input signal, typically
                                   <project>/samples
  --peak-header <file>             File containing the header of --peaks separated by
                                   whitespace or newlines (default: peak columns are named
                                   "_additional_<count>")
  --naming <string>                Naming convention for TF output files ('id', 'name',
                                   'name_id', 'id_name') (default: 'name_id')
  --motif-pvalue <float>           Set p-value threshold for motif scanning (default:
                                   1e-4)
  --bound-pvalue <float>           Set p-value threshold for bound/unbound split (default:
                                   0.001)
  --cluster-threshold <float>      Set the clustering threshold. Motifs below this
                                   threshold will be assigned to one cluster (default:
                                   0.5)
  --pseudo <float>                 Pseudocount for calculating log2fcs (default: estimated
                                   from data)
  --skip-excel                     Skip creation of Excel files to speed up large motif
                                   analyses
  --output-peaks <bed>             Gives the possibility to set the output peak set
                                   differently than the input --peaks. This will limit all
                                   analysis to the regions in --output-peaks. NOTE:
                                   --peaks must still be set to the full peak set!
  --norm-off                       Turn off normalization of footprint scores
  --normalization {condition-quantile,sample-quantile,none}
                                   Signal normalization mode (default: none; --norm-off
                                   maps to none)
  --aggregate-signals [<bigwig> ...]
                                   Corrected cut-site bigWigs used for embedded aggregate
                                   profiles
  --plot-aggregate {sig,all,top,off}
                                   Embed aggregate profiles in HTML reports for
                                   significant, all, top-N, or no motifs (default: sig)
  --plot-aggregate-top-n <int>     Maximum number of motifs to aggregate when --plot-
                                   aggregate sig/top or fallback selection is used
                                   (default: 20)
  --plot-aggregate-motifs <motif> [<motif> ...]
                                   Ordered motif IDs, names, or output prefixes to embed
                                   as aggregate profiles; overrides automatic aggregate-
                                   motif selection
  --default-aggregate-plots <int>  Number of aggregate profiles initially displayed in
                                   interactive reports (default: 4; maximum: 12)
  --aggregate-pvalue-threshold <float>
                                   P-value threshold for --plot-aggregate sig (default:
                                   0.05)
  --aggregate-flank <bp>           Flank around motif centers for embedded aggregate
                                   profiles (default: 100)
  --aggregate-normalization {match,none,sample-quantile,size-factor}
                                   Normalization for embedded aggregate profiles (default:
                                   match --normalization)
  --aggregate-site-set {all,bound}
                                   Motif-site BEDs used for embedded aggregate profiles:
                                   all motif hits or sample-specific bound sites (default:
                                   all)
  --motif-outputs {auto,summary,full}
                                   Per-motif output mode. For match-motifs, auto writes
                                   compact caches plus per-motif BED files; summary writes
                                   only main result/report tables and caches; full writes
                                   per-motif BED and overview files synchronously. For
                                   diff-footprints, auto writes full motif outputs only
                                   when aggregate reports need them.
  --report-label <text>            Optional method label shown under the report subtitle
                                   in interactive HTML reports
  --prefix <prefix>                Prefix for overview files in --outdir folder (default:
                                   motif_matches)

Run arguments:
  --outdir <directory>             Output directory to place motif tables, BED files, and
                                   plots in (default: motif_matches_output)
  --cores <int>                    Number of cores to use for computation (default: all
                                   available cores)
  --sample-workers <int>           Number of input samples to process concurrently in
                                   project/sample-output-root mode (default: auto when
                                   --cores is set)
  --split <int>                    Split of multiprocessing jobs (default: 100)
  --debug                          Creates an additional '_debug.pdf'-file with debug
                                   plots
  --verbosity <int>                Level of output logging (0: silent, 1: errors/warnings,
                                   2: info, 3: stats, 4: debug, 5: spam) (default: 3)
```

</div>

## `diff-footprints`

<div class="fp-api-card" markdown="1">

Compare motif-associated footprint scores across conditions or between
user-defined region sets measured in the same sample(s).

**Example command**

```bash
diff-footprints --sample-table project/metadata/samples.tsv --comparison-table project/metadata/comparisons.tsv \
  --genome hg38.fa.gz --peaks project/peaks/merged_peaks_filtered.bed --motif-db jaspar2026_vertebrates --outdir project
```

**Primary inputs**

- `--sample-table` — samples, conditions, footprint tracks, and reusable motif-result folders.
- `--comparison-table` — condition pairs to compare.
- `--genome` — reference genome FASTA.
- `--peaks` — accessible-region BED file.
- `--motif-db` — built-in motif database name.
- `--outdir` — project directory for statistics, figures, and HTML reports.

**Main outputs**

Each comparison is written below
`{project}/comparisons/{comparison}/`, where `{comparison}` is taken from the
comparison table and `{prefix}` defaults to `diff_footprints`:

| Path | Meaning |
| --- | --- |
| `{prefix}_results.txt` | Tab-separated motif-level differential footprint statistics; change direction is `cond1 - cond2`. |
| `{prefix}_results.xlsx` | Excel copy of the result table unless `--skip-excel` is used. |
| `{prefix}_distances.txt` | Motif distances used for clustering related motifs. |
| `{prefix}_{cond1}_{cond2}.html` | Portable interactive report with volcano, motif, and embedded aggregate-profile views. |
| `{prefix}_replicate_report.tsv` | Long-form per-replicate diagnostic data when replicate reporting is active. |
| `{prefix}_replicate_summary.tsv` | Motif-level replicate agreement summary. |
| `{prefix}_replicate_report.png` | Replicate diagnostic figure. |
| `{prefix}_figures.pdf` and `{prefix}_clusters.pdf` | Optional static summaries written with `--static-plots`. |
| `{motif}/beds/{motif}_{condition}_bound.bed` | Motif instances classified as bound for a condition when full motif outputs are required. |

Region-set analyses use the same result/report patterns and add confidence
intervals, motif prevalence, region counts, per-replicate effects, and matching
balance to the result tables.

For a region-set comparison, use `--comparison-axis regions`, provide two or
more BED files with `--regions`, and name them with `--region-labels`. An
optional `--region-strata-column` preserves accessibility or other matching
strata during resampling. One sample uses a stratified label-permutation test;
two or more biological replicates use a paired empirical-Bayes model.
Use `--plot-aggregate-motifs` to choose an ordered aggregate panel without
limiting the motifs tested, and `--default-aggregate-plots` to set its initial
size.

**Complete options**

```text
usage: diff-footprints [-h] [--signals [<bigwig> ...]] [--peaks <bed>] [--genome <fasta>]
                       [--motifs [<motifs> ...]] [--motif-db <name>] [--list-motif-dbs]
                       [--sample-names [<name> ...]] [--cond-names [<name> ...]]
                       [--comparison-axis {conditions,regions}] [--regions [<bed> ...]]
                       [--region-labels [<name> ...]] [--region-strata-column <int>]
                       [--region-permutations <int>] [--region-bootstrap <int>]
                       [--min-regions-per-set <int>] [--random-seed <int>]
                       [--sample-table <tsv>] [--comparison-table <tsv>]
                       [--layout {custom,project}] [--sample-dirs [<directory> ...]]
                       [--project-dir <directory>] [--peak-header <file>]
                       [--naming <string>] [--motif-pvalue <float>]
                       [--bound-pvalue <float>] [--cluster-threshold <float>]
                       [--pseudo <float>] [--time-series] [--time-course] [--skip-excel]
                       [--output-peaks <bed>] [--norm-off]
                       [--normalization {condition-quantile,sample-quantile,none}]
                       [--replicate-report {auto,on,off}] [--replicate-map <tsv>]
                       [--replicate-report-out <tsv>] [--replicate-summary-out <tsv>]
                       [--replicate-figure-out <figure>]
                       [--aggregate-signals [<bigwig> ...]]
                       [--plot-aggregate {sig,all,top,off}] [--plot-aggregate-top-n <int>]
                       [--plot-aggregate-motifs <motif> [<motif> ...]]
                       [--default-aggregate-plots <int>]
                       [--aggregate-pvalue-threshold <float>] [--aggregate-flank <bp>]
                       [--aggregate-normalization {match,none,sample-quantile,size-factor}]
                       [--aggregate-site-set {all,bound}] [--reuse-existing-results]
                       [--motif-outputs {auto,summary,full}] [--static-plots]
                       [--per-motif-plots] [--skew-report] [--report-label <text>]
                       [--outdir <directory>] [--prefix <prefix>] [--cores <int>]
                       [--split <int>] [--debug] [--verbosity <int>]

__________________________________________________________________________________________

                                 fp-tools diff-footprints
__________________________________________________________________________________________

diff-footprints takes motifs, footprint signals, and genome sequence as input to compare
motif-associated footprint evidence across biological conditions or user-defined region
sets. Region-set mode gives every region equal weight and supports matching strata and
paired biological replicates.

Usage:
diff-footprints --signals <bigwig1> (<bigwig2> (...)) --genome <genome.fasta> --peaks
<peaks.bed> [--motif-db jaspar2026_vertebrates | --motifs <motifs.txt>]
diff-footprints --comparison-axis regions --signals <bigwig1> [<replicate2> ...] --regions
<set1.bed> <set2.bed> --genome <genome.fasta> [--region-labels <set1> <set2>]

Output files:
- <outdir>/<prefix>_results.{txt,xlsx}
- <outdir>/<prefix>_distances.txt
- <outdir>/<prefix>_<condition1>_<condition2>.html
- optional <outdir>/<prefix>_figures.pdf with --static-plots
- optional <outdir>/<prefix>_clusters.pdf with --static-plots
- optional <outdir>/<TF>/plots/<TF>_log2fcs.pdf with --per-motif-plots
- optional <outdir>/<TF>/<TF>_overview.{txt,xlsx} with --motif-outputs full
- <outdir>/<TF>/beds/<TF>_<condition>_bound.bed (per motif-condition pair)
- <outdir>/<TF>/beds/<TF>_<condition>_unbound.bed (per motif-condition pair)

------------------------------------------------------------------------------------------

Required arguments:
  --signals [<bigwig> ...]         Signal per condition (.bigwig format)
  --peaks <bed>                    Peaks.bed containing open chromatin regions across all
                                   conditions (not used with --comparison-axis regions)
  --genome <fasta>                 Genome .fasta file

Optional arguments:
  --motifs [<motifs> ...]          Motif file(s) in pfm/jaspar/meme/transfac format; if
                                   omitted, the built-in JASPAR 2026 vertebrates set is
                                   used
  --motif-db <name>                Built-in motif database to use or add to --motifs
                                   (default when --motifs is omitted:
                                   jaspar2026_vertebrates)
  --list-motif-dbs                 List available built-in motif databases and exit
  --sample-names [<name> ...]      Sample labels for --signals; distinct from --cond-names
                                   and used for per-sample score columns, replicate
                                   reports, and aggregate profiles (default: prefix of
                                   each --signals file)
  --cond-names [<name> ...]        Condition labels for --signals; repeat names to define
                                   biological replicates (default: prefix of each
                                   --signals file)
  --comparison-axis {conditions,regions}
                                   Compare biological conditions or region sets measured
                                   in the same sample(s) (default: conditions)
  --regions [<bed> ...]            Two or more non-overlapping BED files for --comparison-
                                   axis regions
  --region-labels [<name> ...]     Labels for --regions (default: BED filename stems)
  --region-strata-column <int>     Optional 1-based BED column containing matching strata
  --region-permutations <int>      Within-stratum label permutations for a one-sample
                                   region comparison (default: 100000)
  --region-bootstrap <int>         Within-set, within-stratum bootstrap samples for
                                   confidence intervals (default: 1000)
  --min-regions-per-set <int>      Minimum motif-containing regions required in each set
                                   (default: 10)
  --random-seed <int>              Random seed for region-set resampling (default: 1)
  --sample-table <tsv>             Project sample table with sample and condition columns
                                   plus bam/peaks for upstream steps
  --comparison-table <tsv>         Project comparison table with comparison, cond1, and
                                   cond2 columns
  --layout {custom,project}        Use fp-tools standard project output layout under
                                   --outdir (default: project when --sample-table is
                                   provided)
  --sample-dirs [<directory> ...]  Sample output folders containing match_motifs/ outputs
                                   to reuse for differential analysis
  --project-dir <directory>        Parent folder containing sample output folders for
                                   folder-based differential analysis, typically
                                   <project>/samples
  --peak-header <file>             File containing the header of --peaks separated by
                                   whitespace or newlines (default: peak columns are named
                                   "_additional_<count>")
  --naming <string>                Naming convention for TF output files ('id', 'name',
                                   'name_id', 'id_name') (default: 'name_id')
  --motif-pvalue <float>           Set p-value threshold for motif scanning (default:
                                   1e-4)
  --bound-pvalue <float>           Set p-value threshold for bound/unbound split (default:
                                   0.001)
  --cluster-threshold <float>      Set the clustering threshold. Motifs below this
                                   threshold will be assigned to one cluster (default:
                                   0.5)
  --pseudo <float>                 Pseudocount for calculating log2fcs (default: estimated
                                   from data)
  --time-series                    Will only compare signals1<->signals2<->signals3 (...)
                                   in order of input, and skip all-against-all comparison.
  --time-course                    Alias for --time-series; compare adjacent ordered
                                   conditions only.
  --skip-excel                     Skip creation of Excel files to speed up large motif
                                   analyses
  --output-peaks <bed>             Gives the possibility to set the output peak set
                                   differently than the input --peaks. This will limit all
                                   analysis to the regions in --output-peaks. NOTE:
                                   --peaks must still be set to the full peak set!
  --norm-off                       Turn off normalization of footprint scores across
                                   conditions
  --normalization {condition-quantile,sample-quantile,none}
                                   Cross-sample normalization mode (default: none; --norm-
                                   off maps to none)
  --replicate-report {auto,on,off}
                                   Write replicate-aware differential-footprint
                                   diagnostics (default: auto for repeated condition names
                                   or --replicate-map)
  --replicate-map <tsv>            Optional TSV with condition/replicate or
                                   condition/n_replicates columns
  --replicate-report-out <tsv>     Output long-form replicate diagnostic TSV (default:
                                   <outdir>/<prefix>_replicate_report.tsv)
  --replicate-summary-out <tsv>    Output replicate diagnostic summary TSV (default:
                                   <outdir>/<prefix>_replicate_summary.tsv)
  --replicate-figure-out <figure>  Output replicate diagnostic figure (default:
                                   <outdir>/<prefix>_replicate_report.png)
  --aggregate-signals [<bigwig> ...]
                                   Corrected cut-site bigWigs used for embedded aggregate
                                   profiles
  --plot-aggregate {sig,all,top,off}
                                   Embed aggregate profiles in HTML reports for
                                   significant, all, top-N, or no motifs (default: sig)
  --plot-aggregate-top-n <int>     Maximum number of motifs to aggregate when --plot-
                                   aggregate sig/top or fallback selection is used
                                   (default: 20)
  --plot-aggregate-motifs <motif> [<motif> ...]
                                   Ordered motif IDs, names, or output prefixes to embed
                                   as aggregate profiles; overrides automatic aggregate-
                                   motif selection
  --default-aggregate-plots <int>  Number of aggregate profiles initially displayed in
                                   interactive reports (default: 4; maximum: 12)
  --aggregate-pvalue-threshold <float>
                                   P-value threshold for --plot-aggregate sig (default:
                                   0.05)
  --aggregate-flank <bp>           Flank around motif centers for embedded aggregate
                                   profiles (default: 100)
  --aggregate-normalization {match,none,sample-quantile,size-factor}
                                   Normalization for embedded aggregate profiles (default:
                                   match --normalization)
  --aggregate-site-set {all,bound}
                                   Motif-site BEDs used for embedded aggregate profiles:
                                   all motif hits or condition-specific bound sites
                                   (default: all)
  --reuse-existing-results         Regenerate final diff-footprints reports from existing
                                   <prefix>_results.txt and per-motif BEDs without
                                   rescanning motifs
  --motif-outputs {auto,summary,full}
                                   Per-motif output mode. For match-motifs, auto writes
                                   compact caches plus per-motif BED files; summary writes
                                   only main result/report tables and caches; full writes
                                   per-motif BED and overview files synchronously. For
                                   diff-footprints, auto writes full motif outputs only
                                   when aggregate reports need them.
  --static-plots                   Also write static volcano and cluster PDF summaries. By
                                   default diff-footprints writes the interactive HTML
                                   report without these PDFs.
  --per-motif-plots                Also write one diagnostic log2 fold-change PDF per
                                   motif. Disabled by default to keep diff-footprints
                                   fast.
  --skew-report                    Also write the optional skew/shift PDF report. Disabled
                                   by default.
  --report-label <text>            Optional method label shown under the report subtitle
                                   in interactive HTML reports
  --prefix <prefix>                Prefix for overview files in --outdir folder (default:
                                   diff_footprints)

Run arguments:
  --outdir <directory>             Output directory to place motif tables, BED files, and
                                   plots in (default: diff_footprints_output)
  --cores <int>                    Number of cores to use for computation (default: all
                                   available cores)
  --split <int>                    Split of multiprocessing jobs (default: 100)
  --debug                          Creates an additional '_debug.pdf'-file with debug
                                   plots
  --verbosity <int>                Level of output logging (0: silent, 1: errors/warnings,
                                   2: info, 3: stats, 4: debug, 5: spam) (default: 3)
```

</div>

## `normalize-bigwig`

<div class="fp-api-card" markdown="1">

Scale corrected cut-site signals using statistics measured over the same
background regions. Use this optional step when samples require an explicitly
shared signal scale before downstream scoring or plotting.

**Example command**

```bash
normalize-bigwig --sample-table project/metadata/samples.tsv --background project/peaks/merged_peaks_filtered.bed \
  --outdir project --method background-scale --stat q95 --target median
```

**Primary inputs**

- `--sample-table` — project samples whose `{sample}_corrected.bw` files are normalized together.
- `--background` — shared BED intervals used to calculate comparable background statistics.
- `--outdir` — project directory represented by `{project}` below.
- `--method` — transformation; `background-scale` multiplies each signal by a shared-target scale factor.
- `--stat` — within-sample background statistic; the example uses the 95th percentile.
- `--target` — across-sample target for the selected statistic; the example uses the median.

**Main outputs**

| Path | Meaning |
| --- | --- |
| `{project}/samples/{sample}/normalize/{sample}_corrected_q95_scaled.bw` | Q95-scaled bias-corrected cut-site signal bigWig for one sample. |
| `{project}/logs/normalize_q95/normalize_bigwig_qc.tsv` | Background statistics, selected statistic, target, and scale factor for every sample. |
| `{project}/logs/normalize_q95/normalize_bigwig_manifest.tsv` | Sample-to-input/output signal mapping for downstream use. |

In custom layout, default outputs use
`{outdir}/{input_stem}.background_scale_{stat}.bw`, plus the two QC tables in
`{outdir}`. `background-zscore` instead writes standardized signal and uses a
method-specific filename suffix.

**Complete options**

```text
usage: normalize-bigwig [-h] [--bigwigs BIGWIGS [BIGWIGS ...]]
                        [--background BACKGROUND] [--outdir OUTDIR]
                        [--sample-names [SAMPLE_NAMES ...]]
                        [--sample-table SAMPLE_TABLE]
                        [--layout {custom,project}]
                        [--sample-output-root SAMPLE_OUTPUT_ROOT]
                        [--method {background-scale,background-zscore,none}]
                        [--stat STAT] [--target {median,mean}]
                        [--chrom-sizes CHROM_SIZES] [--workers WORKERS]

Normalize input signal bigWigs using robust statistics from shared background
BED regions. For corrected cut-site bigWigs, the recommended method is
background-scale.

options:
  -h, --help            show this help message and exit
  --bigwigs BIGWIGS [BIGWIGS ...]
                        Input bigWig files to normalize together.
  --background BACKGROUND
                        Shared background BED used to estimate sample
                        statistics.
  --outdir OUTDIR       Output directory for normalized bigWig QC tables and
                        default outputs.
  --sample-names [SAMPLE_NAMES ...]
                        Sample labels for --bigwigs when using project layout.
  --sample-table SAMPLE_TABLE
                        Project sample table with sample, condition, bam, and
                        peaks columns.
  --layout {custom,project}
                        Use fp-tools standard project output layout under
                        --outdir (default: project when --sample-table is
                        provided).
  --sample-output-root SAMPLE_OUTPUT_ROOT
                        Sample output root; writes each sample under
                        <root>/<sample>/normalize, typically
                        <project>/samples.
  --method {background-scale,background-zscore,none}
                        Normalization method (default: background-scale).
  --stat STAT           Background statistic used by background-scale
                        (default: q90). Use median, iqr, or quantiles such as
                        q90, q95, q97.5, or q99.
  --target {median,mean}
                        Across-sample target statistic for background-scale
                        (default: median).
  --chrom-sizes CHROM_SIZES
                        Optional chromosome sizes file for output
                        validation/header.
  --workers WORKERS     Number of input signal bigWigs to normalize
                        concurrently (default: all available cores, capped by
                        input count).
```

</div>

## `plot-aggregate`

<div class="fp-api-card" markdown="1">

Plot average signal around motif sites or other genomic regions as a static
figure or interactive HTML report.

Multiple user-defined BED files are supported through `--TFBS`. Multiple
`--regions` BED files can restrict or compare distinct regions of interest.

**Example command**

```bash
plot-aggregate --sample-table project/metadata/samples.tsv --motifs SPIB CEBPB --site-set bound --outdir project
```

**Primary inputs**

- `--sample-table` — samples, conditions, and bias-corrected cut-site signal bigWigs used for the aggregate profiles.
- `--motifs` — motif names or identifiers to plot.
- `--site-set` — motif-site set; the example uses bound sites.
- `--outdir` — project directory containing motif results and receiving plots.

**Main outputs**

- `{project}/reports/plot_aggregate.html` — default project-layout interactive aggregate report with motif-centered signal profiles.
- the exact `--output` path — static PDF/PNG/SVG or interactive HTML in custom layout.
- the exact `--output-txt` path — optional per-position aggregate values.
- the exact `--output-aggregated-signals`, `--output-aggregated-scores`, and `--output-aggregated-stats` paths — optional source tables when requested.
- the exact `--output` path in `--motif-grid` mode — multipage motif-by-comparison PDF built from a review bundle.

When both signal types are available, use footprint score bigWigs for motif
statistics and bias-corrected cut-site signal bigWigs for observed aggregate
profiles; label the chosen signal explicitly in figure captions.

```bash
plot-aggregate \
  --input-html project/reports/review_multi_comparisons/index.html \
  --motif-grid \
  --output project/reports/motif_aggregate_grid.pdf
```

**Complete options**

```text
usage: plot-aggregate [-h] [--TFBS [<bed> ...]] [--signals [<bigwig> ...]]
                      [--match-dir [<directory> ...]] [--sample-dirs [<directory> ...]]
                      [--sample-table <tsv>] [--layout {custom,project}]
                      [--manifest <tsv>] [--input-html [<html> ...]]
                      [--regions [<bed> ...]] [--whitelist [<bed> ...]]
                      [--blacklist [<bed> ...]] [--output] [--outdir <directory>]
                      [--output-txt] [--output-csv] [--output_aggregated_signals]
                      [--output_aggregated_scores] [--multiscale-npz <npz>]
                      [--output-multiscale-aggregate] [--title] [--format {auto,pdf,html}]
                      [--flank] [--motifs [<motif> ...]] [--site-set {bound,all,unbound}]
                      [--top-n <int>] [--default-layout {1x1,1x2,2x2,2x3}]
                      [--hide-summary] [--TFBS-labels [...]] [--signal-labels [...]]
                      [--cond-names [<name> ...]] [--region-labels [...]]
                      [--control-label <label>] [--grid <rows>x<cols>] [--share-y]
                      [--normalize]
                      [--normalization {none,condition-quantile,sample-quantile}]
                      [--normalization-comparison-output] [--output_aggregated_stats]
                      [--show-replicate-sd] [--negate] [--smooth <int>] [--log-transform]
                      [--plot-boundaries] [--signal-on-x] [--remove-outliers <float>]
                      [--motif-grid] [--rows-per-page ROWS_PER_PAGE]
                      [--order-htmls [ORDER_HTMLS ...]] [--fill-missing-profiles]
                      [--recompute-missing-profiles] [--repeat-column-labels {none,row}]
                      [--cores <int>] [--verbosity <int>]

__________________________________________________________________________________________

                                 fp-tools plot-aggregate
__________________________________________________________________________________________

Input / output arguments:
  --TFBS [<bed> ...]                    TFBS sites (*required)
  --signals [<bigwig> ...]              Signals in bigwig format (*required)
  --match-dir [<directory> ...]         match-motifs output directory or directories to
                                        use as the motif-site source
  --sample-dirs [<directory> ...]       Alias for --match-dir in HTML mode; sample or
                                        differential output directories containing motif
                                        BEDs
  --sample-table <tsv>                  Project sample table with sample and condition
                                        columns
  --layout {custom,project}             Use fp-tools standard project output layout under
                                        --outdir (default: project when --sample-table is
                                        provided)
  --manifest <tsv>                      TSV with sample, signal, and match_dir/sample_dir
                                        columns for HTML mode
  --input-html [<html> ...]             Existing aggregate or diff-footprints HTML
                                        payload(s) to merge in HTML mode
  --regions [<bed> ...]                 Regions to overlap with TFBS (optional)
  --whitelist [<bed> ...]               Only plot sites overlapping whitelist (optional)
  --blacklist [<bed> ...]               Exclude sites overlapping blacklist (optional)
  --output                              Path to output plot (default: fp-
                                        tools_aggregate.pdf)
  --outdir <directory>                  Project directory used with --layout project
  --output-txt                          Path to output file for aggregates in .txt-format
                                        (default: None)
  --output-csv                          Path to aggregated signal CSV output (default:
                                        None)
  --output_aggregated_signals           Path to CSV file for per-base aggregated signals
                                        (default: None)
  --output_aggregated_scores            Path to CSV file for aggregated footprint-score
                                        table (default: None)
  --multiscale-npz <npz>                Optional call-footprints --output-multiscale-npz
                                        sidecar to render as a scale-by-position aggregate
                                        figure
  --output-multiscale-aggregate         Path for the optional multiscale aggregate figure
                                        (default: <output stem>_multiscale.<output ext>)

Plot arguments:
  --title                               Title of plot (default: "Aggregated signals")
  --format {auto,pdf,html}              Output format for --output. auto uses the output
                                        file extension (default: auto)
  --flank                               Flanking basepairs (+/-) to show in plot (counted
                                        from middle of the TFBS) (default: 60)
  --motifs [<motif> ...]                Motif prefixes, names, or IDs to plot from
                                        --match-dir
  --site-set {bound,all,unbound}        Motif-site BED set to use from --match-dir
                                        (default: bound)
  --top-n <int>                         Number of motifs to plot from --match-dir when
                                        --motifs is omitted (default: 12)
  --default-layout {1x1,1x2,2x2,2x3}    Initial HTML subplot layout (default: 2x2)
  --hide-summary                        Hide the TF site-count summary sidebar in HTML
                                        mode
  --TFBS-labels [ ...]                  Labels used for each TFBS file (default: prefix of
                                        each --TFBS)
  --signal-labels [ ...]                Labels used for each signal file (default: prefix
                                        of each --signals)
  --cond-names [<name> ...]             Condition names for --signals; repeated names are
                                        averaged as replicates
  --region-labels [ ...]                Labels used for each regions file (default: prefix
                                        of each --regions)
  --control-label <label>               Overlay each non-control signal against this
                                        control signal label (must match one of --signal-
                                        labels)
  --grid <rows>x<cols>                  Explicit grid layout for subplots, e.g. 2x5 or
                                        3x4. Panels fill in order of the input signal
                                        files.
  --share-y                             Share y-axis range across plots
                                        (none/signals/sites/both). Use "--share-y signals"
                                        if bigwig signals have similar ranges. Use "--
                                        share_y sites" if sites per bigwig are comparable,
                                        but bigwigs themselves aren't comparable (default:
                                        none)
  --normalize                           Normalize the aggregate signal(s) to be between
                                        0-1 (default: the true range of values is shown)
  --normalization {none,condition-quantile,sample-quantile}
                                        diff-footprints-compatible quantile normalization
                                        before aggregate plotting (default: none)
  --normalization-comparison-output     Optional paired raw-vs-normalized aggregate figure
  --output_aggregated_stats             Path to CSV file for aggregate mean/SD/stat
                                        summaries (default: None)
  --show-replicate-sd                   Draw replicate SD ribbons when --cond-names
                                        contains repeated condition names
  --negate                              Negate overlap with regions
  --smooth <int>                        Smooth output signal by taking the mean of
                                        <smooth> bp windows (default: 1 (no smooth)
  --log-transform                       Log transform the signals before aggregation
  --plot-boundaries                     Plot TFBS boundaries (Note: estimated from first
                                        region in each --TFBS)
  --signal-on-x                         Show signals on x-axis and TFBSs on y-axis
                                        (default: signal is on y-axis)
  --remove-outliers <float>             Value between 0-1 indicating the percentile of
                                        regions to include, e.g. 0.99 to remove the sites
                                        with 1% highest values (default: 1)
  --motif-grid                          Create a multi-page motif-by-comparison PDF from
                                        one review-multi-comparisons report
  --rows-per-page ROWS_PER_PAGE         Motif rows per page in --motif-grid mode (default:
                                        16)
  --order-htmls [ORDER_HTMLS ...]       Optional review reports used to define one shared
                                        motif order in --motif-grid mode
  --fill-missing-profiles               Fill missing motif profiles from profiles embedded
                                        elsewhere in the review report
  --recompute-missing-profiles          Recompute missing motif profiles from project
                                        bigWigs and motif BEDs
  --repeat-column-labels {none,row}     Repeat comparison labels in every motif row
                                        (default: none)

Run arguments:
  --cores <int>                         Worker processes for recomputing missing motif
                                        profiles
  --verbosity <int>                     Level of output logging (0: silent, 1:
                                        errors/warnings, 2: info, 3: stats, 4: debug, 5:
                                        spam) (default: 3)
```

</div>

## `review-multi-comparisons`

<div class="fp-api-card" markdown="1">

Combine differential-footprint reports as a scalable browser bundle or one
self-contained HTML report.

**Example command**

```bash
review-multi-comparisons --inputs project/comparisons --output-dir project/reports/review_multi_comparisons \
  --default-comparison "HNF4A + FOXA2" "No HNF4A/FOXA2" \
  --default-aggregate-motifs MA1494.2 MA0484.3 MA0047.4 MA0148.5 MA0046.3 MA0153.2 MA0102.5 MA0466.4 \
  --default-aggregate-plots 8 --documentation-url https://oncologylab.github.io/fp-tools/
```

**Primary inputs**

- `--inputs` — report files or directories containing differential reports.
- `--output-dir` — destination for the complete static bundle.
- `--default-comparison` — condition or region pair shown first.
- `--default-aggregate-motifs` — ordered motif panel shown first.
- `--default-aggregate-plots` — initial number of aggregate panels.
- `--documentation-url` — optional link back to the documentation site.

**Main outputs**

`{bundle}` is the `--output-dir`:

| Path | Meaning |
| --- | --- |
| `{bundle}/index.html` | Browser entry point; open or publish this file together with the full bundle. |
| `{bundle}/app.js` and `{bundle}/styles.css` | Local application code and styling. |
| `{bundle}/data/metadata.json` | Comparison index and payload checksums. |
| `{bundle}/data/reports/{comparison}.json.gz` | Compact data for one comparison. |
| `{bundle}/data/profiles/` | Aggregate-profile shards loaded on demand. |
| `{bundle}/data/logos/` | Motif logo assets. |

Project mode defaults to
`{project}/reports/review_multi_comparisons/index.html`. The directory is a
portable unit; copying only `index.html` produces a broken report.

## Standalone output

Use `--output-html` instead of `--output-dir`. Aggregate profiles are optional,
and `--labels` keeps repeated condition pairs distinct. The exact output is the
path passed to `--output-html`; it is one portable HTML file with coordinated
volcano, ranked-motif, logo, and SVG-export views. Aggregate controls appear
only when profiles exist. One **Comparison** list selects the exact input record
in `--labels` order, so repeated condition pairs remain distinct.

```bash
review-multi-comparisons --inputs baseline/report.html dose1/report.html dose2/report.html \
  --labels Baseline "Dose 1" "Dose 2" --output-html review.html
```

**Complete options**

```text
usage: review-multi-comparisons [-h] [--inputs INPUTS [INPUTS ...]]
                                [--labels [LABELS ...]]
                                [--output-dir OUTPUT_DIR | --output-html OUTPUT_HTML]
                                [--outdir OUTDIR] [--layout {custom,project}]
                                [--default-comparison <group1> <group2>]
                                [--default-aggregate-motifs <motif> [<motif> ...]]
                                [--default-aggregate-plots <int>]
                                [--documentation-url <url>]
                                [--fill-missing-aggregate-profiles]
                                [--recompute-missing-aggregate-profiles]
                                [--aggregate-flank AGGREGATE_FLANK]
                                [--cores CORES] [--title TITLE]

Combine diff-footprints reports into a static browser bundle or one self-
contained HTML file.

options:
  -h, --help            show this help message and exit
  --inputs INPUTS [INPUTS ...]
                        diff-footprints HTML files or directories containing
                        diff_footprints_*.html files; directories are searched
                        recursively.
  --labels [LABELS ...]
                        Optional labels, one per resolved input HTML.
  --output-dir OUTPUT_DIR
                        Output directory for index.html, JavaScript, CSS, and
                        compact static data files.
  --output-html OUTPUT_HTML
                        One self-contained HTML report; aggregate profiles are
                        optional.
  --outdir OUTDIR       Project directory used with --layout project.
  --layout {custom,project}
                        Use fp-tools standard project output layout under
                        --outdir (default: project when only --outdir is
                        provided).
  --default-comparison <group1> <group2>
                        Region or condition pair initially shown in the static
                        browser
  --default-aggregate-motifs <motif> [<motif> ...]
                        Ordered motif IDs, names, or output prefixes initially
                        shown
  --default-aggregate-plots <int>
                        Number of aggregate profiles initially shown (default:
                        4; maximum: 12)
  --documentation-url <url>
                        Optional link back to the documentation site
  --fill-missing-aggregate-profiles
                        Fill missing motif aggregate panels from profiles
                        embedded elsewhere in the combined review payload.
  --recompute-missing-aggregate-profiles
                        Recompute still-missing motif aggregate panels from
                        project sample bigWigs and match-motifs BEDs.
  --aggregate-flank AGGREGATE_FLANK
                        Flank used when recomputing missing aggregate
                        profiles, or 'auto' to match the existing report axis
                        (default: auto).
  --cores CORES         Worker processes for --recompute-missing-aggregate-
                        profiles (default: all available cores).
  --title TITLE
```

</div>

## `run-yaml-workflow`

<div class="fp-api-card" markdown="1">

Run one or more fp-tools jobs from a reusable YAML configuration.

**Example command**

```bash
run-yaml-workflow --config examples/gui_configs/diff_footprints_single.yml
```

**Primary inputs**

- `--config` — command-compatible YAML configuration exported by the GUI or written directly.

**Main outputs**

- The exact files documented for each command named in the YAML; YAML does not create a separate analysis format.
- Standard output containing the expanded command lines when `--dry-run` is used.
- `{run_root}/{job_id}/config.yml` and `command.txt` — normalized per-job configuration and exact command.
- `{run_root}/{job_id}/status.json`, `stdout.log`, and `stderr.log` — completion state and captured command output.
- `{run_root}/batch_index.tsv` — one-row-per-job batch status index.

Paths are resolved according to the YAML runner and remain independent of GUI
state. Inspect the dry-run expansion before starting a long workflow.

**Complete options**

```text
usage: run-yaml-workflow [-h] --config CONFIG [--run-root RUN_ROOT]
                         [--only [ONLY ...]] [--dry-run] [--list-jobs]
                         [--fail-fast]

Run fp-tools jobs from a YAML config file.

options:
  -h, --help           show this help message and exit
  --config CONFIG      Path to YAML config.
  --run-root RUN_ROOT  Optional directory for run metadata/logs.
  --only [ONLY ...]    Optional tool filter, e.g. diff-footprints.
  --dry-run            Print expanded commands without running.
  --list-jobs          List expanded jobs and exit.
  --fail-fast          Stop at first failed job.
```

</div>

## `fp-tools-gui`

<div class="fp-api-card" markdown="1">

Launch the browser interface for configuring and running fp-tools commands.

Bulk GUI workflows start from coordinate-sorted BAM/BAI files and matching
peak BED files. The GUI does not perform FASTQ-to-BAM preprocessing.
Missing inputs and unsupported options are reported before a run starts.

The GUI is available through the Python package, the complete container, and
the self-contained desktop downloads on the
[release page](https://github.com/oncologylab/fp-tools/releases).

**Example command**

```bash
fp-tools-gui --host 127.0.0.1 --port 8891 --run-dir project/gui_runs --no-browser
```

**Primary inputs**

- `--host` — interface on which the GUI listens (default: `127.0.0.1`).
- `--port` — fixed browser port.
- `--run-dir` — directory for GUI-managed configurations and runs.
- `--no-browser` — start the server without opening a local browser.

**Main outputs**

- `{run_dir}/{timestamp}_{label}/config.yml` — reusable command-compatible YAML saved for a configured run.
- `{run_dir}/{timestamp}_{label}/status.json`, `launcher_stdout.log`, and `launcher_stderr.log` — launcher state and captured batch-runner output.
- `{run_dir}/{timestamp}_{label}/{job_id}/status.json`, `command.txt`, `stdout.log`, and `stderr.log` — per-job state, exact command, and analysis logs.
- The exact analysis files documented by the selected command; the GUI does not introduce GUI-only scientific outputs.

Files under `{run_dir}` are local run state. A saved YAML remains runnable with
`run-yaml-workflow --config {run_dir}/{timestamp}_{label}/config.yml`.

## Local computer

Open the desktop executable or run `fp-tools-gui`. A browser opens after the
server is ready. If it does not, open the local URL printed in the terminal.

## Remote Linux server

Start fp-tools on the server without exposing a network port:

```bash
fp-tools-gui --no-browser --port 8891
```

On your computer, create an SSH tunnel and keep that terminal open:

```bash
ssh -N -L 8891:127.0.0.1:8891 USER@SERVER
```

Open `http://127.0.0.1:8891`. Binding with `--host 0.0.0.0` is also supported,
but fp-tools does not add authentication; protect direct network access with a
firewall, VPN, or reverse proxy.

**Complete options**

```text
usage: fp-tools-gui [-h] [--host HOST] [--port PORT] [--run-dir RUN_DIR]
                    [--no-browser]

Launch the fp-tools browser interface.

options:
  -h, --help         show this help message and exit
  --host HOST        Bind address (default: 127.0.0.1).
  --port PORT        Optional fixed port (default: first free port from 8891).
  --run-dir RUN_DIR  Directory for GUI-managed runs.
  --no-browser       Do not open a local browser automatically.
```

</div>

## `fp-tools-runtime`

<div class="fp-api-card" markdown="1">

Inspect, install, or repair the private external-tool runtime. Linux provides
raw-read and de novo motif components; macOS and Windows provide the optional
de novo motif component.

**Example command**

```bash
fp-tools-runtime status
```

**Primary inputs**

The `status` action takes no input files.

**Main outputs**

The command reports each runtime component, platform, installation state, and
cache location. `install core` and `install homer` are Linux-only raw-read
components. The MEME Suite component is installed only when requested by de
novo motif discovery.

**Complete options**

```text
usage: fp-tools-runtime [-h] {status,install,repair} ...

Inspect, install, or repair the managed fp-tools runtime.

positional arguments:
  {status,install,repair}
    status              Report managed runtime availability and installation
                        state.
    install             Install a runtime component.
    repair              Repair a runtime component.

options:
  -h, --help            show this help message and exit
```

</div>

## `discover-motifs`

<div class="fp-api-card" markdown="1">

Prepare or run de novo motif discovery from candidate footprint intervals or
an existing FASTA file.

**Example command**

```bash
discover-motifs --candidates project/samples/sample/footprints/sample_candidate_footprints.bed --genome hg38.fa.gz \
  --flank 75 --method streme --known-motif-db jaspar2026_vertebrates --outdir project/de_novo/sample --execute
```

**Primary inputs**

- `--candidates` — candidate-footprint BED intervals.
- `--genome` — reference genome used to extract candidate sequences.
- `--flank` — bases included on each side of a candidate center.
- `--method` — discovery method; the example uses STREME.
- `--known-motif-db` — optional known-motif database for Tomtom matching.
- `--outdir` — directory for candidate FASTA files and discovery results.
- `--execute` — run discovery immediately using the managed MEME Suite runtime.

**Main outputs**

`{outdir}` is the selected discovery directory:

| Path | Meaning |
| --- | --- |
| `{outdir}/candidate_sequences.fa` | Reference sequences extracted around candidate footprint intervals. |
| `{outdir}/run_motif_discovery.sh` | Reproducible MEME/DREME/STREME command plan. |
| `{outdir}/{method}/streme.txt` or the method-equivalent MEME output | De novo motif models when `--execute` is used. |
| `{outdir}/tomtom/tomtom.tsv` | Optional similarity matches to the selected known-motif database. |
| `{outdir}/motif_summary.tsv` and `motif_summary.html` | Summary targets written by the generated plan after discovery and matching complete. |

**Complete options**

```text
usage: discover-motifs [-h] (--fasta FASTA | --candidates CANDIDATES)
                       [--genome GENOME] [--flank FLANK] --outdir OUTDIR
                       [--script SCRIPT] [--method {meme,dreme,streme}]
                       [--known-motifs KNOWN_MOTIFS]
                       [--known-motif-db KNOWN_MOTIF_DB] [--list-motif-dbs]
                       [--extra-args ...] [--execute]
                       [--runtime {auto,managed,system,container}]

Prepare or run a de novo motif discovery command plan.

options:
  -h, --help            show this help message and exit
  --fasta FASTA         Existing candidate FASTA.
  --candidates CANDIDATES
                        Candidate BED from call-footprints --output-bed or
                        another BED-like source.
  --genome GENOME       Genome FASTA, required when --candidates is used.
  --flank FLANK         If >0 with --candidates, export +/- flank bp around
                        each candidate center.
  --outdir OUTDIR       External motif discovery output directory.
  --script SCRIPT       Output shell script path. Defaults to
                        <outdir>/run_motif_discovery.sh.
  --method {meme,dreme,streme}
  --known-motifs KNOWN_MOTIFS
                        Optional known motif database for Tomtom comparison.
  --known-motif-db KNOWN_MOTIF_DB
                        Optional built-in motif database for Tomtom
                        comparison.
  --list-motif-dbs      List available built-in motif databases and exit.
  --extra-args ...      Additional arguments appended to MEME/DREME/STREME.
  --execute             Run the generated script immediately.
  --runtime {auto,managed,system,container}
                        External-tool runtime: auto/managed provisions the
                        pinned fp-tools runtime, system uses PATH, and
                        container uses the complete image (default: auto).
```

</div>

## `summarize-motifs`

<div class="fp-api-card" markdown="1">

Summarize MEME, STREME, DREME, and Tomtom results in a compact report.

**Example command**

```bash
summarize-motifs --meme-txt project/de_novo/sample/streme/streme.txt \
  --tomtom-tsv project/de_novo/sample/tomtom/tomtom.tsv --out-tsv project/de_novo/sample/motif_summary.tsv
```

**Primary inputs**

- `--meme-txt` — MEME-compatible discovery output.
- `--tomtom-tsv` — optional Tomtom known-motif matches.
- `--out-tsv` — compact output table for discovered motifs and matches.

**Main outputs**

- the exact `--out-tsv` path — tab-separated discovered motif IDs, consensus sequences, significance values, and known-database matches when available.
- the exact `--out-html` path — optional portable HTML table containing the same summary and motif logos when available.

The command does not rename the requested output prefix; in the example the
primary file is `project/de_novo/sample/motif_summary.tsv`.

**Complete options**

```text
usage: summarize-motifs [-h] [--meme-txt MEME_TXT] [--tomtom-tsv TOMTOM_TSV]
                        --out-tsv OUT_TSV [--out-html OUT_HTML]
                        [--title TITLE]

Summarize MEME/Tomtom outputs into TSV and HTML reports.

options:
  -h, --help            show this help message and exit
  --meme-txt MEME_TXT   MEME text output, usually meme.txt.
  --tomtom-tsv TOMTOM_TSV
                        Tomtom TSV output, usually tomtom.tsv.
  --out-tsv OUT_TSV     Output motif summary TSV.
  --out-html OUT_HTML   Optional output HTML report.
  --title TITLE
```

</div>

## `pseudobulk-fragments`

<div class="fp-api-card" markdown="1">

Group single-cell ATAC fragments by a cell-annotation column to create
pseudobulk inputs.

**Example command**

```bash
pseudobulk-fragments --fragments pbmc_fragments.tsv.gz --annotations cell_annotations.tsv --group-by cell_type \
  --genome-sizes hg38.chrom.sizes --write-cutsite-bigwigs --outdir project/pseudobulk/fragments
```

**Primary inputs**

- `--fragments` — single-cell fragment TSV or TSV.GZ file.
- `--annotations` — barcode-level cell annotation table.
- `--group-by` — annotation column used to define pseudobulk groups.
- `--genome-sizes` — chromosome sizes used to write signal tracks.
- `--write-cutsite-bigwigs` — write cut-site bigWigs for retained groups.
- `--outdir` — directory for grouped fragments, tracks, and QC outputs.

**Main outputs**

For each sanitized `{group}` under `{outdir}`:

| Path | Meaning |
| --- | --- |
| `{group}.fragments.tsv` or `{group}.fragments.tsv.gz` | Fragments assigned to the group; compressed/indexed form is controlled by the command options. |
| `{group}.fragments.tsv.gz.tbi` | Optional Tabix index for random genomic access. |
| `{group}.cutsites.cpm.bw` | Optional CPM-normalized cut-site signal bigWig written by `--write-cutsite-bigwigs`. |
| `{group}.pseudo_pairs.sorted.bam` and `.bai` | Optional pseudo-paired alignment used by `atac-correct` with read shift `0 0`. |
| `pseudobulk_manifest.tsv` | Per-group paths, cell/fragment counts, and filter status. |
| `fp_tools_manifest.yml` | Machine-readable run settings and retained groups. |
| `pseudobulk_downstream_commands.sh` | Optional generated downstream command examples. |

**Complete options**

```text
usage: pseudobulk-fragments [-h] --fragments FRAGMENTS --annotations
                            ANNOTATIONS --group-by GROUP_BY
                            [--barcode-column BARCODE_COLUMN]
                            [--no-strip-barcode-suffix]
                            [--include-chroms INCLUDE_CHROMS]
                            [--exclude-chroms EXCLUDE_CHROMS]
                            [--min-cells MIN_CELLS]
                            [--min-fragments MIN_FRAGMENTS] --outdir OUTDIR
                            [--compress-output] [--index-output]
                            [--write-cutsite-bigwigs] [--write-pseudo-bams]
                            [--no-cpm-normalize] [--write-downstream-commands]
                            [--genome-sizes GENOME_SIZES] [--cores CORES]

Group single-cell ATAC fragments into pseudobulk fragment files.

options:
  -h, --help            show this help message and exit
  --fragments FRAGMENTS
                        10x-style fragments TSV/TSV.GZ with barcode in column
                        4.
  --annotations ANNOTATIONS
                        Cell annotation TSV or CSV.
  --group-by GROUP_BY   Comma-separated annotation columns to group by, e.g.
                        donor,cell_type.
  --barcode-column BARCODE_COLUMN
                        Annotation barcode column (default: barcode).
  --no-strip-barcode-suffix
                        Require exact barcode matches instead of matching
                        AAAC-1 to AAAC.
  --include-chroms INCLUDE_CHROMS
                        Comma-separated chromosomes to keep, e.g.
                        chr1,chr2,chrX.
  --exclude-chroms EXCLUDE_CHROMS
                        Comma-separated chromosomes to skip, e.g. chrM,chrY.
  --min-cells MIN_CELLS
                        Minimum cells for passes_filters (default: 1).
  --min-fragments MIN_FRAGMENTS
                        Minimum fragments for passes_filters (default: 1).
  --outdir OUTDIR       Output directory.
  --compress-output     Write grouped fragments as .tsv.gz files.
  --index-output        BGZF-compress and tabix-index grouped fragments for
                        random access.
  --write-cutsite-bigwigs
                        Write one sparse cut-site bigWig per kept pseudobulk
                        group.
  --write-pseudo-bams   Write sorted pseudo-paired BAMs for kept groups; use
                        atac-correct --read_shift 0 0 on these BAMs.
  --no-cpm-normalize    Write raw cut counts instead of CPM-normalized bigWig
                        values.
  --write-downstream-commands
                        Write a shell script for BED/BAM/bigWig generation
                        from kept pseudobulk groups.
  --genome-sizes GENOME_SIZES
                        Two-column chromosome sizes file used by generated
                        bedtools/UCSC commands and cut-site bigWigs.
  --cores CORES         Cores for compression, bigWig writing, and generated
                        samtools commands (default: all available cores).
```

</div>

## `find-signature-fp`

<div class="fp-api-card" markdown="1">

Calculate and plot per-cell footprint signatures from completed pseudobulk or
motif analyses.

**Example command**

```bash
find-signature-fp --annotations cell_annotations.tsv --fragments pbmc_fragments.tsv.gz --h5ad pbmc_embedding.h5ad \
  --tf-site-dir marker_motif_sites --all-motif-results project/pseudobulk/pseudobulk_diff_footprints_results.txt \
  --outdir project/pseudobulk/signature_fp
```

**Primary inputs**

- `--annotations` — barcode-level cell annotation table.
- `--fragments` — indexed single-cell fragment file.
- `--h5ad` — single-cell object containing the spectral or UMAP embedding.
- `--tf-site-dir` — motif-site directories from the footprint analysis.
- `--all-motif-results` — completed motif-level differential result table.
- `--outdir` — directory for per-cell scores, heatmaps, and UMAP figures.

**Main outputs**

Under `{outdir}` the default names include:

| Path | Meaning |
| --- | --- |
| `knn_footprint_signature_scores.tsv` | Per-cell KNN-smoothed footprint protection scores for selected TFs. |
| `knn_footprint_orientation_summary.tsv` | Direction/orientation checks used to make marker scores comparable. |
| `chromvar_like_motif_activity_scores.tsv` | Companion accessibility-derived motif activity scores. |
| `knn_footprint_signature_umap.svg` and `.pdf` | Per-marker footprint-signature UMAP panels. |
| `per_cell_footprint_signature_heatmap.svg` and `.pdf` | Selected-marker per-cell heatmap. |
| `single_cell_footprinting_summary.svg` and `.pdf` | Combined heatmap and representative UMAP summary. |
| `all_motif_per_cell_footprint_signature_heatmap.tsv` | Optional all-motif score matrix and metadata when all-motif inputs are supplied. |

Additional top-motif and all-TF review files use their requested output prefix.

**Complete options**

```text
usage: find-signature-fp [-h] --annotations ANNOTATIONS --fragments FRAGMENTS
                         --h5ad H5AD [--tf-site-dir TF_SITE_DIR] --outdir
                         OUTDIR [--markers MARKERS]
                         [--max-sites-per-tf MAX_SITES_PER_TF] [--knn KNN]
                         [--flank FLANK]
                         [--center-half-width CENTER_HALF_WIDTH]
                         [--flank-inner FLANK_INNER]
                         [--flank-outer FLANK_OUTER] [--bin-size BIN_SIZE]
                         [--marker-groups MARKER_GROUPS]
                         [--all-motif-diff-dir ALL_MOTIF_DIFF_DIR]
                         [--all-motif-results ALL_MOTIF_RESULTS]
                         [--all-motif-score-table ALL_MOTIF_SCORE_TABLE]
                         [--marker-score-table MARKER_SCORE_TABLE]
                         [--all-motif-batch-size ALL_MOTIF_BATCH_SIZE]
                         [--max-sites-per-motif MAX_SITES_PER_MOTIF]
                         [--max-motifs MAX_MOTIFS]
                         [--top-motif-signatures-per-cell-type TOP_MOTIF_SIGNATURES_PER_CELL_TYPE]
                         [--top-motif-min-specificity TOP_MOTIF_MIN_SPECIFICITY]
                         [--summary-output-prefix SUMMARY_OUTPUT_PREFIX]
                         [--all-tf-review-prefix ALL_TF_REVIEW_PREFIX]
                         [--all-tf-review-panels-per-page ALL_TF_REVIEW_PANELS_PER_PAGE]
                         [--skip-all-tf-review-pdfs]
                         [--no-create-fragment-index]

Generate per-cell footprint-signature heatmaps and UMAP reports.

options:
  -h, --help            show this help message and exit
  --annotations ANNOTATIONS
                        Cell annotation TSV/CSV with barcode, cell type, and
                        UMAP columns.
  --fragments FRAGMENTS
                        10x-style fragments TSV/TSV.GZ used to count cut sites
                        around motif centers.
  --h5ad H5AD           AnnData file containing the single-cell embedding used
                        for KNN smoothing.
  --tf-site-dir TF_SITE_DIR
                        Optional directory containing marker motif-site BED
                        files named by TF. When omitted, marker sites are
                        taken from --all-motif-diff-dir and --all-motif-
                        results.
  --outdir OUTDIR       Output directory for signature score tables, heatmaps,
                        and UMAP reports.
  --markers MARKERS     Comma-separated marker TFs to score and plot (default:
                        STAT6,FOSB,CEBPA,IRF8,RELA,ZNF683,NR4A1,SMAD3).
  --max-sites-per-tf MAX_SITES_PER_TF
                        Maximum marker motif sites per TF for selected-marker
                        UMAP scoring (default: 1500).
  --knn KNN             Number of nearest neighbors used to smooth per-cell
                        cut-site profiles (default: 75).
  --flank FLANK         Motif-centered half-window in bp for fragment counting
                        (default: 100).
  --center-half-width CENTER_HALF_WIDTH
                        Half-width in bp of the protected center window
                        (default: 10).
  --flank-inner FLANK_INNER
                        Inner flank distance from motif center in bp (default:
                        25).
  --flank-outer FLANK_OUTER
                        Outer flank distance from motif center in bp (default:
                        100).
  --bin-size BIN_SIZE   Bin size for the companion chromVAR-like motif
                        activity score (default: 500).
  --marker-groups MARKER_GROUPS
                        Comma-separated TF:cell_type pairs used to orient KNN
                        marker scores for UMAP review.
  --all-motif-diff-dir ALL_MOTIF_DIFF_DIR
                        Optional differential-footprint output directory
                        containing */beds/*_all.bed files for all-motif per-
                        cell heatmap scoring.
  --all-motif-results ALL_MOTIF_RESULTS
                        Differential-footprint results table used to order and
                        annotate all-motif heatmap rows.
  --all-motif-score-table ALL_MOTIF_SCORE_TABLE
                        Existing all-motif per-cell heatmap TSV to redraw
                        all/top heatmaps without rescoring fragments.
  --marker-score-table MARKER_SCORE_TABLE
                        Existing KNN marker score table used to orient
                        selected marker rows in top heatmaps and summary
                        UMAPs.
  --all-motif-batch-size ALL_MOTIF_BATCH_SIZE
                        Number of motif signatures to score per batch for the
                        all-motif heatmap.
  --max-sites-per-motif MAX_SITES_PER_MOTIF
                        Maximum motif instances per motif for all-motif
                        heatmap scoring; use 0 for all sites.
  --max-motifs MAX_MOTIFS
                        Optional all-motif smoke-test limit.
  --top-motif-signatures-per-cell-type TOP_MOTIF_SIGNATURES_PER_CELL_TYPE
                        Top cell-type-specific all-motif signatures to keep
                        per broad cell type (default: 40).
  --top-motif-min-specificity TOP_MOTIF_MIN_SPECIFICITY
                        Minimum dominant-vs-next cell-type mean z-score
                        difference for top all-motif heatmap rows (default:
                        0.5).
  --summary-output-prefix SUMMARY_OUTPUT_PREFIX
                        Output prefix for the combined heatmap and UMAP
                        summary SVG when all-motif heatmap data are available.
  --all-tf-review-prefix ALL_TF_REVIEW_PREFIX
                        Output prefix for three multi-page all-TF signature
                        review PDFs grouped by dominant broad cell type.
  --all-tf-review-panels-per-page ALL_TF_REVIEW_PANELS_PER_PAGE
                        Number of TF signature UMAP panels per all-TF review
                        PDF page (default: 12).
  --skip-all-tf-review-pdfs
                        Do not write the three all-TF signature review PDFs.
  --no-create-fragment-index
                        Do not create a tabix index for the fragment file when
                        it is missing.
```

</div>

## `sc-footprinting`

<div class="fp-api-card" markdown="1">

Run grouping, bias correction, footprint scoring, motif analysis, and per-cell
signature reporting for single-cell ATAC-seq data.

**Example command**

```bash
sc-footprinting --fragments pbmc_fragments.tsv.gz --annotations cell_annotations.tsv --h5ad cell_embedding.h5ad \
  --group-by cell_type --genome-sizes hg38.chrom.sizes --genome hg38.fa.gz --peaks merged_peaks.bed \
  --motif-db jaspar2026_vertebrates --outdir project/pseudobulk
```

**Primary inputs**

- `--fragments` — single-cell fragment file.
- `--annotations` — barcode-level cell annotation table.
- `--h5ad` — AnnData file containing the cell embedding used for KNN smoothing.
- `--group-by` — annotation column used to define pseudobulk groups.
- `--genome-sizes` — chromosome sizes used to write grouped signal tracks.
- `--genome` — reference genome FASTA.
- `--peaks` — accessible-region BED file.
- `--motif-db` — built-in motif database name.
- `--outdir` — directory for pseudobulk tracks, motif results, and reports.

**Main outputs**

`{outdir}` contains a complete staged workflow:

| Path | Meaning |
| --- | --- |
| `pseudobulk/{group}.fragments.tsv.gz` and `.tbi` | Indexed fragments for each retained cell group. |
| `pseudobulk/{group}.cutsites.cpm.bw` | Group cut-site signal bigWig. |
| `pseudobulk/{group}.pseudo_pairs.sorted.bam` and `.bai` | Pseudo-paired alignment used for bias correction. |
| `atacorrect/{group}/{group}_corrected.bw` | Bias-corrected cut-site signal per group. |
| `footprints/{group}_footprints.bw` | Footprint score signal per group. |
| `diff_footprints/pseudobulk_diff_footprints_results.txt` | Optional motif-level group comparison results. |
| `plots/single_cell_footprinting/` | Per-cell score tables, heatmaps, and UMAP figures from `find-signature-fp`. |
| `pseudobulk_footprint_manifest.tsv` | Group paths and workflow completion state. |
| `pseudobulk_footprint_commands.sh` | Exact generated commands for reproducibility. |
| `logs/{stage}.stdout.log` and `{stage}.stderr.log` | Captured output for each stage. |

**Complete options**

```text
usage: sc-footprinting [-h] --fragments FRAGMENTS --annotations ANNOTATIONS
                       --group-by GROUP_BY --outdir OUTDIR
                       [--genome-sizes GENOME_SIZES] --genome GENOME --peaks
                       PEAKS [--blacklist BLACKLIST]
                       [--barcode-column BARCODE_COLUMN]
                       [--no-strip-barcode-suffix]
                       [--include-chroms INCLUDE_CHROMS]
                       [--exclude-chroms EXCLUDE_CHROMS] [--groups GROUPS]
                       [--min-cells MIN_CELLS] [--min-fragments MIN_FRAGMENTS]
                       [--no-cpm-normalize] [--top-n TOP_N]
                       [--read-shift FWD REV] [--motifs [MOTIFS ...]]
                       [--motif-db MOTIF_DB] [--list-motif-dbs]
                       [--peak-header PEAK_HEADER] [--diff-prefix DIFF_PREFIX]
                       [--diff-normalization {condition-quantile,sample-quantile,none}]
                       [--diff-plot-aggregate {sig,all,top,off}]
                       [--skip-excel | --no-skip-excel]
                       [--tf-site-dir TF_SITE_DIR]
                       [--site-summary SITE_SUMMARY] [--tfs TFS]
                       [--plot-flank PLOT_FLANK] [--plot-script PLOT_SCRIPT]
                       --h5ad SINGLE_CELL_SIGNATURE_H5AD
                       [--single-cell-signature-outdir SINGLE_CELL_SIGNATURE_OUTDIR]
                       [--single-cell-signature-markers SINGLE_CELL_SIGNATURE_MARKERS]
                       [--single-cell-signature-fig-prefix SINGLE_CELL_SIGNATURE_FIG_PREFIX]
                       [--single-cell-signature-all-motif-score-table SINGLE_CELL_SIGNATURE_ALL_MOTIF_SCORE_TABLE]
                       [--single-cell-signature-marker-score-table SINGLE_CELL_SIGNATURE_MARKER_SCORE_TABLE]
                       [--single-cell-signature-top-per-cell-type SINGLE_CELL_SIGNATURE_TOP_PER_CELL_TYPE]
                       [--single-cell-signature-top-min-specificity SINGLE_CELL_SIGNATURE_TOP_MIN_SPECIFICITY]
                       [--single-cell-signature-knn SINGLE_CELL_SIGNATURE_KNN]
                       [--single-cell-signature-max-sites-per-motif SINGLE_CELL_SIGNATURE_MAX_SITES_PER_MOTIF]
                       [--single-cell-signature-max-motifs SINGLE_CELL_SIGNATURE_MAX_MOTIFS]
                       [--cores CORES] [--resume] [--force] [--dry-run]
                       [--fail-fast]

Run the complete pseudobulk and per-cell footprint workflow from single-cell
fragments.

options:
  -h, --help            show this help message and exit
  --fragments FRAGMENTS
                        10x-style fragments TSV/TSV.GZ with barcode in column
                        4.
  --annotations ANNOTATIONS
                        Cell annotation TSV or CSV.
  --group-by GROUP_BY   Comma-separated annotation columns to group by.
  --outdir OUTDIR       Output directory for the full pseudobulk footprint
                        workflow.
  --genome-sizes GENOME_SIZES
                        Two-column chromosome sizes file used for fragment-
                        derived cut-site bigWigs.
  --genome GENOME       Genome FASTA for atac-correct.
  --peaks PEAKS         Peak BED used for atac-correct and footprint scoring.
  --blacklist BLACKLIST
                        Optional blacklist BED for atac-correct.
  --barcode-column BARCODE_COLUMN
                        Annotation barcode column (default: barcode).
  --no-strip-barcode-suffix
                        Require exact barcode matches instead of matching
                        AAAC-1 to AAAC.
  --include-chroms INCLUDE_CHROMS
                        Comma-separated chromosomes to keep.
  --exclude-chroms EXCLUDE_CHROMS
                        Comma-separated chromosomes to skip.
  --groups GROUPS       Comma-separated pseudobulk groups to process after
                        grouping; default processes all retained groups.
  --min-cells MIN_CELLS
                        Minimum cells for passes_filters (default: 1).
  --min-fragments MIN_FRAGMENTS
                        Minimum fragments/reads for passes_filters (default:
                        1).
  --no-cpm-normalize    Write raw cut counts instead of CPM-normalized cut-
                        site bigWigs for fragment input.
  --top-n TOP_N         Optional top N candidate footprints per group.
  --read-shift FWD REV  Override the atac-correct read shift for fragment cut
                        sites (default: 0 0).
  --motifs [MOTIFS ...]
                        Optional motif file(s); when provided, run motif-aware
                        diff-footprints on pseudobulk footprint tracks.
  --motif-db MOTIF_DB   Built-in motif database for motif matching (default:
                        jaspar2026_vertebrates); can be combined with
                        --motifs.
  --list-motif-dbs      List available built-in motif databases and exit.
  --peak-header PEAK_HEADER
                        Optional peak-header file passed to diff-footprints.
  --diff-prefix DIFF_PREFIX
                        Prefix for optional motif-aware diff-footprints
                        outputs.
  --diff-normalization {condition-quantile,sample-quantile,none}
                        Normalization mode for optional motif-aware diff-
                        footprints outputs (default: none).
  --diff-plot-aggregate {sig,all,top,off}
                        Aggregate plot selection for optional motif-aware
                        diff-footprints HTML/PDF outputs.
  --skip-excel, --no-skip-excel
                        Skip Excel files for optional diff-footprints outputs
                        (default: on).
  --tf-site-dir TF_SITE_DIR
                        Optional motif-centered BED directory to plot
                        corrected footprint aggregates.
  --site-summary SITE_SUMMARY
                        Optional motif-centered site summary TSV for plotting.
  --tfs TFS             Comma-separated TFs or 'auto' for plotting (default:
                        auto).
  --plot-flank PLOT_FLANK
                        Flank for optional aggregate plots (default: 100).
  --plot-script PLOT_SCRIPT
                        Plotting script path for optional aggregate plots.
  --h5ad SINGLE_CELL_SIGNATURE_H5AD, --single-cell-signature-h5ad SINGLE_CELL_SIGNATURE_H5AD
                        AnnData file containing the cell embedding used for
                        KNN footprint-signature smoothing.
  --single-cell-signature-outdir SINGLE_CELL_SIGNATURE_OUTDIR
                        Output directory for optional per-cell signature
                        reports (default:
                        <outdir>/plots/single_cell_footprinting).
  --single-cell-signature-markers SINGLE_CELL_SIGNATURE_MARKERS
                        Comma-separated marker TFs for optional per-cell
                        signature UMAPs (default:
                        STAT6,FOSB,CEBPA,IRF8,RELA,ZNF683,NR4A1,SMAD3).
  --single-cell-signature-fig-prefix SINGLE_CELL_SIGNATURE_FIG_PREFIX
                        Output prefix for the combined single-cell footprint-
                        signature SVG (default: single_cell_footprinting).
  --single-cell-signature-all-motif-score-table SINGLE_CELL_SIGNATURE_ALL_MOTIF_SCORE_TABLE
                        Existing all-motif per-cell signature TSV; skips
                        rescoring all motif sites for the signature heatmap.
  --single-cell-signature-marker-score-table SINGLE_CELL_SIGNATURE_MARKER_SCORE_TABLE
                        Existing KNN marker score TSV used for marker rows and
                        UMAP plots.
  --single-cell-signature-top-per-cell-type SINGLE_CELL_SIGNATURE_TOP_PER_CELL_TYPE
                        Top all-motif signatures to keep per cell type in the
                        signature heatmap (default: 40).
  --single-cell-signature-top-min-specificity SINGLE_CELL_SIGNATURE_TOP_MIN_SPECIFICITY
                        Minimum dominant-vs-next cell-type z-score difference
                        for top heatmap rows (default: 0.5).
  --single-cell-signature-knn SINGLE_CELL_SIGNATURE_KNN
                        KNN size for optional per-cell footprint-signature
                        smoothing (default: 75).
  --single-cell-signature-max-sites-per-motif SINGLE_CELL_SIGNATURE_MAX_SITES_PER_MOTIF
                        Maximum motif instances per motif for optional all-
                        motif per-cell heatmap scoring; use 0 for all sites
                        (default: 200).
  --single-cell-signature-max-motifs SINGLE_CELL_SIGNATURE_MAX_MOTIFS
                        Optional smoke-test limit for all-motif per-cell
                        heatmap scoring.
  --cores CORES         Cores for grouping, atac-correct, and footprint
                        scoring (default: 1).
  --resume              Skip atac-correct/call-footprints steps whose expected
                        outputs already exist.
  --force               Run atac-correct/call-footprints even if outputs
                        already exist.
  --dry-run             Write manifests and commands without running atac-
                        correct, call-footprints, motif detection, or plots.
  --fail-fast           Stop after the first failed group command.
```

</div>
