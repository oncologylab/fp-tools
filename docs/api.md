# API Reference

The public interface of `fp-tools` is command-first. Each command below includes the required inputs, main outputs, a practical example, and the complete option reference generated from the current `--help` output.

## Command Overview

| Command | Purpose |
| --- | --- |
| [`atac-correct`](#atac-correct) | Bias-correct ATAC-seq cut-site signal. |
| [`call-footprints`](#call-footprints) | Create footprint score tracks from one or more corrected bigWig signals. |
| [`match-motifs`](#match-motifs) | Scan motifs for one or more footprint tracks and write per-sample bound/unbound motif-site calls. |
| [`diff-footprints`](#diff-footprints) | Compare motif-associated footprint scores across conditions or biological replicates. |
| [`normalize-bigwig`](#normalize-bigwig) | Normalize bigWig tracks with a shared background-region scale estimate. |
| [`plot-aggregate`](#plot-aggregate) | Plot aggregate signal around motif sites or region sets as static output or HTML. |
| [`review-multi-comparisons`](#review-multi-comparisons) | Review multiple differential-footprint HTML reports in one page. |
| [`run-workflow`](#run-workflow) | Run a saved YAML workflow configuration. |
| [`fp-tools-gui`](#fp-tools-gui) | Launch the optional browser GUI for command-compatible workflows. |
| [`motif-discovery`](#motif-discovery) | Prepare or run de novo motif discovery from candidate footprint intervals or FASTA. |
| [`motif-summary`](#motif-summary) | Summarize MEME/STREME/DREME and Tomtom motif discovery outputs. |
| [`fp-tools-score-variants`](#fp-tools-score-variants) | Annotate variants with footprint, sequence, candidate-overlap, and optional motif/model score changes. |
| [`pseudobulk-fragments`](#pseudobulk-fragments) | Group single-cell ATAC fragments into pseudobulk fragment files and optional cut-site tracks. |
| [`find-signature-fp`](#find-signature-fp) | Plot per-cell footprint-signature heatmaps and UMAP reports from completed pseudobulk or motif analysis outputs. |
| [`pseudobulk-footprints`](#pseudobulk-footprints) | Run grouping, correction, footprint scoring, motif reports, aggregate plots, and optional signature reporting for single-cell ATAC pseudobulk analyses. |

## Command Manuals

### `atac-correct`

<div class="fp-api-card" markdown="1">

Bias-correct ATAC-seq cut-site signal.

**Input parameters**

- One or more BAM files from ATAC-seq experiments.
- Reference genome FASTA.
- One shared `merged_peaks.bed`, or multiple per-sample peak BEDs that are merged internally; optional blacklist and q95 scaling regions.

**Output**

- Uncorrected, expected, bias, and corrected bigWig tracks.
- QC PDF and run logs in the output directory.
- For multi-BAM runs, one output subfolder per sample. If multiple peak BEDs are supplied, fp-tools merges them and saves `merged_peaks.bed`.
- With a sample table and `--outdir project`, each sample is written under `<project>/samples/<sample>/atac_correct/`, merged peaks are written under `<project>/peaks/`, and downstream commands can reuse the same sample table. Use `--layout custom` for fully manual paths.

**Example commands**

```bash
atac-correct \
  --bams sample.bam \
  --genome hg38.fa.gz \
  --peaks merged_peaks.bed \
  --blacklist hg38.blacklist.bed \
  --outdir project/samples/sample/atac_correct

atac-correct \
  --sample-table project/metadata/samples.tsv \
  --genome hg38.fa.gz \
  --blacklist hg38.blacklist.bed \
  --outdir project
```

**Options**

This option reference is generated from `atac-correct --help` and lists every accepted option for the command.

```text
usage: atac-correct [-h] [--bams [<bam> ...]] [-g <fasta>] [-p [<bed> ...]]
                    [--regions-in <bed>] [--regions-out <bed>] [--blacklist <bed>]
                    [--extend <int>] [--split-strands] [--norm-off]
                    [--write-tracks [<track> ...]] [--track-off [<track> ...]] [--skip-qc]
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
- optional <outdir>/<prefix>_atacorrect.pdf unless --skip-qc is used

------------------------------------------------------------------------------------------

Required arguments:
  --bams [<bam> ...]               One or more .bam files containing reads to be corrected
  -g <fasta>, --genome <fasta>     A .fasta-file containing whole genomic sequence
  -p [<bed> ...], --peaks [<bed> ...]
                                   One shared merged peak BED, or multiple per-sample
                                   peak BEDs to merge internally

Optional arguments:
  --regions-in <bed>               Input regions for estimating bias (default: regions not
                                   in peaks.bed)
  --regions-out <bed>              Output regions (default: peaks.bed)
  --blacklist <bed>                Blacklisted regions in .bed-format (default: None)
  --extend <int>                   Extend output regions with basepairs
                                   upstream/downstream (default: 100)
  --split-strands                  Write out tracks per strand
  --norm-off                       Switches off normalization based on number of reads
  --write-tracks [<track> ...]     bigWig tracks to write (default: corrected; use all for
                                   corrected, uncorrected, bias, and expected)
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

### `call-footprints`

<div class="fp-api-card" markdown="1">

Create footprint score tracks from one or more corrected bigWig signals.

**Input parameters**

- One or more corrected cut-site bigWigs passed with `--signals`.
- BED regions where footprint scores should be computed.

**Output**

- Footprint score bigWig per input signal.
- Optional BED coordinates for candidate footprint peaks used by de novo motif discovery when `--call-candidates`, `--output-bed`, or `--output-beds` is used.
- With `--sample-output-root`, outputs are written under `<root>/<sample>/footprints/`.

**Example commands**

```bash
call-footprints \
  --signals A_corrected.bw B_corrected.bw \
  --sample-names A B \
  --regions merged_peaks.bed \
  --sample-output-root project/samples
```

**Options**

This option reference is generated from `call-footprints --help` and lists every accepted option for the command.

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
                       [--footprint-kernel {fast,legacy}] [--window <int>]
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
  --footprint-kernel {fast,legacy}      Footprint scoring kernel (default: fast; use
                                        legacy for exact historical floating-point
                                        behavior)

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

### `match-motifs`

<div class="fp-api-card" markdown="1">

Scan motifs for one or more footprint tracks and write per-sample motif binding summaries plus compact reuse caches.

**Input parameters**

- One or more footprint score bigWigs passed with `--signals`.
- Genome FASTA and peak BED.
- A built-in motif database through `--motif-db` or custom motif files through `--motifs`.
- Optional `--sample-names` labels for per-sample columns.

**Output**

- Summary tables with motif binding statistics for each sample.
- Motif-site and background-score caches used by `diff-footprints --sample-dirs`; cached comparisons may create internal per-motif shard caches on first reuse.
- Per-motif folders containing `<motif>_all.bed`, `<motif>_<sample>_bound.bed`, and `<motif>_<sample>_unbound.bed` by default. In project mode these BED folders are materialized in the background after report-ready outputs. Use `--motif-outputs summary` to skip permanent BED folders.
- With project/sample-table input or `--sample-output-root`, fp-tools uses a shared motif scan by default for multi-sample runs and writes one normal sample folder under `<root>/<sample>/match_motifs/`.

**Example commands**

```bash
match-motifs \
  --signals A_footprints.bw B_footprints.bw \
  --sample-names A B \
  --genome hg38.fa.gz \
  --peaks merged_peaks.bed \
  --motif-db jaspar2026_vertebrates \
  --sample-output-root project/samples
```

**Options**

This option reference is generated from `match-motifs --help` and lists every accepted option for the command.

```text
usage: match-motifs [-h] [--signals [<bigwig> ...]] [--peaks <bed>] [--genome <fasta>]
                    [--motifs [<motifs> ...]] [--motif-db <name>] [--list-motif-dbs]
                    [--sample-names [<name> ...]] [--cond-names [<name> ...]]
                    [--sample-table <tsv>] [--layout {custom,project}]
                    [--sample-output-root <directory>] [--peak-header <file>]
                    [--naming <string>] [--motif-pvalue <float>] [--bound-pvalue <float>]
                    [--cluster-threshold <float>] [--pseudo <float>] [--skip-excel]
                    [--output-peaks <bed>] [--norm-off]
                    [--normalization {condition-quantile,sample-quantile,none}]
                    [--aggregate-signals [<bigwig> ...]]
                    [--plot-aggregate {sig,all,top,off}] [--plot-aggregate-top-n <int>]
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
- optional <outdir>/<TF>/<TF>_overview.{txt,xlsx} and <outdir>/<TF>/beds/*.bed with
--motif-outputs full

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
  --plot-aggregate-top-n <int>     Maximum number of motifs to aggregate when --plot-aggregate sig/top
                                   or fallback selection is used (default: 20)
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
                                   diff-footprints, auto writes full motif outputs only when
                                   aggregate reports need them.
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

### `diff-footprints`

<div class="fp-api-card" markdown="1">

Compare motif-associated footprint scores across conditions or biological replicates.

**Input parameters**

- Two or more footprint score bigWigs passed with `--signals`, standardized sample folders from prior `match-motifs` runs passed with `--sample-dirs` / `--project-dir`, or a project-mode sample/comparison table.
- `--cond-names` to define conditions; repeat names for replicates.
- Optional `--sample-names` to override file-derived sample labels without changing conditions.
- Genome FASTA, peak BED, and motif database or motif files.
- Optional corrected cut-site bigWigs in `--aggregate-signals` for embedded aggregate profiles.

**Output**

- Motif-level differential footprint tables.
- Optional per-motif overview tables and bound/unbound BED files when `--motif-outputs full` is used.
- Standalone interactive HTML report with volcano summaries and aggregate profiles.
- Optional static volcano/cluster PDFs with `--static-plots`, optional per-motif diagnostic PDFs with `--per-motif-plots`, and optional skew/shift PDF with `--skew-report`.
- Replicate diagnostic tables when replicate groups are present.

**Example commands**

```bash
diff-footprints \
  --sample-table project/metadata/samples.tsv \
  --comparison-table project/metadata/comparisons.tsv \
  --genome hg38.fa.gz \
  --peaks project/peaks/merged_peaks_filtered.bed \
  --motif-db jaspar2026_vertebrates \
  --outdir project
```

Project-mode comparisons use a generic table with `comparison`, `cond1`, and `cond2` columns. Each comparison reuses cached `match_motifs` folders from the project samples directory, so folder mode skips motif rescanning, motif-site rescoring, and footprint bigWig rereads for the differential table. The first cached comparison can build internal per-motif shard caches; later comparisons reuse those shards.

**Options**

This option reference is generated from `diff-footprints --help` and lists every accepted option for the command.

```text
usage: diff-footprints [-h] [--signals [<bigwig> ...]] [--peaks <bed>] [--genome <fasta>]
                       [--motifs [<motifs> ...]] [--motif-db <name>] [--list-motif-dbs]
                       [--sample-names [<name> ...]] [--cond-names [<name> ...]]
                       [--sample-dirs [<directory> ...]] [--project-dir <directory>]
                       [--peak-header <file>] [--naming <string>] [--motif-pvalue <float>]
                       [--bound-pvalue <float>] [--cluster-threshold <float>]
                       [--pseudo <float>] [--time-series] [--time-course] [--skip-excel]
                       [--output-peaks <bed>] [--norm-off]
                       [--normalization {condition-quantile,sample-quantile,none}]
                       [--replicate-report {auto,on,off}] [--replicate-map <tsv>]
                       [--replicate-report-out <tsv>] [--replicate-summary-out <tsv>]
                       [--replicate-figure-out <figure>]
                       [--aggregate-signals [<bigwig> ...]]
                       [--plot-aggregate {sig,all,top,off}] [--plot-aggregate-top-n <int>]
                       [--aggregate-pvalue-threshold <float>] [--aggregate-flank <bp>]
                       [--aggregate-normalization {match,none,sample-quantile,size-factor}]
                       [--aggregate-site-set {all,bound}] [--reuse-existing-results]
                       [--static-plots] [--per-motif-plots] [--skew-report]
                       [--report-label <text>] [--outdir <directory>] [--prefix <prefix>]
                       [--cores <int>] [--split <int>] [--debug] [--verbosity <int>]

__________________________________________________________________________________________

                                 fp-tools diff-footprints
__________________________________________________________________________________________

diff-footprints takes motifs, footprint signals, and genome sequence as input to infer
motif-associated bound sites and compare footprint evidence across conditions. The method
ranks motifs by signal differences across input conditions and reports motif-level and
site-level results.

Usage:
diff-footprints --signals <bigwig1> (<bigwig2> (...)) --genome <genome.fasta> --peaks
<peaks.bed> [--motif-db jaspar2026_vertebrates | --motifs <motifs.txt>]

Output files:
- <outdir>/<prefix>_results.{txt,xlsx}
- <outdir>/<prefix>_distances.txt
- <outdir>/<prefix>_<condition1>_<condition2>.html
- optional <outdir>/<prefix>_figures.pdf with --static-plots
- optional <outdir>/<prefix>_clusters.pdf with --static-plots
- optional <outdir>/<TF>/plots/<TF>_log2fcs.pdf with --per-motif-plots
- <outdir>/<TF>/<TF>_overview.{txt,xlsx} (per motif)
- <outdir>/<TF>/beds/<TF>_all.bed (per motif)
- <outdir>/<TF>/beds/<TF>_<condition>_bound.bed (per motif-condition pair)
- <outdir>/<TF>/beds/<TF>_<condition>_unbound.bed (per motif-condition pair)

------------------------------------------------------------------------------------------

Required arguments:
  --signals [<bigwig> ...]         Signal per condition (.bigwig format)
  --peaks <bed>                    Peaks.bed containing open chromatin regions across all
                                   conditions
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
  --sample-dirs [<directory> ...]  Sample output folders containing match_motifs/ outputs
                                   to reuse for differential analysis
  --project-dir <directory>        Parent folder containing sample output folders for
                                   folder-based differential analysis
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
  --plot-aggregate-top-n <int>     Maximum number of motifs to aggregate when --plot-aggregate sig/top
                                   or fallback selection is used (default: 20)
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

### `normalize-bigwig`

<div class="fp-api-card" markdown="1">

Normalize bigWig tracks with a shared background-region scale estimate.

**Input parameters**

- A project sample table with `sample`, `condition`, `bam`, and `peaks` columns, or one or more explicit bigWig tracks.
- Shared background BED regions. In project mode this is usually `project/peaks/merged_peaks_filtered.bed`, written by `atac-correct`.
- Output project directory and optional labels.

**Output**

- Per-sample normalized bigWig files under `project/samples/<sample>/normalize/`.
- A summary table of background statistics and scaling factors.
- A manifest of normalized output tracks.

**Example commands**

```bash
normalize-bigwig \
  --sample-table project/metadata/samples.tsv \
  --background project/peaks/merged_peaks_filtered.bed \
  --outdir project \
  --method background-scale \
  --stat q95 \
  --target median
```

Project mode reads corrected tracks from each sample's `atac_correct` folder and writes q95-scaled tracks back into the same sample folder structure for `call-footprints` to reuse.

**Options**

This option reference is generated from `normalize-bigwig --help` and lists every accepted option for the command.

```text
usage: normalize-bigwig [-h] [--bigwigs BIGWIGS [BIGWIGS ...]]
                        [--background BACKGROUND] [--outdir OUTDIR]
                        [--sample-names [SAMPLE_NAMES ...]]
                        [--sample-table SAMPLE_TABLE] [--layout {custom,project}]
                        [--sample-output-root SAMPLE_OUTPUT_ROOT]
                        [--method {background-scale,background-zscore,none}]
                        [--stat STAT] [--target {median,mean}]
                        [--chrom-sizes CHROM_SIZES] [--workers WORKERS]

Normalize bigWig tracks using robust statistics from shared background BED
regions. For corrected cut-site bigWigs, the recommended method is background-
scale.

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
  --workers WORKERS     Number of bigWig tracks to normalize concurrently
                        (default: all available cores, capped by input count).
```

</div>

### `plot-aggregate`

<div class="fp-api-card" markdown="1">

Plot aggregate signal around motif sites or region sets as static output or HTML.

**Input parameters**

- Signal bigWigs passed with `--signals`.
- Either TFBS BED files with `--TFBS` or a `match-motifs` output directory with `--match-dir`.
- Optional motif names/IDs, site set, regions, labels, and normalization mode.

**Output**

- PDF or HTML aggregate plots depending on `--format` or output extension.
- Optional aggregated signal and summary tables.

**Example commands**

```bash
plot-aggregate \
  --sample-table project/metadata/samples.tsv \
  --motifs SPIB CEBPB \
  --site-set bound \
  --outdir project
```

Project mode reads each sample's `match_motifs` folder and corrected cut-site bigWig, then writes the aggregate browser to `project/reports/plot_aggregate.html` unless `--output` is provided.

**Options**

This option reference is generated from `plot-aggregate --help` and lists every accepted option for the command.

```text
usage: plot-aggregate [-h] [--TFBS [<bed> ...]] [--signals [<bigwig> ...]]
                      [--match-dir [<directory> ...]] [--regions [<bed> ...]]
                      [--whitelist [<bed> ...]] [--blacklist [<bed> ...]] [--output]
                      [--output-txt] [--output-csv] [--output_aggregated_signals]
                      [--output_aggregated_scores] [--multiscale-npz <npz>]
                      [--output-multiscale-aggregate] [--title] [--format {auto,pdf,html}]
                      [--flank] [--motifs [<motif> ...]] [--site-set {bound,all,unbound}]
                      [--top-n <int>] [--default-layout {1x1,1x2,2x2,2x3}]
                      [--TFBS-labels [...]] [--signal-labels [...]]
                      [--cond-names [<name> ...]] [--region-labels [...]]
                      [--control-label <label>] [--grid <rows>x<cols>] [--share-y]
                      [--normalize]
                      [--normalization {none,condition-quantile,sample-quantile}]
                      [--normalization-comparison-output] [--output_aggregated_stats]
                      [--show-replicate-sd] [--negate] [--smooth <int>] [--log-transform]
                      [--plot-boundaries] [--signal-on-x] [--remove-outliers <float>]
                      [--verbosity <int>]

__________________________________________________________________________________________

                                 fp-tools plot-aggregate
__________________________________________________________________________________________

Input / output arguments:
  --TFBS [<bed> ...]                    TFBS sites (*required)
  --signals [<bigwig> ...]              Signals in bigwig format (*required)
  --match-dir [<directory> ...]         match-motifs output directory or directories to
                                        use as the motif-site source
  --regions [<bed> ...]                 Regions to overlap with TFBS (optional)
  --whitelist [<bed> ...]               Only plot sites overlapping whitelist (optional)
  --blacklist [<bed> ...]               Exclude sites overlapping blacklist (optional)
  --output                              Path to output plot (default: fp-
                                        tools_aggregate.pdf)
  --output-txt                          Path to output file for aggregates in .txt-format
                                        (default: None)
  --output-csv                          Legacy alias for aggregated signal CSV output
                                        (default: None)
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

Run arguments:
  --verbosity <int>                     Level of output logging (0: silent, 1:
                                        errors/warnings, 2: info, 3: stats, 4: debug, 5:
                                        spam) (default: 3)
```

</div>

### `review-multi-comparisons`

<div class="fp-api-card" markdown="1">

Review multiple diff-footprints HTML reports in one interactive HTML file.

**Input parameters**

- One or more `diff_footprints_*.html` files, or directories containing those files.
- Optional labels for the resolved comparisons.

**Output**

- A standalone HTML report with comparison bar/volcano plot pairs, shared selected motifs, aggregate profiles, and editable SVG export buttons.

**Example commands**

```bash
review-multi-comparisons \
  --outdir project
```

**Options**

This option reference is generated from `review-multi-comparisons --help` and lists every accepted option for the command.

```text
usage: review-multi-comparisons [-h] [--inputs INPUTS [INPUTS ...]]
                                [--labels [LABELS ...]] [--output OUTPUT]
                                [--outdir OUTDIR] [--layout {custom,project}]
                                [--title TITLE]

Review multiple diff-footprints HTML reports in one interactive HTML file.

options:
  -h, --help            show this help message and exit
  --inputs INPUTS [INPUTS ...]
                        diff-footprints HTML files or directories containing
                        diff_footprints_*.html files; directories are searched
                        recursively.
  --labels [LABELS ...]
                        Optional labels, one per resolved input HTML.
  --output OUTPUT       Output standalone review HTML.
  --outdir OUTDIR       Project directory used with --layout project.
  --layout {custom,project}
                        Use fp-tools standard project output layout under
                        --outdir (default: project when only --outdir is
                        provided).
  --title TITLE
```

</div>

### `run-workflow`

<div class="fp-api-card" markdown="1">

Run a saved YAML workflow configuration.

**Input parameters**

- YAML config exported from the GUI or written by hand.
- Optional run root and tool filters.

**Output**

- Runs each configured command or prints the commands with dry-run mode.

**Example commands**

```bash
run-workflow \
  --config examples/gui_configs/diff_footprints.yml
```

**Options**

This option reference is generated from `run-workflow --help` and lists every accepted option for the command.

```text
usage: run-workflow [-h] --config CONFIG [--run-root RUN_ROOT]
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

### `fp-tools-gui`

<div class="fp-api-card" markdown="1">

Launch the optional browser GUI for command-compatible workflows.

**Input parameters**

- Host, port, and optional run directory.

**Output**

- A local or server-hosted Streamlit GUI that writes reusable YAML configs and launches fp-tools commands.

**Example commands**

```bash
fp-tools-gui \
  --host 0.0.0.0 \
  --port 8891
```

**Options**

This option reference is generated from `fp-tools-gui --help` and lists every accepted option for the command.

```text
usage: fp-tools-gui [-h] [--host HOST] [--port PORT] [--run-dir RUN_DIR]

Launch the fp-tools Streamlit GUI.

options:
  -h, --help         show this help message and exit
  --host HOST        Bind address for the GUI server.
  --port PORT        Optional fixed port.
  --run-dir RUN_DIR  Directory for GUI-managed runs.
```

</div>

### `motif-discovery`

<div class="fp-api-card" markdown="1">

Prepare or run de novo motif discovery from candidate footprint intervals or FASTA.

**Input parameters**

- Candidate BED intervals or a FASTA file.
- Reference genome FASTA when extracting candidate-centered sequences.
- Motif discovery method and optional known motif database for comparison.

**Output**

- Candidate FASTA files.
- Runnable MEME/STREME/DREME command plan or executed motif discovery output.
- Optional Tomtom comparison and summary files.

**Example commands**

```bash
motif-discovery \
  --candidates project/samples/sample/footprints/sample_candidate_footprints.bed \
  --genome hg38.fa.gz \
  --flank 75 \
  --method streme \
  --known-motif-db jaspar2026_vertebrates \
  --outdir project/de_novo/sample
```

**Options**

This option reference is generated from `motif-discovery --help` and lists every accepted option for the command.

```text
usage: motif-discovery [-h] (--fasta FASTA | --candidates CANDIDATES)
                       [--genome GENOME] [--flank FLANK] --outdir OUTDIR
                       [--script SCRIPT] [--method {meme,dreme,streme}]
                       [--known-motifs KNOWN_MOTIFS]
                       [--known-motif-db KNOWN_MOTIF_DB] [--list-motif-dbs]
                       [--extra-args ...] [--execute]

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
```

</div>

### `motif-summary`

<div class="fp-api-card" markdown="1">

Summarize MEME/STREME/DREME and Tomtom motif discovery outputs.

**Input parameters**

- Motif discovery text/XML outputs and optional Tomtom TSV.
- Output TSV or HTML paths.

**Output**

- Compact motif summary table and optional HTML report with motif names, consensus strings, and database matches.

**Example commands**

```bash
motif-summary \
  --meme-txt project/de_novo/sample/streme/streme.txt \
  --tomtom-tsv project/de_novo/sample/tomtom/tomtom.tsv \
  --out-tsv project/de_novo/sample/motif_summary.tsv
```

**Options**

This option reference is generated from `motif-summary --help` and lists every accepted option for the command.

```text
usage: motif-summary [-h] [--meme-txt MEME_TXT] [--tomtom-tsv TOMTOM_TSV]
                     --out-tsv OUT_TSV [--out-html OUT_HTML] [--title TITLE]

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

### `fp-tools-score-variants`

<div class="fp-api-card" markdown="1">

Annotate variants with footprint, sequence, candidate-overlap, and optional motif/model score changes.

**Input parameters**

- Variant VCF/BED-like input.
- Reference genome FASTA.
- Optional footprint bigWigs, candidate BEDs, motif files, and trained models.

**Output**

- Variant-level TSV with requested footprint, sequence, motif, or model delta columns.

**Example commands**

```bash
fp-tools-score-variants \
  --variants variants.vcf \
  --genome hg38.fa.gz \
  --out project/variants/variant_scores.tsv
```

**Options**

This option reference is generated from `fp-tools-score-variants --help` and lists every accepted option for the command.

```text
usage: fp-tools-score-variants [-h] --variants VARIANTS --genome GENOME --out
                               OUT [--candidate-scores CANDIDATE_SCORES]
                               [--sequence-flank SEQUENCE_FLANK]
                               [--kmer-size KMER_SIZE] [--motifs [MOTIFS ...]]
                               [--motif-db MOTIF_DB] [--list-motif-dbs]
                               [--motif-flank MOTIF_FLANK]
                               [--tfbs-model TFBS_MODEL]

Annotate variants with genome allele checks and footprint/candidate overlaps.

options:
  -h, --help            show this help message and exit
  --variants VARIANTS   BED-like variants: chrom start end name ref alt.
  --genome GENOME       Genome FASTA, optionally gzipped.
  --out OUT             Output TSV.
  --candidate-scores CANDIDATE_SCORES
                        Optional BED-like scored candidates or footprint
                        intervals.
  --sequence-flank SEQUENCE_FLANK
                        Flanking bases on each side for ref/alt sequence-
                        context delta features.
  --kmer-size KMER_SIZE
                        K-mer size for exact ref/alt disruption features.
  --motifs [MOTIFS ...]
                        Optional JASPAR/MEME motif files for best ref/alt PWM
                        delta scoring.
  --motif-db MOTIF_DB   Optional built-in motif database for best ref/alt PWM
                        delta scoring; can be combined with --motifs.
  --list-motif-dbs      List available built-in motif databases and exit.
  --motif-flank MOTIF_FLANK
                        Flanking bases on each side for motif ref/alt delta
                        scoring.
  --tfbs-model TFBS_MODEL
                        Optional fp-tools tabular TFBS model pickle for
                        ref/alt probability deltas.
```

</div>

### `pseudobulk-fragments`

<div class="fp-api-card" markdown="1">

Group single-cell ATAC fragments into pseudobulk fragment files and optional cut-site tracks.

**Input parameters**

- Fragment TSV/TSV.GZ file.
- Cell annotation table and grouping column.
- Optional genome sizes for cut-site bigWigs.

**Output**

- One pseudobulk fragment file per group.
- Manifest and QC summary tables.
- Optional indexed fragments and CPM-normalized cut-site bigWigs.

**Example commands**

```bash
pseudobulk-fragments \
  --fragments pbmc_fragments.tsv.gz \
  --annotations cell_annotations.tsv \
  --group-by cell_type \
  --genome-sizes hg38.chrom.sizes \
  --write-cutsite-bigwigs \
  --outdir project/pseudobulk/fragments
```

**Options**

This option reference is generated from `pseudobulk-fragments --help` and lists every accepted option for the command.

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
                        atac-correct --bams ... --read_shift 0 0 on these BAMs.
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

### `find-signature-fp`

<div class="fp-api-card" markdown="1">

Plot per-cell footprint-signature heatmaps and UMAP reports from completed pseudobulk or motif analysis outputs.

**Input parameters**

- Cell annotations, fragments, and single-cell embedding H5AD.
- TF site directory and motif result tables from prior analysis.
- Output directory and optional marker groups.

**Output**

- Marker footprint-signature heatmaps.
- UMAP panels colored by cell type and selected motif-associated footprint signatures.
- Source tables for the plotted signatures.

**Example commands**

```bash
find-signature-fp \
  --annotations cell_annotations.tsv \
  --fragments pbmc_fragments.tsv.gz \
  --h5ad pbmc_embedding.h5ad \
  --tf-site-dir marker_motif_sites \
  --all-motif-results project/pseudobulk/pseudobulk_diff_footprints_results.txt \
  --outdir project/pseudobulk/signature_fp
```

**Options**

This option reference is generated from `find-signature-fp --help` and lists every accepted option for the command.

```text
usage: find-signature-fp [-h] --annotations ANNOTATIONS --fragments FRAGMENTS
                         --h5ad H5AD --tf-site-dir TF_SITE_DIR --outdir OUTDIR
                         [--markers MARKERS]
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
                        Directory containing marker motif-site BED files named
                        by TF.
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

### `pseudobulk-footprints`

<div class="fp-api-card" markdown="1">

Run grouping, correction, footprint scoring, motif reports, aggregate plots, and optional signature reporting for single-cell ATAC pseudobulk analyses.

**Input parameters**

- Single-cell fragments or BAM input plus annotation table.
- Grouping column, genome, peaks, and optional motif database.
- Optional TF site directory and H5AD for signature reports.

**Output**

- Pseudobulk fragments and pseudo-BAMs.
- Corrected cut-site and footprint score bigWigs.
- Motif-aware differential reports and aggregate plots.
- Optional single-cell footprint-signature heatmaps and UMAPs.

**Example commands**

```bash
pseudobulk-footprints \
  --fragments pbmc_fragments.tsv.gz \
  --annotations cell_annotations.tsv \
  --group-by cell_type \
  --genome-sizes hg38.chrom.sizes \
  --genome hg38.fa.gz \
  --peaks merged_peaks.bed \
  --motif-db jaspar2026_vertebrates \
  --outdir project/pseudobulk
```

**Options**

This option reference is generated from `pseudobulk-footprints --help` and lists every accepted option for the command.

```text
usage: pseudobulk-footprints [-h] (--fragments FRAGMENTS | --bam BAM)
                             --annotations ANNOTATIONS --group-by GROUP_BY
                             --outdir OUTDIR [--genome-sizes GENOME_SIZES]
                             --genome GENOME --peaks PEAKS
                             [--blacklist BLACKLIST]
                             [--barcode-column BARCODE_COLUMN]
                             [--bam-barcode-tag BAM_BARCODE_TAG]
                             [--no-strip-barcode-suffix]
                             [--include-chroms INCLUDE_CHROMS]
                             [--exclude-chroms EXCLUDE_CHROMS]
                             [--groups GROUPS] [--min-cells MIN_CELLS]
                             [--min-fragments MIN_FRAGMENTS]
                             [--no-cpm-normalize] [--top-n TOP_N]
                             [--read-shift FWD REV] [--motifs [MOTIFS ...]]
                             [--motif-db MOTIF_DB] [--list-motif-dbs]
                             [--peak-header PEAK_HEADER]
                             [--diff-prefix DIFF_PREFIX]
                             [--diff-normalization {condition-quantile,sample-quantile,none}]
                             [--diff-plot-aggregate {sig,all,top,off}]
                             [--skip-excel | --no-skip-excel]
                             [--tf-site-dir TF_SITE_DIR]
                             [--site-summary SITE_SUMMARY] [--tfs TFS]
                             [--plot-flank PLOT_FLANK]
                             [--plot-script PLOT_SCRIPT]
                             [--single-cell-signature-h5ad SINGLE_CELL_SIGNATURE_H5AD]
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

Run a full pseudobulk ATAC footprint workflow from single-cell fragments or
tagged BAMs.

options:
  -h, --help            show this help message and exit
  --fragments FRAGMENTS
                        10x-style fragments TSV/TSV.GZ with barcode in column
                        4.
  --bam BAM             All-cell BAM with cell barcodes stored in a read tag,
                        e.g. CB.
  --annotations ANNOTATIONS
                        Cell annotation TSV or CSV.
  --group-by GROUP_BY   Comma-separated annotation columns to group by.
  --outdir OUTDIR       Output directory for the full pseudobulk footprint
                        workflow.
  --genome-sizes GENOME_SIZES
                        Two-column chromosome sizes file; required for
                        fragment input cut-site bigWigs and pseudo-BAMs.
  --genome GENOME       Genome FASTA for atac-correct.
  --peaks PEAKS         Peak BED used for atac-correct and footprint scoring.
  --blacklist BLACKLIST
                        Optional blacklist BED for atac-correct.
  --barcode-column BARCODE_COLUMN
                        Annotation barcode column (default: barcode).
  --bam-barcode-tag BAM_BARCODE_TAG
                        BAM read tag containing cell barcodes for --bam input
                        (default: CB).
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
  --read-shift FWD REV  Override atac-correct read shift; default is 0 0 for
                        fragment pseudo-BAMs and 4 -5 for tagged BAM input.
  --motifs [MOTIFS ...]
                        Optional motif file(s); when provided, run motif-aware
                        diff-footprints on pseudobulk footprint tracks.
  --motif-db MOTIF_DB   Optional built-in motif database for motif-aware diff-
                        footprints; can be combined with --motifs.
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
  --single-cell-signature-h5ad SINGLE_CELL_SIGNATURE_H5AD
                        Optional h5ad with cell embeddings/counts; with
                        --fragments and --tf-site-dir, write per-cell KNN
                        footprint-signature heatmaps and UMAP reports.
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
