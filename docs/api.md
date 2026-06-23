# API Reference

The public interface of `fp-tools` is command-first. This page lists the full command manuals for the current primary commands.

## Command Overview

| Command | Purpose |
| --- | --- |
| [`atac-correct`](#atac-correct) | Bias-correct ATAC-seq cut-site signal. |
| [`call-footprints`](#call-footprints) | Create footprint score tracks from one or more bigWig signals. |
| [`match-motifs`](#match-motifs) | Scan motifs for one or more footprint tracks and infer bound/unbound motif sites. |
| [`diff-footprints`](#diff-footprints) | Compare motif-associated footprint scores across conditions or replicates. |
| [`normalize-bigwig`](#normalize-bigwig) | Normalize bigWig tracks using shared background regions. |
| [`plot-aggregate`](#plot-aggregate) | Plot aggregate signal around motif sites as PDF/SVG-style output or HTML. |
| [`plot-aggregate-batch`](#plot-aggregate-batch) | Compatibility command for interactive aggregate HTML reports. |
| [`run-workflow`](#run-workflow) | Run a saved YAML workflow config. |
| [`fp-tools-gui`](#fp-tools-gui) | Launch the optional browser GUI. |
| [`motif-discovery`](#motif-discovery) | Prepare or run a de novo motif discovery command plan. |
| [`motif-summary`](#motif-summary) | Summarize MEME/Tomtom motif discovery outputs. |
| [`fp-tools-score-variants`](#fp-tools-score-variants) | Optional utility for annotating variants with footprint, sequence, and motif scores. |
| [`pseudobulk-fragments`](#pseudobulk-fragments) | Group single-cell ATAC fragments into pseudobulk files. |
| [`find-signature-fp`](#find-signature-fp) | Plot per-cell footprint-signature heatmaps and UMAP reports. |
| [`pseudobulk-footprints`](#pseudobulk-footprints) | Run grouping, ATAC correction, footprint scoring, reports, aggregate plots, and optional signature reporting. |

## Command Manuals

### `atac-correct`

Bias-correct ATAC-seq cut-site signal.

```text
usage: atac-correct [-h] [-b <bam>] [-g <fasta>] [-p <bed>] [--regions-in <bed>]
                    [--regions-out <bed>] [--blacklist <bed>] [--extend <int>]
                    [--split-strands] [--norm-off] [--track-off [<track> ...]]
                    [--scale-corrected {auto,none,q95}] [--scale-background <bed>]
                    [--scale-corrected-bigwigs [<bigwig> ...]]
                    [--scale-target {median,mean}] [--scale-chrom-sizes <chrom.sizes>]
                    [--drop-chroms [<chrom> ...]] [--k_flank <int>]
                    [--read_shift <int> <int>] [--bg_shift <int>] [--window <int>]
                    [--score_mat <mat>] [--bias-pkl <obj>] [--prefix <prefix>]
                    [--outdir <directory>] [--cores <int>] [--split <int>]
                    [--verbosity <int>]

__________________________________________________________________________________________

                                  fp-tools atac-correct
__________________________________________________________________________________________

atac-correct corrects ATAC-seq cutsite signal for Tn5 sequence bias.

Usage:
atac-correct --bam <reads.bam> --genome <genome.fa> --peaks <peaks.bed>

Output files:
- <outdir>/<prefix>_uncorrected.bw
- <outdir>/<prefix>_bias.bw
- <outdir>/<prefix>_expected.bw
- <outdir>/<prefix>_corrected.bw
- <outdir>/<prefix>_atacorrect.pdf

------------------------------------------------------------------------------------------

Required arguments:
  -b <bam>, --bam <bam>            A .bam-file containing reads to be corrected
  -g <fasta>, --genome <fasta>     A .fasta-file containing whole genomic sequence
  -p <bed>, --peaks <bed>          A .bed-file containing ATAC peak regions

Optional arguments:
  --regions-in <bed>               Input regions for estimating bias (default: regions not
                                   in peaks.bed)
  --regions-out <bed>              Output regions (default: peaks.bed)
  --blacklist <bed>                Blacklisted regions in .bed-format (default: None)
  --extend <int>                   Extend output regions with basepairs
                                   upstream/downstream (default: 100)
  --split-strands                  Write out tracks per strand
  --norm-off                       Switches off normalization based on number of reads
  --track-off [<track> ...]        Switch off writing of individual .bigwig-tracks
                                   (uncorrected/bias/expected/corrected)
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
  --prefix <prefix>                Prefix for output files (default: same as .bam file)
  --outdir <directory>             Output directory for files (default: current working
                                   directory)
  --cores <int>                    Number of cores to use for computation (default: all
                                   available cores)
  --split <int>                    Split of multiprocessing jobs (default: 100)
  --verbosity <int>                Level of output logging (0: silent, 1: errors/warnings,
                                   2: info, 3: stats, 4: debug, 5: spam) (default: 3)
```

### `call-footprints`

Create footprint score tracks from one or more bigWig signals.

```text
usage: call-footprints [-h] [-s <bigwig>] [--signals [<bigwig> ...]] [-o <bigwig>]
                       [--outputs [<bigwig> ...]] [-r <bed>] [--score <score>]
                       [--absolute] [--extend <int>] [--smooth <int>]
                       [--min-limit <float>] [--max-limit <float>] [--scales [<int> ...]]
                       [--multiscale-summary <method>] [--output-multiscale-npz <npz>]
                       [--output-multiscale-npzs [<npz> ...]] [--output-bed <bed>]
                       [--output-beds [<bed> ...]] [--output-bed-dir <directory>]
                       [--top-n <int>] [--min-score <float>] [--call-width <bp>]
                       [--min-distance <bp>] [--fp-min <int>] [--fp-max <int>]
                       [--flank-min <int>] [--flank-max <int>] [--window <int>]
                       [--outdir <directory>] [--cores <int>] [--split <int>]
                       [--verbosity <int>]

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

Parameters for score == sum:
  --window <int>                        The window for calculation of sum (default: 100)

Run arguments:
  --outdir <directory>                  Output directory used with --signals when
                                        --outputs is not supplied
  --cores <int>                         Number of cores to use for computation (default:
                                        all available cores)
  --split <int>                         Split of multiprocessing jobs (default: 100)
  --verbosity <int>                     Level of output logging (0: silent, 1:
                                        errors/warnings, 2: info, 3: stats, 4: debug, 5:
                                        spam) (default: 3)
```

### `match-motifs`

Scan motifs for one or more footprint tracks and infer bound/unbound motif sites.

```text
usage: match-motifs [-h] [--signals [<bigwig> ...]] [--peaks <bed>] [--genome <fasta>]
                    [--motifs [<motifs> ...]] [--motif-db <name>] [--list-motif-dbs]
                    [--cond-names [<name> ...]] [--peak-header <file>] [--naming <string>]
                    [--motif-pvalue <float>] [--bound-pvalue <float>]
                    [--cluster-threshold <float>] [--pseudo <float>] [--skip-excel]
                    [--output-peaks <bed>] [--norm-off]
                    [--normalization {condition-quantile,sample-quantile,none}]
                    [--aggregate-signals [<bigwig> ...]]
                    [--plot-aggregate {sig,all,top,off}] [--plot-aggregate-top-n <int>]
                    [--aggregate-pvalue-threshold <float>] [--aggregate-flank <bp>]
                    [--aggregate-normalization {match,none,sample-quantile,size-factor}]
                    [--aggregate-site-set {all,bound}] [--report-label <text>]
                    [--outdir <directory>] [--prefix <prefix>] [--cores <int>]
                    [--split <int>] [--debug] [--verbosity <int>]

__________________________________________________________________________________________

                                  fp-tools match-motifs
__________________________________________________________________________________________

match-motifs scans motifs in open chromatin regions for one or more footprint score tracks
and infers sample-specific bound and unbound motif sites.

Usage:
match-motifs --signals <footprints.bw> [<more_footprints.bw> ...] --genome <genome.fasta>
--peaks <peaks.bed> [--motif-db jaspar2026_vertebrates | --motifs <motifs.txt>]

Output files:
- <outdir>/<prefix>_figures.pdf
- <outdir>/<prefix>_results.{txt,xlsx}
- <outdir>/<prefix>_distances.txt
- <outdir>/<TF>/<TF>_overview.{txt,xlsx} (per motif)
- <outdir>/<TF>/beds/<TF>_all.bed (per motif)
- <outdir>/<TF>/beds/<TF>_<sample>_bound.bed (per motif)
- <outdir>/<TF>/beds/<TF>_<sample>_unbound.bed (per motif)

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
  --cond-names [<name> ...]        Optional sample name for --signals (default: prefix of
                                   --signals)
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
  --plot-aggregate-top-n <int>     Number of motifs to aggregate when --plot-aggregate top
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
  --report-label <text>            Optional method label shown under the report subtitle
                                   in interactive HTML reports
  --prefix <prefix>                Prefix for overview files in --outdir folder (default:
                                   motif_matches)

Run arguments:
  --outdir <directory>             Output directory to place motif tables, BED files, and
                                   plots in (default: motif_matches_output)
  --cores <int>                    Number of cores to use for computation (default: all
                                   available cores)
  --split <int>                    Split of multiprocessing jobs (default: 100)
  --debug                          Creates an additional '_debug.pdf'-file with debug
                                   plots
  --verbosity <int>                Level of output logging (0: silent, 1: errors/warnings,
                                   2: info, 3: stats, 4: debug, 5: spam) (default: 3)
```

### `diff-footprints`

Compare motif-associated footprint scores across conditions or replicates.

```text
usage: diff-footprints [-h] [--signals [<bigwig> ...]] [--peaks <bed>] [--genome <fasta>]
                       [--motifs [<motifs> ...]] [--motif-db <name>] [--list-motif-dbs]
                       [--cond-names [<name> ...]] [--peak-header <file>]
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
                       [--aggregate-pvalue-threshold <float>] [--aggregate-flank <bp>]
                       [--aggregate-normalization {match,none,sample-quantile,size-factor}]
                       [--aggregate-site-set {all,bound}] [--reuse-existing-results]
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
- <outdir>/<prefix>_figures.pdf
- <outdir>/<prefix>_results.{txt,xlsx}
- <outdir>/<prefix>_distances.txt
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
  --cond-names [<name> ...]        Names of conditions fitting to --signals (default:
                                   prefix of --signals)
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
  --plot-aggregate-top-n <int>     Number of motifs to aggregate when --plot-aggregate top
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

### `normalize-bigwig`

Normalize bigWig tracks using shared background regions.

```text
usage: normalize-bigwig [-h] --bigwigs BIGWIGS [BIGWIGS ...] --background
                        BACKGROUND --outdir OUTDIR
                        [--method {background-scale,background-zscore,none}]
                        [--stat STAT] [--target {median,mean}]
                        [--chrom-sizes CHROM_SIZES]

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
  --outdir OUTDIR       Output directory for normalized bigWigs and QC tables.
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
```

### `plot-aggregate`

Plot aggregate signal around motif sites as PDF/SVG-style output or HTML.

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

### `plot-aggregate-batch`

Compatibility command for interactive aggregate HTML reports. Prefer `plot-aggregate --format html` for new workflows.

```text
usage: plot-aggregate-batch [-h] [--manifest MANIFEST]
                            [--input-html [INPUT_HTML ...]] --output OUTPUT
                            [--flank FLANK] [--top-n TOP_N]
                            [--motifs [MOTIFS ...]]
                            [--site-set {bound,all,unbound}]
                            [--normalization {none,sample-quantile,condition-quantile}]
                            [--default-layout {1x1,1x2,2x2,2x3}]
                            [--title TITLE] [--hide-summary]

Create an interactive aggregate HTML report from match-motifs or embedded
diff-footprints outputs.

options:
  -h, --help            show this help message and exit
  --manifest MANIFEST   TSV with sample, signal, and match_dir columns.
  --input-html [INPUT_HTML ...]
                        Existing aggregate/diff-footprints HTML report(s) with
                        embedded reportPayloadB64 payloads.
  --output OUTPUT       Output self-contained HTML file.
  --flank FLANK         Flank around motif centers for aggregate profiles
                        (default: 100).
  --top-n TOP_N         Number of motifs to preload from manifest mode
                        (default: 30).
  --motifs [MOTIFS ...]
                        Motif prefixes, names, or IDs to preload from manifest
                        mode.
  --site-set {bound,all,unbound}
                        Motif-site BED set to use from match directories in
                        manifest mode (default: bound).
  --normalization {none,sample-quantile,condition-quantile}
                        Profile scaling for manifest mode (default: none).
  --default-layout {1x1,1x2,2x2,2x3}
                        Initial panel grid layout (default: 2x2).
  --title TITLE
  --hide-summary        Hide the TF site summary sidebar in the HTML report.
```

### `run-workflow`

Run a saved YAML workflow config.

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

### `fp-tools-gui`

Launch the optional browser GUI.

```text
usage: fp-tools-gui [-h] [--host HOST] [--port PORT] [--run-dir RUN_DIR]

Launch the fp-tools Streamlit GUI.

options:
  -h, --help         show this help message and exit
  --host HOST        Bind address for the GUI server.
  --port PORT        Optional fixed port.
  --run-dir RUN_DIR  Directory for GUI-managed runs.
```

### `motif-discovery`

Prepare or run a de novo motif discovery command plan.

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

### `motif-summary`

Summarize MEME/Tomtom motif discovery outputs.

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

### `fp-tools-score-variants`

Annotate variants with footprint, sequence, and motif scores.

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

### `pseudobulk-fragments`

Group single-cell ATAC fragments into pseudobulk files.

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

### `find-signature-fp`

Plot per-cell footprint-signature heatmaps and UMAP reports.

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

### `pseudobulk-footprints`

Run grouping, ATAC correction, footprint scoring, reports, aggregate plots, and optional signature reporting.

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
