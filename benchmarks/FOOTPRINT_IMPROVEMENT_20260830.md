# Footprint improvement experiment, 2026-08-30

## Scope

This time-boxed run used branch `research/footprint-improvement-20260830` and
did not change `main`. The locked study is
`manifests/footprint_detectability_v1.spec.json`. Generated public-data outputs
remain under the ignored `results/footprint_detectability_v1/` directory.

The complete executable design contains 206 signal jobs and 17,430 downstream
evaluation tasks. The local evidence matrix comprised 704,615 ChIP-labeled
motif-site records across K562, HepG2, HCT116, A549, MCF-7, and Panc1. Input
paths, byte sizes, and SHA-256 hashes are recorded in the generated evaluation
summary.

## Frozen candidate

The first candidate combined within-cell/TF percentile ranks of the existing
footprint score and PWM score without using ChIP labels:

`1 - (1 - footprint_rank) * (1 - PWM_rank)`

Candidate parameters were frozen after K562/HepG2 development, then evaluated
once on the MCF-7/A549/HCT116/Panc1 cell-line holdout. Chromosomes 17 and 18
were validation chromosomes; chromosomes 19, 20, 21, 22, and X were test
chromosomes. Confidence intervals used 1,000 chromosome-block bootstrap
replicates.

| Split | Baseline AUROC | Candidate AUROC | Delta | Baseline AUPRC | Candidate AUPRC | Relative AUPRC gain |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Development validation | 0.7975 | 0.8515 | +0.0540 | 0.6490 | 0.7820 | +20.5% |
| Development test | 0.7856 | 0.8380 | +0.0524 | 0.5956 | 0.7294 | +22.5% |
| Locked cell-line test | 0.7397 | 0.7499 | +0.0102 | 0.5611 | 0.5930 | +5.7% |

The candidate failed the prespecified locked-holdout mean AUROC, relative
AUPRC, and strong-positive non-regression gates. Its maximum positive-control
AUROC loss was 0.0647. It is rejected, is not a production default, and must
not be retuned on this holdout.

## Biological diagnosis

The result is TF-family dependent rather than a uniform absence of footprint
information:

- CTCF improved in A549 (+0.0487 AUROC), HCT116 (+0.0926), and MCF-7
  (+0.0346), all supported by chromosome-block bootstrap intervals.
- TCF7L2 improved in Panc1 (+0.0256 AUROC), also with bootstrap support.
- MYC regressed in A549 (-0.0647 AUROC) and MCF-7 (-0.0492). In these tasks,
  the existing footprint signal was strong and the PWM ranking was weak.
- FOXA1 (-0.0069) and GATA3 (+0.0013) were not rescued. Their PWM ranking was
  also weak, so adding motif evidence cannot solve their detectability problem.
- REST was consistently weak by footprint score but strong by PWM evidence;
  fusion substantially rescued REST in the supplemental K562, HepG2, and
  HCT116 checks.
- JUND was context dependent: footprint discrimination was strong in K562 and
  HepG2 but weak in HCT116, arguing against a universal JUND failure.
- The available MAX and GATA1 site tables had constant PWM evidence (AUROC
  0.5). This is a motif-evidence coverage or mapping audit flag, not evidence
  that those TFs lack an ATAC footprint.

These results reject one global evidence-fusion rule. They support the next
hypothesis: estimate a label-free evidence-reliability or abstention state and
combine evidence only when the auxiliary channel is demonstrably informative.
That candidate requires a new independent holdout, naked-DNA bias controls,
and perturbation or orthogonal occupancy validation.

## Ten-million-fragment correction slice

One ENCODE replicate each from K562 (`ENCFF077FBI`) and HepG2 (`ENCFF624SON`)
was deterministically sampled at 10 million fragments with seed 2026. Raw,
PWM-corrected, and DWM-corrected signals were footprint-scored at the same
ChIP-labeled motif centers on chromosomes 17, 18, 19, 20, 21, 22, and X. The
comparison contained nine cell/TF tasks and used only sites finite in all three
arms.

| Correction versus raw | Mean AUROC delta | Mean AUPRC delta | AUROC winner | AUPRC winner |
| --- | ---: | ---: | ---: | ---: |
| PWM | +0.00285 | +0.00903 | 4 of 9 tasks | 3 of 9 tasks |
| DWM | +0.00269 | +0.01052 | 4 of 9 tasks | 6 of 9 tasks |

Raw had the highest AUROC only for K562 JUND. DWM had a positive AUPRC point
estimate in all nine tasks. In 1,000 chromosome-block bootstrap replicates,
both correction models robustly improved CTCF AUROC and AUPRC in both cell
lines. DWM also had a robust K562 REST AUPRC gain, but significantly reduced
HepG2 MAX AUROC. The K562 JUND PWM AUROC point estimate was lower than raw and
had only 2.8% bootstrap probability of a positive delta. Direct DWM-versus-PWM
intervals favored DWM for HepG2 CTCF AUROC/AUPRC and K562 CTCF AUPRC, but
favored PWM for HepG2 JUND and MAX AUROC.

This slice does not support global Tn5 overcorrection, but it does show
TF-family-specific correction tradeoffs. DWM remains the default because this
single-depth slice is insufficient to justify a production change. The next
matrix stage should repeat these comparisons across depths, seeds, replicates,
and naked-DNA controls, and should treat strong correction disagreement as a
possible abstention or reliability feature rather than selecting a model from
the same labels.

## Engineering outcomes

- Added an opt-in, label-free evidence-fusion primitive and a locked evaluator
  that records rejected candidates and input provenance.
- Added fixed-site bigWig evaluation so raw, PWM-corrected, and DWM-corrected
  arms can be compared at identical ChIP-labeled motif centers.
- Removed a redundant full BAM counting pass when a validated fragment count
  is supplied to deterministic downsampling.
- Bound generated correction commands to the selected Python environment. A
  provenance check caught and quarantined partial output from an older global
  fp-tools 0.1.8 executable before it entered the evidence matrix.
- Made `python -m fp_tools.cli` a supported execution path for environment-bound
  `atac-correct` plans.

The production footprint scorer remains unchanged because the scientific
promotion gates did not pass.

## Shape detectability and dual-geometry follow-up

The next development stage addressed two different failure modes separately.
`plot-aggregate` now has an opt-in per-site outer-flank RMS normalization,
site-level 95% confidence bands, quantitative central-depletion diagnostics,
and explicit `strong`, `detectable`, `weak`, `not detected`, and
`underpowered` labels. Its default `--share-y none` behavior now gives panels
independent scales instead of silently forcing one global range. On the locked
chromosome test, the DWM-corrected ChIP-positive aggregates were strong for
K562/HepG2 CTCF and JUND, strong for K562 REST, detectable for both MAX tasks
and HepG2 REST, and underpowered for K562 GATA1 (26 positive sites). This mode
makes failure visible but does not treat aggregate shape as proof of occupancy.

A second opt-in `call-footprints --score hybrid` candidate retained the
existing score and added weight 0.2 from a locally standardized 33 bp central
depletion with symmetric 32 bp shoulders. Parameters were selected using
chromosome 17 and frozen before chromosome 18 validation and the locked
chromosome 19/20/21/22/X test. The executable bigWigs exactly reproduced the
prototype metrics.

| Locked test task | AUROC delta | AUPRC delta |
| --- | ---: | ---: |
| K562 CTCF | +0.0214 | +0.0323 |
| K562 REST | +0.0169 | +0.0146 |
| HepG2 CTCF | +0.0313 | +0.0400 |
| HepG2 REST | +0.0188 | +0.0032 |
| Mean across all nine tasks | +0.0064 | +0.0080 |

Chromosome-block bootstrap intervals supported both CTCF gains and the K562
REST gains; HepG2 REST had 98.3% and 99.7% probabilities of positive AUROC and
AUPRC deltas. The same test showed significant regressions for HepG2 JUND/MAX
and K562 MAX, while K562 GATA1 remained underpowered. The hybrid arm is useful
for testing wide-footprint families but fails the strong-positive
non-regression requirement and is not the default. The next scorer should use
motif-width or label-free geometry reliability to select the wide channel,
then be tested on a new independent cell-line holdout and naked-DNA controls.

As a post hoc routing check, a label-free minimum motif width of 15 bp selected
hybrid scores only for CTCF (15 bp) and REST (20 bp), while retaining the
standard score for GATA1, JUND, and MAX. On the locked chromosome test this
raised mean AUROC by 0.0098 and mean AUPRC by 0.0100 with zero task-level
regressions. This gate was evaluated after viewing the universal-hybrid result,
so it is a hypothesis for a new preregistered experiment, not a validated
promotion result. The full matrix now includes an explicit `fp_tools_hybrid`
method arm so that hypothesis can be tested without overwriting the standard
score.

## Per-TF geometry search with accessibility-matched labels

The universal hybrid was too small and inconsistent to justify production.
The next experiment therefore searched each TF separately across raw, PWM,
and DWM signals; 5--49 bp centers; 4--64 bp shoulders; 0--12 bp gaps; mean or
minimum shoulders; mean, minimum, or lower-quartile centers; five local scale
models; and four asymmetry penalties. Geometry was screened on chromosomes
1--16, candidates were shortlisted on chromosomes 17--18, and chromosomes
19--22/X were not read until the candidates were frozen.

ENCODE discovery found released replicate-aware GRCh38 IDR ChIP peaks for all
21 K562/HepG2 development tasks. Twelve unique JASPAR 2026 motifs produced
cell-specific motif scans within each cell line's released ATAC peaks.
Positive labels required ChIP peak and summit support; ambiguous near-peak
sites were excluded. A fivefold negative pool was optimally matched to
positives on motif score and local raw ATAC coverage. Tasks without common
support or adequate positive counts remain explicitly confounded or
underpowered.

The quick staged search evaluated roughly two thousand hypotheses per TF and
selected different correction/geometry families rather than one global rule.
The frozen matched-site test produced:

| Task | Test AUROC delta | Test AUPRC delta | Candidate test AUROC |
| --- | ---: | ---: | ---: |
| K562 CTCF | +0.1615 | +0.1809 | 0.7732 |
| HepG2 CTCF | +0.1383 | +0.1651 | 0.7531 |
| K562 MEF2A | +0.1439 | +0.1312 | 0.6639 |
| K562 ARID3A | +0.0752 | +0.0579 | 0.5931 |
| HepG2 MYC | +0.0768 | +0.0379 | 0.6075 |
| K562 MEF2D | +0.0593 | +0.0517 | 0.6005 |

The same test rejected universal application: K562 FOXA1, MYC, and TCF7L2;
HepG2 FOXA1, ZNF384, and ZNF558 were flat or worse, and several had inadequate
common support or fewer than 100 positive test sites. Aggregate figures show a
large, visually distinct CTCF depletion and paired shoulders, while many weak
TFs retain overlapping, noisy positive/negative profiles. Thus the principal
problem is a mixture of fixed-kernel mismatch, TF-specific correction effects,
low depth, underpowered occupancy labels, and genuine absence of stable
aggregate shape—not correction alone.

Applying the prespecified minimum of 500 positive sites plus matching-balance
checks leaves only K562 and HepG2 CTCF as strong, balanced, point-gain-passing
results. HepG2/K562 ZNF384 are adequately powered but not detected; HepG2
ZNF362 and K562 MYC remain accessibility-confounded; the remaining tasks are
underpowered on test chromosomes even when their exploratory deltas are large.
This separation prevents attractive small-sample gains from being presented as
validated TF improvements.

Paired chromosome-block bootstrap supported both CTCF gains in both metrics:
HepG2 AUROC 95% interval 0.094--0.180 and AUPRC 0.123--0.214; K562 AUROC
0.088--0.198 and AUPRC 0.119--0.225. The exploratory K562 MEF2A gain also had
positive intervals in both metrics, and HepG2 MYC did likewise, but each had
far fewer than 500 positive test sites. K562 MEF2D supported AUPRC improvement
but its AUROC interval crossed zero. These uncertainty results reinforce the
prespecified distinction between an interesting small-sample signal and a
promotion-eligible result.

These results justify continued depth, randomization, replicate, and new-cell
validation, but not a main-branch scorer change. A future implementation must
route only validated TF/family geometries and abstain when positive-site count,
matching balance, replicate stability, or aggregate-shape reliability fails.

Holding each frozen geometry fixed and changing only raw/PWM/DWM input showed
large TF-specific correction ranges: AUROC ranges were 0.165 for HepG2 MYC,
0.175 for HepG2 ZNF558, 0.164 for K562 FOXA1, and 0.116 for K562 MYC. Raw was
best for both MYC tasks and several zinc-finger tasks; DWM was best for both
CTCF tasks, K562 MEF2A, and K562 MEF2D. HepG2 ZNF558's validation-selected DWM
arm reversed on test and raw was best post hoc. Thus correction can be the
principal problem for particular TF/context pairs, but it is not globally
under- or over-correcting. Correction disagreement and cross-split instability
should become abstention evidence rather than a label-tuned automatic choice.

The exhaustive follow-up expanded the quick screen to roughly six thousand
configurations per TF. It preserved the CTCF result and modestly increased the
exploratory K562 MEF2A and ARID3A point estimates, but did not broadly rescue
the remaining TFs. Some selected configurations transferred worse than the
smaller screen. This is direct evidence that a larger geometry search can add
selection overfit without adding footprint information; exhaustive search is
therefore a development stress test, not a method improvement by itself.

Full-depth DWM signals from all three available replicates per cell line were
then scored with the frozen quick-screen geometries. Replicate averaging raised
CTCF AUROC from 0.761 to 0.803 in HepG2 and from 0.773 to 0.818 in K562. It also
raised HepG2 ZNF362 and ZNF384 by about 0.069 AUROC, although both remained
below a strong-information threshold. K562 MEF2D increased by 0.015. In
contrast, HepG2 MEF2A stayed at chance, HepG2 MYC fell by 0.036 with DWM, and
K562 MEF2A fell by 0.024. More reads and replicate averaging therefore rescue
sampling-limited tasks such as CTCF, but do not repair a wrong correction arm,
an unstable aggregate shape, or absent occupancy information.

No production scorer is promoted from this matrix. A substantial next method
must combine TF/family-specific geometry with correction-reliability and
minimum-information gates, and it must report `not detected` or `underpowered`
instead of manufacturing a confident score for unsupported TF/context pairs.

## Learned-profile stress test

To test whether the hand-built score was hiding a more complex footprint, each
cell/TF task was also evaluated with 228 learned-profile hypotheses. The matrix
crossed raw, PWM, DWM, and three-signal features; full and symmetry-folded
profiles; unscaled and outer-flank RMS-normalized inputs; and matched-template,
regularized logistic, shrinkage-LDA, and ExtraTrees models. Full profiles were
tested in both genomic and motif-strand orientation. Model and feature
choice used chromosomes 17--18 only, followed by refitting on development
chromosomes and one evaluation on chromosomes 19--22/X.

The validation-selected learned models reached AUROC/AUPRC 0.799/0.808 for
HepG2 CTCF and 0.800/0.777 for K562 CTCF. They did not rescue the adequately
powered ZNF384 tasks (0.565/0.554 HepG2 and 0.607/0.592 K562), and HepG2 MEF2A
and MEF2D remained near chance. HepG2 FOXA1 reached 0.778 AUROC and K562 MYC
0.825, but their positive and negative sites retained large accessibility
imbalances (absolute SMD 0.822 and 0.842), so those values cannot be attributed
to footprint shape.

The learned models therefore agree with the geometry, correction-transfer,
depth, and aggregate-profile experiments: CTCF contains strong transferable
shape information, a few low-count tasks are promising but unresolved, and
several TFs lack a stable discriminative profile in the available ATAC data.
The generated complete experiment matrix assigns every cell/TF task a power,
confounding, response, and likely-driver state so unsupported tasks cannot be
silently averaged into an attractive global result.

Strand orientation did not produce a promotion-eligible rescue. HepG2 TBP and
K562 TCF7L2 had large oriented-profile test gains, but orientation had reduced
their validation AUROC; K562/HepG2 ARID3A had small consistent-looking gains
with too few test positives. No task improved by at least 0.03 AUROC on both
validation and test. K562 MEF2D reached 0.654 AUROC with a full-depth,
replicate-mean oriented DWM model, but only 154 positive test sites were
available. Directional footprint shape remains an explicit hypothesis for
larger datasets, not a validated general rescue.

The validation-selected DWM classifiers were also transferred to an
independent 10-million-fragment read sample. CTCF changed by only +0.002 AUROC
in HepG2 and -0.006 in K562. In contrast, absolute seed shifts exceeded 0.05
for HepG2 MEF2D and ZNF558 and for K562 ARID3A, FOXA1, and MEF2A; all five are
underpowered or confounded. This independently confirms that CTCF has a stable
learnable shape while several attractive small-task models are sampling
artifacts or unresolved leads.

The complete raw/PWM/DWM correction comparison across both read seeds reached
the same conclusion. DWM was the CTCF winner for both metrics in both cells
and both seeds; the maximum correction-arm seed shift was 0.0085 AUROC in
HepG2 and 0.0132 in K562. Raw remained the K562 MYC winner. K562 MEF2A kept
DWM as its winner but remained underpowered, while K562 MEF2D switched from
DWM to PWM and had a maximum 0.0836 AUROC seed shift. HepG2 ZNF558 retained
raw as the AUROC winner, but its DWM arm varied by 0.125. These results support
stable TF-specific routing only for CTCF; MEF2D and ZNF558 require abstention
until larger independent data resolve correction instability.

The final aggregate figures expose the same decisions visually. Panel headers
are color-coded and explicitly labeled `DETECTED / RESCUED`, `UNDERPOWERED`,
`ACCESSIBILITY-CONFOUNDED`, or `WEAK / CONTEXT-DEPENDENT`. This prevents a
generic motif-centered depletion, such as the nearly overlapping ZNF384
positive and negative profiles, from being mistaken for occupancy detection.

Finally, DWM-only learned profiles were frozen from the development
chromosomes and transferred independently to all three full-depth biological
replicates. The replicate-mean classifier reached AUROC/AUPRC 0.837/0.853 in
HepG2 CTCF and 0.834/0.812 in K562 CTCF. Individual replicate AUROCs were
0.823--0.835 (standard deviation 0.006) and 0.822--0.830 (standard deviation
0.004), respectively. This is the strongest evidence in the matrix that the
CTCF improvement is a stable shape result rather than a one-sample artifact.
Powered ZNF384 remained reproducibly weak (individual AUROC 0.616--0.622 in
HepG2 and 0.617--0.636 in K562), so additional reads alone do not turn its
aggregate depletion into occupancy discrimination.

## Direct expected-bias diagnosis

The raw, expected-bias, and corrected profiles were next scored with identical
frozen geometry for every available TF under both PWM and DWM correction. This
separates a correction effect from a geometry effect. For CTCF, the expected
bias score opposed occupancy labels (AUROC 0.423 HepG2 and 0.450 K562 under
DWM), and subtracting it improved AUROC by 0.064 and 0.044. CTCF is therefore
not a case of harmful global overcorrection; correction and TF-specific
geometry solve different parts of its error.

The opposite pattern appeared for several TFs. DWM correction changed AUROC by
-0.165 for HepG2 MYC, -0.175 for HepG2 ZNF558, and -0.100 for K562 MYC. The
expected-bias track itself predicted the labels for K562 MYC (AUROC 0.791) and
HepG2 ZNF558 (0.605). Consequently, raw-signal discrimination in those tasks
partly follows sequence-bias structure. Its removal lowering AUROC is not by
itself proof of overcorrection; it may be removal of a false predictive cue.
Naked-DNA Tn5 and orthogonal occupancy data are required to distinguish those
possibilities.

Large corrected-versus-expected score correlations remained in many MEF2,
TBP, FOXA1, MYC, and ZNF384 arms. These correlations are diagnostic flags, not
proof of residual enzyme bias, because expected signal, accessibility, and
motif sequence are not independent. Together with the correction-transfer and
random-seed results they identify the specific next correction improvement:
estimate per-TF reliability from naked-DNA residuals and abstain when PWM/DWM,
read seeds, or biological replicates disagree.

## Untouched WTC11 MAX transfer

The LOG81/+4,-4 shared-strand anchored-FDA route was preregistered for WTC11
MAX before the ENCODE ChIP peak content was downloaded. The frozen evaluation
initially appeared strong (AUROC +0.105 and relative AUPRC +13.5%, with both
biological replicates positive), but the matched cohort failed its covariate
gate: the maximum absolute standardized difference was 1.098. That result is
retained as a failed locked test rather than treated as validation.

An explicitly post-unblinding sensitivity then selected the largest
common-support subset using only motif score, accessibility, GC, and peak
position. The selected 306 sites per class had maximum absolute SMD 0.098.
On this balanced cohort the frozen route improved AUROC by 0.033 and relative
AUPRC by 6.8%; chromosome-block bootstrap probabilities of positive gain were
0.860 and 0.829. It therefore failed the prespecified AUPRC and uncertainty
gates.

A complete post-lock factorial evaluated 284 combinations across conventional
DWM and mitochondrial-trained LOG81 and SELMA10 models, both +4,-4 and +4,-5
shifts, and spline/FDA/GP/hybrid detectors. No configuration passed both point
gates on the balanced cohort. The best AUROC gain was +0.045 (SELMA10/+4,-4,
anchored antisymmetric FDA) with only +5.3% relative AUPRC; the best LOG81
shared-strand FDA reached +0.039 AUROC and +7.1% relative AUPRC. A supervised
functional-PC ceiling also remained below chance after balance. WTC11 MAX is
therefore classified as low-information/assay-limited in this dataset rather
than rescued by a different bias model, shift convention, or functional
detector. No performance PDF is eligible for external sharing from this test.

## Covariate-residualized FDA and null-calibrated calls

Residualizing functional-PC scores against motif score and log accessibility
before the unsupervised mixture improved mean AUROC by only 0.008 across the
complete K562/HepG2 bias-and-shift factorial and reduced mean relative AUPRC
by 23%. Some apparently large TBP, MEF2, and ZNF AUROC gains were accompanied
by substantial precision loss, while ARID3A and CTCF often regressed. This
variant is rejected as a detector improvement.

The first naked-DNA experiment exposed a more fundamental problem with every
two-state detector: a mixture always partitions its input, even when no TF is
bound, so a universal posterior threshold of 0.5 manufactured a bound state.
Candidate false-positive rates reached 42--48% of all sites for several count
models, and the conventional DWM references also reached 12--40%. Posterior
probabilities are therefore useful ranks but are not calibrated binding
probabilities.

A label-free empirical threshold from naked-DNA chromosomes 1--15, using the
upper 2.5% tail, transferred to naked-DNA chromosomes 16--18 with every one of
30 candidate/reference groups below the 5% safety ceiling. The maximum
all-site rates were 2.76% for candidates and 2.55% for DWM; maximum rates among
informative sites were 4.74% and 4.21%. This fixes the artificial-null call
rate but does not change AUROC or AUPRC.

Naked DNA alone did not reproduce broad accessibility structure in cellular
ATAC. A second label-free null was therefore constructed by cyclically moving
each motif-oriented ATAC profile by 35 or 45 bp in either direction. This
preserves per-site coverage and cut distributions while breaking alignment to
the motif center. The dual-null rule takes the stricter per-TF threshold from
naked DNA and shifted cellular ATAC. On the untouched naked-DNA chromosome
panel, all 30 groups again remained below 5%; the shifted-ATAC null became the
limiting threshold for every promoted count detector.

On development validation chromosomes, the dual-null rule reduced candidate
ChIP-negative call rates for MEF2, TBP, and ZNF models from roughly 13--25% to
0.2--1.8%, while retaining modest sensitivity. Continuous-score improvements
over DWM remained large for several tasks, including HepG2 MEF2D (+0.236
AUROC), HepG2 MEF2A (+0.152), K562 ZNF384 (+0.140), and K562 MEF2D (+0.108).
HepG2 FOXA1 remained unsafe at 14.6% and is an explicit abstention case.

These are development-label results, and the corresponding family routes did
not transfer on the already-scored GM12878/IMR-90 holdout. The 2.5% dual-null
choice was also made after observing the first null-calibration behavior.
Consequently it is a mechanistic advance and a preregistered hypothesis for a
new holdout, not a package promotion. A new cell/TF occupancy panel is still
required. No Box performance folder is created from these label-free results.

The shifted-ATAC calibration was then made fourfold cross-fitted: every null
score is produced by a model that did not fit that site's fold. Thresholds
calibrated on offsets -45, -35, +35, and +45 bp transferred to separate -40
and +40 bp offsets. All 15 candidate groups passed the 5% ceiling; the maximum
all-site and informative-site rates were 3.50% and 3.59%. A prespecified
secondary tail of 2.5% was retained because 5% raised the maximum held-out
informative rate to 7.40% and left only 7 of 15 candidates passing. The 1%
tail also passed but discarded more true cellular signal.

Applying the cross-fitted dual thresholds to the full-depth first naked-DNA
replicate again left all candidate groups below 5%, with maxima of 3.00% over
all sites and 3.21% over informative sites. The count-GP MEF2 and TBP routes
made no false calls in that panel. This is useful depth transfer, but it is
not independent replication because the 10M calibration and full-depth panel
come from the same biological library.

The frozen policy was then applied to the second GSE164997 naked-DNA library
without pooling or refitting. The first 10 million input pairs yielded
7,472,170 coordinate-sorted, deduplicated, properly paired fragments after the
same alignment and MAPQ filters. All 15 TF/cell candidate groups passed: the
maximum false-call rate was 2.00% over all motif sites and 3.54% among sites
with nonzero signal. The MEF2, TBP, ZNF, and MYC count routes made no false
calls. The uncalibrated 0.5 mixture cutoff still reached 49.5% over all sites
and 84.3% among informative sites, independently confirming that mixture
posteriors are ranks rather than binding probabilities. Candidate-only scoring
was verified to reproduce all 3,000 prior-replicate probabilities exactly, and
the frozen-threshold applicator reproduced every prior full-depth rate before
the second replicate was opened. This supplies independent biological-library
negative-control evidence. Wilson 95% upper bounds exceeded 5% only for FOXA1
(8.75% in HepG2 and 5.28% in K562); the maximum upper bound for every non-FOX
route was 4.39%. The FOX family is therefore an explicit abstention pending
more independent null depth. This result does not rescue the failed external
occupancy-transfer gate and therefore does not justify package promotion.

Finite-null empirical p-values and Benjamini-Hochberg q-values were added as
an audit output. They generated no false discoveries in naked DNA, but also no
true calls in development ATAC: roughly 800 null scores per TF do not provide
the tail resolution needed for site-level FDR across approximately one
thousand sites. Generalized-Pareto tail extrapolation was tested and rejected;
it removed useful count-model calls and produced unsafe FOXA1 artifacts. The
supported result is therefore conservative TF-level detectability plus
dual-null thresholded candidates, not calibrated site-level FDR.

Separately, the earlier CTCF geometry result now has concise research-only
before/after reports in the yilab Box folder
`Yaoxiang/fp-tools/2026-08-31_CTCF_footprint_improvement_research`. The reports
show the held-out paired ROC/PR curves, aggregate profiles, chromosome-block
bootstrap intervals, and three-biological-replicate transfer. Their scope is
explicitly CTCF-specific; this upload does not change the package default or
override the external-validation requirements for a general method.

A post-hoc, label-free sensitivity analysis fit score-to-covariate regressions
on four chromosomes and residualized the fifth, rotating across chr19--22/X.
After removing motif-score, accessibility, and central-accessibility effects,
the TF-specific geometry still improved AUROC/AUPRC by +0.071/+0.069 in HepG2
and +0.091/+0.112 in K562. The concise Box PDFs and metrics table now include
this result, explicitly labeled post hoc; it strengthens the shape
interpretation without converting the CTCF result into a universal claim.

An exact CTCF naked-DNA check then froze separate candidate and conventional
score thresholds on 1,000 ChIP-negative validation sites per cell line before
applying them to GSE164997 naked-DNA motif sites. The frozen TF geometry called
1/116 finite HepG2 sites (0.86%; Wilson upper 95% bound 4.72%) and 1/123 finite
K562 sites (0.81%; upper bound 4.46%). Both point estimates and uncertainty
bounds therefore pass the prespecified 5% negative-control target. The raw
conventional-score thresholds transferred poorly across samples (101/116 and
111/123 calls), which is evidence of score-scale non-portability rather than a
biological occupancy estimate. The Box reports now show the candidate safety
result. This closes the exact naked-DNA check for the CTCF geometry but does not
replace the pending external-cell transfer gate or justify package-wide
promotion.

## Depth-matrix artifact integrity audit

An incremental 25M refresh exposed a benchmark-infrastructure failure before
its metrics were used: four interrupted HepG2 corrected bigWigs had parseable
headers and large files but zero covered bases. The ablation runner had treated
existence as completion and marked them `skipped_existing`. Three additional
zero-coverage 50M files were active writers at the time of the audit and were
correctly recognized as incomplete once their runner state was considered.

The contaminated seed-2027 evaluation was quarantined. The runner and depth
evaluator now require every bigWig input to have chromosomes and at least one
covered base; the evaluator also validates raw and expected tracks and binds
profile caches to signal size and modification time. Focused runner/depth tests
pass. A repair run is rebuilding the four stale 25M tracks rather than treating
them as biological evidence. Until those repairs complete, only the valid 10M
matrix and the valid seed-2026 25M slice may be interpreted.
