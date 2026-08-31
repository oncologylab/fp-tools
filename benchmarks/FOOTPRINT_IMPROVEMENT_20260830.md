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
cell/TF task was also evaluated with 152 learned-profile hypotheses. The matrix
crossed raw, PWM, DWM, and three-signal features; full and symmetry-folded
profiles; unscaled and outer-flank RMS-normalized inputs; and matched-template,
regularized logistic, shrinkage-LDA, and ExtraTrees models. Model and feature
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
