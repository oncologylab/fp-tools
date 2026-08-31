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
