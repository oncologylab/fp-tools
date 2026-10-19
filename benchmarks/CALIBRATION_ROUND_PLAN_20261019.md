# Calibration follow-up: frozen plan — 19 October 2026

Machine-readable policy: `benchmarks/scripts/calibration_round/protocol.py`.
The freeze (`PLAN.freeze.json`) hashes that policy, the scientific code and
every permitted input before any new score is computed.

## Diagnosis from existing evidence (no new scoring)

The next-round tables already contain within-stratum metrics. Averaged at full
depth over 20 available tasks (equal library, then task), the held nested
cohort shows:

| Score | AUROC, all sites | AUROC, accessibility quartiles 1–4 | AUROC, count strata 10–99 |
| --- | ---: | --- | --- |
| LOG21 count-tail | 0.557 | 0.588 / 0.526 / 0.495 / 0.549 | 0.486 / 0.511 / 0.469 |
| DWM binomial | 0.497 | 0.514 / 0.503 / 0.495 / 0.519 | 0.470 / 0.499 / 0.461 |
| Covariates only (no central shape) | 0.623 | — | — |

Within comparable coverage, no detector separates ChIP-positive from
ChIP-negative motif sites much better than chance. Residual canonical
protection is reproducible only for CTCF (3/3 libraries in each cell);
MEF2D, ZNF362 and TBP show it in at most two of three libraries. The
dominant limitation for ranking is therefore the observations — absent or
unmeasurable protection at most labeled sites, plus coverage confounding —
not the detector's functional form. A more flexible detector cannot create
protection the counts do not contain.

What *is* correctable is the null mean: 7/21 LOG21 central-mean control
intervals exclude zero, and dispersion was modeled before that mean was
corrected. A wrong null mean biases every expected profile, residual, call
and abstention decision, whatever the detector.

## The single candidate

**Verified failure addressed.** LOG21, trained on chromatin mitochondrial
reads, mispredicts the central-versus-shoulder cut share on protein-free
nuclear naked DNA for several tasks.

**Why the existing method cannot address it.** LOG21 is frozen and must not
be refit against occupancy labels. Its coefficients are shared across
positions; a transfer error in the central share is not visible to the
detector, which treats `p` as exact.

**Candidate (`LOG21NC_count_tail`).** Per held chromosome `f`, fit
`logit q = a + b·logit p_LOG21` on naked-DNA replicate-1 full-depth windows
from chr1–15 excluding `f` (deduplicated physical windows, strand
observations, no labels, no chromatin data). Score with the unchanged
strand-conditioned count-tail detector (cap 300, ten decimals) and the
unchanged strict threshold 1.3010299957. Sites whose zero-central-cut
ceiling cannot exceed the threshold are reported as `ABSTAIN_UNREACHABLE`,
never as negative calls.

**Simpler ablation (`LOG21NC1_count_tail`).** Intercept only (`b = 1`).
**Reference.** Unchanged `LOG21_count_tail` (zero added parameters).
**Factorial control.** `DWM_count_tail`: library DWM propensities with the
same detector. This is not native production DWM correction.

**Falsified if** held-chromosome control log-likelihood does not improve,
discrepant control tasks do not fall to ≤3/21, or discrimination drops more
than 0.01 AUROC below the reference.

**Added complexity.** Two scalars per deployment, and a requirement for a
protein-free control library from a comparable Tn5 protocol.

## Populations and endpoints

Original nested cohort only (6,020 sites), six libraries, depths full and
10M (seed 2026). No rematching. All 21 tasks; HepG2 HMGA1 is unavailable
for chromatin metrics but evaluable on controls. All 4,200 control windows
are retained with exclusion reasons. Control endpoints: pooled held-out
log-likelihood, per-task observed/expected central ratio with a 2,000-draw
chromosome-block interval, Pearson dispersion after each mean, calibration
curves by LOG21 `p` and GC, calls and Wilson bounds. Chromatin endpoints:
the next-round fold metrics and aggregation, count-stratified AUROC,
abstention, regional label-level central O/E and the whole-fragment
synthetic challenge re-scored with each null.

## Single advancement screen (full depth)

A calibration, B accuracy (≥0.03 AUROC and ≥10% relative AP versus
`DWM_count_tail`, interval excluding zero), C no loss >0.01 versus the
reference, D CTCF loss ≤0.02, E two families and both cells improved,
F canonical regional protection in two difficult families, G control calls
(Wilson upper ≤5%, paired increase ≤1 point), H direction in ≥5/6 libraries
at both depths. The ablation is preferred unless the primary adds ≥1.92
held-out log-likelihood units.

If only A passes, the result is a development-qualified null correction and
a rejected general replacement. If A fails, the candidate is rejected. Seed
panels, other naked-DNA replicates and the complete-workflow resource
protocol run only after a full pass.
