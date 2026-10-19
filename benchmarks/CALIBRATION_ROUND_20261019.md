# Calibration follow-up — result and decision (19 October 2026)

**Decision: rejected.** Label-free naked-DNA recalibration of the frozen
LOG21 null does not qualify as a general replacement and fails its own
calibration gate. No candidate advances to seed, replicate or resource
panels. Main, released defaults and the manuscript are unchanged; SK-N-SH,
GM23338 labels and naked-DNA replicate 3 stayed closed.

Plan: `CALIBRATION_ROUND_PLAN_20261019.md` (committed `99e135b` before
scoring). Baseline: `CALIBRATION_ROUND_BASELINE_20261019.md`. Model and
interface: `CALIBRATION_ROUND_MODEL_20261019.md`.

## What was tested

One primary candidate (`a + b·logit p` recalibration of the LOG21 central
share, fitted on naked-DNA replicate 1 with each evaluation chromosome held
out), one simpler ablation (intercept only), the unchanged LOG21 count-tail
reference, and a DWM-propensity factorial control with the same detector.
All arms share identical observations and the original 6,020-site nested
cohort (six libraries, full and 10M depth, seed 2026). The reference arm
reproduced every archived next-round score, ceiling and call exactly.

## Findings

**1. Most apparent overdispersion was null-mean error.** On held naked-DNA
chromosomes LOG21 predicts 1,399 central cuts where 1,659 occur (−16%). One
intercept fixes the pooled mean and lowers Pearson dispersion from 1.22 to
1.04; adding the slope gains only 0.31 log-likelihood units, so the ablation
is preferred. Overdispersion modelling in the previous round partly absorbed
this mean error.

**2. Remaining null error is motif-specific.** After recalibration 6/21
control task ratios still exclude 1 (7/21 before), in both directions — CTCF
and MYC over-predicted, K562 FOXA1, MEF2D and TBP under-predicted. GC
tertiles are flat after recalibration (0.98–1.03). No monotone map of `p`
can fix opposite-signed motif errors; that would require changing the
sequence model itself. DWM propensities have 4/21 discrepancies, but its
control model was fitted on the same naked-DNA library (in-sample).

**3. The naked-DNA null does not transfer to chromatin motif sites.** At
ChIP-negative chromatin sites the recalibrated central O/E is 0.34–0.99 by
task (ZNF362/ZNF384 ≈0.34–0.42). Nominal calls at full depth occur in 11.8%
of ChIP-negative and 18.3% of ChIP-positive site/library rows, against ≈0%
on naked DNA. A null that is calibrated on protein-free DNA is therefore not
calibrated in chromatin, and recalibration raises calls on both labels.

**4. Discrimination gains are coverage-driven.** Full-depth mean AUROC:
candidate 0.568, reference 0.557, DWM factorial 0.531, covariates-only
(motif score, accessibility, depth; no central shape) 0.623. Within
coverage strata the candidate is 0.5131 versus 0.5144 for covariates-only and
0.5125 for the reference (DWM factorial 0.504). The +0.011 gain over the
reference shrinks to +0.0005 within strata.

**5. Canonical protection is reproducible only for CTCF.** Positive÷negative
central O/E is 0.69 (HepG2) and 0.67 (K562) for CTCF with intervals below 1.
No difficult family shows it; FOXA1 (both cells), HepG2 MEF2A, K562 ARID3A
and K562 MYC show *more* central cleavage at ChIP-positive sites.

**6. Power is low.** In 138 examined whole-fragment synthetic windows,
removing every central-endpoint fragment is called in 9.5% of resamples
(candidate) and 7.8% (reference); half removal in 3.0% and 1.9%. Windows were
stratified, so these are not population rates.

## Frozen screen (full depth, 20 available of 21 tasks)

| Gate | Result | Evidence |
| --- | --- | --- |
| A calibration | FAIL | held-out log-likelihood −4106.1 vs −4134.4 (better), but 6/21 discrepant (limit 3) |
| B accuracy vs DWM factorial | FAIL | AUROC +0.0373 [0.0204, 0.0535]; relative AP +4.24% (limit 10%) |
| C vs reference | pass | +0.0110 AUROC |
| D CTCF | pass | +0.057 HepG2, +0.027 K562 |
| E families/cells | pass | 7 families, 2 cells (ranking only) |
| F canonical protection | FAIL | 0 difficult families |
| G examined controls | pass | max informative Wilson upper 3.4%; max paired increase 0.6 points |
| H direction | FAIL | 6/6 libraries at full depth, 4/6 at 10M |

The complete promotion ledger is `results/tables/promotion_gate_ledger.tsv`
in the evidence package: accuracy, AP, shape and stability FAIL; untouched
controls and new contexts NOT_TESTED_LOCKED; complete-workflow CPU/memory
and functional uncertainty INCOMPLETE (not run because the screen failed).
The DWM arm here is a raw-propensity factorial control, not the native
fp-tools DWM residual comparator and not TOBIAS.

## What limits footprinting here

In order of evidence: (1) at matched coverage, ChIP-positive motif sites
carry little or no reproducible central depletion except for CTCF;
(2) per-site counts are too low to establish depletion where it might exist;
(3) sequence-only nulls learned on protein-free DNA misplace the chromatin
central mean in motif-specific ways. A different detector cannot fix (1) or
(2). Fixing (3) without learning occupancy needs a label-free chromatin null
reference, which these data do not provide.

**Single justified next diagnostic, if authorized:** fit a sequence-bias
model on nuclear naked-DNA cuts (replicate 1 only, chr1–15, no labels) and
test whether it removes the motif-specific naked-DNA errors on held
chromosomes *before* any chromatin scoring. If it does not, stop sequence-null
work for these factors.

## Verification and evidence

- Plan freeze SHA-256 `b381427a6da07b909ad3d91c7434e6be84c73bbe5ba2e00cae672eb053885e40`;
  `run.py` refuses to start if policy or scientific code differs from it.
  One attempt, exit 0; no failed or repeated attempts.
- All 158 permitted inputs (20 tables/arrays and 138 synthetic-challenge blocks)
  were re-hashed before each use; the next-round
  archive matched its recorded SHA-256. No protected identity was opened.
- Reference parity: the recomputed LOG21 and uniform count-tail arms match
  every archived next-round score, ceiling and call (12 library/depth
  panels, 6,020 sites each arm).
- 13 focused tests pass (`test_core.py`): recalibration recovery and
  likelihood maximum, closed-form intercept, zero-information rejection,
  brute-force tail and ceiling, call/abstain rule, fold exclusion, stratified
  AUROC and the unavailable-task denominator.
- Independent read-only verifier (`verify.py`; generic BFGS refit, explicit
  binomial convolution, rank AUROC, sorted AP): 32 recalibration fits,
  33,200 exact tail recomputations, all calls, 6,300 fold metrics, all
  16 screen summaries with bootstrap intervals, control ratios, gate A/B
  logic and the 21-task ledger — zero failures.
- Full replay from packaged data (`reanalyse.py`, no alignments): all 16
  result tables regenerated and equal to the packaged tables.
- Figures and both slides were rendered and inspected.

Resource evidence: none qualified. The complete-workflow CPU/memory protocol
was intentionally not run because the screen failed; the unrun plan remains
`benchmarks/manifests/end_to_end_resource_protocol_v1.json`. The whole round
ran in a few minutes on 16 workers; that is not a workflow benchmark.

## Evidence archive

`benchmarks/results/fp_tools_calibration_evidence_2026-10-19_eaa31b1.zip`
(ignored output; the filename records the packaging commit).

- Size: 101,661,263 bytes (96.95 MiB), 219 entries including `SHA256SUMS`.
- SHA-256: `ef1e2367c37f610441cf8e719f311538b9b30c954bd3ec079f66b1fd250b07c9`.
- Integrity, single-root, duplicate and path-traversal checks passed; the
  archive was extracted into a new empty directory.
- From that extraction: `python -I -B code/verify.py .` exited 0 (218
  manifested payloads, 32 fits, 33,200 exact tails, 16,086 calls, 6,300 fold
  metrics, zero failures); 13 tests passed; `code/reanalyse.py` regenerated
  all 16 result tables with zero mismatches, writing outside the package.
- Contents exclude alignments, FASTA, bigWig, protected datasets and
  credentials; local paths are replaced by `<LOCAL>`. Verifier log:
  `benchmarks/results/calibration_round_fresh_verifier_v1.log`.
