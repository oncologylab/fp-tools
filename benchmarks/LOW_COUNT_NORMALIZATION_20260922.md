# Low-count normalization diagnosis — 22 September 2026

The examined detector's control exceedances are associated with finite
normalization amplification, not activation of its near-zero fallback.
This motivates a bounded normalization comparison on training data; it does
not establish causality, improvement or safety for a new operator.

## Frozen scope

`scripts/diagnose_low_count_normalization.py` freezes the existing serialized
geometries, thresholds, sequence-bias caches, training panels, replicate-1 BAM
and index, hg38 blacklist, code and runtime before evaluation. It uses only
chromosomes 1–15. Replicate-1 cuts are freshly extracted; only sequence
propensities and label-free coordinates are reused from the older control
cache. The two cell-context tasks use the same naked-DNA replicate, not two
independent biological controls. No fitting, calibration or threshold change
occurs. Final replicate 3 and new promotion labels remain unopened.

Raw-count and conditional-parametric arms are included. The broader
dinucleotide weight matrix (DWM) control reference remains unqualified and
is not silently substituted. The current frozen detector removes a linear
outer-flank trend, divides by outer root mean square (RMS), uses denominator
one only when RMS is at most 1e-8, clips at plus/minus 12 and rounds geometry
scores to ten decimal places. Scores are not binding probabilities.

The audit records these components and frozen-threshold decisions per site.
Count strata are fixed at 0, 1, 2–4, 5–9, 10–19, 20–49, 50–99 and at least
100 cuts. Training scores use full-training fitted models and are descriptive,
not cross-validated accuracy estimates. Thresholds remain those previously
calibrated on chromosomes 16–18; no tuning labels were newly evaluated.

## Results

All 8,106 matched training sites and 4,200 control sites are retained in the
membership/output records, producing 24,612 rows across two arms. All control
coordinates are on training chromosomes. Controls include 1,569 zero-cut sites,
35 blacklist overlaps and 200 sites with an unavailable HepG2 HMGA1 model;
these categories overlap. There are 2,497 informative eligible control sites
with evaluated thresholds in each arm. No invalid profile was observed.

| Group and arm | Evaluated sites | Threshold exceedances | Median cuts | Median multiplier |
| --- | ---: | ---: | ---: | ---: |
| Training negatives, parametric | 3,948 | 235 | 9 | 2.48 |
| Training positives, parametric | 3,948 | 468 | 10 | 2.48 |
| Replicate-1 controls, parametric | 2,497 | 100 | 2 | 4.46 |
| Replicate-1 controls, raw | 2,497 | 15 | 2 | 6.48 |

The pooled control counts are descriptive: tasks share a BAM, can share
windows, and use different models/thresholds. They are not a general safety
pass. Per-task discrepancies remain important:

| Task | Parametric exceedances / eligible sites |
| --- | ---: |
| HepG2 CCCTC-binding factor (CTCF) | 12 / 120 |
| K562 CTCF | 14 / 123 |
| HepG2 MEF2D | 4 / 128 |
| K562 MEF2D | 11 / 134 |

Among informative parametric controls, zero profiles activate the near-zero
fallback. The minimum RMS is 0.0340, far above 1e-8. Exceeding profiles have
median RMS 0.07345 and multiplier 13.61; other profiles have median RMS 0.22987
and multiplier 4.35. Median cuts are two in both groups. Mean clipped fractions
are 0.00657 and 0.00313 respectively.

The association also appears within individual tasks and count strata. At
exactly one cut, HepG2 CTCF has six exceedances among 54 sites: median
multiplier 20.27 for exceedances versus 14.58 otherwise. K562 CTCF has ten
among 48: medians 18.74 versus 14.80. This is still an association with a
component of the score itself, not an independent causal test. Pooled cut
strata alone do not explain all task-specific behavior.

## Verification and reproducibility

Nine focused tests passed in 2.23 seconds. Separate artifact verification
passed. `scripts/verify_low_count_normalization.py` independently recomputed
all 24,612 component rows with least squares, reconstructed geometry weights,
and checked scores and every threshold decision. It also reconstructed
replicate-1 expected profiles using log-sum-exp probabilities and freshly
extracted strand totals, and checked strand-count sums. This does not constitute
a second independent BAM extraction.

Artifacts: `benchmarks/results/low_count_normalization_20260922_v1/`.

| Artifact | SHA-256 |
| --- | --- |
| Configuration | `0d12d6d355074b002dfd2ae94b899893c5c94d859fb676066876352e1b625244` |
| Profile manifest | `f0c99ae415c22eafe18d36374b025fdcd2a92b0d2b7a419f8e458fff4dd8db12` |
| Diagnostic manifest | `f82c0e028013b4ab466a6dfea444fbb73075074959a3801a731bc735b26e27b4` |

Durable job evidence is in `low_count_normalization_jobs_20260922/execute_attempt1/`.
Execution exited zero in 158.245 seconds. This is diagnostic elapsed time,
not complete-workflow central processing unit (CPU) or memory evidence.

## Next comparison

Freeze a fixed unit scale floor, `denominator = max(RMS, 1)`, as a simple
non-amplifying candidate; compare it against the existing rule using the same
training chromosome cross-fitting, matching and geometry/protection rules.
Do not select per-task floors, retune against replicate 2, or alter the archived
failed models and thresholds. This is a new operator requiring new calibration,
exact-operator safety, shape and full-workflow resource validation if it is
worth advancing. A reduced exceedance rate alone would not establish improved
footprint detection. No promotion gate is newly passed; main is unchanged.
