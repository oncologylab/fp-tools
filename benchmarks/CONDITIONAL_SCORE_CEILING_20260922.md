# Count-conditioned score capacity — 22 September 2026

At the current counts, most positive-label sites cannot reach the fixed
conditional-depletion threshold even if their central count were zero. The
constraint is particularly strong for the examined zinc-finger tasks, where
the parametric bias model predicts little central cleavage. This is a limit
of the frozen score under its model and observed strand totals, not a general
information-theoretic limit or proof that accessibility lacks a footprint.

## Frozen definition

For strand totals `n_plus`, `n_minus` within center plus shoulders, and fixed
central probabilities `p_plus`, `p_minus`, the zero-central-count probability is

`P0 = (1 - p_plus)^n_plus * (1 - p_minus)^n_minus`.

The maximum attainable score is `-log10(P0)`, capped at 300 and rounded to ten
decimal places exactly as the evaluated scorer. It is compared with the
unchanged strict threshold `1.3010299957`. Geometry is the corresponding
held-chromosome geometry for training rows and frozen full-training geometry
for controls. No refitting, threshold change or new score selection occurs.
The probability model's independent-cut and fixed-propensity assumptions
remain unverified for biological calibration.

All 32,424 training rows (8,106 sites across two arms and two detectors) and
8,400 examined replicate-1 control rows remain. Unavailable training tasks
remain explicitly represented in the 210 task/group summaries. Counts that
permit the threshold do not guarantee a call: observed central depletion may
be insufficient. Holding total counts fixed while moving cuts out of the
center is a conditional bound, not a biological alternative or a depth forecast.

## Training results

For positive-label sites under selected parametric geometry:

- 3,589/4,053 (88.6%) cannot reach threshold at their observed strand totals.
- 464/4,053 can reach threshold; 105 are called and 359 are not.
- The fixed-geometry parametric comparator has 312 reachable sites and 40 calls.

| Task | Positive-label sites | Threshold reachable | Called | Median maximum score |
| --- | ---: | ---: | ---: | ---: |
| HepG2 CCCTC-binding factor (CTCF) | 300 | 157 | 50 | 1.3852 |
| K562 CTCF | 245 | 122 | 46 | 1.2914 |
| HepG2 MEF2D | 164 | 4 | 0 | 0.1369 |
| K562 MEF2D | 256 | 3 | 0 | 0.0982 |
| HepG2 ZNF362 | 291 | 0 | 0 | 0.0215 |
| HepG2 ZNF384 | 303 | 0 | 0 | 0.0211 |
| K562 ZNF384 | 238 | 0 | 0 | 0.0225 |

HepG2 MYC illustrates the other limitation: 61/100 positive sites can reach
threshold, but only three do. Count capacity alone does not establish actual
protection. Positive labels derive from chromatin immunoprecipitation (ChIP)
and do not guarantee an observable accessibility footprint.

## Separating counts from sequence probabilities

The conventional geometry (center plus/minus seven bases, shoulders 12–40
bases) uses identical counts in the uniform and parametric arms. Uniform
central probability is 15/73 = 20.55%. Among informative positive-label sites,
the median count-weighted parametric central probabilities are 3.06% for
HepG2 ZNF362, 2.93% for HepG2 ZNF384 and 2.89% for K562 ZNF384. Median
conditioned counts are two, two and three respectively. Their median
parametric expected central counts are only 0.0488, 0.0480 and 0.0543.

With that same fixed geometry, zero zinc-finger positives can reach the
parametric threshold, whereas 8/291, 7/303 and 10/238 respectively can reach
the uniform threshold. This isolates a sequence-probability effect from a
geometry or count-support change. It does not prove the parametric prediction
is correct, that correction is excessive, or that adding depth will recover
binding-associated protection. Calibration against actual naked-DNA central
counts is the next required diagnostic.

## What the control result can establish

Among 3,967 eligible replicate-1 sites in each arm, only two parametric sites
could reach threshold with zero central cuts; neither is called. No uniform
site can reach it. The previous zero-call result therefore provides little
challenge to nominal-tail calibration at this depth. It remains a valid
observed zero-call result, not broad evidence of safety at higher counts.
All 4,200 original sites remain recorded, including overlapping blacklist and
unavailable-model exclusions. These are already examined development controls.

## Verification and artifacts

Seven focused tests passed in 3.69 seconds: analytic versus exact zero-count
tail, bounds across every feasible central count in the small-count grid,
degenerate probabilities, zero counts, capped scores and invalid inputs,
alongside the existing probability/geometry tests. Separate artifact
verification passed. Independent probability products and probability-mass
convolution verified every one of the 40,824 rows and all 210 summaries;
observed scores had zero numerical discrepancy. This is not a new independent
BAM extraction.

Artifacts: `benchmarks/results/conditional_score_ceiling_20260922_v1/`.

| Artifact | SHA-256 |
| --- | --- |
| Configuration | `e973222ffb4cd15edecbddde274ee88091e3e046da00606d2a3133f48c2cec4c` |
| Result manifest | `71d7f0dc198c01595b5fee32bbe2c14dcf6149d9fa1981fa4f1878066bd9db60` |
| Job specification | `59c5ff5889b2871e27a37a2a04c247abb557287770846c41f9d0bf0b152dba07` |

The immutable launch/exit receipts are under
`conditional_score_ceiling_jobs_20260922/execute_attempt1/`; execution exited
zero in 28.323 seconds. This does not measure complete-workflow central
processing unit (CPU) or memory performance.

Next freeze a fixed-geometry central/shoulder calibration audit across all
tasks on permitted naked-DNA replicate 1, separating strands and retaining
the existing exclusions. Compare observed central counts with predictions
before interpreting small predicted probabilities as true sequence-related
cleavage deficits. Any later depth experiment needs verified BAM identities,
fixed candidates and training coordinates. Do not lower thresholds or open
untouched controls to force progress. No promotion gate is newly established;
main, held test labels and final controls remain unchanged.
