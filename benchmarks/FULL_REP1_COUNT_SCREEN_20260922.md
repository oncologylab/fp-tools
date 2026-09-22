# Frozen count detector on higher-count replicate-1 controls

The parametric detector makes two nominal calls, both HepG2 CTCF. All 20 tasks
with an available detector meet the descriptive 5% Wilson upper-bound check
on both eligible and positive-count denominators. This is examined development
evidence, not untouched-control safety or promotion evidence.

## Frozen design and results

Apply the exact existing 42 training model records and unchanged rounded
threshold 1.3010299957, with strict greater-than calls. Require valid, nonblacklisted
sites, an available detector and positive conditioned counts. Retain all 4,200
sites and both arms; HepG2 HMGA1 remains unavailable. No model or threshold is
refitted. Full-input preparation differs from the smaller alignment, so this
is not a controlled depth-only comparison.

| Arm | Eligible sites | Informative sites | Calls | Tasks meeting both bounds |
| --- | ---: | ---: | ---: | ---: |
| Conditional parametric | 3,967 | 3,343 | 2 | 20/20 available |
| Conditional uniform | 3,967 | 3,342 | 5 | 19/20 available |

Parametric HepG2 CTCF has 2/162 informative calls, Wilson upper bound 0.043891;
K562 CTCF has 0/184, upper bound 0.020451. Uniform HepG2 CTCF has 4/185,
upper bound 0.054266 (does not meet the bound); K562 has 1/184, upper bound
0.030138. Every other available task has zero calls in both arms. Missing
HepG2 HMGA1 is retained as unavailable, not counted as passing.

The analytic ceiling remains a major limitation. Reapplying the previously
tested `audit_conditional_score_ceiling.score_ceiling` to these saved strand
counts and probabilities yields only 65 eligible parametric sites and 48
uniform sites capable of reaching the unchanged threshold even with no central
cuts. Parametric CTCF contributes 36 HepG2 and 13 K562 reachable sites. No
examined ZNF362/ZNF384 site is reachable. The ceiling is the rounded, capped
negative base-ten logarithm of `(1-p_plus)^n_plus * (1-p_minus)^n_minus`.
It is not a biological power calculation. Thus the low observed call rate
still provides a limited challenge to the tail model and cannot establish
performance at greater per-site counts or strong sensitivity to protection.

The conventional dinucleotide weight matrix (DWM) comparator was not rescored
here: the required paired false-positive increase gate remains unestablished.
Independent-cut assumptions, correlated sites, uncalibrated biological interval
coverage and the already examined control prohibit a promotion claim.

## Verification and provenance

Twelve focused tests passed in 3.85 seconds across orchestration, calibration,
conditional scores and ceilings. The supervised screen exited 0 after
79.27278033597395 seconds, including verification; this is not a whole-workflow
resource measurement. Independent reconstruction verified all 42 models,
8,400 site rows and 42 rate/Wilson summaries, with maximum score discrepancy
zero. Model records were copied exactly from their frozen artifact and checked
for equality before scoring; independent verification rechecks the historical
training selection rather than altering it.

Ignored outputs: `benchmarks/results/full_rep1_count_screen_20260922_v1/`.
Durable receipts: `benchmarks/results/full_rep1_count_screen_jobs_20260922/execute_attempt1/`.

- Parent configuration SHA-256:
  `35918c9d3a6e527dfcdb2eebeb8257c95b7475bed4b53b46a789b2ef272c8b3a`.
- Job specification SHA-256:
  `66b27c141e6ec1d29a225dea73305a84b203fe564f04a5e2b5fac517f51e0eeb`.
- Derived screen configuration SHA-256:
  `aa0384ce8d55fb0fd4fe07029123f90adc6efc48cb6a6442eb7c5dc87255c1d1`.
- Model artifact SHA-256:
  `74ca0dd12ed93e6a4f0976935a442f0662b39cc0d61a4376f33474900c232cdf`.
- Screen lock SHA-256:
  `618f912fe74bb5dfe0d4ec89ead750604036db263610381cfaf019a98ab18d54`.
- Screen result SHA-256:
  `58a0f25efffe789c7011857917ca3f7f3a6b2d18657a59756a10a3827d4d51b1`.

Next audit and freeze higher-depth training-chromosome inputs to evaluate the
unchanged fold-specific detector with more counts, retaining pair membership
and auditing coverage/motif/accessibility balance. Do not choose geometry on
examined test chromosomes or adjust thresholds to observed controls. New
holdout labels and final naked-DNA replicate 3 remain unopened.
