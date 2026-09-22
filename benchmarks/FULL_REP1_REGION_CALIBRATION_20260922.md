# Higher-count replicate-1 calibration

Status: completed with exit 0 and independent verification. No promotion.

The study uses the audited full filtered replicate-1 alignment at exactly the
same 4,200 label-free chromosome 1–15 sites as the sparse-count calibration.
All 21 tasks remain present. Sequence propensities, strand orientation,
geometry (7, 12, 40), exclusion rules, bootstrap weights and paired likelihood
definition are unchanged. No model is fitted, no threshold is changed, and no
new holdout labels or final naked-DNA control are opened.

This is examined development evidence. Full-input deduplication, mitochondrial
removal and read-level blacklist filtering differ from the first-ten-million
input-pair preparation. It is not a controlled depth-only comparison and does
not qualify independent safety. See `REP1_DEPTH_INPUT_AUDIT_20260922.md`.

## Frozen execution

Script: `benchmarks/scripts/calibrate_full_rep1_regions.py`.
It reuses the verified strand-aware extraction, regional analysis and
independent probability/convolution verifier without changing their source.
Verification rejects absent extraction; partial results are preserved.

Ignored result directory:
`benchmarks/results/full_rep1_region_calibration_20260922_v1/`.

- Study configuration SHA-256:
  `33d716bbaf973027be456da09732decd73672e642995d4dc8123220616b2f039`.
- Job specification SHA-256:
  `4ce505c140baf5bc95703b4b7ed0515331f6762062912208405ff9bd32544617`.
- Supervisor launch: process 958690; child 958691, start ticks 120806942,
  boot `664b7c80-18d6-4147-bdf4-2eae53f766e6`.
- Launch receipt and logs:
  `benchmarks/results/full_rep1_region_calibration_jobs_20260922/execute_attempt1/`.

The job freezes input/code hashes before extraction and freezes the derived
calibration manifest before calculating summaries. It then executes the
independent verifier. Successful process exit must be accompanied by verified
profile and calibration manifests and scientific review.

## Validation

Ten focused tests passed in 3.87 seconds: new orchestration checks plus existing
regional calibration, conditional protection and score-ceiling tests.
The orchestration checks cover existing-output preservation, rejection of
changed policy, and refusal to extract missing profiles in verify mode.
The previous checksum-valid sparse-count result was verified without rerunning
its experiment. No full package regression or resource qualification is claimed.

The complete run exited 0 after 851.3603534430731 seconds. The independent
verifier reproduced all 50,400 site-region rows, 378 calibration/predictive
intervals, 63 likelihood comparisons and 5,670 chromosome blocks. The supervisor
receipt itself does not certify science; verification is recorded in stdout.
This elapsed time is not a complete-workflow resource benchmark.

- Profiles manifest SHA-256:
  `5ded0f45025394bdb9e938010564dc9484a0b2abe449884c95eefe9171f6e543`.
- Derived calibration configuration SHA-256:
  `ec18c2a9876bf57042fddcb7d69efb4957b6993dfefe685176bf3b83dc7271e8`.
- Calibration result SHA-256:
  `5767c81e3a6464a3cc48d08e2e7d6b2a7075514af566f0cd0f81a9779d82e7c8`.

## Scientific findings

| Task | Central cuts | Parametric expectation | Independent-null interval |
| --- | ---: | ---: | --- |
| HepG2 ZNF362 | 20 | 17.33 | 10–26 |
| HepG2 ZNF384 | 19 | 16.49 | 9–25 |
| K562 ZNF384 | 16 | 15.45 | 8–23 |
| HepG2 CTCF | 98 | 97.56 | 81–115 |
| K562 CTCF | 95 | 109.19 | 92–127 |

The zinc-finger deficit remains broadly consistent with sequence-related
cleavage, not transcription-factor protection. Both CTCF observed central
counts now fall within the model-predictive intervals. This does not erase the
earlier K562 discrepancy or isolate depth as its cause: counts and preparation
both differ, and the samples are not independent biological replicates.

Seven of 21 pooled central counts exceed their independent-null predictive
interval: HepG2 FOXA1, HMGA1, MEF2D and TBP; K562 FOXA1, MEF2D and TBP.
Three pooled left-shoulder and two right-shoulder counts are also outside.
Per-strand central discrepancies occur in six tasks on each strand. All are
exploratory, unadjusted comparisons under independent-cut assumptions; biological
coverage is uncalibrated. Better likelihood than uniform does not imply correct
calibration. No occupancy-leakage conclusion follows from these discrepancies.

All 21 tasks and all three regions are shown in
`benchmarks/results/full_rep1_calibration_plots_20260922_v1/calibration.pdf`
(also PNG and SVG). The rendered plot was inspected: labels, all tasks, region
panels, legend and the uncertainty limitation are visible without overlap.
`plot_full_rep1_calibration.py` records input/output checksums and preserves
existing output directories. These are calibration plots, not binding-positive
versus binding-negative aggregate protection plots.

The subsequent unchanged-detector screen is documented in
`FULL_REP1_COUNT_SCREEN_20260922.md`. All evidence remains development-only.
