# Frozen shape-operator development screen

The spline tuning gains do **not** establish a safe footprint improvement.
The exact frozen parametric spline failed the prespecified false-positive-rate
(FPR) screen in both CTCF contexts. Retain these failures; do not change model
weights, thresholds, or candidate selection after examining this control.

This is a development screen on previously examined naked-DNA replicate 2.
It is not untouched negative-control evidence. No new holdout labels or final
naked-DNA replicate 3 were opened. No model was refitted. Main and the
production dinucleotide weight matrix (DWM) default were unchanged.

## Method and support

All 11 fixed-normalization artifacts were evaluated with their trained
geometry and hierarchical spline. The ledger retains 21 transcription-factor
tasks × 3 corrections × 2 detectors = 126 rows: 22 have models and calibrated
thresholds; 104 explicitly lack a fitted model. Four of the 22 have no valid
DWM control expectation and consequently cannot produce an informative rate.

Thresholds were calibrated only on chr16–18 negatives with positive observed
counts and valid positive expectation in both correction arms. All thresholds
were frozen before control scoring. The baseline is DWM-expected signed
deviance with conventional geometry on the same corrected observed counts;
it is not a full production-pipeline comparison. Dispersion remains frozen:
HepG2 0; K562 1.6144318781609228.

Absolute rates use candidate-valid sites and, separately, positive-count
sites. Zero-count sites are never called. Relative rates require valid DWM
expectation and matching normalized observed profiles. Missing baseline
coverage cannot produce a relative pass. The table below uses positive-count
support; upper bounds are 95% Wilson bounds.

| Cell | Correction | Detector | Calls / informative sites | Upper bound | Paired increase | FPR screen |
| --- | --- | --- | --- | --- | --- | --- |
| HepG2 | DWM expected | Trained geometry | 0 / 121 | 3.08% | 0 pp | Pass |
| HepG2 | DWM expected | Spline | 28 / 121 | 31.41% | 23.14 pp | Fail |
| HepG2 | Parametric | Trained geometry | 0 / 122 | 3.05% | 0 pp | Pass |
| HepG2 | Parametric | Spline | 5 / 122 | 9.24% | 4.13 pp | Fail |
| K562 | Raw | Trained geometry | 0 / 115 | 3.23% | 0 pp | Pass |
| K562 | Raw | Spline | 0 / 115 | 3.23% | 0 pp | Pass |
| K562 | DWM expected | Trained geometry | 0 / 114 | 3.26% | 0 pp | Pass |
| K562 | DWM expected | Spline | 14 / 114 | 19.56% | 12.28 pp | Fail |
| K562 | Parametric | Trained geometry | 0 / 115 | 3.23% | 0 pp | Pass |
| K562 | Parametric | Spline | 5 / 115 | 9.78% | 4.39 pp | Fail |

All rows in this table concern CTCF. Comparable paired support is 121 in
HepG2 and 114 in K562. Calibration uses 980 and 890 common-support negatives,
respectively. ZNF362/ZNF384 lack comparable DWM control support and cannot
pass the relative screen. No difficult-family rescue is established.

Six rows pass only this FPR screen. They do not establish occupancy-leakage
safety, calibrated functional uncertainty, significant protection-shape
improvement, depth/replicate stability, or whole-workflow resource gates.
Do not select a new candidate retrospectively from the six passing rows.
The previously reported spline area under the receiver operating
characteristic curve and precision-recall curve improvements remain tuning
point estimates with this explicit safety failure attached.

## Reproduction and validation

Output directory: `benchmarks/results/fixed_shape_operator_screen_20260922_v1/`.
The phased runner is `benchmarks/scripts/screen_fixed_shape_operators.py`;
commands and phase ordering are documented in the benchmark README.

| Artifact | SHA-256 |
| --- | --- |
| Configuration | `1074877deb1f4781e7e6115c5325ff02bfb6cd8c17f86dd5def004bfdfcfe5d5` |
| Frozen thresholds | `a863657470f05f9161ad66a16013444b064088a40b18d7d36e50e6f535940965` |
| Screen manifest | `2e335e5c6eca6a037e9d746d9d7f211e4eab70cf998707215607f8dfcec8699a` |

Calibration and control execution have durable exit-0 receipts (19.89 s and
4.32 s). These cached-profile timings are not computational-gate evidence.
Independent `--verify-only` passed. A separate calculation checked all 126
per-site call/support summaries, zero-count exclusions, and informative Wilson
bounds. Eighteen focused screen, normalization, and safety tests passed in
1.61 s. The full 96-run corrected-count reassessment remains a separate job;
its final completion and verification are not asserted here.
