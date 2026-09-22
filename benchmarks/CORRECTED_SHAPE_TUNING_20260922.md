# Corrected-data training and tuning comparison

Models were fitted on chr1–15 and frozen before scoring chr16–18. All 21
cell/transcription-factor tasks remain in the output (189 correction/detector
rows). This is development evidence from replicate 1, 10 million fragments,
seed 2026, not a test, transfer, depth-stability or promotion result.

The comparison uses the same balanced pairs for raw counts, dinucleotide
weight matrix (DWM) expected-signal residuals and frozen parametric residuals.
The three detector arms are conventional geometry, training-selected geometry,
and the existing spline detector with family/global template shrinkage.
Training central enrichment and insufficient power prevent candidate fitting;
their tasks remain in the ledger. Eleven model artifacts were frozen, with
exact prediction parity after serialization. No tuning-driven refit or policy
selection occurred.

## CTCF tuning results

AUROC is area under the receiver operating characteristic curve; AUPRC is area
under the precision-recall curve. These point estimates have no new confidence
intervals or significance claim.

| Cell | Correction | Detector | AUROC | AUPRC |
| --- | --- | --- | ---: | ---: |
| HepG2 | DWM expected | Conventional | 0.6031 | 0.5929 |
| HepG2 | DWM expected | Spline | 0.6382 | 0.6696 |
| HepG2 | Parametric | Conventional | 0.6377 | 0.6413 |
| HepG2 | Parametric | Training-selected geometry | 0.6478 | 0.6841 |
| HepG2 | Parametric | Spline | 0.7063 | 0.6960 |
| K562 | DWM expected | Conventional | 0.6209 | 0.6168 |
| K562 | DWM expected | Spline | 0.6706 | 0.6885 |
| K562 | Parametric | Conventional | 0.6445 | 0.6456 |
| K562 | Parametric | Training-selected geometry | 0.6697 | 0.6690 |
| K562 | Parametric | Spline | 0.6764 | 0.6926 |

There are 358 matched tuning pairs for HepG2 CTCF and 253 for K562 CTCF.
Within the parametric correction, spline ranking gains approximately 0.069 and
0.032 AUROC over conventional geometry. Changing the detector does not change
the input aggregate profiles, so this is not evidence of 25% stronger aggregate
separation. The signed central contrasts of parametric profiles are negative
in both contexts, but the simple point-estimate shape screen does not establish
calibrated functional uncertainty or significant canonical protection.

## Difficult families and limits

The experiment does not broadly rescue difficult families. HepG2 ZNF384's DWM
spline AUROC is approximately 0.486; K562's is 0.456. K562 parametric ZNF384
spline AUROC is 0.594, but its aggregate center is enriched, so this is not a
footprint improvement. HepG2 FOXA1 and ZNF362 also lack protective parametric
aggregate centers. Small positive discrimination changes cannot override that
failure. All other task/correction results and unavailable candidates remain
in `tuning_metrics.tsv`; no favorable-task mean is substituted for the locked
task denominator.

Seven tasks have adequately powered balanced cohorts in both splits. That
does not make all their fitted candidates protective. Many other tasks remain
underpowered, and training shape screens reject additional candidates. The
comparison is against the same residual-scoring operator, not a claim that
the full production DWM pipeline has been outperformed.

No independent negative-control safety is established for these newly learned
geometries or spline scores. Earlier direct-residual safety does not transfer
to them. No functional uncertainty calibration, complete-workflow resource
measurement, new holdout evaluation or general replacement has passed.

## Artifacts and verification

Output directory: `benchmarks/results/corrected_shape_development_20260922_v1/`.

- Configuration SHA-256:
  `0dce35c8f1c78c7042f7b01d226d674c214c42dda7f788f66e5288f60fa01302`.
- Model-freeze SHA-256:
  `a9780ac8006fabc1c945803c0ae0149a4bbcf69be73f95470e546eb48a28b675`.
- Verified tuning-manifest SHA-256:
  `74af0ce9cd66bc9d06af81d3edc1f95f04f9ec7c84da703c6de00e645c97c9ec`.
- Training exited 0 in 5.68 seconds; frozen models verified without refitting.
  Tuning exited 0 in 5.89 seconds; separate `--verify-only` passed. These cached
  profile-stage durations do not establish end-to-end runtime or memory gates.
- Nine focused development/candidate tests passed in 1.58 seconds, including
  exact serialized prediction parity and refusal to score before model freeze.

Per-site scores and per-correction aggregate curve data are checksummed.
Rendered aggregate plots and their visual review remain outstanding. Main,
the published manuscript and locked holdout labels remain unchanged.
