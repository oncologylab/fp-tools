# Corrected-data training and tuning comparison

**Subsequent normalization qualification:** a training-only inference audit
found that legacy raw-spline predictions depended on prediction batch contents
when a site had no outer-flank variation. The frozen results below remain the
recorded full-task-batch calculation; they do not establish batch-independent
raw-spline inference. A separate corrected model copy now freezes the fallback
from training data, without changing learned weights. Its tuning reassessment
and new-operator safety screen remain pending. DWM and parametric training
scores had no material batch dependence in the audited sites; this is not a
substitute for their corrected-artifact tuning and safety verification.

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
All 21 tasks were subsequently rendered across seven pages in
`benchmarks/results/corrected_shape_tuning_plots_20260922_v1/all_task_aggregates.pdf`.
All seven page images were visually reviewed: labels and curves are readable,
underpowered tasks remain explicit, and there is no general visual rescue of
difficult families. CTCF retains the clearest central depletion and shoulders;
the parametric curves do not establish stronger protection than DWM. Axes share
limits within each task, and the central/shoulder windows are shaded. These
are unsmoothed point-estimate contrasts with no confidence bands.

The PDF and seven images match their recorded checksums. PDF SHA-256:
`27a92fbb6bfffb066b28c944c17e869963478caa67ce76ec029d16acd52ed0dd`.
Rendering manifest SHA-256:
`18b83cdd449db1d832972077f7e6d4e6f6a412e0ffdc564a395bb7b3b3dd2ec1`.
Three renderer tests passed in 2.48 seconds, and twenty documentation contracts
passed in 0.42 seconds. Main, the published manuscript and locked holdout labels
remain unchanged.

## Training-frozen normalization correction

The original raw-spline models changed scores between full-task batches and
individual-site inference at 59/512 HepG2 ZNF362, 52/596 HepG2 ZNF384 and
71/560 K562 CTCF training sites. Maximum absolute differences were approximately
1.613, 2.723 and 3.841 score units. This is a reproducibility defect, not new
biological evidence. The other eight model/correction combinations had only
floating-point differences in this check.

`freeze_shape_normalization.py` copies each artifact and adds a fallback scale
computed from its original chr1–15 training profiles. No model is refitted;
all original learned arrays and artifacts are preserved. Every corrected model
reproduces its original full-training-batch scores and agrees between batch
and singleton inference within 1e-10 tolerance. The real stage exited 0 in
20.45 seconds and separate verification passed. Ten related tests passed in
1.62 seconds.

Output: `benchmarks/results/corrected_shape_normalization_20260922_v1/`.
Configuration SHA-256:
`106495fea837880a8b31aa6fa02432950e31ca6c487b3cd08718e67599543b1b`.
Corrected model-freeze SHA-256:
`52986cf07c0861f39ef2a96966056ba7867b4d454ee45c3d475c9e7da2aa1f84`.
This correction establishes deterministic inference for the audited models;
it does not establish improved accuracy, calibrated uncertainty or safety.

### Corrected-artifact tuning reassessment

The frozen model copies were subsequently scored on the identical tuning
sites, without refitting or candidate selection. All 189 metric rows and
21 tasks remain. DWM and parametric spline scores have no material change
(maximum numerical differences below 1e-15), and their metrics are unchanged.
Conventional and training-selected geometry scores are unchanged.

Raw-spline scores change at 47 HepG2 ZNF362, 54 HepG2 ZNF384 and 70 K562 CTCF
tuning sites. AUROC changes from 0.564280 to 0.564297, 0.549280 to 0.549966,
and 0.639941 to 0.639254, respectively. These are implementation-reassessment
differences, not new candidate improvements. All corrected spline predictions
agree between singleton and full-task inference within 1e-10 tolerance.
The aggregate curve artifacts are reused as originally defined per-task
normalized input profiles; detector normalization does not retroactively
replace those descriptive curves.

Output: `benchmarks/results/fixed_shape_tuning_20260922_v1/`. Execution exited
0 in 12.24 seconds and separate verification passed. Nine related tests passed
in 1.47 seconds. Configuration SHA-256:
`16b6b1ea9001d2e0c530bded93844fb9b8073509648ce3671fc376f9d38c1c7a`.
Result SHA-256:
`24096520aaab92947a0351a9208b32575392b5512d0a5aaab6ddf5ac71b5b65d`.

The new-operator false-positive screen is specified in
`benchmarks/manifests/fixed_shape_operator_safety_v1.json`. It includes every
fitted model and retains missing candidates in the 21-task ledger. Thresholds
must use common positive baseline expectation and be frozen before control
scoring. Its execution configuration and results are still pending; no new
operator safety pass is claimed.
