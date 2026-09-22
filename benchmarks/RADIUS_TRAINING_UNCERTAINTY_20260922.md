# Training chromosome-block sensitivity

The wider-background candidate retains a HepG2 zinc-finger discrimination
signal, but central protection remains uncertain. This study does not justify
advancing the unchanged candidate to untouched promotion labels.

## Frozen design

Reused the checksum-verified full-pool comparison: 4,053 matched pairs,
21 transcription-factor tasks and seven fixed methods. Only training
chromosomes 1–15 were evaluated. Drew 15 chromosomes with replacement
2,000 times per cell, with seeds 2026 (HepG2) and 2027 (K562). Each cell's
multiplicities were shared across tasks and methods. Cohorts and parameters
were fixed; there was no rematching or candidate selection. Chromosome
reweighting can change original pair and covariate balance.

Reported pointwise percentile intervals for discrimination gains relative
to dinucleotide weight matrix (DWM) residuals and normalized positive-minus-
negative profile contrasts. Center: offsets −5 through 5; shoulders: −40
through −15 and 15 through 40. Canonical direction requires a negative
center below both shoulders. Its resampling fraction is not a binding
probability or a calibrated posterior probability.

## Results

Area under the receiver operating characteristic curve (AUROC) gain intervals
and relative area under the precision-recall curve (AUPRC) gain intervals
for the radius-500 arm:

| Task | AUROC gain interval | Relative AUPRC gain interval | Canonical fraction |
| --- | ---: | ---: | ---: |
| HepG2 ZNF362 | [0.041923, 0.101381] | [4.56%, 17.29%] | 74.50% |
| HepG2 ZNF384 | [0.012317, 0.105471] | [2.19%, 21.44%] | 76.30% |
| K562 ZNF384 | [−0.039084, 0.077325] | [−9.33%, 12.78%] | 89.00% |
| HepG2 CCCTC-binding factor (CTCF) | [−0.016977, 0.077913] | [−1.39%, 15.74%] | 99.90% |
| K562 CTCF | [−0.024921, 0.062816] | [0.88%, 21.33%] | 96.05% |

HepG2 ZNF362 and ZNF384 center-contrast intervals both cross zero
([−0.095586, 0.047796] and [−0.065925, 0.036460]). Their numerical
discrimination gains cannot be called established footprint improvements.
The earlier exact-cohort functional-separation gains remain below 25%.
HepG2 MEF2D also has positive discrimination intervals, but only 164 pairs
and an uncertain center contrast; it remains underpowered.

All 147 rows were retained. Each of 140 available rows has 2,000 valid
resamples. Seven HepG2 HMGA1 rows have no balanced cohort and explicitly
retain 2,000 unavailable resamples each. No draws were retried.

These are exploratory training diagnostics, without multiplicity correction,
not principal holdout intervals or calibrated simultaneous functional bands.
No new safety, depth/replicate stability, or end-to-end resource gate was
passed. Untouched labels and final naked-DNA controls remain unopened.

## Verification and artifacts

- Durable execution completed with exit 0 in 14.308574 seconds; this is
  diagnostic execution time, not a resource-gate measurement.
- Separate checksum/code/runtime artifact verification passed.
- Independent recomputation verified every interval, valid-resample count
  and canonical fraction across all 147 rows.
- Explicit chromosome replication, independent rank metrics and least-squares
  profile normalization matched 420 saved score/shape draws: three fixed
  draws for each of 140 available task/method rows.
- Four focused tests passed in 2.52 seconds, covering weighted tied ranks,
  missing classes, shape weighting and source evaluation behavior.

Ignored artifacts: `benchmarks/results/radius_training_bootstrap_20260922_v1/`.
Durable receipts: `benchmarks/results/radius_training_bootstrap_jobs_20260922/`.

SHA-256 checksums:

```text
configuration.freeze.json b0131538a58a59c34cfda32622123a0fb6fe086e4462abeda85b3d55834f6b39
bootstrap_manifest.json   676993dd079e045a4d197c4b9160ae25acc7e750aa06201786ea32ee390e80ca
```

Next development priority: explain weak signed protection on the balanced
training cohorts before expanding model complexity. Any new detector must
have a frozen training specification, identical-site correction controls and
exact-operator negative-control testing. Keep the existing resource protocol
separate; it still needs executable full-workflow parity and isolated timing.
