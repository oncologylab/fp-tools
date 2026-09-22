# Fixed-model tuning chromosome sensitivity

The CTCF geometry result is context-dependent. In HepG2, the frozen parametric
trained geometry improves discrimination on each tuning chromosome. In K562,
the pooled gain includes a loss on chr16; its descriptive chromosome-block
interval crosses zero. The small chr16 cohort makes that loss imprecise, but
does not justify excluding it. No models, thresholds, or cohorts were changed.

This analysis uses already examined chr16–18 tuning data only. It is a
post-hoc diagnostic, not calibrated functional uncertainty or the principal
untouched promotion test. All 189 task/correction/detector rows are retained
against two baselines: 378 comparisons, comprising 86 numerical non-self
comparisons, 84 self-comparisons, and 208 unavailable comparisons. All 27
ordered three-chromosome bootstrap draws are enumerated, without Monte Carlo
sampling. Undefined single-class draws make an interval unavailable.

## Separating correction and detector effects

Area under the receiver operating characteristic curve (AUROC) and relative
area under the precision-recall curve (AUPRC) gains below are paired on
identical sites. DWM means dinucleotide weight matrix. Intervals are the
2.5th and 97.5th percentiles of the 27 resamples; only three independent
chromosome blocks are available, so coverage is not established.

| Cell | Parametric detector | Baseline | AUROC gain [descriptive interval] | Relative AUPRC gain |
| --- | --- | --- | --- | --- |
| HepG2 | Conventional | DWM conventional | 0.03460 [0.02305, 0.05484] | 8.17% |
| HepG2 | Trained geometry | Parametric conventional | 0.01003 [0.00403, 0.04002] | 6.67% |
| HepG2 | Trained geometry | DWM conventional | 0.04463 [0.02756, 0.08905] | 15.38% |
| K562 | Conventional | DWM conventional | 0.02368 [0.01988, 0.03494] | 4.66% |
| K562 | Trained geometry | Parametric conventional | 0.02518 [-0.07848, 0.04814] | 3.63% |
| K562 | Trained geometry | DWM conventional | 0.04887 [-0.04622, 0.06801] | 8.46% |

The conventional-detector comparison isolates correction. The
same-correction comparison isolates detector geometry. The combined
comparison changes both. These numerical discrimination gains do not change
the underlying aggregate curves or establish stronger canonical protection.

## Chromosome support and direction

These rows compare parametric trained geometry with DWM conventional scoring.
Matching was performed across each complete tuning task, not independently
within chromosome; per-chromosome class counts therefore differ.

| Cell | Chromosome | Positive / negative sites | AUROC gain |
| --- | --- | --- | --- |
| HepG2 | chr16 | 11 / 25 | 0.15636 |
| HepG2 | chr17 | 230 / 260 | 0.04480 |
| HepG2 | chr18 | 117 / 73 | 0.02143 |
| K562 | chr16 | 6 / 22 | -0.14394 |
| K562 | chr17 | 165 / 177 | 0.06992 |
| K562 | chr18 | 82 / 54 | 0.02687 |

All spline comparisons also remain in the output. Their tuning gains must
be interpreted alongside the exact-operator false-positive failures in
`FIXED_SHAPE_OPERATOR_SCREEN_20260922.md`; this analysis does not rescue them.
The geometry screen passes do not establish cross-context stability,
difficult-family improvement, or complete safety. Nothing is promoted.

## Evidence and reproduction

Runner: `benchmarks/scripts/audit_shape_chromosome_sensitivity.py`.
Output: `benchmarks/results/shape_chromosome_sensitivity_20260922_v1/`.

- Configuration SHA-256:
  `69027a08d0d7892e394fdbf4abce0930052d3af4f6cf870d198bdba34fbba2d2`.
- Result-manifest SHA-256:
  `38765a69ca38f27116b0be16838579a56204965099da0a679da63cbc407f4908`.
- Supervised execution: exit 0, 32.16 seconds; independent `--verify-only`
  passed. This duration is not a whole-workflow resource measurement.
- Five focused tests passed in 1.95 seconds: weighted/explicit resample
  agreement, complete enumeration, missing-class handling, identical-site
  alignment, and rejection of test chromosomes and invalid inputs.
- The combined candidate, development, normalization, reassessment,
  operator-screen, comparison, and sensitivity suites passed: 36 tests in
  2.54 seconds.

Main, production defaults, original model artifacts, and locked new holdouts
remain unchanged.
