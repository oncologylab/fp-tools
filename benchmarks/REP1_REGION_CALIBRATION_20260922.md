# Replicate-1 region calibration: development evidence

This audit supports sequence-related central cleavage deficits for the examined
zinc-finger tasks, but does not establish general calibration or safe footprint
calling. Nothing is eligible for promotion. Main and the production dinucleotide
weight matrix (DWM) method remain unchanged.

## Frozen design

Use the existing freshly extracted naked-DNA replicate-1 counts and sequence-only
parametric propensities at all 4,200 archived training-chromosome sites. Apply the
recorded blacklist exclusions without dropping original rows. Fix geometry at
(7, 12, 40) for every task: 15 central bases and two 29-base shoulders. Condition
on each strand's total across these regions. Compare parametric probabilities
with uniform probabilities (15, 29, 29)/73. No models or thresholds are fitted.

Retain all 21 tasks, including HepG2 HMGA1, for which no occupancy detector was
available. Evaluate strands separately and pooled. The outputs contain 50,400
site-region rows, 378 calibration summaries, 63 paired likelihood comparisons,
and 5,670 chromosome-block rows. Resample the same 15 training chromosomes with
2,000 shared multinomial draws, seed 2026; retain zero-count draws as invalid
rather than redrawing. Report descriptive percentile intervals.

Predictive count intervals convolve independent site-strand binomial counts.
They assume independent cuts and correct propensities; their biological coverage
has not been validated. Neither these intervals nor the block comparisons are
adjusted for multiple comparisons. This is examined development control data,
not independent candidate safety or final negative-control evidence.

## Results

Pooled central counts show that uniform cleavage is a poor expectation for the
three examined zinc-finger tasks:

| Task | Total region cuts | Central observed | Parametric expected | Uniform expected |
| --- | ---: | ---: | ---: | ---: |
| HepG2 ZNF362 | 87 | 3 | 3.75 | 17.88 |
| HepG2 ZNF384 | 77 | 2 | 3.95 | 15.82 |
| K562 ZNF384 | 78 | 6 | 3.65 | 16.03 |

All three observed counts fall inside the parametric predictive interval [1, 8].
Paired parametric log-likelihood gains per cut over uniform are respectively
0.16650 (block interval 0.05494–0.23910), 0.15629 (0.04006–0.24554), and
0.09285 (−0.03112–0.19413). The K562 comparison remains uncertain.
These results support a sequence-related cleavage deficit, not transcription
factor protection or an absence of information in accessibility data.

Calibration is not universal. Four pooled central counts fall outside the
model-based predictive intervals:

| Task | Observed | Parametric expected | Predictive interval |
| --- | ---: | ---: | --- |
| HepG2 HMGA1 | 19 | 10.53 | 5–17 |
| K562 CTCF | 13 | 21.78 | 14–30 |
| K562 FOXA1 | 28 | 16.84 | 10–24 |
| K562 MEF2D | 21 | 11.94 | 6–18 |

These are exploratory discrepancies, not multiplicity-adjusted failures or
proof of occupancy leakage. Central overprediction for K562 CTCF particularly
warrants deeper control calibration before interpreting a protection score.
The earlier zero-call screen barely challenged the scoring threshold because
only two eligible parametric control sites could reach it at these counts.

Across all tasks, 18/21 pooled likelihood gains are positive; only five
descriptive block intervals exclude zero on the positive side. This is not
evidence of general calibration, biological uncertainty coverage, or promotion.

## Reproducibility and verification

Ignored artifacts: `benchmarks/results/rep1_region_calibration_20260922_v1/`.

- Configuration SHA-256:
  `0dd9b495a32d71452f27b59dcb10783511944063666c04a6d3524e5afcca0d44`.
- Result manifest SHA-256:
  `ddff71108c168140e5d66f3d239fc27e507074e4c60e6c7e739f7aceb839b79e`.
- Immutable job specification SHA-256:
  `db22d9d5c3f4c937c5c582ef49c5ebb4e5790a3197d3d7b2562dd6f841527f68`.
- Durable job receipt in `rep1_region_calibration_jobs_20260922/execute_attempt1/`:
  exit 0, elapsed 6.588232700014487 seconds. This is diagnostic execution time,
  not complete-workflow central processing unit (CPU) or memory qualification.
- Seven focused tests passed in 4.10 seconds across region calibration,
  conditional protection scoring, and score-ceiling tests.
- Frozen-artifact verification passed. Independent reconstruction verified all
  50,400 site-region rows, 378 calibration/predictive intervals, 63 likelihood
  comparisons, and 5,670 blocks. It used explicit exponentiation for probabilities
  and one-Bernoulli-at-a-time predictive convolution. No new independent read
  extraction was performed in this stage.

## Next permitted work

Audit higher-depth replicate-1 input provenance before freezing another study.
Larger alignment intermediates and trimmed reads exist locally, but existence
does not establish a completed, checksum-valid final alignment. Do not substitute
an intermediate or change filtering silently. Any new extraction and derived
files must stay in this research worktree. Keep fixed site membership, correction
and detector effects separate, and preserve all task denominators.

Do not tune on previously examined test chromosomes or transfer labels. Keep
SK-N-SH and GM23338 labels and final naked-DNA replicate 3 unopened. End-to-end
resource gates, difficult-family protection, robust sensitivity and independently
qualified safety remain incomplete.
