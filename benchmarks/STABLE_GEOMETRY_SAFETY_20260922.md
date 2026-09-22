# Exact stable-geometry models fail corrected-signal safety

Normalized corrected CTCF detectors fail the previously examined naked-DNA
replicate-2 screen. Do not advance these unchanged operators to untouched
promotion labels. Stable ranks and training discrimination gains did not
establish negative-control safety.

## Frozen models and thresholds

Applied the unchanged stable geometry rule to chromosomes 1–15, preserving
all 21 tasks and seven corrections. The 147 records comprise 40 selected
geometries, 70 underpowered fixed-geometry fallbacks, 30 no-protection
fixed-geometry fallbacks and seven unavailable HepG2 HMGA1 records.

Both cells' corrected CCCTC-binding factor (CTCF) detectors select geometry
(20,25,60); K562 MEF2D selects (3,8,60). HepG2 MEF2D remains underpowered.
The 47,995-byte detector JSON is only the geometry layer, not the complete
deployable bias-model artifact. Independent equivalent-solver checks verified
every selection and full serialized-score parity. Scores use fixed per-site
normalization, ten-decimal rounding and zero-read ties, not probabilities.

Reused the existing development safety design: negative-label sites on tuning
chromosomes 16–18, common valid support with positive observed and expected
totals, at least 200 negatives, and the existing 5% upper order-statistic
threshold. All thresholds were frozen before control scoring. The reference
is the rounded, unnormalized dinucleotide weight matrix (DWM) residual score
with conventional (7,12,40) geometry, not a new full production workflow.
Transferring these thresholds from ATAC negatives to naked DNA does not
guarantee nominal error control; this screen tests that transfer empirically.

Sixty combinations had sufficient calibration data. All 147 rows remain in
the final ledger: 60 screened, 80 missing matching wider-background control
profiles and seven missing models. Missing profiles are not safety passes.
All control data are previously examined replicate 2; final replicate 3 and
SK-N-SH/GM23338 labels remain unopened.

## Results

Informative sites have positive observed totals. Wilson limits are 95% upper
confidence bounds for the false-positive rate (FPR).

| Task and correction | Calls / informative sites | FPR | Wilson upper | Result |
| --- | ---: | ---: | ---: | --- |
| HepG2 CTCF, DWM | 26 / 121 | 21.49% | 29.62% | Fails |
| K562 CTCF, DWM | 24 / 114 | 21.05% | 29.42% | Fails |
| HepG2 CTCF, conditional parametric | 11 / 122 | 9.02% | 15.43% | Fails |
| K562 CTCF, conditional parametric | 12 / 115 | 10.43% | 17.36% | Fails |
| K562 MEF2D, conditional parametric | 5 / 129 | 3.88% | 8.75% | Fails bound; no paired DWM support |

Corrected CTCF paired increases over the conventional DWM reference also
exceed one percentage point. K562 MEF2D/DWM cannot be assessed because its
control artifact has no informative positive-expectation support. Seventeen
DWM task rows have zero informative candidate sites: zero calls there are
not evidence of safety.

Only three rows pass the implemented narrow screen: raw-count CTCF in both
cells (zero calls among 122/115 informative sites) and the raw-count K562
ARID3A fallback. The latter has only one comparable DWM site, making its
relative comparison extremely weak. All other non-CTCF raw-task comparisons
have zero paired DWM support. These passes do not establish bias removal,
occupancy specificity, general accuracy/shape improvement or promotion.

## Preserved failed attempt

The first calibration exited 1 before thresholds or control scoring: the new
wrapper resolved a relative reference against the research worktree instead
of the original data root. That same-named track failed the expected checksum.
Both original tracks match their frozen hashes at the declared data root;
no original data or expected hash changed.

Preserved the failed v1 configuration, runner and durable receipt. The v2
runner explicitly freezes the original read-only data root, with a regression
test for calibration and control loading. This changes path resolution, not
scientific inputs, chromosome membership or the threshold rule.

## Verification

- Fifteen focused tests passed in 1.50 seconds: serialization, zero reads,
  stable ranks, baseline geometry, unavailable controls, safety support/calls
  and data-root resolution.
- All 147 full-training choices and complete serialized scores independently checked.
- All 147 calibration counts and thresholds independently reproduced using
  sorted order statistics and equivalent normalization.
- All 147 control rows independently checked: coordinates, scores, validity,
  call masks, denominators, Wilson bounds using SciPy, paired differences and
  gate decisions. An optional statsmodels import was unavailable; the audit
  used installed SciPy without changing the environment.
- Separate artifact, code, runtime and threshold-identity verification passed.
- Calibration and evaluation exited zero in 24.923284 and 11.048272 seconds.
  These are diagnostic elapsed times, not resource-gate measurements.

## Artifacts and next work

Ignored directories under `benchmarks/results/`:
`stable_geometry_models_20260922_v1/`,
`stable_geometry_screen_20260922_v2/`, and
`stable_geometry_screen_v2_jobs_20260922/`.
The failed v1 screen and receipt directories are preserved.

SHA-256:

```text
models configuration    a856496e40d89790d1c370bad68309e33abc57c27417179fd297404230b1d73e
models.freeze.json      7cfea9cf4c0b9ed7d3bf080aea416d5b88fdfaabc875a0ee0ee96df007a5495d
detectors.json          36360391266f1b5ffb2aaaacad87781b830ff49aa9c87ee424893f527d77aa85
screen v2 configuration 5ac308e27dec74a632afdd37d903f92c6f77d4af6c2397b194ed7faa40092352
thresholds.freeze.json  dd5031ca534ac28096655d11c424f089dd08f2b8ec6e19c6c8513d56fd49ba24
screen_manifest.json    588670f1a9acf789c7ebc5796f1662554f41e8a7254379a0c4995b2a14147a39
```

Next diagnose normalization/depth behavior on permitted training and
replicate-1 negative-control material, and audit missing expected coverage in
the archived DWM controls. Do not tune thresholds or choose a correction to
pass observed replicate-2 results. Preserve failed operators and missing-
support cases. Any successor needs a separately frozen design and exact-
operator safety checks. Full-workflow output parity, resource qualification
and the remaining scientific gates are still incomplete.
