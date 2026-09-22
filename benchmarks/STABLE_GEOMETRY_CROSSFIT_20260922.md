# Numerically stable training geometry comparison

The numerical rank failure in the preceding diagnostic is resolved on the
complete current dataset. This establishes reproducible scores and selections,
not a general footprint improvement or permission to open promotion labels.

## Prespecified numerical successor

Preserved the previous run and its pinned calculation code. The successor
uses the same cohorts, chromosome/pair exclusions, eleven geometries, fitting
support requirement, seven correction arms and three detectors. Only the
numerical policy changed:

- Round geometry scores to ten decimal places before fitting ranks and
  evaluation. The rounding error is at most about 5e-11 score units. This is
  small compared with meaningful normalized profile contrasts and larger
  than the arithmetic differences observed in the failed diagnostic.
- Compare fitting area under the receiver operating characteristic curve
  (AUROC) through exact integer rank-sum numerators. Exact candidate ties
  retain the first geometry; integration roundoff cannot select a different
  geometry from an exact mathematical tie.
- Require the negative center and each shoulder-minus-center contrast to
  exceed a 1e-10 numerical margin. Near-zero arithmetic residuals cannot
  count as canonical protection.

The policy was tested on synthetic sparse and flat profiles and frozen before
the successor evaluated biological data. It was not searched for improved
biological performance. It is a numerical score representation, not a binding
threshold, probability calibration or uncertainty model. Decimal rounding
does not guarantee equivalence for every possible future input at a rounding
boundary; future candidates and datasets must retain numerical parity checks.

## Complete verification

All 2,205 fits, 6,615 fold rows (6,216 evaluable), saved site scores and 441
summary rows passed independent checking. All fitting eligibility decisions,
candidate rank metrics and selected geometries agree under an equivalent
least-squares normalization calculation. All archived held scores agree
exactly under both calculations. Zero fold metrics changed beyond arithmetic
roundoff: maximum AUROC recomputation difference 1.11e-16 and area under the
precision-recall curve (AUPRC) difference zero. The prior run had 287 affected
fold rows. Ten selected geometries changed relative to that unstable run,
all in the raw-count correction arm; the earlier artifacts remain intact.

Ten focused tests passed in 3.15 seconds, including sparse/flat profiles,
individual versus batch scoring, solver parity, tied ranks, finite-input and
missing-class validation, protection margins, exclusion of matched-pair
partners and held-label isolation in the original unchanged split logic.
Separate artifact/code/runtime verification passed. The durable execution
returned zero in 34.751810 seconds; this is diagnostic elapsed time, not a
full-workflow central processing unit or memory qualification.

## Scientific findings remain limited

Mean within-chromosome AUROC, comparing selected geometry against fixed
geometry after the same normalization and correction:

| Task and correction | Fixed | Selected | Change |
| --- | ---: | ---: | ---: |
| HepG2 CCCTC-binding factor (CTCF), dinucleotide weight matrix (DWM) | 0.589444 | 0.656105 | +0.066661 |
| K562 CTCF, DWM | 0.630200 | 0.657641 | +0.027441 |
| K562 MEF2D, DWM | 0.506479 | 0.586154 | +0.079675 |
| HepG2 ZNF362, radius-500 parametric | 0.581932 | 0.576334 | −0.005598 |
| HepG2 ZNF384, radius-500 parametric | 0.584812 | 0.577601 | −0.007211 |
| K562 ZNF384, radius-500 parametric | 0.541540 | 0.530952 | −0.010588 |

K562 MEF2D's DWM mean AUPRC increases from 0.570264 to 0.615092
(7.86% relative), and its selected-geometry protection direction holds in
10/15 chromosomes. HepG2 MEF2D remains underpowered and uses the fixed
fallback. This does not establish improvement across cell contexts. Numerical
gains in HepG2 FOXA1 do not establish the fixed-center canonical protection
gate; changing geometry cannot itself change the underlying aggregate curve.

All 21 tasks remain explicit, including unavailable HepG2 HMGA1. There are
543 protective selections, 1,225 underpowered fallbacks and 437 no-protection
fallbacks. These are training diagnostics with full-training cohort matching,
not fully nested or untouched validation. Equal-chromosome means differ from
pooled-site metrics in earlier reports and must not be substituted for them.
No significance, simultaneous uncertainty, negative-control safety,
depth/seed/replicate stability or computational promotion gate is newly passed.

## Artifacts and next action

Ignored outputs: `benchmarks/results/stable_geometry_crossfit_20260922_v1/`.
Execution evidence: `benchmarks/results/stable_geometry_crossfit_jobs_20260922/`.
Independent full-dataset audit:

```bash
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /home/exouser/projects/fp-tools/.venv/bin/python \
  benchmarks/scripts/verify_stable_geometry_crossfit.py \
  --outdir benchmarks/results/stable_geometry_crossfit_20260922_v1
```

SHA-256:

```text
configuration.freeze.json bbfe92b52e1e84fd05027a224933423839bdb1fc8ade9b696757c77afa832258
crossfit_manifest.json    500dcb12af488b3df68c0f9025d4c945d1ed6a2b05106dae10748146aec969bc
```

Next assess chromosome-wise stability and signed protection for the exact
stable detector, retaining all tasks and correction controls. A deployment
candidate would additionally need a frozen full-training model, exact-operator
negative-control testing and an executable end-to-end output/resource
contract. Do not advance the unchanged zinc-finger route to untouched labels.
