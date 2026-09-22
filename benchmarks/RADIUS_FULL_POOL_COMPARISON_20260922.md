# Seven fixed methods on directly balanced training cohorts

The broader-count-balanced comparison supports a narrow HepG2 zinc-finger
discrimination signal, but not a general footprint improvement. The unchanged
radius-500 policy does not establish the required aggregate-shape gain or
improvement across difficult families and cell contexts. Do not advance it
to untouched promotion labels on the strength of this study.

## Frozen comparison

All seven arms were recomputed on the 4,053 pairs from the direct three-
covariate matching study: raw counts, dinucleotide weight matrix (DWM)
expected-profile residuals, existing conditional parametric residuals,
global parametric radii 100/200/500, and uniform-sequence radius 500.
The same fixed dispersion, conventional detector geometry, normalization
fallback 1 and zero-read tie rule were retained. No model was refitted or
selected. Imported research helpers, direct code, inputs, cohort identities
and runtime were recorded in the execution freeze before scoring.

All 21 tasks and 147 task/method rows are present: 140 have recomputed scores
and curves; seven explicitly lack a balanced HepG2 HMGA1 cohort. Ten tasks
meet the 200-pair requirement; ten are underpowered. Only chromosomes 1–15
were selected and evaluated. Source tables include already examined tuning
rows but these were not selected or scored. Previous cohorts, metrics and
curves were not used as substitutes for this exact-cohort evaluation.

## Results

Radius-500 changes relative to the matched DWM residual comparator:

| Task | Pairs | AUROC gain | Relative AUPRC gain | Functional-separation gain |
| --- | ---: | ---: | ---: | ---: |
| HepG2 ZNF362 | 291 | +0.073051 | +10.74% | +5.20% |
| HepG2 ZNF384 | 303 | +0.061181 | +11.27% | +6.25% |
| K562 ZNF384 | 238 | +0.017019 | +0.66% | +1.53% |
| HepG2 CTCF | 300 | +0.027567 | +7.03% | −0.70% |
| K562 CTCF | 245 | +0.018209 | +10.74% | +12.00% |

AUROC is area under the receiver operating characteristic curve; AUPRC is
area under the precision-recall curve. CTCF is CCCTC-binding factor.
These are training point estimates, not independent validation or principal
chromosome-block confidence intervals. DWM here is the matched expected-
profile residual comparator, not a newly executed full production workflow.

Radius 500 has canonical central-negative point contrasts in the rows above,
but their standardized aggregate-separation changes fall well below 25%.
CTCF discrimination remains below the previous conditional parametric arm
in both cells. HepG2 ZNF362/ZNF384 discrimination survives direct broad-count
matching, while the K562 ZNF384 gain is substantially weaker than in the
initial two-covariate cohort. Matching defines a different population; this
does not isolate a causal explanation for the difference between analyses.

FOXA1 still lacks canonical protection. Other numerical gains often lack
improved shape or sufficient power. MEF2A has 190 HepG2 and 169 K562 pairs,
below the frozen threshold; its numerical trends do not satisfy a second
difficult-family gate. No full-denominator promotion mean is claimed when
one task is unavailable and others are underpowered.

All 21 aggregate panels were visually reviewed across three pages, including
the explicit unavailable-task panel. Several zinc-finger contrasts are shallow
and noisy; the numerical canonical-sign flag alone is not a stable footprint
claim. No calibrated confidence bands or biological safety are established.

## Verification and artifacts

Eight relevant tests passed in 2.45 seconds, including all-task synthetic
evaluation, unavailable rows, zero-read ties, scalar/halo parity and independent
normalization agreement. Supervised execution exited zero in 40.58 seconds.
This cached-profile comparison is not a full-workflow resource benchmark.
Separate configuration/code/input/output verification passed.

The independent reviewer recalculated scores from saved per-site panels,
AUROC from rank sums, and AUPRC from saved score/label pairs. It recalculated
aggregate normalization with a separate least-squares formulation and checked
both group means, central/shoulder contrasts and canonical directions for
all 140 scored rows. It checked the seven unavailable rows and rendered every
task. The generated review manifest remains unchanged; subsequent visual
inspection of its three images is documented here.

Outputs:

- `results/radius_full_pool_comparison_20260922_v1/`: metrics, curves, scores,
  exact site table, per-site profile panels and frozen manifest.
- `results/radius_full_pool_review_20260922_v1/`: independent verification
  manifest, three-page aggregate PDF and three rendered page images.
- `results/radius_full_pool_comparison_jobs_20260922/execute_attempt1/`:
  immutable command specification and launch/exit evidence.

SHA-256:

- Configuration: `d186387c24c291f8cf3a9aab21efa56f21ef045a14883820e9e54fc11ec5fa9b`
- Comparison manifest: `7e8c98ec0bb440a58aad9570c7de1fdb23ca97b51e1a790806ce0c1f85fdec84`
- Review manifest: `f3d1e7a004dcde8b3e6374cc9505a5b8ee048d439195919186f6e019bab223b2`
- Aggregate PDF: `4f9f4884ddf7f5c3832e84435fbde13bd9b28e7f9c9351bde83e5ee7aed285e8`

The next useful diagnostic is chromosome-block uncertainty for the balanced
training gains and signed shape contrasts, retaining every task and method.
Such intervals remain development evidence and cannot satisfy an untouched
promotion gate. Further detector changes must be developed on permitted
training data; exact-operator safety, depth/replicate stability, calibrated
functional uncertainty and full-workflow resource gates remain open.
No new holdout or final negative-control labels were opened. Main and
production scientific defaults remain unchanged.
