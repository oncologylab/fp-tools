# Fixed background-radius comparison: training evidence only

The 500-base-pair background yields a promising training diagnostic for
ZNF362/ZNF384, but does not pass the promotion gates. Further development must
test whether broader accessibility, rather than localized protection, explains
the gain. No policy was selected or evaluated on new holdout labels.

## Design

The frozen comparison used the same 3,949 matched pairs across all 21 tasks
and six arms: raw counts, dinucleotide weight matrix (DWM) expected-profile
residuals, existing conditional parametric residuals, and global parametric
background radii 100, 200 and 500 base pairs. All three global arms use the
same inner radius 50, half-count prior, fixed deviance dispersion and
conventional detector geometry. Only the outer radius changes between them.
The DWM arm is a matched expected-profile residual comparator, not a newly
executed complete production workflow.

Only chromosomes 1–15 were selected and scored. The shared source reader loads
existing tables that also contain previously examined tuning rows, but no
tuning rows were selected or evaluated. Ten tasks have at least 200 pairs;
eleven are underpowered. None was omitted. All original pairs survived common
profile support. Absolute standardized mean differences remain at most 0.1
for motif score and the original 201-position observed-count covariate.
This does not establish balance for the newly used wider background counts.

## Numerical findings

Area under the receiver operating characteristic curve (AUROC) and relative
area under the precision-recall curve (AUPRC) gain for radius 500 versus DWM:

| Cell/task | DWM AUROC | Radius-500 AUROC | Relative AUPRC gain | Canonical contrast |
| --- | ---: | ---: | ---: | --- |
| HepG2 ZNF362 | 0.511536 | 0.608994 | 16.22% | Yes |
| HepG2 ZNF384 | 0.521373 | 0.638541 | 17.02% | Yes |
| K562 ZNF384 | 0.508660 | 0.615442 | 15.68% | Yes |
| HepG2 CTCF | 0.594552 | 0.600334 | 5.19% | Yes, weak central contrast |
| K562 CTCF | 0.618622 | 0.652347 | 14.54% | Yes |

CCCTC-binding factor (CTCF) still has lower discrimination than the previous
conditional parametric residual in both cells (0.606990 and 0.669707).
Radius 200 fails the central-negative contrast check for HepG2 CTCF.
FOXA1 remains inconsistent with canonical protection; MYC and many other
tasks remain weak or underpowered. All task results are preserved.

Across all 21 tasks, the descriptive mean AUROC is 0.522274 for radius 500
versus 0.497995 for DWM, a gain of approximately 0.02428. This includes
underpowered tasks and is not an independent promotion estimate.
Standardized functional-separation increases for the three zinc-finger rows
above are only 7.86%, 22.68% and 10.45%, respectively, below the required 25%.
These related zinc-finger tasks do not establish improvement in two distinct
difficult transcription-factor families. No significance interval, biological
safety, depth/replicate stability or calibrated uncertainty is established.

Every aggregate panel was rendered and visually reviewed across three pages.
The contrasts show a central dip for the zinc-finger tasks at radius 500,
but the curves remain noisy. The point-estimate shape criterion does not
establish a stable biological footprint. Aggregate curves use fixed detrending
and normalization with fallback 1 and clipping at 12; no confidence bands
are implied.

## Verification

The supervised execution exited zero in 34.06 seconds. This is a cached-profile
comparison, not an end-to-end runtime benchmark. Separate frozen-artifact
verification passed. An independent review recomputed all 126 AUROC/AUPRC
rows from saved scores and verified identical coordinates, labels, class
counts and training chromosomes across all six arms.

The prefix-sum radius-100 implementation agrees with the previous scores to
maximum absolute difference 1.547e-13; no task's AUROC or AUPRC changed.
Fifteen relevant tests passed in 2.72 seconds, including a complete synthetic
21-task execution, zero-read ties, pair exclusions, coordinate identities,
scalar/halo numerical parity and independent diagnostic calculations.

Artifacts are under `results/radius_training_comparison_20260922_v1/` and
`results/radius_training_review_20260922_v1/`. The automatic review manifest
leaves visual inspection false; subsequent inspection of all three rendered
pages is recorded here without rewriting the generated manifest.

SHA-256:

- Configuration: `c448b0a18ef1e4ce092eebdd4d6fe9a5861a5f389f13bdddc92c7f4c6db68beb`
- Result manifest: `f168c4c08210179b097e1b3135cdf3a2602036b602a2e4cea79c747c8df1168f`
- Review manifest: `2dd9bb4e7d92eb307b043224918716a67355eb15bce583267093d4911ba6d26b`
- Aggregate PDF: `fc0bd1f2ec7bab591d98deecf03e1d5fc7abc9f7efebc8fd6e69becb5c1f2849`

## Next decision

Before selecting a policy, audit broader-flank accessibility balance and
dependence of the score on those counts on training chromosomes. Preserve
the current cohorts and report any additional matched/stratified analysis as
secondary, including lost power. A uniform-sequence background comparator
can separate the value of the bias model from the new spatial baseline; it
requires a separate freeze before evaluation. Do not use test chromosomes,
new transfer labels or final naked-DNA controls to develop these adjustments.

No method was promoted, uploaded, or integrated into main. Production defaults
and manuscript files remain unchanged.
