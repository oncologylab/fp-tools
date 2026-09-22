# Fixed unit-floor geometry comparison — 22 September 2026

A fixed unit scale floor avoids amplifying sparse profiles but does not provide
a broad footprint improvement. Retain it as a failed general-replacement
candidate and development diagnostic. No safety thresholds were recalibrated,
and no control or new holdout was evaluated in this experiment.

## Frozen design

The candidate detrends outer flanks, divides by `max(RMS, 1)` where RMS means
root mean square, clips at plus/minus 12, and rounds geometry scores to ten
decimal places. The floor is identical for every task and was not fitted.
All 8,106 matched training sites, 21 tasks and seven correction arms remain.
The same 11-geometry grid, central-protection eligibility rule, minimum 200
sites per class, and explicit fallback rules apply.

Each held chromosome among chromosomes 1–15 is excluded from fitting together
with its matched pair mates. Matching itself was previously fixed on the full
training cohort, so this is not a fully nested matching study. Discrimination
is measured within each held chromosome and averaged equally over available
chromosomes; differently fitted scores are never pooled. The unavailable
HepG2 HMGA1 cohort remains explicit. Archived comparator stages were reused.

## Discrimination and shape

AUROC means area under the receiver operating characteristic curve; AUPRC
means area under the precision-recall curve. Across the 20 available tasks,
the conditional-parametric candidate's descriptive mean AUROC gain over the
unnormalized conventional dinucleotide weight matrix (DWM) detector is 0.01558;
the mean per-task relative AUPRC gain is 3.30%. Neither reaches the respective
0.03 or 10% target, and these training results cannot establish promotion.
The 500-base background arm gives 0.02511 and 4.66%, respectively.

| Task / conditional-parametric correction | Prior normalized selected AUROC | Unit-floor selected AUROC | Unit-floor gain over conventional DWM |
| --- | ---: | ---: | ---: |
| HepG2 CCCTC-binding factor (CTCF) | 0.65208 | 0.64902 | 0.06970 |
| K562 CTCF | 0.70119 | 0.69044 | 0.09976 |
| K562 MEF2D | 0.55902 | 0.57052 | 0.08406 |
| HepG2 ZNF362 | 0.55547 | 0.54461 | 0.01700 |
| HepG2 ZNF384 | 0.58059 | 0.58792 | 0.03954 |
| K562 ZNF384 | 0.53105 | 0.47587 | -0.01050 |

The K562 ZNF384 loss is already present with fixed geometry: its parametric
AUROC falls from 0.51482 to 0.48014 when changing normalization. Geometry
selection adds further loss relative to the prior selected detector. It is
therefore not solely a geometry-selection failure.

Unit-floor parametric fixed-center protection holds in 10/15 HepG2 CTCF folds
and 12/15 K562 CTCF folds, but only 5/15 K562 MEF2D, 5/15 HepG2 ZNF362,
6/15 HepG2 ZNF384 and 3/15 K562 ZNF384 folds. The fixed center is plus/minus
five bases, compared with shoulders at 15–40 bases, requiring a negative
central contrast below both shoulders. Selected-geometry protection is also
recorded separately. Numerical discrimination gains are not established
footprint gains in these difficult families.

Of 2,205 fits, 464 select an eligible geometry, 1,225 use the underpowered
fallback and 516 use the no-protection fallback. Full-cohort, unsmoothed
contrast plots retain all 21 tasks and show three correction arms. They have
no uncertainty bands and do not substitute for held-fold evidence.

## Validation and numerical limits

Nine focused tests passed in 2.04 seconds. Frozen-artifact verification passed.
Independent least-squares normalization reproduced all 2,205 geometry choices,
4,410 folds (4,144 evaluable), 294 summaries, 140 contrast curves and 4,410
paired effects. Fold discrimination and paired effects agree within 1e-12.

Exact rounded-score parity is **false**: one HepG2 ZNF362, 500-base background,
chromosome-7 selected-detector score changes from 0.0098870447 to 0.0098870448
under the equivalent solver. It changes neither selection nor fold metrics.
All other saved scores agree exactly. Report this as bounded numerical
sensitivity, not exact parity or proven threshold stability.

The original exact-score verifier and its failure are preserved under
`verification_attempt1/`. The next verifier completed the model/score/shape
checks but rejected a paired-gain difference of 5.55e-17 at a zero absolute
tolerance after CSV reload; that snapshot is in `verification_attempt2/`.
The final verifier explicitly reports score sensitivity and checks metric
parity at 1e-12. Neither scientific outputs nor normalization rules changed.

## Artifacts and decision

Results: `benchmarks/results/unit_floor_geometry_20260922_v1/`.

| Artifact | SHA-256 |
| --- | --- |
| Configuration | `e27164c9a408136c2dbb6253990661be9d0826f12a022b169d288837452a7b64` |
| Result manifest | `fba4eac366f1770d947247bb07d39cf14b94ee39cbc5153f96a7bb17df16a4c9` |
| Job specification | `821f26dc95763dd11069def83c9e297ea152c6e7d2e4a45ac6c5be81ddb6956e` |

Durable launch/exit records are in `unit_floor_geometry_jobs_20260922/execute_attempt1/`.
Execution exited zero in 29.456 seconds; this is not complete-workflow central
processing unit (CPU) or memory qualification. Plots are in
`unit_floor_geometry_plots_20260922_v2/`; v1's overlapping title/legend and
plotting source are preserved. Only the layout changed between plot versions.

Do not advance the unit floor as a general replacement or tune floors per
task. Next specify a count-conditioned central-versus-shoulder candidate,
using explicit strand counts and frozen sequence probabilities, before
evaluation on training data. Its statistical assumptions, low-count behavior,
accuracy and protection shape must be tested rather than inferred from a
nominal probability formula. All safety, uncertainty, replication and
complete-workflow resource gates remain outstanding. Main and untouched
promotion labels/final controls remain unchanged.
