# Secondary broad-count balance sensitivity analysis

The larger-background result remains exploratory. Balancing distant counts
weakens HepG2 ZNF384 and leaves both ZNF384 tasks underpowered. HepG2 ZNF362
retains its original balanced cohort and numerical result. This analysis
does not establish a general footprint improvement.

## Frozen design

The cohort rule was frozen before selecting pairs and did not read footprint
score values or shape outcomes. It uses three covariates: motif score,
log1p central 201-position counts, and log1p counts at absolute offsets
101–500. Existing complete pairs are ordered by Euclidean covariate distance
after within-class pooled-standard-deviation scaling, with stable pair-ID
ties. The largest prefix with absolute standardized mean difference (SMD)
at most 0.1 in all three covariates is retained. There is no new caliper,
rematching or optimization using footprint results. Already balanced cohorts
are preserved in full. This is a prefix rule, not globally optimal subset
selection, and mean balance does not establish distributional balance.

All selection and evaluation is restricted to training chromosomes 1–15.
Original primary cohorts/results remain unchanged. A separate execution
freeze preceded recomputation of seven methods' rank metrics on the selected
sites. All methods use exactly the same secondary site set for each task.

## Power and numerical results

The secondary cohorts retain 2,517 of 3,949 original pairs. Only three tasks
retain at least 200 pairs: HepG2 ARID3A, HepG2 ZNF362, and K562 MEF2D.
Seventeen tasks are underpowered; HepG2 HMGA1 has no balanced prefix and is
retained explicitly with unavailable metrics. All 21 tasks and all 147
task/method rows are reported; 140 rows have scores and seven are unavailable.

Area under the receiver operating characteristic curve (AUROC):

| Task | Original → retained pairs | DWM residual | Radius 500 | Uniform radius 500 |
| --- | ---: | ---: | ---: | ---: |
| HepG2 ZNF362 | 256 → 256 | 0.511536 | 0.608994 | 0.569305 |
| HepG2 ZNF384 | 298 → 89 | 0.575369 | 0.600240 | 0.596831 |
| K562 ZNF384 | 241 → 92 | 0.540938 | 0.609700 | 0.570475 |
| HepG2 CTCF | 299 → 94 | 0.514486 | 0.563603 | 0.529991 |
| K562 CTCF | 280 → 122 | 0.584722 | 0.616971 | 0.597689 |

DWM means dinucleotide weight matrix; the comparator is the matched
expected-profile residual, not a new complete production run. CTCF means
CCCTC-binding factor.

For HepG2 ZNF384, radius-500 area under the precision-recall curve (AUPRC)
is 0.555371 versus DWM 0.583706, and its AUROC is close to the uniform
control. The stronger primary result does not survive this check unchanged.
K562 ZNF384 retains a numerical gain, but its 92 pairs are insufficient for
the prespecified power requirement. HepG2 ZNF362's unchanged result supports
continued investigation of that task, not transfer or a package-wide claim.

The sensitivity analysis changes the sampled population as well as reducing
power; it cannot prove that the original gain was caused by confounding.
No aggregate profiles were recomputed on these subsets. Original aggregate
shape evidence must not be attributed to the smaller cohorts. These are
rank diagnostics only, without significance, uncertainty or safety claims.

## Verification

Six relevant tests passed in 2.59 seconds, covering pair preservation,
covariate-only exclusion, unavailable cohorts, all 147 synthetic evaluation
rows and the accessibility features. Both cohort and score stages passed
separate artifact verification. Independent checks confirmed intact pairs,
all three balance thresholds for every retained task, exact coordinate/label
joins, 140 recomputed metric rows and seven unavailable rows. AUROC was also
recomputed using the rank-sum identity, independently of the scoring runner.

Artifacts:

- `results/radius_balanced_pairs_20260922_v1/`: frozen rule, all original
  sites with selection indicators, retained sites and all 21 diagnostics.
- `results/radius_balanced_scores_20260922_v1/`: separately frozen secondary
  evaluation and all 147 rank-metric rows.

SHA-256:

- Cohort configuration: `2680ae79134ca3870162922b4f2a35ac1d48fa26ac904fe47eb45924c66e46a0`
- Cohort manifest: `c8aa0f656ec23332ca627aa040884aebf27a1e63f4c2c752b160ca7c556c0383`
- Evaluation configuration: `ca24669c3218bc9dff0330b864262f37b91c3389bc4f063b756fcf390552eb85`
- Evaluation manifest: `774204eabf8dc1b8e4538d9d55b1260254e35f72e67b5cfdeb01c741ebe4d388`

## Next development step

The original two-covariate pair assignments severely limit power once distant
counts are considered. Before any further model selection, freeze a direct
three-covariate matching analysis using the already available full training
site pool. Preserve both current primary and secondary evidence, report
support and lost power for all tasks, and do not use scores to construct the
cohorts. Recompute shape only on the resulting exact cohorts; do not borrow
shape conclusions from previous site sets. No new holdout or final negative
control labels should be opened for this development work.

Main, production scientific defaults, releases and the manuscript remain
unchanged. Nothing is promoted or uploaded.
