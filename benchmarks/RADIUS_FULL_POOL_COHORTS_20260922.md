# Direct three-covariate training cohorts

Direct matching from the full existing training pool restores power lost by
trimming the older pairs. These cohorts balance distant accessibility as well
as motif score and central counts. No footprint performance was evaluated
during cohort construction, and earlier primary and secondary results remain
unchanged.

## Frozen construction

The protocol and execution inputs/code were frozen before matching. Only
chromosomes 1–15 were selected. The shared source tables include previously
examined tuning rows; those were not selected or evaluated. Supporting sequence,
alignment counts, dinucleotide weight matrix (DWM) expected profiles and the
existing parametric expectations must all be valid on common support.

Within each task, matching uses motif score, log1p central 201-position counts,
and log1p counts at absolute offsets 101–500. Each covariate is divided by its
standard deviation across all supported sites; constant scales use 1.
Assignment maximizes the number of pairs within Euclidean distance 0.5 and
then minimizes total squared distance. The larger class is not subsampled.
The largest distance-ordered prefix meeting absolute standardized mean
difference (SMD) at most 0.1 for every covariate is retained. This final prefix
is not claimed to be the largest balanced subset among all possible subsets.
No scores or aggregate-shape outcomes enter either step.

## Support and power

All 21 tasks remain explicit. The resulting cohorts contain 4,053 pairs:
ten tasks have at least 200 pairs, ten are underpowered, and HepG2 HMGA1
has no balanced pairs. HMGA1 has only three supported positive sites, so this
is a label-power limitation rather than proof of absent footprint information.

| Task | Trimmed previous pairs | Direct three-covariate pairs |
| --- | ---: | ---: |
| HepG2 CTCF | 94 | 300 |
| K562 CTCF | 122 | 245 |
| HepG2 ZNF362 | 256 | 291 |
| HepG2 ZNF384 | 89 | 303 |
| K562 ZNF384 | 92 | 238 |
| HepG2 MEF2A | 168 | 190 |
| K562 MEF2A | 98 | 169 |

CTCF means CCCTC-binding factor. Both MEF2A tasks remain below the fixed
200-pair threshold. Counts from different cohort definitions must not be
interpreted as biological replication or pooled to manufacture power.
Every retained task meets all three mean-balance thresholds. Distributional
balance, statistical uncertainty and model performance remain separate checks.

## Verification

Fourteen relevant tests passed in 1.43 seconds. The supervised matching run
exited zero in 15.31 seconds. Separate frozen-artifact verification passed.
An independent sparse-graph maximum-bipartite-matching calculation confirmed
the pre-balance maximum-cardinality count for every task. Independent checks
also confirmed unique sites, complete positive/negative pairs, the fixed
distance caliper, membership in supported training sites, and all three
reported SMD values. No completed extraction or earlier study was rerun.

Artifacts are under `results/radius_full_pool_20260922_v1/`: all training sites
with support indicators, selected pair memberships with original and wider
profile indices, and all task diagnostics. Durable execution evidence is under
`results/radius_full_pool_jobs_20260922/execute_attempt1/`.

SHA-256:

- Configuration: `8034ffae8515de0f93af4be697d422dafbc149b8146153c9351d56ce64fffa37`
- Cohort manifest: `680a16cc6d95dfd6ee37155bb9842ede3213695d9e016ab783ef219ba66f755d`

The next comparison requires a separate execution freeze. It must recompute
all seven fixed scoring arms and aggregate curves on these exact sites;
metrics or shape conclusions from previous cohorts cannot be reused as if
membership were unchanged. Preserve the complete denominator and underpowered
tasks. This study does not justify opening untouched transfer labels or final
negative controls, selecting a model for production, or claiming any promotion
or full-workflow computational gate. Main remains unchanged.
