# Covariate-matched ChromBPNet comparison sensitivity

This is post-hoc development evidence on previously examined K562/HepG2
chromosomes. It keeps the direct parametric model, all six comparison methods,
the original scoring geometry, and each cell's recorded dispersion unchanged.
It does not update any candidate, original comparison, safety result or
promotion decision. The broader research goal remains incomplete.

## Frozen diagnostic and provenance

`manifests/chrombpnet_covariate_diagnostic_v1.json` fixes three policies:
motif-score matching; motif plus local accessibility (primary); and those two
plus observed cut total. Matching uses 1:1 positive/negative pairs within a
chromosome, requiring every feature difference to be at most 0.2 training
standard deviations. Assignment maximizes valid pair count before minimizing
squared distance. Coordinate ordering makes ties deterministic.

Means and standard deviations come only from chr1–15. The accessibility and
observed-total variables use `log1p`. Accessibility is the 201-position
normalized raw-track sum used by the archived matching tables. The direct
observed total comes from the BAM-derived profile and is a separate covariate.
The archived training tables lacked accessibility, so the same recorded 10M
uncorrected tracks were used to supply it. Their sums reproduce all 9,842
HepG2 and 6,198 K562 archived test accessibility values within floating-point
tolerance. Track SHA-256 values are locked in the plan.

The first attempt stopped on that missing training column, before matching or
adjusted evaluation. The completed second attempt freezes the configuration,
training scales and every matched cohort before reading model score columns.
Selection uses the known ChIP classes and covariates; it is not label-free or
new holdout evidence. Every declared policy is reported, including failures.
No policy or scoring parameter was selected after seeing adjusted performance.

Original scores are reused. Profiles for the new subsets are reconstructed
from the recorded artifacts, with all reconstructed scores checked against
the original six-method score table. Metric intervals use 1,000 chromosome
block draws; aggregate bands use 500 draws. The five test chromosomes provide
only five independent blocks. Bands are pointwise diagnostic intervals,
not calibrated simultaneous functional uncertainty.

## CTCF after primary matching

| Cell | Pairs | Motif-score SMD | Log-accessibility SMD | DWM AUROC / AUPRC | Parametric AUROC / AUPRC | ChromBPNet bias AUROC / AUPRC |
| --- | ---: | ---: | ---: | --- | --- | --- |
| HepG2 | 634 | −0.0005 | 0.0092 | 0.6109 / 0.5958 | 0.6586 / 0.6715 | 0.6697 / 0.6969 |
| K562 | 609 | 0.0049 | 0.0039 | 0.6500 / 0.6203 | 0.7013 / 0.6869 | 0.7122 / 0.6954 |

HepG2 retains an AUROC gain of 0.0477 over DWM (95% interval 0.0146–0.0765)
and a relative AUPRC gain of 12.70% (7.80–16.89%). Its gains over raw also
have positive intervals: AUROC 0.0154–0.0494; relative AUPRC 3.73–10.44%.
The parametric model remains below ChromBPNet bias in both metrics.

K562 retains an AUROC gain of 0.0513 over DWM (0.0290–0.0665) and a relative
AUPRC gain of 10.73% (4.51–14.97%). Its intervals against raw include zero:
AUROC −0.0064 to 0.0587; relative AUPRC −0.19% to 10.11%. Thus the matched
analysis does not establish a gain over raw in both cells. The K562 interval
against ChromBPNet bias also includes zero; that does not establish equivalence.

Both CTCF tasks retain canonical central protection. HepG2's change is chiefly
stronger shoulders: its center contrast becomes less negative (−0.066 versus
−0.105 for DWM), while center-to-shoulder depletion increases from 0.234 to
0.266. K562 depletion increases from 0.293 to 0.339. The standardized aggregate
separation gains are 16.42% and 13.42%, below 25%; CTCF is a positive control,
not a second difficult TF family.

HepG2 keeps dispersion 0; K562 keeps 1.6144318781609228. They remain separate
recorded experiments, not a common-operator cross-cell validation.

## Prespecified sensitivity policies

| Policy | HepG2 CTCF pairs | HepG2 AUROC gain / relative AUPRC gain over DWM | K562 CTCF pairs | K562 AUROC gain / relative AUPRC gain over DWM |
| --- | ---: | --- | ---: | --- |
| Motif only | 756 | 0.0504 / 13.37% | 750 | 0.0483 / 7.44% |
| Motif + accessibility (primary) | 634 | 0.0477 / 12.70% | 609 | 0.0513 / 10.73% |
| Motif + accessibility + observed total | 579 | 0.0547 / 14.30% | 542 | 0.0665 / 14.58% |

The primary match leaves log-observed-total SMDs of 0.246 and 0.237. The third
policy reduces all three matched-covariate SMDs to at most 0.012 in both cells.
CTCF gains over DWM remain supported there. HepG2 still has positive intervals
against raw, while K562 still has intervals crossing zero. Motif-only matching
also leaves K562's AUROC interval against raw crossing zero. No sensitivity
policy is promoted to replace the declared primary policy.

## Other TFs and general gates

Seven of 21 primary-policy tasks retain at least 200 pairs and maximum absolute
matched-covariate SMD at most 0.1. Fourteen are underpowered. HepG2 HMGA1 has
only one pair, so it remains in the matching/power table but has no aggregate
or method-metric estimate. The remaining 20 tasks yield 360 method/task/policy
rows across all three policies and all six methods.

| Eligible primary task | Pairs | Parametric AUROC | AUROC gain over DWM | Relative AUPRC gain | Canonical protection |
| --- | ---: | ---: | ---: | ---: | --- |
| HepG2 CTCF | 634 | 0.6586 | 0.0477 | 12.70% | Yes |
| HepG2 FOXA1 | 491 | 0.5121 | −0.0062 | 0.37% | No |
| HepG2 ZNF362 | 689 | 0.5663 | 0.0267 | −0.68% | No |
| HepG2 ZNF384 | 658 | 0.5529 | −0.0148 | −5.28% | No |
| K562 CTCF | 609 | 0.7013 | 0.0513 | 10.73% | Yes |
| K562 MYC | 253 | 0.4977 | 0.0001 | 0.08% | No |
| K562 ZNF384 | 690 | 0.5705 | 0.0240 | 1.62% | No |

The within-cell eligible-task mean gains are only 0.0133 AUROC / 1.78% relative
AUPRC in HepG2 and 0.0251 / 4.14% in K562. Neither meets the general mean-gain
gates. Larger unsigned separation in ZNF384 does not recover central protection
and must not be called a footprint improvement. The matched cohorts do not
reveal a broad rescue of the difficult TFs.

Matching restricts the population and can remove real biological effects of
occupancy on accessibility. These results are conditional sensitivity estimates;
they do not identify a causal correction effect or eliminate unmeasured
confounding. The two cell contexts still do not provide new promotion evidence.
Exact-method depth/seed/replicate and resource gates, stronger occupancy-leakage
diagnostics and calibrated functional uncertainty remain incomplete. Final
naked-DNA replicate 3 and SK-N-SH/GM23338 labels remain unopened.

## Artifacts and verification

Output directory:
`benchmarks/results/chrombpnet_covariate_diagnostic_20260908_v2/`.
It includes frozen scales/cohorts, all balance and power statuses, per-cell
profile caches, six completed cell/policy evaluations, metric/bootstrap/profile/
shape tables, and primary-policy aggregate PDFs (11 HepG2 pages, nine K562
pages). All 20 pages were visually checked, including full-page CTCF renders.

Configuration ID:
`9376267c3a2a0b76bb8221aa5ab6079f5bdf32f2fb5cd25fe9889627b68afafa`.
Selection ID:
`0b1c9656e99e2c836ec9afa1c2ab9f5270312bbc876b713b4119b8d985e2a807`.
Diagnostic ID:
`96b92f8d72807b05a1e8e9d5d0d5bdd8b68acd5be1d5e9e3134f0a4e3edb07fd`.
Calculation script SHA-256:
`cf13af449d861c5c52d4148103de8b1b6ba09d99ef14472de9a1979ba21b2f1e`.
The convenience `parametric_sensitivity_summary.tsv` derives gains directly
from the saved per-method metrics and interval/shape flags; its SHA-256 is
`69309f696753a64c0bc8e8cf8631364e34b12a6f32ae173ee219ef08e93e88fc`.

The focused suite passes 54 tests. Checks cover training-only scales,
maximum-cardinality caliper matching, score/row-order independence, chromosome
boundaries, degenerate covariates, raw-track agreement, empty results and
checksum-protected resumption. Independent verification checked the nested
selection outputs, both profile caches, all six evaluations and calculation
code hashes. Reinvoking the completed diagnostic verified and skipped it.

Calculation code is preserved in commit `c39dcd9`. A subsequent verification
fix checks cohort files nested beneath the parent manifest and adds
`--verify-only --code-ref c39dcd9` for historical artifacts after code updates.
The recorded calculation-code SHA-256 values are checked against Git blobs;
the verifier does not rerun matching, scoring, bootstrapping or rendering.
This verification passed on the completed run, and the expanded suite passes
56 tests, including nested-file mutation and wrong/mutable code-reference
rejection. The scientific outputs and their identifiers are unchanged.
