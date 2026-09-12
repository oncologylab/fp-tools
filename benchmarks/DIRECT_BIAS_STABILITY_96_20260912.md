# Complete 96-run direct-bias development stability

All 60 missing conventional DWM corrections and all 96 stability evaluations
are complete and independently verified. CTCF retains useful development
gains in both cells, but the candidate does not pass general promotion gates.
Production DWM and the final manuscript are unchanged.

## Scope and provenance

This is a post-hoc development experiment on already examined chr19–22/X,
not untouched validation. HepG2 and K562 each contribute three biological
replicates, five seeds (2026–2030) at 10M, 25M and 50M, and three deterministic
full-depth endpoints: 48 runs per cell. Seeds are not independent biological
replicates. Candidate coefficients, dispersions, geometry and labels were
not retuned. Expectation support and complete matched pairs were frozen
across all runs before new extraction or scoring.

The four identical-site arms are raw, fixed original DWM, sample/depth-specific
DWM, and frozen parametric bias, each with conventional geometry. This isolates
expectation/correction effects under a fixed detector; it does not establish
an improved detector or causally isolate sequence bias from occupancy.
The completed ChromBPNet reference comparison is recorded separately in
[the reference audit](CHROMBPNET_AUDIT_20260908.md).

Preparation finished at 2026-09-10 20:52:33 UTC; the dependent evaluation
finished at 21:26:54 UTC, both exit 0. Independent read-only checks completed
on September 12. No completed scientific stage was rerun. All 36 archived
extractions were reused; no archived evaluations were reused because common
support changed. The resulting matched CTCF cohorts contain 573 HepG2 and
529 K562 pairs, compared with 577 and 534 in the archived 36-run analysis.

Ignored outputs: `results/direct_bias_stability_96_20260908_v1/`.
The directory contains 7,872 metric rows, 7,872 shape rows, 2,016 coverage
rows, 288 bootstrap rows, per-run profiles and a 21-page aggregate PDF.

| Record | SHA-256 |
| --- | --- |
| Evaluation plan | `b646b58e22ddce161e58894ac3beb06fcb37d1e64739ed432309a0322e0ce7e1` |
| Configuration freeze | `9e870f9a8b66bf3e37f7c330e097ccbecef8c83eb36f4b0985b97e9cb7f1d599` |
| Cohort freeze | `3f9ce350545bde4931a923f74ba92039b80b50c81df7f618d033ec3ab3216ea8` |
| Result manifest | `52fc25ca38d8a83b301ace2980a4947583a5605743cf729a66c6d9c8051da798` |

Result ID: `66909de3c482dac223864055c9684103d22a0774d9f64010c770de0a29496034`.
Preparation configuration and sampling checks are recorded in
[the input report](DIRECT_BIAS_REPLICATE_DEPTHS_20260908.md).

## CTCF results

Full-depth matched results for the frozen parametric candidate follow.
Delta and relative gain use the sample-specific DWM comparator.

| Cell / replicate | AUROC | AUPRC | AUROC gain | Relative AUPRC gain |
| --- | ---: | ---: | ---: | ---: |
| HepG2 1 | 0.735585 | 0.765447 | 0.113362 | 28.55% |
| HepG2 2 | 0.705938 | 0.753633 | 0.100567 | 28.10% |
| HepG2 3 | 0.735805 | 0.764750 | 0.122341 | 26.41% |
| K562 1 | 0.768572 | 0.754478 | 0.150907 | 26.57% |
| K562 2 | 0.761615 | 0.752067 | 0.110867 | 20.20% |
| K562 3 | 0.776184 | 0.762464 | 0.127608 | 23.29% |

All six full-depth endpoints have positive chromosome-bootstrap lower bounds
for AUROC and relative AUPRC gains against raw and both DWM comparators.
Across all 96 matched CTCF runs, both point-estimate metrics improve over raw
and sample DWM. AUROC gains over sample DWM range from 0.021539–0.122341
in HepG2 and 0.022088–0.150907 in K562. Relative AUPRC gains range from
8.10–28.55% and 5.31–26.57%, respectively. Thus improvement direction does
not mean every run exceeds the prespecified magnitude thresholds.

The principal low-depth checks are the six 10M seed-2026 runs. K562 replicate
3's AUROC and relative AUPRC intervals against sample DWM cross zero (lower
bounds -0.014288 and -0.011589). Against raw, both intervals cross zero in
K562 replicates 1 and 3. All six low-depth comparisons against fixed DWM
have positive lower bounds. These intervals use five chromosome blocks;
1,000 bootstrap draws do not create 1,000 independent biological observations.

Canonical protection passes in 46/48 matched HepG2 runs and 48/48 K562 runs.
The exceptions are HepG2 replicate 2 at 10M, seeds 2026 and 2028. Each has
only one significant central point; seed 2028 also has a positive central
contrast (0.003386). All original-support CTCF runs pass the shape rule.
The expanded result therefore qualifies the earlier 36-run stability finding.

## Per-TF and aggregate audit

Each nonempty matched task has 48 runs. Pair counts and candidate AUROC
ranges are shown below; canonical counts alone do not establish power.

| Cell | TF | Pairs | Powered | AUROC range | Canonical runs |
| --- | --- | ---: | --- | --- | ---: |
| HepG2 | ARID3A | 118 | no | 0.451–0.566 | 0 |
| HepG2 | CTCF | 573 | yes | 0.618–0.736 | 46 |
| HepG2 | FOXA1 | 447 | yes | 0.402–0.513 | 0 |
| HepG2 | MEF2A | 54 | no | 0.445–0.616 | 17 |
| HepG2 | MEF2D | 77 | no | 0.484–0.627 | 11 |
| HepG2 | MYC | 73 | no | 0.436–0.561 | 0 |
| HepG2 | TBP | 38 | no | 0.447–0.636 | 11 |
| HepG2 | TCF7L2 | 108 | no | 0.395–0.538 | 0 |
| HepG2 | ZNF362 | 641 | yes | 0.529–0.627 | 0 |
| HepG2 | ZNF384 | 639 | yes | 0.519–0.591 | 0 |
| HepG2 | ZNF558 | 21 | no | 0.401–0.676 | 14 |
| K562 | ARID3A | 51 | no | 0.449–0.626 | 9 |
| K562 | CTCF | 529 | yes | 0.674–0.776 | 48 |
| K562 | FOXA1 | 8 | no | 0.297–0.750 | 23 |
| K562 | MEF2A | 22 | no | 0.519–0.777 | 23 |
| K562 | MEF2D | 65 | no | 0.430–0.620 | 0 |
| K562 | MYC | 225 | yes | 0.478–0.543 | 0 |
| K562 | TBP | 44 | no | 0.381–0.569 | 0 |
| K562 | TCF7L2 | 10 | no | 0.310–0.710 | 15 |
| K562 | ZNF384 | 657 | yes | 0.545–0.597 | 0 |

HepG2 HMGA1 has no retained matched cohort; its original-support task is
underpowered. It remains in the task denominator. No other powered matched
TF has canonical protection. FOXA1 and MYC show central enrichment rather
than protection. ZNF362/ZNF384 show broad shoulders without a protected
center; larger amplitude or occupancy discrimination is not footprint rescue.
Underpowered MEF2, TBP and other tasks cannot certify family-level gains.

All 21 aggregate PDF pages were visually inspected for labels, layout, center,
shoulders, asymmetry and replicate/depth variation. No obvious clipping was
found. Thin curves show runs and thick curves their mean; the mean can conceal
the two failed CTCF runs. CTCF's shoulders strengthen, especially at full depth,
without a uniformly deeper center. Candidate matched shoulders stay positive:
HepG2 left 0.161–0.443, right 0.150–0.410; K562 left 0.158–0.378, right
0.143–0.340. Shoulder asymmetry ranges from -0.111–0.094 and -0.166–0.053.
Functional separation increases over sample DWM by 6.01–72.07% and
6.82–73.53%, respectively; this does not satisfy a two-difficult-family gate.

Coverage uses positive expectation in every run and retains zero-cut sites
and complete matched pairs. CTCF motif-score absolute SMD is below 0.012
and reference accessibility SMD below 0.008. Current observed-count SMD
drifts up to 0.146 in HepG2 and 0.154 in K562. Across all matched tasks,
observed-count SMD ranges from -0.144 to 0.453. Matching does not eliminate
all current-depth accessibility confounding. The separate
[36-run balance sensitivity analysis](DIRECT_BIAS_DEPTH_BALANCE_20260908.md)
is post-hoc diagnostic evidence and must not be substituted for this design.

## Safety and promotion disposition

Existing independent naked-DNA replicate-2 safety is narrow CTCF evidence:
paired candidate calls are 0/121 in HepG2 (Wilson upper bound 3.0771%) and
1/114 in K562 (4.8020%; increase 0.8772 percentage points). Other tasks do
not have sufficient matched DWM coverage to certify safety. These archived
checks were not extended by the ATAC stability analysis. They do not prove
that all residual sequence structure or occupancy leakage is absent.

Final naked-DNA replicate 3 and SK-N-SH/GM23338 labels remain unopened.
The original freeze SHA-256 remains
`7d29af02e99507864e6cc194d4ec8e30df79f95bf4ab55d468c6c30eea70a354`,
with all 13 `labels_read` flags false. The result manifest records no candidate
retuning, no new holdout/control reads and `promotion_eligible: false`.

Two difficult families have not improved; low-depth shape and interval
exceptions remain. Functional bands are pointwise diagnostics, not calibrated
functional uncertainty. The archived small-model/component resource result
does not certify end-to-end CPU or memory gates. Untouched safety, external
generalization and absence of occupancy leakage are not established. No
method is promoted, no holdouts are consumed, and no new Box report is uploaded.

## Independent verification

All 60 sampling targets and 180 BAM/index/metadata files passed the earlier
independent matrix, checksum and alignment checks. All 60 correction stages
subsequently passed frozen-command, input BAM/BAI, seven-output checksum and
four-BigWig readability checks. Historical preparation verification used
code reference `604a4f2` and preserved the frozen configuration.

The live-runtime check initially rejected a changed Python build timestamp:
both runtimes report Python 3.12.3/GCC 13.3.0, but the frozen build date is
June 19 and the current build date August 31, 2026. Declared library versions
and samtools 1.21 match. Supported historical verification passed; it verifies
the recorded computation and does not certify identical current interpreter
bytes. No scientific computation was rerun in the changed runtime.

From the research root, the following read-only verification exited 0 and
reported all 96 runs, frozen inputs, archived reuse and outputs verified:

```bash
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  "$FP_PYTHON" benchmarks/scripts/extend_direct_bias_stability.py \
  --plan benchmarks/manifests/direct_bias_stability_96_v1.json \
  --data-root "$FP_DATA_ROOT" \
  --outdir benchmarks/results/direct_bias_stability_96_20260908_v1 \
  --verify-only --code-ref c9ecd82
```

`FP_PYTHON` identifies the existing project virtualenv interpreter and
`FP_DATA_ROOT` the original data repository. Verification includes the exact
96-run matrix, 192 nested extraction/evaluation stages, source/code identities,
cohort freezes, 36 extraction reuses, zero metric reuses and output hashes.

Documentation-contract tests passed (20 tests, 0.46 seconds); the strict
MkDocs build passed (3.47 seconds). The browser audit passed all 32 pages at
three viewport sizes while emulating a dark system preference.
No package source or interface changed.
No preparation, evaluation or queue process was active at the final process
check. Main remained clean at `214dc1a`; the frozen manuscript PDF, source ZIP
and response DOCX each retained its recorded SHA-256.
