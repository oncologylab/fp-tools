# Direct-bias development stability

All 36 available runs completed on previously examined chr19–22/X. The direct
CPU bias model has a reproducible CTCF-specific result across depths, seeds
and biological replicates. No other powered TF passes the canonical-shape
screen. General promotion remains unsupported; production DWM and the frozen
final manuscript are unchanged.

## Frozen design and coverage qualification

Calculation commit: `c453a37`. Plan: `manifests/direct_bias_stability_v1.json`.
The model, geometry and recorded dispersions were frozen before new scoring:
HepG2 0; K562 1.6144318781609228. These remain separate recorded experiments.
Scores subtract the mean center (±7 bp) from mean shoulders (12–40 bp) on
total-matched signed-deviance residuals. Display profiles separately remove
the outer-flank trend and normalize by its RMS; display normalization does
not replace the scoring operator. Canonical-shape auditing uses the fixed
±5 bp center and 15–40 bp shoulders from the earlier diagnostic.
The candidate is the mitochondrial LOG21 +4/−5 L2 0.0001 five-seed coefficient
ensemble, with NPZ SHA-256
`620f232f634423b8753e4d77f8cf72df87482c09193cb5e7e6b038917c47e89a`.

The available matrix has five seeds (2026–2030) at 10M, 25M and 50M for
replicate 1, plus all three full-depth biological replicates per cell.
Full-depth data have one deterministic endpoint per replicate. Replicates 2
and 3 still lack the depth/seed factorial panel; randomizations of one library
are not independent biological replication.

The existing partial-bias depth matrices used tuning chromosomes. Their
checksum-locked BAM input registries were reused, but their scores and tuning
profiles were not substituted for test evidence. The two existing 10M seed-2026
test count profiles were reused; only 34 missing test count profiles were
extracted. Cached strand-specific sequence propensities reproduce the original
expected cuts and comparison scores. No bias fitting, DWM correction,
ChromBPNet training or completed prediction/evaluation stage was rerun.

Preflight found that original finite support included zero-filled DWM
expectations: 65/1,961 HepG2 and 66/1,975 K562 CTCF sites. A new per-cell stage
froze positive expectation support across both the original full-depth DWM
track and every sample/depth-matched DWM track, before scoring new runs.
It retained 9,361/9,613 HepG2 and 5,972/6,120 K562 original common sites.
CTCF retained 1,890 and 1,894 sites, respectively. The exclusions are explicit
in per-site coverage tables; zero observed-count sites remain in the analysis.

The second cohort starts with the previously frozen motif/accessibility/cut-total
match at 10M seed 2026. Coverage exclusions preserve complete pairs without
rematching: 577 HepG2 and 534 K562 CTCF pairs. The same loci, labels and methods
are used in every run. Comparisons distinguish fixed original DWM expectations
from sample/depth-matched DWM expectations. At HepG2 10M seed 2026, the latter
scores CTCF at 0.6657 AUROC on qualified original support versus 0.6463 for
fixed DWM, showing that comparator estimation also contributes to the result.

## CTCF results

All 36 CTCF point estimates improve over raw, fixed DWM and sample-matched DWM
on both cohorts. All retain canonical protection. On fixed matched sites:

| Cell | Depth / replicate scope | AUROC gain over raw, range |
| --- | --- | ---: |
| HepG2 | 10M, five seeds | 0.0226–0.0307 |
| HepG2 | 25M, five seeds | 0.0357–0.0382 |
| HepG2 | 50M, five seeds | 0.0442–0.0488 |
| HepG2 | Full, three biological replicates | 0.0592–0.0627 |
| K562 | 10M, five seeds | 0.0243–0.0398 |
| K562 | 25M, five seeds | 0.0386–0.0492 |
| K562 | 50M, five seeds | 0.0356–0.0518 |
| K562 | Full, three biological replicates | 0.0575–0.0666 |

The full-depth matched comparison is:

| Replicate | Sample DWM AUROC / AUPRC | Parametric AUROC / AUPRC | AUROC gain over raw, 95% interval | Relative AUPRC gain over raw |
| --- | --- | --- | --- | ---: |
| HepG2 1 | 0.6237 / 0.5963 | 0.7362 / 0.7657 | 0.0592 [0.0396, 0.0689] | 11.62% |
| HepG2 2 | 0.6068 / 0.5892 | 0.7074 / 0.7548 | 0.0617 [0.0374, 0.0794] | 13.48% |
| HepG2 3 | 0.6146 / 0.6055 | 0.7367 / 0.7654 | 0.0627 [0.0464, 0.0707] | 10.80% |
| K562 1 | 0.6197 / 0.5971 | 0.7697 / 0.7553 | 0.0575 [0.0106, 0.0853] | 9.52% |
| K562 2 | 0.6511 / 0.6262 | 0.7610 / 0.7519 | 0.0579 [0.0185, 0.0812] | 9.88% |
| K562 3 | 0.6474 / 0.6176 | 0.7761 / 0.7628 | 0.0666 [0.0137, 0.1025] | 11.57% |

All six full-replicate AUROC and relative AUPRC intervals exclude zero against
raw and both DWM comparators. The principal reference and all full replicates
use 1,000 chromosome-block draws; other depth/seed rows report point estimates.
K562's matched 10M seed-2026 interval against raw still includes zero
(AUROC −0.0084 to 0.0601; relative AUPRC −0.73% to 10.46%). There are only five
independent chromosome blocks. These intervals are diagnostic, unadjusted for
multiple comparisons, and do not establish simultaneous functional coverage.

Full-depth matched standardized aggregate separation increases by 60.6–71.8%
in HepG2 and 56.2–72.4% in K562 versus sample-matched DWM. CTCF is a positive
control, not two difficult families. Both shoulders strengthen and their
asymmetry decreases. HepG2 replicates 1 and 2 have nearly unchanged center
contrast; replicate 3 and all K562 replicates have more negative centers.
The result must not be described as uniformly deeper central depletion.

## Balance and difficult-TF limitations

Motif and reference-accessibility SMDs on the retained CTCF pairs stay below
0.012. Current observed-cut-total SMD exceeds 0.1 in 11/36 matched runs:
HepG2 full replicate 3 and ten K562 runs, including all three full replicates.
The maximum is 0.1344. Matching the 10M reference does not guarantee balance
at other depths or in other libraries. This is conditional sensitivity,
not causal elimination of accessibility confounding.

Only CTCF passes canonical protection among powered tasks in any run. On the
full-depth matched sites, HepG2 FOXA1 has AUROC 0.402–0.428, ZNF362 0.597–0.623,
and ZNF384 0.570–0.590; K562 MYC has 0.497–0.500 and ZNF384 0.581–0.595.
All these tasks remain below raw in AUROC in every full replicate. Their
central enrichment or weak/noncanonical shape is not a footprint improvement.
Twelve original-support tasks and fourteen matched tasks remain underpowered.

Across eligible full-depth task/replicate rows, mean relative AUPRC gains over
sample-matched DWM are 4.90%/3.18% on qualified original support and
6.76%/7.59% on matched support for HepG2/K562. All remain below 10%. Mean
AUROC gains on matched support are 0.0327/0.0439, but passing one numerical
component does not satisfy the general or difficult-family gates.

## Safety and computational scope

The model, direct expectation operator and per-cell dispersions match the
existing independent-library candidate screen exactly. On informative sites
with comparable DWM expectations, HepG2 has 0/121 calls (Wilson upper 3.08%);
K562 has 1/114 (4.80%, increase over DWM 0.88 percentage points). Each excludes
one informative site without DWM expectation; the absolute-candidate screen
also passes with those included: 0/122 and 1/115. These are existing naked-DNA
replicate-2 results, not new null experiments or depth-stratified validation.
Sparse counts and absence of a large motif-residual flag do not prove absence
of occupancy leakage. Final naked-DNA replicate 3 remains untouched.

A separate guarded-load audit verifies that the coefficient model opens only
its NPZ and JSON sidecar: 10,007 bytes total, with 6,816 bytes of coefficient
arrays in memory. It does not open upstream ensemble-member models. This is a
coefficient-artifact measurement; it does not measure pipeline peak memory or
CPU runtime and excludes input data, analysis caches and the Python runtime.
Those computational gates and calibrated uncertainty remain incomplete.
Direct geometry scores are not calibrated binding probabilities; inherited
metric-table calibration fields must not be used as evidence that they are.

The locked holdout manifest still has all 13 `labels_read` flags false and
`holdout_labels_scored: false`; its SHA-256 remains
`7d29af02e99507864e6cc194d4ec8e30df79f95bf4ab55d468c6c30eea70a354`.
No SK-N-SH/GM23338 ChIP label contents were read. No promotion or integration
is authorized by this CTCF-only development result.

## Artifacts and verification

Ignored output directory: `results/direct_bias_stability_20260908_v1/`.
It contains 2,952 metric rows, 192 paired-bootstrap rows, 2,952 signed-shape
rows, 756 coverage rows and a 21-page aggregate PDF. All pages were visually
reviewed. Every frozen input, calculation-code blob and nested output passed
verification with `--verify-only --code-ref c453a37`; no completed stage reran.

- Configuration ID: `4a0836bdf6370dd444bce20501d209236f924dc1939956df877ca06078ce9db7`.
- Result ID: `a4024f20a273e601156961d8a5bb517d4d049a69b164dcfe7ddc162a62ff1e28`.
- Configuration JSON SHA-256: `7998850dd9eece77b668e0df04c8dbfa02a4bbca30adbdcee83170cb97add433`.
- Result JSON SHA-256: `94b4b887bef9431240967c7392970427a40d59470f8b1bc5ad069fa61bc0cfe0`.
- Model-dependency audit SHA-256: `c7e37324a82ec3f4fece336b61d7db64b72724a19fc0fa5054d61548cc0b7abd`.

The convenience `stability_parametric_summary.tsv` joins the candidate metrics,
signed shapes and coverage by cell/sample/depth/seed/cohort/TF, then joins the
sample-DWM functional separation and computes their relative difference.
It has 738 rows and SHA-256
`0e7b439387888e890e58a3f79005fd6f0746a456dd4fdbd4029901e75f56b7a8`.
It is a derivative summary, outside the immutable result manifest.

The final one-page report is in
`results/reports/CTCF_direct_bias_stability_20260908_v2/`. Its renderer requires
all 36 CTCF directions, canonical protection, positive full-replicate metric
intervals and exact-model/dispersion candidate safety before writing a report.
The PDF was visually checked after a layout-only correction; all six report
files passed local checksum and portability checks.

- Report ID: `69008e863b29730b82415325bce59f343a5b7f61fc4ebcad7114b57a4bbc23cb`.
- PDF SHA-256: `68baae527fd76bbb25871962c321d54219a1653a2d36871338d5f53f81ea70f3`.
- Box folder: `yilab:Yaoxiang/fp-tools/2026-09-08_CTCF_direct_bias_stability_coverage/`.
- Remote verification on 2026-09-08: `rclone check --download --one-way`
  reported six matching files and zero differences; the remote listing has
  exactly the PDF, PNG, metrics TSV, README, manifest and SHA256SUMS.

Validation: 104 focused scientific/report/docs/profile tests passed. The
strict MkDocs build passed, and the browser audit passed for 32 pages at
three viewport widths. Research commits and reports do not alter production
defaults, manuscript files or existing Box reports.
