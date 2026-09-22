# Gamma background: identical-site training comparison

The new low-count background produces complete profiles on the existing
matched training cohorts, but does not establish improved footprint detection.
Do not advance this unchanged policy to untouched holdouts. No production
code, default, manuscript, or main-branch file changed.

## Design and scope

Five arms used identical sites and conventional geometry: center within
7 base pairs, shoulders 12–40 base pairs. They were raw counts, dinucleotide
weight matrix (DWM) expected-profile residuals, existing conditional
parametric residuals, global Gamma posterior means with fixed-dispersion
deviance residuals, and global Gamma predictive moment-standardized residuals.
The fixed dispersion remained 0 for HepG2 and 1.6144318781609228 for K562.
The previous two correction arms retain their frozen per-window expectation
rescaling; global expectations do not use that rescaling. Thus comparing
existing and global parametric arms tests expectation construction, while
comparing the two global arms tests residual formulation. Detector geometry
does not change. DWM here is the matched expected-profile residual comparator,
not a newly executed complete production workflow.

Only chromosome 1–15 rows were selected and scored. The existing source tables
contain already examined chromosome 16–18 rows, which are loaded by the shared
reader but excluded from this evaluation. No test, new transfer holdout,
or naked-DNA control labels were evaluated. No model was fitted or selected.
This is training evidence, not independent validation.

All 3,949 original balanced pairs were retained in every arm: 7,898 sites,
21 tasks, and 105 metric rows. Ten tasks had at least 200 pairs; eleven were
underpowered and remain explicit, including HepG2 HMGA1 with only two pairs.
No coverage exclusions changed the previously frozen motif-score and count
balance. Pairing was not repeated after observing outcomes. Scores from
zero-read neighborhoods remain ties. Aggregate plots use a fixed fallback
scale of 1, outer-flank linear detrending/root-mean-square normalization, and
clipping at 12; plots are descriptive, without calibrated confidence bands.

## Findings

For CCCTC-binding factor (CTCF), area under the receiver operating
characteristic curve (AUROC) was:

| Cell | DWM residual | Existing parametric | Gamma fixed dispersion | Gamma predictive |
| --- | ---: | ---: | ---: | ---: |
| HepG2 | 0.594552 | 0.606990 | 0.588511 | 0.580379 |
| K562 | 0.618622 | 0.669707 | 0.632895 | 0.626161 |

CTCF central protection remains visible, but the global model loses
discrimination relative to the existing conditional parametric correction.
The full 21-task descriptive mean AUROC differences from DWM are +0.003937
for Gamma fixed dispersion and +0.003172 for Gamma predictive. These averages
include underpowered tasks and are not promotion-gate estimates.

HepG2 ZNF362 has a small canonical protection contrast and Gamma predictive
AUROC 0.560791 versus DWM 0.511536. However, HepG2 ZNF384 has numerical
discrimination gain without canonical central protection, while K562 ZNF384
has little discrimination gain. FOXA1 and MYC contrasts remain inconsistent
with canonical protection. Underpowered K562 MEF2A is exploratory only.
These results do not demonstrate improvements in two difficult families and
two cell contexts. All 21 contrast panels were visually reviewed; no plots
were dropped because they were noisy or unfavorable.

The separate synthetic stress study already shows inadequate uncertainty
under overdispersion. Neither that issue nor exact-operator biological safety,
replicate/depth stability, or full-workflow resource gates is resolved here.
Finite low-count estimates solve a coverage problem, not the scientific goal.

## Execution, verification and files

The first attempt was deliberately terminated after 84.12 seconds when archive
access was found to decompress full arrays for each site. It produced no
metrics. Its original code snapshot, configuration and termination receipt
are preserved under `results/gamma_training_jobs_20260922/`. Loading each
archive member once fixed the implementation cost without changing scientific
settings. The second attempt used a fresh freeze and finished in 13.79 seconds.
This cached-profile execution time is not an end-to-end resource benchmark.

Thirteen focused tests passed in 2.51 seconds. Checks cover coordinate/strand
identity, paired exclusions, single-decompression access, exact score
serialization, and the Gamma-model and simulation regressions. The separate
comparison `--verify-only` passed. The independent review recalculated all
105 AUROC and area under the precision-recall curve (AUPRC) rows, verified
identical coordinates/labels across arms, and checked chromosome and class
balance. Round-trip floating-point parsing is required: default text parsing
can alter neighboring score values and change rank ties.

Artifacts:

- `results/gamma_training_comparison_20260922_v2/`: frozen configuration,
  all metrics, per-site scores, mean curves and result manifest.
- `results/gamma_training_review_20260922_v2/`: three-page aggregate PDF,
  three rendered images containing all 21 tasks, and independent review
  manifest. The automatic manifest leaves visual review false; the subsequent
  human-visible image inspection is recorded here, without rewriting it.

SHA-256 identities:

- Execution configuration: `d83fd8115756a2fa080d66573b22a806a527182182c68f04c9122a99c5709187`
- Result manifest: `0301edc448a429431f0f600b558821e63c4dad150fa31880d503d61deb551628`
- Review manifest: `e5b613ea9e104e2ca2d4512d28de4fad4b0c110bd1230eb2d705c4280372bc92`
- Aggregate PDF: `aa7446832fc5d503dee8ff0114294cdafd906bc2f899135403921373ecaef272`

No untouched labels were opened, no promotion decision is supported, and
no results were uploaded. Further work must address local-rate estimation
and biological overdispersion using permitted development data before any
new safety or holdout evaluation.
