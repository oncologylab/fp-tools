# Bounded next-round development investigation

Execution began September 22, 2026; evidence packaging continued September 23.
This investigation uses an isolated research branch based on `0569e2f`.
Main, the published manuscript, released defaults and upstream projects were
not changed. Protected new-context labels and final naked-DNA replicate 3
remained closed. The original scientific gates remain unchanged.

The fixed experiment compared signed beta-binomial depletion with a
zero-centered Normal(0,1) prior against the corresponding binomial fallback,
frozen count-tail scorers and separate dinucleotide weight matrix (DWM)
propensity controls. The 21-base log-linear bias model (LOG21) remained the
frozen five-seed mitochondrial ensemble. No bias refitting, geometry search,
new method selection or threshold relaxation was performed.

All twelve initial biological-library/depth endpoints ran: three libraries
per cell context, 10 million complete fragments and full depth, seed 2026.
The 8,106-site historical replay is nonnested. The new held-training-chromosome
cohort has 6,020 sites; all 21 tasks remain recorded, including unavailable
HepG2 HMGA1. Seven of 20 available cohorts satisfy the original 200-pair
eligibility rule; diagnostic averages are not promotion evidence.

## Results and decision

Across 20 available tasks, with equal chromosome, library and task weighting:

| Comparison at full depth | Mean area under receiver operating characteristic curve change | Mean relative average-precision change |
|---|---:|---:|
| LOG21 versus DWM propensities; same binomial detector | +0.022675 | +2.7036% |
| LOG21 binomial versus LOG21 count-tail | -0.037118 | -4.5703% |
| LOG21 beta-binomial versus LOG21 count-tail | -0.043105 | -5.5111% |
| LOG21 beta-binomial versus LOG21 binomial | -0.005987 | -0.9782% |

These are development ranking comparisons, not native end-to-end TOBIAS
results. All 60 available task/library full-depth panels remain imbalanced in
actual coverage after matching on archived covariates. A fitting-fold-only
covariate baseline has mean area under the receiver operating characteristic
curve of approximately 0.622802. Thus the old count-tail ranking advantage
does not isolate biological footprint discrimination. No retrospective
rematching was performed to improve these results.

Changing only the detector leaves identical-site observed, expected and
residual profiles unchanged. The native full-depth DWM residual/strand
comparator remains unqualified. Expanded native control support processed
all 4,200 windows, but two historical expected shapes depend on region layout;
the exact old layout reproduces all 281 overlapping historical profiles.

Seven of 21 exploratory LOG21 control central-mean intervals exclude zero,
without multiplicity correction. Each candidate fails 18 of 21 development
control rate/Wilson/paired-propensity checks. These are examined replicate-1
controls with dependence, not untouched safety evidence. Functional bands and
complete-workflow resources are not qualified.

The sole frozen advancement screen was not met. Do not expand to the full
depth/seed/resource matrices or add candidate methods. A negative result with
these confounding and calibration limitations is the intended deliverable.

## Verification and artifacts

All six chromatin fragment-sharing inventories completed; no complete
fragment linked training windows across fitting/evaluation chromosomes.
An independent coordinate oracle checked retained endpoints. Complete arrays
and per-site data are kept under ignored `benchmarks/results/`.

The portable score replay reproduced all 601,398 score/cohort rows exactly,
including all calls. The independent verifier reproduced 21,548 available
fold metrics, 3,528 control-rate rows, 420 standardized shape rows and 10,080
regional shape summaries. Maximum standardized-shape discrepancy was
1.88e-13. A separate full-depth control refit reproduced all 60 fitted
dispersion and threshold records exactly. Fifteen focused tests passed in
3.183 seconds. These checks establish reproducibility, not biological safety.

The delivery includes a two-slide editable PowerPoint and PDF preview,
complete derived headline arrays, source tables, runnable reanalysis,
read-only verification, provenance and preserved failed attempts. The final
ZIP receipt and fresh-extraction verifier output are stored beside the
archive. Only Linux/source execution and LibreOffice rendering were tested;
native Windows/macOS presentation rendering was not tested.
