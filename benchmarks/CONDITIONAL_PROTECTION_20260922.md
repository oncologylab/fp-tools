# Strand-conditioned protection score — 22 September 2026

Conditioning on strand-specific read counts improves training discrimination
more broadly than the unit normalization floor, but precision-recall gains,
canonical shape, control precision and detection power remain insufficient
for promotion. This is a research scorer, not a production option.

## Scoring operator and assumptions

For each fixed center/shoulder geometry, condition separately on the positive-
and negative-strand counts within the center plus both shoulders. Exclude the
gap and other flanks. Sum frozen sequence propensities to obtain each strand's
central probability, or use uniform propensities for the correction control.
Under independent cuts, central counts follow the sum of two binomials. The
score is the negative base-10 logarithm of the lower-tail probability at the
observed central count, capped at 300 and rounded to ten decimal places.
Zero conditioned counts yield zero score. Scores are not binding probabilities.

The probability formula assumes correct fixed propensities and independent
cuts within and between strands. Paired fragments, overdispersion and imperfect
bias predictions can violate these assumptions. Synthetic validity does not
establish biological calibration. This method changes ranking/calling, not
the underlying corrected signal or aggregate profile.

The same 11 geometries, minimum 200 fitting sites per class, chromosome/pair-
mate exclusions, protection rule and explicit fallbacks apply. Geometry
eligibility uses unit-floor residual contrasts with negative central signal
below both shoulders; eligible geometries are ranked by the new count score's
training discrimination. Uniform and parametric arms each have fixed and
selected-geometry versions, separating correction and detector effects.

## Training evidence

All 8,106 matched sites, 21 tasks and chromosomes 1–15 remain. All geometry
probabilities were finite; no support filtering was needed. The unavailable
HepG2 HMGA1 task remains explicit. Matching was frozen on the full training
cohort and is not nested inside chromosome cross-fitting. Metrics average
equally over evaluable held chromosomes rather than pooling scores from
differently fitted geometries.

AUROC means area under the receiver operating characteristic curve; AUPRC
means area under the precision-recall curve. Descriptive means across the
20 available tasks, relative to the unnormalized conventional dinucleotide
weight matrix (DWM) detector, are:

| Propensities / geometry | Mean AUROC gain | Mean per-task relative AUPRC gain |
| --- | ---: | ---: |
| Parametric / fixed | 0.03907 | 5.00% |
| Parametric / selected | 0.04252 | 5.82% |
| Uniform / fixed | 0.01496 | 0.55% |
| Uniform / selected | 0.01858 | 1.36% |

The parametric AUROC mean exceeds the numerical 0.03 target descriptively;
the AUPRC mean does not reach 10%. Neither is locked promotion evidence.
Selected parametric examples are:

| Task | AUROC | Gain over conventional DWM | Relative AUPRC gain | Fixed-center protection folds |
| --- | ---: | ---: | ---: | ---: |
| HepG2 CCCTC-binding factor (CTCF) | 0.60609 | 0.02678 | 6.39% | 10/15 |
| K562 CTCF | 0.67522 | 0.08454 | 20.54% | 12/15 |
| HepG2 MEF2D | 0.57211 | 0.09164 | 10.99% | 5/15 |
| K562 MEF2D | 0.60128 | 0.11483 | 15.32% | 5/15 |
| HepG2 ZNF362 | 0.57034 | 0.04273 | 4.24% | 5/15 |
| HepG2 ZNF384 | 0.57392 | 0.02554 | 7.38% | 6/15 |
| K562 ZNF384 | 0.53745 | 0.05108 | 0.30% | 3/15 |

HepG2 MEF2D has only 164 matched pairs and uses the underpowered fallback.
Fixed-center protection tests plus/minus five bases against shoulders at
15–40 bases; selected-geometry protection is also recorded separately.
These numerical gains do not establish a difficult-family footprint rescue.
Of 630 fits, 105 select a protective geometry, 350 are underpowered and 175
have no eligible protective contrast.

## Replicate-1 null audit and power

All 42 full-training model records were frozen before control scoring. The
nominal threshold was fixed at rounded `-log10(0.05) = 1.3010299957`, with
strict greater-than calls and positive conditioned counts required. It was
not calibrated or adjusted to the control data. Replicate 1 is already
examined development data, not independent safety evidence. All 4,200
control sites remain in each arm's ledger, including 35 blacklist overlaps
and unavailable models.

There are zero calls among 1,471 informative parametric task-sites and 1,441
uniform task-sites. These task-sites share a BAM and can share windows;
pooled counts do not prove general safety. Only 9/20 available parametric
tasks and 8/20 uniform tasks meet both all-eligible and informative-site
Wilson upper bounds of at most 5%. Parametric examples:

| Task | Calls / informative sites | Wilson upper |
| --- | ---: | ---: |
| HepG2 CTCF | 0 / 72 | 5.07% — insufficient precision |
| K562 CTCF | 0 / 95 | 3.89% |
| K562 MEF2D | 0 / 101 | 3.66% |
| K562 ZNF384 | 0 / 63 | 5.75% — insufficient precision |

Wilson limits are descriptive binomial intervals, not validated chromosome-
block uncertainty for dependent windows. No paired-DWM safety claim is made.
The maximum observed control score is only 0.6818 for parametric and 0.9787
for uniform propensities, below the nominal threshold.

The same fixed threshold applied to held training scores shows low sensitivity
at the existing development depth: parametric selected geometry calls 105/4,053
positive-label sites (2.59%) and 28/4,053 negatives (0.69%). CTCF positives
contribute 50/300 HepG2 and 46/245 K562 calls. Both MEF2D tasks and the examined
ZNF362/ZNF384 tasks have zero positive-label calls. Labels come from chromatin
immunoprecipitation (ChIP); a label is not proof of an observable accessibility
footprint. Conservatism must not be mistaken for useful detection power.

## Verification and artifacts

Ten focused tests passed in 2.84 seconds, including exhaustive small-count
probability enumeration, nominal-tail control under the stated independent
null, rounded-threshold behavior, strand exchange, region accounting, zero
counts, invalid support and geometry fallbacks. Independent direct probability
convolution exactly reproduced all 178,332 training grid scores with zero
numerical discrepancy; all 630 fits, 1,260 folds and 84 summaries were checked.
The control verifier independently checked full-model choices, all 8,400 site
rows, probabilities, counts, exclusions, calls and 42 rate/Wilson summaries.
There is no new independent BAM extraction in this experiment.

Training artifacts: `benchmarks/results/conditional_protection_20260922_v1/`.
Control artifacts: `benchmarks/results/conditional_protection_rep1_20260922_v1/`.

| Artifact | SHA-256 |
| --- | --- |
| Training configuration | `cbf826133caa4d7af4ff0bd8f8cd6fc60f957aff878b7fee662249f5e11d3f37` |
| Training result | `ea02b6bdab9550d5e6604273bb171f5f5cf0f5685f867fe91c62b18023a764a1` |
| Control configuration | `3c4a6f953d2911c3cfd0727252888a11084391e5320b8d7c91979239c7c065e2` |
| Full models | `4e1aeded40b4bf84629c1b9ce0382c682642554ddc2c5181d3fd0a73610a7592` |
| Pre-evaluation control lock | `c11eb7e1d30f29186781fe204383d89d6140283831148c368c7ee6d6c76498ea` |
| Control result | `2b008a855ecabfd7f732db28e981df884022531ec9baa66f464f6eba66fa37ae` |

Both jobs have immutable launch/exit receipts under their corresponding
`conditional_protection*_jobs_20260922/execute_attempt1/` directories.
Training execution took 60.653 seconds and the model/control audit 13.891
seconds. These are diagnostic elapsed times, not complete-workflow central
processing unit (CPU) or memory measurements.

Next freeze an analytic score-ceiling audit: at each observed strand total,
compute the most extreme possible central depletion (zero central cuts) and
determine whether the unchanged threshold is reachable. This is a conditional
information limit, not a power estimate under a biological alternative. Use
it to decide whether more permitted development depth is informative, rather
than lowering thresholds to manufacture calls. Shape, safety, uncertainty,
replication and complete-workflow resource gates remain unresolved. Main,
SK-N-SH/GM23338 labels and final naked-DNA replicate 3 remain unchanged.
