# fp-tools Development Plan

Last updated: 2026-09-23

## Next-major-version research goal

**Bounded next-round experiments and verified evidence delivery completed:**
`benchmarks/scripts/next_round/` implements the requested fixed-geometry,
zero-centred binomial and beta-binomial shrinkage comparison. The experiment
policy and permitted input hashes were frozen before new scoring. Independent
quadrature fixtures (40 model/count combinations), exact zero-information and
support checks, calibration order-statistic checks, and independent CIGAR/pair
fixtures passed. All 21 tasks remain in the ledger. The original 8,106-site
replay is explicitly nonnested; a separately frozen training-chromosome
matching procedure retained 6,020 held sites. New fragment accounting and
same-source depth arrays completed for all six chromatin libraries and the
replicate-1 naked-DNA control. All twelve initial 10M/full scoring endpoints
completed. Fifteen focused regression tests passed in 3.183 s. Independent
summary checks reproduced 3,528 control-rate rows, 420 shape rows (maximum
absolute error 1.88e-13), and 147 full-screen task/comparator records.
On 20 available tasks, both shrinkage candidates lose ranking performance
against the existing LOG21 count-tail detector at full depth; the original
21-task denominator and missing HMGA1 remain explicit. Do not advance the
expensive depth/seed or resource matrices. Native full-depth DWM residual and
strand qualification, functional uncertainty, and final controls remain
incomplete or untested. All six fragment-sharing inventories completed with
no shared training windows across fitting/evaluation chromosomes. The portable
replay reproduced 601,398 score/cohort rows with zero score differences and
zero call disagreements. A separate control replay reproduced all 60 full-depth
dispersion/threshold records. The final ZIP was safely extracted into a new
directory and its bundled read-only verifier exited zero, checking all 1,016
manifested payloads, all 21 task records, 30,870 method/run/fold records and the
scientific summaries. The archive contains an editable two-slide PowerPoint
and visually checked PDF preview. Scientific qualification remains partial;
no method is promoted. Archive identity and checks are recorded in
`benchmarks/NEXT_ROUND_20260922.md`.

The frozen LOG21 bias compared with DWM propensities under the same binomial
detector gains mean AUROC 0.022675 and relative average precision 2.7036% across
20 available tasks. This is a development factorial comparison, not a native
TOBIAS benchmark or the qualified original promotion comparator. Both new
shrinkage scorers rank below the existing LOG21 count-tail scorer. Actual
full-depth coverage is imbalanced in all 60 available task/library panels,
and a covariate-only baseline is discriminative: these rank losses do not prove
intrinsic biological inferiority. The cohort was not rematched after scoring.
The two-slide PowerPoint preserves all task results, these qualifications and
the no-advancement decision. See `benchmarks/NEXT_ROUND_20260922.md`.

A September 23 presentation supplement adds the explicit LOG21 feature and
likelihood diagram, fixed-region schematic, and receiver operating
characteristic/precision–recall curves. The four arms share identical held
sites; 2,520 fold/method metrics reproduce exactly. Curves retain equal
chromosome/library/task weighting over 20 available tasks and explicitly
record the unavailable 21st task. Stepwise precision–recall averaging avoids
the area distortion of straight-line interpolation. The supplement changes
no scientific outputs or archived evidence.

For this round, SK-N-SH/GM23338 labels and final naked-DNA replicate 3 stay
closed regardless of development results. No main, release, or manuscript
changes are authorized by this investigation.

**Higher-count replicate-1 calibration and count screen completed:** independent
verification reproduced all 50,400 regional rows/378 summaries, then all 42
frozen models, 8,400 screen rows and 42 rate/Wilson summaries with zero score
discrepancy. Twelve focused tests passed in 3.85 s. Parametric calls are 2/3,343
informative sites; 20/20 available tasks meet both descriptive bounds, but only
65/3,967 eligible sites can reach threshold even with zero central cuts.
Zinc-finger central counts agree broadly with sequence expectations, while
seven pooled central-count discrepancies remain. Both jobs exited 0; no rerun
is needed. See `benchmarks/FULL_REP1_REGION_CALIBRATION_20260922.md` and
`benchmarks/FULL_REP1_COUNT_SCREEN_20260922.md`. No independent safety,
depth-only effect or promotion is established. Next freeze higher-depth
training-input evaluation of unchanged fold-specific detectors, with identical
sites and coverage/motif/accessibility balance audits. Main and untouched
controls/labels remain unchanged.

**Higher-depth replicate-1 input audit completed:** the final filtered alignment
and index pass quick checks, fresh full-file statistics and before/after input
hash checks. It contains 78,358,765 mapped records versus 14,957,194 in the
first-ten-million-pair alignment. Historical summary fields have different
counting definitions and must not supply depth denominators. Full preparation
also applies read-level blacklist and mitochondrial filtering; the files are
not a controlled depth-only comparison. See
`benchmarks/REP1_DEPTH_INPUT_AUDIT_20260922.md`. Next freeze higher-count
development calibration at the same sites with unchanged model, regions and
thresholds, explicitly retaining the preparation limitation. No new profiles
or scientific scores were opened in this input audit.

**Fixed-region replicate-1 calibration completed:** all 21 tasks and 4,200
original sites are retained. Independent verification reproduced 50,400
site-region rows, 378 calibration summaries, 63 likelihood comparisons and
5,670 chromosome blocks; seven focused tests passed in 4.10 s. Zinc-finger
central counts agree broadly with low parametric expectations, supporting
sequence-related cleavage deficits rather than uniform cleavage. Exploratory
central-count discrepancies remain for K562 CTCF, FOXA1 and MEF2D and HepG2
HMGA1. These unadjusted, sparse-count development results do not establish
general calibration, biological uncertainty coverage or promotion safety.
See `benchmarks/REP1_REGION_CALIBRATION_20260922.md` for frozen hashes and all
limitations. Next audit higher-depth replicate-1 provenance before freezing
new extraction; larger intermediates alone are not verified final inputs.
Main, untouched labels and final naked-DNA controls remain unchanged.

**Conditional score ceiling explains much of low sensitivity:** the frozen
audit retains all 40,824 rows and 210 task/group summaries. Seven focused tests
passed in 3.69 s; independent probability-product and convolution checks
verified all rows/summaries with zero observed-score discrepancy. Of 4,053
positive-label training sites, 3,589 (88.6%) cannot reach the unchanged nominal
threshold under the selected parametric model even with zero central cuts.
Only 464 can reach it, of which 105 are called. All positive ZNF362/ZNF384
sites are unreachable in the examined tasks. With fixed geometry and identical
counts, their parametric expected central probability is about 3% versus
20.55% under uniform propensities. This is model-conditional, not proof of an
absence of footprint information or a sequencing-depth recommendation.
Only two eligible parametric replicate-1 control sites can reach threshold;
none can in the uniform arm, limiting what zero control calls establish.
See `benchmarks/CONDITIONAL_SCORE_CEILING_20260922.md`. Next freeze a fixed-geometry
observed-versus-expected central-count calibration audit on permitted naked-DNA
replicate 1 before attributing the zinc-finger limit entirely to depth. Keep
strand/shoulder effects separate and use all tasks; do not select an analysis
only for the weakest family or change thresholds. Further depth work must use
verified BAM identities and training coordinates, never examined test labels.

**Count-conditioned candidate improves ranks but remains power-limited:** a
strand-conditioned binomial depletion score was frozen and tested with uniform
and fixed parametric propensities on the unchanged 8,106 matched training sites.
Ten focused tests passed in 2.84 s, including exact small-count null enumeration.
Independent convolution exactly reproduced all 178,332 grid scores, 630 fits,
1,260 folds and 84 summaries. Conditional-parametric selected geometry has
descriptive mean AUROC gain 0.04252 and relative AUPRC gain 5.82% over conventional
DWM across 20 available tasks; shape remains inconsistent and no promotion
follows. Frozen full-training models give zero nominal calls on examined
replicate 1, but only 9/20 available parametric tasks meet the informative
Wilson bound because of limited counts. Independent audit verified all 42
models, 8,400 control rows and rate/bound summaries. At the nominal threshold,
held training positives yield only 105/4,053 calls (negatives 28/4,053); MEF2D
and examined zinc-finger tasks have no positive-label calls. See
`benchmarks/CONDITIONAL_PROTECTION_20260922.md`. Next freeze the analytic
count-conditioned score ceiling (zero central cuts at fixed strand totals)
to identify sites incapable of reaching the unchanged threshold. Use that
diagnosis to decide whether additional permitted development depth can help;
do not lower thresholds or claim zero calls prove broad safety. Main, untouched
holdouts, final controls and all prior failures remain unchanged.

**Unit-floor candidate does not broadly rescue performance:** frozen training
cross-fitting retains all 21 tasks/seven correction arms and produces 2,205
fits, 4,410 folds and 294 summaries. Conditional-parametric mean gains over
conventional DWM are 0.01558 AUROC and 3.30% relative AUPRC across 20 available
tasks, below targets and not promotion evidence. CTCF remains strongest;
K562 ZNF384 loses 0.05518 AUROC against the former normalized detector and
fixed-center protection remains weak in difficult families. Nine tests passed
in 2.04 s. Independent verification reproduced all geometry choices and fold
metrics; one score differs by 1e-10 under equivalent solvers, so exact score
parity is explicitly false. Failed verification snapshots are preserved.
See `benchmarks/UNIT_FLOOR_GEOMETRY_20260922.md`. Do not advance this as a general
replacement or tune per-task floors. Next specify and freeze a count-conditioned
central-versus-shoulder scoring candidate using training data, with explicit
strand counts and sequence probabilities, before evaluating it. Assess whether
it can suppress low-count noise without discarding real protection; any null
calibration claim requires its assumptions and exact operator to be checked.
Untouched controls, promotion labels, main and all prior failures remain intact.

**Low-count normalization diagnostic complete:** all 4,200 replicate-1 control
sites are on training chromosomes; the membership ledger retains 35 blacklist
overlaps, 1,569 zero-count sites and 200 unavailable-model sites (categories
overlap). Existing thresholds evaluate 2,497 informative eligible sites per
available arm. Conditional parametric scores exceed threshold at 100 sites;
CTCF contributes 12/120 HepG2 and 14/123 K562, and K562 MEF2D 11/134.
No informative parametric profile uses the zero-variation fallback; minimum
scale is 0.034, far above 1e-8. Called profiles have median normalization
multiplier 13.61 versus 4.35 for noncalled profiles. This is a descriptive
association, not proof of causality or a safety pass. Nine focused tests passed
in 2.23 s; independent least-squares verification checked all 24,612 rows,
scores and threshold decisions, plus equivalent expected-profile construction.
See `benchmarks/LOW_COUNT_NORMALIZATION_20260922.md`. Next freeze a training-only
comparison of a fixed unit scale floor against current normalization, using
the same chromosome cross-fitting, matching, geometry and protection rule.
Do not tune floors per task or against replicate 2, and do not change existing
models/thresholds. This tests a concrete non-amplifying candidate; it does not
relax any scientific, safety or full-workflow resource gate.

**Blacklist and output-layout effects separated:** the all-site blacklist
intervention leaves every eligible observed profile and expected shape
unchanged, removes signal at 12 blacklisted sites, and retains 2,880 positive
expected profiles among all 4,200 rows. All 269 historical positive windows
remain, with the same three layout-induced shape differences. Twelve focused
tests passed in 2.44 s; independent interval reconstruction verified 25,200
profiles and both comparison masks. See `benchmarks/BLACKLISTED_CONTROL_LAYOUTS_20260922.md`.
The expanded reference remains unqualified; do not waive its overlap gate or
retune candidates on replicate 2. The exclusion/layout diagnosis is now
resolved sufficiently to move to the outstanding low-count normalization
diagnostic using training and replicate-1 controls. Freeze counts, normalization
components and fixed-model outputs before inspecting their relationships;
preserve all failed operators. Full-workflow resource and promotion gates
remain outstanding, and untouched labels/final controls stay unopened.

**Blacklist intervention resolves historical coverage discrepancies:** adding
the checksum-locked hg38 blacklist to the matched historical-region command
reproduces all 269 positive expected windows, with zero support, expected-shape
or observed-shape discrepancies. All 4,200 sites remain, including 35 explicitly
blacklisted motifs/windows. Eleven focused tests passed in 2.11 s; independent
interval reconstruction verified 8,400 profiles and all exclusion flags.
See `benchmarks/CONTROL_BLACKLIST_20260922.md`. This establishes a sufficient
reconstruction, not complete provenance for the archived command. Next freeze
a blacklist-aware all-site counterpart to distinguish exclusion effects from
output-layout effects. The three existing layout differences remain real;
preserve the failed expanded-reference qualification and do not rescore or
retune failed candidates. Replicate 2 remains examined development evidence.

**Current checkpoint — historical-region isolation verified:** all 269 archived
positive expected profiles reproduce their shape with the historical output
layout. Expanding the layout changes three shapes under otherwise identical
settings. Two additional historical coverage discrepancies overlap the verified
hg38 blacklist; this is a provenance clue, not proof of the original command.
Nine focused tests passed in 2.00 s and immutable artifact verification passed.
See `benchmarks/CONTROL_REGION_LAYOUT_20260922.md`. Next freeze a matched
blacklist-aware reconstruction, retaining excluded-site denominators and all
previous artifacts. Do not rescore against the unqualified expanded reference.
After reference qualification, diagnose low-count normalization on training
and replicate-1 controls, then require exact-operator safety before broader
candidate evaluation. Complete-workflow CPU/memory gates remain unmeasured.
Main and untouched promotion controls remain outside this work.

**All-site control reference constructed but overlap qualification failed:**
all 4,200 observed profiles match corrected BAM extraction; positive expected
coverage increases from 269 to 2,892 windows. Three old positive windows
(two CTCF, one FOXA1) change expected shape where merged output regions change
the correction-window layout. Phase and boundary context have not been
separately isolated. Eight tests passed in 1.91 s; model
count arrays are unchanged. Do not rescreen using the unqualified reference.
See `benchmarks/ALL_SITE_CONTROL_REFERENCE_20260922.md`. Next freeze a
historical-region-only rerun with the same model, implementation and disabled
global scaling to isolate region-layout effects. Preserve all failed
operators, thresholds and prior references. No promotion gate is passed.

**Negative-control DWM scope limitation verified:** all 17,311 expected-track
bases fall within historical CTCF regions plus/minus 100. Of 4,200 control
windows, 3,931 have no track entries; positive expectation exists in 267 CTCF
and only two non-CTCF windows. One non-CTCF window has mismatched observed
shape. Archived coordinate-valid flags do not establish expected support.
Twelve tests passed in 1.42 s; all site classifications, summaries and track
containment independently checked. See `benchmarks/CONTROL_TRACK_SCOPE_20260922.md`.
Next freeze an all-site control-reference construction with the recorded
inputs and existing DWM bias model, requiring overlap parity before any
rescreening. Do not alter failed candidates or thresholds. Low-count
normalization diagnosis on training/replicate-1 controls remains outstanding.

**Full stable models and exact-operator safety screen complete:** 147 model
records frozen; 60 calibrated/scored combinations, 80 missing wider control
profiles and seven missing models retained. Corrected CTCF fails replicate-2
safety: normalized DWM calls about 21% of informative sites and conditional
parametric calls 9–10%. K562 MEF2D/parametric has Wilson upper 8.75%, also
failing; its DWM control has no informative support. Only raw CTCF in both
cells and a raw K562 ARID3A fallback pass the narrow screen; ARID3A has only
one paired DWM site. No promotion follows. Fifteen tests passed in 1.50 s;
all model, calibration and control rows independently verified. A failed
data-root resolution attempt was preserved and corrected in a separate v2
runner; original data hashes match. See `benchmarks/STABLE_GEOMETRY_SAFETY_20260922.md`.
Next audit normalization/depth behavior on training and replicate-1 controls,
and missing expected coverage in archived DWM control profiles. Do not retune
thresholds against observed replicate-2 failures. Main, final controls and
untouched promotion labels remain unchanged.

**Stable detector chromosome-direction audit verified:** all 8,820 paired
effects, 2,205 fixed-center audits and 588 summaries independently checked;
12 tests passed in 2.60 s and all three 21-task maps visually reviewed.
K562 MEF2D/DWM improves discrimination on 12/15 chromosomes, but fixed-center
protection passes only 4/15; do not call it an established footprint gain.
CTCF has stronger agreement; zinc-finger gains remain inconsistent or worse.
See `benchmarks/STABLE_GEOMETRY_DIRECTION_20260922.md`. Next freeze the
unchanged full-training geometry rule and exact model artifacts, preserve all
tasks/fallbacks, then audit their exact operators on permitted development
negative controls. Interpret the narrow MEF2D shape explicitly; do not borrow
safety from an earlier operator. Main and untouched promotion labels remain
unchanged. Full-workflow resource qualification remains incomplete.

**Stable numerical successor verified:** score rounding, exact rank-sum
selection and a numerical protection margin were frozen before the repeat.
All 2,205 fits and saved held scores agree under equivalent normalization
solvers; zero fold metrics change beyond arithmetic roundoff (previously 287).
Ten tests passed in 3.15 s; full independent verification covers 6,615 fold
rows and 441 summaries. Geometry helps CTCF and K562 MEF2D diagnostically,
but does not rescue the zinc-finger tasks. No scientific promotion gate is
newly passed. See `benchmarks/STABLE_GEOMETRY_CROSSFIT_20260922.md`.
Next quantify chromosome-wise stability and signed protection of this exact
detector before any full-training candidate freeze or negative-control test.

**Balanced-cohort geometry cross-fitting complete, numerical limitation found:**
all 2,205 fits, 6,615 fold rows and 441 summaries verified; seven tests passed
in 2.96 s. Equivalent normalization perturbs scores by at most 2.91e-13 but
changes rank metrics in 287 rows (maximum chromosome AUROC change 0.175).
Do not accept this as a numerically qualified candidate evaluation. Preserve
the run; next freeze a numerical tie policy for fitting and scoring, with
sparse/flat-profile, batching and equivalent-solver checks, before a successor
run. Do not select numerical tolerances to improve biological results.
Details: `benchmarks/BALANCED_GEOMETRY_CROSSFIT_20260922.md`. No untouched
labels were opened, no promotion gate passed, and main remains unchanged.

**Training chromosome-block sensitivity verified:** 2,000 fixed draws per
cell across all 147 rows. HepG2 ZNF362/ZNF384 discrimination intervals exclude
zero, but central-protection intervals cross zero; K562 ZNF384 discrimination
also remains uncertain. No general footprint improvement or promotion gate
is established. Four tests passed in 2.52 s; artifact verification, all 147
interval/count checks and 420 independent explicit-replication score/shape
checks passed. See `benchmarks/RADIUS_TRAINING_UNCERTAINTY_20260922.md`.
Next diagnose weak signed protection on balanced training cohorts before
adding model complexity; freeze any new detector and test its exact operator
on permitted negative controls. Full-workflow resource qualification remains
incomplete. Main and unopened promotion labels remain untouched.

**Full-pool seven-arm comparison verified and reviewed:** all 147 task/method
rows retained; 140 score/curve rows independently checked against saved panels
and seven unavailable rows confirmed. All 21 aggregate panels visually checked.
Radius 500 retains HepG2 ZNF362/ZNF384 AUROC gains of 0.0731/0.0612, but K562
ZNF384 gains only 0.0170. Their functional-separation gains are 5.20%, 6.25%
and 1.53%, far below 25%. Do not advance the unchanged policy to untouched
promotion labels. Eight tests passed in 2.45 s; separate artifact verification
passed. Exact evidence is in `benchmarks/RADIUS_FULL_POOL_COMPARISON_20260922.md`.
Next freeze chromosome-block uncertainty diagnostics for balanced training
scores and signed shape; they cannot substitute for untouched validation.

**Direct three-covariate cohorts complete and verified:** matching the complete
supported training pool restores 4,053 balanced pairs and ten powered tasks.
HepG2/K562 CTCF retain 300/245 pairs, HepG2 ZNF362 291, and HepG2/K562 ZNF384
303/238. MEF2A remains underpowered (190/169). All 21 tasks remain explicit.
Fourteen tests passed in 1.43 s; separate artifact verification and independent
sparse-graph maximum-cardinality, pair/caliper/support/balance checks passed.
No footprint scores or shape outcomes were used for construction. Next freeze
and recompute all seven fixed arms and exact-cohort aggregates; do not reuse
earlier metrics or curves as if site membership were unchanged. Evidence is in
`benchmarks/RADIUS_FULL_POOL_COHORTS_20260922.md`.

**Secondary broad-count balance sensitivity completed:** the covariate-only
frozen prefix rule retains 2,517/3,949 pairs, leaving only three powered tasks,
17 underpowered tasks and one unavailable task. All 147 method/task rows remain
explicit. HepG2 ZNF384 weakens (radius-500 AUPRC 0.5554 versus DWM 0.5837);
K562 ZNF384 retains a numerical gain with only 92 pairs. HepG2 ZNF362 keeps
all 256 original pairs and its earlier result. No secondary shape claim is
made. Six tests passed, both artifact verifiers passed, and independent pair,
balance and rank-metric checks passed. Exact evidence is in
`benchmarks/RADIUS_SECONDARY_BALANCE_20260922.md`. Next freeze direct
three-covariate matching from the full existing training pool to address power
without selecting on scores; preserve both earlier analyses and all tasks.

**Radius accessibility/control audit complete:** all original pairs and tasks
retained. Counts at offsets 101–500 are imbalanced for CTCF (SMD −0.746/−0.522)
and ZNF384 (−0.210/−0.225), despite balanced central counts. The parametric
radius-500 arm still exceeds its uniform-sequence control by AUROC 0.0397
for HepG2 ZNF362, 0.0567 for HepG2 ZNF384 and 0.0470 for K562 ZNF384.
This does not eliminate broader-accessibility confounding. Seven relevant
tests passed; separate artifact verification and independent recomputation
of 105 associations, 21 metrics and 735 correlations passed. No candidate
selection, rematching or holdout evaluation occurred. Next freeze a secondary
broad-count-balanced analysis, selecting pairs only by covariates, retaining
the primary evidence and reporting lost power. Exact evidence is in
`benchmarks/RADIUS_ACCESSIBILITY_AUDIT_20260922.md`.

**Radius comparison completed and reviewed:** all 126 metric rows independently
recomputed, all 21 aggregate panels visually checked, all 3,949 original pairs
retained. Radius 500 shows promising training discrimination and canonical
point contrasts for HepG2 ZNF362/ZNF384 and K562 ZNF384. It is not promotion
evidence: broader accessibility remains unaudited, functional-separation gains
are below 25%, multiple difficult families are not established, and safety/
uncertainty/stability/resource gates remain open. CTCF remains below the earlier
conditional parametric comparator. The 100-radius implementation's maximum
score difference is 1.547e-13 with no changed AUROC/AUPRC. Next freeze a
training-only broad-accessibility and uniform-sequence comparator audit before
candidate selection. See `benchmarks/RADIUS_TRAINING_COMPARISON_20260922.md`.

**Wider extraction and independent audit complete:** both cells finished
in 641.24 s; separate extraction and audit verification passed. All 126
diagnostic rows were independently recomputed. Empty strand backgrounds
decrease markedly with wider flanks, but CTCF left/right rate disagreement
persists. This is label-free training evidence, not better discrimination or
uncertainty calibration. Exact results and hashes are in
`benchmarks/WIDE_BACKGROUND_DIAGNOSTIC_20260922.md`.

**Historical fixed-radius comparison launch (now verified above):** the runner preserves the
three reference arms and compares radii 100/200/500 on complete matched pairs,
with fixed dispersion and geometry. It checks archived site hashes, observed
counts, motif/count balance and numeric/rank parity against the earlier
100-radius result. Fifteen relevant tests passed in 2.72 s, including a full
synthetic 21-task, 126-row execution. Configuration SHA-256 is
`c448b0a18ef1e4ce092eebdd4d6fe9a5861a5f389f13bdddc92c7f4c6db68beb`.
Execution evidence is under
`benchmarks/results/radius_training_jobs_20260922/execute_attempt1/`.
No scientific completion or candidate selection is implied by launch.

**Fixed-radius scoring primitive ready:** `radius_background.py` implements
the same half-count background model at prespecified outer radii 100, 200 and
500 with prefix sums. The fixed deviance dispersion and detector geometry
remain unchanged. Fourteen relevant tests passed in 1.62 s, including scalar
agreement, 100-radius parity, arbitrary log-bias offsets, complete scoring halo
parity, missing sequence and zero-count ties. The training comparison protocol
is `benchmarks/manifests/radius_training_comparison_v1.json`; its real-data
execution configuration is not yet frozen, and no radius has been evaluated
or selected by discrimination. Preserve all 21 tasks and complete common-
support pairs, audit any balance changes, and report floating-point rank
sensitivity against the earlier 100-radius result.

**HepG2 wider profiles independently verified:** the stage contains 20,019
sites, all alignment-valid, with shared central counts and bias values matching
the prior extraction. K562 extraction child 677733 was confirmed live at
7 min 30 s. Full extraction and common-support diagnostic results remain
pending; do not rerun the completed HepG2 stage.

**Independent wider-background auditor ready:** `audit_wide_background.py`
recalculates flank counts with prefix sums and exposures in log space, checks
every archived diagnostic row, and reports all 126 task/strand/radius rows on
common sequence support. It rejects changed coordinates, missing task rows and
changed checksums. Twenty-nine related tests passed in 2.41 s. The extraction
child 677733 was rechecked live at 3 min 41 s; no scientific result is inferred
from elapsed time. After extraction completes, run its independent
`--verify-only` before freezing and executing this audit. Left/right disagreement
is a descriptive diagnostic: it can reflect sampling noise, sequence-model
error or real regulatory heterogeneity, and is not a calibrated uncertainty
measure. Overlapping site windows are not independent biological replicates.

**Wider-context training diagnostic launched:** a new label-free extraction
uses unchanged training coordinates, alignments, bias model and cut convention,
with 600-base-pair context on each side. The frozen diagnostic reports
zero-background frequency and exposure-adjusted left/right rate disagreement
for outer radii 100, 200 and 500 (inner radius 50), retaining all tasks.
No radius is selected and no labels or footprint scores are evaluated.
Twelve related tests passed in 1.53 s. Execution configuration SHA-256 is
`50483e6776d99a162a84b49aece5db1ff78abc5ea14fad32e468f93542ae7937`.
Durable launch/exit evidence is under
`benchmarks/results/wide_training_jobs_20260922/execute_attempt1/`;
supervisor 677716 and child 677733 (start ticks 119307628) were confirmed live
at launch. Recheck process identity and output verification before resuming.
The frozen protocol is `benchmarks/manifests/wide_background_training_v1.json`.
This tests background noise/heterogeneity, not scientific improvement or
resource qualification. A later accuracy evaluation needs a separate freeze.

**Gamma background training comparison complete; do not advance this policy:**
all 3,949 matched pairs across 21 tasks remained available in five correction/
residual arms with fixed conventional geometry. All 105 saved metric rows were
independently recomputed and all 21 aggregate plots visually reviewed.
The new predictive residual reduced CTCF training discrimination relative to
the existing conditional parametric residual in both cells. Full-denominator
mean AUROC gain over the matched DWM residual was only 0.00317 (descriptive,
including underpowered tasks). Coverage recovery is not accuracy improvement.
Thirteen focused tests passed in 2.51 s. No tuning/test/control evaluation or
full-resource qualification follows from this result. Preserve the negative
result and avoid spending untouched labels on this unchanged policy.
See `benchmarks/GAMMA_TRAINING_COMPARISON_20260922.md` for exact evidence.

**Low-count model simulation verified:** the separately frozen Gamma–Poisson
prototype completed all 50 synthetic cases (50,000 draws each). Nine focused
tests passed in 1.56 s, and independent artifact verification passed. Nominal
95% scalar predictive coverage was 94.964–100% under the Poisson working
model but fell to 15.018% in the strongest overdispersion stress case.
This is not calibrated biological or functional uncertainty. No real-data
scores, labels, controls, or promotion gates were evaluated by this stage.
The next step requires a separate execution freeze for identical-site
training-chromosome comparisons, retaining all 21 tasks and separating
expectation changes from residual and detector changes. Do not retune this
prior using examined tests or controls. See
`benchmarks/GAMMA_POISSON_SIMULATION_20260922.md`.

**Global training parity verified; coverage remains limiting:** both cells
completed with zero observed/expected/validity overlap conflicts and matching
scalar/halo checks. Independent profile and coverage verification passed.
Conventional scoring support is 6,453/20,019 HepG2 and 5,248/14,987 K562 sites;
complete 201-position support is only 3,861 and 3,246. Central exclusions are
empty strand-specific backgrounds, not missing sequence. All 21 tasks remain
explicit. No labels or real discrimination scores were evaluated. Preserve
these exclusions and freeze any revised low-count model separately; do not
reuse prior safety/accuracy claims. Exact evidence is in
`benchmarks/GLOBAL_PARAMETRIC_TRAINING_PARITY_20260922.md`.

**Global scoring primitive tested:** `global_parametric_scoring.py` converts
genomic expectations directly to signed-deviance residuals and the existing
conventional geometry, without window-total rescaling. A 140-base-pair halo
preserves complete count-to-footprint results across chunk boundaries for
Poisson and negative-binomial residuals. Missing expectation remains missing
evidence. Fifteen related tests passed (1.43 s). This is synthetic numerical
verification only; no real-data scores, safety, or full-pipeline claim exists.

**Historical extraction launch (now completed and verified above):** the execution configuration was frozen
at SHA-256 `d903eec175130993afe35e6f3c64701e25399d1e29f6d2b98f8750825849b328`.
Supervisor 620697 tracked child 620699 (start ticks 119031082). The job used
the existing 10-million-fragment replicate-1 BAMs and all recorded training
coordinates: 20,019 HepG2 and 14,987 K562. Its scope remains numerical/output
parity and coverage, with no labels read or discrimination scoring.

**Historical training-profile runner implementation:** the initial-stage runner
extracts 401-position strand counts and fresh sequence propensities on
label-free chr1–15 coordinates. It audits scalar/halo agreement, genomic
overlaps, archived observed counts, and missing coverage before any scoring.
Nine focused tests passed (1.50 s). At that stage the supervised input freeze
was active; numerical compatibility and coverage have since been established
as reported above.
Protocol, model sidecar, reference/index, source alignment/index, code, and
runtime identities are recorded before extraction. No tuning, control, or
test scoring is authorized by this initial stage.

**Corrected 96-run comparison complete:** all inputs/outputs were independently
verified, all 21 aggregate pages visually reviewed, and the unchanged ledger
calculation rebuilt for corrected counts. CTCF full-depth mean AUROC gains
versus sample-matched DWM are 0.111969 (HepG2) and 0.130023 (K562); extraction
changes each run's CTCF AUROC by less than 0.0008. Low-depth failures remain:
two HepG2 shape exceptions and a K562 replicate-3 10-million-read interval
crossing zero. All 21 tasks and 746 shape exceptions remain explicit. This is
post-hoc reassessment on already examined test chromosomes, not untouched
evidence or data for retuning. Forty-one relevant tests passed (4.04 s).
See `benchmarks/CORRECTED_COUNT_REASSESSMENT_20260922.md` for exact hashes,
the complete scope, comparisons, and outstanding gates.

**Globally consistent prototype:** a research-only primitive now assigns one
expectation per genomic position from same-strand counts and sequence bias in
the symmetric 50–100 base-pair outer flanks. Ten primitive/output-contract
tests passed (1.38 s), including scalar, halo-chunk, overlap, strand, and
missing-support checks. The initial training-only protocol is locked in
`global_parametric_track_development_v1.json`; no real-data evaluation or
complete output pipeline exists yet. Earlier conditional-profile accuracy,
safety, and timing claims do not transfer to this new operator.

**Corrected 96-run verification passed:** independent full input/output
verification exited 0 after 937.31 s. The corrected result manifest SHA-256
is `71235831ce3559b5e5c81a5b30e4564879cf7afb11593cc69c221e9eb3d8b979`.
The historical-versus-corrected comparison and aggregate rendering are now
running with their input/script identities frozen in a supervised job.

**Full-workflow output compatibility audited:** the current conditional
expectation cannot be flattened exactly into one genomic track. Training-only
overlap checks found expected conflicts in 1,498/1,558 HepG2 pairs and 777/814
K562 pairs, with zero observed-cut disagreements. Five tests and independent
real-result verification passed. This is window-specific normalization, not
a new extraction error. Do not silently average profiles or call site-array
timing a production-equivalent resource benchmark. A globally consistent
operator needs its own frozen scientific evaluation before equivalent-output
resource measurements. See `benchmarks/SITE_CONDITIONED_OUTPUT_CONTRACT_20260922.md`.

**Tuning sensitivity verified:** all 378 comparisons are complete with
independent result verification. Parametric trained geometry versus DWM
conventional scoring improves pooled CTCF AUROC by 0.04463 in HepG2 and 0.04887
in K562, but the K562 three-block descriptive interval crosses zero
[-0.04622, 0.06801]. K562 chr16 has only six positives and shows a loss; it
remains in the analysis. HepG2 gains are positive on all three chromosomes.
This does not establish calibrated uncertainty or a general footprint
improvement. See `benchmarks/SHAPE_CHROMOSOME_SENSITIVITY_20260922.md`.

**Corrected 96-run execution completed:** the supervised reassessment exited
0 after 5,738.28 s. Independent full input/output verification is running as
child PID 584758 (supervisor 584757); no final corrected scientific comparison
is claimed before that check finishes. Verification specification SHA-256:
`c86e2fbf5110b4814f6cfe320a4d11e48a373202147f0e16c1cce46db3af2515`.

**Tuning chromosome sensitivity:** a new research-only helper preserves all
189 metric rows against two explicit baselines, enumerates all 27 ordered
three-chromosome resamples, and records missing-class draws without dropping
them. Five focused tests passed (1.95 s). Configuration is frozen at
`69027a08d0d7892e394fdbf4abce0930052d3af4f6cf870d198bdba34fbba2d2`;
execution is supervised. This is post-hoc tuning sensitivity, not calibrated
uncertainty, model selection, or untouched promotion evidence.

**Exact spline operator fails the control screen:** all 126 ledger rows are
complete and independently verified. Of 22 fitted operator rows, six pass
only the previously examined replicate-2 FPR screen; 104 missing-model rows
remain explicit. Parametric CTCF splines call 5/122 HepG2 and 5/115 K562
informative control sites, with Wilson upper bounds 9.24% and 9.78% and paired
increases above one percentage point. Their tuning gains therefore do not
establish a safe improvement. Trained CTCF geometries pass this narrow screen;
no difficult-family rescue or promotion is established. Preserve all models
and thresholds without control-driven retuning or retrospective selection.
See `benchmarks/FIXED_SHAPE_OPERATOR_SCREEN_20260922.md` for exact hashes,
support, validation, and limitations. Main remains clean at `636db5f`.

**Fixed-operator screen implemented:** `screen_fixed_shape_operators.py`
preserves all 126 task/operator rows, excludes missing expected coverage from
calibration, and freezes thresholds before the previously examined replicate-2
screen. Missing models, insufficient calibration, and unavailable baseline
support cannot pass. Zero-count sites are never called. Eighteen focused
screen/normalization/safety tests passed (1.61 s). The execution configuration
is frozen at SHA-256
`1074877deb1f4781e7e6115c5325ff02bfb6cd8c17f86dd5def004bfdfcfe5d5`;
supervised calibration is running. No control result is claimed yet.

**Fixed-normalization tuning verified:** corrected artifacts were reassessed
on identical tuning sites with all 189 metric rows preserved. DWM and
parametric spline metrics are unchanged; raw-spline metrics change slightly
(K562 CTCF AUROC 0.639941 to 0.639254). All corrected spline scores pass the
singleton/batch check. No refit or candidate selection occurred. Nine related
tests and independent real-result verification passed. The next operator
safety protocol is locked in `fixed_shape_operator_safety_v1.json`: evaluate
all fitted candidates, use common positive calibration expectation, and freeze
thresholds before examining replicate-2 scores. Its execution remains pending.

**Prediction consistency correction:** before screening the new operators on
controls, a training-only audit found batch-dependent raw-spline predictions.
The cause was recomputing a zero-outer-variation fallback from the prediction
batch. A research-only helper now freezes that fallback from training data.
All 11 corrected copies preserve learned weights and full-training scores,
and pass batch/singleton agreement checks. Original artifacts are unchanged.
Ten related tests and independent real-artifact verification passed. Exact
hashes and counts are in `benchmarks/CORRECTED_SHAPE_TUNING_20260922.md`.
Reassess the frozen copies on tuning data, then calibrate the new operators
on valid common expected support before the previously examined replicate-2
screen. Neither this consistency fix nor prior direct-operator safety certifies
the newly learned geometries/spline scores.

**Corrected tuning plot review completed:** all 21 tasks are rendered on seven
pages, with underpowered tasks and lack of confidence bands explicit. All seven
page images were visually reviewed and eight output checksums verified. CTCF
has the clearest depletion and shoulders; no general difficult-family rescue
or stronger parametric protection over DWM is established. Three rendering
tests and twenty documentation contracts passed. Exact plot hashes and
interpretation are recorded in `benchmarks/CORRECTED_SHAPE_TUNING_20260922.md`.

**Corrected candidate training/tuning completed:** 11 model artifacts were
fitted on chr1–15, frozen and verified before chr16–18 scoring. All 21 tasks
remain in 189 correction/detector rows. The parametric CTCF spline improves
tuning AUROC over the same correction's conventional geometry from 0.6377 to
0.7063 in HepG2 and 0.6445 to 0.6764 in K562. This is point-estimate development
discrimination, not evidence of stronger aggregate protection or promotion.
Difficult-family shapes/performance remain inadequate. No policy was selected
or refitted after tuning. Nine focused tests and independent result
verification passed. See `benchmarks/CORRECTED_SHAPE_TUNING_20260922.md` for
all scopes and exact hashes. Aggregate plotting/review, new-operator safety,
calibrated uncertainty and end-to-end resource measurements remain outstanding.

**Balanced cohort v3 verified:** a covariate-only rule takes the largest
balanced prefix of v2 pairs ordered by distance and then pair ID. This keeps
the 0.5 caliper, 0.1 balance bound, 200-per-class requirement, original expected
support and paired memberships. It is not a global optimal balanced-subset
claim and does not inspect any footprint score. Twenty-four related focused
tests passed in 1.52 seconds. The real execution exited 0 in 5.73 seconds and
independent verification passed. Freeze SHA-256:
`3295079e67d0681e8f8d9b4d2e25e85099bebaac2679f58020834eb5b903db46`.
Output: `benchmarks/results/corrected_shape_cohorts_20260922_v3/`; manifest
SHA-256 `9b19e5644b4c7cb60ee2d5477b156e927953a60810a0e87f236cd5ce19a568b3`.

Seventeen of 42 task/split rows are eligible and 25 underpowered. Seven tasks
have both eligible training and tuning cohorts: HepG2 ARID3A, CTCF, FOXA1,
ZNF362, ZNF384; K562 CTCF and ZNF384. CTCF train/tune pairs are 299/358 in
HepG2 and 280/253 in K562; ZNF384 pairs are 298/270 and 241/231. All other
tasks remain in the denominator. These are replicate-1, 10M, seed-2026
development cohorts only. Next freeze and fit detector artifacts on training
chromosomes, then evaluate tuning profiles without refitting. No candidate
performance, protection improvement, safety or resource pass is established.

**Caliper assignment diagnosis:** exhaustive synthetic checks show that forcing
a full minimum-cost assignment before filtering can discard every valid pair.
The new constrained matcher maximizes in-caliper cardinality, then minimizes
distance. Nineteen matching/cohort/candidate tests passed in 1.47 seconds.
Version 2 preserves version 1 and changes only this assignment objective; it
retains the 0.5 pair-distance limit, 0.1 absolute standardized mean-difference
limit, 200-per-class power requirement, original support and chromosome splits.
Its freeze SHA-256 is
`1b05603cc09b3c69f63c60af20e89c5a8b5fe3e753d5d35aeb99eaac421dff2a`.
Execution exited 0 in 8.24 seconds and independent verification passed.
Result: `benchmarks/results/corrected_shape_cohorts_20260922_v2/`, manifest
SHA-256 `39300d43e6cee4da314da08999729ffd6c48202dc83667957d42c6411cec023d`.

CTCF train/tune pair counts increase to 567/651 in HepG2 and 505/595 in K562,
but all four remain covariate-imbalanced. Across all 42 rows, 25 are unbalanced
and 17 underpowered; none is eligible. This diagnoses an assignment problem
and remaining systematic covariate imbalance, not absence of footprint signal.
No detector has been fitted or tuning score inspected. Balance must become an
explicit covariate-only subset-selection constraint before candidate fitting;
do not weaken the current gates or choose pairs using footprint outcomes.

**Corrected cohort preparation completed:** the first shared-support matching
protocol completed in 31.99 seconds and passed independent verification.
Configuration SHA-256:
`738101be0f93fe665cf628fb4edec5637d8aee302575a8bd3c72648ae9033e70`;
cohort manifest SHA-256:
`1593a42b7d9ca3ec90b6508a2dbb88627a0732481ff319da14350fe155e7ceef`.
Outputs are in `benchmarks/results/corrected_shape_cohorts_20260922_v1/`.
All 42 task/split rows are retained: 40 are underpowered and two unbalanced;
none is eligible for candidate fitting under this protocol. Ten focused
cohort/candidate tests passed in 1.45 seconds. No detector was fitted or scored.

For CTCF, common expected support remains large (HepG2 train/tune: 1,864/2,088
sites; K562: 1,857/2,089), but optimal full assignment followed by the fixed
distance caliper retains only 15/24 and 18/21 pairs. Motif score and observed
coverage differ substantially between classes. Diagnose the assignment and
caliper interaction before attributing the loss solely to absent biological
information; do not relax power/balance gates or overwrite this result.
The archived site metadata also labels chr16 as train. This new preparation
uses the requested chromosome coordinates directly (chr1–15 train; chr16–18
tune), rather than trusting that legacy column. Previously examined tuning
data are not described as untouched evidence.

**Next corrected-data candidate comparison:**
`benchmarks/manifests/corrected_shape_development_v1.json` fixes a training/tuning
factorial comparison of raw counts, DWM-expected residuals and frozen-parametric
residuals with conventional geometry, training-selected symmetric geometry and
the existing spline detector with family/global shrinkage. The new
`corrected_shape_candidates.py` fits only chr1–15 and rejects training templates
with central enrichment. Five focused tests passed, including rejection of
tuning/test inputs, preservation of underpowered tasks, family-prior isolation
and unchanged model coefficients during prediction. Matched cohort preparation,
execution/model freezes and the real tuning comparison remain to be implemented;
no new candidate-performance result is claimed from these synthetic tests.

**Verified corrected control result:** all four safety profiles and both
threshold analyses have completed and passed independent verification. The
fixed-threshold candidate call counts are unchanged across all 21 tasks.
CTCF control profiles and calls are unchanged, but the new calibration audit
finds zero baseline expectation in 9/989 HepG2 and 5/895 K562 calibration
negatives. The frozen protocol therefore withholds relative-safety passes;
numerically favorable CTCF control rates alone cannot resolve that gap. See
`benchmarks/CORRECTED_COUNT_SAFETY_20260922.md` for exact evidence and hashes.
The 96-run corrected-count reassessment is still running. No promotion or
new holdout access is justified.

**Corrected safety protocol:**
`benchmarks/manifests/corrected_count_safety_protocol_v1.json` fixes the next
negative-control reassessment before reading corrected control scores. Both
calibration and naked-DNA replicate-2 profiles used the affected extractor:
the archived inputs contain 36,931/26,647 development sites and 2,400/1,800
control sites for HepG2/K562, respectively. Sequence log-bias arrays are
available, so fresh observed counts can be combined with the same model
without fitting or new labels. Fixed historical thresholds isolate extraction
effects; independently recalibrated thresholds use only original tuning
negatives and the original rule. Both analyses are mandatory, with the latter
designated the primary development screen before outcomes are known.
Zero baseline expectation, profile disagreements and changed denominators
must remain explicit. This protocol is not an execution freeze or safety
result. `benchmarks/scripts/rebuild_corrected_safety_profiles.py` now provides
separate freeze, extraction and verification stages for the four archived
development/control inputs. It preserves site order, strand orientation,
original validity and fresh validity separately, and refuses changed inputs
or incomplete output reuse. Five focused tests passed, including invalid
sequence rows and corruption detection. The real input freeze/extraction and
corrected safety scorer remain outstanding.

The subsequent real four-input freeze completed successfully. Configuration
SHA-256: `dc053ee9d0e0d6154cabe0e3915a6a78ba39887ffe28e2850a5059db4aa3a9e9`.
Fresh extraction is now supervised with receipts under
`benchmarks/results/corrected_safety_jobs_20260922/execute_attempt1/`, with
that configuration hash checked before execution. Outputs are under
`benchmarks/results/corrected_safety_profiles_20260922_v1/`. Eighteen combined
safety-profile, corrected-count and extraction regression tests passed in
3.56 seconds. Completion and independent verification remain pending; these
profiles are not themselves safety scores or a passed safety gate.

`benchmarks/scripts/reassess_corrected_safety.py` implements the separate
comparison freeze, calibration-only and control-evaluation stages. The fixed
threshold analysis is explicitly diagnostic; the corrected calibration is the
prespecified primary development screen. It retains all sites, records changes
in baseline-comparable membership, reports residual diagnostics and curves,
and prevents unavailable calibration expectation from establishing relative
safety. Tests cover threshold freezing before control access, both analyses,
missing baseline support, invalid-site retention and output corruption/partial
output preservation. Real scoring awaits completed profile verification.

**Current execution, September 22:** the corrected-count freeze completed
successfully (exit 0; 761.65 seconds), with configuration SHA-256
`59dd9a562d31e49ebd7bd2303d1fc458412f8170dad174d7292696ccc13d7b32`.
The unchanged-candidate 96-run reassessment was launched at 08:06:38 UTC
under the durable experiment supervisor. Its command checks that exact
configuration hash before executing, and the runner checks frozen calculation
code and runtime. Job specification SHA-256:
`bcf41ee29fbdc7635648eab950e59ef2a7cf57f1d0ee6ac017a0a6a83bac8f74`.
Receipts are in
`benchmarks/results/corrected_count_jobs_20260922/execute_attempt1/`;
scientific outputs are in
`benchmarks/results/corrected_count_reassessment_20260922_v1/`.
Launch is not completion: require the exit receipt, independent
`--verify-only` verification, paired comparison, and aggregate-plot review
before interpreting results. This remains post-hoc implementation reassessment,
not candidate tuning or untouched validation. Main remains clean at `636db5f`.

**September 22 extraction qualification:** the research strand-profile reader
had duplicate-window and soft-clipping count defects. They are now fixed and
confirmed against an independent training-only oracle. Historical stability,
reference and safety computations remain preserved, but their scientific
interpretation requires corrected-count reassessment; checksums alone do not
establish extraction correctness. See
`benchmarks/STRAND_EXTRACTION_AUDIT_20260922.md` before reusing their conclusions.

Work exclusively in the research worktree on
`research/footprint-improvement-20260830`. Improve footprint detection and
interpretable protection across transcription factor (TF) families while
measuring complete-workflow CPU time and peak memory. Main, released packages,
the published manuscript and the production dinucleotide weight matrix (DWM)
default remain unchanged. No integration or publication is authorized here.

September 22 status audit: research HEAD is `8a4a958`, with 18 local commits
after remote research HEAD `888c6e2`. Main is clean at `636db5f`. No relevant
experiment process was found. The completed HepG2 ChromBPNet reference and
96-run study supersede the older queued-stage handoff. The three primary
96-run manifests and 198 referenced stage/summary files match their recorded
checksums. This is a manifest/stage-summary audit, not a fresh verification of
every nested BAM/profile input or a scientific rerun. The original holdout
freeze matches its recorded hash and all 13 `labels_read` flags remain false;
no holdout label contents were opened.

The completed 96-run matrix is post-hoc development evidence on already
examined chromosomes. CTCF point gains over raw and sample-specific DWM are
positive throughout, but canonical protection passes 94/96 runs and some
low-depth bootstrap intervals cross zero. No other powered matched TF has
canonical protection. Independent exact-operator negative-control support is
limited to CTCF. Sequence inference uses about 47% of DWM CPU time, but this
excludes extraction, fitting, normalization, scoring and output. End-to-end
resource gates, calibrated functional uncertainty, occupancy-leakage exclusion,
two difficult families and untouched validation remain unestablished.

Proceed in this order:

1. Build a per-TF gate ledger from the existing reports, preserving original
   task denominators, coverage exclusions and evidence splits. Reverify nested
   inputs and the recorded runtime before reuse; do not rerun valid stages.
2. Freeze a fair complete-workflow resource protocol with identical inputs,
   outputs, core budgets and cache conditions. Measure fitting separately from
   reuse, CPU and elapsed time, and absolute peak memory in isolated processes.
   Require numerical/shape equivalence for implementation-only optimizations.
3. Diagnose the low-depth CTCF exceptions and difficult-family failures.
   Existing test results are descriptive only. Fit new candidates on chr1–15
   and select on chr16–18; never retune using examined test or transfer labels.
   Separate correction from scoring geometry in a frozen factorial comparison.
   Prefer simple parametric or hierarchical-spline models before more complex
   Gaussian-process models. Treat insufficient labels as a power limitation.
4. Freeze each candidate's code, data, model, thresholds and evaluation before
   scoring. Audit signed central protection, shoulders, asymmetry, coverage,
   motif/accessibility balance and exact-operator negative-control behavior.
   Require the existing scientific, uncertainty and computational gates in full;
   narrow CTCF gains and faster sequence inference cannot waive them.
5. Keep SK-N-SH/GM23338 labels and final naked-DNA replicate 3 locked until the
   prespecified prerequisites pass. Record failed candidates without promotion.
   Any later production integration requires separate authorization.

The exact current scientific results and limitations are in
`benchmarks/DIRECT_BIAS_STABILITY_96_20260912.md`; component resource results
are in `benchmarks/DIRECT_BIAS_INFERENCE_20260908.md`. No new scientific
experiment was launched during this status-and-goal setup.

### September 22 gate-ledger implementation

The new benchmark-only ledger validates the checksummed summary tables and
preserves all 21 tasks and 48 runs per task. HepG2 HMGA1 has one matched site
per class in the coverage table and no scored rows; it is retained explicitly
instead of disappearing from the denominator. Principal interval availability,
power, canonical protection and unresolved resource/uncertainty/safety gates
are separate fields. No new scientific score or candidate selection is made.

The ledger reproduces 46/48 canonical CTCF runs in HepG2 and 48/48 in K562.
The two failures are HepG2 replicate 2, 10M, seeds 2026 and 2028: both have one
significant central point, with central contrasts -0.054764 and +0.003386.
Both retain positive shoulders. K562 has 5/6 principal sample-DWM comparisons
with both lower bounds positive; HepG2 has 6/6. These observations describe
examined test data and cannot be used to tune a new candidate.

Forty-two focused ledger, stability-extension and report tests passed in
26.28 seconds. The real ledger and its exact immutable rerun passed, with
SHA-256 `33d65e9025ccb23ef09d640bf2e8757107547ba62656297f9c5c360693ed4d88`.
Output: `benchmarks/results/direct_bias_gate_ledger_20260922_v1/gate_ledger.json`.
The full historical nested-input verifier subsequently exited 0 and verified
all 96 runs, frozen inputs, archived reuse and outputs without rerunning a
stage (`--verify-only --code-ref c9ecd82`). The browser audit also passed all
32 pages at three viewport sizes with a dark-system preference. Historical
verification checks recorded code; it does not certify an identical current
scientific runtime for a new experiment.
No scientific computation, holdout access or production change occurred.

### Complete-workflow resource evidence contract

The research promotion checker previously accepted a generic runtime/memory
table without end-to-end scope or matched workload identities. Twelve new
regressions failed before correction. Resource evidence now requires total
process-tree CPU time, absolute process-tree peak memory, fitting plus
application, frozen benchmark/workload identities, matching core/cache settings
and all six repetitions per workload. It rejects nonfinite/negative values,
duplicate or unmatched runs and missing repetition pairs. Worst workload
ratios prevent a fast case from hiding a slow one. Wall time remains recorded
but cannot substitute for CPU time. Legacy component tables fail this gate.
These checks validate table consistency, not the authenticity of declarations;
the benchmark producer must still verify its frozen matrix and artifacts.

The premeasurement protocol is
`benchmarks/manifests/end_to_end_resource_protocol_v1.json`, SHA-256
`eeff0c5ecf2bee1741c8719205a8072b9f736e93d08fb0ef230c9df5b38adc2b`.
It specifies both development cells, all three biological replicates at 10M
and full depth, six timing repetitions and a shared four-core budget. It
includes fitting, input reading, correction, scoring and equivalent output
products. Frozen-model reuse and same-operator correction comparisons remain
separate diagnostics. This locks the protocol, not an executable benchmark:
exact command/input/code/runtime hashes, output-contract parity and validated
process-tree instrumentation are required before measurement. No runtime
ratio or computational promotion pass is claimed. The host's cgroup root is
not writable by the current user; a suitable isolated memory-accounting path
must be verified rather than substituting sampled or parent-only peaks.
After adding the six-repetition guard, all 31 focused resource, promotion and
ledger tests passed in 1.00 second. All 20 documentation-contract tests and the
strict MkDocs build passed. No timing or scientific model run was started.

### Process-tree measurement qualification

The cgroup root is not writable, but the existing user systemd manager can
create a dedicated delegated service. A benchmark-only synthetic probe now
moves its supervisor outside a fresh measured subtree, places the command in
that subtree before exec and reads the kernel's `cpu.stat` and `memory.peak`
after the command and descendants exit. It rejects leftover descendants and
failed commands, retains their measurements and cleans up only its own units.
The measured memory is kernel-accounted cgroup memory, including charged file
cache and kernel memory, not resident-set size. Both benchmark arms must use
this same metric and the frozen cache policy.

The qualified four-case probe used no biological data or model. Two sequential
32-MiB child allocations used 0.498529 CPU seconds and peaked at 42.209280 MB;
two simultaneous allocations used 0.497108 CPU seconds and peaked at 79.167488
MB. An exit-7 command and a parent that left a child running were both rejected.
All probe services were removed. These results distinguish descendant CPU and
simultaneous memory from parent-only or cumulative measurements; they are not
scientific speed or memory improvements. All 23 focused instrumentation and
resource-gate tests passed in 0.47 seconds.

Evidence: `benchmarks/results/process_tree_probe_20260922_v2/probe.json`,
SHA-256 `74fa8f93de827bffc2a948a0da66c7db80c3fd504d771d134c29912e46c2fce0`.
The earlier successful three-case probe remains preserved separately. The
complete-workflow execution freeze still requires scientific command adapters,
equivalent output products, numerical parity, input identities and the exact
runtime. No candidate timing, fitting, holdout reading or promotion occurred.

### Corrected research cut extraction

The BAM-to-score parity work exposed a long alignment being counted twice
when fetched by disjoint windows, and missed soft-clipped cuts outside queried
alignment spans. Three boundary regressions failed before the fix. The reader
now traverses each requested chromosome once and assigns each cut once per
containing site, preserving legitimate overlap and positional row order.
Forty focused regressions passed. No production code was changed.

A frozen, label-free chromosome-1 audit compared old and corrected extraction
on both development cells against independent sorted-cut counting. Exact
oracle parity passed at 2,381 HepG2 and 1,943 K562 sites. Counts changed at 35
and 15 sites, removing 173 and 45 duplicate assignments respectively. This
does not quantify the impact on test metrics or safety. Archived artifacts
remain unchanged, and no model was retuned. The new audit report records input,
code and result hashes and the requirements for a separately frozen corrected
reassessment. Runtime measurement must use the corrected reader.

### Fixed-candidate corrected-count reassessment

`reassess_corrected_strand_counts.py` freezes a new 96-run computation using
the original matrix, exact site/label order, complete matched pairs, model
propensities, dispersions, scoring geometry and bootstrap settings. Both cells'
contexts are frozen before corrected extraction. Every BAM is re-extracted;
no archived observed counts are reused. Existing valid DWM expectations are
reused with verified coordinate order. New outputs live separately from the
historical experiment and retain explicit post-hoc/non-promotion status.

Fifty-one focused reassessment, extraction and stability tests passed in
28.70 seconds. Fixtures cover all 96 runs, verification without recomputation,
corrupted output rejection, preserved incomplete artifacts, changed settings
and matrix rejection, and prohibitions on retuning, holdout reads and old-count
reuse. The new runner requires a configuration freeze before execution and
binds resumed stages to their configuration, run identity and checksums.
It does not establish a new method or close the independent-control gates.
The supervised freeze was launched with durable receipts under
`benchmarks/results/corrected_count_jobs_20260922/freeze_attempt1/` and a new
scientific output root `benchmarks/results/corrected_count_reassessment_20260922_v1/`.
Its job specification SHA-256 is
`3cdaec83fa33c2a613ba5af8e084933caacfcd58a1243c099c853b07842e35e4`.
Check the exit receipt before launching execution; launch evidence alone does
not certify the freeze. At the launch audit, library versions matched the
historical run and only the Python 3.12.3 build timestamp differed; the new
configuration records the current runtime. All 20 documentation-contract
tests and the strict MkDocs build passed after adding the runner.

The before/after reporter `compare_corrected_count_reassessment.py` is separate
from the pinned calculation. It requires the full corrected run matrix and
paired metric/shape/coverage/bootstrap identities, rejects changed class counts
or eligibility, retains unscored tasks via the coverage denominator, and emits
paired tables plus corrected aggregate plots. Six focused tests passed in
1.85 seconds. Reported reader-effect differences are descriptive; they are not
confidence intervals for the implementation change or a new method result.

## Current Baseline

fp-tools is a command-first Python 3.11–3.13 package for bulk and pseudobulk
ATAC-seq footprinting, motif analysis, differential reports, and browser/YAML
wrappers. Scientific workflow logic belongs in
`src/fp_tools/tools/` and shared helpers in `src/fp_tools/utils/`; CLI and GUI
layers remain thin.

The supported console commands are declared in `pyproject.toml`. Public docs,
examples, package metadata, and smoke tests must use that declaration as the
source of truth. Removed TOBIAS-style console aliases must not be restored.

The current project layout, YAML compatibility, output filenames, and command
names are stable public contracts. Internal cleanup should preserve them unless
a separately planned breaking release explicitly changes the contract.

## Completed Capabilities

- Linux CLI/container raw-read ATAC preparation with modern and legacy
  processing profiles.
- A managed, versioned external-tool runtime downloaded on first use, with
  Linux raw-read components, cross-platform optional MEME components, and the
  complete Linux container retained as an explicit backend.
- Cross-platform one-command bulk analysis from coordinate-sorted BAM/BAI and
  matching peak BED files through portable multi-comparison HTML reports.
- Separate Linux CLI/container FASTQ preparation through `prepare-atac`, which
  writes the BAM/BAI and peak BED sample table accepted by the downstream bulk
  workflow.
- Tn5 bias correction, background scaling, footprint scoring, candidate calls,
  known-motif matching, de novo motif preparation, and aggregate plotting.
- Replicate-aware differential footprint reports, including per-sample motif
  matrices and empirical-Bayes residual-variance moderation with biological
  samples as the inferential units. Bundle and self-contained multi-comparison
  reviews share one browser implementation; standalone reports preserve input
  order with an exact-record selector and remove aggregate-only controls when
  profiles are unavailable. The shared browser also supports aggregate-neutral
  volcano highlighting, lightweight user-selected TF labels, and a two-mode
  ranked-motif waterfall with reciprocal score/significance labels, directional
  blue/red color strength, a stable square volcano viewport, and comparison
  labels embedded in exported waterfall, volcano, and combined-panel SVGs.
- Region-set differential footprinting for one sample or paired biological
  replicates, with equal region weighting, optional matching strata,
  resampling-based uncertainty, paired empirical-Bayes inference, and the same
  portable interactive report format.
- A compact HepG2 HNF4A/FOXA2 region-set example with four mutually exclusive,
  baseline-matched region classes, six pairwise comparisons, three biological
  replicates, complete motif statistics, and curated all-site aggregate views.
- Pseudobulk fragment/BAM generation and single-cell signature reporting.
- Streamlit GUI whose saved YAML runs through `run-yaml-workflow`.
- Local tests, package builds, GitHub CI, GitHub Pages deployment, and manual
  PyPI publication.
- A complete amd64/arm64 container with the external genomics toolchain and a
  browser GUI as its default entry point.
- BAM-first self-contained GUI executables for Windows x64 and Apple Silicon,
  with branded operating-system icons, a bundled native application window,
  and frozen-safe command and workflow dispatch. Every GUI starts from BAM/BAI
  plus peak BED files. Linux and Intel macOS use the browser-based Python
  package; the complete container remains available on amd64 and arm64.
- The Apple Silicon download remains an explicitly unsigned preview. Its DMG
  includes the same one-line quarantine-removal instruction shown on the
  installation page, and CI verifies the ad-hoc bundle signature, ARM64
  executable, DMG contents, quarantine removal, and native-window launch.
- Branded GitHub and MkDocs presentation with responsive embedded report and
  GUI demos plus automated desktop/mobile browser audits.
- Responsive GUI validation output with readable wrapped paths and locally
  scrolling YAML previews, verified in the native Windows and macOS audits.
  Validation uses concise field names, preserves whole words where possible,
  and failed desktop audits retain screenshots as CI artifacts for diagnosis.
- Responsive GUI navigation and forms across every destination: the mobile
  sidebar has named open/close controls, the persistent sidebar header stays
  compact, workspace details remain available on demand, loaded YAML refreshes
  all form controls, and main-panel labels use explicit high-contrast colors.
  Frozen Windows and Apple Silicon audits exercise every page at desktop and
  narrow widths and retain screenshots for visual review.
- The documentation, browser GUI, static demos, and native Windows/macOS
  windows use one explicit light theme regardless of the operating-system
  preference. Browser audits emulate a dark system preference to guard this
  presentation contract.
- A task-oriented GUI home page, neutral first-open forms, centralized field
  labels, compact advanced options, clearly disabled launch controls, and
  responsive form/run panels. Every route now reuses the Home page's light
  hero-card, typography, borders, controls, and action colors; form controls
  remain visibly bounded on white cards, and the form/run layout is compact at
  desktop widths before stacking cleanly at narrower breakpoints. Fresh desktop
  installations no longer open with repository-relative example paths
  presented as user inputs.
- Browser-verified plot controls and SVG exports for aggregate and
  aggregate-free reports, including explicit subplot-bound checks.
- Command-aware logging across shared analysis engines, so `match-motifs` and
  `diff-footprints` retain their public names in direct and wrapper-run logs.
- A manifest-pinned, storage-conscious ENCODE workflow for 17 biological
  replicates from seven cancer cell lines. Its resumable runner uses a
  pair-specific union of released IDR-thresholded peaks, peak-q95 scaling, and
  the same differential workflow as the preserved three-replicate
  K562-HepG2 report. The dependency-free static browser has two directional
  selectors and reads compact canonical report payloads.
- Benchmark-only footprint detectability auditing that ranks motifs within
  source analyses, collapses repeated pairwise/control rows to independent
  biological contexts, joins exact-name RNA evidence, records hashed inputs,
  and emits both machine-readable tables and a searchable HTML report. A
  separate matched-label classifier applies prespecified correction, scoring,
  coverage, ambiguity, and information-limit diagnoses without treating low
  aggregate scores as evidence of TF absence.
- A locked footprint-improvement study with 35 ENCODE cell/TF tasks,
  development/holdout and chromosome partitions, summit-supported motif-site
  labels, matched negative controls, deterministic nested depth subsets,
  raw/PWM/DWM/reused-bias correction arms, genomic-block uncertainty,
  depth/correction summaries, one-time holdout promotion gates, and tiered
  nutrient-stress replication. The external nutrient resources are pinned to
  GSE144833 and GSE137034 (including GSE137031/GSE137032).
- A label-free percentile-rank evidence-fusion primitive and locked evaluator.
  Its first footprint/PWM soft-OR candidate passed the available K562/HepG2
  development slice but failed the one-time cell-line holdout because it did
  not meet mean-gain or strong-positive non-regression gates. It remains
  experimental; the production footprint score is unchanged.
- An opt-in aggregate-shape detectability mode with per-site outer-flank RMS
  normalization, site-level 95% confidence bands, quantitative depletion
  diagnostics, explicit underpowered/not-detected states, and functional
  per-panel or grouped y-axis scaling. This mode improves visualization and
  failure identification without changing the production footprint score or
  presenting aggregate shape as proof of occupancy.
- An opt-in dual-geometry `call-footprints --score hybrid` arm that adds a
  low-weight, locally standardized wide symmetric depletion channel. It
  reproducibly improved CTCF and REST on locked K562/HepG2 chromosomes but
  regressed some JUND/MAX tasks, so it remains experimental and cannot replace
  the production `footprint` default.
- A staged per-TF footprint research engine that searches correction arm,
  center/shoulder geometry, local normalization, center statistic, and
  asymmetry penalty on train chromosomes, freezes candidates on validation
  chromosomes, and evaluates identical finite sites only after test unlock.
  Reproducible ENCODE peak discovery, cell-specific motif scans, summit-based
  labels, optimal motif/accessibility matching, matched legacy comparisons,
  and bound/unbound aggregate figures support the experiment. The first
  10-million-fragment test produced large CTCF gains and selective MEF2A,
  MEF2D, ARID3A, and MYC gains, but also clear TF/context regressions and
  underpowered tasks. It remains research-only and does not change the default
  scorer.
- A second locked research specification for a fast ChromBPNet-inspired
  alternative. The implementation now includes conditional parametric
  sequence-bias models, pooled sample adaptation, calibrated Poisson/NB
  residuals, weighted functional PCA, spline and 25-knot Matérn
  Gaussian-process-equivalent footprint mixtures, an exact-GP reference,
  hybrid FDA-GP scoring, and replicate-level functional differential tests.
  Model files use checksummed NPZ plus JSON rather than executable pickle.
  These APIs remain research-only until the GM12878/IMR-90, naked-DNA,
  uncertainty, runtime, and non-regression gates pass.
- A label-free parametric-bias benchmark now extracts GC-matched low-signal
  nonpeak or mitochondrial control windows, tests +4/-5 and +4/-4, fits
  SELMA-equivalent adjacent-simplex and 81-bp interaction models, and compares
  sample, pooled, and pooled-plus-adapted estimates. Validation likelihood must
  improve over the uniform conditional control before any model advances. The
  functional benchmark can consume a frozen NPZ model directly and derive
  total-preserving expected motif-site profiles without producing genome-wide
  intermediate tracks.
- The functional research path now has a strand-aware profile artifact. It
  preserves plus/minus observed and expected cuts, swaps strands when reverse
  motifs are orientation-aligned, and exposes combined, shared-strand, and
  antisymmetric signed-deviance channels. This enables genuine TF-specific
  asymmetry experiments rather than inferring asymmetry from a combined track.
- Functional evaluation now produces bootstrap-banded aggregate profiles,
  explicit depth/width/shoulder/asymmetry/periodicity descriptors, blinded
  multi-correction PDF panels, and unsupervised phenotype clusters learned
  from normalized residual curves rather than motif-family names. Numeric
  discrimination is therefore reviewed alongside visible footprint shape.
- Replicate-aware differential functional orchestration now consumes frozen
  strand-profile artifacts and reports the complete condition-difference
  curve, pointwise and simultaneous uncertainty, global permutation testing,
  BH-adjusted significance, replicate consistency, and changes in footprint
  depth, width, shoulders, asymmetry, and periodicity. Nutrient-stress data
  remain locked until model selection is complete.
- Bias-model leakage diagnostics now score every eligible development JASPAR
  motif against nearby local-sequence controls without reading ChIP labels.
  They report bootstrap-banded response functions, broad center/flank effects,
  and cross-cell/model concordance. A sequence response is a review flag—not
  automatic evidence of occupancy leakage—until naked-DNA and held-out TF
  behavior distinguish intrinsic Tn5 preference from overlearned chromatin.
- A fail-closed functional promotion auditor now implements the complete
  preregistered decision: locked external metrics, functional separation,
  positive-control non-regression, clustered bootstrap evidence, naked-DNA
  false positives, motif-response review, stability, CPU/memory/model size,
  uncertainty coverage, and GP-versus-spline justification. Development-only
  evidence or any missing table cannot produce a promotion pass.
- The frozen parametric factorization experiment now has checksum-locked
  SK-N-SH and GM23338 holdout manifests, LOG21/LOG41/LOG81 and SELMA10
  coefficient grids, bounded sample bias-strength calibration, flank-only
  accessibility splines, hierarchical latent footprint profiles, five
  residual ablations, and immutable resumable stage manifests. Independent
  enzyme controls are compared with paired block bootstraps against a safe
  TOBIAS-style all-pairs DWM reference. Candidate-screening and final
  naked-DNA libraries are extracted without fitting, and motif-centered
  observed-minus-predicted residuals have an explicit leakage gate. All new
  artifacts use checksummed NPZ plus JSON and remain research-only. Continuous
  residual safety is evaluated with cutoffs frozen from chr16--18 matched
  ChIP-negative sites before they are applied to label-free naked-DNA
  replicate 2; both finite-site and nonzero-cut false-positive rates must pass
  before a residual can enter the immutable configuration freeze.
- The first frozen conditional-multinomial checkpoint conditions explicitly on
  each motif site's cut total and learns only profile shape. Across the nine
  internally eligible K562/HepG2 tasks it improved mean AUROC by 0.0139,
  relative AUPRC by 4.38%, and functional separation by 28.7% over DWM, while
  all 14 evaluated task thresholds passed naked-DNA replicate-2 safety. The
  gain is not general: the same candidate was 0.0184 AUROC and 0.73% relative
  AUPRC below the exact raw-signal guardrail on average. CTCF is the clear
  internal exception, improving AUROC/AUPRC over raw by 0.0647/15.95% in HepG2
  and 0.0421/9.43% in K562 with visibly stronger canonical aggregate
  protection. FOXA1's smaller numerical gain does not yet have a convincing
  canonical aggregate shape. A count/FDA rank-consensus screen also failed the
  10% AUPRC gate, and FDA-heavy mixtures produced naked-DNA false positives.
  These results remain internal research evidence: full-depth confirmation,
  control-source qualification, the official ChromBPNet comparison, and
  locked external validation are still required before any promotion.
- Independent control-source qualification now favors mitochondrial-only
  training over combined mitochondrial plus naked-DNA training for both cut
  conventions. The retained internal +4/-5 reference remains the five-seed
  mitochondrial LOG21 ensemble; adding the naked-DNA library did not improve
  its paired held-out control likelihood. A validation-frozen partial-bias
  experiment identifies one reproducible TF-specific exception to the broader
  failure: CTCF shrinkage improves both AUROC and AUPRC over raw signal in
  HepG2 and K562 at 10M, 25M, and 50M fragments for every seed, and across all
  three biological replicates in each cell. Mean full-replicate AUROC gains
  over raw are 0.0410 in HepG2 and 0.0378 in K562. The same policy makes zero
  calls in the independent naked-DNA replicate-2 CTCF panels, with a 1.88%
  Wilson upper bound. For other eligible TFs, large gains over DWM generally
  return to raw-signal performance rather than add information; depth reports
  therefore distinguish genuine gains over raw from recovery after DWM
  overcorrection. This is a CTCF-specific internal research route, not evidence
  for a general correction replacement, and the package default is unchanged.
- A separate integrity-checked detector ledger now prevents TF-specific shape
  gains from being conflated with bias-correction gains. The globally frozen
  `count_only` operator passes independent naked-DNA replicate-2 safety for all
  14 evaluable task panels. FOXA1 has a reproducible occupancy-classification
  gain: on the exact common chromosome-test support it improves AUROC by
  0.0624 and relative AUPRC by 5.92% over raw signal, both metrics improve in
  all three HepG2 biological replicates, and 25M/50M seeds are stable. Visual
  and numeric shape auditing rejects this as a footprint improvement, however:
  the candidate has central enrichment rather than protection (depletion
  -0.0423). It is therefore classified as an occupancy-signal diagnostic, not
  a FOX footprint model. The detector also improves CTCF over DWM, but its
  AUROC comparison with raw signal changes with support/depth and therefore
  does not supersede the simpler CTCF shrinkage result. The other seven
  internally eligible tasks fail the raw guardrail, and 12 tasks remain
  underpowered. No general detector or correction promotion follows.
- A five-family secondary screen now compares validation-frozen count, FDA,
  anchored FDA, residualized FDA, and hybrid FDA-GP policies against raw and
  DWM scores on exact common test support. It joins independent naked-DNA
  safety and requires central protection, a negative ChIP-positive-versus-
  matched-negative center contrast, and pointwise aggregate evidence before a
  discrimination gain is called a footprint gain. Only HepG2 CTCF passes the
  strong shape guard: count improves AUROC by 0.0725 and relative AUPRC by
  15.90% over raw (bootstrap lower bounds 0.0406 and 8.96%), while anchored
  FDA improves AUROC by 0.0545 and relative AUPRC by 3.76% (lower bounds
  0.0408 and 0.80%). K562 CTCF remains support-sensitive, and no difficult TF
  family produces a safety-passing significant shape-guarded gain over raw.
  This confirms that the current FDA/GP variants do not yet solve hard-TF
  footprint detection and prevents accessibility-driven occupancy gains from
  being reported as footprints.
- A canonical-protection conditional model now constrains the learned center
  to depletion and the shoulders to enrichment. Validation-frozen models were
  evaluated for all 14 K562/HepG2 task contexts and passed independent
  naked-DNA replicate-2 safety. On exact chromosome-test support, CTCF gains
  over raw were 0.0589 AUROC/14.01% relative AUPRC in HepG2 and 0.0508/8.94%
  in K562 with visibly canonical protection. The model beats DWM on both
  metrics in all six biological replicates, but its full-replicate AUROC is
  below raw in all six; no other TF has a defensible shape-guarded gain. It is
  therefore a useful detector ablation, not a promotion candidate.
- Both complete 120-fit mitochondrial and naked-DNA sequence-bias grids are
  now collated into five-seed geometric ensembles. The naked-DNA grid ranks
  LOG21, +4/-5, L2 0.01 first on held-out conditional likelihood; cross-source
  qualification still retains the mitochondrial LOG21 ensemble for the
  internal footprint experiment. Resumption is fail-closed: NPZ/JSON
  checksums and the full source, shift, feature, penalty, depth, sample, seed,
  window, and cut identity must match before an existing fit is rescored.
- The official ChromBPNet reference is pinned to commit
  `09938fdb4397ec0006510e5251e48920a505d4de` and a digest-locked TensorFlow
  2.8 GPU container. On identical K562 test sites and conventional geometry,
  frozen CPU parametric bias scores CTCF at AUROC 0.7421/AUPRC 0.7180, close
  to ChromBPNet bias at 0.7509/0.7229 and above DWM at 0.6810/0.6455. Its
  paired gains over DWM have positive chromosome-block intervals; its AUPRC
  difference from ChromBPNet includes zero. Neither bias method is general:
  both regress eligible MYC, while the parametric model gives only a small
  ZNF384 gain over DWM and remains below raw. The checksum-locked HepG2 bias
  model, regulatory training, predictions, and same-site comparison are now
  checksum-verified complete. The HepG2 comparison used dispersion zero,
  whereas the recorded K562 comparison used 1.6144318781609228; their results
  must remain separate experiments rather than a pooled common-operator claim.
- A completed-reference audit now verifies comparison inputs, outputs and
  canonical identities and reports all 21 development TF tasks with signed
  protection, shoulders, asymmetry, label-dependent site retention, motif-score
  balance and accessibility balance. Only CTCF has significant parametric gains
  over both DWM and raw with canonical protection in both recorded comparisons.
  HepG2 CTCF scores 0.7149 AUROC/0.7388 AUPRC versus DWM 0.6477/0.6250 and
  ChromBPNet bias 0.7253/0.7570. The parametric mean gains across eligible tasks
  are only 0.0112 AUROC/1.99% relative AUPRC in HepG2 and 0.0053/1.25% in
  K562, below the general promotion gates. CTCF motif scores remain imbalanced
  (standardized differences 0.482 HepG2 and 0.463 K562). The older deep-reference
  report cited factorization-residual naked-DNA safety, which does not certify
  the direct parametric scoring operator. A new direct-operator audit freezes
  code, inputs and validation-only thresholds before scoring existing replicate-2
  profiles. Its coverage qualification reuses those scores and thresholds:
  only CTCF has sufficient archived DWM control coverage. On DWM-comparable
  informative sites, HepG2 has 0/121 false calls (Wilson upper 3.08%); K562
  has 1/114 (upper 4.80%, increase over DWM 0.88 percentage points). Each
  excludes one informative site lacking a DWM expectation; the absolute
  candidate FPR also passes with those sites retained. This is retrospective
  development candidate screening, not final untouched safety. Sparse counts
  and the absence of a large residual flag do not establish no occupancy
  leakage. Other TFs cannot pass relative safety from this CTCF-region DWM
  archive. See
  `benchmarks/CHROMBPNET_AUDIT_20260908.md` for the evidence and remaining work.
- The one-page direct-bias renderer now rejects legacy factorization safety,
  requires exact comparison identity plus coverage-qualified direct safety,
  and gates on significant paired gains over DWM/raw and canonical protection.
  It attaches safety only to the candidate, reports the actual cell/dispersion
  and covariate imbalance, and preserves historical reports. The new HepG2
  CTCF report is development evidence with candidate-screen safety; no general
  correction, final untouched safety or package promotion is claimed. Its six
  files were uploaded to the new Box folder
  `Yaoxiang/fp-tools/2026-09-08_HepG2_CTCF_direct_bias_candidate_screen`
  on 2026-09-08; `rclone check --download --one-way` reported six matches and
  zero differences. Historical report folders were preserved.
- A frozen post-hoc covariate sensitivity now covers all 21 development tasks
  with three declared matching policies and unchanged scores/dispersions. The
  primary motif/accessibility match retains 634 HepG2 and 609 K562 CTCF pairs
  with covariate SMDs below 0.01. Parametric gains over DWM remain significant
  in both cells (AUROC +0.0477/+0.0513, relative AUPRC +12.70%/+10.73%) with
  canonical protection, but K562's intervals against raw include zero. A third
  fixed policy also balances direct observed cut total and keeps that pattern.
  Fourteen tasks become underpowered; balanced FOXA1/MYC remain near chance
  and powered ZNF384 tasks remain weak without canonical protection. General
  mean-gain and difficult-family gates still fail. This is sensitivity on
  examined chromosomes, not new holdout evidence or a new candidate. See
  `benchmarks/CHROMBPNET_COVARIATE_DIAGNOSTIC_20260908.md` for complete results.
  Historical verification now checks nested cohort artifacts and can validate
  recorded calculation-code blobs against their immutable Git commit without
  rerunning completed experiments.
- Direct-bias stability now has a fixed research plan for the 36 available
  development BAMs: five seeds at 10M/25M/50M for replicate 1 plus all three
  full-depth biological replicates per cell. Existing partial-bias matrices
  use tuning chromosomes and cannot supply the missing test-site evidence.
  Preflight also found zero DWM expectations at 65 HepG2 and 66 K562 CTCF
  sites in the historical finite-support comparison. The new diagnostic
  freezes positive expected coverage across all runs, preserves complete
  existing matched pairs, and keeps zero observed-count sites. It separates
  the original full-depth DWM expectation from sample/depth-matched DWM.
  All 36 runs are complete and verified from calculation commit `c453a37`.
  CTCF retains canonical protection and positive point gains over raw and both
  DWM comparators throughout. All six full biological replicates have positive
  AUROC/AUPRC bootstrap intervals against raw; matched AUROC gains are
  0.0592–0.0627 in HepG2 and 0.0575–0.0666 in K562. No other powered TF passes
  the canonical-shape screen. Current cut-total SMD drifts above 0.1 in 11/36
  matched CTCF runs (maximum 0.134), so fixed matching does not eliminate
  accessibility confounding. General relative-AUPRC and difficult-family gates
  still fail; full runtime/memory, uncertainty and untouched safety remain
  incomplete. See `benchmarks/DIRECT_BIAS_STABILITY_20260908.md` for the full
  evidence, expected-coverage exclusions and limitations.
- A concise CTCF stability report now requires all declared run directions,
  full-replicate paired intervals, canonical protection and matching direct
  safety identities before rendering. Its 36-row metrics TSV retains balance
  drift and clearly distinguishes candidate-screen safety from untouched final
  safety. The six report files were placed in the new Box folder
  `Yaoxiang/fp-tools/2026-09-08_CTCF_direct_bias_stability_coverage`; download
  verification reported six matches and zero differences on 2026-09-08.
- A separate frozen sequence-inference benchmark compares the archived
  production Cython DWM scorer with the unchanged parametric model on 512
  coordinate-selected loci per development cell. Six paired repetitions use
  alternating method order and fresh single-CPU processes. It validates the
  candidate against archived strand propensities and measures absolute Linux
  worker peak memory. Its component ratios cannot certify the full workflow
  runtime/memory gates; fitting, BAM/FASTA I/O, local correction, detector
  scoring and track output remain outside this measurement.
  All 24 workers completed: candidate/DWM mean CPU ratios are 0.47385 in
  HepG2 and 0.47654 in K562, with maximum worker-memory ratios 0.99247 and
  0.99202. Candidate propensities reproduce the archived values; every method
  is numerically stable across timing repetitions. See
  `benchmarks/DIRECT_BIAS_INFERENCE_20260908.md` for the frozen protocol,
  checksums, validation and explicit disposition of the remaining gates.
- Nested depth-matrix resume now verifies source identity, sampling metadata,
  BAM and index checksums. Previously it silently skipped any existing BAM,
  including interrupted outputs. Regression fixtures reproduced six failures
  before the fix; the corrected implementation also rejects a source that
  differs from an explicitly frozen checksum. Historical data are preserved;
  this stricter path prepares missing inputs without rewriting old manifests.
- The missing depth/seed inputs for HepG2/K562 replicates 2 and 3 now have a
  separate frozen preparation plan. It declares all 60 missing combinations,
  the original pair-preserving hash selection and the same refitted DWM
  command used for replicate 1. A native fragment counter is verified against
  the original predicate for all 4,096 SAM flag combinations, including MAPQ
  zero and unpaired reads. Preparation checks source and nested output hashes,
  pins the production source tree/runtime, limits concurrency to two jobs and
  uses a kernel lock. This stage prepares data only; subsequent evaluation
  must freeze support across all 96 runs and preserve complete matched pairs.
  The four fragment-count stages are complete: HepG2 replicates 2/3 contain
  96,480,203/105,920,129 usable fragments, and K562 replicates 2/3 contain
  76,393,585/114,471,886. All support the requested 50M depth. All four sampling
  stages completed on 2026-09-09 UTC: independent checks passed for the exact
  60-target matrix and all 180 BAM/index/metadata files, including checksums,
  parent identities, metadata and indexed alignment counts. Conventional DWM
  correction and the 96-run scientific evaluation completed on September 10.
  Independent verification passed for all 60 corrections and all 96 runs on
  September 12; no completed scientific stage was rerun.
  See `benchmarks/DIRECT_BIAS_REPLICATE_DEPTHS_20260908.md` for the frozen
  identities, actual selected-depth ranges and remaining requirements.
- The first sampling process terminated with exit 143 before completing any
  sampling stage; the dependent queue correctly failed without starting an
  evaluation. The termination cause is unknown. All four source BAM/BAI pairs
  and count stages reverified, and 30 interrupted BAMs (54,098,858,454 bytes)
  were preserved with checksums and unchanged file identities. No partial BAM
  was accepted as complete. A separate logged-job supervisor now pins the
  command specification, records supervisor/child identities and terminal
  status, and forwards stop signals to its own child group. Detached launch
  survives the short launcher process; receipt failure stops an unrecorded
  child. Resume retains the original preparation/evaluator code, source tree,
  runtime, inputs, scientific settings and completed counts.
- The complete 96-run evaluator now extends the archived 36-run calculation
  with the 60 prepared replicate/depth inputs. It verifies nested sampling and
  correction provenance, pins the existing scientific code and runtime, and
  freezes both cells' all-run expectation coverage before new extraction or
  scoring. All 36 archived profiles are reused after coordinate/label and
  DWM-track checks. Metrics are reused only on identical original and matched
  cohorts; changed coverage requires a fresh comparison on the new support.
  The plan pins the evaluator checksum, and queued execution also pins the
  plan checksum. A small 96-run fixture exercises unchanged and changed
  coverage, complete pairs, zero-cut sites, reuse, verification and nested
  output corruption. The real calculation is complete: 36 archived
  extractions were reused and all 96 evaluations used the newly frozen common
  support. Matched CTCF AUROC/AUPRC point gains over raw and sample DWM are
  positive in all runs, but protection passes only 94/96 and some low-depth
  bootstrap intervals cross zero. No other powered matched TF has canonical
  protection. All 21 aggregate pages were reviewed. The result does not pass
  promotion gates; holdouts and final naked-DNA controls remain unopened.
  See `benchmarks/DIRECT_BIAS_STABILITY_96_20260912.md` for verified hashes,
  per-TF results, balance limitations and the historical-runtime audit.
- A separate current-depth balance diagnostic addresses the observed cut-total
  drift in the 36-run archive. Its fixed secondary rule retains original
  partners only within the existing 0.2-SD caliper on current log1p cut totals,
  using the original chr1–15 per-cell/TF scales. All selections precede new
  scores; underpowered and zero-pair outcomes remain explicit. Primary results
  and the queued 96-run evaluation are unchanged. This is a post-hoc
  sensitivity check, not a new model or causal validation. The completed
  36-run calculation balances all measured CTCF covariates below 0.04 SMD.
  All six full-depth CTCF runs retain significant gains over both DWM
  comparators, but HepG2 replicate 1's AUROC interval against raw includes
  zero, and several lower-depth raw/protection checks fail. FOXA1/MYC lose
  power after filtering; excluding them must not change the task denominator
  used to claim a promotion gate. CTCF remains the only powered canonical
  task. See `benchmarks/DIRECT_BIAS_DEPTH_BALANCE_20260908.md` for the complete
  result, 21-page figure audit, frozen identities and limitations.
- The significant, safety-passing internal CTCF partial-bias result has a
  concise one-page report, metrics, checksums, and limitations in Box under
  `Yaoxiang/fp-tools/2026-09-02_CTCF_frozen_parametric_bias_shrinkage`.
  This is explicitly a research artifact; it does not change the package
  default and is not evidence for a general TF correction.
- The K562 deep-reference comparison has a second checksum-verified Box report
  under
  `Yaoxiang/fp-tools/2026-09-02_K562_CTCF_parametric_vs_ChromBPNet`.
  It contains a one-page same-site aggregate and metric comparison, its source
  metrics, and explicit limitations; the report does not claim general
  promotion from the CTCF result.

## Current ENCODE Resource

- All 21 pairwise comparisons are available in the static browser, and every
  aggregate profile uses the complete motif-site set rather than bound-only
  sites.
- The preserved K562-HepG2 comparison retains its reference scientific
  payload, and the browser loads compact report data, profile shards, and motif
  logos on demand.
- The public report page embeds this browser with K562 versus HepG2 as the
  default comparison.

## Near-Term Priorities

1. Keep the seven-line ENCODE cancer resource reproducible and
   storage-conscious. Preserve all 1,019 motifs, 17 biological replicates, and
   21 prespecified contrasts; resource membership remains independent of
   observed differential results, and the K562-HepG2 comparison must retain its
   exact preserved scientific payload.
2. Keep README, MkDocs, examples, CLI help, package metadata, and the static
   cancer-cell-line browser synchronized.
3. Add focused regression tests for every user-visible bug or command-contract
   change.
4. Reduce avoidable warnings and multiprocessing fragility without changing
   scientific output.
5. Improve GUI cancellation and output previews only through existing command
   and YAML interfaces.
6. Keep benchmark claims proportional to validated public data and recorded
   metrics.
7. Keep public pages free of broken assets, browser errors, responsive
   overflow, inaccessible navigation, and broken live embeds.
8. Keep region-set examples outcome-independent: define groups from external
   annotations, match baseline signal before testing, report all motifs, and
   use explicit display motifs only to configure the initial browser view.
9. Maintain BAM/peak scientific I/O parity on Windows, macOS, and Linux, plus
   desktop-app parity on Windows x64 and Apple Silicon. Test Linux raw-read and
   cross-platform MEME runtime artifacts separately, and keep the complete
   multi-architecture Linux container as a reproducible alternative backend.
10. Continue the locked footprint detectability matrix beyond the completed
    K562/HepG2 10-million-fragment, seed-2026 raw/PWM/DWM correction slice.
    Add the remaining depths, seeds, replicates, and naked-DNA controls before
    changing the default correction. The first label-free
    footprint/PWM fusion candidate has been rejected after its one-time
    MCF-7/A549/HCT116/Panc1 evaluation and must not be retuned on that holdout.
    Keep
    nutrient data as a prospective application rather than model-training data,
    and do not promote a new scorer until held-out performance, calibration,
    naked-DNA false-discovery control, and strong-positive non-regression gates
    pass.
11. Run the `footprint_functional_v1` stages in order: establish the cut-site
    convention and control-trained parametric bias models; compare residual
    formulations; freeze spline/FDA/GP/hybrid models on K562/HepG2; use the
    already examined GM12878/IMR-90 only as transfer diagnostics. Keep the
    current DWM correction and footprint score as defaults unless every
    scientific and CPU performance gate passes.
12. Complete `frozen_parametric_factorization_v1` without opening its new
    holdout labels out of order: finish the control-source grid, retain at most
    two frozen ensembles, select one residual on chr16--18, open chr19--22/X
    only after the configuration hash is written, then run depth/replicate,
    independent naked-DNA, diagnostic transfer, ChromBPNet reference, and the
    new SK-N-SH/GM23338 holdout stages. Nutrient-stress analysis and Box
    before/after reports follow only after a candidate passes its safety and
    significance gates. Do not merge this research branch into `main`.

## Deferred or Experimental Work

- Broad supervised TFBS prediction and calibration.
- Motif-removal, motif-relaxed, or motif-free recovery claims.
- Variant-scoring case studies beyond the existing optional utility.
- Footprint competition or nucleosome-decomposition claims.
- Cross-study RNA/ATAC causal interpretation where study and cell state are
  confounded.
- Detectability-aware occupancy probabilities, abstention statuses, and
  TF-family hierarchical models. The benchmark schemas are available, but
  these remain experimental until the matched-label and perturbation gates in
  the footprint detectability study are satisfied.

These areas may have scaffolding under `benchmarks/` or internal tool modules,
but they are not promoted as primary package workflows until their datasets,
metrics, validation, tests, and documentation are complete.

## Maintenance Gates

Before pushing broad changes, run the full unittest suite, all console-script
help checks, YAML dry runs, `pip check`, strict MkDocs, release artifact checks,
and relevant focused regressions. Release, GitHub Actions, Pages, and PyPI
instructions belong only in `RELEASE_CHECKLIST.md`.

The GUI remains a thin wrapper included in the standard installation. Future
GUI changes must keep direct CLI use primary, retain reusable YAML, avoid
hosted-service assumptions, enforce the BAM/BAI plus peak BED starting point,
and use the current static demo and layout as the visual baseline.
