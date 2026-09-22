# Balanced-cohort geometry diagnostic and numerical limitation

The simple geometry comparison is complete, but numerical rank instability
prevents accepting it as a qualified candidate evaluation. Preserve this run
and address near-tied scores in a separately frozen successor. No method is
selected or promoted from these results.

## Design frozen before evaluation

Used the existing 4,053 matched training pairs, all 21 transcription-factor
tasks, seven correction arms, and the existing eleven geometries. Each of
chromosomes 1–15 was held out in turn. Fitting excluded that chromosome and
both members of every matched pair touching it. Minimum fitting support was
200 sites per class. Among geometries with a negative fitting center below
both shoulders, selected the greatest fitting area under the receiver
operating characteristic curve (AUROC); exact ties retained the first grid
entry. Insufficient support or absence of protection retained the conventional
geometry with an explicit fallback reason.

Compared three detectors within each correction arm: fixed geometry on
unnormalized profiles, fixed geometry after per-site normalization, and
selected geometry after the same normalization. Normalization used fixed
outer detrending, root mean square scaling, fallback 1 and clipping at 12.
Zero-read 201-position windows had score zero. Evaluation was within each
held chromosome, with equal-chromosome mean metrics; differently fitted
scores were not pooled. Cohorts were originally matched using the full
training pool, so this is not fully nested validation. No tuning, examined
test, transfer, or untouched promotion labels entered this stage.

## Completion and preliminary observations

All 2,205 fitting records and 6,615 evaluation rows were retained. There were
543 supported protective selections, 1,225 underpowered fallbacks and 437
no-protection fallbacks. A different geometry was chosen in 528 fits.
Of 441 task/arm/detector summaries, 357 have 15 evaluable chromosomes, 42
have 14, 21 have 13, and 21 have none (HepG2 HMGA1).

For CCCTC-binding factor (CTCF), the dinucleotide weight matrix (DWM)
residual arm selected center/shoulder geometry (20,25,60) in all HepG2
folds. K562 selected (20,25,60) in eight folds and (15,20,60) in five,
with two underpowered fallbacks. K562 MEF2D selected (3,8,60) in fourteen
DWM folds, with one underpowered fallback. These reinforce geometry as a
development question, not a proven general improvement.

For radius-500 residuals, selected geometry did not improve the original
mean within-chromosome AUROC relative to normalized fixed geometry in
HepG2 ZNF362 (0.576334 versus 0.581932), HepG2 ZNF384 (0.577601 versus
0.584812), or K562 ZNF384 (0.530952 versus 0.541540). These are provisional
arithmetic-specific diagnostic values. Changing a detector does not change
the underlying aggregate profile, and a changed geometry's protection
criterion is not the original fixed-center shape gate.

## Numerical failure found by independent verification

An equivalent least-squares normalization implementation reproduced scores
within 2.90878432451791e-13, but changed AUROC or area under the
precision-recall curve (AUPRC) in 287 evaluable rows: 148 normalized-fixed
and 139 normalized-selected rows. No unnormalized-fixed rows changed.
The largest chromosome-level absolute changes were AUROC 0.175 and AUPRC
0.067619. The largest absolute equal-chromosome mean changes were AUROC
0.0139413 and AUPRC 0.0067127. The extreme chromosome case was K562 TBP,
raw-count correction arm, chromosome 14, normalized-fixed detector; its
maximum score perturbation was only 2.775557561562892e-17.

Tiny floating-point distinctions are being treated as real ranks. Archived
metrics are reproducible from archived scores, but this is not numerical
robustness. The sensitivity recomputes held scores for the archived geometries;
it does not establish robustness of fitting selection under alternative
arithmetic. A successor must apply a prespecified numerical tie policy both
during selection and evaluation, with synthetic sparse/flat-profile tests,
batching checks and equivalent-solver checks. Do not tune that policy to
increase biological performance. Preserve this failed diagnostic unchanged.

## Validation and provenance

- Seven focused tests passed in 2.96 seconds, covering chromosome/pair
  exclusions, held-label isolation, protection/fallback rules, weighted
  bootstrap ranks and source evaluation behavior.
- Durable execution exited 0 in 44.864883 seconds. This is diagnostic
  elapsed time, not end-to-end central processing unit or memory evidence.
- Separate checksum, code and runtime verification passed.
- Independent audit checked all 2,205 fitting choices against rank-sum
  metrics and the prescribed fitting arithmetic; all saved scores and
  coordinates; all 6,615 evaluation rows (6,216 evaluable); all available
  signed contrasts; and all 441 summary rows.
- Independent least-squares normalization exposed the numerical limitation
  above. It is recorded as a failure, not hidden by artifact verification.

Artifacts are ignored under `benchmarks/results/`:

```text
balanced_geometry_crossfit_20260922_v1/configuration.freeze.json
SHA-256 63ed40e61ce84bd7f83cac3f9ed94ec1e6d7884fc8af7574d9abe36fd9b48949
balanced_geometry_crossfit_20260922_v1/crossfit_manifest.json
SHA-256 38f3e03a55dcb273e47bda8dc17643746f17a661790763de1113cec93893b34a
balanced_geometry_solver_sensitivity_20260922_v1.tsv
SHA-256 5cc0b75a9a6e395fafb15ed6b3e243b0e7f12c383ef1ab6981a6ef1218c9dd1a
```

Durable execution evidence is in
`balanced_geometry_crossfit_jobs_20260922/execute_attempt1/`.
The independent audit is reproducible with
`benchmarks/scripts/verify_balanced_geometry_crossfit.py --outdir
benchmarks/results/balanced_geometry_crossfit_20260922_v1`.
An optional fresh `--sensitivity-output` path saves the diagnostic table.

The research goal remains active. Negative-control qualification, calibrated
uncertainty, depth/seed/replicate stability and whole-workflow resource gates
remain incomplete. Main and the production default are unchanged.
