# Completed ChromBPNet development reference audit

Verified on 2026-09-08, starting from research commit `888c6e2`.
This is development/diagnostic evidence. No promotion decision can pass from
these results alone. Production defaults and the final manuscript are unchanged.

## Completed stages and provenance

HepG2 regulatory training, all three prediction tracks, and the same-site
evaluation are complete. The 19 regulatory outputs, six prediction outputs,
five comparison outputs, recorded inputs and controller identities passed
SHA-256 verification. No completed stage was rerun. The external source tree
was clean at pinned commit `09938fdb4397ec0006510e5251e48920a505d4de`; its folds
are chr1–15 train, chr16–18 validation and chr19–22/X test.

The new report also verifies all inputs, outputs and canonical identities of
the completed K562 comparison. Input/result identities:

| Cell | Input ID | Result ID | Dispersion |
| --- | --- | --- | --- |
| HepG2 | `af3f30f68d78620c67f0b8a76d082e35de489814aa3f8aa7c15f8cd71e2b3a13` | `7634745f6c1eef990db1bbc0f8268a32274f930fe3980ea7a54f710f73a87404` | 0 |
| K562 | `ddea6ec5b1b0d61cba1484e066ad7b18e8c49fb08d18a79b79681ae26cc08797` | `79c814e45acc42a500ebc93a67f57112bc7dca0ae7e52d6b9a9c7949319a659b` | 1.6144318781609228 |

Those dispersions were already recorded before this audit. They are not a
common-operator cross-cell experiment and were not changed after viewing results.
Observed profiles use the frozen 10M-fragment, seed-2026 development artifacts.
These chromosomes were previously examined for other candidate models.

## What the comparison measures

The raw arm scores observed cuts with fixed center/shoulder geometry. DWM,
parametric and ChromBPNet bias arms total-match their expected profiles to each
observed motif-site profile, form signed-deviance residuals, and apply the same
geometry within a cell. The DWM arm is therefore a controlled expected-bias
reference, not a direct run of the production corrected-track footprint scorer.

ChromBPNet regulatory and no-bias arms instead score total-matched predicted
profiles. Their gains cannot be assigned solely to bias correction or a new
footprint detector. Residual aggregate panels use outer-flank normalization;
raw and predicted panels retain their stated cut units. The audit does not
overlay different units or equate unsigned separation with protection.

## Findings

On identical HepG2 CTCF sites (997 positive, 964 negative):

| Method | AUROC | AUPRC |
| --- | ---: | ---: |
| Raw geometry | 0.6817 | 0.6909 |
| DWM expected + deviance | 0.6477 | 0.6250 |
| Frozen parametric bias + deviance | 0.7149 | 0.7388 |
| ChromBPNet bias + deviance | 0.7253 | 0.7570 |
| ChromBPNet regulatory prediction | 0.7522 | 0.7580 |
| ChromBPNet no-bias prediction | 0.8392 | 0.8350 |

HepG2 parametric-minus-DWM AUROC has a chromosome-block 95% interval
0.0384–0.0857; relative AUPRC gain has an interval 12.69–22.87%. Against raw,
the intervals are 0.0109–0.0496 AUROC and 2.35–12.71% relative AUPRC. The
parametric model is below ChromBPNet bias here: AUROC interval
−0.0169 to −0.0011; relative AUPRC interval −3.51% to −0.35%.

Only CTCF shows both significant metric gains over DWM and raw and canonical
central protection in the two recorded comparisons. The parametric CTCF
standardized functional-separation gains over DWM are 16.77% in HepG2 and
11.60% in K562, below 25%. Neither is a difficult-family result.

Across the six eligible HepG2 tasks, mean parametric gains are only 0.0112
AUROC and 1.99% relative AUPRC; across three eligible K562 tasks they are
0.0053 and 1.25%. Six tasks in each cell are underpowered. HepG2 FOXA1 and
ZNF384 regress with parametric bias; FOXA1 has central enrichment. Larger
unsigned separation in FOXA1, MYC or ZNF384 does not establish protection.
No general replacement passes the prespecified mean-gain or difficult-family
requirements.

The common-site filter retains 99.7% of HepG2 CTCF positives but 96.4% of
negatives; K562 retains 100% and 97.5%. Retained CTCF motif-score standardized
differences are 0.482 in HepG2 and 0.463 in K562. Accessibility differences are
0.094 and 0.036. HepG2 FOXA1 and K562 MYC retain substantial accessibility
imbalance (0.805 and 0.832). These are explicit confounding diagnostics;
the audit does not rematch sites or retune a candidate on examined labels.

## Safety correction and remaining requirements

The older one-page K562 report cited
`development_tune/mito_log21_l2_1e-4_shift_4_-5/residual_safety_rep2/`.
That safety evaluator uses the factorization model's `expected_unbound`,
including fitted sample bias strength and a per-site flank spline. The direct
comparison instead uses total-matched `plus_expected + minus_expected`.
The cited factorization FPR cannot certify this direct scoring operator.
Historical artifacts remain unchanged; a new safety-passing report requires
exact-method verification. Partial-bias and factorization safety remain
evidence only for the methods actually tested.

Next: test the recorded direct operator on existing replicate-2 naked-DNA
profiles with validation-only frozen thresholds and audit residual motif
structure; preserve the distinction between candidate screening and final
untouched safety. Assess whether CTCF gains survive a separately frozen
diagnostic for motif-score/accessibility imbalance. No new SK-N-SH/GM23338
labels or final naked-DNA replicate-3 results were read by this audit.
The general development gates already fail, so this result does not justify
spending untouched holdouts or promoting a package method.

## Artifacts and verification

Ignored output directory in the research worktree:
`benchmarks/results/chrombpnet_development_audit_20260908_v2/`.
It contains 12-page HepG2 and nine-page K562 aggregate PDFs, all 126 method/task
metric rows, signed shape descriptors, site-balance diagnostics, method
decisions, limitations and output checksums. All 21 pages were visually checked.

Audit ID: `ab3aed5d523c9c0d2c99d7d5a4e2e7788fe859dc6a55f91c7733704f3d74c3ab`.
Per-TF metrics SHA-256:
`81b05c3e831f466c578f500064cf85cfd9b7dd74b203f28b2f865fcce11ce672`.

Focused verification: 27 tests passed, covering the audit, reference evaluator,
reference runner, checksum queue, existing renderer and residual-safety helper.
The audit tests cover coordinate identity, mismatched labels, mutation rejection,
unopened-cell rejection, signed protection and degenerate covariate balance.
Strict MkDocs build passed; the Playwright documentation audit passed all 32
pages at three viewport sizes with a dark system preference. Final PDF pages
were pixel-identical to the 21 pages visually reviewed during the audit.
