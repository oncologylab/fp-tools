# Verified corrected-count 96-run reassessment

The corrected extractor preserves the narrow CTCF conclusion and the broader
failure to rescue difficult transcription-factor families. It does not
establish a general replacement or satisfy promotion gates. All historical
outputs remain unchanged.

This is a post-hoc implementation reassessment of the **previously examined
chr19–22/X** fixed-candidate study. It is not the newer chr16–18 geometry/spline
tuning study, an untouched evaluation, or permission to retune against these
test results. No new labels were opened and no candidate was refitted.

## Design and verified scope

All 96 cell/replicate/depth/seed runs completed, using fresh corrected cut
extraction, unchanged frozen sequence coefficients, site order, matched
cohorts, scoring operator, dispersion, and principal bootstrap identities.
Each of 21 tasks retains 48 expected runs within its cell. Twenty matched
tasks are scorable; HepG2 HMGA1 remains explicitly unscored. Seven tasks are
powered under the original requirement of at least 200 sites per class.

The candidate is frozen parametric bias with conventional signed-deviance
geometry. The comparator below uses the sample/depth-matched dinucleotide
weight matrix (DWM) expectation with the same residual geometry. It is not
a whole-production-pipeline or learned-spline comparison.

Full-depth values average three biological replicates. Seeds are not treated
as independent biological replicates. AUROC means area under the receiver
operating characteristic curve; AUPRC means area under the precision-recall
curve. Every task remains in `per_tf_summary.tsv.gz` and the gate ledger.

| Powered task | Full-depth mean AUROC gain | Relative AUPRC gain | Canonical protection / 48 runs |
| --- | --- | --- | --- |
| HepG2 CTCF | 0.111969 | 27.69% | 46 |
| HepG2 FOXA1 | -0.049606 | -6.65% | 0 |
| HepG2 ZNF362 | 0.051473 | 4.85% | 0 |
| HepG2 ZNF384 | 0.014600 | 0.84% | 0 |
| K562 CTCF | 0.130023 | 23.36% | 48 |
| K562 MYC | -0.045643 | -7.72% | 0 |
| K562 ZNF384 | 0.046494 | 6.80% | 0 |

The table is a powered-task description, not a selected-denominator mean-gain
gate. Numerical occupancy discrimination without canonical protection is
not a footprint improvement.

## Extraction effect and remaining failures

The largest absolute CTCF AUROC change across runs is 0.000752 in HepG2 and
0.000729 in K562. Full-depth mean gains changed from 0.112090 to 0.111969 and
from 0.129794 to 0.130023, respectively. Their shape-pass counts are unchanged.
Small-cohort metrics can change more: the largest matched-task AUROC change
is 0.020000 for K562 TCF7L2 and the largest AUPRC change is 0.061111 for HepG2
ZNF558. Neither task is powered. These are numerical extraction effects, not
new estimates of statistical significance.

The same HepG2 CTCF replicate-2 10-million-read runs remain shape exceptions:
seed 2026 has central contrast -0.056385 but fails the prespecified protection
criterion; seed 2028 has central contrast +0.001540. Both retain positive
shoulders. K562 CTCF replicate-3 at 10 million reads, seed 2026, has a paired
chromosome-block AUROC gain interval [-0.014312, 0.056015] and relative AUPRC
gain interval [-0.011738, 0.125445]. It cannot pass the positive-interval gate.

All six CTCF full-depth replicate intervals remain positive versus the
sample-matched DWM scorer. This does not erase the low-depth failure.
The maximum absolute current-count standardized mean difference across
CTCF runs is 0.15145 in HepG2 and 0.15863 in K562. Fixed cohort membership was
preserved to isolate extraction; this reassessment did not rematch away
depth-dependent balance limitations.

The corrected ledger retains 21 tasks and 746 shape exceptions across
scorable matched task/runs. Difficult-family improvement, untouched safety,
occupancy-leakage exclusion, calibrated functional uncertainty, complete
output compatibility, and whole-workflow computational gates remain
unestablished. Separate corrected safety and spline reports retain their
failures. Production defaults remain unchanged.

## Verification and artifacts

- All 96 cases: supervised exit 0 after 5,738.28 s.
- Independent full input/output verification: exit 0 after 937.31 s.
- Comparison: exit 0; six output checksums verified independently.
- `pdfinfo` confirms 21 A4 pages. All 21 rendered pages were visually reviewed
  for layout and profile interpretation. CTCF shows the clearest protection;
  difficult-family profiles remain weak, enriched, or noisy. HMGA1's absent
  matched panels reflect its unscored support, not a missing result.
- Forty-one extraction, ledger, comparison, and new expectation/contract
  tests passed in 4.04 s. These do not certify native desktop builds.

| Artifact | SHA-256 |
| --- | --- |
| Corrected result manifest | `71235831ce3559b5e5c81a5b30e4564879cf7afb11593cc69c221e9eb3d8b979` |
| Comparison manifest | `9e49fec131010eecfd7f2c0961f01618ee4a33d2f1cd226842689543888f9c2a` |
| Aggregate PDF | `8156344d3943e1a4242312765ba82485c7c5b4fdbb87a31395ef234feffbf007` |
| Corrected gate ledger | `4f405fca8756c255fba29abcf611a1e8ae0267e4b2386d0409e7cc18b1fc394b` |

Results are in `benchmarks/results/corrected_count_reassessment_20260922_v1/`,
`corrected_count_comparison_20260922_v1/`, and
`corrected_count_gate_ledger_20260922_v1/`. The latter two paths are also under
`benchmarks/results/`. Nothing was uploaded, promoted, or merged into main.
