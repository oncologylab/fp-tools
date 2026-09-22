# Site-conditioned profiles and the full-workflow resource contract

The frozen research expectation profiles cannot be represented exactly by
one genomic expected-signal track. Overlapping windows assign different
expected values to the same genomic base, despite identical observed cuts.
This is a property of the conditional model, not an extraction defect.

For each strand and motif window, `_expected_from_scores` multiplies
sequence-derived probabilities by that window's observed cut total. Its
normalization denominator also depends on the window. The corrected residual
profiles inherit this dependence. Production `atac-correct` instead writes
`uncorrected`, `bias`, `expected`, and `corrected` bigWig tracks by default;
`call-footprints` consumes a genomic signal track.

Averaging, prioritizing, or selecting centers when exporting the site
profiles would define a different operator. The full-workflow resource
benchmark must not silently do this or substitute site-array timings for
equivalent production products. The existing resource protocol is unchanged.
This finding does not invalidate site-based comparisons or show that a
globally consistent candidate is impossible.

## Training-only evidence

The audit uses corrected, checksum-locked chr1–15 profiles. It selects valid
sites, the first sorted transcription-factor site at each distinct center,
and adjacent overlapping centers. Selection never uses labels, scores, or
disagreement. Motif orientation is restored to genomic orientation before
comparison. Expected agreement uses relative tolerance 1e-8 and absolute
tolerance 1e-10; observed counts require exact equality.

| Cell | Training sites | Overlapping pairs | Pairs with expected conflicts | Pair-position conflicts | Observed conflicts |
| --- | --- | --- | --- | --- | --- |
| HepG2 | 20,019 | 1,558 | 1,498 | 191,008 | 0 |
| K562 | 14,987 | 814 | 777 | 94,956 | 0 |

Pair-position counts are not unique genomic base counts. This subset
disproves exact single-track representability; it does not enumerate all
overlaps or certify other products. No tuning, test, transfer, or naked-DNA
control scores were evaluated. No model or production command changed.

## Required next work

Define a globally consistent correction and its products using permitted
training data, freeze it as a new research operator, and verify coordinate,
strand, coverage, and numerical behavior before tuning evaluation. Preserve
the current conditional candidate separately. A numerically different
operator needs its own scientific and exact-control evidence; existing
results cannot be transferred. Then freeze complete resource commands and
equivalent outputs, including fitting and serialization costs.

## Verification

Runner: `benchmarks/scripts/audit_site_conditioned_tracks.py`.
Output: `benchmarks/results/site_conditioned_tracks_20260922_v1/`.

- Configuration SHA-256:
  `aec5f02fe8b8ba3df10545e6fbb82458cfdaa81b0ac510e7e72a89059038ec28`.
- Result-manifest SHA-256:
  `e49740561b9442beecbfa17d0bc5053d96f00af2f1d5d83ac876f552da570c72`.
- Supervised execution: exit 0, 3.56 s; independent `--verify-only` passed.
- Five tests passed in 1.45 s: a constructive two-window conflict,
  reverse-strand restoration, invalid/nonoverlapping support, consistent
  track control, and malformed-input rejection.

No computational or scientific promotion gate is claimed. Main is unchanged.
