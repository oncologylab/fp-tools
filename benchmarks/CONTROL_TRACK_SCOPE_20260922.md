# Archived negative-control DWM reference has CTCF-only scope

The missing dinucleotide weight matrix (DWM) comparator coverage is a reference
scope limitation, not evidence of safe zero scores across transcription
factors. All 17,311 recorded expected-track bases lie within the historical
CCCTC-binding factor (CTCF) regions expanded by 100 bases. The artifact metadata
also names the earlier `ctcf_exact_geometry_safety` analysis.

## Frozen coverage-only audit

Inspected the already examined naked-DNA replicate-2 artifacts for all 4,200
sites and 21 tasks. Pinned the profile manifests, expected and observed tracks,
historical CTCF BED, corrected observed profiles, code and runtime before
calculating coverage. No candidate or threshold was changed or rescored.

The legacy sparse-track loader marks any readable in-bounds genomic window
as valid and substitutes zero for missing entries. Thus its `valid` flag
describes coordinate readability, not positive expected coverage. The recent
safety screen separately required positive expectation and therefore retained
these comparisons as unavailable; it did not count them as safety passes.

## Results

| Category | Site windows |
| --- | ---: |
| Coordinate-valid in the archived artifact | 4,200 |
| No expected-track entries | 3,931 |
| Recorded entries, but zero expectation throughout | 0 |
| Positive expectation at one or more positions | 269 |
| All 201 positions explicitly recorded | 0 |

Expected tracks are sparse: lack of dense entries does not by itself imply
an invalid prediction. Nor does an empty sparse window alone distinguish
an unprocessed region from a processed region with no output signal. The
historical region provenance and independent interval-containment check
establish the narrow scope of this particular reference.

Of 400 CTCF windows, HepG2 has 137 with positive expectation and K562 has
130. All 400 coordinates match the historical CTCF BED. After requiring
positive observed totals and matching observed-profile shape, 121/114 are
comparable. Two CTCF windows have positive fresh observed totals but no
expected-track entries; they remain unsupported, not presumed safe.

Only two of 3,800 non-CTCF windows have positive expectation: one K562 ARID3A
window with 38 recorded bases and one K562 FOXA1 window with 32. The former
has matching observed shape; the latter does not. No non-CTCF coordinates
match the historical CTCF BED exactly, but these two windows intersect its
expanded track support. This explains the one-site ARID3A comparison in the
previous screen. There are 2,301 positive-observed non-CTCF windows, so absent
expectation cannot generally be explained by absence of observed insertions.

Across all tasks, 2,301 windows have positive fresh observed counts but no
expected entries. The corrected CTCF false-positive failures remain valid
for their supported sites. Nothing here reverses those failures or qualifies
the unsupported tasks.

## Verification

- Twelve focused tests passed in 1.42 seconds, including explicit separation
  of absent entries, recorded zeros and positive expectation, malformed input
  rejection, data-root handling and existing safety support/call tests.
- Separate input/code/runtime/output verification passed.
- Archived observed and expected arrays exactly match direct track extraction
  after applying the archived zero-fill and motif-orientation rules.
- Independent interval-based extraction verified all 4,200 finite/positive
  base counts and categories and all 21 task summaries.
- Interval-derived expected totals agree with the stored float32 sums at
  relative tolerance 1e-6 and absolute tolerance 1e-8; maximum absolute
  difference was 4.8408285e-5. An initial float64-level tolerance was too strict
  for the archived float32 summation. Coverage classifications agree exactly.
- Every recorded expected-track interval is contained within the historical
  BED plus/minus 100; zero recorded bases fall outside it. This establishes
  track support, not recovery of an unrecorded historical command line.

Ignored artifacts:
`benchmarks/results/control_track_scope_20260922_v1/` contains per-site and
per-task tables and an immutable manifest. SHA-256:

```text
configuration.freeze.json 84265db9da3c51678e88538ea05e42ca0a772078e4bc8980534ec0ed128401ab
coverage_manifest.json    46540ec11fedde98e8579311d61f49485493c4f7c73203bcdb5c525da6a2a882
```

## Next work

Freeze a separate all-site control-reference construction using the recorded
BAM, reference genome and existing DWM bias model, with an explicit output
region contract. Check observed cuts and expected-profile shape on the old
overlap before any rescreening; changed support must not silently replace the
old evidence. Preserve all failed detectors, thresholds and original tracks.
This is a reference-coverage repair, not candidate retuning or a safety pass.
Investigating normalization at low counts on training and replicate-1 controls
remains outstanding. Main and final untouched controls remain unchanged.
