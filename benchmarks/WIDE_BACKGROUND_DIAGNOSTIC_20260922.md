# Wider local background: training-only diagnostic

Increasing background width reduces empty counts but does not eliminate
left/right rate disagreement. These are label-free development diagnostics,
not evidence of improved footprint detection or calibrated uncertainty.

## Frozen design and verification

The unchanged replicate-1, 10-million-fragment, seed-2026 alignments and frozen
parametric bias model were evaluated at all existing training coordinates on
chromosomes 1–15: 20,019 HepG2 and 14,987 K562 sites. Extraction covers 600 base
pairs on each side. Backgrounds include absolute offsets 50–100, 50–200, or
50–500, with both endpoints included. No labels were read, no model was refit,
and no window was selected.

Extraction exited zero after 641.24 seconds. Separate verification passed for
both stages and the final manifest. Shared central counts and bias values
agree with the earlier extraction. This is extraction time, not a full-workflow
computational benchmark.

The independent audit recomputed all 126 task/strand/radius summaries using
prefix-sum counts and log-space exposure sums. It retained all 21 tasks and
compared radii on common sequence support. Its separate artifact verification
also passed. The extraction and diagnostic stages were not rerun.

## Results

The table counts site–strand observations, not unique reads or independent
biological replicates. Overlapping windows can reuse the same reads.

| Cell | Common site–strand observations | Empty at radius 100 | Empty at radius 200 | Empty at radius 500 |
| --- | ---: | ---: | ---: | ---: |
| HepG2 | 40,038 | 11,722 | 4,876 | 1,464 |
| K562 | 29,972 | 8,655 | 3,747 | 1,152 |

HepG2 retains all 40,038 available site–strand observations; K562 retains
29,972 of 29,974 on common sequence support. The same denominator is used
at all three radii. All task-level exclusions remain in the output table.

For CCCTC-binding factor (CTCF), median strand-background counts increase
from 2 to 10 in HepG2 and from 2 to 8 in K562. However, median absolute log2
left/right exposure-adjusted rate disagreement remains approximately 1.0–1.14
across strands and radii. Its 90th percentile remains approximately 2.59–2.89.
Thus more background counts do not make the two flanks agree.

The reported rates use an explicit half-count prior for zero counts.
Disagreement can arise from sampling noise, sequence-model error, or actual
regulatory heterogeneity. It does not identify which cause dominates, and
must not be called evidence of bias overcorrection or transcription-factor
binding. No hypothesis-test significance or biological replication is claimed.

## Next evaluation and limitations

The separate fixed-radius training comparison holds the bias model, deviance
dispersion, conventional detector geometry and original matched pairs fixed.
It must retain complete pairs on common support, report changed power and
covariate balance, and separate numerical discrimination from canonical
central protection. Its execution freeze is separate from this diagnostic.
No new transfer labels, test-chromosome evaluation, or negative-control
replicates are authorized by this initial comparison.

## Artifact identities

Outputs are under `results/wide_training_profiles_20260922_v1/` and
`results/wide_background_audit_20260922_v1/`. Launch and exit records are
under `results/wide_training_jobs_20260922/execute_attempt1/`.

SHA-256:

- Extraction configuration: `50483e6776d99a162a84b49aece5db1ff78abc5ea14fad32e468f93542ae7937`
- Extraction manifest: `e2600056d9b24d660931de57a00e0a9c9bbf67878a2fba330825a12a9dded7ef`
- Audit configuration: `2695b20e5eaf51d5ad6ff03dee4ec1ab79bae1d56a30893f37b172a43dcb7ef6`
- Audit manifest: `3ca2c9e7315ed07c31efed0bef64f40572287e4404eee3096d7c5c08057c68fa`

Main, released packages, production defaults, and the manuscript remain
unchanged. No result establishes a promotion gate.
