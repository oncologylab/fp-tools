# Global expectation: training parity and coverage

The new genomic expectation resolves the overlap inconsistency of the
conditional-window representation. It does not yet establish better
footprint detection. Coverage at 10 million fragments is a substantial
limitation and must remain visible before any discrimination evaluation.

This initial stage uses only chr1–15 coordinates from HepG2 and K562
replicate 1 at 10 million fragments, seed 2026. Label columns were not
deserialized. No model was fitted, geometry selected, discrimination score
evaluated, or control/holdout opened. Existing coefficients remain frozen.

## Operator and numerical verification

Each strand's expected count at a genomic position is its sequence propensity
multiplied by the observed-to-propensity ratio in symmetric 50–100 base-pair
outer flanks. Counts and propensities were freshly extracted over 401
positions; the central 201 positions have a complete background halo.
Invalid expectation remains NaN with an explicit validity mask.

The runner checked an independent scalar calculation at every site's center
on both strands and smaller halo windows for the first 64 coordinate-order
sites per cell. Observed central profiles agree exactly with the archived
corrected extraction after restoring genomic orientation. Adjacent overlapping
windows, including duplicate centers, have no observed, expected, or validity
conflicts.

| Cell | Training sites | Checked overlap pairs | Shared supported pair-positions | Shared unavailable pair-positions | Conflicts |
| --- | --- | --- | --- | --- | --- |
| HepG2 | 20,019 | 1,918 | 199,280 | 73,124 | 0 |
| K562 | 14,987 | 1,181 | 133,413 | 40,091 | 0 |

Pair-positions are not unique genomic bases. This verifies the recorded
profiles and overlap subset, not a complete exported genomic track pipeline.
The separate scoring primitive also passes synthetic count-to-residual-to-
footprint parity with a 140-base-pair halo. It uses the conventional geometry
without reintroducing per-window total rescaling. Its real-data discrimination
and safety have not been evaluated.

## Coverage limitation

An independent prefix-sum calculation reconstructs the expected-validity
mask from raw counts and sequence availability and matches every saved
position mask. The conventional geometry requires only its fixed central
and shoulder positions; full aggregate profiles require all 201 positions.
Neither support definition was selected using labels or performance.

| Cell | Positive observed-count sites | Supported motif centers | Complete conventional scoring support | Complete 201-position profiles |
| --- | --- | --- | --- | --- |
| HepG2 | 18,052 | 11,709 | 6,453 | 3,861 |
| K562 | 13,421 | 8,984 | 5,248 | 3,246 |

At the central base, all alignments and sequence windows are available.
The missing expectation is explained by empty plus/minus-strand background
windows: HepG2 has 5,812 empty plus, 5,910 empty minus, and 3,412 empty both;
K562 has 4,358, 4,299, and 2,654, respectively. These reason counts overlap.
Positive cuts somewhere in a profile do not guarantee a supported local
background on both strands.

All 21 transcription-factor tasks remain in `tasks.tsv`, with offset-specific
support in `positions.tsv`. CTCF has 675/711 sites supporting the conventional
geometry in HepG2/K562, but only 377/400 complete profiles. HepG2 CTCF therefore
cannot meet the 200-sites-per-class requirement on complete profiles even
before labels or matching are considered. Class balance, accessibility
balance, and power for other comparisons have not been established.

Do not fill unavailable expectations with zero, silently discard exclusions,
or transfer earlier conditional-profile accuracy/safety claims. Any revised
low-count background model must be frozen separately and developed using
permitted training evidence before tuning or controls are evaluated.

## Artifacts and validation

Profiles: `benchmarks/results/global_parametric_training_profiles_20260922_v1/`.
Coverage: `benchmarks/results/global_parametric_training_coverage_20260922_v1/`.

| Artifact | SHA-256 |
| --- | --- |
| Profile configuration | `d903eec175130993afe35e6f3c64701e25399d1e29f6d2b98f8750825849b328` |
| Profile result manifest | `a2f3401aaf928ff3f8379449d5f8a66710077c15c992687990be4c16547e202d` |
| Coverage configuration | `eaf66c19ef609b1a676549d5b04e9ceb3c39f5345bbc49d71de27a3118224648` |
| Coverage result manifest | `f22f9c93852d3ea142f14894f04af788c3cb98adee2f84b1c07dda3322b2d883` |

Profile execution has a durable exit-0 receipt after 376.52 s; independent
full verification exited 0 after 6.08 s. These durations are not evidence
for the complete-workflow computational gate. Coverage execution and its
independent `--verify-only` completed successfully. Nine initial profile/
expectation tests passed in 1.50 s; fifteen scoring/profile/expectation tests
passed in 1.43 s; nine coverage/scoring tests passed in 1.38 s. These suites
overlap and their counts must not be added together.
The final combined run of all four suites passed all 18 tests in 1.51 s.

Main, the production default, original artifacts, and locked holdouts remain
unchanged. No promotion is warranted.
