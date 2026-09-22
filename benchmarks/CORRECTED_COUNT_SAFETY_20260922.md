# Corrected-count development safety reassessment

This is a post-hoc reassessment of previously examined naked-DNA replicate 2.
It is not final untouched safety, candidate optimization, or promotion evidence.
No new holdout labels or final control replicate were read. Production and the
published manuscript are unchanged.

The corrected reader reextracted the original development and control sites,
retaining the frozen sequence model, scoring operator, dispersion, cut shift,
mapping filters and site order. All four profile stages completed and passed
independent verification. The comparison configuration was frozen, followed by
both threshold analyses, before corrected control scoring. The completed
comparison passed independent verification.

## Results and interpretation

There are 21 cell/transcription-factor tasks, each with 200 original control
sites. Corrected counting changed seven HepG2 sites and three K562 sites:
motif-window cut assignments changed from 2,845 to 2,835 and from 2,212 to 2,209,
respectively. These are window assignments, not unique library fragment counts.
No CTCF control profile changed. With the historical thresholds, candidate
false-positive counts and valid-site denominators are unchanged in all 21 tasks.

| CTCF context | Informative sites | Candidate calls | DWM-comparable sites | Paired Wilson upper 95% bound | Paired increase over DWM |
| --- | ---: | ---: | ---: | ---: | ---: |
| HepG2 | 122 | 0 | 121 | 3.0771% | 0 percentage points |
| K562 | 115 | 1 | 114 | 4.8020% | 0.8772 percentage points |

DWM means dinucleotide weight matrix. The table is identical for fixed and
recalibrated thresholds. Both CTCF absolute false-positive screens pass and
neither context has a flagged central motif residual. The HepG2 parametric
threshold changes slightly, from approximately 0.211875 to 0.211788; this does
not change any control call. K562's threshold is unchanged at this precision.

**No task receives a qualified relative-safety pass.** The original calibration
rule includes informative tuning negatives with zero DWM expectation. This
affects 9 of 989 CTCF negatives in HepG2 and 5 of 895 in K562. The prespecified
corrected protocol therefore withholds the relative-safety claim, even though
the CTCF control rates remain numerically favorable. This is an additional
baseline-calibration qualification, not numerical deterioration caused by
correcting the read counts. Most other tasks also lack sufficient comparable
DWM control support; zero-filled baseline profiles cannot establish safety.

A future development candidate must define and freeze valid baseline
calibration support before evaluating its thresholds. Do not retroactively
change this comparison or choose the more favorable threshold analysis. Final
untouched controls remain locked until all prerequisites pass.

## Reproduction and validation

- Profile outputs: `benchmarks/results/corrected_safety_profiles_20260922_v1/`.
  Profile execution exited 0 after 371.15 seconds; `--verify-only` passed.
- Safety outputs: `benchmarks/results/corrected_safety_reassessment_20260922_v1/`.
  Configuration SHA-256:
  `0b2fab46058990d8536fa0897f8699dbe15c77681f449312c68c259c74364452`.
- Threshold freeze SHA-256:
  `a9e196371d0878d619e8ba18b81827099e9c329c77de8928e2b166783f5f6ee7`.
- Verified safety manifest SHA-256:
  `4fd462c4852b88776784bec0e039b9cd62bfd2cdbbbd4001e5ec6338b1f0143a`.
- Safety scoring exited 0 after 4.56 seconds; separate `--verify-only` passed.
  These execution durations are not complete-workflow resource benchmarks.
- Twenty-five focused safety/profile tests passed in 2.09 seconds; twenty
  documentation contracts passed in 0.38 seconds; strict MkDocs passed in
  2.83 seconds. No native desktop validation applies to these benchmark scripts.

Per-site counts and validity/membership changes, both threshold tables, residual
estimates and aggregate curve data are retained with checksums. This report
does not certify plot review, functional uncertainty, occupancy-leakage
exclusion, broader-family improvement or complete-workflow resource gates.
