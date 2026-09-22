# Blacklist and output-layout effects — 22 September 2026

The all-site blacklist intervention changes only excluded profiles. It does
not explain the three expected-shape changes caused by the expanded output
layout. Those are reproducible properties of the current correction operator
under different output-region layouts. No candidate or threshold was changed.

This is diagnostic evidence from already examined naked-DNA replicate 2.
It is not an independent safety result or evidence of a footprint improvement.

## Frozen construction

`scripts/compare_blacklisted_control_layouts.py` copies the original all-site
command, preserves its complete region list, and adds the checksum-locked hg38
blacklist with a fresh output destination. The saved dinucleotide weight matrix
(DWM), count arrays, read shifts, numerical settings and disabled global scaling
remain fixed. Historical and expanded constructions are compared with the
blacklist held constant; the expanded construction is also compared with and
without the blacklist. The original archive and failed references are preserved.

Artifacts are under `benchmarks/results/blacklisted_control_layouts_20260922_v1/`.

| Artifact | SHA-256 |
| --- | --- |
| Configuration | `224034f981e08ce2f9ddf3d5b5db2352f3a111ed3432ff1ec7e9ade02264d250` |
| Review manifest | `076f5b2d8fb922d7779d7ee18ca8613ad22f5b9c430afd187cf2e8a0b60f3f7f` |
| Job specification | `f4043b6e4df598b985523426e738d409eda0115e72829b1c98b1a0cbae1a40db` |

The durable launch/exit records are under
`blacklisted_control_layouts_jobs_20260922/execute_attempt1/`.
Execution exited zero in 52.690 seconds, including input verification and
review; correction itself reported 22.815 seconds. Neither timing establishes
complete-workflow central processing unit (CPU) or memory performance.

## Results and denominator

All 4,200 rows remain; 35 motifs/profile windows are explicitly blacklisted.
There are 2,880 positive expected profiles, compared with 2,892 without the
blacklist. No blacklisted window has positive expectation after exclusion.
Observed profiles change at 12 blacklisted sites and at zero eligible sites.
Expected proportions change at zero eligible previously positive sites when
adding the blacklist. Empty profiles are not evidence of biological negatives.

All 269 historical positive windows remain covered. The same three normalized
expected profiles differ between historical and expanded output layouts:

| Task | Motif interval | Maximum absolute proportion difference |
| --- | --- | ---: |
| HepG2 CTCF | chr10:28991925–28991940 | 0.028374 |
| K562 CTCF | chr9:98416493–98416508 | 0.006280 |
| K562 FOXA1 | chr9:98416381–98416389 | 0.056100 |

Window phase and other boundary effects remain components of the layout
intervention; this does not isolate them individually. Shape comparisons use
relative tolerance 1e-5 and absolute tolerance 1e-7, without claiming equality
of globally scaled expected counts.

## Verification and decision

Twelve focused tests passed in 2.44 seconds. Independent direct bigWig interval
reconstruction checked 25,200 profiles across six tracks, both expected-shape
comparison masks and observed-count equality. Both saved profile archives
match the independently reconstructed arrays exactly. Frozen-artifact
verification is recorded separately from the successful process exit.

The original expanded-reference overlap gate remains failed, and the new
manifest explicitly sets `qualified_for_rescreen: false`. Do not silently
substitute these tracks into a previous safety claim. The reference diagnosis
now distinguishes exclusion from layout effects; further repetition of the
same intervention is unnecessary.

Proceed to low-count normalization diagnosis on training and replicate-1
controls, freezing the descriptive components before evaluation. Any future
revised correction/reference contract needs its own frozen paired analysis;
do not borrow safety from an older operator. No promotion gate is newly passed.
Main and final untouched controls remain unchanged.
