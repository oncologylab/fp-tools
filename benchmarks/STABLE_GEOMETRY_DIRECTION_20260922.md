# Chromosome direction and fixed-center protection

The stable detector's numerical gains are not a general footprint improvement.
CCCTC-binding factor (CTCF) has the clearest agreement between detector gains
and fixed-center protection. K562 MEF2D has a repeatable discrimination gain
across training chromosomes, but much weaker agreement with the fixed-center
shape criterion. Zinc-finger tasks remain inconsistent or worse.

## Frozen diagnostic

Retained all 21 tasks, seven correction arms and fifteen training chromosomes.
Compared four prespecified effects on identical held-chromosome sites:

- Normalization: normalized fixed geometry versus unnormalized fixed geometry.
- Detector: selected versus fixed geometry with the same normalization and correction.
- Correction: each correction versus dinucleotide weight matrix (DWM), both
  using their selected normalized detector.
- Combined: selected normalized detector versus unnormalized fixed DWM.

The correction comparison holds the detector-selection procedure fixed;
fitted geometries can differ between corrections. It therefore does not isolate
only expected-profile changes with identical detector coefficients. The fixed
detector comparison in the earlier full-pool study remains the complementary
control for that question.

Each effect has 147 task/arm summaries, including unavailable tasks. Summaries
give equal-chromosome means, signs, extremes and the range of means obtained
by omitting each chromosome. These omission ranges are descriptive sensitivity
ranges, not confidence intervals. Fitting sets overlap across folds, and
matching used the full training pool; no independent-validation claim is made.

Recomputed a common protection criterion from the saved profiles: center
offsets −5 through 5 must be negative and below both shoulders (−40 through
−15 and 15 through 40), by more than the fixed 1e-10 numerical margin.
This criterion does not change when a detector selects a different geometry.
No score-selected subsets were used to make a stronger-looking aggregate.

## Detector-effect findings

AUROC means area under the receiver operating characteristic curve.
All rows below have fifteen evaluable chromosomes.

| Task and correction | Mean AUROC gain | Positive / negative / tied chromosomes | Fixed-center protection | Gain and protection together |
| --- | ---: | ---: | ---: | ---: |
| HepG2 CTCF, DWM | +0.066661 | 13 / 2 / 0 | 12 | 11 |
| K562 CTCF, DWM | +0.027441 | 9 / 4 / 2 | 11 | 6 |
| K562 CTCF, conditional parametric | +0.054818 | 10 / 3 / 2 | 12 | 7 |
| K562 MEF2D, DWM | +0.079675 | 12 / 2 / 1 | 4 | 3 |
| K562 MEF2D, conditional parametric | +0.058215 | 13 / 1 / 1 | 5 | 4 |
| HepG2 ZNF362, radius 500 | −0.005597 | 6 / 7 / 2 | 6 | 3 |
| HepG2 ZNF384, radius 500 | −0.007211 | 3 / 5 / 7 | 8 | 1 |
| K562 ZNF384, radius 500 | −0.010588 | 2 / 6 / 7 | 5 | 0 |

The DWM detector-gain omission ranges stay positive for HepG2 CTCF
[0.053333, 0.080842], K562 CTCF [0.019071, 0.037437] and K562 MEF2D
[0.072866, 0.092509]. Thus none of these mean gains depends solely on
one chromosome. This does not establish significance or stability across
biological replicates and sequencing depths.

K562 MEF2D passes its selected-geometry protection criterion in ten chromosomes
in the prior report, versus four under the fixed-center criterion here.
The selected geometry often has a narrower center. This difference warrants
an exact-shape and negative-control audit; it neither proves the absence of a
narrow footprint nor establishes the locked aggregate-improvement gate.
HepG2 MEF2D remains underpowered and uses the conventional fallback.

The mean of per-chromosome relative area under the precision-recall curve
(AUPRC) gains is 8.51% for K562 MEF2D/DWM. This differs from the 7.86% ratio
of mean AUPRCs in the previous report because averaging and division do not
commute. Neither statistic substitutes for the prespecified promotion metric.

## Verification and artifacts

- Twelve focused tests passed in 2.60 seconds.
- Durable execution exited zero in 16.429472 seconds; this is diagnostic
  elapsed time, not computational-gate evidence.
- Separate checksum/code/runtime verification passed.
- Independent recomputation checked all 8,820 paired effects, 2,205
  chromosome/arm fixed-center audits and 588 task/effect summaries, including
  missing-class rows, signs and explicit leave-one-out means. Shape checking
  used the equivalent least-squares normalizer.
- All three maps were visually reviewed, covering all 21 tasks for DWM,
  conditional parametric and radius-500 correction. Other arms remain in
  the complete tables. Color scales differ by panel and are labeled; dots
  indicate fixed-center protection, not statistical significance.

Ignored output directory:
`benchmarks/results/stable_geometry_direction_20260922_v1/`.
Execution evidence:
`benchmarks/results/stable_geometry_direction_jobs_20260922/`.

SHA-256:

```text
configuration.freeze.json       009ee65deac8a2a63b815d14003658c90f4546cb8994b66f3f7c342850d6b343
direction_manifest.json         cc0457fce3eef49914c1c696a9ff5ad184a97a6b189bc9953d90485473ef61e0
all_task_detector_direction.pdf c3d61b1c087e68cf57099155552819415aca11d99c490bb6898b94b8e4e3972b
```

The immutable manifest retains its pre-review visual flag; the completed
visual review is recorded here. No new holdout labels or final control were
opened. No candidate was promoted, uploaded or integrated into main.

Next freeze full-training geometry artifacts under the unchanged stable rule,
retain every task and fallback, and qualify their exact scoring operators on
permitted development negative controls. Use the existing safety protocol as
the starting point; do not transfer safety from a different correction or
detector. The narrower MEF2D profile needs explicit shape interpretation.
Full-workflow output parity and resource qualification remain incomplete.
