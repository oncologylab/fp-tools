# Direct-bias sequence inference resources

The frozen parametric model used 47.4% of production DWM's sequence-inference
CPU time in HepG2 and 47.7% in K562 on identical loci. Absolute worker peak
memory was nearly equal, about 226–228 MiB. This is a component measurement;
the full workflow runtime and memory promotion gates remain unmeasured.
No scientific model, detector, production default or manuscript changed.

## Frozen measurement

Calculation commit: `feb705f`. Plan:
`manifests/direct_bias_inference_benchmark_v1.json`.

The plan selects 512 sites per cell by the lowest seeded coordinate SHA-256
from the existing all-run expectation-coverage cohort. Selection uses cell,
TF and genomic coordinates; it does not use labels or numerical model scores.
These are previously examined chr19–22/X development sites. The site tables,
sequences, archived reference propensities, model files and calculation code
were frozen before timing.

Both scorers predict the same 201 genomic positions on each strand. The
candidate is the unchanged mitochondrial LOG21 +4/−5 L2 0.0001 coefficient
ensemble. The comparator uses the actual compiled production DWM scorer and
the two checksum-pinned local 10M seed-2026 replicate-1 DWM model artifacts.
It does not substitute the Python control-reference implementation. Each DWM
sequence is trimmed to produce exactly the requested output positions.
The respective model contexts remain 21 and 25 bp.

There are six repetitions per cell and method, with alternating method order.
All 24 workers have distinct process IDs and are pinned to logical CPU 0, with
the numerical-library thread environment set to one. Both import the same
modules and load the identical pre-extracted sequence bundle before timing.
Each warms up on the first 16 frozen sites; the measured pass then covers all
512 sites. All repetitions are retained.

The shared Linux x86-64 host reports an AMD EPYC-Milan processor. The frozen
environment records Python and numerical-library versions, thread settings
and the native extension checksum. Unrelated host workloads were present;
these are measurements on this host, not a dedicated-hardware guarantee.
CPU time is reported alongside elapsed time to expose scheduling effects.

## Results

Means use all six repetitions. Peak memory is the largest absolute Linux
`VmHWM` among each method's six fresh workers, including common imports,
sequence inputs, model, predictions and output serialization.

| Cell | Method | Mean elapsed seconds | Mean CPU seconds | Maximum worker MiB |
| --- | --- | ---: | ---: | ---: |
| HepG2 | Production DWM sequence | 2.8340 | 2.8330 | 228.12 |
| HepG2 | Frozen parametric sequence | 1.3433 | 1.3424 | 226.40 |
| K562 | Production DWM sequence | 2.8151 | 2.8140 | 228.09 |
| K562 | Frozen parametric sequence | 1.3414 | 1.3410 | 226.27 |

| Cell | Candidate / DWM elapsed ratio | CPU ratio | Maximum worker memory ratio |
| --- | ---: | ---: | ---: |
| HepG2 | 0.47398 | 0.47385 | 0.99247 |
| K562 | 0.47652 | 0.47654 | 0.99202 |

Parametric elapsed ranges were 1.3337–1.3588 seconds in HepG2 and
1.3338–1.3523 in K562. DWM ranges were 2.8244–2.8437 and 2.8016–2.8433,
respectively. Elapsed and CPU means agree closely. The common process baseline
dominates these memory measurements; the small memory difference does not
establish a practical reduction in full-analysis peak memory.

Every run produced all 205,824 finite strand/position values. Each method's
numerical output is identical across its repetitions. Every parametric run
also reproduces its archived strand propensities after reversing the saved
motif orientation, within absolute and relative tolerance 1e-12. DWM and
parametric values have their respective production log-score units; equality
between methods is neither required nor claimed.

## Scope and advancement decision

The timed component excludes fitting, BAM extraction, FASTA I/O, DWM local
regression, expected-count normalization, residual calculation, detector
scoring and track output. Model loading is recorded separately. Its ratios
cannot be entered as end-to-end promotion evidence. The benchmark uses
absolute fresh-process memory, avoiding both cumulative child peaks and
differences of process high-water marks in older engineering/fit helpers.
Those older measurements are not used to certify this candidate.

The [36-run scientific study](DIRECT_BIAS_STABILITY_20260908.md) remains the
source of discrimination, protection, balance and candidate-screen safety
evidence. The current gate disposition is:

| Requirement | Current evidence and disposition |
| --- | --- |
| Mean AUROC gain ≥0.03 | Qualified matched full-depth development means exceed 0.03; this is not untouched holdout evidence. |
| Relative AUPRC gain ≥10% | Fails the general development comparison: matched full-depth means are 6.76% in HepG2 and 7.59% in K562. |
| Protection/separation in two difficult families | Fails: CTCF is the only powered task with canonical protection in all runs, and is a positive control. |
| Improvement in two new cell holdouts | Not evaluated; SK-N-SH and GM23338 labels remain locked. |
| CTCF non-regression and principal intervals | Supported in the examined development comparisons; all six full replicates have positive paired intervals against raw and both DWM comparators. |
| Independent naked-DNA candidate screening | Passes for the exact CTCF direct operator on informative replicate-2 support; sparse controls do not certify all TFs or depths. |
| Untouched final safety and occupancy leakage | Final replicate 3 remains unopened; absence of a large residual flag does not prove absence of occupancy leakage. |
| Depth, seed and biological-replicate stability | Positive CTCF direction in the available 36 runs; replicates 2–3 still lack the depth/seed factorial panel. Cut-total balance drifts in 11 runs. |
| CPU ≤2× and peak memory ≤1.5× DWM | This inference component is within both numerical ratios; complete workflow measurements are still missing. |
| Model artifact ≤25 MB | The prior guarded-load audit measured 10,007 bytes for the NPZ and required JSON sidecar. This does not include data or process memory. |
| Calibrated functional uncertainty; GP/hybrid benefit | Not established. Pointwise aggregate intervals and geometry scores do not establish calibrated functional coverage or binding probabilities. |

General advancement remains unsupported by the existing development failures.
Faster sequence inference does not rescue the difficult-family result or
justify opening new labels. No new Box upload is warranted by this resource
diagnostic alone; the existing narrow CTCF reports retain their stated scope.

## Artifacts and verification

Ignored output directory: `results/direct_bias_inference_20260908_v1/`.
It contains the frozen input bundles and site lists, 24 worker measurements
and prediction arrays, the complete measurement TSV, component ratios and
their checksum manifests. The configuration records 15 source input files.

- Configuration ID: `3143dfee8f643d687690c9e9e2b90ac8a47637cb1e4c1345861eed6f48d06b9a`.
- Configuration JSON SHA-256: `05973ed4c19ba375d32875ad783f9b800575bfa9ffb25f22e32dac6e9a14eee4`.
- Result ID: `7db527cf9338c6d42d240cd2303ea9afe2c05813b232f88d670c311d8471b353`.
- Result JSON SHA-256: `e48ee22355d6a2fc23670416d08dcecf4195b46201002763bf09ff6dfc1f8d6c`.
- Measurement TSV SHA-256: `b6cfd2150de7717216cb1ad3d8b92fff5d5edaf6a7608d9fc81287d3914112d0`.
- Ratio TSV SHA-256: `0a559843a73f322abceb15cee77e1d4c50458a41649189c6ac6ac2b2dc837533`.

All inputs, frozen code and nested worker outputs passed `--verify-only
--code-ref feb705f` without repeating a timing run. Fifty focused benchmark,
strand-profile, parametric-model and documentation tests passed, including
13 new resource-benchmark tests. The strict MkDocs build and the 32-page
browser audit at three widths passed before measurement.
The final documentation edits also passed all 20 documentation-contract tests
and another strict build. After measurement, `main` remained clean at
`214dc1a`, all four frozen manuscript package checksums matched the handoff,
and all 13 holdout `labels_read` flags remained false. The holdout manifest
SHA-256 remains `7d29af02e99507864e6cc194d4ec8e30df79f95bf4ab55d468c6c30eea70a354`.
