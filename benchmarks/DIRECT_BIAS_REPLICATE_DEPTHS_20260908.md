# Completing the direct-bias replicate/depth panel

The existing stability experiment covers 36 runs: all depths and seeds for
replicate 1, plus the full-depth endpoint for all three biological replicates
in HepG2 and K562. The new input-preparation plan covers the remaining 60
combinations for replicates 2 and 3. All sampled BAMs are complete and
independently verified. All 60 DWM corrections and the expanded 96-run
evaluation completed on September 10 and independently verified on September
12. See [the final stability report](DIRECT_BIAS_STABILITY_96_20260912.md)
for results, checksums and limitations. No promotion is claimed. The sections
below retain the preparation protocol and its earlier validation history.

## Frozen scope

Plan: `manifests/direct_bias_replicate_depths_v1.json`.
Preparation code: `604a4f2`; sampler-integrity fix: `66e64bc`.
Ignored output directory: `results/direct_bias_replicate_depth_inputs_20260908_v1/`.

The four full-depth BAMs and indexes match their existing frozen registry.
The genome, genome index, merged peaks, blacklist and archived correction
commands are checksum-verified. The preparation configuration pins code,
the production source tree, native extensions, Python/library versions and
a two-job concurrency limit. A kernel lock prevents duplicate controllers.

- Configuration ID: `e0546cdd1edf72d66f6eac15c94131a59bdc579634d3b2427e71ab9bb2bfb497`.
- Configuration SHA-256: `bede0ce8d5d43c3c86bd0491831783fbd31ae5ce91c7df425ecd111a656d40e7`.

The candidate coefficients, separate recorded dispersions, geometry, labels
and chromosome splits are unchanged. Preparation uses no ChIP labels and
cannot refit the candidate. SK-N-SH/GM23338 labels and final naked-DNA
replicate 3 remain outside its scope.

## Verified fragment counts

The original downsampler counts usable unpaired alignments or usable first
mates. Unmapped, secondary, supplementary, QC-failed and duplicate alignments
are excluded. No new numeric MAPQ threshold is added to this denominator;
MAPQ ≥30 remains the later profile-extraction rule.

Counting uses the pinned `pysam` SAMtools filter in an isolated process. Its
predicate matches the original Python function for all 4,096 combinations
of the defined SAM flag bits, including unpaired reads and MAPQ zero. The
filter uses the documented [SAMtools flag-expression semantics](https://www.htslib.org/doc/samtools.html#FILTER_EXPRESSIONS).

| Cell / replicate | Usable fragments | Count-stage ID |
| --- | ---: | --- |
| HepG2 2 | 96,480,203 | `655b1896111304c8e80dfeb10a4fcc3af43f9b732aaedea227b0202dd5fcbeab` |
| HepG2 3 | 105,920,129 | `a2766fcd51eb2ec82b97c2a054c26b0d32317db2f98b754b6c855649f8e1ea67` |
| K562 2 | 76,393,585 | `a898f10fee69d703dfb9c40098ce85218e6e3e41b8daf687dfc9bfce544ab5de` |
| K562 3 | 114,471,886 | `b3dca2ed44df3c2b6d245a739d9406c9d9ee0337521743e782d8b5f78dc6592e` |

Every library supports 10M, 25M and 50M targets. The query-name/seed BLAKE2b
rule preserves mates and nested depth subsets. These are expected target
depths under hash sampling; each completed subset records its actual selected
fragment count. A requested target above available fragments fails rather
than silently becoming a full-depth replicate.

## Verified sampling matrix

All four sampling stages completed on 2026-09-09 UTC under the unchanged
configuration above. Each contains 15 targets: 10M, 25M and 50M fragments
at seeds 2026–2030. Independent read-only verification passed for all 60
targets and 180 BAM, index and metadata files.

The checks covered every output checksum, the exact depth/seed matrix,
parent count/configuration identities, frozen preparation code, sampling
metadata identities and recorded source identity. BAM quickchecks passed;
all headers declared coordinate order and all indexes were readable.
Indexed mapped-alignment totals equaled twice the recorded selected-fragment
counts, with zero unmapped alignments. Completed count and sampling stages
were preserved.

The following SHA-256 values bind each sample's `sampling.json` in the
ignored output directory above, including all 45 nested output records:

| Sample | Sampling-manifest SHA-256 |
| --- | --- |
| HepG2 replicate 2 | `861ca1a4a923acec56c81bce79eef4a89c32904cf866d0ea367868700f5a1f7a` |
| HepG2 replicate 3 | `c400a9a39a4d5043027caf1da9108777dc0cce51f2de00a60a7813a738b33cac` |
| K562 replicate 2 | `8ff1ffaf5f451e2de5709e4c030c6387fea9ba3cbf786c6e8eb3c192bed57cab` |
| K562 replicate 3 | `b793cc7bdda2aa2c5fe34119e82540a5514af1d6b2f1b831076a4e16e41f80ca` |

Actual selected depths vary under the frozen hash-sampling rule. Across the
20 subsets at each nominal depth, the recorded counts span:

| Nominal depth | Minimum selected fragments | Maximum selected fragments |
| --- | ---: | ---: |
| 10M | 9,996,258 | 10,004,550 |
| 25M | 24,992,477 | 25,008,077 |
| 50M | 49,991,996 | 50,006,069 |

This verifies input integrity. All 60 conventional DWM corrections and the
final preparation registry must pass their own checks before the queued
96-run evaluation can calculate expanded footprint metrics.

## Remaining preparation and evaluation

With sampling complete, the controller fits conventional DWM separately
for each sample/depth/seed. It uses the recorded replicate-1 correction command,
including all four output tracks and four CPU cores. It does not substitute
the full-depth reused-bias comparator for a newly fitted low-depth comparator.
Each correction must finish successfully and provide nonempty tracks, its
bias artifact and diagnostic PDF before its checksum manifest is written.

Resume now validates the source hash, sampling identity, BAM and index hashes.
Six regression failures demonstrated the previous existence-only behavior.
The fix preserves complete valid targets and rejects changed or interrupted
ones without overwriting them. Old experiment manifests and already completed
scientific outputs are unchanged. Temporary BAM files are not completion
evidence.

The 96-run evaluator now enforces positive DWM expectation support across
every run and freezes both cells' cohort identities before new extraction
or scoring. Coverage exclusions retain complete existing matched pairs;
zero observed counts remain in the cohort. The previous 36 observed/expected
profile arrays are reused on the new common support, with new extraction only
for the missing 60 runs. Metrics are reused only if both site order and cohort
memberships are unchanged; changed support requires a new comparison.

Evaluator: `scripts/extend_direct_bias_stability.py`.
Plan: `manifests/direct_bias_stability_96_v1.json`.
Queue: `manifests/direct_bias_stability_96_queue_v1.json`.
The queue waits for the existing preparation to finish and binds the
evaluator's plan checksum. The plan in turn pins the evaluator code and
the original 36-run configuration/result checksums. Preparation source,
count, sampling, BAM/BAI, DWM command and nested output identities must
verify before any scientific calculation begins.

- Evaluator SHA-256: `07a8f9da64c4e3e028d90e2f0022ece79fa0d2fed3ba0cc3c1a3affe2a3414d4`.
- Evaluation-plan SHA-256: `b646b58e22ddce161e58894ac3beb06fcb37d1e64739ed432309a0322e0ce7e1`.
- Source 36-run configuration SHA-256: `7998850dd9eece77b668e0df04c8dbfa02a4bbca30adbdcee83170cb97add433`.
- Source 36-run result SHA-256: `94b4b887bef9431240967c7392970427a40d59470f8b1bc5ad069fa61bc0cfe0`.

This extension preserves the original mitochondrial LOG21 ensemble and
conventional geometry, including HepG2 dispersion 0 and K562 dispersion
1.6144318781609228. The four arms remain raw, fixed original DWM, sample/depth
DWM and frozen parametric expectation, all compared on the same sites.
Principal metric intervals cover the six full-depth endpoints and the six
10M seed-2026 runs. Repeated sampling seeds do not create new biological
replicates. Functional bands remain pointwise diagnostics, not a calibrated
functional uncertainty certificate.

These are development stability diagnostics on previously examined chr19–22/X.
Input preparation does not resolve the difficult-family, uncertainty,
untouched-safety or full workflow resource gates. Production DWM remains the
default and the final manuscript remains frozen.

## Interrupted sampling and recovery

The first preparation process returned exit 143 after its last partial-file
write at 2026-09-08 22:21:11 UTC. Its process was absent on inspection. The
dependent queue then returned exit 1 because the watched preparation had
ended without the required registry. No sampling stage or evaluation had
completed. The termination cause could not be established from the available
evidence.

All four source BAM/BAI pairs and count-stage identities reverified. The 30
partial BAMs, totaling 54,098,858,454 bytes, had no BGZF end-of-file marker.
They were checksum-inventoried and atomically moved into
`results/direct_bias_replicate_depth_inputs_20260908_interrupted_v1/`.
Device, inode, byte size and modification time were unchanged after each
move. All partial data and completed count manifests were preserved; no
partial BAM was accepted as a completed subset.

- Inventory ID: `c093024f991dafe04294e689b9b85a2975d7c12c7607dfdae684ddd5d5bd5ae7`.
- Preservation ID: `2c53229839f497bc425c06b7f771ab4be87e2cd9fc524fd3c1498f93a3f93c4c`.

Resumption reuses the four counts and the unchanged preparation/evaluation
commands. The separate `run_logged_experiment.py` supervisor pins a command
specification, preserves its exact bytes and records process identities,
stdout and terminal status. A new attempt directory keeps earlier evidence
intact. Detaching the supervisor from its short launcher is tested; it does
not establish why the earlier process stopped or prevent external shutdown.
A successful process exit still requires the scientific manifest verifiers.

## Validation

All 101 focused evaluator, stability, sampler, preparation, queue, study and
documentation-contract tests passed in 29.34 seconds. The small 96-run fixture
retains 36 archived extractions, reuses 18 evaluations on unchanged support,
recomputes the other cell's 18 evaluations after a coverage change, and
extracts/evaluates the 60 new runs. It checks both support freezes before any
new run, complete matched pairs, zero-cut retention, empty underpowered
interval tables, completed-run reuse and rejection of changed nested outputs.
The fixture's two-page A4 PDF was visually checked. These fixture numbers
verify execution behavior and are not biological results.

A read-only check of the real archived data verified all 36 extraction
coordinate identities and their sample DWM profiles on 9,361 HepG2 and
5,972 K562 sites. All 12 frozen scientific source files and the original
runtime match; no real extraction or scoring was repeated. The strict
MkDocs build passed in 4.16 seconds, and the browser audit passed all 32
pages at three widths. The four count stages remain complete; the recovery
above preserves them. The subsequent real sampling verification is recorded
above. The subsequent full preparation and 96-run evaluation are complete;
their independent verification and results are recorded in the final report
linked above.
