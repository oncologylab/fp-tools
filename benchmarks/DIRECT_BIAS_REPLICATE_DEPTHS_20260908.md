# Completing the direct-bias replicate/depth panel

The existing stability experiment covers 36 runs: all depths and seeds for
replicate 1, plus the full-depth endpoint for all three biological replicates
in HepG2 and K562. The new input-preparation plan covers the remaining 60
combinations for replicates 2 and 3. These inputs are still being prepared;
no expanded scientific result or promotion is claimed.

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

## Remaining preparation and evaluation

The controller creates missing subsets, then fits conventional DWM separately
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

After all 60 new input sets pass verification, the 96-run evaluation must
freeze positive DWM expectation support across every run before scoring.
Coverage exclusions must retain complete existing matched pairs; zero observed
counts must remain in the cohort. The previous 36 observed/expected profile
arrays should be reused on the new common support, with new extraction only
for the missing 60 runs. Any changed cohort requires a new comparison rather
than combining metrics from different site sets.

These are development stability diagnostics on previously examined chr19–22/X.
Input preparation does not resolve the difficult-family, uncertainty,
untouched-safety or full workflow resource gates. Production DWM remains the
default and the final manuscript remains frozen.

## Validation

All 60 focused sampler, preparation, study and documentation-contract tests
passed. The strict MkDocs build passed, and the browser audit passed all 32
pages at three widths. The four count stages completed with their source
identities intact. The full preparation/evaluation matrix remains incomplete.
