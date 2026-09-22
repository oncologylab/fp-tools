# Higher-depth replicate-1 input audit

The completed `NakedDNA_rep1.filtered.bam` is available and passes file/index
consistency checks. It can support a newly frozen development-control study,
but it is not a filtering-matched depth extension of `first10m.bam`.
No new profiles, labels or footprint scores were evaluated in this audit.

## Evidence

The historical preparation log links the full alignment to SRR13450940
(naked-DNA replicate 1), with 55,862,627 input read pairs and 55,530,574 pairs
after fastp filtering. The alignment header records Bowtie2 2.5.5 and samtools
1.24. Both alignments use the same reference dictionary and coordinate sorting,
proper-pair filtering, minimum mapping quality 30 and duplicate removal.
The full preparation additionally removes mitochondrial alignments and applies
read-level blacklist intersection filtering. The smaller preparation takes the
first ten million trimmed pairs and deduplicates after aligning that subset.
Deduplication at different input depths is another potential difference.

Fresh full-file flag statistics and independent index counts agree:

| Alignment | Mapped records | First-mate records | Second-mate records |
| --- | ---: | ---: | ---: |
| Full filtered | 78,358,765 | 39,179,382 | 39,179,383 |
| First 10 million input pairs | 14,957,194 | 7,478,597 | 7,478,597 |

Both contain zero secondary, supplementary or duplicate-flagged records.
Quick checks pass. Fresh full-file statistics also reproduce the archived
flag-statistic categories. All frozen input hashes were rechecked after the
audit. This verifies the current files, not a fresh rebuild from raw reads.
Raw-read checksums and the original reference-index build were not independently
reconstructed. Proper-pair flags do not prove that both mates remain present
after read-level filtering; the complete retained pair count was not verified.

Do not confuse historical quality-control summary fields with alignment counts.
The recorded `usable_reads = 39,179,382` follows the preparation code's division
of paired alignment records by two. The recorded `fragment_count = 19,581,308`
comes from a helper requiring both first-mate status and positive template
length, excluding the opposite template orientation. It must not define the
depth denominator for this research. Neither historical files nor production
code were changed.

## Frozen artifacts and validation

Research script: `benchmarks/scripts/audit_rep1_depth_inputs.py`.
Ignored artifacts: `benchmarks/results/rep1_depth_input_audit_20260922_v1/`.

- Configuration SHA-256:
  `b482720f9db9743b48edce88a15fcb7972ce8e82c33f90444f4b248f8348c280`.
- Result SHA-256:
  `38306f2cef853ffa9d9eb4695a48dbb44a1154b4ae18b6bf78452c18dceef42d`.
- Full alignment SHA-256:
  `a04cb6abc7bd06f00f0c813e31ebba774e24481d2d7075cba036e6e2638afd66`.
- Alignment index SHA-256:
  `7e0e25282ead41ce92aa51292549fa79ebbc68e2f584d63772335a3c297d3d4a`.

The audit command completed with exit 0. Its checks cover both entire alignment
files and index totals, source stability, coordinate sorting, dictionary parity,
and absence of flagged excluded alignment classes. Script syntax validation and
comparison against archived full-file statistics passed. No scientific scoring
or production interfaces changed, and no full package regression suite was run.

## Next experiment

Freeze extraction at the same 4,200 development-control sites on chromosomes
1–15, reusing only the existing sequence propensities. Retain all exclusions and
task denominators, fixed region geometry and unchanged score thresholds. Treat
this as higher-count development calibration under the recorded filtering,
not a causal depth comparison or independent safety qualification. A depth-only
comparison requires matched preprocessing or a separately frozen subsampling
design from one common source. Preserve the previous sparse-count evidence.
Do not inspect new holdout labels or final naked-DNA replicate 3.
