# Historical blacklist intervention — 22 September 2026

Adding the verified hg38 blacklist resolves the two historical coverage
discrepancies in the preceding region-layout diagnostic. This is an intervention
on already examined naked-DNA replicate 2, not independent candidate safety
evidence. No candidate was scored and no thresholds changed.

## Frozen comparison

`scripts/audit_control_blacklist.py` copies the historical-region command,
adds only the blacklist and changes the output destination. All scientific
inputs, the saved dinucleotide weight matrix (DWM) model, disabled global
scaling and computational settings remain fixed. Model count and background
count arrays are checked for exact identity after execution.

Artifacts remain under `benchmarks/results/control_blacklist_20260922_v1/`:

| Artifact | SHA-256 |
| --- | --- |
| Configuration | `c775afca04f848f7fc22b8f07feeee8c4a493f51906fc0ba443df7213b240fe6` |
| Review manifest | `bad4ad902eba36bf9251e60fe54e1b08c6e30a42752195bad028716454a901f9` |
| Per-site table | `7a6d174f9798133124649b6915d4c49e9d015fdd10414a5088131744c210c70b` |
| hg38 blacklist | `31c69342df43bbc19dd8ef2886611a8150cb53e90f341de14d30e742c9251737` |

The configuration includes input, implementation and runtime identities.
The durable supervisor specification, launch and exit receipts are under
`control_blacklist_jobs_20260922/execute_attempt1/`. Execution exited zero in
24.035 seconds. This is diagnostic elapsed time, not a complete-workflow
central processing unit (CPU) or peak-memory qualification.

## Results

All 4,200 original rows remain. Both archive and reconstruction have 269 positive
expected profiles. There are zero expected-support changes, zero normalized
expected-shape failures and zero observed-shape failures. Comparisons use
relative tolerance 1e-5 and absolute tolerance 1e-7 after total normalization;
this is shape parity, not equality of globally scaled expected counts.

Thirty-five motifs overlap the blacklist; the same 35 profile windows overlap.
These are explicit annotation flags, not biological negative calls. This
includes HepG2 CTCF chr15:20624328–20624343 and K562 CTCF
chr4:11808–11823, the two discrepancies in the no-blacklist construction.
No rows were dropped or silently reassigned to the negative class.

This establishes that the blacklist is sufficient to recover archived support
and shape under the tested historical layout. It does not establish the full
original command provenance. It also does not eliminate the three previously
measured shape changes induced by expanding output regions.

## Validation and next step

- Eleven focused tests passed in 2.11 seconds, covering command isolation,
  half-open interval boundaries, separate motif/window flags and existing
  layout, reference, scope and stable-screen helpers.
- Separate frozen-artifact verification passed.
- Independent direct bigWig interval reconstruction checked all 8,400 observed
  and expected profiles against archived support and normalized shape.
- Independent dataframe interval intersections checked all 4,200 motif and
  window exclusion-flag pairs.

Next freeze a blacklist-aware all-site counterpart, keeping all previous
constructions and failures intact. Compare exclusion and region-layout effects
before deciding whether a new reference can support an explicitly revised
diagnostic protocol. Do not silently replace the failed overlap gate. Low-count
normalization diagnosis on training/replicate-1 controls and complete-workflow
resource qualification remain outstanding. Untouched promotion labels and
final naked-DNA replicate 3 remain unopened; no promotion gate is newly passed.
