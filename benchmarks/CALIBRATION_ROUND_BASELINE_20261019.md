# Calibration follow-up: baseline before new work — 19 October 2026

Branch `research/footprint-calibration-20261019`, worktree
`fp-tools-calibration`, created from `60ff647`
(`research/footprint-next-round-20260922`). That tip contains the original
research history: `0569e2f` is an ancestor, and the eight later commits are the
bounded next-round investigation. Both earlier checkouts were clean at their
recorded commits; main was clean at `1d8631e`. No research process was running.
The host had 32 cores, about 93 GiB available memory and 867 GB free disk.

The next-round archive
`fp_tools_next_round_evidence_2026-09-23_97ec5a4.zip` (522,301,298 bytes)
re-hashes to the recorded SHA-256
`236aa5e744611a81d263cd90d52d30ef4ed80ae4d9315c7ef01076dad94cf0b6`.
Its unpacked source folder `next_round_20260922_v1/` is the governing
derived-data source for this round. `FP_TOOLS_NEXT_ROUND_CODEX.md` was not found
in the three checkouts; the packaged `NEXT_DEEP_RESEARCH_PROMPT.md`,
`EXPERIMENT_PLAN.freeze.json` and `LIMITATIONS_AND_BLOCKERS.md` were read instead.

## Already complete — do not repeat

| Work | Governing evidence | Status |
| --- | --- | --- |
| Strand-extraction correction and independent oracle | `STRAND_EXTRACTION_AUDIT_20260922.md` | Done for the corrected extractor; historical 96-run, ChromBPNet and transfer results remain unrequalified |
| Strand-conditioned count-tail detector (LOG21, uniform) | `CONDITIONAL_PROTECTION_20260922.md` | Retrospective nonnested replay: +0.04252 AUROC, +5.82% relative AP versus the fp-tools DWM residual comparator |
| Attainable-significance ceiling | `CONDITIONAL_SCORE_CEILING_20260922.md` | 3,589/4,053 (88.6%) positive-label sites cannot reach 1.3010299957 even with zero central cuts |
| Full replicate-1 region calibration and count screen | `FULL_REP1_REGION_CALIBRATION_20260922.md`, `FULL_REP1_COUNT_SCREEN_20260922.md` | 7/21 pooled central counts outside the independent-null interval; 65/3,967 control sites can reach threshold |
| Signed binomial and beta-binomial shrinkage | `NEXT_ROUND_20260922.md` | Rejected by the frozen screen (−0.0371 and −0.0431 AUROC versus LOG21 count-tail at full depth) |
| Complete-fragment inventories, 10M/full arrays for six chromatin libraries and naked-DNA replicate 1 | next-round `provenance/`, `arrays/` | Done; reuse, do not re-extract |
| LOG21 versus DWM propensities under one binomial detector | next-round `full_screen_summary.tsv` | +0.022675 AUROC, +2.70% AP (development factorial only) |
| Sequence-inference component timing | `DIRECT_BIAS_INFERENCE_20260908.md` | LOG21 ≈0.47× DWM CPU for that component only; not whole-workflow |

## Unresolved, and which evidence governs it

1. **Coverage confounding of the ranking estimand.** After matching on archived
   replicate-1 10M covariates, actual full-depth coverage is imbalanced in
   60/60 available task/library panels, and a fitting-fold-only covariate
   model reaches mean AUROC 0.6228. A count-tail significance grows with
   counts, so its AUROC mixes depletion with accessibility. The shrinkage
   losses therefore do not show which detector better finds depletion.
   Governing data: next-round `per_site_scores.tsv.gz`,
   `counts_and_expectations.tsv.gz`, `pool_library_covariates.tsv.gz`.
2. **Null mean calibration.** 7/21 LOG21 central-mean intervals on examined
   naked-DNA replicate 1 exclude zero; DWM propensities were not calibrated
   the same way. Variance has been modeled (beta-binomial) without first
   correcting the mean. Governing data: `control_counts.full.tsv.gz`.
3. **Information.** Most positives cannot reach the threshold; zinc-finger
   sites have expected central counts near 0.05. This is an attainable bound,
   not power under a biological alternative.
4. **Comparator qualification.** Native full-depth DWM residual/strand output
   is unqualified; the archived DWM residual exists only for the nonnested
   replay. This round cannot qualify the native production comparator without
   new full-depth production DWM runs and remains explicit about it.
5. **Biological protection.** A detector-only change cannot alter the
   identical-site aggregate profile. Whether canonical protection exists in
   difficult families is a property of the observations and expectations.

## Inputs and access boundaries

Permitted and protected inputs are listed in
`benchmarks/scripts/calibration_round/inputs.json`. This round reads only the
next-round derived package tables, arrays and frozen LOG21 file, verified by
SHA-256 before use. No BAM, FASTA or bigWig is opened unless a later frozen
plan says so. SK-N-SH ENCSR587TRP and GM23338 ENCSR485TLP labels, naked-DNA
replicate 3, and nutrient-stress data are denied by identity and resolved path.
Replicate 2 naked DNA is examined data and is not used here.

All 21 tasks stay in every ledger, including unavailable HepG2 HMGA1.
Historical scorers, the cap of 300, rounding to ten decimals and the strict
threshold 1.3010299957 are unchanged; new rules get new names.
