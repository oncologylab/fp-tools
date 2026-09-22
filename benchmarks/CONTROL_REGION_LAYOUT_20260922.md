# Historical control-region isolation — 22 September 2026

This is a development reference diagnostic on already examined naked-DNA
replicate 2, not untouched safety or promotion evidence. No candidate was scored.
The production dinucleotide weight matrix (DWM) default remains unchanged.

## Paired construction

`isolate_control_region_layout.py` reuses the frozen expanded-reference command,
changing only its output-region file to the historical peak file and its output
destination/prefix. Bias-model count arrays remain identical; global scaling
is disabled in both constructions. Completed artifacts are preserved under
`benchmarks/results/control_region_layout_20260922_v1/`.

- Configuration SHA-256:
  `7913ce31f44d55fcf1ccc282b0bc162d5e137e2d5127c1db0fcde770601f0f8e`
- Review manifest SHA-256:
  `d1a6378c67d692ac634b5e5e742997ff047f828566039a2cb7528e937d105fd3`

All 4,200 original site rows remain. All 269 archived positive expected
profiles reproduce their normalized shape. Expanded output regions change
three shapes: HepG2 CTCF chr10:28991925–28991940, K562 CTCF
chr9:98416493–98416508 and K562 FOXA1 chr9:98416381–98416389.
The intervention isolates output-region layout; correction-window phase and
other boundary effects have not been separately isolated.

## Remaining coverage provenance

Two sites have observed and expected signal in the new historical-region run
but were empty in the archive. Both overlap the existing hg38 blacklist:

| Site | Blacklist interval |
| --- | --- |
| HepG2 CTCF chr15:20624328–20624343 | chr15:20005600–22606300 |
| K562 CTCF chr4:11808–11823 | chr4:0–69200 |

The blacklist SHA-256 is
`31c69342df43bbc19dd8ef2886611a8150cb53e90f341de14d30e742c9251737`,
also recorded in `manifests/direct_bias_replicate_depths_v1.json`.
The new constructions omit a blacklist. Overlap is consistent with a historical
exclusion, but does not establish which original command was used. A separately
frozen matched blacklist-aware rerun is required; preserve exclusion flags and
the original denominator. Do not silently reinterpret excluded windows as
biological negatives or waive the expanded reference's failed qualification.

## Verification and next work

Nine focused tests passed in 2.00 seconds across region-layout, all-site control,
track-scope and stable-screen helpers. The separate `--stage verify` passed.
An independent read of all 4,200 table rows reproduced the 269/0/3 positive,
historical-shape-failure and expanded-shape-failure counts. Direct interval
intersection confirmed both blacklist overlaps. These checks are not an
independent reconstruction of every numerical profile.

Reference qualification, low-count normalization diagnosis on training and
replicate-1 controls, and complete-workflow central processing unit (CPU)/memory
qualification remain outstanding. Do not change failed thresholds based on
replicate-2 results. Final naked-DNA replicate 3 and new promotion labels stay
unopened. No promotion or main-branch integration is justified.
