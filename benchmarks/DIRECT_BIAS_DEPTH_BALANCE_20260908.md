# Direct-bias sensitivity to current-depth cut balance

CTCF retains a full-depth advantage over both DWM comparators after the
current-cut balance filter. The comparison with raw signal and the lower-depth
checks are less consistent. CTCF remains the only powered task that passes the
canonical-protection screen. This completed secondary analysis provides no
general replacement or promotion evidence.

## Fixed protocol and evidence scope

Calculation commit: `e8986ed`.
Plan: `manifests/direct_bias_depth_balance_v1.json`.
Output directory: `results/direct_bias_depth_balance_20260908_v1/`.

The [primary 36-run study](DIRECT_BIAS_STABILITY_20260908.md) matched motif
score, reference accessibility and reference cut totals. Current cut-total
SMD exceeded 0.1 in 11/36 CTCF runs. This secondary diagnostic keeps each
original pair only when its current log1p cut totals differ by at most 0.2
training standard deviations. The per-cell/TF scales already fitted on
chr1–15 are reused. There is one caliper, no new partners, no model selection,
and no adjustment after inspecting the secondary results. All 36 pair masks
were written before new scoring; zero-cut pairs can remain when the same
symmetric rule passes.

The source is the already examined chr19–22/X development archive: five
nested seeds at 10M/25M/50M for replicate 1, and one full-depth endpoint for
each of three biological replicates in HepG2 and K562. Counts and expectations
are reused without BAM extraction or model fitting. The mitochondrial LOG21
ensemble, strand-total-preserving expectation, signed-deviance operator,
conventional geometry, and dispersions (HepG2 0; K562 1.6144318781609228) are
unchanged. Every within-run comparison uses identical sites with positive
fixed and sample DWM expectations. Primary metrics and files are preserved.

This is post-hoc conditioning on ATAC counts. It changes the retained sites
and can change task difficulty or introduce selection effects. It does not
identify a causal correction effect or turn the development data into an
untouched validation set.

## Full-depth CTCF

The filter retains 44–54% of the original full-depth CTCF pairs. All measured
covariate SMDs are below 0.04 in every CTCF run; current cut-total SMD now
ranges from −0.0244 to 0.0369, compared with a previous maximum of 0.1344.

| Sample | Retained pairs | Cut-total SMD, original → filtered | Parametric AUROC / AUPRC | AUROC gain over raw, 95% interval |
| --- | ---: | --- | --- | --- |
| HepG2 rep1 | 293 / 577 | 0.0516 → 0.0118 | 0.6984 / 0.7310 | 0.0235 [−0.0128, 0.0434] |
| HepG2 rep2 | 256 / 577 | 0.0259 → −0.0035 | 0.7402 / 0.7648 | 0.0728 [0.0368, 0.0934] |
| HepG2 rep3 | 259 / 577 | 0.1081 → 0.0205 | 0.7214 / 0.7415 | 0.0654 [0.0510, 0.0808] |
| K562 rep1 | 289 / 534 | 0.1167 → 0.0146 | 0.7672 / 0.7536 | 0.0627 [0.0149, 0.0833] |
| K562 rep2 | 260 / 534 | 0.1344 → 0.0078 | 0.7657 / 0.7583 | 0.0651 [0.0194, 0.0956] |
| K562 rep3 | 266 / 534 | 0.1148 → 0.0059 | 0.7744 / 0.7690 | 0.0703 [0.0109, 0.1048] |

All six full-depth runs retain positive AUROC and relative AUPRC intervals
against both fixed and sample DWM. Gains over sample DWM are 0.0845–0.1636
AUROC and 21.2–29.3% relative AUPRC. All six pass the canonical-protection
screen. Standardized aggregate separation is 60.6–80.4% stronger than sample
DWM within the filtered cohorts; the plots show the central contrast and
shoulders on shared vertical axes. These percentages must not be interpreted
as an additional improvement over the primary study, which uses other sites.

HepG2 replicate 1's AUROC interval against raw includes zero. Its AUPRC
interval remains positive. The other five full-depth runs have positive
intervals for both metrics against raw. The change in HepG2 replicate 1 can
reflect both the selected site population and the smaller sample; it cannot
be attributed solely to removal of confounding.

## Depth and power limitations

The filtered lower-depth results do not support an all-run robustness claim:

- HepG2 10M seed 2027 retains 198 pairs, below the 200-per-class requirement.
- AUROC point gains over raw are negative in HepG2 10M seed 2028 (−0.0051)
  and 25M seed 2030 (−0.0026).
- Relative AUPRC point gains over raw are negative in HepG2 10M seed 2030
  (−1.62%) and K562 10M seed 2030 (−0.47%).
- The pointwise central-protection screen fails in HepG2 10M seed 2028,
  25M seeds 2028/2029, 50M seed 2026, and K562 10M seed 2027. Mean centers
  remain negative and mean depletion remains positive; these five failures
  concern insufficient pointwise exclusion of zero, not proof of absent
  footprint information.

The screen passes in 31/36 CTCF runs, including the underpowered 198-pair
run. K562's 10M seed-2026 metric interval against raw still includes zero.
Metric intervals use 1,000 chromosome-block draws for all full-depth endpoints
and the two 10M seed-2026 references; other runs have point estimates.
Aggregate bands use 500 draws. Only five chromosome blocks are available;
these unadjusted pointwise diagnostics do not certify simultaneous functional
coverage or calibrated binding probabilities.

## Other TFs and the task denominator

On filtered full-depth sites, HepG2 ZNF362 remains near AUROC 0.604–0.608,
HepG2 ZNF384 near 0.533–0.569, and K562 ZNF384 near 0.573–0.585. Some
numerical gains over raw become positive, but none passes canonical protection.
These are not new footprint improvements or difficult-family rescues.

HepG2 FOXA1 retains only 136–187 pairs and K562 MYC only 95–106 at full depth;
both become underpowered. Dropping them changes the task set. Means over only
the remaining powered task/replicate rows give relative AUPRC gains of
11.24%/16.10% over sample DWM for HepG2/K562. On the original eligible task
set, the descriptive means are 7.97%/10.51%, including the now-underpowered
rows. Neither calculation establishes the locked mean-gain gate. A changed
subset, lost power, or a smaller set of TFs cannot replace the original
promotion study.

All 21 source cell/TF tasks remain in the 378-row retention table. It includes
288 underpowered task/run rows and 16 zero-pair rows. No new difficult family
has acquired powered canonical protection.

## Integrity, artifacts and disposition

The frozen source graph was checked before selection and again after
completion by the figure command's full historical-code verifier. All 36
nested evaluations passed. The combined tables contain 2,876 method/task/run
metric and shape rows and 150 bootstrap comparisons, with primary and
secondary roles explicitly labeled. The all-TF PDF has 21 A4 pages; every
page was visually reviewed, including non-estimable and underpowered tasks.
The renderer also creates HepG2 and K562 CTCF PNG previews.

- Configuration ID: `6337d9417e839966901c65878dc7efbf3c1f532f242f12d65a47dc92bd3235a1`.
- Configuration SHA-256: `c7efa722036c5fe5a2d4377330c9f7864da982456dd40a7d657671622da579c9`.
- Joint selection ID: `37eed82dbe58847cd5423c20ae48af2b2805953c50a0c231c2e6baf409966e05`.
- Joint selection SHA-256: `e10689af7c4f5c8f4832e2aee737de3333d92eb74c737a3ace081e6f33b88d56`.
- Result ID: `94e2814a9b32cbb45355e8898a3c6a5c211366c920ddba4cfc1e2f0b7c4a4e8e`.
- Result SHA-256: `99974a0c45b0cadabb62ddcdf23ed752bd54e7d2f35d0846d6f5393e93546160`.
- Figure ID: `aa2d5d0e5a8b72fe8ec405aef67ef8570694d4e438e9fa530d086bd42c66b7a5`.
- PDF SHA-256: `b62f36eda0020d370c969dce12af989ed42deca851efe023a9bfbc8437b37911`.

The 72 focused calculation/stability/documentation checks and the renderer
fixture test passed. The strict documentation build and 32-page browser audit
at three widths also passed. The original 36-run results and the separate
96-run queue remain unchanged. This diagnostic does not repeat naked-DNA
safety, establish absence of occupancy leakage, or measure complete-workflow
resources. SK-N-SH/GM23338 labels and final naked-DNA replicate 3 remain
unopened. Production DWM and the frozen manuscript are unchanged.
