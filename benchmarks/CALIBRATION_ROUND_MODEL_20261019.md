# What the evaluated method is — and is not

Research documentation for the October 2026 calibration follow-up. None of
this is a production option; `atac-correct`, `call-footprints` and the
released DWM default are unchanged.

## Structure

For one motif window and one strand `s` (motif orientation):

1. **Observations.** Strand-specific Tn5 cut endpoints, two per complete
   fragment, counted in three fixed regions: centre −7..+7 (15 bp) and
   shoulders −40..−12 and +12..+40 (58 bp). `n_s` = cuts in the three
   regions; `c_s` = central cuts. Everything else in the 201-bp window is
   ignored by the detector.
2. **Sequence bias (LOG21, frozen).** A log-linear model reads the 21-bp
   sequence context around every candidate cut position and returns a
   relative cleavage propensity at one-base spacing. Its 21-bp *context* is
   not an output *resolution*: it predicts one value per base. The null
   central share `p_s` is the propensity mass in the centre divided by the
   mass in centre plus shoulders.
3. **Null-mean recalibration (evaluated candidate).**
   `q_s = expit(a + b·logit p_s)`; the ablation fixes `b = 1`. `a` and `b`
   are fitted by maximum likelihood on protein-free naked-DNA windows only.
4. **Detector (unchanged).** Conditional on `n_+` and `n_−`, the central
   count is a sum of two binomials. Score =
   `−log10 P(Bin(n_+,q_+) + Bin(n_−,q_−) ≤ c_+ + c_−)`, capped at 300 and
   rounded to ten decimals.
5. **Information and calls.** Ceiling = score at zero central cuts. A site
   is `CALL` if score > 1.3010299957 (strict), `NO_CALL` if the ceiling
   exceeds the threshold but the score does not, `ABSTAIN_UNREACHABLE` if
   even zero central cuts could not cross the threshold, and
   `NO_INFORMATION` if `n = 0`. Abstention is not a negative call.

Effect, uncertainty and calling remain distinct: the observed/expected
central ratio is the effect estimate, the ceiling states how much evidence
the counts can carry, and the call is a nominal tail test. The score is a
ranking statistic, not a binding probability.

## Training data

| Component | Source | Labels used |
| --- | --- | --- |
| LOG21 | Five-seed ensemble on chromatin-library mitochondrial reads (K562/HepG2, shifts +4/−5, L2 0.0001); 7,519-byte NPZ | No |
| Recalibration `a, b` | Naked-DNA replicate 1 (GSE164997 SRR13450940), full depth, 5,286 strand observations from deduplicated windows on chr1–15; each evaluation fold leaves its chromosome out | No |

Deployable record (all chr1–15): `a = 0.1171, b = 0.9381` (primary);
`a = 0.2184, b = 1` (ablation). These were never used for evaluation.

## What a user would need

A protein-free control library from a comparable Tn5 protocol, the frozen
LOG21 file and the two scalars. Retraining would be required for a different
enzyme lot or protocol, a different genome build, or a library with a
different fragment-end convention (the +4/−5 shift). No automatic retraining
rule is claimed.

## Failure conditions found in this round

- A global intercept fixes the pooled naked-DNA mean but not motif-specific
  errors: 6/21 task ratios still exclude 1 after recalibration, in opposite
  directions (CTCF and MYC over-predicted; K562 FOXA1, MEF2D and TBP
  under-predicted).
- The naked-DNA null does not transfer to chromatin motif sites. At
  ChIP-negative sites the recalibrated observed/expected central ratio is
  0.34–0.99 by task, and nominal full-depth calls occur at 11.8% of
  ChIP-negative site/library rows versus 18.3% of ChIP-positive rows.
- Information is scarce: under the candidate, 37% of ChIP-positive and 55% of
  ChIP-negative site/library rows abstain at full depth, and over 90% at 10M
  fragments. A whole-fragment
  synthetic removal of every central-endpoint fragment is called in only
  about 10% of resamples.
