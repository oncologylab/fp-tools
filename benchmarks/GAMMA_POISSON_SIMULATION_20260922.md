# Gamma–Poisson low-count background: synthetic development check

This research-only model addresses undefined maximum-likelihood expectations
when a strand has no outer-flank reads. It has not established better footprint
detection, biological safety, or calibrated functional uncertainty. Main and
the production dinucleotide weight matrix (DWM) default remain unchanged.

## Frozen model and scope

For each strand, let N be the count in offsets with absolute distance 50–100
base pairs, B the summed bias exposure there, and b the center exposure.
The Poisson Jeffreys prior gives a Gamma posterior with shape N + 1/2 and
rate B. The predicted center count has mean `(N + 1/2) * b/B` and variance
`mean + (N + 1/2) * (b/B)^2`. This is the Gamma–Poisson negative-binomial
parameterization documented in the
[Stan reference](https://mc-stan.org/docs/2_23/functions-reference/negative-binomial-distribution.html);
the prior and spatial assumptions are choices of this research prototype.

Strand means and variances add under a working independence assumption.
Their sum is not generally negative binomial when their exposure ratios
differ. Scoring uses a moment-standardized residual and the unchanged
conventional center/shoulder geometry. Missing sequence remains unavailable;
zero-read 201-position neighborhoods produce score ties, not informative
calls. Background counts and prior fractions remain explicit. A finite
prior-supported expectation is not measured signal.

The protocol `manifests/gamma_poisson_background_v1.json` was frozen before
simulation. Seed 2026 generated 50,000 trials for each of 50 cases: five
background means, five exposure ratios, and two target-count generators.
The background is Poisson in both generators. The stress generator makes
center counts negative binomial with variance `mean + mean^2`.
No real data or labels were read, no prior was retuned, and no model was
selected using these results. This stage is ineligible for promotion.

## Results and limitations

All 50 cases are retained in `coverage.tsv` under
`results/gamma_poisson_simulation_20260922_v1/`.

| Generator | Observed coverage of nominal 95% intervals |
| --- | --- |
| Poisson working model | 94.964–100% |
| Overdispersed center stress | 15.018–100% |

The minimum Poisson coverage has a pointwise 95% Wilson interval of
94.769–95.152%. Discrete low-count intervals can be conservative. In the
strongest stress case (background mean 100, exposure ratio 10), coverage
was 15.018%, with interval 14.708–15.334%. This failure is retained, not
repaired by changing the prior after evaluation.

At true background mean zero, the predicted mean remains positive because
of the prior; the prior-only fraction is one. For example, exposure ratio
10 gives predicted mean 5 despite true target mean zero. Such predictions
must not be described as evidence of binding or measured accessibility.

Scalar count coverage is neither a simultaneous functional confidence band
nor motif-call false-positive control. These simulations do not measure
spatial correlation, strand dependence, bias-model error, or biological
heterogeneity. No accuracy, safety, or computational gate transfers from
previous operators to this one.

## Verification and artifact identities

Nine focused tests passed in 1.56 seconds. They cover independent numerical
Gamma–Poisson integration, predictive probability mass, scalar posterior
calculations, chunk parity, scale invariance, missing sequence, zero-count
ties, invalid inputs, and deterministic simulation stress behavior.
The separate `--verify-only` invocation exited zero and verified frozen
configuration, code, runtime, and output identities. Completed stages were
not rerun.

SHA-256 identities:

- Configuration: `9380697b3a5ffd360beae75ff11d5bbabe30fd4c54c3b54fdd3a739459c20d6d`
- Result manifest: `9dd0ded2a98bfb7778513bbafcb5070be01e0af2a496e0c16b5d8004d547f2b1`
- Coverage table: `967e7fa3f8bc1dbd9fafaa003523a75607f0b838b2ae6d4d33d9d764d70fea85`

The next permitted evaluation needs a separate configuration/code/input
freeze and training chromosomes 1–15 only. It must retain identical sites,
complete matched pairs and the full 21-task denominator, separating
expectation construction, residual formulation and detector geometry.
Untouched transfer labels and final naked-DNA controls remain unopened.
