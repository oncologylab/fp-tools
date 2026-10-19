# Calibration reporting correction

This correction changes presentation, not the archived calibration experiment.
The rejected candidate and every historical numerical table remain unchanged.

- Failure applies to the tested models, counts, geometry and cohorts. It does
  not prove that the observations rather than any possible detector are the
  sole remaining limitation.
- ChIP-negative sites are not known unbound sites. Their 11.8% call rate is
  a diagnostic label-negative call rate, not demonstrated false-positive risk.
- Opposite-sign motif errors show that the tested map failed, not that no
  monotone mapping could ever work.
- CTCF 0.69/0.67 are positive/negative **observed-to-expected central-cut ratios**,
  not raw cleavage ratios. They do not establish rescue in difficult families.
- Genome-build compatibility and model transfer are distinct. A genome-build
  change alone does not establish a need to retrain sequence-bias coefficients.

ROC and precision–recall curves now average within held chromosome, then
biological library, then available task, matching the headline scalar metrics.
No optimistic precision envelope is used. Average precision is calculated
separately from the original predictions; area under a displayed mean curve is
not its estimate. All 21 tasks remain visible; unavailable tasks are not zeros.

Synthetic-power plots use realized central depletion from normalized whole-
fragment resampling. The perturbation parameter r is retained separately. It is
not the marginal fraction of central cleavage removed. Undefined effects at
zero-source-center windows remain in the exported ledger and are excluded from
both axes' means, not silently assigned a depletion.

Regenerated outputs: `benchmarks/results/calibration_reporting_20261019_v2/`.
Original source: calibration round at 690978f; no new fitting or scoring.
Four focused plotting tests pass. All 19 copied historical tables match their
source SHA-256 hashes. Primary equal-task AUROC/AP remain 0.568187/0.630054;
DWM-propensity same-detector comparator remains 0.530931/0.605755. This is not a
direct end-to-end TOBIAS benchmark and is not new scientific qualification.
