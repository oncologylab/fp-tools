# Nuclear sequence-bias follow-up: failed chromatin screen

The frozen primary candidate is rejected for footprint improvement. This is a
completed bounded experiment with explicitly incomplete promotion qualification,
not a claim that no footprint detector could ever improve these observations.
Main, released scientific defaults, manuscript and historical artifacts were
not changed. No protected labels or naked-DNA replicate 2/3 were opened.

## Execution

The primary freeze preceded training; conditional chromatin code and input
identities were frozen at `ff23cdf` before scoring. The run used 39,179,095
complete naked-DNA replicate-1 fragments; 78,358,190 independent coordinate
checks had zero mismatches. All 150 nuclear LOG21 fits and 30 ensembles finished.
The primary regularization remained 0.01; 0.0001 stayed an ablation.

The control component passed: 2/21 task mean discrepancies, and positive paired
held-chromosome likelihood intervals against both frozen mitochondrial LOG21
and intercept-only calibration. The unchanged 6,020-site training cohort was
then evaluated across all six biological libraries at full and 10M depth, seed
2026. All 21 tasks remain in the ledger; HepG2 HMGA1 remains unavailable.

## Scientific result

| Comparison with same-detector DWM propensities | AUROC gain | Relative average-precision gain | Paired AUROC interval |
|---|---:|---:|---|
| Full-depth primary nuclear model | 0.012710 | 1.130% | [-0.006426, 0.030941] |
| 10M primary nuclear model | -0.001267 | -0.179% | [-0.011679, 0.009284] |

AUROC means area under the receiver operating characteristic curve. Average
precision retains the original implementation and equal chromosome, library,
then task weighting; it is not the area under the displayed mean PR curve.
Full-depth accuracy requirements of 0.03 AUROC and 10% relative average
precision both fail; the principal interval includes zero. The primary model
also loses against frozen mitochondrial and intercept-only LOG21. Its mean
within-count-stratum AUROC is 0.481796, versus 0.504191 for DWM propensities.
There is no supported general footprint rescue.

The control result does not establish chromatin-null calibration. Source,
preparation, depth and window selection differ from mitochondrial training,
preventing attribution uniquely to nuclear versus mitochondrial DNA. The
synthetic fragment challenge is descriptive: among 68 defined source windows
with complete synthetic central depletion, primary mean call sensitivity is
13.95%. This is not population power, a biological binding simulation or a
sequencing recommendation; zero-source-center windows have undefined depletion.

Full-depth residual sign diagnostics reproduce canonical signs in all three
libraries for CTCF in both cells and for HepG2 ZNF362. These are model-dependent
shape diagnostics, not proof of a 25% difficult-family rescue. Raw observed cuts
are unchanged. All-task raw/expected/residual figures and source tables are
provided, including failures and unavailable tasks. Additional count/expected-
count-stratified shape panels and functional-band calibration remain incomplete.

Native DWM strand and complete-pair qualification remains incomplete, despite
processed-region support for all 4,200 controls. The comparison above is not an
end-to-end native TOBIAS benchmark. The expensive depth/seed and end-to-end
resource matrices were intentionally stopped after the failed accuracy screen.
Untouched safety and functional uncertainty remain unqualified.

## Evidence and checks

Archive: `benchmarks/results/fp_tools_next_round_evidence_2026-10-19_1473f2d.zip`

Size: 374,477,584 bytes. SHA-256:
`8953b12dd95b4d90aadf4e6b525fc6afdab9da02dbc4737f0f930bbd7fe1cc7f`.

The archive includes complete derived primary arrays, all-site counts/scores,
all-task tables, model coefficients, freezes, provenance, control simulations,
ROC/precision-recall curves, all-task profile figures, corrected historical
slides, an independent verifier and explicit reanalysis entry points. No raw
BAM/genome payloads or protected contents are included. A preliminary package
was retained locally after finding an absolute workspace path in its copied
historical receipt; the delivered package sanitizes it and records both hashes.

Fourteen focused nuclear tests passed. Independent calculations on the original
outputs reproduce 12,040 sequence-prediction sites, 72,240 new exact tails/calls
and 8,820 defined rank/AP rows. Twelve paired summaries reproduce using the
frozen aggregation contract; that replay is not an independent implementation
of matching or full-refit uncertainty. Fresh-extraction verification receipts
are under `benchmarks/results/nuclear_evidence_fresh/`.

Fresh-extraction checks completed successfully: all 773 payload checksums,
8,400 control tails/calls, 1,380 power rows, 15 fitting-fold exclusions,
12,040 chromatin sequence predictions, 72,240 new chromatin tails/calls and
8,820 rank/AP rows. Explicit aggregation replay matched all 12 screen rows,
252 paired-task rows and 294 task-metric rows. The representative raw figure
reproduced 402 profile values and was visually inspected. All original code
hashes in the main, diagnostic and chromatin freezes resolve to preserved
packaged source. The privacy scan found no personal absolute paths in the
delivered text/code payloads. These checks establish reproducibility within
the documented scope, not scientific promotion.

From the extracted top-level folder:

```bash
python -B code/verify.py .
python -B code/verify_chromatin.py . --output ../chromatin_replay
python -B code/reanalyse.py . --outdir ../control_replay
python -B code/replot_chromatin.py . --output ../representative_chromatin.png
```

Stop this candidate. No promotion, integration, upload or publication is
authorized by these results. Any future round needs a distinct falsifiable
question rather than new detector tuning on the same examined labels.
