"""Write narrowly scoped reports from completed, inspectable result tables."""
import argparse
from pathlib import Path
import json
from verify import read_table


def execute(out):
    def write(name,text):(out/name).write_text(text.strip()+'\n')
    baseline=(out/'BASELINE_STATE.md').read_text().split('\n## This executed round')[0]
    write('BASELINE_STATE.md',baseline+'''

## This executed round

The original source at 0569e2f governs archived extraction/scoring context; it is not a new validation of historical ChromBPNet or transfer results.
The initial plan was frozen in the new research worktree before new scores. Source commits 8791f7b and 4a274f0 record the signed model and algebraic near-zero-dispersion repair; 32764aa records the extraction/cohort/initial-panel execution implementation.
Original replay arrays and old thresholds remain unchanged. Archived control fits and model fixtures were checked against the numerical repair's affected domain; archived fits were above that domain and their refit differences were below 5.2e-11. They were retained with their original source identity.
The first new 10M calibration failed numerically; the first full calibration was stopped, both attempts retained. Both affected calibrations then completed from the repaired likelihood, and all 12 new full/10M scoring endpoints used the repaired code.
The first support audit treated sparse zeros as missing; the second mixed normalized/raw observation units. Both remain visible. Version 3 corrects support and units; an additional fixed-layout replay isolates the two remaining expected-shape differences.
The final metadata audit separated the historical nominal 10M label from verified complete-fragment depth. All 245,042 nuisance-covariate rows, count features and fitted depth constants were unchanged; balance-table bytes were identical. The original metadata table and repair receipt are retained. This was not a scientific refit.
The current source and packaged byte identities are in ORIGIN_INDEX.json. Per-stage freezes and launch/exit receipts govern exact commands; success receipts alone do not qualify the science.
The old full-replicate-1 reports' proposed next steps are historical: this round actually executed the new full-depth panel, and its failed advancement screen stops further expansion.
''')
    screen=read_table(out/'tables/full_screen_summary.tsv')
    rows=['| Method | Comparator | Mean AUROC change | Mean relative average-precision change | Available / original tasks |',
          '|---|---|---:|---:|---:|']
    for r in screen.itertuples():
        rows.append(f'| {r.method} | {r.reference} | {r.mean_auroc_gain:+.6f} | {r.mean_relative_auprc_gain:+.2%} | {r.available_tasks}/21 |')
    summary='\n'.join(rows)
    controls=read_table(out/'tables/control_mean_variance.tsv.gz')
    center=controls[controls.depth.eq('full')&controls.arm.eq('LOG21')&controls.strand.eq('combined')&controls.region.eq('center')&controls.stratum_type.eq('all')]
    mismatched=int(((center.mean_error_low>0)|(center.mean_error_high<0)).sum())
    paired=read_table(out/'tables/paired_control_rates.tsv')
    failed=paired[paired.depth.eq('full')&paired.development_check.eq('FAIL')]
    inventory=read_table(out/'tables/fragment_inventory.tsv')
    nuisance=read_table(out/'tables/covariate_baseline_task_metrics.tsv')
    nuisance_auc=nuisance[nuisance.depth.eq('full')].auroc.mean()
    balance=read_table(out/'tables/actual_library_covariate_balance.tsv')
    current=balance[balance.depth.eq('full')&balance.phase.eq('after')&balance.covariate.eq('actual_whole201')&balance.status.eq('AVAILABLE')]
    imbalanced=int(current.standardized_mean_difference.abs().gt(.1).sum())
    write('README_FIRST.md',f'''
# Bounded development evidence: read this first

This is an executed research round and a **partial scientific qualification**, not a release or a promoted footprint detector.
Main, released defaults and the published manuscript were not changed. New holdout labels and final naked-DNA replicate 3 stayed closed.

Start with EXECUTIVE_FINDINGS.md, LIMITATIONS_AND_BLOCKERS.md and BASELINE_STATE.md.
figures/fp_tools_development_summary.pptx is an editable two-slide briefing; the adjacent PDF is its Linux-rendered preview. Slide notes identify comparators, sources and limitations.
EXPERIMENT_PLAN.freeze.json records the pre-score policy; later implementation freezes and preserved failures are in provenance/ and logs/.
ORIGINAL_PROMOTION_REQUIREMENTS.md retains the original gates. No general scientific PASS is asserted.

The complete 21-task ledger includes unavailable HepG2 HMGA1. The historical 8,106-site replay is nonnested retrospective development.
The new 6,020-site held-training-chromosome cohort uses fitting-fold-only scales and a frozen held-sampling rule.
The union of replay and new cohort is 11,064 sites; arrays retain every observation used for headline metrics.
All 35,006 original pool rows remain in pool.tsv.gz. Optional unselected pool profile arrays are omitted only as documented in ORIGIN_INDEX.json.

From a fresh extraction, with Python and code/requirements.txt installed:

```bash
python -B code/verify.py . --figure-check
python -B code/reanalyse.py . --outdir ../fresh_reanalysis --workers 4
python -B code/recalibrate.py . --outdir ../fresh_calibration --depth full --refit-dispersion --workers 4
```

The first command is read-only and never fits, repairs or regenerates missing scientific output.
The second is the separate run entry point: it recomputes frozen site scores, calls and metrics from complete derived data, with no raw alignment files.
The third separately refits control dispersion and replays the calibration thresholds on fitting chromosomes. It does not alter the packaged scientific evidence.
Use code/package.py only to create a fresh ZIP from an existing verified folder; it does not perform scientific calculations.
Dependencies and original environment identities are included. Source execution code is included separately under code/source/; raw-data extraction requires the authorized original checkout and inputs, which are not distributed.

AUROC means area under the receiver operating characteristic curve. AUPRC is computed as average precision, the packet's area-under-precision-recall implementation.
DWM means dinucleotide weight matrix. LOG21 is the frozen 21-base log-linear sequence-bias model.
TF means transcription factor; CTCF means CCCTC-binding factor; FPR means false-positive rate. GC means guanine-cytosine content; CPU means central processing unit.
The DWM residual arm is an fp-tools expected-profile/residual comparator, **not a native end-to-end TOBIAS benchmark**.

Read the data dictionary before pooling rows: the same site can appear in two cohorts, several methods, libraries or synthetic repetitions.
These rows are not independent biological observations. Model posterior quantities are conditional on a working likelihood and are not binding probabilities.
''')
    write('EXECUTIVE_FINDINGS.md',f'''
# Executive findings

**Decision: stop this candidate search and do not advance either shrinkage method to the expensive robustness/resource matrix.**
The complete initial panel ran: six biological libraries, each at 10 million complete fragments and full depth, seed 2026.
The archived original-setting replay also ran. No protected holdout/control contents were used.

The proposed beta-binomial and simpler binomial models were implemented with the frozen zero-centered Normal(0,1) prior and deterministic numerical integration.
One numerical cancellation bug near zero dispersion was repaired algebraically; affected attempts were preserved and rerun, with unchanged hypotheses, cohorts and thresholds.

## Full-depth development comparison

Equal chromosome, then equal biological library, then equal available task. Missing tasks are explicitly excluded from these diagnostic means, never treated as zero or pass.
The original denominator remains 21. A NaN against the archived DWM arm means that the requisite native full-depth comparator is unavailable, not zero gain.

{summary}

These are ranking diagnostics. A detector-only change leaves the observed, expected and residual profiles on identical sites unchanged, so it cannot establish the required 25% stronger underlying protection.
All-task principal panels use all held matched sites; no favorable score-selected subset is substituted. Some actual profiles can show canonical protection, but the new score did not create it.

**Important confounding limit:** matching was frozen on archived replicate-1 count covariates. It did not remain balanced for actual full-depth coverage: {imbalanced}/{len(current)} available task/library panels have an absolute standardized mean difference above 0.1. The fitting-fold-only motif/accessibility/depth baseline has mean full-depth AUROC {nuisance_auc:.6f} on the complete held cohort. Thus the old significance score's ranking advantage is not proof of superior biological footprint detection; count confounding may contribute. The cohort was not rematched after seeing this result. Complete actual-library pool covariates, balance and stratified results are retained.

## Mechanistic conclusions

1. Signed shrinkage and an extra dispersion parameter do not establish a general rescue. Read the complete table above rather than selecting one favorable transcription factor (TF), library or depth.
2. Information limitation is real in low-count settings, but the zero-central-count ceiling is only an attainable-significance bound. Model-based and whole-fragment synthetic experiments separately measure sensitivity and reveal dependence on expected central counts. They do not prove absence of a biological footprint or justify a sequencing recommendation.
3. Null mean calibration remains unresolved: {mismatched}/21 full-depth LOG21 task-panel central-mean discrepancies have exploratory chromosome-block intervals excluding zero. This count is not multiplicity-adjusted, and the panels reuse one examined naked-DNA library. Wider beta-binomial intervals do not repair a shifted mean.
4. {len(failed)} full-depth candidate/task development-control checks fail at least one informative-site rate, descriptive Wilson upper-bound, or same-operator paired-DWM criterion. These are examined controls and conditional on available support, not untouched safety evidence. The original native-DWM safety comparison is incomplete.
5. The remaining limitation is a combination of sparse local information, imperfect null expectations and task-dependent profile biology. This experiment does not identify a universal missing detector. Stop adding model complexity here; if a later round is authorized, investigate null mean transfer and finish native comparator/strand qualification before any new holdout use.

## What changed from the old packet

New complete-fragment extraction/inventory, nested held-chromosome matching, both signed-shrinkage candidates, chromosome-held control calibration, six-library initial depth panel, fragment-based challenge, and all-task shape/confounding evidence.
The previous full-replicate-1 region/count reports were already completed and were not presented as future work.
There are {int(inventory.complete_fragments.sum()):,} complete source fragments across the six chromatin libraries plus the single control library; this total is an accounting result, not independent safety replication.

The beta-binomial primary has not justified production complexity over the fallback. Neither is promoted. Functional uncertainty calibration, full depth/seed robustness, qualified complete-workflow resources and untouched controls remain unmet or untested.
''')
    write('LIMITATIONS_AND_BLOCKERS.md','''
# Limitations and blockers

- All evidence is development or historical replay. Training chromosomes 1–15 only; no tuning expansion, new test-label evaluation or fresh transfer confirmation.
- The replay's matching is nonnested. The new evaluation isolates fitting by whole chromosome, including pooled tasks; held labels enter only the frozen matching/sampling rule. Overlapping training folds and prior development choices still limit inference. Fragment-sharing inventories record cross-chromosome checks.
- Matching uses archived count covariates. Actual full-depth coverage balance exceeds the original 0.1 standardized-mean-difference limit in all 60 available task/library panels. This is a measured confounding limitation, not repaired by changing the evaluation cohort after scoring. Covariate-only discrimination is substantial; rank losses do not prove an intrinsically worse biological detector.
- The archived DWM residual comparator exists only for the original replay. Raw DWM sequence propensities with the new scoring operator are a separate factorial control, not a reconstruction of production correction. The all-site native control track has processed-region support, but expanded region layout changes two historical expected shapes and native strand expectations remain unqualified. No native full-depth residual comparator was invented.
- HMGA1 in HepG2 has no balanced evaluation cohort. Other underpowered tasks retain their rows and pre-score eligibility flags. Available-task diagnostic means cannot certify the original promotion denominator.
- Replicate-1 naked DNA is already examined development data. Replicate 2 was not used again; final replicate 3 stayed closed. Reusing a control library across tasks/cells does not create biological replication. Site-overlap dependence limits Wilson interpretation.
- Beta-binomial strand marginals form a working composite likelihood. A marginal dispersion parameter does not remove fragment-end dependence. Its chromatin transfer is not established. The shape residual's historical dispersion parameter is separate from the new count-model rho.
- Empirical challenge windows are deliberately stratified and few; their average null call fraction is not a population-wide false-positive rate. Resampling 2,000 times does not create 2,000 independent source libraries. Whole-fragment perturbations have no single common theta truth when strand marginals shift differently; scalar coverage is explicitly undefined there.
- Control null audits include count, expected-count, sequence/GC and native-support strata. The separate template-length audit is descriptive and does not qualify a length-dependent calibration model. Complete endpoint spans/orientations are preserved for empirical challenges; per-site recorded template lengths and cut-endpoint spans are distinct quantities.
- Model-fixture interval coverage is conditional on its generating model. Functional bands were not qualified. No site interval or posterior is advertised as a binding probability.
- The fragment-link augmentation independently recovered complete-fragment/depth accounting after the first extraction had omitted per-site fragment summaries. Control 10M window membership is exactly replayed, but the global control cutoff is bounded rather than recovered outside those windows.
- Historical ChromBPNet, transfer and 96-run stability results were not independently requalified. No direct native TOBIAS run was performed.
- Formal resource measurements did not run because the advancement screen was not met, production-equivalent candidate outputs are incomplete, and absolute simultaneous process-tree memory instrumentation was not qualified. Job elapsed times and pure model sizes do not replace these gates.
- The whole pre-matching pool table is retained, but nonselected 201-position arrays may be omitted from the compact package with source hashes. Complete arrays, counts and scores for every headline site are retained. Raw BAM/FASTA and local-path traces are intentionally absent.
- Some source/preparation scripts require original raw inputs and native extensions. Portable derived-data reanalysis and the independent verifier do not.
''')
    write('DATA_DICTIONARY.md','''
# Data dictionary and interpretation

Stable site_id identifies cell/task/motif coordinates; physical parent_locus_id identifies chromosome:center. Different tasks can share a physical window.
TFBS (transcription factor binding site) is the historical column prefix for motif locations, not an assertion of binding. Coordinates use the zero-based, half-open browser extensible data (BED) convention. center is integer motif midpoint. Array columns are -100 through +100 inclusive in motif orientation; minus motifs reverse positions and swap genomic strands.
pair_id joins positive/negative matched evaluation mates within one task/chromosome. overlap_group and shared_fragment_component identify related windows; fold is the held chromosome across every pooled task.
chip_label is permitted historical training-set chromatin immunoprecipitation occupancy labeling, not proof of a canonical footprint. Unavailable labels/tasks are not imputed.

counts_and_expectations.tsv.gz contains every site/library/depth and control record. Individual counts.* and control_counts.* tables are the immutable stage outputs used by the verifier.
pool_library_covariates.tsv.gz retains actual whole-window and strand counts for every original fitting-pool site/library/depth, including sites whose optional full profiles are omitted. These covariates permit nuisance-baseline refitting without raw alignments.
plus/minus_center: positions -7..7. left: -40..-12. right: 12..40. total: these three disjoint regions, not the entire 201-base array.
unique_fragments_* and cut_assignments_* refer to the complete 201-base window. They must not be confused with regional cut totals or mapped alignment records.
LOG21/DWM_*_p are central probabilities normalized over center+shoulders. Uniform p=15/73. DWM here is pure sequence propensity, not a native expected-track output.
Expected central count is strand total multiplied by its null probability. Validity flags distinguish unsupported expectations from zero-count sites.
depth=archived_10M retains historical preparation; it is not a causal depth comparison with new filtered sources. New 10M and full arrays share one complete-pair source and deterministic seed 2026.
For archived_10M, complete_library_fragments is unknown. depth_covariate_constant retains the historical nominal value 10000000 solely to reproduce the original nuisance calculation; it is constant within a library and centers out. It is not a verified fragment inventory. New-library complete depths are independently inventoried.
The initial 10M/full panel uses its predeclared depth-specific control dispersion fits and calling thresholds. It is not the deferred locked-apply-at-depth robustness test.

per_site_scores.tsv.gz repeats sites by cohort/method. score is frozen count-tail significance or posterior mean signed theta. Positive theta denotes depletion; negative theta denotes enrichment.
rho is the held-chromosome control-fitted dispersion. mean/median/sd/low/high/probability_positive are conditional working-model quantities. prior_only/no_information prevent calls when total counts are zero.
threshold uses each method's recorded policy; call uses strict >, not >=. The old count-tail cap of 300, rounding to 10 decimals and threshold 1.3010299957 remain unchanged.
attainable_ceiling.tsv.gz is the zero-center maximum old significance, not biological power. The fixed geometry and all original strata are retained.

calibration.tsv rates use supported/eligible development-control rows; paired_control_rates.tsv separately records all original 200 sites/task, informative eligibility, calls and paired denominators. Excluded control outcomes are not silently called safe.
Wilson intervals are descriptive under dependence. paired_control_rates uses same-operator DWM shrinkage, never the historical native-DWM policy.
model_power.tsv and model_simulations contain every synthetic count draw and interval. empirical_power.tsv/power.tsv and empirical_simulations arrays retain complete fragment challenges, block identities, realized depletion and draw-level counts/calls.

fold_metrics.tsv uses tie-aware AUROC and average precision, equal chromosome weighting. task_metrics and paired_deltas retain undefined tasks/folds; full_screen_* adds equal biological-library then equal-task weighting.
Missing/undefined values are empty fields/NaN, never zero. Matched-cohort average precision is not genome-wide precision. Bootstrap intervals condition on fixed predictions, not complete refitting.

shape_curves.tsv.xz includes complete source curves by task/library/depth, label, strand, weighting and count/expected strata, including empty strata. It is losslessly recompressed from the original gzip table. Unstratified principal rows have empty stratum fields.
Raw/expected units are cuts per site per base; residual units are signed deviance residuals. Shape metrics use fixed outside-flank detrending/scaling, not the plotted raw units.
Canonical diagnostic center is -5..5 and shoulders -40..-15,15..40. Scoring geometry is separate. Standardized functional separation is the unchanged root-mean-square (RMS) standardized label-curve difference within +/-50, with the source variance-floor rule.
principal_regional_shapes.tsv gives signed center/left/right contrasts and asymmetry for both prespecified geometries, separately by label, strand and weighting in the explicit plotted units. It does not replace or renormalize the standardized separation gate.
detector_shape_identity.tsv records the algebraic zero profile change between identical-track detector arms. It is not evidence that a new method improved measured protection.

arrays are non-pickle NumPy NPZ files; every observation/propensity array has a site_id alignment. Models contain frozen coefficients and per-fold calibration records, not raw protected inputs.
For libraries 2/3, identical LOG21 arrays are losslessly shared through a LOG21_reference pointing to the corresponding library 1 propensity file. Array equality and site order are checked before deduplication; code/verify.py resolves these local references. No expectation values are rounded or reconstructed by model fitting.
ORIGIN_INDEX.json separates original-file hashes from sanitized or subset packaged-file hashes. SHA256SUMS authenticates the packaged bytes only. Local absolute paths are replaced by logical root placeholders.
''')
    write('NEXT_DEEP_RESEARCH_PROMPT.md','''
# Independent review request

Review this executed negative/partial research round, beginning with README_FIRST.md, EXECUTIVE_FINDINGS.md, the frozen plan and blockers. Recompute tables with code/verify.py; use code/reanalyse.py for a separate score replay.
Do not request a generic literature review or search for a positive model. Challenge the evidence and its limits:

1. Are whole-chromosome matching, fitting, pooling and control calibration genuinely isolated, including duplicate loci and shared fragments? Does held-label matched sampling support only the stated conditional estimand?
2. Does zero-centered shrinkage reduce information or induce prior-dominated protection? Inspect zero-count flags, prior-width sensitivity and the independent numerical fixtures.
3. Is dispersion identified under the working likelihood, and does any improvement justify it over rho=0? Separate null mean error from variance and fragment dependence; do not assume mitochondrial-to-nuclear transfer.
4. Does conservative calibration conceal weak power? Compare exact attainable ceilings, model fixtures and whole-fragment perturbations using their actual independent source-block counts. Do not read stratified-window averages as population FPR.
5. Are claimed ranking differences independent of count/motif/accessibility confounding? Specifically challenge the measured full-depth coverage imbalance and strong covariate baseline. Recalculate every available task, retaining unavailable HMGA1 and the original 21-task denominator.
6. Are canonical central protection and both shoulders present in the prespecified held populations? A detector-only ranking change cannot strengthen an unchanged aggregate. Is an honest no-rescue conclusion warranted?
7. Are DWM support/layout/strand limitations correctly bounded, and are archived residual, raw propensity and production-workflow comparisons kept separate? This is not a direct TOBIAS benchmark.
8. Are functional uncertainty and full-workflow resource claims appropriately left incomplete? What single falsifiable next diagnostic, if any, is justified without another broad search?

Protected new-context labels and final naked-DNA replicate 3 must remain closed; do not suggest retuning against examined test/transfer labels. Preserve all original promotion gates. Recommend stopping if the complete evidence does not justify further model complexity.
''')
    reports={
      'A_INPUT_QUALIFICATION.md':'''Hypothesis: exact observations and qualified references allow paired comparison. Complete read-name/read-group pairs were inventoried and both retained endpoints checked against an independent coordinate oracle that reads the alignment operations (CIGAR). Full and exact 10M arrays were extracted for all six libraries and control replicate 1. Fixture edge cases and control old/new extraction parity are recorded. The original pre-matching pool exists; new held matching has 6,020 sites and replay has 8,106. Native control support covers 4,200 processed windows: 4,020 positive expectations and 180 valid zeros. Old-layout replay matches 281/281 historical profiles; expanded layout matches 279/281. Native strand/reference qualification remains INCOMPLETE. Retain source-specific limitations; this is not a safety qualification.''',
      'B_INFORMATION_CALIBRATION.md':f'''Hypothesis: non-calls can reflect insufficient information rather than reliable safety. Executed 210 model-based scenarios with 2,000 evaluation and separate calibration draws, and 138 nonoverlapping empirical source windows with 2,000 whole-fragment resamples per perturbation. Full source-block/dependence identities and unavailable cases are retained. Count/expected-count ceilings, strand/shoulder/GC/support null means and variances are available. {mismatched}/21 LOG21 central-mean task intervals exclude zero before multiplicity correction. Neither small call fractions nor working-model interval coverage qualify real biological calibration. No fragment perturbation is described as a binding model.''',
      'C_CANDIDATE_COMPARISON.md':'''Hypothesis: signed shrinkage improves protection ranking, and beta-binomial dispersion adds value beyond rho=0. Both fixed candidates and uniform/LOG21 old count arms plus DWM-propensity factorial controls ran on identical arrays. Archived native-DWM residual replay remains separate. Initial panel: six libraries, 10M/full, seed 2026. See full_screen_summary/task_deltas, complete fold/site tables, confound diagnostics and all-task shape curves. Actual full-depth coverage remains imbalanced in all 60 available task/library panels; ranking differences do not isolate footprint discrimination. Missing native full-depth residual and calibrated functional uncertainty prevent qualification. No detector-only track gain exists. Decision: no advancement and no new methods.''',
      'D_DEPTH_REPLICATION.md':'''All three biological libraries per cell were included in the initial 10M/full seed 2026 panel, using exact complete-pair subsampling. This is useful replication/depth development evidence. The full 25M/50M and 2027–2030 matrix was intentionally not run because the sole frozen full-depth advancement screen was not met. No favorable depth, task, replicate or seed was substituted. depth_replication.tsv retains every expected case.''',
      'E_RESOURCES_PACKAGE.md':'''The unchanged end_to_end_resource_protocol_v1.json remains binding. The formal six-repeat, four-core, isolated fit-and-apply matrix was intentionally not run after the unmet scientific screen. Output equivalence and absolute simultaneous process-tree memory are incomplete; exact deployable execution configuration cannot honestly be supplied for a score-only prototype. resource_results.tsv preserves all cases and blockers. The package contains complete headline observations/counts/scores and independent read-only verification. ZIP checksum/path integrity establishes transport integrity, not scientific qualification.'''}
    for name,body in reports.items():write('experiments/'+name,'# '+name.replace('_',' ').replace('.md','')+'\n\n'+body)
    write('provenance/PRIMARY_READING.md','''
# Bounded primary-source context

No upstream code was installed or substituted for frozen local implementations.
msCentipede: Raj et al., PLOS ONE2015, DOI10.1371/journal.pone.0138030; paper and official README inspected for variation across sites/replicates. Official repository pin ddd39e5cb1d38152b4e2c681f32a56f42b10766e: https://github.com/rajanil/msCentipede/tree/ddd39e5cb1d38152b4e2c681f32a56f42b10766e . This two-region model is not msCentipede.
PRINT/seq2PRINT: Hu et al., Nature2025, DOI10.1038/s41586-024-08443-4. Metadata/context checked; full-text access encountered site/recaptcha limitations in this round. Do not claim an independent algorithm reproduction or new protein-addition validation.
TraceBind: Lin et al., Nature Communications2026, DOI10.1038/s41467-026-73164-3; official README at9592bfcdc65e0dd755d4b7b430281c1d3794fb3e inspected for mitochondrial adaptation and coverage-stratified thresholding: https://github.com/lyx-lin/TraceBind/tree/9592bfcdc65e0dd755d4b7b430281c1d3794fb3e . No new bias fitting or dynamic geometry was adopted.
ChromBPNet: official README at packet pin09938fdb4397ec0006510e5251e48920a505d4de inspected: https://github.com/kundajelab/chrombpnet/tree/09938fdb4397ec0006510e5251e48920a505d4de . Bias, regulatory and no-bias predictions remain distinct historical products; no new neural-network experiment ran.
''')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,required=True);execute(p.parse_args().outdir.resolve())
