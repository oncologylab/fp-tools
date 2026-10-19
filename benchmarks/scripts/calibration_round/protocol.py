"""Frozen policy for the October 2026 calibration follow-up.

This file is hashed into the plan freeze before any new scoring. Changing it
after the freeze requires a new attempt directory; earlier attempts are kept.
"""

SOURCE = 'next_round_20260922_v1'
CHROMOSOMES = [f'chr{i}' for i in range(1, 16)]
LIBRARIES = [f'{cell}_rep{r}' for cell in ['HepG2', 'K562'] for r in (1, 2, 3)]
DEPTHS = ['full', '10M']
THRESHOLD = 1.3010299957
CAP = 300.
DECIMALS = 10
GEOMETRY = (7, 12, 40)

# Arms share identical observations, sites and detector; only the null
# central probability differs. DWM here is a raw sequence-propensity factorial
# control, not native production DWM correction.
ARMS = {
    'uniform_count_tail': 'uniform central share 15/73 (existing arm, recomputed for parity)',
    'DWM_count_tail': 'library DWM propensities, unchanged count-tail detector (factorial control)',
    'LOG21_count_tail': 'frozen LOG21 propensities, unchanged count-tail detector (reference; existing arm)',
    'LOG21NC1_count_tail': 'ablation: LOG21 logit shifted by one naked-DNA intercept per held fold',
    'LOG21NC_count_tail': 'primary: LOG21 logit recalibrated by naked-DNA intercept and slope per held fold',
}
PRIMARY = 'LOG21NC_count_tail'
ABLATION = 'LOG21NC1_count_tail'
REFERENCE = 'LOG21_count_tail'
COMPARATOR = 'DWM_count_tail'

POLICY = {
    'schema': 'fp-tools-calibration-round-v1',
    'question': 'Does a label-free naked-DNA correction of the LOG21 null central mean fix control '
                'miscalibration without losing discrimination, and does it change the information '
                'available at chromatin sites?',
    'addresses': 'verified 7/21 LOG21 control central-mean discrepancies; variance was modeled before the mean',
    'simplest_alternative': 'unchanged LOG21 count-tail; intercept-only ablation tests whether a global shift suffices',
    'falsified_if': ['held-chromosome control log-likelihood does not improve over LOG21',
                     'control discrepant-task count is not reduced to at most 3/21',
                     'discrimination falls more than 0.01 AUROC below LOG21 count-tail'],
    'added_complexity': 'two scalars per deployment; requires a naked-DNA (or equivalent protein-free) control library',
    'recalibration': {
        'model': 'central_count_s ~ Binomial(n_s, expit(a + b*logit(p_s))) per strand observation',
        'ablation': 'b fixed at 1',
        'p_clip': [1e-6, 1 - 1e-6],
        'training_rows': 'naked-DNA replicate 1 full depth; eligible; LOG21 valid on both strands; '
                         'deduplicated physical window representative; strand total > 0',
        'labels_used': False, 'chromatin_data_used': False,
        'folds': 'leave one chromosome out over chr1-15; chromatin and control sites on chromosome f use fit excluding f',
        'deployable_record': 'additional all-chr1-15 fit, recorded but never used for evaluation',
        'optimizer': 'Newton-Raphson on the binomial log-likelihood, tolerance 1e-12, maximum 100 iterations',
    },
    'detector': 'conditional_protection_score.depletion_score unchanged; cap 300; rounded to 10 decimals',
    'calls': {'threshold': THRESHOLD, 'rule': 'strict greater-than', 'informative': 'conditioned total > 0',
              'abstain': 'attainable ceiling (zero central cuts) <= threshold; abstention is not a negative call',
              'rule_name': 'count_tail_nominal_abstain_v1'},
    'populations': {
        'control': 'all 4,200 naked-DNA replicate-1 task windows; exclusions retained with reasons',
        'chromatin': 'nested held-training-chromosome cohort (6,020 sites), six libraries, depths full and 10M, seed 2026',
        'cohort_changes': 'none; original nested cohort only; no rematching',
        'tasks': 21, 'unavailable': 'HepG2 HMGA1 retained as unavailable',
    },
    'control_endpoints': {
        'heldout_loglik': 'sum over held folds of binomial log-likelihood, deduplicated windows, both cells pooled',
        'task_mean': 'per task ratio sum observed central / sum expected central over common-support windows '
                     '(eligible; LOG21 and DWM valid on both strands); 2000 synchronized chromosome-block '
                     'bootstrap draws, seed 2026, percentile 95% interval on the ratio',
        'discrepant': 'ratio interval excludes 1 (equivalently log ratio excludes 0)',
        'variance': 'Pearson chi-square / df after each arm mean, held folds, deduplicated windows; diagnostic only',
        'curves': 'observed versus expected central share by fixed LOG21 p bins and GC tertiles',
        'calls': 'per task calls, informative and eligible Wilson upper bounds; paired call-rate difference against DWM_count_tail',
    },
    'chromatin_endpoints': {
        'metrics': 'sklearn roc_auc_score and average_precision_score per held chromosome',
        'aggregation': 'equal chromosome, then equal biological library, then equal available task (next-round screen)',
        'relative_ap': 'per task mean AP difference divided by mean reference AP, then mean over tasks',
        'bootstrap': 'synchronized chromosome multinomial weights, 2000 draws, seed 2026, >=1900 valid draws',
        'stratified': 'within fold and count-stratum AUROC, reported as diagnostic',
        'covariate_baseline': 'existing fitting-fold-only covariate predictions, reused without refitting',
        'information': 'abstention fraction and expected central count among positives',
        'shape': 'regional central log(O/E) by label for each arm, full depth, three libraries pooled per cell; '
                 'protection effect = positive minus negative; chromosome-block interval',
        'power': 'existing whole-fragment synthetic challenge draws re-scored with each arm; diagnostic, not population power',
    },
    'difficult_families': {'MEF2': ['MEF2A', 'MEF2D'], 'ZNF': ['ZNF362', 'ZNF384', 'ZNF558'],
                           'FOX': ['FOXA1'], 'TBP': ['TBP'], 'bHLH': ['MYC'], 'TCF': ['TCF7L2'],
                           'ARID': ['ARID3A'], 'HMGA': ['HMGA1']},
    'positive_control': ['CTCF'],
    'advancement_screen': {
        'depth': 'full', 'comparator': COMPARATOR,
        'A_calibration': 'primary held-out loglik > reference, and discrepant tasks <= 3/21 and fewer than reference',
        'B_accuracy': 'mean AUROC gain >= 0.03 and mean relative AP gain >= 0.10 versus comparator; AUROC interval excludes zero',
        'C_reference': 'primary mean AUROC not below LOG21_count_tail by more than 0.01',
        'D_positive_control': 'no CTCF task AUROC loss > 0.02 versus comparator',
        'E_contexts': 'AUROC gain > 0 versus comparator in tasks from >= 2 families and in both cells',
        'F_shape': 'canonical protection (positive central log O/E < 0 and positive-minus-negative < 0, both intervals '
                   'excluding zero) in >= 2 difficult families',
        'G_controls': 'every task (all 21 evaluable on controls) informative Wilson upper <= 0.05 and paired call '
                      'increase over comparator <= 0.01; examined replicate 1, not untouched safety',
        'H_direction': 'AUROC gain over comparator positive in >= 5/6 libraries at both depths',
        'ablation_preference': 'prefer ablation if primary minus ablation held-out loglik < 1.92',
        'outcomes': {'all': 'DEVELOPMENT_CRITERION_MET; still NOT_TESTED_LOCKED for untouched controls and new contexts',
                     'A_only': 'calibration component development-qualified; general replacement rejected',
                     'A_fails': 'calibration candidate rejected'},
        'missing_evidence': 'not met',
    },
    'not_run_unless_advanced': ['seed panel', 'replicate-2/3 naked DNA', 'complete-workflow resource protocol'],
    'protected': ['ENCSR587TRP', 'ENCSR485TLP', 'SK-N-SH', 'GM23338', 'NakedDNA_rep3', 'naked_rep3',
                  'naked_dna_10m_rep3', 'SRR13450942', 'nutrient'],
    'bootstrap_seed': 2026, 'bootstrap_draws': 2000,
}
