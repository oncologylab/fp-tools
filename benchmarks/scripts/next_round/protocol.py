"""Frozen scientific policy. Changes after scoring require a preserved attempt."""
POLICY = {
    'schema': 'fp-tools-next-round-v1',
    'hypotheses': ['signed shrinkage improves ranking over significance',
                   'control overdispersion adds useful information beyond binomial shrinkage',
                   'low-count non-calls can reflect inadequate attainable significance',
                   'ranking changes do not change identical-site measured protection'],
    'training_chromosomes': [f'chr{i}' for i in range(1,16)],
    'tuning_allowed': False, 'protected_inputs_allowed': False,
    'task_denominator': 21, 'geometry': [7,12,40], 'canonical_geometry': [5,15,40],
    'bias': {'source': 'five-seed mitochondrial LOG21', 'shift':[4,-5],
             'l2':.0001, 'sha256':'620f232f634423b8753e4d77f8cf72df87482c09193cb5e7e6b038917c47e89a',
             'refit':False},
    'models': {'primary':'bias-offset beta-binomial with theta ~ Normal(0,1)',
               'fallback':'same model with rho=0', 'ranking':'posterior mean signed theta',
               'strand_likelihood':'working composite, dependence not solved',
               'rho_range':[0,.95], 'rho_fitting':'unique rep1 control windows on fitting chromosomes, theta=0',
               'unidentified':'profile support width > .25 or upper boundary or numerical failure',
               'prior_sd_sensitivity':[.5,2.], 'sensitivity_selection_allowed':False},
    'integration': {'method':'mode-split adaptive quadrature', 'independent_oracle':'768-node Gauss-Legendre using SciPy PMFs',
                    'mean_probability_atol':1e-8, 'interval_atol':1e-6, 'quadrature_tolerance':2e-11},
    'calls': {'old_threshold':1.3010299957, 'old_cap':300.,'old_decimals':10,
              'new_rule':'ceil(.95*(n+1)) order statistic; strict greater-than',
              'minimum_unique_informative_controls':200, 'pooling':'one library/depth/propensity, across tasks deduplicated by genomic window',
              'fitting':'leave held chromosome out', 'zero_information_call':False},
    'cohorts': {'replay':'exact original 8106 sites, nonnested retrospective development',
                'primary':'original pre-match pool, held-chromosome matching with fitting-fold-only scaling',
                'covariates':['motif_score','log1p_central_201','log1p_outer_101_500'],
                'caliper':.5, 'assignment':'maximum cardinality then minimum squared distance',
                'balance':'largest distance-sorted prefix with every absolute pooled-SD SMD <= .1',
                'minimum_matched_pairs_for_task_gate':200,
                'anchor':'recorded rep1 10M counts, fixed membership across new libraries/depths',
                'grouping':'whole chromosome held across all tasks and contexts; overlapping windows recorded',
                'missing_task':'HepG2 HMGA1 retained, never counted as zero or pass'},
    'endpoints': {'primary_depth':'10M seed2026', 'initial_panel':'all three biological libraries per context at full and 10M where complete-pair qualified',
                  'only_advancement_screen':'full depth: equal held chromosome then equal biological library then equal task',
                  'incomplete_panel':'no advancement; available endpoints are diagnostic only',
                  'arms':['archived_DWM_residual','uniform_count_tail','LOG21_count_tail','LOG21_binomial',
                          'LOG21_beta_binomial','DWM_binomial','DWM_beta_binomial'],
                  'unsupported_DWM':'no zero filling; mark arm incomplete and preserve other diagnostic arms'},
    'bootstrap': {'draws':2000,'seed':2026,'unit':'synchronized chromosome blocks',
                  'undefined':'no redraw; fewer than1900 valid draws is incomplete',
                  'scope':'fixed predictions, not full-refit uncertainty'},
    'model_simulation': {'repetitions':2000,'seed':2026,'n_pairs':[[0,0],[1,0],[2,2],[5,5],[10,10],[25,25],[50,50]],
                         'p_pairs':[[.2,.2],[1e-4,1e-4],[.05,.4]], 'theta':[0.,.5,1.,2.],
                         'rho':[0.,.1], 'misspecification':'null mean shift log odds +.5 and strand asymmetry +.5/-.5; correlated endpoints in empirical challenge'},
    'empirical_simulation': {'source':'naked DNA replicate1 only; complete fragments',
                             'window_selection':'coordinate-hash order, nonoverlap, max10 per count/expected stratum',
                             'r':[1.,.75,.5,.25,0.], 'repetitions':2000,'seed':2026,
                             'weight':'r**m per whole fragment, m central endpoints; 0**0=1',
                             'normalization':'normalize within window, sample K source fragments with replacement',
                             'empty_weight_population':'unavailable', 'truth':'realized marginal depletion, not binding biology'},
    'count_strata':['0','1','2-4','5-9','10-19','20-49','50-99','100+'],
    'expected_strata':['<0.1','0.1-<0.5','0.5-<1','1-<3','3+'],
    'shape': {'principal_population':'all held matched sites, not score-selected',
              'representations':['raw','sequence expectation','signed deviance residual'],
              'summaries':['site equal','count weighted','strand specific','combined','count and expected count strata'],
              'standardized_definition':'unchanged standardized_functional_separation on fixed outside-flank normalization',
              'detector_identity':'same counts/expectations/cohort imply exactly unchanged aggregate profile'},
    'advancement': {'requirements':'all original scientific numeric, shape, context, positive-control and principal interval gates, plus qualified development FPR/Wilson/paired DWM',
                    'missing_evidence':'not met', 'fallback_preferred':True,
                    'primary_extra_requirement':'at least .01 mean AUROC over fallback, positive paired interval, no AUPRC or calibration worsening',
                    'both_fail':'stop search; package negative and incomplete evidence'},
    'resources': {'elapsed_hours':72,'package_reserve_hours':4,'science_max_cores':32,
                  'maximum_ram_GiB':64,'maximum_scratch_GiB':500,'formal_benchmark_cores':4,
                  'formal_protocol':'end_to_end_resource_protocol_v1.json unchanged',
                  'formal_matrix_if_screen_not_met':'intentionally not run'},
    'promotion_possible_this_round':False,
}

PROMOTION = '''Mean AUROC gain at least 0.03 over conventional DWM; relative AUPRC gain at least 10%;
at least 25% stronger standardized aggregate separation with canonical protection in at least two difficult TF families;
improvement in at least two TF families and two cell contexts; no CTCF or other positive-control AUROC loss greater than 0.02;
principal chromosome-block bootstrap interval excludes zero; untouched naked-DNA false-positive rate at most 5%,
Wilson upper bound at most 5%, and no increase over DWM greater than one percentage point;
no evidence the bias component learned TF occupancy; stable direction across depths, seeds and biological replicates;
complete-workflow CPU runtime at most 2x DWM, peak memory at most 1.5x DWM, total deployed model at most 25 MB;
functional uncertainty coverage must be calibrated. GP/hybrid must materially outperform the simpler hierarchical spline
or provide substantially better calibrated uncertainty within the runtime gate. Otherwise retain the simpler method
or keep FDA/GP as research diagnostics.

AUROC: area under the receiver operating characteristic curve. AUPRC: average precision implementation of area
under the precision-recall curve. DWM: dinucleotide weight matrix. TF: transcription factor.
CTCF: CCCTC-binding factor. CPU: central processing unit. GP: Gaussian process. FDA: functional data analysis.
The archived DWM arm is an fp-tools comparator, not an end-to-end TOBIAS benchmark.
Locked SK-N-SH/GM23338 labels and final naked-DNA replicate3 remain closed throughout this round.
'''
