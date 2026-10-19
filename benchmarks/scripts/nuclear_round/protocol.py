"""Prespecified control-first nuclear LOG21 experiment; no outcome-driven search."""
CHROMS = [f'chr{i}' for i in range(1,16)]
SEEDS = list(range(2026,2031))
POLICY = dict(
    schema='nuclear-bias-round-v1', primary_l2=.01, ablation_l2=.0001,
    chromosomes=CHROMS, seeds=SEEDS, shifts=[4,-5], context=21, pair_distances=[1,5,9],
    windows_per_fold=3000, window_size=200, margin=41,
    sampling='smallest blake2b coordinate hashes among genome-grid windows with ACGT sequence and no halo blacklist overlap; no count or label selection',
    fitting='leave one training chromosome out, five seeds per fold; equal coefficient ensemble',
    extraction='both members of unique reciprocal primary MAPQ>=30 nonduplicate nonQCfail read-name/read-group pairs; independent CIGAR oracle',
    optimizer=dict(epochs=100,batch_windows=64,learning_rate=.03,tolerance=1e-7,patience=12),
    detector=dict(geometry=[7,12,40],cap=300,round_decimals=10,threshold=1.3010299957,comparison='strict greater than'),
    screen=dict(bootstrap_draws=2000,seed=2026,minimum_valid_draws=1900,
        conditional_likelihood='paired chromosome sum, unique controls, lower95%>0 versus mitochondrial and intercept-only',
        task_means='at most 3/21 intervals excluding observed/expected=1',
        controls='all tasks informative Wilson upper <=.05 and same-site DWM call increase <=.01; native DWM qualification separate',
        primary_failure='stop; do not substitute ablation or retune'),
    protected_contents_allowed=False, tuning_allowed=False, chromatin_until_control_pass=False,
    resource_budget=dict(max_workers=32,max_memory_gib=64,max_new_scratch_gib=500,elapsed_hours=72,packaging_reserve_hours=4))
