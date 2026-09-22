from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks/scripts"))
import audit_shape_chromosome_sensitivity as audit


def fixture():
    return pd.DataFrame(dict(cell="K562", tf="TEST", TFBS_chr=np.repeat(audit.BLOCKS, 4),
        TFBS_start=np.arange(12), TFBS_end=np.arange(12)+1, TFBS_strand="+",
        chip_label=[0, 0, 1, 1]*3, score=[.1, .4, .3, .8]*3))


def test_weighted_resample_matches_explicit_block_duplication():
    candidate = fixture()
    baseline = candidate.copy()
    baseline.score = baseline.score.iloc[::-1].to_numpy()
    paired = audit.paired_scores(candidate, baseline)
    weights = np.repeat([2, 1, 0], 4)
    weighted = audit.delta(paired, weights)
    repeated = paired.iloc[np.repeat(np.arange(12), weights)]
    explicit = audit.delta(repeated, np.ones(len(repeated)))
    assert weighted == pytest.approx(explicit)


def test_all_ordered_draws_are_kept_and_identical_methods_have_zero_delta():
    f = fixture()
    result = audit.resamples(audit.paired_scores(f, f))
    assert len(result) == 27 and result.blocks.nunique() == 27
    assert result.valid.all()
    assert not result.auroc_gain.any() and not result.auprc_gain.any()


def test_missing_class_draws_remain_explicit():
    f = fixture()
    f.loc[f.TFBS_chr.eq('chr16'), 'chip_label'] = 0
    result = audit.resamples(audit.paired_scores(f, f))
    assert len(result) == 27
    assert (~result.valid).sum() == 1
    assert result.loc[~result.valid, 'auroc_gain'].isna().all()


def test_reordered_sites_align_but_changed_labels_or_support_fail():
    f = fixture()
    assert audit.delta(audit.paired_scores(f, f.iloc[::-1]), np.ones(12))['auroc_gain'] == 0
    changed = f.copy()
    changed.loc[0, 'chip_label'] = 1
    with pytest.raises(ValueError, match='support differs'):
        audit.paired_scores(f, changed)
    with pytest.raises(ValueError, match='support differs'):
        audit.paired_scores(f, f.iloc[1:])


def test_test_chromosomes_duplicate_sites_and_invalid_weights_reject():
    f = fixture()
    with pytest.raises(ValueError, match='duplicate'):
        audit.paired_scores(pd.concat([f, f]), f)
    changed = f.copy()
    changed.loc[0, 'TFBS_chr'] = 'chr19'
    with pytest.raises(ValueError, match='tuning'):
        audit.paired_scores(changed, changed)
    with pytest.raises(ValueError, match='weights'):
        audit.delta(audit.paired_scores(f, f), -np.ones(12))
