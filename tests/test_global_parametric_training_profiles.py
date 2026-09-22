from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks/scripts"))
import build_global_parametric_training_profiles as builder


def sites():
    return pd.DataFrame(dict(cell=['HepG2']*3, tf=['A', 'B', 'C'], TFBS_chr=['chr1', 'chr1', 'chr16'],
        TFBS_start=[400, 450, 600], TFBS_end=[401, 451, 601], TFBS_strand=['+', '-', '+'],
        chip_label=['do not deserialize labels']*3))


def test_training_coordinate_reader_does_not_retain_labels_or_tuning(tmp_path):
    path = tmp_path/'sites.tsv'
    sites().to_csv(path, sep='\t', index=False)
    frame, selection = builder.training_sites(path, 'HepG2', ['chr1'])
    assert len(frame) == 2 and 'chip_label' not in frame
    assert selection.tolist() == [True, True, False]
    with pytest.raises(ValueError, match='coordinate'):
        builder.training_sites(path, 'K562', ['chr1'])


def test_real_sequence_inference_is_genomic_and_overlap_consistent(tmp_path):
    import pysam
    from fp_tools.tools.parametric_bias import BiasFeatureSpec, ConditionalSequenceBiasModel
    path = tmp_path/'genome.fa'
    path.write_text('>chr1\n'+'ACGTTGCA'*150+'\n')
    pysam.faidx(str(path))
    model = ConditionalSequenceBiasModel(BiasFeatureSpec.selma10())
    model.main[4, 0] = 1
    frame = sites().iloc[:2].drop(columns='chip_label')
    plus, minus = builder.infer_log_bias(frame, model, path)
    assert plus.shape == minus.shape == (2, 401)
    np.testing.assert_allclose(plus[0, 50:], plus[1, :351], rtol=1e-12)
    np.testing.assert_allclose(minus[0, 50:], minus[1, :351], rtol=1e-12)
    counts = np.ones_like(plus)
    expected, valid = builder.strand_expectation(counts, plus)
    result = builder.overlap_summary(frame, counts[:, 100:301], expected, valid)
    assert result['shared_valid_positions'] == 151
    assert result['expected_conflicts'] == result['observed_conflicts'] == 0


def test_overlap_assertions_detect_count_validity_and_expectation_conflicts():
    frame = sites().iloc[:2]
    observed, expected, valid = np.ones((2, 201)), np.ones((2, 201)), np.ones((2, 201), bool)
    for target in ('counts', 'expectation', 'validity'):
        obs, exp, mask = observed.copy(), expected.copy(), valid.copy()
        if target == 'counts': obs[1, 0] = 2
        elif target == 'expectation': exp[1, 0] = 2
        else: mask[1, 0] = False
        with pytest.raises(AssertionError):
            builder.overlap_summary(frame, obs, exp, mask)


def test_missing_sequence_and_zero_background_remain_unavailable():
    counts = np.ones((2, 401))
    bias = np.zeros_like(counts)
    bias[0, 200] = np.nan
    counts[1] = 0
    expected, valid = builder.strand_expectation(counts, bias)
    assert not valid.any() and np.isnan(expected).all()
    with pytest.raises(ValueError, match='401'):
        builder.strand_expectation(counts[:, :201], bias[:, :201])
