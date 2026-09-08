from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks/scripts"))


def module():
    import evaluate_chrombpnet_covariate_diagnostic
    return evaluate_chrombpnet_covariate_diagnostic


def sites(features, labels, chromosomes=None):
    n = len(labels)
    return pd.DataFrame({"cell": ["K562"] * n, "tf": ["CTCF"] * n,
                         "TFBS_chr": chromosomes or ["chr19"] * n,
                         "TFBS_start": np.arange(n) * 20, "TFBS_end": np.arange(n) * 20 + 10,
                         "TFBS_strand": ["+"] * n, "chip_label": labels,
                         "motif_score": features})


def test_matching_uses_training_scales_only():
    m = module()
    train = sites([1., 3., 999., 999.], [0, 1, 0, 1], ["chr1", "chr2", "chr16", "chr19"])
    scale = m.fit_scales(train, ["motif_score"], ["chr1", "chr2"])
    assert scale["motif_score"]["mean"] == 2
    assert scale["motif_score"]["sd"] == pytest.approx(np.sqrt(2))
    assert scale["training_sites"] == 2


def test_assignment_maximizes_valid_pairs_before_minimizing_distance():
    m = module()
    # Positive 0 can match both negatives; positive .18 can only match .1.
    frame = sites([0., .18, .1, -.19], [1, 1, 0, 0])
    scale = {"motif_score": {"mean": 0., "sd": 1.}}
    matched = m.match_sites(frame, ["motif_score"], scale, caliper=.2)
    assert len(matched) == 4
    assert matched.groupby("pair_id").chip_label.nunique().eq(2).all()
    assert matched.match_maximum_standardized_difference.le(.2).all()


def test_matching_is_independent_of_model_scores_and_row_order():
    m = module()
    frame = sites([0., .18, .1, -.19, 50.], [1, 1, 0, 0, 0])
    scale = {"motif_score": {"mean": 0., "sd": 1.}}
    frame["candidate_score"] = np.arange(5)
    first = m.match_sites(frame, ["motif_score"], scale, caliper=.2)
    frame["candidate_score"] = -frame.candidate_score
    second = m.match_sites(frame.sample(frac=1, random_state=55), ["motif_score"], scale, caliper=.2)
    pd.testing.assert_frame_equal(first[m.SITE_KEYS + ["pair_id"]], second[m.SITE_KEYS + ["pair_id"]])


def test_pairing_never_crosses_chromosomes_or_accepts_nonfinite_covariates():
    m = module()
    frame = sites([0., 0., 0., np.nan], [1, 0, 0, 1], ["chr19", "chr20", "chr19", "chr19"])
    scale = {"motif_score": {"mean": 0., "sd": 1.}}
    matched = m.match_sites(frame, ["motif_score"], scale, caliper=.2)
    assert len(matched) == 2
    assert matched.TFBS_chr.eq("chr19").all()


def test_degenerate_training_covariate_requires_exact_agreement():
    m = module()
    frame = sites([3., 3., 4., 4.], [1, 0, 1, 0])
    scale = {"motif_score": {"mean": 3., "sd": 0.}}
    matched = m.match_sites(frame, ["motif_score"], scale, caliper=.2)
    assert len(matched) == 2
    assert matched.motif_score.eq(3).all()


def test_unopened_cell_plan_is_rejected_before_input_read(tmp_path):
    import json
    m = module()
    plan = json.loads((Path(m.__file__).parents[1] / "manifests/chrombpnet_covariate_diagnostic_v1.json").read_text())
    plan["cells"] = ["SK-N-SH"]
    with pytest.raises(ValueError, match="development"):
        m.resolve_plan(plan, tmp_path)


def test_no_power_result_keeps_readable_tables_and_verified_resume(tmp_path):
    m = module()
    selected = sites([1., 1.], [1, 0])
    selected["pair_id"] = "pair"
    plan = {"primary_policy": "motif_accessibility"}
    result = m.evaluate_policy("K562", "motif_only", selected, selected, selected, {},
                               pd.DataFrame(), plan, tmp_path, tmp_path)
    assert pd.read_csv(tmp_path / "metrics.tsv", sep="\t").empty
    assert pd.read_csv(tmp_path / "bootstrap.tsv", sep="\t").empty
    assert not result["promotion_eligible"]
    again = m.evaluate_policy("K562", "motif_only", selected, selected, selected, {},
                               pd.DataFrame(), plan, tmp_path, tmp_path)
    assert again == result
    (tmp_path / "metrics.tsv").write_text("changed\n")
    with pytest.raises(ValueError, match="checksum"):
        m.evaluate_policy("K562", "motif_only", selected, selected, selected, {},
                           pd.DataFrame(), plan, tmp_path, tmp_path)


def test_raw_track_supplies_training_covariate_and_must_reproduce_archived_values():
    m = module()
    frame = sites([1., 2.], [0, 1])
    observed = np.asarray([[1., 2.], [3., 4.]])
    extracted = np.asarray([6., 14.])
    result = m.add_covariates(frame, observed, accessibility=extracted)
    np.testing.assert_allclose(result.log_accessibility, np.log1p(extracted))
    np.testing.assert_allclose(result.observed_total, [3., 7.])
    frame["accessibility"] = extracted
    assert m.add_covariates(frame, observed, accessibility=extracted).accessibility.equals(frame.accessibility)
    with pytest.raises(ValueError, match="archived accessibility"):
        m.add_covariates(frame, observed, accessibility=extracted * 2)
