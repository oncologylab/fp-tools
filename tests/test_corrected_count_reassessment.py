import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks/scripts"))
import reassess_corrected_strand_counts as m


def fixture(tmp_path, monkeypatch):
    root = tmp_path
    output = root / "new"
    output.mkdir()
    bam = root / "dummy.bam"
    bam.write_bytes(b"fixture input")
    digest = m.file_sha256(bam)
    matrix = [{"cell": cell, "sample": f"{cell}_rep{rep}", "depth": depth, "seed": seed,
               **{key: str(bam) for key in ("bam", "bai", "dwm_expected")},
               **{key + "_sha256": digest for key in ("bam", "bai", "dwm_expected")}}
              for cell in ("HepG2", "K562") for rep in (1, 2, 3)
              for depth in ("10m", "25m", "50m", "full")
              for seed in ([2026] if depth == "full" else range(2026, 2031))]
    contexts, memberships = {}, {}
    for cell in ("HepG2", "K562"):
        sites = pd.DataFrame({"cell": [cell] * 4, "tf": ["CTCF"] * 4, "TFBS_chr": ["chr19"] * 4,
                              "TFBS_start": [1000, 2000, 3000, 4000], "TFBS_end": [1010, 2010, 3010, 4010],
                              "TFBS_strand": ["+", "-", "+", "-"], "chip_label": [0, 0, 1, 1]})
        cohorts = {"original_support_with_expected_coverage": list(range(4)),
                   "motif_accessibility_observed_total": list(range(4))}
        memberships[cell] = {"site_and_cohort_identity": m.extension.cohort_identity(sites, cohorts)}
        sp = output / f"{cell}.tsv"
        sites.to_csv(sp, sep="\t", index=False)
        ap = output / f"{cell}.npz"
        np.savez(ap, selected=np.arange(4), site_hash=m.base.site_hashes(sites),
                 plus_log_bias=np.zeros((4, 201)), minus_log_bias=np.zeros((4, 201)))
        tp = output / f"{cell}.tracks.npz"
        np.savez(tp, site_hash=m.base.site_hashes(sites), fixed_dwm=np.ones((4, 201)),
                 **{f"run_{i}": np.ones((4, 201)) for i in range(48)})
        frame = pd.DataFrame(matrix).query("cell == @cell").sort_values(m.base.RUN_KEYS)
        contexts[cell] = {"sites": m.record(sp), "arrays": m.record(ap), "tracks": m.record(tp),
                          "cohorts": cohorts, "runs": frame[m.base.RUN_KEYS].to_dict("records")}
    settings = {"dispersion": {"HepG2": 0., "K562": 1.6}}
    config = {"schema": m.SCHEMA, "evidence": "post_hoc_fixed_candidate_corrected_count_reassessment",
              "source_result_sha256": m.SOURCE_HASHES["stability_manifest.json"],
              "scientific_settings": settings, "matrix": matrix, "contexts": contexts,
              "inputs": [], "code": [], "runtime": m.extension.runtime(),
              **{key: False for key in ("candidate_retuned", "new_holdout_labels_read",
                 "final_untouched_control_read", "promotion_eligible", "archived_observed_profiles_reused")}}
    config["configuration_id"] = m.canonical_id(config)
    m.write_immutable(output / "configuration.freeze.json", config)
    monkeypatch.setattr(m, "source_documents", lambda: {
        "configuration.freeze.json": {"scientific_settings": settings, "matrix": matrix},
        "cohorts.freeze.json": {"memberships": memberships}})
    calls = []

    def extraction(sites, bam, arrays):
        calls.append(str(bam))
        return np.ones((4, 201)), np.ones((4, 201))

    def evaluate(row, sites, cohorts, methods, observed, plan, directory, identity, root):
        outputs = {}
        run = {k: getattr(row, k) for k in m.base.RUN_KEYS}
        for name in ("metrics", "bootstrap", "profiles", "shapes", "coverage"):
            p = directory / f"{name}.tsv.gz"
            pd.DataFrame([{**run, "tf": "CTCF", "n_sites": len(sites)}]).to_csv(p, sep="\t", index=False)
            outputs[name] = m.record(p)
        stage = {"configuration_id": identity, "run": run, "outputs": outputs}
        stage["stage_id"] = m.canonical_id(stage)
        m.write_immutable(directory / "evaluation.json", stage)

    monkeypatch.setattr(m, "corrected_arrays", extraction)
    monkeypatch.setattr(m.base, "evaluate_run", evaluate)
    return root, output, config, calls


def test_all_96_reextract_and_completed_outputs_verify_without_recomputation(tmp_path, monkeypatch):
    root, output, _, calls = fixture(tmp_path, monkeypatch)
    m.execute(root, output)
    assert len(calls) == 96
    m.execute(root, output, verify_only=True)
    assert len(calls) == 96
    result = json.loads((output / "reassessment_manifest.json").read_text())
    assert len(result["stages"]) == 96 and result["promotion_eligible"] is False
    profile = output / "HepG2_rep1/10m/seed_2026/profiles.npz"
    profile.write_bytes(b"corrupt")
    with pytest.raises(ValueError, match="checksum"):
        m.execute(root, output, verify_only=True)


@pytest.mark.parametrize("field", ["candidate_retuned", "new_holdout_labels_read", "archived_observed_profiles_reused"])
def test_retuning_and_old_count_reuse_rejected(tmp_path, monkeypatch, field):
    _, _, config, _ = fixture(tmp_path, monkeypatch)
    config[field] = True
    with pytest.raises(ValueError, match="scope"):
        m.validate_configuration(config)


def test_incomplete_extraction_is_preserved(tmp_path, monkeypatch):
    root, output, _, calls = fixture(tmp_path, monkeypatch)
    path = output / "HepG2_rep1/10m/seed_2026/profiles.npz"
    path.parent.mkdir(parents=True)
    path.write_bytes(b"partial")
    with pytest.raises(ValueError, match="incomplete extraction"):
        m.execute(root, output)
    assert not calls and path.read_bytes() == b"partial"


@pytest.mark.parametrize("change", ["dispersion", "matrix"])
def test_frozen_settings_cannot_be_changed(tmp_path, monkeypatch, change):
    _, _, config, _ = fixture(tmp_path, monkeypatch)
    altered = json.loads(json.dumps(config))
    if change == "dispersion":
        altered["scientific_settings"]["dispersion"]["K562"] = 0.
    else:
        altered["matrix"] = altered["matrix"][:-1]
    with pytest.raises(ValueError, match="settings or matrix"):
        m.validate_configuration(altered)
