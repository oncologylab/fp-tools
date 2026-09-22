import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks/scripts"))
import rebuild_corrected_safety_profiles as m


def inputs(cell="HepG2", role="validation"):
    sites = pd.DataFrame(dict(cell=[cell] * 2, tf=["CTCF"] * 2,
        TFBS_chr=["chr16"] * 2, TFBS_start=[1000, 2000], TFBS_end=[1010, 2010], TFBS_strand=["+", "-"]))
    if role == "validation":
        sites["chip_label"] = [0, 1]
    arrays = dict(site_hash=m.base.site_hashes(sites), valid=np.array([True, False]),
                  plus_log_bias=np.zeros((2, 201)), minus_log_bias=np.zeros((2, 201)))
    return sites, arrays


def test_scope_and_site_order():
    sites, arrays = inputs()
    m.validate_sites(sites, arrays, "HepG2", "validation")
    with pytest.raises(ValueError, match="labels"):
        m.validate_sites(sites, arrays, "HepG2", "naked")
    sites.loc[0, "TFBS_chr"] = "chr19"
    with pytest.raises(ValueError, match="non-development"):
        m.validate_sites(sites, arrays, "HepG2", "validation")
    sites, arrays = inputs()
    arrays["site_hash"] = arrays["site_hash"][::-1]
    with pytest.raises(ValueError, match="order"):
        m.validate_sites(sites, arrays, "HepG2", "validation")


def test_rebuild_preserves_strand_totals_orientation_and_invalid_rows(monkeypatch):
    sites, arrays = inputs()
    plus, minus = np.zeros((2, 201)), np.zeros((2, 201))
    plus[:, 20], minus[:, 30] = 3, 7
    arrays["plus_log_bias"][1] = np.nan
    arrays["minus_log_bias"][1] = np.nan
    monkeypatch.setattr(m.base, "extract_strand_cut_profiles", lambda *a, **k: (plus, minus, np.ones(2, dtype=bool)))
    out = m.rebuild(sites, arrays, "fixture")
    assert out["plus_observed"][1, 170] == 7
    assert out["minus_observed"][1, 180] == 3
    np.testing.assert_allclose(out["expected"][0].sum(), 10)
    assert out["valid"].tolist() == [True, False]
    assert out["archived_valid"].tolist() == [True, False]
    assert len(out["observed"]) == 2


def fixture(tmp_path, monkeypatch):
    protocol = tmp_path / "protocol.json"
    protocol.write_text("{}")
    monkeypatch.setattr(m, "PROTOCOL_SHA", m.record(protocol)["sha256"])
    output = tmp_path / "outputs"
    output.mkdir()
    items = []
    for cell in ("HepG2", "K562"):
        for role in ("validation", "naked"):
            sites, arrays = inputs(cell, role)
            p = tmp_path / f"{cell}.{role}.tsv"
            sites.to_csv(p, sep="\t", index=False)
            z = p.with_suffix(".npz")
            np.savez(z, **arrays)
            items.append(dict(cell=cell, role=role, sites=m.record(p), arrays=m.record(z),
                              **{k: m.record(protocol) for k in ("bam", "bai", "manifest")}))
    config = dict(schema="fp-tools-corrected-safety-profiles-v1", runtime=m.extension.runtime(),
                  protocol=m.record(protocol), source_configuration=m.record(protocol), items=items,
                  code=[], profiles_scored=False, candidate_fitted=False, final_control_read=False)
    config["configuration_id"] = m.canonical_id(config)
    m.write_immutable(output / "configuration.freeze.json", config)
    monkeypatch.setattr(m.base, "extract_strand_cut_profiles", lambda *a, **k: (np.ones((2, 201)), np.ones((2, 201)), np.ones(2, dtype=bool)))
    return output


def test_all_four_inputs_resume_and_verify_without_extraction(tmp_path, monkeypatch):
    output = fixture(tmp_path, monkeypatch)
    m.execute(tmp_path, output)
    monkeypatch.setattr(m.base, "extract_strand_cut_profiles", lambda *a, **k: pytest.fail("unexpected rerun"))
    m.execute(tmp_path, output, verify_only=True)
    m.execute(tmp_path, output)
    assert len(json.loads((output / "profiles_manifest.json").read_text())["stages"]) == 4
    with (output / "HepG2.naked.npz").open("ab") as f:
        f.write(b"changed")
    with pytest.raises(ValueError):
        m.execute(tmp_path, output, verify_only=True)


def test_partial_extraction_is_preserved(tmp_path, monkeypatch):
    output = fixture(tmp_path, monkeypatch)
    p = output / "HepG2.validation.npz"
    p.write_bytes(b"partial")
    with pytest.raises(ValueError, match="partial output"):
        m.execute(tmp_path, output)
    assert p.read_bytes() == b"partial"


def test_missing_stage_fails_verification(tmp_path, monkeypatch):
    output = fixture(tmp_path, monkeypatch)
    with pytest.raises(ValueError, match="missing corrected"):
        m.execute(tmp_path, output, verify_only=True)
