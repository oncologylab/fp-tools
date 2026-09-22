import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks/scripts"))
import reassess_corrected_safety as m


def fixture(tmp_path, monkeypatch, missing_baseline=False):
    n = 220
    sites = pd.DataFrame(dict(cell=["HepG2"]*n, tf=["CTCF"]*n, TFBS_chr=["chr16"]*n,
        TFBS_start=np.arange(n)*300+1000, TFBS_end=np.arange(n)*300+1010, TFBS_strand=["+"]*n, chip_label=[0]*n))
    control = sites.drop(columns="chip_label")
    arrays = dict(observed=np.ones((n, 201)), expected=np.ones((n, 201)), valid=np.ones(n, dtype=bool))
    baseline = np.ones((n, 201))
    if missing_baseline:
        baseline[0] = 0
    old = pd.DataFrame([dict(cell="HepG2", tf="CTCF", method=k, threshold=1.0, validation_negatives=n, threshold_eligible=True) for k in m.safety.METHODS])
    path = tmp_path / "source"
    path.write_text("fixture")
    protocol = json.loads(m.builder.PROTOCOL.read_text())
    protocol["source_thresholds"] = m.safety.record(path)
    config = dict(configuration_id="fixture", protocol=protocol, corrected_stages={"HepG2":dict(validation="validation",naked="naked")},
                  source_configuration=dict(configuration_id="old", cells={"HepG2":dict(dwm_track=m.safety.record(path), dispersion=0.0, naked="naked", naked_dwm="dwm")}))
    monkeypatch.setattr(m, "preflight", lambda *a:config)
    monkeypatch.setattr(m, "read_profiles", lambda entry, root:(sites if entry=="validation" else control, arrays))
    monkeypatch.setattr(m.safety, "extract_prediction_profiles", lambda *a, **k:(baseline,np.ones(n,dtype=bool)))
    monkeypatch.setattr(m.safety, "read_thresholds", lambda *a:({},old))
    def load(entry, *a, **k):
        assert k["control"]
        return control, arrays["observed"], baseline if entry=="dwm" else arrays["expected"], arrays["valid"], {}
    monkeypatch.setattr(m.safety,"load_artifact",load)
    monkeypatch.setattr(m.safety,"motif_residual_effect",lambda *a,**k:{"motif_residual_flag":False})
    out = tmp_path/"output"
    out.mkdir()
    return out, arrays


def test_calibration_does_not_score_control_and_both_analyses_complete(tmp_path,monkeypatch):
    out, _ = fixture(tmp_path,monkeypatch)
    load = m.safety.load_artifact
    monkeypatch.setattr(m.safety,"load_artifact",lambda *a,**k:pytest.fail("control loaded during calibration"))
    m.calibrate(tmp_path,out)
    m.calibrate(tmp_path,out)
    threshold_bytes=(out/"thresholds.tsv").read_bytes()
    monkeypatch.setattr(m.safety,"load_artifact",load)
    m.evaluate(tmp_path,out)
    table=pd.read_csv(out/"corrected_safety.tsv",sep="\t")
    assert set(table.analysis)==set(m.ANALYSES)
    assert table.total_sites.tolist()==[220,220]
    assert table.passes_candidate_screen.all()
    assert (out/"thresholds.tsv").read_bytes()==threshold_bytes
    m.evaluate(tmp_path,out,verify_only=True)
    with (out/"corrected_sites.tsv").open("a") as f:f.write("corrupt")
    with pytest.raises(ValueError):m.evaluate(tmp_path,out,verify_only=True)


def test_missing_calibration_baseline_prevents_relative_safety_pass(tmp_path,monkeypatch):
    out,_=fixture(tmp_path,monkeypatch,missing_baseline=True)
    m.calibrate(tmp_path,out)
    m.evaluate(tmp_path,out)
    table=pd.read_csv(out/"corrected_safety.tsv",sep="\t")
    assert not table.baseline_calibration_qualified.any()
    assert not table.passes_candidate_screen.any()
    assert table.absolute_fpr_screen_pass.all()


def test_control_requires_frozen_thresholds(tmp_path,monkeypatch):
    out,_=fixture(tmp_path,monkeypatch)
    monkeypatch.setattr(m,"read_profiles",lambda *a:pytest.fail("control accessed before threshold freeze"))
    with pytest.raises(FileNotFoundError):m.evaluate(tmp_path,out)


def test_invalid_sites_are_retained_in_output_and_denominator(tmp_path,monkeypatch):
    out,arrays=fixture(tmp_path,monkeypatch)
    arrays["valid"][0]=False
    m.calibrate(tmp_path,out)
    m.evaluate(tmp_path,out)
    table=pd.read_csv(out/"corrected_safety.tsv",sep="\t")
    assert table.invalid_sites.tolist()==[1,1]
    assert table.support.tolist()==[219,219]
    assert len(pd.read_csv(out/"corrected_sites.tsv",sep="\t"))==220


def test_partial_scoring_is_not_overwritten(tmp_path,monkeypatch):
    out,_=fixture(tmp_path,monkeypatch)
    m.calibrate(tmp_path,out)
    p=out/"corrected_sites.tsv";p.write_text("partial")
    with pytest.raises(ValueError,match="partial"):m.evaluate(tmp_path,out)
    assert p.read_text()=="partial"
