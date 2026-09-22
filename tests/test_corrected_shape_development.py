import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"benchmarks/scripts"))
import run_corrected_shape_development as m


def fixture(tmp_path,monkeypatch):
    rng=np.random.default_rng(2026)
    y=np.array([0]*20+[1]*20)
    sites=pd.DataFrame(dict(cell=["HepG2"]*40,tf=["A"]*40,motif_family=["family"]*40,
        TFBS_chr=["chr1"]*40,TFBS_start=np.arange(40)*300+1000,TFBS_end=np.arange(40)*300+1010,
        TFBS_strand=["+"]*40,chip_label=y))
    profiles=rng.normal(0,.1,(40,201));profiles[np.ix_(y==1,np.abs(m.candidates.POSITIONS)<=7)]-=2
    p=json.loads(m.cohorts.parent.base.PROTOCOL.read_text());p["minimum_sites_per_class"]=10
    config=dict(configuration_id="fixture",protocol=p)
    ledger=pd.DataFrame([dict(cell="HepG2",tf="A",status="eligible",matched_pairs=20)])
    calls=[]
    def load(config,root,split):
        calls.append(split)
        s=sites.copy()
        s["TFBS_chr"]="chr1" if split=="train" else "chr16"
        return s,{c:profiles.copy() for c in p["corrections"]},ledger
    monkeypatch.setattr(m,"preflight",lambda *a:config)
    monkeypatch.setattr(m,"load_split",load)
    return config,sites,profiles,calls


def test_frozen_models_precede_tuning_and_resume_without_refit(tmp_path,monkeypatch):
    _,_,_,calls=fixture(tmp_path,monkeypatch)
    m.fit(tmp_path,tmp_path)
    assert calls==["train"]
    frozen=(tmp_path/"models.freeze.json").read_bytes()
    m.fit(tmp_path,tmp_path)
    assert calls==["train"]
    m.evaluate(tmp_path,tmp_path)
    assert calls==["train","tune"]
    assert (tmp_path/"models.freeze.json").read_bytes()==frozen
    t=pd.read_csv(tmp_path/"tuning_metrics.tsv",sep="\t")
    assert len(t)==9 and t.auroc.notna().all()
    m.evaluate(tmp_path,tmp_path,verify_only=True)
    with (tmp_path/"model_000.npz").open("ab") as f:f.write(b"changed")
    with pytest.raises(ValueError):m.evaluate(tmp_path,tmp_path,verify_only=True)


def test_tuning_refuses_missing_model_freeze(tmp_path,monkeypatch):
    _,_,_,calls=fixture(tmp_path,monkeypatch)
    with pytest.raises(FileNotFoundError):m.evaluate(tmp_path,tmp_path)
    assert calls==[]


def test_artifact_scores_exactly_reproduce_fitted_scores(tmp_path,monkeypatch):
    _,sites,profiles,_=fixture(tmp_path,monkeypatch)
    models,_=m.candidates.fit_candidates(sites,profiles,minimum_per_class=10)
    model=models[("HepG2","A")]
    path=tmp_path/"frozen.npz";m.save_model(path,model)
    before=m.candidates.score_candidates(model,profiles)
    after=m.frozen_scores(path,profiles)
    for key in before:np.testing.assert_array_equal(before[key],after[key])
    with pytest.raises(ValueError):m.save_model(path,model)


def test_training_rejects_tuning_profiles_even_if_loader_is_wrong(tmp_path,monkeypatch):
    _,sites,profiles,_=fixture(tmp_path,monkeypatch)
    sites.TFBS_chr="chr16"
    monkeypatch.setattr(m,"load_split",lambda *a:(sites,{"raw":profiles},pd.DataFrame()))
    with pytest.raises(ValueError,match="training chromosomes"):m.fit(tmp_path,tmp_path)
