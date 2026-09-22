#!/usr/bin/env python3
"""Reassess frozen normalization on identical tuning sites without refitting."""

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

sys.path.insert(0,str(Path(__file__).resolve().parent))
import freeze_shape_normalization as fixed
from audit_direct_bias_safety import record,write_immutable
from audit_chrombpnet_comparison import canonical_id,checked_path

ROOT=fixed.ROOT
MODELS=ROOT/"benchmarks/results/corrected_shape_normalization_20260922_v1"
MODELS_HASH="52986cf07c0861f39ef2a96966056ba7867b4d454ee45c3d475c9e7da2aa1f84"
TUNING_HASH="74af0ce9cd66bc9d06af81d3edc1f95f04f9ec7c84da703c6de00e645c97c9ec"
KEYS=["cell","tf","correction","detector"]
SITES=["cell","tf","TFBS_chr","TFBS_start","TFBS_end","TFBS_strand","chip_label"]


def replace_scores(metrics,stored,sites,values,cell,tf,correction):
    audit=[]
    for detector,score in values.items():
        mask=stored.cell.eq(cell)&stored.tf.eq(tf)&stored.correction.eq(correction)&stored.detector.eq(detector)
        old=stored.loc[mask]
        if old.empty:
            raise ValueError("a corrected model cannot introduce a new detector")
        pd.testing.assert_frame_equal(old[SITES].reset_index(drop=True),sites[SITES].reset_index(drop=True),check_dtype=False)
        score=np.asarray(score,dtype=float)
        if score.shape!=(len(old),) or not np.isfinite(score).all():raise ValueError("invalid reassessed scores")
        if detector!="hierarchical_spline":
            np.testing.assert_allclose(old.score,score,rtol=1e-12,atol=1e-12)
            continue
        target=metrics.cell.eq(cell)&metrics.tf.eq(tf)&metrics.correction.eq(correction)&metrics.detector.eq(detector)
        if int(target.sum())!=1:raise ValueError("metric task identity missing or duplicated")
        labels=sites.chip_label.to_numpy(int)
        auc=roc_auc_score(labels,score);ap=average_precision_score(labels,score)
        prior=metrics.loc[target].iloc[0]
        delta=np.abs(old.score.to_numpy()-score)
        audit.append(dict(cell=cell,tf=tf,correction=correction,detector=detector,sites=len(score),
            changed_scores=int((delta>1e-10).sum()),maximum_score_change=float(delta.max()),
            previous_auroc=prior.auroc,corrected_auroc=auc,previous_auprc=prior.auprc,corrected_auprc=ap))
        metrics.loc[target,"auroc"],metrics.loc[target,"auprc"]=auc,ap
        stored.loc[mask,"score"]=score
    return audit


def freeze(root,out):
    checked_path({"path":str(MODELS/"models.freeze.json"),"sha256":MODELS_HASH},ROOT)
    fixed.execute(root,MODELS,verify_only=True)
    tuning=fixed.SOURCE/"tuning_manifest.json"
    checked_path({"path":str(tuning),"sha256":TUNING_HASH},ROOT)
    fixed.parent.evaluate(root,fixed.SOURCE,verify_only=True)
    if out.exists():raise ValueError("use a fresh output directory")
    parent_config=json.loads((MODELS/"configuration.freeze.json").read_text())
    models=json.loads((MODELS/"models.freeze.json").read_text())
    result=json.loads(tuning.read_text())
    config=dict(models=record(MODELS/"models.freeze.json"),tuning=record(tuning),
        parent_configuration=record(fixed.SOURCE/"configuration.freeze.json"),
        inputs=parent_config["inputs"]+list(models["outputs"].values())+list(result["outputs"].values()),
        code=parent_config["code"]+[record(Path(__file__))],runtime=fixed.parent.base.source.builder.extension.runtime(),
        model_refitted=False,candidate_selected=False,test_labels_read=False,promotion_eligible=False)
    config["configuration_id"]=canonical_id(config)
    out.mkdir(parents=True);write_immutable(out/"configuration.freeze.json",config)
    print("Frozen identical-site normalization reassessment",flush=True)


def execute(root,out,verify_only=False):
    config=fixed.parent.base.source.builder.extension.read_document(out/"configuration.freeze.json","configuration_id")
    if (config["runtime"]!=fixed.parent.base.source.builder.extension.runtime()
            or config["models"]["sha256"]!=MODELS_HASH or config["tuning"]["sha256"]!=TUNING_HASH
            or any(config[k] is not False for k in ("model_refitted","candidate_selected","test_labels_read","promotion_eligible"))):
        raise ValueError("reassessment scope or runtime changed")
    for item in config["inputs"]+config["code"]+[config["models"],config["tuning"],config["parent_configuration"]]:checked_path(item,root)
    path=out/"reassessment_manifest.json"
    if path.exists():
        fixed.parent.base.source.builder.base.verify_stage(path,config["configuration_id"],root)
        print("Verified frozen-normalization tuning reassessment",flush=True);return
    if verify_only or list(out.glob("*.tsv")):raise ValueError("missing result or partial outputs preserved")
    parent_config=fixed.parent.preflight(root,fixed.SOURCE)
    sites,panels,_=fixed.parent.load_split(parent_config,root,"tune")
    original=json.loads(checked_path(config["tuning"],root).read_text())
    models=json.loads(checked_path(config["models"],root).read_text())
    metrics=pd.read_csv(checked_path(original["outputs"]["metrics"],root),sep="\t")
    stored=pd.read_csv(checked_path(original["outputs"]["scores"],root),sep="\t")
    if metrics.duplicated(KEYS).any():raise ValueError("duplicated metric task")
    rows=[]
    for item in models["models"]:
        ix=np.flatnonzero((sites.cell.eq(item["cell"])&sites.tf.eq(item["tf"])).to_numpy())
        model=checked_path(item["artifact"],root);values=panels[item["correction"]][ix]
        batch=fixed.scores(model,values)
        if "hierarchical_spline" in batch:
            single=np.array([fixed.scores(model,row[None,:])["hierarchical_spline"][0] for row in values])
            np.testing.assert_allclose(batch["hierarchical_spline"],single,rtol=1e-10,atol=1e-10)
        rows.extend(replace_scores(metrics,stored,sites.iloc[ix],batch,item["cell"],item["tf"],item["correction"]))
    outputs={"curves":original["outputs"]["curves"]}
    for name,table in (("metrics",metrics),("scores",stored),("changes",pd.DataFrame(rows))):
        p=out/(name+".tsv");table.to_csv(p,sep="\t",index=False);outputs[name]=record(p)
    result=dict(configuration_id=config["configuration_id"],outputs=outputs,models=config["models"],
                evidence="fixed_model_implementation_reassessment_on_tuning",model_refitted=False,
                candidate_selected=False,test_labels_read=False,promotion_eligible=False,
                aggregate_curves_scope="original_per_task_normalized_input_profiles_unchanged")
    result["stage_id"]=canonical_id(result);write_immutable(path,result)
    print("Completed fixed-normalization tuning reassessment",flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root",type=Path,required=True);parser.add_argument("--outdir",type=Path,required=True)
    mode=parser.add_mutually_exclusive_group();mode.add_argument("--freeze-only",action="store_true");mode.add_argument("--verify-only",action="store_true")
    args=parser.parse_args();out=args.outdir.resolve()
    if not out.is_relative_to(ROOT/"benchmarks/results"):parser.error("use research results")
    if args.freeze_only:freeze(args.data_root.resolve(),out)
    else:
        with fixed.parent.base.source.builder.extension.exclusive_run(out):execute(args.data_root.resolve(),out,args.verify_only)
