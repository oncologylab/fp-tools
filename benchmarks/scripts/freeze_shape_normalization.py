#!/usr/bin/env python3
"""Add training-derived normalization to copies of frozen shape model artifacts."""

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0,str(Path(__file__).resolve().parent))
import run_corrected_shape_development as parent
import fixed_shape_normalization as normalization
from audit_direct_bias_safety import record,write_immutable
from audit_chrombpnet_comparison import checked_path,canonical_id

ROOT=parent.ROOT
SOURCE=ROOT/"benchmarks/results/corrected_shape_development_20260922_v1"
SOURCE_HASH="a9780ac8006fabc1c945803c0ae0149a4bbcf69be73f95470e546eb48a28b675"


def scores(path,profiles):
    values=np.asarray(profiles,dtype=float)
    result=parent.candidates.score_candidates(None,values)
    with np.load(path,allow_pickle=False) as data:
        result["trained_geometry"]=values@parent.candidates.geometry_weights(data["geometry"])
        if "discriminant" in data:
            if "normalization_fallback" not in data:
                raise ValueError("spline requires training-frozen normalization")
            scale=float(data["score_scale"])
            if not np.isfinite(scale) or scale<=0:
                raise ValueError("invalid frozen score scale")
            normalized=normalization.normalize_fixed(values,parent.candidates.POSITIONS,float(data["normalization_fallback"]))
            result["hierarchical_spline"]=(normalized@data["discriminant"]-float(data["midpoint"]))/scale
    return result


def copy_model(source,path,training_profiles):
    if path.exists():
        raise ValueError("do not overwrite model artifacts")
    with np.load(source,allow_pickle=False) as data:
        arrays={k:data[k] for k in data.files}
    if "normalization_fallback" in arrays:
        raise ValueError("source model already has a frozen fallback")
    arrays["normalization_fallback"]=np.array(normalization.training_fallback(training_profiles,parent.candidates.POSITIONS))
    np.savez_compressed(path,**arrays)
    with np.load(path,allow_pickle=False) as restored:
        for name,values in arrays.items():
            np.testing.assert_array_equal(values,restored[name])
    before=parent.frozen_scores(source,training_profiles)
    after=scores(path,training_profiles)
    for key in before:
        np.testing.assert_allclose(before[key],after[key],rtol=1e-10,atol=1e-10)
    if "hierarchical_spline" in after:
        singleton=np.array([scores(path,row[None,:])["hierarchical_spline"][0] for row in training_profiles])
        np.testing.assert_allclose(after["hierarchical_spline"],singleton,rtol=1e-10,atol=1e-10)
        legacy_singleton=np.array([parent.frozen_scores(source,row[None,:])["hierarchical_spline"][0] for row in training_profiles])
        difference=np.abs(before["hierarchical_spline"]-legacy_singleton)
    else:
        difference=np.zeros(len(training_profiles))
    return dict(training_sites=len(training_profiles),fallback=float(arrays["normalization_fallback"]),
                legacy_batch_dependent_sites=int((difference>1e-10).sum()),
                legacy_maximum_score_difference=float(difference.max()),fixed_batch_invariant=True,
                learned_weights_unchanged=True)


def freeze(root,out):
    checked_path({"path":str(SOURCE/"models.freeze.json"),"sha256":SOURCE_HASH},ROOT)
    old=parent.preflight(root,SOURCE)
    models=parent.verify_models(old,root,SOURCE)
    if out.exists():raise ValueError("use a fresh output directory")
    config=dict(parent_configuration=record(SOURCE/"configuration.freeze.json"),parent_models=record(SOURCE/"models.freeze.json"),
                inputs=old["inputs"]+list(models["outputs"].values()),
                code=old["code"]+[record(Path(__file__)),record(Path(normalization.__file__))],
                runtime=parent.base.source.builder.extension.runtime(),training_chromosomes=old["protocol"]["train_chromosomes"],
                model_refitted=False,tuning_scored=False,new_holdout_labels_read=False,promotion_eligible=False)
    config["configuration_id"]=canonical_id(config)
    out.mkdir(parents=True);write_immutable(out/"configuration.freeze.json",config)
    print("Frozen training-only model normalization correction",flush=True)


def execute(root,out,verify_only=False):
    config=parent.base.source.builder.extension.read_document(out/"configuration.freeze.json","configuration_id")
    if (config["runtime"]!=parent.base.source.builder.extension.runtime() or config["parent_models"]["sha256"]!=SOURCE_HASH
            or any(config[k] is not False for k in ("model_refitted","tuning_scored","new_holdout_labels_read","promotion_eligible"))):
        raise ValueError("normalization correction scope changed")
    for entry in config["inputs"]+config["code"]+[config["parent_configuration"],config["parent_models"]]:checked_path(entry,root)
    path=out/"models.freeze.json"
    if path.exists():
        done=parent.base.source.builder.base.verify_stage(path,config["configuration_id"],root)
        if done.get("model_refitted") is not False or done.get("tuning_scored") is not False:
            raise ValueError("completed model scope changed")
        for entry in done["models"]:checked_path(entry["artifact"],root)
        print("Verified training-frozen normalization models",flush=True);return
    if verify_only or list(out.glob("model_*.npz")) or (out/"normalization_audit.tsv").exists():
        raise ValueError("missing result or partial model artifacts preserved")
    old=parent.preflight(root,SOURCE);models=parent.verify_models(old,root,SOURCE)
    if config["training_chromosomes"]!=old["protocol"]["train_chromosomes"]:
        raise ValueError("training chromosomes changed")
    sites,panels,_=parent.load_split(old,root,"train")
    if not sites.TFBS_chr.isin(config["training_chromosomes"]).all():raise ValueError("non-training input")
    rows,new_models,outputs=[],[],{}
    for index,item in enumerate(models["models"]):
        ix=sites.cell.eq(item["cell"])&sites.tf.eq(item["tf"])
        p=out/f"model_{index:03d}.npz"
        row=copy_model(checked_path(item["artifact"],root),p,panels[item["correction"]][ix])
        row.update(cell=item["cell"],tf=item["tf"],correction=item["correction"]);rows.append(row)
        new_models.append(dict(item,artifact=record(p)));outputs[p.name]=record(p)
    p=out/"normalization_audit.tsv";pd.DataFrame(rows).to_csv(p,sep="\t",index=False);outputs["audit"]=record(p)
    result=dict(configuration_id=config["configuration_id"],models=new_models,outputs=outputs,
                parent_models=config["parent_models"],model_refitted=False,tuning_scored=False,promotion_eligible=False)
    result["stage_id"]=canonical_id(result);write_immutable(path,result)
    print("Frozen normalization for all 11 model artifacts without refitting",flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root",type=Path,required=True);parser.add_argument("--outdir",type=Path,required=True)
    mode=parser.add_mutually_exclusive_group();mode.add_argument("--freeze-only",action="store_true");mode.add_argument("--verify-only",action="store_true")
    args=parser.parse_args();out=args.outdir.resolve()
    if not out.is_relative_to(ROOT/"benchmarks/results"):parser.error("use research results")
    if args.freeze_only:freeze(args.data_root.resolve(),out)
    else:
        with parent.base.source.builder.extension.exclusive_run(out):execute(args.data_root.resolve(),out,args.verify_only)
