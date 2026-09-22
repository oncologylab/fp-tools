#!/usr/bin/env python3
"""Freeze, fit on chr1-15, then compare frozen detectors on chr16-18 only."""

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
import balance_corrected_shape_cohorts as cohorts
import corrected_shape_candidates as candidates
from audit_direct_bias_safety import record, write_immutable
from audit_chrombpnet_comparison import canonical_id, checked_path
from fp_tools.tools.functional_footprints import standardized_functional_separation

ROOT = cohorts.ROOT
SOURCE = ROOT / "benchmarks/results/corrected_shape_cohorts_20260922_v3"
SOURCE_HASH = "9b19e5644b4c7cb60ee2d5477b156e927953a60810a0e87f236cd5ce19a568b3"
base = cohorts.parent.base


def save_model(path, model):
    arrays = {"geometry": np.asarray(model["geometry"], dtype=int)}
    detector = model["spline"]
    if detector is not None:
        arrays.update(discriminant=detector.discriminant_, midpoint=np.array(detector.midpoint_),
                      score_scale=np.array(detector.score_scale_), template=detector.footprint_template_,
                      template_standard_error=detector.template_standard_error_, pooled_variance=detector.pooled_variance_)
    if path.exists():
        raise ValueError("model artifact must not be overwritten")
    np.savez_compressed(path, **arrays)


def frozen_scores(path, profiles):
    values = np.asarray(profiles, dtype=float)
    if values.ndim != 2 or values.shape[1] != 201 or not np.isfinite(values).all():
        raise ValueError("finite 201-position profiles required")
    result = candidates.score_candidates(None, values)
    if path is None:
        return result
    with np.load(path, allow_pickle=False) as data:
        result["trained_geometry"] = values @ candidates.geometry_weights(data["geometry"])
        if "discriminant" in data:
            scale = float(data["score_scale"])
            if data["discriminant"].shape != (201,) or not np.isfinite(scale) or scale <= 0:
                raise ValueError("invalid frozen spline parameters")
            result["hierarchical_spline"] = (candidates.normalize_functional_profiles(values, candidates.POSITIONS) @ data["discriminant"] - float(data["midpoint"])) / scale
    return result


def freeze(root, out):
    checked_path({"path": str(SOURCE / "cohorts_manifest.json"), "sha256": SOURCE_HASH}, ROOT)
    cohorts.execute(root, SOURCE, verify_only=True)
    if out.exists():
        raise ValueError("freeze requires a fresh directory")
    old = json.loads((SOURCE / "configuration.freeze.json").read_text())
    result = json.loads((SOURCE / "cohorts_manifest.json").read_text())
    source_config = json.loads(checked_path(old["protocol"]["dwm_source_configuration"], ROOT).read_text())
    paths = [Path(__file__), Path(candidates.__file__), ROOT / "src/fp_tools/tools/functional_footprints.py",
             ROOT / "benchmarks/scripts/evaluate_parametric_factorization.py",
             ROOT / "benchmarks/scripts/evaluate_chrombpnet_reference.py"]
    config = dict(protocol=old["protocol"], cohort_configuration=record(SOURCE / "configuration.freeze.json"),
                  cohort_manifest=record(SOURCE / "cohorts_manifest.json"),
                  inputs=old["inputs"]+list(result["outputs"].values()),
                  code=old["code"]+[record(p) for p in paths], runtime=base.source.builder.extension.runtime(),
                  dispersion={c:d["dispersion"] for c,d in source_config["cells"].items()},
                  evidence="training_and_tuning_candidate_development", test_labels_read=False,
                  new_holdout_labels_read=False, final_control_read=False, promotion_eligible=False)
    config["configuration_id"] = canonical_id(config)
    out.mkdir(parents=True)
    write_immutable(out / "configuration.freeze.json", config)
    print("Frozen training/tuning candidate experiment", flush=True)


def preflight(root, out):
    config = base.source.builder.extension.read_document(out / "configuration.freeze.json", "configuration_id")
    if (config["runtime"] != base.source.builder.extension.runtime()
            or config.get("evidence") != "training_and_tuning_candidate_development"
            or config["cohort_manifest"]["sha256"] != SOURCE_HASH
            or any(config.get(k) is not False for k in ("test_labels_read", "new_holdout_labels_read", "final_control_read", "promotion_eligible"))):
        raise ValueError("candidate experiment scope or runtime changed")
    for item in config["inputs"]+config["code"]+[config["cohort_manifest"],config["cohort_configuration"]]:
        checked_path(item, root)
    previous = json.loads(checked_path(config["cohort_configuration"],root).read_text())
    if config["protocol"] != previous["protocol"] or config["protocol"]["geometry_grid"] != [list(g) for g in candidates.GEOMETRIES]:
        raise ValueError("candidate or cohort protocol changed")
    return config


def load_split(config, root, split):
    if split not in ("train","tune"):
        raise ValueError("only training and tuning splits allowed")
    manifest = json.loads(checked_path(config["cohort_manifest"],root).read_text())
    tables, panels, ledgers = [], {c:[] for c in config["protocol"]["corrections"]}, []
    for entry in config["protocol"]["profile_inputs"]:
        cell = entry["cell"]
        sites, arrays = base.source.read_profiles(entry, ROOT)
        pairs = pd.read_csv(checked_path(manifest["outputs"][cell+"_pairs"],root),sep="\t")
        diag = pd.read_csv(checked_path(manifest["outputs"][cell+"_diagnostics"],root),sep="\t")
        ledgers.append(diag.loc[diag.split.eq(split)].copy())
        selected = pairs.loc[pairs.split.eq(split)].site_index.to_numpy(int)
        subset = sites.iloc[selected].copy().reset_index(drop=True)
        permitted = config["protocol"]["train_chromosomes" if split=="train" else "tune_chromosomes"]
        if not subset.TFBS_chr.isin(permitted).all():
            raise ValueError("cohort crosses chromosome split")
        with np.load(checked_path(manifest["outputs"][cell+"_dwm"],root),allow_pickle=False) as data:
            if not np.array_equal(data["site_hash"],arrays["site_hash"]) or not data["support"][selected].all():
                raise ValueError("expected support identity changed")
            dwm = data["expected"][selected]
        observed = arrays["observed"][selected]
        tables.append(subset)
        panels["raw"].append(observed)
        for name, expected in (("DWM_expected",dwm),("parametric_direct",arrays["expected"][selected])):
            panels[name].append(base.source.safety.direct_scores(observed,expected,config["dispersion"][cell])[2])
    return pd.concat(tables,ignore_index=True), {c:np.concatenate(v) for c,v in panels.items()}, pd.concat(ledgers,ignore_index=True)


def fit(root, out):
    config = preflight(root,out)
    manifest = out / "models.freeze.json"
    if manifest.exists():
        verify_models(config,root,out)
        print("Verified frozen training models; no refit",flush=True)
        return
    if list(out.glob("model_*.npz")) or (out/"training.tsv").exists():
        raise ValueError("partial training artifacts preserved")
    sites, panels, ledger = load_split(config,root,"train")
    models, rows, outputs = [], [], {}
    for correction, profiles in panels.items():
        fitted, summary = candidates.fit_candidates(sites,profiles,
            minimum_per_class=config["protocol"]["minimum_sites_per_class"], prior_strength=config["protocol"]["spline"]["prior_strength"])
        for task in ledger.itertuples(index=False):
            selected = summary.loc[summary.cell.eq(task.cell)&summary.tf.eq(task.tf)]
            row = selected.iloc[0].to_dict() if len(selected) else dict(cell=task.cell,tf=task.tf,status="no_training_pairs")
            row.update(correction=correction, cohort_status=task.status)
            key=(task.cell,task.tf)
            if key in fitted:
                if task.status!="eligible":
                    raise ValueError("model fitted on ineligible cohort")
                path=out/f"model_{len(models):03d}.npz"
                save_model(path,fitted[key])
                ix=sites.cell.eq(task.cell)&sites.tf.eq(task.tf)
                before=candidates.score_candidates(fitted[key],profiles[ix])
                after=frozen_scores(path,profiles[ix])
                for name in before:
                    np.testing.assert_allclose(before[name],after[name],rtol=0,atol=0)
                model=dict(cell=task.cell,tf=task.tf,correction=correction,artifact=record(path))
                models.append(model);outputs[path.name]=record(path)
            rows.append(row)
    path=out/"training.tsv"
    pd.DataFrame(rows).to_csv(path,sep="\t",index=False)
    outputs["training"]=record(path)
    result=dict(configuration_id=config["configuration_id"],models=models,outputs=outputs,
                training_chromosomes=config["protocol"]["train_chromosomes"],tuning_scored=False,promotion_eligible=False)
    result["stage_id"]=canonical_id(result)
    write_immutable(manifest,result)
    print(f"Frozen {len(models)} training model artifacts before tuning evaluation",flush=True)


def verify_models(config,root,out):
    result=base.source.builder.base.verify_stage(out/"models.freeze.json",config["configuration_id"],root)
    if (result.get("tuning_scored") is not False or result.get("promotion_eligible") is not False
            or result["training_chromosomes"]!=config["protocol"]["train_chromosomes"]):
        raise ValueError("models were not frozen before tuning")
    keys=[(m["cell"],m["tf"],m["correction"]) for m in result["models"]]
    if len(keys)!=len(set(keys)):
        raise ValueError("duplicate frozen model")
    for m in result["models"]:
        checked_path(m["artifact"],root)
    return result


def evaluate(root,out,verify_only=False):
    config=preflight(root,out)
    frozen=verify_models(config,root,out)
    manifest=out/"tuning_manifest.json"
    if manifest.exists():
        done=base.source.builder.base.verify_stage(manifest,config["configuration_id"],root)
        if done["models"]!=record(out/"models.freeze.json") or done.get("promotion_eligible") is not False:
            raise ValueError("tuning result uses different models or scope")
        print("Verified completed tuning comparison",flush=True)
        return
    if verify_only or list(out.glob("tuning_*.tsv")):
        raise ValueError("missing result or partial tuning output preserved")
    lookup={(m["cell"],m["tf"],m["correction"]):m["artifact"] for m in frozen["models"]}
    sites,panels,ledger=load_split(config,root,"tune")
    rows,curves,scores_out=[],[],[]
    for task in ledger.itertuples(index=False):
        ix=np.flatnonzero((sites.cell.eq(task.cell)&sites.tf.eq(task.tf)).to_numpy())
        y=sites.iloc[ix].chip_label.to_numpy(int)
        for correction,profiles in panels.items():
            model=lookup.get((task.cell,task.tf,correction))
            scores=frozen_scores(checked_path(model,root) if model else None,profiles[ix]) if len(ix) else {}
            norm=candidates.normalize_functional_profiles(profiles[ix],candidates.POSITIONS) if len(ix) else np.empty((0,201))
            shape={}
            if len(np.unique(y))==2:
                delta=norm[y==1].mean(axis=0)-norm[y==0].mean(axis=0)
                center=delta[np.abs(candidates.POSITIONS)<=5].mean()
                left=delta[(candidates.POSITIONS>=-40)&(candidates.POSITIONS<=-15)].mean()
                right=delta[(candidates.POSITIONS>=15)&(candidates.POSITIONS<=40)].mean()
                shape=dict(center_contrast=center,left_shoulder=left,right_shoulder=right,
                    asymmetry=(right-left)/(abs(right)+abs(left)) if abs(right)+abs(left)>0 else 0,
                    point_estimate_protection=candidates.protection_direction(delta),
                    functional_separation=standardized_functional_separation(norm,y,candidates.POSITIONS))
                curves.append(pd.DataFrame(dict(cell=task.cell,tf=task.tf,correction=correction,position=candidates.POSITIONS,
                    positive_mean=norm[y==1].mean(axis=0),negative_mean=norm[y==0].mean(axis=0),positive_minus_negative=delta)))
            for detector in config["protocol"]["detectors"]:
                row=dict(cell=task.cell,tf=task.tf,correction=correction,detector=detector,
                         cohort_status=task.status,matched_pairs=task.matched_pairs,
                         status=task.status if detector in scores else "no_frozen_candidate",**shape)
                if detector in scores and len(np.unique(y))==2:
                    row.update(auroc=roc_auc_score(y,scores[detector]),auprc=average_precision_score(y,scores[detector]))
                    frame=sites.iloc[ix][["cell","tf","TFBS_chr","TFBS_start","TFBS_end","TFBS_strand","chip_label"]].copy()
                    frame["correction"],frame["detector"],frame["score"]=correction,detector,scores[detector]
                    scores_out.append(frame)
                rows.append(row)
    outputs={}
    for name,table in (("metrics",pd.DataFrame(rows)),("curves",pd.concat(curves)),("scores",pd.concat(scores_out))):
        path=out/f"tuning_{name}.tsv"
        table.to_csv(path,sep="\t",index=False);outputs[name]=record(path)
    result=dict(configuration_id=config["configuration_id"],outputs=outputs,models=record(out/"models.freeze.json"),
                evidence="training_and_tuning_candidate_development",promotion_eligible=False,
                uncertainty_calibrated=False,candidate_selected=False,test_labels_read=False)
    result["stage_id"]=canonical_id(result)
    write_immutable(manifest,result)
    print("Completed tuning comparison without refitting or candidate selection",flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root",type=Path,required=True)
    parser.add_argument("--outdir",type=Path,required=True)
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument("--freeze-only",action="store_true")
    mode.add_argument("--fit-only",action="store_true")
    mode.add_argument("--verify-only",action="store_true")
    args=parser.parse_args();out=args.outdir.resolve()
    if not out.is_relative_to(ROOT/"benchmarks/results"):
        parser.error("use research results")
    if args.freeze_only:
        freeze(args.data_root.resolve(),out)
    else:
        with base.source.builder.extension.exclusive_run(out):
            if args.fit_only:fit(args.data_root.resolve(),out)
            else:evaluate(args.data_root.resolve(),out,args.verify_only)
