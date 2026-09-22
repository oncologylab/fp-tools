#!/usr/bin/env python3
"""Screen frozen shape operators on previously examined naked-DNA replicate 2."""

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import freeze_shape_normalization as fixed
import reassess_corrected_safety as parent
from audit_chrombpnet_comparison import canonical_id, checked_path

ROOT = fixed.ROOT
PROTOCOL = ROOT / "benchmarks/manifests/fixed_shape_operator_safety_v1.json"
PROTOCOL_SHA = "938f16fc8c96c1e8a0d51fbe579bd8b48363b69a22cd372794528f01a7a174dc"
safety, builder = parent.safety, parent.builder
KEYS = ["cell", "tf", "correction", "detector"]


def support(arrays, dwm, dwm_valid):
    """Calibration requires real expectation support in both correction arms."""
    return (arrays["valid"] & dwm_valid & (arrays["observed"].sum(axis=1) > 0)
            & (arrays["expected"].sum(axis=1) > 0) & (dwm.sum(axis=1) > 0))


def rate_summary(candidate, baseline, totals, valid, comparable, threshold, baseline_threshold, eligible):
    candidate, baseline, totals = map(np.asarray, (candidate, baseline, totals))
    valid = np.asarray(valid) & np.isfinite(candidate) & np.isfinite(totals)
    informative = valid & (totals > 0)
    paired = informative & comparable & np.isfinite(baseline)
    calls = informative & (candidate > threshold) & bool(eligible)
    baseline_calls = paired & (baseline > baseline_threshold) & bool(eligible)
    row = dict(candidate_valid_sites=int(valid.sum()), informative_sites=int(informative.sum()),
               baseline_comparable_sites=int(paired.sum()), false_positive_calls=int(calls.sum()),
               paired_candidate_calls=int(calls[paired].sum()), paired_baseline_calls=int(baseline_calls.sum()))
    passes = bool(eligible and np.isfinite(threshold) and np.isfinite(baseline_threshold))
    for name, mask in (("", valid), ("informative_", informative)):
        n = int(mask.sum())
        rate = calls.sum() / n if n else np.nan
        upper = safety.wilson_interval(int(calls.sum()), n)[1]
        row[name + "false_positive_rate"] = rate
        row[name + "wilson_upper_bound"] = upper
        passes &= bool(n and rate <= .05 and upper <= .05)
    n = int(paired.sum())
    increase = (int(calls[paired].sum()) - int(baseline_calls.sum())) / n if n else np.nan
    row["paired_increase_over_DWM"] = increase
    row["passes_fpr_screen"] = bool(passes and n and increase <= .01)
    row["threshold_evaluated"] = bool(eligible and np.isfinite(threshold) and np.isfinite(baseline_threshold))
    if not row["threshold_evaluated"]:
        for name in ("false_positive_rate", "wilson_upper_bound", "informative_false_positive_rate",
                     "informative_wilson_upper_bound", "paired_increase_over_DWM"):
            row[name] = np.nan
    return row, calls, baseline_calls


def freeze(root, out):
    protocol = json.loads(checked_path(dict(path=str(PROTOCOL), sha256=PROTOCOL_SHA), ROOT).read_text())
    if out.exists():
        raise ValueError("use a fresh output directory")
    sources = {k: checked_path(protocol[k], ROOT) for k in
               ("models", "corrected_profile_configuration", "direct_safety_configuration", "corrected_tuning_result")}
    direct = parent.preflight(root, sources["direct_safety_configuration"].parent)
    fixed.execute(root, sources["models"].parent, verify_only=True)
    models = json.loads(sources["models"].read_text())["models"]
    model_config = json.loads((sources["models"].parent / "configuration.freeze.json").read_text())
    paths = {e["path"] for e in direct["code"] + model_config["code"]}
    paths.add(str(Path(__file__)))
    inputs = direct["inputs"] + [safety.record(p) for p in sources.values()]
    inputs += [m["artifact"] for m in models]
    config = dict(protocol=protocol, models=models, direct=direct, inputs=inputs,
                  code=[safety.record(Path(p)) for p in sorted(paths)], runtime=builder.extension.runtime(),
                  promotion_eligible=False, new_holdout_labels_read=False, model_refitted=False)
    config["configuration_id"] = canonical_id(config)
    out.mkdir(parents=True)
    safety.write_immutable(out / "configuration.freeze.json", config)
    print("Frozen all model, input, code and runtime identities", flush=True)


def preflight(root, out):
    config = builder.extension.read_document(out / "configuration.freeze.json", "configuration_id")
    protocol = json.loads(checked_path(dict(path=str(PROTOCOL), sha256=PROTOCOL_SHA), ROOT).read_text())
    if (config["protocol"] != protocol or config["runtime"] != builder.extension.runtime()
            or any(config[k] is not False for k in ("promotion_eligible", "new_holdout_labels_read", "model_refitted"))):
        raise ValueError("frozen scope or runtime changed")
    for item in config["inputs"] + config["code"]:
        checked_path(item, root)
    return config


def load_cell(config, root, cell, control):
    direct = config["direct"]
    source = direct["source_configuration"]["cells"][cell]
    sites, arrays = parent.read_profiles(direct["corrected_stages"][cell]["naked" if control else "validation"], root)
    if control:
        if any("label" in c.lower() or "chip" in c.lower() for c in sites.columns):
            raise ValueError("control must not contain labels")
        ds, observed, dwm, dv, _ = safety.load_artifact(source["naked_dwm"], root, cell, control=True)
        order = pd.Index(safety.site_hashes(ds)).get_indexer(safety.site_hashes(sites))
        if np.any(order < 0):
            raise ValueError("control site support missing")
        observed, dwm, dv = observed[order], dwm[order], dv[order]
        comparable, _ = parent.coverage.comparable_profiles(arrays["observed"], observed, dwm)
        comparable &= arrays["valid"] & dv
    else:
        if not sites.TFBS_chr.isin([f"chr{i}" for i in range(1, 19)]).all():
            raise ValueError("non-development calibration input")
        dwm, dv = safety.extract_prediction_profiles(sites, checked_path(source["dwm_track"], root), 100, require_dense=False)
        dwm = np.nan_to_num(dwm)
        comparable = support(arrays, dwm, dv)
    panels = {"raw": arrays["observed"]}
    baseline, _, panels["DWM_expected"] = safety.direct_scores(arrays["observed"], dwm, source["dispersion"])
    panels["parametric_direct"] = safety.direct_scores(arrays["observed"], arrays["expected"], source["dispersion"])[2]
    valid = {"raw": arrays["valid"],
             "DWM_expected": arrays["valid"] & dv & (dwm.sum(axis=1) > 0),
             "parametric_direct": arrays["valid"] & ((arrays["expected"].sum(axis=1) > 0)
                                                      | (arrays["observed"].sum(axis=1) == 0))}
    return sites, arrays["observed"].sum(axis=1), panels, baseline, valid, comparable


def operators(config, root, cell, tf, panels, ix):
    models = {(m["cell"], m["tf"], m["correction"]): m for m in config["models"]}
    for correction in config["protocol"]["corrections"]:
        model = models.get((cell, tf, correction))
        scores = fixed.scores(checked_path(model["artifact"], root), panels[correction][ix]) if model else {}
        for detector in config["protocol"]["detectors"]:
            yield correction, detector, scores.get(detector)


def calibrate(root, out):
    config = preflight(root, out)
    if (out / "thresholds.freeze.json").exists():
        builder.base.verify_stage(out / "thresholds.freeze.json", config["configuration_id"], root)
        print("Verified frozen thresholds; no recalibration", flush=True)
        return
    if (out / "thresholds.tsv").exists():
        raise ValueError("partial calibration preserved")
    rows = []
    for cell in config["protocol"]["cells"]:
        sites, totals, panels, baseline, _, common = load_cell(config, root, cell, False)
        for tf, task in sites.groupby("tf"):
            ix = task.index.to_numpy()
            negatives = task.TFBS_chr.isin(["chr16", "chr17", "chr18"]).to_numpy() & task.chip_label.eq(0).to_numpy()
            mask = negatives & common[ix] & np.isfinite(baseline[ix])
            for correction, detector, scores in operators(config, root, cell, tf, panels, ix):
                usable = mask & np.isfinite(scores) if scores is not None else np.zeros(len(ix), bool)
                n = int(usable.sum())
                eligible = n >= config["protocol"]["minimum_calibration_negatives"]
                rows.append(dict(cell=cell, tf=tf, correction=correction, detector=detector,
                    calibration_negatives=n, calibration_exclusions=int(negatives.sum())-n,
                    threshold_eligible=eligible, model_available=scores is not None,
                    candidate_threshold=safety.conformal_upper_threshold(scores[usable], .05) if eligible else np.nan,
                    baseline_threshold=safety.conformal_upper_threshold(baseline[ix][usable], .05) if eligible else np.nan))
    table = pd.DataFrame(rows)
    if len(table) != 126 or table.duplicated(KEYS).any():
        raise ValueError("calibration task denominator changed")
    path = out / "thresholds.tsv"
    table.to_csv(path, sep="\t", index=False)
    receipt = dict(configuration_id=config["configuration_id"], outputs={"thresholds": safety.record(path)}, control_profiles_scored=False)
    receipt["stage_id"] = canonical_id(receipt)
    safety.write_immutable(out / "thresholds.freeze.json", receipt)
    print("All 126 threshold ledger rows frozen before control scoring", flush=True)


def evaluate(root, out, verify_only=False):
    config = preflight(root, out)
    frozen = builder.base.verify_stage(out / "thresholds.freeze.json", config["configuration_id"], root)
    if frozen.get("control_profiles_scored") is not False:
        raise ValueError("threshold provenance changed")
    path = out / "screen_manifest.json"
    if path.exists():
        result = builder.base.verify_stage(path, config["configuration_id"], root)
        if (result["thresholds"] != safety.record(out / "thresholds.freeze.json")
                or result.get("promotion_eligible") is not False
                or result.get("evidence") != config["protocol"]["evidence"]):
            raise ValueError("result scope or threshold identity changed")
        print("Verified completed operator false-positive screen", flush=True)
        return
    if verify_only or list(out.glob("screen_*.tsv")):
        raise ValueError("missing result or partial outputs preserved")
    thresholds = pd.read_csv(checked_path(frozen["outputs"]["thresholds"], root), sep="\t").set_index(KEYS)
    rows, site_rows = [], []
    for cell in config["protocol"]["cells"]:
        sites, totals, panels, baseline, valid, comparable = load_cell(config, root, cell, True)
        for tf, task in sites.groupby("tf"):
            ix = task.index.to_numpy()
            for correction, detector, scores in operators(config, root, cell, tf, panels, ix):
                threshold = thresholds.loc[(cell, tf, correction, detector)]
                eligible = bool(threshold.threshold_eligible)
                values = scores if scores is not None else np.full(len(ix), np.nan)
                row, calls, bcalls = rate_summary(values, baseline[ix], totals[ix], valid[correction][ix],
                    comparable[ix], threshold.candidate_threshold, threshold.baseline_threshold, eligible)
                row.update(dict(zip(KEYS, (cell, tf, correction, detector))), **threshold.to_dict(), total_sites=len(ix),
                    unavailable_candidate_sites=int((~(valid[correction][ix] & np.isfinite(values))).sum()),
                    status="missing_model" if scores is None else ("screened" if eligible else "insufficient_calibration"))
                rows.append(row)
                frame = task[safety.SITE_KEYS].copy()
                frame["correction"], frame["detector"] = correction, detector
                frame["candidate_score"], frame["baseline_score"] = values, baseline[ix]
                frame["total_cuts"] = totals[ix]
                frame["candidate_valid"] = valid[correction][ix] & np.isfinite(values)
                frame["baseline_comparable"] = comparable[ix]
                frame["candidate_threshold"], frame["baseline_threshold"] = threshold.candidate_threshold, threshold.baseline_threshold
                frame["candidate_call"], frame["baseline_call"] = calls, bcalls
                frame["threshold_eligible"] = eligible
                site_rows.append(frame)
    if len(rows) != 126:
        raise ValueError("control task denominator changed")
    outputs = {}
    for name, table in (("summary", pd.DataFrame(rows)), ("sites", pd.concat(site_rows, ignore_index=True))):
        p = out / f"screen_{name}.tsv"
        table.to_csv(p, sep="\t", index=False)
        outputs[name] = safety.record(p)
    result = dict(configuration_id=config["configuration_id"], outputs=outputs,
                  thresholds=safety.record(out / "thresholds.freeze.json"), evidence=config["protocol"]["evidence"],
                  promotion_eligible=False, model_refitted=False, new_holdout_labels_read=False)
    result["stage_id"] = canonical_id(result)
    safety.write_immutable(path, result)
    print("Completed all 126 operator ledger rows on already examined replicate 2", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--freeze-only", action="store_true")
    modes.add_argument("--calibrate-only", action="store_true")
    modes.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    out, root = args.outdir.resolve(), args.data_root.resolve()
    if not out.is_relative_to(ROOT / "benchmarks/results"):
        parser.error("use research benchmark results")
    if args.freeze_only:
        freeze(root, out)
    else:
        with builder.extension.exclusive_run(out):
            if args.calibrate_only:
                calibrate(root, out)
            else:
                evaluate(root, out, args.verify_only)
