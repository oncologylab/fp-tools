#!/usr/bin/env python3
"""Compare corrected control-2 safety with fixed and recalibrated thresholds."""

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_direct_bias_safety as safety
import audit_direct_bias_safety_coverage as coverage
import rebuild_corrected_safety_profiles as builder
from audit_chrombpnet_comparison import canonical_id, checked_path, check_identity

ANALYSES = ("fixed_threshold_extraction_effect", "corrected_calibration")


def read_profiles(entry, root):
    doc = json.loads(checked_path(entry, root).read_text())
    check_identity(doc, "stage_id")
    sites = pd.read_csv(checked_path(doc["outputs"]["sites"], root), sep="\t")
    with np.load(checked_path(doc["outputs"]["profiles"], root), allow_pickle=False) as data:
        arrays = {k: data[k] for k in data.files}
    if not np.array_equal(arrays["site_hash"], safety.site_hashes(sites)):
        raise ValueError("corrected profile site order changed")
    return sites, arrays


def calibrate_corrected(sites, arrays, dwm, dwm_valid, dispersion, protocol):
    valid = arrays["valid"] & dwm_valid
    totals = arrays["observed"].sum(axis=1)
    frames = []
    for method, expected in zip(safety.METHODS, (arrays["expected"], dwm)):
        scores = safety.direct_scores(arrays["observed"], expected, dispersion)[0]
        frame = safety.calibrate(sites, scores, valid, totals,
            target=protocol["target_validation_negative_call_rate"],
            minimum_negatives=protocol["minimum_validation_negatives"])
        frame["method"] = method
        missing = {}
        for tf, task in sites.groupby("tf"):
            ix = task.index.to_numpy()
            member = task.TFBS_chr.isin(protocol["calibration_chromosomes"]).to_numpy() & task.chip_label.eq(0).to_numpy()
            eligible = member & valid[ix] & (totals[ix] > 0) & np.isfinite(scores[ix])
            missing[tf] = int((eligible & (dwm[ix].sum(axis=1) <= 0)).sum())
        frame["zero_dwm_expectation_calibration_negatives"] = frame.tf.map(missing)
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)


def qualify(candidate, baseline, totals, common, disagreement, cutoff, *, residual_flag, calibration_qualified):
    row = coverage.qualify_row(candidate, baseline, totals, common, disagreement,
        candidate_threshold=cutoff[0].threshold, baseline_threshold=cutoff[1].threshold,
        threshold_eligible=all(bool(c.threshold_eligible) for c in cutoff), motif_residual_flag=residual_flag)
    row["baseline_calibration_qualified"] = bool(calibration_qualified)
    row["passes_candidate_screen"] &= bool(calibration_qualified)
    return row


def freeze(root, profiles, output):
    builder.execute(root, profiles, verify_only=True)
    if output.exists():
        raise ValueError("use a fresh safety output directory")
    protocol = json.loads(checked_path({"path": str(builder.PROTOCOL), "sha256": builder.PROTOCOL_SHA}, ROOT).read_text())
    source = json.loads(checked_path(protocol["source_configuration"], ROOT).read_text())
    inputs = [safety.record(builder.PROTOCOL)]
    for key in ("source_configuration", "source_thresholds", "source_result"):
        inputs.append(safety.record(checked_path(protocol[key], ROOT)))
    threshold_doc = json.loads(checked_path(protocol["source_thresholds"], ROOT).read_text())
    inputs.append(safety.record(checked_path(threshold_doc["thresholds"], root)))
    stages = {}
    for cell in protocol["cells"]:
        stages[cell] = {}
        for role in ("validation", "naked"):
            p = profiles / f"{cell}.{role}.json"
            stages[cell][role] = safety.record(p)
            doc = json.loads(p.read_text())
            inputs.append(safety.record(p))
            inputs.extend(safety.record(checked_path(e, root)) for e in doc["outputs"].values())
        for role in ("naked", "naked_dwm"):
            p = checked_path(source["cells"][cell][role], root)
            doc = json.loads(p.read_text())
            inputs.append(safety.record(p))
            for name, digest in (("sites", "sites_sha256"), ("profiles_npz", "profiles_sha256")):
                inputs.append(safety.record(checked_path({"path": doc[name], "sha256": doc[digest]}, root)))
        inputs.append(safety.record(checked_path(source["cells"][cell]["dwm_track"], root)))
    build_config = json.loads((profiles / "configuration.freeze.json").read_text())
    paths = {Path(e["path"]) for e in build_config["code"]}
    paths.update([Path(__file__), Path(coverage.__file__), Path(safety.__file__)])
    config = dict(schema="fp-tools-corrected-safety-reassessment-v1", protocol=protocol,
                  corrected_stages=stages, inputs=inputs, source_configuration=source,
                  profile_configuration=safety.record(profiles / "configuration.freeze.json"),
                  profile_manifest=safety.record(profiles / "profiles_manifest.json"),
                  code=[safety.record(p) for p in sorted(paths)], runtime=builder.extension.runtime(),
                  promotion_eligible=False, final_untouched_control_read=False, candidate_retuned=False)
    config["configuration_id"] = canonical_id(config)
    output.mkdir(parents=True)
    safety.write_immutable(output / "configuration.freeze.json", config)
    print("Frozen corrected safety comparison inputs", flush=True)


def preflight(root, output):
    config = builder.extension.read_document(output / "configuration.freeze.json", "configuration_id")
    if (config["schema"] != "fp-tools-corrected-safety-reassessment-v1" or config["runtime"] != builder.extension.runtime()
            or any(config[k] is not False for k in ("promotion_eligible", "final_untouched_control_read", "candidate_retuned"))):
        raise ValueError("safety reassessment scope or runtime changed")
    protocol = json.loads(checked_path({"path": str(builder.PROTOCOL), "sha256": builder.PROTOCOL_SHA}, ROOT).read_text())
    if config["protocol"] != protocol:
        raise ValueError("safety protocol changed")
    for item in config["inputs"] + config["code"] + [config["profile_configuration"], config["profile_manifest"]]:
        checked_path(item, root)
    return config


def calibrate(root, output):
    config = preflight(root, output)
    source, protocol = config["source_configuration"], config["protocol"]
    frames = []
    for cell, stages in config["corrected_stages"].items():
        sites, arrays = read_profiles(stages["validation"], root)
        dwm, valid = safety.extract_prediction_profiles(sites, checked_path(source["cells"][cell]["dwm_track"], root), 100, require_dense=False)
        frame = calibrate_corrected(sites, arrays, np.nan_to_num(dwm), valid, source["cells"][cell]["dispersion"], protocol)
        frame["cell"] = cell
        frames.append(frame)
    corrected = pd.concat(frames, ignore_index=True)
    _, old = safety.read_thresholds(checked_path(protocol["source_thresholds"], ROOT), source["configuration_id"], root)
    keys = ["cell", "tf", "method"]
    if set(map(tuple, old[keys].values)) != set(map(tuple, corrected[keys].values)):
        raise ValueError("calibration task denominator changed")
    fixed = old.merge(corrected[keys + ["zero_dwm_expectation_calibration_negatives"]], on=keys, validate="one_to_one")
    fixed["analysis"], corrected["analysis"] = ANALYSES
    table = pd.concat([fixed, corrected], ignore_index=True)
    path = output / "thresholds.tsv"
    receipt = output / "thresholds.freeze.json"
    if receipt.exists():
        doc = builder.base.verify_stage(receipt, config["configuration_id"], root)
        pd.testing.assert_frame_equal(pd.read_csv(checked_path(doc["outputs"]["thresholds"], root), sep="\t"), table, check_dtype=False)
        return
    if path.exists():
        raise ValueError("partial thresholds preserved")
    table.to_csv(path, sep="\t", index=False)
    doc = dict(configuration_id=config["configuration_id"], outputs={"thresholds": safety.record(path)}, control_profiles_scored=False)
    doc["stage_id"] = canonical_id(doc)
    safety.write_immutable(receipt, doc)
    print("Both threshold analyses frozen before control scoring", flush=True)


def evaluate(root, output, verify_only=False):
    config = preflight(root, output)
    receipt = builder.base.verify_stage(output / "thresholds.freeze.json", config["configuration_id"], root)
    if receipt.get("control_profiles_scored") is not False:
        raise ValueError("thresholds were not frozen before control scoring")
    result_path = output / "safety_manifest.json"
    if result_path.exists():
        result = builder.base.verify_stage(result_path, config["configuration_id"], root)
        checked_path(result["thresholds"], root)
        if (result["thresholds"] != safety.record(output / "thresholds.freeze.json")
                or result.get("evidence") != "post_hoc_corrected_implementation_candidate_screen"
                or any(result.get(k) is not False for k in ("promotion_eligible", "final_untouched_control_read", "candidate_retuned"))):
            raise ValueError("completed safety scope or thresholds changed")
        print("Verified completed corrected safety reassessment", flush=True)
        return
    if verify_only:
        raise ValueError("missing safety result")
    if list(output.glob("corrected_*.tsv")):
        raise ValueError("partial scoring output preserved")
    thresholds = pd.read_csv(checked_path(receipt["outputs"]["thresholds"], root), sep="\t").set_index(["analysis", "cell", "tf", "method"])
    rows, sites_out, curves, residuals = [], [], [], []
    source, protocol = config["source_configuration"], config["protocol"]
    for cell, stages in config["corrected_stages"].items():
        sites, arrays = read_profiles(stages["naked"], root)
        old_sites, old_observed, _, old_valid, _ = safety.load_artifact(source["cells"][cell]["naked"], root, cell, control=True)
        dwm_sites, dwm_observed, dwm_expected, dwm_valid, _ = safety.load_artifact(source["cells"][cell]["naked_dwm"], root, cell, control=True)
        if not np.array_equal(safety.site_hashes(old_sites), safety.site_hashes(sites)):
            raise ValueError("historical control support changed")
        order = pd.Index(safety.site_hashes(dwm_sites)).get_indexer(safety.site_hashes(sites))
        if np.any(order < 0):
            raise ValueError("DWM control support missing")
        dwm_observed, dwm_expected, dwm_valid = dwm_observed[order], dwm_expected[order], dwm_valid[order]
        observed, valid = arrays["observed"], arrays["valid"]
        common, disagreement = coverage.comparable_profiles(observed, dwm_observed, dwm_expected)
        old_common, _ = coverage.comparable_profiles(old_observed, dwm_observed, dwm_expected)
        common &= valid & dwm_valid
        old_common &= old_valid & dwm_valid
        totals = observed.sum(axis=1)
        methods = {name: safety.direct_scores(observed, expected, source["cells"][cell]["dispersion"])
                   for name, expected in zip(safety.METHODS, (arrays["expected"], dwm_expected))}
        site_table = sites.copy()
        site_table["valid"], site_table["archived_valid"] = valid, old_valid
        site_table["dwm_comparable"], site_table["archived_dwm_comparable"] = common, old_common
        site_table["total_cuts"], site_table["archived_total_cuts"] = totals, old_observed.sum(axis=1)
        for name, values in methods.items():
            site_table[name + "_score"] = values[0]
        sites_out.append(site_table)
        for tf, task in sites.groupby("tf"):
            all_ix = task.index.to_numpy()
            ix = all_ix[valid[all_ix]]
            informative = ix[totals[ix] > 0]
            residual_flag = True
            if len(informative) >= 2:
                effect = safety.motif_residual_effect(observed[informative], methods["parametric_direct"][1][informative], np.arange(-100, 101),
                    bootstraps=protocol["residual_bootstraps"], seed=protocol["seed"], threshold=protocol["motif_residual_absolute_center_flank_limit"])
                residual_flag = bool(effect["motif_residual_flag"])
                residuals.append(dict(cell=cell, tf=tf, informative_sites=len(informative), **effect))
                for method, values in methods.items():
                    curves.append(pd.DataFrame(dict(cell=cell, tf=tf, method=method, position=np.arange(-100, 101),
                        mean_observed=observed[informative].mean(axis=0), mean_expected=values[1][informative].mean(axis=0), mean_deviance=values[2][informative].mean(axis=0))))
            for analysis in ANALYSES:
                cutoff = [thresholds.loc[(analysis, cell, tf, m)] for m in safety.METHODS]
                row = qualify(methods["parametric_direct"][0][ix], methods["DWM_expected"][0][ix], totals[ix], common[ix], disagreement[ix], cutoff,
                    residual_flag=residual_flag, calibration_qualified=all(c.zero_dwm_expectation_calibration_negatives == 0 for c in cutoff))
                row.update(cell=cell, tf=tf, analysis=analysis, total_sites=len(all_ix), invalid_sites=len(all_ix)-len(ix),
                    decision_role="diagnostic_only" if analysis == ANALYSES[0] else "primary_development_candidate_screen",
                    archived_comparable=int(old_common[all_ix].sum()), comparable_gained=int((common & ~old_common)[all_ix].sum()),
                    comparable_lost=int((old_common & ~common)[all_ix].sum()), changed_count_sites=int(np.any(observed[all_ix] != old_observed[all_ix], axis=1).sum()))
                rows.append(row)
    outputs = {}
    for name, table in (("safety", pd.DataFrame(rows)), ("sites", pd.concat(sites_out)),
                        ("residuals", pd.DataFrame(residuals)), ("curves", pd.concat(curves) if curves else pd.DataFrame())):
        path = output / f"corrected_{name}.tsv"
        table.to_csv(path, sep="\t", index=False)
        outputs[name] = safety.record(path)
    result = dict(configuration_id=config["configuration_id"], outputs=outputs, thresholds=safety.record(output / "thresholds.freeze.json"),
                  evidence="post_hoc_corrected_implementation_candidate_screen", promotion_eligible=False,
                  final_untouched_control_read=False, candidate_retuned=False)
    result["stage_id"] = canonical_id(result)
    safety.write_immutable(result_path, result)
    print("Completed corrected replicate-2 safety reassessment", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--profiles", type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--freeze-only", action="store_true")
    modes.add_argument("--calibrate-only", action="store_true")
    modes.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    output = args.outdir.resolve()
    if not output.is_relative_to(ROOT / "benchmarks/results"):
        parser.error("output must stay in research benchmark results")
    if args.freeze_only:
        if args.profiles is None:
            parser.error("freeze requires --profiles")
        freeze(args.data_root.resolve(), args.profiles.resolve(), output)
    else:
        with builder.extension.exclusive_run(output):
            if args.calibrate_only:
                calibrate(args.data_root.resolve(), output)
            else:
                evaluate(args.data_root.resolve(), output, args.verify_only)
