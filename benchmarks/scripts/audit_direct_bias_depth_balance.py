#!/usr/bin/env python3
"""Sensitivity of archived direct-bias results to current-depth cut balance.

Keep original partners, apply the existing training-scale 0.2-SD caliper to
current log1p cut totals, and freeze all 36 selections before any new scores.
This post-hoc diagnostic preserves the primary study and cannot establish a
causal effect, change a candidate, or open new holdout/control labels.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd

REPOSITORY = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_direct_bias_stability as base
from audit_chrombpnet_comparison import canonical_id, checked_path, data_path
from audit_direct_bias_safety import record, write_immutable
from extend_direct_bias_stability import exclusive_run, read_document, runtime

SCHEMA = "fp-tools-direct-bias-depth-balance-v1"
SECONDARY = "existing_pairs_current_depth_caliper"


def validate_plan(plan):
    expected = {"schema": "fp-tools-direct-bias-depth-balance-plan-v1",
                "evidence": "post_hoc_development_current_depth_balance_sensitivity",
                "cells": ["HepG2", "K562"], "evaluation_chromosomes": base.TEST_CHROMOSOMES,
                "candidate_retuning_allowed": False, "new_holdout_labels_allowed": False,
                "final_untouched_control_allowed": False, "available_matrix": "the_complete_archived_36_run_study",
                "source_cohort": "motif_accessibility_observed_total", "secondary_cohort": SECONDARY,
                "selection": "retain_existing_complete_pairs_if_absolute_current_log1p_cut_total_difference_is_within_frozen_training_caliper",
                "caliper_standard_deviations": 0.2, "scales": "unchanged_original_chr1_15_per_cell_and_tf",
                "maximum_absolute_smd": 0.1, "zero_cut_pairs": "retain_if_the_same_symmetric_caliper_passes",
                "model_scores_consulted_for_selection": False, "rematching_allowed": False,
                "all_run_selections_frozen_before_scoring": True, "all_source_tfs_reported": True,
                "primary_results_recomputed": False,
                "candidate_and_metric_settings": "unchanged_from_the_36_run_source_configuration"}
    if any(plan.get(key) != value for key, value in expected.items()):
        raise ValueError("only the fixed post-hoc current-depth balance diagnostic is allowed")


def qualify_pairs(sites, matched):
    """Reuse complete original partners on the archived all-run DWM support."""
    if matched.duplicated(base.SITE_KEYS).any() or matched.groupby("pair_id").size().ne(2).any():
        raise ValueError("original matching must contain unique sites and complete pairs")
    present = pd.MultiIndex.from_frame(matched[base.SITE_KEYS]).isin(pd.MultiIndex.from_frame(sites[base.SITE_KEYS]))
    pair_present = pd.Series(present, index=matched.index).groupby(matched.pair_id).transform("all")
    selected = matched.loc[pair_present].copy()
    selected["site_index"] = base.align_sites(sites, selected)
    for _, pair in selected.groupby("pair_id"):
        if (set(pair.chip_label) != {0, 1} or pair.tf.nunique() != 1 or pair.TFBS_chr.nunique() != 1
                or pair.cell.nunique() != 1):
            raise ValueError("each original pair needs opposite labels in one TF and chromosome")
    return selected


def select_pairs(sites, pairs, observed, scales, caliper=0.2):
    if caliper != 0.2 or observed.shape != (len(sites), 201) or not np.isfinite(observed).all() or (observed < 0).any():
        raise ValueError("selection requires the fixed caliper and finite nonnegative archived cut profiles")
    totals = np.log1p(observed.sum(axis=1))
    rows = []
    for pair_id, pair in pairs.groupby("pair_id", sort=True):
        tf = pair.tf.iloc[0]
        sd = float(scales[tf]["log_observed_total"]["sd"])
        if not np.isfinite(sd) or sd < 0:
            raise ValueError("training covariate scale is invalid")
        positive = int(pair.loc[pair.chip_label.eq(1), "site_index"].iloc[0])
        negative = int(pair.loc[pair.chip_label.eq(0), "site_index"].iloc[0])
        gap = abs(float(totals[positive] - totals[negative]))
        standardized_gap = gap / sd if sd else (0. if gap <= 1e-12 else float("inf"))
        rows.append({"tf": tf, "pair_id": pair_id, "chromosome": pair.TFBS_chr.iloc[0],
                     "positive_index": positive, "negative_index": negative,
                     "positive_log_total": totals[positive], "negative_log_total": totals[negative],
                     "training_sd": sd, "standardized_gap": standardized_gap,
                     "selected": standardized_gap <= caliper + 1e-12})
    return pd.DataFrame(rows, columns=["tf", "pair_id", "chromosome", "positive_index", "negative_index",
                                       "positive_log_total", "negative_log_total", "training_sd", "standardized_gap", "selected"])


def member_indexes(pairs, selected):
    return pairs.loc[pairs.pair_id.isin(selected.loc[selected.selected, "pair_id"]), "site_index"].to_numpy(dtype=int)


def retention_rows(sites, pairs, selected, observed, minimum):
    totals = observed.sum(axis=1)
    rows = []
    for tf, task in sites.groupby("tf", sort=True):
        before = pairs.loc[pairs.tf.eq(tf), "site_index"].to_numpy(dtype=int)
        keep_pairs = selected.loc[selected.tf.eq(tf)]
        after = member_indexes(pairs, keep_pairs)
        row = {"tf": tf, "original_pairs": len(before) // 2, "retained_pairs": len(after) // 2,
               "status": "eligible" if len(after) // 2 >= minimum else "underpowered"}
        for label, indexes in (("original", before), ("retained", after)):
            subset = sites.iloc[indexes]
            positive = subset.chip_label.to_numpy(dtype=int) == 1
            values = {"observed_total": np.log1p(totals[indexes]), "motif_score": subset.motif_score.to_numpy(),
                      "reference_log_accessibility": np.log1p(subset.accessibility.to_numpy())}
            assessable = min(int(positive.sum()), int((~positive).sum())) >= 2
            for name, data in values.items():
                row[f"{name}_smd_{label}"] = (base.standardized_difference(data[positive], data[~positive])
                                               if assessable else float("nan"))
            row[f"zero_cut_sites_{label}"] = int((totals[indexes] == 0).sum())
        balance = [row[name + "_smd_retained"] for name in ("observed_total", "motif_score", "reference_log_accessibility")]
        row["all_measured_covariates_balanced"] = bool(np.isfinite(balance).all() and max(map(abs, balance)) <= .1)
        rows.append(row)
    return pd.DataFrame(rows)


def load_cell(cell, source, result, matched, root):
    support = base.verify_stage(checked_path(result["supports"][cell], root), source["configuration_id"], root)
    coverage = pd.read_csv(checked_path(support["outputs"]["coverage"], root), sep="\t")
    present = coverage.all_run_dwm_present.to_numpy(dtype=bool)
    reference = source["references"][cell]["profile_manifest"]
    original = pd.read_csv(data_path(reference["sites"], root), sep="\t")
    sites = original.iloc[base.align_sites(original, coverage.loc[present])].reset_index(drop=True)
    if set(sites.cell) != {cell} or not set(sites.TFBS_chr).issubset(base.TEST_CHROMOSOMES):
        raise ValueError("source sites are outside the examined development chromosomes")
    strict = matched.loc[matched.cell.eq(cell) & matched.policy.eq("motif_accessibility_observed_total")]
    pairs = qualify_pairs(sites, strict)
    with np.load(checked_path(support["outputs"]["expected_tracks"], root), allow_pickle=False) as data:
        if not np.array_equal(data["site_hash"], base.site_hashes(coverage)):
            raise ValueError("source DWM profiles do not match frozen coverage coordinates")
        fixed_dwm = data["fixed_dwm"][present]
    return sites, pairs, fixed_dwm


def archived_arrays(stage, sites, root, *, counts_only=False):
    with np.load(checked_path(stage["outputs"]["profiles"], root), allow_pickle=False) as data:
        if not np.array_equal(data["site_hash"], base.site_hashes(sites)):
            raise ValueError("archived counts have different site identities")
        if counts_only:
            return data["observed"]
        return data["observed"], data["expected"], data["sample_dwm"]


def freeze_configuration(plan, plan_path, root, outdir):
    source_path = checked_path(plan["source_configuration"], root)
    source = read_document(source_path, "configuration_id")
    result_path = checked_path(plan["source_result"], root)
    result = base.verify_result(result_path.parent, root)
    if result["configuration_id"] != source["configuration_id"] or runtime() != source["runtime"]:
        raise ValueError("source study or scientific runtime differs from the frozen calculation")
    selection = read_document(checked_path(plan["matched_selection"], root), "selection_id")
    scales_path = checked_path(plan["training_scales"], root)
    scales = json.loads(scales_path.read_text())
    if (selection["outputs"]["scales"]["sha256"] != plan["training_scales"]["sha256"]
            or scales.get("training_chromosomes") != [f"chr{i}" for i in range(1, 16)]
            or scales.get("model_scores_consulted") is not False
            or scales["configuration_id"] != selection["configuration_id"]
            or source["settings"]["matched_cohorts"]["sha256"] != plan["matched_selection"]["sha256"]):
        raise ValueError("matching scales or cohorts differ from the original training freeze")
    entries = source["input_records"] + [source["plan"], *[plan[key] for key in
              ("source_configuration", "source_result", "matched_selection", "training_scales")]]
    entries += list(selection["outputs"].values()) + list(result["outputs"].values())
    for entry in list(result["supports"].values()) + [value for item in result["stages"] for value in item.values()]:
        entries.append(entry)
        entries.extend(json.loads(checked_path(entry, root).read_text())["outputs"].values())
    inputs = {}
    for entry in entries:
        path = str(data_path(entry["path"], root).resolve())
        if path in inputs and inputs[path] != entry["sha256"]:
            raise ValueError("one source path has conflicting frozen checksums")
        inputs[path] = entry["sha256"]
    # base.verify_result verified the source graph. Check additional selection
    # leaves as well; no BAM extraction or completed scientific stage is rerun.
    for entry in selection["outputs"].values():
        checked_path(entry, root)
    code = source["code"] + [record(Path(__file__)), record(Path(__file__).with_name("extend_direct_bias_stability.py"))]
    config = {"schema": SCHEMA, "settings": plan, "plan": record(plan_path), "source_configuration": plan["source_configuration"],
              "source_configuration_id": source["configuration_id"], "source_result": plan["source_result"],
              "scientific_settings": source["settings"], "code": code, "runtime": runtime(),
              "input_records": [{"path": path, "sha256": digest} for path, digest in sorted(inputs.items())]}
    config["configuration_id"] = canonical_id(config)
    write_immutable(outdir / "configuration.freeze.json", config)
    print(f"Froze current-depth balance diagnostic: {config['configuration_id']}", flush=True)
    return config, source, result, selection, scales


def freeze_cohorts(config, source, result, selection, scales, outdir, root):
    matched = pd.read_csv(checked_path(selection["outputs"]["cohorts"], root), sep="\t")
    cells = {cell: load_cell(cell, source, result, matched, root) for cell in config["settings"]["cells"]}
    final = outdir / "cohorts.freeze.json"
    if final.exists():
        done = base.verify_stage(final, config["configuration_id"], root)
        return cells, done
    pair_rows, retention, runs = [], [], []
    for item in result["stages"]:
        stage = base.verify_stage(checked_path(item["extraction"], root), source["configuration_id"], root)
        run = stage["run"]
        sites, pairs, _ = cells[run["cell"]]
        observed = archived_arrays(stage, sites, root, counts_only=True)
        selected = select_pairs(sites, pairs, observed, scales["scales"][run["cell"]])
        pair_rows.append(selected.assign(**run))
        retention.append(retention_rows(sites, pairs, selected, observed, source["settings"]["minimum_sites_per_class"]).assign(**run))
        runs.append(run)
    base.validate_matrix(pd.DataFrame(runs))
    outputs = {}
    for name, frames in (("pairs", pair_rows), ("retention", retention)):
        path = outdir / f"{name}.tsv.gz"
        pd.concat(frames, ignore_index=True).to_csv(path, sep="\t", index=False, compression={"method": "gzip", "mtime": 0})
        outputs[name] = record(path)
    done = {"configuration_id": config["configuration_id"], "outputs": outputs, "runs": runs,
            "model_scores_consulted": False, "new_evaluation_started": False, "partners_changed": False}
    done["stage_id"] = canonical_id(done)
    write_immutable(final, done)
    print("Froze all 36 pair selections before new evaluation.", flush=True)
    return cells, done


def select_run(frame, run):
    selected = np.ones(len(frame), dtype=bool)
    for name, value in run.items():
        selected &= frame[name].eq(value).to_numpy()
    return frame.loc[selected]


def empty_evaluation(run, result, retention, outdir, configuration_id, root):
    """Keep explicit coverage even when no TF has two retained pairs."""
    outputs = {}
    for name in ("metrics", "bootstrap", "profiles", "shapes", "coverage"):
        frame = pd.read_csv(checked_path(result["outputs"][name], root), sep="\t", nrows=0)
        if name == "coverage":
            frame = retention[["tf", "status"]].copy()
            frame["n_positive"] = retention.retained_pairs.to_numpy()
            frame["n_negative"] = retention.retained_pairs.to_numpy()
            frame["zero_cut_sites"] = retention.zero_cut_sites_retained.to_numpy()
            for field in ("observed_total", "motif_score", "reference_log_accessibility"):
                frame[field + "_smd"] = retention[field + "_smd_retained"].to_numpy()
            frame = frame.assign(**run, cohort=SECONDARY)
        path = outdir / f"{name}.tsv.gz"
        frame.to_csv(path, sep="\t", index=False, compression={"method": "gzip", "mtime": 0})
        outputs[name] = record(path)
    done = {"configuration_id": configuration_id, "run": run, "outputs": outputs,
            "bootstrap_run": run["depth"] == "full" or (run["depth"] == "10m" and run["seed"] == 2026),
            "status": "no_task_with_two_retained_pairs"}
    done["stage_id"] = canonical_id(done)
    write_immutable(outdir / "evaluation.json", done)
    return done


def evaluate(config, source, result, cells, cohorts, outdir, root):
    # Reading the complete immutable joint selection precedes any scoring.
    base.verify_stage(outdir / "cohorts.freeze.json", config["configuration_id"], root)
    pair_table = pd.read_csv(checked_path(cohorts["outputs"]["pairs"], root), sep="\t")
    retention = pd.read_csv(checked_path(cohorts["outputs"]["retention"], root), sep="\t")
    completed = []
    for item in result["stages"]:
        stage = read_document(checked_path(item["extraction"], root), "stage_id")
        run = stage["run"]
        row = SimpleNamespace(**run)
        directory = outdir / row.sample / row.depth / f"seed_{row.seed}"
        directory.mkdir(parents=True, exist_ok=True)
        manifest = directory / "evaluation.json"
        if manifest.is_file():
            base.verify_stage(manifest, config["configuration_id"], root, expected_run=run)
        else:
            kept = select_run(retention, run)
            if kept.retained_pairs.max() < 2:
                empty_evaluation(run, result, kept, directory, config["configuration_id"], root)
            else:
                sites, pairs, fixed_dwm = cells[row.cell]
                observed, expected, sample_dwm = archived_arrays(stage, sites, root)
                methods = base.score_panel(observed, expected, fixed_dwm, sample_dwm, source["settings"]["dispersion"][row.cell])
                indexes = member_indexes(pairs, select_run(pair_table, run))
                base.evaluate_run(row, sites, {SECONDARY: indexes}, methods, observed, source["settings"],
                                  directory, config["configuration_id"], root)
        completed.append({"run": run, "source_extraction": item["extraction"], "evaluation": record(manifest)})
    outputs = {"retention": cohorts["outputs"]["retention"]}
    for name in ("metrics", "bootstrap", "profiles", "shapes", "coverage"):
        frames = []
        original = pd.read_csv(checked_path(result["outputs"][name], root), sep="\t")
        original = original.loc[original.cohort.eq(config["settings"]["source_cohort"])].copy()
        original["evidence_role"] = "preserved_primary"
        frames.append(original)
        for item in completed:
            stage = base.verify_stage(checked_path(item["evaluation"], root), config["configuration_id"], root, expected_run=item["run"])
            table = pd.read_csv(checked_path(stage["outputs"][name], root), sep="\t")
            if not table.empty:
                frames.append(table.assign(evidence_role="post_hoc_current_depth_sensitivity"))
        path = outdir / f"balance_{name}.tsv.gz"
        pd.concat(frames, ignore_index=True).to_csv(path, sep="\t", index=False, compression={"method": "gzip", "mtime": 0})
        outputs[name] = record(path)
    done = {"schema": SCHEMA, "configuration_id": config["configuration_id"], "outputs": outputs,
            "cohorts": record(outdir / "cohorts.freeze.json"), "stages": completed, "primary_results_recomputed": False,
            "candidate_retuned": False, "new_holdout_labels_read": False, "final_untouched_control_read": False,
            "promotion_eligible": False, "causal_effect_established": False, "evidence": config["settings"]["evidence"]}
    done["stage_id"] = canonical_id(done)
    write_immutable(outdir / "balance_manifest.json", done)
    return done


def verify_result(outdir, root, code_ref=None):
    config = read_document(outdir / "configuration.freeze.json", "configuration_id")
    validate_plan(config["settings"])
    if code_ref is not None and not re.fullmatch(r"[0-9a-fA-F]{7,40}", code_ref):
        raise ValueError("--code-ref must be an immutable Git commit hash")
    for entry in config["code"]:
        if code_ref is None:
            checked_path(entry, root)
        else:
            name = Path(entry["path"]).relative_to(REPOSITORY).as_posix()
            if sha256(subprocess.check_output(["git", "show", f"{code_ref}:{name}"], cwd=REPOSITORY)).hexdigest() != entry["sha256"]:
                raise ValueError("historical diagnostic code checksum mismatch")
    for entry in config["input_records"] + [config["plan"]]:
        checked_path(entry, root)
    source = read_document(checked_path(config["source_configuration"], root), "configuration_id")
    if config["scientific_settings"] != source["settings"] or config["source_configuration_id"] != source["configuration_id"]:
        raise ValueError("depth-balance diagnostic changed the archived scientific settings")
    result = base.verify_stage(outdir / "balance_manifest.json", config["configuration_id"], root)
    if (result.get("schema") != SCHEMA or result.get("evidence") != config["settings"]["evidence"]
            or any(result.get(name) is not False for name in ("promotion_eligible", "candidate_retuned", "new_holdout_labels_read",
                   "final_untouched_control_read", "primary_results_recomputed", "causal_effect_established"))):
        raise ValueError("unexpected depth-balance evidence scope")
    cohorts = base.verify_stage(checked_path(result["cohorts"], root), config["configuration_id"], root)
    if any(cohorts.get(name) is not False for name in ("model_scores_consulted", "new_evaluation_started", "partners_changed")):
        raise ValueError("cohort freeze did not preserve the fixed selection protocol")
    base.validate_matrix(pd.DataFrame(cohorts["runs"]))
    runs = []
    for item in result["stages"]:
        run = item["run"]
        base.verify_stage(checked_path(item["source_extraction"], root), config["source_configuration_id"], root, expected_run=run)
        base.verify_stage(checked_path(item["evaluation"], root), config["configuration_id"], root, expected_run=run)
        runs.append(run)
    base.validate_matrix(pd.DataFrame(runs))
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--freeze-only", action="store_true")
    parser.add_argument("--verify-only", action="store_true")
    parser.add_argument("--code-ref")
    args = parser.parse_args(argv)
    if args.code_ref is not None and not args.verify_only:
        parser.error("--code-ref requires --verify-only")
    plan = json.loads(args.plan.read_text())
    validate_plan(plan)
    outdir = args.outdir.resolve()
    if args.verify_only or (outdir / "balance_manifest.json").is_file():
        config = read_document(outdir / "configuration.freeze.json", "configuration_id")
        if record(args.plan)["sha256"] != config["plan"]["sha256"]:
            raise ValueError("requested plan differs from completed depth-balance diagnostic")
        verify_result(outdir, args.data_root, args.code_ref)
        print("Verified the complete 36-run sensitivity and preserved primary evidence; no stage rerun.")
        return 0
    with exclusive_run(outdir):
        config, source, result, selection, scales = freeze_configuration(plan, args.plan, args.data_root, outdir)
        cells, cohorts = freeze_cohorts(config, source, result, selection, scales, outdir, args.data_root)
        if not args.freeze_only:
            done = evaluate(config, source, result, cells, cohorts, outdir, args.data_root)
            print(f"Completed current-depth balance diagnostic: {done['stage_id']}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
