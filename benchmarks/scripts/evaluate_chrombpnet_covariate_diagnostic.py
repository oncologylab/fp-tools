#!/usr/bin/env python3
"""Frozen covariate-matching sensitivity for previously examined comparisons.

The candidate, scoring geometry and per-cell dispersion remain unchanged.
Fit covariate scales on chr1-15, freeze within-chromosome matched cohorts
without consulting model scores, then summarize all prespecified policies.
This diagnostic cannot create new holdout evidence or authorize promotion.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
from importlib.metadata import version
import json
from pathlib import Path
import re
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))
sys.path.insert(0, str(SCRIPT_DIR))

from audit_chrombpnet_comparison import (
    DWM, METHODS, PARAMETRIC, PINNED_COMMIT, SITE_KEYS, canonical_id, check_identity,
    checked_path, shape_audit, standardized_difference,
)
from audit_direct_bias_safety import one_record, record, write_immutable
from build_strand_functional_profiles import site_hashes
from evaluate_chrombpnet_reference import extract_prediction_profiles, metric_row, score_methods
from evaluate_frozen_functional_policy import aggregate_curve
from evaluate_parametric_factorization import block_bootstrap_delta

SCHEMA = "fp-tools-chrombpnet-covariate-diagnostic-v1"
FEATURES = ("motif_score", "log_accessibility", "log_observed_total")
TRAIN_CHROMOSOMES = [f"chr{i}" for i in range(1, 16)]
TEST_CHROMOSOMES = ["chr19", "chr20", "chr21", "chr22", "chrX"]


def fit_scales(frame, features, chromosomes):
    train = frame.loc[frame.TFBS_chr.isin(chromosomes), list(features)].to_numpy(dtype=float)
    train = train[np.isfinite(train).all(axis=1)]
    if len(train) < 2:
        raise ValueError("at least two finite training sites are required for covariate scales")
    return {"training_sites": len(train), **{feature: {"mean": float(train[:, i].mean()),
                                                     "sd": float(train[:, i].std(ddof=1))}
                                            for i, feature in enumerate(features)}}


def match_sites(frame, features, scales, *, caliper):
    """Find the most valid 1:1 pairs, then minimize distance; ignore score columns."""
    if caliper <= 0 or not np.isfinite(caliper) or frame.duplicated(SITE_KEYS).any():
        raise ValueError("matching requires a positive caliper and unique loci")
    work = frame.loc[frame.chip_label.isin([0, 1])].copy()
    z = work[list(features)].to_numpy(dtype=float, copy=True)
    for i, feature in enumerate(features):
        center, sd = scales[feature]["mean"], scales[feature]["sd"]
        if sd == 0:
            z[:, i] = np.where(np.isclose(z[:, i], center, rtol=0, atol=1e-12), 0, np.nan)
        else:
            z[:, i] = (z[:, i] - center) / sd
    work[["_z_" + name for name in features]] = z
    work = work.loc[np.isfinite(z).all(axis=1)]
    rows = []
    for chromosome, group in work.groupby("TFBS_chr", sort=True):
        groups = [group.loc[group.chip_label.eq(label)].sort_values(SITE_KEYS).reset_index(drop=True)
                  for label in (1, 0)]
        if any(g.empty for g in groups):
            continue
        values = [g[["_z_" + name for name in features]].to_numpy(dtype=float) for g in groups]
        delta = np.abs(values[0][:, None, :] - values[1][None, :, :])
        allowed = delta.max(axis=2) <= caliper + 1e-12
        # Each allowed cost <= 1. One forbidden edge costs more than all allowed
        # edges combined, so cardinality is optimized before squared distance.
        distance = np.mean((delta / caliper) ** 2, axis=2)
        cost = np.where(allowed, distance, min(len(g) for g in groups) + 1.0)
        positives, negatives = linear_sum_assignment(cost)
        for positive, negative in zip(positives, negatives):
            if not allowed[positive, negative]:
                continue
            pair = pd.concat([groups[0].iloc[[positive]], groups[1].iloc[[negative]]], ignore_index=True)
            identity = "|".join(str(value) for value in pair[SITE_KEYS].to_numpy().ravel())
            pair["pair_id"] = sha256(identity.encode()).hexdigest()[:20]
            pair["match_maximum_standardized_difference"] = float(delta[positive, negative].max())
            rows.append(pair)
    columns = list(frame.columns) + ["pair_id", "match_maximum_standardized_difference"]
    if not rows:
        return pd.DataFrame(columns=columns)
    return pd.concat(rows, ignore_index=True)[columns].sort_values(
        ["TFBS_chr", "pair_id", "chip_label"], kind="mergesort").reset_index(drop=True)


def add_covariates(sites, observed, *, accessibility=None):
    result = sites.copy()
    if accessibility is not None:
        if "accessibility" in result:
            if not np.allclose(result.accessibility, accessibility, rtol=1e-6, atol=1e-6):
                raise ValueError("archived accessibility differs from the recorded raw track")
        else:
            result["accessibility"] = accessibility
    if "accessibility" not in result:
        raise ValueError("accessibility must come from the recorded raw track or archived covariate")
    if len(result) != len(observed) or np.any(observed < 0) or np.any(result.accessibility < 0):
        raise ValueError("covariates require aligned nonnegative observed cuts and accessibility")
    result["log_accessibility"] = np.log1p(result.accessibility.to_numpy(dtype=float))
    result["observed_total"] = observed.sum(axis=1)
    result["log_observed_total"] = np.log1p(result.observed_total.to_numpy(dtype=float))
    return result


def extract_accessibility(sites, track, root):
    values, valid = extract_prediction_profiles(sites, checked_path(track, root), 100, require_dense=False)
    if not valid.all():
        raise ValueError("accessibility track lacks a recorded development chromosome/window")
    # Match the float32 cache used by the original accessibility-matching helper.
    return np.abs(values.astype(np.float32)).sum(axis=1)


def read_profiles(item, root, cell, *, training=False):
    document = json.loads(checked_path(item, root).read_text())
    metadata = document["metadata"]
    if (document.get("schema") != "fp-tools-strand-functional-profiles-v1"
            or metadata.get("labels_used") is not False or metadata.get("flank") != 100
            or metadata.get("read_shift") != [4, -5]):
        raise ValueError("unexpected development profile construction")
    sites = pd.read_csv(checked_path({"path": document["sites"], "sha256": document["sites_sha256"]}, root), sep="\t")
    allowed = TRAIN_CHROMOSOMES + ["chr16", "chr17", "chr18"] if training else TEST_CHROMOSOMES
    if (set(sites.cell) != {cell} or sites.duplicated(SITE_KEYS).any()
            or not set(sites.TFBS_chr).issubset(allowed) or not sites.chip_label.isin([0, 1]).all()):
        raise ValueError("profile sites are outside the recorded development scope")
    array_path = checked_path({"path": document["profiles_npz"], "sha256": document["profiles_sha256"]}, root)
    with np.load(array_path, allow_pickle=False) as arrays:
        if not np.array_equal(site_hashes(sites), arrays["site_hash"]):
            raise ValueError("profile coordinate order mismatch")
        observed = arrays["plus_observed"] + arrays["minus_observed"]
        expected = None if training else arrays["plus_expected"] + arrays["minus_expected"]
        valid = arrays["valid"].astype(bool)
    return sites, observed, expected, valid, metadata


def resolve_plan(plan, root):
    if (plan.get("schema") != "fp-tools-chrombpnet-covariate-diagnostic-plan-v1"
            or sorted(plan.get("cells", [])) != ["HepG2", "K562"]
            or plan.get("training_chromosomes") != TRAIN_CHROMOSOMES
            or plan.get("evaluation_chromosomes") != TEST_CHROMOSOMES
            or plan.get("new_holdout_labels_allowed") is not False
            or plan.get("candidate_retuning_allowed") is not False):
        raise ValueError("diagnostic accepts only the fixed, previously examined development scope")
    if (plan["primary_policy"] not in plan["policies"]
            or any(not names or not set(names).issubset(FEATURES) for names in plan["policies"].values())
            or set(plan.get("accessibility_tracks", {})) != set(plan["cells"])
            or plan["minimum_sites_per_class"] < 200 or plan["bootstrap_iterations"] < 1000):
        raise ValueError("invalid fixed matching policies or power/bootstrap requirements")
    source = json.loads(checked_path(plan["source_plan"], root).read_text())
    registry = json.loads(checked_path(source["validation_input_registry"], root).read_text())
    cells = {}
    for item in source["comparisons"]:
        cell = item["cell"]
        if cell not in plan["cells"] or cell in cells:
            raise ValueError("comparison cell differs from development plan")
        comparison_path = checked_path(item, root)
        comparison = json.loads(comparison_path.read_text())
        freeze_path = comparison_path.with_name("chrombpnet_comparison_inputs.freeze.json")
        freeze = json.loads(freeze_path.read_text())
        check_identity(comparison, "result_id")
        check_identity(freeze, "input_id")
        if (freeze["input_id"] != comparison["input_id"] or freeze["cells"] != [cell]
                or freeze["pinned_chrombpnet_commit"] != PINNED_COMMIT
                or freeze["common_finite_support"] is not True
                or freeze["test_labels_previously_opened_for_other_models"] is not True):
            raise ValueError("comparison is not the recorded common-support development experiment")
        for entry in freeze["inputs"]:
            checked_path(entry, root)
        checked_path(comparison["outputs"]["site_scores"], root)
        checked_path(plan["accessibility_tracks"][cell], root)
        cells[cell] = {"comparison": item, "comparison_result_id": comparison["result_id"],
                       "freeze": record(freeze_path), "scores": comparison["outputs"]["site_scores"],
                       "test": one_record(freeze["inputs"], lambda r: r["purpose"] == f"{cell}-test-manifest"),
                       "training": one_record(registry["validation_inputs"], lambda r: Path(r["path"]).name == f"{cell}.json"),
                       "accessibility": plan["accessibility_tracks"][cell],
                       "tracks": {name: one_record(freeze["inputs"], lambda r: r["purpose"] == f"{cell}-{name}")
                                  for name in ("dwm", "bias", "regulatory", "nobias")},
                       "dispersion": freeze["dispersion"]}
    if set(cells) != set(plan["cells"]):
        raise ValueError("source plan lacks a development cell")
    return cells


def matched_diagnostics(before, after, features, minimum, maximum_smd):
    result = {"original_positive": int(before.chip_label.eq(1).sum()),
              "original_negative": int(before.chip_label.eq(0).sum()), "matched_pairs": len(after) // 2}
    for name in FEATURES:
        for suffix, frame in (("before", before), ("after", after)):
            result[f"{name}_smd_{suffix}"] = standardized_difference(
                frame.loc[frame.chip_label.eq(1), name], frame.loc[frame.chip_label.eq(0), name])
    values = [abs(result[f"{name}_smd_after"]) for name in features]
    result["maximum_matched_covariate_smd"] = max(values) if np.isfinite(values).all() else np.nan
    result["balance_pass"] = bool(np.isfinite(result["maximum_matched_covariate_smd"])
                                  and result["maximum_matched_covariate_smd"] <= maximum_smd)
    result["power_pass"] = result["matched_pairs"] >= minimum
    result["status"] = ("underpowered" if not result["power_pass"] else
                        "eligible_balanced" if result["balance_pass"] else "balance_failed")
    return result


def verify_outputs(path, identity, root):
    document = json.loads(path.read_text())
    check_identity(document, identity)
    for entry in document["outputs"].values():
        target = checked_path(entry, root)
        if target.suffix == ".json":
            child = json.loads(target.read_text())
            for child_identity in ("selection_id", "cache_id", "evaluation_id"):
                if child_identity in child and "outputs" in child:
                    verify_outputs(target, child_identity, root)
                    break
    return document


def verify_historical_result(result_path, root, *, code_ref=None):
    """Verify frozen inputs/outputs without rebuilding a configuration or scoring.

    An explicit immutable Git revision can supply calculation code after a
    verification-only code update. Every blob must match its recorded SHA-256.
    """
    if code_ref is not None and not re.fullmatch(r"[0-9a-fA-F]{7,40}", code_ref):
        raise ValueError("--code-ref must be an immutable Git commit hash")
    done = json.loads(result_path.read_text())
    if (done.get("schema") != SCHEMA or done.get("new_holdout_labels_read") is not False
            or done.get("candidate_retuned") is not False):
        raise ValueError("verification accepts only the frozen development diagnostic")
    config = json.loads(result_path.with_name("configuration.freeze.json").read_text())
    check_identity(config, "configuration_id")
    if (config["configuration_id"] != done["configuration_id"]
            or not config["inputs"] or not set(config["inputs"]).issubset({"HepG2", "K562"})):
        raise ValueError("completed result and development configuration differ")
    for item in config["code"]:
        if code_ref is None:
            checked_path(item, root)
        else:
            relative = Path(item["path"]).relative_to(REPOSITORY).as_posix()
            data = subprocess.run(["git", "show", f"{code_ref}:{relative}"], cwd=REPOSITORY,
                                  capture_output=True, check=True).stdout
            if sha256(data).hexdigest() != item["sha256"]:
                raise ValueError(f"Git calculation code checksum mismatch: {relative}")
    checked_path(config["plan"], root)
    verified = set()

    def verify_inputs(node):
        if isinstance(node, dict) and {"path", "sha256"}.issubset(node):
            key = (node["path"], node["sha256"])
            if key in verified:
                return
            path = checked_path(node, root)
            verified.add(key)
            if path.suffix == ".json":
                metadata = json.loads(path.read_text())
                if metadata.get("schema") == "fp-tools-strand-functional-profiles-v1":
                    for field, hash_field in (("sites", "sites_sha256"), ("profiles_npz", "profiles_sha256")):
                        verify_inputs({"path": metadata[field], "sha256": metadata[hash_field]})
        elif isinstance(node, dict):
            for value in node.values():
                verify_inputs(value)
        elif isinstance(node, list):
            for value in node:
                verify_inputs(value)

    verify_inputs(config["inputs"])
    done = verify_outputs(result_path, "diagnostic_id", root)
    selection = verify_outputs(checked_path(done["outputs"]["selection"], root), "selection_id", root)
    if (selection["selection_id"] != done["selection_id"]
            or selection["configuration_id"] != config["configuration_id"]):
        raise ValueError("completed cohort identity differs from diagnostic")
    checked_caches = set()
    for item in done["evaluations"]:
        evaluation = verify_outputs(checked_path(item, root), "evaluation_id", root)
        cache_path = result_path.parent / evaluation["cell"] / "profile_cache.json"
        if cache_path not in checked_caches:
            verify_outputs(cache_path, "cache_id", root)
            checked_caches.add(cache_path)
    return done


def load_selected_profiles(inputs, root, selected, outdir):
    """Reconstruct the new subsets' profiles and verify their stored scores."""
    path = outdir / "profile_cache.json"
    selection_hash = canonical_id(selected[SITE_KEYS].to_dict("records"))
    if path.exists():
        done = verify_outputs(path, "cache_id", root)
        if done["selection_hash"] != selection_hash:
            raise ValueError("completed profile cache selection differs")
        with np.load(checked_path(done["outputs"]["profiles"], root), allow_pickle=False) as data:
            return {name: data[name] for name in METHODS}
    cell = str(selected.cell.iloc[0])
    sites, observed, expected, valid, _ = read_profiles(inputs["test"], root, cell)
    order = pd.Index(site_hashes(sites)).get_indexer(site_hashes(selected))
    if np.any(order < 0) or not valid[order].all():
        raise ValueError("selected cohort is outside valid profile support")
    values = {}
    for name, entry in inputs["tracks"].items():
        values[name], finite = extract_prediction_profiles(selected, checked_path(entry, root), 100, require_dense=name != "dwm")
        if not finite.all():
            raise ValueError("selected cohort is outside recorded common prediction support")
    methods = score_methods(observed[order], values["dwm"], expected[order], values["bias"],
                            values["regulatory"], values["nobias"], np.arange(-100, 101), dispersion=inputs["dispersion"])
    stored = pd.read_csv(checked_path(inputs["scores"], root), sep="\t")
    original_order = pd.Index(site_hashes(stored)).get_indexer(site_hashes(selected))
    if np.any(original_order < 0):
        raise ValueError("selected loci missing from stored scores")
    for name, (scores, _) in methods.items():
        np.testing.assert_allclose(scores, stored.iloc[original_order][name], rtol=1e-6, atol=1e-7,
                                   err_msg=f"reconstructed operator differs from recorded {name}")
    output = outdir / "selected_profiles.npz"
    np.savez_compressed(output, **{name: value[1] for name, value in methods.items()})
    done = {"selection_hash": selection_hash, "scores_verified_against_original": True,
            "outputs": {"profiles": record(output)}}
    done["cache_id"] = canonical_id(done)
    write_immutable(path, done)
    return {name: value[1] for name, value in methods.items()}


def render_profiles(profiles, metrics, diagnostics, policy, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(path, metadata={"CreationDate": None, "ModDate": None}) as pdf:
        for (cell, tf), task in profiles.groupby(["cell", "tf"], sort=True):
            fig, axes = plt.subplots(2, 3, figsize=(11.69, 8.27))
            fig.suptitle(f"{cell} {tf} · {policy}\nPost-hoc covariate-matched development sensitivity", fontsize=13)
            for ax, (method, (label, unit, color)) in zip(axes.flat, METHODS.items()):
                curve = task[task.method.eq(method)].sort_values("position")
                row = metrics[metrics.tf.eq(tf) & metrics.method.eq(method)].iloc[0]
                ax.plot(curve.position, curve.positive_minus_negative, color=color, linewidth=1.3)
                ax.fill_between(curve.position, curve.lower_95, curve.upper_95, color=color, alpha=.18)
                ax.axhline(0, color="#888888", linewidth=.6)
                ax.axvspan(-5, 5, color="#cccccc", alpha=.2)
                ax.set(title=f"{label}\nAUROC {row.auroc:.3f} · AUPRC {row.auprc:.3f}",
                       xlabel="Position relative to motif center (bp)", ylabel=f"Positive − negative ({unit})", xlim=(-100, 100))
                ax.title.set_fontsize(9)
                ax.tick_params(labelsize=8)
                ax.spines[["top", "right"]].set_visible(False)
            d = diagnostics[diagnostics.tf.eq(tf)].iloc[0]
            fig.text(.06, .045, f"Matched pairs: {d.matched_pairs}; {d.status}; maximum matched-covariate SMD {d.maximum_matched_covariate_smd:.3f}.\n"
                     "Same stored scores and fixed geometry. Pointwise chromosome-bootstrap bands; prediction panels have different units.\n"
                     "Conditional sensitivity, not new holdout evidence, causal identification, or a promotion result.", fontsize=8)
            fig.subplots_adjust(left=.08, right=.98, top=.83, bottom=.17, wspace=.38, hspace=.48)
            pdf.savefig(fig); plt.close(fig)


def evaluate_policy(cell, policy, selected, scores, profile_sites, profiles, diagnostics, plan, outdir, root):
    manifest_path = outdir / "evaluation_manifest.json"
    if manifest_path.exists():
        return verify_outputs(manifest_path, "evaluation_id", root)
    metric_rows, bootstrap_rows, profile_rows = [], [], []
    for tf, task in selected.groupby("tf", sort=True):
        if len(task) < 4:
            continue
        indexes = pd.Index(site_hashes(profile_sites)).get_indexer(site_hashes(task))
        score_indexes = pd.Index(site_hashes(scores)).get_indexer(site_hashes(task))
        if np.any(indexes < 0) or np.any(score_indexes < 0):
            raise ValueError("matched cohort does not align to scored/profile support")
        stored = scores.iloc[score_indexes].reset_index(drop=True)
        if not np.array_equal(stored.chip_label, task.chip_label):
            raise ValueError("matched cohort labels differ from stored comparison")
        labels = task.chip_label.to_numpy(dtype=int)
        d = diagnostics[diagnostics.tf.eq(tf)].iloc[0]
        for method in METHODS:
            values = profiles[method][indexes]
            row = metric_row(task, stored[method].to_numpy(), values, np.arange(-100, 101),
                             method=method, minimum_sites_per_class=plan["minimum_sites_per_class"])
            pairs = pd.DataFrame({"pair": task.pair_id.to_numpy(), "label": labels, "score": stored[method].to_numpy()}).pivot(
                index="pair", columns="label", values="score")
            row.update(policy=policy, status=d.status,
                       matched_pair_concordance=float(((pairs[1] > pairs[0]) + .5 * (pairs[1] == pairs[0])).mean()))
            metric_rows.append(row)
            curve = aggregate_curve(values, labels, task.TFBS_chr.to_numpy(),
                                    iterations=plan["aggregate_bootstrap_iterations"], seed=plan["seed"])
            profile_rows.extend({"cell": cell, "tf": tf, "policy": policy, "method": method, "position": p,
                                 "positive_mean": curve["positive_mean"][i], "negative_mean": curve["negative_mean"][i],
                                 "positive_minus_negative": curve["difference"][i],
                                 "lower_95": curve["lower_95"][i], "upper_95": curve["upper_95"][i]}
                                for i, p in enumerate(range(-100, 101)))
        if d.power_pass:
            for baseline in (DWM, "raw_geometry", "ChromBPNet_bias_conventional_geometry"):
                bootstrap_rows.append({"cell": cell, "tf": tf, "policy": policy, "method": PARAMETRIC,
                                       "baseline": baseline, "status": d.status,
                                       **block_bootstrap_delta(task, stored[PARAMETRIC].to_numpy(), stored[baseline].to_numpy(),
                                                               iterations=plan["bootstrap_iterations"], seed=plan["seed"])})
    metrics, bootstrap, curves = pd.DataFrame(metric_rows), pd.DataFrame(bootstrap_rows), pd.DataFrame(profile_rows)
    if metrics.empty:
        metrics = pd.DataFrame(columns=["cell", "tf", "policy", "method", "status", "auroc", "auprc"])
    if bootstrap.empty:
        bootstrap = pd.DataFrame(columns=["cell", "tf", "policy", "method", "baseline", "status",
                                          "auroc_gain_lower_95", "auroc_gain_upper_95",
                                          "relative_auprc_gain_lower_95", "relative_auprc_gain_upper_95"])
    if curves.empty:
        curves = pd.DataFrame(columns=["cell", "tf", "policy", "method", "position",
                                      "positive_mean", "negative_mean", "positive_minus_negative", "lower_95", "upper_95"])
    outputs = {}
    for name, table in {"metrics": metrics, "bootstrap": bootstrap, "profiles": curves,
                        "shape": shape_audit(curves) if not curves.empty else pd.DataFrame(
                            columns=["cell", "tf", "method", "canonical_protection"])}.items():
        path = outdir / f"{name}.tsv"
        table.to_csv(path, sep="\t", index=False)
        outputs[name] = record(path)
    if policy == plan["primary_policy"] and not curves.empty:
        pdf = outdir / "matched_aggregate_profiles.pdf"
        render_profiles(curves, metrics, diagnostics, policy, pdf)
        outputs["pdf"] = record(pdf)
    done = {"cell": cell, "policy": policy, "outputs": outputs, "promotion_eligible": False}
    done["evaluation_id"] = canonical_id(done)
    write_immutable(manifest_path, done)
    print(f"Completed {cell} {policy}: {len(metrics)} method/task rows", flush=True)
    return done


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--verify-only", action="store_true", help="Verify existing artifacts without computing any stage")
    parser.add_argument("--code-ref", help="Calculation-code commit hash, only with --verify-only")
    args = parser.parse_args(argv)
    if args.verify_only:
        verify_historical_result(args.outdir / "diagnostic_manifest.json", args.data_root, code_ref=args.code_ref)
        print("Verified historical inputs, calculation code and nested outputs; no stage rerun.")
        return 0
    if args.code_ref is not None:
        parser.error("--code-ref requires --verify-only")
    plan = json.loads(args.plan.read_text())
    inputs = resolve_plan(plan, args.data_root)
    code_names = [Path(__file__), SCRIPT_DIR / "audit_chrombpnet_comparison.py", SCRIPT_DIR / "audit_direct_bias_safety.py",
                  SCRIPT_DIR / "build_strand_functional_profiles.py", SCRIPT_DIR / "evaluate_chrombpnet_reference.py",
                  SCRIPT_DIR / "evaluate_frozen_functional_policy.py", SCRIPT_DIR / "evaluate_parametric_factorization.py",
                  SCRIPT_DIR / "evaluate_functional_footprints.py", REPOSITORY / "src/fp_tools/tools/parametric_bias.py",
                  REPOSITORY / "src/fp_tools/tools/functional_footprints.py"]
    configuration = {"plan": record(args.plan), "inputs": inputs, "code": [record(path) for path in code_names],
                     "runtime": {"python": sys.version, **{name: version(name) for name in ("numpy", "pandas", "scipy", "scikit-learn", "pyBigWig")}}}
    configuration["configuration_id"] = canonical_id(configuration)
    args.outdir.mkdir(parents=True, exist_ok=True)
    write_immutable(args.outdir / "configuration.freeze.json", configuration)
    result_path = args.outdir / "diagnostic_manifest.json"
    if result_path.exists():
        done = verify_outputs(result_path, "diagnostic_id", args.data_root)
        if done["configuration_id"] != configuration["configuration_id"]:
            raise ValueError("completed diagnostic configuration mismatch")
        for item in done["evaluations"]:
            verify_outputs(checked_path(item, args.data_root), "evaluation_id", args.data_root)
        print("Verified completed covariate diagnostic; no stage rerun.")
        return 0
    selection_path = args.outdir / "matched_cohorts.freeze.json"
    if selection_path.exists():
        selection = verify_outputs(selection_path, "selection_id", args.data_root)
        if selection["configuration_id"] != configuration["configuration_id"]:
            raise ValueError("matched cohort configuration mismatch")
        cohorts = pd.read_csv(checked_path(selection["outputs"]["cohorts"], args.data_root), sep="\t")
        diagnostics = pd.read_csv(checked_path(selection["outputs"]["diagnostics"], args.data_root), sep="\t")
    else:
        scales = {}
        for cell, source in inputs.items():
            train, observed, _, valid, train_metadata = read_profiles(source["training"], args.data_root, cell, training=True)
            accessibility = extract_accessibility(train, source["accessibility"], args.data_root)
            train = add_covariates(train, observed, accessibility=accessibility).loc[valid]
            scales[cell] = {tf: fit_scales(task, FEATURES, TRAIN_CHROMOSOMES) for tf, task in train.groupby("tf", sort=True)}
            test_metadata = json.loads(checked_path(source["test"], args.data_root).read_text())["metadata"]
            for key in ("bam_sha256", "bias_model_sha256", "read_shift", "genome_sha256", "minimum_mapq", "keep_duplicates"):
                if train_metadata[key] != test_metadata[key]:
                    raise ValueError(f"training covariate source differs from recorded comparison: {key}")
        scale_freeze = {"configuration_id": configuration["configuration_id"], "scales": scales,
                        "training_chromosomes": TRAIN_CHROMOSOMES, "model_scores_consulted": False}
        write_immutable(args.outdir / "covariate_scales.freeze.json", scale_freeze)
        cohort_parts, diagnostics_rows = [], []
        for cell, source in inputs.items():
            sites, observed, _, valid, _ = read_profiles(source["test"], args.data_root, cell)
            accessibility = extract_accessibility(sites, source["accessibility"], args.data_root)
            sites = add_covariates(sites, observed, accessibility=accessibility)
            # Read only coordinates and known labels; no model score column is
            # consulted until every policy's selected cohort has been frozen.
            support = pd.read_csv(checked_path(source["scores"], args.data_root), sep="\t", usecols=SITE_KEYS + ["chip_label"])
            cohort = support.merge(sites, on=SITE_KEYS, how="left", validate="one_to_one", suffixes=("_score", ""), indicator=True)
            if not cohort._merge.eq("both").all() or not cohort.chip_label.eq(cohort.chip_label_score).all():
                raise ValueError("common support fails coordinate/label join")
            if not valid[pd.Index(site_hashes(sites)).get_indexer(site_hashes(cohort))].all():
                raise ValueError("original score support includes an invalid profile")
            cohort = cohort.drop(columns=["_merge", "chip_label_score"])
            for policy, features in plan["policies"].items():
                for tf, task in cohort.groupby("tf", sort=True):
                    selected = match_sites(task, features, scales[cell][tf], caliper=plan["caliper_standard_deviations_per_feature"])
                    selected["policy"] = policy
                    cohort_parts.append(selected)
                    diagnostics_rows.append({"cell": cell, "tf": tf, "policy": policy,
                                             **matched_diagnostics(task, selected, features, plan["minimum_sites_per_class"], plan["maximum_absolute_smd"])})
        cohorts, diagnostics = pd.concat(cohort_parts, ignore_index=True), pd.DataFrame(diagnostics_rows)
        outputs = {"scales": record(args.outdir / "covariate_scales.freeze.json")}
        for name, table in (("cohorts", cohorts), ("diagnostics", diagnostics)):
            path = args.outdir / f"matched_{name}.tsv"
            table.to_csv(path, sep="\t", index=False)
            outputs[name] = record(path)
        selection = {"configuration_id": configuration["configuration_id"], "outputs": outputs,
                     "model_scores_consulted": False, "candidate_retuned": False}
        selection["selection_id"] = canonical_id(selection)
        write_immutable(selection_path, selection)
    print("All matched cohorts frozen; evaluating stored scores and exact reconstructed profiles.", flush=True)
    print(diagnostics.loc[diagnostics.policy.eq(plan["primary_policy"]),
                          ["cell", "tf", "matched_pairs", "status"]].to_string(index=False), flush=True)
    evaluations = []
    for cell, source in inputs.items():
        cell_dir = args.outdir / cell
        cell_dir.mkdir(exist_ok=True)
        profile_sites = cohorts.loc[cohorts.cell.eq(cell)].drop_duplicates(SITE_KEYS).sort_values(SITE_KEYS).reset_index(drop=True)
        if profile_sites.empty:
            continue
        profiles = load_selected_profiles(source, args.data_root, profile_sites, cell_dir)
        scores = pd.read_csv(checked_path(source["scores"], args.data_root), sep="\t")
        for policy in plan["policies"]:
            outdir = cell_dir / policy
            outdir.mkdir(exist_ok=True)
            selected = cohorts.loc[cohorts.cell.eq(cell) & cohorts.policy.eq(policy)].copy()
            diagnostics_part = diagnostics.loc[diagnostics.cell.eq(cell) & diagnostics.policy.eq(policy)]
            evaluate_policy(cell, policy, selected, scores, profile_sites, profiles, diagnostics_part, plan, outdir, args.data_root)
            evaluations.append(record(outdir / "evaluation_manifest.json"))
    done = {"schema": SCHEMA, "configuration_id": configuration["configuration_id"], "selection_id": selection["selection_id"],
            "outputs": {"selection": record(selection_path)}, "evaluations": evaluations,
            "evidence": "post_hoc_development_sensitivity", "promotion_eligible": False,
            "new_holdout_labels_read": False, "candidate_retuned": False}
    done["diagnostic_id"] = canonical_id(done)
    write_immutable(result_path, done)
    print(result_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
