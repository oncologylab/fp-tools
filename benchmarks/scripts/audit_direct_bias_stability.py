#!/usr/bin/env python3
"""Apply the recorded direct-bias operator to missing development test profiles.

Use fixed sites, cached sequence propensities and checksum-locked existing BAMs.
Retain zero-count sites, distinguish fixed from sample-matched DWM expectations,
and never reinterpret tuning matrices as independent test evidence. No fitting
or new holdout/control labels are permitted by this diagnostic.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
from importlib.metadata import version
import json
import os
from pathlib import Path
import re
import subprocess
import sys

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))
sys.path.insert(0, str(SCRIPT_DIR))

from audit_chrombpnet_comparison import (
    DWM, PARAMETRIC, SITE_KEYS, canonical_id, check_identity, checked_path,
    data_path, shape_audit, standardized_difference,
)
from audit_direct_bias_safety import direct_scores, record, write_immutable
from build_strand_functional_profiles import extract_strand_cut_profiles, orient_strand_log_bias, site_hashes
from evaluate_chrombpnet_covariate_diagnostic import TEST_CHROMOSOMES, resolve_plan, verify_outputs
from evaluate_chrombpnet_reference import extract_prediction_profiles, metric_row, score_methods as reference_score_methods
from evaluate_frozen_functional_depth_matrix import expected_from_cached_log_bias
from evaluate_frozen_functional_policy import aggregate_curve
from evaluate_parametric_factorization import block_bootstrap_delta, geometry_score
from fp_tools.tools.functional_footprints import normalize_functional_profiles

SCHEMA = "fp-tools-direct-bias-stability-v1"
RAW = "raw_geometry"
SAMPLE_DWM = "DWM_sample_matched_geometry"
METHODS = (RAW, DWM, SAMPLE_DWM, PARAMETRIC)
RUN_KEYS = ["cell", "sample", "depth", "seed"]
POSITIONS = np.arange(-100, 101)


def validate_plan(plan):
    if (plan.get("schema") != "fp-tools-direct-bias-stability-plan-v1"
            or plan.get("cells") != ["HepG2", "K562"]
            or plan.get("evaluation_chromosomes") != TEST_CHROMOSOMES
            or any(plan.get(key) is not False for key in
                   ("candidate_retuning_allowed", "new_holdout_labels_allowed", "final_untouched_control_allowed"))):
        raise ValueError("stability accepts only the frozen, previously examined development scope")
    expected = {
        "depths": ["10m", "25m", "50m", "full"], "seeds": list(range(2026, 2031)),
        "replicates": [1, 2, 3], "cohorts": ["original_support_with_expected_coverage", "motif_accessibility_observed_total"],
        "cohort_membership": "fixed_original_or_matched_sites_intersected_with_all_run_expected_coverage",
        "missing_expected_sites": "exclude_before_scoring_and_retain_complete_matched_pairs",
        "zero_cut_sites": "retain_as_ties_on_fixed_support_and_report_counts",
        "bootstrap_runs": "all_full_replicates_and_10m_seed_2026_reference",
        "dwm_comparators": ["original_full_depth_replicate_1", "sample_and_depth_matched"],
        "available_matrix": "replicate_1_depth_seed_matrix_plus_all_full_replicates",
        "full_depth_randomization": "one_deterministic_endpoint_per_replicate",
    }
    if (any(plan.get(key) != value for key, value in expected.items())
            or plan.get("minimum_sites_per_class", 0) < 200
            or plan.get("bootstrap_iterations", 0) < 1000
            or plan.get("aggregate_bootstrap_iterations", 0) < 500):
        raise ValueError("stability plan changed the fixed cohort, matrix or evidence requirements")


def validate_matrix(frame):
    expected = {(cell, f"{cell}_rep1", depth, seed)
                for cell in ("HepG2", "K562") for depth in ("10m", "25m", "50m")
                for seed in range(2026, 2031)}
    expected |= {(cell, f"{cell}_rep{rep}", "full", 2026)
                 for cell in ("HepG2", "K562") for rep in (1, 2, 3)}
    if (frame.duplicated(RUN_KEYS).any()
            or set(frame[RUN_KEYS].itertuples(index=False, name=None)) != expected):
        raise ValueError("input matrix must contain exactly the 36 prespecified development runs")


def align_sites(source, selected):
    if source.duplicated(SITE_KEYS).any() or selected.duplicated(SITE_KEYS).any():
        raise ValueError("site coordinates must be unique within a cell and TF")
    index = pd.MultiIndex.from_frame(source[SITE_KEYS]).get_indexer(pd.MultiIndex.from_frame(selected[SITE_KEYS]))
    if np.any(index < 0) or not np.array_equal(source.iloc[index].chip_label, selected.chip_label):
        raise ValueError("fixed cohort loci or labels differ from archived support")
    return index


def expected_for_observed(plus, minus, plus_log_bias, minus_log_bias):
    """Inputs are already motif-oriented; preserve separate strand totals."""
    first, valid_first = expected_from_cached_log_bias(plus, plus_log_bias)
    second, valid_second = expected_from_cached_log_bias(minus, minus_log_bias)
    if not (valid_first & valid_second).all():
        raise ValueError("frozen candidate lacks a strand propensity on fixed support")
    return first + second


def score_panel(observed, expected, fixed_dwm, sample_dwm, dispersion):
    result = {RAW: (geometry_score(observed, POSITIONS), observed)}
    for method, prediction in ((DWM, fixed_dwm), (SAMPLE_DWM, sample_dwm), (PARAMETRIC, expected)):
        if not np.all(np.asarray(prediction).sum(axis=1)[observed.sum(axis=1) > 0] > 0):
            raise ValueError(f"{method} has no expectation for an informative fixed site")
        score, _, residual = direct_scores(observed, prediction, dispersion)
        result[method] = (score, normalize_functional_profiles(residual, POSITIONS))
    return result


def verify_stage(path, configuration_id, root, *, expected_run=None):
    done = json.loads(path.read_text())
    check_identity(done, "stage_id")
    if done["configuration_id"] != configuration_id:
        raise ValueError("completed stage belongs to a different frozen configuration")
    if expected_run is not None and done.get("run") != expected_run:
        raise ValueError("completed stage belongs to a different sample, depth or seed")
    for entry in done["outputs"].values():
        checked_path(entry, root)
    return done


def freeze_inputs(plan, root):
    """Verify source identities and all transitive data inputs before new scoring."""
    records = {}

    def check(item):
        key = (item["path"], item["sha256"])
        if key not in records:
            records[key] = record(checked_path(item, root))
        return Path(records[key]["path"])

    source_plan = json.loads(check(plan["source_plan"]).read_text())
    cells = resolve_plan(source_plan, root)
    matrices = []
    for item in plan["input_registries"]:
        document = json.loads(check(item).read_text())
        check_identity(document, "depth_input_id")
        if document.get("models_refitted_by_depth") is not False:
            raise ValueError("archived matrix does not certify fixed models")
        matrices.append(pd.read_csv(check(document["signals"]), sep="\t"))
    matrix = pd.concat(matrices, ignore_index=True)
    validate_matrix(matrix)
    for row in matrix.itertuples(index=False):
        for field in ("bam", "bai", "dwm_expected"):
            check({"path": getattr(row, field), "sha256": getattr(row, field + "_sha256")})
        print(f"Verified inputs: {row.cell} {row.sample} {row.depth} seed {row.seed}", flush=True)
    selection_path = check(plan["matched_cohorts"])
    selection = verify_outputs(selection_path, "selection_id", root)
    matched = pd.read_csv(check(selection["outputs"]["cohorts"]), sep="\t")
    for entry in selection["outputs"].values():
        check(entry)
    references = {}
    for cell, inputs in cells.items():
        if inputs["dispersion"] != plan["dispersion"][cell]:
            raise ValueError("dispersion differs from the exact recorded comparison")
        for field in ("comparison", "freeze", "scores", "test"):
            check(inputs[field])
        check(inputs["tracks"]["dwm"])
        document = json.loads(check(inputs["test"]).read_text())
        metadata = document["metadata"]
        if (document.get("schema") != "fp-tools-strand-functional-profiles-v1"
                or metadata.get("labels_used") is not False or metadata.get("flank") != 100
                or metadata.get("read_shift") != [4, -5] or metadata.get("minimum_mapq") != 30
                or metadata.get("keep_duplicates") is not False
                or metadata.get("bias_model_sha256") != plan["bias_model_sha256"]):
            raise ValueError("test profiles changed the frozen model or read extraction contract")
        for field, hash_field in (("sites", "sites_sha256"), ("profiles_npz", "profiles_sha256")):
            check({"path": document[field], "sha256": document[hash_field]})
        check({"path": metadata["bias_model"], "sha256": metadata["bias_model_sha256"]})
        # Freeze the sidecar as well as the coefficient NPZ, without refitting.
        model_sidecar = data_path(metadata["bias_model"], root).with_suffix(".json")
        if model_sidecar.is_file():
            check(record(model_sidecar))
        references[cell] = {"inputs": inputs, "profile_manifest": document}
    return matrix, matched, references, list(records.values())


def load_reference(cell, reference, matched, root):
    document = reference["profile_manifest"]
    source = pd.read_csv(data_path(document["sites"], root), sep="\t")
    sites = pd.read_csv(data_path(reference["inputs"]["scores"]["path"], root), sep="\t")
    if (set(source.cell) != {cell} or not set(source.TFBS_chr).issubset(TEST_CHROMOSOMES)
            or not source.chip_label.isin([0, 1]).all()):
        raise ValueError("source profiles are outside the declared development scope")
    order = align_sites(source, sites)
    for name in ("motif_score", "accessibility"):
        sites[name] = source.iloc[order][name].to_numpy()
    with np.load(data_path(document["profiles_npz"], root), allow_pickle=False) as data:
        if not np.array_equal(data["site_hash"], site_hashes(source)) or not data["valid"][order].all():
            raise ValueError("source profile identities or validity differ from archived support")
        arrays = {key: data[key][order] for key in
                  ("plus_observed", "minus_observed", "plus_expected", "minus_expected", "plus_log_bias", "minus_log_bias")}
    fixed_dwm, valid = extract_prediction_profiles(
        sites, data_path(reference["inputs"]["tracks"]["dwm"]["path"], root), 100, require_dense=False)
    if not valid.all():
        raise ValueError("fixed DWM track lacks a recorded development window")
    strict = matched.loc[matched.cell.eq(cell) & matched.policy.eq("motif_accessibility_observed_total")].copy()
    cohorts = {"original_common_support": np.arange(len(sites)),
               "motif_accessibility_observed_total": align_sites(sites, strict)}
    expected = expected_for_observed(arrays["plus_observed"], arrays["minus_observed"],
                                     arrays["plus_log_bias"], arrays["minus_log_bias"])
    np.testing.assert_allclose(expected, arrays["plus_expected"] + arrays["minus_expected"], rtol=1e-10, atol=1e-10,
                               err_msg="cached strand propensity does not reproduce original expected cuts")
    # Reproduce the historical scores, including their zero-filled expectation
    # intervals. The new support stage excludes those before any new scoring.
    methods = reference_score_methods(arrays["plus_observed"] + arrays["minus_observed"], fixed_dwm, expected,
                                      fixed_dwm, fixed_dwm, fixed_dwm, POSITIONS,
                                      dispersion=reference["inputs"]["dispersion"])
    for method in (RAW, DWM, PARAMETRIC):
        np.testing.assert_allclose(methods[method][0], sites[method], rtol=1e-6, atol=1e-7,
                                   err_msg=f"operator differs from original {cell} {method}")
    return sites, arrays, fixed_dwm, cohorts


def qualify_support(sites, strict, present):
    """Filter using expectation availability only; keep matched pairs intact."""
    present = np.asarray(present, dtype=bool)
    if present.shape != (len(sites),) or strict.groupby("pair_id").size().ne(2).any():
        raise ValueError("coverage mask or frozen matched pairs are malformed")
    strict = strict.copy()
    strict["_present"] = present[align_sites(sites, strict)]
    complete_pairs = strict.groupby("pair_id")._present.transform("all")
    retained = sites.loc[present].reset_index(drop=True)
    return np.flatnonzero(present), {
        "original_support_with_expected_coverage": np.arange(len(retained)),
        "motif_accessibility_observed_total": align_sites(retained, strict.loc[complete_pairs]),
    }


def freeze_support(cell, frame, sites, arrays, fixed_dwm, matched, outdir, configuration_id, root):
    path = outdir / "support.json"
    if path.exists():
        done = verify_stage(path, configuration_id, root)
        if done.get("cell") != cell:
            raise ValueError("frozen support belongs to a different development cell")
        with np.load(checked_path(done["outputs"]["expected_tracks"], root), allow_pickle=False) as data:
            if not np.array_equal(data["site_hash"], site_hashes(sites)):
                raise ValueError("frozen expected-track support has changed coordinates")
            tracks = {key: data[key] for key in data.files if key != "site_hash"}
        coverage = pd.read_csv(checked_path(done["outputs"]["coverage"], root), sep="\t")
        if not np.array_equal(site_hashes(coverage), site_hashes(sites)):
            raise ValueError("frozen coverage table order differs")
    else:
        tracks = {"fixed_dwm": fixed_dwm}
        coverage = sites[SITE_KEYS + ["chip_label"]].copy()
        coverage["fixed_dwm_present"] = np.isfinite(fixed_dwm).all(axis=1) & (fixed_dwm.sum(axis=1) > 0)
        for index, row in enumerate(frame.itertuples(index=False)):
            values, valid = extract_prediction_profiles(sites, data_path(row.dwm_expected, root), 100, require_dense=False)
            tracks[f"run_{index}"] = values
            coverage[f"run_{index}_present"] = valid & np.isfinite(values).all(axis=1) & (values.sum(axis=1) > 0)
        coverage["all_run_dwm_present"] = coverage.filter(regex="_present$").all(axis=1)
        outputs = {}
        target = outdir / "expected_tracks.npz"
        np.savez_compressed(target, **tracks, site_hash=site_hashes(sites))
        outputs["expected_tracks"] = record(target)
        target = outdir / "expected_coverage.tsv.gz"
        coverage.to_csv(target, sep="\t", index=False, compression={"method": "gzip", "mtime": 0})
        outputs["coverage"] = record(target)
        done = {"configuration_id": configuration_id, "cell": cell, "outputs": outputs,
                "runs": frame[RUN_KEYS].to_dict("records"), "selection_uses_new_scores": False}
        done["stage_id"] = canonical_id(done)
        write_immutable(path, done)
    if done["runs"] != frame[RUN_KEYS].to_dict("records"):
        raise ValueError("frozen expectation tracks differ from requested run order")
    strict = matched.loc[matched.cell.eq(cell) & matched.policy.eq("motif_accessibility_observed_total")]
    selected, cohorts = qualify_support(sites, strict, coverage.all_run_dwm_present.to_numpy())
    print(f"Froze {cell} expected coverage: {len(selected)}/{len(sites)} original sites; "
          f"{len(cohorts['motif_accessibility_observed_total'])} matched sites", flush=True)
    return (sites.iloc[selected].reset_index(drop=True), {key: value[selected] for key, value in arrays.items()},
            {key: value[selected] for key, value in tracks.items()}, cohorts)


def extract_run(row, sites, arrays, reference, sample_dwm, outdir, configuration_id, root):
    path = outdir / "extraction.json"
    if path.exists():
        done = verify_stage(path, configuration_id, root,
                            expected_run={key: getattr(row, key) for key in RUN_KEYS})
        with np.load(checked_path(done["outputs"]["profiles"], root), allow_pickle=False) as data:
            if not np.array_equal(data["site_hash"], site_hashes(sites)):
                raise ValueError("resumed extraction site order differs")
            return data["observed"], data["expected"], data["sample_dwm"]
    reference_bam = reference["profile_manifest"]["metadata"]["bam_sha256"]
    reused = row.bam_sha256 == reference_bam
    if reused:
        plus, minus = arrays["plus_observed"], arrays["minus_observed"]
    else:
        plus, minus, valid = extract_strand_cut_profiles(
            sites.reset_index(drop=True), data_path(row.bam, root), flank=100, read_shift=(4, -5),
            minimum_mapq=30, keep_duplicates=False)
        if not valid.all():
            raise ValueError("BAM lacks a fixed development window")
        # The strand swap/reversal is an involution and applies to cut arrays too.
        plus, minus, _ = orient_strand_log_bias(plus, minus, sites.TFBS_strand.tolist())
    expected = expected_for_observed(plus, minus, arrays["plus_log_bias"], arrays["minus_log_bias"])
    observed = plus + minus
    target = outdir / "profiles.npz"
    np.savez_compressed(target, observed=observed, expected=expected, sample_dwm=sample_dwm, site_hash=site_hashes(sites))
    done = {"configuration_id": configuration_id, "run": {key: getattr(row, key) for key in RUN_KEYS},
            "site_identity": canonical_id(sites[SITE_KEYS].to_dict("records")),
            "reference_observed_reused": reused, "outputs": {"profiles": record(target)}}
    done["stage_id"] = canonical_id(done)
    write_immutable(path, done)
    print(f"Extracted {row.sample} {row.depth} {row.seed}; original counts reused={reused}", flush=True)
    return observed, expected, sample_dwm


def evaluate_run(row, sites, cohorts, methods, observed, plan, outdir, configuration_id, root):
    path = outdir / "evaluation.json"
    if path.exists():
        return verify_stage(path, configuration_id, root,
                            expected_run={key: getattr(row, key) for key in RUN_KEYS})
    metric_rows, bootstrap_rows, curves, coverage = [], [], [], []
    run = {key: getattr(row, key) for key in RUN_KEYS}
    bootstrap_run = row.depth == "full" or (row.depth == "10m" and row.seed == 2026)
    for cohort, membership in cohorts.items():
        for tf, task in sites.iloc[membership].groupby("tf", sort=True):
            indexes = task.index.to_numpy(dtype=int)
            labels = task.chip_label.to_numpy(dtype=int)
            totals = observed[indexes].sum(axis=1)
            n_positive, n_negative = int((labels == 1).sum()), int((labels == 0).sum())
            coverage.append({**run, "cohort": cohort, "tf": tf, "n_positive": n_positive, "n_negative": n_negative,
                             "zero_cut_sites": int((totals == 0).sum()),
                             "observed_total_smd": standardized_difference(np.log1p(totals[labels == 1]), np.log1p(totals[labels == 0])),
                             "motif_score_smd": standardized_difference(task.loc[labels == 1, "motif_score"], task.loc[labels == 0, "motif_score"]),
                             "reference_log_accessibility_smd": standardized_difference(
                                 np.log1p(task.loc[labels == 1, "accessibility"]), np.log1p(task.loc[labels == 0, "accessibility"])),
                             "status": "eligible" if min(n_positive, n_negative) >= plan["minimum_sites_per_class"] else "underpowered"})
            if min(n_positive, n_negative) < 2:
                continue
            for method, (score, profiles) in methods.items():
                metric_rows.append({**metric_row(task, score[indexes], profiles[indexes], POSITIONS,
                                                 method=method, minimum_sites_per_class=plan["minimum_sites_per_class"]),
                                    **run, "cohort": cohort})
                curve = aggregate_curve(profiles[indexes], labels, task.TFBS_chr.to_numpy(),
                                        iterations=plan["aggregate_bootstrap_iterations"], seed=plan["seed"])
                curves.extend({**run, "cohort": cohort, "tf": tf, "method": method, "position": position,
                               "positive_mean": curve["positive_mean"][i], "negative_mean": curve["negative_mean"][i],
                               "positive_minus_negative": curve["difference"][i],
                               "lower_95": curve["lower_95"][i], "upper_95": curve["upper_95"][i]}
                              for i, position in enumerate(range(-100, 101)))
            if bootstrap_run and min(n_positive, n_negative) >= plan["minimum_sites_per_class"]:
                for baseline in (DWM, SAMPLE_DWM, RAW):
                    bootstrap_rows.append({**run, "cohort": cohort, "tf": tf, "method": PARAMETRIC, "baseline": baseline,
                                           **block_bootstrap_delta(task, methods[PARAMETRIC][0][indexes], methods[baseline][0][indexes],
                                                                   iterations=plan["bootstrap_iterations"], seed=plan["seed"])})
    metrics = pd.DataFrame(metric_rows)
    for method, suffix in ((DWM, "fixed_dwm"), (SAMPLE_DWM, "sample_dwm"), (RAW, "raw")):
        baseline = metrics.loc[metrics.method.eq(method), ["cohort", "tf", "auroc", "auprc"]].rename(
            columns={"auroc": f"{suffix}_auroc", "auprc": f"{suffix}_auprc"})
        metrics = metrics.merge(baseline, on=["cohort", "tf"], validate="many_to_one")
        metrics[f"auroc_gain_over_{suffix}"] = metrics.auroc - metrics[f"{suffix}_auroc"]
        metrics[f"relative_auprc_gain_over_{suffix}"] = metrics.auprc / metrics[f"{suffix}_auprc"].clip(lower=1e-8) - 1
    profiles = pd.DataFrame(curves)
    shapes = pd.concat([shape_audit(group).assign(cohort=cohort) for cohort, group in profiles.groupby("cohort")], ignore_index=True)
    for key, value in run.items():
        shapes[key] = value
    bootstrap = pd.DataFrame(bootstrap_rows) if bootstrap_rows else pd.DataFrame(columns=RUN_KEYS + ["cohort", "tf", "method", "baseline"])
    outputs = {}
    for name, table in (("metrics", metrics), ("bootstrap", bootstrap), ("profiles", profiles), ("shapes", shapes),
                        ("coverage", pd.DataFrame(coverage))):
        target = outdir / f"{name}.tsv.gz"
        table.to_csv(target, sep="\t", index=False, compression={"method": "gzip", "mtime": 0})
        outputs[name] = record(target)
    done = {"configuration_id": configuration_id, "run": run, "outputs": outputs, "bootstrap_run": bootstrap_run}
    done["stage_id"] = canonical_id(done)
    write_immutable(path, done)
    print(f"Evaluated {row.sample} {row.depth} {row.seed}: {len(metrics)} method/task/cohort rows", flush=True)
    return done


def render_stability(metrics, profiles, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    colors = {DWM: "#777777", SAMPLE_DWM: "#7B3294", PARAMETRIC: "#D55E00"}
    names = {DWM: "Fixed DWM", SAMPLE_DWM: "Sample DWM", PARAMETRIC: "Parametric"}
    with PdfPages(output, metadata={"CreationDate": None, "ModDate": None}) as pdf:
        for (cell, tf), task in profiles.groupby(["cell", "tf"], sort=True):
            fig, axes = plt.subplots(2, 4, figsize=(11.69, 8.27))
            fig.suptitle(f"{cell} {tf}: direct-bias development stability\nFixed original and covariate-matched sites; no refitting", fontsize=13)
            for row_index, cohort in enumerate(("original_support_with_expected_coverage", "motif_accessibility_observed_total")):
                for column, depth in enumerate(("10m", "25m", "50m", "full")):
                    ax = axes[row_index, column]
                    panel = task.loc[task.cohort.eq(cohort) & task.depth.eq(depth)]
                    if panel.empty:
                        ax.set_axis_off(); continue
                    for method, color in colors.items():
                        data = panel.loc[panel.method.eq(method)]
                        for _, group in data.groupby(["sample", "seed"]):
                            ax.plot(group.position, group.positive_minus_negative, color=color, alpha=.12, linewidth=.6)
                        mean = data.groupby("position").positive_minus_negative.mean()
                        ax.plot(mean.index, mean.values, color=color, linewidth=1.1, label=names[method])
                    selected = metrics.loc[metrics.cell.eq(cell) & metrics.tf.eq(tf) & metrics.cohort.eq(cohort)
                                           & metrics.depth.eq(depth) & metrics.method.eq(PARAMETRIC)]
                    ax.set(title=f"{depth}: {'coverage subset' if row_index == 0 else 'matched'}; {selected.status.iloc[0]}\n"
                                 f"ΔAUROC vs raw {selected.auroc_gain_over_raw.min():+.3f} to {selected.auroc_gain_over_raw.max():+.3f}",
                           xlabel="Motif offset (bp)", xlim=(-100, 100))
                    ax.title.set_fontsize(8); ax.tick_params(labelsize=7)
                    ax.axhline(0, color="#888888", linewidth=.5); ax.axvspan(-5, 5, alpha=.15, color="#aaaaaa")
                    ax.spines[["top", "right"]].set_visible(False)
                    if column == 0:
                        ax.set_ylabel("Positive − negative\n(normalized deviance)", fontsize=8)
            axes[0, 0].legend(fontsize=7, frameon=False)
            fig.text(.065, .045, "All sites have positive fixed and sample DWM expectations in every run; complete matched pairs are retained. Thin curves: runs; thick curves: their mean.\n"
                     "Five nested seeds of replicate 1 at 10M/25M/50M; three biological replicates at full. Raw scores are a guardrail; raw count curves have different units.\n"
                     "Examined chr19–22/X. Per-run metrics and pointwise bootstrap bands are in tables. Seed means are not independent replication or promotion evidence.", fontsize=7.5)
            fig.subplots_adjust(left=.07, right=.98, top=.82, bottom=.16, wspace=.35, hspace=.5)
            pdf.savefig(fig); plt.close(fig)


def verify_result(outdir, root, code_ref=None):
    config = json.loads((outdir / "configuration.freeze.json").read_text())
    check_identity(config, "configuration_id")
    validate_plan(config["settings"])
    if code_ref is not None and not re.fullmatch(r"[0-9a-fA-F]{7,40}", code_ref):
        raise ValueError("--code-ref must be an immutable Git commit hash")
    for item in config["code"]:
        if code_ref is None:
            checked_path(item, root)
        else:
            relative = Path(item["path"]).relative_to(REPOSITORY).as_posix()
            blob = subprocess.run(["git", "show", f"{code_ref}:{relative}"], cwd=REPOSITORY, capture_output=True, check=True).stdout
            if sha256(blob).hexdigest() != item["sha256"]:
                raise ValueError("historical calculation code checksum mismatch")
    checked_path(config["plan"], root)
    for item in config["input_records"]:
        checked_path(item, root)
    result = verify_stage(outdir / "stability_manifest.json", config["configuration_id"], root)
    if result.get("schema") != SCHEMA or result.get("promotion_eligible") is not False:
        raise ValueError("unexpected stability evidence scope")
    if len(result["stages"]) != 36:
        raise ValueError("completed result lacks a prespecified run")
    if set(result["supports"]) != {"HepG2", "K562"}:
        raise ValueError("completed result lacks frozen per-cell expected coverage")
    for item in result["supports"].values():
        verify_stage(checked_path(item, root), config["configuration_id"], root)
    runs = []
    for item in result["stages"]:
        run = None
        for stage in item.values():
            done = verify_stage(checked_path(stage, root), config["configuration_id"], root, expected_run=run)
            run = done["run"]
        runs.append(run)
    validate_matrix(pd.DataFrame(runs))
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--verify-only", action="store_true")
    parser.add_argument("--code-ref", help="Immutable calculation commit, for verification only")
    args = parser.parse_args(argv)
    if args.code_ref is not None and not args.verify_only:
        parser.error("--code-ref requires --verify-only")
    if args.verify_only or (args.outdir / "stability_manifest.json").is_file():
        frozen = json.loads((args.outdir / "configuration.freeze.json").read_text())
        if record(args.plan)["sha256"] != frozen["plan"]["sha256"]:
            raise ValueError("requested plan differs from completed stability configuration")
        verify_result(args.outdir, args.data_root, args.code_ref)
        print("Verified all frozen stability inputs, code and outputs; no stage rerun.")
        return 0
    plan = json.loads(args.plan.read_text())
    validate_plan(plan)
    matrix, matched, references, inputs = freeze_inputs(plan, args.data_root)
    code = [Path(__file__), *[SCRIPT_DIR / (name + ".py") for name in
            ("audit_chrombpnet_comparison", "audit_direct_bias_safety", "build_strand_functional_profiles",
             "evaluate_chrombpnet_covariate_diagnostic", "evaluate_chrombpnet_reference", "evaluate_frozen_functional_depth_matrix",
             "evaluate_frozen_functional_policy", "evaluate_parametric_factorization", "evaluate_functional_footprints")],
            *[REPOSITORY / "src/fp_tools/tools" / (name + ".py") for name in ("parametric_bias", "functional_footprints")]]
    config = {"settings": plan, "plan": record(args.plan), "input_records": inputs, "references": references,
              "code": [record(path) for path in code],
              "runtime": {"python": sys.version,
                          "thread_environment": {name: os.environ.get(name) for name in
                                                 ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")},
                          **{name: version(name) for name in
                          ("numpy", "pandas", "scipy", "scikit-learn", "pysam", "pyBigWig", "matplotlib")}}}
    config["configuration_id"] = canonical_id(config)
    args.outdir.mkdir(parents=True, exist_ok=True)
    write_immutable(args.outdir / "configuration.freeze.json", config)
    print(f"Froze configuration {config['configuration_id']}", flush=True)
    completed, supports = [], {}
    for cell, frame in matrix.groupby("cell", sort=True):
        frame = frame.sort_values(["sample", "depth", "seed"]).reset_index(drop=True)
        sites, arrays, fixed_dwm, cohorts = load_reference(cell, references[cell], matched, args.data_root)
        support_dir = args.outdir / cell
        support_dir.mkdir(parents=True, exist_ok=True)
        sites, arrays, tracks, cohorts = freeze_support(cell, frame, sites, arrays, fixed_dwm, matched,
                                                       support_dir, config["configuration_id"], args.data_root)
        supports[cell] = record(support_dir / "support.json")
        for index, row in enumerate(frame.itertuples(index=False)):
            outdir = args.outdir / row.sample / row.depth / f"seed_{row.seed}"
            outdir.mkdir(parents=True, exist_ok=True)
            if (outdir / "evaluation.json").is_file():
                run = {key: getattr(row, key) for key in RUN_KEYS}
                verify_stage(outdir / "extraction.json", config["configuration_id"], args.data_root, expected_run=run)
                verify_stage(outdir / "evaluation.json", config["configuration_id"], args.data_root, expected_run=run)
            else:
                observed, expected, sample_dwm = extract_run(row, sites, arrays, references[cell], tracks[f"run_{index}"],
                                                           outdir, config["configuration_id"], args.data_root)
                methods = score_panel(observed, expected, tracks["fixed_dwm"], sample_dwm, plan["dispersion"][cell])
                evaluate_run(row, sites, cohorts, methods, observed, plan, outdir, config["configuration_id"], args.data_root)
            completed.append({name: record(outdir / (name + ".json")) for name in ("extraction", "evaluation")})
    outputs, tables = {}, {}
    for name in ("metrics", "bootstrap", "profiles", "shapes", "coverage"):
        frames = []
        for item in completed:
            stage = json.loads(Path(item["evaluation"]["path"]).read_text())
            frame = pd.read_csv(checked_path(stage["outputs"][name], args.data_root), sep="\t")
            if not frame.empty:
                frames.append(frame)
        tables[name] = pd.concat(frames, ignore_index=True)
        target = args.outdir / f"stability_{name}.tsv.gz"
        tables[name].to_csv(target, sep="\t", index=False, compression={"method": "gzip", "mtime": 0})
        outputs[name] = record(target)
    pdf = args.outdir / "stability_aggregate_profiles.pdf"
    render_stability(tables["metrics"], tables["profiles"], pdf)
    outputs["pdf"] = record(pdf)
    done = {"schema": SCHEMA, "configuration_id": config["configuration_id"], "outputs": outputs,
            "stages": completed, "supports": supports, "new_holdout_labels_read": False, "candidate_retuned": False,
            "promotion_eligible": False, "evidence": plan["evidence"]}
    done["stage_id"] = canonical_id(done)
    write_immutable(args.outdir / "stability_manifest.json", done)
    print(f"Completed 36 frozen development runs; result {done['stage_id']}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
