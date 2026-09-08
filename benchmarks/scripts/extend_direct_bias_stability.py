#!/usr/bin/env python3
"""Evaluate the complete 96-run development matrix without retuning.

Reuse all 36 archived observed/expected matrices. Freeze both cells' common
expectation support before extracting the 60 missing profiles, and reuse old
metrics only if site order and both cohort memberships are unchanged.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
from importlib.metadata import version
from contextlib import contextmanager
import json
import os
from pathlib import Path
import re
import subprocess
import sys

import numpy as np
import pandas as pd

REPOSITORY = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_direct_bias_stability as base
from audit_chrombpnet_comparison import canonical_id, check_identity, checked_path, data_path, file_sha256
from audit_direct_bias_safety import record, write_immutable

SCHEMA = "fp-tools-direct-bias-stability-96-v1"
THREADS = {name: "1" for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")}
INPUT_FIELDS = ("bam", "bai", "dwm_expected")


def read_document(path, field):
    document = json.loads(Path(path).read_text())
    check_identity(document, field)
    return document


def validate_plan(plan):
    fixed = {"schema": "fp-tools-direct-bias-stability-96-plan-v1",
             "evidence": "post_hoc_development_complete_replicate_depth_stability",
             "cells": ["HepG2", "K562"], "replicates": [1, 2, 3],
             "depths": ["10m", "25m", "50m", "full"], "seeds": list(range(2026, 2031)),
             "full_depth_randomization": "one_deterministic_endpoint_per_replicate",
             "evaluation_chromosomes": base.TEST_CHROMOSOMES,
             "candidate_retuning_allowed": False, "new_holdout_labels_allowed": False,
             "final_untouched_control_allowed": False,
             "candidate_and_metric_settings": "unchanged_from_the_36_run_source_configuration",
             "support_policy": "positive_expected_coverage_across_all_96_runs_and_complete_existing_matched_pairs",
             "support_freeze_order": "both_cells_before_any_new_run_extraction_or_scoring",
             "archived_profiles": "reuse_all_36_on_coordinate_verified_new_support",
             "archived_metrics": "reuse_only_if_sites_and_both_cohort_memberships_are_identical",
             "zero_count_sites": "retain_as_ties"}
    if any(plan.get(key) != value for key, value in fixed.items()):
        raise ValueError("only the fixed 96-run development extension is allowed")


def run_key(row):
    return tuple(getattr(row, key) for key in base.RUN_KEYS)


def validate_matrix(frame):
    expected = {(cell, f"{cell}_rep{rep}", depth, seed)
                for cell in ("HepG2", "K562") for rep in (1, 2, 3)
                for depth in ("10m", "25m", "50m") for seed in range(2026, 2031)}
    expected |= {(cell, f"{cell}_rep{rep}", "full", 2026)
                 for cell in ("HepG2", "K562") for rep in (1, 2, 3)}
    if frame.duplicated(base.RUN_KEYS).any() or set(frame[base.RUN_KEYS].itertuples(index=False, name=None)) != expected:
        raise ValueError("matrix must contain exactly all 96 development runs")


def validate_merge(old, new):
    base.validate_matrix(old)
    if len(new) != 60:
        raise ValueError("prepared extension must contain exactly the missing 60 runs")
    matrix = pd.concat([old, new], ignore_index=True)
    validate_matrix(matrix)
    return matrix.sort_values(base.RUN_KEYS).reset_index(drop=True)


def bind_prepared_rows(frame, stages):
    """Bind every table row to the exact correction and sampling inputs."""
    by_run = {}
    for stage in stages:
        sample, depth, seed = stage["sample"], int(stage["depth"]), int(stage["seed"])
        if depth not in (10_000_000, 25_000_000, 50_000_000) or seed not in range(2026, 2031):
            raise ValueError("prepared correction has an undeclared depth or seed")
        key = (sample.rsplit("_", 1)[0], sample, f"{depth // 1_000_000}m", seed)
        if key in by_run:
            raise ValueError("duplicate prepared correction stage")
        by_run[key] = stage
    if len(by_run) != 60 or set(by_run) != set(frame[base.RUN_KEYS].itertuples(index=False, name=None)):
        raise ValueError("prepared table and correction-stage identities differ")
    for row in frame.itertuples(index=False):
        stage = by_run[run_key(row)]
        expected = [entry for name, entry in stage["outputs"].items() if name.endswith("_expected.bw")]
        if (len(expected) != 1 or expected[0]["sha256"] != row.dwm_expected_sha256
                or stage["input_bam"]["sha256"] != row.bam_sha256):
            raise ValueError("prepared table substituted a correction track or BAM")


def verify_archived_stage(item, configuration_id, check, expected_run=None):
    stage = read_document(check(item), "stage_id")
    if stage["configuration_id"] != configuration_id or (expected_run is not None and stage.get("run") != expected_run):
        raise ValueError("archived stage identity differs from the source study")
    for entry in stage["outputs"].values():
        check(entry)
    return stage


def freeze_inputs(plan, root):
    # The completed manifest is a prerequisite, not inferred from partial BAMs.
    prepared_path = data_path(plan["prepared_inputs"]["path"], root)
    if not prepared_path.is_file():
        raise ValueError("the 60-run input preparation is incomplete; no evaluation may start")
    records = {}

    def check(entry):
        path = data_path(entry["path"], root).resolve()
        key = str(path)
        if key in records and records[key]["sha256"] != entry["sha256"]:
            raise ValueError("one input path has conflicting frozen checksums")
        if key not in records:
            checked_path(entry, root)
            records[key] = {"path": key, "sha256": entry["sha256"]}
        return path

    old_config = read_document(check(plan["source_configuration"]), "configuration_id")
    old_result = read_document(check(plan["source_result"]), "stage_id")
    base.validate_plan(old_config["settings"])
    if (old_result.get("schema") != base.SCHEMA or old_result["configuration_id"] != old_config["configuration_id"]
            or any(old_result.get(key) is not False for key in
                   ("promotion_eligible", "candidate_retuned", "new_holdout_labels_read"))
            or len(old_result["stages"]) != 36 or set(old_result["supports"]) != {"HepG2", "K562"}):
        raise ValueError("source is not the complete recorded 36-run development experiment")
    for entry in old_config["code"] + old_config["input_records"] + [old_config["plan"]]:
        check(entry)
    for entry in old_result["outputs"].values():
        check(entry)
    old_frames = []
    for entry in old_config["settings"]["input_registries"]:
        registry = read_document(check(entry), "depth_input_id")
        old_frames.append(pd.read_csv(check(registry["signals"]), sep="\t"))
    old_frame = pd.concat(old_frames, ignore_index=True)
    base.validate_matrix(old_frame)
    archives = []
    for item in old_result["stages"]:
        extraction = verify_archived_stage(item["extraction"], old_config["configuration_id"], check)
        verify_archived_stage(item["evaluation"], old_config["configuration_id"], check, extraction["run"])
        archives.append({"run": extraction["run"], **item})
    base.validate_matrix(pd.DataFrame([entry["run"] for entry in archives]))
    for cell, entry in old_result["supports"].items():
        support = verify_archived_stage(entry, old_config["configuration_id"], check)
        if support.get("cell") != cell:
            raise ValueError("source support belongs to the wrong cell")
    return freeze_prepared(plan, root, check, records, old_config, old_result, old_frame, archives)


def preparation_stage(entry, config, sample, kind, check, depth=None, seed=None):
    stage = read_document(check(entry), "stage_id")
    if (stage.get("configuration_id") != config["configuration_id"] or stage.get("sample") != sample
            or stage.get("stage") != kind or stage.get("depth") != depth or stage.get("seed") != seed):
        raise ValueError("preparation stage belongs to a different frozen input request")
    for output in stage.get("outputs", {}).values():
        check(output)
    return stage


def bind_sampling_row(row, sampling, counts, source, check):
    depth = int(row.depth.removesuffix("m")) * 1_000_000
    key = f"{depth}:{row.seed}"
    for field in ("bam", "bai"):
        entry = sampling["outputs"][key + ":" + field]
        if entry["sha256"] != getattr(row, field + "_sha256"):
            raise ValueError("prepared table substituted a sampled BAM or BAI")
    metadata = read_document(check(sampling["outputs"][key + ":metadata"]), "metadata_id")
    fixed = {"schema": "fp-tools-nested-depth-matrix-v1", "source_bam_sha256": source["bam"]["sha256"],
             "output_bam_sha256": row.bam_sha256, "output_bai_sha256": row.bai_sha256,
             "available_fragments": counts["available_fragments"], "target_fragments": depth, "seed": row.seed,
             "sampling_fraction": depth / counts["available_fragments"],
             "selection_rule": "blake2b(seed,query_name)<sampling_fraction", "nested_with_same_seed": True}
    if any(metadata.get(key) != value for key, value in fixed.items()):
        raise ValueError("sampled BAM metadata differs from its source, count or requested depth/seed")


def bind_correction_command(stage, preparation, manifest_path):
    prefix = f"{stage['sample']}.{stage['depth'] // 1_000_000}m.s{stage['seed']}"
    replacements = {"{bam}": stage["input_bam"]["path"], "{prefix}": prefix, "{outdir}": str(manifest_path.parent)}
    command = [replacements.get(value, value) for value in preparation["correction_template"]]
    if stage.get("command") != command:
        raise ValueError("prepared correction changed the frozen DWM command or its inputs")


def freeze_prepared(plan, root, check, records, old_config, old_result, old_frame, archives):
    config_path = check(plan["preparation_configuration"])
    prepared_path = data_path(plan["prepared_inputs"]["path"], root)
    prepared = read_document(check(record(prepared_path)), "stage_id")
    preparation = read_document(config_path, "configuration_id")
    if (preparation["configuration_id"] != plan["prepared_inputs"]["configuration_id"]
            or prepared.get("configuration_id") != preparation["configuration_id"]
            or prepared.get("schema") != "fp-tools-direct-bias-replicate-depth-inputs-v1"
            or prepared.get("labels_scored") is not False or prepared.get("candidate_refitted") is not False
            or len(prepared["correction_manifests"]) != 60):
        raise ValueError("prerequisite is not the complete fixed-candidate 60-run preparation")
    for entry in ([preparation["plan_file"]] + preparation["sources"] + preparation["native_extensions"]
                  + list(preparation["common"].values())):
        check(entry)
    for name, digest in preparation["code"].items():
        check({"path": str(REPOSITORY / name), "sha256": digest})
    new_frame = pd.read_csv(check(prepared["outputs"]["signals"]), sep="\t")
    matrix = validate_merge(old_frame, new_frame)
    counts, sampling = {}, {}
    for sample, source in preparation["samples"].items():
        for name in ("bam", "bai"):
            check(source[name])
        sample_path = config_path.parent / sample / "sampling.json"
        sampling[sample] = preparation_stage(record(sample_path), preparation, sample, "sampling", check)
        counts[sample] = preparation_stage(sampling[sample]["count_manifest"], preparation, sample, "fragment_count", check)
        count = counts[sample]
        if (count["source_bam"]["sha256"] != source["bam"]["sha256"] or count["available_fragments"] < 50_000_000
                or count["count_rule"] != preparation["plan"]["fragment_count_rule"]):
            raise ValueError("sampled BAMs are not bound to the frozen source fragment count")
    corrections = []
    for entry in prepared["correction_manifests"]:
        manifest_path = check(entry)
        stage = read_document(manifest_path, "stage_id")
        bind_correction_command(stage, preparation, manifest_path)
        corrections.append(preparation_stage(entry, preparation, stage["sample"], "correction", check,
                                             stage["depth"], stage["seed"]))
    bind_prepared_rows(new_frame, corrections)
    for row in matrix.itertuples(index=False):
        for name in INPUT_FIELDS:
            check({"path": getattr(row, name), "sha256": getattr(row, name + "_sha256")})
    for row in new_frame.itertuples(index=False):
        bind_sampling_row(row, sampling[row.sample], counts[row.sample], preparation["samples"][row.sample], check)
    selection = read_document(check(old_config["settings"]["matched_cohorts"]), "selection_id")
    matched = pd.read_csv(check(selection["outputs"]["cohorts"]), sep="\t")
    return matrix, matched, old_config, old_result, archives, sorted(records.values(), key=lambda entry: entry["path"])


def cohort_identity(sites, cohorts):
    return canonical_id({"sites": sites[base.SITE_KEYS + ["chip_label"]].to_dict("records"),
                         "cohorts": {name: np.asarray(index, dtype=int).tolist() for name, index in cohorts.items()}})


def archived_cohorts(sites, matched, support, root):
    coverage = pd.read_csv(checked_path(support["outputs"]["coverage"], root), sep="\t")
    order = base.align_sites(coverage, sites)
    if len(coverage) != len(sites):
        raise ValueError("archived coverage is not the complete original reference support")
    mask = coverage.iloc[order].all_run_dwm_present.to_numpy(dtype=bool)
    selected, cohorts = base.qualify_support(sites, matched, mask)
    return sites.iloc[selected].reset_index(drop=True), cohorts


def freeze_supports(matrix, matched, source_config, source_result, outdir, configuration_id, root):
    """Complete both support stages and their joint identity before any run."""
    contexts, supports, memberships = {}, {}, {}
    for cell, frame in matrix.groupby("cell", sort=True):
        frame = frame.sort_values(base.RUN_KEYS).reset_index(drop=True)
        reference = source_config["references"][cell]
        sites, arrays, fixed_dwm, _ = base.load_reference(cell, reference, matched, root)
        old_support = base.verify_stage(checked_path(source_result["supports"][cell], root),
                                        source_config["configuration_id"], root)
        strict = matched.loc[matched.cell.eq(cell) & matched.policy.eq("motif_accessibility_observed_total")]
        old_sites, old_cohorts = archived_cohorts(sites, strict, old_support, root)
        directory = outdir / cell
        directory.mkdir(parents=True, exist_ok=True)
        sites, arrays, tracks, cohorts = base.freeze_support(cell, frame, sites, arrays, fixed_dwm, matched,
                                                            directory, configuration_id, root)
        base.align_sites(old_sites, sites)  # The all-96 intersection cannot gain old unsupported sites.
        old_identity, new_identity = cohort_identity(old_sites, old_cohorts), cohort_identity(sites, cohorts)
        memberships[cell] = {"site_and_cohort_identity": new_identity, "archived_identity": old_identity,
                             "archived_metrics_reusable": old_identity == new_identity,
                             "site_count": len(sites), "matched_site_count": len(cohorts["motif_accessibility_observed_total"])}
        contexts[cell] = (frame, sites, arrays, tracks, cohorts, old_sites)
        supports[cell] = record(directory / "support.json")
    if set(supports) != {"HepG2", "K562"}:
        raise ValueError("both development supports must be frozen before any run")
    joint = {"configuration_id": configuration_id, "supports": supports, "memberships": memberships,
             "new_run_extraction_started": False, "new_run_scoring_started": False,
             "selection_uses_new_scores": False}
    joint["stage_id"] = canonical_id(joint)
    write_immutable(outdir / "cohorts.freeze.json", joint)
    return contexts, joint


def subset_archived_profiles(archive, old_sites, sites, sample_dwm, root):
    extraction = read_document(checked_path(archive["extraction"], root), "stage_id")
    order = base.align_sites(old_sites, sites)
    with np.load(checked_path(extraction["outputs"]["profiles"], root), allow_pickle=False) as data:
        if not np.array_equal(data["site_hash"], base.site_hashes(old_sites)):
            raise ValueError("archived extraction does not have the frozen old site order")
        observed, expected, archived_dwm = (data[name][order] for name in ("observed", "expected", "sample_dwm"))
    if not np.array_equal(archived_dwm, sample_dwm):
        raise ValueError("archived DWM profiles differ from the checksum-bound current track")
    return observed, expected, archived_dwm


def reuse_extraction(row, archive, old_sites, sites, sample_dwm, outdir, configuration_id, root):
    path = outdir / "extraction.json"
    run = {key: getattr(row, key) for key in base.RUN_KEYS}
    if path.exists():
        stage = base.verify_stage(path, configuration_id, root, expected_run=run)
        if stage.get("reused_from") != archive["extraction"]:
            raise ValueError("resumed profiles have a different archived source")
        with np.load(checked_path(stage["outputs"]["profiles"], root), allow_pickle=False) as data:
            if not np.array_equal(data["site_hash"], base.site_hashes(sites)):
                raise ValueError("resumed archived extraction has different sites")
            if not np.array_equal(data["sample_dwm"], sample_dwm):
                raise ValueError("resumed archived extraction has different DWM profiles")
            return data["observed"], data["expected"], data["sample_dwm"]
    observed, expected, sample_dwm = subset_archived_profiles(archive, old_sites, sites, sample_dwm, root)
    target = outdir / "profiles.npz"
    np.savez_compressed(target, observed=observed, expected=expected, sample_dwm=sample_dwm,
                        site_hash=base.site_hashes(sites))
    stage = {"configuration_id": configuration_id, "run": run,
             "site_identity": canonical_id(sites[base.SITE_KEYS].to_dict("records")),
             "reused_from": archive["extraction"], "outputs": {"profiles": record(target)}}
    stage["stage_id"] = canonical_id(stage)
    write_immutable(path, stage)
    print(f"Reused archived profiles: {row.sample} {row.depth} {row.seed}", flush=True)
    return observed, expected, sample_dwm


def reuse_evaluation(row, archive, membership, source_configuration_id, outdir, configuration_id, root):
    if (not membership["archived_metrics_reusable"]
            or membership["site_and_cohort_identity"] != membership["archived_identity"]):
        raise ValueError("archived metrics cannot be reused after a site or cohort change")
    run = {key: getattr(row, key) for key in base.RUN_KEYS}
    previous = base.verify_stage(checked_path(archive["evaluation"], root), source_configuration_id, root, expected_run=run)
    stage = {"configuration_id": configuration_id, "run": run, "outputs": previous["outputs"],
             "bootstrap_run": previous["bootstrap_run"], "reused_from": archive["evaluation"],
             "site_and_cohort_identity": membership["site_and_cohort_identity"]}
    stage["stage_id"] = canonical_id(stage)
    write_immutable(outdir / "evaluation.json", stage)
    return stage


@contextmanager
def exclusive_run(outdir):
    import fcntl
    outdir.mkdir(parents=True, exist_ok=True)
    with (outdir / ".evaluation.lock").open("a+") as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise ValueError("another 96-run evaluator holds this output directory") from exc
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def runtime():
    return {"python": sys.version, "thread_environment": {name: os.environ.get(name) for name in THREADS},
            **{name: version(name) for name in ("numpy", "pandas", "scipy", "scikit-learn", "pysam", "pyBigWig", "matplotlib")}}


def execute(plan, plan_path, root, outdir):
    if any(os.environ.get(name) != value for name, value in THREADS.items()):
        raise ValueError("the frozen evaluation requires one BLAS, OMP and MKL thread")
    if checked_path(plan["evaluator_code"], root).resolve() != Path(__file__).resolve():
        raise ValueError("plan must pin this evaluator's own code")
    matrix, matched, source_config, source_result, archives, inputs = freeze_inputs(plan, root)
    if runtime() != source_config["runtime"]:
        raise ValueError("evaluation runtime differs from the archived 36-run scientific calculation")
    preparation = read_document(data_path(plan["preparation_configuration"]["path"], root), "configuration_id")
    code = source_config["code"] + [{"path": str(REPOSITORY / name), "sha256": digest}
                                   for name, digest in preparation["code"].items()] + [record(Path(__file__))]
    code_paths = {entry["path"] for entry in code}
    config = {"schema": SCHEMA, "settings": plan, "plan": record(plan_path),
              "input_records": [entry for entry in inputs if entry["path"] not in code_paths],
              "scientific_settings": source_config["settings"], "matrix": matrix.to_dict("records"),
              "code": code, "runtime": runtime(),
              "source_configuration_id": source_config["configuration_id"], "archives": archives}
    config["configuration_id"] = canonical_id(config)
    write_immutable(outdir / "configuration.freeze.json", config)
    print(f"Froze 96-run configuration: {config['configuration_id']}", flush=True)
    contexts, joint = freeze_supports(matrix, matched, source_config, source_result, outdir, config["configuration_id"], root)
    archive_map = {tuple(item["run"][key] for key in base.RUN_KEYS): item for item in archives}
    completed = []
    for cell, (frame, sites, arrays, tracks, cohorts, old_sites) in contexts.items():
        for index, row in enumerate(frame.itertuples(index=False)):
            directory = outdir / row.sample / row.depth / f"seed_{row.seed}"
            directory.mkdir(parents=True, exist_ok=True)
            run = {key: getattr(row, key) for key in base.RUN_KEYS}
            archive = archive_map.get(run_key(row))
            if (directory / "evaluation.json").exists():
                base.verify_stage(directory / "extraction.json", config["configuration_id"], root, expected_run=run)
                base.verify_stage(directory / "evaluation.json", config["configuration_id"], root, expected_run=run)
            else:
                if archive:
                    observed, expected, sample_dwm = reuse_extraction(row, archive, old_sites, sites, tracks[f"run_{index}"],
                                                                     directory, config["configuration_id"], root)
                else:
                    observed, expected, sample_dwm = base.extract_run(row, sites, arrays, source_config["references"][cell],
                                                                      tracks[f"run_{index}"], directory, config["configuration_id"], root)
                if archive and joint["memberships"][cell]["archived_metrics_reusable"]:
                    reuse_evaluation(row, archive, joint["memberships"][cell], source_config["configuration_id"],
                                     directory, config["configuration_id"], root)
                else:
                    methods = base.score_panel(observed, expected, tracks["fixed_dwm"], sample_dwm,
                                               source_config["settings"]["dispersion"][cell])
                    base.evaluate_run(row, sites, cohorts, methods, observed, source_config["settings"],
                                      directory, config["configuration_id"], root)
            completed.append({name: record(directory / (name + ".json")) for name in ("extraction", "evaluation")})
    return finish(config, joint, completed, outdir, root)


def render_stability(metrics, profiles, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    methods = {base.DWM: ("#777777", "Fixed DWM"), base.SAMPLE_DWM: ("#7B3294", "Sample DWM"),
               base.PARAMETRIC: ("#D55E00", "Parametric")}
    with PdfPages(output, metadata={"CreationDate": None, "ModDate": None}) as pdf:
        for (cell, tf), task in profiles.groupby(["cell", "tf"], sort=True):
            fig, axes = plt.subplots(2, 4, figsize=(11.69, 8.27))
            fig.suptitle(f"{cell} {tf}: complete development replicate/depth stability\n"
                         "Fixed candidate and geometry; support frozen across all 96 runs", fontsize=13)
            for row_index, cohort in enumerate(("original_support_with_expected_coverage", "motif_accessibility_observed_total")):
                for column, depth in enumerate(("10m", "25m", "50m", "full")):
                    ax = axes[row_index, column]
                    panel = task.loc[task.cohort.eq(cohort) & task.depth.eq(depth)]
                    if panel.empty:
                        ax.set_axis_off()
                        continue
                    for method, (color, label) in methods.items():
                        data = panel.loc[panel.method.eq(method)]
                        for _, group in data.groupby(["sample", "seed"]):
                            ax.plot(group.position, group.positive_minus_negative, color=color, alpha=.10, linewidth=.6)
                        mean = data.groupby("position").positive_minus_negative.mean()
                        ax.plot(mean.index, mean.values, color=color, linewidth=1.1, label=label)
                    selected = metrics.loc[metrics.cell.eq(cell) & metrics.tf.eq(tf) & metrics.cohort.eq(cohort)
                                           & metrics.depth.eq(depth) & metrics.method.eq(base.PARAMETRIC)]
                    ax.set(title=f"{depth}: {'coverage subset' if row_index == 0 else 'matched'}; {selected.status.iloc[0]}\n"
                                 f"ΔAUROC vs raw {selected.auroc_gain_over_raw.min():+.3f} to {selected.auroc_gain_over_raw.max():+.3f}",
                           xlabel="Motif offset (bp)", xlim=(-100, 100))
                    ax.title.set_fontsize(8)
                    ax.tick_params(labelsize=7)
                    ax.axhline(0, color="#888888", linewidth=.5)
                    ax.axvspan(-5, 5, alpha=.15, color="#aaaaaa")
                    ax.spines[["top", "right"]].set_visible(False)
                    if column == 0:
                        ax.set_ylabel("Positive − negative\n(normalized deviance)", fontsize=8)
            axes[0, 0].legend(fontsize=7, frameon=False)
            fig.text(.065, .045,
                     "Positive fixed and sample DWM expectation in every run; complete matched pairs retained. Thin curves: runs; thick curves: their mean.\n"
                     "Three biological replicates at every depth; five nested seeds per replicate at 10M/25M/50M and three deterministic full-depth endpoints.\n"
                     "Examined chr19–22/X. Per-run pointwise bands are in tables. Seed means are not independent replication or promotion evidence; raw scores are a guardrail.",
                     fontsize=7.5)
            fig.subplots_adjust(left=.07, right=.98, top=.82, bottom=.16, wspace=.35, hspace=.5)
            pdf.savefig(fig)
            plt.close(fig)


def verify_completed_stages(stages, config, joint, root):
    if len(stages) != 96:
        raise ValueError("completed result lacks the exact 96-run matrix")
    archive_map = {tuple(item["run"][key] for key in base.RUN_KEYS): item for item in config["archives"]}
    runs, reused_metrics = [], 0
    for item in stages:
        if set(item) != {"extraction", "evaluation"}:
            raise ValueError("each run requires extraction and evaluation manifests")
        extraction = base.verify_stage(checked_path(item["extraction"], root), config["configuration_id"], root)
        run = extraction["run"]
        evaluation = base.verify_stage(checked_path(item["evaluation"], root), config["configuration_id"], root, expected_run=run)
        archive = archive_map.get(tuple(run[key] for key in base.RUN_KEYS))
        if archive:
            if extraction.get("reused_from") != archive["extraction"]:
                raise ValueError("one of the 36 archived profiles was not reused from its frozen source")
            checked_path(extraction["reused_from"], root)
        elif "reused_from" in extraction:
            raise ValueError("a missing replicate profile incorrectly claims archived extraction")
        if "reused_from" in evaluation:
            membership = joint["memberships"][run["cell"]]
            if (not archive or evaluation["reused_from"] != archive["evaluation"]
                    or not membership["archived_metrics_reusable"]
                    or membership["archived_identity"] != membership["site_and_cohort_identity"]
                    or evaluation.get("site_and_cohort_identity") != membership["site_and_cohort_identity"]):
                raise ValueError("archived metrics were reused for a different site/cohort identity")
            previous = base.verify_stage(checked_path(evaluation["reused_from"], root), config["source_configuration_id"], root,
                                         expected_run=run)
            if evaluation["outputs"] != previous["outputs"]:
                raise ValueError("reused evaluation changed an archived output")
            reused_metrics += 1
        runs.append(run)
    validate_matrix(pd.DataFrame(runs))
    return reused_metrics


def finish(config, joint, completed, outdir, root):
    reused_metrics = verify_completed_stages(completed, config, joint, root)
    outputs, tables = {}, {}
    for name in ("metrics", "bootstrap", "profiles", "shapes", "coverage"):
        frames = []
        for item in completed:
            stage = read_document(Path(item["evaluation"]["path"]), "stage_id")
            frame = pd.read_csv(checked_path(stage["outputs"][name], root), sep="\t")
            if not frame.empty:
                frames.append(frame)
        tables[name] = pd.concat(frames, ignore_index=True) if frames else frame.iloc[:0].copy()
        target = outdir / f"stability_{name}.tsv.gz"
        tables[name].to_csv(target, sep="\t", index=False, compression={"method": "gzip", "mtime": 0})
        outputs[name] = record(target)
    pdf = outdir / "stability_aggregate_profiles.pdf"
    render_stability(tables["metrics"], tables["profiles"], pdf)
    outputs["pdf"] = record(pdf)
    done = {"schema": SCHEMA, "configuration_id": config["configuration_id"], "outputs": outputs,
            "stages": completed, "supports": joint["supports"], "cohorts": record(outdir / "cohorts.freeze.json"),
            "archived_extractions_reused": 36, "archived_evaluations_reused": reused_metrics,
            "new_holdout_labels_read": False, "final_untouched_control_read": False, "candidate_retuned": False,
            "promotion_eligible": False, "evidence": config["settings"]["evidence"]}
    done["stage_id"] = canonical_id(done)
    write_immutable(outdir / "stability_manifest.json", done)
    print(f"Completed 96 frozen development runs; reused 36 extractions and {reused_metrics} evaluations; "
          f"result {done['stage_id']}", flush=True)
    return done


def verify_result(outdir, root, code_ref=None):
    config = read_document(outdir / "configuration.freeze.json", "configuration_id")
    validate_plan(config["settings"])
    validate_matrix(pd.DataFrame(config["matrix"]))
    base.validate_matrix(pd.DataFrame([item["run"] for item in config["archives"]]))
    if code_ref is not None and not re.fullmatch(r"[0-9a-fA-F]{7,40}", code_ref):
        raise ValueError("--code-ref must be an immutable Git commit hash")
    for item in config["code"]:
        if code_ref is None:
            checked_path(item, root)
        else:
            relative = Path(item["path"]).relative_to(REPOSITORY).as_posix()
            blob = subprocess.check_output(["git", "show", f"{code_ref}:{relative}"], cwd=REPOSITORY)
            if sha256(blob).hexdigest() != item["sha256"]:
                raise ValueError("historical calculation code checksum mismatch")
    checked_path(config["plan"], root)
    for item in config["input_records"]:
        checked_path(item, root)
    source = read_document(data_path(config["settings"]["source_configuration"]["path"], root), "configuration_id")
    if config["scientific_settings"] != source["settings"] or config["source_configuration_id"] != source["configuration_id"]:
        raise ValueError("extension changed the archived candidate or scientific settings")
    result = base.verify_stage(outdir / "stability_manifest.json", config["configuration_id"], root)
    if (result.get("schema") != SCHEMA or result.get("evidence") != config["settings"]["evidence"]
            or any(result.get(key) is not False for key in
                   ("promotion_eligible", "candidate_retuned", "new_holdout_labels_read", "final_untouched_control_read"))
            or result.get("archived_extractions_reused") != 36):
        raise ValueError("unexpected completed development evidence scope")
    joint = read_document(checked_path(result["cohorts"], root), "stage_id")
    if (joint.get("configuration_id") != config["configuration_id"] or joint["supports"] != result["supports"]
            or set(joint["supports"]) != {"HepG2", "K562"}
            or any(joint.get(key) is not False for key in
                   ("new_run_extraction_started", "new_run_scoring_started", "selection_uses_new_scores"))):
        raise ValueError("completed result lacks the joint support freeze before extraction/scoring")
    for cell, entry in joint["supports"].items():
        support = base.verify_stage(checked_path(entry, root), config["configuration_id"], root)
        expected = pd.DataFrame(config["matrix"]).query("cell == @cell")[base.RUN_KEYS].to_dict("records")
        if support.get("cell") != cell or support["runs"] != expected:
            raise ValueError("completed support has a different cell or replicate/depth/seed matrix")
    reused_metrics = verify_completed_stages(result["stages"], config, joint, root)
    if reused_metrics != result["archived_evaluations_reused"]:
        raise ValueError("completed result misstates the number of reused evaluations")
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--plan-sha256", help="Required plan checksum when launched by a frozen queue")
    parser.add_argument("--verify-only", action="store_true")
    parser.add_argument("--code-ref", help="Immutable calculation commit, for verification only")
    args = parser.parse_args(argv)
    if args.code_ref is not None and not args.verify_only:
        parser.error("--code-ref requires --verify-only")
    if args.plan_sha256 is not None and file_sha256(args.plan) != args.plan_sha256:
        raise ValueError("queued plan checksum changed before evaluation")
    plan = json.loads(args.plan.read_text())
    validate_plan(plan)
    outdir = args.outdir.resolve()
    if args.verify_only or (outdir / "stability_manifest.json").is_file():
        frozen = read_document(outdir / "configuration.freeze.json", "configuration_id")
        if record(args.plan)["sha256"] != frozen["plan"]["sha256"]:
            raise ValueError("requested plan differs from the completed extension")
        verify_result(outdir, args.data_root, args.code_ref)
        print("Verified all 96 runs, frozen inputs, archived reuse and outputs; no stage rerun.")
    else:
        with exclusive_run(outdir):
            execute(plan, args.plan, args.data_root, outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
