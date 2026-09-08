#!/usr/bin/env python3
"""Prepare only the missing 60 development replicate/depth/seed inputs.

Freeze source BAMs and the recorded correction command before preparation.
Count usable fragments with the original rule, create nested pair-preserving
subsets, and fit the sample/depth-matched conventional DWM comparator. This
script cannot score labels, fit a candidate or open a final negative control.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from hashlib import sha256
from importlib.metadata import version
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import time

import pandas as pd
import pysam

REPOSITORY = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from downsample_bam_depth_matrix import build_targets, downsample_matrix, file_sha256
from downsample_bam_fragments import fragment_representative
from run_footprint_ablation_plan import output_is_complete

SCHEMA = "fp-tools-direct-bias-replicate-depth-inputs-v1"
SAMPLES = ["HepG2_rep2", "HepG2_rep3", "K562_rep2", "K562_rep3"]
CODE = ["benchmarks/scripts/prepare_direct_bias_replicate_depths.py",
        "benchmarks/scripts/downsample_bam_depth_matrix.py",
        "benchmarks/scripts/downsample_bam_fragments.py",
        "benchmarks/scripts/run_footprint_ablation_plan.py"]
THREADS = {name: "1" for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")}
COUNT_EXPRESSION = "(flag & 3844) == 0 && (!(flag & 1) || (flag & 64))"


def identity(document, field):
    return sha256(json.dumps({k: v for k, v in document.items() if k != field},
                             sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def read_json(path, field):
    document = json.loads(Path(path).read_text())
    if document.get(field) != identity(document, field):
        raise ValueError(f"invalid {field}: {path}")
    return document


def write_json(path, document, field):
    document = dict(document)
    document[field] = identity(document, field)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as handle:
        handle.write(json.dumps(document, sort_keys=True, indent=2, allow_nan=False) + "\n")
    return document


def checked(entry, root=Path("/")):
    path = Path(entry["path"])
    path = path if path.is_absolute() else root / path
    if file_sha256(path) != entry["sha256"]:
        raise ValueError(f"input/output checksum mismatch: {path}")
    return path.resolve()


def record(path):
    return {"path": str(Path(path).resolve()), "sha256": file_sha256(Path(path))}


def checked_record(entry, root):
    return {"path": str(checked(entry, root)), "sha256": entry["sha256"]}


def validate_plan(plan):
    fixed = {"schema": "fp-tools-direct-bias-replicate-depth-plan-v1",
             "scope": "missing_development_replicate_depth_inputs", "samples": SAMPLES,
             "depths": [10_000_000, 25_000_000, 50_000_000], "seeds": list(range(2026, 2031)),
             "sampling_rule": "blake2b(seed,query_name)<target_fragments/usable_fragments",
             "fragment_count_rule": "downsample_bam_fragments.fragment_representative",
             "fragment_count_implementation": "pinned_pysam_view_filter_verified_against_all_4096_SAM_flag_combinations",
             "depths_above_available": "reject_without_capping_or_relabeling_full_depth",
             "new_holdout_labels_allowed": False, "final_untouched_control_allowed": False,
             "candidate_retuning_allowed": False, "dwm_cores": 4, "maximum_parallel_jobs": 2,
             "candidate_model_sha256": "620f232f634423b8753e4d77f8cf72df87482c09193cb5e7e6b038917c47e89a",
             "dwm_comparator": "refit_on_each_sample_depth_seed_using_recorded_replicate_1_command",
             "evaluation_chromosomes": ["chr19", "chr20", "chr21", "chr22", "chrX"],
             "profile_extraction_minimum_mapq": 30, "profile_read_shift": [4, -5],
             "evaluation_policy": "freeze_96_run_expected_coverage_and_retain_complete_existing_matched_pairs_before_scoring"}
    if any(plan.get(key) != value for key, value in fixed.items()):
        raise ValueError("only the fixed missing development input matrix is allowed")


def source_rows(frame):
    selected = frame[frame["sample"].isin(SAMPLES)].copy()
    if (len(selected) != 4 or selected["sample"].duplicated().any()
            or set(selected["sample"]) != set(SAMPLES)
            or not selected.depth.eq("full").all() or not selected.seed.eq(2026).all()
            or any(row.cell != row.sample.rsplit("_", 1)[0] for row in selected.itertuples(index=False))):
        raise ValueError("source registry must have exactly the four full development replicates")
    return selected.sort_values("sample").reset_index(drop=True)


def correction_template(commands, root, plan, common):
    templates = []
    for cell in ("HepG2", "K562"):
        job = f"correct:{cell}_rep1.10m.s2026:fp_tools_dwm"
        match = commands[commands.job_id.eq(job)]
        if len(match) != 1:
            raise ValueError("reference correction command is missing or duplicated")
        args = shlex.split(match.iloc[0].command)
        if args[1:3] != ["-m", "fp_tools.cli"] or "--bias-pkl" in args:
            raise ValueError("reference command is not a refitted conventional DWM")
        args[0] = sys.executable
        for flag, name in (("--genome", "genome"), ("--peaks", "peaks"), ("--blacklist", "blacklist")):
            value = Path(args[args.index(flag)+1])
            value = value if value.is_absolute() else root / value
            if value.resolve() != Path(common[name]["path"]):
                raise ValueError("reference correction inputs differ from the frozen plan")
            args[args.index(flag)+1] = common[name]["path"]
        if (args[args.index("--score_mat")+1] != "DWM" or args[args.index("--write-tracks")+1] != "all"
                or int(args[args.index("--cores")+1]) != plan["dwm_cores"]):
            raise ValueError("reference correction model, output or CPU settings differ")
        for flag, token in (("--bams", "{bam}"), ("--prefix", "{prefix}"), ("--outdir", "{outdir}")):
            args[args.index(flag)+1] = token
        templates.append(args)
    if templates[0] != templates[1]:
        raise ValueError("development cells used different DWM correction options")
    return templates[0]


def freeze(plan, plan_path, root, outdir):
    registry_path = checked(plan["source_registry"], root)
    registry = read_json(registry_path, "depth_input_id")
    if registry.get("models_refitted_by_depth") is not False:
        raise ValueError("source registry does not identify the fixed candidate")
    table_path = checked(registry["signals"], root)
    rows = source_rows(pd.read_csv(table_path, sep="\t"))
    common = {name: checked_record(plan[name], root) for name in ("genome", "genome_index", "peaks", "blacklist")}
    commands_path = checked(plan["replicate_1_commands"], root)
    full_commands_path = checked(plan["full_replicate_commands"], root)
    template = correction_template(pd.read_csv(commands_path, sep="\t"), root, plan, common)
    samples = {}
    for row in rows.itertuples(index=False):
        samples[row.sample] = {"cell": row.cell}
        for name in ("bam", "bai"):
            samples[row.sample][name] = checked_record({"path": getattr(row, name), "sha256": getattr(row, name+"_sha256")}, root)
        with pysam.AlignmentFile(samples[row.sample]["bam"]["path"], "rb") as bam:
            if bam.header.to_dict().get("HD", {}).get("SO") != "coordinate" or not bam.has_index():
                raise ValueError("source must be coordinate sorted with an index")
        print(f"Verified frozen source: {row.sample}", flush=True)
    subprocess.run(["git", "diff", "--quiet", "HEAD", "--", "src/fp_tools"], cwd=REPOSITORY, check=True)
    outdir.mkdir(parents=True, exist_ok=False)
    return write_json(outdir / "configuration.freeze.json", {
        "schema": SCHEMA, "plan": plan, "plan_file": record(plan_path), "samples": samples,
        "common": common, "correction_template": template,
        "sources": [record(p) for p in (registry_path, table_path, commands_path, full_commands_path)],
        "code": {name: file_sha256(REPOSITORY / name) for name in CODE},
        "source_tree": subprocess.check_output(["git", "rev-parse", "HEAD:src/fp_tools"], cwd=REPOSITORY, text=True).strip(),
        "native_extensions": [record(path) for path in sorted((REPOSITORY / "src/fp_tools").rglob("*.so"))],
        "python": sys.version, "samtools_version": pysam.__samtools_version__,
        "versions": {name: version(name) for name in ("pysam", "numpy", "pandas", "scipy", "matplotlib")},
        "threads": THREADS, "labels_scored": False,
    }, "configuration_id")


def verify_configuration(config, plan_path, code_ref=None):
    if file_sha256(plan_path) != config["plan_file"]["sha256"]:
        raise ValueError("plan differs from the frozen input configuration")
    for name, expected in config["code"].items():
        actual = (sha256(subprocess.check_output(["git", "show", f"{code_ref}:{name}"], cwd=REPOSITORY)).hexdigest()
                  if code_ref else file_sha256(REPOSITORY / name))
        if actual != expected:
            raise ValueError(f"frozen preparation code changed: {name}")
    ref = code_ref or "HEAD"
    if subprocess.check_output(["git", "rev-parse", f"{ref}:src/fp_tools"], cwd=REPOSITORY, text=True).strip() != config["source_tree"]:
        raise ValueError("frozen production source tree changed")
    if not code_ref:
        subprocess.run(["git", "diff", "--quiet", "HEAD", "--", "src/fp_tools"], cwd=REPOSITORY, check=True)
        if (sys.version != config["python"] or pysam.__samtools_version__ != config["samtools_version"]
                or any(version(name) != value for name, value in config["versions"].items())):
            raise ValueError("preparation runtime differs from the frozen environment")
    for entry in config["sources"] + config["native_extensions"]:
        checked(entry)


def verify_stage(path, config, sample, stage, *, depth=None, seed=None):
    done = read_json(path, "stage_id")
    if (done.get("configuration_id") != config["configuration_id"] or done.get("sample") != sample
            or done.get("stage") != stage or done.get("depth") != depth or done.get("seed") != seed):
        raise ValueError("completed input stage belongs to a different frozen request")
    for entry in done.get("outputs", {}).values():
        checked(entry)
    if stage == "sampling":
        count_path = checked(done["count_manifest"])
        verify_stage(count_path, config, sample, "fragment_count")
    return done


@contextmanager
def exclusive_run(outdir):
    import fcntl
    with (outdir / ".preparation.lock").open("a+") as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError("another input preparation process holds this output lock") from exc
        handle.seek(0)
        handle.truncate()
        handle.write(str(os.getpid()) + "\n")
        handle.flush()
        yield


def native_fragment_count(path):
    # Isolate pysam's samtools dispatcher because stdout capture is process-wide.
    # This expression is the original fragment_representative flag predicate;
    # exhaustive flag fixtures verify it, including unpaired reads and MAPQ 0.
    code = "import pysam, sys; print(pysam.view('-c', '-e', sys.argv[1], sys.argv[2]), end='')"
    return int(subprocess.check_output([sys.executable, "-c", code, COUNT_EXPRESSION, str(path)],
                                       text=True, env={**os.environ, **THREADS}).strip())


def count_sample(config, sample, outdir):
    path = outdir / sample / "fragment_count.json"
    bam_path = checked(config["samples"][sample]["bam"])
    checked(config["samples"][sample]["bai"])
    if path.exists():
        return verify_stage(path, config, sample, "fragment_count")
    print(f"Counting usable fragments: {sample}", flush=True)
    started = time.perf_counter()
    fragments = native_fragment_count(bam_path)
    if fragments < max(config["plan"]["depths"]):
        raise ValueError(f"{sample} has only {fragments} usable fragments; cannot construct the declared depths")
    result = write_json(path, {"schema": SCHEMA, "configuration_id": config["configuration_id"],
                              "sample": sample, "stage": "fragment_count",
                              "available_fragments": fragments, "source_bam": config["samples"][sample]["bam"],
                              "count_rule": config["plan"]["fragment_count_rule"], "filter_expression": COUNT_EXPRESSION,
                              "wall_seconds": time.perf_counter() - started}, "stage_id")
    print(f"Count complete: {sample} {fragments:,} fragments", flush=True)
    return result


def sample_matrix(config, sample, outdir):
    manifest = outdir / sample / "sampling.json"
    count_path = outdir / sample / "fragment_count.json"
    counts = verify_stage(count_path, config, sample, "fragment_count")
    if manifest.exists():
        return verify_stage(manifest, config, sample, "sampling")
    targets = build_targets(outdir, sample, config["plan"]["depths"], config["plan"]["seeds"])
    print(f"Creating missing depth matrix: {sample}", flush=True)
    rows = downsample_matrix(config["samples"][sample]["bam"]["path"], targets,
                             available_fragments=counts["available_fragments"], index_threads=2,
                             expected_source_sha256=config["samples"][sample]["bam"]["sha256"])
    outputs = {}
    for target in targets:
        for kind, path in (("bam", target.output_bam), ("bai", Path(str(target.output_bam)+".bai")),
                           ("metadata", target.metadata_path)):
            outputs[f"{target.target_fragments}:{target.seed}:{kind}"] = record(path)
    return write_json(manifest, {"schema": SCHEMA, "configuration_id": config["configuration_id"],
                                 "sample": sample, "stage": "sampling", "count_manifest": record(count_path),
                                 "outputs": outputs, "records": rows}, "stage_id")


def correct_run(config, sample, depth, seed, outdir):
    label = f"{depth // 1_000_000}m"
    prefix = f"{sample}.{label}.s{seed}"
    directory = outdir / "signals" / sample / label / f"seed_{seed}" / "fp_tools_dwm"
    manifest = directory / "correction.json"
    if manifest.exists():
        return verify_stage(manifest, config, sample, "correction", depth=depth, seed=seed)
    sampling = read_json(outdir / sample / "sampling.json", "stage_id")
    if sampling["configuration_id"] != config["configuration_id"] or sampling["sample"] != sample:
        raise ValueError("sampling manifest differs from correction inputs")
    key = f"{depth}:{seed}"
    bam = checked(sampling["outputs"][key+":bam"])
    checked(sampling["outputs"][key+":bai"])
    if directory.exists():
        raise ValueError(f"partial correction requires review before resume: {directory}")
    directory.mkdir(parents=True)
    substitutions = {"{bam}": str(bam), "{prefix}": prefix, "{outdir}": str(directory)}
    command = [substitutions.get(value, value) for value in config["correction_template"]]
    print(f"Correcting conventional DWM: {prefix}", flush=True)
    start = time.perf_counter()
    with (directory / "command.log").open("x") as log:
        subprocess.run(command, cwd=REPOSITORY, env={**os.environ, **THREADS, "PYTHONPATH": str(REPOSITORY / "src")},
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    elapsed = time.perf_counter() - start
    paths = [directory / f"{prefix}_{track}.bw" for track in ("uncorrected", "bias", "expected", "corrected")]
    if not all(output_is_complete(path) for path in paths):
        raise ValueError("correction returned without complete nonempty tracks")
    paths += [directory / f"{prefix}_AtacBias.pickle", directory / f"{prefix}_atacorrect.pdf", directory / "command.log"]
    result = write_json(manifest, {"schema": SCHEMA, "configuration_id": config["configuration_id"],
                                  "sample": sample, "stage": "correction", "depth": depth, "seed": seed,
                                  "command": command, "wall_seconds": elapsed,
                                  "input_bam": sampling["outputs"][key+":bam"],
                                  "outputs": {path.name: record(path) for path in paths}}, "stage_id")
    print(f"Correction complete: {prefix} ({elapsed:.1f}s)", flush=True)
    return result


def run(config, outdir, phase):
    samples, workers = config["plan"]["samples"], config["plan"]["maximum_parallel_jobs"]
    if phase in ("count", "all"):
        with ThreadPoolExecutor(max_workers=workers) as pool:
            list(pool.map(lambda sample: count_sample(config, sample, outdir), samples))
    if phase == "count":
        return
    if phase in ("sample", "all"):
        with ThreadPoolExecutor(max_workers=workers) as pool:
            list(pool.map(lambda sample: sample_matrix(config, sample, outdir), samples))
    if phase == "sample":
        return
    for entry in config["common"].values():
        checked(entry)
    jobs = [(sample, depth, seed) for depth in config["plan"]["depths"] for seed in config["plan"]["seeds"] for sample in samples]
    if phase != "verify":
        with ThreadPoolExecutor(max_workers=workers) as pool:
            list(pool.map(lambda job: correct_run(config, *job, outdir), jobs))
    rows, manifests = [], []
    for sample in samples:
        checked(config["samples"][sample]["bam"])
        checked(config["samples"][sample]["bai"])
        verify_stage(outdir / sample / "fragment_count.json", config, sample, "fragment_count")
        verify_stage(outdir / sample / "sampling.json", config, sample, "sampling")
    for sample, depth, seed in jobs:
        label = f"{depth // 1_000_000}m"
        prefix = f"{sample}.{label}.s{seed}"
        directory = outdir / "signals" / sample / label / f"seed_{seed}"
        manifest = directory / "fp_tools_dwm/correction.json"
        done = verify_stage(manifest, config, sample, "correction", depth=depth, seed=seed)
        manifests.append(record(manifest))
        row = {"sample": sample, "cell": config["samples"][sample]["cell"], "depth": label, "seed": seed}
        for name, path in (("bam", directory / f"{prefix}.bam"), ("bai", directory / f"{prefix}.bam.bai"),
                           ("dwm_expected", directory / "fp_tools_dwm" / f"{prefix}_expected.bw")):
            entry = record(path)
            row[name], row[name+"_sha256"] = entry["path"], entry["sha256"]
        rows.append(row)
    final = outdir / "prepared_inputs.json"
    if final.exists():
        done = read_json(final, "stage_id")
        if done["configuration_id"] != config["configuration_id"] or done["correction_manifests"] != manifests:
            raise ValueError("prepared input registry differs from the complete verified matrix")
        checked(done["outputs"]["signals"])
    elif phase == "verify":
        raise ValueError("final prepared input registry is missing")
    else:
        table = outdir / "depth_inputs.tsv"
        pd.DataFrame(rows).to_csv(table, sep="\t", index=False)
        write_json(final, {"schema": SCHEMA, "configuration_id": config["configuration_id"],
                           "correction_manifests": manifests, "outputs": {"signals": record(table)},
                           "labels_scored": False, "candidate_refitted": False}, "stage_id")
    print("Verified all 60 prepared input runs.", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--phase", choices=["freeze", "count", "sample", "correct", "all", "verify"], default="all")
    parser.add_argument("--code-ref")
    args = parser.parse_args()
    if args.code_ref and args.phase != "verify":
        parser.error("--code-ref is available only for verification")
    plan = json.loads(args.plan.read_text())
    validate_plan(plan)
    args.outdir = args.outdir.resolve()
    config_path = args.outdir / "configuration.freeze.json"
    if config_path.exists():
        config = read_json(config_path, "configuration_id")
        verify_configuration(config, args.plan, args.code_ref)
    elif args.phase not in ("freeze", "all"):
        raise ValueError("configuration must be frozen before any preparation stage")
    else:
        config = freeze(plan, args.plan, args.data_root, args.outdir)
    if args.phase == "freeze":
        print("Frozen missing replicate/depth input preparation: " + config["configuration_id"])
        return
    with exclusive_run(args.outdir):
        run(config, args.outdir, args.phase)


if __name__ == "__main__":
    main()
