#!/usr/bin/env python3
"""Measure frozen sequence inference in fresh, single-CPU Linux processes.

This is a component diagnostic, not an end-to-end computational promotion test.
Both methods score the same 201 genomic positions on both strands. The DWM is
the archived production Cython scorer, not the slower Python control reference.
No labels, cut counts, model fitting, correction or detector selection are used.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
from importlib.metadata import version
import json
import os
from pathlib import Path
import pickle
import platform
import subprocess
import sys
import time

import numpy as np
import pandas as pd

REPOSITORY = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

# Keep both workers' import baseline equal, including the archived pickle class.
from fp_tools.tools.atacorrect_functions import AtacBias  # noqa: F401,E402
from fp_tools.tools.parametric_bias import ConditionalSequenceBiasModel, strand_log_bias  # noqa: E402
from fp_tools.utils.fasta import open_fasta  # noqa: E402
from fp_tools.utils import sequences as native_sequences  # noqa: E402

SCHEMA = "fp-tools-direct-bias-inference-v1"
METHODS = ("production_dwm_sequence", "frozen_parametric_sequence")
KEYS = ["cell", "tf", "TFBS_chr", "TFBS_start", "TFBS_end", "TFBS_strand"]
CHROMOSOMES = ["chr19", "chr20", "chr21", "chr22", "chrX"]
CODE = ["benchmarks/scripts/benchmark_direct_bias_inference.py",
        "src/fp_tools/tools/parametric_bias.py", "src/fp_tools/tools/atacorrect_functions.py",
        "src/fp_tools/utils/sequences.pyx", "src/fp_tools/utils/fasta.py"]
THREADS = {key: "1" for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                               "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")}


def digest_file(path):
    result = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def identity(document, field):
    return sha256(json.dumps({k: v for k, v in document.items() if k != field},
                             sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def checked(entry, root):
    path = Path(entry["path"])
    if not path.is_absolute():
        path = (REPOSITORY if entry.get("root") == "research" else root) / path
    if digest_file(path) != entry["sha256"]:
        raise ValueError(f"checksum mismatch: {path}")
    return path.resolve()


def record(path):
    return {"path": str(Path(path).resolve()), "sha256": digest_file(path)}


def write_json(path, document, field):
    document = dict(document)
    document[field] = identity(document, field)
    with Path(path).open("x") as handle:
        handle.write(json.dumps(document, indent=2, sort_keys=True, allow_nan=False) + "\n")
    return document


def read_json(path, field):
    document = json.loads(Path(path).read_text())
    if document.get(field) != identity(document, field):
        raise ValueError(f"invalid {field}: {path}")
    return document


def validate_plan(plan):
    fixed = {"schema": "fp-tools-direct-bias-inference-plan-v1",
             "scope": "sequence_inference_component_only", "cells": ["HepG2", "K562"],
             "chromosomes": CHROMOSOMES, "sites_per_cell": 512, "flank": 100, "margin": 41,
             "repetitions": 6, "warmup_sites": 16, "seed": 2026,
             "selection": "lowest_coordinate_sha256_among_all_run_expected_coverage",
             "methods": list(METHODS), "process_order": "alternate_method_order_each_repetition",
             "model_training_allowed": False, "new_holdout_labels_allowed": False,
             "final_untouched_control_allowed": False, "promotion_evidence": False,
             "threads": 1, "cpu_affinity_index": 0}
    if any(plan.get(k) != v for k, v in fixed.items()):
        raise ValueError("only the frozen development component diagnostic is supported")


def select_sites(frame, count, seed):
    """Selection consumes coordinates and an already frozen coverage flag only."""
    if frame.duplicated(KEYS).any() or not set(frame.TFBS_chr).issubset(CHROMOSOMES):
        raise ValueError("duplicate sites or chromosomes outside examined development support")
    selected = frame.loc[frame.all_run_dwm_present.eq(True), KEYS].copy()
    selected["selection_hash"] = [sha256(json.dumps([seed, *row], separators=(",", ":")).encode()).hexdigest()
                                  for row in selected.itertuples(index=False, name=None)]
    if len(selected) < count:
        raise ValueError("insufficient frozen coverage-qualified sites")
    return selected.sort_values("selection_hash", kind="stable").head(count).reset_index(drop=True)


def sequence_scores(method, model, sequences, width, margin):
    output = np.empty((len(sequences), 2, width), dtype=np.float64)
    for index, sequence in enumerate(sequences):
        if method == METHODS[1]:
            plus, minus, _ = strand_log_bias(model, str(sequence), margin + np.arange(width))
        elif method == METHODS[0]:
            half = int(model.bias["forward"].length) // 2
            if model.bias["reverse"].length != 2 * half + 1 or half > margin:
                raise ValueError("DWM strands have incompatible context lengths")
            # Trim to exactly width output positions, avoiding extra DWM work.
            codes = native_sequences.nuc_to_num(str(sequence)[margin-half:margin+width+half])
            reverse = np.asarray([1, 0, 3, 2, 4], dtype=np.int64)[codes[::-1]].copy()
            plus = model.bias["forward"].score_sequence(codes)[half:half+width]
            minus = model.bias["reverse"].score_sequence(reverse)[half:half+width][::-1]
        else:
            raise ValueError("unknown sequence inference method")
        output[index, 0], output[index, 1] = plus, minus
    output[~np.isfinite(output)] = np.nan
    return output


def array_identity(values):
    values = np.asarray(values, dtype="<f8").copy()
    values[~np.isfinite(values)] = np.nan
    return sha256(str(values.shape).encode() + values.tobytes()).hexdigest()


def peak_rss_mib():
    # VmHWM belongs to this exec's address space; never subtract high-water marks
    # or use a parent's cumulative RUSAGE_CHILDREN across multiple workers.
    fields = Path("/proc/self/status").read_text().splitlines()
    for line in fields:
        if line.startswith("VmHWM:"):
            return int(line.split()[1]) / 1024.0
    raise RuntimeError("Linux VmHWM is required for the declared measurement")


def worker(config_path, cell, method, repetition, outdir):
    config = read_json(config_path, "configuration_id")
    plan = config["plan"]
    validate_plan(plan)
    if cell not in plan["cells"] or method not in METHODS or repetition not in range(plan["repetitions"]):
        raise ValueError("worker identity differs from the frozen methods or cells")
    os.sched_setaffinity(0, {config["runtime"]["cpu"]})
    if any(os.environ.get(k) != v for k, v in THREADS.items()):
        raise ValueError("worker thread environment differs from the frozen single-CPU plan")
    inputs = config["cells"][cell]
    root = Path("/")
    bundle = checked(inputs["bundle"], root)
    with np.load(bundle, allow_pickle=False) as values:
        sequences = values["sequences"].copy()
    started = time.perf_counter()
    if method == METHODS[0]:
        # This checksum-pinned artifact was generated locally by the recorded
        # production run. Do not accept downloaded/unreviewed pickle inputs.
        with checked(inputs["dwm_model"], root).open("rb") as handle:
            model = pickle.load(handle)
        if not isinstance(model, AtacBias) or model.stype != "DWM":
            raise ValueError("archived model is not the production DWM")
    else:
        checked(config["parametric_sidecar"], root)
        model = ConditionalSequenceBiasModel.load(checked(config["parametric_model"], root))
    load_seconds = time.perf_counter() - started
    width, margin = 2 * plan["flank"] + 1, plan["margin"]
    sequence_scores(method, model, sequences[:plan["warmup_sites"]], width, margin)
    wall_start, cpu_start = time.perf_counter(), time.process_time()
    prediction = sequence_scores(method, model, sequences, width, margin)
    cpu_seconds, wall_seconds = time.process_time() - cpu_start, time.perf_counter() - wall_start
    prediction_id = array_identity(prediction)
    outdir.mkdir(parents=True, exist_ok=False)
    scores = outdir / "sequence_scores.npz"
    np.savez_compressed(scores, scores=prediction)
    output_record = record(scores)
    memory = peak_rss_mib()
    return write_json(outdir / "worker.json", {
        "schema": SCHEMA, "configuration_id": config["configuration_id"],
        "cell": cell, "method": method, "repetition": repetition,
        "pid": os.getpid(), "cpu": config["runtime"]["cpu"],
        "sites": len(sequences), "positions_per_strand": width,
        "inference_wall_seconds": wall_seconds, "inference_cpu_seconds": cpu_seconds,
        "model_load_seconds": load_seconds, "worker_peak_rss_mib": memory,
        "finite_values": int(np.isfinite(prediction).sum()), "prediction_id": prediction_id,
        "output": output_record, "promotion_evidence": False,
    }, "worker_id")


def freeze(plan, plan_path, root, outdir):
    inputs = {}

    def use(entry):
        path = checked(entry, root)
        inputs[str(path)] = {"path": str(path), "sha256": entry["sha256"]}
        return path

    source_path = use(plan["source_configuration"])
    source = read_json(source_path, "configuration_id")
    cells = {}
    genome = use(plan["genome"])
    use(plan["genome_index"])
    model = use(plan["parametric_model"])
    sidecar = use(plan["parametric_sidecar"])
    outdir.mkdir(parents=True, exist_ok=False)
    for cell in plan["cells"]:
        support = read_json(use(plan["supports"][cell]), "stage_id")
        if support["configuration_id"] != source["configuration_id"] or support["cell"] != cell:
            raise ValueError("support differs from the frozen stability configuration")
        # Explicit usecols excludes ChIP labels even though these are examined data.
        coverage = pd.read_csv(use(support["outputs"]["coverage"]), sep="\t",
                               usecols=KEYS + ["all_run_dwm_present"])
        selected = select_sites(coverage, plan["sites_per_cell"], plan["seed"])
        if set(selected.cell) != {cell}:
            raise ValueError("site cell differs from its frozen input")
        ref = source["references"][cell]["profile_manifest"]
        if (ref["metadata"]["bias_model_sha256"] != plan["parametric_model"]["sha256"]
                or ref["metadata"]["genome_sha256"] != plan["genome"]["sha256"]):
            raise ValueError("reference profiles used a different bias model or genome")
        sites = pd.read_csv(use({"path": ref["sites"], "sha256": ref["sites_sha256"]}), sep="\t", usecols=KEYS)
        if sites.duplicated(KEYS).any():
            raise ValueError("archived reference loci are not unique")
        indexes = pd.MultiIndex.from_frame(sites[KEYS]).get_indexer(pd.MultiIndex.from_frame(selected[KEYS]))
        if (indexes < 0).any():
            raise ValueError("selected loci are missing from the archived reference")
        with np.load(use({"path": ref["profiles_npz"], "sha256": ref["profiles_sha256"]}), allow_pickle=False) as values:
            archived = np.stack([values["plus_log_bias"][indexes], values["minus_log_bias"][indexes]], axis=1)
        reverse = selected.TFBS_strand.eq("-").to_numpy()
        archived[reverse] = archived[reverse, ::-1, ::-1]
        sequences = []
        with open_fasta(genome) as fasta:
            for row in selected.itertuples(index=False):
                center = (int(row.TFBS_start) + int(row.TFBS_end)) // 2
                sequence = fasta.fetch(row.TFBS_chr, center-plan["flank"]-plan["margin"],
                                       center+plan["flank"]+plan["margin"]+1).upper()
                if len(sequence) != 2 * (plan["flank"] + plan["margin"]) + 1:
                    raise ValueError("reference sequence is incomplete")
                sequences.append(sequence)
        bundle = outdir / f"{cell}.sequence_inputs.npz"
        np.savez_compressed(bundle, sequences=np.asarray(sequences), archived_parametric=archived)
        selection_path = outdir / f"{cell}.sites.tsv"
        selected.to_csv(selection_path, sep="\t", index=False)
        dwm_model = use(plan["dwm_models"][cell])
        cells[cell] = {"bundle": record(bundle), "selection": record(selection_path),
                       "dwm_model": record(dwm_model)}
    cpu_model = next((line.split(":", 1)[1].strip() for line in Path("/proc/cpuinfo").read_text().splitlines()
                      if line.startswith("model name")), "unavailable")
    runtime = {"python": sys.version, "platform": platform.system(), "machine": platform.machine(),
               "cpu_model": cpu_model, "cpu": sorted(os.sched_getaffinity(0))[plan["cpu_affinity_index"]],
               "versions": {name: version(name) for name in ("numpy", "scipy", "pandas", "pysam", "matplotlib")},
               "native_extension": record(native_sequences.__file__), "threads": THREADS}
    return write_json(outdir / "configuration.freeze.json", {
        "schema": SCHEMA, "plan": plan, "plan_file": record(plan_path), "cells": cells,
        "inputs": list(inputs.values()), "parametric_model": record(model),
        "parametric_sidecar": record(sidecar), "runtime": runtime,
        "code": {name: digest_file(REPOSITORY / name) for name in CODE},
        "promotion_evidence": False,
    }, "configuration_id")


def verify_configuration(config, plan_path, code_ref=None):
    if config["plan_file"]["sha256"] != digest_file(plan_path):
        raise ValueError("plan differs from frozen configuration")
    for entry in config["inputs"]:
        checked(entry, Path("/"))
    checked(config["runtime"]["native_extension"], Path("/"))
    for name, digest in config["code"].items():
        actual = (sha256(subprocess.check_output(["git", "show", f"{code_ref}:{name}"], cwd=REPOSITORY)).hexdigest()
                  if code_ref else digest_file(REPOSITORY / name))
        if actual != digest:
            raise ValueError(f"frozen code changed: {name}")
    for inputs in config["cells"].values():
        for entry in inputs.values():
            checked(entry, Path("/"))


def verify_worker(path, config, cell, method, repetition):
    done = read_json(path, "worker_id")
    if (done["configuration_id"] != config["configuration_id"] or done["cell"] != cell
            or done["method"] != method or done["promotion_evidence"] is not False
            or done["sites"] != config["plan"]["sites_per_cell"] or done["repetition"] != repetition):
        raise ValueError("worker differs from the frozen benchmark")
    with np.load(checked(done["output"], Path("/")), allow_pickle=False) as values:
        prediction = values["scores"]
    if prediction.shape != (done["sites"], 2, 201) or array_identity(prediction) != done["prediction_id"]:
        raise ValueError("worker prediction identity or dimensions differ")
    if method == METHODS[1]:
        with np.load(config["cells"][cell]["bundle"]["path"], allow_pickle=False) as values:
            np.testing.assert_allclose(prediction, values["archived_parametric"], rtol=1e-12, atol=1e-12, equal_nan=True)
    for key in ("inference_wall_seconds", "inference_cpu_seconds", "worker_peak_rss_mib"):
        if not np.isfinite(done[key]) or done[key] <= 0:
            raise ValueError("invalid component resource measurement")
    return done


def summarize(rows):
    frame = pd.DataFrame(rows)
    result = []
    for cell, group in frame.groupby("cell", sort=True):
        if (len(group) != 12 or group.duplicated(["method", "repetition"]).any()
                or set(group.method) != set(METHODS)
                or any(set(part.repetition) != set(range(6)) for _, part in group.groupby("method"))):
            raise ValueError("resource summary requires all six paired repetitions per cell")
        if (group.groupby("method").prediction_id.nunique() != 1).any():
            raise ValueError("predictions changed between identical timing repetitions")
        baseline = group[group.method.eq(METHODS[0])]
        candidate = group[group.method.eq(METHODS[1])]
        result.append({"cell": cell, "scope": "sequence_inference_component_only",
                       "component_mean_wall_ratio": candidate.inference_wall_seconds.mean() / baseline.inference_wall_seconds.mean(),
                       "component_mean_cpu_ratio": candidate.inference_cpu_seconds.mean() / baseline.inference_cpu_seconds.mean(),
                       "component_max_rss_ratio": candidate.worker_peak_rss_mib.max() / baseline.worker_peak_rss_mib.max(),
                       "promotion_evidence": False, "end_to_end_runtime_measured": False,
                       "end_to_end_memory_measured": False})
    return result


def run(config, config_path, outdir, verify_only):
    rows, stages = [], []
    for cell in config["plan"]["cells"]:
        for repetition in range(config["plan"]["repetitions"]):
            for method in METHODS[::1 if repetition % 2 == 0 else -1]:
                stage_dir = outdir / cell / f"repeat_{repetition}" / method
                stage = stage_dir / "worker.json"
                if not stage.exists():
                    if verify_only or stage_dir.exists():
                        raise ValueError(f"missing or partial benchmark worker: {stage}")
                    print(f"Measuring {cell} repetition {repetition}: {method}", flush=True)
                    subprocess.run([sys.executable, str(Path(__file__).resolve()), "--worker-config", str(config_path),
                                    "--cell", cell, "--method", method, "--repetition", str(repetition),
                                    "--outdir", str(stage_dir)],
                                   env={**os.environ, **THREADS}, check=True)
                done = verify_worker(stage, config, cell, method, repetition)
                rows.append({k: v for k, v in done.items() if k not in ("output", "schema")})
                stages.append(record(stage))
    summary = summarize(rows)
    final = outdir / "inference_manifest.json"
    if final.exists():
        done = read_json(final, "result_id")
        if done["configuration_id"] != config["configuration_id"] or done["workers"] != stages:
            raise ValueError("completed resource result differs from the verified matrix")
        for entry in done["outputs"].values():
            checked(entry, Path("/"))
        return done
    if verify_only:
        raise ValueError("benchmark result manifest is missing")
    metrics_path, summary_path = outdir / "component_measurements.tsv", outdir / "component_ratios.tsv"
    pd.DataFrame(rows).to_csv(metrics_path, sep="\t", index=False)
    pd.DataFrame(summary).to_csv(summary_path, sep="\t", index=False)
    return write_json(final, {"schema": SCHEMA, "configuration_id": config["configuration_id"],
                             "workers": stages, "outputs": {"measurements": record(metrics_path),
                                                              "ratios": record(summary_path)},
                             "promotion_evidence": False}, "result_id")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path)
    parser.add_argument("--data-root", type=Path)
    parser.add_argument("--outdir", type=Path, required=True)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--freeze-only", action="store_true")
    mode.add_argument("--verify-only", action="store_true")
    parser.add_argument("--code-ref")
    parser.add_argument("--worker-config", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--cell", help=argparse.SUPPRESS)
    parser.add_argument("--method", help=argparse.SUPPRESS)
    parser.add_argument("--repetition", type=int, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker_config:
        worker(args.worker_config, args.cell, args.method, args.repetition, args.outdir)
        return
    if not args.plan or not args.data_root or (args.code_ref and not args.verify_only):
        parser.error("--plan and --data-root are required; --code-ref is verification-only")
    plan = json.loads(args.plan.read_text())
    validate_plan(plan)
    config_path = args.outdir / "configuration.freeze.json"
    if config_path.exists():
        config = read_json(config_path, "configuration_id")
        verify_configuration(config, args.plan, args.code_ref)
    elif args.verify_only:
        raise ValueError("frozen benchmark configuration is missing")
    else:
        config = freeze(plan, args.plan, args.data_root, args.outdir)
    if args.freeze_only:
        print("Frozen sequence inputs, archived models, code and component measurement scope.")
        return
    if not args.verify_only:
        if sys.version != config["runtime"]["python"] or any(version(k) != v for k, v in config["runtime"]["versions"].items()):
            raise ValueError("runtime differs from the frozen benchmark")
    result = run(config, config_path, args.outdir, args.verify_only)
    print(("Verified without rerunning: " if args.verify_only else "Completed component diagnostic: ") + result["result_id"])


if __name__ == "__main__":
    main()
