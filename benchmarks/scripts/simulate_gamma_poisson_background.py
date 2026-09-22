#!/usr/bin/env python3
"""Frozen scalar predictive coverage simulation, including model misspecification."""

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gamma_poisson_global_expectation as model
import build_global_parametric_training_profiles as shared
from audit_direct_bias_safety import record, write_immutable, wilson_interval
from audit_chrombpnet_comparison import canonical_id, checked_path

ROOT = shared.ROOT
PROTOCOL = ROOT / "benchmarks/manifests/gamma_poisson_background_v1.json"


def case(background_mean, ratio, generator, *, replicates, rng, level=.95):
    if (not isinstance(replicates, int) or replicates < 1 or not np.isfinite(background_mean)
            or background_mean < 0 or not np.isfinite(ratio) or ratio <= 0):
        raise ValueError("invalid synthetic scenario")
    counts = rng.poisson(background_mean, replicates)
    target_mean = background_mean * ratio
    if generator == "poisson_working_model":
        observed = rng.poisson(target_mean, replicates)
    elif generator == "overdispersed_center_stress_phi1":
        observed = rng.negative_binomial(1, 1/(1+target_mean), replicates)
    else:
        raise ValueError("unknown generator")
    unique, inverse = np.unique(counts, return_inverse=True)
    low, high = model.predictive_interval(unique + model.PRIOR_SHAPE, ratio, level)
    lower, upper = low[inverse], high[inverse]
    covered = (observed >= lower) & (observed <= upper)
    interval = wilson_interval(int(covered.sum()), replicates)
    return dict(background_mean=background_mean, exposure_ratio=ratio, generator=generator,
        simulations=replicates, covered=int(covered.sum()), coverage=float(covered.mean()),
        coverage_lower_95=float(interval[0]), coverage_upper_95=float(interval[1]),
        below_interval=int((observed < lower).sum()), above_interval=int((observed > upper).sum()),
        mean_interval_width=float((upper-lower).mean()), prior_only_fraction=float((counts == 0).mean()),
        true_target_mean=target_mean, mean_predicted_count=float(((counts+.5)*ratio).mean()),
        nominal_coverage=level, functional_uncertainty_gate_eligible=False)


def freeze(out):
    if out.exists():
        raise ValueError("use a fresh output directory")
    protocol = json.loads(PROTOCOL.read_text())
    config = dict(protocol=protocol, inputs=[record(PROTOCOL)],
        code=[record(Path(__file__)), record(Path(model.__file__)),
              record(ROOT / "benchmarks/scripts/global_parametric_scoring.py")],
        runtime=shared.parent.extension.runtime(), real_data_read=False,
        prior_retuned=False, promotion_eligible=False)
    config["configuration_id"] = canonical_id(config)
    out.mkdir(parents=True)
    write_immutable(out / "configuration.freeze.json", config)


def execute(out, verify_only=False):
    config = shared.parent.extension.read_document(out / "configuration.freeze.json", "configuration_id")
    if (config["runtime"] != shared.parent.extension.runtime()
            or any(config[k] is not False for k in ("real_data_read", "prior_retuned", "promotion_eligible"))):
        raise ValueError("simulation scope or runtime changed")
    for item in config["inputs"] + config["code"]:
        checked_path(item, ROOT)
    if config["protocol"] != json.loads(PROTOCOL.read_text()):
        raise ValueError("protocol changed")
    result_path = out / "simulation_manifest.json"
    if result_path.exists():
        shared.parent.base.verify_stage(result_path, config["configuration_id"], ROOT)
        print("Verified scalar predictive simulation", flush=True)
        return
    if verify_only or (out / "coverage.tsv").exists():
        raise ValueError("missing result or partial output preserved")
    settings = config["protocol"]["synthetic"]
    rng = np.random.default_rng(settings["seed"])
    rows = [case(mean, ratio, generator, replicates=settings["replicates_per_case"], rng=rng,
                 level=settings["nominal_coverage"])
            for generator in settings["generators"] for mean in settings["background_means"]
            for ratio in settings["center_to_background_exposure_ratios"]]
    target = out / "coverage.tsv"
    pd.DataFrame(rows).to_csv(target, sep="\t", index=False)
    result = dict(configuration_id=config["configuration_id"], outputs=dict(coverage=record(target)),
        real_data_read=False, prior_retuned=False, promotion_eligible=False,
        evidence="scalar_synthetic_predictive_coverage_with_explicit_misspecification_stress")
    result["stage_id"] = canonical_id(result)
    write_immutable(result_path, result)
    print(f"Completed all {len(rows)} cases; no functional or biological calibration claim", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--freeze-only", action="store_true")
    modes.add_argument("--verify-only", action="store_true")
    args = parser.parse_args(); out = args.outdir.resolve()
    if not out.is_relative_to(ROOT / "benchmarks/results"):
        parser.error("use research results")
    if args.freeze_only:
        freeze(out)
    else:
        with shared.parent.extension.exclusive_run(out):
            execute(out, args.verify_only)
