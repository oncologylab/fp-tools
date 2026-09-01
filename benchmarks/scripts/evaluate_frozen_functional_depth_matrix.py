#!/usr/bin/env python3
"""Evaluate frozen parametric functional models across BAM depths and seeds."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict
from hashlib import sha256
import json
from pathlib import Path
import sys
from time import perf_counter

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))
sys.path.insert(0, str(SCRIPT_DIR))

from build_strand_functional_profiles import (  # noqa: E402
    extract_strand_cut_profiles,
    orient_strand_log_bias,
    site_hashes,
    write_profiles,
)
from evaluate_frozen_functional_policy import (  # noqa: E402
    aggregate_curve,
    candidate_score_and_profile,
    metric_record,
    validate_policy,
)
from evaluate_functional_footprints import (  # noqa: E402
    chromosome_split,
    load_or_extract_profiles,
)
from evaluate_parametric_factorization import (  # noqa: E402
    load_safe_configuration,
    residual_score,
)
from evaluate_strand_functional_templates import PROFILE_ARRAYS  # noqa: E402
from evaluate_strand_label_free_models import file_sha256  # noqa: E402
from fp_tools.tools.functional_footprints import (  # noqa: E402
    construct_strand_functional_profiles,
    normalize_functional_profiles,
    orient_profiles,
    profile_descriptors,
)
from fp_tools.tools.parametric_factorization import (  # noqa: E402
    FrozenParametricFactorization,
)


SCHEMA = "fp-tools-frozen-functional-depth-matrix-v1"


def expected_from_cached_log_bias(
    observed: np.ndarray,
    log_bias: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Apply the frozen sequence propensity while preserving each site total."""

    observed = np.asarray(observed, dtype=float)
    log_bias = np.asarray(log_bias, dtype=float)
    if observed.ndim != 2 or log_bias.shape != observed.shape:
        raise ValueError("observed counts and cached log bias must have equal 2D shape")
    finite = np.isfinite(log_bias)
    valid = finite.any(axis=1)
    row_maximum = np.max(np.where(finite, log_bias, -np.inf), axis=1)
    centered = np.where(
        finite,
        log_bias - np.where(valid, row_maximum, 0.0)[:, None],
        -np.inf,
    )
    propensity = np.exp(centered)
    denominator = propensity.sum(axis=1)
    total = observed.sum(axis=1)
    expected = np.divide(
        propensity * total[:, None],
        denominator[:, None],
        out=np.zeros_like(propensity),
        where=denominator[:, None] > 0,
    )
    return expected, valid


def scale_expected_to_observed(
    observed: np.ndarray,
    expected: np.ndarray,
) -> np.ndarray:
    """Remove track-normalization differences while preserving expected shape."""

    observed = np.asarray(observed, dtype=float)
    expected = np.maximum(np.asarray(expected, dtype=float), 0.0)
    if observed.ndim != 2 or expected.shape != observed.shape:
        raise ValueError("observed and expected profiles must have equal 2D shape")
    expected_total = expected.sum(axis=1)
    observed_total = observed.sum(axis=1)
    return np.divide(
        expected * observed_total[:, None],
        expected_total[:, None],
        out=np.zeros_like(expected),
        where=expected_total[:, None] > 0,
    )


def parse_sample(value: str) -> tuple[str, str]:
    fields = value.split(",", 1)
    if len(fields) != 2 or not all(fields):
        raise argparse.ArgumentTypeError("sample must use CELL,SAMPLE")
    return fields[0], fields[1]


def parse_reference(value: str) -> tuple[str, str, Path]:
    fields = value.split(",", 2)
    if len(fields) != 3 or not all(fields):
        raise argparse.ArgumentTypeError("reference must use MODEL,CELL,JSON")
    return fields[0], fields[1], Path(fields[2])


def discover_signals(
    root: Path,
    samples: dict[str, str],
    depths: list[str],
    seeds: list[int],
    *,
    allow_incomplete: bool,
) -> pd.DataFrame:
    rows = []
    for cell, sample in sorted(samples.items()):
        for depth in depths:
            for seed in seeds:
                directory = root / sample / depth / f"seed_{seed}"
                bams = sorted(directory.glob("*.bam"))
                expected = sorted((directory / "fp_tools_dwm").glob("*_expected.bw"))
                if len(bams) != 1 or len(expected) != 1:
                    if allow_incomplete:
                        continue
                    raise ValueError(
                        f"expected one BAM and DWM expected BigWig in {directory}"
                    )
                bam = bams[0]
                bai = Path(str(bam) + ".bai")
                if not bai.is_file():
                    alternative = bam.with_suffix(".bai")
                    bai = alternative if alternative.is_file() else bai
                if not all(path.is_file() and path.stat().st_size > 0 for path in (bam, bai, expected[0])):
                    if allow_incomplete:
                        continue
                    raise ValueError(f"incomplete depth input in {directory}")
                rows.append(
                    {
                        "cell": cell,
                        "sample": sample,
                        "depth": depth,
                        "seed": seed,
                        "bam": str(bam),
                        "bai": str(bai),
                        "dwm_expected": str(expected[0]),
                    }
                )
    return pd.DataFrame(rows)


def checksum_signals(signals: pd.DataFrame) -> pd.DataFrame:
    output = signals.copy()
    for column in ("bam", "bai", "dwm_expected"):
        output[f"{column}_sha256"] = [file_sha256(path) for path in output[column]]
    return output


def load_reference(
    path: Path,
    *,
    cell: str,
    study: dict,
    policy_tfs: set[str],
) -> tuple[pd.DataFrame, dict[str, np.ndarray], dict]:
    document = json.loads(path.read_text(encoding="utf-8"))
    if document.get("schema") != "fp-tools-strand-functional-profiles-v1":
        raise ValueError(f"unsupported reference profile artifact: {path}")
    if document.get("metadata", {}).get("labels_used") is not False:
        raise ValueError(f"reference profiles do not certify label-free construction: {path}")
    metadata_cell = document.get("metadata", {}).get("cell")
    if metadata_cell is not None and str(metadata_cell) != cell:
        raise ValueError(f"reference cell does not match {cell}: {path}")
    profiles_path = Path(document["profiles_npz"])
    sites_path = Path(document["sites"])
    if file_sha256(profiles_path) != document["profiles_sha256"]:
        raise ValueError(f"reference profile checksum mismatch: {profiles_path}")
    if file_sha256(sites_path) != document["sites_sha256"]:
        raise ValueError(f"reference site checksum mismatch: {sites_path}")
    sites = pd.read_csv(sites_path, sep="\t").reset_index(drop=True)
    if "cell" not in sites or set(sites["cell"].astype(str)) != {cell}:
        raise ValueError(f"reference sites do not exclusively contain {cell}: {path}")
    split = sites["TFBS_chr"].map(
        lambda chromosome: chromosome_split(str(chromosome), study)
    )
    selected = split.eq("validation") & sites["tf"].astype(str).isin(policy_tfs)
    indexes = np.flatnonzero(selected.to_numpy())
    sites = sites.iloc[indexes].reset_index(drop=True)
    sites["chromosome_split"] = "validation"
    with np.load(profiles_path, allow_pickle=False) as source:
        required = set(PROFILE_ARRAYS + ("valid", "site_hash", "plus_log_bias", "minus_log_bias"))
        missing = sorted(required.difference(source.files))
        if missing:
            raise ValueError(
                f"reference artifact is missing arrays {', '.join(missing)}: {path}"
            )
        if not np.array_equal(
            np.asarray(source["site_hash"], dtype=np.uint64),
            site_hashes(pd.read_csv(sites_path, sep="\t")),
        ):
            raise ValueError(f"reference site order mismatch: {path}")
        arrays = {
            name: np.asarray(source[name])[indexes]
            for name in required
        }
    return sites, arrays, document


def artifact_prefix(outdir: Path, row) -> Path:
    return (
        outdir
        / "parametric_profiles"
        / str(row.cell)
        / str(row.depth)
        / f"seed_{int(row.seed)}"
        / "profiles"
    )


def build_depth_artifact(
    row,
    *,
    outdir: Path,
    sites: pd.DataFrame,
    reference_arrays: dict[str, np.ndarray],
    reference_document: dict,
    flank: int,
    minimum_mapq: int,
) -> dict:
    started = perf_counter()
    prefix = artifact_prefix(outdir, row)
    json_path = Path(str(prefix) + ".json")
    cache_path = (
        outdir
        / "dwm_profile_cache"
        / str(row.cell)
        / str(row.depth)
        / f"seed_{int(row.seed)}.npz"
    )
    if json_path.is_file():
        document = json.loads(json_path.read_text(encoding="utf-8"))
        if (
            document.get("metadata", {}).get("bam_sha256") != row.bam_sha256
            or document.get("metadata", {}).get("reference_profiles_sha256")
            != reference_document["profiles_sha256"]
        ):
            raise ValueError(f"resumable depth artifact changed inputs: {json_path}")
        if file_sha256(document["profiles_npz"]) != document["profiles_sha256"]:
            raise ValueError(f"resumable depth profile checksum mismatch: {json_path}")
        if file_sha256(document["sites"]) != document["sites_sha256"]:
            raise ValueError(f"resumable depth site checksum mismatch: {json_path}")
    else:
        strands = sites["TFBS_strand"].astype(str).tolist()
        genomic_plus_bias, genomic_minus_bias, _combined = orient_strand_log_bias(
            reference_arrays["plus_log_bias"],
            reference_arrays["minus_log_bias"],
            strands,
        )
        plus, minus, bam_valid = extract_strand_cut_profiles(
            sites,
            row.bam,
            flank=flank,
            read_shift=(4, -5),
            minimum_mapq=minimum_mapq,
            keep_duplicates=False,
        )
        plus_expected, plus_bias_valid = expected_from_cached_log_bias(
            plus, genomic_plus_bias
        )
        minus_expected, minus_bias_valid = expected_from_cached_log_bias(
            minus, genomic_minus_bias
        )
        bias_valid = plus_bias_valid & minus_bias_valid
        profiles = construct_strand_functional_profiles(
            plus,
            minus,
            plus_expected,
            minus_expected,
            strands,
            dispersion=0.0,
        )
        with np.errstate(invalid="ignore"):
            combined_log_bias = (
                np.logaddexp(
                    reference_arrays["plus_log_bias"],
                    reference_arrays["minus_log_bias"],
                )
                - np.log(2.0)
            )
        write_profiles(
            prefix,
            sites,
            profiles,
            bam_valid & bias_valid,
            {
                "cell": row.cell,
                "sample": row.sample,
                "depth": row.depth,
                "seed": int(row.seed),
                "bam": row.bam,
                "bam_sha256": row.bam_sha256,
                "reference_profiles": reference_document["profiles_npz"],
                "reference_profiles_sha256": reference_document["profiles_sha256"],
                "bias_model": reference_document["metadata"]["bias_model"],
                "bias_model_sha256": reference_document["metadata"][
                    "bias_model_sha256"
                ],
                "read_shift": [4, -5],
                "flank": flank,
                "minimum_mapq": minimum_mapq,
                "keep_duplicates": False,
                "labels_used": False,
                "chromosome_split": "validation",
            },
            log_bias=(
                reference_arrays["plus_log_bias"],
                reference_arrays["minus_log_bias"],
                combined_log_bias,
            ),
        )
    load_or_extract_profiles(
        sites,
        row.dwm_expected,
        cache_path,
        flank,
    )
    return {
        "cell": row.cell,
        "sample": row.sample,
        "depth": row.depth,
        "seed": int(row.seed),
        "artifact": str(json_path),
        "dwm_cache": str(cache_path),
        "build_seconds": perf_counter() - started,
    }


def load_built_artifact(
    record: dict,
) -> tuple[pd.DataFrame, dict[str, np.ndarray], np.ndarray, np.ndarray]:
    document = json.loads(Path(record["artifact"]).read_text(encoding="utf-8"))
    with np.load(document["profiles_npz"], allow_pickle=False) as source:
        arrays = {
            name: np.asarray(source[name])
            for name in PROFILE_ARRAYS + ("valid", "site_hash")
        }
    sites = pd.read_csv(document["sites"], sep="\t").reset_index(drop=True)
    if not np.array_equal(site_hashes(sites), arrays["site_hash"]):
        raise ValueError(f"built profile site order mismatch: {record['artifact']}")
    with np.load(record["dwm_cache"], allow_pickle=False) as source:
        if not np.array_equal(
            np.asarray(source["site_hash"], dtype=np.uint64), arrays["site_hash"]
        ):
            raise ValueError(f"DWM depth cache site order mismatch: {record['dwm_cache']}")
        dwm = np.asarray(source["profiles"], dtype=float)
        dwm_valid = np.asarray(source["valid"], dtype=bool)
    dwm = orient_profiles(dwm, sites["TFBS_strand"].astype(str))
    return sites, arrays, dwm, dwm_valid


def summarize_metrics(metrics: pd.DataFrame) -> pd.DataFrame:
    keys = [
        "cell",
        "tf",
        "motif_family",
        "candidate_id",
        "method",
        "depth",
    ]
    summary = (
        metrics.groupby(keys, sort=True)
        .agg(
            seeds=("seed", "nunique"),
            auroc_mean=("auroc", "mean"),
            auroc_sd=("auroc", "std"),
            auprc_mean=("auprc", "mean"),
            auprc_sd=("auprc", "std"),
            auroc_gain_over_dwm_mean=("auroc_gain_over_dwm", "mean"),
            auroc_gain_over_dwm_sd=("auroc_gain_over_dwm", "std"),
            relative_auprc_gain_over_dwm_mean=(
                "relative_auprc_gain_over_dwm",
                "mean",
            ),
            relative_auprc_gain_over_dwm_sd=(
                "relative_auprc_gain_over_dwm",
                "std",
            ),
            auroc_gain_over_raw_mean=("auroc_gain_over_raw", "mean"),
            auroc_gain_over_raw_sd=("auroc_gain_over_raw", "std"),
            relative_auprc_gain_over_raw_mean=(
                "relative_auprc_gain_over_raw",
                "mean",
            ),
            relative_auprc_gain_over_raw_sd=(
                "relative_auprc_gain_over_raw",
                "std",
            ),
            functional_separation_mean=("functional_separation", "mean"),
            functional_separation_sd=("functional_separation", "std"),
            functional_separation_relative_change_over_dwm_mean=(
                "functional_separation_relative_change_over_dwm",
                "mean",
            ),
            functional_separation_relative_change_over_raw_mean=(
                "functional_separation_relative_change_over_raw",
                "mean",
            ),
            brier_mean=("brier", "mean"),
            calibration_error_mean=("calibration_error", "mean"),
            prediction_seconds_mean=("prediction_seconds", "mean"),
        )
        .reset_index()
    )
    direction = (
        metrics.groupby(keys, sort=True)["auroc_gain_over_dwm"]
        .apply(lambda values: float(np.mean(values > 0)))
        .rename("auroc_gain_positive_fraction")
        .reset_index()
    )
    return summary.merge(direction, on=keys, validate="one_to_one")


def add_depth_baseline_deltas(task_rows: list[dict]) -> list[dict]:
    """Attach DWM and raw guardrail deltas to one TF/depth/seed result set."""

    by_method = {str(row["method"]): row for row in task_rows}
    try:
        dwm = by_method["DWM_conventional_geometry"]
        raw = by_method["raw_geometry"]
    except KeyError as exc:
        raise ValueError("depth metrics require DWM and raw baselines") from exc
    for metric in task_rows:
        metric["auroc_gain_over_dwm"] = metric["auroc"] - dwm["auroc"]
        metric["relative_auprc_gain_over_dwm"] = (
            metric["auprc"] - dwm["auprc"]
        ) / max(dwm["auprc"], 1e-8)
        metric["functional_separation_relative_change_over_dwm"] = (
            metric["functional_separation"]
            / max(dwm["functional_separation"], 1e-8)
            - 1.0
        )
        metric["auroc_gain_over_raw"] = metric["auroc"] - raw["auroc"]
        metric["relative_auprc_gain_over_raw"] = (
            metric["auprc"] - raw["auprc"]
        ) / max(raw["auprc"], 1e-8)
        metric["functional_separation_relative_change_over_raw"] = (
            metric["functional_separation"]
            / max(raw["functional_separation"], 1e-8)
            - 1.0
        )
    return task_rows


def classify_depth(summary: pd.DataFrame) -> pd.DataFrame:
    candidates = summary[summary["method"].str.startswith("frozen_")]
    rows = []
    for (cell, tf, family), group in candidates.groupby(
        ["cell", "tf", "motif_family"], sort=True
    ):
        values = group.set_index("depth")
        low = values.loc["10m"] if "10m" in values.index else None
        high_name = next(
            (name for name in ("full", "50m", "25m") if name in values.index),
            None,
        )
        if low is None or high_name is None:
            continue
        high = values.loc[high_name]
        depth_gain = float(high.auroc_mean - low.auroc_mean)
        if float(high.auroc_mean) >= 0.65:
            status = "detectable_at_high_depth"
        elif depth_gain >= 0.03:
            status = "power_limited"
        elif float(high.auroc_gain_over_dwm_mean) >= 0.03:
            status = "model_improved_but_assay_weak"
        else:
            status = "shape_or_assay_limited"
        rows.append(
            {
                "cell": cell,
                "tf": tf,
                "motif_family": family,
                "low_depth": "10m",
                "high_depth": high_name,
                "low_auroc": float(low.auroc_mean),
                "high_auroc": float(high.auroc_mean),
                "depth_auroc_gain": depth_gain,
                "high_auroc_gain_over_dwm": float(high.auroc_gain_over_dwm_mean),
                "high_relative_auprc_gain_over_dwm": float(
                    high.relative_auprc_gain_over_dwm_mean
                ),
                "high_seed_positive_fraction": float(
                    high.auroc_gain_positive_fraction
                ),
                "classification": status,
            }
        )
    return pd.DataFrame(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--policy", type=Path, required=True)
    parser.add_argument("--reference", action="append", type=parse_reference, required=True)
    parser.add_argument("--reference-configuration", type=Path, required=True)
    parser.add_argument("--signals-root", type=Path, required=True)
    parser.add_argument("--sample", action="append", type=parse_sample, required=True)
    parser.add_argument("--depth", action="append", default=[])
    parser.add_argument("--seed", action="append", type=int, default=[])
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--minimum-mapq", type=int, default=30)
    parser.add_argument("--minimum-sites-per-class", type=int, default=100)
    parser.add_argument("--allow-incomplete", action="store_true")
    args = parser.parse_args(argv)

    policy, models = validate_policy(args.policy)
    study = json.loads(Path(policy["study"]["path"]).read_text(encoding="utf-8"))
    configuration = load_safe_configuration(args.reference_configuration)
    factorization = FrozenParametricFactorization.load(
        configuration["factorization_model"]["path"]
    )
    dispersion = float(factorization.total_dispersion_)
    samples = dict(args.sample)
    references = {(model, cell): path for model, cell, path in args.reference}
    policy_keys = {
        (str(record["bias_configuration"]), str(record["cell"]))
        for record, _candidate, _model in models
    }
    if set(references) != policy_keys or set(samples) != {cell for _model, cell in policy_keys}:
        raise ValueError("samples/references do not exactly match policy cells")
    depths = args.depth or ["10m", "25m", "50m"]
    seeds = args.seed or list(range(2026, 2031))
    signals = discover_signals(
        args.signals_root,
        samples,
        depths,
        seeds,
        allow_incomplete=args.allow_incomplete,
    )
    if signals.empty:
        raise ValueError("no complete depth signals were found")
    signals = checksum_signals(signals)
    args.outdir.mkdir(parents=True, exist_ok=True)
    signal_path = args.outdir / "depth_inputs.tsv"
    signals.to_csv(signal_path, sep="\t", index=False)

    reference_data = {}
    for key, path in sorted(references.items()):
        tfs = {
            str(record["tf"])
            for record, _candidate, _model in models
            if str(record["bias_configuration"]) == key[0]
            and str(record["cell"]) == key[1]
        }
        reference_data[key] = load_reference(
            path,
            cell=key[1],
            study=study,
            policy_tfs=tfs,
        )
    freeze = {
        "schema": "fp-tools-frozen-functional-depth-inputs-v1",
        "policy": {"path": str(args.policy), "sha256": file_sha256(args.policy)},
        "policy_id": policy["policy_id"],
        "reference_configuration": {
            "path": str(args.reference_configuration),
            "sha256": file_sha256(args.reference_configuration),
        },
        "signals": {"path": str(signal_path), "sha256": file_sha256(signal_path)},
        "references": [
            {"model": key[0], "cell": key[1], "path": str(path), "sha256": file_sha256(path)}
            for key, path in sorted(references.items())
        ],
        "validation_labels_used_for_original_selection": True,
        "models_refitted_by_depth": False,
    }
    canonical = json.dumps(freeze, sort_keys=True, separators=(",", ":"))
    freeze["depth_input_id"] = sha256(canonical.encode()).hexdigest()
    freeze_path = args.outdir / "depth_inputs.freeze.json"
    rendered = json.dumps(freeze, indent=2, sort_keys=True) + "\n"
    if freeze_path.exists() and freeze_path.read_text(encoding="utf-8") != rendered:
        raise ValueError(f"immutable depth input freeze differs: {freeze_path}")
    freeze_path.write_text(rendered, encoding="utf-8")

    built = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {}
        for row in signals.itertuples(index=False):
            key = next(key for key in references if key[1] == str(row.cell))
            sites, arrays, document = reference_data[key]
            future = executor.submit(
                build_depth_artifact,
                row,
                outdir=args.outdir,
                sites=sites,
                reference_arrays=arrays,
                reference_document=document,
                flank=int(study["profile_flank_bp"]),
                minimum_mapq=args.minimum_mapq,
            )
            futures[future] = (row.cell, row.depth, row.seed)
        for future in as_completed(futures):
            record = future.result()
            built.append(record)
            print(
                f"built {record['cell']} depth={record['depth']} seed={record['seed']}",
                flush=True,
            )

    metrics_rows = []
    curve_rows = []
    with threadpool_limits(limits=1):
        for built_record in sorted(
            built, key=lambda row: (row["depth"], row["seed"], row["cell"])
        ):
            cell = str(built_record["cell"])
            sites, arrays, dwm_expected, dwm_valid = load_built_artifact(built_record)
            positions = np.arange(arrays["plus_observed"].shape[1], dtype=float)
            positions -= arrays["plus_observed"].shape[1] // 2
            valid = (
                arrays["valid"].astype(bool)
                & dwm_valid
                & np.isfinite(dwm_expected).all(axis=1)
            )
            for record, candidate, model in models:
                if str(record["cell"]) != cell:
                    continue
                selected = sites["tf"].astype(str).eq(str(record["tf"])).to_numpy() & valid
                indexes = np.flatnonzero(selected)
                labels = sites.iloc[indexes]["chip_label"].to_numpy(dtype=int)
                if len(indexes) == 0 or np.unique(labels).size != 2:
                    continue
                observed = arrays["plus_observed"][indexes] + arrays["minus_observed"][indexes]
                started = perf_counter()
                candidate_score, candidate_profiles, fitted_profile = candidate_score_and_profile(
                    candidate,
                    model,
                    arrays,
                    indexes,
                    positions,
                    motif_score=sites.iloc[indexes]["motif_score"].to_numpy(dtype=float),
                )
                candidate_seconds = perf_counter() - started
                started = perf_counter()
                dwm_score, dwm_profiles = residual_score(
                    observed,
                    scale_expected_to_observed(
                        observed,
                        dwm_expected[indexes],
                    ),
                    positions,
                    "deviance",
                    dispersion,
                )
                dwm_seconds = perf_counter() - started
                raw_score = geometry_score(observed, positions)
                direct_expected = arrays["plus_expected"][indexes] + arrays["minus_expected"][indexes]
                started = perf_counter()
                direct_score, direct_profiles = residual_score(
                    observed,
                    direct_expected,
                    positions,
                    "deviance",
                    dispersion,
                )
                direct_seconds = perf_counter() - started
                methods = (
                    (
                        "DWM_conventional_geometry",
                        dwm_score,
                        normalize_functional_profiles(dwm_profiles, positions),
                        dwm_seconds,
                        None,
                    ),
                    (
                        "raw_geometry",
                        raw_score,
                        observed,
                        0.0,
                        None,
                    ),
                    (
                        "LOG21_direct_geometry",
                        direct_score,
                        normalize_functional_profiles(direct_profiles, positions),
                        direct_seconds,
                        None,
                    ),
                    (
                        f"frozen_{candidate.candidate_id}",
                        candidate_score,
                        candidate_profiles,
                        candidate_seconds,
                        fitted_profile,
                    ),
                )
                task_rows = []
                for method, score, profiles, seconds, model_profile in methods:
                    metric = metric_record(
                        record=record,
                        candidate=candidate,
                        method=method,
                        labels=labels,
                        score=score,
                        profiles=profiles,
                        positions=positions,
                        prediction_seconds=seconds,
                        minimum_sites_per_class=args.minimum_sites_per_class,
                        split_name="validation",
                    )
                    metric.update(
                        {
                            "sample": built_record["sample"],
                            "depth": built_record["depth"],
                            "seed": built_record["seed"],
                        }
                    )
                    task_rows.append(metric)
                    curve = aggregate_curve(
                        profiles,
                        labels,
                        sites.iloc[indexes]["TFBS_chr"].astype(str).to_numpy(),
                        iterations=0,
                        seed=args.seed[0] if args.seed else 2026,
                    )
                    descriptors = asdict(
                        profile_descriptors(curve["difference"], positions)
                    )
                    for offset, position in enumerate(positions.astype(int)):
                        curve_rows.append(
                            {
                                "cell": cell,
                                "tf": record["tf"],
                                "motif_family": record["motif_family"],
                                "candidate_id": candidate.candidate_id,
                                "method": method,
                                "sample": built_record["sample"],
                                "depth": built_record["depth"],
                                "seed": built_record["seed"],
                                "position": position,
                                "positive_mean": curve["positive_mean"][offset],
                                "negative_mean": curve["negative_mean"][offset],
                                "positive_minus_negative": curve["difference"][offset],
                                "frozen_model_profile": (
                                    np.nan
                                    if model_profile is None
                                    else model_profile[offset]
                                ),
                                **descriptors,
                            }
                        )
                metrics_rows.extend(add_depth_baseline_deltas(task_rows))

    metrics = pd.DataFrame(metrics_rows)
    profiles = pd.DataFrame(curve_rows)
    summary = summarize_metrics(metrics)
    classification = classify_depth(summary)
    metrics_path = args.outdir / "frozen_functional_depth_metrics.tsv.gz"
    profiles_path = args.outdir / "frozen_functional_depth_profiles.tsv.gz"
    summary_path = args.outdir / "frozen_functional_depth_summary.tsv"
    classification_path = args.outdir / "frozen_functional_depth_classification.tsv"
    artifacts_path = args.outdir / "frozen_functional_depth_artifacts.tsv"
    metrics.to_csv(metrics_path, sep="\t", index=False)
    profiles.to_csv(profiles_path, sep="\t", index=False)
    summary.to_csv(summary_path, sep="\t", index=False)
    classification.to_csv(classification_path, sep="\t", index=False)
    pd.DataFrame(built).to_csv(artifacts_path, sep="\t", index=False)
    manifest = {
        "schema": SCHEMA,
        "policy_id": policy["policy_id"],
        "depth_input_id": freeze["depth_input_id"],
        "models_refitted_by_depth": False,
        "raw_signal_guardrail": True,
        "validation_labels_used_for_original_selection": True,
        "depths": depths,
        "seeds": seeds,
        "dispersion": dispersion,
        "outputs": {
            "metrics": {"path": str(metrics_path), "sha256": file_sha256(metrics_path)},
            "profiles": {"path": str(profiles_path), "sha256": file_sha256(profiles_path)},
            "summary": {"path": str(summary_path), "sha256": file_sha256(summary_path)},
            "classification": {
                "path": str(classification_path),
                "sha256": file_sha256(classification_path),
            },
            "artifacts": {
                "path": str(artifacts_path),
                "sha256": file_sha256(artifacts_path),
            },
        },
    }
    (args.outdir / "frozen_functional_depth_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(classification.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
