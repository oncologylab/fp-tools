#!/usr/bin/env python3
"""Qualify an existing direct-bias safety screen using archived DWM coverage.

Reuse all scores and frozen thresholds. Sparse DWM tracks may be absent outside
the regions originally processed; zero-filled arrays do not establish coverage.
This stage never recalibrates thresholds or reruns a completed scoring stage.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from audit_direct_bias_safety import (
    METHODS, OPERATOR, SCHEMA, SITE_KEYS, canonical_id, check_identity, checked_path,
    load_artifact, read_thresholds, record, safety_row, site_hashes, write_immutable,
)

SCHEMA_COVERAGE = "fp-tools-direct-bias-covered-candidate-safety-v1"


def profile_indexes(sites, valid, scores):
    retained = np.flatnonzero(valid)
    if scores.duplicated(SITE_KEYS).any() or len(scores) != len(retained):
        raise ValueError("stored score support differs from valid control profiles")
    order = pd.Index(site_hashes(sites.iloc[retained])).get_indexer(site_hashes(scores))
    if np.any(order < 0):
        raise ValueError("stored scores and covered profiles use different motif loci")
    return retained[order]


def comparable_profiles(observed, dwm_observed, dwm_expected):
    """Require positive DWM expectation and matching observed count shape.

    DWM uncorrected tracks may use a constant library normalization. Comparing
    proportions avoids treating that harmless scale change as different reads.
    """
    arrays = [np.asarray(a, dtype=float) for a in (observed, dwm_observed, dwm_expected)]
    if any(a.shape != arrays[0].shape or a.ndim != 2 or not np.isfinite(a).all()
           or np.any(a < 0) for a in arrays):
        raise ValueError("coverage audit requires aligned finite nonnegative profiles")
    total = [a.sum(axis=1) for a in arrays]
    informative = total[0] > 0
    fractions = [np.divide(a, t[:, None], out=np.zeros_like(a), where=t[:, None] > 0)
                 for a, t in zip(arrays[:2], total[:2])]
    equal_shape = np.all(np.isclose(fractions[0], fractions[1], rtol=1e-5, atol=1e-7), axis=1)
    comparable = informative & (total[1] > 0) & (total[2] > 0) & equal_shape
    disagreement = informative & (total[1] > 0) & ~equal_shape
    return comparable, disagreement


def qualify_row(candidate, baseline, totals, comparable, disagreement, *, candidate_threshold,
                baseline_threshold, threshold_eligible, motif_residual_flag):
    common = np.asarray(comparable, dtype=bool)
    if common.shape != np.asarray(candidate).shape or np.asarray(disagreement).shape != common.shape:
        raise ValueError("coverage flags must match score support")
    arguments = dict(candidate_threshold=candidate_threshold, baseline_threshold=baseline_threshold,
                     threshold_eligible=threshold_eligible)
    full = safety_row(candidate, baseline, totals, **arguments)
    paired = safety_row(np.asarray(candidate)[common], np.asarray(baseline)[common],
                        np.asarray(totals)[common], **arguments)
    # A filled zero baseline must not make a relative safety test easier to pass.
    absolute_pass = bool(threshold_eligible and full["support"] > 0 and full["informative_support"] > 0
                         and all(full[prefix + "false_positive_rate"] <= .05
                                 and full[prefix + "false_positive_rate_upper_95"] <= .05
                                 for prefix in ("", "informative_")))
    result = {key: value for key, value in full.items()
              if key not in {"dwm_calls", "passes_candidate_screen"}
              and not key.endswith("increase_over_dwm")}
    result.update({"paired_" + key: value for key, value in paired.items()
                   if key in {"support", "false_positive_calls", "dwm_calls", "false_positive_rate",
                              "false_positive_rate_upper_95", "false_positive_rate_increase_over_dwm"}})
    result.update(dwm_unavailable_informative_sites=int(full["informative_support"] - common.sum()),
                  observed_shape_disagreements=int(np.asarray(disagreement).sum()),
                  motif_residual_flag=bool(motif_residual_flag), absolute_fpr_screen_pass=absolute_pass,
                  passes_candidate_screen=bool(absolute_pass and paired["passes_candidate_screen"]
                                               and not np.asarray(disagreement).any() and not motif_residual_flag))
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--safety-manifest", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)
    source = json.loads(args.safety_manifest.read_text())
    check_identity(source, "safety_id")
    if (source.get("schema") != SCHEMA or source.get("operator") != OPERATOR
            or source.get("replicate") != "NakedDNA_rep2" or source.get("new_holdout_labels_read") is not False
            or source.get("final_untouched_safety") is not False):
        raise ValueError("coverage audit accepts only the direct replicate-2 candidate screen")
    config = json.loads(checked_path(source["configuration"], args.data_root).read_text())
    check_identity(config, "configuration_id")
    if (config["configuration_id"] != source["configuration_id"]
            or not config["cells"] or not set(config["cells"]).issubset({"K562", "HepG2"})):
        raise ValueError("safety configuration does not match development scope")
    for item in config["code"]:
        checked_path(item, args.data_root)
    frozen, thresholds = read_thresholds(checked_path(source["thresholds"], args.data_root),
                                         config["configuration_id"], args.data_root)
    if frozen["threshold_id"] != source["threshold_id"]:
        raise ValueError("safety threshold identity mismatch")
    tables = {name: pd.read_csv(checked_path(item, args.data_root), sep="\t")
              for name, item in source["outputs"].items()}
    freeze = {"source_safety": record(args.safety_manifest), "code": record(Path(__file__)),
              "scoring_code": record(Path(__file__).with_name("audit_direct_bias_safety.py")),
              "same_thresholds_and_scores": True}
    freeze["coverage_input_id"] = canonical_id(freeze)
    args.outdir.mkdir(parents=True, exist_ok=True)
    write_immutable(args.outdir / "coverage_inputs.freeze.json", freeze)
    manifest_path = args.outdir / "coverage_manifest.json"
    if manifest_path.exists():
        done = json.loads(manifest_path.read_text())
        check_identity(done, "coverage_id")
        if done["coverage_input_id"] != freeze["coverage_input_id"]:
            raise ValueError("completed coverage audit input mismatch")
        for item in done["outputs"].values():
            checked_path(item, args.data_root)
        print("Verified completed coverage qualification; no stage rerun.")
        return 0
    lookup = thresholds.set_index(["cell", "tf", "method"])
    rows, site_frames = [], []
    for cell, inputs in config["cells"].items():
        sites, observed, _, valid, _ = load_artifact(inputs["naked"], args.data_root, cell, control=True)
        dwm_sites, dwm_observed, dwm_expected, dwm_valid, _ = load_artifact(inputs["naked_dwm"], args.data_root, cell, control=True)
        order = pd.Index(site_hashes(dwm_sites)).get_indexer(site_hashes(sites))
        if np.any(order < 0):
            raise ValueError("DWM coverage lacks the candidate's motif loci")
        dwm_observed, dwm_expected = dwm_observed[order], dwm_expected[order]
        valid &= dwm_valid[order]
        scores = tables["site_scores"].query("cell == @cell").reset_index(drop=True)
        indexes = profile_indexes(sites, valid, scores)
        np.testing.assert_allclose(scores.total_cuts, observed[indexes].sum(axis=1), rtol=0, atol=0)
        common, disagreement = comparable_profiles(observed[indexes], dwm_observed[indexes], dwm_expected[indexes])
        scores["dwm_comparable"] = common
        scores["observed_shape_disagreement"] = disagreement
        site_frames.append(scores)
        for tf, task in scores.groupby("tf", sort=True):
            ix = task.index.to_numpy(dtype=int)
            cutoffs = [lookup.loc[(cell, tf, name)] for name in METHODS]
            residual = tables["motif_residuals"].query("cell == @cell and tf == @tf and method == 'parametric_direct'")
            if len(residual) != 1:
                raise ValueError("candidate motif residual diagnostic is missing or ambiguous")
            row = qualify_row(task.parametric_direct_score.to_numpy(), task.DWM_expected_score.to_numpy(),
                              task.total_cuts.to_numpy(), common[ix], disagreement[ix],
                              candidate_threshold=cutoffs[0].threshold, baseline_threshold=cutoffs[1].threshold,
                              threshold_eligible=all(bool(r.threshold_eligible) for r in cutoffs),
                              motif_residual_flag=bool(residual.motif_residual_flag.iloc[0]))
            row.update(cell=cell, tf=tf, dispersion=inputs["dispersion"],
                       candidate_threshold=cutoffs[0].threshold, baseline_threshold=cutoffs[1].threshold)
            rows.append(row)
    outputs = {}
    for name, table in {"qualified_safety": pd.DataFrame(rows), "covered_site_scores": pd.concat(site_frames)}.items():
        path = args.outdir / f"{name}.tsv"
        table.to_csv(path, sep="\t", index=False)
        outputs[name] = record(path)
    result = {"schema": SCHEMA_COVERAGE, "coverage_input_id": freeze["coverage_input_id"],
              "source_safety": record(args.safety_manifest), "configuration_id": source["configuration_id"],
              "threshold_id": source["threshold_id"], "operator": OPERATOR, "outputs": outputs,
              "evidence": "post_hoc_coverage_qualified_development_candidate_screen", "replicate": "NakedDNA_rep2",
              "new_holdout_labels_read": False, "final_untouched_safety": False, "promotion_eligible": False}
    result["coverage_id"] = canonical_id(result)
    write_immutable(manifest_path, result)
    print(pd.DataFrame(rows)[["cell", "tf", "informative_support", "false_positive_calls", "paired_support",
                             "paired_false_positive_calls", "paired_dwm_calls", "passes_candidate_screen"]].to_string(index=False))
    print(manifest_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
