#!/usr/bin/env python3
"""Validate the locked fp-tools footprint-detectability study specification."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REQUIRED_TASK_FIELDS = {
    "cell",
    "tf",
    "motif_id",
    "motif_family",
    "role",
    "split",
}
ALLOWED_SPLITS = {"development", "locked_holdout"}
ALLOWED_ROLES = {
    "positive_control",
    "difficult",
    "weak_shape",
    "motif_ambiguity",
    "transfer",
}
REQUIRED_DIAGNOSTIC_THRESHOLDS = {
    "minimum_positive_sites",
    "correction_delta_auroc",
    "scorer_delta_auroc",
    "scorer_relative_delta_auprc",
    "information_limit_auroc",
    "weak_auroc",
    "detectable_auroc",
}


def load_spec(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def validate_spec(spec: dict) -> list[str]:
    errors: list[str] = []
    for key in (
        "schema_version",
        "study_id",
        "assembly",
        "random_seed",
        "development_cells",
        "locked_holdout_cells",
        "chromosome_split",
        "depth_fragments",
        "diagnostic_thresholds",
        "promotion_gates",
        "negative_control_resources",
        "tasks",
        "nutrient_application",
    ):
        if key not in spec:
            errors.append(f"missing required top-level field: {key}")

    if errors:
        return errors
    if spec["schema_version"] != 1:
        errors.append("schema_version must be 1")
    if spec["assembly"] != "GRCh38":
        errors.append("assembly must be GRCh38 for v1")
    if not isinstance(spec["random_seed"], int):
        errors.append("random_seed must be an integer")
    negative_controls = spec["negative_control_resources"]
    if not negative_controls.get("naked_dna_accession") or not negative_controls.get("runs"):
        errors.append("negative_control_resources must lock naked-DNA accession and runs")

    chromosome_split = spec["chromosome_split"]
    expected_split_names = {"train", "validation", "test"}
    if set(chromosome_split) != expected_split_names:
        errors.append("chromosome_split must contain train, validation, and test")
    else:
        seen: set[str] = set()
        for split_name in ("train", "validation", "test"):
            chromosomes = chromosome_split[split_name]
            duplicates = seen.intersection(chromosomes)
            if duplicates:
                errors.append(
                    f"chromosomes assigned to multiple splits: {', '.join(sorted(duplicates))}"
                )
            seen.update(chromosomes)

    depth_values = spec["depth_fragments"]
    if not depth_values or depth_values[-1] != "full":
        errors.append("depth_fragments must end with 'full'")
    numeric_depths = [value for value in depth_values if value != "full"]
    if any(not isinstance(value, int) or value <= 0 for value in numeric_depths):
        errors.append("numeric depth_fragments must be positive integers")
    if numeric_depths != sorted(set(numeric_depths)):
        errors.append("numeric depth_fragments must be unique and increasing")

    missing_thresholds = REQUIRED_DIAGNOSTIC_THRESHOLDS.difference(
        spec["diagnostic_thresholds"]
    )
    if missing_thresholds:
        errors.append(
            "diagnostic_thresholds missing: " + ", ".join(sorted(missing_thresholds))
        )

    development_cells = set(spec["development_cells"])
    holdout_cells = set(spec["locked_holdout_cells"])
    overlap = development_cells.intersection(holdout_cells)
    if overlap:
        errors.append(
            "development and holdout cells overlap: " + ", ".join(sorted(overlap))
        )

    task_keys: set[tuple[str, str, str]] = set()
    for index, task in enumerate(spec["tasks"]):
        missing = REQUIRED_TASK_FIELDS.difference(task)
        if missing:
            errors.append(f"task {index} missing: {', '.join(sorted(missing))}")
            continue
        key = (task["cell"], task["tf"], task["motif_id"])
        if key in task_keys:
            errors.append(f"duplicate task: {'/'.join(key)}")
        task_keys.add(key)
        if task["split"] not in ALLOWED_SPLITS:
            errors.append(f"task {'/'.join(key)} has invalid split {task['split']}")
        if task["role"] not in ALLOWED_ROLES:
            errors.append(f"task {'/'.join(key)} has invalid role {task['role']}")
        expected_cells = development_cells if task["split"] == "development" else holdout_cells
        if task["cell"] not in expected_cells:
            errors.append(
                f"task {'/'.join(key)} cell is inconsistent with split {task['split']}"
            )

    for cell_group, split_name in (
        (development_cells, "development"),
        (holdout_cells, "locked_holdout"),
    ):
        represented = {
            task["cell"] for task in spec["tasks"] if task.get("split") == split_name
        }
        missing_cells = cell_group.difference(represented)
        if missing_cells:
            errors.append(
                f"{split_name} cells without tasks: {', '.join(sorted(missing_cells))}"
            )

    nutrient = spec["nutrient_application"]
    for key in (
        "discovery_cell",
        "transfer_cells",
        "primary_contrast",
        "external_pdac_accession",
        "external_mechanistic_accession",
        "external_resources",
        "candidate_rules",
    ):
        if key not in nutrient:
            errors.append(f"nutrient_application missing: {key}")
    resources = nutrient.get("external_resources", {})
    for accession_key in ("external_pdac_accession", "external_mechanistic_accession"):
        accession = nutrient.get(accession_key)
        if accession and accession not in resources:
            errors.append(f"external_resources missing locked accession: {accession}")
    pdac = resources.get(nutrient.get("external_pdac_accession"), {})
    for assay in ("atac_runs", "rna_runs"):
        states = pdac.get(assay, {})
        if set(states) != {"non_adapted", "adapted", "reverse_adapted"}:
            errors.append(f"external PDAC {assay} must lock all three adaptation states")
            continue
        runs = [run for state_runs in states.values() for run in state_runs]
        if len(runs) != len(set(runs)):
            errors.append(f"external PDAC {assay} contains duplicate run accessions")
    return errors


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--spec",
        type=Path,
        default=Path("benchmarks/manifests/footprint_detectability_v1.spec.json"),
    )
    args = parser.parse_args(argv)
    errors = validate_spec(load_spec(args.spec))
    if errors:
        for error in errors:
            print(error)
        return 1
    print(f"Validated footprint study specification: {args.spec}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
