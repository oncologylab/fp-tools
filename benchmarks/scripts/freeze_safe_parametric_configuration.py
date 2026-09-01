#!/usr/bin/env python3
"""Create an immutable factorization freeze after residual safety screening."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(SCRIPT_DIR))

from evaluate_factorization_residual_safety import SCHEMA as SAFETY_SCHEMA  # noqa: E402
from evaluate_parametric_factorization import (  # noqa: E402
    RESIDUAL_TIE_ORDER,
    sha256_file,
)


SCHEMA = "fp-tools-parametric-factorization-configuration-freeze-v2"
PROVISIONAL_SCHEMA = "fp-tools-parametric-factorization-configuration-freeze-v1"


def _truth(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "1", "yes"}


def verify_record(record: dict[str, str], *, purpose: str) -> Path:
    path = Path(record["path"])
    if not path.is_file():
        raise ValueError(f"missing {purpose}: {path}")
    if sha256_file(path) != record["sha256"]:
        raise ValueError(f"checksum mismatch for {purpose}: {path}")
    return path


def select_safe_residual(
    residual_selection: pd.DataFrame,
    residual_safety: pd.DataFrame,
) -> tuple[str, pd.DataFrame]:
    required_selection = {
        "residual",
        "mean_relative_auprc_gain",
        "standard_error",
        "passes_ctcf_gate",
    }
    required_safety = {"residual", "passes_naked_dna_safety"}
    missing = required_selection.difference(residual_selection.columns)
    if missing:
        raise ValueError(
            "residual selection is missing columns: " + ", ".join(sorted(missing))
        )
    missing = required_safety.difference(residual_safety.columns)
    if missing:
        raise ValueError(
            "residual safety is missing columns: " + ", ".join(sorted(missing))
        )
    if residual_selection["residual"].duplicated().any():
        raise ValueError("residual selection contains duplicate residuals")
    if residual_safety["residual"].duplicated().any():
        raise ValueError("residual safety contains duplicate residuals")
    merged = residual_selection.merge(
        residual_safety,
        on="residual",
        how="left",
        validate="one_to_one",
    )
    if merged["passes_naked_dna_safety"].isna().any():
        missing_residuals = merged.loc[
            merged["passes_naked_dna_safety"].isna(), "residual"
        ].astype(str)
        raise ValueError(
            "residual safety is incomplete: " + ", ".join(missing_residuals)
        )
    merged["passes_ctcf_gate"] = merged["passes_ctcf_gate"].map(_truth)
    merged["passes_naked_dna_safety"] = merged[
        "passes_naked_dna_safety"
    ].map(_truth)
    merged["eligible_for_freeze"] = (
        merged["passes_ctcf_gate"] & merged["passes_naked_dna_safety"]
    )
    eligible = merged[merged["eligible_for_freeze"]].copy()
    if eligible.empty:
        raise RuntimeError(
            "no residual passed both CTCF non-regression and naked-DNA safety"
        )
    best = eligible.sort_values(
        "mean_relative_auprc_gain", ascending=False, kind="mergesort"
    ).iloc[0]
    threshold = float(best["mean_relative_auprc_gain"] - best["standard_error"])
    eligible_within = eligible[
        eligible["mean_relative_auprc_gain"] >= threshold
    ]
    order = {name: index for index, name in enumerate(RESIDUAL_TIE_ORDER)}
    unknown = set(eligible_within["residual"].astype(str)).difference(order)
    if unknown:
        raise ValueError("unknown residuals in freeze: " + ", ".join(sorted(unknown)))
    selected = min(
        eligible_within["residual"].astype(str), key=lambda value: order[value]
    )
    merged["within_one_se_after_safety"] = (
        merged["eligible_for_freeze"]
        & (merged["mean_relative_auprc_gain"] >= threshold)
    )
    merged["selected_after_safety"] = merged["residual"].astype(str).eq(selected)
    return selected, merged.sort_values(
        ["selected_after_safety", "eligible_for_freeze", "mean_relative_auprc_gain"],
        ascending=[False, False, False],
        kind="mergesort",
    )


def build_configuration(
    provisional_path: Path,
    selection_path: Path,
    safety_path: Path,
    output: Path,
) -> dict:
    provisional = json.loads(provisional_path.read_text(encoding="utf-8"))
    if provisional.get("schema") != PROVISIONAL_SCHEMA:
        raise ValueError("unsupported provisional configuration schema")
    if provisional.get("test_labels_opened") is not False:
        raise ValueError("provisional configuration has already opened test labels")
    for record in (
        provisional["factorization_model"],
        provisional["factorization_model_metadata"],
        provisional["bias_calibration"],
        provisional["bias_calibration_metadata"],
        provisional["study"],
        *provisional["inputs"],
    ):
        verify_record(record, purpose="provisional input")
    safety = json.loads(safety_path.read_text(encoding="utf-8"))
    if safety.get("schema") != SAFETY_SCHEMA:
        raise ValueError("unsupported residual-safety schema")
    if safety.get("naked_dna_labels_used") is not False:
        raise ValueError("residual safety does not certify label-free naked DNA")
    for key in (
        "factorization_model",
        "factorization_model_metadata",
        "thresholds",
        "detail",
        "summary",
        "site_scores",
    ):
        verify_record(safety[key], purpose=f"residual safety {key}")
    for record in safety["validation_inputs"] + safety["naked_dna_inputs"]:
        verify_record(record, purpose="residual safety input")
    if safety["factorization_model"] != provisional["factorization_model"]:
        raise ValueError("residual safety used a different factorization model")
    if safety["factorization_model_metadata"] != provisional[
        "factorization_model_metadata"
    ]:
        raise ValueError("residual safety used different model metadata")
    selection = pd.read_csv(selection_path, sep="\t")
    safety_summary = pd.read_csv(safety["summary"]["path"], sep="\t")
    selected, audited = select_safe_residual(selection, safety_summary)
    audited_path = output.parent / "factorization_residual_selection_safe.tsv"
    audited_path.parent.mkdir(parents=True, exist_ok=True)
    audited.to_csv(audited_path, sep="\t", index=False)
    if selected not in set(safety.get("passing_residuals", [])):
        raise ValueError("selected residual is absent from the safety manifest pass list")
    document = {
        **{
            key: provisional[key]
            for key in (
                "factorization_model",
                "factorization_model_metadata",
                "bias_calibration",
                "bias_calibration_metadata",
                "study",
                "inputs",
            )
        },
        "schema": SCHEMA,
        "base_configuration_id": provisional["configuration_id"],
        "provisional_configuration": {
            "path": str(provisional_path),
            "sha256": sha256_file(provisional_path),
        },
        "residual_selection": {
            "path": str(selection_path),
            "sha256": sha256_file(selection_path),
        },
        "residual_safety": {
            "path": str(safety_path),
            "sha256": sha256_file(safety_path),
            "safety_id": safety["safety_id"],
        },
        "safe_selection_audit": {
            "path": str(audited_path),
            "sha256": sha256_file(audited_path),
        },
        "selected_residual": selected,
        "safety_qualified": True,
        "selection_rule": (
            "best mean difficult-task relative AUPRC among CTCF- and naked-DNA-"
            "eligible residuals; fixed simplicity order within one standard error"
        ),
        "test_labels_opened": False,
    }
    canonical = json.dumps(document, sort_keys=True, separators=(",", ":"))
    document["configuration_id"] = sha256(canonical.encode()).hexdigest()
    rendered = json.dumps(document, indent=2, sort_keys=True) + "\n"
    if output.exists() and output.read_text(encoding="utf-8") != rendered:
        raise ValueError(f"refusing to replace a different immutable freeze: {output}")
    output.write_text(rendered, encoding="utf-8")
    return document


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--provisional-configuration", type=Path, required=True)
    parser.add_argument("--residual-selection", type=Path, required=True)
    parser.add_argument("--residual-safety", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    build_configuration(
        args.provisional_configuration,
        args.residual_selection,
        args.residual_safety,
        args.out,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
