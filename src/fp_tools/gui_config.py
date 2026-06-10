"""Shared config helpers for optional GUI and YAML-driven batch execution.

This module defines the normalized config shape used by the GUI and the
optional ``fp-tools-run --config ...`` path. Direct CLI commands remain the
primary interface and do not depend on this layer.
"""

from __future__ import annotations

from collections.abc import Mapping
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

CONFIG_VERSION = 1

TOOL_ALIASES = {
    "atac-correct": "ATACorrect",
    "atacorrect": "ATACorrect",
    "score-footprints": "FootprintScores",
    "footprintscores": "FootprintScores",
    "scorebigwig": "FootprintScores",
    "score-bigwig": "FootprintScores",
    "detect-tf-binding": "BINDetect",
    "bindetect": "BINDetect",
    "plot-aggregate": "PlotAggregate",
    "plotaggregate": "PlotAggregate",
}

RESERVED_KEYS = {
    "sample_id",
    "comparison_id",
    "job_id",
    "tool",
    "label",
    "name",
    "description",
}

LIST_FLAGS = {
    "signals",
    "cond_names",
    "cond-names",
    "tfbs",
    "regions",
}

REQUIRED_FIELDS = {
    "ATACorrect": ("bam", "genome", "peaks"),
    "FootprintScores": ("signal", "regions", "output"),
    "BINDetect": ("motifs", "signals", "genome", "peaks"),
    "PlotAggregate": ("TFBS", "signals", "output"),
}

FLAG_NAME_MAP = {
    "peak_header": "--peak-header",
    "cond_names": "--cond-names",
    "output_txt": "--output-txt",
    "output_csv": "--output-csv",
    "output_aggregated_signals": "--output_aggregated_signals",
    "output_aggregated_scores": "--output_aggregated_scores",
    "control_label": "--control-label",
    "TFBS_labels": "--TFBS-labels",
    "signal_labels": "--signal-labels",
    "region_labels": "--region-labels",
    "share_y": "--share-y",
    "log_transform": "--log-transform",
    "plot_boundaries": "--plot-boundaries",
    "signal_on_x": "--signal-on-x",
    "remove_outliers": "--remove-outliers",
}


@dataclass
class JobSpec:
    job_id: str
    section: str
    tool: str
    params: dict[str, Any]
    command: list[str]


def canonical_tool_name(name: str) -> str:
    key = str(name).strip().lower()
    if key not in TOOL_ALIASES:
        raise ValueError(f"Unsupported tool in config: {name}")
    return TOOL_ALIASES[key]


def load_yaml_config(path: str | Path) -> dict[str, Any]:
    with Path(path).expanduser().open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, Mapping):
        raise ValueError("Top-level YAML config must be a mapping.")
    return dict(data)


def dump_yaml_config(config: Mapping[str, Any], path: str | Path) -> None:
    path = Path(path).expanduser()
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(
            dict(config),
            handle,
            sort_keys=False,
            default_flow_style=False,
            allow_unicode=False,
        )


def make_single_config(tool: str, params: Mapping[str, Any], job_id: str = "run") -> dict[str, Any]:
    return {
        "version": CONFIG_VERSION,
        "run_mode": "single",
        "defaults": {},
        "samples": [
            {
                "sample_id": job_id,
                "tool": canonical_tool_name(tool),
                **dict(params),
            }
        ],
        "comparisons": [],
    }


def normalize_config(config: Mapping[str, Any]) -> dict[str, Any]:
    raw = deepcopy(dict(config))
    if "tool" in raw and "samples" not in raw and "comparisons" not in raw:
        tool = raw.pop("tool")
        job_id = raw.pop("job_id", "run")
        return make_single_config(tool, raw, job_id=job_id)

    version = int(raw.get("version", CONFIG_VERSION))
    defaults = raw.get("defaults", {}) or {}
    samples = raw.get("samples", []) or []
    comparisons = raw.get("comparisons", []) or []
    run_mode = raw.get("run_mode", "batch" if len(samples) + len(comparisons) > 1 else "single")
    run_root = raw.get("run_root")

    if not isinstance(defaults, Mapping):
        raise ValueError("'defaults' must be a mapping when present.")
    if not isinstance(samples, list):
        raise ValueError("'samples' must be a list when present.")
    if not isinstance(comparisons, list):
        raise ValueError("'comparisons' must be a list when present.")
    if not samples and not comparisons:
        raise ValueError("Config must contain at least one item in 'samples' or 'comparisons'.")

    return {
        "version": version,
        "run_mode": run_mode,
        "run_root": run_root,
        "defaults": dict(defaults),
        "samples": [dict(item) for item in samples],
        "comparisons": [dict(item) for item in comparisons],
    }


def expand_jobs(config: Mapping[str, Any], only_tools: set[str] | None = None) -> list[JobSpec]:
    normalized = normalize_config(config)
    jobs: list[JobSpec] = []
    defaults = normalized["defaults"]

    for section in ("samples", "comparisons"):
        for idx, item in enumerate(normalized[section], start=1):
            if not isinstance(item, Mapping):
                raise ValueError(f"Each item in '{section}' must be a mapping.")
            merged = dict(defaults)
            merged.update(dict(item))
            tool = canonical_tool_name(str(merged.get("tool", "")))
            if only_tools and tool not in only_tools:
                continue

            job_id = str(
                merged.get("sample_id")
                or merged.get("comparison_id")
                or merged.get("job_id")
                or f"{tool.lower()}_{idx:03d}"
            )
            params = {k: v for k, v in merged.items() if k not in RESERVED_KEYS}
            command = build_cli_command(tool, params)
            jobs.append(JobSpec(job_id=job_id, section=section, tool=tool, params=params, command=command))
    return jobs


def build_cli_command(tool: str, params: Mapping[str, Any]) -> list[str]:
    command = [tool]
    extras = list(params.get("extra_args", []) or [])

    ordered_keys = [key for key in params.keys() if key != "extra_args"]
    for key in ordered_keys:
        value = params[key]
        if value is None or value == "":
            continue
        flag = _key_to_flag(key)
        if isinstance(value, bool):
            if value:
                command.append(flag)
            continue
        if isinstance(value, Mapping):
            raise ValueError(f"Nested mapping for '{key}' is not supported in YAML CLI configs.")
        if isinstance(value, list):
            if not value:
                continue
            command.append(flag)
            command.extend(str(item) for item in value)
            continue
        command.extend([flag, str(value)])

    command.extend(str(arg) for arg in extras)
    return command


def config_to_yaml_text(config: Mapping[str, Any]) -> str:
    normalized = normalize_config(config)
    return yaml.safe_dump(normalized, sort_keys=False, default_flow_style=False, allow_unicode=False)


def parse_yaml_text(text: str) -> dict[str, Any]:
    data = yaml.safe_load(text) or {}
    if not isinstance(data, Mapping):
        raise ValueError("YAML text must define a mapping.")
    return normalize_config(dict(data))


def validate_config(config: Mapping[str, Any]) -> list[str]:
    normalized = normalize_config(config)
    errors: list[str] = []

    for section in ("samples", "comparisons"):
        for idx, item in enumerate(normalized[section], start=1):
            tool = canonical_tool_name(str(item.get("tool", "")))
            required = REQUIRED_FIELDS.get(tool, ())
            job_name = str(
                item.get("sample_id")
                or item.get("comparison_id")
                or item.get("job_id")
                or f"{tool.lower()}_{idx:03d}"
            )
            for field in required:
                value = item.get(field)
                if isinstance(value, list):
                    if not [v for v in value if str(v).strip()]:
                        errors.append(f"{job_name}: missing required field '{field}'")
                elif str(value or "").strip() == "":
                    errors.append(f"{job_name}: missing required field '{field}'")

    return errors


def _key_to_flag(key: str) -> str:
    key = str(key).strip()
    if key.startswith("--"):
        return key
    if key in FLAG_NAME_MAP:
        return FLAG_NAME_MAP[key]
    return f"--{key.replace('_', '-')}"
