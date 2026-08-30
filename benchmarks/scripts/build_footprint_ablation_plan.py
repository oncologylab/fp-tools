#!/usr/bin/env python3
"""Build executable depth/correction and locked method-evaluation plans."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shlex
import sys

import pandas as pd


REQUIRED_SAMPLE_COLUMNS = ("sample", "cell", "bam", "peaks", "fragments")


def load_study(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def load_samples(path: Path, check_paths: bool = False) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t")
    missing = [column for column in REQUIRED_SAMPLE_COLUMNS if column not in frame]
    if missing:
        raise ValueError(f"sample manifest is missing columns: {', '.join(missing)}")
    if frame["sample"].duplicated().any():
        raise ValueError("sample manifest contains duplicate sample names")
    frame["fragments"] = pd.to_numeric(frame["fragments"], errors="raise").astype(int)
    if (frame["fragments"] <= 0).any():
        raise ValueError("sample fragment counts must be positive")
    if check_paths:
        absent = [
            str(value)
            for column in ("bam", "peaks")
            for value in frame[column]
            if not Path(value).exists()
        ]
        if absent:
            raise FileNotFoundError("input paths do not exist: " + ", ".join(absent[:5]))
    return frame


def command_text(arguments: list[str]) -> str:
    return shlex.join(arguments)


def depth_label(depth: int | str) -> str:
    return "full" if depth == "full" else f"{int(depth) // 1_000_000}m"


def build_signal_plan(
    study: dict,
    samples: pd.DataFrame,
    genome: Path,
    blacklist: Path,
    outdir: Path,
    python: str = sys.executable,
    cores: int = 1,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    downsampler = Path(__file__).with_name("downsample_bam_fragments.py")
    base_seed = int(study["random_seed"])
    randomizations = int(study["depth_randomizations"])
    corrections = list(study["native_corrections"])
    for sample_row in samples.sort_values(["cell", "sample"]).itertuples(index=False):
        available = int(sample_row.fragments)
        for depth in study["depth_fragments"]:
            if depth != "full" and int(depth) > available:
                continue
            seeds = [base_seed + index for index in range(randomizations)] if depth != "full" else [base_seed]
            for seed in seeds:
                label = depth_label(depth)
                subset_id = f"{sample_row.sample}.{label}.s{seed}"
                subset_dir = outdir / "signals" / str(sample_row.sample) / label / f"seed_{seed}"
                if depth == "full":
                    bam = Path(str(sample_row.bam))
                    dependency = ""
                else:
                    bam = subset_dir / f"{subset_id}.bam"
                    dependency = f"downsample:{subset_id}"
                    downsample_command = [
                        python,
                        str(downsampler),
                        "--bam",
                        str(sample_row.bam),
                        "--out",
                        str(bam),
                        "--target-fragments",
                        str(depth),
                        "--seed",
                        str(seed),
                    ]
                    rows.append(
                        {
                            "job_id": dependency,
                            "stage": "downsample",
                            "sample": sample_row.sample,
                            "cell": sample_row.cell,
                            "depth": depth,
                            "seed": seed,
                            "correction": "",
                            "depends_on": "",
                            "expected_output": str(bam),
                            "command": command_text(downsample_command),
                        }
                    )

                dwm_job = f"correct:{subset_id}:fp_tools_dwm"
                dwm_dir = subset_dir / "fp_tools_dwm"
                dwm_command = [
                    "atac-correct",
                    "--bams",
                    str(bam),
                    "--genome",
                    str(genome),
                    "--peaks",
                    str(sample_row.peaks),
                    "--blacklist",
                    str(blacklist),
                    "--outdir",
                    str(dwm_dir),
                    "--prefix",
                    subset_id,
                    "--score_mat",
                    "DWM",
                    "--write-tracks",
                    "all",
                    "--cores",
                    str(cores),
                ]
                rows.append(
                    {
                        "job_id": dwm_job,
                        "stage": "correction",
                        "sample": sample_row.sample,
                        "cell": sample_row.cell,
                        "depth": depth,
                        "seed": seed,
                        "correction": "fp_tools_dwm",
                        "depends_on": dependency,
                        "expected_output": str(dwm_dir / f"{subset_id}_corrected.bw"),
                        "command": command_text(dwm_command),
                    }
                )
                if "raw" in corrections:
                    rows.append(
                        {
                            "job_id": f"signal:{subset_id}:raw",
                            "stage": "derived_signal",
                            "sample": sample_row.sample,
                            "cell": sample_row.cell,
                            "depth": depth,
                            "seed": seed,
                            "correction": "raw",
                            "depends_on": dwm_job,
                            "expected_output": str(dwm_dir / f"{subset_id}_uncorrected.bw"),
                            "command": "",
                        }
                    )
                if "fp_tools_pwm" in corrections:
                    pwm_dir = subset_dir / "fp_tools_pwm"
                    pwm_command = dwm_command.copy()
                    pwm_command[pwm_command.index(str(dwm_dir))] = str(pwm_dir)
                    pwm_command[pwm_command.index("DWM")] = "PWM"
                    rows.append(
                        {
                            "job_id": f"correct:{subset_id}:fp_tools_pwm",
                            "stage": "correction",
                            "sample": sample_row.sample,
                            "cell": sample_row.cell,
                            "depth": depth,
                            "seed": seed,
                            "correction": "fp_tools_pwm",
                            "depends_on": dependency,
                            "expected_output": str(pwm_dir / f"{subset_id}_corrected.bw"),
                            "command": command_text(pwm_command),
                        }
                    )
                if "fp_tools_reused_bias" in corrections and depth != "full":
                    reuse_dir = subset_dir / "fp_tools_reused_bias"
                    full_id = f"{sample_row.sample}.full.s{base_seed}"
                    full_bias = (
                        outdir
                        / "signals"
                        / str(sample_row.sample)
                        / "full"
                        / f"seed_{base_seed}"
                        / "fp_tools_dwm"
                        / f"{full_id}_AtacBias.pickle"
                    )
                    reuse_command = dwm_command.copy()
                    reuse_command[reuse_command.index(str(dwm_dir))] = str(reuse_dir)
                    reuse_command.extend(["--bias-pkl", str(full_bias)])
                    rows.append(
                        {
                            "job_id": f"correct:{subset_id}:fp_tools_reused_bias",
                            "stage": "correction",
                            "sample": sample_row.sample,
                            "cell": sample_row.cell,
                            "depth": depth,
                            "seed": seed,
                            "correction": "fp_tools_reused_bias",
                            "depends_on": ";".join(filter(None, [dependency, f"correct:{full_id}:fp_tools_dwm"])),
                            "expected_output": str(reuse_dir / f"{subset_id}_corrected.bw"),
                            "command": command_text(reuse_command),
                        }
                    )
    return pd.DataFrame(rows)


def build_evaluation_plan(study: dict, signal_plan: pd.DataFrame) -> pd.DataFrame:
    signals = signal_plan[signal_plan["stage"].isin(["correction", "derived_signal"])].copy()
    signals = signals.rename(columns={"expected_output": "signal"})
    rows: list[dict[str, object]] = []
    for task in study["tasks"]:
        matching = signals[signals["cell"] == task["cell"]]
        for signal in matching.itertuples(index=False):
            for method in study["whole_methods"]:
                rows.append(
                    {
                        "evaluation_id": f"{task['cell']}:{task['tf']}:{signal.sample}:{signal.depth}:s{signal.seed}:{signal.correction}:{method}",
                        **task,
                        "sample": signal.sample,
                        "depth": signal.depth,
                        "seed": signal.seed,
                        "correction": signal.correction,
                        "method": method,
                        "signal": signal.signal,
                        "depends_on": signal.job_id,
                    }
                )
    return pd.DataFrame(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, default=Path("benchmarks/manifests/footprint_detectability_v1.spec.json"))
    parser.add_argument("--samples", type=Path, required=True)
    parser.add_argument("--genome", type=Path, required=True)
    parser.add_argument("--blacklist", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--cores", type=int, default=1)
    parser.add_argument("--check-paths", action="store_true")
    args = parser.parse_args(argv)

    study = load_study(args.study)
    samples = load_samples(args.samples, check_paths=args.check_paths)
    signal_plan = build_signal_plan(
        study,
        samples,
        args.genome,
        args.blacklist,
        args.outdir,
        python=args.python,
        cores=args.cores,
    )
    evaluation_plan = build_evaluation_plan(study, signal_plan)
    args.outdir.mkdir(parents=True, exist_ok=True)
    signal_path = args.outdir / "ablation_commands.tsv"
    evaluation_path = args.outdir / "evaluation_matrix.tsv"
    signal_plan.to_csv(signal_path, sep="\t", index=False)
    evaluation_plan.to_csv(evaluation_path, sep="\t", index=False)
    print(f"wrote {len(signal_plan):,} signal jobs to {signal_path}")
    print(f"wrote {len(evaluation_plan):,} evaluation tasks to {evaluation_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
