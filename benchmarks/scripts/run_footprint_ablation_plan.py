#!/usr/bin/env python3
"""Execute or inspect a footprint ablation command plan with dependency checks."""

from __future__ import annotations

import argparse
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from datetime import datetime, timezone
from pathlib import Path
import re
import shlex
import subprocess

import numpy as np
import pandas as pd


REQUIRED_COLUMNS = ["job_id", "stage", "depends_on", "expected_output", "command"]


def order_plan(frame: pd.DataFrame, order: str) -> pd.DataFrame:
    """Order jobs for either strict plan fidelity or breadth-first coverage."""

    if order == "plan":
        return frame.reset_index(drop=True)
    if order != "breadth-first":
        raise ValueError(f"unknown scheduling order: {order}")
    required = {"sample", "depth", "seed", "stage"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(
            "breadth-first scheduling requires columns: " + ", ".join(sorted(missing))
        )
    ordered = frame.copy()
    numeric_depth = pd.to_numeric(ordered["depth"], errors="coerce")
    finite_depths = numeric_depth[np.isfinite(numeric_depth)]
    full_rank = float(finite_depths.max() + 1) if len(finite_depths) else 1.0
    ordered["_depth_order"] = numeric_depth.fillna(full_rank)
    ordered["_seed_order"] = pd.to_numeric(ordered["seed"], errors="coerce").fillna(-1)
    ordered["_stage_order"] = ordered["stage"].map(
        {"downsample": 0, "correction": 1, "derived_signal": 2}
    ).fillna(3)
    ordered["_plan_order"] = np.arange(len(ordered))
    ordered = ordered.sort_values(
        [
            "_depth_order",
            "_seed_order",
            "sample",
            "_stage_order",
            "correction",
            "_plan_order",
        ],
        kind="mergesort",
    )
    return ordered.drop(
        columns=["_depth_order", "_seed_order", "_stage_order", "_plan_order"]
    ).reset_index(drop=True)


def dependencies(value: object) -> list[str]:
    if pd.isna(value) or not str(value).strip():
        return []
    return [item for item in str(value).split(";") if item]


def validate_plan(frame: pd.DataFrame) -> None:
    missing = [column for column in REQUIRED_COLUMNS if column not in frame]
    if missing:
        raise ValueError(f"ablation plan is missing columns: {', '.join(missing)}")
    if frame["job_id"].duplicated().any():
        raise ValueError("ablation plan contains duplicate job IDs")
    known = set(frame["job_id"])
    unknown = sorted(
        {
            dependency
            for value in frame["depends_on"]
            for dependency in dependencies(value)
            if dependency not in known
        }
    )
    if unknown:
        raise ValueError("ablation plan has unknown dependencies: " + ", ".join(unknown[:5]))


def execute_plan(
    frame: pd.DataFrame,
    status_path: Path,
    selected_jobs: set[str] | None = None,
    dry_run: bool = False,
    workers: int = 1,
    log_dir: Path | None = None,
) -> pd.DataFrame:
    validate_plan(frame)
    if workers < 1:
        raise ValueError("workers must be positive")
    selected_jobs = set(frame["job_id"]) if selected_jobs is None else selected_jobs
    known = frame.set_index("job_id")
    pending = [job for job in frame["job_id"] if job in selected_jobs]
    completed = {
        job
        for job, row in known.iterrows()
        if Path(str(row.expected_output)).exists()
    }
    records: list[dict[str, object]] = []
    running = {}
    executor = ThreadPoolExecutor(max_workers=workers)

    def record(job: str, row: pd.Series, state: str, timestamp: str, job_log: Path | None) -> None:
        completed.add(job)
        if job in pending:
            pending.remove(job)
        records.append(
            {
                "job_id": job,
                "stage": row.stage,
                "state": state,
                "expected_output": str(row.expected_output),
                "log": str(job_log) if job_log is not None else "",
                "timestamp_utc": timestamp,
            }
        )
        status_path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(records).to_csv(status_path, sep="\t", index=False)
        print(f"{state}\t{job}", flush=True)

    def run_command(job: str, row: pd.Series) -> tuple[str, Path | None]:
        if not str(row.command).strip() or str(row.command) == "nan":
            raise RuntimeError(
                f"derived job {job} is missing expected output {row.expected_output}"
            )
        job_log = None
        if log_dir is not None:
            log_dir.mkdir(parents=True, exist_ok=True)
            safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", job)
            job_log = log_dir / f"{safe_name}.log"
            with job_log.open("w", encoding="utf-8") as handle:
                subprocess.run(
                    shlex.split(str(row.command)),
                    check=True,
                    stdout=handle,
                    stderr=subprocess.STDOUT,
                )
        else:
            subprocess.run(shlex.split(str(row.command)), check=True)
        output = Path(str(row.expected_output))
        if not output.exists():
            raise RuntimeError(f"job {job} completed without expected output {output}")
        return "completed", job_log

    try:
        while pending or running:
            progressed = False
            for job in pending.copy():
                if len(running) >= workers:
                    break
                row = known.loc[job]
                required = dependencies(row.depends_on)
                if not all(dependency in completed for dependency in required):
                    continue
                output = Path(str(row.expected_output))
                timestamp = datetime.now(timezone.utc).isoformat()
                if output.exists():
                    record(job, row, "skipped_existing", timestamp, None)
                elif dry_run:
                    record(job, row, "planned", timestamp, None)
                else:
                    future = executor.submit(run_command, job, row)
                    running[future] = (job, row, timestamp)
                    pending.remove(job)
                progressed = True

            if running:
                finished, _ = wait(tuple(running), return_when=FIRST_COMPLETED)
                for future in finished:
                    job, row, timestamp = running.pop(future)
                    try:
                        state, job_log = future.result()
                    except Exception as error:
                        for active in running:
                            active.cancel()
                        raise RuntimeError(f"ablation job failed: {job}") from error
                    record(job, row, state, timestamp, job_log)
                    progressed = True
            if not progressed and pending:
                blocked = ", ".join(pending[:5])
                raise RuntimeError(
                    f"no runnable jobs remain; dependencies are absent or were not selected: {blocked}"
                )
    finally:
        executor.shutdown(wait=True, cancel_futures=True)
    return pd.DataFrame(records)


def filter_jobs(
    frame: pd.DataFrame,
    *,
    job_ids: list[str],
    stages: list[str],
    samples: list[str],
    depths: list[str],
    corrections: list[str],
) -> set[str] | None:
    """Select jobs by matrix fields, then close the selection over dependencies."""

    if not any((job_ids, stages, samples, depths, corrections)):
        return None
    mask = pd.Series(True, index=frame.index)
    if stages:
        mask &= frame["stage"].astype(str).isin(stages)
    if samples:
        mask &= frame["sample"].astype(str).isin(samples)
    if depths:
        mask &= frame["depth"].astype(str).isin(depths)
    if corrections:
        mask &= frame["correction"].astype(str).isin(corrections)
    selected = set(frame.loc[mask, "job_id"]).union(job_ids)
    known = frame.set_index("job_id")
    required = list(selected)
    while required:
        job = required.pop()
        if job not in known.index:
            raise ValueError(f"unknown selected job: {job}")
        for dependency in dependencies(known.loc[job].depends_on):
            if dependency not in selected:
                selected.add(dependency)
                required.append(dependency)
    return selected


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--status", type=Path)
    parser.add_argument("--job-id", action="append", default=[])
    parser.add_argument("--stage", action="append", default=[])
    parser.add_argument("--sample", action="append", default=[])
    parser.add_argument("--depth", action="append", default=[])
    parser.add_argument("--correction", action="append", default=[])
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument(
        "--order",
        choices=("plan", "breadth-first"),
        default="plan",
        help="Breadth-first alternates cells/seeds at lower depths before deeper jobs.",
    )
    parser.add_argument("--log-dir", type=Path)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    frame = pd.read_csv(args.plan, sep="\t", keep_default_na=False)
    selected = filter_jobs(
        frame,
        job_ids=args.job_id,
        stages=args.stage,
        samples=args.sample,
        depths=args.depth,
        corrections=args.correction,
    )
    status = args.status or args.plan.with_name("ablation_status.tsv")
    execute_plan(
        order_plan(frame, args.order),
        status,
        selected_jobs=selected,
        dry_run=args.dry_run,
        workers=args.workers,
        log_dir=args.log_dir,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
