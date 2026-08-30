#!/usr/bin/env python3
"""Execute or inspect a footprint ablation command plan with dependency checks."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path
import shlex
import subprocess

import pandas as pd


REQUIRED_COLUMNS = ["job_id", "stage", "depends_on", "expected_output", "command"]


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
) -> pd.DataFrame:
    validate_plan(frame)
    selected_jobs = selected_jobs or set(frame["job_id"])
    known = frame.set_index("job_id")
    pending = [job for job in frame["job_id"] if job in selected_jobs]
    completed = {
        job
        for job, row in known.iterrows()
        if Path(str(row.expected_output)).exists()
    }
    records: list[dict[str, object]] = []
    while pending:
        progressed = False
        for job in pending.copy():
            row = known.loc[job]
            required = dependencies(row.depends_on)
            if not all(dependency in completed for dependency in required):
                continue
            output = Path(str(row.expected_output))
            timestamp = datetime.now(timezone.utc).isoformat()
            if output.exists():
                state = "skipped_existing"
            elif dry_run:
                state = "planned"
            elif not str(row.command).strip() or str(row.command) == "nan":
                raise RuntimeError(f"derived job {job} is missing expected output {output}")
            else:
                subprocess.run(shlex.split(str(row.command)), check=True)
                if not output.exists():
                    raise RuntimeError(f"job {job} completed without expected output {output}")
                state = "completed"
            completed.add(job)
            pending.remove(job)
            progressed = True
            records.append(
                {
                    "job_id": job,
                    "stage": row.stage,
                    "state": state,
                    "expected_output": str(output),
                    "timestamp_utc": timestamp,
                }
            )
            status_path.parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame(records).to_csv(status_path, sep="\t", index=False)
            print(f"{state}\t{job}")
        if not progressed:
            blocked = ", ".join(pending[:5])
            raise RuntimeError(
                f"no runnable jobs remain; dependencies are absent or were not selected: {blocked}"
            )
    return pd.DataFrame(records)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--status", type=Path)
    parser.add_argument("--job-id", action="append", default=[])
    parser.add_argument("--stage", action="append", default=[])
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    frame = pd.read_csv(args.plan, sep="\t", keep_default_na=False)
    selected: set[str] | None = None
    if args.job_id or args.stage:
        selected = set(args.job_id)
        if args.stage:
            selected.update(frame.loc[frame["stage"].isin(args.stage), "job_id"])
        required = list(selected)
        while required:
            job = required.pop()
            row = frame.set_index("job_id").loc[job]
            for dependency in dependencies(row.depends_on):
                if dependency not in selected:
                    selected.add(dependency)
                    required.append(dependency)
    status = args.status or args.plan.with_name("ablation_status.tsv")
    execute_plan(frame, status, selected_jobs=selected, dry_run=args.dry_run)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
