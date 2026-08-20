#!/usr/bin/env python
"""Run the complete bulk ATAC-seq footprinting workflow for explicit comparisons."""

from __future__ import annotations

import argparse
import json
import os
import shlex
import subprocess
import sys
from pathlib import Path

from fp_tools.platform_support import require_raw_read_preparation_support
from fp_tools.runtime import RuntimeProvisionError, add_runtime_argument, prepare_command_runtime
from fp_tools.utils.project_layout import (
    analysis_peaks_path,
    comparison_dir,
    corrected_bigwig_path,
    footprint_bigwig_path,
    match_motifs_dir,
    read_comparison_table,
    read_sample_table,
)
from fp_tools.utils.subprocess_commands import resolve_fp_tools_subprocess


def _quote(command: list[str]) -> str:
    return " ".join(shlex.quote(str(part)) for part in command)


def _run(command: list[str], stdout_path: Path, stderr_path: Path) -> int:
    resolved = resolve_fp_tools_subprocess(command)
    env = os.environ.copy()
    env.setdefault("MPLCONFIGDIR", str(stdout_path.parent / ".mplconfig"))
    env.setdefault("XDG_CACHE_HOME", str(stdout_path.parent / ".cache"))
    Path(env["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    Path(env["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)
    with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open("w", encoding="utf-8") as stderr:
        result = subprocess.run(resolved, stdout=stdout, stderr=stderr, text=True, env=env, check=False)
    return int(result.returncode)


def _all_nonempty(paths: list[Path]) -> bool:
    return bool(paths) and all(path.is_file() and path.stat().st_size > 0 for path in paths)


def _resolve_review_format(args: argparse.Namespace) -> str:
    plot_aggregate = str(getattr(args, "plot_aggregate", "all"))
    requested = str(getattr(args, "review_format", "auto"))
    resolved = "standalone" if requested == "auto" and plot_aggregate == "off" else requested
    if resolved == "auto":
        resolved = "bundle"
    if resolved == "bundle" and plot_aggregate == "off":
        raise ValueError(
            "--review-format bundle requires aggregate profiles; use "
            "--review-format standalone, auto, or none with --plot-aggregate off"
        )
    return resolved


def _stage_complete(
    stage: str,
    project: Path,
    samples,
    comparisons,
    review_format: str = "bundle",
) -> bool:
    if stage == "atac-correct":
        return analysis_peaks_path(project).is_file() and _all_nonempty(
            [corrected_bigwig_path(project, row.sample) for row in samples]
        )
    if stage == "call-footprints":
        return _all_nonempty([footprint_bigwig_path(project, row.sample) for row in samples])
    if stage == "match-motifs":
        return all(match_motifs_dir(project, row.sample).is_dir() for row in samples)
    if stage == "diff-footprints":
        return all(
            any(comparison_dir(project, row.comparison).glob("diff_footprints_*.html"))
            for row in comparisons
        )
    if stage == "review-multi-comparisons":
        if review_format == "standalone":
            return _all_nonempty([project / "reports" / "review_multi_comparisons.html"])
        return (project / "reports" / "review_multi_comparisons" / "index.html").is_file()
    return False


def build_commands(args: argparse.Namespace) -> list[tuple[str, list[str]]]:
    project = Path(args.outdir).expanduser().resolve()
    review_format = _resolve_review_format(args)
    shared = ["--sample-table", str(args.sample_table), "--layout", "project", "--outdir", str(project)]
    motif_args = ["--motif-db", str(args.motif_db)] if args.motif_db else []
    if args.motifs:
        motif_args.extend(["--motifs", *[str(path) for path in args.motifs]])
    atac = [
        "atac-correct",
        *shared,
        "--genome",
        str(args.genome),
        "--cores",
        str(args.cores),
    ]
    if args.blacklist:
        atac.extend(["--blacklist", str(args.blacklist)])
    footprints = ["call-footprints", *shared, "--cores", str(args.cores)]
    motifs = [
        "match-motifs",
        *shared,
        "--genome",
        str(args.genome),
        *motif_args,
        "--cores",
        str(args.cores),
    ]
    differential = [
        "diff-footprints",
        *shared,
        "--comparison-table",
        str(args.comparison_table),
        "--genome",
        str(args.genome),
        *motif_args,
        "--normalization",
        str(args.normalization),
        "--plot-aggregate",
        str(args.plot_aggregate),
        "--aggregate-site-set",
        "all",
        "--cores",
        str(args.cores),
    ]
    commands = [
        ("atac-correct", atac),
        ("call-footprints", footprints),
        ("match-motifs", motifs),
        ("diff-footprints", differential),
    ]
    if review_format != "none":
        review = [
            "review-multi-comparisons",
            "--inputs",
            str(project / "comparisons"),
        ]
        if review_format == "standalone":
            review.extend(
                [
                    "--output-html",
                    str(project / "reports" / "review_multi_comparisons.html"),
                ]
            )
        else:
            review.extend(
                [
                    "--output-dir",
                    str(project / "reports" / "review_multi_comparisons"),
                ]
            )
        commands.append(("review-multi-comparisons", review))
    return commands


def build_prepare_command(args: argparse.Namespace) -> list[str]:
    """Build the raw-read preparation stage used by the end-to-end wrapper."""

    command = [
        "prepare-atac",
        "--samples",
        str(args.reads_table),
        "--genome",
        str(args.genome),
        "--outdir",
        str(Path(args.outdir).expanduser().resolve()),
        "--cores",
        str(args.cores),
        "--runtime",
        "system",
    ]
    scalar_options = {
        "config": "--config",
        "profile": "--profile",
        "reference_dir": "--reference-dir",
        "fasta": "--fasta",
        "bowtie2_index": "--bowtie2-index",
        "blacklist": "--blacklist",
        "tss": "--tss",
        "macs_genome_size": "--macs-genome-size",
        "max_parallel_samples": "--max-parallel-samples",
        "memory_gb": "--memory-gb",
    }
    for attribute, flag in scalar_options.items():
        value = getattr(args, attribute, None)
        if value is not None:
            command.extend([flag, str(value)])
    if args.keep_intermediates:
        command.append("--keep-intermediates")
    if args.force:
        command.append("--no-resume")
    if args.fail_fast:
        command.append("--fail-fast")
    return command


def _prepared_inputs(project: Path) -> tuple[Path, str, str | None]:
    sample_table = project / "metadata" / "samples.tsv"
    reference_path = project / "metadata" / "reference.json"
    if not sample_table.is_file() or not reference_path.is_file():
        raise ValueError(
            "Raw-read preparation did not produce metadata/samples.tsv and metadata/reference.json"
        )
    reference = json.loads(reference_path.read_text(encoding="utf-8"))
    genome = str(reference.get("fasta") or "")
    if not genome:
        raise ValueError("The prepared reference metadata does not contain a genome FASTA")
    blacklist = str(reference.get("blacklist") or "") or None
    return sample_table, genome, blacklist


def _planned_prepared_inputs(args: argparse.Namespace, project: Path) -> argparse.Namespace:
    """Resolve predictable downstream paths without creating or downloading data."""

    planned = argparse.Namespace(**vars(args))
    planned.sample_table = str(project / "metadata" / "samples.tsv")
    reference_root = Path(
        args.reference_dir or Path.home() / ".cache" / "fp-tools" / "references"
    ).expanduser()
    planned.genome = str(
        Path(args.fasta).expanduser().resolve()
        if args.fasta
        else reference_root / str(args.genome) / f"{args.genome}.fa"
    )
    if args.blacklist:
        planned.blacklist = str(Path(args.blacklist).expanduser().resolve())
    elif args.genome in {"hg38", "mm10"}:
        planned.blacklist = str(
            reference_root / str(args.genome) / f"{args.genome}.blacklist.bed"
        )
    else:
        planned.blacklist = None
    return planned


def run_bulk_footprinting(args: argparse.Namespace) -> int:
    project = Path(args.outdir).expanduser().resolve()
    project.mkdir(parents=True, exist_ok=True)
    log_dir = project / "logs" / "bulk_footprinting"
    log_dir.mkdir(parents=True, exist_ok=True)
    comparisons = read_comparison_table(args.comparison_table)
    prepare_command = None
    if args.reads_table:
        prepare_command = build_prepare_command(args)
        if args.dry_run:
            print(f"[prepare-atac] {_quote([*prepare_command, '--dry-run'])}")
            for label, command in build_commands(_planned_prepared_inputs(args, project)):
                print(f"[{label}] {_quote(command)}")
            return 0
        if not (args.resume and not args.force and (project / "metadata" / "samples.tsv").is_file()):
            code = _run(
                prepare_command,
                log_dir / "prepare-atac.stdout.log",
                log_dir / "prepare-atac.stderr.log",
            )
            if code:
                print(f"prepare-atac failed with exit code {code}; see {log_dir}", file=sys.stderr)
                return code
        else:
            print("[resume] prepare-atac: complete")
        sample_table, genome, blacklist = _prepared_inputs(project)
        args.sample_table = str(sample_table)
        args.genome = genome
        args.blacklist = blacklist

    samples = read_sample_table(args.sample_table)
    conditions = {row.condition for row in samples}
    unknown = sorted({value for row in comparisons for value in (row.cond1, row.cond2)} - conditions)
    if unknown:
        raise ValueError(f"Comparison table references unknown conditions: {', '.join(unknown)}")

    review_format = _resolve_review_format(args)
    commands = build_commands(args)
    command_file = log_dir / "bulk_footprinting_commands.sh"
    lines = ["#!/usr/bin/env bash", "set -euo pipefail", "", "# Generated by bulk-footprinting.", ""]
    logged_commands = ([('prepare-atac', prepare_command)] if prepare_command else []) + commands
    for label, command in logged_commands:
        lines.extend([f"# {label}", _quote(command), ""])
    command_file.write_text("\n".join(lines), encoding="utf-8")
    command_file.chmod(0o755)

    exit_code = 0
    for label, command in commands:
        if args.resume and not args.force and _stage_complete(
            label,
            project,
            samples,
            comparisons,
            review_format=review_format,
        ):
            print(f"[resume] {label}: complete")
            continue
        if args.dry_run:
            print(f"[{label}] {_quote(command)}")
            continue
        code = _run(command, log_dir / f"{label}.stdout.log", log_dir / f"{label}.stderr.log")
        if code:
            print(f"{label} failed with exit code {code}; see {log_dir}", file=sys.stderr)
            exit_code = code
            if args.fail_fast:
                break
            break
    print(f"Wrote command log to {command_file}")
    return exit_code


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="bulk-footprinting", description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--sample-table", help="TSV with sample, condition, BAM, and peak BED columns.")
    source.add_argument(
        "--reads-table",
        help=(
            "Linux CLI/container only: TSV/CSV with local FASTQ paths or public run "
            "accessions; runs preparation before footprinting."
        ),
    )
    parser.add_argument("--comparison-table", required=True, help="TSV with comparison, cond1, and cond2 columns.")
    parser.add_argument("--genome", required=True, help="Reference FASTA for aligned input, or hg38/mm10/custom label for raw reads.")
    parser.add_argument("--blacklist", help="Optional blacklist BED used during bias correction.")
    parser.add_argument("--config", help="Optional prepare-atac YAML for raw reads.")
    parser.add_argument("--profile", choices=["modern", "homer-atac"], default="modern", help="Raw-read processing profile (default: modern).")
    parser.add_argument("--reference-dir", help="Reference cache root for raw reads.")
    parser.add_argument("--fasta", help="Custom reference FASTA for raw reads.")
    parser.add_argument("--bowtie2-index", help="Existing Bowtie2 index prefix for raw reads.")
    parser.add_argument("--tss", help="Optional TSS BED for raw-read QC.")
    parser.add_argument("--macs-genome-size", help="MACS3 genome size for a custom genome.")
    parser.add_argument("--max-parallel-samples", type=int, help="Maximum raw-read samples processed concurrently.")
    parser.add_argument("--memory-gb", type=float, help="Total memory budget for raw-read processing.")
    parser.add_argument("--keep-intermediates", action="store_true", help="Keep raw-read intermediate files.")
    parser.add_argument("--motifs", nargs="*", help="Optional motif files.")
    parser.add_argument("--motif-db", default="jaspar2026_vertebrates", help="Built-in motif database (default: jaspar2026_vertebrates).")
    parser.add_argument("--normalization", choices=["none", "condition-quantile", "sample-quantile"], default="none", help="Differential-stage normalization (default: none).")
    parser.add_argument("--plot-aggregate", choices=["sig", "all", "top", "off"], default="all", help="Aggregate profiles generated by diff-footprints (default: all).")
    parser.add_argument("--review-format", choices=["auto", "bundle", "standalone", "none"], default="auto", help="Combined-review output; auto uses standalone when aggregation is off and bundle otherwise (default: auto).")
    parser.add_argument("--outdir", required=True, help="Project output directory.")
    parser.add_argument("--cores", type=int, default=1, help="Total worker cores passed to each stage (default: 1).")
    parser.add_argument("--resume", action="store_true", help="Skip stages whose expected outputs are complete.")
    parser.add_argument("--force", action="store_true", help="Rerun stages even when outputs already exist.")
    parser.add_argument("--dry-run", action="store_true", help="Validate inputs and print the commands without running them.")
    parser.add_argument("--fail-fast", action="store_true", help="Stop at the first failed stage.")
    add_runtime_argument(parser)
    return parser


def _require_input_file(value: str | None, flag: str) -> None:
    if not value:
        return
    path = Path(value).expanduser()
    if not path.is_file():
        raise ValueError(f"{flag} file does not exist: {value}")


def _preflight_input_paths(args: argparse.Namespace) -> None:
    """Reject missing workflow inputs before runtime provisioning or writes."""

    _require_input_file(args.reads_table or args.sample_table, "--reads-table" if args.reads_table else "--sample-table")
    _require_input_file(args.comparison_table, "--comparison-table")
    if args.reads_table:
        for value, flag in (
            (args.config, "--config"),
            (args.fasta, "--fasta"),
            (args.blacklist, "--blacklist"),
            (args.tss, "--tss"),
        ):
            _require_input_file(value, flag)
    else:
        _require_input_file(args.genome, "--genome")
        _require_input_file(args.blacklist, "--blacklist")
        for motif in args.motifs or []:
            _require_input_file(motif, "--motifs")


def main(argv: list[str] | None = None) -> int:
    raw_args = list(sys.argv[1:] if argv is None else argv)
    parser = build_parser()
    args = parser.parse_args(raw_args)
    if args.resume and args.force:
        parser.error("--resume and --force are mutually exclusive")
    if args.reads_table:
        try:
            require_raw_read_preparation_support()
        except ValueError as exc:
            parser.error(str(exc))
    elif any(
        argument.split("=", 1)[0]
        in {
            "--config",
            "--profile",
            "--reference-dir",
            "--fasta",
            "--bowtie2-index",
            "--tss",
            "--macs-genome-size",
            "--max-parallel-samples",
            "--memory-gb",
            "--keep-intermediates",
        }
        for argument in raw_args
    ):
        parser.error("raw-read preparation options require --reads-table")
    try:
        _preflight_input_paths(args)
        if args.reads_table and not args.dry_run:
            delegated = prepare_command_runtime(
                "bulk-footprinting",
                raw_args,
                args.runtime,
                "core",
                {
                    "--sample-table",
                    "--reads-table",
                    "--comparison-table",
                    "--genome",
                    "--blacklist",
                    "--config",
                    "--reference-dir",
                    "--fasta",
                    "--bowtie2-index",
                    "--tss",
                    "--outdir",
                    "--motifs",
                },
            )
            if delegated is not None:
                return delegated
            if args.profile == "homer-atac":
                prepare_command_runtime("bulk-footprinting", raw_args, args.runtime, "homer", set())
        return run_bulk_footprinting(args)
    except (ValueError, RuntimeProvisionError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    raise SystemExit(main())
