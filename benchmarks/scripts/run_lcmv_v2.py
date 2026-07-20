#!/usr/bin/env python3
"""Run the restartable LCMV CD8 ATAC/RNA v2 expansion on a local server."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SELECTION = ROOT / "benchmarks/manifests/compact/lcmv_cd8_libraries_v2.tsv"
COMPARISONS = ROOT / "benchmarks/manifests/compact/lcmv_cd8_comparisons_v2.tsv"
RAW = ROOT / "data/public/raw/lcmv_cd8_bulk"
META = RAW / "metadata/v2"
ADDITIONS = RAW / "atac_v2_additions"
V1 = ROOT / "data/public/processed/lcmv_cd8_bulk_fp_rna"
PROJECT = ROOT / "data/public/processed/lcmv_cd8_bulk_fp_rna_v2"
GENOME = ROOT / "data/public/raw/gse192390_mm10_reference/genome/mm10_no_alt_analysis_set_ENCODE.fasta"
BLACKLIST = ROOT / "data/public/raw/gse192390_mm10_reference/blacklist/ENCFF547MET_mm10.bed"
INDEX = RAW / "reference/mm10/bowtie2/mm10"
TSS = RAW / "reference/mm10/mm10.tss.bed"
GTF = RAW / "reference/mm10/gencode_m25/gencode.vM25.annotation.gtf"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in fields} for row in rows)


def fingerprint(paths: list[Path], values: list[str]) -> str:
    digest = hashlib.sha256()
    for path in paths:
        digest.update(str(path.resolve()).encode())
        digest.update(path.read_bytes())
    for value in values:
        digest.update(value.encode())
    return digest.hexdigest()


class Runner:
    def __init__(self, args: argparse.Namespace):
        self.args = args
        self.state = PROJECT / "state"
        self.state.mkdir(parents=True, exist_ok=True)
        self.env = os.environ.copy()
        self.env["PATH"] = f"{args.atac_tool_bin}:{self.env.get('PATH', '')}"
        self.python = Path(sys.executable)
        self.fp_bin = ROOT / ".venv/bin"

    def command(self, command: list[object], *, env: dict[str, str] | None = None) -> None:
        argv = [str(value) for value in command]
        print("+", " ".join(argv), flush=True)
        if not self.args.dry_run:
            subprocess.run(argv, cwd=ROOT, env=env or self.env, check=True)

    def run_stage(self, name: str, inputs: list[Path], action) -> None:
        token = fingerprint(inputs, [name, str(self.args.cores), str(self.args.memory_gb)])
        state_path = self.state / f"{name}.json"
        if self.args.resume and state_path.is_file():
            prior = json.loads(state_path.read_text())
            if prior.get("fingerprint") == token and prior.get("status") == "complete":
                print(f"Skipping verified stage: {name}")
                return
        action()
        if not self.args.dry_run:
            state_path.write_text(
                json.dumps({"stage": name, "fingerprint": token, "status": "complete"}, indent=2) + "\n",
                encoding="utf-8",
            )

    def guard_space(self) -> None:
        free_gib = shutil.disk_usage(ROOT).free / 2**30
        if free_gib < self.args.min_free_gb:
            raise RuntimeError(
                f"Only {free_gib:.1f} GiB free; v2 requires at least {self.args.min_free_gb:.1f} GiB"
            )

    def resolve(self) -> None:
        self.command(
            [
                self.python,
                ROOT / "benchmarks/scripts/resolve_lcmv_cd8_collection.py",
                "--selection", SELECTION,
                "--root", RAW,
                "--metadata-dir", META,
            ]
        )

    def download(self) -> None:
        self.guard_space()
        self.command(
            [
                self.python,
                ROOT / "benchmarks/scripts/download_manifest.py",
                "--manifest", META / "download_manifest.tsv",
                "--report", META / "download_report.tsv",
                "--downloader", "auto",
            ]
        )

    def prepare_atac(self) -> None:
        self.guard_space()
        resolved = read_tsv(META / "resolved_runs.tsv")
        existing = {
            row["sample"]
            for row in read_tsv(RAW / "metadata/all_atac_samples.tsv")
        }
        new = [
            {
                **row,
                "sample": row["run_accession"],
                "condition": row["harmonized_condition"],
            }
            for row in resolved
            if row["assay"] == "ATAC" and row["run_accession"] not in existing
        ]
        fields = [
            "run_accession", "sample", "condition", "replicate", "fastq_1", "fastq_2",
            "fastq_1_md5", "fastq_2_md5", "author", "year", "series", "subseries",
            "gsm_accession", "assay", "collection", "infection", "day", "tissue",
            "original_condition", "harmonized_condition", "broader_state", "technical_partition",
            "condition_pair_id", "rna_match_status", "include_in_primary_paired_analysis",
        ]
        sample_table = META / "new_atac_samples.tsv"
        if not self.args.dry_run:
            write_tsv(sample_table, new, fields)
        self.command(
            [
                self.fp_bin / "prepare-atac",
                "--samples", sample_table,
                "--outdir", ADDITIONS,
                "--profile", "legacy-atac",
                "--config", RAW / "atac_primary/metadata/resolved_settings.yml",
                "--fasta", GENOME,
                "--bowtie2-index", INDEX,
                "--blacklist", BLACKLIST,
                "--tss", TSS,
                "--cores", self.args.cores,
                "--memory-gb", self.args.memory_gb,
                "--fail-fast",
            ]
        )

    def build_metadata(self) -> None:
        self.command(
            [
                self.python,
                ROOT / "benchmarks/scripts/build_lcmv_cd8_downstream.py",
                "--resolved-runs", META / "resolved_runs.tsv",
                "--project", PROJECT,
                "--comparisons", COMPARISONS,
                "--atac-outputs", RAW / "metadata/all_atac_samples.tsv",
                "--new-atac-root", ADDITIONS,
            ]
        )

    def merge_technical(self) -> None:
        samtools = Path(self.args.atac_tool_bin) / "samtools"
        bedtools = Path(self.args.atac_tool_bin) / "bedtools"
        inventory = read_tsv(PROJECT / "metadata/multimodal_inventory.tsv")
        for row in inventory:
            srrs = row["ATAC_SRR"].split(",")
            if len(srrs) == 1:
                continue
            out_bam = Path(row["ATAC_BAM"])
            out_bed = Path(row["ATAC_peaks"])
            out_bam.parent.mkdir(parents=True, exist_ok=True)
            resolved = read_tsv(META / "resolved_runs.tsv")
            by_run = {item["run_accession"]: item for item in resolved}
            existing_outputs = {item["sample"]: item for item in read_tsv(RAW / "metadata/all_atac_samples.tsv")}
            bams, beds = [], []
            for run in srrs:
                output = existing_outputs.get(run)
                if output is None:
                    base = ADDITIONS / "samples" / run
                    output = {
                        "bam": str(base / "alignment" / f"{run}.filtered.bam"),
                        "peaks": str(base / "peaks" / f"{run}.narrowPeak"),
                    }
                if run not in by_run:
                    raise ValueError(f"Unknown technical partition {run}")
                bams.append(output["bam"])
                beds.append(output["peaks"])
            self.command([samtools, "merge", "-f", "-@", self.args.cores, out_bam, *bams])
            self.command([samtools, "index", "-@", self.args.cores, out_bam])
            if not self.args.dry_run:
                with (out_bed.parent / "technical_peaks.bed").open("wb") as target:
                    for path in beds:
                        target.write(Path(path).read_bytes())
            if not self.args.dry_run:
                sorted_bed = out_bed.parent / "technical_peaks.sorted.bed"
                subprocess.run(
                    [str(bedtools), "sort", "-i", str(out_bed.parent / "technical_peaks.bed")],
                    stdout=sorted_bed.open("wb"), check=True, env=self.env,
                )
                subprocess.run(
                    [str(bedtools), "merge", "-i", str(sorted_bed)],
                    stdout=out_bed.open("wb"), check=True, env=self.env,
                )

    def atac_primary(self) -> None:
        samples = PROJECT / "metadata/samples.tsv"
        peaks = PROJECT / "peaks/merged_peaks_filtered.bed"
        self.command([self.fp_bin / "atac-correct", "--sample-table", samples, "--outdir", PROJECT, "--genome", GENOME, "--blacklist", BLACKLIST, "--cores", self.args.cores])
        self.command([self.fp_bin / "normalize-bigwig", "--sample-table", samples, "--outdir", PROJECT, "--background", peaks, "--method", "background-scale", "--stat", "q95", "--target", "median", "--workers", self.args.cores])
        self.command([self.fp_bin / "call-footprints", "--sample-table", samples, "--outdir", PROJECT, "--regions", peaks, "--cores", self.args.cores])
        self.command([self.fp_bin / "match-motifs", "--sample-table", samples, "--outdir", PROJECT, "--genome", GENOME, "--peaks", peaks, "--motif-db", "jaspar2026_vertebrates", "--cores", self.args.cores])
        self.command([self.fp_bin / "diff-footprints", "--sample-table", samples, "--comparison-table", PROJECT / "metadata/comparisons.tsv", "--outdir", PROJECT, "--genome", GENOME, "--peaks", peaks, "--motif-db", "jaspar2026_vertebrates", "--cores", self.args.cores])
        self.command([self.fp_bin / "review-multi-comparisons", "--outdir", PROJECT, "--title", "LCMV CD8 v2 paired within-study footprint comparisons"])

    def atac_supporting(self) -> None:
        support = PROJECT / "supporting_assay_only"
        samples = PROJECT / "metadata/supporting_samples.tsv"
        combined = PROJECT / "metadata/supporting_analysis_samples.tsv"
        peaks = PROJECT / "peaks/merged_peaks_filtered.bed"
        self.command([self.fp_bin / "atac-correct", "--sample-table", samples, "--outdir", support, "--genome", GENOME, "--blacklist", BLACKLIST, "--cores", self.args.cores])
        self.command([self.fp_bin / "normalize-bigwig", "--sample-table", samples, "--outdir", support, "--background", peaks, "--method", "background-scale", "--stat", "q95", "--target", "median", "--workers", self.args.cores])
        self.command([self.fp_bin / "call-footprints", "--sample-table", samples, "--outdir", support, "--regions", peaks, "--cores", self.args.cores])
        self.command([self.fp_bin / "match-motifs", "--sample-table", samples, "--outdir", support, "--genome", GENOME, "--peaks", peaks, "--motif-db", "jaspar2026_vertebrates", "--cores", self.args.cores])
        if not self.args.dry_run:
            (support / "samples").mkdir(parents=True, exist_ok=True)
            for row in read_tsv(PROJECT / "metadata/samples.tsv"):
                link = support / "samples" / row["sample"]
                if not link.exists():
                    link.symlink_to(PROJECT / "samples" / row["sample"], target_is_directory=True)
            if not (support / "peaks").exists():
                (support / "peaks").symlink_to(PROJECT / "peaks", target_is_directory=True)
        self.command([self.fp_bin / "diff-footprints", "--sample-table", combined, "--comparison-table", PROJECT / "metadata/supporting_atac_comparisons.tsv", "--outdir", support, "--genome", GENOME, "--peaks", peaks, "--motif-db", "jaspar2026_vertebrates", "--cores", self.args.cores])
        self.command([self.fp_bin / "review-multi-comparisons", "--outdir", support, "--title", "LCMV CD8 v2 supporting ATAC-only comparisons"])

    def initialize_rna(self) -> list[dict[str, str]]:
        rows = read_tsv(PROJECT / "rna/metadata/samples.tsv")
        old_samples = {row["sample"] for row in read_tsv(V1 / "rna/metadata/samples.tsv")}
        for layer in ("uniform_kallisto", "paper_specific"):
            root = PROJECT / "rna" / layer
            root.mkdir(parents=True, exist_ok=True)
            for sample in sorted(old_samples):
                link = root / sample
                if not link.exists():
                    link.symlink_to(V1 / "rna" / layer / sample, target_is_directory=True)
        return [row for row in rows if row["sample"] not in old_samples]

    def stage_reads(self, row: dict[str, str]) -> tuple[Path, Path | None]:
        sample = row["sample"]
        root = PROJECT / "rna/staged_fastq" / sample
        root.mkdir(parents=True, exist_ok=True)
        read1 = root / "read1.fastq.gz"
        read2 = root / "read2.fastq.gz" if row["library_layout"] == "PAIRED" else None
        sources1 = [Path(value) for value in row["fastq_1"].split(";") if value]
        sources2 = [Path(value) for value in row["fastq_2"].split(";") if value]
        for output, sources in ((read1, sources1), (read2, sources2)):
            if output is None or output.exists():
                continue
            if len(sources) == 1:
                output.symlink_to(sources[0])
            else:
                with output.open("wb") as target:
                    for source in sources:
                        with source.open("rb") as handle:
                            shutil.copyfileobj(handle, target, 8 * 1024 * 1024)
        return read1, read2

    def rna(self) -> None:
        current = Path(self.args.rna_tool_bin)
        legacy = Path(self.args.legacy_rna_tool_bin)
        new = self.initialize_rna()
        k21 = RAW / "reference/mm10/gencode_m25/kallisto_0.51_k21/transcripts.idx"
        for row in new:
            sample = row["sample"]
            read1, read2 = self.stage_reads(row)
            uniform = PROJECT / "rna/uniform_kallisto" / sample
            uniform.mkdir(parents=True, exist_ok=True)
            command = [current / "kallisto", "quant", "-i", k21, "-o", uniform, "-t", self.args.cores, "-b", "0"]
            if read2 is None:
                command += ["--single", "-l", "200", "-s", "20", read1]
            else:
                command += [read1, read2]
            if not (uniform / "run_info.json").is_file():
                self.command(command)
            paper = PROJECT / "rna/paper_specific" / sample
            paper.mkdir(parents=True, exist_ok=True)
            tophat = paper / "tophat"
            command = [legacy / "tophat2", "-p", self.args.cores, "-G", GTF, "-o", tophat, RAW / "reference/mm10/bowtie2/mm10", read1]
            if read2 is not None:
                command.append(read2)
            if not (tophat / "accepted_hits.bam").is_file():
                self.command(command)
            bam = tophat / "accepted_hits.bam"
            order = "pos"
            if read2 is not None and not (tophat / "accepted_hits.name_sorted.bam").is_file():
                name_bam = tophat / "accepted_hits.name_sorted.bam"
                self.command([legacy / "samtools", "sort", "-n", "-@", self.args.cores, "-m", "1G", "-o", name_bam, bam])
                bam, order = name_bam, "name"
            elif read2 is not None:
                bam, order = tophat / "accepted_hits.name_sorted.bam", "name"
            counts = paper / "htseq_counts.tsv"
            print("+", legacy / "htseq-count", "...", counts)
            if not self.args.dry_run and not (paper / "htseq_counts.done").is_file():
                with counts.open("wb") as output:
                    subprocess.run([legacy / "htseq-count", "-f", "bam", "-r", order, "-s", "no", "-m", "union", bam, GTF], stdout=output, check=True, env=self.env)
                (paper / "htseq_counts.done").write_text("complete\n", encoding="utf-8")
        self.command([self.python, ROOT / "benchmarks/scripts/summarize_lcmv_rna.py", "--project", PROJECT, "--gtf", GTF])
        self.command([current / "Rscript", ROOT / "benchmarks/scripts/analyze_lcmv_rna.R", PROJECT])

    def validate(self) -> None:
        command: list[object] = [
            self.python, ROOT / "benchmarks/scripts/validate_lcmv_outputs.py",
            "--raw-metadata", META, "--project", PROJECT,
            "--selection", SELECTION, "--comparisons", COMPARISONS,
            "--atac-qc", RAW / "atac_primary/reports/qc_summary.tsv",
            "--atac-qc", RAW / "atac_supplemental/reports/qc_summary.tsv",
            "--atac-qc", ADDITIONS / "reports/qc_summary.tsv",
            "--verify-checksums",
        ]
        self.command(command)

    def summarize(self) -> None:
        self.command(
            [
                self.python,
                ROOT / "benchmarks/scripts/summarize_lcmv_v2_outputs.py",
                "--raw", RAW,
                "--project", PROJECT,
                "--checksums",
            ]
        )


STAGES = (
    "resolve", "download", "prepare_atac", "build_metadata", "merge_technical",
    "atac_primary", "atac_supporting", "rna", "validate", "summarize",
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", action="append", choices=STAGES, help="Run only this stage; repeat as needed.")
    parser.add_argument("--cores", type=int, default=max(1, os.cpu_count() or 1))
    parser.add_argument("--memory-gb", type=float, default=128.0)
    parser.add_argument("--min-free-gb", type=float, default=500.0)
    parser.add_argument("--atac-tool-bin", type=Path, default=Path("/home/exouser/miniforge3/envs/fp-tools-atac/bin"))
    parser.add_argument("--rna-tool-bin", type=Path, default=Path("/home/exouser/miniforge3/envs/lcmv-rna/bin"))
    parser.add_argument("--legacy-rna-tool-bin", type=Path, default=Path("/home/exouser/miniforge3/envs/lcmv-rna-legacy/bin"))
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--no-resume", dest="resume", action="store_false")
    parser.set_defaults(resume=True)
    args = parser.parse_args()
    runner = Runner(args)
    selected = args.stage or list(STAGES)
    stage_inputs = {
        "resolve": [SELECTION],
        "download": [META / "download_manifest.tsv"],
        "prepare_atac": [META / "resolved_runs.tsv"],
        "build_metadata": [META / "resolved_runs.tsv", COMPARISONS],
        "merge_technical": [PROJECT / "metadata/multimodal_inventory.tsv"],
        "atac_primary": [PROJECT / "metadata/samples.tsv", PROJECT / "metadata/comparisons.tsv"],
        "atac_supporting": [PROJECT / "metadata/supporting_samples.tsv", PROJECT / "metadata/supporting_atac_comparisons.tsv"],
        "rna": [PROJECT / "rna/metadata/samples.tsv", PROJECT / "metadata/rna_comparisons.tsv"],
        "validate": [SELECTION, COMPARISONS, META / "resolved_runs.tsv"],
        "summarize": [PROJECT / "validation/audit_summary.json"],
    }
    for name in selected:
        runner.run_stage(name, stage_inputs[name], getattr(runner, name))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
