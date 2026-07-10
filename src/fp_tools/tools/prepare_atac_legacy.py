"""ATAC-only compatibility runner for the historical fp-tools preprocessing route."""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path
from typing import Any

import pysam


def _concatenate(inputs: list[Path], output: Path) -> Path:
    if len(inputs) == 1:
        return inputs[0]
    with output.open("wb") as target:
        for source in inputs:
            with source.open("rb") as handle:
                shutil.copyfileobj(handle, target, length=8 * 1024 * 1024)
    return output


def _remove_xs(input_bam: Path, output_bam: Path) -> int:
    kept = 0
    with (
        pysam.AlignmentFile(input_bam, "rb") as source,
        pysam.AlignmentFile(output_bam, "wb", template=source) as target,
    ):
        for read in source.fetch(until_eof=True):
            if not read.has_tag("XS"):
                target.write(read)
                kept += 1
    return kept


def _ensure_bigwig(output: Path, chrom_sizes: Path) -> None:
    """Convert HOMER's fallback bedGraph when its UCSC converter is unavailable."""
    from fp_tools.tools.prepare_atac import _bedgraph_to_bigwig

    with output.open("rb") as handle:
        magic = handle.read(4)
    if magic in {b"&\xfc\x8f\x88", b"\x88\x8f\xfc&"}:
        return
    bedgraph = output.with_suffix(".homer.bedgraph")
    sorted_bedgraph = output.with_suffix(".sorted.bedgraph")
    temporary_bigwig = output.with_suffix(".tmp.bw")
    output.replace(bedgraph)
    try:
        with sorted_bedgraph.open("w", encoding="utf-8") as handle:
            subprocess.run(
                ["bedtools", "sort", "-faidx", str(chrom_sizes), "-i", str(bedgraph)],
                check=True,
                stdout=handle,
            )
        _bedgraph_to_bigwig(sorted_bedgraph, chrom_sizes, temporary_bigwig)
        temporary_bigwig.replace(output)
        bedgraph.unlink(missing_ok=True)
    except Exception:
        temporary_bigwig.unlink(missing_ok=True)
        if not output.exists() and bedgraph.exists():
            bedgraph.replace(output)
        raise
    finally:
        sorted_bedgraph.unlink(missing_ok=True)


def process_legacy_sample(
    sample, root: Path, reference, settings: dict[str, Any], resume: bool = True
) -> dict[str, str]:
    """Run the paired-end legacy ATAC route while retaining modern project metadata."""
    from fp_tools.tools.prepare_atac import (
        MITO_CHROMS,
        _bam_count,
        _fingerprint,
        _fragment_metrics,
        _frip,
        _relative_link,
        _run,
        _run_pipeline,
        _tss_enrichment,
        materialize_run_fastqs,
    )

    sample_root = root / "samples" / sample.sample
    alignment = sample_root / "alignment"
    peaks_dir = sample_root / "peaks"
    tracks = sample_root / "tracks"
    qc = sample_root / "qc"
    work = sample_root / ".work"
    fastq_dir = root / "fastq"
    for directory in (alignment, peaks_dir, tracks, qc, work):
        directory.mkdir(parents=True, exist_ok=True)

    fastqs = [materialize_run_fastqs(run, fastq_dir, settings) for run in sample.runs]
    if not fastqs or any(r2 is None for _, r2 in fastqs):
        raise ValueError(
            f"legacy-atac requires paired-end FASTQs for sample {sample.sample}"
        )

    final_bam = alignment / f"{sample.sample}.filtered.bam"
    final_bai = Path(str(final_bam) + ".bai")
    peak_bed = peaks_dir / f"{sample.sample}.narrowPeak"
    homer_txt = peaks_dir / f"{sample.sample}.rp10m.narrowpeaks.txt"
    track_bw = tracks / f"{sample.sample}.rp10m.bw"
    state_path = sample_root / "state.json"
    fingerprint = _fingerprint(sample, reference, settings, fastqs)
    if (
        resume
        and state_path.exists()
        and all(
            path.exists()
            for path in (final_bam, final_bai, peak_bed, homer_txt, track_bw)
        )
    ):
        state = json.loads(state_path.read_text(encoding="utf-8"))
        if (
            state.get("fingerprint") == fingerprint
            and state.get("status") == "complete"
        ):
            return {
                "sample": sample.sample,
                "condition": sample.condition,
                "bam": str(final_bam),
                "peaks": str(peak_bed),
                "bigwig": str(track_bw),
                "status": "cached",
            }

    cores = max(1, int(settings["resources"]["cores"]))
    threads = str(cores)
    log = qc / "commands.log"
    combined_r1 = _concatenate(
        [r1 for r1, _ in fastqs], work / f"{sample.sample}.R1.fastq.gz"
    )
    combined_r2 = _concatenate(
        [r2 for _, r2 in fastqs if r2], work / f"{sample.sample}.R2.fastq.gz"
    )

    if settings["trim"].get("enabled"):
        trim_command = [
            "trim_galore",
            "--paired",
            "--fastqc",
            "--cores",
            threads,
            "-o",
            str(work),
            str(combined_r1),
            str(combined_r2),
            "--basename",
            sample.sample,
        ]
        trim_command.extend(
            str(value) for value in settings["trim"].get("extra_args", [])
        )
        _run(trim_command, log)
        trimmed_r1 = work / f"{sample.sample}_val_1.fq.gz"
        trimmed_r2 = work / f"{sample.sample}_val_2.fq.gz"
    else:
        trimmed_r1, trimmed_r2 = combined_r1, combined_r2

    initial = work / f"{sample.sample}.{reference.assembly}_alignment.bam"
    bowtie = [
        "bowtie2",
        "-p",
        threads,
        "-x",
        str(reference.index_prefix),
        "-1",
        str(trimmed_r1),
        "-2",
        str(trimmed_r2),
        *[str(value) for value in settings["align"].get("extra_args", [])],
        "-X",
        str(settings["align"]["max_insert"]),
        "--interleaved",
        "-",
    ]
    _run_pipeline(
        bowtie, ["samtools", "view", "-@", threads, "-bh", "-o", str(initial), "-"], log
    )
    with (qc / "mapped.txt").open("w", encoding="utf-8") as output:
        _run(
            ["samtools", "view", "-@", threads, "-c", "-F", "12", str(initial)],
            log,
            stdout=output,
        )

    sorted_bam = work / "coordinate_sorted.bam"
    rg_bam = work / "read_groups.bam"
    dupmark = work / "dupmark.bam"
    dupmetrics = qc / "dupmetrics.txt"
    _run(
        [
            "samtools",
            "sort",
            "-@",
            threads,
            "-m",
            "384M",
            "-o",
            str(sorted_bam),
            str(initial),
        ],
        log,
    )
    _run(
        [
            "samtools",
            "addreplacerg",
            "-r",
            "@RG\tID:RG1\tSM:SampleName\tPL:Illumina\tLB:Library.fa",
            "-o",
            str(rg_bam),
            str(sorted_bam),
        ],
        log,
    )
    _run(
        [
            "picard",
            "-Xmx8g",
            "MarkDuplicates",
            f"I={rg_bam}",
            f"O={dupmark}",
            f"M={dupmetrics}",
            f"TMP_DIR={work}",
        ],
        log,
    )

    nodup = work / "nodup.bam"
    unique = work / "unique.bam"
    unique_sorted = work / "unique.sorted.bam"
    _run(
        [
            "samtools",
            "view",
            "-@",
            threads,
            "-b",
            "-F",
            "1024",
            "-o",
            str(nodup),
            str(dupmark),
        ],
        log,
    )
    unique_count = _remove_xs(nodup, unique)
    (qc / "no2nd.txt").write_text(f"{unique_count}\n", encoding="utf-8")
    _run(
        [
            "samtools",
            "sort",
            "-@",
            threads,
            "-m",
            "384M",
            "-o",
            str(unique_sorted),
            str(unique),
        ],
        log,
    )
    if reference.blacklist and reference.blacklist.exists():
        blacklist_unsorted = work / "blacklist_filtered.bam"
        with blacklist_unsorted.open("wb") as output:
            _run(
                [
                    "bedtools",
                    "intersect",
                    "-v",
                    "-a",
                    str(unique_sorted),
                    "-b",
                    str(reference.blacklist),
                    "-ubam",
                ],
                log,
                stdout=output,
            )
        _run(
            [
                "samtools",
                "sort",
                "-@",
                threads,
                "-m",
                "384M",
                "-o",
                str(final_bam),
                str(blacklist_unsorted),
            ],
            log,
        )
    else:
        shutil.copy2(unique_sorted, final_bam)
    _run(["samtools", "index", "-@", threads, str(final_bam)], log)
    (qc / "filtered.txt").write_text(f"{_bam_count(final_bam)}\n", encoding="utf-8")
    with (qc / "flagstat.tsv").open("w", encoding="utf-8") as output:
        _run(
            ["samtools", "flagstat", "-@", threads, "-O", "tsv", str(final_bam)],
            log,
            stdout=output,
        )

    homer_bam = work / "homer.no_chrEBV.bam"
    with (
        pysam.AlignmentFile(final_bam, "rb") as source,
        pysam.AlignmentFile(homer_bam, "wb", template=source) as target,
    ):
        for read in source.fetch(until_eof=True):
            if read.reference_name != "chrEBV":
                target.write(read)
    tag_directory = work / "TagDirectory"
    _run(["makeTagDirectory", str(tag_directory), "-sspe", str(homer_bam)], log)
    _run(
        [
            "makeUCSCfile",
            str(tag_directory),
            "-bigWig",
            str(reference.chrom_sizes),
            "-o",
            str(track_bw),
        ],
        log,
    )
    _ensure_bigwig(track_bw, reference.chrom_sizes)
    _run(
        [
            "findPeaks",
            str(tag_directory),
            "-style",
            str(settings["peaks"]["homer_style"]),
            "-L",
            str(settings["peaks"]["homer_local_fold"]),
            "-localSize",
            str(settings["peaks"]["homer_local_size"]),
            "-fdr",
            str(settings["peaks"]["qvalue"]),
            "-o",
            str(homer_txt),
            *[str(value) for value in settings["peaks"].get("extra_args", [])],
        ],
        log,
    )
    with peak_bed.open("w", encoding="utf-8") as output:
        _run(["pos2bed.pl", str(homer_txt)], log, stdout=output)
    filtered_peak = peaks_dir / f"{sample.sample}.narrowPeak.filtered.bed"
    with (
        peak_bed.open(encoding="utf-8") as source,
        filtered_peak.open("w", encoding="utf-8") as target,
    ):
        for line in source:
            fields = line.rstrip().split("\t")
            if (
                len(fields) >= 3
                and fields[1].isdigit()
                and fields[2].isdigit()
                and fields[0] not in MITO_CHROMS
            ):
                target.write(line)

    metrics = {
        "sample": sample.sample,
        "condition": sample.condition,
        "profile": "legacy-atac",
        "paired_end": True,
        "usable_reads": max(1, _bam_count(final_bam) // 2),
        "peaks": sum(
            1
            for line in peak_bed.open(encoding="utf-8")
            if len(line.rstrip().split("\t")) >= 3
            and line.rstrip().split("\t")[1].isdigit()
            and line.rstrip().split("\t")[2].isdigit()
        ),
        "frip": round(_frip(final_bam, peak_bed, log), 6),
        "tss_enrichment": _tss_enrichment(final_bam, reference.tss),
        "mitochondrial_fraction": None,
    }
    idxstats = pysam.idxstats(str(final_bam)).splitlines()
    mapped_by_chrom = {
        fields[0]: int(fields[2])
        for line in idxstats
        if len(fields := line.split("\t")) >= 3
    }
    total_mapped = sum(mapped_by_chrom.values())
    metrics["mitochondrial_fraction"] = (
        round(
            sum(mapped_by_chrom.get(chrom, 0) for chrom in MITO_CHROMS) / total_mapped,
            6,
        )
        if total_mapped
        else None
    )
    metrics.update(_fragment_metrics(final_bam, qc / "fragment_lengths.tsv"))
    metrics["warnings"] = []
    (qc / "metrics.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")

    legacy = root / sample.sample
    suffix = reference.assembly
    _relative_link(final_bam, legacy / f"{sample.sample}.{suffix}.filtered.bam")
    _relative_link(final_bai, legacy / f"{sample.sample}.{suffix}.filtered.bam.bai")
    _relative_link(track_bw, legacy / f"{sample.sample}.{suffix}.rp10m.bw")
    _relative_link(
        homer_txt, legacy / f"{sample.sample}.{suffix}.rp10m.narrowpeaks.txt"
    )
    _relative_link(peak_bed, legacy / f"{sample.sample}.{suffix}.rp10m.narrowpeaks.bed")
    _relative_link(
        filtered_peak,
        legacy / f"{sample.sample}.{suffix}.rp10m.narrowpeaks.filtered.bed",
    )
    state_path.write_text(
        json.dumps(
            {"status": "complete", "fingerprint": fingerprint, "metrics": metrics},
            indent=2,
        ),
        encoding="utf-8",
    )

    if not settings.get("cleanup", {}).get("keep_intermediates"):
        shutil.rmtree(work, ignore_errors=True)
    return {
        "sample": sample.sample,
        "condition": sample.condition,
        "bam": str(final_bam),
        "peaks": str(peak_bed),
        "bigwig": str(track_bw),
        "status": "complete",
    }
