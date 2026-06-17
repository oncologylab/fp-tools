"""Fixed motif-site differential backends for diff-footprints."""

from __future__ import annotations

import base64
import concurrent.futures
import os
import re
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pyBigWig
import pysam
import scipy.stats

from fp_tools.tools.bindetect_functions import build_bindetect_aggregate_payload, plot_interactive_bindetect
from fp_tools.utils.logger import FpToolsLogger
from fp_tools.utils.motifs import MotifList
from fp_tools.utils.utilities import check_cores, make_directory


_COUNT_INTERVAL_INDEX = None
_COUNT_CHROMS = None
_COUNT_MOTIFS = None


@dataclass(frozen=True)
class SiteRecord:
    motif_prefix: str
    chrom: str
    start: int
    end: int
    name: str
    score: str
    strand: str

    @property
    def center(self) -> int:
        return int((self.start + self.end) // 2)


def benjamini_hochberg(pvalues) -> np.ndarray:
    pvals = np.asarray(pvalues, dtype=float)
    qvals = np.full(pvals.shape, np.nan, dtype=float)
    finite = np.isfinite(pvals)
    if not finite.any():
        return qvals
    clipped = np.clip(pvals[finite], 0.0, 1.0)
    order = np.argsort(clipped)
    ranked = clipped[order]
    adjusted = ranked * float(len(ranked)) / np.arange(1, len(ranked) + 1, dtype=float)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    unsorted = np.empty_like(adjusted)
    unsorted[order] = np.clip(adjusted, 0.0, 1.0)
    qvals[finite] = unsorted
    return qvals


def _as_float(value, default=np.nan) -> float:
    try:
        out = float(value)
    except Exception:
        return default
    return out if np.isfinite(out) else default


def _safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_") or "motif"


def _read_result_metadata(reference_dirs: list[str]) -> dict[str, dict[str, str]]:
    metadata: dict[str, dict[str, str]] = {}
    for ref_dir in reference_dirs:
        for table in Path(ref_dir).glob("*_results.txt"):
            try:
                df = pd.read_csv(table, sep="\t", dtype=str)
            except Exception:
                continue
            if "output_prefix" not in df.columns:
                continue
            for _, row in df.iterrows():
                prefix = str(row["output_prefix"])
                metadata.setdefault(prefix, {})
                for col in ("name", "motif_id", "cluster"):
                    if col in row and pd.notna(row[col]):
                        metadata[prefix][col] = str(row[col])
    return metadata


def load_motif_site_reference(reference_dirs: list[str], site_set: str = "bound-union", window: int = 100) -> tuple[dict[str, list[SiteRecord]], dict[str, dict[str, str]]]:
    if not reference_dirs:
        raise ValueError("--site-reference-dirs is required for fixed-site differential methods")
    if window <= 0:
        raise ValueError("--site-window must be > 0")

    metadata = _read_result_metadata(reference_dirs)
    by_motif: dict[str, dict[tuple[str, int, str], SiteRecord]] = {}
    for ref_dir in reference_dirs:
        patterns = ["*/beds/*_all.bed"] if site_set == "all" else ["*/beds/*_bound.bed"]
        for pattern in patterns:
            for bed in Path(ref_dir).glob(pattern):
                motif_prefix = bed.parent.parent.name
                with bed.open("r", encoding="utf-8") as handle:
                    for line in handle:
                        if not line.strip() or line.startswith(("#", "track", "browser")):
                            continue
                        fields = line.rstrip("\n").split("\t")
                        if len(fields) < 3:
                            continue
                        try:
                            start = int(fields[1])
                            end = int(fields[2])
                        except ValueError:
                            continue
                        chrom = fields[0]
                        name = fields[3] if len(fields) > 3 and fields[3] else motif_prefix
                        score = fields[4] if len(fields) > 4 and fields[4] else "0"
                        strand = fields[5] if len(fields) > 5 and fields[5] in {"+", "-"} else "."
                        center = int((start + end) // 2)
                        half = window // 2
                        win_start = max(0, center - half)
                        win_end = win_start + window
                        record = SiteRecord(motif_prefix, chrom, win_start, win_end, name, score, strand)
                        key = (chrom, center, strand)
                        by_motif.setdefault(motif_prefix, {})[key] = record

    reference = {motif: sorted(records.values(), key=lambda r: (r.chrom, r.start, r.end, r.strand)) for motif, records in by_motif.items()}
    reference = {motif: records for motif, records in reference.items() if records}
    if not reference:
        raise ValueError(f"No motif sites found in --site-reference-dirs using --reference-site-set {site_set}")
    return reference, metadata


def write_reference_beds(reference: dict[str, list[SiteRecord]], outdir: str, conditions: list[str]) -> None:
    for motif, records in reference.items():
        bed_dir = Path(outdir) / motif / "beds"
        bed_dir.mkdir(parents=True, exist_ok=True)
        lines = [
            "\t".join([r.chrom, str(r.start), str(r.end), r.name, r.score, r.strand])
            for r in records
        ]
        text = "\n".join(lines) + ("\n" if lines else "")
        all_path = bed_dir / f"{motif}_all.bed"
        all_path.write_text(text, encoding="utf-8")
        for condition in conditions:
            condition_path = bed_dir / f"{motif}_{condition}_bound.bed"
            if condition_path.exists():
                condition_path.unlink()
            try:
                os.link(all_path, condition_path)
            except OSError:
                condition_path.write_text(text, encoding="utf-8")


def _cut_positions_by_chrom(bam_path: str, chroms: set[str], read_shift: tuple[int, int], min_mapq: int) -> dict[str, np.ndarray]:
    positions = {chrom: [] for chrom in chroms}
    shift_fwd, shift_rev = read_shift
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for chrom in sorted(chroms):
            if chrom not in bam.references:
                continue
            for read in bam.fetch(chrom):
                if read.is_unmapped or read.is_duplicate or read.is_secondary or read.is_supplementary:
                    continue
                if read.mapping_quality < min_mapq:
                    continue
                pos = read.reference_end - 1 + shift_rev if read.is_reverse else read.reference_start + shift_fwd
                if pos >= 0:
                    positions[chrom].append(int(pos))
    return {chrom: np.asarray(sorted(vals), dtype=np.int64) for chrom, vals in positions.items()}


def _build_interval_index(reference: dict[str, list[SiteRecord]]) -> dict[str, dict[str, tuple[np.ndarray, np.ndarray]]]:
    interval_index: dict[str, dict[str, tuple[np.ndarray, np.ndarray]]] = {}
    for motif, records in reference.items():
        by_chrom: dict[str, list[SiteRecord]] = {}
        for record in records:
            by_chrom.setdefault(record.chrom, []).append(record)
        interval_index[motif] = {}
        for chrom, chrom_records in by_chrom.items():
            starts = np.asarray([r.start for r in chrom_records], dtype=np.int64)
            ends = np.asarray([r.end for r in chrom_records], dtype=np.int64)
            interval_index[motif][chrom] = (starts, ends)
    return interval_index


def _build_chrom_interval_index(reference: dict[str, list[SiteRecord]]) -> tuple[list[str], dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]]:
    motifs = sorted(reference)
    motif_to_idx = {motif: idx for idx, motif in enumerate(motifs)}
    rows_by_chrom: dict[str, list[tuple[int, int, int]]] = {}
    for motif in motifs:
        motif_idx = motif_to_idx[motif]
        for record in reference[motif]:
            rows_by_chrom.setdefault(record.chrom, []).append((record.start, record.end, motif_idx))
    chrom_index: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for chrom, rows in rows_by_chrom.items():
        starts = np.asarray([row[0] for row in rows], dtype=np.int64)
        ends = np.asarray([row[1] for row in rows], dtype=np.int64)
        motif_idx = np.asarray([row[2] for row in rows], dtype=np.int64)
        chrom_index[chrom] = (starts, ends, motif_idx)
    return motifs, chrom_index


def _count_one_sample_from_interval_index(bam_path: str, sample_name: str, read_shift, min_mapq: int) -> tuple[str, dict[str, int]]:
    cuts = _cut_positions_by_chrom(bam_path, _COUNT_CHROMS or set(), tuple(read_shift), min_mapq)
    motifs = _COUNT_MOTIFS or []
    totals = np.zeros(len(motifs), dtype=np.int64)
    for chrom, (starts, ends, motif_idx) in (_COUNT_INTERVAL_INDEX or {}).items():
        arr = cuts.get(chrom)
        if arr is None or arr.size == 0:
            continue
        interval_counts = np.searchsorted(arr, ends, side="left") - np.searchsorted(arr, starts, side="left")
        totals += np.bincount(motif_idx, weights=interval_counts, minlength=len(motifs)).astype(np.int64)
    sample_counts = {motif: int(totals[idx]) for idx, motif in enumerate(motifs)}
    return sample_name, sample_counts


def build_cutcount_matrix(reference: dict[str, list[SiteRecord]], bam_paths: list[str], sample_names: list[str], read_shift=(4, -5), min_mapq: int = 30, cores: int = 1, logger=None) -> pd.DataFrame:
    if len(bam_paths) != len(sample_names):
        raise ValueError("--count-bams must match sample count")
    global _COUNT_INTERVAL_INDEX, _COUNT_CHROMS, _COUNT_MOTIFS
    _COUNT_MOTIFS, _COUNT_INTERVAL_INDEX = _build_chrom_interval_index(reference)
    _COUNT_CHROMS = {record.chrom for records in reference.values() for record in records}
    matrix = pd.DataFrame(0, index=sorted(reference), columns=sample_names, dtype=np.int64)
    tasks = list(zip(bam_paths, sample_names))
    worker_count = min(max(1, int(cores)), len(tasks))
    if logger is not None:
        logger.info(f"Counting shifted Tn5 insertions for {len(tasks)} BAMs using {worker_count} worker(s)")
    if worker_count == 1:
        results = []
        for bam_path, sample_name in tasks:
            if logger is not None:
                logger.info(f"Counting shifted Tn5 insertions for {sample_name}")
            results.append(_count_one_sample_from_interval_index(bam_path, sample_name, read_shift, min_mapq))
            if logger is not None:
                logger.info(f"Finished raw cut counting for {sample_name}")
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=worker_count) as executor:
            futures = [
                executor.submit(_count_one_sample_from_interval_index, bam_path, sample_name, read_shift, min_mapq)
                for bam_path, sample_name in tasks
            ]
            results = []
            for future in concurrent.futures.as_completed(futures):
                sample_name, sample_counts = future.result()
                results.append((sample_name, sample_counts))
                if logger is not None:
                    logger.info(f"Finished raw cut counting for {sample_name}")
    for sample_name, sample_counts in results:
        for motif, count in sample_counts.items():
            matrix.at[motif, sample_name] = count
    return matrix


def _bigwig_window_mean(bw, record: SiteRecord) -> float:
    chroms = bw.chroms()
    if record.chrom not in chroms:
        return np.nan
    end = min(record.end, chroms[record.chrom])
    if end <= record.start:
        return np.nan
    value = bw.stats(record.chrom, record.start, end, type="mean")[0]
    return _as_float(value)


def build_score_matrix(reference: dict[str, list[SiteRecord]], score_signals: list[str], sample_names: list[str]) -> pd.DataFrame:
    if len(score_signals) != len(sample_names):
        raise ValueError("--score-signals must match sample count")
    matrix = pd.DataFrame(np.nan, index=sorted(reference), columns=sample_names, dtype=float)
    for signal, sample_name in zip(score_signals, sample_names):
        with pyBigWig.open(signal) as bw:
            for motif, records in reference.items():
                vals = [_bigwig_window_mean(bw, record) for record in records]
                finite = np.asarray([v for v in vals if np.isfinite(v)], dtype=float)
                matrix.at[motif, sample_name] = float(np.mean(finite)) if finite.size else np.nan
    return matrix


def build_score_matrix_from_reference_beds(reference: dict[str, list[SiteRecord]], score_reference_dir: str, sample_names: list[str], condition_names: list[str]) -> pd.DataFrame:
    """Build motif x sample score matrix from existing *_all.bed sample-score columns."""

    matrix = pd.DataFrame(np.nan, index=sorted(reference), columns=sample_names, dtype=float)
    n_samples = len(sample_names)
    n_conditions = len(condition_names)
    for motif, records in reference.items():
        wanted = {(r.chrom, r.center, r.strand) for r in records}
        sums = np.zeros(n_samples, dtype=float)
        counts = 0
        bed = Path(score_reference_dir) / motif / "beds" / f"{motif}_all.bed"
        if not bed.exists():
            continue
        with bed.open("r", encoding="utf-8") as handle:
            for line in handle:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 6 + n_samples:
                    continue
                try:
                    center = (int(fields[1]) + int(fields[2])) // 2
                except ValueError:
                    continue
                strand = fields[5] if fields[5] in {"+", "-"} else "."
                if (fields[0], center, strand) not in wanted:
                    continue
                sample_start = len(fields) - n_samples - (2 * n_conditions)
                if sample_start < 6:
                    continue
                vals = np.asarray([_as_float(v) for v in fields[sample_start:sample_start + n_samples]], dtype=float)
                if np.isfinite(vals).all():
                    sums += vals
                    counts += 1
        if counts > 0:
            matrix.loc[motif, sample_names] = sums / float(counts)
    return matrix


def run_pydeseq2(counts_by_motif: pd.DataFrame, sample_names: list[str], sample_to_condition: dict[str, str], comparison: tuple[str, str], cores: int) -> pd.DataFrame:
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except Exception as exc:
        raise ImportError('PyDESeq2 is required for --method deseq2-cutcount. Install with: pip install "fp-tools-bio[deseq2]"') from exc

    count_table = counts_by_motif.T.loc[sample_names].astype(int)
    metadata = pd.DataFrame({"condition": [sample_to_condition[s] for s in sample_names]}, index=sample_names)
    dds = DeseqDataSet(
        counts=count_table,
        metadata=metadata,
        design="~condition",
        ref_level=["condition", comparison[1]],
        min_replicates=2,
        n_cpus=max(1, int(cores)),
        quiet=True,
    )
    dds.deseq2()
    stats = DeseqStats(dds, contrast=["condition", comparison[0], comparison[1]], n_cpus=max(1, int(cores)), quiet=True)
    stats.summary()
    result = stats.results_df.copy()
    result.index.name = "output_prefix"
    return result.reset_index()


def moderated_footprint_score(score_by_motif: pd.DataFrame, sample_names: list[str], condition_samples: dict[str, list[str]], comparison: tuple[str, str]) -> pd.DataFrame:
    cond1, cond2 = comparison
    rows = []
    group1 = condition_samples[cond1]
    group2 = condition_samples[cond2]
    df_resid = max(1, len(group1) + len(group2) - 2)
    variances = []
    raw_rows = []
    for motif, values in score_by_motif.iterrows():
        x = pd.to_numeric(values[group1], errors="coerce").to_numpy(dtype=float)
        y = pd.to_numeric(values[group2], errors="coerce").to_numpy(dtype=float)
        x = x[np.isfinite(x)]
        y = y[np.isfinite(y)]
        if len(x) == 0 or len(y) == 0:
            pooled_var = np.nan
        else:
            vx = np.var(x, ddof=1) if len(x) > 1 else 0.0
            vy = np.var(y, ddof=1) if len(y) > 1 else 0.0
            pooled_var = ((len(x) - 1) * vx + (len(y) - 1) * vy) / max(1, len(x) + len(y) - 2)
        raw_rows.append((motif, x, y, pooled_var))
        if np.isfinite(pooled_var) and pooled_var > 0:
            variances.append(pooled_var)

    prior_var = float(np.median(variances)) if variances else 1.0
    if prior_var <= 0 or not np.isfinite(prior_var):
        prior_var = 1.0
    prior_df = 4.0
    total_df = df_resid + prior_df

    for motif, x, y, pooled_var in raw_rows:
        if len(x) == 0 or len(y) == 0:
            effect = np.nan
            pvalue = 1.0
            tstat = np.nan
            se = np.nan
            welch_p = 1.0
        else:
            effect = float(np.mean(x) - np.mean(y))
            if not np.isfinite(pooled_var):
                pooled_var = prior_var
            mod_var = (prior_df * prior_var + df_resid * max(float(pooled_var), 0.0)) / total_df
            se = float(np.sqrt(mod_var * (1.0 / len(x) + 1.0 / len(y))))
            tstat = effect / se if se > 0 else 0.0
            pvalue = float(2.0 * scipy.stats.t.sf(abs(tstat), df=total_df))
            welch = scipy.stats.ttest_ind(x, y, equal_var=False, nan_policy="omit")
            welch_p = float(welch.pvalue) if np.isfinite(welch.pvalue) else 1.0
        rows.append(
            {
                "output_prefix": motif,
                "footprint_score_delta": effect,
                "moderated_t": tstat,
                "moderated_se": se,
                "moderated_df": total_df,
                "pvalue": pvalue,
                "welch_pvalue": welch_p,
                f"{cond1}_mean_score": float(np.mean(x)) if len(x) else np.nan,
                f"{cond2}_mean_score": float(np.mean(y)) if len(y) else np.nan,
            }
        )
    result = pd.DataFrame(rows)
    result["padj"] = benjamini_hochberg(result["pvalue"].to_numpy(dtype=float))
    return result


def _load_logo_lookup(motif_paths, naming: str) -> dict[str, object]:
    lookup = {}
    if not motif_paths:
        return lookup
    try:
        motifs = MotifList()
        for path in motif_paths:
            motifs += MotifList().from_file(path)
        for motif in motifs:
            motif.set_prefix(naming)
            lookup[motif.prefix] = motif
    except Exception:
        return {}
    return lookup


def _result_to_info_table(results: pd.DataFrame, metadata: dict[str, dict[str, str]], reference: dict[str, list[SiteRecord]], comparison: tuple[str, str], method: str) -> pd.DataFrame:
    cond1, cond2 = comparison
    base = f"{cond1}_{cond2}"
    rows = []
    for _, row in results.iterrows():
        prefix = str(row["output_prefix"])
        meta = metadata.get(prefix, {})
        if method == "deseq2-cutcount":
            change = _as_float(row.get("log2FoldChange"), 0.0)
            pvalue = _as_float(row.get("pvalue"), 1.0)
            qvalue = _as_float(row.get("padj"), 1.0)
        else:
            change = _as_float(row.get("footprint_score_delta"), 0.0)
            pvalue = _as_float(row.get("pvalue"), 1.0)
            qvalue = _as_float(row.get("padj"), 1.0)
        rows.append(
            {
                "output_prefix": prefix,
                "name": meta.get("name", prefix.rsplit("_", 1)[0] if "_" in prefix else prefix),
                "motif_id": meta.get("motif_id", ""),
                "cluster": meta.get("cluster", prefix),
                "total_tfbs": len(reference.get(prefix, [])),
                f"{cond1}_bound": len(reference.get(prefix, [])),
                f"{cond2}_bound": len(reference.get(prefix, [])),
                f"{base}_change": round(float(change), 5),
                f"{base}_pvalue": f"{max(float(pvalue), 1e-308):.5E}",
                f"{base}_qvalue_bh": f"{float(qvalue) if np.isfinite(qvalue) else 1.0:.5E}",
                f"{base}_significant_fdr05": bool(np.isfinite(qvalue) and qvalue <= 0.05),
                f"{base}_mean_delta_fp": round(float(change), 5),
                f"{base}_mean_log2fc": round(float(change), 5),
            }
        )
    info_table = pd.DataFrame(rows)
    changes = pd.to_numeric(info_table[f"{base}_change"], errors="coerce").fillna(0.0)
    pvals = pd.to_numeric(info_table[f"{base}_pvalue"], errors="coerce").fillna(1.0)
    pval_min = np.percentile(pvals[pvals > 0], 5) if (pvals > 0).any() else 1.0
    lo, hi = np.percentile(changes, [5, 95]) if len(changes) else (0.0, 0.0)
    info_table[f"{base}_highlighted"] = (changes < lo) | (changes > hi) | (pvals < pval_min)
    return info_table


def _motif_records_for_html(info_table: pd.DataFrame, comparison: tuple[str, str], motif_lookup: dict[str, object], outdir: str) -> list[SimpleNamespace]:
    cond1, cond2 = comparison
    base = f"{cond1}_{cond2}"
    records = []
    for _, row in info_table.iterrows():
        prefix = str(row["output_prefix"])
        change = _as_float(row[f"{base}_change"], 0.0)
        pvalue = max(_as_float(row[f"{base}_pvalue"], 1.0), 1e-308)
        highlighted = str(row.get(f"{base}_highlighted", "")).lower() in {"true", "1", "yes"}
        group = f"{cond2}_up" if highlighted and change < 0 else (f"{cond1}_up" if highlighted else "n.s.")
        logo_path = Path(outdir) / prefix / f"{prefix}.png"
        logo = ""
        if logo_path.exists():
            logo = base64.b64encode(logo_path.read_bytes()).decode("utf-8")
        motif = motif_lookup.get(prefix)
        records.append(
            SimpleNamespace(
                prefix=prefix,
                name=str(row.get("name", prefix)),
                id=str(row.get("motif_id", "")),
                change=change,
                pvalue=pvalue,
                logpvalue=-np.log10(pvalue),
                highlighted=highlighted,
                group=group,
                base=logo,
                counts=getattr(motif, "counts", None),
            )
        )
    return records


def run_fixed_site_differential(args) -> None:
    logger = FpToolsLogger("diff-footprints", args.verbosity)
    logger.begin()
    args.cores = check_cores(args.cores, logger)
    make_directory(args.outdir)

    if args.method == "deseq2-cutcount":
        if not args.count_bams:
            raise ValueError("--count-bams is required for --method deseq2-cutcount")
        args.signals = list(args.count_bams)
    elif args.method == "footprint-score":
        if not args.score_signals and not getattr(args, "score_reference_dir", None):
            raise ValueError("--score-signals or --score-reference-dir is required for --method footprint-score")
        args.signals = list(args.score_signals or args.aggregate_signals or [])
    else:
        raise ValueError(f"Unsupported fixed-site method: {args.method}")

    from fp_tools.tools.bindetect import _prepare_condition_metadata

    _prepare_condition_metadata(args)
    if len(args.cond_names) != 2:
        raise ValueError(f"--method {args.method} currently supports exactly two conditions")

    logger.info("Loading fixed motif-site reference")
    reference, metadata = load_motif_site_reference(args.site_reference_dirs or [], args.reference_site_set, args.site_window)
    logger.info(f"Loaded {len(reference)} motifs with {sum(len(records) for records in reference.values()):,} fixed motif-site windows")
    logger.info("Writing fixed-site BED reference for aggregate reports")
    write_reference_beds(reference, args.outdir, args.cond_names)
    reference_rows = [
        {"output_prefix": motif, "n_sites": len(records)}
        for motif, records in sorted(reference.items())
    ]
    pd.DataFrame(reference_rows).to_csv(Path(args.outdir) / f"{args.prefix}_fixed_site_reference_summary.tsv", sep="\t", index=False)

    sample_names = list(args.sample_names)
    comparison = args.comparisons[0]
    if args.method == "deseq2-cutcount":
        matrix = build_cutcount_matrix(reference, list(args.count_bams), sample_names, tuple(args.read_shift), int(args.min_mapq), args.cores, logger)
        matrix.to_csv(Path(args.outdir) / f"{args.prefix}_motif_cutcount_matrix.tsv", sep="\t")
        logger.info("Running PyDESeq2 on motif-site cut counts")
        method_results = run_pydeseq2(matrix, sample_names, args.sample_to_condition, comparison, args.cores)
        method_results.to_csv(Path(args.outdir) / f"{args.prefix}_deseq2_native_results.tsv", sep="\t", index=False)
        label = f"Method: motif-site raw Tn5 cut-count DESeq2; reference: {args.reference_site_set}; window: {args.site_window} bp"
    else:
        if getattr(args, "score_reference_dir", None):
            logger.info("Building footprint-score matrix from existing motif-site score columns")
            matrix = build_score_matrix_from_reference_beds(reference, args.score_reference_dir, sample_names, args.cond_names)
        else:
            logger.info("Building footprint-score matrix from footprint-score bigWigs")
            matrix = build_score_matrix(reference, list(args.score_signals), sample_names)
        matrix.to_csv(Path(args.outdir) / f"{args.prefix}_motif_footprint_score_matrix.tsv", sep="\t")
        logger.info("Running empirical-Bayes moderated test on motif footprint scores")
        method_results = moderated_footprint_score(matrix, sample_names, args.condition_samples, comparison)
        method_results.to_csv(Path(args.outdir) / f"{args.prefix}_footprint_score_ebayes_results.tsv", sep="\t", index=False)
        label = f"Method: motif-site footprint-score empirical Bayes; reference: {args.reference_site_set}; window: {args.site_window} bp"

    info_table = _result_to_info_table(method_results, metadata, reference, comparison, args.method)
    results_path = Path(args.outdir) / f"{args.prefix}_results.txt"
    info_table.to_csv(results_path, sep="\t", index=False)

    motif_lookup = _load_logo_lookup(getattr(args, "motifs", None), args.naming)
    motif_records = _motif_records_for_html(info_table, comparison, motif_lookup, args.outdir)
    aggregate_data = None
    if getattr(args, "aggregate_signals", None) and getattr(args, "plot_aggregate", "off") != "off":
        aggregate_data = build_bindetect_aggregate_payload(motif_records, info_table, comparison, args)
    html_out = Path(args.outdir) / f"{args.prefix}_{comparison[0]}_{comparison[1]}.html"
    plot_interactive_bindetect(
        motif_records,
        list(comparison),
        str(html_out),
        aggregate_data=aggregate_data,
        title="Differential footprint report",
        report_label=getattr(args, "report_label", None) or label,
    )
    logger.info(f"Wrote {results_path}")
    logger.info(f"Wrote {html_out}")
    logger.end()
