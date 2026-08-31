#!/usr/bin/env python3
"""Build a combined ChIP-labelled motif-site matrix for ENCODE TF tasks."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from build_footprint_site_labels import label_sites, propensity_match, read_peaks, read_sites


def motif_site_file(root: Path, motif_id: str) -> Path:
    matches = sorted(root.glob(f"*_{motif_id}/beds/*_{motif_id}_all.bed"))
    if len(matches) != 1:
        raise ValueError(f"expected one all-site BED for {motif_id}; found {len(matches)}")
    return matches[0]


def cell_motif_root(root: Path, cell: str) -> Path:
    candidate = root / cell
    return candidate if candidate.is_dir() else root


def build_matrix(
    study: dict[str, object],
    selected_peaks: pd.DataFrame,
    motif_root: Path,
    split: str,
    positive_distance: int,
    negative_distance: int,
    negative_ratio: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    natural_parts = []
    matched_parts = []
    audit_parts = []
    peaks_by_task = selected_peaks.set_index(["cell", "tf"])
    for task in study["tasks"]:
        if task["split"] != split:
            continue
        cell, tf, motif_id = str(task["cell"]), str(task["tf"]), str(task["motif_id"])
        if (cell, tf) not in peaks_by_task.index:
            raise ValueError(f"selected peak manifest lacks {cell}/{tf}")
        peak_row = peaks_by_task.loc[(cell, tf)]
        sites = read_sites(motif_site_file(cell_motif_root(motif_root, cell), motif_id))
        labelled = label_sites(
            sites,
            read_peaks(Path(str(peak_row.local_path))),
            positive_summit_distance=positive_distance,
            negative_peak_distance=negative_distance,
        )
        labelled.insert(0, "cell", cell)
        labelled.insert(1, "tf", tf)
        labelled.insert(2, "motif", motif_id)
        labelled.insert(3, "chip_accession", str(peak_row.file_accession))
        labelled["motif_family"] = str(task["motif_family"])
        labelled["role"] = str(task["role"])
        labelled["study_split"] = split
        natural = labelled[labelled["label"].isin([0, 1])].copy()
        positive = natural[natural["label"] == 1]
        negative = natural[natural["label"] == 0]
        maximum_positives = len(negative) // negative_ratio
        matching_input = natural
        if len(positive) > maximum_positives:
            positive = positive.sample(n=maximum_positives, random_state=seed)
            matching_input = pd.concat([positive, negative], ignore_index=True)
        matched = propensity_match(
            matching_input,
            ["motif_score"],
            negative_ratio=negative_ratio,
            seed=seed,
        )
        natural_parts.append(natural)
        matched_parts.append(matched)
        counts = labelled["label_reason"].value_counts().rename_axis("label_reason").reset_index(name="sites")
        counts.insert(0, "cell", cell)
        counts.insert(1, "tf", tf)
        counts["chip_accession"] = str(peak_row.file_accession)
        audit_parts.append(counts)
    natural = pd.concat(natural_parts, ignore_index=True)
    matched = pd.concat(matched_parts, ignore_index=True)
    audit = pd.concat(audit_parts, ignore_index=True)
    rename = {
        "chrom": "TFBS_chr", "start": "TFBS_start", "end": "TFBS_end",
        "strand": "TFBS_strand", "label": "chip_label",
    }
    return natural.rename(columns=rename), matched.rename(columns=rename), audit


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, required=True)
    parser.add_argument("--selected-peaks", type=Path, required=True)
    parser.add_argument("--motif-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--split", default="development")
    parser.add_argument("--positive-summit-distance", type=int, default=100)
    parser.add_argument("--negative-peak-distance", type=int, default=500)
    parser.add_argument("--negative-ratio", type=int, default=1)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args(argv)

    study = json.loads(args.study.read_text(encoding="utf-8"))
    selected = pd.read_csv(args.selected_peaks, sep="\t")
    natural, matched, audit = build_matrix(
        study, selected, args.motif_root, args.split,
        args.positive_summit_distance, args.negative_peak_distance,
        args.negative_ratio, args.seed,
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    natural.to_csv(args.outdir / "encode_tf_site_labels.tsv.gz", sep="\t", index=False, compression="gzip")
    matched.to_csv(args.outdir / "encode_tf_site_labels_motif_matched.tsv.gz", sep="\t", index=False, compression="gzip")
    audit.to_csv(args.outdir / "encode_tf_site_label_audit.tsv", sep="\t", index=False)
    print(audit.pivot_table(index=["cell", "tf"], columns="label_reason", values="sites", fill_value=0).to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
