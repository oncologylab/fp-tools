#!/usr/bin/env python3
"""Build and audit the frozen genomic expectation on training coordinates only."""

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rebuild_corrected_safety_profiles as parent
import build_strand_functional_profiles as strands
import global_parametric_expectation as operator
from audit_chrombpnet_comparison import canonical_id, checked_path
from audit_direct_bias_safety import record, write_immutable

ROOT = parent.ROOT
PROTOCOL = ROOT / "benchmarks/manifests/global_parametric_track_development_v1.json"
PROTOCOL_SHA = "5ce5e9afe15493b57b49ff6683f6f29d801f34728675e9f923ed567084bd499f"
SOURCE = ROOT / "benchmarks/results/corrected_safety_profiles_20260922_v1"
SOURCE_SHA = "dc053ee9d0e0d6154cabe0e3915a6a78ba39887ffe28e2850a5059db4aa3a9e9"
SITE_KEYS = parent.base.SITE_KEYS


def training_sites(path, cell, chromosomes):
    # Do not deserialize any ChIP label column, including tuning labels.
    sites = pd.read_csv(path, sep="\t", usecols=SITE_KEYS)
    if set(sites.cell) != {cell} or sites.duplicated(SITE_KEYS).any():
        raise ValueError("invalid source coordinate support")
    selection = sites.TFBS_chr.isin(chromosomes).to_numpy()
    selected = sites.loc[selection].reset_index(drop=True)
    if selected.empty or any("label" in c.lower() or "chip" in c.lower() for c in selected.columns):
        raise ValueError("training stage requires label-free coordinates")
    return selected, selection


def infer_log_bias(sites, model, genome):
    width, flank = 401, 200
    margin = max(41, model.feature_spec.context_length // 2 + 1)
    plus = np.full((len(sites), width), np.nan)
    minus = np.full_like(plus, np.nan)
    with strands.open_fasta(genome) as fasta:
        lengths = dict(zip(fasta.references, fasta.lengths))
        for index, row in enumerate(sites.itertuples(index=False)):
            center = (int(row.TFBS_start) + int(row.TFBS_end)) // 2
            start, end = center-flank-margin, center+flank+margin+1
            if row.TFBS_chr not in lengths or start < 0 or end > lengths[row.TFBS_chr]:
                continue
            first, second, valid = strands.strand_log_bias(model,
                fasta.fetch(row.TFBS_chr, start, end).upper(), margin+np.arange(width))
            plus[index, valid], minus[index, valid] = first[valid], second[valid]
    return plus, minus


def strand_expectation(observed, log_bias):
    if observed.shape != log_bias.shape or observed.ndim != 2 or observed.shape[1] != 401:
        raise ValueError("401-position genomic strand arrays required")
    expected = np.full((len(observed), 201), np.nan)
    valid = np.zeros(expected.shape, bool)
    for index, (counts, bias) in enumerate(zip(observed, log_bias)):
        result, available = operator.expected_track(counts, bias)
        expected[index], valid[index] = result[100:301], available[100:301]
        # Independent scalar calculation at every site's central base.
        flank = np.r_[np.arange(100, 151), np.arange(250, 301)]
        supported = np.isfinite(bias[100:301]).all() and counts[flank].sum() > 0
        if bool(available[200]) != bool(supported):
            raise ValueError("scalar support differs")
        if supported:
            offset = bias[100:301].max()
            scalar = np.exp(bias[200]-offset) * counts[flank].sum() / np.exp(bias[flank]-offset).sum()
            np.testing.assert_allclose(result[200], scalar, rtol=1e-10, atol=1e-10)
        # Fixed coordinate-order subset; no label or score selection.
        if index < 64:
            for center in (120, 200, 280):
                sub, mask = operator.expected_track(counts[center-100:center+101], bias[center-100:center+101])
                if mask[100] != available[center]:
                    raise ValueError("halo chunk support differs")
                np.testing.assert_allclose(sub[100], result[center], rtol=1e-10, atol=1e-10, equal_nan=True)
    return expected, valid


def overlap_summary(sites, observed, expected, valid):
    rows = sites.reset_index(drop=True).copy()
    rows["row"] = np.arange(len(rows))
    rows["center"] = (rows.TFBS_start + rows.TFBS_end) // 2
    pairs, shared, missing = 0, 0, 0
    for _, group in rows.groupby("TFBS_chr"):
        records = list(group.sort_values(["center", "tf", "TFBS_strand"], kind="stable").itertuples(index=False))
        for left, right in zip(records, records[1:]):
            shift = int(right.center-left.center)
            if shift >= 201:
                continue
            a, b = slice(shift, 201), slice(0, 201-shift)
            np.testing.assert_array_equal(observed[left.row, a], observed[right.row, b])
            np.testing.assert_array_equal(valid[left.row, a], valid[right.row, b])
            mask = valid[left.row, a]
            np.testing.assert_allclose(expected[left.row, a][mask], expected[right.row, b][mask], rtol=1e-8, atol=1e-10)
            pairs += 1
            shared += int(mask.sum())
            missing += int((~mask).sum())
    return dict(overlap_pairs=pairs, shared_valid_positions=shared, shared_unavailable_positions=missing,
                observed_conflicts=0, expected_conflicts=0, validity_conflicts=0)


def freeze(root, out):
    protocol = json.loads(checked_path(dict(path=str(PROTOCOL), sha256=PROTOCOL_SHA), ROOT).read_text())
    old = json.loads(checked_path(dict(path=str(SOURCE / "configuration.freeze.json"), sha256=SOURCE_SHA), ROOT).read_text())
    if out.exists():
        raise ValueError("use a fresh output directory")
    inputs, cells = {}, {}

    def use(entry):
        p = checked_path(entry, root)
        inputs[str(p)] = record(p)
        return record(p)

    use(record(PROTOCOL)); use(record(SOURCE / "configuration.freeze.json"))
    for item in old["items"]:
        if item["role"] != "validation":
            continue
        cell = item["cell"]
        doc = json.loads(checked_path(use(item["manifest"]), root).read_text())
        meta = doc["metadata"]
        if any(meta[k] != protocol[k] for k in ("read_shift", "minimum_mapq", "keep_duplicates", "bias_model_sha256")):
            raise ValueError("scientific input convention differs")
        model = use(dict(path=meta["bias_model"], sha256=meta["bias_model_sha256"]))
        genome = use(dict(path=meta["genome"], sha256=meta["genome_sha256"]))
        sidecar = use(record(Path(model["path"]).with_suffix(".json")))
        fai = use(record(Path(genome["path"] + ".fai")))
        corrected = use(record(SOURCE / f"{cell}.validation.json"))
        stage = json.loads(Path(corrected["path"]).read_text())
        if stage["outputs"]["sites"] != item["sites"]:
            raise ValueError("corrected and source site identities differ")
        cells[cell] = dict(bam=use(item["bam"]), bai=use(item["bai"]), model=model, sidecar=sidecar,
            genome=genome, fai=fai, sites=use(item["sites"]), corrected_profiles=use(stage["outputs"]["profiles"]))
        sites, _ = training_sites(cells[cell]["sites"]["path"], cell, protocol["initial_evaluation_chromosomes"])
        cells[cell]["training_site_hashes"] = strands.site_hashes(sites).tolist()
        print(f"Verified training inputs: {cell}, {len(sites)} label-free sites", flush=True)
    if sorted(cells) != protocol["cells"]:
        raise ValueError("development cells differ")
    paths = {e["path"] for e in old["code"]}
    paths.update([str(Path(__file__)), operator.__file__, strands.__file__])
    config = dict(protocol=protocol, cells=cells, inputs=list(inputs.values()),
        code=[record(Path(p)) for p in sorted(paths)], runtime=parent.extension.runtime(),
        profiles_scored=False, model_refitted=False, labels_read=False, promotion_eligible=False)
    config["configuration_id"] = canonical_id(config)
    out.mkdir(parents=True)
    write_immutable(out / "configuration.freeze.json", config)
    print("Frozen complete training-only execution configuration", flush=True)


def execute(root, out, verify_only=False):
    config = parent.extension.read_document(out / "configuration.freeze.json", "configuration_id")
    protocol = json.loads(checked_path(dict(path=str(PROTOCOL), sha256=PROTOCOL_SHA), ROOT).read_text())
    if (config["protocol"] != protocol or config["runtime"] != parent.extension.runtime()
            or any(config[k] is not False for k in ("profiles_scored", "model_refitted", "labels_read", "promotion_eligible"))):
        raise ValueError("frozen scope or runtime differs")
    for item in config["inputs"] + config["code"]:
        checked_path(item, root)
    stages = []
    for cell, item in config["cells"].items():
        receipt, target = out / f"{cell}.json", out / f"{cell}.npz"
        if receipt.exists():
            parent.base.verify_stage(receipt, config["configuration_id"], root)
        elif verify_only or target.exists() or (out / f"{cell}.sites.tsv.gz").exists():
            raise ValueError("missing result or partial output preserved")
        else:
            sites, selection = training_sites(item["sites"]["path"], cell, protocol["initial_evaluation_chromosomes"])
            if strands.site_hashes(sites).tolist() != item["training_site_hashes"]:
                raise ValueError("training coordinates changed")
            plus, minus, bam_valid = strands.extract_strand_cut_profiles(sites, item["bam"]["path"],
                flank=200, read_shift=(4, -5), minimum_mapq=30, keep_duplicates=False)
            model = strands.ConditionalSequenceBiasModel.load(item["model"]["path"])
            first, second = infer_log_bias(sites, model, item["genome"]["path"])
            plus_expected, plus_valid = strand_expectation(plus, first)
            minus_expected, minus_valid = strand_expectation(minus, second)
            observed = plus[:, 100:301] + minus[:, 100:301]
            expected = plus_expected + minus_expected
            valid = plus_valid & minus_valid & bam_valid[:, None]
            expected[~valid] = np.nan
            with np.load(item["corrected_profiles"]["path"], allow_pickle=False) as old:
                np.testing.assert_array_equal(old["site_hash"][selection], strands.site_hashes(sites))
                original = old["observed"][selection].copy()
            reverse = sites.TFBS_strand.eq("-").to_numpy()
            original[reverse] = original[reverse, ::-1]
            np.testing.assert_array_equal(observed[bam_valid], original[bam_valid])
            audit = overlap_summary(sites.loc[bam_valid], observed[bam_valid], expected[bam_valid], valid[bam_valid])
            audit.update(cell=cell, sites=len(sites), bam_valid_sites=int(bam_valid.sum()),
                fully_supported_sites=int(valid.all(axis=1).sum()), supported_positions=int(valid.sum()),
                scalar_centers_checked_per_strand=len(sites), halo_sites_checked_per_strand=min(64, len(sites)))
            np.savez_compressed(target, observed=observed, expected=expected, position_valid=valid,
                plus_observed=plus, minus_observed=minus, plus_log_bias=first, minus_log_bias=second,
                plus_expected=plus_expected, minus_expected=minus_expected, bam_valid=bam_valid,
                valid=valid.all(axis=1), site_hash=strands.site_hashes(sites))
            site_path = out / f"{cell}.sites.tsv.gz"
            sites.to_csv(site_path, sep="\t", index=False, compression={"method": "gzip", "mtime": 0})
            stage = dict(configuration_id=config["configuration_id"], cell=cell, audit=audit,
                outputs=dict(profiles=record(target), sites=record(site_path)), profiles_scored=False, labels_read=False)
            stage["stage_id"] = canonical_id(stage)
            write_immutable(receipt, stage)
        stages.append(record(receipt))
        print("Verified global training profiles:", cell, flush=True)
    result = dict(configuration_id=config["configuration_id"], stages=stages, outputs={},
        profiles_scored=False, model_refitted=False, labels_read=False, promotion_eligible=False)
    result["stage_id"] = canonical_id(result)
    path = out / "profiles_manifest.json"
    if verify_only:
        if parent.base.verify_stage(path, config["configuration_id"], root) != result:
            raise ValueError("completed stage identities differ")
    else:
        write_immutable(path, result)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--freeze-only", action="store_true")
    modes.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    out = args.outdir.resolve()
    if not out.is_relative_to(ROOT / "benchmarks/results"):
        parser.error("use research results")
    if args.freeze_only:
        freeze(args.data_root.resolve(), out)
    else:
        with parent.extension.exclusive_run(out):
            execute(args.data_root.resolve(), out, args.verify_only)
