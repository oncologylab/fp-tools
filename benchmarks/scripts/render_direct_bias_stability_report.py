#!/usr/bin/env python3
"""Render a narrow, safety-linked CTCF report from completed development runs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))
sys.path.insert(0, str(SCRIPT_DIR))

from audit_chrombpnet_comparison import canonical_id, check_identity, checked_path, file_sha256
from audit_direct_bias_safety import read_thresholds, record, write_immutable
from audit_direct_bias_stability import METHODS, PARAMETRIC, RAW, DWM, SAMPLE_DWM, RUN_KEYS, SCHEMA, validate_plan, verify_stage
from render_chrombpnet_parametric_report import validate_safety_link

COHORT = "motif_accessibility_observed_total"
KEYS = RUN_KEYS + ["cohort", "tf"]
SAFETY_FIELDS = ["false_positive_calls", "informative_support", "informative_false_positive_rate_upper_95",
                 "paired_false_positive_calls", "paired_support", "paired_false_positive_rate_upper_95",
                 "paired_false_positive_rate_increase_over_dwm"]


def ctcf(frame):
    return frame.loc[frame.tf.eq("CTCF") & frame.cohort.eq(COHORT)].copy()


def select_evidence(metrics, shapes, bootstrap, coverage, safety):
    full = ctcf(metrics).loc[lambda f: f.depth.eq("full")]
    expected = {(cell, f"{cell}_rep{rep}", "full", 2026) for cell in ("HepG2", "K562") for rep in (1, 2, 3)}
    if (full.duplicated(RUN_KEYS + ["method"]).any() or len(full) != 24
            or set(full[RUN_KEYS].itertuples(index=False, name=None)) != expected
            or not full.groupby(RUN_KEYS).method.apply(lambda x: set(x) == set(METHODS)).all()
            or not full.status.eq("eligible").all()):
        raise ValueError("report requires all six eligible full biological replicates and four methods")
    candidates = ctcf(metrics).loc[lambda f: f.method.eq(PARAMETRIC)].copy()
    if len(candidates) != 36 or candidates.duplicated(RUN_KEYS).any():
        raise ValueError("report requires all 36 CTCF point-estimate runs")
    for baseline in ("fixed_dwm", "sample_dwm", "raw"):
        if not (candidates[f"auroc_gain_over_{baseline}"].gt(0).all()
                and candidates[f"relative_auprc_gain_over_{baseline}"].gt(0).all()):
            raise ValueError("CTCF direction is not stable against every prespecified baseline")
    shape = ctcf(shapes).loc[lambda f: f.method.eq(PARAMETRIC)]
    if (len(shape) != 36 or shape.duplicated(RUN_KEYS).any()
            or not shape.canonical_protection.eq(True).all()):
        raise ValueError("report requires canonical CTCF protection across the available matrix")
    candidates = candidates.merge(shape, on=KEYS + ["method"], validate="one_to_one")
    balance = ctcf(coverage).drop(columns=["n_positive", "n_negative", "status"], errors="ignore")
    if len(balance) != 36 or not np.isfinite(balance.observed_total_smd).all():
        raise ValueError("report requires every run's observed-count balance audit")
    candidates = candidates.merge(balance, on=KEYS, validate="one_to_one")
    # Matching was frozen at 10M. Drift is reported, never silently repaired.
    candidates["current_cut_balance_exceeds_0p1"] = candidates.observed_total_smd.abs().gt(.1)
    paired = ctcf(bootstrap)
    for baseline, suffix in ((DWM, "fixed_dwm"), (SAMPLE_DWM, "sample_dwm"), (RAW, "raw")):
        selected = paired.loc[paired.baseline.eq(baseline)]
        full_pairs = selected.loc[selected.depth.eq("full")]
        if (set(full_pairs[RUN_KEYS].itertuples(index=False, name=None)) != expected
                or full_pairs.duplicated(RUN_KEYS).any()
                or not full_pairs.auroc_gain_lower_95.gt(0).all()
                or not full_pairs.relative_auprc_gain_lower_95.gt(0).all()):
            raise ValueError("report requires significant full-replicate gains over DWM and raw")
        fields = [f"{metric}_{bound}_95" for metric in ("auroc_gain", "relative_auprc_gain") for bound in ("lower", "upper")]
        candidates = candidates.merge(selected[KEYS + fields].rename(columns={key: f"{key}_vs_{suffix}" for key in fields}),
                                      on=KEYS, how="left", validate="one_to_one")
    safe = safety.loc[safety.tf.eq("CTCF")]
    if (len(safe) != 2 or set(safe.cell) != {"HepG2", "K562"} or safe.cell.duplicated().any()
            or not safe.passes_candidate_screen.eq(True).all()
            or not safe.informative_false_positive_rate_upper_95.le(.05).all()
            or not safe.paired_false_positive_rate_upper_95.le(.05).all()
            or not safe.paired_false_positive_rate_increase_over_dwm.le(.01).all()):
        raise ValueError("report requires passing exact-operator candidate safety in both cells")
    candidates = candidates.merge(safe[["cell"] + SAFETY_FIELDS], on="cell", validate="many_to_one")
    candidates["safety_scope"] = "existing_naked_rep2_candidate_screen_not_final_untouched_safety"
    return candidates.sort_values(RUN_KEYS).reset_index(drop=True)


def render(summary, profiles, pdf, png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(8.27, 11.69))
    fig.text(.07, .958, "CTCF: direct-bias stability", fontsize=20, weight="bold")
    fig.text(.07, .925, "Frozen CPU model • examined development chromosomes • research only", fontsize=10)
    fig.text(.07, .897, "Before/after: sample-matched DWM vs frozen parametric bias, with identical deviance scoring.\n"
             "Full-depth biological replicates; sites and matched pairs stay fixed from the 10M reference.", fontsize=9)
    selected = ctcf(profiles).loc[lambda f: f.depth.eq("full")]
    for column, cell in enumerate(("HepG2", "K562")):
        ax = fig.add_axes([.13 + .45 * column, .615, .32, .23])
        for method, label, color in ((SAMPLE_DWM, "Before: sample DWM", "#777777"),
                                     (PARAMETRIC, "After: parametric", "#D55E00")):
            group = selected.loc[selected.cell.eq(cell) & selected.method.eq(method)]
            for _, part in group.groupby("sample"):
                ax.plot(part.position, part.positive_minus_negative, color=color, alpha=.2, linewidth=.6)
            average = group.groupby("position").positive_minus_negative.mean()
            ax.plot(average.index, average.values, color=color, linewidth=1.4, label=label)
        pairs = int(summary.loc[summary.cell.eq(cell), "n_sites"].iloc[0]) // 2
        ax.set(title=f"{cell} CTCF · {pairs} fixed pairs", xlim=(-60, 60), xlabel="Motif offset (bp)")
        ax.axhline(0, color="#888888", linewidth=.5); ax.axvspan(-5, 5, color="#aaaaaa", alpha=.15)
        ax.spines[["top", "right"]].set_visible(False); ax.tick_params(labelsize=8)
        ax.title.set_fontsize(11)
        if column == 0:
            ax.set_ylabel("Positive − negative\n(normalized deviance)", fontsize=9)
            ax.legend(fontsize=8, frameon=False)
    fig.text(.07, .552, "Thin curves show three biological replicates; thick curves show their mean.\n"
             "Central protection persists, with stronger shoulders; the change is not uniformly deeper central depletion.", fontsize=8.5)
    full = summary.loc[summary.depth.eq("full")].sort_values(["cell", "sample"])
    rows = []
    for row in full.itertuples(index=False):
        rows.append([row.sample.replace("_rep", " #"), f"{row.sample_dwm_auroc:.3f} / {row.sample_dwm_auprc:.3f}",
                     f"{row.auroc:.3f} / {row.auprc:.3f}",
                     f"{row.auroc_gain_over_raw:+.3f} [{row.auroc_gain_lower_95_vs_raw:.3f}, {row.auroc_gain_upper_95_vs_raw:.3f}]",
                     f"{100 * row.relative_auprc_gain_over_raw:+.1f}%"])
    ax = fig.add_axes([.065, .342, .88, .178]); ax.axis("off")
    table = ax.table(cellText=rows, colLabels=["Full replicate", "DWM AUC / AP", "CPU AUC / AP", "ΔAUC vs raw [95% CI]", "ΔAP vs raw"],
                     colWidths=[.16, .18, .18, .32, .16], cellLoc="center", bbox=[0, 0, 1, 1])
    table.auto_set_font_size(False); table.set_fontsize(7.8)
    for (row, col), item in table.get_celld().items():
        item.set_linewidth(.3); item.set_edgecolor("#cccccc")
        if row == 0:
            item.set_facecolor("#ececec"); item.set_text_props(weight="bold")
    fig.text(.07, .314, "All six full-replicate AUC and relative AP gains over raw and DWM have positive chromosome-block intervals.\n"
             "All 36 CTCF point estimates improve; 10M/25M/50M seeds share one biological library per cell.", fontsize=8.5)
    safe = summary.groupby("cell").first()
    safety_text = "Candidate screen, existing naked-DNA replicate 2 (informative DWM-comparable sites):\n"
    safety_text += ";  ".join(f"{cell}: {int(safe.loc[cell, 'paired_false_positive_calls'])}/{int(safe.loc[cell, 'paired_support'])} calls, "
                             f"Wilson upper {100 * safe.loc[cell, 'paired_false_positive_rate_upper_95']:.2f}%" for cell in ("HepG2", "K562"))
    fig.text(.07, .267, safety_text, fontsize=9)
    drift = int(summary.current_cut_balance_exceeds_0p1.sum())
    maximum = float(summary.observed_total_smd.abs().max())
    fig.text(.07, .17, f"Limits: current cut-total balance exceeds SMD 0.1 in {drift}/36 matched runs (maximum {maximum:.3f}).\n"
             "Matching controls the 10M reference covariates; it does not remove every accessibility cue at other depths.\n"
             "Only CTCF passes the canonical-shape screen among powered TFs. General relative-AP and difficult-family gates fail.\n"
             "Sparse naked-DNA screening does not establish no occupancy leakage or final untouched safety.\n"
             "New-cell labels and final naked-DNA replicate 3 remain unopened; CPU/memory and uncertainty gates remain incomplete.",
             fontsize=8.5, linespacing=1.5)
    fig.text(.07, .065, "Recorded dispersions: HepG2 0; K562 1.61443. Separate development experiments; no cross-cell retuning.\n"
             "AUC = AUROC; AP = AUPRC; AP gains are relative. Five chromosome blocks; intervals are diagnostic.\n"
             "Production DWM defaults remain unchanged. Full metrics, provenance and limitations accompany this PDF.", fontsize=8)
    fig.savefig(pdf, metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(png, dpi=180)
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stability-manifest", type=Path, required=True)
    parser.add_argument("--safety-manifest", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)
    config_path = args.stability_manifest.with_name("configuration.freeze.json")
    config = json.loads(config_path.read_text()); check_identity(config, "configuration_id"); validate_plan(config["settings"])
    result = verify_stage(args.stability_manifest, config["configuration_id"], args.data_root)
    if (result.get("schema") != SCHEMA or result.get("candidate_retuned") is not False
            or result.get("new_holdout_labels_read") is not False or result.get("promotion_eligible") is not False):
        raise ValueError("report accepts only the frozen development stability result")
    for item in result["supports"].values():
        verify_stage(checked_path(item, args.data_root), config["configuration_id"], args.data_root)
    for item in result["stages"]:
        for entry in item.values():
            verify_stage(checked_path(entry, args.data_root), config["configuration_id"], args.data_root)
    safety = json.loads(args.safety_manifest.read_text())
    source = json.loads(checked_path(safety["source_safety"], args.data_root).read_text())
    safe_config = json.loads(checked_path(source["configuration"], args.data_root).read_text())
    thresholds, _ = read_thresholds(checked_path(source["thresholds"], args.data_root), safe_config["configuration_id"], args.data_root)
    if thresholds["threshold_id"] != source["threshold_id"]:
        raise ValueError("candidate safety threshold identity differs")
    safe_rows = pd.read_csv(checked_path(safety["outputs"]["qualified_safety"], args.data_root), sep="\t")
    for cell in ("HepG2", "K562"):
        comparison_id = config["references"][cell]["inputs"]["comparison_result_id"]
        validate_safety_link(safety, source, safe_config, {"result_id": comparison_id}, cell=cell)
        if (safe_config["cells"][cell]["bias_model_sha256"] != config["settings"]["bias_model_sha256"]
                or safe_config["cells"][cell]["dispersion"] != config["settings"]["dispersion"][cell]
                or not safe_rows.loc[safe_rows.cell.eq(cell) & safe_rows.tf.eq("CTCF"), "dispersion"].eq(config["settings"]["dispersion"][cell]).all()):
            raise ValueError("candidate safety uses a different model or dispersion")
    tables = {key: pd.read_csv(checked_path(result["outputs"][key], args.data_root), sep="\t")
              for key in ("metrics", "shapes", "bootstrap", "coverage", "profiles")}
    summary = select_evidence(tables["metrics"], tables["shapes"], tables["bootstrap"], tables["coverage"], safe_rows)
    if args.outdir.exists() and any(args.outdir.iterdir()):
        raise ValueError("use a new report directory; historical reports are immutable")
    args.outdir.mkdir(parents=True, exist_ok=True)
    pdf, png, metrics, readme = [args.outdir / name for name in ("CTCF_direct_bias_stability.pdf", "CTCF_direct_bias_stability.png", "metrics.tsv", "README.md")]
    summary.to_csv(metrics, sep="\t", index=False)
    render(summary, tables["profiles"], pdf, png)
    readme.write_text(
        "# CTCF direct-bias development stability\n\n"
        "This report compares sample/depth-matched DWM expected cuts with the frozen CPU parametric expectation using identical total-matched signed-deviance scoring and fixed geometry. "
        "It is not a direct production corrected-track comparison or a new detector. HepG2 dispersion is 0; K562 is 1.6144318781609228. "
        "They remain separate recorded experiments. The candidate and parameters were not retuned.\n\n"
        "The frozen matrix contains five seeds (2026–2030) at 10M/25M/50M in replicate 1, plus three full-depth biological replicates per cell. "
        "The 36 point estimates are not 36 independent biological replicates. Replicates 2 and 3 still lack the depth/seed factorial panel. "
        "All loci are previously examined chr19–22/X; new SK-N-SH/GM23338 labels and final naked-DNA replicate 3 remain unopened.\n\n"
        "Original finite support contained zero-filled DWM expectations, including 65 HepG2 and 66 K562 CTCF sites. "
        "The new fixed support requires positive expectation in every run. Existing 10M motif/accessibility/observed-total matches retain complete pairs without rematching: 577 HepG2 and 534 K562 pairs. "
        "Zero observed-count sites remain. Current cut-total SMD exceeds 0.1 in 11/36 matched CTCF runs, reaching 0.134; motif and reference-accessibility SMDs remain below 0.012. "
        "These are conditional sensitivity results, not causal removal of accessibility confounding.\n\n"
        "All 36 CTCF point estimates improve against raw, fixed DWM and sample-matched DWM and retain canonical protection. "
        "All six full-replicate AUC and relative AP gains have positive 1,000-draw chromosome-block intervals. "
        "Only the 10M seed-2026 reference and full replicates have metric intervals; other seed/depth rows contain point estimates. "
        "K562's matched 10M reference interval against raw still includes zero. Five chromosomes provide only five independent blocks; intervals are not multiplicity-adjusted or calibrated simultaneous functional bands. "
        "Plots show replicate curves and their mean, not a pooled confidence band. Stronger shoulders contribute; center depletion does not increase uniformly.\n\n"
        "The exact unchanged direct operator passes the existing naked-DNA replicate-2 candidate screen: HepG2 0/121 paired informative calls (Wilson upper 3.08%); K562 1/114 (4.80%, increase over DWM 0.88 percentage points). "
        "One additional informative site per cell lacks a DWM expectation and remains in the passing absolute-candidate calculation (HepG2 0/122; K562 1/115). "
        "Sparse control counts and absence of a large residual flag do not prove no occupancy leakage. This is not final untouched safety or depth-stratified null validation.\n\n"
        "Only CTCF has canonical protection among powered tasks. Full-depth mean relative AP gains remain below 10%, and no second difficult TF family passes. "
        "CPU runtime, pipeline peak memory and calibrated uncertainty remain unestablished. No method is promoted; production defaults remain unchanged. "
        "Historical reports are preserved; this report explicitly adds expected-track coverage and replicate/depth evidence.\n\n"
        f"Stability result: `{result['stage_id']}`. Configuration: `{config['configuration_id']}`. Calculation code: `c453a37`. "
        f"Safety coverage: `{safety['coverage_id']}`. Bias NPZ: `{config['settings']['bias_model_sha256']}`.\n",
        encoding="utf-8")
    sources = {"stability": {"file": args.stability_manifest.name, "sha256": file_sha256(args.stability_manifest)},
               "configuration": {"file": config_path.name, "sha256": file_sha256(config_path)},
               "candidate_safety": {"file": args.safety_manifest.name, "sha256": file_sha256(args.safety_manifest)}}
    manifest = {"schema": "fp-tools-direct-bias-stability-report-v1", "stability_id": result["stage_id"],
                "configuration_id": config["configuration_id"], "safety_coverage_id": safety["coverage_id"],
                "sources": sources, "renderer_sha256": file_sha256(Path(__file__)), "promotion_eligible": False,
                "outputs": {path.name: {"bytes": path.stat().st_size, "sha256": file_sha256(path)} for path in (pdf, png, metrics, readme)}}
    manifest["report_id"] = canonical_id(manifest)
    write_immutable(args.outdir / "manifest.json", manifest)
    (args.outdir / "SHA256SUMS").write_text("".join(f"{file_sha256(path)}  {path.name}\n" for path in sorted(args.outdir.iterdir()) if path.is_file()))
    print(f"Rendered safety-linked CTCF stability report {manifest['report_id']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
