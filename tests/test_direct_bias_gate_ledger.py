from pathlib import Path
import sys
import csv
import gzip
import json

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks/scripts"))
import build_direct_bias_gate_ledger as ledger


def fixture():
    tables = {n: [] for n in ("metrics", "coverage", "shapes", "bootstrap")}
    for k in sorted(ledger.expected_keys()):
        row = {**dict(zip(ledger.KEYS, k)), "cohort": ledger.COHORT}
        small = row["tf"] == "HMGA1"
        counts = {"n_positive": 1 if small else 200, "n_negative": 1 if small else 200}
        tables["coverage"].append({**row, **counts, "status": "underpowered" if small else "eligible",
                                   "observed_total_smd": "" if small else "0.05"})
        if small:
            continue
        for method, value in ((ledger.CANDIDATE, .7), (ledger.BASELINE, .6)):
            tables["metrics"].append({**row, **counts, "method": method,
                                      "auroc": value, "auprc": value, "functional_separation": value,
                                      "auroc_gain_over_sample_dwm": .1,
                                      "relative_auprc_gain_over_sample_dwm": .7 / .6 - 1})
        tables["shapes"].append({**row, "method": ledger.CANDIDATE, "canonical_protection": "True",
                                 "center_contrast": -.1, "left_shoulder": .1,
                                 "right_shoulder": .1, "significant_center_points": 3})
        if row["depth"] == "full" or (row["depth"], row["seed"]) == ("10m", "2026"):
            tables["bootstrap"].append({**row, "method": ledger.CANDIDATE, "baseline": ledger.BASELINE,
                                        "auroc_gain_lower_95": .01, "relative_auprc_gain_lower_95": .01})
    return tables


def test_unscored_tasks_and_missing_safety_are_preserved():
    rows, exceptions = ledger.summarize(fixture())
    assert len(rows) == 21 and not exceptions
    small = next(r for r in rows if r["tf"] == "HMGA1")
    assert small["expected_runs"] == 48 and small["scored_runs"] == 0
    assert small["auroc_gain_range"] == [None, None]
    assert all(not r["promotion_eligible"] for r in rows)
    assert all(r["safety_status"] == "not_assessed_by_this_ledger" for r in rows)


def test_shape_failure_and_crossing_interval_remain_visible():
    tables = fixture()
    tables["shapes"][0]["canonical_protection"] = "False"
    tables["bootstrap"][0]["auroc_gain_lower_95"] = -.01
    rows, exceptions = ledger.summarize(tables)
    assert len(exceptions) == 1
    affected = next(r for r in rows if (r["cell"], r["tf"]) == ("HepG2", "ARID3A"))
    assert affected["canonical_protection_runs"] == 47
    assert affected["principal_intervals_both_positive"] == 5


@pytest.mark.parametrize("failure", ["task", "duplicate", "shape", "interval", "counts", "gain"])
def test_reject_incomplete_or_inconsistent_comparisons(failure):
    tables = fixture()
    if failure == "task":
        tables["coverage"] = [r for r in tables["coverage"] if r["tf"] != "HMGA1"]
    elif failure == "duplicate":
        tables["metrics"].append(tables["metrics"][0])
    elif failure == "shape":
        tables["shapes"].pop()
    elif failure == "interval":
        tables["bootstrap"].pop()
    elif failure == "counts":
        tables["metrics"][0]["n_positive"] = 201
    else:
        tables["metrics"][0]["auroc_gain_over_sample_dwm"] = .4
    with pytest.raises(ValueError):
        ledger.summarize(tables)


def test_changed_manifest_rejected_before_output_creation(tmp_path):
    path = tmp_path / "manifest.json"
    path.write_text("{}")
    with pytest.raises(ValueError, match="checksum"):
        ledger.build(path, "0" * 64, tmp_path / "output")
    assert not (tmp_path / "output").exists()


def recorded_fixture(tmp_path):
    outputs = {}
    for name, rows in fixture().items():
        path = tmp_path / f"{name}.tsv.gz"
        with gzip.open(path, "wt") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        outputs[name] = {"path": str(path), "sha256": ledger.digest(path)}
    document = {"schema": "fp-tools-direct-bias-stability-96-v1",
                "evidence": "post_hoc_development_complete_replicate_depth_stability",
                "stage_id": "fixture", "outputs": outputs,
                **{k: False for k in ("candidate_retuned", "new_holdout_labels_read",
                                      "final_untouched_control_read", "promotion_eligible")}}
    path = tmp_path / "manifest.json"
    path.write_text(json.dumps(document))
    return path, document


def test_exact_rerun_is_immutable_and_changed_output_is_rejected(tmp_path):
    path, _ = recorded_fixture(tmp_path)
    out = tmp_path / "output"
    result = ledger.build(path, ledger.digest(path), out)
    assert ledger.build(path, ledger.digest(path), out) == result
    (out / "gate_ledger.json").write_text("changed")
    with pytest.raises(ValueError, match="overwrite"):
        ledger.build(path, ledger.digest(path), out)


@pytest.mark.parametrize("mutation", ["table", "scope"])
def test_input_or_evidence_mutation_rejected(tmp_path, mutation):
    path, document = recorded_fixture(tmp_path)
    if mutation == "table":
        Path(document["outputs"]["metrics"]["path"]).write_bytes(b"changed")
    else:
        document["new_holdout_labels_read"] = True
        path.write_text(json.dumps(document))
    with pytest.raises(ValueError):
        ledger.build(path, ledger.digest(path), tmp_path / "output")
    assert not (tmp_path / "output").exists()
