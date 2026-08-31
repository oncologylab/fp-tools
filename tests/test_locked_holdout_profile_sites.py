import importlib.util
from pathlib import Path
import sys

import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "benchmarks" / "scripts" / "build_locked_holdout_profile_sites.py"
spec = importlib.util.spec_from_file_location("build_locked_holdout_profile_sites", SCRIPT)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)


def fixture_sites():
    rows = []
    for tf in ("A", "B"):
        for split, count in (("train", 8), ("validation", 3), ("test", 5)):
            for index in range(count):
                rows.append(
                    {
                        "cell": "Cell",
                        "tf": tf,
                        "motif_id": f"M{tf}",
                        "motif_family": tf,
                        "TFBS_chr": {"train": "chr1", "validation": "chr16", "test": "chr19"}[split],
                        "TFBS_start": index * 20,
                        "TFBS_end": index * 20 + 10,
                        "TFBS_strand": "+" if index % 2 else "-",
                        "motif_score": float(index),
                        "peak_start": index * 20,
                        "peak_end": index * 20 + 100,
                        "chromosome_split": split,
                    }
                )
    return pd.DataFrame(rows)


def test_selection_caps_train_and_retains_every_test_site_deterministically():
    frame = fixture_sites()
    first, counts = module.select_profile_sites(
        frame, maximum_train_per_tf=4, seed=2026
    )
    second, _ = module.select_profile_sites(
        frame, maximum_train_per_tf=4, seed=2026
    )
    pd.testing.assert_frame_equal(first, second)
    assert set(first["chromosome_split"]) == {"train", "test"}
    assert (counts["selected_train"] == 4).all()
    assert (counts["selected_test"] == 5).all()
    assert (counts["selected_validation"] == 0).all()


def test_validation_refuses_label_columns():
    frame = fixture_sites()
    frame["chip_label"] = 0
    with pytest.raises(ValueError, match="forbidden label columns"):
        module.validate_label_free(frame, Path("fixture.tsv"))


def test_route_filter_keeps_only_selected_bias_tasks(tmp_path):
    frame = fixture_sites()
    sites_path = tmp_path / "sites.tsv"
    frame.to_csv(sites_path, sep="\t", index=False)
    routes = pd.DataFrame(
        {
            "cell": ["Cell", "Cell"],
            "tf": ["A", "B"],
            "bias_configuration": ["NEW", "DWM"],
        }
    )
    route_path = tmp_path / "routes.tsv"
    routes.to_csv(route_path, sep="\t", index=False)
    output = tmp_path / "selected.tsv.gz"
    assert module.main(
        [
            "--sites", str(sites_path),
            "--out", str(output),
            "--counts-out", str(tmp_path / "counts.tsv"),
            "--manifest-out", str(tmp_path / "manifest.json"),
            "--routes", str(route_path),
            "--bias-configuration", "NEW",
            "--maximum-train-per-tf", "4",
        ]
    ) == 0
    filtered = pd.read_csv(output, sep="\t")
    assert set(filtered["tf"]) == {"A"}
    assert filtered["chromosome_split"].value_counts().to_dict() == {
        "test": 5,
        "train": 4,
    }
