from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parents[1] / "benchmarks" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from evaluate_differential_functional_profiles import (  # noqa: E402
    benjamini_hochberg,
    evaluate_differential_profiles,
)


def test_bh_adjustment_is_monotone_in_rank() -> None:
    pvalues = np.asarray([0.04, 0.001, 0.02, np.nan])
    adjusted = benjamini_hochberg(pvalues)
    order = np.argsort(pvalues[:3])
    assert np.all(np.diff(adjusted[:3][order]) >= 0)
    assert np.isnan(adjusted[3])


def test_differential_profile_manifest_runs_replicate_level_test(tmp_path: Path) -> None:
    rng = np.random.default_rng(14)
    flank = 20
    width = flank * 2 + 1
    x = np.arange(-flank, flank + 1)
    shape = -1.2 * np.exp(-0.5 * np.square(x / 4.0))
    rows = []
    for condition in ("stress", "control"):
        for replicate in ("r1", "r2", "r3"):
            sites = pd.DataFrame(
                {
                    "TFBS_chr": ["chr1"] * 80,
                    "TFBS_start": np.arange(80) * 100 + 1000,
                    "TFBS_end": np.arange(80) * 100 + 1010,
                    "TFBS_strand": ["+"] * 80,
                    "tf": ["TFX"] * 80,
                }
            )
            profile = rng.normal(scale=0.4, size=(80, width))
            if condition == "stress":
                profile += shape
            npz = tmp_path / f"{condition}.{replicate}.npz"
            np.savez_compressed(npz, combined_residual=profile, valid=np.ones(80, dtype=bool))
            site_path = tmp_path / f"{condition}.{replicate}.sites.tsv"
            sites.to_csv(site_path, sep="\t", index=False)
            rows.append(
                {
                    "sample": f"{condition}_{replicate}",
                    "condition": condition,
                    "replicate": replicate,
                    "profiles_npz": npz,
                    "sites_tsv": site_path,
                }
            )
    summary, curves = evaluate_differential_profiles(
        pd.DataFrame(rows),
        [("stress", "control")],
        channel="combined",
        flank=flank,
        bootstraps=100,
        seed=3,
    )
    assert len(summary) == 1
    assert len(curves) == width
    assert summary.loc[0, "change_depletion"] > 0.5
    center = curves.loc[curves["position"] == 0, "difference"].iloc[0]
    assert center < -0.5
