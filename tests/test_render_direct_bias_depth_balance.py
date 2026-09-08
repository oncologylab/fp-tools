from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks/scripts"))
import render_direct_bias_depth_balance as m


def test_all_tf_figures_include_unestimable_tasks_and_both_evidence_roles(tmp_path):
    metrics, profiles, retention = [], [], []
    for cell, tf in (("HepG2", "CTCF"), ("K562", "CTCF"), ("K562", "SPARSE")):
        for depth in ("10m", "25m", "50m", "full"):
            run = {"cell": cell, "tf": tf, "depth": depth, "sample": cell + "_rep1", "seed": 2026}
            retention.append({**run, "original_pairs": 4 if tf == "CTCF" else 0,
                              "retained_pairs": 3 if tf == "CTCF" else 0})
            if tf == "SPARSE":
                continue
            for cohort in ("motif_accessibility_observed_total", m.audit.SECONDARY):
                metrics.append({**run, "cohort": cohort, "method": m.audit.base.PARAMETRIC, "auroc_gain_over_raw": .01})
                for method in (m.audit.base.DWM, m.audit.base.SAMPLE_DWM, m.audit.base.PARAMETRIC):
                    for position in range(-100, 101):
                        profiles.append({**run, "cohort": cohort, "method": method, "position": position,
                                         "positive_minus_negative": -.1 * np.exp(-(position / 10) ** 2)})
    outputs, pages = m.render(pd.DataFrame(metrics), pd.DataFrame(profiles), pd.DataFrame(retention), tmp_path)
    assert pages == [{"cell": "HepG2", "tf": "CTCF"}, {"cell": "K562", "tf": "CTCF"}, {"cell": "K562", "tf": "SPARSE"}]
    assert set(outputs) == {"current_depth_balance_all_tf.pdf", "HepG2_CTCF_current_depth_balance.png", "K562_CTCF_current_depth_balance.png"}
    pdf = (tmp_path / "current_depth_balance_all_tf.pdf").read_bytes()
    assert pdf.startswith(b"%PDF") and len(re.findall(rb"/Type\s*/Page\b", pdf)) == 3
    for entry in outputs.values():
        assert m.record(Path(entry["path"])) == entry
