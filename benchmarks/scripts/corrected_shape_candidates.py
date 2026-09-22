"""Training-only, interpretable detector candidates for corrected research data.

This module does not select on tuning/test outcomes or estimate probabilities.
Template direction is a development screen, not calibrated shape uncertainty.
"""

from pathlib import Path
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from fp_tools.tools.functional_footprints import FunctionalTemplateDetector, normalize_functional_profiles

POSITIONS = np.arange(-100, 101)
TRAIN = {f"chr{i}" for i in range(1, 16)}
# Conventional first; exact ties therefore retain the conventional geometry.
GEOMETRIES = [(7, 12, 40)] + [(c, c + 5, e) for c in (3, 5, 10, 15, 20) for e in (40, 60)]


def geometry_weights(geometry):
    center, shoulder_start, shoulder_end = geometry
    if not 0 < center < shoulder_start < shoulder_end <= 100:
        raise ValueError("invalid footprint geometry")
    inner = np.abs(POSITIONS) <= center
    outer = (np.abs(POSITIONS) >= shoulder_start) & (np.abs(POSITIONS) <= shoulder_end)
    return outer.astype(float) / outer.sum() - inner.astype(float) / inner.sum()


def protection_direction(template):
    values = np.asarray(template, dtype=float)
    if values.shape != POSITIONS.shape or not np.isfinite(values).all():
        raise ValueError("template must contain 201 finite positions")
    center = values[np.abs(POSITIONS) <= 5].mean()
    left = values[(POSITIONS >= -40) & (POSITIONS <= -15)].mean()
    right = values[(POSITIONS >= 15) & (POSITIONS <= 40)].mean()
    return bool(center < 0 and left > center and right > center)


def fit_candidates(sites, profiles, *, minimum_per_class=200, prior_strength=250.0):
    """Fit each task on training chromosomes; preserve failed task records.

    Input sites must already use the same frozen matching/coverage policy for
    all correction arms. Family priors average other eligible task contrasts
    with equal task weights. Global priors are used only without family peers.
    """
    if (sites.empty or not set(sites.TFBS_chr).issubset(TRAIN)
            or not set(sites.cell).issubset({"HepG2", "K562"})
            or not sites.chip_label.isin([0, 1]).all()):
        raise ValueError("candidate fitting accepts training chromosomes and development cells only")
    values = np.asarray(profiles, dtype=float)
    if values.shape != (len(sites), 201) or not np.isfinite(values).all():
        raise ValueError("training profiles must match finite site rows")
    if minimum_per_class < 2 or prior_strength < 0:
        raise ValueError("invalid training support or prior strength")
    sites = sites.reset_index(drop=True)
    normalized = normalize_functional_profiles(values, POSITIONS)
    tasks, contrasts = {}, {}
    for key, task in sites.groupby(["cell", "tf"], sort=True):
        ix = task.index.to_numpy()
        y = task.chip_label.to_numpy(dtype=int)
        families = task.motif_family.unique()
        if len(families) != 1:
            raise ValueError("task has ambiguous motif family assignment")
        tasks[key] = dict(indexes=ix, labels=y, family=str(families[0]))
        if min(int((y == 0).sum()), int((y == 1).sum())) >= minimum_per_class:
            contrasts[key] = normalized[ix][y == 1].mean(axis=0) - normalized[ix][y == 0].mean(axis=0)
    models, rows = {}, []
    for key, task in tasks.items():
        ix, y = task["indexes"], task["labels"]
        row = dict(cell=key[0], tf=key[1], motif_family=task["family"],
                   train_positive=int((y == 1).sum()), train_negative=int((y == 0).sum()))
        if key not in contrasts:
            rows.append(dict(row, status="underpowered", geometry=None, spline_available=False))
            continue
        if not protection_direction(contrasts[key]):
            rows.append(dict(row, status="no_training_central_protection", geometry=None, spline_available=False))
            continue
        scores = [roc_auc_score(y, values[ix] @ geometry_weights(g)) for g in GEOMETRIES]
        geometry = GEOMETRIES[int(np.argmax(scores))]
        # Exclude the same TF in both contexts: its other-cell data must not
        # masquerade as independent family support.
        others = [k for k in contrasts if k[1] != key[1]]
        family = [k for k in others if tasks[k]["family"] == task["family"]]
        peers = family or others
        prior = np.mean([contrasts[k] for k in peers], axis=0) if peers else None
        detector = FunctionalTemplateDetector(POSITIONS, smoother="spline", spline_penalty=10.0,
                                               window_limit=50.0, variance_shrinkage=0.5, variance_floor=0.05)
        try:
            detector.fit(values[ix], y, prior_template=prior, prior_strength=prior_strength if peers else 0.0)
        except ValueError as exc:
            if "non-zero functional template" not in str(exc):
                raise
            detector = None
        spline_ok = detector is not None and protection_direction(detector.footprint_template_)
        models[key] = dict(geometry=geometry, spline=detector if spline_ok else None)
        rows.append(dict(row, status="trained", geometry=list(geometry), spline_available=spline_ok,
                         spline_status="trained" if spline_ok else "no_protective_spline_template",
                         prior_scope="family" if family else "global" if peers else "none",
                         prior_tasks=[list(k) for k in peers], train_geometry_auroc=max(scores)))
    return models, pd.DataFrame(rows)


def score_candidates(model, profiles):
    values = np.asarray(profiles, dtype=float)
    if values.ndim != 2 or values.shape[1] != 201 or not np.isfinite(values).all():
        raise ValueError("candidate scoring requires finite 201-position profiles")
    scores = {"conventional": values @ geometry_weights(GEOMETRIES[0])}
    if model is not None:
        scores["trained_geometry"] = values @ geometry_weights(model["geometry"])
        if model["spline"] is not None:
            scores["hierarchical_spline"] = model["spline"].decision_function(values)
    return scores
