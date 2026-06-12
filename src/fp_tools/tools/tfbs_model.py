#!/usr/bin/env python
"""Motif-centric supervised TFBS model helpers and CLIs."""

from __future__ import annotations

import argparse
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

MODEL_VERSION = "fp-tools-tfbs-tabular-v1"


def infer_feature_columns(frame: pd.DataFrame, label_column: str | None = None, id_columns: list[str] | None = None) -> list[str]:
    """Infer numeric feature columns from a table."""

    excluded = set(id_columns or [])
    if label_column is not None:
        excluded.add(label_column)
    features = [
        column for column in frame.columns
        if column not in excluded and pd.api.types.is_numeric_dtype(frame[column])
    ]
    if not features:
        raise ValueError("No numeric feature columns were found. Use --feature-columns to specify them explicitly.")
    return features


def _read_table(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def _write_metrics(path: str | Path, metrics: dict[str, float | int | str]) -> None:
    with Path(path).open("w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        for key, value in metrics.items():
            handle.write(f"{key}\t{value}\n")


def _binary_metrics(y_true: np.ndarray, probabilities: np.ndarray) -> dict[str, float | int]:
    metrics: dict[str, float | int] = {
        "n": int(len(y_true)),
        "positives": int(np.sum(y_true == 1)),
        "negatives": int(np.sum(y_true == 0)),
    }
    if len(np.unique(y_true)) == 2:
        metrics["auroc"] = float(roc_auc_score(y_true, probabilities))
        metrics["auprc"] = float(average_precision_score(y_true, probabilities))
        metrics["brier"] = float(brier_score_loss(y_true, probabilities))
    return metrics


def train_tabular_model(
    train_table: str | Path,
    model_out: str | Path,
    label_column: str = "label",
    feature_columns: list[str] | None = None,
    metrics_out: str | Path | None = None,
    seed: int = 2026,
) -> dict[str, float | int | str]:
    """Train a compact motif-centric tabular classifier and save it as pickle."""

    frame = _read_table(train_table)
    if label_column not in frame.columns:
        raise ValueError(f"Label column '{label_column}' was not found in {train_table}")
    if feature_columns is None:
        feature_columns = infer_feature_columns(frame, label_column=label_column)

    y = frame[label_column].astype(int).to_numpy()
    if len(np.unique(y)) != 2:
        raise ValueError("Training labels must contain both positive and negative examples.")

    model = Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            ("classifier", LogisticRegression(max_iter=1000, class_weight="balanced", random_state=seed)),
        ]
    )
    model.fit(frame[feature_columns], y)
    probabilities = model.predict_proba(frame[feature_columns])[:, 1]
    metrics: dict[str, float | int | str] = {
        "model_version": MODEL_VERSION,
        "feature_columns": ",".join(feature_columns),
        **_binary_metrics(y, probabilities),
    }

    model_path = Path(model_out)
    model_path.parent.mkdir(parents=True, exist_ok=True)
    with model_path.open("wb") as handle:
        pickle.dump(
            {
                "version": MODEL_VERSION,
                "label_column": label_column,
                "feature_columns": feature_columns,
                "model": model,
            },
            handle,
        )

    if metrics_out is not None:
        Path(metrics_out).parent.mkdir(parents=True, exist_ok=True)
        _write_metrics(metrics_out, metrics)
    return metrics


def predict_tabular_model(
    model_path: str | Path,
    feature_table: str | Path,
    output: str | Path,
    probability_column: str = "binding_probability",
) -> pd.DataFrame:
    """Apply a saved motif-centric tabular classifier to a feature table."""

    with Path(model_path).open("rb") as handle:
        bundle = pickle.load(handle)
    if bundle.get("version") != MODEL_VERSION:
        raise ValueError(f"Unsupported model version: {bundle.get('version')}")

    frame = _read_table(feature_table)
    feature_columns = bundle["feature_columns"]
    missing = [column for column in feature_columns if column not in frame.columns]
    if missing:
        raise ValueError(f"Missing feature columns in {feature_table}: {missing}")

    probabilities = bundle["model"].predict_proba(frame[feature_columns])[:, 1]
    out_frame = frame.copy()
    out_frame[probability_column] = probabilities
    out_path = Path(output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_frame.to_csv(out_path, sep="\t", index=False)
    return out_frame


def train_main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Train a motif-centric supervised TFBS model from a feature TSV.")
    parser.add_argument("--train-table", required=True, help="Training feature TSV with a binary label column.")
    parser.add_argument("--model-out", required=True, help="Output pickle model path.")
    parser.add_argument("--label-column", default="label", help="Binary label column name (default: label).")
    parser.add_argument("--feature-columns", nargs="*", default=None, help="Feature columns to use. Defaults to all numeric non-label columns.")
    parser.add_argument("--metrics-out", default=None, help="Optional metrics TSV output path.")
    parser.add_argument("--seed", type=int, default=2026, help="Random seed (default: 2026).")
    args = parser.parse_args(argv)

    metrics = train_tabular_model(
        args.train_table,
        args.model_out,
        label_column=args.label_column,
        feature_columns=args.feature_columns,
        metrics_out=args.metrics_out,
        seed=args.seed,
    )
    for key, value in metrics.items():
        print(f"{key}\t{value}")
    return 0


def predict_main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Predict TFBS binding probabilities from a saved fp-tools tabular model.")
    parser.add_argument("--model", required=True, help="Model pickle from the TFBS model training utility.")
    parser.add_argument("--features", required=True, help="Feature TSV to score.")
    parser.add_argument("--out", required=True, help="Output TSV with appended binding probability column.")
    parser.add_argument("--probability-column", default="binding_probability", help="Output probability column name.")
    args = parser.parse_args(argv)

    predict_tabular_model(args.model, args.features, args.out, probability_column=args.probability_column)
    print(args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(train_main())
