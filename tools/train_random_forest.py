#!/usr/bin/env python3
"""Train a random forest model from MashScreen CSV output."""

from __future__ import annotations

import argparse
import csv
import pickle
from pathlib import Path

import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, confusion_matrix, precision_recall_fscore_support
from sklearn.model_selection import GridSearchCV, StratifiedKFold, cross_val_predict


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Train a random forest serotype classifier using a MashScreen matrix CSV. "
            "The first column must contain labels in the format '<class>-<sample_id>'."
        )
    )
    parser.add_argument(
        "--input-csv",
        default="/workspace/data/VRD_Jonathan/training_genomes_mashscreen.csv",
        help="Path to the MashScreen training CSV.",
    )
    parser.add_argument(
        "--model-out",
        default="model/model.rfm",
        help="Path where the trained model pickle will be written.",
    )
    parser.add_argument(
        "--max-cv-folds",
        type=int,
        default=5,
        help="Maximum number of stratified cross-validation folds.",
    )
    parser.add_argument(
        "--random-state",
        type=int,
        default=42,
        help="Random seed for reproducible training.",
    )
    return parser.parse_args()


def load_training_data(csv_path: Path) -> tuple[np.ndarray, np.ndarray, list[str], list[str]]:
    with csv_path.open("r", newline="") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError("Training CSV is empty.") from exc

        if len(header) < 2:
            raise ValueError("Training CSV must have at least two columns (label + features).")

        all_feature_cols = [str(c) for c in header[1:]]
        rows = [row for row in reader if row]

    if not rows:
        raise ValueError("Training CSV has no data rows.")

    labels_all = [str(row[0]).split("-", 1)[0] for row in rows]
    label_counts: dict[str, int] = {}
    for label in labels_all:
        label_counts[label] = label_counts.get(label, 0) + 1

    excluded_zero_sample = [c for c in all_feature_cols if label_counts.get(c, 0) == 0]
    retained_classes = [c for c in all_feature_cols if c not in excluded_zero_sample]
    retained_index = {cls: idx for idx, cls in enumerate(all_feature_cols)}

    X_rows: list[list[float]] = []
    y_list: list[str] = []
    dropped_row_count = 0
    for row in rows:
        label = str(row[0]).split("-", 1)[0]
        if label not in retained_classes:
            dropped_row_count += 1
            continue

        features = []
        for cls in retained_classes:
            col_idx = retained_index[cls] + 1  # +1 because row[0] is label column
            try:
                features.append(float(row[col_idx]))
            except (IndexError, ValueError) as exc:
                raise ValueError(
                    f"Invalid feature value for class column '{cls}' in row: {row[:2]}"
                ) from exc

        X_rows.append(features)
        y_list.append(label)

    if dropped_row_count > 0:
        print(
            "Dropped rows with labels not represented in retained class columns: "
            f"{dropped_row_count}"
        )

    X = np.asarray(X_rows, dtype=float)
    y = np.asarray(y_list, dtype=str)

    if X.shape[0] == 0:
        raise ValueError("No training rows remain after filtering classes.")
    if len(np.unique(y)) < 2:
        raise ValueError("At least two classes are required to train a classifier.")

    return X, y, excluded_zero_sample, retained_classes


def choose_cv_folds(y: np.ndarray, max_cv_folds: int) -> int:
    classes, counts = np.unique(y, return_counts=True)
    min_count = int(counts.min())

    if min_count < 2:
        problematic = [str(c) for c, n in zip(classes, counts) if n < 2]
        raise ValueError(
            "Cannot run stratified cross-validation while ensuring every class appears in "
            f"every test fold. Classes with <2 samples: {problematic}"
        )

    cv_folds = min(max_cv_folds, min_count)
    if cv_folds < 2:
        raise ValueError("Need at least 2 folds for cross-validation.")
    return cv_folds


def per_class_accuracy(y_true: np.ndarray, y_pred: np.ndarray, classes: list[str]) -> dict[str, float]:
    cm = confusion_matrix(y_true, y_pred, labels=classes)
    total = cm.sum()
    per_class_acc: dict[str, float] = {}
    for idx, cls in enumerate(classes):
        tp = cm[idx, idx]
        fn = cm[idx, :].sum() - tp
        fp = cm[:, idx].sum() - tp
        tn = total - tp - fn - fp
        per_class_acc[cls] = (tp + tn) / total if total else 0.0
    return per_class_acc


def train_model(
    X: np.ndarray,
    y: np.ndarray,
    cv_folds: int,
    random_state: int,
) -> tuple[GridSearchCV, np.ndarray]:
    base_model = RandomForestClassifier(
        n_estimators=1000,
        random_state=random_state,
        n_jobs=-1,
    )

    param_grid = {
        "n_estimators": [1000, 1200],
        "max_depth": [None, 30],
        "min_samples_split": [2, 4],
        "min_samples_leaf": [1, 2],
        "max_features": ["sqrt", "log2"],
    }

    splitter = StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=random_state)
    grid = GridSearchCV(
        estimator=base_model,
        param_grid=param_grid,
        scoring="accuracy",
        cv=splitter,
        n_jobs=-1,
        verbose=1,
        refit=True,
    )
    grid.fit(X, y)

    y_pred_cv = cross_val_predict(
        grid.best_estimator_,
        X,
        y,
        cv=splitter,
        n_jobs=-1,
    )

    return grid, y_pred_cv


def main() -> None:
    args = parse_args()

    input_csv = Path(args.input_csv)
    if not input_csv.exists():
        raise FileNotFoundError(f"Input CSV not found: {input_csv}")

    model_out = Path(args.model_out)
    if not model_out.is_absolute():
        model_out = (Path(__file__).resolve().parents[1] / model_out).resolve()

    X, y, excluded_zero_sample, retained_classes = load_training_data(input_csv)

    print(f"Loaded training rows: {X.shape[0]}")
    print(f"Loaded features: {X.shape[1]}")
    print(f"Classes retained: {len(retained_classes)}")
    if excluded_zero_sample:
        print(
            "Excluded classes with zero samples: "
            + ", ".join(sorted(excluded_zero_sample))
        )
    else:
        print("Excluded classes with zero samples: none")

    cv_folds = choose_cv_folds(y, args.max_cv_folds)
    print(f"Using stratified {cv_folds}-fold cross-validation")

    grid, y_pred_cv = train_model(X, y, cv_folds, args.random_state)

    classes = sorted(np.unique(y).tolist())
    precision, recall, _, _ = precision_recall_fscore_support(
        y,
        y_pred_cv,
        labels=classes,
        zero_division=0,
    )
    class_acc = per_class_accuracy(y, y_pred_cv, classes)
    overall_acc = accuracy_score(y, y_pred_cv)

    print("Best hyperparameters:")
    for k, v in grid.best_params_.items():
        print(f"  {k}: {v}")
    print(f"Best mean CV accuracy (grid score): {grid.best_score_:.6f}")
    print(f"Overall CV accuracy (cross_val_predict): {overall_acc:.6f}")

    print("Per-class metrics (from cross-validated predictions):")
    print("class,accuracy,precision,recall")
    for idx, cls in enumerate(classes):
        print(f"{cls},{class_acc[cls]:.6f},{precision[idx]:.6f},{recall[idx]:.6f}")

    model_out.parent.mkdir(parents=True, exist_ok=True)
    with model_out.open("wb") as handle:
        pickle.dump(grid.best_estimator_, handle)

    print(f"Saved best model to: {model_out}")


if __name__ == "__main__":
    main()