#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import (
    r2_score,
    mean_absolute_error,
    classification_report,
    confusion_matrix,
    ConfusionMatrixDisplay,
)


def first_existing_column(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None


def normalize_tier_label(x):
    x = str(x).strip()
    mapping = {
        "Elite": "Elite candidate",
        "High": "High candidate",
        "Moderate": "Moderate candidate",
        "Low": "Low priority",
    }
    x = mapping.get(x, x)
    x = x.replace("Cautionary ", "")
    return x


def main():
    ap = argparse.ArgumentParser(
        description="Functional-module-only Random Forest interpretation for LAB-Score."
    )
    ap.add_argument("--scores", required=True, help="Input LAB_score_v1.tsv file")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--threads", type=int, default=-1)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    scores = Path(args.scores)
    outdir = Path(args.outdir)
    figdir = outdir / "figures"

    outdir.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(scores, sep="\t")

    feature_cols = [
        "CAZyme_total",
        "CAZyme_total.1",
        "acid_energy.acid_energy",
        "bileresistance.bileresistance",
        "gutpersistence.gutpersistence",
        "osmoticstress.osmoticstress",
        "heatstress.heatstress",
        "coldstress.coldstress",
        "gaba.gaba",
        "vitamins.vitamins",
        "immunomodulation.immunomodulation",
        "bacteriocins.bacteriocin",
        "cellenvelope_eps.cellenvelope_eps",
        "adhesion_surface.adhesion_surface",
        "adhesion_biofilm.adhesion_biofilm",
        "carbohydrate.carbohydrate",
        "metabolism.metabolism",
        "defense_crispr.defense_crispr",
        "antipath_qs.antipath_qs",
        "alkalinestress.alkalinestress",
    ]

    # avoid duplicate CAZyme column if both exist
    if "CAZyme_total" in df.columns and "CAZyme_total.1" in df.columns:
        feature_cols.remove("CAZyme_total.1")

    feature_cols = [c for c in feature_cols if c in df.columns]

    if "LAB_score_v1" not in df.columns:
        raise ValueError("Missing required column: LAB_score_v1")

    if len(feature_cols) < 3:
        raise ValueError(f"Too few functional features found: {feature_cols}")

    X = df[feature_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    y_score = pd.to_numeric(df["LAB_score_v1"], errors="coerce").fillna(0)

    # ------------------------------------------------------------
    # Regression model
    # ------------------------------------------------------------
    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y_score,
        test_size=0.2,
        random_state=args.seed,
    )

    rf_reg = RandomForestRegressor(
        n_estimators=500,
        random_state=args.seed,
        n_jobs=args.threads,
        max_features="sqrt",
    )

    rf_reg.fit(X_train, y_train)
    pred = rf_reg.predict(X_test)

    reg_r2 = r2_score(y_test, pred)
    reg_mae = mean_absolute_error(y_test, pred)

    reg_metrics = pd.DataFrame(
        {
            "metric": ["R2", "MAE"],
            "value": [reg_r2, reg_mae],
        }
    )

    reg_metrics.to_csv(
        outdir / "RF_functional_only_regression_metrics.tsv",
        sep="\t",
        index=False,
    )

    reg_imp = pd.DataFrame(
        {
            "feature": feature_cols,
            "importance": rf_reg.feature_importances_,
        }
    ).sort_values("importance", ascending=False)

    reg_imp.to_csv(
        outdir / "RF_functional_only_regression_feature_importance.tsv",
        sep="\t",
        index=False,
    )

    # ------------------------------------------------------------
    # Classification model
    # ------------------------------------------------------------
    tier_col = first_existing_column(
        df,
        ["Candidate_tier_v1_1", "Candidate_tier", "Priority_class"],
    )

    if tier_col is None:
        raise ValueError("Missing candidate-tier column.")

    clf_df = df.copy()
    clf_df["Candidate_group_simple"] = clf_df[tier_col].map(normalize_tier_label)

    # Remove critical review from biological candidate-group classification
    clf_df = clf_df[clf_df["Candidate_group_simple"] != "Critical safety review"].copy()

    # Keep classes with enough samples
    class_counts = clf_df["Candidate_group_simple"].value_counts()
    valid_classes = class_counts[class_counts >= 5].index.tolist()
    clf_df = clf_df[clf_df["Candidate_group_simple"].isin(valid_classes)].copy()

    Xc = clf_df[feature_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    yc = clf_df["Candidate_group_simple"]

    if yc.nunique() < 2:
        raise ValueError(f"Too few classes for classification: {yc.value_counts().to_dict()}")

    Xc_train, Xc_test, yc_train, yc_test = train_test_split(
        Xc,
        yc,
        test_size=0.2,
        random_state=args.seed,
        stratify=yc,
    )

    rf_clf = RandomForestClassifier(
        n_estimators=500,
        random_state=args.seed,
        n_jobs=args.threads,
        class_weight="balanced",
        max_features="sqrt",
    )

    rf_clf.fit(Xc_train, yc_train)
    yc_pred = rf_clf.predict(Xc_test)

    report = pd.DataFrame(
        classification_report(yc_test, yc_pred, output_dict=True)
    ).T

    report.to_csv(
        outdir / "RF_functional_only_candidate_group_classification_report.tsv",
        sep="\t",
    )

    clf_acc = float(report.loc["accuracy", "precision"])
    weighted_f1 = float(report.loc["weighted avg", "f1-score"])

    model_summary = pd.DataFrame(
        [
            {
                "Model": "Functional modules only",
                "Regression_R2": reg_r2,
                "Regression_MAE": reg_mae,
                "Classification_accuracy": clf_acc,
                "Weighted_F1": weighted_f1,
                "n_genomes": len(df),
                "n_classification_genomes": len(clf_df),
                "n_features": len(feature_cols),
            }
        ]
    )

    model_summary.to_csv(
        outdir / "RF_functional_only_model_summary.tsv",
        sep="\t",
        index=False,
    )

    clf_imp = pd.DataFrame(
        {
            "feature": feature_cols,
            "importance": rf_clf.feature_importances_,
        }
    ).sort_values("importance", ascending=False)

    clf_imp.to_csv(
        outdir / "RF_functional_only_candidate_group_feature_importance.tsv",
        sep="\t",
        index=False,
    )

    cm = confusion_matrix(yc_test, yc_pred, labels=rf_clf.classes_)
    cm_df = pd.DataFrame(cm, index=rf_clf.classes_, columns=rf_clf.classes_)
    cm_df.to_csv(
        outdir / "RF_functional_only_candidate_group_confusion_matrix.tsv",
        sep="\t",
    )

    # ------------------------------------------------------------
    # Figures
    # ------------------------------------------------------------
    top = clf_imp.head(20).sort_values("importance")

    plt.figure(figsize=(8, 7))
    plt.barh(top["feature"], top["importance"])
    plt.xlabel("Random Forest importance")
    plt.ylabel("Functional module")
    plt.title("Functional modules driving LAB-Score candidate groups")
    plt.tight_layout()
    plt.savefig(figdir / "Figure_functional_only_candidate_group_feature_importance.png", dpi=300)
    plt.savefig(figdir / "Figure_functional_only_candidate_group_feature_importance.pdf")
    plt.close()

    disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=rf_clf.classes_)
    disp.plot(values_format="d", xticks_rotation=45)
    plt.title("Functional-only model: LAB-Score candidate groups")
    plt.tight_layout()
    plt.savefig(figdir / "Figure_functional_only_candidate_group_confusion_matrix.png", dpi=300)
    plt.savefig(figdir / "Figure_functional_only_candidate_group_confusion_matrix.pdf")
    plt.close()

    print("\nSaved functional-only ML outputs to:")
    print(outdir)

    print("\nModel summary:")
    print(model_summary.to_string(index=False))

    print("\nTop 20 functional modules:")
    print(clf_imp.head(20).to_string(index=False))


if __name__ == "__main__":
    main()
