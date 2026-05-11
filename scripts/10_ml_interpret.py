#!/usr/bin/env python3
"""
10_ml_interpret.py
==================
Machine learning interpretation of LAB-Score v1.1.
Handles small datasets gracefully (n < 10).
"""

import argparse, os, warnings
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import classification_report, confusion_matrix, r2_score, mean_absolute_error

warnings.filterwarnings("ignore")

ap = argparse.ArgumentParser()
ap.add_argument("--scores", required=True)
ap.add_argument("--outdir", required=True)
args = ap.parse_args()
os.makedirs(args.outdir, exist_ok=True)

# SHAP is intentionally not imported here. It is slow and is not used by this script.
SHAP = False

df = pd.read_csv(args.scores, sep="\t")
print(f"[ML] Loaded {len(df)} genomes")

# ── Minimum sample check ──────────────────────────────────────────────────────
MIN_SAMPLES = 10
if len(df) < MIN_SAMPLES:
    print(f"[ML] WARNING: Only {len(df)} genomes — ML requires at least {MIN_SAMPLES} samples.")
    print(f"[ML] Skipping ML models. Add more genomes to get meaningful ML results.")
    print(f"[ML] Writing placeholder outputs...")

    # Write empty performance summary
    pd.DataFrame([
        {"Model": "RF_Regression_FUNC", "Metric": "Note",
         "Value": f"Skipped — only {len(df)} genomes (need ≥{MIN_SAMPLES})"},
    ]).to_csv(f"{args.outdir}/model_performance_summary.tsv", sep="\t", index=False)

    # Write Spearman correlations (works with any n)
    FUNC_FEATURES = [
        "acid_energy.acid_energy","bileresistance.bileresistance",
        "gutpersistence.gutpersistence","osmoticstress.osmoticstress",
        "heatstress.heatstress","coldstress.coldstress","gaba.gaba",
        "vitamins.vitamins","immunomodulation.immunomodulation",
        "bacteriocins.bacteriocin","cellenvelope_eps.cellenvelope_eps",
        "adhesion_surface.adhesion_surface","adhesion_biofilm.adhesion_biofilm",
        "carbohydrate.carbohydrate","metabolism.metabolism",
        "defense_crispr.defense_crispr","antipath_qs.antipath_qs",
        "alkalinestress.alkalinestress","CAZyme_total",
    ]
    FUNC_FEATURES = [f for f in FUNC_FEATURES if f in df.columns]

    spear_rows = []
    for feat in FUNC_FEATURES:
        try:
            r, p = stats.spearmanr(df[feat].fillna(0), df["LAB_score_v1"])
            spear_rows.append({"feature": feat, "spearman_r": round(r,4), "p_value": round(p,6)})
        except Exception:
            pass

    if spear_rows:
        pd.DataFrame(spear_rows).sort_values("spearman_r", ascending=False)\
          .to_csv(f"{args.outdir}/module_spearman_correlation.tsv", sep="\t", index=False)
        print(f"[ML] Spearman correlations saved (n={len(df)}).")

    # Write placeholder feature importance tables
    for fname in [
        "RF_regression_all_feature_importance.tsv",
        "RF_regression_functional_feature_importance.tsv",
        "RF_classification_all_feature_importance.tsv",
        "RF_classification_functional_feature_importance.tsv",
    ]:
        pd.DataFrame(columns=["feature","importance"])\
          .to_csv(f"{args.outdir}/{fname}", sep="\t", index=False)

    print(f"[ML] Done. Add ≥{MIN_SAMPLES} genomes for full ML analysis.")
    exit(0)

# ── Define feature sets ───────────────────────────────────────────────────────
ALL_FEATURES = [
    "GI_survival_raw","Functional_raw","Fermentation_raw","CAZyme_total",
    "acid_energy.acid_energy","bileresistance.bileresistance",
    "gutpersistence.gutpersistence","osmoticstress.osmoticstress",
    "heatstress.heatstress","coldstress.coldstress","gaba.gaba",
    "vitamins.vitamins","immunomodulation.immunomodulation",
    "bacteriocins.bacteriocin","cellenvelope_eps.cellenvelope_eps",
    "adhesion_surface.adhesion_surface","adhesion_biofilm.adhesion_biofilm",
    "carbohydrate.carbohydrate","metabolism.metabolism",
    "defense_crispr.defense_crispr","antipath_qs.antipath_qs",
    "alkalinestress.alkalinestress",
    "amrfinder_hits","resfinder_hits","vfdb_hits",
    "biogenic_amines","hemolysin_markers",
]
FUNC_FEATURES = [
    "acid_energy.acid_energy","bileresistance.bileresistance",
    "gutpersistence.gutpersistence","osmoticstress.osmoticstress",
    "heatstress.heatstress","coldstress.coldstress","gaba.gaba",
    "vitamins.vitamins","immunomodulation.immunomodulation",
    "bacteriocins.bacteriocin","cellenvelope_eps.cellenvelope_eps",
    "adhesion_surface.adhesion_surface","adhesion_biofilm.adhesion_biofilm",
    "carbohydrate.carbohydrate","metabolism.metabolism",
    "defense_crispr.defense_crispr","antipath_qs.antipath_qs",
    "alkalinestress.alkalinestress","CAZyme_total",
]

# Filter to existing columns only
ALL_FEATURES  = [f for f in ALL_FEATURES  if f in df.columns]
FUNC_FEATURES = [f for f in FUNC_FEATURES if f in df.columns]

# Clean data — only use columns that exist
df_clean = df.dropna(subset=["LAB_score_v1"]).copy()
for col in ALL_FEATURES + FUNC_FEATURES:
    if col not in df_clean.columns:
        df_clean[col] = 0
df_clean[ALL_FEATURES]  = df_clean[ALL_FEATURES].fillna(0)
df_clean[FUNC_FEATURES] = df_clean[FUNC_FEATURES].fillna(0)

y_reg = df_clean["LAB_score_v1"]
X_all = df_clean[ALL_FEATURES]
X_fnc = df_clean[FUNC_FEATURES]

le    = LabelEncoder()
y_cls = le.fit_transform(df_clean["Priority_class"].fillna("Low"))

print(f"[ML] Feature sets: ALL={len(ALL_FEATURES)}, FUNC_ONLY={len(FUNC_FEATURES)}")
print(f"[ML] Samples for ML: {len(df_clean)}")

# ── Adaptive train/test split ─────────────────────────────────────────────────
TEST_SIZE = min(0.2, max(2/len(df_clean), 0.1))
CV_FOLDS  = min(5, len(df_clean) - 1)

def train_rf_regression(X, y, label=""):
    Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=TEST_SIZE, random_state=42)
    rf = RandomForestRegressor(n_estimators=200, max_features="sqrt", random_state=42, n_jobs=-1)
    rf.fit(Xtr, ytr)
    ypred = rf.predict(Xte)
    r2  = r2_score(yte, ypred) if len(yte) > 1 else float('nan')
    mae = mean_absolute_error(yte, ypred)
    cv  = cross_val_score(rf, X, y, cv=CV_FOLDS, scoring="r2", n_jobs=-1) if len(X) >= CV_FOLDS else np.array([r2])
    print(f"\n[ML] Regression ({label}): R²={r2:.4f}  MAE={mae:.4f}  CV_R²={cv.mean():.4f}±{cv.std():.4f}")
    fi = pd.Series(rf.feature_importances_, index=X.columns).sort_values(ascending=False)
    return rf, r2, mae, cv, fi, (Xte, yte, ypred)

def train_rf_classification(X, y, label=""):
    # Need at least 2 samples per class for stratified split
    unique, counts = np.unique(y, return_counts=True)
    can_stratify = all(c >= 2 for c in counts)
    Xtr, Xte, ytr, yte = train_test_split(
        X, y, test_size=TEST_SIZE, random_state=42,
        stratify=y if can_stratify else None
    )
    rf = RandomForestClassifier(n_estimators=200, max_features="sqrt",
                                class_weight="balanced", random_state=42, n_jobs=-1)
    rf.fit(Xtr, ytr)
    ypred = rf.predict(Xte)
    acc = (ypred == yte).mean()
    cv  = cross_val_score(rf, X, y, cv=CV_FOLDS, scoring="f1_weighted", n_jobs=-1) if len(X) >= CV_FOLDS else np.array([acc])
    print(f"\n[ML] Classification ({label}): Accuracy={acc:.4f}  CV_F1={cv.mean():.4f}±{cv.std():.4f}")
    fi = pd.Series(rf.feature_importances_, index=X.columns).sort_values(ascending=False)
    cm = confusion_matrix(yte, ypred)
    return rf, acc, cv, fi, cm, (Xte, yte, ypred)

# Run models
print("\n" + "="*60)
rf_reg_all, r2_all, mae_all, cv_reg_all, fi_reg_all, reg_pred_all = train_rf_regression(X_all, y_reg, "ALL")
rf_reg_fn,  r2_fn,  mae_fn,  cv_reg_fn,  fi_reg_fn,  reg_pred_fn  = train_rf_regression(X_fnc, y_reg, "FUNC")
rf_cls_all, acc_all, cv_cls_all, fi_cls_all, cm_all, cls_pred_all  = train_rf_classification(X_all, y_cls, "ALL")
rf_cls_fn,  acc_fn,  cv_cls_fn,  fi_cls_fn,  cm_fn,  cls_pred_fn   = train_rf_classification(X_fnc, y_cls, "FUNC")

# Save tables
for fi, fname in [
    (fi_reg_all, "RF_regression_all_feature_importance.tsv"),
    (fi_reg_fn,  "RF_regression_functional_feature_importance.tsv"),
    (fi_cls_all, "RF_classification_all_feature_importance.tsv"),
    (fi_cls_fn,  "RF_classification_functional_feature_importance.tsv"),
]:
    fi.reset_index().rename(columns={"index":"feature", 0:"importance"})\
      .to_csv(f"{args.outdir}/{fname}", sep="\t", index=False)

pd.DataFrame([
    {"Model":"RF_Regression_ALL",  "Metric":"R2",       "Value": r2_all},
    {"Model":"RF_Regression_ALL",  "Metric":"MAE",      "Value": mae_all},
    {"Model":"RF_Regression_FUNC", "Metric":"R2",       "Value": r2_fn},
    {"Model":"RF_Regression_FUNC", "Metric":"MAE",      "Value": mae_fn},
    {"Model":"RF_Class_ALL",       "Metric":"Accuracy", "Value": acc_all},
    {"Model":"RF_Class_FUNC",      "Metric":"Accuracy", "Value": acc_fn},
]).to_csv(f"{args.outdir}/model_performance_summary.tsv", sep="\t", index=False)

# Spearman correlations
spear_rows = []
for feat in FUNC_FEATURES:
    try:
        r, p = stats.spearmanr(df_clean[feat], df_clean["LAB_score_v1"])
        spear_rows.append({"feature": feat, "spearman_r": round(r,4), "p_value": round(p,6)})
    except Exception:
        pass
pd.DataFrame(spear_rows).sort_values("spearman_r", ascending=False)\
  .to_csv(f"{args.outdir}/module_spearman_correlation.tsv", sep="\t", index=False)

# Figures
COLORS = {"dark": "#1B4332", "mid": "#52B788", "light": "#B7E4C7"}
fig, axes = plt.subplots(1, 2, figsize=(18, 8))
for ax, fi, title in [
    (axes[0], fi_reg_fn.head(16), f"Regression (Functional-only)\nR² = {r2_fn:.3f}"),
    (axes[1], fi_cls_fn.head(16), f"Classification (Functional-only)\nAccuracy = {acc_fn:.3f}"),
]:
    fi_plot = fi.sort_values(ascending=True)
    colors  = ["#2D6A4F" if v >= fi_plot.median() else "#95D5B2" for v in fi_plot.values]
    ax.barh(fi_plot.index, fi_plot.values, color=colors, edgecolor="white", height=0.7)
    ax.axvline(fi_plot.median(), color="#F4A261", ls="--", lw=1.5, label="Median")
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_xlabel("Feature Importance")
    ax.legend(fontsize=9)
    ax.grid(axis="x", alpha=0.3)
plt.suptitle("Random Forest Feature Importance — LAB-Score v1.1", fontsize=14, fontweight="bold")
plt.tight_layout()
fig.savefig(f"{args.outdir}/fig_ml_feature_importance.png", dpi=300, bbox_inches="tight")
plt.close()

# Actual vs predicted
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
for ax, (Xte, yte, ypred), r2, label in [
    (axes[0], reg_pred_all, r2_all, f"All features (R²={r2_all:.3f})"),
    (axes[1], reg_pred_fn,  r2_fn,  f"Functional only (R²={r2_fn:.3f})"),
]:
    ax.scatter(yte, ypred, alpha=0.6, color="#52B788", edgecolors="#2D6A4F", s=60)
    if len(yte) > 1:
        lo, hi = min(yte.min(), ypred.min()), max(yte.max(), ypred.max())
        ax.plot([lo, hi], [lo, hi], "r--", lw=1.5)
    ax.set_xlabel("Actual LAB-Score"); ax.set_ylabel("Predicted LAB-Score")
    ax.set_title(label, fontsize=12, fontweight="bold")
    ax.grid(alpha=0.3)
plt.suptitle("Actual vs Predicted", fontsize=13, fontweight="bold")
plt.tight_layout()
fig.savefig(f"{args.outdir}/fig_ml_actual_vs_predicted.png", dpi=300, bbox_inches="tight")
plt.close()

# Confusion matrix
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
for ax, cm, title in [
    (axes[0], cm_all, f"All features (Acc={acc_all:.3f})"),
    (axes[1], cm_fn,  f"Functional only (Acc={acc_fn:.3f})"),
]:
    sns.heatmap(cm, annot=True, fmt="d", cmap="YlGn",
                xticklabels=le.classes_, yticklabels=le.classes_,
                ax=ax, linewidths=0.5, cbar=False)
    ax.set_xlabel("Predicted"); ax.set_ylabel("Actual")
    ax.set_title(title, fontsize=12, fontweight="bold")
plt.suptitle("Classification Confusion Matrix", fontsize=13, fontweight="bold")
plt.tight_layout()
fig.savefig(f"{args.outdir}/fig_ml_confusion_matrix.png", dpi=300, bbox_inches="tight")
plt.close()

print(f"\n[ML] All outputs saved to: {args.outdir}")
