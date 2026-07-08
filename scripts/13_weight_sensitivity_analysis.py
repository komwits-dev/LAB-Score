#!/usr/bin/env python3
"""
13_weight_sensitivity_analysis.py
=================================
Robustness analysis for LAB-Score v1.2 component weights.

Compares the main safety-weighted LAB-Score against alternative weighting
schemes and reports rank correlation, top-candidate overlap, and class counts.
"""

import argparse
import os
import pandas as pd
from scipy.stats import spearmanr

ap = argparse.ArgumentParser(description="LAB-Score v1.2 weight sensitivity analysis")
ap.add_argument("--scores", required=True, help="LAB_score_v1.tsv from 09_calculate_lab_score.py")
ap.add_argument("--outdir", required=True, help="Output directory")
ap.add_argument("--top-n", type=int, default=100, help="Top-N overlap to evaluate [100]")
args = ap.parse_args()
os.makedirs(args.outdir, exist_ok=True)

df = pd.read_csv(args.scores, sep="\t")

required = ["Safety_score", "GI_survival_score", "Functional_score", "Fermentation_score"]
missing = [c for c in required if c not in df.columns]
if missing:
    raise SystemExit(f"[sensitivity] Missing required columns: {missing}")

if "Accession" not in df.columns:
    for alt in ["accession", "genome", "Strain", "strain"]:
        if alt in df.columns:
            df["Accession"] = df[alt].astype(str)
            break
if "Accession" not in df.columns:
    df["Accession"] = [f"genome_{i+1}" for i in range(len(df))]

weight_sets = {
    "main_safety_heavy_0.45_0.25_0.20_0.10": {
        "Safety_score": 0.45,
        "GI_survival_score": 0.25,
        "Functional_score": 0.20,
        "Fermentation_score": 0.10,
    },
    "balanced_equal_0.25_each": {
        "Safety_score": 0.25,
        "GI_survival_score": 0.25,
        "Functional_score": 0.25,
        "Fermentation_score": 0.25,
    },
    "function_heavy_0.30_0.20_0.35_0.15": {
        "Safety_score": 0.30,
        "GI_survival_score": 0.20,
        "Functional_score": 0.35,
        "Fermentation_score": 0.15,
    },
    "GI_survival_heavy_0.30_0.35_0.25_0.10": {
        "Safety_score": 0.30,
        "GI_survival_score": 0.35,
        "Functional_score": 0.25,
        "Fermentation_score": 0.10,
    },
    "fermentation_heavy_0.30_0.20_0.20_0.30": {
        "Safety_score": 0.30,
        "GI_survival_score": 0.20,
        "Functional_score": 0.20,
        "Fermentation_score": 0.30,
    },
}

for name, weights in weight_sets.items():
    df[f"LAB_score_{name}"] = sum(df[col].astype(float) * w for col, w in weights.items())

def class_from_score(score):
    if score >= 95:
        return "Elite"
    if score >= 80:
        return "High"
    if score >= 60:
        return "Moderate"
    return "Low"

main_col = "LAB_score_main_safety_heavy_0.45_0.25_0.20_0.10"
if "LAB_score_v1" in df.columns:
    # Use the pipeline-calculated score as the official reference.
    main_reference = df["LAB_score_v1"].astype(float)
else:
    main_reference = df[main_col].astype(float)

top_main = set(df.assign(_main=main_reference).sort_values("_main", ascending=False).head(args.top_n)["Accession"])

rows = []
for name in weight_sets:
    col = f"LAB_score_{name}"
    class_col = f"Priority_{name}"
    df[class_col] = df[col].apply(class_from_score)
    rho, p = spearmanr(main_reference, df[col])
    top_alt = set(df.sort_values(col, ascending=False).head(args.top_n)["Accession"])
    overlap = len(top_main & top_alt)
    counts = df[class_col].value_counts().to_dict()
    rows.append({
        "weighting_scheme": name,
        "spearman_vs_main_LAB_score": rho,
        "p_value": p,
        f"top{args.top_n}_overlap_with_main": overlap,
        f"top{args.top_n}_overlap_percent": round(overlap / args.top_n * 100, 2),
        "n_elite": counts.get("Elite", 0),
        "n_high": counts.get("High", 0),
        "n_moderate": counts.get("Moderate", 0),
        "n_low": counts.get("Low", 0),
    })

summary = pd.DataFrame(rows)
summary.to_csv(os.path.join(args.outdir, "weight_sensitivity_summary.tsv"), sep="\t", index=False)

id_cols = [c for c in ["Accession", "genome", "full_species", "genus", "species", "strain"] if c in df.columns]
score_cols = [f"LAB_score_{x}" for x in weight_sets]
priority_cols = [f"Priority_{x}" for x in weight_sets]
df[id_cols + (["LAB_score_v1"] if "LAB_score_v1" in df.columns else []) + score_cols + priority_cols].to_csv(
    os.path.join(args.outdir, "weight_sensitivity_scores.tsv"), sep="\t", index=False
)

print("[sensitivity] Saved:", os.path.join(args.outdir, "weight_sensitivity_summary.tsv"))
print(summary.to_string(index=False))
