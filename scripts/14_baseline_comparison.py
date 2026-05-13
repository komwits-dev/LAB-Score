#!/usr/bin/env python3
"""
14_baseline_comparison.py
=========================
Compares LAB-Score v1.1 with simpler single-component or reduced-component
baseline rankings.
"""

import argparse
import os
import pandas as pd
from scipy.stats import spearmanr

ap = argparse.ArgumentParser(description="LAB-Score baseline comparison")
ap.add_argument("--scores", required=True, help="LAB_score_v1.tsv from 09_calculate_lab_score.py")
ap.add_argument("--outdir", required=True, help="Output directory")
ap.add_argument("--top-n", type=int, default=100, help="Top-N overlap to evaluate [100]")
args = ap.parse_args()
os.makedirs(args.outdir, exist_ok=True)

df = pd.read_csv(args.scores, sep="\t")

if "LAB_score_v1" not in df.columns:
    raise SystemExit("[baseline] Missing LAB_score_v1 column")

if "Accession" not in df.columns:
    for alt in ["accession", "genome", "Strain", "strain"]:
        if alt in df.columns:
            df["Accession"] = df[alt].astype(str)
            break
if "Accession" not in df.columns:
    df["Accession"] = [f"genome_{i+1}" for i in range(len(df))]

for col in ["Safety_score", "GI_survival_score", "Functional_score", "Fermentation_score"]:
    if col not in df.columns:
        raise SystemExit(f"[baseline] Missing required column: {col}")

# Simple baseline scores
df["Safety_only_score"] = df["Safety_score"].astype(float)
df["GI_only_score"] = df["GI_survival_score"].astype(float)
df["Functional_only_score"] = df["Functional_score"].astype(float)
df["Fermentation_only_score"] = df["Fermentation_score"].astype(float)
df["Function_plus_fermentation_score"] = 0.50 * df["Functional_score"].astype(float) + 0.50 * df["Fermentation_score"].astype(float)
df["No_safety_score"] = 0.40 * df["GI_survival_score"].astype(float) + 0.40 * df["Functional_score"].astype(float) + 0.20 * df["Fermentation_score"].astype(float)

baseline_cols = [
    "Safety_only_score",
    "GI_only_score",
    "Functional_only_score",
    "Fermentation_only_score",
    "Function_plus_fermentation_score",
    "No_safety_score",
]

main_col = "LAB_score_v1"
top_main = set(df.sort_values(main_col, ascending=False).head(args.top_n)["Accession"])
if "Priority_class" in df.columns:
    elite_main = set(df[df["Priority_class"] == "Elite"]["Accession"])
else:
    elite_main = set(df[df[main_col] >= 85]["Accession"])

rows = []
for col in baseline_cols:
    rho, p = spearmanr(df[main_col], df[col])
    top_base = set(df.sort_values(col, ascending=False).head(args.top_n)["Accession"])
    rows.append({
        "baseline": col,
        "spearman_vs_LAB_score": rho,
        "p_value": p,
        f"top{args.top_n}_overlap_with_LAB_score": len(top_main & top_base),
        f"top{args.top_n}_overlap_percent": round(len(top_main & top_base) / args.top_n * 100, 2),
        f"elite_recovered_in_baseline_top{args.top_n}": len(elite_main & top_base),
        "n_main_elite": len(elite_main),
    })

summary = pd.DataFrame(rows)
summary.to_csv(os.path.join(args.outdir, "baseline_comparison_summary.tsv"), sep="\t", index=False)

id_cols = [c for c in ["Accession", "genome", "full_species", "genus", "species", "strain", "LAB_score_v1", "Priority_class", "Candidate_tier_v1_1"] if c in df.columns]
df[id_cols + baseline_cols].to_csv(os.path.join(args.outdir, "baseline_scores.tsv"), sep="\t", index=False)

print("[baseline] Saved:", os.path.join(args.outdir, "baseline_comparison_summary.tsv"))
print(summary.to_string(index=False))
