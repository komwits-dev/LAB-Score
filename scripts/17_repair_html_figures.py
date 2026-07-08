#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUTDIR = Path("results_from_existing_prokka")
SCORE = OUTDIR / "09_scores/LAB_score_v1.tsv"
FIGDIR = OUTDIR / "11_figures"
FIGDIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(SCORE, sep="\t")

# -----------------------------
# Figure 2: LAB-score distribution
# -----------------------------
plt.figure(figsize=(10, 6))

# Background threshold zones
plt.axvspan(df["LAB_score_v1"].min(), 50, alpha=0.12)
plt.axvspan(60, 80, alpha=0.12)
plt.axvspan(80, 95, alpha=0.12)
plt.axvspan(85, df["LAB_score_v1"].max() + 2, alpha=0.12)

plt.hist(df["LAB_score_v1"], bins=40, edgecolor="black", linewidth=0.4)

for x in [60, 80, 95]:
    plt.axvline(x, linestyle="--", linewidth=1.5)

plt.text(40, plt.ylim()[1]*0.92, "Low\n(<60)", ha="center", va="top", fontsize=11, fontweight="bold")
plt.text(60, plt.ylim()[1]*0.92, "Moderate\n(60–79.9)", ha="center", va="top", fontsize=11, fontweight="bold")
plt.text(77.5, plt.ylim()[1]*0.92, "High\n(80–94.9)", ha="center", va="top", fontsize=11, fontweight="bold")
plt.text(88, plt.ylim()[1]*0.92, "Elite\n(≥95)", ha="center", va="top", fontsize=11, fontweight="bold")

plt.xlabel("LAB-Score v1")
plt.ylabel("Number of genomes")
plt.title("LAB-Score distribution with prioritization thresholds", fontsize=15, fontweight="bold")
plt.text(0.98, 0.03, f"n = {len(df):,} genomes", transform=plt.gca().transAxes, ha="right", fontsize=11)
plt.grid(axis="y", linestyle="--", alpha=0.35)
plt.tight_layout()
plt.savefig(FIGDIR / "fig2_score_distribution.png", dpi=300)
plt.close()

# -----------------------------
# Figure 3: candidate tier distribution
# -----------------------------
tier_col = "Candidate_tier_v1_1" if "Candidate_tier_v1_1" in df.columns else "Priority_class"
safety_col = "Refined_safety_status" if "Refined_safety_status" in df.columns else "flag"

tier_counts = df[tier_col].value_counts().reset_index()
tier_counts.columns = ["Candidate tier", "n"]

preferred = [
    "Cautionary Elite candidate",
    "Elite candidate",
    "Cautionary High candidate",
    "High candidate",
    "Cautionary Moderate candidate",
    "Moderate candidate",
    "Cautionary Low priority",
    "Low priority",
    "Critical safety review",
]

tier_counts["order"] = tier_counts["Candidate tier"].apply(
    lambda x: preferred.index(x) if x in preferred else 99
)
tier_counts = tier_counts.sort_values(["order", "n"], ascending=[True, False]).drop(columns="order")

safety_counts = df[safety_col].value_counts().reset_index()
safety_counts.columns = ["Safety status", "n"]

fig, axes = plt.subplots(1, 2, figsize=(13, 6), gridspec_kw={"width_ratios": [1.7, 1]})

# Panel A
y = np.arange(len(tier_counts))
axes[0].barh(y, tier_counts["n"])
axes[0].set_yticks(y)
axes[0].set_yticklabels(tier_counts["Candidate tier"])
axes[0].invert_yaxis()
axes[0].set_xlabel("Number of genomes")
axes[0].set_title("Candidate tier distribution", fontweight="bold")

for i, v in enumerate(tier_counts["n"]):
    pct = v / len(df) * 100
    axes[0].text(v + max(tier_counts["n"])*0.01, i, f"{v:,} ({pct:.1f}%)", va="center", fontsize=9)

axes[0].grid(axis="x", linestyle="--", alpha=0.35)

# Panel B
axes[1].pie(
    safety_counts["n"],
    labels=[f"{a}\n{b:,} ({b/len(df)*100:.1f}%)" for a, b in zip(safety_counts["Safety status"], safety_counts["n"])],
    startangle=90,
    textprops={"fontsize": 9}
)
axes[1].set_title("Refined safety status", fontweight="bold")

fig.suptitle("Refined safety-gated candidate tier distribution", fontsize=15, fontweight="bold")
plt.tight_layout()
plt.savefig(FIGDIR / "fig3_candidate_tiers.png", dpi=300)
plt.close()

# -----------------------------
# Figure 5: functional feature importance if available
# -----------------------------
possible_importance = [
    OUTDIR / "12_ml_functional_only/RF_functional_only_candidate_group_feature_importance.tsv",
    OUTDIR / "10_ml/RF_functional_only_candidate_group_feature_importance.tsv",
    OUTDIR / "10_ml/RF_candidate_group_feature_importance.tsv",
]

imp_file = None
for p in possible_importance:
    if p.exists():
        imp_file = p
        break

if imp_file:
    imp = pd.read_csv(imp_file, sep="\t")
    imp = imp.head(20).copy()
    imp = imp.sort_values("importance")

    plt.figure(figsize=(8, 7))
    plt.barh(imp["feature"], imp["importance"])
    plt.xlabel("Random Forest importance")
    plt.ylabel("Functional module")
    plt.title("Functional modules driving LAB-Score candidate groups", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(FIGDIR / "fig5_functional_drivers.png", dpi=300)
    plt.close()
    print("Saved Figure 5 from:", imp_file)
else:
    print("No feature-importance file found for Figure 5.")

print("Saved figures to:", FIGDIR)
