#!/usr/bin/env python3
"""
09_calculate_lab_score.py
=========================
Calculates LAB-Score v1.2 from the master matrix.

Formula:
    LAB_score_v1 = 0.45 × Safety_score
                 + 0.25 × GI_survival_score
                 + 0.20 × Functional_score
                 + 0.10 × Fermentation_score

Priority classes:     ≥95=Elite  80-94=High  60-79=Moderate  <60=Low
Refined safety gate:  Critical / Cautionary / Pass
"""

import argparse, os
import pandas as pd
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--matrix", required=True)
ap.add_argument("--outdir", required=True)
args = ap.parse_args()
os.makedirs(args.outdir, exist_ok=True)

df = pd.read_csv(args.matrix, sep="\t")
print(f"[score] Input: {df.shape[0]} genomes × {df.shape[1]} features")

# ── Percentile rank helper ─────────────────────────────────────────────────────
def pct_rank(series):
    return series.rank(pct=True, method="average") * 100

def safe_col(df, col, default=0):
    return df[col] if col in df.columns else pd.Series(default, index=df.index)

# ═══════════════════════════════════════════════════════════════════════════════
# 1. Safety Score
# ═══════════════════════════════════════════════════════════════════════════════
# Penalty-based: start at 100, deduct for each safety concern
critical_amr_present = (
    safe_col(
        df,
        "amrfinder_critical_hits",
        default=safe_col(df, "amrfinder_hits").astype(float),
    ).astype(float) > 0
)
caution_amr_present = (
    safe_col(df, "amrfinder_caution_hits").astype(float) > 0
) & (~critical_amr_present)

df["Safety_penalty"] = (
    critical_amr_present.astype(int) * 40 +
    caution_amr_present.astype(int) * 15 +
    (safe_col(df,"vfdb_hits").astype(float)         > 0).astype(int) * 30 +
    (safe_col(df,"resfinder_hits").astype(float)    > 0).astype(int) * 20 +
    (safe_col(df,"hemolysin_markers").astype(float) > 0).astype(int) * 15 +
    (safe_col(df,"biogenic_amines").astype(float)   > 0).astype(int) * 15
)
df["Safety_score"] = (100 - df["Safety_penalty"]).clip(lower=0)

# ═══════════════════════════════════════════════════════════════════════════════
# 2. GI Survival Score
# ═══════════════════════════════════════════════════════════════════════════════
GI_COLS = [
    "acid_energy.acid_energy",
    "bileresistance.bileresistance",
    "gutpersistence.gutpersistence",
    "osmoticstress.osmoticstress",
    "heatstress.heatstress",
    "coldstress.coldstress",
]
GI_COLS = [c for c in GI_COLS if c in df.columns]
df["GI_survival_raw"]   = df[GI_COLS].sum(axis=1)
df["GI_survival_score"] = pct_rank(df["GI_survival_raw"])

# ═══════════════════════════════════════════════════════════════════════════════
# 3. Functional Score (probiotic-associated)
# ═══════════════════════════════════════════════════════════════════════════════
FUNC_COLS = [
    "gaba.gaba",
    "vitamins.vitamins",
    "immunomodulation.immunomodulation",
    "bacteriocins.bacteriocins",
    "cellenvelope_eps.cellenvelope_eps",
    "adhesion_surface.adhesion_surface",
    "adhesion_biofilm.adhesion_biofilm",
    "defense_crispr.defense_crispr",
    "antipath_qs.antipath_qs",
    "alkalinestress.alkalinestress",
]
FUNC_COLS = [c for c in FUNC_COLS if c in df.columns]
df["Functional_raw"]   = df[FUNC_COLS].sum(axis=1)
df["Functional_score"] = pct_rank(df["Functional_raw"])

# ═══════════════════════════════════════════════════════════════════════════════
# 4. Fermentation Score
# ═══════════════════════════════════════════════════════════════════════════════
FERM_COLS = [
    "carbohydrate.carbohydrate",
    "metabolism.metabolism",
    "CAZyme_total",
]
FERM_COLS = [c for c in FERM_COLS if c in df.columns]
df["Fermentation_raw"]   = df[FERM_COLS].sum(axis=1)
df["Fermentation_score"] = pct_rank(df["Fermentation_raw"])

# ═══════════════════════════════════════════════════════════════════════════════
# 5. Composite LAB-Score v1.2
# ═══════════════════════════════════════════════════════════════════════════════
df["LAB_score_v1"] = (
    0.45 * df["Safety_score"] +
    0.25 * df["GI_survival_score"] +
    0.20 * df["Functional_score"] +
    0.10 * df["Fermentation_score"]
).round(6)

# ═══════════════════════════════════════════════════════════════════════════════
# 6. Priority classes (raw, before safety gate)
# ═══════════════════════════════════════════════════════════════════════════════
def assign_class(score):
    if score >= 95: return "Elite"
    elif score >= 80: return "High"
    elif score >= 60: return "Moderate"
    else: return "Low"

df["Priority_class"] = df["LAB_score_v1"].apply(assign_class)

# ═══════════════════════════════════════════════════════════════════════════════
# 7. Refined safety-gated candidate tiers (v1.1)
# ═══════════════════════════════════════════════════════════════════════════════
def refined_tier(row):
    flag = str(row.get("flag", "")).strip()
    pc   = row["Priority_class"]

    # Keep low-scoring genomes as "Low priority" rather than "Low candidate".
    # This makes the final v1.1 tier names consistent with the manuscript text.
    if flag == "Critical":
        return "Critical", "Critical safety review"
    elif flag == "Caution":
        if pc == "Low":
            return "Caution", "Cautionary Low priority"
        return "Caution", f"Cautionary {pc} candidate"
    else:
        if pc == "Low":
            return "Pass", "Low priority"
        return "Pass", f"{pc} candidate"

tiers = df.apply(refined_tier, axis=1)
df["Refined_safety_status"] = [t[0] for t in tiers]
df["Candidate_tier_v1_1"]   = [t[1] for t in tiers]

# ═══════════════════════════════════════════════════════════════════════════════
# 8. Select and order output columns
# ═══════════════════════════════════════════════════════════════════════════════
ID_COLS = ["genome","accession","genus","species","full_species","strain"]
ID_COLS = [c for c in ID_COLS if c in df.columns]

SCORE_COLS = [
    "Safety_score","GI_survival_score","Functional_score","Fermentation_score",
    "LAB_score_v1","Priority_class",
    "Refined_safety_status","Candidate_tier_v1_1",
    "Safety_penalty","flag","reasons",
    "amrfinder_hits","amrfinder_critical_hits","amrfinder_caution_hits",
    "amrfinder_critical_symbols","amrfinder_caution_symbols",
    "vfdb_hits","amrfinder_stress_hits","resfinder_hits",
    "biogenic_amines","biogenic_amine_details",
    "nonspecific_decarboxylase_annotations","nonspecific_decarboxylase_details",
    "gaba_related_annotations","gaba_related_details",
    "hemolysin_markers","hemolysin_marker_details",
    "hemolysin_like_annotations","hemolysin_like_details",
    "GI_survival_raw","Functional_raw","Fermentation_raw",
    "CAZyme_total",
]
SCORE_COLS = [c for c in SCORE_COLS if c in df.columns]

GENOME_QC = ["Completeness","Contamination","contigs","total_length","gc_pct",
             "assembly_level","Contigs","TotalLength","GC%"]
GENOME_QC = [c for c in GENOME_QC if c in df.columns]

DETAIL_COLS = GI_COLS + FUNC_COLS + FERM_COLS
DETAIL_COLS = list(dict.fromkeys(DETAIL_COLS))

out_cols = ID_COLS + SCORE_COLS + GENOME_QC + DETAIL_COLS
out_cols = [c for c in out_cols if c in df.columns]

out = df[out_cols].sort_values("LAB_score_v1", ascending=False).reset_index(drop=True)

# ═══════════════════════════════════════════════════════════════════════════════
# 9. Summary statistics
# ═══════════════════════════════════════════════════════════════════════════════
print(f"\n[score] LAB-Score range: {out['LAB_score_v1'].min():.2f} – {out['LAB_score_v1'].max():.2f}")
print(f"[score] Mean: {out['LAB_score_v1'].mean():.2f}  Median: {out['LAB_score_v1'].median():.2f}")

print("\n[score] Priority class counts:")
print(out["Priority_class"].value_counts().to_string())

print("\n[score] Candidate tier counts (v1.1):")
print(out["Candidate_tier_v1_1"].value_counts().to_string())

print("\n[score] Top 10 strains:")
show_cols = [c for c in ["genome","full_species","LAB_score_v1","Priority_class","Candidate_tier_v1_1"] if c in out.columns]
print(out[show_cols].head(10).to_string(index=False))

# ═══════════════════════════════════════════════════════════════════════════════
# 10. Save outputs
# ═══════════════════════════════════════════════════════════════════════════════
# Main scores
out.to_csv(f"{args.outdir}/LAB_score_v1.tsv", sep="\t", index=False)
print(f"\n[score] Saved: {args.outdir}/LAB_score_v1.tsv")

# Species mean scores
sp_group_cols = [c for c in ["genus","full_species","species"] if c in out.columns]
if sp_group_cols:
    sp_col = sp_group_cols[0] if "full_species" not in sp_group_cols else "full_species"
    sp_df = out.groupby(sp_col).agg(
        N                 = ("LAB_score_v1","count"),
        mean_LAB_score    = ("LAB_score_v1","mean"),
        median_LAB_score  = ("LAB_score_v1","median"),
        max_LAB_score     = ("LAB_score_v1","max"),
        min_LAB_score     = ("LAB_score_v1","min"),
        sd_LAB_score      = ("LAB_score_v1","std"),
        Elite             = ("Priority_class", lambda x: (x=="Elite").sum()),
        High              = ("Priority_class", lambda x: (x=="High").sum()),
        Moderate          = ("Priority_class", lambda x: (x=="Moderate").sum()),
        Low               = ("Priority_class", lambda x: (x=="Low").sum()),
        Critical          = ("Refined_safety_status", lambda x: (x=="Critical").sum()),
    ).reset_index().sort_values("mean_LAB_score", ascending=False)
    sp_df.to_csv(f"{args.outdir}/species_mean_scores.tsv", sep="\t", index=False)
    print(f"[score] Saved: {args.outdir}/species_mean_scores.tsv")

# Priority class counts
out["Priority_class"].value_counts().reset_index().rename(
    columns={"index":"Priority_class","Priority_class":"Count"}
).to_csv(f"{args.outdir}/priority_class_counts.tsv", sep="\t", index=False)

# Top 100 non-critical
top100 = out[out["Refined_safety_status"] != "Critical"].head(100)
top100.to_csv(f"{args.outdir}/top100_candidates.tsv", sep="\t", index=False)
print(f"[score] Saved: {args.outdir}/top100_candidates.tsv  (n={len(top100)})")

# Score distribution summary
dist_stats = out["LAB_score_v1"].describe(percentiles=[0.10,0.25,0.50,0.75,0.90])
dist_stats.to_csv(f"{args.outdir}/LAB_score_distribution_summary.tsv", sep="\t", header=["value"])
print(f"[score] Saved: {args.outdir}/LAB_score_distribution_summary.tsv")

print("\n[score] All outputs saved.")
