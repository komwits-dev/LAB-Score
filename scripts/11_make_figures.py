#!/usr/bin/env python3
"""
11_make_figures.py — Publication figures for LAB-Score v1.1/v2.0

Fig 2: Score distribution + priority class pie
Fig 3: Refined safety-gated candidate tiers + safety status
Fig 4: Species-level LAB-Score and candidate-tier summary
Fig 5: ML feature importance

Fixes in this version:
  * Fig. 4 no longer becomes blank when species_mean_scores.tsv is empty or has only one row.
  * Species summaries are rebuilt directly from LAB_score_v1.tsv, so tier counts are always available.
  * If species metadata were not resolved, Fig. 4 falls back to genome-level top-candidate and tier summaries.
  * The final tier names include "Low priority" and "Cautionary Low priority".
  * Each figure is exported as PNG, PDF, and TIFF.
"""

import argparse, os, warnings, textwrap
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy.stats import gaussian_kde

warnings.filterwarnings("ignore")

ap = argparse.ArgumentParser()
ap.add_argument("--scores",  required=True)
ap.add_argument("--species", required=True)
ap.add_argument("--ml",      required=True)
ap.add_argument("--outdir",  required=True)
args = ap.parse_args()
os.makedirs(args.outdir, exist_ok=True)

# ── Theme ─────────────────────────────────────────────────────────────────────
PALETTE = {
    "Elite":    "#1B4332",
    "High":     "#2D6A4F",
    "Moderate": "#74C69D",
    "Low":      "#D8F3DC",
    "Critical": "#C1121F",
    "Caution":  "#E07C24",
    "Pass":     "#2D6A4F",
}

TIER_COLORS = {
    "Elite candidate":               "#1B4332",
    "High candidate":                "#2D6A4F",
    "Moderate candidate":            "#52B788",
    "Low priority":                  "#B7E4C7",
    "Low candidate":                 "#B7E4C7",  # backward-compatible with old score files
    "Cautionary Elite candidate":    "#95D5B2",
    "Cautionary High candidate":     "#E07C24",
    "Cautionary Moderate candidate": "#F4A261",
    "Cautionary Low priority":       "#FFDDD2",
    "Cautionary Low candidate":      "#FFDDD2",  # backward-compatible with old score files
    "Critical safety review":        "#C1121F",
}
TIER_ORDER = [
    "Critical safety review",
    "Cautionary Low priority",
    "Cautionary Low candidate",
    "Low priority",
    "Low candidate",
    "Cautionary Moderate candidate",
    "Moderate candidate",
    "Cautionary High candidate",
    "High candidate",
    "Cautionary Elite candidate",
    "Elite candidate",
]

plt.rcParams.update({
    "font.family":       "DejaVu Sans",
    "font.size":         11,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "figure.dpi":        300,
    "savefig.dpi":       600,
})


def save_figure(fig, basename):
    """Save publication-ready outputs in three formats. TIFF uses LZW compression."""
    fig.savefig(os.path.join(args.outdir, f"{basename}.png"),  dpi=600, bbox_inches="tight")
    fig.savefig(os.path.join(args.outdir, f"{basename}.pdf"),  dpi=600, bbox_inches="tight")
    fig.savefig(os.path.join(args.outdir, f"{basename}.tiff"), dpi=600, bbox_inches="tight",
                pil_kwargs={"compression": "tiff_lzw"})


def wrap_labels(labels, width=26):
    return ["\n".join(textwrap.wrap(str(x), width=width, break_long_words=False)) for x in labels]


def clean_species_column(df):
    """Return a robust species label series from the score table."""
    if "full_species" in df.columns:
        sp = df["full_species"].astype(str).str.strip()
    elif {"genus", "species"}.issubset(df.columns):
        sp = (df["genus"].astype(str).str.strip() + " " + df["species"].astype(str).str.strip())
    elif "species" in df.columns:
        sp = df["species"].astype(str).str.strip()
    else:
        sp = pd.Series("Unknown sp.", index=df.index)
    sp = sp.replace({"": "Unknown sp.", "nan": "Unknown sp.", "None": "Unknown sp."})
    return sp


def build_species_summary(score_df):
    """Build species summary directly from LAB_score_v1.tsv."""
    tmp = score_df.copy()
    tmp["_species_label"] = clean_species_column(tmp)

    if "Candidate_tier_v1_1" not in tmp.columns:
        tmp["Candidate_tier_v1_1"] = tmp.get("Priority_class", "Unknown")

    base = tmp.groupby("_species_label").agg(
        N                 = ("LAB_score_v1", "count"),
        mean_LAB_score    = ("LAB_score_v1", "mean"),
        median_LAB_score  = ("LAB_score_v1", "median"),
        max_LAB_score     = ("LAB_score_v1", "max"),
        min_LAB_score     = ("LAB_score_v1", "min"),
        sd_LAB_score      = ("LAB_score_v1", "std"),
    ).reset_index().rename(columns={"_species_label": "full_species"})

    tier_counts = pd.crosstab(tmp["_species_label"], tmp["Candidate_tier_v1_1"]).reset_index().rename(columns={"_species_label": "full_species"})
    out = base.merge(tier_counts, on="full_species", how="left")
    for t in TIER_ORDER:
        if t not in out.columns:
            out[t] = 0
    out["sd_LAB_score"] = out["sd_LAB_score"].fillna(0)
    return out


def species_metadata_is_resolved(score_df):
    sp = clean_species_column(score_df)
    unique = set(sp.dropna().unique())
    bad = {"Unknown sp.", "Unknown", "sp.", "Unknown nan", "nan sp."}
    if len(unique) == 0:
        return False
    if len(unique - bad) == 0:
        return False
    return True


def load_fi(fname):
    try:
        return pd.read_csv(os.path.join(args.ml, fname), sep="\t")
    except Exception:
        return pd.DataFrame(columns=["feature","importance"])


df = pd.read_csv(args.scores, sep="\t")
print(f"[figures] Loaded {len(df)} genomes")

try:
    sp_df_from_file = pd.read_csv(args.species, sep="\t")
except Exception:
    sp_df_from_file = pd.DataFrame()

fi_reg  = load_fi("RF_regression_functional_feature_importance.tsv")
fi_cls  = load_fi("RF_classification_functional_feature_importance.tsv")

# ═══════════════════════════════════════════════════════════════════════════════
# Fig 2 — LAB-Score Distribution
# ═══════════════════════════════════════════════════════════════════════════════
fig = plt.figure(figsize=(16, 6))
gs  = gridspec.GridSpec(1, 3, width_ratios=[3, 0.05, 1.2], figure=fig)
ax_dist = fig.add_subplot(gs[0])
ax_pie  = fig.add_subplot(gs[2])

scores = df["LAB_score_v1"].dropna().values
x_grid = np.linspace(max(0, scores.min()-2), min(100, scores.max()+2), 500)

regions = [
    (0,  50, "#FFDDD2", "Low\n(<50)"),
    (50, 70, "#D8F3DC", "Moderate\n(50–69)"),
    (70, 85, "#95D5B2", "High\n(70–84)"),
    (85,100, "#2D6A4F", "Elite\n(≥85)"),
]
for lo, hi, col, label in regions:
    ax_dist.axvspan(lo, hi, color=col, alpha=0.35, zorder=1)
    ax_dist.text((lo + hi) / 2, 0.001, label, ha="center", va="bottom",
                 fontsize=8.5, color="#333", fontweight="bold", zorder=5)

ax_dist.hist(scores, bins=min(30, len(scores)), density=True,
             color="#52B788", edgecolor="white", alpha=0.6, zorder=2)

if len(scores) > 3:
    try:
        kde = gaussian_kde(scores, bw_method=0.3)
        ax_dist.plot(x_grid, kde(x_grid), color="#1B4332", lw=2.5, zorder=3, label="KDE")
        ax_dist.legend(fontsize=9)
    except Exception:
        pass

for val, col in [(50,"#E07C24"),(70,"#2D6A4F"),(85,"#1B4332")]:
    ax_dist.axvline(val, color=col, lw=1.8, ls="--", zorder=4)

ax_dist.set_xlabel("LAB-Score v1.1")
ax_dist.set_ylabel("Density")
ax_dist.set_title(f"LAB-Score Distribution (n = {len(scores):,} genomes)", fontweight="bold")

class_counts = df["Priority_class"].value_counts()
pie_order  = [c for c in ["Elite","High","Moderate","Low"] if c in class_counts.index]
pie_vals   = [class_counts[c] for c in pie_order]
pie_colors = [PALETTE[c] for c in pie_order]

wedges, _, autotexts = ax_pie.pie(
    pie_vals, labels=None, autopct="%1.1f%%",
    colors=pie_colors, startangle=140,
    wedgeprops=dict(edgecolor="white", linewidth=2),
    pctdistance=0.75,
)
for at in autotexts:
    at.set_fontsize(9)

legend_patches = [mpatches.Patch(color=PALETTE[c], label=f"{c} (n={class_counts[c]:,})") for c in pie_order]
ax_pie.legend(handles=legend_patches, loc="lower center", bbox_to_anchor=(0.5, -0.25), fontsize=9, frameon=False)
ax_pie.set_title("Priority Classes", fontweight="bold", fontsize=12)

plt.suptitle("Figure 2. LAB-Score v1.1 Distribution and Priority Classification",
             fontsize=13, fontweight="bold", y=1.02)
plt.tight_layout()
save_figure(fig, "fig2_score_distribution")
plt.close()
print("[figures] Fig 2 saved")

# ═══════════════════════════════════════════════════════════════════════════════
# Fig 3 — Candidate Tiers
# ═══════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(17, 6), gridspec_kw={"width_ratios":[2.25,1]})

tier_counts = df["Candidate_tier_v1_1"].value_counts()
plot_tiers  = [t for t in TIER_ORDER if t in tier_counts.index]
# Add any unexpected categories at the end instead of dropping them.
plot_tiers += [t for t in tier_counts.index if t not in plot_tiers]
vals        = [tier_counts[t] for t in plot_tiers]
colors_bar  = [TIER_COLORS.get(t, "#999999") for t in plot_tiers]

ax = axes[0]
bars = ax.barh(plot_tiers, vals, color=colors_bar, edgecolor="white", height=0.65, zorder=2)
for bar, val in zip(bars, vals):
    ax.text(val + max(vals)*0.01, bar.get_y() + bar.get_height()/2,
            f"n = {val:,}", va="center", fontsize=10)
ax.set_xlabel("Number of Genomes")
ax.set_title("Candidate Tier Distribution", fontweight="bold")
ax.set_xlim(0, max(vals) * 1.25)
ax.grid(axis="x", alpha=0.4, zorder=1)
ax.set_axisbelow(True)

ax2 = axes[1]
safety_counts = df["Refined_safety_status"].value_counts()
s_order  = [s for s in ["Pass","Caution","Critical"] if s in safety_counts.index]
s_colors = [PALETTE.get(s,"#999999") for s in s_order]
s_vals   = [safety_counts[s] for s in s_order]

wedges2, _, autotexts2 = ax2.pie(
    s_vals, labels=None, autopct="%1.1f%%",
    colors=s_colors, startangle=90,
    wedgeprops=dict(edgecolor="white", linewidth=2),
    pctdistance=0.75,
)
for at in autotexts2:
    at.set_fontsize(10)

safety_patches = [mpatches.Patch(color=c, label=f"{s} (n={v:,})") for s, c, v in zip(s_order, s_colors, s_vals)]
ax2.legend(handles=safety_patches, loc="lower center", bbox_to_anchor=(0.5, -0.20), fontsize=9, frameon=False)
ax2.set_title("Refined Safety Status", fontweight="bold")

plt.suptitle("Figure 3. Refined Safety-Gated Candidate Tier Classification",
             fontsize=13, fontweight="bold", y=1.02)
plt.tight_layout()
save_figure(fig, "fig3_candidate_tiers")
plt.close()
print("[figures] Fig 3 saved")

# ═══════════════════════════════════════════════════════════════════════════════
# Fig 4 — Species summary / within-species fallback / genome fallback
# ═══════════════════════════════════════════════════════════════════════════════
def safe_label_col(score_df):
    """Create compact genome/strain labels for within-species Fig. 4."""
    if "accession" in score_df.columns:
        base = score_df["accession"].astype(str)
    elif "genome" in score_df.columns:
        base = score_df["genome"].astype(str)
    else:
        base = pd.Series([f"Genome {i+1}" for i in range(len(score_df))], index=score_df.index)

    if "strain" in score_df.columns:
        strain = score_df["strain"].astype(str).str.strip()
        strain = strain.replace({"": "", "nan": "", "None": ""})
        return [f"{b}\n{s}" if s else str(b) for b, s in zip(base, strain)]
    return base.astype(str).tolist()


def italic_species_name(name):
    """Return a matplotlib title-safe italic species name when possible."""
    name = str(name).strip()
    parts = name.split()
    if len(parts) >= 2 and parts[0].lower() not in {"unknown", "nan", "none"}:
        genus, species = parts[0], parts[1]
        rest = " ".join(parts[2:])
        out = rf"$\it{{{genus}}}$ $\it{{{species}}}$"
        if rest:
            out += f" {rest}"
        return out
    return name


def draw_single_species_fig4(score_df, sp_summary):
    """Fig. 4 mode for a valid dataset containing only one species."""
    species_name = sp_summary["full_species"].iloc[0] if len(sp_summary) else "single species"
    plot = score_df.copy().sort_values("LAB_score_v1", ascending=True)
    if len(plot) > 30:
        # Keep the figure readable. If the single-species dataset is very large,
        # show the strongest 30 genomes and summarize all tiers in panel B.
        plot = score_df.sort_values("LAB_score_v1", ascending=False).head(30).sort_values("LAB_score_v1", ascending=True)

    labels = wrap_labels(safe_label_col(plot), width=26)
    colors = [TIER_COLORS.get(t, "#999999") for t in plot.get("Candidate_tier_v1_1", pd.Series("", index=plot.index))]
    fig_h = max(7.5, 0.34 * len(plot) + 3.0)

    fig, axes = plt.subplots(
        1, 2,
        figsize=(18, fig_h),
        gridspec_kw={"width_ratios": [1.45, 1.0]}
    )

    # Panel A — genome/strain ranking within the species
    ax = axes[0]
    y = np.arange(len(plot))
    ax.barh(y, plot["LAB_score_v1"], color=colors, edgecolor="white", height=0.72)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7.5)
    ax.set_xlim(0, 103)
    ax.set_xlabel("LAB-Score v1.1")
    ax.set_title(f"A. Within-species genome ranking\n{italic_species_name(species_name)}", fontweight="bold")
    for cutoff in [50, 70, 85]:
        ax.axvline(cutoff, color="#666666", linestyle="--", linewidth=1, alpha=0.55)
    for yi, score in enumerate(plot["LAB_score_v1"]):
        ax.text(score + 0.8, yi, f"{score:.1f}", va="center", fontsize=7.5)
    ax.grid(axis="x", alpha=0.25)
    ax.set_axisbelow(True)

    # Panel B — all candidate tiers in the single-species dataset
    ax2 = axes[1]
    tier_counts2 = score_df["Candidate_tier_v1_1"].value_counts() if "Candidate_tier_v1_1" in score_df.columns else pd.Series(dtype=int)
    plot_tiers2 = [t for t in TIER_ORDER if t in tier_counts2.index]
    plot_tiers2 += [t for t in tier_counts2.index if t not in plot_tiers2]

    if len(plot_tiers2) == 0:
        ax2.text(0.5, 0.5, "Candidate-tier data not available", ha="center", va="center", transform=ax2.transAxes)
        ax2.axis("off")
    else:
        vals2 = [int(tier_counts2[t]) for t in plot_tiers2]
        colors2 = [TIER_COLORS.get(t, "#999999") for t in plot_tiers2]
        bars = ax2.barh(plot_tiers2, vals2, color=colors2, edgecolor="white", height=0.65)
        for bar, val in zip(bars, vals2):
            ax2.text(val + max(vals2) * 0.02, bar.get_y() + bar.get_height()/2, f"n={val:,}", va="center", fontsize=9)
        ax2.set_xlabel("Number of genomes")
        ax2.set_title("B. Safety-gated candidate tiers", fontweight="bold")
        ax2.set_xlim(0, max(vals2) * 1.35)
        ax2.grid(axis="x", alpha=0.25)
        ax2.set_axisbelow(True)

    fig.suptitle(
        f"Figure 4. Within-species LAB-Score Summary for {italic_species_name(species_name)}",
        fontsize=14,
        fontweight="bold",
        y=0.98,
    )
    fig.text(
        0.5, 0.012,
        "This dataset contains one resolved species; therefore Fig. 4 summarizes genome-level prioritization within that species rather than between-species comparisons.",
        ha="center",
        fontsize=9,
        color="#555555",
    )
    plt.tight_layout(rect=[0, 0.04, 1, 0.94])
    save_figure(fig, "fig4_species_summary")
    plt.close()
    print("[figures] Fig 4 saved — within-species mode")


def draw_multi_species_fig4(sp_summary):
    """Fig. 4 mode for datasets with two or more resolved species."""
    sp_plot = sp_summary[sp_summary["N"] >= 2].copy()
    if sp_plot.empty:
        sp_plot = sp_summary.copy()
    sp_plot = sp_plot.sort_values(["mean_LAB_score", "N"], ascending=[False, False]).head(20)

    fig, axes = plt.subplots(1, 2, figsize=(18, max(8, 0.42*len(sp_plot) + 3)), gridspec_kw={"width_ratios":[1.05,1.35]})

    y = np.arange(len(sp_plot))
    ylabels = wrap_labels(sp_plot["full_species"], 28)

    # Panel A: score range with mean point and N labels.
    ax = axes[0]
    ax.hlines(y, sp_plot["min_LAB_score"], sp_plot["max_LAB_score"], color="#95D5B2", lw=5, alpha=0.85, zorder=1)
    size = 60 + (sp_plot["N"] / sp_plot["N"].max()) * 260
    ax.scatter(sp_plot["mean_LAB_score"], y, s=size, color="#1B4332", edgecolor="white", linewidth=1.2, zorder=3)
    for yi, row in enumerate(sp_plot.itertuples(index=False)):
        ax.text(row.max_LAB_score + 1.0, yi, f"n={int(row.N)}", va="center", fontsize=8.5, color="#333333")
    for v in [50, 70, 85]:
        ax.axvline(v, color="#999999", lw=1.0, ls="--", alpha=0.55)
    ax.set_yticks(y)
    ax.set_yticklabels(ylabels, fontsize=8.5)
    ax.invert_yaxis()
    ax.set_xlim(0, 103)
    ax.set_xlabel("LAB-Score v1.1")
    ax.set_title("A. Top species by mean LAB-Score\nline = min–max; dot = mean; size = n", fontweight="bold")
    ax.grid(axis="x", alpha=0.25)

    # Panel B: stacked candidate-tier counts showing every final category that exists.
    ax2 = axes[1]
    sp_s = sp_plot.set_index("full_species")
    tier_cols = [t for t in TIER_ORDER if t in sp_s.columns and sp_s[t].sum() > 0]
    tier_cols += [t for t in sp_s.columns if t in tier_counts.index and t not in tier_cols and sp_s[t].sum() > 0]
    bottom = np.zeros(len(sp_s))
    for col in tier_cols:
        vals2 = sp_s[col].fillna(0).values.astype(float)
        ax2.barh(y, vals2, left=bottom, color=TIER_COLORS.get(col, "#999999"),
                 label=col, edgecolor="white", height=0.72)
        bottom += vals2
    for yi, total in enumerate(bottom):
        ax2.text(total + max(bottom.max()*0.01, 0.15), yi, f"{int(total):,}", va="center", fontsize=8.5)
    ax2.set_yticks(y)
    ax2.set_yticklabels(ylabels, fontsize=8.5)
    ax2.invert_yaxis()
    ax2.set_xlabel("Number of genomes")
    ax2.set_title("B. Final safety-gated candidate tiers by species", fontweight="bold")
    ax2.grid(axis="x", alpha=0.25)
    ax2.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False, fontsize=8.5)

    plt.suptitle("Figure 4. Species-level LAB-Score and Safety-Gated Candidate Tier Summary",
                 fontsize=14, fontweight="bold", y=1.01)
    plt.tight_layout()
    save_figure(fig, "fig4_species_summary")
    plt.close()
    print("[figures] Fig 4 saved — species-level mode")


def draw_unresolved_species_fig4(score_df):
    """Fig. 4 mode when species metadata are unresolved."""
    top = score_df.sort_values("LAB_score_v1", ascending=False).head(25).copy()
    labels = wrap_labels(safe_label_col(top), 30)
    colors = [TIER_COLORS.get(t, "#999999") for t in top.get("Candidate_tier_v1_1", pd.Series("", index=top.index))]

    fig, axes = plt.subplots(1, 2, figsize=(18, 9), gridspec_kw={"width_ratios":[1.3,1]})
    ax = axes[0]
    y = np.arange(len(top))
    ax.barh(y, top["LAB_score_v1"], color=colors, edgecolor="white", height=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7.5)
    ax.invert_yaxis()
    for v in [50,70,85]:
        ax.axvline(v, color="#999999", lw=1.0, ls="--", alpha=0.6)
    for yi, score in enumerate(top["LAB_score_v1"]):
        ax.text(score + 0.8, yi, f"{score:.1f}", va="center", fontsize=8)
    ax.set_xlim(0, 103)
    ax.set_xlabel("LAB-Score v1.1")
    ax.set_title("A. Top-ranked genomes by LAB-Score", fontweight="bold")
    ax.grid(axis="x", alpha=0.25)

    ax2 = axes[1]
    tier_counts2 = score_df["Candidate_tier_v1_1"].value_counts() if "Candidate_tier_v1_1" in score_df.columns else pd.Series(dtype=int)
    plot_tiers2 = [t for t in TIER_ORDER if t in tier_counts2.index]
    plot_tiers2 += [t for t in tier_counts2.index if t not in plot_tiers2]
    vals2 = [int(tier_counts2[t]) for t in plot_tiers2]
    colors2 = [TIER_COLORS.get(t, "#999999") for t in plot_tiers2]
    if len(vals2):
        bars = ax2.barh(plot_tiers2, vals2, color=colors2, edgecolor="white", height=0.65)
        for bar, val in zip(bars, vals2):
            ax2.text(val + max(vals2)*0.01, bar.get_y()+bar.get_height()/2, f"n={val:,}", va="center", fontsize=9)
        ax2.set_xlim(0, max(vals2) * 1.25)
    else:
        ax2.text(0.5, 0.5, "Candidate-tier data not available", ha="center", va="center", transform=ax2.transAxes)
    ax2.set_xlabel("Number of genomes")
    ax2.set_title("B. Final candidate-tier distribution", fontweight="bold")
    ax2.grid(axis="x", alpha=0.25)

    note = "Species metadata were not resolved; Fig. 4 uses a genome-level fallback summary."
    fig.text(0.5, 0.015, note, ha="center", fontsize=9.5, color="#666666")
    plt.suptitle("Figure 4. Genome-level LAB-Score Summary When Species Metadata Are Unresolved",
                 fontsize=14, fontweight="bold", y=0.98)
    plt.tight_layout(rect=[0, 0.04, 1, 0.95])
    save_figure(fig, "fig4_species_summary")
    plt.close()
    print("[figures] Fig 4 saved — unresolved-species fallback mode")


species_resolved = species_metadata_is_resolved(df)
sp_summary = build_species_summary(df)
n_species = sp_summary.shape[0]

# Save the rebuilt species summary for inspection/reproducibility before plotting.
sp_summary.sort_values(["mean_LAB_score", "N"], ascending=[False, False]).to_csv(
    os.path.join(args.outdir, "fig4_rebuilt_species_summary.tsv"), sep="\t", index=False
)

if species_resolved and n_species >= 2:
    draw_multi_species_fig4(sp_summary)
elif species_resolved and n_species == 1:
    draw_single_species_fig4(df, sp_summary)
else:
    draw_unresolved_species_fig4(df)

# ═══════════════════════════════════════════════════════════════════════════════
# Fig 5 — ML Feature Importance
# ═══════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(18, 8))

for ax, fi_df, title in [
    (axes[0], fi_reg, "Regression\n(Functional modules → LAB-Score)"),
    (axes[1], fi_cls, "Classification\n(Functional modules → Priority class)"),
]:
    if fi_df.empty or len(fi_df) == 0:
        ax.text(0.5, 0.5, "ML data not available\n(need ≥10 genomes)",
                ha="center", va="center", transform=ax.transAxes, fontsize=11, color="#666")
        ax.axis("off")
        continue

    fi_plot = fi_df.head(16).set_index("feature")["importance"].sort_values(ascending=True)
    med     = fi_plot.median()
    colors  = ["#1B4332" if v >= fi_plot.quantile(0.75) else
               "#52B788" if v >= med else "#B7E4C7"
               for v in fi_plot.values]

    ax.barh(fi_plot.index, fi_plot.values, color=colors, edgecolor="white", height=0.7, zorder=2)
    ax.axvline(med, color="#E07C24", ls="--", lw=1.5, label="Median")
    for i, val in enumerate(fi_plot.values):
        ax.text(val + fi_plot.max()*0.01, i, f"{val:.3f}", va="center", fontsize=8.5)

    ax.set_title(title, fontweight="bold", fontsize=12)
    ax.set_xlabel("Feature Importance (Mean Decrease Impurity)")
    ax.legend(fontsize=9, frameon=False)
    ax.grid(axis="x", alpha=0.3, zorder=1)
    ax.set_axisbelow(True)

plt.suptitle("Figure 5. Machine Learning Feature Importance — LAB-Score v1.1",
             fontsize=13, fontweight="bold", y=1.02)
plt.tight_layout()
save_figure(fig, "fig5_ml_feature_importance")
plt.close()
print("[figures] Fig 5 saved")

print(f"\n[figures] All figures saved to: {args.outdir}")
