#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd
import html


def read_table(path, sep="\t"):
    path = Path(path)
    if path.exists():
        return pd.read_csv(path, sep=sep)
    return None


def first_existing_column(df, candidates):
    if df is None:
        return None
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
    return mapping.get(x, x)


def harmonize_safety_status(x):
    x = str(x).strip()
    mapping = {
        "Pass": "No safety marker detected",
        "Review": "Safety review",
        "Critical": "Critical safety review",
        "Caution": "Cautionary safety review",
    }
    return mapping.get(x, x)


def percent(n, total):
    if total == 0:
        return 0
    return round(n / total * 100, 2)


def make_candidate_tier_summary(score):
    tier_col = first_existing_column(
        score,
        ["Candidate_tier_v1_1", "Candidate_tier", "candidate_tier", "Priority_class", "priority_class"],
    )

    if tier_col is None:
        return pd.DataFrame(columns=["Candidate_tier", "n_genomes", "percentage"])

    tmp = score.copy()
    tmp["Candidate_tier"] = tmp[tier_col].map(normalize_tier_label)

    out = (
        tmp["Candidate_tier"]
        .value_counts()
        .rename_axis("Candidate_tier")
        .reset_index(name="n_genomes")
    )
    out["percentage"] = (out["n_genomes"] / out["n_genomes"].sum() * 100).round(2)

    preferred_order = [
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

    out["order"] = out["Candidate_tier"].apply(
        lambda x: preferred_order.index(x) if x in preferred_order else 99
    )
    out = out.sort_values(["order", "n_genomes"], ascending=[True, False])
    out = out.drop(columns=["order"])

    return out


def make_safety_summary(score):
    status_col = first_existing_column(
        score,
        ["Refined_safety_status", "Safety_gate", "flag", "safety_status"],
    )

    if status_col is None:
        return pd.DataFrame(columns=["Refined_safety_status", "n_genomes", "percentage"])

    tmp = score.copy()
    tmp["Refined_safety_status"] = tmp[status_col].map(harmonize_safety_status)

    out = (
        tmp["Refined_safety_status"]
        .value_counts()
        .rename_axis("Refined_safety_status")
        .reset_index(name="n_genomes")
    )
    out["percentage"] = (out["n_genomes"] / out["n_genomes"].sum() * 100).round(2)

    preferred_order = [
        "No safety marker detected",
        "Cautionary safety review",
        "Critical safety review",
        "Safety review",
    ]
    out["order"] = out["Refined_safety_status"].apply(
        lambda x: preferred_order.index(x) if x in preferred_order else 99
    )
    out = out.sort_values(["order", "n_genomes"], ascending=[True, False])
    out = out.drop(columns=["order"])

    return out


def make_species_summary(score):
    genus_col = first_existing_column(score, ["genus", "Genus"])
    species_col = first_existing_column(score, ["species", "Species_x", "Species"])
    full_col = first_existing_column(score, ["full_species", "Species_full"])

    if full_col is None:
        if genus_col is not None and species_col is not None:
            score = score.copy()
            score["full_species_auto"] = score[genus_col].astype(str) + " " + score[species_col].astype(str)
            full_col = "full_species_auto"
        else:
            return None

    tier_col = first_existing_column(score, ["Candidate_tier_v1_1", "Candidate_tier", "Priority_class"])
    if tier_col is None:
        return None

    tmp = score.copy()
    tmp["Candidate_tier_clean"] = tmp[tier_col].map(normalize_tier_label)

    out = (
        tmp.groupby(full_col)
        .agg(
            n_genomes=("LAB_score_v1", "count"),
            mean_LAB_score=("LAB_score_v1", "mean"),
            median_LAB_score=("LAB_score_v1", "median"),
            max_LAB_score=("LAB_score_v1", "max"),
            n_elite_tier=("Candidate_tier_clean", lambda x: x.astype(str).str.contains("Elite", case=False).sum()),
            n_high_tier=("Candidate_tier_clean", lambda x: x.astype(str).str.contains("High", case=False).sum()),
            n_critical_review=("Candidate_tier_clean", lambda x: (x == "Critical safety review").sum()),
        )
        .reset_index()
        .rename(columns={full_col: "Species"})
        .sort_values(["n_elite_tier", "n_high_tier", "mean_LAB_score"], ascending=[False, False, False])
    )

    for c in ["mean_LAB_score", "median_LAB_score", "max_LAB_score"]:
        out[c] = out[c].round(3)

    return out


def df_html(df, max_rows=50):
    if df is None or len(df) == 0:
        return "<p class='warn'>No table available.</p>"

    show = df.head(max_rows).copy()

    for c in show.columns:
        if pd.api.types.is_float_dtype(show[c]):
            show[c] = show[c].round(3)

    return show.to_html(index=False, classes="tbl", border=0)


def find_first_existing(paths):
    for p in paths:
        p = Path(p)
        if p.exists():
            return p
    return None


def img_html(path, title):
    if path is None or not Path(path).exists():
        return f"<div class='fig missing'>Missing figure: {html.escape(title)}</div>"
    return f"""
    <div class='fig'>
      <h3>{html.escape(title)}</h3>
      <img src='{html.escape(str(path))}'>
    </div>
    """


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pipeline-outdir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--title", default="LAB-Score v1.1 Summary Report")
    args = ap.parse_args()

    root = Path(args.pipeline_outdir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    score_path = root / "09_scores" / "LAB_score_v1.tsv"
    score = read_table(score_path)

    if score is None:
        raise FileNotFoundError(f"Missing score table: {score_path}")

    n_genomes = len(score)

    tier_col = first_existing_column(score, ["Candidate_tier_v1_1", "Candidate_tier", "Priority_class"])
    tier_series = score[tier_col].map(normalize_tier_label) if tier_col else pd.Series([], dtype=str)

    safety_col = first_existing_column(score, ["Refined_safety_status", "Safety_gate", "flag"])
    safety_series = score[safety_col].map(harmonize_safety_status) if safety_col else pd.Series([], dtype=str)

    n_elite_tier = int(tier_series.astype(str).str.contains("Elite", case=False).sum())
    n_high_tier = int(tier_series.astype(str).str.contains("High", case=False).sum())
    n_no_safety = int((safety_series == "No safety marker detected").sum())
    n_critical = int((tier_series == "Critical safety review").sum())

    tier_summary = make_candidate_tier_summary(score)
    safety_summary = make_safety_summary(score)
    species_summary = make_species_summary(score)

    # Existing ML and robustness outputs
    ml_summary = read_table(root / "12_ml_functional_only" / "RF_functional_only_model_summary.tsv")
    if ml_summary is None:
        ml_summary = read_table(root / "10_ml" / "ML_model_performance_summary.tsv")
    if ml_summary is None:
        ml_summary = read_table(root / "06_results_v1_1" / "Table5_ML_model_performance_summary.tsv")

    func_importance = read_table(root / "12_ml_functional_only" / "RF_functional_only_candidate_group_feature_importance.tsv")
    if func_importance is None:
        func_importance = read_table(root / "10_ml" / "RF_functional_only_candidate_group_feature_importance.tsv")
    if func_importance is None:
        func_importance = read_table(root / "06_results_v1_1" / "Table6_top20_functional_modules_driving_candidate_groups.tsv")

    sensitivity = read_table(root / "13_sensitivity" / "weight_sensitivity_summary.tsv")
    baseline = read_table(root / "14_baselines" / "baseline_comparison_summary.tsv")
    manifest = read_table(root / "15_release" / "release_manifest.tsv")

    rf_acc = "NA"
    if ml_summary is not None:
        acc_col = first_existing_column(ml_summary, ["Classification_accuracy", "accuracy"])
        if acc_col is not None:
            try:
                rf_acc = f"{float(ml_summary[acc_col].iloc[0]):.3f}"
            except Exception:
                rf_acc = "NA"

    # Figures: look for current pipeline figure names
    fig2 = find_first_existing([
        root / "11_figures" / "fig2_score_distribution.png",
        root / "11_figures" / "Figure_LAB_score_distribution_thresholds.png",
        root / "04_figures" / "Figure_LAB_score_distribution_thresholds.png",
    ])

    fig3 = find_first_existing([
        root / "11_figures" / "fig3_candidate_tiers.png",
        root / "11_figures" / "Figure_v1_1_refined_candidate_tiers.png",
        root / "04_figures" / "Figure_safety_gated_candidate_tiers.png",
    ])

    fig5 = find_first_existing([
        root / "11_figures" / "fig5_functional_drivers.png",
        root / "12_ml_functional_only" / "figures" / "Figure_functional_only_candidate_group_feature_importance.png",
        root / "12_ml_functional_only" / "Figure_functional_only_candidate_group_feature_importance.png",
        root / "10_ml" / "Figure_functional_only_candidate_group_feature_importance.png",
        root / "04_figures_ml_functional_only" / "Figure_functional_only_candidate_group_feature_importance.png",
    ])

    css = """
    body { margin:0; font-family:Arial, Helvetica, sans-serif; background:#f4f7fb; color:#1f2933; }
    header { padding:34px 44px; color:white; background:linear-gradient(135deg,#10375c,#2274a5,#1b998b); }
    header h1 { margin:0; font-size:34px; }
    header p { margin:8px 0 0 0; font-size:16px; opacity:0.95; }
    main { padding:24px 42px 50px 42px; }
    .cards { display:grid; grid-template-columns:repeat(6,1fr); gap:14px; margin-bottom:24px; }
    .card { background:white; border:1px solid #d9e2ec; border-radius:14px; padding:16px; box-shadow:0 2px 7px rgba(0,0,0,0.06); }
    .card .label { color:#8492a6; font-size:12px; font-weight:700; text-transform:uppercase; letter-spacing:0.6px; }
    .card .value { margin-top:8px; color:#10375c; font-size:28px; font-weight:800; }
    .tabs { display:flex; flex-wrap:wrap; gap:8px; margin-bottom:18px; }
    .tabbtn { border:1px solid #d9e2ec; background:white; color:#10375c; padding:11px 15px; border-radius:999px; font-weight:700; cursor:pointer; }
    .tabbtn.active { background:#10375c; color:white; }
    .tab { display:none; background:white; border:1px solid #d9e2ec; border-radius:16px; padding:24px; box-shadow:0 2px 7px rgba(0,0,0,0.05); }
    .tab.active { display:block; }
    h2 { color:#10375c; margin-top:0; }
    h3 { color:#10375c; }
    .note { background:#eef7ff; border-left:5px solid #2274a5; padding:13px 16px; border-radius:10px; margin:14px 0 20px 0; }
    .warn { color:#b42318; font-weight:700; }
    .tbl { width:100%; border-collapse:collapse; font-size:13px; margin-top:14px; }
    .tbl th { background:#10375c; color:white; padding:9px; text-align:left; }
    .tbl td { padding:8px 9px; border-bottom:1px solid #d9e2ec; }
    .tbl tr:nth-child(even) { background:#f8fafc; }
    .figgrid { display:grid; grid-template-columns:repeat(2,1fr); gap:18px; }
    .fig { background:white; border:1px solid #d9e2ec; border-radius:14px; padding:14px; text-align:center; }
    .fig img { max-width:100%; border-radius:8px; }
    .missing { color:#b42318; }
    footer { text-align:center; color:#8492a6; margin-top:28px; font-size:13px; }
    @media (max-width:1100px) { .cards { grid-template-columns:repeat(2,1fr); } .figgrid { grid-template-columns:1fr; } }
    """

    out_html = outdir / "LAB_SCORE_v1_1_summary_report.html"

    html_text = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>{html.escape(args.title)}</title>
<style>{css}</style>
</head>
<body>
<header>
  <h1>{html.escape(args.title)}</h1>
  <p>Genome-based, safety-aware prioritization report for lactic acid bacteria candidates</p>
</header>

<main>
  <div class="cards">
    <div class="card"><div class="label">Genomes scored</div><div class="value">{n_genomes}</div></div>
    <div class="card"><div class="label">Elite-tier genomes</div><div class="value">{n_elite_tier}</div></div>
    <div class="card"><div class="label">High-tier genomes</div><div class="value">{n_high_tier}</div></div>
    <div class="card"><div class="label">No safety marker / Pass</div><div class="value">{n_no_safety}</div></div>
    <div class="card"><div class="label">Critical review</div><div class="value">{n_critical}</div></div>
    <div class="card"><div class="label">Functional RF accuracy</div><div class="value">{rf_acc}</div></div>
  </div>

  <div class="tabs">
    <button class="tabbtn active" onclick="openTab(event,'overview')">Overview</button>
    <button class="tabbtn" onclick="openTab(event,'tiers')">Candidate tiers</button>
    <button class="tabbtn" onclick="openTab(event,'safety')">Safety status</button>
    <button class="tabbtn" onclick="openTab(event,'top')">Top candidates</button>
    <button class="tabbtn" onclick="openTab(event,'species')">Species summary</button>
    <button class="tabbtn" onclick="openTab(event,'ml')">Functional drivers / ML</button>
    <button class="tabbtn" onclick="openTab(event,'sensitivity')">Sensitivity</button>
    <button class="tabbtn" onclick="openTab(event,'baseline')">Baselines</button>
    <button class="tabbtn" onclick="openTab(event,'files')">Output files</button>
  </div>

  <section id="overview" class="tab active">
    <h2>Overview</h2>
    <p>LAB-Score v1.1 integrates safety screening, gastrointestinal survival modules, probiotic-associated functions, and fermentation-related features into a continuous score and refined candidate tier.</p>
    <div class="note"><strong>Interpretation:</strong> LAB-Score is a decision-support framework for selecting candidate genomes for downstream experimental validation. It is not definitive proof of probiotic efficacy or safety.</div>
    <div class="figgrid">
      {img_html(fig2, "Figure 2. LAB-Score distribution")}
      {img_html(fig3, "Figure 3. Refined candidate tiers")}
      {img_html(fig5, "Figure 5. Functional drivers")}
    </div>
  </section>

  <section id="tiers" class="tab">
    <h2>Candidate tiers</h2>
    <p>Counts of refined LAB-Score candidate tiers detected in the score table.</p>
    {df_html(tier_summary, max_rows=30)}
    <h3>Score table preview</h3>
    {df_html(score, max_rows=10)}
  </section>

  <section id="safety" class="tab">
    <h2>Safety status</h2>
    <p>Refined safety status based on detected safety-associated markers.</p>
    {df_html(safety_summary, max_rows=20)}
  </section>

  <section id="top" class="tab">
    <h2>Top candidates</h2>
    <p>Top genomes sorted by LAB_score_v1.</p>
    {df_html(score.sort_values("LAB_score_v1", ascending=False), max_rows=50)}
  </section>

  <section id="species" class="tab">
    <h2>Species summary</h2>
    <p>Species-level aggregation of LAB-Score and candidate-tier composition.</p>
    {df_html(species_summary, max_rows=50)}
  </section>

  <section id="ml" class="tab">
    <h2>Functional drivers / ML</h2>
    <p>Functional-module-only Random Forest analysis provides an interpretability layer for LAB-Score candidate grouping.</p>
    <h3>ML model performance</h3>
    {df_html(ml_summary, max_rows=20)}
    <h3>Functional feature importance</h3>
    {df_html(func_importance, max_rows=30)}
  </section>

  <section id="sensitivity" class="tab">
    <h2>Weight sensitivity analysis</h2>
    {df_html(sensitivity, max_rows=30)}
  </section>

  <section id="baseline" class="tab">
    <h2>Baseline comparison</h2>
    {df_html(baseline, max_rows=30)}
  </section>

  <section id="files" class="tab">
    <h2>Output files</h2>
    {df_html(manifest, max_rows=50)}
    <p>Main score table: <code>{html.escape(str(score_path))}</code></p>
  </section>

  <footer>LAB-Score v1.1 HTML report generated from pipeline outputs.</footer>
</main>

<script>
function openTab(evt, name) {{
  let tabs = document.getElementsByClassName("tab");
  for (let i=0; i<tabs.length; i++) tabs[i].classList.remove("active");
  let buttons = document.getElementsByClassName("tabbtn");
  for (let i=0; i<buttons.length; i++) buttons[i].classList.remove("active");
  document.getElementById(name).classList.add("active");
  evt.currentTarget.classList.add("active");
}}
</script>
</body>
</html>
"""

    out_html.write_text(html_text, encoding="utf-8")

    print("[html] Saved:", out_html)
    print("[html] Genomes:", n_genomes)
    print("[html] Elite-tier genomes:", n_elite_tier)
    print("[html] High-tier genomes:", n_high_tier)
    print("[html] No safety marker:", n_no_safety)
    print("[html] Critical review:", n_critical)


if __name__ == "__main__":
    main()
