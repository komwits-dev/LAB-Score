#!/usr/bin/env python3
"""
16_make_html_summary_report_portable.py

Portable manuscript-style LAB-Score v1.1 HTML summary report.

Key fixes:
- Copies figures into 16_html_report/assets so the report remains usable after transfer.
- Reads ML outputs from 10_ml.
- Separates raw score tiers from safety-gated candidate labels.
- Uses compact, interpretable tables rather than all master-matrix columns.
- Shows explicit warnings for missing outputs and unusually high caution rates.
"""

import argparse
import html
import shutil
from pathlib import Path

import pandas as pd


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pipeline-outdir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--title", default="LAB-Score v1.1 Summary Report")
    return ap.parse_args()


def read_tsv(base: Path, rel: str):
    path = base / rel
    if not path.is_file():
        return None
    try:
        return pd.read_csv(path, sep="\t")
    except Exception as exc:
        print(f"[html] WARNING: could not read {path}: {exc}")
        return None


def table_html(df, columns=None, max_rows=30, empty_message="No table available."):
    if df is None or df.empty:
        return f"<p class='warn'>{html.escape(empty_message)}</p>"

    show = df.copy()
    if columns:
        existing = [c for c in columns if c in show.columns]
        if existing:
            show = show[existing]

    show = show.head(max_rows).copy()
    for col in show.columns:
        if pd.api.types.is_float_dtype(show[col]):
            show[col] = show[col].map(
                lambda value: "" if pd.isna(value) else f"{value:.3f}"
            )
        else:
            show[col] = show[col].fillna("")

    return (
        "<div class='table-wrap'>"
        + show.to_html(index=False, classes="tbl", border=0, escape=True)
        + "</div>"
    )


def count_equals(df, col, value):
    if df is None or col not in df.columns:
        return 0
    return int(df[col].astype(str).eq(value).sum())


def count_contains(df, col, text):
    if df is None or col not in df.columns:
        return 0
    return int(df[col].astype(str).str.contains(text, case=False, na=False).sum())


def copy_figure(base: Path, outdir: Path, title: str, candidates):
    source = None
    for rel in candidates:
        path = base / rel
        if path.is_file():
            source = path
            break

    if source is None:
        expected = " or ".join(str(base / rel) for rel in candidates)
        return (
            "<div class='fig missing'>"
            f"<h3>{html.escape(title)}</h3>"
            f"<p>Missing figure: {html.escape(expected)}</p>"
            "</div>"
        )

    assets = outdir / "assets"
    assets.mkdir(parents=True, exist_ok=True)
    destination = assets / source.name
    shutil.copy2(source, destination)

    return (
        "<div class='fig'>"
        f"<h3>{html.escape(title)}</h3>"
        f"<img src='assets/{html.escape(destination.name)}' "
        f"alt='{html.escape(title)}'>"
        "</div>"
    )


def metric_value(df, model, metric):
    if df is None or not {"Model", "Metric", "Value"}.issubset(df.columns):
        return "NA"
    selected = df.loc[
        df["Model"].astype(str).eq(model)
        & df["Metric"].astype(str).eq(metric),
        "Value",
    ]
    if selected.empty:
        return "NA"
    try:
        return f"{float(selected.iloc[0]):.3f}"
    except Exception:
        return html.escape(str(selected.iloc[0]))


args = parse_args()
base = Path(args.pipeline_outdir).expanduser().resolve()
outdir = Path(args.outdir).expanduser().resolve()
outdir.mkdir(parents=True, exist_ok=True)
out_html = outdir / "LAB_SCORE_v1_1_summary_report.html"

scores = read_tsv(base, "09_scores/LAB_score_v1.tsv")
priority = read_tsv(base, "09_scores/priority_class_counts.tsv")
species_file = read_tsv(base, "09_scores/species_mean_scores.tsv")
top100 = read_tsv(base, "09_scores/top100_candidates.tsv")

# Build a robust species summary directly from the final score table.
# This avoids dependence on column naming differences in species_mean_scores.tsv.
species = None
if scores is not None and not scores.empty:
    species_col = None
    for candidate in ["full_species", "Species", "species"]:
        if candidate in scores.columns:
            species_col = candidate
            break

    if species_col is not None and "LAB_score_v1" in scores.columns:
        work = scores.copy()
        work[species_col] = work[species_col].fillna("Unknown sp.").astype(str)

        species = (
            work.groupby(species_col, dropna=False)
            .agg(
                n_genomes=("LAB_score_v1", "size"),
                mean_LAB_score=("LAB_score_v1", "mean"),
                median_LAB_score=("LAB_score_v1", "median"),
                max_LAB_score=("LAB_score_v1", "max"),
            )
            .reset_index()
            .rename(columns={species_col: "Species"})
        )

        if "Priority_class" in work.columns:
            elite_counts = (
                work.assign(_elite=work["Priority_class"].astype(str).eq("Elite"))
                .groupby(species_col)["_elite"].sum()
            )
            high_counts = (
                work.assign(_high=work["Priority_class"].astype(str).eq("High"))
                .groupby(species_col)["_high"].sum()
            )
            species["n_elite_tier"] = species["Species"].map(elite_counts).fillna(0).astype(int)
            species["n_high_tier"] = species["Species"].map(high_counts).fillna(0).astype(int)
        else:
            species["n_elite_tier"] = 0
            species["n_high_tier"] = 0

        if "Refined_safety_status" in work.columns:
            critical_counts = (
                work.assign(
                    _critical=work["Refined_safety_status"]
                    .astype(str)
                    .str.contains("critical", case=False, na=False)
                )
                .groupby(species_col)["_critical"].sum()
            )
            species["n_critical_review"] = (
                species["Species"].map(critical_counts).fillna(0).astype(int)
            )
        else:
            species["n_critical_review"] = 0

        species = species.sort_values(
            ["mean_LAB_score", "n_genomes"],
            ascending=[False, False],
        ).reset_index(drop=True)

# Fall back to the precomputed species table only when a score-derived
# summary cannot be created.
if species is None:
    species = species_file
ml_perf = read_tsv(base, "10_ml/model_performance_summary.tsv")
ml_feat = read_tsv(base, "10_ml/RF_classification_functional_feature_importance.tsv")
sensitivity = read_tsv(base, "13_sensitivity/weight_sensitivity_summary.tsv")
baselines = read_tsv(base, "14_baselines/baseline_comparison_summary.tsv")
manifest = read_tsv(base, "15_release/LAB_SCORE_v1_1_release_manifest.tsv")

n_genomes = len(scores) if scores is not None else "NA"
raw_elite = count_equals(scores, "Priority_class", "Elite")
raw_high = count_equals(scores, "Priority_class", "High")

pass_n = 0
caution_n = 0
critical_n = 0
if scores is not None and "Refined_safety_status" in scores.columns:
    safety = scores["Refined_safety_status"].astype(str)
    pass_n = int(safety.isin(["Pass", "No safety marker detected"]).sum())
    caution_n = int(safety.isin(["Caution", "Cautionary safety review"]).sum())
    critical_n = int(safety.isin(["Critical", "Critical safety review"]).sum())

safety_clear_high_elite = 0
if scores is not None and "Candidate_tier_v1_1" in scores.columns:
    tier = scores["Candidate_tier_v1_1"].astype(str)
    safety_clear_high_elite = int(
        tier.isin(["Elite candidate", "High candidate"]).sum()
    )

func_acc = metric_value(ml_perf, "RF_Class_FUNC", "Accuracy")
func_r2 = metric_value(ml_perf, "RF_Regression_FUNC", "R2")

tier_counts = None
if scores is not None and "Candidate_tier_v1_1" in scores.columns:
    tier_counts = (
        scores["Candidate_tier_v1_1"]
        .fillna("Unclassified")
        .value_counts()
        .rename_axis("Candidate_tier_v1_1")
        .reset_index(name="n_genomes")
    )
    tier_counts["percentage"] = (
        100.0 * tier_counts["n_genomes"] / len(scores)
    )

safety_counts = None
if scores is not None and "Refined_safety_status" in scores.columns:
    safety_counts = (
        scores["Refined_safety_status"]
        .fillna("Unclassified")
        .value_counts()
        .rename_axis("Refined_safety_status")
        .reset_index(name="n_genomes")
    )
    safety_counts["percentage"] = (
        100.0 * safety_counts["n_genomes"] / len(scores)
    )

caution_warning = ""
if isinstance(n_genomes, int) and n_genomes > 0 and caution_n / n_genomes >= 0.80:
    caution_warning = (
        "<div class='warning'><strong>Safety-screen review recommended:</strong> "
        f"{caution_n:,} of {n_genomes:,} genomes "
        f"({100*caution_n/n_genomes:.1f}%) were assigned a cautionary status. "
        "Before manuscript interpretation, inspect whether broad annotation terms "
        "such as glutamate decarboxylase, generic PLP-dependent decarboxylase, or "
        "hemolysin-III-family proteins are being treated as adverse markers.</div>"
    )

fig_overview = "".join(
    [
        copy_figure(
            base,
            outdir,
            "Figure 1. LAB-Score workflow",
            [
                "11_figures/fig1_workflow_overview.svg",
                "11_figures/lab_score_workflow.svg",
                "11_figures/fig1_workflow_overview.png",
                "11_figures/lab_score_workflow_overview_diagram.png",
            ],
        ),
        copy_figure(
            base,
            outdir,
            "Figure 2. LAB-Score distribution",
            ["11_figures/fig2_score_distribution.png"],
        ),
        copy_figure(
            base,
            outdir,
            "Figure 3. Refined candidate tiers",
            ["11_figures/fig3_candidate_tiers.png"],
        ),
        copy_figure(
            base,
            outdir,
            "Figure 4. Species summary",
            ["11_figures/fig4_species_summary.png"],
        ),
        copy_figure(
            base,
            outdir,
            "Figure 5. Functional drivers",
            [
                "11_figures/fig5_ml_feature_importance.png",
                "10_ml/fig_ml_feature_importance.png",
            ],
        ),
    ]
)

fig_ml = "".join(
    [
        copy_figure(
            base,
            outdir,
            "ML feature importance",
            ["10_ml/fig_ml_feature_importance.png"],
        ),
        copy_figure(
            base,
            outdir,
            "ML actual versus predicted",
            ["10_ml/fig_ml_actual_vs_predicted.png"],
        ),
        copy_figure(
            base,
            outdir,
            "ML confusion matrix",
            ["10_ml/fig_ml_confusion_matrix.png"],
        ),
    ]
)

score_preview_columns = [
    "accession",
    "full_species",
    "strain",
    "LAB_score_v1",
    "Priority_class",
    "Refined_safety_status",
    "Candidate_tier_v1_1",
    "Safety_score",
    "GI_survival_score",
    "Functional_score",
    "Fermentation_score",
    "amrfinder_hits",
    "vfdb_hits",
    "resfinder_hits",
    "biogenic_amines",
    "hemolysin_markers",
    "reasons",
]

top_columns = [
    "accession",
    "full_species",
    "strain",
    "LAB_score_v1",
    "Priority_class",
    "Refined_safety_status",
    "Candidate_tier_v1_1",
    "Safety_score",
    "GI_survival_score",
    "Functional_score",
    "Fermentation_score",
    "reasons",
]

species_columns = [
    "Species",
    "n_genomes",
    "mean_LAB_score",
    "median_LAB_score",
    "max_LAB_score",
    "n_elite_tier",
    "n_high_tier",
    "n_critical_review",
]

file_rows = []
for rel in [
    "09_scores/LAB_score_v1.tsv",
    "09_scores/top100_candidates.tsv",
    "09_scores/species_mean_scores.tsv",
    "10_ml/model_performance_summary.tsv",
    "10_ml/RF_classification_functional_feature_importance.tsv",
    "13_sensitivity/weight_sensitivity_summary.tsv",
    "14_baselines/baseline_comparison_summary.tsv",
    "15_release/LAB_SCORE_v1_1_release_manifest.tsv",
    "12_report/LAB_score_report.html",
]:
    p = base / rel
    file_rows.append(
        {
            "status": "Available" if p.exists() else "Missing",
            "path": rel,
            "size_bytes": p.stat().st_size if p.is_file() else "",
        }
    )
files_df = pd.DataFrame(file_rows)

html_text = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(args.title)}</title>
<style>
:root {{
  --navy:#10375c; --blue:#2274a5; --teal:#1b998b; --orange:#d97706;
  --red:#b42318; --bg:#f4f7fb; --card:#ffffff; --text:#1f2933;
  --muted:#64748b; --border:#d9e2ec;
}}
* {{ box-sizing:border-box; }}
body {{ margin:0; font-family:Arial,Helvetica,sans-serif; background:var(--bg); color:var(--text); }}
header {{ padding:32px 40px; color:white; background:linear-gradient(135deg,var(--navy),var(--blue),var(--teal)); }}
header h1 {{ margin:0; font-size:34px; }}
header p {{ margin:8px 0 0; font-size:16px; }}
main {{ padding:24px 36px 48px; }}
.cards {{ display:grid; grid-template-columns:repeat(6,minmax(150px,1fr)); gap:14px; margin-bottom:22px; }}
.card {{ background:var(--card); border:1px solid var(--border); border-radius:14px; padding:16px; box-shadow:0 2px 7px rgba(0,0,0,.05); }}
.label {{ color:var(--muted); font-size:12px; font-weight:700; text-transform:uppercase; letter-spacing:.5px; }}
.value {{ margin-top:8px; color:var(--navy); font-size:28px; font-weight:800; }}
.tabs {{ display:flex; flex-wrap:wrap; gap:8px; margin-bottom:18px; }}
.tabbtn {{ border:1px solid var(--border); background:white; color:var(--navy); padding:11px 15px; border-radius:999px; font-weight:700; cursor:pointer; }}
.tabbtn.active {{ background:var(--navy); color:white; }}
.tab {{ display:none; background:white; border:1px solid var(--border); border-radius:16px; padding:24px; box-shadow:0 2px 7px rgba(0,0,0,.05); }}
.tab.active {{ display:block; }}
h2,h3 {{ color:var(--navy); }}
h2 {{ margin-top:0; }}
.note {{ background:#eef7ff; border-left:5px solid var(--blue); padding:13px 16px; border-radius:10px; margin:14px 0 20px; }}
.warning {{ background:#fff7e6; border-left:5px solid var(--orange); padding:13px 16px; border-radius:10px; margin:14px 0 20px; }}
.warn,.missing {{ color:var(--red); font-weight:700; }}
.figgrid {{ display:grid; grid-template-columns:repeat(2,minmax(300px,1fr)); gap:18px; }}
.fig {{ background:white; border:1px solid var(--border); border-radius:14px; padding:14px; text-align:center; }}
.fig img {{ width:100%; height:auto; border-radius:8px; }}
.table-wrap {{ width:100%; overflow-x:auto; margin-top:14px; }}
.tbl {{ width:100%; border-collapse:collapse; font-size:13px; white-space:nowrap; }}
.tbl th {{ background:var(--navy); color:white; padding:9px; text-align:left; }}
.tbl td {{ padding:8px 9px; border-bottom:1px solid var(--border); }}
.tbl tr:nth-child(even) {{ background:#f8fafc; }}
footer {{ text-align:center; color:var(--muted); margin-top:28px; font-size:13px; }}
@media (max-width:1200px) {{ .cards {{ grid-template-columns:repeat(3,1fr); }} }}
@media (max-width:800px) {{ main {{ padding:18px; }} .cards {{ grid-template-columns:repeat(2,1fr); }} .figgrid {{ grid-template-columns:1fr; }} }}
</style>
</head>
<body>
<header>
  <h1>{html.escape(args.title)}</h1>
  <p>Genome-based, safety-aware prioritization report for lactic acid bacteria candidates</p>
</header>
<main>
  <div class="cards">
    <div class="card"><div class="label">Genomes scored</div><div class="value">{n_genomes}</div></div>
    <div class="card"><div class="label">Raw Elite score tier</div><div class="value">{raw_elite}</div></div>
    <div class="card"><div class="label">Raw High score tier</div><div class="value">{raw_high}</div></div>
    <div class="card"><div class="label">No safety marker / Pass</div><div class="value">{pass_n}</div></div>
    <div class="card"><div class="label">Safety-clear High + Elite</div><div class="value">{safety_clear_high_elite}</div></div>
    <div class="card"><div class="label">Critical review</div><div class="value">{critical_n}</div></div>
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
    <p>LAB-Score integrates safety screening, gastrointestinal-survival modules,
    probiotic-associated functions, and fermentation-related features into a continuous
    score and safety-gated candidate designation.</p>
    <div class="note"><strong>Interpretation:</strong> LAB-Score is a decision-support
    framework for selecting genomes for downstream experimental validation. It is not
    definitive proof of probiotic efficacy or safety.</div>
    {caution_warning}
    <div class="figgrid">{fig_overview}</div>
  </section>

  <section id="tiers" class="tab">
    <h2>Candidate tiers</h2>
    <p>Raw score tiers and safety-gated candidate labels must not be interpreted as the same quantity.</p>
    {table_html(tier_counts, max_rows=30)}
    <h3>Compact score preview</h3>
    {table_html(scores, score_preview_columns, max_rows=25)}
  </section>

  <section id="safety" class="tab">
    <h2>Safety status</h2>
    {table_html(safety_counts, max_rows=20)}
    {caution_warning}
  </section>

  <section id="top" class="tab">
    <h2>Top candidates</h2>
    <p>Highest LAB-Score genomes. Safety status and reasons must be reviewed alongside rank.</p>
    {table_html(top100, top_columns, max_rows=50)}
  </section>

  <section id="species" class="tab">
    <h2>Species summary</h2>
    {table_html(species, species_columns, max_rows=60)}
  </section>

  <section id="ml" class="tab">
    <h2>Functional drivers / ML</h2>
    <p>Functional-only Random Forest results are presented as an interpretability layer,
    not as independent phenotype validation.</p>
    <p><strong>Functional RF classification accuracy:</strong> {func_acc}
    &nbsp;&nbsp; <strong>Functional RF regression R²:</strong> {func_r2}</p>
    <h3>Model performance</h3>
    {table_html(ml_perf, max_rows=20)}
    <h3>Functional feature importance</h3>
    {table_html(ml_feat, ["feature", "importance"], max_rows=30)}
    <div class="figgrid">{fig_ml}</div>
  </section>

  <section id="sensitivity" class="tab">
    <h2>Weight sensitivity analysis</h2>
    {table_html(sensitivity, max_rows=30)}
  </section>

  <section id="baseline" class="tab">
    <h2>Baseline comparison</h2>
    {table_html(baselines, max_rows=30)}
  </section>

  <section id="files" class="tab">
    <h2>Output files</h2>
    {table_html(files_df, max_rows=30)}
    <h3>Release manifest</h3>
    {table_html(manifest, max_rows=50)}
  </section>

  <footer>LAB-Score v1.1 portable HTML report generated from pipeline outputs.</footer>
</main>
<script>
function openTab(evt, name) {{
  const tabs = document.getElementsByClassName("tab");
  for (let i=0; i<tabs.length; i++) tabs[i].classList.remove("active");
  const buttons = document.getElementsByClassName("tabbtn");
  for (let i=0; i<buttons.length; i++) buttons[i].classList.remove("active");
  document.getElementById(name).classList.add("active");
  evt.currentTarget.classList.add("active");
}}
</script>
</body>
</html>
"""

out_html.write_text(html_text, encoding="utf-8")
print(f"[html] Saved: {out_html}")
print(f"[html] Assets: {outdir / 'assets'}")
