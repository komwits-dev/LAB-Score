#!/usr/bin/env python3
"""
15_make_release_manifest.py
===========================
Creates a TSV index of key LAB-Score v1.1 pipeline outputs for reproducibility.
"""

import argparse
from pathlib import Path
import pandas as pd

ap = argparse.ArgumentParser(description="Create LAB-Score release manifest")
ap.add_argument("--pipeline-outdir", required=True, help="Pipeline output directory, e.g. lab_score_out")
ap.add_argument("--outdir", required=True, help="Output directory for manifest")
args = ap.parse_args()

base = Path(args.pipeline_outdir)
outdir = Path(args.outdir)
outdir.mkdir(parents=True, exist_ok=True)

files = [
    ("Master matrix", "08_master_matrix/master_matrix.tsv"),
    ("Main LAB-Score table", "09_scores/LAB_score_v1.tsv"),
    ("Species mean scores", "09_scores/species_mean_scores.tsv"),
    ("Priority class counts", "09_scores/priority_class_counts.tsv"),
    ("Top 100 candidates", "09_scores/top100_candidates.tsv"),
    ("Score distribution summary", "09_scores/LAB_score_distribution_summary.tsv"),
    ("ML model performance", "10_ml/model_performance_summary.tsv"),
    ("ML functional feature importance", "10_ml/RF_classification_functional_feature_importance.tsv"),
    ("Figure outputs", "11_figures"),
    ("Interactive strain report", "12_report/LAB_score_report.html"),
    ("Weight sensitivity summary", "13_sensitivity/weight_sensitivity_summary.tsv"),
    ("Baseline comparison summary", "14_baselines/baseline_comparison_summary.tsv"),
    ("Manuscript HTML summary report", "16_html_report/LAB_SCORE_v1_1_summary_report.html"),
]

rows = []
for desc, rel in files:
    path = base / rel
    rows.append({
        "description": desc,
        "relative_path": rel,
        "path": str(path),
        "exists": path.exists(),
        "size_bytes": path.stat().st_size if path.exists() and path.is_file() else None,
    })

df = pd.DataFrame(rows)
out = outdir / "LAB_SCORE_v1_1_release_manifest.tsv"
df.to_csv(out, sep="\t", index=False)
print("[manifest] Saved:", out)
print(df.to_string(index=False))
