# LAB-Score v1.1 Manuscript Outputs Quick Start

This guide runs the manuscript-ready robustness/reporting modules after you already have `LAB_score_v1.tsv`.

```bash
# From the LAB-Score repository root
python scripts/13_weight_sensitivity_analysis.py \
  --scores lab_score_out/09_scores/LAB_score_v1.tsv \
  --outdir lab_score_out/13_sensitivity

python scripts/14_baseline_comparison.py \
  --scores lab_score_out/09_scores/LAB_score_v1.tsv \
  --outdir lab_score_out/14_baselines

python scripts/15_make_release_manifest.py \
  --pipeline-outdir lab_score_out \
  --outdir lab_score_out/15_release

python scripts/16_make_html_summary_report.py \
  --pipeline-outdir lab_score_out \
  --outdir lab_score_out/16_html_report
```

Open the manuscript summary report:

```bash
xdg-open lab_score_out/16_html_report/LAB_SCORE_v1_1_summary_report.html
```

Key outputs:

| Output | Description |
|---|---|
| `13_sensitivity/weight_sensitivity_summary.tsv` | Robustness to alternative component weights |
| `14_baselines/baseline_comparison_summary.tsv` | Comparison with simplified baseline rankings |
| `15_release/LAB_SCORE_v1_1_release_manifest.tsv` | File manifest for reproducibility |
| `16_html_report/LAB_SCORE_v1_1_summary_report.html` | Tabbed manuscript-style HTML summary report |

Interpretation note: LAB-Score is a genome-based decision-support framework for prioritizing strains for downstream validation. It is not definitive proof of probiotic efficacy or safety.
