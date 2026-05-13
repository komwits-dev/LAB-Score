# LAB-Score

<p align="center">
  <img src="docs/images/lab_score_workflow.svg" alt="LAB-Score workflow" width="100%">
</p>

<p align="center">
  <b>A safety-gated genomic prioritization framework for lactic acid bacteria candidates</b>
</p>

<p align="center">
  <img alt="Version" src="https://img.shields.io/badge/version-v3.6-2D6A4F">
  <img alt="Score" src="https://img.shields.io/badge/scoring-LAB--Score%20v1.1-0096C7">
  <img alt="Python" src="https://img.shields.io/badge/python-3.10-blue">
  <img alt="License" src="https://img.shields.io/badge/license-MIT-green">
</p>

---

## Overview

**LAB-Score** is a genome-based framework for prioritizing lactic acid bacteria (LAB) candidates for probiotic, fermentation, and functional-food applications.

The pipeline integrates:

- genome annotation or existing Prokka annotation import;
- metadata and species resolution;
- safety screening;
- antimicrobial resistance and virulence marker detection;
- functional trait panel scanning;
- dbCAN-based CAZyme profiling;
- genome quality integration;
- LAB-Score v1.1 calculation;
- machine-learning interpretation;
- interactive radar-based strain comparison and exportable reports;
- weight sensitivity analysis and baseline comparison;
- manuscript-style tabbed HTML summary report.

> **Pipeline version:** v3.6  
> **Scoring model:** LAB-Score v1.1  
> **Recommended use:** candidate prioritization before experimental validation.

---

## Key features

- Accepts **raw genome assemblies** or **existing Prokka annotation folders**
- Supports nested Prokka layouts such as `species/accession/*.gff`
- Performs safety-gated LAB candidate classification
- Screens functional modules relevant to LAB applications
- Integrates dbCAN CAZyme annotation
- Produces genome-level and species-level score summaries
- Generates publication-ready figures and tables
- Includes an interactive HTML report
- Supports radar-based strain comparison
- Exports selected strain information as TSV, JSON, PNG, SVG, and HTML reports
- Includes weight sensitivity analysis and baseline ranking comparison
- Generates a manuscript-style tabbed HTML summary report

---

## Interactive report preview

LAB-Score v3.6 generates an interactive HTML report for exploring large-scale LAB genome prioritization results.

### Overview dashboard

The overview dashboard summarizes the number of analyzed genomes, candidate-tier distribution, LAB-Score distribution, safety-gated classes, and top species.

![LAB-Score overview dashboard](docs/images/overview_dashboard.png)

### Strain-level report and export panel

Users can click any strain to view a detailed report containing LAB-Score breakdown, safety analysis, species identification, genome quality, CAZyme profile, and functional trait information. Each strain report can be exported as TSV, JSON, radar figure, or HTML report.

![LAB-Score strain detail modal](docs/images/strain_detail_modal.png)

### Radar-based strain comparison

The Radar Compare tab allows users to search, filter, and select strains for side-by-side comparison across functional modules. Selected strains can be exported as figures, comparison tables, or strain-level reports.

![LAB-Score radar comparison](docs/images/radar_compare.png)

### Gene presence / abundance heatmap

The Gene Presence tab provides an interactive heatmap for comparing functional panel abundance across selected genomes or species.

![LAB-Score gene presence heatmap](docs/images/gene_presence_heatmap.png)

---

## LAB-Score v1.1 refined safety-gated logic

LAB-Score uses a **refined safety-gated strategy**. Genomes with favorable functional profiles are prioritized only after safety screening, while genomes with cautionary or critical markers are retained but clearly separated for downstream interpretation.

### Refined safety status

| Safety status | Definition |
|---|---|
| No safety marker detected / Pass | No AMR, virulence, hemolysin-associated, or biogenic amine-associated marker detected |
| Cautionary safety review / Caution | Putative hemolysin-associated or biogenic amine-associated marker detected |
| Critical safety review / Critical | AMR gene or virulence-associated marker detected |

### Candidate-tier thresholds

| LAB-Score range | Tier |
|---|---|
| ≥85 | Elite candidate |
| 70–84.9 | High candidate |
| 50–69.9 | Moderate candidate |
| <50 | Low priority |

If a genome has cautionary markers, the score-based tier receives a cautionary prefix, for example `Cautionary High candidate`. Genomes with AMR or virulence-associated markers are assigned to `Critical safety review`, regardless of the score.

> LAB-Score is a genome-based decision-support framework for candidate prioritization. It does not replace phenotypic safety testing or experimental probiotic/fermentation validation.

---

## Manuscript v1.1 release outputs

The current release adds manuscript-scale robustness and reporting modules for LAB-Score v1.1.

| Output | Description |
|---|---|
| `09_scores/LAB_score_v1.tsv` | Genome-level LAB-Score results |
| `09_scores/top100_candidates.tsv` | Top 100 non-critical candidate genomes |
| `09_scores/species_mean_scores.tsv` | Species-level LAB-Score summary |
| `10_ml/model_performance_summary.tsv` | Random Forest model performance |
| `10_ml/RF_classification_functional_feature_importance.tsv` | Functional-module feature importance |
| `13_sensitivity/weight_sensitivity_summary.tsv` | Robustness across alternative scoring weights |
| `14_baselines/baseline_comparison_summary.tsv` | Comparison with simplified baseline rankings |
| `15_release/LAB_SCORE_v1_1_release_manifest.tsv` | File manifest for reproducibility |
| `16_html_report/LAB_SCORE_v1_1_summary_report.html` | Tabbed manuscript-style HTML summary report |

### Manuscript figures

The `docs/manuscript/figures/` folder contains manuscript-ready figure drafts:

| Figure | Description |
|---|---|
| Figure 1 | LAB-Score workflow |
| Figure 2 | LAB-Score distribution with thresholds |
| Figure 3 | Refined safety-gated candidate tier distribution |
| Figure 4 | Composition of LAB-Score v1.1 candidate tiers |
| Figure 5 | Functional modules driving LAB-Score candidate prioritization |

### Robustness analyses

The v3.6 release adds two robustness analyses:

1. **Weight sensitivity analysis** — compares the main safety-heavy LAB-Score v1.1 formula with equal, function-heavy, GI-survival-heavy, and fermentation-heavy weighting schemes.
2. **Baseline comparison** — compares integrated LAB-Score rankings against safety-only, GI-only, functional-only, fermentation-only, function-plus-fermentation, and no-safety baseline rankings.

---

## Repository structure

```text
LAB-Score/
├── README.md
├── QUICK_START.md
├── QUICK_START_RADAR.md
├── QUICK_START_MANUSCRIPT.md
├── VERSION.txt
├── LICENSE
├── .gitignore
├── install.sh
├── run_all.sh
├── docs/
│   └── images/
│       ├── lab_score_workflow.svg
│       ├── overview_dashboard.png
│       ├── strain_detail_modal.png
│       ├── radar_compare.png
│       └── gene_presence_heatmap.png
└── scripts/
    ├── 00_prepare_prokka_input.sh
    ├── 00_resolve_metadata.py
    ├── 00b_species_verify.py
    ├── 01_run_prokka.sh
    ├── 02_run_amrfinder.sh
    ├── 03_run_resfinder.sh
    ├── 04_run_dbcan.sh
    ├── 05_run_checkm.sh
    ├── 06_safety_screen.py
    ├── 07_panel_scan.py
    ├── 08_build_master_matrix.py
    ├── 09_calculate_lab_score.py
    ├── 10_ml_interpret.py
    ├── 11_make_figures.py
    ├── 12_interactive_report.py
    ├── 13_weight_sensitivity_analysis.py
    ├── 14_baseline_comparison.py
    ├── 15_make_release_manifest.py
    └── 16_make_html_summary_report.py
```

---

## What does `install.sh` do?

The `install.sh` script prepares the software environment required to run LAB-Score. It was designed to avoid common dependency conflicts in bioinformatics workflows by separating the main LAB-Score environment from the dbCAN environment.

### Environment design

LAB-Score uses two conda environments:

| Environment | Purpose |
|---|---|
| `lab_score_clean` | Main LAB-Score workflow, Prokka, AMRFinderPlus, scoring, figures, and interactive report |
| `dbcan_clean` | Dedicated environment for run_dbCAN, DIAMOND, HMMER, Prodigal, and CAZyme annotation |

This separation is important because run_dbCAN can create Python dependency conflicts when installed together with other bioinformatics tools.

### Main functions of `install.sh`

`install.sh` performs the following steps:

1. Checks whether conda is available.
2. Creates the main LAB-Score conda environment.
3. Creates a separate dbCAN conda environment.
4. Installs core Python packages such as `pandas`, `numpy`, `scikit-learn`, and `matplotlib`.
5. Installs Prokka and required Perl dependencies.
6. Installs AMRFinderPlus and updates the AMRFinderPlus database.
7. Installs FastANI for optional species verification.
8. Installs run_dbCAN, DIAMOND, HMMER, Prodigal, and supporting dbCAN dependencies.
9. Downloads or checks the dbCAN database.
10. Performs a final prerequisite check.

### Recommended installation command

Replace `/path/to/dbcan_database` with the directory where you want to store or already have the dbCAN database.

```bash
bash install.sh \
  --env-name lab_score_clean \
  --dbcan-env-name dbcan_clean \
  --dbcan-db-dir /path/to/dbcan_database \
  --skip-optional
```

Example using a project-local database directory:

```bash
bash install.sh \
  --env-name lab_score_clean \
  --dbcan-env-name dbcan_clean \
  --dbcan-db-dir ./databases/dbcan \
  --skip-optional
```

### Check installation only

```bash
bash run_all.sh \
  --check-only \
  --env lab_score_clean \
  --dbcan-env-name dbcan_clean \
  --dbcan-db-dir /path/to/dbcan_database
```

For the project-local example above:

```bash
bash run_all.sh \
  --check-only \
  --env lab_score_clean \
  --dbcan-env-name dbcan_clean \
  --dbcan-db-dir ./databases/dbcan
```

### Useful installation options

| Option | Meaning |
|---|---|
| `--env-name` | Name of the main LAB-Score conda environment |
| `--dbcan-env-name` | Name of the separate dbCAN conda environment |
| `--dbcan-db-dir` | Directory for the dbCAN database |
| `--skip-optional` | Skip optional tools such as CheckM, ResFinder, or SHAP |
| `--skip-dbcan-db` | Do not download the dbCAN database |
| `--recreate-env` | Delete and recreate the main LAB-Score environment |
| `--recreate-dbcan-env` | Delete and recreate the dbCAN environment |
| `--check` | Check prerequisites only without installing |

### When should I use `--recreate-dbcan-env`?

Use this option if run_dbCAN is broken or if the system accidentally uses an old run_dbCAN from another conda environment.

```bash
bash install.sh \
  --env-name lab_score_clean \
  --dbcan-env-name dbcan_clean \
  --dbcan-db-dir /path/to/dbcan_database \
  --recreate-dbcan-env \
  --skip-optional
```

---

## Quick start

### Option 1: Run from genome assemblies

```bash
bash run_all.sh \
  -i genomes \
  -o results \
  -t 32 \
  --env lab_score_clean \
  --dbcan-env-name dbcan_clean \
  --dbcan-db-dir ./databases/dbcan
```

### Option 2: Run from existing Prokka results

```bash
bash run_all.sh \
  -p annotations_prokka_GCA \
  -o results_from_existing_prokka \
  -t 32 \
  --env lab_score_clean \
  --dbcan-env-name dbcan_clean \
  --dbcan-db-dir ./databases/dbcan \
  --skip-install \
  -n
```

The `-p` option accepts a Prokka result directory. Nested layouts are supported, for example:

```text
annotations_prokka_GCA/
├── Lactiplantibacillus_plantarum/
│   ├── GCA_000000001.1/
│   │   ├── GCA_000000001.1.gff
│   │   ├── GCA_000000001.1.faa
│   │   └── GCA_000000001.1.fna
```

---

## Main outputs

| Output | Description |
|---|---|
| `09_scores/LAB_score_v1.tsv` | Genome-level LAB-Score results and candidate tiers |
| `08_master_matrix/master_matrix.tsv` | Integrated feature matrix used for scoring |
| `09_scores/species_mean_scores.tsv` | Species-level LAB-Score summary when multiple species are available |
| `11_figures/` | Publication-ready summary figures |
| `12_report/LAB_score_report.html` | Interactive strain-level report with radar comparison |
| `08_master_matrix/dbcan_cazy_summary.tsv` | CAZyme family summary from dbCAN outputs |
| `10_ml/feature_importance.tsv` | Machine-learning feature interpretation |

---

## Interactive radar report

Generate only the interactive report:

```bash
python scripts/12_interactive_report.py \
  --scores results_from_existing_prokka/09_scores/LAB_score_v1.tsv \
  --matrix results_from_existing_prokka/08_master_matrix/master_matrix.tsv \
  --metadata results_from_existing_prokka/00_metadata/genome_metadata.tsv \
  --outdir results_from_existing_prokka/12_report
```

Open:

```text
results_from_existing_prokka/12_report/LAB_score_report.html
```

---

## Important interpretation note

LAB-Score is a **genome-based prioritization framework**, not a replacement for phenotypic validation.

Candidate genomes assigned to elite or high tiers represent strains with favorable predicted profiles after safety gating. Their functional potential should be confirmed using targeted phenotypic assays, such as antimicrobial susceptibility testing, acid and bile tolerance, adhesion assays, bacteriocin activity, carbohydrate utilization, exopolysaccharide production, GABA production, or other application-specific experiments.

---


## Citation

If you use LAB-Score, please cite this repository and the associated manuscript when available.

```text
LAB-Score: a safety-gated genomic framework for prioritizing lactic acid bacteria candidates with interactive strain-level interpretation.
```

---

## License

This project is released under the MIT License.

---

## Running robustness/report modules only

If you already have `LAB_score_v1.tsv`, you can run only the new robustness and report modules:

```bash
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

Open the summary report:

```bash
xdg-open lab_score_out/16_html_report/LAB_SCORE_v1_1_summary_report.html
```
