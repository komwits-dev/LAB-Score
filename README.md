# LAB-Score Pipeline v3.3

**Important v3.3 change:** the installer now uses a separate `dbcan_clean` conda environment for dbCAN/run_dbcan and a main `lab_score_clean` environment for the rest of the pipeline. This avoids mixed pip/conda `ClobberError` problems seen in older dirty environments.

See `QUICK_START.md` for recommended commands.


End-to-end pipeline: genome FASTA → species metadata → annotation → safety screening → functional panels → LAB-Score v1.1 → ML interpretation → publication figures.



The easiest fix: **delete your current Workflow section** and replace it with this exact block.

Paste this into `README.md`:

````markdown
## Workflow

```mermaid
flowchart TD
    A["Input genomes or existing Prokka results"] --> B["Metadata resolution"]
    A --> C["Genome annotation or imported Prokka"]

    C --> D["Safety screening"]
    C --> E["Functional trait detection"]
    C --> F["CAZyme profiling with dbCAN"]
    C --> G["Genome quality assessment"]

    B --> H["Master feature matrix"]
    D --> H
    E --> H
    F --> H
    G --> H

    H --> I["LAB-Score v1.1 calculation"]
    I --> J["Safety-gated candidate tier assignment"]
    I --> K["Machine-learning interpretation"]
    I --> L["Interactive radar-based report"]

    J --> M["Elite candidate"]
    J --> N["High candidate"]
    J --> O["Moderate candidate"]
    J --> P["Low priority"]
    J --> Q["Cautionary candidate"]
    J --> R["Critical safety review"]

---

## Pipeline overview

```
FASTA genomes
    │
    ▼
Step 00 ── resolve_metadata.py  ─→ species/genus/strain (NCBI API or user TSV or filename)
Step 01 ── run_prokka.sh        ─→ genome annotation (.gff .faa)
Step 02 ── run_amrfinder.sh     ─→ AMR + virulence genes
Step 03 ── run_resfinder.sh     ─→ acquired resistance genes
Step 04 ── run_dbcan.sh         ─→ CAZyme families
Step 05 ── run_checkm.sh        ─→ genome completeness/contamination  [OPTIONAL]
Step 06 ── safety_screen.py     ─→ biogenic amines + hemolysin + AMR merge
Step 07 ── panel_scan.py        ─→ 18 functional probiotic panels
Step 08 ── build_master_matrix  ─→ master_matrix.tsv (all features)
Step 09 ── calculate_lab_score  ─→ LAB_score_v1.tsv + species summary + top100
Step 10 ── ml_interpret.py      ─→ RF regression + classification + SHAP
Step 11 ── make_figures.py      ─→ Fig 2–5 (publication quality, 300 dpi)
```

---

## Quick start

### 1. Install dependencies

```bash
# Create conda environment
conda create -n lab_score -c conda-forge -c bioconda \
    python=3.10 prokka amrfinderplus dbcan \
    pandas numpy scikit-learn matplotlib seaborn scipy biopython \
    -y

conda activate lab_score

# Optional: SHAP for ML interpretation
pip install shap

# Optional: ResFinder
pip install resfinder

# Optional: CheckM
pip install checkm-genome
```

### 2. Set up your genomes

```bash
mkdir genomes/
cp /path/to/your/genomes/*.fna genomes/
```

Supported filename formats:
- NCBI accessions: `GCF_000001405.40.fna`
- Genus_species_strain: `Lactobacillus_acidophilus_NCFM.fna`
- Any format: falls back to filename parsing

### 3. Run the pipeline

```bash
# Basic run (auto-detect species from NCBI for GCF accessions)
bash run_all.sh -i genomes/ -o results/ -t 16 -n

# With your own metadata file
bash run_all.sh -i genomes/ -o results/ -t 16 -m metadata.tsv

# With CheckM genome QC
bash run_all.sh -i genomes/ -o results/ -t 16 -n -q

# Full options
bash run_all.sh -h
```

---

## Options

| Flag | Description | Default |
|------|-------------|---------|
| `-i` | Input genome directory (required) | — |
| `-o` | Output directory | `lab_score_out` |
| `-t` | CPU threads | `8` |
| `-e` | Conda environment name | `lab_score` |
| `-m` | Metadata TSV file (optional) | auto |
| `-q` | Run CheckM QC (optional) | off |
| `-n` | Auto-fetch species from NCBI API | off |

---

## Metadata TSV format (optional)

If you provide `-m metadata.tsv`, it should have these columns (tab-separated):

```
genome          genus           species         strain          assembly_level
GCF_000001.1    Lactobacillus   acidophilus     NCFM            Complete
GCF_000002.1    Limosilactobacillus reuteri     DSM20016        Chromosome
```

Any unrecognised genomes will fall back to NCBI API (if `-n`) or filename parsing.

---

## Outputs

```
results/
├── 00_metadata/
│   └── genome_metadata.tsv          ← species/genus/strain per genome
├── 01_prokka/
│   └── {genome}/                    ← .gff .faa .ffn etc.
├── 02_amrfinder/
│   └── {genome}.tsv                 ← AMR + virulence hits
├── 03_resfinder/
│   └── resfinder_summary.tsv        ← acquired resistance summary
├── 04_dbcan/
│   └── {genome}/overview.txt        ← CAZyme annotation
├── 05_checkm/
│   └── checkm_summary.tsv           ← completeness/contamination
├── 06_safety/
│   └── safety_summary.tsv           ← per-genome safety flags
├── 07_panels/
│   └── all_panels.tsv               ← 18 functional panel counts
├── 08_master_matrix/
│   └── master_matrix.tsv            ← ALL features merged
├── 09_scores/
│   ├── LAB_score_v1.tsv             ★ MAIN OUTPUT
│   ├── species_mean_scores.tsv      ← species-level summary
│   ├── top100_candidates.tsv        ← top 100 non-critical
│   └── LAB_score_distribution_summary.tsv
├── 10_ml/
│   ├── RF_regression_*_feature_importance.tsv
│   ├── RF_classification_*_feature_importance.tsv
│   ├── model_performance_summary.tsv
│   ├── module_spearman_correlation.tsv
│   ├── fig_ml_feature_importance.png
│   ├── fig_ml_actual_vs_predicted.png
│   ├── fig_ml_confusion_matrix.png
│   ├── fig_ml_spearman_heatmap.png
│   └── fig_ml_shap_summary.png      ← if shap installed
└── 11_figures/
    ├── fig2_score_distribution.png  ★
    ├── fig3_candidate_tiers.png     ★
    ├── fig4_species_summary.png     ★
    └── fig5_ml_feature_importance.png ★
```

---

## LAB-Score v1.1 formula

```
LAB_score_v1 = 0.45 × Safety_score
             + 0.25 × GI_survival_score
             + 0.20 × Functional_score
             + 0.10 × Fermentation_score
```

### Component definitions

| Component | Source panels | Weight |
|-----------|--------------|--------|
| Safety_score | AMRFinder + ResFinder + VFDB + hemolysin + biogenic amines | 0.45 |
| GI_survival_score | acid_energy, bile_resistance, gut_persistence, osmotic_stress, heat_stress, cold_stress, GABA | 0.25 |
| Functional_score | vitamins, immunomodulation, bacteriocins, EPS, adhesion, CRISPR, QS, alkaline | 0.20 |
| Fermentation_score | carbohydrate, metabolism, CAZymes | 0.10 |

### Priority classes

| LAB-Score | Class |
|-----------|-------|
| ≥ 85 | Elite |
| 70–84.9 | High |
| 50–69.9 | Moderate |
| < 50 | Low |

### Refined safety gate (v1.1)

| Safety marker | Status prefix |
|---------------|---------------|
| AMR gene / virulence gene / ResFinder hit | Critical safety review |
| Hemolysin / biogenic amine marker | Cautionary |
| None | (no prefix) |

---

## Functional panels (18 total)

| Panel | Key genes/functions |
|-------|-------------------|
| acid_energy | ATP synthase (atpA-I), gadA/B/C |
| adhesion_biofilm | biofilm, mub, MapA, mucin-binding |
| adhesion_surface | S-layer, sortase, pilus, LPXTG |
| alkalinestress | Na+/H+ antiporter, nhaA/B/C |
| antipath_qs | luxS, quorum sensing, agrABCD |
| bacteriocins | nisin, lantibiotics, pediocin, sakacin |
| bileresistance | BSH, bile salt hydrolase |
| carbohydrate | PTS, lacZ, amylase, glucosidase |
| cellenvelope_eps | EPS biosynthesis, teichoic acid |
| coldstress | cspA-E, cold shock proteins |
| defense_crispr | cas1-9, CRISPR, R/M systems |
| gaba | gadA/B, GABA transaminase |
| gutpersistence | mucin binding, NADH oxidase |
| heatstress | groEL/ES, dnaK/J, clpB/P |
| immunomodulation | LTA, flagellin, peptidoglycan hydrolase |
| metabolism | ldh, ackA, pta, pyruvate kinase |
| osmoticstress | opuA/B/C, glycine betaine, ectoine |
| vitamins | riboflavin, folate, cobalamin, thiamine |

---

## Citation

> [Your manuscript] — *in preparation*


## v2.4 update: using existing Prokka results

This version can start from either raw genome FASTA files or an existing Prokka output directory.

### Run from genome FASTA files

```bash
bash run_all.sh -i genomes/ -o results/ -t 32 -n
```

### Run from existing Prokka output only

Use this when you already have folders such as `results/01_prokka/sample/*.gff`, `*.faa`, and `*.fna`.
The pipeline will skip Prokka, import the existing Prokka folders, and use the Prokka `.fna` files for metadata and ResFinder.

```bash
bash run_all.sh -p previous_results/01_prokka -o results_from_prokka -t 32 -n
```

### Recommended mode: existing Prokka + original genome FASTA

This is the safest mode because ResFinder and CheckM use the original genome FASTA files, while safety and functional scans reuse the existing Prokka annotations.

```bash
bash run_all.sh \
  -i genomes/ \
  -p previous_results/01_prokka \
  -o results_from_prokka \
  -t 32 \
  -n
```

The `-p` option expects one subdirectory per genome, for example:

```text
01_prokka/
├── GCA_000239955.2_ASM23995v2_genomic/
│   ├── GCA_000239955.2.gff
│   ├── GCA_000239955.2.faa
│   └── GCA_000239955.2.fna
└── GCA_001436895.1_ASM143689v1_genomic/
    ├── GCA_001436895.1.gff
    ├── GCA_001436895.1.faa
    └── GCA_001436895.1.fna
```

Downstream sample IDs are taken from the Prokka **folder name**, not from the internal Prokka prefix, so outputs remain aligned across AMRFinder, dbCAN, safety screening, panel scanning, and the master matrix.
"# LAB-Score" 
