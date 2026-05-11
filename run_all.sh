#!/usr/bin/env bash
# =============================================================================
# LAB-Score Pipeline v3.3
# run_all.sh — Full pipeline from FASTA and/or existing/nested Prokka results
#              to LAB-Score + ML interpretation + publication-ready figures.
#
# By default, this version checks and installs missing prerequisite tools before
# running. Use --skip-install when the environment is already prepared.
#
# Usage:
#   bash run_all.sh [OPTIONS]
#
# Required input: provide at least one of:
#   -i, --input DIR              Input genome FASTA directory
#   -p, --prokka DIR             Existing Prokka result directory; skips Prokka
#
# Main options:
#   -o, --outdir DIR             Output directory [default: lab_score_out]
#   -t, --threads INT            Threads [default: 8]
#   -e, --env NAME               Main conda env name [default: lab_score_clean]
#   --dbcan-env-name NAME       Separate dbCAN conda env [default: dbcan_clean]
#   -m, --metadata FILE          Metadata TSV with genome/accession/species fields
#   -q, --checkm                 Run CheckM QC [optional]
#   -n, --ncbi                   NCBI auto-fetch species metadata [kept for compatibility]
  --ref-dir DIR                Directory of LAB reference genomes for FastANI
  --download-refs              Auto-download LAB reference genomes from NCBI
#
# Installation/preflight options:
#   --skip-install               Do not run installer before pipeline
#   --check-only                 Check prerequisites only, then exit
#   --dbcan-db-dir DIR           dbCAN database directory
#   --skip-dbcan-db              Do not download dbCAN database during install
#   --skip-optional              Skip optional ResFinder/CheckM/SHAP installation
#   --force-dbcan-reinstall      Alias for --recreate-dbcan-env during install
#   --recreate-env              Recreate main conda env during install
#   --recreate-dbcan-env        Recreate separate dbCAN env during install
#
# Examples:
#   bash run_all.sh -i genomes/ -o results/ -t 32 -n
#   bash run_all.sh -p annotations_prokka_GCA -o results_from_prokka -t 32 -n
#   bash run_all.sh -p annotations_prokka_GCA -o results -t 32 --dbcan-db-dir /media/mecob/komwit/db/dbcan
#   bash run_all.sh -p annotations_prokka_GCA -o results -t 32 --skip-install
# =============================================================================

set -euo pipefail

INDIR=""
PROKKA_IN=""
OUTDIR="lab_score_out"
THREADS=8
CONDA_ENV="lab_score_clean"
DBCAN_CONDA_ENV="dbcan_clean"
METADATA=""
RUN_CHECKM=false
NCBI_FETCH=false
REF_DIR=""
DOWNLOAD_REFS=false
AUTO_INSTALL=true
CHECK_ONLY=false
DBCAN_DB_DIR_ARG=""
SKIP_DBCAN_DB=false
SKIP_OPTIONAL=false
FORCE_DBCAN_REINSTALL=false
RECREATE_ENV=false
RECREATE_DBCAN_ENV=false
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
  grep '^#' "$0" | grep -v '^#!/' | sed 's/^# \?//'
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input) INDIR="${2:?Missing input directory}"; shift 2 ;;
    -p|--prokka) PROKKA_IN="${2:?Missing Prokka directory}"; shift 2 ;;
    -o|--outdir) OUTDIR="${2:?Missing output directory}"; shift 2 ;;
    -t|--threads) THREADS="${2:?Missing thread count}"; shift 2 ;;
    -e|--env) CONDA_ENV="${2:?Missing conda env name}"; shift 2 ;;
    --dbcan-env-name) DBCAN_CONDA_ENV="${2:?Missing dbCAN conda env name}"; shift 2 ;;
    -m|--metadata) METADATA="${2:?Missing metadata file}"; shift 2 ;;
    -q|--checkm) RUN_CHECKM=true; shift ;;
    -n|--ncbi) NCBI_FETCH=true; shift ;;
    --ref-dir) REF_DIR="${2:?Missing ref dir}"; shift 2 ;;
    --download-refs) DOWNLOAD_REFS=true; shift ;;
    --skip-install) AUTO_INSTALL=false; shift ;;
    --check-only) CHECK_ONLY=true; shift ;;
    --dbcan-db-dir) DBCAN_DB_DIR_ARG="${2:?Missing dbCAN database directory}"; shift 2 ;;
    --skip-dbcan-db) SKIP_DBCAN_DB=true; shift ;;
    --skip-optional) SKIP_OPTIONAL=true; shift ;;
    --force-dbcan-reinstall) FORCE_DBCAN_REINSTALL=true; RECREATE_DBCAN_ENV=true; shift ;;
    --recreate-env) RECREATE_ENV=true; shift ;;
    --recreate-dbcan-env) RECREATE_DBCAN_ENV=true; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[ERROR] Unknown option: $1"; usage; exit 1 ;;
  esac
done

if [[ -z "$INDIR" && -z "$PROKKA_IN" && "$CHECK_ONLY" != true ]]; then
  echo "[ERROR] Provide either genome FASTA input (-i) or existing Prokka output (-p)."
  exit 1
fi
if [[ -n "$INDIR" && ! -d "$INDIR" ]]; then
  echo "[ERROR] Input genome directory not found: $INDIR"
  exit 1
fi
if [[ -n "$PROKKA_IN" && ! -d "$PROKKA_IN" ]]; then
  echo "[ERROR] Existing Prokka directory not found: $PROKKA_IN"
  exit 1
fi
if [[ -n "$METADATA" && ! -f "$METADATA" ]]; then
  echo "[ERROR] Metadata file not found: $METADATA"
  exit 1
fi

# -----------------------------------------------------------------------------
# Step -1: install/check prerequisites before running
# -----------------------------------------------------------------------------
INSTALL_ARGS=(--env-name "$CONDA_ENV" --dbcan-env-name "$DBCAN_CONDA_ENV")
if [[ -n "$DBCAN_DB_DIR_ARG" ]]; then
  INSTALL_ARGS+=(--dbcan-db-dir "$DBCAN_DB_DIR_ARG")
  export DBCAN_DB_DIR="$DBCAN_DB_DIR_ARG"
fi
if [[ "$SKIP_DBCAN_DB" == true ]]; then
  INSTALL_ARGS+=(--skip-dbcan-db)
fi
if [[ "$SKIP_OPTIONAL" == true ]]; then
  INSTALL_ARGS+=(--skip-optional)
fi
if [[ "$FORCE_DBCAN_REINSTALL" == true ]]; then
  INSTALL_ARGS+=(--force-dbcan-reinstall)
fi
if [[ "$RECREATE_ENV" == true ]]; then
  INSTALL_ARGS+=(--recreate-env)
fi
if [[ "$RECREATE_DBCAN_ENV" == true ]]; then
  INSTALL_ARGS+=(--recreate-dbcan-env)
fi

if [[ "$CHECK_ONLY" == true ]]; then
  bash "$SCRIPT_DIR/install.sh" --check "${INSTALL_ARGS[@]}"
  exit $?
fi

if [[ "$AUTO_INSTALL" == true ]]; then
  echo "======================================================================"
  echo "  LAB-Score Pipeline v3.3 prerequisite setup"
  echo "======================================================================"
  bash "$SCRIPT_DIR/install.sh" "${INSTALL_ARGS[@]}"
else
  echo "[Preflight] Auto-install skipped by user. Checking required prerequisites..."
  bash "$SCRIPT_DIR/install.sh" --check "${INSTALL_ARGS[@]}"
fi

# Activate conda environment after install/check.
activate_env() {
  local env=$1
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true
  conda activate "$env"
}
activate_env "$CONDA_ENV"
# Prevent base/old Miniconda Perl variables from breaking Prokka.
unset PERL5LIB PERL_LOCAL_LIB_ROOT PERL_MB_OPT PERL_MM_OPT
export CONDA_ENV
export DBCAN_CONDA_ENV

if [[ -z "${DBCAN_DB_DIR:-}" ]]; then
  DBCAN_ENV_PREFIX="$(conda env list | awk -v env="$DBCAN_CONDA_ENV" '$1==env {print $NF; exit}')"
  export DBCAN_DB_DIR="${DBCAN_ENV_PREFIX}/db"
fi

mkdir -p \
  "$OUTDIR/00_inputs/genomes" \
  "$OUTDIR/00_metadata" \
  "$OUTDIR/01_prokka" \
  "$OUTDIR/02_amrfinder" \
  "$OUTDIR/03_resfinder" \
  "$OUTDIR/04_dbcan" \
  "$OUTDIR/05_checkm" \
  "$OUTDIR/06_safety" \
  "$OUTDIR/07_panels" \
  "$OUTDIR/08_master_matrix" \
  "$OUTDIR/09_scores" \
  "$OUTDIR/10_ml" \
  "$OUTDIR/11_figures" \
  "$OUTDIR/12_report" \
  "$OUTDIR/logs"

LOG="$OUTDIR/logs/pipeline_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

echo "======================================================================"
echo "  LAB-Score Pipeline v3.3"
echo "  $(date)"
echo "======================================================================"
echo "  Genome input:     ${INDIR:-none}"
echo "  Existing Prokka:  ${PROKKA_IN:-none}"
echo "  Output:           $OUTDIR"
echo "  Threads:          $THREADS"
echo "  Main env:         $CONDA_ENV"
echo "  dbCAN env:        $DBCAN_CONDA_ENV"
echo "  Metadata:         ${METADATA:-auto}"
echo "  CheckM:           $RUN_CHECKM"
echo "  NCBI:             $NCBI_FETCH"
echo "  dbCAN DB:         ${DBCAN_DB_DIR}"
echo "  Auto-install:     $AUTO_INSTALL"
echo "======================================================================"

WORK_GENOMES="$INDIR"

if [[ -n "$PROKKA_IN" ]]; then
  echo ""
  echo "[Prepare] Importing existing Prokka outputs..."
  bash "$SCRIPT_DIR/scripts/00_prepare_prokka_input.sh" \
    "$PROKKA_IN" "$OUTDIR/01_prokka" "$OUTDIR/00_inputs/genomes"

  if [[ -z "$INDIR" ]]; then
    WORK_GENOMES="$OUTDIR/00_inputs/genomes"
  fi
  echo "[Prepare] Done."
fi

if [[ -z "$WORK_GENOMES" || ! -d "$WORK_GENOMES" ]]; then
  echo "[ERROR] No usable genome FASTA directory available. Provide -i or ensure Prokka .fna files exist in -p."
  exit 1
fi

NGENOMES=$(find "$WORK_GENOMES" -maxdepth 1 \( -name "*.fna" -o -name "*.fa" -o -name "*.fasta" \) 2>/dev/null | wc -l)
if [[ "$NGENOMES" -eq 0 ]]; then
  echo "[ERROR] No genome FASTA files (.fna/.fa/.fasta) found in working genome directory: $WORK_GENOMES"
  echo "        Use -i /path/to/genomes, or provide Prokka results containing *.fna files."
  exit 1
fi

echo "[Prepare] Working genome FASTA directory: $WORK_GENOMES ($NGENOMES genomes)"

echo ""
echo "[Step 0] Resolving species metadata..."
python "$SCRIPT_DIR/scripts/00_resolve_metadata.py" \
  --indir    "$WORK_GENOMES" \
  --outdir   "$OUTDIR/00_metadata" \
  --metadata "${METADATA:-}" \
  --ncbi
echo "[Step 0] Done."

echo ""
echo "[Step 0b] Verifying species identity with FastANI..."
python "$SCRIPT_DIR/scripts/00b_species_verify.py" \
  --indir    "$WORK_GENOMES" \
  --metadata "$OUTDIR/00_metadata/genome_metadata.tsv" \
  --outdir   "$OUTDIR/00_metadata" \
  --threads  "$THREADS" \
  ${REF_DIR:+--ref_dir "$REF_DIR"} \
  $([[ "$DOWNLOAD_REFS" == true ]] && echo "--download_refs" || true)
echo "[Step 0b] Done."

if [[ -n "$PROKKA_IN" ]]; then
  echo ""
  echo "[Step 1] Prokka skipped — using existing Prokka results from: $PROKKA_IN"
  PROKKA_COUNT=$(find "$OUTDIR/01_prokka" -mindepth 2 -maxdepth 2 -name "*.gff" | wc -l)
  echo "[Step 1] Imported Prokka GFF files: $PROKKA_COUNT"
else
  echo ""
  echo "[Step 1] Running Prokka annotation..."
  bash "$SCRIPT_DIR/scripts/01_run_prokka.sh" \
    "$WORK_GENOMES" "$OUTDIR/01_prokka" "$THREADS"
  echo "[Step 1] Done."
fi

echo ""
echo "[Step 2] Running AMRFinder..."
bash "$SCRIPT_DIR/scripts/02_run_amrfinder.sh" \
  "$OUTDIR/01_prokka" "$OUTDIR/02_amrfinder" "$THREADS"
echo "[Step 2] Done."

echo ""
echo "[Step 3] Running ResFinder..."
bash "$SCRIPT_DIR/scripts/03_run_resfinder.sh" \
  "$WORK_GENOMES" "$OUTDIR/03_resfinder" "$THREADS"
echo "[Step 3] Done."

echo ""
echo "[Step 4] Running dbCAN..."
bash "$SCRIPT_DIR/scripts/04_run_dbcan.sh" \
  "$OUTDIR/01_prokka" "$OUTDIR/04_dbcan" "$THREADS"
echo "[Step 4] Done."

if [[ "$RUN_CHECKM" == true ]]; then
  echo ""
  echo "[Step 5] Running CheckM..."
  bash "$SCRIPT_DIR/scripts/05_run_checkm.sh" \
    "$WORK_GENOMES" "$OUTDIR/05_checkm" "$THREADS"
  echo "[Step 5] Done."
else
  echo ""
  echo "[Step 5] CheckM skipped (use -q to enable)."
  echo -e "genome\tCompleteness\tContamination" > "$OUTDIR/05_checkm/checkm_summary.tsv"
fi

echo ""
echo "[Step 6] Running safety screening..."
python "$SCRIPT_DIR/scripts/06_safety_screen.py" \
  --prokka_dir   "$OUTDIR/01_prokka" \
  --amrfinder    "$OUTDIR/02_amrfinder" \
  --resfinder    "$OUTDIR/03_resfinder" \
  --outdir       "$OUTDIR/06_safety"
echo "[Step 6] Done."

echo ""
echo "[Step 7] Scanning functional panels..."
python "$SCRIPT_DIR/scripts/07_panel_scan.py" \
  --prokka_dir "$OUTDIR/01_prokka" \
  --outdir     "$OUTDIR/07_panels"
echo "[Step 7] Done."

echo ""
echo "[Step 8] Building master matrix..."
python "$SCRIPT_DIR/scripts/08_build_master_matrix.py" \
  --metadata    "$OUTDIR/00_metadata/genome_metadata.tsv" \
  --safety      "$OUTDIR/06_safety/safety_summary.tsv" \
  --panels      "$OUTDIR/07_panels/all_panels.tsv" \
  --dbcan       "$OUTDIR/04_dbcan" \
  --checkm      "$OUTDIR/05_checkm/checkm_summary.tsv" \
  --outdir      "$OUTDIR/08_master_matrix"
echo "[Step 8] Done."

echo ""
echo "[Step 9] Calculating LAB-Score v1.1..."
python "$SCRIPT_DIR/scripts/09_calculate_lab_score.py" \
  --matrix  "$OUTDIR/08_master_matrix/master_matrix.tsv" \
  --outdir  "$OUTDIR/09_scores"
echo "[Step 9] Done."

echo ""
echo "[Step 10] Running ML interpretation..."
python "$SCRIPT_DIR/scripts/10_ml_interpret.py" \
  --scores  "$OUTDIR/09_scores/LAB_score_v1.tsv" \
  --outdir  "$OUTDIR/10_ml"
echo "[Step 10] Done."

echo ""
echo "[Step 11] Generating figures..."
python "$SCRIPT_DIR/scripts/11_make_figures.py" \
  --scores    "$OUTDIR/09_scores/LAB_score_v1.tsv" \
  --species   "$OUTDIR/09_scores/species_mean_scores.tsv" \
  --ml        "$OUTDIR/10_ml" \
  --outdir    "$OUTDIR/11_figures"
echo "[Step 11] Done."
echo ""
echo "[Step 12] Generating interactive HTML report..."
python "$SCRIPT_DIR/scripts/12_interactive_report.py" \
  --scores   "$OUTDIR/09_scores/LAB_score_v1.tsv" \
  --matrix   "$OUTDIR/08_master_matrix/master_matrix.tsv" \
  --metadata "$OUTDIR/00_metadata/genome_metadata.tsv" \
  --outdir   "$OUTDIR/12_report"
echo "[Step 12] Done."


echo ""
echo "======================================================================"
echo "  LAB-Score Pipeline COMPLETE"
echo "  $(date)"
echo "======================================================================"
echo ""
echo "  Key outputs:"
echo "    Master matrix : $OUTDIR/08_master_matrix/master_matrix.tsv"
echo "    LAB Scores    : $OUTDIR/09_scores/LAB_score_v1.tsv"
echo "    Top 100       : $OUTDIR/09_scores/top100_candidates.tsv"
echo "    Figures       : $OUTDIR/11_figures/"
echo "    ML results    : $OUTDIR/10_ml/"
echo "    dbCAN status  : $OUTDIR/04_dbcan/dbcan_status.tsv"
echo "    Full log      : $LOG
    Interactive report: $OUTDIR/12_report/LAB_score_report.html
    Species verify: $OUTDIR/00_metadata/genome_metadata.tsv (ani_species, species_status columns)"
echo ""
python - "$OUTDIR/09_scores/LAB_score_v1.tsv" <<'PY'
import pandas as pd, sys
try:
    df = pd.read_csv(sys.argv[1], sep="\t")
    print("  Priority class summary:")
    col = "Priority_class" if "Priority_class" in df.columns else "Candidate_tier_v1_1"
    for cls, n in df[col].value_counts().items():
        print(f"    {cls:<32} : {n:>5} genomes")
    print(f"\n  Total genomes scored: {len(df)}")
except Exception as e:
    print(f"  Summary unavailable: {e}")
PY
