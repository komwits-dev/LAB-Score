#!/usr/bin/env bash
# 04_run_dbcan.sh — robust dbCAN runner for LAB-Score v3.3
# Usage: bash 04_run_dbcan.sh <prokka_dir> <outdir> <threads>
#
# Accepts either:
#   1) flattened Prokka folders:     prokka_dir/sample/*.faa
#   2) nested Prokka folders:        prokka_dir/species/sample/*.faa
#   3) any recursive directory tree containing Prokka .faa files
#
# Notes:
#   - Uses dbCAN CAZyme_annotation in protein mode on Prokka .faa files.
#   - Does NOT silently treat dbCAN failures as biological zeroes. Failures are
#     written to dbcan_status.tsv so empty results can be diagnosed.
#   - Optional environment variable: DBCAN_DB_DIR=/path/to/dbCAN/db

set -uo pipefail

PROKKA_DIR=${1:?"Prokka directory required"}
OUTDIR=${2:?"output directory required"}
THREADS=${3:-8}
CONDA_ENV=${CONDA_ENV:-lab_score_clean}
DBCAN_CONDA_ENV=${DBCAN_CONDA_ENV:-dbcan_clean}

JOBS=$(( THREADS > 4 ? 4 : THREADS ))
[[ $JOBS -lt 1 ]] && JOBS=1
CPUS_PER=$(( THREADS / JOBS ))
[[ $CPUS_PER -lt 1 ]] && CPUS_PER=1

mkdir -p "$OUTDIR" logs/dbcan

# -----------------------------------------------------------------------------
# Locate dbCAN database directory
# Priority:
#   1. DBCAN_DB_DIR environment variable
#   2. current active conda prefix /db
#   3. lab_score conda env /db
# -----------------------------------------------------------------------------
DB_DIR="${DBCAN_DB_DIR:-}"
if [[ -z "$DB_DIR" && -n "${CONDA_PREFIX:-}" && -d "${CONDA_PREFIX}/db" ]]; then
  DB_DIR="${CONDA_PREFIX}/db"
fi
if [[ -z "$DB_DIR" ]]; then
  ENV_PATH=$(conda env list 2>/dev/null | awk -v env="$DBCAN_CONDA_ENV" '$1==env {print $NF}' | head -1 || true)
  if [[ -z "$ENV_PATH" ]]; then
    ENV_PATH="$(conda info --base 2>/dev/null)/envs/${DBCAN_CONDA_ENV}"
  fi
  DB_DIR="${ENV_PATH}/db"
fi

STATUS="$OUTDIR/dbcan_status.tsv"
MANIFEST="$OUTDIR/dbcan_input_manifest.tsv"
echo -e "sample\tfaa\tstatus\tnote" > "$STATUS"
echo -e "sample\tsource_dir\tfaa" > "$MANIFEST"

write_empty_dbcan() {
  local sampledir="$1"
  mkdir -p "$sampledir"
  # Keep a standard empty dbCAN-like overview for downstream parser compatibility.
  echo -e "Gene ID\tEC#\tHMMER\tDIAMOND\tdbCAN_sub\t#ofTools" > "$sampledir/overview.txt"
}

sanitize_sample() {
  local s="$1"
  s="${s// /_}"
  s="${s//\//_}"
  s="${s//:/_}"
  echo "$s"
}

# Build manifest recursively from all .faa-containing directories.
declare -A seen_samples
while IFS= read -r -d '' faa; do
  sample_dir=$(dirname "$faa")
  base=$(basename "$sample_dir")
  sample=$(sanitize_sample "$base")

  if [[ -n "${seen_samples[$sample]:-}" && "${seen_samples[$sample]}" != "$sample_dir" ]]; then
    parent_safe=$(sanitize_sample "$(basename "$(dirname "$sample_dir")")")
    sample="${parent_safe}__${sample}"
  fi
  if [[ -n "${seen_samples[$sample]:-}" && "${seen_samples[$sample]}" != "$sample_dir" ]]; then
    echo "[dbcan] WARNING: duplicate sample name skipped: $sample_dir" >&2
    continue
  fi
  seen_samples[$sample]="$sample_dir"
  echo -e "${sample}\t${sample_dir}\t${faa}" >> "$MANIFEST"
done < <(find -L "$PROKKA_DIR" -type f -name "*.faa" -print0 | sort -z)

TOTAL=$(tail -n +2 "$MANIFEST" | wc -l)
echo "[dbcan] Input Prokka directory: $PROKKA_DIR"
echo "[dbcan] Detected protein FASTA files: $TOTAL"
echo "[dbcan] Processing with $JOBS parallel jobs; $CPUS_PER threads/job"
echo "[dbcan] Using database: $DB_DIR"
echo "[dbcan] Using conda env: $DBCAN_CONDA_ENV"

DBCAN_ENV_PREFIX=$(conda env list | awk -v env="$DBCAN_CONDA_ENV" '$1==env {print $NF; exit}')
RUN_DBCAN="$DBCAN_ENV_PREFIX/bin/run_dbcan"
if [[ -z "$DBCAN_ENV_PREFIX" || ! -x "$RUN_DBCAN" ]]; then
  echo "[dbcan] ERROR: Cannot find absolute run_dbcan in conda env: $DBCAN_CONDA_ENV" >&2
  echo "[dbcan] Expected: $RUN_DBCAN" >&2
  exit 1
fi
export RUN_DBCAN

if [[ "$TOTAL" -eq 0 ]]; then
  echo "[dbcan] ERROR: No .faa files found recursively in: $PROKKA_DIR" >&2
  echo "[dbcan] Check that existing Prokka folders contain .faa files." >&2
  exit 1
fi

if ! conda run -n "$DBCAN_CONDA_ENV" "$RUN_DBCAN" --help >/dev/null 2>&1; then
  echo "[dbcan] ERROR: run_dbcan is broken in conda env: $DBCAN_CONDA_ENV" >&2
  echo "[dbcan] Absolute path tested: $RUN_DBCAN" >&2
  echo "[dbcan] Install/activate dbCAN or run the pipeline install script." >&2
  exit 1
fi

if [[ ! -d "$DB_DIR" ]]; then
  echo "[dbcan] ERROR: dbCAN database directory not found: $DB_DIR" >&2
  echo "[dbcan] Fix by running one of these:" >&2
  echo "  conda run -n $DBCAN_CONDA_ENV $RUN_DBCAN database --db_dir $DB_DIR" >&2
  echo "  export DBCAN_DB_DIR=/path/to/dbcan/db" >&2
  exit 1
fi

run_dbcan_sample() {
  local sample="$1"
  local faa="$2"
  local sampledir="$OUTDIR/$sample"
  local log="logs/dbcan/${sample}.log"
  mkdir -p "$sampledir"

  if [[ ! -s "$faa" ]]; then
    echo "[dbcan] WARNING: Empty/missing FAA for $sample"
    write_empty_dbcan "$sampledir"
    echo -e "${sample}\t${faa}\tNO_FAA\tMissing or empty protein FASTA" >> "$STATUS"
    return 0
  fi

  if [[ -s "$sampledir/overview.txt" || -s "$sampledir/dbCAN_hmm_results.tsv" ]]; then
    echo "[dbcan] Skipping $sample (existing result)"
    echo -e "${sample}\t${faa}\tSKIPPED\tExisting result present" >> "$STATUS"
    return 0
  fi

  echo "[dbcan] Running: $sample"

  # dbCAN v5 syntax. Current docs use comma-separated methods; some older
  # builds accept repeated --methods flags. Try documented syntax first, then
  # fall back to repeated flags, then finally default methods.
  if conda run -n "$DBCAN_CONDA_ENV" "$RUN_DBCAN" CAZyme_annotation \
      --input_raw_data "$faa" \
      --output_dir "$sampledir" \
      --db_dir "$DB_DIR" \
      --mode protein \
      --methods diamond,hmm,dbCANsub \
      --threads "$CPUS_PER" \
      > "$log" 2>&1; then
    :
  else
    echo "[dbcan] Comma-separated method command failed for $sample; retrying with repeated --methods flags" >> "$log"
    if conda run -n "$DBCAN_CONDA_ENV" "$RUN_DBCAN" CAZyme_annotation \
        --input_raw_data "$faa" \
        --output_dir "$sampledir" \
        --db_dir "$DB_DIR" \
        --mode protein \
        --methods hmm --methods diamond --methods dbCANsub \
        --threads "$CPUS_PER" \
        >> "$log" 2>&1; then
      :
    else
      echo "[dbcan] Repeated method command failed for $sample; retrying with dbCAN default methods" >> "$log"
      conda run -n "$DBCAN_CONDA_ENV" "$RUN_DBCAN" CAZyme_annotation \
        --input_raw_data "$faa" \
        --output_dir "$sampledir" \
        --db_dir "$DB_DIR" \
        --mode protein \
        --threads "$CPUS_PER" \
        >> "$log" 2>&1
    fi
  fi

  if [[ -s "$sampledir/overview.txt" ]]; then
    nlines=$(wc -l < "$sampledir/overview.txt")
    if [[ "$nlines" -gt 1 ]]; then
      echo "[dbcan] Done: $sample ($((nlines-1)) CAZyme rows in overview)"
      echo -e "${sample}\t${faa}\tOK\t${nlines} overview lines" >> "$STATUS"
    else
      echo "[dbcan] WARNING: $sample produced overview header only"
      echo -e "${sample}\t${faa}\tEMPTY_RESULT\tOverview header only; check log" >> "$STATUS"
    fi
  elif [[ -s "$sampledir/dbCAN_hmm_results.tsv" || -s "$sampledir/dbCANsub_hmm_results.tsv" || -s "$sampledir/diamond.out" ]]; then
    echo "[dbcan] Done: $sample (no overview, detailed files present)"
    echo -e "${sample}\t${faa}\tOK_NO_OVERVIEW\tDetailed dbCAN result files present" >> "$STATUS"
  else
    echo "[dbcan] FAILED: $sample — see $log"
    write_empty_dbcan "$sampledir"
    echo -e "${sample}\t${faa}\tFAILED\tNo dbCAN output files; check log" >> "$STATUS"
  fi
}

export -f run_dbcan_sample write_empty_dbcan
export OUTDIR DB_DIR CPUS_PER CONDA_ENV DBCAN_CONDA_ENV STATUS RUN_DBCAN

tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r sample source_dir faa; do
  printf '%s\t%s\n' "$sample" "$faa"
done | xargs -P "$JOBS" -n 2 bash -c 'run_dbcan_sample "$1" "$2"' _

DONE=$(find "$OUTDIR" -mindepth 1 -maxdepth 1 -type d | wc -l)
echo "[dbcan] Output directories: $DONE / $TOTAL"
echo "[dbcan] Status summary:"
tail -n +2 "$STATUS" | cut -f3 | sort | uniq -c || true
echo "[dbcan] Status file: $STATUS"
