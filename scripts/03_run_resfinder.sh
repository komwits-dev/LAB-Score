#!/usr/bin/env bash
# 03_run_resfinder.sh
# Runs ResFinder on all genome FASTA files
# Usage: bash 03_run_resfinder.sh <indir> <outdir> <threads>

set -euo pipefail
INDIR=$1
OUTDIR=$2
THREADS=$3
JOBS=$(( THREADS > 4 ? 4 : THREADS ))

mkdir -p "$OUTDIR" logs/resfinder

# Check ResFinder is available
if ! command -v resfinder &>/dev/null && ! python -c "import resfinder" &>/dev/null; then
  echo "[resfinder] WARNING: ResFinder not found. Creating empty results."
  # Write empty summary so pipeline continues
  echo -e "genome\tresfinder_hits\tresfinder_genes" > "$OUTDIR/resfinder_summary.tsv"
  GENOMES=$(find "$INDIR" -maxdepth 1 \( -name "*.fna" -o -name "*.fa" -o -name "*.fasta" \) | sort)
  while IFS= read -r f; do
    sample=$(basename "$f" | sed 's/\.\(fna\|fa\|fasta\)$//')
    echo -e "$sample\t0\t" >> "$OUTDIR/resfinder_summary.tsv"
  done <<< "$GENOMES"
  exit 0
fi

run_resfinder() {
  local fasta=$1
  local sample
  sample=$(basename "$fasta" | sed 's/\.\(fna\|fa\|fasta\)$//')
  local sampledir="$OUTDIR/$sample"

  if [[ -f "$sampledir/ResFinder_results_tab.txt" ]]; then
    echo "[resfinder] Skipping $sample (already done)"
    return 0
  fi

  mkdir -p "$sampledir"
  echo "[resfinder] Running: $sample"

  python -m resfinder \
    -i  "$fasta" \
    -o  "$sampledir" \
    -s  "Other" \
    --acquired \
    --db_path_res "$(python -c 'import resfinder; import os; print(os.path.dirname(resfinder.__file__))')/db_resfinder" \
    > "logs/resfinder/${sample}.log" 2>&1 \
  && echo "[resfinder] Done: $sample" \
  || {
    echo "[resfinder] FAILED: $sample"
    echo "" > "$sampledir/ResFinder_results_tab.txt"
  }
}

export -f run_resfinder
export OUTDIR

GENOMES=$(find "$INDIR" -maxdepth 1 \( -name "*.fna" -o -name "*.fa" -o -name "*.fasta" \) | sort)
TOTAL=$(echo "$GENOMES" | wc -l)
echo "[resfinder] Processing $TOTAL genomes with $JOBS parallel jobs"

echo "$GENOMES" | xargs -P "$JOBS" -I{} bash -c 'run_resfinder "$@"' _ {}

# Summarise all results
echo -e "genome\tresfinder_hits\tresfinder_genes" > "$OUTDIR/resfinder_summary.tsv"
while IFS= read -r f; do
  sample=$(basename "$f" | sed 's/\.\(fna\|fa\|fasta\)$//')
  res_file="$OUTDIR/$sample/ResFinder_results_tab.txt"
  if [[ -f "$res_file" ]]; then
    hits=$(grep -v "^#" "$res_file" 2>/dev/null | grep -v "^$" | wc -l || echo 0)
    genes=$(grep -v "^#" "$res_file" 2>/dev/null | awk -F'\t' '{print $1}' | sort -u | tr '\n' ',' | sed 's/,$//' || echo "")
  else
    hits=0; genes=""
  fi
  echo -e "$sample\t$hits\t$genes"
done <<< "$GENOMES" >> "$OUTDIR/resfinder_summary.tsv"

echo "[resfinder] Summary saved: $OUTDIR/resfinder_summary.tsv"
