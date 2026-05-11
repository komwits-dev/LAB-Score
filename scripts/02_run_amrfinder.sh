#!/usr/bin/env bash
# 02_run_amrfinder.sh
# Runs AMRFinderPlus on Prokka protein FASTA files.
# Sample IDs are taken from Prokka sample directory names, not internal Prokka file prefixes.

set -uo pipefail
PROKKA_DIR=${1:?"Prokka directory required"}
OUTDIR=${2:?"output directory required"}
THREADS=${3:-8}
JOBS=$(( THREADS > 8 ? 8 : THREADS ))
[[ $JOBS -lt 1 ]] && JOBS=1

mkdir -p "$OUTDIR" logs/amrfinder

AMR_CMD=""
for cmd in amrfinder amrfinderplus ncbi-amrfinderplus; do
  if command -v "$cmd" &>/dev/null; then
    AMR_CMD="$cmd"
    break
  fi
done

mapfile -t SAMPLE_DIRS < <(find "$PROKKA_DIR" -mindepth 1 -maxdepth 1 -type d | sort)

write_empty_amr() {
  local out="$1"
  echo -e "Name\tProtein identifier\tContig id\tStart\tStop\tStrand\tGene symbol\tSequence name\tScope\tElement type\tElement subtype\tClass\tSubclass\tMethod\tTarget length\tReference sequence length\t% Coverage of reference sequence\t% Identity to reference sequence\tAlignment length\tAccession of closest sequence\tName of closest sequence\tHMM id\tHMM description" > "$out"
}

if [[ -z "$AMR_CMD" ]]; then
  echo "[amrfinder] WARNING: AMRFinder not found. Writing empty results."
  for sample_dir in "${SAMPLE_DIRS[@]}"; do
    sample=$(basename "$sample_dir")
    write_empty_amr "$OUTDIR/${sample}.tsv"
  done
  exit 0
fi

echo "[amrfinder] Using command: $AMR_CMD"

if [[ ! -f "$OUTDIR/.db_updated" ]]; then
  echo "[amrfinder] Updating database..."
  $AMR_CMD --update 2>/dev/null && touch "$OUTDIR/.db_updated" || echo "[amrfinder] DB update skipped"
fi

run_amrfinder_sample() {
  local sample_dir=$1
  local sample
  sample=$(basename "$sample_dir")
  local faa
  faa=$(find -L "$sample_dir" -maxdepth 1 -name "*.faa" | sort | head -n 1 || true)
  local out="$OUTDIR/${sample}.tsv"

  if [[ -z "$faa" ]]; then
    echo "[amrfinder] WARNING: No .faa found for $sample — writing empty result"
    write_empty_amr "$out"
    return 0
  fi

  if [[ -f "$out" && -s "$out" ]]; then
    echo "[amrfinder] Skipping $sample (already done)"
    return 0
  fi

  echo "[amrfinder] Running: $sample"
  $AMR_CMD --protein "$faa" --threads 4 --plus --output "$out" \
    > "logs/amrfinder/${sample}.log" 2>&1 \
  && echo "[amrfinder] Done: $sample" \
  || {
    echo "[amrfinder] FAILED: $sample — writing empty result"
    write_empty_amr "$out"
  }
}

export -f run_amrfinder_sample write_empty_amr
export OUTDIR AMR_CMD

TOTAL=${#SAMPLE_DIRS[@]}
echo "[amrfinder] Processing $TOTAL genomes"
printf '%s\0' "${SAMPLE_DIRS[@]}" | xargs -0 -P "$JOBS" -I{} bash -c 'run_amrfinder_sample "$@"' _ {}
DONE=$(find "$OUTDIR" -name "*.tsv" | wc -l)
echo "[amrfinder] Completed: $DONE / $TOTAL"
