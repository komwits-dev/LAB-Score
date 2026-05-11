#!/usr/bin/env bash
# 01_run_prokka.sh
set -uo pipefail
# Prevent base/old Miniconda Perl variables from breaking Prokka.
unset PERL5LIB PERL_LOCAL_LIB_ROOT PERL_MB_OPT PERL_MM_OPT

INDIR=$1
OUTDIR=$2
THREADS=$3
JOBS=$(( THREADS > 4 ? 4 : THREADS ))
CPUS_PER=$(( THREADS / JOBS ))
[[ $CPUS_PER -lt 1 ]] && CPUS_PER=1

mkdir -p "$OUTDIR" logs/prokka

run_prokka() {
  local fasta=$1
  local full_name
  full_name=$(basename "$fasta" | sed 's/\.\(fna\|fa\|fasta\)$//')
  local prefix
  prefix=$(echo "$full_name" | cut -d'_' -f1-2 | cut -c1-20)
  local outdir="$OUTDIR/$full_name"

  if [[ -n "$(find "$outdir" -name "*.gff" 2>/dev/null)" ]]; then
    echo "[prokka] Skipping $full_name (already done)"
    return 0
  fi

  mkdir -p "$outdir"
  echo "[prokka] Annotating: $full_name"
  prokka \
    --outdir  "$outdir" \
    --prefix  "$prefix" \
    --cpus    "$CPUS_PER" \
    --force   \
    --quiet   \
    "$fasta" \
    > "logs/prokka/${full_name}.log" 2>&1 \
  && echo "[prokka] Done: $full_name" \
  || echo "[prokka] FAILED: $full_name — see logs/prokka/${full_name}.log"
}

export -f run_prokka
export OUTDIR CPUS_PER

GENOMES=$(find "$INDIR" -maxdepth 1 \( -name "*.fna" -o -name "*.fa" -o -name "*.fasta" \) | sort)
TOTAL=$(echo "$GENOMES" | wc -l)
echo "[prokka] Processing $TOTAL genomes with $JOBS parallel jobs ($CPUS_PER threads each)"
echo "$GENOMES" | xargs -P "$JOBS" -I{} bash -c 'run_prokka "$@"' _ {}
DONE=$(find "$OUTDIR" -name "*.gff" | wc -l)
echo "[prokka] Completed: $DONE / $TOTAL"
