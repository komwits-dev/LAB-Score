#!/usr/bin/env bash
# 05_run_checkm.sh  (OPTIONAL - only runs if -q flag passed to run_all.sh)
# Usage: bash 05_run_checkm.sh <indir> <outdir> <threads>

set -euo pipefail
INDIR=$1
OUTDIR=$2
THREADS=$3

mkdir -p "$OUTDIR" logs/checkm

if ! command -v checkm &>/dev/null; then
  echo "[checkm] CheckM not found. Writing empty result."
  echo -e "genome\tCompleteness\tContamination\tStrain heterogeneity" > "$OUTDIR/checkm_summary.tsv"
  exit 0
fi

echo "[checkm] Running lineage workflow on all genomes..."
checkm lineage_wf \
  -t "$THREADS" \
  -x fna \
  --tab_table \
  -f "$OUTDIR/checkm_full.tsv" \
  "$INDIR" \
  "$OUTDIR/checkm_wd" \
  > logs/checkm/checkm.log 2>&1

# Parse to simple summary
python - <<'PY'
import pandas as pd, sys
try:
    df = pd.read_csv("$OUTDIR/checkm_full.tsv", sep="\t", comment="[")
    # CheckM output column names vary by version - handle both
    col_map = {}
    for c in df.columns:
        cl = c.lower()
        if "completeness" in cl:   col_map[c] = "Completeness"
        if "contamination" in cl:  col_map[c] = "Contamination"
        if "bin id" in cl or "genome" in cl: col_map[c] = "genome"
    df = df.rename(columns=col_map)
    df[["genome","Completeness","Contamination"]].to_csv(
        "$OUTDIR/checkm_summary.tsv", sep="\t", index=False)
    print(f"[checkm] Saved summary for {len(df)} genomes.")
except Exception as e:
    print(f"[checkm] Parse error: {e}")
    pd.DataFrame(columns=["genome","Completeness","Contamination"]).to_csv(
        "$OUTDIR/checkm_summary.tsv", sep="\t", index=False)
PY

echo "[checkm] Done."
