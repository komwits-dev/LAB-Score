#!/usr/bin/env bash
# =============================================================================
# 00_prepare_prokka_input.sh
# Import existing Prokka annotations into LAB-Score directory layout.
#
# Usage:
#   bash 00_prepare_prokka_input.sh <existing_prokka_dir> <out_prokka_dir> <out_genome_dir>
#
# Accepted Prokka layouts:
#   1) Direct per-sample folders:
#      existing_prokka_dir/sample_1/*.gff *.faa *.fna
#      existing_prokka_dir/sample_2/*.gff *.faa *.fna
#
#   2) Nested species/accession folders:
#      existing_prokka_dir/Species_name/GCA_xxxxx/*.gff *.faa *.fna
#      existing_prokka_dir/Species_name/GCF_xxxxx/*.gff *.faa *.fna
#
#   3) Flat files:
#      existing_prokka_dir/sample_1.gff sample_1.faa sample_1.fna
#
# The script recursively finds every directory containing at least one Prokka
# output file and symlinks that sample into the LAB-Score output layout.
# =============================================================================

set -euo pipefail

PROKKA_IN=${1:?"existing Prokka directory required"}
PROKKA_OUT=${2:?"output Prokka directory required"}
GENOME_OUT=${3:?"output genome FASTA directory required"}

if [[ ! -d "$PROKKA_IN" ]]; then
  echo "[prepare-prokka] ERROR: Prokka directory not found: $PROKKA_IN" >&2
  exit 1
fi

mkdir -p "$PROKKA_OUT" "$GENOME_OUT"
PROKKA_IN_ABS=$(cd "$PROKKA_IN" && pwd)
PROKKA_OUT_ABS=$(cd "$PROKKA_OUT" && pwd)
GENOME_OUT_ABS=$(cd "$GENOME_OUT" && pwd)

MANIFEST="$PROKKA_OUT_ABS/prokka_import_manifest.tsv"
echo -e "sample\tsource_dir\tgff\tfaa\tfna\tstatus" > "$MANIFEST"

count=0
missing_fna=0
missing_gff=0
missing_faa=0
skipped_duplicate=0

declare -A seen_samples

sanitize_sample() {
  local s="$1"
  s="${s// /_}"
  s="${s//\//_}"
  s="${s//:/_}"
  echo "$s"
}

import_sample_dir() {
  local sample_dir="$1"
  local base sample gff faa fna status out_sample f parent_safe

  base=$(basename "$sample_dir")
  sample=$(sanitize_sample "$base")

  gff=$(find -L "$sample_dir" -maxdepth 1 -type f -name "*.gff" | sort | head -n 1 || true)
  faa=$(find -L "$sample_dir" -maxdepth 1 -type f -name "*.faa" | sort | head -n 1 || true)
  fna=$(find -L "$sample_dir" -maxdepth 1 -type f \( -name "*.fna" -o -name "*.fa" -o -name "*.fasta" \) | sort | head -n 1 || true)

  if [[ -z "$gff" && -z "$faa" && -z "$fna" ]]; then
    return 0
  fi

  # Avoid collisions if two nested directories have the same basename.
  if [[ -n "${seen_samples[$sample]:-}" && "${seen_samples[$sample]}" != "$sample_dir" ]]; then
    parent_safe=$(sanitize_sample "$(basename "$(dirname "$sample_dir")")")
    sample="${parent_safe}__${sample}"
    if [[ -n "${seen_samples[$sample]:-}" ]]; then
      skipped_duplicate=$((skipped_duplicate + 1))
      echo "[prepare-prokka] WARNING: duplicate sample name skipped: $sample_dir" >&2
      return 0
    fi
  fi
  seen_samples[$sample]="$sample_dir"

  status="OK"
  count=$((count + 1))
  out_sample="$PROKKA_OUT_ABS/$sample"
  mkdir -p "$out_sample"

  while IFS= read -r -d '' f; do
    ln -sfn "$f" "$out_sample/$(basename "$f")"
  done < <(find -L "$sample_dir" -maxdepth 1 -type f -print0)

  if [[ -z "$gff" ]]; then
    echo "[prepare-prokka] WARNING: no GFF found for $sample" >&2
    missing_gff=$((missing_gff + 1))
    status="MISSING_GFF"
  fi
  if [[ -z "$faa" ]]; then
    echo "[prepare-prokka] WARNING: no FAA found for $sample" >&2
    missing_faa=$((missing_faa + 1))
    [[ "$status" == "OK" ]] && status="MISSING_FAA" || status="${status};MISSING_FAA"
  fi
  if [[ -z "$fna" ]]; then
    echo "[prepare-prokka] WARNING: no FNA/FASTA found for $sample; ResFinder/metadata may need -i genome FASTA input" >&2
    missing_fna=$((missing_fna + 1))
    [[ "$status" == "OK" ]] && status="MISSING_FNA" || status="${status};MISSING_FNA"
  else
    # Link genome FASTA using the sample ID, so downstream metadata, ResFinder,
    # and Prokka-derived IDs stay aligned.
    ln -sfn "$fna" "$GENOME_OUT_ABS/$sample.fna"
  fi

  echo -e "$sample\t$sample_dir\t${gff:-NA}\t${faa:-NA}\t${fna:-NA}\t$status" >> "$MANIFEST"
}

# Recursively import any directory containing Prokka output files.
# This supports nested layouts such as annotations_prokka_GCA/Species/GCA_*/.
while IFS= read -r -d '' sample_dir; do
  import_sample_dir "$sample_dir"
done < <(find -L "$PROKKA_IN_ABS" -mindepth 1 -type d -print0 | sort -z)

# If no directories were imported, try flat Prokka layout in PROKKA_IN itself.
if [[ "$count" -eq 0 ]]; then
  echo "[prepare-prokka] No per-sample directories detected; trying flat Prokka layout..."
  while IFS= read -r -d '' gff; do
    stem=$(basename "$gff" .gff)
    sample=$(sanitize_sample "$stem")
    out_sample="$PROKKA_OUT_ABS/$sample"
    mkdir -p "$out_sample"

    for ext in gff faa fna fa fasta ffn gbk txt tsv sqn tbl err log fsa; do
      for f in "$PROKKA_IN_ABS/$stem"."$ext"; do
        [[ -f "$f" ]] && ln -sfn "$f" "$out_sample/$(basename "$f")"
      done
    done

    faa="$PROKKA_IN_ABS/$stem.faa"
    fna="$PROKKA_IN_ABS/$stem.fna"
    status="OK"

    if [[ ! -f "$faa" ]]; then
      missing_faa=$((missing_faa + 1))
      status="MISSING_FAA"
    fi
    if [[ -f "$fna" ]]; then
      ln -sfn "$fna" "$GENOME_OUT_ABS/$sample.fna"
    else
      missing_fna=$((missing_fna + 1))
      [[ "$status" == "OK" ]] && status="MISSING_FNA" || status="${status};MISSING_FNA"
    fi

    echo -e "$sample\t$PROKKA_IN_ABS\t$gff\t$([[ -f "$faa" ]] && echo "$faa" || echo NA)\t$([[ -f "$fna" ]] && echo "$fna" || echo NA)\t$status" >> "$MANIFEST"
    count=$((count + 1))
  done < <(find -L "$PROKKA_IN_ABS" -maxdepth 1 -type f -name "*.gff" -print0 | sort -z)
fi

if [[ "$count" -eq 0 ]]; then
  echo "[prepare-prokka] ERROR: No usable Prokka outputs found in $PROKKA_IN" >&2
  echo "[prepare-prokka] Expected Prokka directories or flat files containing .gff/.faa/.fna files." >&2
  exit 1
fi

echo "[prepare-prokka] Imported Prokka samples: $count"
echo "[prepare-prokka] Prokka output dir: $PROKKA_OUT_ABS"
echo "[prepare-prokka] Genome FASTA dir:  $GENOME_OUT_ABS"
echo "[prepare-prokka] Manifest:          $MANIFEST"
echo "[prepare-prokka] Missing GFF: $missing_gff | Missing FAA: $missing_faa | Missing FNA: $missing_fna | Duplicate skipped: $skipped_duplicate"
