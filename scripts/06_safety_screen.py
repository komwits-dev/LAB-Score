#!/usr/bin/env python3
"""
06_safety_screen.py — Comprehensive safety screening
Combines AMRFinder + ResFinder + GFF biogenic amine + hemolysin scan
"""

import argparse, os, re
import pandas as pd
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("--prokka_dir",  required=True)
ap.add_argument("--amrfinder",   required=True)
ap.add_argument("--resfinder",   required=True)
ap.add_argument("--outdir",      required=True)
args = ap.parse_args()
os.makedirs(args.outdir, exist_ok=True)

BA_PATTERNS = re.compile(
    r"histidine.decarboxylase|tyrosine.decarboxylase|ornithine.decarboxylase|"
    r"lysine.decarboxylase|arginine.decarboxylase|glutamate.decarboxylase|"
    r"\bhdc[A-Za-z]?\b|\btdc[A-Za-z]?\b|\bodc[A-Za-z]?\b|\badc[A-Za-z]?\b|"
    r"amino.acid.decarboxylase|pyridoxal.phosphate.dependent.decarboxylase",
    re.IGNORECASE
)
HLY_PATTERNS = re.compile(
    r"\bhemolysin\b|\bhaemolysin\b|\bcytolysin\b|\bhly[A-Za-z]?\b|\bvhh\b|"
    r"pore.forming.toxin|beta.hemolysin|alpha.hemolysin",
    re.IGNORECASE
)
AMR_TYPES       = {"AMR", "STRESS", "POINT"}
VIRULENCE_TYPES = {"VIRULENCE"}

def parse_gff_safety(gff_path):
    ba_hits = hly_hits = 0
    try:
        with open(gff_path) as fh:
            for line in fh:
                if line.startswith(("#", ">")) or not line.strip():
                    continue
                parts = line.split("\t")
                if len(parts) < 9:
                    continue
                attrs = parts[8]
                product = gene = ""
                for field in attrs.split(";"):
                    fl = field.lower()
                    if fl.startswith("product="):
                        product = field.split("=", 1)[1].strip()
                    elif fl.startswith("gene="):
                        gene = field.split("=", 1)[1].strip()
                text = f"{product} {gene}"
                if BA_PATTERNS.search(text):  ba_hits  += 1
                if HLY_PATTERNS.search(text): hly_hits += 1
    except Exception as e:
        print(f"[safety] GFF error {gff_path}: {e}")
    return ba_hits, hly_hits

def parse_amrfinder(tsv_path):
    amr_hits = vf_hits = 0
    try:
        df = pd.read_csv(tsv_path, sep="\t")
        for col in ["Element type", "Element_type"]:
            if col in df.columns:
                amr_hits = int(df[col].isin(AMR_TYPES).sum())
                vf_hits  = int(df[col].isin(VIRULENCE_TYPES).sum())
                break
    except Exception:
        pass
    return amr_hits, vf_hits

def parse_resfinder_summary(tsv_path):
    result = {}
    try:
        df = pd.read_csv(tsv_path, sep="\t", dtype=str)
        for _, row in df.iterrows():
            g = str(row.get("genome", "")).strip()
            result[g] = {
                "resfinder_hits":  (lambda x: 0 if (x is None or (isinstance(x, float) and x != x)) else int(float(str(x))))(row.get("resfinder_hits", 0)),
                "resfinder_genes": str(row.get("resfinder_genes", "") or ""),
            }
    except Exception as e:
        print(f"[safety] ResFinder parse error: {e}")
    return result

# ── Discover genomes ──────────────────────────────────────────────────────────
prokka_dir = Path(args.prokka_dir)
samples = sorted([d.name for d in prokka_dir.iterdir() if d.is_dir()]) \
          if prokka_dir.exists() else []
print(f"[safety] Processing {len(samples)} genomes")

resfinder_summary = parse_resfinder_summary(
    os.path.join(args.resfinder, "resfinder_summary.tsv")
)

rows = []
for sample in samples:
    # Find GFF — prefix may differ from folder name
    gff_files = list(Path(args.prokka_dir, sample).glob("*.gff"))
    gff_path  = gff_files[0] if gff_files else None

    # Find AMRFinder result
    amr_path = Path(args.amrfinder) / f"{sample}.tsv"

    ba_hits, hly_hits = parse_gff_safety(str(gff_path)) if gff_path else (0, 0)
    amr_hits, vf_hits = parse_amrfinder(str(amr_path))  if amr_path.exists() else (0, 0)
    res = resfinder_summary.get(sample, {"resfinder_hits": 0, "resfinder_genes": ""})

    if amr_hits > 0 or vf_hits > 0 or res["resfinder_hits"] > 0:
        flag = "Critical"
        reasons = []
        if amr_hits > 0:              reasons.append("AMR_gene_detected")
        if vf_hits  > 0:              reasons.append("Virulence_gene_detected")
        if res["resfinder_hits"] > 0: reasons.append("ResFinder_hit")
        reason_str = ";".join(reasons)
    elif hly_hits > 0 or ba_hits > 0:
        flag = "Caution"
        reasons = []
        if hly_hits > 0: reasons.append("Hemolysin_marker")
        if ba_hits  > 0: reasons.append("Biogenic_amine_gene")
        reason_str = ";".join(reasons)
    else:
        flag = "Pass"
        reason_str = ""

    rows.append({
        "genome":            sample,
        "amrfinder_hits":    amr_hits,
        "vfdb_hits":         vf_hits,
        "resfinder_hits":    res["resfinder_hits"],
        "resfinder_genes":   res["resfinder_genes"],
        "biogenic_amines":   ba_hits,
        "hemolysin_markers": hly_hits,
        "flag":              flag,
        "reasons":           reason_str,
    })

# ── Handle empty case gracefully ──────────────────────────────────────────────
if not rows:
    print("[safety] WARNING: No genomes found — writing empty safety summary.")
    df_out = pd.DataFrame(columns=[
        "genome","amrfinder_hits","vfdb_hits","resfinder_hits",
        "resfinder_genes","biogenic_amines","hemolysin_markers","flag","reasons"
    ])
else:
    df_out = pd.DataFrame(rows)
    print("[safety] Flag summary:")
    print(df_out["flag"].value_counts().to_string())

out_path = os.path.join(args.outdir, "safety_summary.tsv")
df_out.to_csv(out_path, sep="\t", index=False)
print(f"[safety] Saved: {out_path}")
