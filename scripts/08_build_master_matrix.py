#!/usr/bin/env python3
"""
08_build_master_matrix.py
=========================
Merges all per-genome results into a single master matrix.

v2.6 update:
- Robust dbCAN parsing from overview.txt, dbCAN_hmm_results.tsv,
  dbCANsub_hmm_results.tsv, and diamond.out.
- CAZyme_total now counts all detected CAZy families, not only a small fixed list.
- dbcan_status.tsv is merged when available so failed/empty runs are visible.
"""

import argparse
import os
import re
from pathlib import Path
from collections import defaultdict
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--metadata", required=True)
ap.add_argument("--safety",   required=True)
ap.add_argument("--panels",   required=True)
ap.add_argument("--dbcan",    required=True)
ap.add_argument("--checkm",   required=True)
ap.add_argument("--outdir",   required=True)
args = ap.parse_args()
os.makedirs(args.outdir, exist_ok=True)

# Common CAZy family prefixes. We track all families found dynamically.
FAMILY_RE = re.compile(r"\b(?:GH|GT|PL|CE|AA|CBM)\d+(?:[_-][A-Za-z0-9]+)?\b", re.IGNORECASE)
BASE_FAMILY_RE = re.compile(r"^(GH|GT|PL|CE|AA|CBM)(\d+)", re.IGNORECASE)

# Keep these columns even when absent, because some downstream analyses may expect them.
CORE_CAZY_FAMILIES = [
    "AA0","AA10","CBM13","CBM32","CBM48","CBM50","CBM51",
    "CE4","CE7","CE8","CE12","CE14",
    "GH1","GH2","GH3","GH13","GH15","GH18","GH20","GH23",
    "GH24","GH25","GH27","GH28","GH29","GH31","GH32","GH33",
    "GH35","GH36","GH38","GH39","GH42","GH43","GH51","GH53",
    "GH65","GH66","GH67","GH68","GH70","GH73","GH78","GH84",
    "GH85","GH88","GH91","GH92","GH93","GH94","GH105","GH106",
    "GT1","GT2","GT4","GT11","GT14","GT28","GT30","GT35",
    "GT47","GT51","GT57","GT58","GT83","GT92","PL8",
]


def normalize_family(token: str) -> str | None:
    """Convert GH13_e1/GH13-foo/GH13 to base family GH13."""
    if token is None:
        return None
    token = str(token).strip()
    m = BASE_FAMILY_RE.search(token)
    if not m:
        return None
    return f"{m.group(1).upper()}{m.group(2)}"


def extract_families_from_text(text: str) -> list[str]:
    if not isinstance(text, str):
        text = str(text)
    fams = []
    for m in FAMILY_RE.finditer(text):
        fam = normalize_family(m.group(0))
        if fam:
            fams.append(fam)
    return fams


def parse_overview(overview_path: Path) -> tuple[dict, int]:
    """Parse dbCAN overview.txt. Count each gene-family pair once."""
    counts = defaultdict(int)
    total_gene_family = 0
    if not overview_path.exists() or overview_path.stat().st_size == 0:
        return counts, 0

    try:
        df = pd.read_csv(overview_path, sep="\t", dtype=str)
    except Exception:
        return counts, 0

    if df.empty:
        return counts, 0

    # Columns expected from dbCAN include Gene ID, EC#, HMMER, DIAMOND,
    # dbCAN_sub, #ofTools. We scan all non-ID columns to be robust.
    ignore = {"gene id", "geneid", "gene", "ec#", "ec", "#oftools", "numoftools"}
    scan_cols = [c for c in df.columns if re.sub(r"\s+", "", c.lower()) not in ignore]
    if not scan_cols:
        scan_cols = list(df.columns)

    for _, row in df.iterrows():
        row_fams = set()
        for col in scan_cols:
            val = row.get(col, "")
            if pd.isna(val):
                continue
            sval = str(val)
            if sval.strip() in {"", "-", "NA", "nan", "None"}:
                continue
            for fam in extract_families_from_text(sval):
                row_fams.add(fam)
        for fam in row_fams:
            counts[fam] += 1
            total_gene_family += 1

    return counts, total_gene_family


def parse_raw_result_file(path: Path) -> tuple[dict, int]:
    """Fallback parser for detailed dbCAN outputs. Counts family occurrences by line."""
    counts = defaultdict(int)
    total = 0
    if not path.exists() or path.stat().st_size == 0:
        return counts, 0

    try:
        with open(path, "r", errors="ignore") as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                # Avoid counting header-only files.
                if "Gene ID" in line and ("HMMER" in line or "DIAMOND" in line):
                    continue
                fams = set(extract_families_from_text(line))
                for fam in fams:
                    counts[fam] += 1
                    total += 1
    except Exception:
        pass
    return counts, total


def merge_counts(a: dict, b: dict) -> dict:
    out = defaultdict(int)
    for d in (a, b):
        for k, v in d.items():
            out[k] += int(v)
    return out


def parse_dbcan_sample(sample_dir: Path) -> tuple[dict, int, str]:
    """Return family counts, total, and parse source."""
    overview = sample_dir / "overview.txt"
    counts, total = parse_overview(overview)
    if total > 0:
        return counts, total, "overview.txt"

    # Fallback: parse detailed output files if overview is missing/empty.
    fallback_files = [
        sample_dir / "dbCAN_hmm_results.tsv",
        sample_dir / "dbCANsub_hmm_results.tsv",
        sample_dir / "diamond.out",
    ]
    source = []
    merged = defaultdict(int)
    total2 = 0
    for f in fallback_files:
        c, t = parse_raw_result_file(f)
        if t > 0:
            source.append(f.name)
            merged = merge_counts(merged, c)
            total2 += t
    return merged, total2, "+".join(source) if source else "none"


# ── Load all data ─────────────────────────────────────────────────────────────
print("[matrix] Loading metadata...")
meta = pd.read_csv(args.metadata, sep="\t", dtype=str)

print("[matrix] Loading safety results...")
safety = pd.read_csv(args.safety, sep="\t", dtype=str)

print("[matrix] Loading panel results...")
panels = pd.read_csv(args.panels, sep="\t")

print("[matrix] Loading CheckM results...")
try:
    checkm = pd.read_csv(args.checkm, sep="\t")
except Exception:
    checkm = pd.DataFrame(columns=["genome", "Completeness", "Contamination"])

# ── Parse dbCAN ───────────────────────────────────────────────────────────────
print("[matrix] Parsing dbCAN results...")
dbcan_dir = Path(args.dbcan)
cazy_rows = []
all_families = set(CORE_CAZY_FAMILIES)

if dbcan_dir.exists():
    for sample_dir in sorted(dbcan_dir.iterdir()):
        if not sample_dir.is_dir():
            continue
        sample = sample_dir.name
        counts, total, source = parse_dbcan_sample(sample_dir)
        all_families.update(counts.keys())
        row = {"genome": sample, "CAZyme_total": int(total), "dbcan_parse_source": source}
        for fam, val in counts.items():
            row[fam] = int(val)
        cazy_rows.append(row)

if cazy_rows:
    cazy = pd.DataFrame(cazy_rows)
else:
    cazy = pd.DataFrame(columns=["genome", "CAZyme_total", "dbcan_parse_source"])

# Ensure family columns exist and fill missing family values with zero.
for fam in sorted(all_families, key=lambda x: (re.match(r"[A-Z]+", x).group(0), int(re.search(r"\d+", x).group(0)) if re.search(r"\d+", x) else 0)):
    if fam not in cazy.columns:
        cazy[fam] = 0
family_cols = [c for c in cazy.columns if BASE_FAMILY_RE.match(str(c))]
if family_cols:
    cazy[family_cols] = cazy[family_cols].fillna(0).astype(int)
if "CAZyme_total" in cazy.columns:
    cazy["CAZyme_total"] = cazy["CAZyme_total"].fillna(0).astype(int)

# Merge dbCAN execution status if available.
status_path = dbcan_dir / "dbcan_status.tsv"
if status_path.exists():
    try:
        status = pd.read_csv(status_path, sep="\t", dtype=str)
        status = status.rename(columns={"status": "dbcan_status", "note": "dbcan_note"})
        if "sample" in status.columns:
            status = status.rename(columns={"sample": "genome"})
        keep = [c for c in ["genome", "dbcan_status", "dbcan_note"] if c in status.columns]
        status = status[keep].drop_duplicates("genome", keep="last")
        cazy = cazy.merge(status, on="genome", how="left")
    except Exception as e:
        print(f"[matrix] WARNING: could not read dbCAN status: {e}")
else:
    cazy["dbcan_status"] = ""
    cazy["dbcan_note"] = ""

# Diagnostic export for CAZyme results.
cazy_summary_path = os.path.join(args.outdir, "dbcan_cazy_summary.tsv")
cazy.to_csv(cazy_summary_path, sep="\t", index=False)
print(f"[matrix] dbCAN parsed rows: {len(cazy)}; saved summary: {cazy_summary_path}")
if len(cazy) > 0:
    n_positive = int((cazy["CAZyme_total"] > 0).sum())
    print(f"[matrix] Genomes with CAZyme_total > 0: {n_positive} / {len(cazy)}")

# ── Merge ─────────────────────────────────────────────────────────────────────
print("[matrix] Merging all results...")

# Start with panels (has all genomes)
merged = panels.copy()

# Add metadata
merged = merged.merge(meta, on="genome", how="left")

# Add safety
merged = merged.merge(safety, on="genome", how="left")

# Add CAZymes
merged = merged.merge(cazy, on="genome", how="left")

# Add CheckM
if "genome" in checkm.columns:
    keep_cols = [c for c in ["genome", "Completeness", "Contamination"] if c in checkm.columns]
    merged = merged.merge(checkm[keep_cols], on="genome", how="left")
else:
    merged["Completeness"] = None
    merged["Contamination"] = None

# Fill numeric NAs with 0
num_cols = merged.select_dtypes(include=["number"]).columns
merged[num_cols] = merged[num_cols].fillna(0)

# Fill string NAs
str_cols = merged.select_dtypes(include=["object"]).columns
merged[str_cols] = merged[str_cols].fillna("")

print(f"[matrix] Final shape: {merged.shape}")
print(f"[matrix] Columns ({len(merged.columns)}):")
print("  ", list(merged.columns[:20]), "...")

out_path = os.path.join(args.outdir, "master_matrix.tsv")
merged.to_csv(out_path, sep="\t", index=False)
print(f"[matrix] Saved: {out_path}")
