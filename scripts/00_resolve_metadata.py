#!/usr/bin/env python3
"""
00_resolve_metadata.py
======================
Resolves species, genus, and strain for each genome.

Priority order:
  1. User-supplied metadata TSV  (--metadata)
  2. NCBI Datasets API           (auto for GCF_/GCA_ accessions)
  3. Parse from filename          (fallback for any naming scheme)

Output:
  {outdir}/genome_metadata.tsv
"""

import argparse, os, re, json, time, sys
import pandas as pd
from pathlib import Path
from urllib.request import urlopen, Request
from urllib.error import URLError

ap = argparse.ArgumentParser()
ap.add_argument("--indir",    required=True)
ap.add_argument("--outdir",   required=True)
ap.add_argument("--metadata", default="")
ap.add_argument("--ncbi",     action="store_true")
args = ap.parse_args()

os.makedirs(args.outdir, exist_ok=True)

# ── Discover genomes ──────────────────────────────────────────────────────────
EXTS = [".fna", ".fa", ".fasta"]
genome_files = sorted([p for p in Path(args.indir).iterdir() if p.suffix in EXTS])
genomes = [p.stem for p in genome_files]
print(f"[metadata] Found {len(genomes)} genome files.")

# ── Load user metadata ────────────────────────────────────────────────────────
user_meta = {}
if args.metadata and os.path.exists(args.metadata):
    um = pd.read_csv(args.metadata, sep="\t", dtype=str)
    um.columns = [c.lower().strip() for c in um.columns]
    for _, row in um.iterrows():
        rowd = row.to_dict()
        # Allow matching either by genome filename stem or by accession.
        for key_col in ["genome", "accession"]:
            key = str(row.get(key_col, "") or "").strip()
            if key:
                user_meta[key] = rowd
    print(f"[metadata] Loaded {len(user_meta)} metadata keys from user metadata.")

# ── GCA/GCF accession regex ───────────────────────────────────────────────────
GCF_RE = re.compile(r"(GC[FA]_\d+\.\d+)")

# ── NCBI Datasets API fetch ───────────────────────────────────────────────────
def fetch_ncbi_metadata(accessions):
    """Fetch species info for GCF/GCA accessions from NCBI Datasets API."""
    results = {}
    BATCH   = 200
    base_url = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession"

    for i in range(0, len(accessions), BATCH):
        batch = accessions[i:i+BATCH]
        ids   = ",".join(batch)
        url   = f"{base_url}/{ids}/dataset_report?page_size={BATCH}"
        try:
            req  = Request(url, headers={"Accept": "application/json",
                                         "User-Agent": "LAB-Score-Pipeline/2.0"})
            resp = urlopen(req, timeout=30)
            data = json.loads(resp.read())
            for item in data.get("reports", []):
                acc  = item.get("accession", "")
                info = item.get("assembly_info", {})
                org  = item.get("organism", {})
                seq  = item.get("assembly_stats", {})
                name = org.get("organism_name", "")
                parts = name.split(" ", 2)
                results[acc] = {
                    "genus":          parts[0] if len(parts) > 0 else "Unknown",
                    "species":        parts[1] if len(parts) > 1 else "sp.",
                    "full_species":   " ".join(parts[:2]) if len(parts) >= 2 else name or "Unknown sp.",
                    "strain":         org.get("infraspecific_names", {}).get("strain", acc),
                    "assembly_level": info.get("assembly_level", ""),
                    "total_length":   seq.get("total_sequence_length", 0),
                    "gc_pct":         seq.get("gc_percent", 0),
                    "contigs":        seq.get("number_of_contigs", 0),
                }
            time.sleep(0.35)  # NCBI rate limit
        except Exception as e:
            print(f"[metadata] NCBI fetch error (batch {i//BATCH+1}): {e}", file=sys.stderr)
    return results

# ── Filename-based parsing (fallback) ────────────────────────────────────────
def parse_from_filename(stem):
    m = GCF_RE.search(stem)
    accession = m.group(1) if m else stem
    parts = stem.replace("-", "_").split("_")
    if len(parts) >= 2 and parts[0][0].isupper() and len(parts[0]) > 2 and parts[1][0].islower():
        genus   = parts[0]
        species = parts[1]
        strain  = "_".join(parts[2:]) if len(parts) > 2 else stem
    else:
        genus   = "Unknown"
        species = "sp."
        strain  = stem
    return {
        "genus":          genus,
        "species":        species,
        "full_species":   f"{genus} {species}",
        "strain":         strain,
        "accession":      accession,
        "assembly_level": "",
        "total_length":   0,
        "gc_pct":         0,
        "contigs":        0,
    }

# ── Auto-detect GCF/GCA accessions for NCBI lookup ───────────────────────────
# Always try NCBI for GCF/GCA accessions unless user explicitly opted out
ncbi_candidates = {}
for g in genomes:
    m = GCF_RE.search(g)
    if m:
        ncbi_candidates[g] = m.group(1)

ncbi_meta = {}
if ncbi_candidates:
    print(f"[metadata] Detected {len(ncbi_candidates)} GCF/GCA accessions — fetching from NCBI...")
    try:
        ncbi_meta = fetch_ncbi_metadata(list(set(ncbi_candidates.values())))
        print(f"[metadata] NCBI returned {len(ncbi_meta)} records.")
    except Exception as e:
        print(f"[metadata] NCBI fetch failed: {e} — falling back to filename parsing.")
elif args.ncbi:
    print(f"[metadata] --ncbi flag set but no GCF/GCA accessions detected in filenames.")

# ── Build final metadata table ────────────────────────────────────────────────
rows = []
resolved_ncbi = 0
resolved_user = 0
resolved_file = 0

for stem in genomes:
    row = {"genome": stem}

    # 1. User metadata (highest priority)
    accession_guess = ncbi_candidates.get(stem, parse_from_filename(stem).get("accession", stem))

    if stem in user_meta or accession_guess in user_meta:
        u = user_meta.get(stem, user_meta.get(accession_guess))
        row.update({
            "accession":      u.get("accession", stem),
            "genus":          u.get("genus", "Unknown"),
            "species":        u.get("species", "sp."),
            "full_species":   u.get("full_species", f"{u.get('genus','Unknown')} {u.get('species','sp.')}"),
            "strain":         u.get("strain", stem),
            "assembly_level": u.get("assembly_level", ""),
            "total_length":   u.get("total_length", 0),
            "gc_pct":         u.get("gc_pct", 0),
            "contigs":        u.get("contigs", 0),
        })
        resolved_user += 1

    # 2. NCBI metadata (auto for GCF/GCA)
    elif stem in ncbi_candidates and ncbi_candidates[stem] in ncbi_meta:
        acc = ncbi_candidates[stem]
        row.update(ncbi_meta[acc])
        row["accession"] = acc
        resolved_ncbi += 1

    # 3. Filename parsing (fallback)
    else:
        row.update(parse_from_filename(stem))
        if ncbi_candidates.get(stem):
            # Had a GCF accession but NCBI didn't return it
            row["accession"] = ncbi_candidates[stem]
        resolved_file += 1

    rows.append(row)

meta_df = pd.DataFrame(rows)

# Ensure all expected columns exist
for col in ["genome","accession","genus","species","full_species","strain",
            "assembly_level","total_length","gc_pct","contigs"]:
    if col not in meta_df.columns:
        meta_df[col] = ""

# Clean malformed or incomplete species metadata.
def _clean_text(x):
    if pd.isna(x):
        return ""
    return str(x).strip()

for c in ["genus", "species", "full_species", "strain", "accession"]:
    if c in meta_df.columns:
        meta_df[c] = meta_df[c].map(_clean_text)

for i, r in meta_df.iterrows():
    genus = r.get("genus", "") or "Unknown"
    species = r.get("species", "") or "sp."
    full = r.get("full_species", "") or ""
    if not full or full.lower() in {"nan", "none"}:
        full = f"{genus} {species}".strip()
    parts = full.split()
    if len(parts) >= 2 and genus != "Unknown":
        if species.startswith(genus) or species == "sp." or species == "":
            species = parts[1]
    meta_df.at[i, "genus"] = genus
    meta_df.at[i, "species"] = species
    meta_df.at[i, "full_species"] = f"{genus} {species}".strip() if genus != "Unknown" else full

out_path = f"{args.outdir}/genome_metadata.tsv"
meta_df.to_csv(out_path, sep="\t", index=False)

print(f"[metadata] Resolved via: NCBI={resolved_ncbi}  User={resolved_user}  Filename={resolved_file}")
print(f"[metadata] Saved: {out_path}")
print(f"[metadata] Species resolved:")
sp_counts = meta_df["full_species"].value_counts()
for sp, n in sp_counts.head(10).items():
    print(f"  {n:>4}  {sp}")
if len(sp_counts) > 10:
    print(f"  ... and {len(sp_counts)-10} more species")
