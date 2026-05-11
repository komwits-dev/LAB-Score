#!/usr/bin/env python3
"""
00b_species_verify.py
=====================
Genome-based species verification using FastANI.

Compares each genome against a curated LAB reference database
and flags mismatches with NCBI-assigned names.

Outputs 3 new columns added to genome_metadata.tsv:
  ani_species     — top FastANI species match
  ani_identity    — ANI % identity (≥95% = same species)
  species_status  — CONFIRMED / MISMATCH / NO_NCBI / LOW_ANI / UNRESOLVED

Usage:
  python 00b_species_verify.py
    --indir       genomes/
    --metadata    results/00_metadata/genome_metadata.tsv
    --ref_dir     /path/to/lab_ref_genomes/   (or --download_refs)
    --outdir      results/00_metadata/
    --threads     8
    --download_refs   (auto-download LAB reference genomes from NCBI)
"""

import argparse, os, re, subprocess, sys, json, time
import pandas as pd
from pathlib import Path
from urllib.request import urlopen, Request

ap = argparse.ArgumentParser()
ap.add_argument("--indir",         required=True,  help="Input genome FASTA directory")
ap.add_argument("--metadata",      required=True,  help="genome_metadata.tsv from step 00")
ap.add_argument("--outdir",        required=True,  help="Output directory")
ap.add_argument("--ref_dir",       default="",     help="Directory of LAB reference genomes")
ap.add_argument("--threads",       default=8,      type=int)
ap.add_argument("--ani_threshold", default=95.0,   type=float, help="ANI threshold for same species")
ap.add_argument("--download_refs", action="store_true", help="Auto-download LAB reference genomes")
args = ap.parse_args()

os.makedirs(args.outdir, exist_ok=True)

# ── Check FastANI is available ────────────────────────────────────────────────
def check_fastani():
    try:
        result = subprocess.run(["fastANI", "--version"],
                                capture_output=True, text=True, timeout=10)
        print(f"[species_verify] FastANI: {result.stdout.strip() or result.stderr.strip()}")
        return True
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False

HAS_FASTANI = check_fastani()
if not HAS_FASTANI:
    print("[species_verify] WARNING: FastANI not found.")
    print("[species_verify] Install with: conda install -c bioconda fastani")
    print("[species_verify] Skipping genome-based verification — only NCBI names will be used.")

# ── LAB reference species (type strains) ─────────────────────────────────────
# Top LAB species to download as references
LAB_REFERENCES = {
    "Lactobacillus_acidophilus":           "GCF_000011985.1",
    "Lactobacillus_delbrueckii":           "GCF_000056065.1",
    "Lactobacillus_helveticus":            "GCF_000026205.1",
    "Limosilactobacillus_reuteri":         "GCF_000016825.1",
    "Ligilactobacillus_salivarius":        "GCF_000026065.1",
    "Lacticaseibacillus_casei":            "GCF_000309565.1",
    "Lacticaseibacillus_rhamnosus":        "GCF_000026505.1",
    "Lactiplantibacillus_plantarum":       "GCF_000203855.1",
    "Levilactobacillus_brevis":            "GCF_000006165.2",
    "Pediococcus_acidilactici":            "GCF_000026185.1",
    "Pediococcus_pentosaceus":             "GCF_000014205.1",
    "Enterococcus_faecalis":               "GCF_000007785.1",
    "Enterococcus_faecium":                "GCF_000174395.1",
    "Lactococcus_lactis":                  "GCF_000006865.1",
    "Streptococcus_thermophilus":          "GCF_000011825.1",
    "Leuconostoc_mesenteroides":           "GCF_000159175.1",
    "Oenococcus_oeni":                     "GCF_000010045.1",
    "Apilactobacillus_kunkeei":            "GCF_000741865.1",
    "Fructilactobacillus_sanfranciscensis":"GCF_000621405.1",
    "Companilactobacillus_crispatus":      "GCF_000165905.1",
    "Limosilactobac_fermentum":            "GCF_000010045.1",
    "Lactobacillus_gasseri":               "GCF_000014425.1",
    "Lactobacillus_johnsonii":             "GCF_000010045.1",
    "Lactobacillus_crispatus":             "GCF_000165905.1",
}

def download_reference_genomes(ref_dir: str):
    """Download LAB type strain genomes from NCBI."""
    os.makedirs(ref_dir, exist_ok=True)
    downloaded = 0

    print(f"[species_verify] Downloading {len(LAB_REFERENCES)} LAB reference genomes...")

    for species, accession in LAB_REFERENCES.items():
        out_file = Path(ref_dir) / f"{species}.fna"
        if out_file.exists() and out_file.stat().st_size > 1000:
            print(f"[species_verify]   Skip (exists): {species}")
            continue

        # Use NCBI datasets to download
        try:
            url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/download?include_annotation_type=GENOME_FASTA"
            print(f"[species_verify]   Downloading: {species} ({accession})")
            req  = Request(url, headers={"Accept": "application/zip"})
            resp = urlopen(req, timeout=60)

            # Save zip and extract
            zip_path = Path(ref_dir) / f"{species}.zip"
            with open(zip_path, "wb") as f:
                f.write(resp.read())

            # Extract FASTA
            import zipfile
            with zipfile.ZipFile(zip_path) as z:
                for name in z.namelist():
                    if name.endswith(".fna") and "genomic" in name:
                        with z.open(name) as src, open(out_file, "wb") as dst:
                            dst.write(src.read())
                        break

            zip_path.unlink()
            downloaded += 1
            time.sleep(0.5)

        except Exception as e:
            print(f"[species_verify]   FAILED: {species}: {e}")

    print(f"[species_verify] Downloaded {downloaded} reference genomes to {ref_dir}")
    return ref_dir

# ── Get or set ref_dir ────────────────────────────────────────────────────────
ref_dir = args.ref_dir
if not ref_dir:
    ref_dir = os.path.join(args.outdir, "lab_references")

if args.download_refs or not os.path.exists(ref_dir) or \
   len(list(Path(ref_dir).glob("*.fna"))) < 5:
    if args.download_refs:
        download_reference_genomes(ref_dir)
    else:
        print(f"[species_verify] No reference genomes found at: {ref_dir}")
        print(f"[species_verify] Use --download_refs to auto-download, or --ref_dir to specify path")

ref_genomes = list(Path(ref_dir).glob("*.fna")) if os.path.exists(ref_dir) else []
print(f"[species_verify] Reference genomes available: {len(ref_genomes)}")

# ── Load metadata ─────────────────────────────────────────────────────────────
meta = pd.read_csv(args.metadata, sep="\t", dtype=str)
print(f"[species_verify] Loaded metadata: {len(meta)} genomes")

# ── Run FastANI ───────────────────────────────────────────────────────────────
EXTS = [".fna", ".fa", ".fasta"]
query_genomes = sorted([
    p for p in Path(args.indir).iterdir()
    if p.suffix in EXTS
])

ani_results = {}  # genome_stem -> {ani_species, ani_identity}

if HAS_FASTANI and ref_genomes and query_genomes:
    # Write query and reference lists
    query_list = Path(args.outdir) / "fastani_query.txt"
    ref_list   = Path(args.outdir) / "fastani_ref.txt"
    ani_out    = Path(args.outdir) / "fastani_output.txt"

    with open(query_list, "w") as f:
        for g in query_genomes:
            f.write(str(g) + "\n")

    with open(ref_list, "w") as f:
        for r in ref_genomes:
            f.write(str(r) + "\n")

    print(f"[species_verify] Running FastANI: {len(query_genomes)} queries × {len(ref_genomes)} refs...")

    cmd = [
        "fastANI",
        "--ql",     str(query_list),
        "--rl",     str(ref_list),
        "-o",       str(ani_out),
        "--threads", str(args.threads),
        "--minFraction", "0.2",
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        if result.returncode != 0:
            print(f"[species_verify] FastANI stderr: {result.stderr[:500]}")
        else:
            print(f"[species_verify] FastANI completed.")
    except subprocess.TimeoutExpired:
        print("[species_verify] FastANI timed out after 1 hour.")
    except Exception as e:
        print(f"[species_verify] FastANI error: {e}")

    # Parse FastANI output: query ref ANI fragments total_frags
    if ani_out.exists():
        best = {}  # query_stem -> (ani, ref_species)
        with open(ani_out) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 3:
                    continue
                query_path = parts[0]
                ref_path   = parts[1]
                ani_val    = float(parts[2])

                query_stem = Path(query_path).stem
                ref_stem   = Path(ref_path).stem

                if query_stem not in best or ani_val > best[query_stem][0]:
                    best[query_stem] = (ani_val, ref_stem)

        for stem, (ani_val, ref_species) in best.items():
            # Clean up species name from filename
            species_name = ref_species.replace("_", " ").strip()
            ani_results[stem] = {
                "ani_species":  species_name,
                "ani_identity": round(ani_val, 2),
            }

        print(f"[species_verify] ANI results parsed: {len(ani_results)} genomes")

# ── Add verification columns to metadata ─────────────────────────────────────
def determine_status(row, ani_threshold):
    ncbi_sp = str(row.get("full_species", "")).strip()
    ani_sp  = str(row.get("ani_species",  "")).strip()
    ani_id  = row.get("ani_identity", 0)

    if not ani_sp or ani_sp == "nan":
        if not HAS_FASTANI:
            return "NO_FASTANI"
        elif not ref_genomes:
            return "NO_REFS"
        else:
            return "UNRESOLVED"

    try:
        ani_id = float(ani_id)
    except (ValueError, TypeError):
        ani_id = 0

    if ani_id < args.ani_threshold:
        return "LOW_ANI"

    if not ncbi_sp or ncbi_sp in ("Unknown sp.", "nan", ""):
        return "NO_NCBI"

    # Compare genus + species (first 2 words)
    ncbi_parts = ncbi_sp.lower().split()[:2]
    ani_parts  = ani_sp.lower().split()[:2]

    if ncbi_parts == ani_parts:
        return "CONFIRMED"
    elif ncbi_parts[0] == ani_parts[0]:  # Same genus, different species
        return "GENUS_MATCH"
    else:
        return "MISMATCH"

# Apply results
meta["ani_species"]  = meta["genome"].map(lambda g: ani_results.get(g, {}).get("ani_species",  ""))
meta["ani_identity"] = meta["genome"].map(lambda g: ani_results.get(g, {}).get("ani_identity", ""))
meta["species_status"] = meta.apply(lambda r: determine_status(r, args.ani_threshold), axis=1)

# Use ANI species when NCBI is unknown
meta["resolved_species"] = meta.apply(
    lambda r: r["ani_species"] if (
        str(r.get("full_species","")) in ("Unknown sp.", "", "nan")
        and str(r.get("ani_species","")) not in ("", "nan")
    ) else r.get("full_species", "Unknown sp."),
    axis=1
)

# Save updated metadata
out_path = os.path.join(args.outdir, "genome_metadata.tsv")
meta.to_csv(out_path, sep="\t", index=False)

# Summary
print(f"\n[species_verify] Species status summary:")
print(meta["species_status"].value_counts().to_string())

if "MISMATCH" in meta["species_status"].values:
    print(f"\n[species_verify] ⚠ MISMATCHES DETECTED:")
    mismatches = meta[meta["species_status"] == "MISMATCH"][
        ["genome","full_species","ani_species","ani_identity"]
    ]
    print(mismatches.to_string(index=False))
    mismatches.to_csv(f"{args.outdir}/species_mismatches.tsv", sep="\t", index=False)

if "GENUS_MATCH" in meta["species_status"].values:
    print(f"\n[species_verify] ℹ GENUS MATCH (possible subspecies/strain):")
    genus_match = meta[meta["species_status"] == "GENUS_MATCH"][
        ["genome","full_species","ani_species","ani_identity"]
    ]
    print(genus_match.to_string(index=False))

print(f"\n[species_verify] Updated metadata saved: {out_path}")
