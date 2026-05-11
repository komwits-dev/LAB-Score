#!/usr/bin/env python3
"""
07_panel_scan.py
================
Scans Prokka GFF files for all 18 functional probiotic panels.
Uses curated gene tokens matching your original scan_gff_*.py scripts.

Output:
  {outdir}/all_panels.tsv  — one row per genome, one column per panel
"""

import argparse, os, re
import pandas as pd
from pathlib import Path
from collections import defaultdict
from urllib.parse import unquote

ap = argparse.ArgumentParser()
ap.add_argument("--prokka_dir", required=True)
ap.add_argument("--outdir",     required=True)
args = ap.parse_args()
os.makedirs(args.outdir, exist_ok=True)

# ═══════════════════════════════════════════════════════════════════════════════
# Curated panel gene tokens
# ═══════════════════════════════════════════════════════════════════════════════
PANELS = {

"acid_energy": [
    r"atpA", r"atpB", r"atpC", r"atpD", r"atpE", r"atpF", r"atpG", r"atpH", r"atpI",
    r"F.type.ATP.synthase", r"ATP.synthase", r"H\+.transporting",
    r"proton.translocating.ATP",
    r"gadA", r"gadB", r"gadC",
    r"glutamate.decarboxylase", r"glutamate.GABA.antiporter",
    r"acid.resistance", r"acid.tolerance", r"proton.motive",
],

"adhesion_biofilm": [
    r"biofilm", r"auto.aggregation", r"aggregation.promoting",
    r"mub", r"MapA", r"Msa",
    r"mucus.binding", r"mucin.binding",
    r"collagen.binding", r"fibronectin.binding",
],

"adhesion_surface": [
    r"S.layer", r"surface.layer", r"slpA", r"slpB",
    r"sortase", r"LPXTG", r"cell.wall.anchor",
    r"pilus", r"pilin", r"fimbriae", r"fimbrial",
    r"surface.protein", r"surface.exposed",
],

"alkalinestress": [
    r"alkali.stress", r"alkali.tolerance",
    r"Na.*H.*antiporter", r"sodium.proton.antiporter",
    r"nhaA", r"nhaB", r"nhaC",
    r"mrpA", r"mrpB", r"mrpC",
    r"cation.proton.antiporter", r"pH.homeostasis",
],

"antipath_qs": [
    r"quorum.sensing", r"quorum-sensing",
    r"luxS", r"S-ribosylhomocysteine",
    r"autoinducer", r"auto.inducer",
    r"agrA", r"agrB", r"agrC", r"agrD",
    r"comX", r"comS", r"pheromone",
],

"bacteriocins": [
    r"bacteriocin", r"nisin", r"lantibiotic",
    r"sakacin", r"plantaricin", r"gassericin",
    r"acidocin", r"enterocin", r"pediocin", r"lactocin",
    r"antimicrobial.peptide", r"bacteriocin.immunity",
],

"bileresistance": [
    r"bile.salt.hydrolase", r"cholylglycine.hydrolase",
    r"bsh", r"bile.acid", r"bile.resistance",
    r"bile.efflux", r"bile.salt.transporter",
],

"carbohydrate": [
    r"phosphotransferase.system", r"PTS.system",
    r"sugar.phosphotransferase", r"EIIA", r"EIIB", r"EIIC",
    r"lacA", r"lacB", r"lacF", r"lacG", r"lacZ", r"lacS",
    r"amylase", r"glucoamylase",
    r"alpha.glucosidase", r"beta.glucosidase",
    r"fructosidase", r"xylosidase",
    r"arabinofuranosidase", r"galactosidase",
    r"sucrase", r"maltase", r"trehalase",
],

"cellenvelope_eps": [
    r"exopolysaccharide", r"EPS.biosynthesis",
    r"epsA", r"epsB", r"epsC", r"epsD", r"epsE",
    r"capsular.polysaccharide", r"glycosyltransferase",
    r"UDP.galactose", r"UDP.glucose",
    r"polysaccharide.biosynthesis",
    r"teichoic.acid", r"lipoteichoic.acid",
],

"coldstress": [
    r"cold.shock", r"cold.inducible", r"cold.adaptation",
    r"cspA", r"cspB", r"cspC", r"cspD", r"cspE",
    r"RNA.helicase", r"DEAD.box.helicase",
    r"cryoprotect",
],

"defense_crispr": [
    r"CRISPR", r"CRISPR.associated",
    r"cas1", r"cas2", r"cas3", r"cas9",
    r"type.I.restriction", r"type.II.restriction",
    r"restriction.modification", r"methyltransferase",
    r"phage.resistance",
],

"gaba": [
    r"glutamate.decarboxylase",
    r"gadA", r"gadB",
    r"GABA.transaminase", r"4.aminobutyrate",
    r"gamma.aminobutyric",
    r"succinate.semialdehyde", r"gabT", r"gabD",
    r"GABA.permease",
],

"gutpersistence": [
    r"mucus.binding", r"mucin.binding", r"mucin.degradation",
    r"gut.colonization", r"intestinal.persistence",
    r"intestinal.epithelium", r"epithelial.adhesion",
    r"oxygen.tolerance", r"NADH.oxidase",
    r"nox", r"noxE",
],

"heatstress": [
    r"heat.shock", r"heat.inducible",
    r"groEL", r"groES", r"dnaK", r"dnaJ", r"grpE",
    r"clpB", r"clpP", r"clpX", r"htpG", r"htpX",
    r"Lon.protease", r"molecular.chaperone",
    r"small.heat.shock", r"hsp",
],

"immunomodulation": [
    r"immunomodulat", r"immunostimulat",
    r"lipoteichoic.acid",
    r"peptidoglycan.hydrolase", r"N.acetylmuramidase",
    r"D.alanylation", r"D.alanine",
    r"flagellin", r"fliC",
    r"muropeptide",
],

"metabolism": [
    r"lactate.dehydrogenase", r"ldhA", r"ldhD",
    r"acetate.kinase", r"ackA",
    r"phosphotransacetylase", r"pta",
    r"pyruvate.kinase", r"pyk",
    r"acetaldehyde.dehydrogenase",
    r"alcohol.dehydrogenase",
    r"diacetyl", r"acetoin",
    r"mixed.acid.fermentation",
],

"osmoticstress": [
    r"osmotic.stress", r"osmostress", r"osmolyte",
    r"betaine.transporter", r"glycine.betaine",
    r"opuA", r"opuB", r"opuC",
    r"carnitine.transporter",
    r"ectoine", r"hydroxyectoine",
    r"compatible.solute", r"proline.transporter",
],

"vitamins": [
    r"riboflavin", r"riboflavin.synthase",
    r"ribB", r"ribC", r"ribD", r"ribE",
    r"folate", r"folate.biosynthesis",
    r"folP", r"folB", r"folC",
    r"cobalamin", r"vitamin.B12",
    r"thiamine", r"pyridoxal", r"pyridoxine",
    r"pantothenate", r"biotin", r"menaquinone",
],

}

# Compile all panels
COMPILED_PANELS = {
    name: re.compile("|".join(patterns), re.IGNORECASE)
    for name, patterns in PANELS.items()
}

# ── Per-genome GFF scanning ───────────────────────────────────────────────────
def parse_gff_attributes(attr_string):
    """Parse GFF3 attributes, handling URL-encoding and semicolon separation."""
    attrs = {}
    for field in attr_string.rstrip('\n').split(';'):
        field = field.strip()
        if '=' in field:
            key, _, val = field.partition('=')
            # URL-decode (Prokka sometimes encodes spaces as %20)
            attrs[key.lower()] = unquote(val).strip()
    return attrs

def scan_gff(gff_path: str) -> dict:
    counts = defaultdict(int)
    try:
        with open(gff_path) as fh:
            for line in fh:
                if line.startswith(("#", ">")) or not line.strip():
                    continue
                parts = line.split("\t")
                if len(parts) < 9 or parts[2] in ("region", "source"):
                    continue
                attrs   = parse_gff_attributes(parts[8])
                product = attrs.get("product", "")
                gene    = attrs.get("gene", "")
                name    = attrs.get("name", "")
                text    = f"{product} {gene} {name}"
                for panel, pat in COMPILED_PANELS.items():
                    if pat.search(text):
                        counts[panel] += 1
    except Exception as e:
        print(f"[panels] Error scanning {gff_path}: {e}")
    return counts

# ── Process all genomes ───────────────────────────────────────────────────────
prokka_dir = Path(args.prokka_dir)
samples    = sorted([d.name for d in prokka_dir.iterdir() if d.is_dir()])
print(f"[panels] Scanning {len(samples)} genomes across {len(PANELS)} panels")

rows = []
for sample in samples:
    # Find GFF — prefix may differ from folder name
    gff_files = list(Path(prokka_dir, sample).glob("*.gff"))
    gff       = gff_files[0] if gff_files else None

    if gff:
        counts = scan_gff(str(gff))
    else:
        print(f"[panels] WARNING: No GFF found for {sample}")
        counts = {}

    row = {"genome": sample}
    for panel in PANELS:
        col = f"{panel}.{panel}"
        row[col] = counts.get(panel, 0)
    row["TOTAL_all_markers"] = sum(counts.get(p, 0) for p in PANELS)
    rows.append(row)

df = pd.DataFrame(rows)
print(f"[panels] Output shape: {df.shape}")

print("[panels] Panel marker count summary:")
for panel in PANELS:
    col = f"{panel}.{panel}"
    if col in df.columns:
        print(f"  {panel:<25} mean={df[col].mean():.1f}  max={df[col].max()}")

out_path = os.path.join(args.outdir, "all_panels.tsv")
df.to_csv(out_path, sep="\t", index=False)
print(f"[panels] Saved: {out_path}")
