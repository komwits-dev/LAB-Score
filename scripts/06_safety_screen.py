#!/usr/bin/env python3
"""
06_safety_screen.py — context-aware LAB safety screening

Conservative interpretation rules:
1. Prokka annotation names alone are not treated as proof of haemolytic toxicity.
   - tlyA, hlyB, ply/pneumolysin, haemolysin-III, putative haemolysin, and
     isolated cytolysin-like annotations are reported for review only.
   - A GFF-derived cautionary cytolysin signal requires an operon-like cluster:
     at least one cytolysin structural gene plus at least two cyl accessory genes
     on the same contig within 30 kb.

2. Biogenic-amine caution requires pathway context rather than an isolated
   decarboxylase annotation.
   - histamine: hdcA + hdcP
   - tyramine: tyrDC/tdcA + tyrP
   - putrescine (ODC): odc + potE
   - cadaverine: cadA/ldc + cadB
   - agmatine deiminase route: aguA + aguB + aguC
   All partners must occur on the same contig within the configured window.

3. gadA/gadB/gadC and glutamate decarboxylase/GABA transport genes are recorded
   as GABA-related and are never penalized as adverse biogenic-amine markers.

4. AMRFinderPlus evidence is stratified.
   - Core AMR/POINT hits supported by EXACT/ALLELE/BLAST evidence are critical.
   - Lower-confidence core hits are critical only for specific clinically
     interpretable determinants (for example bla, tet, erm, aminoglycoside-
     modifying enzymes, dfr, poxtA/optrA/cfr, fex, cat, lnu, mef, or msr).
   - Broad HMM-family calls such as lsa, isolated vanR/vat-family calls, and
     plus-scope AMR records are cautionary.
   - Plus-scope STRESS records are informational only.
   - VIRULENCE records remain critical alerts.

Legacy downstream columns are retained:
- biogenic_amines = number of context-supported adverse BA pathways
- hemolysin_markers = number of context-supported cytolysin clusters
"""

import argparse
import os
import re
from pathlib import Path
from urllib.parse import unquote

import pandas as pd

BA_WINDOW = 20_000
AGU_WINDOW = 30_000
CYL_WINDOW = 30_000

AMR_HIGH_METHODS = {
    "EXACT", "EXACTP", "ALLELE", "ALLELEP", "BLAST", "BLASTP"
}
AMR_LOWER_METHODS = {"HMM", "PARTIAL", "PARTIALP"}

# A lower-confidence match is promoted to Critical only when its symbol
# represents a specific, clinically interpretable acquired-resistance
# determinant. Broad family/regulatory calls such as lsa, vat, vanR and vanS
# remain Caution unless supported by high-confidence evidence.
LOWER_CONFIDENCE_CRITICAL_PATTERNS = [
    re.compile(r"^bla", re.IGNORECASE),
    re.compile(r"^tet(?:\(|[A-Za-z0-9])", re.IGNORECASE),
    re.compile(r"^erm(?:\(|[A-Za-z0-9])", re.IGNORECASE),
    re.compile(r"^(aac|aad|ant|aph)", re.IGNORECASE),
    re.compile(r"^dfr", re.IGNORECASE),
    re.compile(r"^(poxtA|optrA|cfr)", re.IGNORECASE),
    re.compile(r"^fex[A-Za-z0-9]*$", re.IGNORECASE),
    re.compile(r"^cat[A-Za-z0-9]*$", re.IGNORECASE),
    re.compile(r"^lnu(?:\(|[A-Za-z0-9])", re.IGNORECASE),
    re.compile(r"^mef(?:\(|[A-Za-z0-9])", re.IGNORECASE),
    re.compile(r"^msr(?:\(|[A-Za-z0-9])", re.IGNORECASE),
]

ap = argparse.ArgumentParser()
ap.add_argument("--prokka_dir", required=True)
ap.add_argument("--amrfinder", required=True)
ap.add_argument("--resfinder", required=True)
ap.add_argument("--outdir", required=True)
args = ap.parse_args()
os.makedirs(args.outdir, exist_ok=True)


def parse_attrs(attr_string):
    attrs = {}
    for field in attr_string.rstrip("\n").split(";"):
        field = field.strip()
        if "=" not in field:
            continue
        key, value = field.split("=", 1)
        attrs[key.lower()] = unquote(value).strip()
    return attrs


def norm_gene(value):
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def short_detail(feature):
    gene = feature["gene"] or feature["name"]
    product = feature["product"]
    label = f"{gene}:{product}" if gene and product else gene or product
    return f"{feature['contig']}:{feature['start']}-{feature['end']}:{label}"


def dedupe_join(values):
    seen = set()
    out = []
    for value in values:
        value = str(value).strip()
        if value and value not in seen:
            seen.add(value)
            out.append(value)
    return "|".join(out)


def feature_distance(a, b):
    if a["contig"] != b["contig"]:
        return None
    if a["end"] < b["start"]:
        return b["start"] - a["end"]
    if b["end"] < a["start"]:
        return a["start"] - b["end"]
    return 0


def near(a, b, window):
    distance = feature_distance(a, b)
    return distance is not None and distance <= window


def gene_is(feature, names):
    observed = {
        norm_gene(feature["gene"]),
        norm_gene(feature["name"]),
    }
    wanted = {norm_gene(x) for x in names}
    return bool(observed & wanted)


def product_has(feature, pattern):
    return bool(re.search(pattern, feature["product"], re.IGNORECASE))


def load_gff_features(gff_path):
    features = []
    try:
        with open(gff_path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith(("#", ">")) or not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9 or parts[2] in {"region", "source"}:
                    continue
                try:
                    start = int(parts[3])
                    end = int(parts[4])
                except ValueError:
                    continue
                attrs = parse_attrs(parts[8])
                features.append(
                    {
                        "contig": parts[0],
                        "start": start,
                        "end": end,
                        "type": parts[2],
                        "gene": attrs.get("gene", ""),
                        "name": attrs.get("name", ""),
                        "product": attrs.get("product", ""),
                    }
                )
    except Exception as exc:
        print(f"[safety] GFF error {gff_path}: {exc}")
    return features


def find_supported_pairs(left, right, pathway, window):
    hits = []
    for a in left:
        partners = [b for b in right if near(a, b, window)]
        if partners:
            best = min(partners, key=lambda b: feature_distance(a, b))
            hits.append(
                f"{pathway}[{short_detail(a)} + {short_detail(best)}]"
            )
    return hits


def find_agmatine_clusters(agu_a, agu_b, agu_c):
    hits = []
    for a in agu_a:
        bs = [b for b in agu_b if near(a, b, AGU_WINDOW)]
        cs = [c for c in agu_c if near(a, c, AGU_WINDOW)]
        if not bs or not cs:
            continue
        for b in bs:
            for c in cs:
                if b["contig"] != c["contig"]:
                    continue
                cluster_start = min(a["start"], b["start"], c["start"])
                cluster_end = max(a["end"], b["end"], c["end"])
                if cluster_end - cluster_start <= AGU_WINDOW:
                    hits.append(
                        "agmatine_deiminase["
                        + " + ".join(
                            [short_detail(a), short_detail(b), short_detail(c)]
                        )
                        + "]"
                    )
                    break
            if hits:
                break
    return hits


def find_cytolysin_clusters(structural, accessory):
    hits = []
    for structural_gene in structural:
        nearby = [
            f
            for f in accessory
            if near(structural_gene, f, CYL_WINDOW)
        ]
        distinct_accessory = {}
        for f in nearby:
            key = norm_gene(f["gene"] or f["name"])
            if key:
                distinct_accessory[key] = f

        if len(distinct_accessory) >= 2:
            selected = list(distinct_accessory.values())
            hits.append(
                "cytolysin_operon["
                + " + ".join(
                    [short_detail(structural_gene)]
                    + [short_detail(f) for f in selected]
                )
                + "]"
            )
    return hits


def parse_gff_safety(gff_path):
    features = load_gff_features(gff_path)

    gaba = []
    broad_decarboxylase = []
    hly_like = []

    hdc_a, hdc_p = [], []
    tyr_dc, tyr_p = [], []
    odc, pot_e = [], []
    ldc, cad_b = [], []
    agu_a, agu_b, agu_c = [], [], []

    cyl_structural = []
    cyl_accessory = []

    for f in features:
        gene = norm_gene(f["gene"] or f["name"])
        product = f["product"]
        text = f"{f['gene']} {f['name']} {product}"

        # GABA system: report, never penalize.
        if (
            gene in {"gada", "gadb", "gadc"}
            or re.search(r"\bglutamate[\s_-]+decarboxylase\b", text, re.I)
            or re.search(r"\bglutamate[/\s_-]+GABA[\s_-]+antiporter\b", text, re.I)
            or re.search(r"\bgamma[\s_-]+aminobutyr", text, re.I)
        ):
            gaba.append(short_detail(f))
            continue

        # Specific BA genes and partners.
        if gene in {"hdca"} or product_has(f, r"\bhistidine[\s_-]+decarboxylase\b"):
            hdc_a.append(f)
        elif gene in {"hdcp"} or product_has(f, r"histidine[/\s_-]+histamine.*antiporter"):
            hdc_p.append(f)

        if gene in {"tyrdc", "tdca"} or product_has(f, r"\btyrosine[\s_-]+decarboxylase\b"):
            tyr_dc.append(f)
        elif gene in {"tyrp"} or product_has(f, r"tyrosine[/\s_-]+tyramine.*antiporter"):
            tyr_p.append(f)

        if gene in {"odci", "odca", "odc"} or product_has(f, r"\bornithine[\s_-]+decarboxylase\b"):
            odc.append(f)
        elif gene in {"pote"} or product_has(f, r"ornithine[/\s_-]+putrescine.*antiporter"):
            pot_e.append(f)

        if gene in {"cada", "ldcc", "ldci", "ldc"} or product_has(f, r"\blysine[\s_-]+decarboxylase\b"):
            ldc.append(f)
        elif gene in {"cadb"} or product_has(f, r"lysine[/\s_-]+cadaverine.*antiporter"):
            cad_b.append(f)

        if gene in {"agua"} or product_has(f, r"\bagmatine[\s_-]+deiminase\b"):
            agu_a.append(f)
        elif gene in {"agub"} or product_has(f, r"\bputrescine[\s_-]+carbamoyltransferase\b"):
            agu_b.append(f)
        elif gene in {"aguc"} or product_has(f, r"\bcarbamate[\s_-]+kinase\b"):
            agu_c.append(f)

        if re.search(
            r"\bamino[\s_-]+acid[\s_-]+decarboxylase\b|"
            r"\bpyridoxal[\s_-]+phosphate[\s_-]+dependent[\s_-]+decarboxylase\b|"
            r"\bPLP[\s_-]+dependent[\s_-]+decarboxylase\b",
            text,
            re.I,
        ):
            broad_decarboxylase.append(short_detail(f))

        # Context-supported Enterococcus-type cytolysin only.
        if (
            gene in {"cylll", "cylls"}
            or product_has(f, r"\bcytolysin\b.*\b(large|small)\b.*\bsubunit\b")
        ):
            cyl_structural.append(f)
        elif gene in {"cylm", "cylb", "cyla", "cyli"}:
            cyl_accessory.append(f)

        # All annotation-only haemolysin/toxin names are review-only unless
        # they participate in the supported cyl cluster above.
        if re.search(
            r"\bha?emolysin\b|"
            r"\bcytolysin\b|"
            r"\bpneumolysin\b|"
            r"\blisteriolysin\b|"
            r"\bstreptolysin\b|"
            r"\bperfringolysin\b|"
            r"\btlyA\b|"
            r"\bhlyA\b|"
            r"\bhlyB\b|"
            r"\bply\b",
            text,
            re.I,
        ):
            hly_like.append(short_detail(f))

    ba_pathways = []
    ba_pathways += find_supported_pairs(hdc_a, hdc_p, "histamine", BA_WINDOW)
    ba_pathways += find_supported_pairs(tyr_dc, tyr_p, "tyramine", BA_WINDOW)
    ba_pathways += find_supported_pairs(odc, pot_e, "putrescine_ODC", BA_WINDOW)
    ba_pathways += find_supported_pairs(ldc, cad_b, "cadaverine", BA_WINDOW)
    ba_pathways += find_agmatine_clusters(agu_a, agu_b, agu_c)

    cytolysin_clusters = find_cytolysin_clusters(cyl_structural, cyl_accessory)

    # Remove confirmed cluster members from the review-only list where possible.
    ba_pathway_types = sorted(
        {
            hit.split("[", 1)[0].strip()
            for hit in ba_pathways
            if str(hit).strip()
        }
    )

    return {
        # Number of distinct supported pathway types, not duplicate gene copies.
        "biogenic_amines": len(ba_pathway_types),
        "biogenic_amine_pathway_types": ";".join(ba_pathway_types),
        "biogenic_amine_cluster_count": len(ba_pathways),
        "biogenic_amine_details": dedupe_join(ba_pathways),
        "nonspecific_decarboxylase_annotations": len(broad_decarboxylase),
        "nonspecific_decarboxylase_details": dedupe_join(broad_decarboxylase),
        "gaba_related_annotations": len(gaba),
        "gaba_related_details": dedupe_join(gaba),
        "hemolysin_markers": len(cytolysin_clusters),
        "hemolysin_marker_details": dedupe_join(cytolysin_clusters),
        "hemolysin_like_annotations": len(hly_like),
        "hemolysin_like_details": dedupe_join(hly_like),
    }


def lower_confidence_is_specific(symbol, name):
    target = f"{symbol} {name}".strip()
    return any(pattern.search(target) for pattern in LOWER_CONFIDENCE_CRITICAL_PATTERNS)


def parse_amrfinder(tsv_path):
    """
    Parse AMRFinderPlus v4 output and stratify evidence.

    Returned categories:
    - critical AMR: core + high-confidence method, or a specific clinically
      interpretable determinant detected by a lower-confidence method
    - caution AMR: broad/uncertain core HMM or partial match, or plus-scope AMR
    - stress: informational only
    - virulence: critical
    """
    result = {
        "amr_hits": 0,
        "critical_amr_hits": 0,
        "caution_amr_hits": 0,
        "virulence_hits": 0,
        "stress_hits": 0,
        "amr_symbols": "",
        "critical_amr_symbols": "",
        "caution_amr_symbols": "",
        "virulence_symbols": "",
        "stress_symbols": "",
        "amr_details": "",
        "critical_amr_details": "",
        "caution_amr_details": "",
        "virulence_details": "",
        "stress_details": "",
    }

    try:
        df = pd.read_csv(tsv_path, sep="\t", dtype=str)
    except Exception as exc:
        print(f"[safety] AMRFinder parse warning {tsv_path}: {exc}")
        return result

    if df.empty:
        return result

    def normalize_header(value):
        value = str(value).strip().lstrip("#").lower()
        return re.sub(r"[^a-z0-9]+", "", value)

    normalized = {normalize_header(col): col for col in df.columns}

    type_col = normalized.get("type") or normalized.get("elementtype")
    symbol_col = (
        normalized.get("elementsymbol")
        or normalized.get("genesymbol")
        or normalized.get("symbol")
    )
    name_col = (
        normalized.get("elementname")
        or normalized.get("sequencename")
        or normalized.get("name")
    )
    scope_col = normalized.get("scope")
    method_col = normalized.get("method")
    class_col = normalized.get("class")
    subclass_col = normalized.get("subclass")

    if type_col is None:
        print(
            f"[safety] AMRFinder warning: no Type/Element type column in {tsv_path}"
        )
        return result

    buckets = {
        "critical_amr": [],
        "caution_amr": [],
        "virulence": [],
        "stress": [],
    }

    for _, row in df.iterrows():
        raw_type = str(row.get(type_col, "")).strip().upper()
        tokens = {
            token.strip()
            for token in re.split(r"[,;/|]+", raw_type)
            if token.strip()
        }

        symbol = (
            ""
            if symbol_col is None or pd.isna(row.get(symbol_col))
            else str(row.get(symbol_col, "")).strip()
        )
        name = (
            ""
            if name_col is None or pd.isna(row.get(name_col))
            else str(row.get(name_col, "")).strip()
        )
        scope = (
            ""
            if scope_col is None or pd.isna(row.get(scope_col))
            else str(row.get(scope_col, "")).strip().lower()
        )
        method = (
            ""
            if method_col is None or pd.isna(row.get(method_col))
            else str(row.get(method_col, "")).strip().upper()
        )
        drug_class = (
            ""
            if class_col is None or pd.isna(row.get(class_col))
            else str(row.get(class_col, "")).strip()
        )
        subclass = (
            ""
            if subclass_col is None or pd.isna(row.get(subclass_col))
            else str(row.get(subclass_col, "")).strip()
        )

        detail_parts = []
        if symbol:
            detail_parts.append(symbol)
        if name:
            detail_parts.append(name)
        if scope:
            detail_parts.append(f"scope={scope}")
        if method:
            detail_parts.append(f"method={method}")
        if drug_class:
            detail_parts.append(f"class={drug_class}")
        if subclass:
            detail_parts.append(f"subclass={subclass}")

        item = {
            "symbol": symbol,
            "detail": "; ".join(detail_parts),
        }

        if "VIRULENCE" in tokens:
            buckets["virulence"].append(item)
            continue

        if "STRESS" in tokens and not (tokens & {"AMR", "POINT"}):
            buckets["stress"].append(item)
            continue

        if not (tokens & {"AMR", "POINT"}):
            continue

        # Plus-scope AMR records are supplementary and remain cautionary.
        if scope == "plus":
            buckets["caution_amr"].append(item)
        # Core records with strong sequence/allele evidence are critical.
        elif scope == "core" and method in AMR_HIGH_METHODS:
            buckets["critical_amr"].append(item)
        # Specific clinically interpretable lower-confidence determinants
        # remain critical; broad family/regulatory matches remain cautionary.
        elif scope == "core" and method in AMR_LOWER_METHODS:
            if lower_confidence_is_specific(symbol, name):
                buckets["critical_amr"].append(item)
            else:
                buckets["caution_amr"].append(item)
        else:
            # Unknown scope/method: retain conservatively as caution.
            buckets["caution_amr"].append(item)

    critical_rows = buckets["critical_amr"]
    caution_rows = buckets["caution_amr"]
    virulence_rows = buckets["virulence"]
    stress_rows = buckets["stress"]
    all_amr_rows = critical_rows + caution_rows

    result["amr_hits"] = len(all_amr_rows)
    result["critical_amr_hits"] = len(critical_rows)
    result["caution_amr_hits"] = len(caution_rows)
    result["virulence_hits"] = len(virulence_rows)
    result["stress_hits"] = len(stress_rows)

    result["amr_symbols"] = dedupe_join(item["symbol"] for item in all_amr_rows)
    result["critical_amr_symbols"] = dedupe_join(
        item["symbol"] for item in critical_rows
    )
    result["caution_amr_symbols"] = dedupe_join(
        item["symbol"] for item in caution_rows
    )
    result["virulence_symbols"] = dedupe_join(
        item["symbol"] for item in virulence_rows
    )
    result["stress_symbols"] = dedupe_join(
        item["symbol"] for item in stress_rows
    )

    result["amr_details"] = dedupe_join(item["detail"] for item in all_amr_rows)
    result["critical_amr_details"] = dedupe_join(
        item["detail"] for item in critical_rows
    )
    result["caution_amr_details"] = dedupe_join(
        item["detail"] for item in caution_rows
    )
    result["virulence_details"] = dedupe_join(
        item["detail"] for item in virulence_rows
    )
    result["stress_details"] = dedupe_join(
        item["detail"] for item in stress_rows
    )

    return result


def parse_resfinder_summary(tsv_path):
    result = {}
    try:
        df = pd.read_csv(tsv_path, sep="\t", dtype=str)
        for _, row in df.iterrows():
            genome = str(row.get("genome", "")).strip()
            raw_hits = row.get("resfinder_hits", 0)
            try:
                hits = 0 if pd.isna(raw_hits) else int(float(str(raw_hits)))
            except Exception:
                hits = 0
            result[genome] = {
                "resfinder_hits": hits,
                "resfinder_genes": (
                    ""
                    if pd.isna(row.get("resfinder_genes"))
                    else str(row.get("resfinder_genes", ""))
                ),
            }
    except Exception as exc:
        print(f"[safety] ResFinder parse warning: {exc}")
    return result


prokka_dir = Path(args.prokka_dir)
samples = (
    sorted(d.name for d in prokka_dir.iterdir() if d.is_dir())
    if prokka_dir.exists()
    else []
)
print(f"[safety] Processing {len(samples)} genomes")

resfinder_summary = parse_resfinder_summary(
    os.path.join(args.resfinder, "resfinder_summary.tsv")
)

rows = []
for sample in samples:
    gff_files = list((prokka_dir / sample).glob("*.gff"))
    gff_path = gff_files[0] if gff_files else None
    amr_path = Path(args.amrfinder) / f"{sample}.tsv"

    default_gff = {
        "biogenic_amines": 0,
        "biogenic_amine_pathway_types": "",
        "biogenic_amine_cluster_count": 0,
        "biogenic_amine_details": "",
        "nonspecific_decarboxylase_annotations": 0,
        "nonspecific_decarboxylase_details": "",
        "gaba_related_annotations": 0,
        "gaba_related_details": "",
        "hemolysin_markers": 0,
        "hemolysin_marker_details": "",
        "hemolysin_like_annotations": 0,
        "hemolysin_like_details": "",
    }
    gff_result = parse_gff_safety(gff_path) if gff_path else default_gff

    amr_result = (
        parse_amrfinder(amr_path)
        if amr_path.exists()
        else {
            "amr_hits": 0,
            "critical_amr_hits": 0,
            "caution_amr_hits": 0,
            "virulence_hits": 0,
            "stress_hits": 0,
            "amr_symbols": "",
            "critical_amr_symbols": "",
            "caution_amr_symbols": "",
            "virulence_symbols": "",
            "stress_symbols": "",
            "amr_details": "",
            "critical_amr_details": "",
            "caution_amr_details": "",
            "virulence_details": "",
            "stress_details": "",
        }
    )
    amr_hits = amr_result["amr_hits"]
    critical_amr_hits = amr_result["critical_amr_hits"]
    caution_amr_hits = amr_result["caution_amr_hits"]
    vf_hits = amr_result["virulence_hits"]
    stress_hits = amr_result["stress_hits"]
    res = resfinder_summary.get(
        sample,
        {"resfinder_hits": 0, "resfinder_genes": ""},
    )

    if critical_amr_hits > 0 or vf_hits > 0 or res["resfinder_hits"] > 0:
        flag = "Critical"
        reasons = []
        if critical_amr_hits > 0:
            reasons.append("High_confidence_or_specific_AMR_determinant")
        if vf_hits > 0:
            reasons.append("Virulence_gene_detected")
        if res["resfinder_hits"] > 0:
            reasons.append("ResFinder_hit")
    elif (
        caution_amr_hits > 0
        or gff_result["hemolysin_markers"] > 0
        or gff_result["biogenic_amines"] > 0
    ):
        flag = "Caution"
        reasons = []
        if caution_amr_hits > 0:
            reasons.append("Uncertain_or_plus_scope_AMR_evidence")
        if gff_result["hemolysin_markers"] > 0:
            reasons.append("Context_supported_cytolysin_operon")
        if gff_result["biogenic_amines"] > 0:
            reasons.append("Context_supported_biogenic_amine_pathway")
    else:
        flag = "Pass"
        reasons = []

    rows.append(
        {
            "genome": sample,
            "amrfinder_hits": amr_hits,
            "amrfinder_critical_hits": critical_amr_hits,
            "amrfinder_caution_hits": caution_amr_hits,
            "amrfinder_amr_symbols": amr_result["amr_symbols"],
            "amrfinder_critical_symbols": amr_result["critical_amr_symbols"],
            "amrfinder_caution_symbols": amr_result["caution_amr_symbols"],
            "amrfinder_amr_details": amr_result["amr_details"],
            "amrfinder_critical_details": amr_result["critical_amr_details"],
            "amrfinder_caution_details": amr_result["caution_amr_details"],
            "vfdb_hits": vf_hits,
            "amrfinder_virulence_symbols": amr_result["virulence_symbols"],
            "amrfinder_virulence_details": amr_result["virulence_details"],
            "amrfinder_stress_hits": stress_hits,
            "amrfinder_stress_symbols": amr_result["stress_symbols"],
            "amrfinder_stress_details": amr_result["stress_details"],
            "resfinder_hits": res["resfinder_hits"],
            "resfinder_genes": res["resfinder_genes"],
            **gff_result,
            "flag": flag,
            "reasons": ";".join(reasons),
        }
    )

columns = [
    "genome",
    "amrfinder_hits",
    "amrfinder_critical_hits",
    "amrfinder_caution_hits",
    "amrfinder_amr_symbols",
    "amrfinder_critical_symbols",
    "amrfinder_caution_symbols",
    "amrfinder_amr_details",
    "amrfinder_critical_details",
    "amrfinder_caution_details",
    "vfdb_hits",
    "amrfinder_virulence_symbols",
    "amrfinder_virulence_details",
    "amrfinder_stress_hits",
    "amrfinder_stress_symbols",
    "amrfinder_stress_details",
    "resfinder_hits",
    "resfinder_genes",
    "biogenic_amines",
    "biogenic_amine_pathway_types",
    "biogenic_amine_cluster_count",
    "biogenic_amine_details",
    "nonspecific_decarboxylase_annotations",
    "nonspecific_decarboxylase_details",
    "gaba_related_annotations",
    "gaba_related_details",
    "hemolysin_markers",
    "hemolysin_marker_details",
    "hemolysin_like_annotations",
    "hemolysin_like_details",
    "flag",
    "reasons",
]

df_out = pd.DataFrame(rows, columns=columns)
print("[safety] Flag summary:")
print(df_out["flag"].value_counts(dropna=False).to_string())

out_path = Path(args.outdir) / "safety_summary.tsv"
df_out.to_csv(out_path, sep="\t", index=False)
print(f"[safety] Saved: {out_path}")
