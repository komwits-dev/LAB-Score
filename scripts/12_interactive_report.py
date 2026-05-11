#!/usr/bin/env python3
"""
12_interactive_report.py — LAB-Score v1.1 Interactive HTML Report
"""

import argparse, os, json, subprocess
import pandas as pd
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--scores",   required=True)
ap.add_argument("--matrix",   required=True)
ap.add_argument("--metadata", required=True)
ap.add_argument("--outdir",   required=True)
args = ap.parse_args()
os.makedirs(args.outdir, exist_ok=True)

df   = pd.read_csv(args.scores,   sep="\t")
mat  = pd.read_csv(args.matrix,   sep="\t")
meta = pd.read_csv(args.metadata, sep="\t", dtype=str)
print(f"[report] {len(df)} genomes loaded")

# Detect tool versions
def get_ver(cmd):
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=5)
        out = (r.stdout + r.stderr).strip().split("\n")[0]
        return out[:60]
    except Exception:
        return "not detected"

tool_versions = {
    "Prokka":        get_ver(["prokka", "--version"]),
    "AMRFinderPlus": get_ver(["amrfinder", "--version"]),
    "dbCAN":         get_ver(["run_dbcan", "--version"]),
    "Diamond":       get_ver(["diamond", "--version"]),
    "HMMER":         get_ver(["hmmscan", "-h"]),
    "FastANI":       get_ver(["fastANI", "--version"]),
    "Python":        get_ver(["python", "--version"]),
}

PANELS = ["acid_energy","bileresistance","gutpersistence","osmoticstress",
          "heatstress","coldstress","gaba","vitamins","immunomodulation",
          "bacteriocins","cellenvelope_eps","adhesion_surface","adhesion_biofilm",
          "defense_crispr","antipath_qs","alkalinestress","carbohydrate","metabolism"]

PANEL_LABELS = {
    "acid_energy.acid_energy":"Acid/Energy","bileresistance.bileresistance":"Bile Resistance",
    "gutpersistence.gutpersistence":"Gut Persistence","osmoticstress.osmoticstress":"Osmotic Stress",
    "heatstress.heatstress":"Heat Stress","coldstress.coldstress":"Cold Stress",
    "gaba.gaba":"GABA","vitamins.vitamins":"Vitamins",
    "immunomodulation.immunomodulation":"Immunomodulation","bacteriocins.bacteriocin":"Bacteriocins",
    "cellenvelope_eps.cellenvelope_eps":"EPS/Cell Env.","adhesion_surface.adhesion_surface":"Adhesion Surface",
    "adhesion_biofilm.adhesion_biofilm":"Adhesion Biofilm","defense_crispr.defense_crispr":"CRISPR Defense",
    "antipath_qs.antipath_qs":"Anti-pathogen QS","alkalinestress.alkalinestress":"Alkali Stress",
    "carbohydrate.carbohydrate":"Carbohydrate","metabolism.metabolism":"Fermentation"
}

PANEL_COLS_FINAL = [c for c in list(PANEL_LABELS.keys()) if c in df.columns or c in mat.columns]

for col in PANEL_COLS_FINAL:
    if col not in df.columns and col in mat.columns and "genome" in mat.columns:
        df = df.merge(mat[["genome", col]], on="genome", how="left")

def safe(v):
    if v is None or (isinstance(v, float) and np.isnan(v)): return None
    if isinstance(v, (np.integer, int)) and not isinstance(v, bool): return int(v)
    if isinstance(v, (np.floating, float)): return round(float(v), 4)
    return str(v)

all_cols = list(df.columns)
rows_data = [{c: safe(row.get(c)) for c in all_cols} for _, row in df.iterrows()]

tier_counts = df["Candidate_tier_v1_1"].value_counts().to_dict() if "Candidate_tier_v1_1" in df.columns else {}
pc_counts   = df["Priority_class"].value_counts().to_dict() if "Priority_class" in df.columns else {}
sp_col      = "full_species" if "full_species" in df.columns else "genome"
sp_counts   = df[sp_col].value_counts().head(20).to_dict()
score_hist  = df["LAB_score_v1"].dropna().tolist() if "LAB_score_v1" in df.columns else []
panel_maxes = {pc: float(df[pc].max()) if pc in df.columns and df[pc].max() > 0 else 1 for pc in PANEL_COLS_FINAL}

# Radar-ready compact data for all strains. The browser uses this to search/select
# strains for interactive comparison without restricting to only the top candidates.
radar_data = []
for _, row in df.iterrows():
    sp  = str(row.get("full_species", ""))
    gn  = str(row.get("genome", ""))
    st  = str(row.get("strain", ""))
    nm_base = sp if sp not in ("", "nan", "Unknown sp.") else gn[:35]
    nm = (nm_base + (f" | {st}" if st not in ("", "nan", "None") else ""))[:90]
    vals = [float(row.get(pc, 0) or 0) for pc in PANEL_COLS_FINAL]
    radar_data.append({
        "name": nm,
        "genome": gn,
        "species": sp,
        "strain": st,
        "score": safe(row.get("LAB_score_v1")),
        "tier":  str(row.get("Candidate_tier_v1_1", "")),
        "safety": str(row.get("Refined_safety_status", "")),
        "values": vals
    })

DATA = {
    "rows": rows_data, "all_cols": all_cols,
    "tier_counts": tier_counts, "pc_counts": pc_counts,
    "sp_counts": sp_counts, "score_hist": score_hist,
    "radar_data": radar_data,
    "panel_cols": PANEL_COLS_FINAL,
    "panel_labels": [PANEL_LABELS.get(c, c) for c in PANEL_COLS_FINAL],
    "panel_maxes": panel_maxes,
    "n_genomes": len(df),
    "tool_versions": tool_versions,
    "tier_colors": {
        "Elite candidate":"#1B4332","High candidate":"#2D6A4F",
        "Moderate candidate":"#52B788","Low priority":"#D8F3DC",
        "Cautionary Elite candidate":"#74C69D",
        "Cautionary High candidate":"#E07C24",
        "Cautionary Moderate candidate":"#F4A261",
        "Cautionary Low candidate":"#FFDDD2",
        "Critical safety review":"#C1121F",
    }
}

json_str = json.dumps(DATA, allow_nan=False, default=str)

# ── Embed or fallback Plotly ──────────────────────────────────────────────────
import urllib.request as _ur
_plotly_path = os.path.join(os.path.dirname(os.path.abspath(args.outdir)), "plotly.min.js")
_plotly_js = ""
for _p in [_plotly_path,
           os.path.join(os.path.dirname(__file__), "plotly.min.js"),
           os.path.join(args.outdir, "plotly.min.js")]:
    if os.path.exists(_p):
        with open(_p) as _f: _plotly_js = _f.read()
        print(f"[report] Embedding Plotly from {_p}")
        break
if not _plotly_js:
    print("[report] Downloading Plotly.js (~3MB) for offline embedding...")
    try:
        _dl = os.path.join(args.outdir, "plotly.min.js")
        _ur.urlretrieve("https://cdn.plot.ly/plotly-2.27.0.min.js", _dl)
        with open(_dl) as _f: _plotly_js = _f.read()
        print("[report] Plotly downloaded and embedded.")
    except Exception as _e:
        print(f"[report] WARNING: Could not download Plotly ({_e}). Using CDN — requires internet.")
        _plotly_js = "/* Plotly CDN fallback */\nconst _s=document.createElement('script');_s.src='https://cdn.plot.ly/plotly-2.27.0.min.js';document.head.appendChild(_s);" 

# Replace placeholder
HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>LAB-Score v1.1 — Interactive Report</title>
<script>{PLOTLY_JS}</script>
<style>
:root{
  --gd:#1B4332;--gm:#2D6A4F;--gl:#52B788;--gp:#D8F3DC;
  --am:#E07C24;--rd:#C1121F;--bg:#f5f7f5;--cb:#fff;
  --bd:#d0e4d0;--tx:#1a2e1a;--tl:#4a6a4a;
}
*{box-sizing:border-box;margin:0;padding:0}
body{font-family:'Segoe UI',system-ui,sans-serif;background:var(--bg);color:var(--tx);}

nav{background:var(--gd);color:#fff;padding:0 24px;height:54px;display:flex;align-items:center;
  position:sticky;top:0;z-index:200;box-shadow:0 2px 8px rgba(0,0,0,.3)}
nav h1{font-size:1.05rem;font-weight:700}
.nsub{font-size:.75rem;opacity:.7;margin-left:10px}
.ntabs{margin-left:auto;display:flex;gap:4px}
.ntab{padding:5px 14px;border-radius:6px;cursor:pointer;font-size:.8rem;color:rgba(255,255,255,.8);
  border:1px solid transparent;transition:all .2s;white-space:nowrap}
.ntab:hover,.ntab.on{background:rgba(255,255,255,.18);color:#fff;border-color:rgba(255,255,255,.3)}

.page{display:none;padding:20px;max-width:1440px;margin:0 auto}
.page.on{display:block}

.card{background:var(--cb);border:1px solid var(--bd);border-radius:12px;
  padding:18px;margin-bottom:18px;box-shadow:0 1px 4px rgba(0,0,0,.05)}
.card h2{font-size:.95rem;font-weight:700;color:var(--gd);margin-bottom:12px}
.card p{font-size:.8rem;color:var(--tl);margin-bottom:8px;line-height:1.6}

.srow{display:flex;gap:14px;margin-bottom:18px;flex-wrap:wrap}
.scard{flex:1;min-width:110px;background:var(--cb);border:1px solid var(--bd);
  border-radius:10px;padding:14px;text-align:center}
.scard .v{font-size:1.9rem;font-weight:800;color:var(--gd)}
.scard .l{font-size:.68rem;color:var(--tl);margin-top:3px;text-transform:uppercase;letter-spacing:.5px}

.g2{display:grid;grid-template-columns:1fr 1fr;gap:18px}
.g3{display:grid;grid-template-columns:1fr 1fr 1fr;gap:18px}
@media(max-width:900px){.g2,.g3{grid-template-columns:1fr}}

/* Table */
.sw{display:flex;gap:10px;margin-bottom:14px;align-items:center;flex-wrap:wrap}
input[type=text]{flex:1;min-width:180px;padding:7px 11px;border:1px solid var(--bd);
  border-radius:7px;font-size:.83rem;outline:none}
input[type=text]:focus{border-color:var(--gm);box-shadow:0 0 0 3px rgba(82,183,136,.15)}
select{padding:7px 10px;border:1px solid var(--bd);border-radius:7px;font-size:.82rem;
  background:#fff;cursor:pointer}
.rcount{font-size:.75rem;color:var(--tl);white-space:nowrap}

.tw{overflow-x:auto;max-height:600px;overflow-y:auto}
table{width:100%;border-collapse:collapse;font-size:.79rem}
thead tr{background:var(--gd);color:#fff;position:sticky;top:0;z-index:10}
th{padding:9px 11px;text-align:left;font-weight:600;cursor:pointer;white-space:nowrap;user-select:none}
th:hover{background:var(--gm)}
th.asc::after{content:' ▲';font-size:.65em}
th.desc::after{content:' ▼';font-size:.65em}
tbody tr{border-bottom:1px solid var(--bd);cursor:pointer;transition:background .12s}
tbody tr:hover{background:#f0f7f0}
tbody tr.sel{background:#dff0df!important;outline:2px solid var(--gl) inset}
td{padding:7px 11px;white-space:nowrap}

.badge{display:inline-block;padding:2px 8px;border-radius:20px;font-size:.69rem;font-weight:700;white-space:nowrap}
.be{background:#1B4332;color:#fff}.bh{background:#2D6A4F;color:#fff}
.bm{background:#95D5B2;color:#1B4332}.bl{background:#D8F3DC;color:#2D6A4F}
.bc{background:#C1121F;color:#fff}.bca{background:#F4A261;color:#fff}
.bp{background:#D8F3DC;color:#1B4332}

.sb{display:flex;align-items:center;gap:6px}
.sbo{height:7px;border-radius:4px;background:#e0e0e0;width:70px}
.sbi{height:100%;border-radius:4px}

/* MODAL */
#overlay{position:fixed;inset:0;background:rgba(0,0,0,.5);z-index:300;
  display:none;align-items:flex-start;justify-content:center;padding:20px;overflow-y:auto}
#overlay.on{display:flex}
#modal{background:#fff;border-radius:14px;width:min(1100px,98vw);
  margin:auto;box-shadow:0 20px 60px rgba(0,0,0,.3);display:flex;flex-direction:column}
.mhead{background:var(--gd);color:#fff;padding:16px 22px;border-radius:14px 14px 0 0;
  display:flex;align-items:center;gap:10px;position:sticky;top:0;z-index:10}
.mhead h2{font-size:1rem;font-weight:700;flex:1}
.mbtns{display:flex;gap:8px}
.mbtn{padding:6px 14px;border-radius:7px;border:1px solid rgba(255,255,255,.4);
  cursor:pointer;font-size:.78rem;font-weight:600;color:#fff;background:rgba(255,255,255,.15);
  transition:all .2s}
.mbtn:hover{background:rgba(255,255,255,.3)}
.mbtn.export{background:#52B788;border-color:#52B788}
.mbtn.export:hover{background:#40a070}
.mclose{cursor:pointer;font-size:1.4rem;opacity:.8;padding:0 6px;border-radius:6px}
.mclose:hover{background:rgba(255,255,255,.2)}
.mbody{padding:20px;display:grid;grid-template-columns:1fr 1fr;gap:18px}
@media(max-width:750px){.mbody{grid-template-columns:1fr}}
.msec{background:#f8faf8;border:1px solid var(--bd);border-radius:10px;padding:14px}
.msec h3{font-size:.78rem;font-weight:700;color:var(--gd);margin-bottom:10px;
  text-transform:uppercase;letter-spacing:.6px}
.mrow{display:flex;justify-content:space-between;align-items:center;
  padding:4px 0;border-bottom:1px solid #eef2ee;font-size:.79rem}
.mrow:last-child{border-bottom:none}
.mk{color:var(--tl);font-weight:500}.mv{font-weight:600;text-align:right;max-width:58%;word-break:break-word}
.alert{padding:8px 12px;border-radius:8px;font-size:.78rem;font-weight:600;margin-bottom:6px;display:flex;align-items:center;gap:6px}
.ared{background:#fee;border:1px solid #fcc;color:#900}
.aamb{background:#fff3e0;border:1px solid #ffe0b2;color:#7a4000}
.agrn{background:#e8f5e9;border:1px solid #c8e6c9;color:#1b5e20}
.pbar-row{display:flex;align-items:center;gap:8px;padding:3px 0;font-size:.76rem}
.pbar-lbl{width:145px;color:var(--tl);text-align:right;flex-shrink:0;font-size:.74rem}
.pbar-out{flex:1;height:10px;background:#e8e8e8;border-radius:5px;overflow:hidden}
.pbar-in{height:100%;border-radius:5px}
.pbar-v{width:28px;text-align:right;font-weight:700;font-size:.74rem}
.mfull{grid-column:span 2}


/* Radar compare */
.radar-grid{display:grid;grid-template-columns:360px 1fr;gap:18px;align-items:start}
@media(max-width:1050px){.radar-grid{grid-template-columns:1fr}}
.radar-controls{display:flex;gap:8px;flex-wrap:wrap;margin-bottom:10px}
.radar-btn{border:1px solid var(--bd);background:#fff;color:var(--gd);border-radius:7px;padding:7px 10px;cursor:pointer;font-size:.76rem;font-weight:700}
.radar-btn:hover{background:#f0f7f0;border-color:var(--gl)}
.radar-btn.primary{background:var(--gm);color:#fff;border-color:var(--gm)}
.radar-btn.warn{background:#fff3e0;color:#7a4000;border-color:#ffe0b2}
.radar-list{border:1px solid var(--bd);border-radius:10px;max-height:590px;overflow:auto;background:#fff}
.radar-item{padding:9px 10px;border-bottom:1px solid #eef2ee;cursor:pointer;display:grid;grid-template-columns:22px 1fr auto;gap:8px;align-items:center}
.radar-item:last-child{border-bottom:none}
.radar-item:hover{background:#f4faf4}
.radar-item.sel{background:#e4f4e4;outline:2px solid var(--gl);outline-offset:-2px}
.radar-title{font-weight:700;color:var(--gd);font-size:.78rem;line-height:1.25;word-break:break-word}
.radar-sub{font-size:.69rem;color:var(--tl);line-height:1.3;margin-top:3px;word-break:break-word}
.radar-view{border:0;background:#1B4332;color:#fff;border-radius:6px;padding:4px 7px;font-size:.68rem;cursor:pointer}
.radar-chipbar{display:flex;gap:6px;flex-wrap:wrap;margin:8px 0 12px}
.radar-chip{background:#e8f5e9;color:#1B4332;border:1px solid #c8e6c9;border-radius:14px;padding:3px 8px;font-size:.72rem;font-weight:700}
.radar-export{display:flex;gap:8px;flex-wrap:wrap;margin:10px 0}
.radar-note{font-size:.76rem;color:var(--tl);line-height:1.5;margin-bottom:10px}
#radar-selected-table{max-height:260px;overflow:auto;margin-top:10px;border:1px solid var(--bd);border-radius:9px}
#radar-selected-table table{font-size:.74rem}
#radar-selected-table tbody tr{cursor:pointer}

/* About page */
.step-card{background:var(--cb);border:1px solid var(--bd);border-radius:10px;
  padding:16px;margin-bottom:12px;display:flex;gap:14px}
.step-num{width:34px;height:34px;border-radius:50%;background:var(--gd);color:#fff;
  font-weight:800;font-size:.9rem;display:flex;align-items:center;justify-content:center;flex-shrink:0}
.step-body h3{font-size:.88rem;font-weight:700;color:var(--gd);margin-bottom:4px}
.step-body p{font-size:.78rem;color:var(--tl);line-height:1.55}
.step-body .tag{display:inline-block;background:var(--gp);color:var(--gd);
  padding:1px 7px;border-radius:10px;font-size:.68rem;font-weight:600;margin-right:4px;margin-top:4px}
.tool-row{display:flex;align-items:center;justify-content:space-between;
  padding:8px 0;border-bottom:1px solid var(--bd);font-size:.8rem}
.tool-row:last-child{border-bottom:none}
.tool-name{font-weight:700;color:var(--gd)}
.tool-ver{color:var(--tl);font-size:.73rem;max-width:60%;text-align:right;word-break:break-word}
.formula-box{background:#f0f7f0;border:1px solid var(--gl);border-radius:8px;
  padding:14px 18px;font-family:monospace;font-size:.82rem;line-height:1.9;color:var(--gd)}
.hctrl{display:flex;gap:10px;align-items:center;margin-bottom:12px;flex-wrap:wrap;font-size:.8rem}

footer{text-align:center;padding:20px;font-size:.72rem;color:var(--tl);
  border-top:1px solid var(--bd);margin-top:30px}
</style>
</head>
<body>

<nav>
  <h1>🧫 LAB-Score v1.1</h1>
  <span class="nsub">Interactive Genomic Analysis Report</span>
  <div class="ntabs">
    <div class="ntab on"  data-page="overview">📊 Overview</div>
    <div class="ntab"     data-page="strains">🔍 Strains</div>
    <div class="ntab"     data-page="radar">🎯 Radar Compare</div>
    <div class="ntab"     data-page="heatmap">🧬 Gene Presence</div>
    <div class="ntab"     data-page="about">ℹ️ About / Pipeline</div>
  </div>
</nav>

<!-- STRAIN DETAIL MODAL -->
<div id="overlay" onclick="closeModal(event)">
  <div id="modal">
    <div class="mhead">
      <h2 id="modal-title">Strain Details</h2>
      <div class="mbtns">
        <button class="mbtn export" onclick="exportStrain()">⬇ Export TSV</button>
        <button class="mbtn" onclick="exportStrainJSON()">⬇ Export JSON</button>
        <button class="mbtn" onclick="exportCurrentRadarPNG()">⬇ Export Radar</button>
        <button class="mbtn" onclick="exportStrainReportHTML()">⬇ Export Report</button>
      </div>
      <span class="mclose" onclick="closeModal()">✕</span>
    </div>
    <div class="mbody" id="modal-body"></div>
  </div>
</div>

<!-- PAGE: OVERVIEW -->
<div id="page-overview" class="page on">
  <div class="srow" id="stat-row"></div>
  <div class="g2">
    <div class="card"><h2>📊 LAB-Score Distribution</h2><div id="ch-hist" style="height:270px"></div></div>
    <div class="card"><h2>🏆 Priority Classes</h2><div id="ch-pie" style="height:270px"></div></div>
  </div>
  <div class="g2">
    <div class="card"><h2>🛡️ Safety-Gated Candidate Tiers</h2><div id="ch-tiers" style="height:310px"></div></div>
    <div class="card"><h2>🦠 Top Species</h2><div id="ch-sp" style="height:310px"></div></div>
  </div>
  <div class="card"><h2>🔬 Functional Module Overview</h2><div id="ch-modules" style="height:290px"></div></div>
</div>

<!-- PAGE: STRAINS -->
<div id="page-strains" class="page">
  <div class="card">
    <h2>🔍 Strain Browser <span style="font-weight:400;font-size:.8rem;color:var(--tl)">— click any row for full strain report + export</span></h2>
    <div class="sw">
      <input type="text" id="srch" placeholder="Search genome, species, tier..." oninput="filterTbl()">
      <select id="ftier" onchange="filterTbl()"><option value="">All tiers</option></select>
      <select id="fsafe" onchange="filterTbl()">
        <option value="">All safety</option>
        <option value="Pass">Pass</option>
        <option value="Caution">Caution</option>
        <option value="Critical">Critical</option>
      </select>
      <span class="rcount" id="rcount"></span>
    </div>
    <div class="tw"><table><thead><tr id="thead"></tr></thead><tbody id="tbody"></tbody></table></div>
  </div>
</div>


<!-- PAGE: RADAR COMPARE -->
<div id="page-radar" class="page">
  <div class="card">
    <h2>🎯 Radar Compare <span style="font-weight:400;font-size:.8rem;color:var(--tl)">— select any strains, click View for full report, export figure/table/report</span></h2>
    <p class="radar-note">Use the search box and filters to select strains for side-by-side functional module comparison. Values are normalized as percentage of the dataset maximum for each module, so profiles are comparable across strains.</p>
    <div class="radar-grid">
      <div>
        <div class="radar-controls">
          <input type="text" id="radar-search" placeholder="Search genome, species, strain, tier..." oninput="renderRadarSelector()">
          <select id="radar-tier" onchange="renderRadarSelector()"><option value="">All tiers</option></select>
          <select id="radar-safety" onchange="renderRadarSelector()">
            <option value="">All safety</option>
            <option value="Pass">Pass</option>
            <option value="Caution">Caution</option>
            <option value="Critical">Critical</option>
          </select>
        </div>
        <div class="radar-controls">
          <button class="radar-btn primary" onclick="radarPreset('top5')">Top 5</button>
          <button class="radar-btn" onclick="radarPreset('top10')">Top 10</button>
          <button class="radar-btn" onclick="radarPreset('elite')">Elite</button>
          <button class="radar-btn warn" onclick="clearRadarSelection()">Clear</button>
        </div>
        <div class="rcount" id="radar-count"></div>
        <div id="radar-list" class="radar-list"></div>
      </div>
      <div>
        <div class="radar-chipbar" id="radar-chips"></div>
        <div class="radar-export">
          <button class="radar-btn primary" onclick="exportRadarPNG()">⬇ Export figure PNG</button>
          <button class="radar-btn" onclick="exportRadarSVG()">⬇ Export figure SVG</button>
          <button class="radar-btn" onclick="exportRadarTable()">⬇ Export comparison table</button>
          <button class="radar-btn" onclick="exportRadarReport()">⬇ Export selected-strain report</button>
        </div>
        <div id="ch-radar-compare" style="height:560px"></div>
        <div id="radar-selected-table"></div>
      </div>
    </div>
  </div>
</div>

<!-- PAGE: GENE PRESENCE HEATMAP -->
<div id="page-heatmap" class="page">
  <div class="card">
    <h2>🧬 Gene Presence / Abundance Heatmap</h2>
    <div class="hctrl">
      <label><b>Mode:</b></label>
      <select id="hmode" onchange="renderHeatmap()">
        <option value="count">Raw counts</option>
        <option value="presence">Presence / Absence</option>
        <option value="zscore">Z-score</option>
      </select>
      <label><b>Sort:</b></label>
      <select id="hsort" onchange="renderHeatmap()">
        <option value="score">LAB-Score</option>
        <option value="tier">Tier</option>
        <option value="species">Species</option>
      </select>
      <label><b>Show:</b></label>
      <select id="hn" onchange="renderHeatmap()">
        <option value="50">Top 50</option>
        <option value="100">Top 100</option>
        <option value="200">Top 200</option>
        <option value="0">All</option>
      </select>
    </div>
    <div id="ch-hmap" style="height:640px"></div>
  </div>
  <div class="g3">
    <div class="card"><h2>🔋 GI Survival Modules</h2><div id="ch-gi" style="height:260px"></div></div>
    <div class="card"><h2>⚡ Functional Modules</h2><div id="ch-fu" style="height:260px"></div></div>
    <div class="card"><h2>🍞 Fermentation Modules</h2><div id="ch-fe" style="height:260px"></div></div>
  </div>
</div>

<!-- PAGE: ABOUT / PIPELINE -->
<div id="page-about" class="page">

  <div class="g2">
    <div>
      <div class="card">
        <h2>🧫 About LAB-Score v1.1</h2>
        <p>LAB-Score v1.1 is a genome-based scoring framework for prioritizing <i>Lactobacillaceae</i> strains as probiotic or fermentation candidates. It integrates safety screening, GI survival capacity, probiotic-associated functional markers, and fermentation traits into a single interpretable composite score.</p>
        <p>The pipeline was developed to enable high-throughput, reproducible evaluation of large collections of formerly <i>Lactobacillus</i> genomes (now reorganized into 23+ genera), supporting downstream experimental validation decisions.</p>
        <p style="margin-top:8px"><b>Reference:</b> LAB-Score: An interpretable genome-based scoring framework for prioritizing lactic acid bacteria candidates. <i>In preparation.</i></p>
      </div>

      <div class="card">
        <h2>📐 Scoring Formula</h2>
        <div class="formula-box">
LAB_score_v1 =<br>
&nbsp;&nbsp;&nbsp;0.45 × Safety_score<br>
&nbsp;&nbsp;+ 0.25 × GI_survival_score<br>
&nbsp;&nbsp;+ 0.20 × Functional_score<br>
&nbsp;&nbsp;+ 0.10 × Fermentation_score
        </div>
        <p style="margin-top:10px"><b>Priority thresholds:</b></p>
        <div style="margin-top:6px">
          <div style="display:flex;gap:8px;flex-wrap:wrap">
            <span class="badge be">≥85 Elite</span>
            <span class="badge bh">70–84 High</span>
            <span class="badge bm">50–69 Moderate</span>
            <span class="badge bl">&lt;50 Low</span>
          </div>
        </div>
        <p style="margin-top:10px"><b>Safety gate:</b> Genomes with AMR genes, virulence factors, or ResFinder hits receive a "Critical" or "Cautionary" prefix regardless of score.</p>
      </div>

      <div class="card">
        <h2>🔧 Tool Versions</h2>
        <div id="tool-table"></div>
      </div>
    </div>

    <div>
      <div class="card">
        <h2>🔄 Pipeline Steps</h2>

        <div class="step-card">
          <div class="step-num">0</div>
          <div class="step-body">
            <h3>Species Metadata Resolution</h3>
            <p>Automatically fetches species names from NCBI Datasets API for GCF/GCA accessions. Falls back to filename parsing for non-accession filenames.</p>
            <span class="tag">NCBI API</span><span class="tag">GCF/GCA</span>
          </div>
        </div>

        <div class="step-card">
          <div class="step-num">0b</div>
          <div class="step-body">
            <h3>Species Verification (FastANI)</h3>
            <p>Genome-based species identity verification using Average Nucleotide Identity (ANI ≥95% = same species). Flags mismatches between NCBI metadata and genomic evidence.</p>
            <span class="tag">FastANI</span><span class="tag">ANI ≥95%</span>
          </div>
        </div>

        <div class="step-card">
          <div class="step-num">1</div>
          <div class="step-body">
            <h3>Genome Annotation (Prokka)</h3>
            <p>Rapid prokaryotic genome annotation using Prokka. Generates GFF, FAA (protein), and FFN (nucleotide) files for downstream analysis.</p>
            <span class="tag">Prokka</span><span class="tag">Prodigal</span><span class="tag">HMMER</span>
          </div>
        </div>

        <div class="step-card">
          <div class="step-num">2</div>
          <div class="step-body">
            <h3>AMR &amp; Virulence Screening (AMRFinder)</h3>
            <p>Identifies acquired antimicrobial resistance (AMR) genes, virulence factors, and stress response genes using NCBI AMRFinderPlus with the curated reference database.</p>
            <span class="tag">AMRFinderPlus</span><span class="tag">NCBI database</span>
          </div>
        </div>

        <div class="step-card">
          <div class="step-num">3</div>
          <div class="step-body">
            <h3>Acquired Resistance (ResFinder)</h3>
            <p>Screens for acquired antibiotic resistance genes using the ResFinder database. Results contribute to the Critical safety flag.</p>
            <span class="tag">ResFinder</span><span class="tag">Optional</span>
          </div>
        </div>

        <div class="step-card">
          <div class="step-num">4</div>
          <div class="step-body">
            <h3>CAZyme Annotation (dbCAN)</h3>
            <p>Identifies carbohydrate-active enzymes (CAZymes) using dbCAN. CAZyme diversity contributes to the fermentation score, reflecting carbohydrate utilization capacity.</p>
            <span class="tag">dbCAN v5</span><span class="tag">Diamond</span><span class="tag">HMMER</span>
          </div>
        </div>

        <div class="step-card">
          <div class="step-num">5</div>
          <div class="step-body">
            <h3>Genome Quality (CheckM) — Optional</h3>
            <p>Estimates genome completeness and contamination using lineage-specific marker gene sets. Used to filter low-quality assemblies.</p>
            <span class="tag">CheckM</span><span class="tag">Optional (-q)</span>
          </div>
        </div>

        <div class="step-card">
          <div class="step-num">6</div>
          <div class="step-body">
            <h3>Safety Screening</h3>
            <p>GFF-based scan for biogenic amine decarboxylase genes and hemolysin markers. Combined with AMRFinder and ResFinder results to assign safety flags: Pass / Caution / Critical.</p>
            <span class="tag">GFF scan</span><span class="tag">Regex patterns</span>
          </div>
        </div>

        <div class="step-card">
          <div class="step-num">7</div>
          <div class="step-body">
            <h3>Functional Panel Scanning (18 panels)</h3>
            <p>Scans Prokka GFF annotation for 18 curated probiotic-associated functional marker panels including acid/energy metabolism, bile resistance, GABA production, vitamin biosynthesis, immunomodulation, adhesion, bacteriocins, EPS, CRISPR defense, and stress tolerance.</p>
            <span class="tag">Acid/Energy</span><span class="tag">Bile resistance</span><span class="tag">GABA</span>
            <span class="tag">Vitamins</span><span class="tag">Immunomod.</span><span class="tag">Bacteriocins</span>
            <span class="tag">EPS</span><span class="tag">CRISPR</span><span class="tag">+10 more</span>
          </div>
        </div>

        <div class="step-card">
          <div class="step-num">8–9</div>
          <div class="step-body">
            <h3>Master Matrix &amp; LAB-Score Calculation</h3>
            <p>All features merged into a master matrix. LAB-Score v1.1 calculated using the weighted formula. Refined safety gate applied to assign final candidate tiers.</p>
            <span class="tag">pandas</span><span class="tag">numpy</span>
          </div>
        </div>

        <div class="step-card">
          <div class="step-num">10</div>
          <div class="step-body">
            <h3>Machine Learning Interpretation</h3>
            <p>Random Forest regression (predict LAB-Score) and classification (predict priority class) using functional-only features. Provides feature importance and SHAP values to identify which modules drive candidate ranking.</p>
            <span class="tag">Random Forest</span><span class="tag">scikit-learn</span><span class="tag">SHAP</span>
          </div>
        </div>

        <div class="step-card">
          <div class="step-num">11</div>
          <div class="step-body">
            <h3>Publication Figures</h3>
            <p>Generates Fig 2 (score distribution), Fig 3 (candidate tiers), Fig 4 (species summary), and Fig 5 (ML feature importance) at 300 DPI for publication.</p>
            <span class="tag">matplotlib</span><span class="tag">seaborn</span><span class="tag">300 DPI</span>
          </div>
        </div>

        <div class="step-card">
          <div class="step-num">12</div>
          <div class="step-body">
            <h3>Interactive HTML Report</h3>
            <p>This self-contained interactive report. No server required — open in any web browser.</p>
            <span class="tag">Plotly.js</span><span class="tag">Single HTML file</span>
          </div>
        </div>

      </div><!-- end steps card -->
    </div>
  </div><!-- end g2 -->
</div><!-- end about page -->

<footer>LAB-Score v1.1 Interactive Report — Generated by LAB-Score Pipeline v3.5</footer>

<script>
const D = {json_str};
const TC = D.tier_colors;
const PL = D.panel_labels;
const PC = D.panel_cols;
const PM = D.panel_maxes;

function showPage(name) {
  document.querySelectorAll('.page').forEach(p => p.classList.remove('on'));
  document.querySelectorAll('.ntab').forEach(t => t.classList.remove('on'));
  const pg = document.getElementById('page-' + name);
  if (pg) pg.classList.add('on');
  const tab = document.querySelector('[data-page="' + name + '"]');
  if (tab) tab.classList.add('on');
  if (name === 'heatmap') setTimeout(renderHeatmap, 50);
  if (name === 'radar') setTimeout(renderRadar, 50);
}

// Tab click handler — attached after DOM ready
document.addEventListener('DOMContentLoaded', function() {
  document.querySelectorAll('.ntab[data-page]').forEach(function(tab) {
    tab.addEventListener('click', function() {
      showPage(this.getAttribute('data-page'));
    });
  });
});

function tierBadge(t) {
  if (!t) return '';
  const tl = t.toLowerCase();
  let c = 'bm';
  if (tl.includes('critical')) c = 'bc';
  else if (tl.includes('caution')) c = 'bca';
  else if (tl.includes('elite')) c = 'be';
  else if (tl.includes('high')) c = 'bh';
  else if (tl.includes('low')) c = 'bl';
  return `<span class="badge ${c}">${t}</span>`;
}
function safetyBadge(f) {
  if (!f) return '';
  return `<span class="badge ${f==='Critical'?'bc':f==='Caution'?'bca':'bp'}">${f}</span>`;
}
function scorebar(v) {
  const p = Math.min(100, Math.max(0, v||0));
  const col = p>=85?'#1B4332':p>=70?'#2D6A4F':p>=50?'#52B788':'#B7E4C7';
  return `<div class="sb"><span style="width:34px;font-weight:700;font-size:.8rem;color:${col}">${p.toFixed(1)}</span>
    <div class="sbo"><div class="sbi" style="width:${p}%;background:${col}"></div></div></div>`;
}

// ── Overview ──────────────────────────────────────────────────────────────────
function renderOverview() {
  const elite=D.pc_counts['Elite']||0, high=D.pc_counts['High']||0;
  const crit=D.tier_counts['Critical safety review']||0;
  const cauts=Object.entries(D.tier_counts).filter(([k])=>k.includes('Cautionary')).reduce((a,[,v])=>a+v,0);
  document.getElementById('stat-row').innerHTML=[
    [D.n_genomes,'Total Genomes'],[elite,'Elite'],[high,'High'],
    [D.pc_counts['Moderate']||0,'Moderate'],[crit,'Critical Review'],[cauts,'Cautionary'],
  ].map(([v,l])=>`<div class="scard"><div class="v">${v}</div><div class="l">${l}</div></div>`).join('');

  const cfg={responsive:true,displayModeBar:false};
  const lay=(o={})=>Object.assign({margin:{t:10,b:40,l:40,r:20},paper_bgcolor:'transparent',plot_bgcolor:'transparent',font:{size:11}},o);

  Plotly.newPlot('ch-hist',[{x:D.score_hist,type:'histogram',nbinsx:50,
    marker:{color:D.score_hist.map(v=>v>=85?'#1B4332':v>=70?'#2D6A4F':v>=50?'#52B788':'#B7E4C7'),line:{color:'#fff',width:.5}},
    hovertemplate:'Score: %{x:.1f}<br>Count: %{y}<extra></extra>'}],
    lay({shapes:[50,70,85].map((x,i)=>{return{type:'line',x0:x,x1:x,y0:0,y1:1,yref:'paper',
      line:{color:['#E07C24','#2D6A4F','#1B4332'][i],dash:'dash',width:1.8}}}),
      xaxis:{title:'LAB-Score',gridcolor:'#eee'},yaxis:{title:'Count',gridcolor:'#eee'}}),cfg);

  const pkeys=Object.keys(D.pc_counts);
  Plotly.newPlot('ch-pie',[{labels:pkeys,values:Object.values(D.pc_counts),type:'pie',
    marker:{colors:pkeys.map(k=>k==='Elite'?'#1B4332':k==='High'?'#2D6A4F':k==='Moderate'?'#52B788':'#D8F3DC'),line:{color:'#fff',width:2}},
    textinfo:'label+percent',hovertemplate:'%{label}<br>n=%{value}<extra></extra>'}],
    lay({margin:{t:10,b:10,l:10,r:10},showlegend:false}),cfg);

  const tk=Object.keys(D.tier_counts).sort((a,b)=>D.tier_counts[a]-D.tier_counts[b]);
  Plotly.newPlot('ch-tiers',[{y:tk,x:tk.map(k=>D.tier_counts[k]),type:'bar',orientation:'h',
    marker:{color:tk.map(k=>TC[k]||'#aaa')},text:tk.map(k=>D.tier_counts[k]),textposition:'outside',
    hovertemplate:'%{y}<br>n=%{x}<extra></extra>'}],
    lay({margin:{t:10,b:40,l:230,r:60},xaxis:{title:'Genomes',gridcolor:'#eee'},yaxis:{automargin:true}}),cfg);

  const sk=Object.keys(D.sp_counts).slice(0,15).reverse();
  Plotly.newPlot('ch-sp',[{y:sk,x:sk.map(k=>D.sp_counts[k]),type:'bar',orientation:'h',
    marker:{color:'#52B788',line:{color:'#fff',width:1}},text:sk.map(k=>D.sp_counts[k]),textposition:'outside',
    hovertemplate:'%{y}<br>n=%{x}<extra></extra>'}],
    lay({margin:{t:10,b:40,l:200,r:60},xaxis:{title:'Genomes',gridcolor:'#eee'},yaxis:{automargin:true}}),cfg);

  const means=PL.map((_,i)=>{const c=PC[i];return c?D.rows.reduce((a,r)=>a+(r[c]||0),0)/D.rows.length:0});
  const ord=[...means.keys()].sort((a,b)=>means[b]-means[a]);
  Plotly.newPlot('ch-modules',[{x:ord.map(i=>PL[i]),y:ord.map(i=>means[i]),type:'bar',
    marker:{color:ord.map(i=>means[i]>means.reduce((a,b)=>a+b,0)/means.length?'#2D6A4F':'#95D5B2')},
    hovertemplate:'%{x}<br>Mean: %{y:.2f}<extra></extra>'}],
    lay({margin:{t:10,b:110,l:50,r:10},xaxis:{tickangle:-45,automargin:true},yaxis:{title:'Mean count',gridcolor:'#eee'}}),cfg);
}

// ── About page tool versions ──────────────────────────────────────────────────
function renderAbout() {
  const tv = D.tool_versions || {};
  document.getElementById('tool-table').innerHTML =
    Object.entries(tv).map(([name,ver])=>
      `<div class="tool-row"><span class="tool-name">${name}</span><span class="tool-ver">${ver||'—'}</span></div>`
    ).join('');
}

// ── Strain table ──────────────────────────────────────────────────────────────
const TCOLS=[
  {k:'genome',lb:'Genome'},{k:'full_species',lb:'Species'},
  {k:'LAB_score_v1',lb:'LAB-Score'},{k:'Candidate_tier_v1_1',lb:'Tier'},
  {k:'Refined_safety_status',lb:'Safety'},{k:'Safety_score',lb:'Safety Sc.'},
  {k:'GI_survival_score',lb:'GI Score'},{k:'Functional_score',lb:'Func. Score'},
  {k:'Fermentation_score',lb:'Ferm. Score'},{k:'amrfinder_hits',lb:'AMR'},
  {k:'vfdb_hits',lb:'VF'},{k:'biogenic_amines',lb:'Bio.Amine'},
  {k:'hemolysin_markers',lb:'Hemolysin'},{k:'ani_species',lb:'ANI Species'},
  {k:'species_status',lb:'Sp.Status'},
];
let sortK='LAB_score_v1',sortD=-1,fRows=[];

function initTable(){
  const tiers=[...new Set(D.rows.map(r=>r.Candidate_tier_v1_1).filter(Boolean))].sort();
  const sel=document.getElementById('ftier');
  tiers.forEach(t=>{const o=document.createElement('option');o.value=o.textContent=t;sel.appendChild(o)});
  document.getElementById('thead').innerHTML=TCOLS.map(c=>
    `<th onclick="sortTbl('${c.k}')" id="th-${c.k}">${c.lb}</th>`).join('');
  filterTbl();
}

function filterTbl(){
  const q=document.getElementById('srch').value.toLowerCase();
  const ti=document.getElementById('ftier').value;
  const sa=document.getElementById('fsafe').value;
  fRows=D.rows.filter(r=>{
    const txt=[r.genome,r.full_species,r.Candidate_tier_v1_1,r.Refined_safety_status,r.ani_species]
      .filter(Boolean).join(' ').toLowerCase();
    return(!q||txt.includes(q))&&(!ti||r.Candidate_tier_v1_1===ti)&&(!sa||r.Refined_safety_status===sa);
  });
  sortTbl(sortK,true);
}

function sortTbl(k,keep=false){
  if(!keep){sortD=k===sortK?-sortD:-1;sortK=k}
  document.querySelectorAll('th').forEach(t=>t.classList.remove('asc','desc'));
  const th=document.getElementById('th-'+k);
  if(th) th.classList.add(sortD>0?'asc':'desc');
  fRows.sort((a,b)=>{const av=a[k],bv=b[k];if(av==null)return 1;if(bv==null)return-1;
    return typeof av==='number'?sortD*(av-bv):sortD*String(av).localeCompare(String(bv));
  });
  renderTbody();
}

let selIdx=-1;
function renderTbody(){
  document.getElementById('rcount').textContent=`${fRows.length} / ${D.n_genomes} genomes`;
  document.getElementById('tbody').innerHTML=fRows.map(r=>{
    const gi=D.rows.indexOf(r);
    return `<tr onclick="openModal(${gi})" class="${gi===selIdx?'sel':''}">${
      TCOLS.map(c=>{
        const v=r[c.k];
        if(c.k==='LAB_score_v1') return `<td>${scorebar(v)}</td>`;
        if(c.k==='Candidate_tier_v1_1') return `<td>${tierBadge(v)}</td>`;
        if(c.k==='Refined_safety_status') return `<td>${safetyBadge(v)}</td>`;
        if(c.k==='species_status'){
          const col=v==='CONFIRMED'?'#1B4332':v==='MISMATCH'?'#C1121F':v==='GENUS_MATCH'?'#E07C24':'#888';
          return `<td><span style="color:${col};font-weight:600;font-size:.75rem">${v||'—'}</span></td>`;
        }
        if(typeof v==='number') return `<td style="text-align:right">${v.toFixed(1)}</td>`;
        return `<td>${v!=null&&v!=='nan'?v:'—'}</td>`;
      }).join('')
    }</tr>`;
  }).join('');
}


// ── RADAR COMPARE TAB ────────────────────────────────────────────────────────
let radarSelected = new Set();
const RADAR_MAX_SELECTED = 12;
const radarPalette = ['#1B4332','#2D6A4F','#52B788','#74C69D','#E07C24','#F4A261','#C1121F','#6C757D','#4361EE','#7209B7','#B5179E','#4895EF'];

function strainLabel(r, maxLen=70) {
  const sp = r.full_species && r.full_species !== 'Unknown sp.' && r.full_species !== 'nan' ? r.full_species : r.genome;
  const st = r.strain && r.strain !== 'nan' && r.strain !== 'None' ? ` | ${r.strain}` : '';
  let x = `${sp}${st}`;
  return x.length > maxLen ? x.slice(0, maxLen-1) + '…' : x;
}

function initRadar(){
  const tierSel = document.getElementById('radar-tier');
  if(!tierSel) return;
  if(tierSel.options.length <= 1){
    [...new Set(D.rows.map(r=>r.Candidate_tier_v1_1).filter(Boolean))].sort().forEach(t=>{
      const o=document.createElement('option');o.value=o.textContent=t;tierSel.appendChild(o);
    });
  }
  if(radarSelected.size===0) radarPreset('top5', true);
  renderRadarSelector();
  renderRadar();
}

function radarFilteredRows(){
  const q=(document.getElementById('radar-search')?.value||'').toLowerCase();
  const ti=document.getElementById('radar-tier')?.value||'';
  const sa=document.getElementById('radar-safety')?.value||'';
  return D.rows.map((r,idx)=>({r,idx})).filter(({r})=>{
    const txt=[r.genome,r.full_species,r.strain,r.Candidate_tier_v1_1,r.Refined_safety_status,r.ani_species]
      .filter(Boolean).join(' ').toLowerCase();
    return (!q||txt.includes(q)) && (!ti||r.Candidate_tier_v1_1===ti) && (!sa||r.Refined_safety_status===sa);
  }).sort((a,b)=>(b.r.LAB_score_v1||0)-(a.r.LAB_score_v1||0));
}

function renderRadarSelector(){
  const list=document.getElementById('radar-list');
  if(!list) return;
  const rows=radarFilteredRows();
  document.getElementById('radar-count').textContent = `${rows.length} matching strains • ${radarSelected.size} selected`;
  list.innerHTML = rows.slice(0,600).map(({r,idx})=>{
    const checked = radarSelected.has(idx) ? 'checked' : '';
    const sel = radarSelected.has(idx) ? 'sel' : '';
    const score = r.LAB_score_v1!=null ? Number(r.LAB_score_v1).toFixed(1) : '—';
    return `<div class="radar-item ${sel}" onclick="toggleRadarSelect(${idx})">
      <input type="checkbox" ${checked} onclick="event.stopPropagation(); toggleRadarSelect(${idx})">
      <div><div class="radar-title">${strainLabel(r,72)}</div>
      <div class="radar-sub">${r.genome||'—'} • Score ${score} • ${r.Candidate_tier_v1_1||'—'} • ${r.Refined_safety_status||'—'}</div></div>
      <button class="radar-view" onclick="event.stopPropagation(); openModal(${idx})">View</button>
    </div>`;
  }).join('') + (rows.length>600?'<div class="radar-sub" style="padding:10px">Showing first 600 matches. Refine search to narrow results.</div>':'');
  renderRadarChips();
}

function toggleRadarSelect(idx){
  if(radarSelected.has(idx)) radarSelected.delete(idx);
  else {
    if(radarSelected.size >= RADAR_MAX_SELECTED){ alert(`Please select no more than ${RADAR_MAX_SELECTED} strains for a readable radar plot.`); return; }
    radarSelected.add(idx);
  }
  renderRadarSelector();
  renderRadar();
}

function radarPreset(kind, silent=false){
  radarSelected.clear();
  let rows=D.rows.map((r,idx)=>({r,idx}));
  if(kind==='elite') rows=rows.filter(x=>String(x.r.Candidate_tier_v1_1||'').includes('Elite'));
  rows=rows.filter(x=>x.r.Refined_safety_status!=='Critical').sort((a,b)=>(b.r.LAB_score_v1||0)-(a.r.LAB_score_v1||0));
  const n = kind==='top10' ? 10 : kind==='elite' ? Math.min(10, rows.length) : 5;
  rows.slice(0,n).forEach(x=>radarSelected.add(x.idx));
  if(!silent){renderRadarSelector();renderRadar();}
}
function clearRadarSelection(){radarSelected.clear();renderRadarSelector();renderRadar();}
function selectedRadarRows(){ return [...radarSelected].map(idx=>({idx,r:D.rows[idx]})).filter(x=>x.r); }

function renderRadarChips(){
  const el=document.getElementById('radar-chips'); if(!el) return;
  const rows=selectedRadarRows();
  el.innerHTML = rows.length ? rows.map(({r,idx})=>`<span class="radar-chip" onclick="openModal(${idx})">${strainLabel(r,42)} ✦</span>`).join('') : '<span class="radar-sub">No strain selected.</span>';
}

function normalizedPanelValues(r){
  return PC.map(pc=>{
    const maxv = PM[pc] || 1;
    const v = Number(r[pc] || 0);
    return Math.max(0, Math.min(100, (v / maxv) * 100));
  });
}

function renderRadar(){
  const div=document.getElementById('ch-radar-compare'); if(!div) return;
  const rows=selectedRadarRows();
  renderRadarChips();
  renderRadarSelectedTable(rows);
  if(!rows.length){
    Plotly.purge(div); div.innerHTML='<div class="radar-note">Select one or more strains on the left to draw radar profiles.</div>'; return;
  }
  const traces=rows.map(({r,idx},i)=>{
    const vals=normalizedPanelValues(r);
    const theta=[...PL,PL[0]];
    return {type:'scatterpolar',mode:'lines+markers',name:strainLabel(r,38),r:[...vals,vals[0]],theta,
      fill:'toself',opacity:0.78,line:{color:radarPalette[i%radarPalette.length],width:2.4},marker:{size:5},
      customdata:Array(theta.length).fill(idx),
      hovertemplate:'%{fullData.name}<br>%{theta}: %{r:.1f}% of max<extra></extra>'};
  });
  Plotly.newPlot(div,traces,{polar:{radialaxis:{visible:true,range:[0,100],ticksuffix:'%',gridcolor:'#ddd'},angularaxis:{tickfont:{size:10},rotation:90}},
    margin:{t:30,b:40,l:70,r:70},paper_bgcolor:'transparent',plot_bgcolor:'transparent',legend:{orientation:'h',y:-0.18},
    title:{text:`Functional module radar (${rows.length} selected strain${rows.length>1?'s':''})`,font:{size:13,color:'#1B4332'}}},
    {responsive:true,displaylogo:false,modeBarButtonsToAdd:[{name:'Download comparison TSV',icon:Plotly.Icons.disk,click:exportRadarTable}]});
  div.on('plotly_click', function(ev){
    const idx = ev?.points?.[0]?.customdata;
    if(idx!=null) openModal(idx);
  });
}

function renderRadarSelectedTable(rows){
  const el=document.getElementById('radar-selected-table'); if(!el) return;
  if(!rows.length){el.innerHTML=''; return;}
  const cols=['genome','full_species','strain','LAB_score_v1','Candidate_tier_v1_1','Refined_safety_status','Safety_score','GI_survival_score','Functional_score','Fermentation_score'];
  el.innerHTML=`<table><thead><tr>${cols.map(c=>`<th>${c}</th>`).join('')}</tr></thead><tbody>${rows.map(({r,idx})=>`<tr onclick="openModal(${idx})">${cols.map(c=>`<td>${r[c]??''}</td>`).join('')}</tr>`).join('')}</tbody></table>`;
}

function exportRadarPNG(){ Plotly.downloadImage('ch-radar-compare',{format:'png',filename:'LAB_score_radar_compare',width:1400,height:1000,scale:2}); }
function exportRadarSVG(){ Plotly.downloadImage('ch-radar-compare',{format:'svg',filename:'LAB_score_radar_compare',width:1400,height:1000}); }
function exportRadarTable(){
  const rows=selectedRadarRows(); if(!rows.length) return alert('No strain selected.');
  const cols=['genome','accession','full_species','genus','species','strain','LAB_score_v1','Priority_class','Candidate_tier_v1_1','Refined_safety_status','Safety_score','GI_survival_score','Functional_score','Fermentation_score','amrfinder_hits','vfdb_hits','resfinder_hits','biogenic_amines','hemolysin_markers',...PC];
  const tsv=[cols.join('\\t'),...rows.map(({r})=>cols.map(c=>r[c]==null?'':String(r[c]).replace(/\\t/g,' ')).join('\\t'))].join('\\n');
  downloadFile('LAB_score_radar_selected_strains.tsv',tsv,'text/tab-separated-values');
}

function htmlEscape(x){ return String(x??'').replace(/[&<>\"]/g, s=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[s])); }
function makeStrainReportHTML(rows){
  const cols=['genome','accession','full_species','strain','LAB_score_v1','Priority_class','Candidate_tier_v1_1','Refined_safety_status','Safety_score','GI_survival_score','Functional_score','Fermentation_score','amrfinder_hits','vfdb_hits','resfinder_hits','biogenic_amines','hemolysin_markers'];
  const panelHeader=PC.map((c,i)=>`<th>${htmlEscape(PL[i])}</th>`).join('');
  const cards=rows.map(({r})=>`<section><h2>${htmlEscape(strainLabel(r,120))}</h2><table>${cols.map(c=>`<tr><th>${htmlEscape(c)}</th><td>${htmlEscape(r[c])}</td></tr>`).join('')}</table><h3>Functional modules</h3><table><tr>${panelHeader}</tr><tr>${PC.map(c=>`<td>${htmlEscape(r[c]??0)}</td>`).join('')}</tr></table></section>`).join('\\n');
  return `<!DOCTYPE html><html><head><meta charset="UTF-8"><title>LAB-Score selected strain report</title><style>body{font-family:Segoe UI,Arial,sans-serif;margin:28px;color:#1a2e1a}h1{color:#1B4332}h2{margin-top:24px;color:#2D6A4F}table{border-collapse:collapse;width:100%;margin:8px 0 16px;font-size:13px}th,td{border:1px solid #cfe3cf;padding:6px 8px;text-align:left}th{background:#e8f5e9}.note{color:#4a6a4a;font-size:13px}</style></head><body><h1>LAB-Score selected strain report</h1><p class="note">Generated from the interactive LAB-Score report. Selected strains: ${rows.length}</p>${cards}</body></html>`;
}
function exportRadarReport(){
  const rows=selectedRadarRows(); if(!rows.length) return alert('No strain selected.');
  downloadFile('LAB_score_selected_strain_report.html', makeStrainReportHTML(rows),'text/html');
}

// ── MODAL ─────────────────────────────────────────────────────────────────────
let currentRow=null;
function openModal(idx){
  selIdx=idx;
  currentRow=D.rows[idx];
  const r=currentRow;
  renderTbody();

  const sp=r.full_species&&r.full_species!=='Unknown sp.'&&r.full_species!=='nan'?r.full_species:r.genome;
  document.getElementById('modal-title').textContent=`📋 ${sp}`;

  let alerts='';
  if(r.Refined_safety_status==='Critical')
    alerts+=`<div class="alert ared">🚨 CRITICAL — ${r.reasons||'AMR/Virulence detected'}</div>`;
  else if(r.Refined_safety_status==='Caution')
    alerts+=`<div class="alert aamb">⚠️ CAUTION — ${r.reasons||'Hemolysin/biogenic amine marker'}</div>`;
  else
    alerts+=`<div class="alert agrn">✅ No critical safety markers detected</div>`;

  if((r.amrfinder_hits||0)>0) alerts+=`<div class="alert ared">🦠 AMR genes: ${r.amrfinder_hits} hit(s)</div>`;
  if((r.vfdb_hits||0)>0)      alerts+=`<div class="alert ared">⚡ Virulence factors: ${r.vfdb_hits} hit(s)</div>`;
  if((r.resfinder_hits||0)>0) alerts+=`<div class="alert aamb">💊 ResFinder: ${r.resfinder_hits} hit(s)</div>`;
  if((r.biogenic_amines||0)>0)alerts+=`<div class="alert aamb">⚗️ Biogenic amine genes: ${r.biogenic_amines}</div>`;
  if((r.hemolysin_markers||0)>0)alerts+=`<div class="alert aamb">🩸 Hemolysin markers: ${r.hemolysin_markers}</div>`;

  let spAlert='';
  const spStatus=r.species_status||'—';
  const spCol=spStatus==='CONFIRMED'?'#1B4332':spStatus==='MISMATCH'?'#C1121F':spStatus==='GENUS_MATCH'?'#E07C24':'#888';
  if(spStatus==='MISMATCH')
    spAlert=`<div class="alert ared">⚠️ Species mismatch: NCBI="${r.full_species}" vs ANI="${r.ani_species}"</div>`;
  else if(spStatus==='GENUS_MATCH')
    spAlert=`<div class="alert aamb">ℹ️ Same genus, different species-level match by ANI</div>`;

  const row=([k,v])=>`<div class="mrow"><span class="mk">${k}</span><span class="mv">${v??'—'}</span></div>`;

  const scoreSec=`<div class="msec">
    <h3>🏆 LAB-Score Breakdown</h3>
    ${[['LAB-Score v1.1',(r.LAB_score_v1||0).toFixed(2)],['Priority Class',r.Priority_class],
      ['Candidate Tier',r.Candidate_tier_v1_1],['Safety Status',r.Refined_safety_status],
      ['Safety Score (×0.45)',(r.Safety_score||0).toFixed(2)],
      ['GI Survival Score (×0.25)',(r.GI_survival_score||0).toFixed(2)],
      ['Functional Score (×0.20)',(r.Functional_score||0).toFixed(2)],
      ['Fermentation Score (×0.10)',(r.Fermentation_score||0).toFixed(2)],
    ].map(row).join('')}
  </div>`;

  const safeSec=`<div class="msec">
    <h3>🛡️ Safety Analysis</h3>
    ${alerts}
    ${[['AMRFinder hits',r.amrfinder_hits||0],['Virulence factors',r.vfdb_hits||0],
      ['ResFinder hits',r.resfinder_hits||0],['ResFinder genes',r.resfinder_genes||'None'],
      ['Biogenic amines',r.biogenic_amines||0],['Hemolysin markers',r.hemolysin_markers||0],
      ['Safety reasons',r.reasons||'None'],
    ].map(row).join('')}
  </div>`;

  const spSec=`<div class="msec">
    <h3>🦠 Species Identification</h3>
    ${spAlert}
    ${[['Genome',r.genome],['NCBI Species',r.full_species||'Unknown'],
      ['Genus',r.genus||'—'],['Strain',r.strain||'—'],
      ['ANI Species',r.ani_species||'—'],
      ['ANI Identity',r.ani_identity!=null?r.ani_identity+'%':'—'],
      ['Species Status',`<span style="color:${spCol};font-weight:700">${spStatus}</span>`],
      ['Assembly Level',r.assembly_level||'—'],
    ].map(row).join('')}
  </div>`;

  const qcSec=`<div class="msec">
    <h3>🔬 Genome Quality</h3>
    ${[['Completeness',r.Completeness!=null?r.Completeness+'%':r.completeness!=null?r.completeness+'%':'—'],
      ['Contamination',r.Contamination!=null?r.Contamination+'%':r.contamination!=null?r.contamination+'%':'—'],
      ['Total Length (bp)',r.TotalLength||r.total_length||'—'],
      ['Contigs',r.Contigs||r.contigs||'—'],
      ['GC%',r['GC%']||r.gc_pct||'—'],
    ].map(row).join('')}
  </div>`;

  const cazSec=`<div class="msec">
    <h3>🍄 CAZyme &amp; Fermentation</h3>
    ${[['Total CAZymes',r.CAZyme_total||0],
      ['Carbohydrate markers',r['carbohydrate.carbohydrate']||0],
      ['Metabolism markers',r['metabolism.metabolism']||0],
      ['EPS/Cell env.',r['cellenvelope_eps.cellenvelope_eps']||0],
      ['Bacteriocins',r['bacteriocins.bacteriocin']||r['bacteriocins.bacteriocins']||0],
    ].map(row).join('')}
  </div>`;

  const maxVals=PC.map(pc=>PM[pc]||1);
  const panelBars=PC.map((pc,i)=>{
    const v=r[pc]||0;
    const pct=Math.min(100,(v/maxVals[i])*100);
    const col=pct>66?'#1B4332':pct>33?'#52B788':'#95D5B2';
    return `<div class="pbar-row">
      <span class="pbar-lbl">${PL[i]}</span>
      <div class="pbar-out"><div class="pbar-in" style="width:${pct}%;background:${col}"></div></div>
      <span class="pbar-v">${v}</span>
    </div>`;
  }).join('');

  const panelSec=`<div class="msec mfull">
    <h3>🧬 Functional Module Marker Counts</h3>
    ${panelBars}
  </div>`;

  const radarSec=`<div class="msec mfull">
    <h3>🎯 Functional Profile Radar</h3>
    <div id="modal-radar" style="height:420px"></div>
  </div>`;

  document.getElementById('modal-body').innerHTML =
    scoreSec + safeSec + spSec + qcSec + cazSec + panelSec + radarSec;

  // Radar
  setTimeout(()=>{
    const vals=PC.map((pc,i)=>((r[pc]||0)/maxVals[i])*100);
    Plotly.newPlot('modal-radar',[{
      type:'scatterpolar', r:[...vals,vals[0]], theta:[...PL,PL[0]],
      fill:'toself', fillcolor:'rgba(82,183,136,0.2)',
      line:{color:'#2D6A4F',width:2.5},
      hovertemplate:'%{theta}<br>%{r:.1f}% of max<extra></extra>'
    }],{
      polar:{radialaxis:{visible:true,range:[0,100],ticksuffix:'%',tickfont:{size:9},gridcolor:'#ddd'},
              angularaxis:{tickfont:{size:10},rotation:90}},
      margin:{t:30,b:30,l:60,r:60},paper_bgcolor:'transparent',showlegend:false,
      title:{text:'Functional profile (% of dataset maximum)',font:{size:11,color:'#1B4332'}}
    },{responsive:true,displayModeBar:false});
  },120);

  document.getElementById('overlay').classList.add('on');
  document.body.style.overflow='hidden';
}

function closeModal(e){
  if(e&&e.target!==document.getElementById('overlay')&&!e.target.classList.contains('mclose')) return;
  selIdx=-1;
  document.getElementById('overlay').classList.remove('on');
  document.body.style.overflow='';
  renderTbody();
}

// ── EXPORT ───────────────────────────────────────────────────────────────────
function exportStrain(){
  if(!currentRow) return;
  const r=currentRow;
  const keys=Object.keys(r);
  const vals=keys.map(k=>{const v=r[k];return v==null?'':String(v)});
  const tsv=keys.join('\\t')+'\\n'+vals.join('\\t');
  const sp=(r.full_species&&r.full_species!=='nan'?r.full_species:r.genome).replace(/[^a-zA-Z0-9_]/g,'_');
  downloadFile(`LAB_score_${sp}.tsv`, tsv, 'text/tab-separated-values');
}

function exportStrainJSON(){
  if(!currentRow) return;
  const r=currentRow;
  const sp=(r.full_species&&r.full_species!=='nan'?r.full_species:r.genome).replace(/[^a-zA-Z0-9_]/g,'_');
  downloadFile(`LAB_score_${sp}.json`, JSON.stringify(r, null, 2), 'application/json');
}

function exportCurrentRadarPNG(){
  if(!currentRow) return;
  const base=(currentRow.full_species&&currentRow.full_species!=='nan'?currentRow.full_species:currentRow.genome).replace(/[^a-zA-Z0-9_]/g,'_');
  const div=document.getElementById('modal-radar');
  if(!div) return;
  Plotly.downloadImage('modal-radar',{format:'png',filename:`LAB_score_${base}_radar`,width:1200,height:900,scale:2});
}

function exportStrainReportHTML(){
  if(!currentRow) return;
  const idx=D.rows.indexOf(currentRow);
  const base=(currentRow.full_species&&currentRow.full_species!=='nan'?currentRow.full_species:currentRow.genome).replace(/[^a-zA-Z0-9_]/g,'_');
  downloadFile(`LAB_score_${base}_report.html`, makeStrainReportHTML([{idx:idx,r:currentRow}]), 'text/html');
}

function downloadFile(filename, content, mime){
  const a=document.createElement('a');
  a.href=URL.createObjectURL(new Blob([content],{type:mime}));
  a.download=filename;
  a.click();
  URL.revokeObjectURL(a.href);
}

// ── Heatmap ───────────────────────────────────────────────────────────────────
function renderHeatmap(){
  if(!document.getElementById('page-heatmap').classList.contains('on')) return;
  const mode=document.getElementById('hmode').value;
  const sort=document.getElementById('hsort').value;
  const nv=parseInt(document.getElementById('hn').value)||0;
  let rows=[...D.rows];
  if(sort==='score') rows.sort((a,b)=>(b.LAB_score_v1||0)-(a.LAB_score_v1||0));
  else if(sort==='tier') rows.sort((a,b)=>String(a.Candidate_tier_v1_1||'').localeCompare(String(b.Candidate_tier_v1_1||'')));
  else rows.sort((a,b)=>String(a.full_species||'').localeCompare(String(b.full_species||'')));
  if(nv>0) rows=rows.slice(0,nv);
  const xlabels=rows.map(r=>{
    const s=r.full_species&&r.full_species!=='Unknown sp.'&&r.full_species!=='nan'?r.full_species:r.genome;
    return s.substring(0,28);
  });
  let z=PC.map(pc=>rows.map(r=>r[pc]||0));
  if(mode==='presence'){
    const thr=PC.map(pc=>{const vs=D.rows.map(r=>r[pc]||0);return Math.max(1,vs.reduce((a,b)=>a+b,0)/vs.length*.5)});
    z=z.map((row,i)=>row.map(v=>v>=thr[i]?1:0));
  } else if(mode==='zscore'){
    z=z.map(row=>{const m=row.reduce((a,b)=>a+b,0)/row.length;
      const s=Math.sqrt(row.reduce((a,b)=>a+(b-m)**2,0)/row.length)||1;
      return row.map(v=>(v-m)/s);});
  }
  const cscale=mode==='presence'?[[0,'#f5f5f5'],[1,'#1B4332']]:
    mode==='zscore'?[[0,'#C1121F'],[.5,'#f8faf8'],[1,'#1B4332']]:
    [[0,'#f5f5f5'],[.4,'#95D5B2'],[.7,'#2D6A4F'],[1,'#1B4332']];
  Plotly.newPlot('ch-hmap',[{type:'heatmap',z,x:xlabels,y:PL,colorscale:cscale,
    hovertemplate:'Strain: %{x}<br>Panel: %{y}<br>Value: %{z}<extra></extra>'}],
    {margin:{t:10,b:130,l:165,r:60},paper_bgcolor:'transparent',
      xaxis:{tickangle:-50,tickfont:{size:8},automargin:true},
      yaxis:{automargin:true,tickfont:{size:10}}},{responsive:true});
  function catChart(id,pnames){
    const cols=pnames.map(p=>`${p}.${p}`).filter(c=>PC.includes(c));
    const lbs=cols.map(c=>PL[PC.indexOf(c)]);
    const ms=cols.map(c=>rows.reduce((a,r)=>a+(r[c]||0),0)/rows.length);
    Plotly.newPlot(id,[{x:lbs,y:ms,type:'bar',
      marker:{color:'#52B788',line:{color:'#fff',width:1}},
      hovertemplate:'%{x}<br>Mean: %{y:.2f}<extra></extra>'}],
      {margin:{t:10,b:80,l:40,r:10},paper_bgcolor:'transparent',plot_bgcolor:'transparent',
        xaxis:{tickangle:-35,automargin:true,tickfont:{size:8.5}},
        yaxis:{title:'Mean count',gridcolor:'#eee'}},{responsive:true,displayModeBar:false});
  }
  catChart('ch-gi',['acid_energy','bileresistance','gutpersistence','osmoticstress','heatstress','coldstress','gaba']);
  catChart('ch-fu',['vitamins','immunomodulation','bacteriocins','cellenvelope_eps','adhesion_surface','adhesion_biofilm','defense_crispr']);
  catChart('ch-fe',['carbohydrate','metabolism','antipath_qs','alkalinestress']);
}

// ── Init ──────────────────────────────────────────────────────────────────────
renderOverview();
renderAbout();
initTable();
initRadar();
</script>
</body>
</html>"""
HTML = HTML_TEMPLATE.replace('{PLOTLY_JS}', _plotly_js).replace('{json_str}', json_str)

out = os.path.join(args.outdir, "LAB_score_report.html")
with open(out, "w", encoding="utf-8") as f:
    f.write(HTML)
print(f"[report] Saved: {out}  ({os.path.getsize(out)//1024} KB)")
