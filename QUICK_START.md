# LAB-Score Pipeline v3.4 — Quick Start

## Install
```bash
bash install.sh \
  --env-name lab_score_clean \
  --dbcan-env-name dbcan_clean \
  --dbcan-db-dir /path/to/db/dbcan \
  --skip-optional
```

## Run (from FASTA)
```bash
conda activate lab_score_clean
bash run_all.sh \
  -i genomes/ \
  -o results/ \
  -t 32 \
  --dbcan-db-dir /path/to/db/dbcan
```

## Run (from existing Prokka)
```bash
bash run_all.sh \
  -p annotations_prokka_GCA \
  -o results/ \
  -t 32 \
  --dbcan-db-dir /path/to/db/dbcan
```

## Species verification options
```bash
# Auto-download LAB reference genomes (first run)
bash run_all.sh -i genomes/ -o results/ -t 32 --download-refs

# Use your own reference genomes
bash run_all.sh -i genomes/ -o results/ -t 32 --ref-dir /path/to/lab_refs/
```

## Pipeline steps
| Step | Description |
|------|-------------|
| 0    | Species metadata (NCBI API auto-lookup) |
| 0b   | **Species verification (FastANI)** — flags misidentified genomes |
| 1    | Prokka annotation |
| 2    | AMRFinder (AMR + virulence) |
| 3    | ResFinder (acquired resistance) |
| 4    | dbCAN (CAZymes) |
| 5    | CheckM QC (optional, -q) |
| 6    | Safety screening |
| 7    | Functional panel scan (18 panels) |
| 8    | Master matrix |
| 9    | LAB-Score v1.1 calculation |
| 10   | ML interpretation |
| 11   | Publication figures |

## Species status codes (Step 0b output)
| Status | Meaning |
|--------|---------|
| CONFIRMED | NCBI name matches ANI top hit (≥95%) ✅ |
| MISMATCH | NCBI name differs from genome-based ID ⚠️ |
| GENUS_MATCH | Same genus, different species |
| LOW_ANI | <95% ANI — possible novel strain |
| NO_NCBI | No NCBI name — ANI provides ID |
| NO_FASTANI | FastANI not installed (skipped) |

## Check installation
```bash
bash run_all.sh --check-only \
  --env lab_score_clean \
  --dbcan-env-name dbcan_clean \
  --dbcan-db-dir /path/to/db/dbcan
```
