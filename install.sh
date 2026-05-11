#!/usr/bin/env bash
# =============================================================================
# LAB-Score Pipeline v3.3 — prerequisite installer/checker
# =============================================================================
# Safer installer: keeps dbCAN in a separate clean conda environment to avoid
# mixed pip/conda collisions such as numba ClobberError in the main environment.
#
# Recommended:
#   bash install.sh --env-name lab_score_clean --dbcan-env-name dbcan_clean \
#     --dbcan-db-dir /media/mecob/komwit/db/dbcan --skip-optional
# =============================================================================

set -euo pipefail

RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; BLUE='\033[0;34m'; NC='\033[0m'
ok()   { echo -e "${GREEN}[✓]${NC} $1"; }
warn() { echo -e "${YELLOW}[!]${NC} $1"; }
info() { echo -e "${BLUE}[i]${NC} $1"; }
fail() { echo -e "${RED}[✗]${NC} $1"; }
header() { echo -e "\n${BLUE}══════════════════════════════════════════${NC}\n${BLUE}  $1${NC}\n${BLUE}══════════════════════════════════════════${NC}"; }

CHECK_ONLY=false
ENV_NAME="lab_score_clean"
DBCAN_ENV_NAME="dbcan_clean"
DBCAN_DB_DIR_USER=""
SKIP_DBCAN_DB=false
USE_AWS_S3=true
SKIP_OPTIONAL=false
RECREATE_ENV=false
RECREATE_DBCAN_ENV=false

usage() {
  cat <<EOF
LAB-Score Pipeline v3.3 installer

Options:
  --check                       Check tools/databases only; do not install
  -e, --env-name NAME            Main conda env [default: lab_score_clean]
  --dbcan-env-name NAME          Separate dbCAN env [default: dbcan_clean]
  --dbcan-db-dir PATH            dbCAN database directory [default: <dbcan_env_prefix>/db]
  --skip-dbcan-db                Do not download/update dbCAN database
  --no-aws-s3                    Do not use dbCAN AWS S3 download option
  --skip-optional                Skip optional ResFinder/CheckM/SHAP installs
  --recreate-env                 Delete and recreate the main env first
  --recreate-dbcan-env           Delete and recreate the dbCAN env first
  --force-dbcan-reinstall        Alias for --recreate-dbcan-env
  -h, --help                     Show this help

Notes:
  - The main env runs Prokka, AMRFinderPlus, scoring, ML, and figures.
  - The dbCAN env runs run_dbcan only. This avoids dirty main-env collisions.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --check) CHECK_ONLY=true; shift ;;
    -e|--env-name) ENV_NAME="${2:?Missing env name}"; shift 2 ;;
    --dbcan-env-name) DBCAN_ENV_NAME="${2:?Missing dbCAN env name}"; shift 2 ;;
    --dbcan-db-dir) DBCAN_DB_DIR_USER="${2:?Missing dbCAN db path}"; shift 2 ;;
    --skip-dbcan-db) SKIP_DBCAN_DB=true; shift ;;
    --no-aws-s3) USE_AWS_S3=false; shift ;;
    --skip-optional) SKIP_OPTIONAL=true; shift ;;
    --recreate-env) RECREATE_ENV=true; shift ;;
    --recreate-dbcan-env|--force-dbcan-reinstall) RECREATE_DBCAN_ENV=true; shift ;;
    -h|--help) usage; exit 0 ;;
    *) fail "Unknown option: $1"; usage; exit 1 ;;
  esac
done

header "LAB-Score v3.3 prerequisite setup"

if ! command -v conda >/dev/null 2>&1; then
  fail "Conda was not found. Please install Miniconda/Anaconda first."
  exit 1
fi
ok "Conda: $(conda --version)"

CONDA_BASE="$(conda info --base)"
# shellcheck disable=SC1090
source "$CONDA_BASE/etc/profile.d/conda.sh" 2>/dev/null || true

env_exists() {
  local env="$1"
  conda env list | awk '{print $1}' | grep -Fxq "$env"
}

env_prefix() {
  local env="$1"
  conda env list | awk -v env="$env" '$1==env {print $NF; exit}'
}

remove_env_if_requested() {
  local env="$1" recreate="$2"
  if [[ "$recreate" == true && "$(env_exists "$env" && echo yes || echo no)" == yes ]]; then
    if [[ "$CHECK_ONLY" == true ]]; then
      warn "Would recreate env but --check was used: $env"
    else
      warn "Removing existing conda env: $env"
      conda env remove -n "$env" -y
    fi
  fi
}

create_env_if_missing() {
  local env="$1" python_version="${2:-3.10}"
  if ! env_exists "$env"; then
    if [[ "$CHECK_ONLY" == true ]]; then
      fail "Conda environment not found: $env"
      return 1
    fi
    header "Creating conda environment: $env"
    conda create -n "$env" -y -c conda-forge "python=${python_version}"
  fi
}

run_in_env() {
  local env="$1"; shift
  conda run -n "$env" --no-capture-output "$@"
}

# Avoid leakage from another Miniconda/base Perl into Prokka.
# This was observed as XML/Simple.pm lookup in the wrong @INC path.
clean_perl_exports() {
  unset PERL5LIB PERL_LOCAL_LIB_ROOT PERL_MB_OPT PERL_MM_OPT
}

cmd_in_env() {
  local env="$1" cmd="$2"
  local prefix
  prefix="$(env_prefix "$env")"
  [[ -n "$prefix" && -x "$prefix/bin/$cmd" ]]
}

python_import_in_env() {
  local env="$1" module="$2"
  conda run -n "$env" python -c "import ${module}" >/dev/null 2>&1
}

check_tool() {
  local env="$1" label="$2" cmd="$3" optional="${4:-required}"
  local prefix
  prefix="$(env_prefix "$env")"
  if [[ -z "$prefix" || ! -x "$prefix/bin/$cmd" ]]; then
    if [[ "$optional" == "optional" ]]; then
      warn "$label [$env]: not found (optional)"
      return 0
    else
      fail "$label [$env]: NOT FOUND"
      return 1
    fi
  fi

  # Run real usability checks, not only executable existence.
  # Some tools do not support --version; Prokka can be broken by leaked Perl env vars.
  local check_cmd
  case "$cmd" in
    python)
      check_cmd='python --version 2>&1 | head -1'
      ;;
    prokka)
      check_cmd='unset PERL5LIB PERL_LOCAL_LIB_ROOT PERL_MB_OPT PERL_MM_OPT; prokka --version 2>&1 | head -1'
      ;;
    hmmscan)
      check_cmd='hmmscan -h 2>&1 | head -1'
      ;;
    hmmsearch)
      check_cmd='hmmsearch -h 2>&1 | head -1'
      ;;
    diamond)
      check_cmd='diamond version 2>&1 | head -1'
      ;;
    amrfinder)
      check_cmd='amrfinder --version 2>&1 | head -1'
      ;;
    *)
      check_cmd="$cmd --version 2>&1 | head -1"
      ;;
  esac

  local out
  if out=$(conda run -n "$env" --no-capture-output bash -c "$check_cmd" 2>&1); then
    out=$(echo "$out" | head -1)
    [[ -z "$out" ]] && out="usable"
    ok "$label [$env]: $out"
    return 0
  fi

  if [[ "$optional" == "optional" ]]; then
    warn "$label [$env]: found but failed usability check (optional)"
    echo "$out" | head -5
    return 0
  else
    fail "$label [$env]: found but FAILED usability check"
    echo "$out" | head -12
    return 1
  fi
}

check_dbcan_cli() {
  local env="$DBCAN_ENV_NAME" prefix
  prefix="$(env_prefix "$env")"
  [[ -n "$prefix" ]] || return 1
  [[ -x "$prefix/bin/python" ]] || return 1
  [[ -x "$prefix/bin/diamond" ]] || return 1

  # Check if run_dbcan wrapper exists and is executable
  if [[ ! -x "$prefix/bin/run_dbcan" ]]; then
    return 1
  fi

  # Try multiple dbcan module import paths (GitHub vs PyPI vs bioconda versions differ)
  local can_import=false
  for mod in     "import dbcan_cli; import dbcan_cli.run_dbcan"     "from dbcan.main import cli"     "from dbcan import cli"     "import dbcan"
  do
    if conda run -n "$env" --no-capture-output "$prefix/bin/python" -c "$mod" >/dev/null 2>&1; then
      can_import=true
      break
    fi
  done
  [[ "$can_import" == true ]] || return 1

  # Test the wrapper itself
  conda run -n "$env" --no-capture-output "$prefix/bin/run_dbcan" --help >/dev/null 2>&1
}

ensure_run_dbcan_wrapper() {
  local env="$DBCAN_ENV_NAME" prefix
  prefix="$(env_prefix "$env")"
  [[ -n "$prefix" ]] || return 1

  # Detect which dbcan module/entry point is available and build correct wrapper
  local import_line="" entry_line=""

  if conda run -n "$env" --no-capture-output "$prefix/bin/python" -c       "import dbcan_cli; import dbcan_cli.run_dbcan" >/dev/null 2>&1; then
    import_line="from dbcan_cli.run_dbcan import cli_main"
    entry_line="sys.exit(cli_main())"
  elif conda run -n "$env" --no-capture-output "$prefix/bin/python" -c       "from dbcan.main import cli" >/dev/null 2>&1; then
    import_line="from dbcan.main import cli"
    entry_line="sys.exit(cli())"
  elif conda run -n "$env" --no-capture-output "$prefix/bin/python" -c       "from dbcan import cli" >/dev/null 2>&1; then
    import_line="from dbcan import cli"
    entry_line="sys.exit(cli())"
  else
    warn "Cannot detect dbcan entry point — trying generic wrapper"
    import_line="import dbcan"
    entry_line="sys.exit(0)"
  fi

  info "Building run_dbcan wrapper using: $import_line"
  cat > "$prefix/bin/run_dbcan" <<EOF_RUNDBCAN
#!$prefix/bin/python
import re
import sys
$import_line
if __name__ == '__main__':
    sys.argv[0] = re.sub(r'(-script\.pyw|\.exe)?$', '', sys.argv[0])
    $entry_line
EOF_RUNDBCAN
  chmod +x "$prefix/bin/run_dbcan"
}

check_dbcan_db() {
  local dmnd_count hmm_count
  dmnd_count=$(find "$DB_DIR" -type f -name "*.dmnd" 2>/dev/null | wc -l || true)
  hmm_count=$(find "$DB_DIR" -type f \( -name "*.hmm" -o -name "*HMM*" -o -name "*hmm*" \) 2>/dev/null | wc -l || true)
  if [[ "$dmnd_count" -gt 0 && "$hmm_count" -gt 0 ]]; then
    ok "dbCAN database: present ($dmnd_count DIAMOND db, $hmm_count HMM-like files)"
    return 0
  fi
  warn "dbCAN database appears incomplete or missing in: $DB_DIR"
  return 1
}

remove_env_if_requested "$ENV_NAME" "$RECREATE_ENV"
remove_env_if_requested "$DBCAN_ENV_NAME" "$RECREATE_DBCAN_ENV"

create_env_if_missing "$ENV_NAME" 3.10 || exit 1
create_env_if_missing "$DBCAN_ENV_NAME" 3.10 || exit 1

ENV_PREFIX="$(env_prefix "$ENV_NAME")"
DBCAN_ENV_PREFIX="$(env_prefix "$DBCAN_ENV_NAME")"
if [[ -z "$ENV_PREFIX" || -z "$DBCAN_ENV_PREFIX" ]]; then
  fail "Cannot resolve conda environment paths."
  exit 1
fi
DB_DIR="${DBCAN_DB_DIR_USER:-$DBCAN_ENV_PREFIX/db}"

info "Main conda env    : $ENV_NAME"
info "Main env path     : $ENV_PREFIX"
info "dbCAN conda env   : $DBCAN_ENV_NAME"
info "dbCAN env path    : $DBCAN_ENV_PREFIX"
info "dbCAN database    : $DB_DIR"

if [[ "$CHECK_ONLY" == true ]]; then
  header "Checking tools"
  FAILS=0
  check_tool "$ENV_NAME" "Python" python || ((FAILS++))
  python_import_in_env "$ENV_NAME" pandas && ok "pandas [$ENV_NAME]" || { fail "pandas [$ENV_NAME]: NOT FOUND"; ((FAILS++)); }
  python_import_in_env "$ENV_NAME" numpy && ok "numpy [$ENV_NAME]" || { fail "numpy [$ENV_NAME]: NOT FOUND"; ((FAILS++)); }
  python_import_in_env "$ENV_NAME" sklearn && ok "scikit-learn [$ENV_NAME]" || { fail "scikit-learn [$ENV_NAME]: NOT FOUND"; ((FAILS++)); }
  python_import_in_env "$ENV_NAME" matplotlib && ok "matplotlib [$ENV_NAME]" || { fail "matplotlib [$ENV_NAME]: NOT FOUND"; ((FAILS++)); }
  check_tool "$ENV_NAME" "Prokka" prokka || ((FAILS++))
  check_tool "$ENV_NAME" "AMRFinderPlus" amrfinder || ((FAILS++))
  check_tool "$ENV_NAME" "FastANI" fastANI optional || true
  if check_dbcan_cli; then
    ok "run_dbCAN [$DBCAN_ENV_NAME]: usable"
  else
    fail "run_dbCAN [$DBCAN_ENV_NAME]: broken or NOT FOUND"
    ((FAILS++))
  fi
  check_tool "$DBCAN_ENV_NAME" "DIAMOND" diamond || ((FAILS++))
  check_tool "$DBCAN_ENV_NAME" "HMMER/hmmscan" hmmscan || ((FAILS++))
  check_tool "$ENV_NAME" "ResFinder" resfinder optional || true
  check_tool "$ENV_NAME" "CheckM" checkm optional || true
  python_import_in_env "$ENV_NAME" shap && ok "SHAP [$ENV_NAME]" || warn "SHAP [$ENV_NAME]: not found (optional)"
  check_dbcan_db || ((FAILS++))
  echo ""
  if [[ "$FAILS" -gt 0 ]]; then
    fail "$FAILS required checks failed."
    echo "Recommended clean install:"
    echo "  bash install.sh --env-name lab_score_clean --dbcan-env-name dbcan_clean --dbcan-db-dir '$DB_DIR' --recreate-dbcan-env --skip-optional"
    exit 1
  fi
  ok "All required prerequisites are ready."
  exit 0
fi

header "Installing main LAB-Score packages"
conda install -n "$ENV_NAME" -y -c conda-forge -c bioconda \
  pandas numpy scipy scikit-learn matplotlib seaborn biopython tqdm requests joblib \
  prokka ncbi-amrfinderplus blast prodigal perl-xml-simple perl-bioperl-core \
  || { fail "Main package installation failed. If this is an old/dirty env, rerun with --env-name lab_score_clean --recreate-env"; exit 1; }
ok "Main LAB-Score packages ready"

header "Updating AMRFinderPlus database"
run_in_env "$ENV_NAME" amrfinder --update || warn "AMRFinder database update failed; AMRFinder will retry during the pipeline."


header "Installing FastANI (species verification)"
if cmd_in_env "$ENV_NAME" fastANI; then
  ok "FastANI already installed"
else
  info "Installing FastANI via bioconda..."
  conda install -n "$ENV_NAME" -y -c bioconda -c conda-forge fastani \
    && ok "FastANI installed" \
    || warn "FastANI installation failed — species verification (Step 0b) will be skipped"
fi

header "Installing dbCAN environment"
# Install low-level binaries with conda in the separate env. If this env is dirty,
# --recreate-dbcan-env is the intended fix.
conda install -n "$DBCAN_ENV_NAME" -y -c conda-forge -c bioconda \
  python=3.10 diamond hmmer prodigal blast biopython requests tqdm pandas openpyxl pyhmmer pyrodigal \
  || { fail "dbCAN dependency installation failed. Rerun with --recreate-dbcan-env."; exit 1; }

# Remove stale/broken console wrappers and install run_dbcan in dbCAN env only.
rm -f "$DBCAN_ENV_PREFIX/bin/run_dbcan" "$DBCAN_ENV_PREFIX/bin/run_dbcan.py" 2>/dev/null || true
run_in_env "$DBCAN_ENV_NAME" python -m pip install --upgrade pip setuptools wheel || true
run_in_env "$DBCAN_ENV_NAME" python -m pip uninstall -y dbcan run-dbcan-new run-dbcan 2>/dev/null || true

warn "Installing run_dbcan from official GitHub source inside $DBCAN_ENV_NAME."
if run_in_env "$DBCAN_ENV_NAME" python -m pip install --no-cache-dir --force-reinstall "git+https://github.com/bcb-unl/run_dbcan.git"; then
  ok "Installed run_dbcan from GitHub"
else
  warn "GitHub installation failed; trying PyPI fallback package: dbcan"
  if run_in_env "$DBCAN_ENV_NAME" python -m pip install --no-cache-dir --force-reinstall dbcan; then
    ok "Installed dbcan from PyPI fallback"
  else
    warn "PyPI dbcan failed; trying run-dbcan-new fallback"
    run_in_env "$DBCAN_ENV_NAME" python -m pip install --no-cache-dir --force-reinstall run-dbcan-new \
      || { fail "run_dbcan Python package installation failed"; exit 1; }
  fi
fi

# Rebuild the wrapper using the absolute Python inside dbcan_clean, then check.
ensure_run_dbcan_wrapper || true
if ! check_dbcan_cli; then
  fail "run_dbcan is still not usable in $DBCAN_ENV_NAME."
  echo "Debug:"
  echo "  dbCAN env python: $DBCAN_ENV_PREFIX/bin/python"
  echo "  expected run_dbcan: $DBCAN_ENV_PREFIX/bin/run_dbcan"
  echo "  PATH-selected run_dbcan, if any: $(command -v run_dbcan || true)"
  echo "  Expected wrapper head:"
  head -30 "$DBCAN_ENV_PREFIX/bin/run_dbcan" 2>/dev/null || true
  echo "  Python module probe:"
  conda run -n "$DBCAN_ENV_NAME" --no-capture-output "$DBCAN_ENV_PREFIX/bin/python" - <<PY_DEBUG || true
import sys, pkgutil
print("Python:", sys.executable)
dbcan_mods = [m.name for m in pkgutil.iter_modules() if "dbcan" in m.name.lower() or "cazy" in m.name.lower()]
print("dbcan-related modules:", dbcan_mods)
for mod_path in ["dbcan_cli.run_dbcan", "dbcan.main", "dbcan"]:
    try:
        exec(f"import {mod_path.split('.')[0]}")
        print(f"OK: import {mod_path.split('.')[0]}")
    except Exception as e:
        print(f"FAIL: {mod_path}: {e}")
PY_DEBUG
  exit 1
fi
ok "run_dbCAN ready in $DBCAN_ENV_NAME"

mkdir -p "$DB_DIR"
if [[ "$SKIP_DBCAN_DB" == true ]]; then
  warn "Skipping dbCAN database download by request."
elif check_dbcan_db; then
  ok "Using existing dbCAN database: $DB_DIR"
else
  header "Downloading dbCAN database"
  info "This may take time and needs internet access. Directory: $DB_DIR"
  DB_ARGS=(database --db_dir "$DB_DIR")
  if [[ "$USE_AWS_S3" == true ]]; then
    DB_ARGS+=(--aws_s3)
  fi
  DB_ARGS+=(--no-cgc)

  if run_in_env "$DBCAN_ENV_NAME" "$DBCAN_ENV_PREFIX/bin/run_dbcan" "${DB_ARGS[@]}"; then
    ok "dbCAN database downloaded"
  else
    warn "Primary dbCAN database command failed; retrying without AWS/no-cgc flags."
    run_in_env "$DBCAN_ENV_NAME" "$DBCAN_ENV_PREFIX/bin/run_dbcan" database --db_dir "$DB_DIR" \
      && ok "dbCAN database downloaded by fallback command" \
      || { fail "dbCAN database download failed. Set --dbcan-db-dir to an existing dbCAN database or rerun later."; exit 1; }
  fi
  check_dbcan_db || warn "dbCAN database was downloaded but expected files were not detected; inspect: $DB_DIR"
fi

if [[ "$SKIP_OPTIONAL" != true ]]; then
  header "Installing optional ML/QC helpers"
  if python_import_in_env "$ENV_NAME" shap; then
    ok "SHAP already installed"
  else
    run_in_env "$ENV_NAME" pip install shap || warn "SHAP installation failed; ML SHAP plots will be skipped."
  fi

  if cmd_in_env "$ENV_NAME" checkm; then
    ok "CheckM already installed"
  else
    conda install -n "$ENV_NAME" -y -c conda-forge -c bioconda checkm-genome \
      && ok "CheckM installed" \
      || warn "CheckM installation failed; run without -q or install manually."
  fi

  if python_import_in_env "$ENV_NAME" resfinder || cmd_in_env "$ENV_NAME" resfinder; then
    ok "ResFinder already available"
  else
    run_in_env "$ENV_NAME" pip install resfinder || warn "ResFinder installation failed; pipeline will continue without it."
  fi
else
  warn "Optional tools skipped."
fi

header "Final prerequisite check"
FAILS=0
check_tool "$ENV_NAME" "Prokka" prokka || ((FAILS++))
check_tool "$ENV_NAME" "AMRFinderPlus" amrfinder || ((FAILS++))
check_tool "$ENV_NAME" "FastANI" fastANI optional || true
if check_dbcan_cli; then
  ok "run_dbCAN [$DBCAN_ENV_NAME]: usable"
else
  fail "run_dbCAN [$DBCAN_ENV_NAME]: broken or NOT FOUND"
  ((FAILS++))
fi
check_tool "$DBCAN_ENV_NAME" "DIAMOND" diamond || ((FAILS++))
check_tool "$DBCAN_ENV_NAME" "HMMER/hmmscan" hmmscan || ((FAILS++))
python_import_in_env "$ENV_NAME" pandas && ok "pandas [$ENV_NAME]" || { fail "pandas [$ENV_NAME]: NOT FOUND"; ((FAILS++)); }
python_import_in_env "$ENV_NAME" sklearn && ok "scikit-learn [$ENV_NAME]" || { fail "scikit-learn [$ENV_NAME]: NOT FOUND"; ((FAILS++)); }
if [[ "$SKIP_DBCAN_DB" != true ]]; then
  check_dbcan_db || ((FAILS++))
fi

echo ""
if [[ "$FAILS" -gt 0 ]]; then
  fail "$FAILS required prerequisites are still missing. See messages above."
  exit 1
fi

ok "All required prerequisites are ready."
echo ""
echo "Use these before manual dbCAN runs:"
echo "  conda activate $DBCAN_ENV_NAME"
echo "  export DBCAN_DB_DIR='$DB_DIR'"
echo ""
echo "Run pipeline:"
echo "  conda activate $ENV_NAME"
echo "  bash run_all.sh -p /path/to/annotations_prokka_GCA -o results -t 32 --skip-install --dbcan-env-name $DBCAN_ENV_NAME --dbcan-db-dir '$DB_DIR'"
