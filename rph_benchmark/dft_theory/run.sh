#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

DEFAULTS_CONFIG="$PROJECT_ROOT/config/defaults.yaml"
METHODS_SP_CONFIG="$SCRIPT_DIR/config/methods_sp.yaml"
METHODS_GEO_CONFIG="$SCRIPT_DIR/config/methods_geo.yaml"
CASES_CONFIG="$SCRIPT_DIR/config/benchmark_cases.yaml"
OUTPUT_BASE="$PROJECT_ROOT/Output/benchmark_dft_theory/experiments"
CONFSEARCH_ROOT="$PROJECT_ROOT/Output/benchmark"
REACTION_PROFILE='[4+3]_default'
STAGE="all"
SESSION_DIR=""
RUN_ALL=0
FORCE_CLEAN=0
RESUME=0
DRY_RUN=0
PARALLEL=1

RX_ID_INPUTS=()
RX_IDS=()
SP_METHOD_INPUTS=()
SP_POINT_INPUTS=()
GEO_METHOD_INPUTS=()
GEO_SP_READER="SP-1"

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m'

info()  { echo -e "${CYAN}[INFO]${NC} $*"; }
ok()    { echo -e "${GREEN}[  OK]${NC} $*"; }
warn()  { echo -e "${YELLOW}[WARN]${NC} $*"; }
fail()  { echo -e "${RED}[FAIL]${NC} $*"; }

usage() {
    cat <<'EOF'
Usage:
  bash benchmark/dft_theory/run.sh --all --stage all
  bash benchmark/dft_theory/run.sh --rx-id 1,3,8 --stage migrate
  bash benchmark/dft_theory/run.sh --rx-id 1 --stage sp --methods SP-REF,SP-1

Options:
  --stage migrate|baseline|sp|geo|shermo|report|all
  --rx-id ID[,ID...]            Select reactions by rx_id
  --all                         Run all reactions from benchmark_cases.yaml
  --session-dir DIR             Reuse an existing benchmark session directory
  --output DIR                  Output base directory for new sessions
  --confsearch-root DIR         Root directory of confsearch benchmark outputs
  --reaction-profile KEY        Reaction profile key passed into baseline RPH runs
  --methods LIST                Comma-separated SP method ids
  --geo-methods LIST            Comma-separated GEO method ids
  --sp-reader METHOD            SP reader for GEO stage (default: SP-1)
  --parallel N                  Reserved for future parallel dispatch
  --force-clean                 Remove stage output directories before rerunning
  --resume                      Resume an existing session (requires --session-dir)
  --dry-run                     Reserved, currently unsupported
  -h, --help                    Show this help message
EOF
}

require_file() {
    local path="$1"
    if [ ! -f "$path" ]; then
        fail "Required file not found: $path"
        exit 1
    fi
}

split_csv_to_array() {
    local raw="$1"
    local -n target_ref="$2"
    local item
    IFS=',' read -r -a target_ref <<< "$raw"
    for item in "${target_ref[@]}"; do
        if [ -z "$item" ]; then
            fail "Invalid empty value in comma-separated list: $raw"
            exit 1
        fi
    done
}

append_rx_ids() {
    local raw="$1"
    local parsed=()
    local item
    split_csv_to_array "$raw" parsed
    for item in "${parsed[@]}"; do
        RX_ID_INPUTS+=("$item")
    done
}

dedupe_into() {
    local -n output_ref="$1"
    shift
    declare -A seen=()
    local value
    output_ref=()
    for value in "$@"; do
        if [ -n "$value" ] && [ -z "${seen[$value]+x}" ]; then
            seen[$value]=1
            output_ref+=("$value")
        fi
    done
}

load_all_rx_ids() {
    mapfile -t RX_IDS < <(python3 - "$CASES_CONFIG" <<'PYEOF'
import sys
import yaml

with open(sys.argv[1], encoding="utf-8") as f:
    payload = yaml.safe_load(f) or {}
for item in payload.get("cases", []):
    if isinstance(item, dict):
        value = str(item.get("rx_id", "")).strip()
        if value:
            print(value)
PYEOF
)
}

validate_selection() {
    if [ "$RUN_ALL" -eq 1 ] && [ "${#RX_ID_INPUTS[@]}" -gt 0 ]; then
        fail "Use either --all or --rx-id, not both"
        exit 1
    fi
    if [ "$RUN_ALL" -eq 1 ]; then
        load_all_rx_ids
    else
        if [ "${#RX_ID_INPUTS[@]}" -eq 0 ]; then
            fail "One of --rx-id or --all is required"
            exit 1
        fi
        dedupe_into RX_IDS "${RX_ID_INPUTS[@]}"
    fi
    if [ "${#RX_IDS[@]}" -eq 0 ]; then
        fail "No rx_id values resolved"
        exit 1
    fi
}

validate_stage() {
    case "$STAGE" in
        migrate|baseline|sp|geo|shermo|report|all)
            ;;
        *)
            fail "Unsupported stage: $STAGE"
            exit 1
            ;;
    esac
}

resolve_session_dir() {
    if [ -n "$SESSION_DIR" ]; then
        SESSION_DIR="$(python3 - "$SESSION_DIR" <<'PYEOF'
import os
import sys
print(os.path.abspath(os.path.expanduser(sys.argv[1])))
PYEOF
)"
        mkdir -p "$SESSION_DIR"
        return
    fi
    local stamp
    stamp="$(date '+%Y%m%d_%H%M%S')"
    SESSION_DIR="$OUTPUT_BASE/session_${stamp}"
    mkdir -p "$SESSION_DIR"
}

append_optional_list_arg() {
    local flag="$1"
    shift
    local -a values=("$@")
    if [ "${#values[@]}" -gt 0 ]; then
        printf '%s\0' "$flag" "${values[@]}"
    fi
}

run_migrate() {
    local -a cmd=(python3 "$SCRIPT_DIR/stages/migrate.py" --session-dir "$SESSION_DIR" --cases-config "$CASES_CONFIG" --confsearch-root "$CONFSEARCH_ROOT" --rx-id)
    cmd+=("${RX_IDS[@]}")
    info "Running migrate stage"
    "${cmd[@]}"
}

run_baseline() {
    local -a cmd=(python3 "$SCRIPT_DIR/stages/baseline.py" --session-dir "$SESSION_DIR" --cases-config "$CASES_CONFIG" --defaults-config "$DEFAULTS_CONFIG" --reaction-profile "$REACTION_PROFILE" --rx-id)
    cmd+=("${RX_IDS[@]}")
    if [ "$FORCE_CLEAN" -eq 1 ]; then
        cmd+=(--force-clean)
    fi
    info "Running baseline stage"
    "${cmd[@]}"
}

run_sp() {
    local -a cmd=(python3 "$SCRIPT_DIR/stages/sp_benchmark.py" --session-dir "$SESSION_DIR" --mode all --methods-config "$METHODS_SP_CONFIG" --defaults-config "$DEFAULTS_CONFIG" --rx-id)
    cmd+=("${RX_IDS[@]}")
    if [ "$FORCE_CLEAN" -eq 1 ]; then
        cmd+=(--force-clean)
    fi
    if [ "${#SP_METHOD_INPUTS[@]}" -gt 0 ]; then
        cmd+=(--methods)
        cmd+=("${SP_METHOD_INPUTS[@]}")
    fi
    if [ "${#SP_POINT_INPUTS[@]}" -gt 0 ]; then
        cmd+=(--points)
        cmd+=("${SP_POINT_INPUTS[@]}")
    fi
    info "Running SP benchmark stage"
    "${cmd[@]}"
}

run_geo() {
    local -a cmd=(python3 "$SCRIPT_DIR/stages/geo_benchmark.py" --session-dir "$SESSION_DIR" --mode all --methods-config "$METHODS_GEO_CONFIG" --defaults-config "$DEFAULTS_CONFIG" --sp-reader "$GEO_SP_READER" --rx-id)
    cmd+=("${RX_IDS[@]}")
    if [ "$FORCE_CLEAN" -eq 1 ]; then
        cmd+=(--force-clean)
    fi
    if [ "${#GEO_METHOD_INPUTS[@]}" -gt 0 ]; then
        cmd+=(--methods)
        cmd+=("${GEO_METHOD_INPUTS[@]}")
    fi
    info "Running GEO benchmark stage"
    "${cmd[@]}"
}

run_report() {
    local -a cmd=(python3 "$SCRIPT_DIR/evaluate.py" --session-dir "$SESSION_DIR" --phase final --rx-id)
    cmd+=("${RX_IDS[@]}")
    info "Running final report stage"
    "${cmd[@]}"
}

run_shermo() {
    local -a cmd=(python3 "$SCRIPT_DIR/stages/shermo_correction.py" --session-dir "$SESSION_DIR" --mode all --rx-id)
    cmd+=("${RX_IDS[@]}")
    if [ "$FORCE_CLEAN" -eq 1 ]; then
        cmd+=(--force-clean)
    fi
    if [ "${#SP_METHOD_INPUTS[@]}" -gt 0 ]; then
        cmd+=(--methods)
        cmd+=("${SP_METHOD_INPUTS[@]}")
    fi
    info "Running Shermo correction stage"
    "${cmd[@]}"
}

while [ $# -gt 0 ]; do
    case "$1" in
        --rx-id)
            [ $# -ge 2 ] || { fail "Missing value for --rx-id"; exit 1; }
            append_rx_ids "$2"
            shift 2
            ;;
        --all)
            RUN_ALL=1
            shift
            ;;
        --stage)
            [ $# -ge 2 ] || { fail "Missing value for --stage"; exit 1; }
            STAGE="$2"
            shift 2
            ;;
        --session-dir)
            [ $# -ge 2 ] || { fail "Missing value for --session-dir"; exit 1; }
            SESSION_DIR="$2"
            shift 2
            ;;
        --output)
            [ $# -ge 2 ] || { fail "Missing value for --output"; exit 1; }
            OUTPUT_BASE="$2"
            shift 2
            ;;
        --confsearch-root)
            [ $# -ge 2 ] || { fail "Missing value for --confsearch-root"; exit 1; }
            CONFSEARCH_ROOT="$2"
            shift 2
            ;;
        --reaction-profile)
            [ $# -ge 2 ] || { fail "Missing value for --reaction-profile"; exit 1; }
            REACTION_PROFILE="$2"
            shift 2
            ;;
        --methods)
            [ $# -ge 2 ] || { fail "Missing value for --methods"; exit 1; }
            split_csv_to_array "$2" SP_METHOD_INPUTS
            shift 2
            ;;
        --geo-methods)
            [ $# -ge 2 ] || { fail "Missing value for --geo-methods"; exit 1; }
            split_csv_to_array "$2" GEO_METHOD_INPUTS
            shift 2
            ;;
        --sp-points)
            [ $# -ge 2 ] || { fail "Missing value for --sp-points"; exit 1; }
            split_csv_to_array "$2" SP_POINT_INPUTS
            shift 2
            ;;
        --sp-reader)
            [ $# -ge 2 ] || { fail "Missing value for --sp-reader"; exit 1; }
            GEO_SP_READER="$2"
            shift 2
            ;;
        --parallel)
            [ $# -ge 2 ] || { fail "Missing value for --parallel"; exit 1; }
            PARALLEL="$2"
            shift 2
            ;;
        --force-clean)
            FORCE_CLEAN=1
            shift
            ;;
        --resume)
            RESUME=1
            shift
            ;;
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            fail "Unknown argument: $1"
            exit 1
            ;;
    esac
done

require_file "$DEFAULTS_CONFIG"
require_file "$METHODS_SP_CONFIG"
require_file "$METHODS_GEO_CONFIG"
require_file "$CASES_CONFIG"
validate_stage
validate_selection

if [ "$DRY_RUN" -eq 1 ]; then
    fail "--dry-run is not implemented yet"
    exit 1
fi

if [ "$RESUME" -eq 1 ] && [ -z "$SESSION_DIR" ]; then
    fail "--resume requires --session-dir"
    exit 1
fi

resolve_session_dir

info "DFT benchmark runner"
info "Session dir:      $SESSION_DIR"
info "Stage:            $STAGE"
info "RX IDs:           ${RX_IDS[*]}"
info "Confsearch root:  $CONFSEARCH_ROOT"
info "Reaction profile: $REACTION_PROFILE"
info "Parallel:         $PARALLEL"
echo ""

case "$STAGE" in
    migrate)
        run_migrate
        ;;
    baseline)
        run_baseline
        ;;
    sp)
        run_sp
        ;;
    geo)
        run_geo
        ;;
    shermo)
        run_shermo
        ;;
    report)
        run_report
        ;;
    all)
        run_migrate
        run_baseline
        run_sp
        run_geo
        run_shermo
        run_report
        ;;
esac

ok "Benchmark stage completed"
