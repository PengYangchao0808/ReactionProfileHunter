#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
CONFIG="$PROJECT_ROOT/config/defaults.yaml"
DATASET="$PROJECT_ROOT/data/reaxys_cleaned.csv"
OUTPUT_BASE="$PROJECT_ROOT/Output/benchmark"
REACTION_TYPE='[4+3]_default'
FORCE_CLEAN=1
RUN_ALL=0

DEFAULT_PROTOCOLS=(ext full lite zero)
PROTOCOLS=("${DEFAULT_PROTOCOLS[@]}")
RX_ID_INPUTS=()
RX_IDS=()

declare -A RX_PROTO_STATUS=()
declare -A RX_PASS_COUNT=()
declare -A RX_FAIL_COUNT=()

TOTAL_ERRORS=0

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
  bash benchmark/confsearch/run.sh --rx-id 1
  bash benchmark/confsearch/run.sh --rx-id 1,3,8
  bash benchmark/confsearch/run.sh --all
  bash benchmark/confsearch/run.sh --rx-id 1 --protocols ext,lite
  bash benchmark/confsearch/run.sh --rx-id 1 --reaction-type [4+3]_default

Arguments:
  --rx-id ID           Select reactions by rx_id from CSV (repeatable or comma-separated)
  --all                Run all reactions from CSV
  --protocols LIST     Comma-separated protocols (default: ext,full,lite,zero)
  --reaction-type TYPE Override reaction type (default: [4+3]_default)
  --output DIR         Override output base directory (default: PROJECT_ROOT/Output/benchmark)
  --force-clean        Delete existing output dirs before running (default)
  --no-clean           Keep existing output dirs (resume mode)
  -h, --help           Show this help message
EOF
}

timestamp() {
    date '+%Y-%m-%d %H:%M:%S'
}

cleanup() {
    local exit_code=$?
    if [ -d "$OUTPUT_BASE" ]; then
        find "$OUTPUT_BASE" -name ".bench_config.yaml" -delete 2>/dev/null || true
    fi
    exit "$exit_code"
}

trap cleanup EXIT

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

dedupe_rx_ids() {
    local seen_key
    declare -A seen=()
    local rx_id
    RX_IDS=()
    for rx_id in "$@"; do
        seen_key="$rx_id"
        if [ -n "$rx_id" ] && [ -z "${seen[$seen_key]+x}" ]; then
            seen[$seen_key]=1
            RX_IDS+=("$rx_id")
        fi
    done
}

load_all_rx_ids() {
    mapfile -t RX_IDS < <(python3 -c "
import csv
with open('$DATASET', newline='') as f:
    reader = csv.reader(f)
    next(reader, None)
    for row in reader:
        if row and row[0].strip():
            print(row[0].strip())
" )
}

get_reaction_info() {
    local rx_id="$1"
    python3 -c "
import csv, json, sys
with open('$DATASET', newline='') as f:
    for row in csv.DictReader(f):
        if row.get('rx_id','').strip() == '$rx_id':
            info = {
                'product_smiles': row.get('product_smiles_main','').strip(),
                'precursor_smiles': row.get('precursor_smiles','').strip(),
            }
            print(json.dumps(info))
            sys.exit(0)
sys.exit(1)
"
}

validate_protocol_name() {
    local proto="$1"
    case "$proto" in
        ext|full|lite|zero)
            ;;
        *)
            fail "Unsupported protocol: $proto"
            exit 1
            ;;
    esac
}

generate_temp_config() {
    local proto="$1"
    local outdir="$2"
    local tmp_config="$outdir/.bench_config.yaml"

    python3 - "$tmp_config" "$CONFIG" "$DATASET" "$outdir" "$proto" <<'PYEOF'
import sys, yaml
tmp_config_path = sys.argv[1]
defaults_path = sys.argv[2]
dataset_path = sys.argv[3]
output_root = sys.argv[4]
proto = sys.argv[5]

with open(defaults_path) as f:
    cfg = yaml.safe_load(f)

cfg['run']['source'] = 'dataset'
cfg['run']['dataset']['path'] = dataset_path
cfg['run']['dataset']['delimiter'] = ','
cfg['run']['output_root'] = output_root
cfg['run']['max_tasks'] = 0
cfg['run']['filter_ids'] = []
cfg['run']['resume'] = False
cfg['step1']['protocol'] = proto

with open(tmp_config_path, 'w') as f:
    yaml.dump(cfg, f, default_flow_style=False, allow_unicode=True)
PYEOF
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
            usage
            exit 1
        fi
        dedupe_rx_ids "${RX_ID_INPUTS[@]}"
    fi

    if [ "${#RX_IDS[@]}" -eq 0 ]; then
        fail "No rx_id values resolved from dataset"
        exit 1
    fi
}

validate_rx_ids_exist() {
    local rx_id info product_smiles
    for rx_id in "${RX_IDS[@]}"; do
        if ! info=$(get_reaction_info "$rx_id" 2>/dev/null); then
            fail "rx_id not found in dataset: $rx_id"
            exit 1
        fi
        product_smiles=$(INFO_JSON="$info" python3 -c "import json, os; print(json.loads(os.environ['INFO_JSON']).get('product_smiles', ''))")
        if [ -z "$product_smiles" ]; then
            fail "Empty product_smiles_main for rx_id: $rx_id"
            exit 1
        fi
    done
}

validate_outputs() {
    local outdir="$1"
    local rx_id="$2"
    local proto="$3"
    local errors=0
    local provenance=""
    local proto_in_json=""
    local product_xyz=""
    local s2_dir=""
    local s3_dir=""
    local s4_dir=""

    provenance=$(find "$outdir" -name "provenance.json" -path "*/S1_ConfGeneration/*" -print -quit 2>/dev/null || true)
    if [ -z "$provenance" ]; then
        fail "rx${rx_id}/${proto}: provenance.json not found"
        errors=$((errors + 1))
    else
        ok "rx${rx_id}/${proto}: provenance.json found"
        proto_in_json=$(python3 -c "
import json
with open('$provenance') as f:
    data = json.load(f)
print(data.get('protocol', 'MISSING'))
" 2>/dev/null || echo "PARSE_ERROR")
        if [ "$proto_in_json" = "$proto" ]; then
            ok "rx${rx_id}/${proto}: provenance.protocol = ${proto_in_json}"
        else
            fail "rx${rx_id}/${proto}: provenance.protocol = ${proto_in_json}, expected ${proto}"
            errors=$((errors + 1))
        fi
    fi

    product_xyz=$(find "$outdir" \( -name "product_min.xyz" -o -name "product_global_min.xyz" \) -print -quit 2>/dev/null || true)
    if [ -n "$product_xyz" ]; then
        ok "rx${rx_id}/${proto}: product structure found"
    else
        fail "rx${rx_id}/${proto}: product_min.xyz/product_global_min.xyz not found"
        errors=$((errors + 1))
    fi

    local precursor_xyz=""
    precursor_xyz=$(find "$outdir" -name "precursor_global_min.xyz" -print -quit 2>/dev/null || true)
    if [ -n "$precursor_xyz" ]; then
        ok "rx${rx_id}/${proto}: precursor structure found"
    else
        warn "rx${rx_id}/${proto}: precursor structure not found"
    fi

    s2_dir=$(find "$outdir" -maxdepth 2 -type d -name "S2_Retro" -print -quit 2>/dev/null || true)
    s3_dir=$(find "$outdir" -maxdepth 2 -type d -name "S3_*" -print -quit 2>/dev/null || true)
    s4_dir=$(find "$outdir" -maxdepth 2 -type d -name "S4_Data" -print -quit 2>/dev/null || true)
    if [ -z "$s2_dir" ] && [ -z "$s3_dir" ] && [ -z "$s4_dir" ]; then
        ok "rx${rx_id}/${proto}: S2/S3/S4 absent as expected"
    else
        fail "rx${rx_id}/${proto}: unexpected directories present: S2=${s2_dir:-N/A} S3=${s3_dir:-N/A} S4=${s4_dir:-N/A}"
        errors=$((errors + 1))
    fi

    return "$errors"
}

run_one() {
    local rx_id="$1"
    local proto="$2"
    local outdir="$OUTPUT_BASE/rx${rx_id}/${proto}"
    local start_epoch end_epoch elapsed status

    info "=========================================="
    info "[$(timestamp)] Start rx_id=${rx_id} protocol=${proto} (dataset mode)"
    info "Output: ${outdir}"
    info "=========================================="

    if [ -d "$outdir" ]; then
        if [ "$FORCE_CLEAN" -eq 1 ]; then
            warn "Existing output found, removing: $outdir"
            rm -rf "$outdir"
            ok "Cleaned previous output"
        else
            warn "Existing output found, keeping for resume: $outdir"
        fi
    fi

    mkdir -p "$outdir"
    generate_temp_config "$proto" "$outdir"
    ok "Generated temp config for protocol=${proto}"

    start_epoch=$(date +%s)

    if (
        cd "$PROJECT_ROOT"
        python -m rph_core \
            --config "$outdir/.bench_config.yaml" \
            --rx-id "$rx_id" \
            --reaction-type "$REACTION_TYPE" \
            --skip-steps s2,s3,s4 \
            --log-level INFO
    ); then
        status=0
    else
        status=$?
        fail "rx${rx_id}/${proto}: run failed with exit code ${status}"
    fi

    end_epoch=$(date +%s)
    elapsed=$((end_epoch - start_epoch))

    if [ "$status" -eq 0 ]; then
        if validate_outputs "$outdir" "$rx_id" "$proto"; then
            ok "[$(timestamp)] End rx_id=${rx_id} protocol=${proto} (${elapsed}s)"
            return 0
        fi
    fi

    fail "[$(timestamp)] End rx_id=${rx_id} protocol=${proto} (${elapsed}s)"
    return 1
}

run_evaluation() {
    require_file "$SCRIPT_DIR/evaluate.py"
    info "Running post-run evaluation ..."
    if python3 "$SCRIPT_DIR/evaluate.py" \
        --rx-ids "${RX_IDS[@]}" \
        --protocols "${PROTOCOLS[@]}" \
        --output-base "$OUTPUT_BASE"; then
        ok "Evaluation completed"
    else
        fail "Evaluation failed"
        TOTAL_ERRORS=$((TOTAL_ERRORS + 1))
    fi
}

print_summary() {
    local rx_id passed failed total
    echo ""
    info "=========================================="
    info " Benchmark summary"
    info "=========================================="
    printf '%-12s %-8s %-8s %-8s\n' "rx_id" "passed" "failed" "total"
    for rx_id in "${RX_IDS[@]}"; do
        passed=${RX_PASS_COUNT[$rx_id]:-0}
        failed=${RX_FAIL_COUNT[$rx_id]:-0}
        total=$((passed + failed))
        printf '%-12s %-8s %-8s %-8s\n' "$rx_id" "$passed" "$failed" "$total"
    done
    echo ""
    if [ "$TOTAL_ERRORS" -eq 0 ]; then
        ok "Total errors: 0"
    else
        fail "Total errors: ${TOTAL_ERRORS}"
    fi
    for rx_id in "${RX_IDS[@]}"; do
        info "Evaluation report path (rx${rx_id}): ${OUTPUT_BASE}/rx${rx_id}/evaluation/summary.md"
    done
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
        --protocols)
            [ $# -ge 2 ] || { fail "Missing value for --protocols"; exit 1; }
            split_csv_to_array "$2" PROTOCOLS
            shift 2
            ;;
        --reaction-type)
            [ $# -ge 2 ] || { fail "Missing value for --reaction-type"; exit 1; }
            REACTION_TYPE="$2"
            shift 2
            ;;
        --output)
            [ $# -ge 2 ] || { fail "Missing value for --output"; exit 1; }
            OUTPUT_BASE="$2"
            shift 2
            ;;
        --force-clean)
            FORCE_CLEAN=1
            shift
            ;;
        --no-clean)
            FORCE_CLEAN=0
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            fail "Unknown argument: $1"
            usage
            exit 1
            ;;
    esac
done

require_file "$CONFIG"
require_file "$DATASET"
validate_selection

for proto in "${PROTOCOLS[@]}"; do
    validate_protocol_name "$proto"
done

validate_rx_ids_exist
mkdir -p "$OUTPUT_BASE"

info "ReactionProfileHunter benchmark runner (dataset mode)"
info "Dataset:        ${DATASET}"
info "Reaction type:  ${REACTION_TYPE}"
info "Protocols:      ${PROTOCOLS[*]}"
info "RX IDs:         ${RX_IDS[*]}"
info "Output base:    ${OUTPUT_BASE}"
info "Force clean:    ${FORCE_CLEAN}"
echo ""

for rx_id in "${RX_IDS[@]}"; do
    RX_PASS_COUNT[$rx_id]=0
    RX_FAIL_COUNT[$rx_id]=0
    info "Resolved rx_id=${rx_id} (dataset mode: product + precursor)"
    for proto in "${PROTOCOLS[@]}"; do
        if run_one "$rx_id" "$proto"; then
            RX_PROTO_STATUS["${rx_id}:${proto}"]=PASS
            RX_PASS_COUNT[$rx_id]=$((RX_PASS_COUNT[$rx_id] + 1))
        else
            RX_PROTO_STATUS["${rx_id}:${proto}"]=FAIL
            RX_FAIL_COUNT[$rx_id]=$((RX_FAIL_COUNT[$rx_id] + 1))
            TOTAL_ERRORS=$((TOTAL_ERRORS + 1))
        fi
        echo ""
    done
done

run_evaluation
print_summary

if [ "$TOTAL_ERRORS" -ne 0 ]; then
    exit 1
fi
