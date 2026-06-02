#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
RUNNER="$SCRIPT_DIR/run.sh"

SESSION_DIR="${1:-/tmp/rph_dft_bench/session_smoke}"
CONFSEARCH_ROOT="${CONFSEARCH_ROOT:-$PROJECT_ROOT/Output/benchmark}"
RX_IDS="${RX_IDS:-1}"
SP_METHODS="${SP_METHODS:-SP-REF,SP-1}"
GEO_METHODS="${GEO_METHODS:-GEO-1,GEO-3}"

echo "[smoke] session_dir=$SESSION_DIR"
echo "[smoke] confsearch_root=$CONFSEARCH_ROOT"
echo "[smoke] rx_ids=$RX_IDS"
echo "[smoke] sp_methods=$SP_METHODS"
echo "[smoke] geo_methods=$GEO_METHODS"

bash "$RUNNER" --stage migrate --rx-id "$RX_IDS" --session-dir "$SESSION_DIR" --confsearch-root "$CONFSEARCH_ROOT"
bash "$RUNNER" --stage baseline --rx-id "$RX_IDS" --session-dir "$SESSION_DIR" --force-clean
bash "$RUNNER" --stage sp --rx-id "$RX_IDS" --session-dir "$SESSION_DIR" --methods "$SP_METHODS" --force-clean
bash "$RUNNER" --stage geo --rx-id "$RX_IDS" --session-dir "$SESSION_DIR" --geo-methods "$GEO_METHODS" --force-clean
bash "$RUNNER" --stage report --rx-id "$RX_IDS" --session-dir "$SESSION_DIR"

echo "[smoke] completed"
