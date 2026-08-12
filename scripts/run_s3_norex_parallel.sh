#!/usr/bin/env bash
# S3 primary-only batch (rescue disabled) for S2-seed robustness testing.
#
# Runs rx_id 2..20 (rx_id=1 already complete 7/7) through S3 with the rescue
# matrix OFF (refinement.common.rescue.matrix.enabled=false in defaults.yaml):
#   - failures are classified F1-F6 and recorded (S3.4.0_diagnosis)
#   - no rescue methods run -> shortest wall time per reaction
#   - NO fail-fast: every reaction runs to completion so the S2 point-selection
#     robustness can be scored across the whole cluster.
#
# Machine target: 32 cores / 54 GB -> PARALLEL=2 reactions (16 cores / 32 GB each).
# Usage:
#   bash scripts/run_s3_norex_parallel.sh [CSV] [OUTPUT_ROOT] [PARALLEL]
set -u -o pipefail

csv_path="${1:-data/reaxys_cleaned.csv}"
output_root="${2:-RPH_Test_Results/s2_gfn2_benchmark_alpb_v2}"
parallel="${3:-2}"

if [[ ! -f "$csv_path" ]]; then
  echo "[error] CSV not found: $csv_path" >&2
  exit 2
fi

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

mapfile -t reaction_ids < <(
  python3 - "$csv_path" <<'PY'
import csv
import sys

with open(sys.argv[1], encoding="utf-8-sig", newline="") as handle:
    reader = csv.DictReader(handle)
    if "rx_id" not in (reader.fieldnames or []):
        raise SystemExit("CSV has no rx_id column")
    for row in reader:
        value = (row.get("rx_id") or "").strip()
        if value:
            print(value)
PY
)

summary="${TMPDIR:-/tmp}/s3_norex_summary.txt"
: > "$summary"

run_one() {
  local rx="$1"
  local w="$output_root/rx_id_$rx"
  if ! find -L "$w/S2_PEB" -maxdepth 2 -name manifest.json -print -quit 2>/dev/null | grep -q .; then
    echo "rx_id=$rx NO_S2 (skip)" >> "$summary"
    return
  fi
  if python3 - "$w/S3_LowLevel/manifest.json" <<'PY'
import json
import os
import sys

path = sys.argv[1]
if not os.path.isfile(path):
    sys.exit(1)
meta = json.load(open(path, encoding="utf-8"))
failed = [s for s in meta.get("structures", []) if s.get("status") != "complete"]
sys.exit(0 if not failed else 1)
PY
  then
    echo "rx_id=$rx S3_ALREADY_COMPLETE (skip)" >> "$summary"
    return
  fi
  mkdir -p "$w"
  local log="$w/rph_s3_norex.log"
  bin/rph_run \
    --csv "$csv_path" \
    --rx-id "$rx" \
    --output "$w" \
    --resume-policy use-existing-upstream \
    --start-from s3 \
    --stop-after s3 \
    --ui classic \
    > "$log" 2>&1
  echo "rx_id=$rx rc=$?" >> "$summary"
}

echo "[batch] S3 primary-only (rescue disabled), parallel=$parallel, reactions: ${reaction_ids[*]}"
started=$(date '+%Y-%m-%d %H:%M:%S')
echo "[batch] start=$started"

for rx in "${reaction_ids[@]}"; do
  # throttle: wait until fewer than `parallel` jobs are alive
  while (( $(jobs -rp | wc -l) >= parallel )); do
    sleep 20
  done
  run_one "$rx" &
done
wait

finished=$(date '+%Y-%m-%d %H:%M:%S')
echo "[batch] finish=$finished"
echo
echo "=== per-reaction rc ==="
sort -t= -k2 -n "$summary"
