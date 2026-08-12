#!/usr/bin/env bash
# Run S0 -> S2 only for reaction IDs that do not yet have an accepted S2 run.
# Invoke from WSL: bash scripts/run_s0_s2_uncomputed.sh
set -u -o pipefail

csv_path="${1:-data/reaxys_cleaned.csv}"
output_root="${2:-RPH_Test_Results/s0_s2_batch}"
# Override when the baseline changes, for example:
# RPH_COMPLETED_IDS="1 2 3 8 15" bash scripts/run_s0_s2_uncomputed.sh
completed_ids="${RPH_COMPLETED_IDS:-1 2 3 8 15}"
# Keep the normal live dashboard for an interactive WSL terminal.  Users can
# override this for CI/log capture, e.g. RPH_UI=classic bash scripts/....
ui_mode="${RPH_UI:-auto}"

if [[ ! -f "$csv_path" ]]; then
  echo "CSV not found: $csv_path" >&2
  exit 2
fi

mkdir -p "$output_root"
mapfile -t reaction_ids < <(
  python - "$csv_path" <<'PY'
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

is_completed_baseline() {
  local rx_id="$1"
  for completed_id in $completed_ids; do
    [[ "$rx_id" == "$completed_id" ]] && return 0
  done
  return 1
}

failures=()
for rx_id in "${reaction_ids[@]}"; do
  if is_completed_baseline "$rx_id"; then
    echo "[skip baseline] rx_id=$rx_id"
    continue
  fi

  work_dir="$output_root/rx_id_$rx_id"
  if [[ -f "$work_dir/S2_PEB/manifest.json" ]]; then
    echo "[skip complete] rx_id=$rx_id s2_manifest=$work_dir/S2_PEB/manifest.json"
    continue
  fi

  echo "[run] rx_id=$rx_id output=$work_dir"
  bin/rph_run \
    --csv "$csv_path" \
    --rx-id "$rx_id" \
    --output "$work_dir" \
    --stop-after s2 \
    --ui "$ui_mode"
  run_status=$?
  if (( run_status != 0 )); then
    echo "[failed] rx_id=$rx_id; continuing with the next reaction" >&2
    failures+=("$rx_id")
  fi
done

if (( ${#failures[@]} )); then
  echo "S0-S2 batch finished with failures: ${failures[*]}" >&2
  exit 1
fi

echo "S0-S2 batch complete. Outputs: $output_root"
