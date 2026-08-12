#!/usr/bin/env bash
# Recompute S2 for every reaction in the shared batch tree while reusing the
# trusted S0/S1 artifacts.  Invoke from WSL: bash scripts/refresh_s2_batch.sh
set -u -o pipefail

csv_path="${1:-data/reaxys_cleaned.csv}"
output_root="${2:-RPH_Test_Results/s0_s2_batch}"
ui_mode="${RPH_UI:-auto}"

if [[ ! -f "$csv_path" ]]; then
  echo "CSV not found: $csv_path" >&2
  exit 2
fi

bash scripts/prepare_s0_s2_batch_baselines.sh "$output_root"

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

failures=()
for rx_id in "${reaction_ids[@]}"; do
  work_dir="$output_root/rx_id_$rx_id"
  if [[ ! -f "$work_dir/S0_Mechanism/mechanism.json" ]] || ! find -L "$work_dir/S1_ConfSearch" -type f -name manifest.json -print -quit 2>/dev/null | grep -q .; then
    echo "[failed] rx_id=$rx_id has no reusable S0/S1 baseline: $work_dir" >&2
    failures+=("$rx_id")
    continue
  fi

  echo "[refresh S2] rx_id=$rx_id output=$work_dir"
  bin/rph_run \
    --csv "$csv_path" \
    --rx-id "$rx_id" \
    --output "$work_dir" \
    --refresh-stages s2 \
    --stop-after s2 \
    --ui "$ui_mode"
  run_status=$?
  if (( run_status != 0 )); then
    echo "[failed] rx_id=$rx_id; continuing with the next reaction" >&2
    failures+=("$rx_id")
  fi
done

if (( ${#failures[@]} )); then
  echo "S2 refresh finished with failures: ${failures[*]}" >&2
  exit 1
fi

echo "All S2 refreshes complete. Outputs: $output_root"
