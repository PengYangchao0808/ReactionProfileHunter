#!/usr/bin/env bash
# Serial S3 batch with fail-fast for the GFN2-xTB/ALPB benchmark set.
#
# Requirements (enforced by config/defaults.yaml):
#   1. SERIAL   : one reaction at a time (no concurrency)
#   2. RESOURCES: 16 cores / 32 GB per job (refinement.s3.resources)
#   3. FAIL-FAST: every S3 rescue ladder is disabled (refinement.common.*)
#      and the batch STOPS at the first failed structure so the root cause
#      can be investigated without burning compute on remaining reactions.
#
# S2 outputs are reused as-is (no --refresh-stages s2).
#
# Usage:
#   bash scripts/run_s3_serial.sh [CSV] [OUTPUT_ROOT]
#
# Examples:
#   bash scripts/run_s3_serial.sh
#   bash scripts/run_s3_serial.sh data/reaxys_cleaned.csv RPH_Test_Results/s2_gfn2_benchmark_alpb_v2
set -u -o pipefail

csv_path="${1:-data/reaxys_cleaned.csv}"
output_root="${2:-RPH_Test_Results/s2_gfn2_benchmark_alpb_v2}"
ui_mode="${RPH_UI:-classic}"

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
if (( ${#reaction_ids[@]} == 0 )); then
  echo "[error] no rx_id entries found in $csv_path" >&2
  exit 2
fi

declare -A results

check_s3_manifest() {
  # Exit 0 = all structures complete; exit 1 = failed structures; exit 2 = no manifest.
  python3 - "$1" <<'PY'
import json
import os
import sys

work = sys.argv[1]
manifest = os.path.join(work, "S3_LowLevel", "manifest.json")
if not os.path.isfile(manifest):
    print("NO_S3_MANIFEST", file=sys.stderr)
    sys.exit(2)

with open(manifest, encoding="utf-8") as handle:
    meta = json.load(handle)

structures = meta.get("structures") or []
failed = [s for s in structures if s.get("status") != "complete"]
print(f"structures={len(structures)} complete={sum(1 for s in structures if s.get('status') == 'complete')} failed={len(failed)}")
for s in failed:
    print(f"  FAILED {s.get('id')} status={s.get('status')} error={s.get('error') or s.get('opt_error') or ''}")
sys.exit(1 if failed else 0)
PY
}

for rx_id in "${reaction_ids[@]}"; do
  work_dir="$output_root/rx_id_$rx_id"
  log="$work_dir/rph_s3_rerun.log"

  if ! find -L "$work_dir/S2_PEB" -maxdepth 2 -name manifest.json -print -quit 2>/dev/null | grep -q .; then
    echo "[error] rx_id=$rx_id has no completed S2 manifest; S3 requires it" >&2
    echo "        run S2 first (scripts/run_s2_parallel_4core.sh) or fix the S2 outputs" >&2
    exit 3
  fi

  mkdir -p "$work_dir"
  started=$(date '+%Y-%m-%d %H:%M:%S')
  echo "[S3 serial] rx_id=$rx_id start=$started"
  bin/rph_run \
    --csv "$csv_path" \
    --rx-id "$rx_id" \
    --output "$work_dir" \
    --resume-policy use-existing-upstream \
    --start-from s3 \
    --stop-after s3 \
    --ui "$ui_mode" \
    > "$log" 2>&1
  rc=$?
  finished=$(date '+%Y-%m-%d %H:%M:%S')
  echo "[S3 serial] rx_id=$rx_id finish=$finished exit=$rc"

  if (( rc != 0 )); then
    echo
    echo "[FAIL-FAST] rx_id=$rx_id pipeline exited with rc=$rc" >&2
    echo "--- tail $log ---" >&2
    tail -n 30 "$log" >&2
    results[$rx_id]="pipeline_rc=$rc"
    break
  fi

  check_s3_manifest "$work_dir"
  mrc=$?
  if (( mrc != 0 )); then
    echo
    echo "[FAIL-FAST] rx_id=$rx_id has failed structures (manifest check rc=$mrc)" >&2
    echo "--- tail $log ---" >&2
    tail -n 20 "$log" >&2
    results[$rx_id]="manifest_rc=$mrc"
    break
  fi
  results[$rx_id]="ok"
done

echo
echo "=== S3 serial summary ==="
for rx_id in "${reaction_ids[@]}"; do
  printf "  rx_id=%-4s %s\n" "$rx_id" "${results[$rx_id]:-not_run}"
done
if [[ " ${results[*]} " != *" ok "* ]] || [[ " ${results[*]} " == *"not_run"* ]]; then
  echo
  echo "[stopped] batch halted before all reactions ran (fail-fast)." >&2
  exit 1
fi
echo "[ok] all ${#reaction_ids[@]} reactions finished S3"
