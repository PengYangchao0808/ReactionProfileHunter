#!/usr/bin/env bash
# Parallel S2 recompute for the GFN2-xTB/ALPB benchmark set.
#
# Runs every reaction's S2 stage with a 4-core per-reaction budget and up to
# N reactions concurrently.  The per-reaction core budget is enforced by
# config/defaults.yaml (single source of truth):
#   step2.orca_gfn2_scan.relaxed_scan.nproc     -> GFN2-xTB relaxed scan (4)
#   step2.energy_refinement.parallel_jobs       -> B97-3c SP refinement (4)
#
# Usage:
#   bash scripts/run_s2_parallel_4core.sh [CSV] [OUTPUT_ROOT] [PARALLEL_JOBS]
#
# Examples:
#   bash scripts/run_s2_parallel_4core.sh
#       # -> CSV=data/reaxys_cleaned.csv, root=RPH_Test_Results/s2_gfn2_benchmark_alpb_v2, 4-way parallel
#   bash scripts/run_s2_parallel_4core.sh data/reaxys_cleaned.csv /some/root 8
#       # -> 8-way parallel saturates all 32 cores (8 x 4 cores)
#
# Notes:
#   - S0/S1 assets are reused; only S2 is refreshed (--refresh-stages s2).
#   - Per-reaction console output goes to <root>/rx_id_N/rph_s2_rerun.log.
#   - If the checkpoint hash gate rejects a previously completed S2 (config
#     drift after changing defaults.yaml), delete that reaction's stale
#     S2_PEB and pipeline.state, then re-run; the script prints the failing
#     log tail to make the reason visible.
set -u -o pipefail

csv_path="${1:-data/reaxys_cleaned.csv}"
output_root="${2:-RPH_Test_Results/s2_gfn2_benchmark_alpb_v2}"
parallel_jobs="${3:-4}"
ui_mode="${RPH_UI:-classic}"

if ! [[ "$parallel_jobs" =~ ^[0-9]+$ ]] || (( parallel_jobs < 1 )); then
  echo "[error] PARALLEL_JOBS must be a positive integer: '$parallel_jobs'" >&2
  exit 2
fi
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

declare -A pid_to_rx
declare -A rx_to_rc
running=0

# Reap finished jobs until fewer than parallel_jobs remain running.
# Uses a plain integer counter, never ${#assoc[@]} (empty assoc arrays are
# "unbound" under `set -u` in bash 5.x).
enforce_slots() {
  while (( running >= parallel_jobs )); do
    for pid in "${!pid_to_rx[@]}"; do
      if ! kill -0 "$pid" 2>/dev/null; then
        wait "$pid"; rc=$?
        rx_to_rc[${pid_to_rx[$pid]}]=$rc
        unset "pid_to_rx[$pid]"
        (( running-- ))
      fi
    done
    (( running >= parallel_jobs )) && sleep 5
  done
}

for rx_id in "${reaction_ids[@]}"; do
  work_dir="$output_root/rx_id_$rx_id"
  mkdir -p "$work_dir"
  if [[ ! -f "$work_dir/S0_Mechanism/mechanism.json" ]] \
     || ! find -L "$work_dir/S1_ConfSearch" -type f -name manifest.json -print -quit 2>/dev/null | grep -q .; then
    echo "[skip] rx_id=$rx_id has no reusable S0/S1 baseline"
    rx_to_rc[$rx_id]=99
    continue
  fi

  enforce_slots
  echo "[launch] rx_id=$rx_id ($(( running + 1 )) concurrent, max=$parallel_jobs)"
  bin/rph_run \
    --csv "$csv_path" \
    --rx-id "$rx_id" \
    --output "$work_dir" \
    --refresh-stages s2 \
    --stop-after s2 \
    --ui "$ui_mode" \
    > "$work_dir/rph_s2_rerun.log" 2>&1 &
  pid_to_rx[$!]="$rx_id"
  (( running++ ))
done

for pid in "${!pid_to_rx[@]}"; do
  wait "$pid"; rc=$?
  rx_to_rc[${pid_to_rx[$pid]}]=$rc
  unset "pid_to_rx[$pid]"
done

echo
echo "=== S2 parallel summary ==="
failures=()
for rx_id in "${reaction_ids[@]}"; do
  rc="${rx_to_rc[$rx_id]:-99}"
  printf "  rx_id=%-4s exit=%s\n" "$rx_id" "$rc"
  if (( rc != 0 )); then failures+=("$rx_id"); fi
done

if (( ${#failures[@]} )); then
  echo
  echo "[failed] reactions: ${failures[*]}" >&2
  for rx_id in "${failures[@]}"; do
    log="$output_root/rx_id_$rx_id/rph_s2_rerun.log"
    if [[ -f "$log" ]]; then
      echo "--- tail $log ---" >&2
      tail -n 20 "$log" >&2
    fi
  done
fi

echo
echo "=== collecting GFN2 robustness summary -> $output_root/SUMMARY.csv ==="
python3 - "$output_root" <<'PY'
import csv
import json
import os
import sys

root = sys.argv[1]
rows = []
for dir_name in sorted(os.listdir(root)):
    if not dir_name.startswith("rx_id_"):
        continue
    peb = os.path.join(root, dir_name, "S2_PEB")
    if not os.path.isdir(peb):
        continue
    for variant in sorted(os.listdir(peb)):
        manifest = os.path.join(peb, variant, "manifest.json")
        if not os.path.isfile(manifest):
            continue
        try:
            meta = json.load(open(manifest))
        except (OSError, ValueError):
            continue
        drift = ""
        topology_state = ""
        decision = ""
        scan_profile = os.path.join(peb, variant, "scan_profile.json")
        if os.path.isfile(scan_profile):
            try:
                profile = json.load(open(scan_profile))
                drift = ",".join(str(i) for i in (profile.get("excluded_frames") or []))
                tq = profile.get("trajectory_quality") or {}
                topology_state = str(tq.get("topology_state") or "")
                decision = str((tq.get("topology_rescue_decision") or {}).get("decision") or "")
            except (OSError, ValueError):
                pass
        rows.append({
            "rx_id": dir_name,
            "variant": variant,
            "status": meta.get("status", ""),
            "s2_state": meta.get("s2_state", ""),
            "seed_evidence": meta.get("seed_evidence", ""),
            "ts_seed_index": meta.get("ts_seed_index", ""),
            "int_seed_index": meta.get("int_seed_index", ""),
            "excluded_frames": drift,
            "topology_state": topology_state,
            "rescue_decision": decision,
        })

summary_path = os.path.join(root, "SUMMARY.csv")
with open(summary_path, "w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()) if rows else ["rx_id"])
    writer.writeheader()
    writer.writerows(rows)
print(f"wrote {len(rows)} variant records -> {summary_path}")

seeded = sum(1 for row in rows if row["s2_state"] == "gfn2_seeded")
rescue_seeded = sum(1 for row in rows if row["s2_state"] == "rescue_seeded")
unresolved = sum(1 for row in rows if row["s2_state"] == "unresolved")
drifted = sum(1 for row in rows if row["excluded_frames"])
print(f"gfn2_seeded={seeded} rescue_seeded={rescue_seeded} unresolved={unresolved} drifted={drifted}")
PY

if (( ${#failures[@]} )); then
  exit 1
fi
echo "[ok] all ${#reaction_ids[@]} reactions finished S2"
