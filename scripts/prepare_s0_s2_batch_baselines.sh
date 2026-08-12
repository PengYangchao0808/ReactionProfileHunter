#!/usr/bin/env bash
# Bring the already validated S0/S1 assets from the earlier robustness runs
# into the common batch tree.  Relative S1 paths and selected geometries are
# retained through directory links, avoiding an unnecessary 10+ GiB duplicate
# of CREST/SP scratch files.  S2 and run-state artifacts are deliberately
# excluded so every reaction is recomputed by the current S2 implementation.
# Invoke from WSL: bash scripts/prepare_s0_s2_batch_baselines.sh
set -euo pipefail

output_root="${1:-RPH_Test_Results/s0_s2_batch}"
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

declare -A source_runs=(
  [1]="RPH_Test_Results/rx_id_1_test3"
  [2]="RPH_Test_Results/rx_id_2_test"
  [3]="RPH_Test_Results/robustness_rx_3"
  [8]="RPH_Test_Results/robustness_rx_8"
  [15]="RPH_Test_Results/robustness_rx_15"
)

has_s1_manifest() {
  find -L "$1/S1_ConfSearch" -type f -name manifest.json -print -quit 2>/dev/null | grep -q .
}

for rx_id in 1 2 3 8 15; do
  source_dir="${source_runs[$rx_id]}"
  target_dir="$output_root/rx_id_$rx_id"

  if [[ ! -f "$source_dir/S0_Mechanism/mechanism.json" ]] || ! has_s1_manifest "$source_dir"; then
    echo "[error] source rx_id=$rx_id lacks a reusable S0 or S1 manifest: $source_dir" >&2
    exit 2
  fi

  if [[ -e "$target_dir/S0_Mechanism" || -e "$target_dir/S1_ConfSearch" ]]; then
    if [[ -f "$target_dir/S0_Mechanism/mechanism.json" ]] && has_s1_manifest "$target_dir"; then
      echo "[reuse] rx_id=$rx_id already has S0/S1 in $target_dir"
      continue
    fi
    echo "[error] target rx_id=$rx_id contains incomplete S0/S1 assets: $target_dir" >&2
    echo "        Resolve it manually rather than overwriting a partial run." >&2
    exit 3
  fi

  mkdir -p "$target_dir"
  ln -s "$(realpath "$source_dir/S0_Mechanism")" "$target_dir/S0_Mechanism"
  ln -s "$(realpath "$source_dir/S1_ConfSearch")" "$target_dir/S1_ConfSearch"
  echo "[linked] rx_id=$rx_id source=$source_dir"
done

echo "S0/S1 baselines are ready under $output_root.  No S2 or pipeline state was copied."
