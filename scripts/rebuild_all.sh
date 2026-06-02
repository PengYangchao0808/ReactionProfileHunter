#!/bin/bash
# ============================================================
# RPH 数据重建: Phase 1 → Phase 2 → Phase 3
#
# Usage:
#   bash scripts/rebuild_all.sh              # 预览
#   bash scripts/rebuild_all.sh --execute    # 执行全部
#   bash scripts/rebuild_all.sh --phase 1   # 仅 Phase 1
#   bash scripts/rebuild_all.sh --phase 2   # 仅 Phase 2
#   bash scripts/rebuild_all.sh --phase 3   # 仅 Phase 3
# ============================================================
set -euo pipefail

DIE() { echo "ERROR: $*" >&2; exit 1; }
OK()  { echo "[OK] $*"; }

DRY=true
PHASE=""

for arg in "$@"; do
    case "$arg" in
        --execute) DRY=false ;;
        --phase)   shift; PHASE="$1" ;;
        1|2|3)     PHASE="$arg" ;;
    esac
done

RUN_PHASE1() { [[ -z "$PHASE" || "$PHASE" == "1" ]]; }
RUN_PHASE2() { [[ -z "$PHASE" || "$PHASE" == "2" ]]; }
RUN_PHASE3() { [[ -z "$PHASE" || "$PHASE" == "3" ]]; }

EXTRA=""
$DRY || EXTRA="--execute"

echo "========================================"
echo " RPH 数据重建: Phase 1→2→3"
echo " 模式: $($DRY && echo '预览' || echo '执行')"
echo "========================================"
echo ""

# ── Phase 1: S4 特征提取 ──
if RUN_PHASE1; then
    echo "── Phase 1: S4 特征提取 (branch 级, dG=NaN) ──"
    bash scripts/rerun_s4_features.sh $EXTRA
    OK "Phase 1 done"
    echo ""
fi

# ── Phase 2: 条件热力学回填 ──
if RUN_PHASE2; then
    echo "── Phase 2: 条件热力学回填 (condition 级, dG @条件温度) ──"
    python scripts/backfill_condition_thermo.py \
        --rph-root rph_output_backup \
        --config config/defaults.yaml
    OK "Phase 2 done"
    echo ""
fi

# ── Phase 3: ML 数据集构建 ──
if RUN_PHASE3; then
    echo "── Phase 3: ML 数据集构建 ──"
    python -m rph_features.cli build-dataset \
        --rph-root rph_output_backup \
        --output data/processed_rph_v2.2 \
        --barrier-threshold 10.0
    OK "Phase 3 done"
    echo ""
fi

echo "========================================"
echo " 全部完成"
echo " 产出:"
echo "   branches/*/S4_Data/features_raw.csv   (Phase 1)"
echo "   conditions/*/branches/*/condition_features_mlr.csv (Phase 2)"
echo "   data/processed_rph_v2.2/              (Phase 3)"
echo "========================================"
