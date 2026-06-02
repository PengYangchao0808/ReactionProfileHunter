#!/bin/bash
# ============================================================
# 清理所有包含错误 Shermo 数据的文件（不删除任何 QC 计算数据）
#
# 保留 (计算数据): *.log *.out *.fchk *.chk *.xyz *.gjf *.inp
#                  sp_matrix_metadata.json s3_resume.json
#                  conformer_state.json conformer_energies.json
#                  condition_manifest.json reaction_manifest.json
#                  所有目录结构
#
# 删除 (Shermo 派生数据):
#   - 所有 *_Shermo.sum (S1/S3/condition 三处)
#   - shermo_summary.json (S1 汇总)
#   - thermo_calculation.json (条件热力学, 含 SP→Gibbs 伪装 Bug)
#   - merged_features.csv / merged_features.json (条件特征合并)
#   - features_raw.csv / features_mlr.csv (branch S4, 含错误 dG)
#   - deployable_features.csv / qa_metadata.csv / feature_meta.json (S4)
#   - condition_features_mlr.csv (回填测试残留)
#
# Usage:
#   bash scripts/cleanup_old_shermo.sh              # 仅预览
#   bash scripts/cleanup_old_shermo.sh --execute    # 执行删除
# ============================================================

set -euo pipefail
ROOT="rph_output_backup"
DRY=true

for arg in "$@"; do
    case "$arg" in
        --execute) DRY=false ;;
        --root)    shift; ROOT="$1" ;;
    esac
done

if $DRY; then
    echo ">>> 预览模式: 将列出所有文件（加 --execute 执行删除）<<<"
else
    echo ">>> 执行模式: 正在删除... <<<"
fi
echo ""

# ── 删除 Shermo 输出文件 ──
echo "── Shermo .sum 文件 ──"

echo "  S1 finalDFT (conf_*_Shermo.sum @298K):"
find "$ROOT" -name "conf_*_Shermo.sum" -path "*/finalDFT/*" | while read f; do
    echo "    $f"
    $DRY || rm -v "$f"
done

echo "  S3 shermo/ (ts/reactant_Shermo.sum @298K):"
find "$ROOT" -name "*_Shermo.sum" -path "*/S3_TS/shermo/*" | while read f; do
    echo "    $f"
    $DRY || rm -v "$f"
done
# 删除 S3_TS/shermo 空目录
if ! $DRY; then
    find "$ROOT" -type d -name "shermo" -path "*/S3_TS/shermo" -empty -delete 2>/dev/null || true
fi

echo "  conditions/ (旧 condition 级 Shermo @条件温度):"
find "$ROOT" -name "*_Shermo.sum" -path "*/conditions/*" | while read f; do
    echo "    $f"
    $DRY || rm -v "$f"
done

# ── 删除 Shermo 派生数据文件 ──
echo ""
echo "── Shermo 派生数据 ──"

for pattern in \
    "shermo_summary.json" \
    "thermo_calculation.json" \
    "merged_features.csv" \
    "merged_features.json" \
    "condition_features_mlr.csv" \
; do
    echo "  $pattern:"
    find "$ROOT" -name "$pattern" | while read f; do
        echo "    $f"
        $DRY || rm -v "$f"
    done
done

# ── 删除 S4 特征表（含错误 dG） ──
echo ""
echo "── S4 特征表 (branches/*/S4_Data/) ──"

for pattern in \
    "features_raw.csv" \
    "features_mlr.csv" \
    "feature_meta.json" \
    "deployable_features.csv" \
    "qa_metadata.csv" \
; do
    echo "  $pattern:"
    find "$ROOT" -name "$pattern" -path "*/S4_Data/*" | while read f; do
        echo "    $f"
        $DRY || rm -v "$f"
    done
done

# ── 删除 S4 状态文件 ──
echo ""
echo "── S4 状态文件 ──"
find "$ROOT" -name ".rph_step_status.json" -path "*/S4_Data/*" | while read f; do
    echo "    $f"
    $DRY || rm -v "$f"
done

# ── 汇总 ──
echo ""
if $DRY; then
    count=$(find "$ROOT" \( \
        -name "*_Shermo.sum" \
        -o -name "shermo_summary.json" \
        -o -name "thermo_calculation.json" \
        -o -name "merged_features.csv" \
        -o -name "merged_features.json" \
        -o -name "condition_features_mlr.csv" \
        -o -name "features_raw.csv" \
        -o -name "features_mlr.csv" \
        -o -name "feature_meta.json" \
        -o -name "deployable_features.csv" \
        -o -name "qa_metadata.csv" \
        -o -name ".rph_step_status.json" \
    \) | wc -l)
    echo "=== 预览完成: 将删除 $count 个文件 ==="
    echo "=== 加 --execute 执行删除 ==="
else
    echo "=== 清理完成 ==="
    echo ""
    echo "保留的计算数据 (未删除):"
    echo "  *.log *.out *.fchk *.chk *.xyz *.gjf *.inp"
    echo "  sp_matrix_metadata.json s3_resume.json"
    echo "  conformer_state.json conformer_energies.json conformer_thermo.csv"
    echo "  condition_manifest.json reaction_manifest.json"
    echo "  所有目录结构"
fi
