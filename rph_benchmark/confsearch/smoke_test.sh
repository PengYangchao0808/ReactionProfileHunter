#!/usr/bin/env bash
# ============================================================================
# smoke_test.sh — S1 协议 (ext / full / lite / zero) 快速冒烟测试
# ============================================================================
# 定位: 单反应快速验证脚本，与 confsearch/run.sh 正式 benchmark 互补。
# 用法:
#   bash benchmark/confsearch/smoke_test.sh          # 依次跑 ext → full → lite → zero
#   bash benchmark/confsearch/smoke_test.sh lite     # 只跑指定协议
#   bash benchmark/confsearch/smoke_test.sh ext zero # 跑 ext 和 zero
#
# 前置条件:
#   - Gaussian / ORCA / xTB / CREST 已在 defaults.yaml 中正确配置
#   - Python 环境可用
# ============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# 配置
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
CONFIG="$PROJECT_ROOT/config/defaults.yaml"
CLI="python -m rph_core"

# CSV 第一条数据 (rx_id=1)
PRODUCT_SMILES='CC(C)(C)OC(=O)N1CCC[C@@]23C=C[C@@H](CC(=O)[C@H]12)O3'
REACTION_TYPE='[4+3]_default'
OUTPUT_BASE="$PROJECT_ROOT/Output/s1_protocol_test"

ALL_PROTOCOLS=(ext full lite zero)
FORCE_CLEAN=${FORCE_CLEAN:-1}

# 如果命令行指定了协议，只跑指定的
if [ $# -gt 0 ]; then
    ALL_PROTOCOLS=("$@")
fi

# ---------------------------------------------------------------------------
# 颜色辅助
# ---------------------------------------------------------------------------
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

info()  { echo -e "${CYAN}[INFO]${NC} $*"; }
ok()    { echo -e "${GREEN}[  OK]${NC} $*"; }
warn()  { echo -e "${YELLOW}[WARN]${NC} $*"; }
fail()  { echo -e "${RED}[FAIL]${NC} $*"; }

# ---------------------------------------------------------------------------
# 切换协议: 用 sed 修改 defaults.yaml 中 step1.protocol 行
# ---------------------------------------------------------------------------
set_protocol() {
    local proto="$1"
    sed -i "s/^  protocol:.*/  protocol: ${proto}/" "$CONFIG"
    # 验证写入成功
    local actual
    actual=$(grep '^  protocol:' "$CONFIG" | head -1 | awk '{print $2}')
    if [ "$actual" != "$proto" ]; then
        fail "协议切换失败: 期望 ${proto}, 实际 ${actual}"
        exit 1
    fi
}

# ---------------------------------------------------------------------------
# 运行一轮 S1 测试
# ---------------------------------------------------------------------------
run_one() {
    local proto="$1"
    local outdir="${OUTPUT_BASE}/${proto}"

    info "=========================================="
    info " 协议: ${proto}  输出: ${outdir}"
    info "=========================================="

    # 1. 切协议
    info "切换 step1.protocol → ${proto} ..."
    set_protocol "$proto"
    ok "defaults.yaml 已更新"

    local expected_handoff_mode=""
    local expected_effective_modes=""
    case "$proto" in
        ext|default)
            expected_handoff_mode="optimize_all_candidates"
            expected_effective_modes="optimize_all_candidates"
            ;;
        full)
            expected_handoff_mode="optimize_all_survivors_within_window"
            expected_effective_modes="optimize_all_survivors_within_window"
            ;;
        lite)
            expected_handoff_mode="optimize_rank1"
            expected_effective_modes="optimize_rank1,optimize_top2_if_gap_small"
            ;;
        zero)
            expected_handoff_mode="optimize_rank1"
            expected_effective_modes="optimize_rank1,optimize_all_within_0p5_kcal"
            ;;
        *) expected_handoff_mode="" ;;
    esac

    if [ -d "$outdir" ]; then
        if [ "$FORCE_CLEAN" = "1" ]; then
            warn "输出目录已存在，按 FORCE_CLEAN=1 自动删除: ${outdir}"
            rm -rf "$outdir"
            ok "已删除旧目录"
        else
            warn "输出目录已存在: ${outdir}"
            warn "保留现有目录并尝试断点续算 (FORCE_CLEAN=0)"
        fi
    fi

    # 3. 运行
    info "启动 S1-only 运行 ..."
    echo ""
    cd "$PROJECT_ROOT"
    $CLI \
        --smiles "$PRODUCT_SMILES" \
        --output "$outdir" \
        --reaction-type "$REACTION_TYPE" \
        --skip-steps s2,s3,s4 \
        --log-level INFO \
    || {
        fail "协议 ${proto} 运行失败 (退出码 $?)"
        return 1
    }
    echo ""

    # 4. 验证
    info "验证输出 ..."
    local errors=0

    # 4a. provenance.json 存在
    local provenance
    provenance=$(find "$outdir" -name "provenance.json" -path "*/S1_ConfGeneration/*" 2>/dev/null | head -1)
    if [ -z "$provenance" ]; then
        fail "provenance.json 未找到"
        errors=$((errors + 1))
    else
        ok "provenance.json: $provenance"

        # 4b. 协议字段正确
        local proto_in_json
        proto_in_json=$(python3 -c "
import json, sys
with open('$provenance') as f:
    d = json.load(f)
print(d.get('protocol', 'MISSING'))
" 2>/dev/null || echo "PARSE_ERROR")

        if [ "$proto_in_json" = "$proto" ]; then
            ok "provenance.protocol = ${proto_in_json}"
        else
            fail "provenance.protocol = ${proto_in_json}, 期望 ${proto}"
            errors=$((errors + 1))
        fi

        # 4c. 关键字段摘要
        info "--- provenance 摘要 ---"
        python3 -c "
import json
with open('$provenance') as f:
    d = json.load(f)

fields = [
    ('protocol',                  'protocol'),
    ('has_geometry_optimization', 'has_geom_opt'),
]
cs = d.get('conformer_search', {})
fields += [
    ('two_stage_enabled', 'two_stage'),
    ('ngeom_default',     'ngeom_def'),
    ('ngeom_max',         'ngeom_max'),
]
fos = d.get('final_opt_sp', {})
fields += [
    ('freq_requested',     'freq'),
    ('final_sp_requested', 'final_sp'),
    ('selection_mode',     'sel_mode'),
    ('selected_candidate_count', 'sel_count'),
]
funnel = d.get('funnel', {})
handoff = d.get('handoff', {})
fields += [
    ('search_mode', 'funnel_mode'),
    ('ranking_basis', 'funnel_basis'),
    ('candidate_total', 'funnel_total'),
    ('survivor_count', 'funnel_surv'),
    ('mode', 'handoff_mode'),
    ('mode_requested', 'handoff_req'),
    ('mode_effective', 'handoff_eff'),
    ('selected_candidate_count', 'handoff_sel'),
]

for key, label in fields:
    src = d if key in d else cs if key in cs else fos if key in fos else funnel if key in funnel else handoff
    val = src.get(key, 'N/A')
    print(f'  {label:16s} = {val}')
" 2>/dev/null || warn "provenance 解析失败"
        info "--- end ---"

        local has_funnel has_handoff
        has_funnel=$(python3 -c "import json; d=json.load(open('$provenance')); print('ok' if isinstance(d.get('funnel'), dict) else 'bad')" 2>/dev/null || echo "bad")
        has_handoff=$(python3 -c "import json; d=json.load(open('$provenance')); print('ok' if isinstance(d.get('handoff'), dict) else 'bad')" 2>/dev/null || echo "bad")
        if [ "$has_funnel" = "ok" ] && [ "$has_handoff" = "ok" ]; then
            ok "funnel/handoff 元数据存在"
        else
            fail "funnel/handoff 元数据缺失"
            errors=$((errors + 1))
        fi

        if [ -n "$expected_handoff_mode" ]; then
            local handoff_mode_requested handoff_mode_effective
            handoff_mode_requested=$(python3 -c "import json; d=json.load(open('$provenance')); h=d.get('handoff',{}); print(h.get('mode_requested', h.get('mode','MISSING')))" 2>/dev/null || echo "PARSE_ERROR")
            handoff_mode_effective=$(python3 -c "import json; d=json.load(open('$provenance')); h=d.get('handoff',{}); print(h.get('mode_effective', h.get('mode','MISSING')))" 2>/dev/null || echo "PARSE_ERROR")
            if [ "$handoff_mode_requested" = "$expected_handoff_mode" ]; then
                ok "handoff.mode_requested = ${handoff_mode_requested}"
            else
                fail "handoff.mode_requested = ${handoff_mode_requested}, 期望 ${expected_handoff_mode}"
                errors=$((errors + 1))
            fi
            if python3 - <<'PY' "$handoff_mode_effective" "$expected_effective_modes"
import sys
actual = sys.argv[1]
allowed = [item for item in sys.argv[2].split(',') if item]
raise SystemExit(0 if actual in allowed else 1)
PY
            then
                ok "handoff.mode_effective = ${handoff_mode_effective}"
            else
                fail "handoff.mode_effective = ${handoff_mode_effective}, 允许值 ${expected_effective_modes}"
                errors=$((errors + 1))
            fi
        fi
    fi

    # 4d. product_min.xyz 存在
    local product_xyz
    product_xyz=$(find "$outdir" \( -name "product_min.xyz" -o -name "product_global_min.xyz" \) 2>/dev/null | head -1)
    if [ -n "$product_xyz" ]; then
        ok "产物 XYZ: $product_xyz"
    else
        fail "产物 XYZ 未找到"
        errors=$((errors + 1))
    fi

    # 4e. S2/S3/S4 目录不存在 (确认 skip 生效)
    local s2_dir s3_dir s4_dir
    s2_dir=$(find "$outdir" -maxdepth 2 -type d -name "S2_Retro" 2>/dev/null | head -1)
    s3_dir=$(find "$outdir" -maxdepth 2 -type d -name "S3_*" 2>/dev/null | head -1)
    s4_dir=$(find "$outdir" -maxdepth 2 -type d -name "S4_Data" 2>/dev/null | head -1)
    if [ -z "$s2_dir" ] && [ -z "$s3_dir" ] && [ -z "$s4_dir" ]; then
        ok "S2/S3/S4 未创建 (skip 生效)"
    else
        fail "发现不应存在的目录: S2=$s2_dir S3=$s3_dir S4=$s4_dir"
        errors=$((errors + 1))
    fi

    # 5. 结果
    if [ "$errors" -eq 0 ]; then
        ok "协议 ${proto} 全部验证通过 ✓"
    else
        fail "协议 ${proto} 有 ${errors} 个验证失败"
    fi

    return $errors
}

# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------
info "ReactionProfileHunter S1 协议测试"
info "产物 SMILES: ${PRODUCT_SMILES}"
info "反应类型:    ${REACTION_TYPE}"
info "测试协议:    ${ALL_PROTOCOLS[*]}"
info "项目根目录:  ${PROJECT_ROOT}"
echo ""

# 备份原始配置
local_orig_protocol=$(grep '^  protocol:' "$CONFIG" | head -1 | awk '{print $2}')
info "当前协议: ${local_orig_protocol} (测试结束后会恢复)"
echo ""

TOTAL_ERRORS=0

for proto in "${ALL_PROTOCOLS[@]}"; do
    if run_one "$proto"; then
        :
    else
        TOTAL_ERRORS=$((TOTAL_ERRORS + 1))
    fi
    echo ""
    echo ""
done

# 恢复原始协议
info "恢复 step1.protocol → ${local_orig_protocol} ..."
set_protocol "$local_orig_protocol"
ok "配置已恢复"

# ---------------------------------------------------------------------------
# 总结
# ---------------------------------------------------------------------------
echo ""
info "=========================================="
info " 测试总结"
info "=========================================="
if [ "$TOTAL_ERRORS" -eq 0 ]; then
    ok "全部 ${#ALL_PROTOCOLS[@]} 个协议测试通过 ✓"
else
    fail "${TOTAL_ERRORS} / ${#ALL_PROTOCOLS[@]} 个协议测试失败"
    exit 1
fi
