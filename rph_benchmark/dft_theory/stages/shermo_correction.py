#!/usr/bin/env python3
"""
Shermo 热力学校正阶段 (Shermo Thermodynamic Correction Stage)
=============================================================

对 Phase 1 SP benchmark 结果应用 Shermo 热力学校正，
计算 ΔG‡, ΔH‡, ΔG_rxn, ΔH_rxn 等热力学量。

核心原理:
  G_correction = G_sum(Shermo输出) − E_el(Shermo输入)
  校正量仅依赖于振动频率，与电子能方法无关。
  因此只需提取一次校正量，即可应用于所有 benchmark 方法。

工作流:
  1. 从 baseline RPH 输出解析已有 Shermo .sum 文件（提取校正量）
  2. 若无 .sum 文件，从频率 log + baseline SP 能量运行 Shermo
  3. 对每个 benchmark 方法: G = E_SP + G_correction
  4. 计算 ΔG‡, ΔH‡, ΔG_rxn, ΔH_rxn
  5. 评估 MAE vs 参考方法，选 Winner，生成报告

用法:
  # 计算 + 评估
  python benchmark/dft_theory/stages/shermo_correction.py \
      --session-dir .../session_bl_fixed \
      --mode all --rx-id 1 3 8 15

  # 仅计算
  python benchmark/dft_theory/stages/shermo_correction.py \
      --session-dir .../session_bl_fixed --mode run --rx-id 1 3 8 15

  # 仅评估（使用已有 shermo_result.json）
  python benchmark/dft_theory/stages/shermo_correction.py \
      --session-dir .../session_bl_fixed --mode evaluate --rx-id 1 3 8 15
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

# ---------------------------------------------------------------------------
# 项目内依赖 (复用 RPH 核心模块)
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from rph_core.utils.thermo import parse_shermo_sum
from rph_core.utils.shermo_runner import run_shermo
from rph_core.utils.constants import HARTREE_TO_KCAL

from benchmark.dft_theory.lib.state_manager import (
    append_warning,
    init_rx_manifest,
    read_benchmark_manifest,
    read_rx_manifest,
    record_task_status,
    resolve_rx_dir,
    resolve_rx_phase_report_paths,
    set_phase_winner,
    set_stage_status,
    write_benchmark_manifest,
    write_rx_manifest,
    DEFAULT_BENCHMARK_STAGES,
    BASELINE_SCHEMA_VERSION,
)

# ---------------------------------------------------------------------------
# 常量
# ---------------------------------------------------------------------------
SHERMO_RESULT_SCHEMA_VERSION = "benchmark_shermo_v1"
SHERMO_CANONICAL_POINTS = ("precursor", "intermediate", "ts", "product")
SHERMO_PHASE_KEY = "shermo_correction"
DEFAULT_REFERENCE_METHOD = "SP-REF"
DEFAULT_METHOD_IDS = ["SP-REF", "SP-1", "SP-2", "SP-3", "SP-4", "SP-5", "SP-6"]


# ---------------------------------------------------------------------------
# 数据结构
# ---------------------------------------------------------------------------

@dataclass
class ShermoCorrectionData:
    """单个驻点的热力学校正量 (Hartree)"""
    point_label: str
    g_correction: float
    h_correction: float
    u_correction: float
    s_total: Optional[float] = None       # cal/mol·K
    source_sum_file: Optional[Path] = None
    e_input: float = 0.0


@dataclass
class ThermoMethodResult:
    """单个 benchmark 方法在单个反应上的热力学结果 (kcal/mol)"""
    method_id: str
    delta_G_activation: float
    delta_H_activation: float
    delta_G_reaction: float
    delta_H_reaction: float
    g_ts: float
    g_intermediate: float
    g_product: float
    g_precursor: float
    h_ts: float
    h_intermediate: float
    h_product: float
    h_precursor: float


# ---------------------------------------------------------------------------
# Step 0: 路径解析
# ---------------------------------------------------------------------------

def _resolve_baseline_dir(session_dir: Path, rx_id: str) -> Optional[Path]:
    """解析 baseline 输出目录。

    优先使用 session 中的 bl_fixed，其次使用最近修改的 bl_* 目录。
    """
    baselines_root = session_dir.parent.parent / "baselines"
    rx_baselines = baselines_root / f"rx{rx_id}"

    # 优先 bl_fixed
    fixed = rx_baselines / "bl_fixed"
    if fixed.exists():
        return fixed

    # 退而求其次
    candidates = sorted(rx_baselines.glob("bl_*"), key=lambda p: p.stat().st_mtime, reverse=True)
    return candidates[0] if candidates else None


def _resolve_freq_outputs(baseline_dir: Path) -> Dict[str, Path]:
    """查找 baseline 中各驻点的频率输出文件。

    Returns:
        {"precursor": Path, "intermediate": Path, "ts": Path, "product": Path}
        缺失项不包含在返回字典中。
    """
    paths: Dict[str, Path] = {}

    # S1 频率输出 (product + precursor)
    s1 = baseline_dir / "S1_ConfGeneration"
    if s1.exists():
        # product: 取 conf_000 的优化 log（含频率）
        product_logs = sorted((s1 / "product" / "dft").glob("conf_000.log"))
        if product_logs:
            paths["product"] = product_logs[0]
        # precursor
        precursor_logs = sorted((s1 / "precursor" / "dft").glob("conf_000.log"))
        if precursor_logs:
            paths["precursor"] = precursor_logs[0]

    # S3 频率输出 (TS + intermediate)
    s3 = baseline_dir / "S3_TS"
    if s3.exists():
        # TS: berny optimization (含 freq=calcfc)
        ts_log = s3 / "ts_opt" / "berny" / "ts_guess.log"
        if ts_log.exists():
            paths["ts"] = ts_log
        # intermediate: standard optimization (含 freq)
        int_log = s3 / "S3_intermediate_opt" / "standard" / "intermediate.log"
        if int_log.exists():
            paths["intermediate"] = int_log

    return paths


def _resolve_shermo_sum_files(baseline_dir: Path) -> Dict[str, Path]:
    """查找 baseline 中已有的 Shermo .sum 文件。

    按优先级搜索: S3/shermo/ → S1/*/dft/
    """
    paths: Dict[str, Path] = {}

    # S3 Shermo 目录 (TS, reactant/intermediate)
    s3_shermo = baseline_dir / "S3_TS" / "shermo"
    if s3_shermo.exists():
        ts_sum = s3_shermo / "ts_Shermo.sum"
        if ts_sum.exists():
            paths["ts"] = ts_sum
        react_sum = s3_shermo / "reactant_Shermo.sum"
        if react_sum.exists():
            paths["intermediate"] = react_sum

    # S1 product + precursor
    s1 = baseline_dir / "S1_ConfGeneration"
    if s1.exists():
        for label, subdir in [("product", "product"), ("precursor", "precursor")]:
            dft_dir = s1 / subdir / "dft"
            if dft_dir.exists():
                sums = sorted(dft_dir.glob("conf_000_Shermo.sum"))
                if sums:
                    paths[label] = sums[0]

    return paths


def _read_e_input_from_sum(sum_file: Path) -> float:
    """从 Shermo .sum 文件读取输入的电子能量 (Electronic energy line)."""
    content = sum_file.read_text(errors="ignore")
    for line in content.splitlines():
        if "Electronic energy:" in line:
            tokens = line.strip().split()
            for token in reversed(tokens):
                try:
                    return float(token)
                except ValueError:
                    continue
    raise ValueError(f"Cannot read Electronic energy from {sum_file}")


# ---------------------------------------------------------------------------
# Step 1: 校正量提取
# ---------------------------------------------------------------------------

def _extract_corrections_from_sum(
    sum_file: Path,
    point_label: str,
) -> ShermoCorrectionData:
    """从已有 Shermo .sum 文件提取热力学校正量。

    Args:
        sum_file: Shermo .sum 文件路径
        point_label: 驻点标签

    Returns:
        ShermoCorrectionData 包含 G/H/U 校正量
    """
    record = parse_shermo_sum(sum_file)
    e_input = _read_e_input_from_sum(sum_file)

    return ShermoCorrectionData(
        point_label=point_label,
        g_correction=record.g_sum_hartree - e_input,
        h_correction=record.h_sum_hartree - e_input,
        u_correction=record.u_sum_hartree - e_input,
        s_total=record.s_cal_mol_k,
        source_sum_file=sum_file,
        e_input=e_input,
    )


def _run_shermo_for_point(
    freq_log: Path,
    sp_energy_hartree: float,
    output_dir: Path,
    config: dict,
    point_label: str,
) -> ShermoCorrectionData:
    """从频率 log 运行 Shermo 获取热力学校正量。

    用于无已有 .sum 文件时的 fallback。

    Args:
        freq_log: 频率计算 log 文件
        sp_energy_hartree: baseline SP 能量
        output_dir: 输出目录（存放 .sum）
        config: RPH 全局配置
        point_label: 驻点标签

    Returns:
        ShermoCorrectionData
    """
    from rph_core.utils.config_loader import resolve_executable_config

    shermo_bin_raw = config.get("executables", {}).get("shermo", {}).get("path", "Shermo")
    shermo_bin = Path(resolve_executable_config("shermo", config)["path"]
                      if callable(resolve_executable_config) else shermo_bin_raw)

    thermo_cfg = config.get("thermo", {})
    temperature_k = thermo_cfg.get("temperature_k", 298.15)
    pressure_atm = thermo_cfg.get("pressure_atm", 1.0)
    scl_zpe = thermo_cfg.get("scl_zpe", 0.9905)
    ilowfreq = thermo_cfg.get("ilowfreq", 2)
    imagreal = thermo_cfg.get("imagreal")

    output_dir.mkdir(parents=True, exist_ok=True)
    sum_file = output_dir / f"{point_label}_Shermo.sum"

    result = run_shermo(
        shermo_bin=Path(shermo_bin),
        freq_output=freq_log,
        sp_energy=sp_energy_hartree,
        output_file=sum_file,
        temperature_k=temperature_k,
        pressure_atm=pressure_atm,
        scl_zpe=scl_zpe,
        ilowfreq=ilowfreq,
        imagreal=imagreal,
    )

    return ShermoCorrectionData(
        point_label=point_label,
        g_correction=result.g_sum - sp_energy_hartree,
        h_correction=result.h_sum - sp_energy_hartree,
        u_correction=result.u_sum - sp_energy_hartree,
        s_total=result.s_total,
        source_sum_file=sum_file,
        e_input=sp_energy_hartree,
    )


def _resolve_all_corrections(
    baseline_dir: Path,
    config: dict,
    shermo_output_dir: Path,
) -> Dict[str, ShermoCorrectionData]:
    """获取所有驻点的热力学校正量 (优先已有 .sum，fallback 运行 Shermo)。

    Returns:
        {"precursor": ..., "intermediate": ..., "ts": ..., "product": ...}
    """
    corrections: Dict[str, ShermoCorrectionData] = {}
    sum_files = _resolve_shermo_sum_files(baseline_dir)
    freq_files = _resolve_freq_outputs(baseline_dir)

    for point_label in SHERMO_CANONICAL_POINTS:
        sum_file = sum_files.get(point_label)
        if sum_file:
            # 从已有 .sum 提取
            try:
                corrections[point_label] = _extract_corrections_from_sum(sum_file, point_label)
            except Exception as exc:
                print(f"  [WARN] 解析 {sum_file} 失败: {exc}，尝试运行 Shermo")
                sum_file = None

        if point_label not in corrections:
            freq_log = freq_files.get(point_label)
            if freq_log:
                # Fallback: 运行 Shermo
                # 使用 baseline SP 能量作为 e_input（任意值均可，校正量与之无关）
                e_sp = _read_baseline_energy(baseline_dir, point_label)
                try:
                    corrections[point_label] = _run_shermo_for_point(
                        freq_log=freq_log,
                        sp_energy_hartree=e_sp,
                        output_dir=shermo_output_dir,
                        config=config,
                        point_label=point_label,
                    )
                except Exception as exc:
                    print(f"  [ERROR] 运行 Shermo 失败 ({point_label}): {exc}")
            else:
                print(f"  [WARN] 缺少 {point_label} 的频率文件，跳过")

    return corrections


def _read_baseline_energy(baseline_dir: Path, point_label: str) -> float:
    """从 baseline thermo_reference.json 读取 SP 能量。

    用于 Shermo fallback 的 e_input。
    """
    tr_path = baseline_dir / "thermo_reference.json"
    if tr_path.exists():
        tr = json.loads(tr_path.read_text())
        key_map = {
            "ts": "e_ts",
            "intermediate": "e_reactant",
            "product": "e_product",
        }
        key = key_map.get(point_label)
        if key and key in tr:
            return float(tr[key])

    # Fallback: 从 S3 SP output 读取
    s3 = baseline_dir / "S3_TS"
    sp_map = {
        "ts": s3 / "ts_opt" / "L2_SP" / "ts_guess_ts_*.out",
        "intermediate": s3 / "S3_intermediate_opt" / "L2_SP" / "intermediate_opt_*.out",
    }
    pattern = sp_map.get(point_label)
    if pattern:
        candidates = sorted(pattern.parent.glob(pattern.name))
        if candidates:
            return _extract_energy_from_orca_out(candidates[0])

    return 0.0


def _extract_energy_from_orca_out(out_file: Path) -> float:
    """从 ORCA SP 输出提取最终能量。"""
    content = out_file.read_text(errors="ignore")
    for line in reversed(content.splitlines()):
        if "FINAL SINGLE POINT ENERGY" in line:
            tokens = line.strip().split()
            try:
                return float(tokens[-1])
            except (ValueError, IndexError):
                pass
    return 0.0


# ---------------------------------------------------------------------------
# Step 2: 热力学量计算 (G_method = E_SP_method + G_correction)
# ---------------------------------------------------------------------------

def _read_sp_energies(session_dir: Path, rx_id: str, method_id: str) -> Optional[Dict[str, float]]:
    """从 sp_result.json 读取 SP 能量。"""
    result_file = resolve_rx_dir(session_dir, rx_id) / "sp" / method_id / "sp_result.json"
    if not result_file.exists():
        return None

    data = json.loads(result_file.read_text())
    return {
        "precursor": data.get("e_precursor_hartree"),
        "intermediate": data.get("e_intermediate_hartree"),
        "ts": data.get("e_ts_hartree"),
        "product": data.get("e_product_hartree"),
    }


def _compute_thermo_for_method(
    method_id: str,
    sp_energies: Dict[str, Optional[float]],
    corrections: Dict[str, ShermoCorrectionData],
) -> ThermoMethodResult:
    """为单个方法计算完整热力学量。

    G_X = E_SP_X + G_correction_X  (Hartree → kcal/mol via HARTREE_TO_KCAL)
    """

    def _g(label: str) -> Optional[float]:
        e = sp_energies.get(label)
        c = corrections.get(label)
        if e is not None and c is not None:
            return (e + c.g_correction) * HARTREE_TO_KCAL
        return None

    def _h(label: str) -> Optional[float]:
        e = sp_energies.get(label)
        c = corrections.get(label)
        if e is not None and c is not None:
            return (e + c.h_correction) * HARTREE_TO_KCAL
        return None

    g_ts = _g("ts")
    g_int = _g("intermediate")
    g_prod = _g("product")
    g_prec = _g("precursor")
    h_ts = _h("ts")
    h_int = _h("intermediate")
    h_prod = _h("product")
    h_prec = _h("precursor")

    return ThermoMethodResult(
        method_id=method_id,
        delta_G_activation=(g_ts - g_int) if (g_ts is not None and g_int is not None) else float("nan"),
        delta_H_activation=(h_ts - h_int) if (h_ts is not None and h_int is not None) else float("nan"),
        delta_G_reaction=(g_prod - g_int) if (g_prod is not None and g_int is not None) else float("nan"),
        delta_H_reaction=(h_prod - h_int) if (h_prod is not None and h_int is not None) else float("nan"),
        g_ts=g_ts or float("nan"),
        g_intermediate=g_int or float("nan"),
        g_product=g_prod or float("nan"),
        g_precursor=g_prec or float("nan"),
        h_ts=h_ts or float("nan"),
        h_intermediate=h_int or float("nan"),
        h_product=h_prod or float("nan"),
        h_precursor=h_prec or float("nan"),
    )


# ---------------------------------------------------------------------------
# Step 3: 运行阶段
# ---------------------------------------------------------------------------

def _run_shermo_correction_for_rx(
    session_dir: Path,
    rx_id: str,
    method_ids: List[str],
    config: dict,
) -> Dict[str, Any]:
    """为单个反应运行 Shermo 热力学校正。

    Returns:
        结果字典，写入 shermo_result.json
    """
    bl_dir = _resolve_baseline_dir(session_dir, rx_id)
    if not bl_dir:
        return {"status": "failed", "error": f"No baseline directory found for rx{rx_id}"}

    shermo_out = resolve_rx_dir(session_dir, rx_id) / "shermo"
    shermo_out.mkdir(parents=True, exist_ok=True)

    # 1. 获取热力学校正量
    corrections = _resolve_all_corrections(bl_dir, config, shermo_out)
    missing = [p for p in SHERMO_CANONICAL_POINTS if p not in corrections]
    if missing:
        return {"status": "partial_failed", "error": f"Missing corrections: {missing}"}

    # 2. 对每个方法计算热力学量
    method_results: Dict[str, Dict[str, Any]] = {}
    for mid in method_ids:
        sp_energies = _read_sp_energies(session_dir, rx_id, mid)
        if sp_energies is None:
            print(f"  [SKIP] rx{rx_id}/{mid}: sp_result.json not found")
            continue

        thermo = _compute_thermo_for_method(mid, sp_energies, corrections)
        method_results[mid] = {
            "method_id": thermo.method_id,
            "delta_G_activation_kcal": thermo.delta_G_activation,
            "delta_H_activation_kcal": thermo.delta_H_activation,
            "delta_G_reaction_kcal": thermo.delta_G_reaction,
            "delta_H_reaction_kcal": thermo.delta_H_reaction,
            "g_ts_kcal": thermo.g_ts,
            "g_intermediate_kcal": thermo.g_intermediate,
            "g_product_kcal": thermo.g_product,
            "g_precursor_kcal": thermo.g_precursor,
            "h_ts_kcal": thermo.h_ts,
            "h_intermediate_kcal": thermo.h_intermediate,
            "h_product_kcal": thermo.h_product,
            "h_precursor_kcal": thermo.h_precursor,
        }

    # 3. 构建结果文件
    result = {
        "schema_version": SHERMO_RESULT_SCHEMA_VERSION,
        "rx_id": rx_id,
        "baseline_dir": str(bl_dir),
        "corrections": {
            label: {
                "g_correction_hartree": c.g_correction,
                "h_correction_hartree": c.h_correction,
                "u_correction_hartree": c.u_correction,
                "s_total_cal": c.s_total,
                "source_sum_file": str(c.source_sum_file) if c.source_sum_file else None,
            }
            for label, c in corrections.items()
        },
        "method_results": method_results,
        "status": "completed",
        "created_at": datetime.now(timezone.utc).isoformat(),
    }

    result_path = shermo_out / "shermo_result.json"
    result_path.write_text(json.dumps(result, indent=2, ensure_ascii=False))

    # 4. 更新 manifest
    record_task_status(
        session_dir, rx_id, "shermo", "correction",
        status="completed",
        run_dir=shermo_out,
        extra={"result_json": str(result_path)},
    )

    return result


def run_shermo_correction_stage(
    session_dir: Path,
    rx_ids: List[str],
    methods: Optional[List[str]] = None,
    config: Optional[dict] = None,
) -> Dict[str, Any]:
    """Shermo 热力学校正阶段 — 计算入口。

    对每个反应读取已有 Shermo .sum 文件提取热力学校正量，
    然后对每个 benchmark 方法计算 G/H 校正后的热力学量。
    """
    if methods is None:
        methods = DEFAULT_METHOD_IDS
    if config is None:
        config = _load_default_config()

    set_stage_status(session_dir, SHERMO_PHASE_KEY, "running")

    results: Dict[str, Any] = {}
    for rx_id in rx_ids:
        print(f"\n{'='*60}")
        print(f"  Shermo Correction — rx{rx_id}")
        print(f"{'='*60}")
        results[rx_id] = _run_shermo_correction_for_rx(session_dir, rx_id, methods, config)

    # 汇总状态
    completed = sum(1 for v in results.values() if v.get("status") == "completed")
    failed = sum(1 for v in results.values() if v.get("status") != "completed")

    if failed == 0:
        set_stage_status(session_dir, SHERMO_PHASE_KEY, "computed")
    else:
        set_stage_status(session_dir, SHERMO_PHASE_KEY, "partial_failed")

    return {
        "stage": "shermo_correction",
        "total": len(rx_ids),
        "completed": completed,
        "failed": failed,
        "results": results,
    }


# ---------------------------------------------------------------------------
# Step 4: 评估阶段
# ---------------------------------------------------------------------------

def _read_shermo_result(session_dir: Path, rx_id: str) -> Optional[Dict[str, Any]]:
    """读取 shermo_result.json。"""
    path = resolve_rx_dir(session_dir, rx_id) / "shermo" / "shermo_result.json"
    if not path.exists():
        return None
    return json.loads(path.read_text())


def _annotate_shermo_scores(
    method_results: Dict[str, Dict[str, Any]],
    reference_id: str,
) -> Dict[str, Dict[str, float]]:
    """计算各方法 vs 参考方法的 MAE (ΔG‡ + ΔG_rxn)。"""
    ref = method_results.get(reference_id)
    if not ref:
        return {}

    scores: Dict[str, Dict[str, float]] = {}
    for mid, m in method_results.items():
        if mid == reference_id:
            continue
        mae_act = abs(m["delta_G_activation_kcal"] - ref["delta_G_activation_kcal"])
        mae_rxn = abs(m["delta_G_reaction_kcal"] - ref["delta_G_reaction_kcal"])
        scores[mid] = {
            "mae_activation": round(mae_act, 4),
            "mae_reaction": round(mae_rxn, 4),
            "mae_combined": round((mae_act + mae_rxn) / 2, 4),
        }
    return scores


def _choose_shermo_winner(
    scores: Dict[str, Dict[str, float]],
    method_ids: List[str],
) -> Optional[str]:
    """选择 MAE 最小的方法作为 Winner。"""
    candidates = {mid: s["mae_combined"] for mid, s in scores.items()
                  if mid in method_ids}
    if not candidates:
        return None
    return min(candidates, key=candidates.get)


def _write_shermo_rx_report(
    session_dir: Path,
    rx_id: str,
    method_ids: List[str],
    shermo_result: Dict[str, Any],
    scores: Dict[str, Dict[str, float]],
    winner: str,
    reference_id: str,
) -> Tuple[Path, Path]:
    """写 per-rx Shermo 报告 (JSON + Markdown)。"""
    reports = resolve_rx_phase_report_paths(session_dir, rx_id, "shermo")
    reports["json"].parent.mkdir(parents=True, exist_ok=True)

    corrections = shermo_result.get("corrections", {})
    method_results = shermo_result.get("method_results", {})

    # JSON
    json_data = {
        "rx_id": rx_id,
        "phase": "shermo",
        "reference_method": reference_id,
        "winner": winner,
        "corrections": corrections,
        "scores": scores,
        "method_results": method_results,
        "created_at": datetime.now(timezone.utc).isoformat(),
    }
    reports["json"].write_text(json.dumps(json_data, indent=2, ensure_ascii=False))

    # Markdown
    lines = [
        f"# Shermo 热力学校正报告 — rx{rx_id}",
        "",
        f"**Winner**: {winner} | **参考方法**: {reference_id}",
        "",
        "## 热力学校正量 (Hartree)",
        "",
        "| 驻点 | G_correction | H_correction | U_correction | S_total (cal/mol·K) |",
        "|------|-------------|-------------|-------------|---------------------|",
    ]
    for label in SHERMO_CANONICAL_POINTS:
        c = corrections.get(label, {})
        if c:
            lines.append(
                f"| {label} | {c.get('g_correction_hartree', 0):+.6f} "
                f"| {c.get('h_correction_hartree', 0):+.6f} "
                f"| {c.get('u_correction_hartree', 0):+.6f} "
                f"| {c.get('s_total_cal', 0):.2f} |"
            )

    lines.extend([
        "",
        "## 热力学量 (kcal/mol)",
        "",
        "| Method | ΔG‡ | ΔH‡ | ΔG_rxn | ΔH_rxn | MAE (combined) |",
        "|--------|-----|-----|--------|--------|----------------|",
    ])
    for mid in method_ids:
        m = method_results.get(mid, {})
        if m:
            s = scores.get(mid, {})
            mae_str = f"{s.get('mae_combined', 0):.2f}" if s else "—"
            lines.append(
                f"| {mid} "
                f"| {m.get('delta_G_activation_kcal', float('nan')):+.2f} "
                f"| {m.get('delta_H_activation_kcal', float('nan')):+.2f} "
                f"| {m.get('delta_G_reaction_kcal', float('nan')):+.2f} "
                f"| {m.get('delta_H_reaction_kcal', float('nan')):+.2f} "
                f"| {mae_str} |"
            )

    reports["md"].write_text("\n".join(lines) + "\n")
    return reports["json"], reports["md"]


def _write_shermo_global_report(
    session_dir: Path,
    rx_ids: List[str],
    method_ids: List[str],
    votes: Dict[str, int],
    global_winner: str,
    reference_id: str,
    all_scores: Dict[str, Dict[str, Dict[str, float]]],
) -> None:
    """写全局 Shermo 报告。"""
    reports_dir = session_dir / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)

    # JSON
    json_path = reports_dir / "shermo_global_summary.json"
    json_data = {
        "phase": "shermo",
        "reference_method": reference_id,
        "winner": global_winner,
        "votes": votes,
        "created_at": datetime.now(timezone.utc).isoformat(),
    }
    json_path.write_text(json.dumps(json_data, indent=2, ensure_ascii=False))

    # Markdown
    md_path = reports_dir / "shermo_global_summary.md"
    lines = [
        "# Shermo 热力学校正 — 全局报告",
        "",
        f"**全局 Winner**: {global_winner} | **参考方法**: {reference_id}",
        "",
        "## 投票结果",
        "",
        "| 方法 | 得票 |",
        "|------|------|",
    ]
    for mid in sorted(votes, key=votes.get, reverse=True):
        lines.append(f"| {mid} | {votes[mid]}/{len(rx_ids)} |")

    # 综合 MAE 排名 (跨反应平均)
    lines.extend([
        "",
        "## 综合 MAE 排名 (vs SP-REF, ΔG‡ + ΔG_rxn)",
        "",
        "| 方法 | avg MAE_activation | avg MAE_reaction | avg MAE_combined |",
        "|------|--------------------|------------------|------------------|",
    ])

    avg_scores: Dict[str, Dict[str, float]] = {}
    for mid in method_ids:
        if mid == reference_id:
            continue
        act_vals = []
        rxn_vals = []
        for rx_id in rx_ids:
            s = all_scores.get(rx_id, {}).get(mid, {})
            if s:
                act_vals.append(s["mae_activation"])
                rxn_vals.append(s["mae_reaction"])
        if act_vals:
            avg_scores[mid] = {
                "avg_mae_activation": sum(act_vals) / len(act_vals),
                "avg_mae_reaction": sum(rxn_vals) / len(rxn_vals),
            }

    for mid in sorted(avg_scores, key=lambda m: (
            avg_scores[m]["avg_mae_activation"] + avg_scores[m]["avg_mae_reaction"]) / 2):
        s = avg_scores[mid]
        combined = (s["avg_mae_activation"] + s["avg_mae_reaction"]) / 2
        lines.append(
            f"| {mid} | {s['avg_mae_activation']:.2f} | {s['avg_mae_reaction']:.2f} | {combined:.2f} |"
        )

    md_path.write_text("\n".join(lines) + "\n")

    # 更新全局 manifest
    bm = read_benchmark_manifest(session_dir)
    bm["shermo_winner"] = global_winner
    reports = bm.get("reports", {})
    reports["shermo_global_summary_json"] = str(json_path)
    reports["shermo_global_summary_md"] = str(md_path)
    bm["reports"] = reports
    write_benchmark_manifest(session_dir, bm)


def evaluate_shermo_correction_stage(
    session_dir: Path,
    rx_ids: List[str],
    methods: Optional[List[str]] = None,
    reference_method: str = DEFAULT_REFERENCE_METHOD,
) -> Dict[str, Any]:
    """Shermo 热力学校正阶段 — 评估入口。

    投票选出各反应和全局 Winner，生成 per-rx + global 报告。
    """
    if methods is None:
        methods = DEFAULT_METHOD_IDS

    set_stage_status(session_dir, SHERMO_PHASE_KEY, "running")

    votes: Dict[str, int] = {}
    all_scores: Dict[str, Dict[str, Dict[str, float]]] = {}

    for rx_id in rx_ids:
        print(f"\n  Evaluating Shermo — rx{rx_id}")

        shermo_result = _read_shermo_result(session_dir, rx_id)
        if shermo_result is None:
            print(f"    [WARN] No shermo_result.json found")
            continue

        method_results = shermo_result.get("method_results", {})
        if reference_method not in method_results:
            print(f"    [WARN] Reference method {reference_method} not in results")
            continue

        # 评分
        scores = _annotate_shermo_scores(method_results, reference_method)
        all_scores[rx_id] = scores

        if not scores:
            print(f"    [WARN] No comparable methods")
            continue

        # 选 Winner
        winner = _choose_shermo_winner(scores, methods)
        if winner is None:
            continue

        print(f"    Winner: {winner}")

        # 更新 per-rx manifest
        set_phase_winner(session_dir, rx_id, "shermo", winner,
                        extra={"reference_method": reference_method})

        # 写 per-rx 报告
        json_rpt, md_rpt = _write_shermo_rx_report(
            session_dir, rx_id, methods, shermo_result, scores, winner, reference_method
        )

        # 更新 rx manifest 报告路径
        rx_manifest = init_rx_manifest(session_dir, rx_id)
        rx_reports = rx_manifest.get("reports", {})
        rx_reports["shermo_summary_json"] = str(json_rpt)
        rx_reports["shermo_summary_md"] = str(md_rpt)
        rx_manifest["reports"] = rx_reports
        write_rx_manifest(session_dir, rx_id, rx_manifest)

        votes[winner] = votes.get(winner, 0) + 1

    # 全局投票
    if votes:
        global_winner = max(votes, key=votes.get)
        _write_shermo_global_report(
            session_dir, rx_ids, methods, votes, global_winner,
            reference_method, all_scores,
        )
        set_stage_status(session_dir, SHERMO_PHASE_KEY, "done")
        print(f"\n  Global Shermo Winner: {global_winner}")
    else:
        set_stage_status(session_dir, SHERMO_PHASE_KEY, "failed")

    return {
        "stage": "shermo_evaluate",
        "winner": max(votes, key=votes.get) if votes else None,
        "votes": votes,
    }


# ---------------------------------------------------------------------------
# 辅助: 加载配置
# ---------------------------------------------------------------------------

def _load_default_config() -> dict:
    """加载 RPH 全局配置 (defaults.yaml)。"""
    import yaml
    defaults_path = REPO_ROOT / "config" / "defaults.yaml"
    if defaults_path.exists():
        with open(defaults_path) as f:
            return yaml.safe_load(f) or {}
    return {}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Shermo Thermodynamic Correction Stage",
    )
    parser.add_argument("--session-dir", required=True, type=Path, help="Benchmark session directory")
    parser.add_argument("--mode", default="all", choices=["run", "evaluate", "all"],
                       help="Stage mode (default: all)")
    parser.add_argument("--rx-id", nargs="+", required=True, help="Reaction IDs to process")
    parser.add_argument("--methods", nargs="*", default=None, help="Method IDs (default: all 7 SP methods)")
    parser.add_argument("--reference", default=DEFAULT_REFERENCE_METHOD,
                       help=f"Reference method for MAE (default: {DEFAULT_REFERENCE_METHOD})")
    parser.add_argument("--defaults-config", type=Path, help="Path to defaults.yaml")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> None:
    args = _parse_args(argv)

    config = _load_default_config()
    if args.defaults_config and Path(args.defaults_config).exists():
        import yaml
        with open(args.defaults_config) as f:
            config = yaml.safe_load(f) or {}

    method_ids = args.methods if args.methods else DEFAULT_METHOD_IDS
    rx_ids = list(args.rx_id)

    print(f"Shermo Correction Stage — mode={args.mode}, rx_ids={rx_ids}")
    print(f"  Session: {args.session_dir}")

    if args.mode in ("run", "all"):
        run_shermo_correction_stage(
            session_dir=Path(args.session_dir),
            rx_ids=rx_ids,
            methods=method_ids,
            config=config,
        )

    if args.mode in ("evaluate", "all"):
        evaluate_shermo_correction_stage(
            session_dir=Path(args.session_dir),
            rx_ids=rx_ids,
            methods=method_ids,
            reference_method=args.reference,
        )

    print("\nDone.")


if __name__ == "__main__":
    main()
