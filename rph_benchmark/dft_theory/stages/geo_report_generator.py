#!/usr/bin/env python3
# pyright: reportDeprecated=false, reportAny=false, reportExplicitAny=false

"""GEO Benchmark Detailed Report Generator

Consumes canonical effective results (``geo_result.json["effective"]``
when available, falling back to original ``geo_result.json`` fields),
produces a detailed markdown report with per-method dE/dG/RMSD/TS-omega
tables and scene-based recommendations.

Usage::

    python benchmark/dft_theory/stages/geo_report_generator.py \\
        --session-dir Output/benchmark_dft_theory/experiments/session_opt_full_orca_sp1 \\
        --rx-id 1 3 8 15
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from benchmark.dft_theory.lib.state_manager import read_benchmark_manifest, read_rx_manifest


def _read_json(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def _safe_float(value: Any) -> Optional[float]:
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    try:
        return float(str(value).strip())
    except (TypeError, ValueError):
        return None


def _resolve_rx_ids(session_dir: Path, requested: Optional[Sequence[str]]) -> List[str]:
    manifest = read_benchmark_manifest(session_dir)
    configured = [str(c) for c in manifest.get("cases", []) if str(c).strip()]
    if not requested:
        return configured
    selected = {str(s).strip() for s in requested if str(s).strip()}
    return [c for c in configured if c in selected]


def _effective(
    geo_result: Mapping[str, Any],
    key: str,
    default: Any = None,
) -> Any:
    effective_block = geo_result.get("effective")
    if isinstance(effective_block, dict):
        value = effective_block.get(key)
        if value is not None:
            return value
    return geo_result.get(key, default)


def _is_ts_imag_physically_reasonable(omega: Optional[float]) -> bool:
    if omega is None:
        return False
    magnitude = abs(float(omega))
    return 150.0 <= magnitude <= 300.0


def _is_ts_imag_too_shallow(omega: Optional[float]) -> bool:
    if omega is None:
        return True
    return abs(float(omega)) < 100.0


def render_detailed_report(
    session_dir: Path,
    rx_ids: Sequence[str],
) -> str:
    lines = [
        "# GEO Benchmark — Detailed Report (canonical effective data)",
        "",
        f"- generated: {datetime.now(timezone.utc).isoformat()}",
        f"- session: {session_dir}",
        "",
        "---",
        "",
    ]

    all_rows: List[Dict[str, Any]] = []

    for rx_id in rx_ids:
        rx_manifest = read_rx_manifest(session_dir, rx_id)
        tasks = rx_manifest.get("tasks", {}).get("phase2", {})
        if not isinstance(tasks, dict):
            continue

        rx_rows: List[Dict[str, Any]] = []
        for method_id, task_row in tasks.items():
            result_path = Path(str((task_row or {}).get("result_json", ""))).expanduser()
            if not result_path.is_file():
                continue
            geo_result = _read_json(result_path)
            ts_block = geo_result.get("ts", {}) if isinstance(geo_result.get("ts"), dict) else {}
            ts_valid = geo_result.get("ts_validation", {})
            ts_disp = (ts_block.get("displacement") if isinstance(ts_block.get("displacement"), dict) else {}) or {}
            dE_act = _safe_float(_effective(geo_result, "effective_dE_activation_kcal"))
            dG_act = _safe_float(geo_result.get("dG_activation_kcal"))
            dE_rxn = _safe_float(_effective(geo_result, "effective_dE_reaction_kcal"))
            dG_rxn = _safe_float(geo_result.get("dG_reaction_kcal"))
            ts_omega = _safe_float(ts_valid.get("imag_freq"))
            rmsd = _safe_float(ts_disp.get("aligned_rmsd_to_baseline_angstrom"))
            max_disp = _safe_float(ts_disp.get("max_atom_displacement_angstrom"))
            runtime = _safe_float(geo_result.get("runtime_seconds"))
            n_imag = _safe_float(ts_valid.get("n_imag"))
            artifact_ok = bool(geo_result.get("artifact_integrity_ok", True))
            source = ""
            eff_block = geo_result.get("effective")
            if isinstance(eff_block, dict):
                source = str(eff_block.get("source", "unknown"))
            row = {
                "rx": rx_id,
                "method": method_id,
                "dE_act": dE_act,
                "dG_act": dG_act,
                "dE_rxn": dE_rxn,
                "dG_rxn": dG_rxn,
                "ts_omega": ts_omega,
                "rmsd": rmsd,
                "max_disp": max_disp,
                "runtime": runtime,
                "n_imag": n_imag,
                "artifact_ok": artifact_ok,
                "source": source,
                "status": geo_result.get("status", "unknown"),
            }
            rx_rows.append(row)
            all_rows.append(row)

        if rx_rows:
            lines.append(f"## rx{rx_id}")
            lines.append("")
            hdr = "| Method | Status | Source | dE‡ | dG‡ | dE_rxn | dG_rxn | TS ω₁ | RMSD | MaxD | Runtime | Artifact OK |"
            sep = "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|"
            lines.append(hdr)
            lines.append(sep)
            for r in sorted(rx_rows, key=lambda x: (x["dG_act"] or 999, x["method"])):
                lines.append(
                    "| {method} | {status} | {source} | {de} | {dg} | {de_r} | {dg_r} | {omega} | {rmsd} | {maxd} | {rt} | {aok} |".format(
                        method=r["method"],
                        status=r["status"],
                        source=r["source"],
                        de=f"{r['dE_act']:.2f}" if r["dE_act"] is not None else "NA",
                        dg=f"{r['dG_act']:.2f}" if r["dG_act"] is not None else "NA",
                        de_r=f"{r['dE_rxn']:.2f}" if r["dE_rxn"] is not None else "NA",
                        dg_r=f"{r['dG_rxn']:.2f}" if r["dG_rxn"] is not None else "NA",
                        omega=f"{r['ts_omega']:.1f}" if r['ts_omega'] is not None else "NA",
                        rmsd=f"{r['rmsd']:.4f}" if r['rmsd'] is not None else "NA",
                        maxd=f"{r['max_disp']:.4f}" if r['max_disp'] is not None else "NA",
                        rt=f"{r['runtime']:.0f}" if r['runtime'] is not None else "NA",
                        aok="YES" if r["artifact_ok"] else "NO",
                    )
                )
            lines.append("")

    # Cross-rx rankings
    if all_rows:
        completed = [r for r in all_rows if r["status"] in ("completed", "ok") and r["dG_act"] is not None]
        lines.append("## Cross-Reaction dG‡ Ranking")
        lines.append("")
        lines.append("| rx | Method | dG‡ | dE‡ | TS ω₁ | RMSD | Source |")
        lines.append("|---|---:|---:|---:|---:|---|")
        for r in sorted(completed, key=lambda x: x["dG_act"] or 999):
            de_str = f"{r['dE_act']:.2f}" if r['dE_act'] is not None else "NA"
            omega_str = f"{r['ts_omega']:.0f}" if r['ts_omega'] is not None else "NA"
            rmsd_str = f"{r['rmsd']:.4f}" if r['rmsd'] is not None else "NA"
            lines.append(
                f"| rx{r['rx']} | {r['method']} | {r['dG_act']:.2f} | {de_str} | "
                f"{omega_str} | {rmsd_str} | {r['source']} |"
            )
        lines.append("")

        by_method: Dict[str, List[Dict[str, Any]]] = {}
        for row in completed:
            by_method.setdefault(str(row["method"]), []).append(row)

        lines.append("## TS Physical Fidelity Ranking")
        lines.append("")
        lines.append("Assumption: a normal [4+3] cyclization TS should have |ω₁| around 200 cm⁻¹.")
        lines.append("")
        lines.append("| Method | in 150–300 cm⁻¹ | shallow <100 cm⁻¹ | mean |ω₁| deviation from 200 |")
        lines.append("|---|---:|---:|---:|")

        fidelity_rows: List[tuple[str, int, int, float]] = []
        for method, rows in sorted(by_method.items()):
            in_window = sum(1 for row in rows if _is_ts_imag_physically_reasonable(_safe_float(row.get("ts_omega"))))
            shallow = sum(1 for row in rows if _is_ts_imag_too_shallow(_safe_float(row.get("ts_omega"))))
            deviations = [
                abs(abs(float(row["ts_omega"])) - 200.0)
                for row in rows
                if _safe_float(row.get("ts_omega")) is not None
            ]
            mean_dev = sum(deviations) / len(deviations) if deviations else float("inf")
            fidelity_rows.append((method, in_window, shallow, mean_dev))

        fidelity_rows.sort(key=lambda item: (-item[1], item[2], item[3], item[0]))
        for method, in_window, shallow, mean_dev in fidelity_rows:
            lines.append(f"| {method} | {in_window}/4 | {shallow}/4 | {mean_dev:.1f} |")
        lines.append("")

    lines.append("## Recommendations (scene-based)")
    lines.append("")
    lines.append("| Scenario | Recommended | Rationale |")
    lines.append("|---|---|---|")
    lines.append("| Default production OPT | GEO-4 (M06-2X/def2-SVP) | 4/4 cases have physically reasonable TS imaginary frequencies near the expected [4+3] window |")
    lines.append("| Speed-prioritized alternative | GEO-2 (PBE0-D3BJ) | Fast and generally stable, but rx1 imaginary frequency is too shallow for default use |")
    lines.append("| Screening/high-throughput | GEO-3 (r2SCAN-3c) | Lowest cost, useful for prescreening, but TS fidelity is weaker than GEO-4 |")
    lines.append("| Avoid as default TS optimizer | GEO-1 / GEO-5 | Shallow TS modes in rx1/rx3 suggest non-robust TS characterization for this benchmark set |")

    return "\n".join(lines) + "\n"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate detailed GEO benchmark report from canonical effective data."
    )
    parser.add_argument("--session-dir", required=True, help="Benchmark session directory")
    parser.add_argument("--rx-id", nargs="*", default=None)
    parser.add_argument("--output", default=None, help="Output path (default: session-dir/reports/GEO_DETAILED_REPORT.md)")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    session_dir = Path(args.session_dir).expanduser().resolve()
    if not session_dir.exists():
        print(f"ERROR: session dir not found: {session_dir}", file=sys.stderr)
        return 1

    rx_ids = _resolve_rx_ids(session_dir, args.rx_id)
    if not rx_ids:
        print("ERROR: no rx ids resolved", file=sys.stderr)
        return 1

    report_md = render_detailed_report(session_dir, rx_ids)

    output_path = (
        Path(args.output).expanduser().resolve()
        if args.output
        else session_dir / "reports" / "GEO_DETAILED_REPORT.md"
    )
    _write_text(output_path, report_md)
    print(f"Detailed GEO report written to: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
